!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file NewmarkModified.f90
!!\brief Matrix-free "modified equation" Newmark, regional variant for
!! SEM3D (see docs/superpowers/specs/2026-07-28-newmark-modified-sem3d-design.md
!! and SEM2D/SRC/NewmarkModified.F90 for the reference derivation). Only
!! sdom/fdom (non-PML, non-DG) DOFs flagged %modified (core+buffer, see
!! SEM3D/SRC/irons_dtcrit.F90's select_modified_region) receive the order-m
!! correction; spmldom/fpmldom/sdomdg elements and solid-fluid-interface
!! DOFs keep running through the classic Newmark call below, untouched.
!! Mono-rank only (enforced at startup by irons_dtcrit.F90).
!!
!! depl(n+1) = 2*depl(n) - depl(n-1) + dt^2*Minv*(Fext-K*depl(n))
!!           + sum_{k=2}^m c_k * A^k * depl(n),   A = -Minv*K,
!!           c_k = 2*dt^(2k)/(2k)!
!! The first line is exactly what the base (classic, order-1) Newmark call
!! already computes. This routine: (1) snapshots depl(n) for every sdom/fdom
!! DOF, (2) calls the UNMODIFIED classic Newmark for the base step, (3) adds
!! the k=2..m correction terms via k-1 extra matrix-free applications of A
!! (reusing calcul_forces_iso/aniso/fluid, restricted each iteration to the
!! shrinking block_active_*(k) set -- see irons_dtcrit.F90's module-level
!! comment on elem_hop for why the halo must shrink with k), and (4) resyncs
!! Veloc = (Depla_new - depl(n))/dt at every corrected DOF so the next
!! step's classic Newmark call sees a consistent (Depla,Veloc) pair again.
!<
#include "index.h"
module snewmark_modified
    use constants
    use sdomain
    use mtimestep, only : Newmark
    use m_irons_dtcrit, only : block_active_sdom, block_active_fdom, &
        dof_touched_sdom, dof_touched_fdom, dof_corrected_sdom, dof_corrected_fdom, fact2k
    implicit none

    ! un = snapshot of depl(n) taken before the base step (ALL DOFs, not just
    ! the corrected set -- the shrinking-halo gather at early iterations
    ! needs a valid seed up to hop<=modified_order-1, wider than the final
    ! corrected set at hop=0). w/wout = current/next matrix-free iterate
    ! (A^(k-1)*depl(n) -> A^k*depl(n)); corr = accumulated sum_{k=2}^m c_k*A^k*depl(n).
    real(fpp), dimension(:,:), allocatable, save :: sdom_un, sdom_w, sdom_wout, sdom_corr
    real(fpp), dimension(:), allocatable, save :: fdom_un, fdom_w, fdom_wout, fdom_corr
    logical, save :: scratch_ready = .false.

contains

    !> Allocate the module-level scratch on first use (after allocate_domain
    !! has sized champs(0)%Depla/Veloc/Phi/VelPhi, so the shapes match).
    subroutine ensure_scratch(Tdomain)
        implicit none
        type(domain), intent(in) :: Tdomain
        if (scratch_ready) return
        allocate(sdom_un(0:Tdomain%sdom%nglltot-1, 0:2))
        allocate(sdom_w(0:Tdomain%sdom%nglltot-1, 0:2))
        allocate(sdom_wout(0:Tdomain%sdom%nglltot-1, 0:2))
        allocate(sdom_corr(0:Tdomain%sdom%nglltot-1, 0:2))
        allocate(fdom_un(0:Tdomain%fdom%nglltot-1))
        allocate(fdom_w(0:Tdomain%fdom%nglltot-1))
        allocate(fdom_wout(0:Tdomain%fdom%nglltot-1))
        allocate(fdom_corr(0:Tdomain%fdom%nglltot-1))
        scratch_ready = .true.
    end subroutine ensure_scratch

    subroutine NewmarkModified(Tdomain, ntime)
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: ntime
        integer :: k, mm, idx
        real(fpp) :: dt, ck

        if (.not. Tdomain%TimeD%velocity_scheme) return
        call ensure_scratch(Tdomain)
        mm = Tdomain%TimeD%modified_order
        dt = Tdomain%TimeD%dtmin

        ! 1) Snapshot depl(n) for every sdom/fdom DOF before the base step
        ! overwrites Depla/Phi.
        if (Tdomain%sdom%nglltot > 0) sdom_un = Tdomain%sdom%champs(0)%Depla
        if (Tdomain%fdom%nglltot > 0) fdom_un = Tdomain%fdom%champs(0)%Phi

        ! 2) Base order-1 step: reuse the existing, unmodified classic
        ! Newmark -- this is exactly the k=1 term of the recursion.
        call Newmark(Tdomain, ntime)

        if (mm < 2) return

        ! 3) k=2..m correction: seed w := depl(n), chain A*w -> A^2*w -> ...
        if (Tdomain%sdom%nglltot > 0) then
            sdom_w = sdom_un
            sdom_corr = 0._fpp
        end if
        if (Tdomain%fdom%nglltot > 0) then
            fdom_w = fdom_un
            fdom_corr = 0._fpp
        end if

        do k = 1, mm
            call apply_a_regional(Tdomain, k)

            if (k >= 2) then
                ck = 2._fpp*dt**(2*k)/fact2k(k)
                do idx = 0, Tdomain%sdom%nglltot - 1
                    if (dof_corrected_sdom(idx)) sdom_corr(idx,:) = sdom_corr(idx,:) + ck*sdom_wout(idx,:)
                end do
                do idx = 0, Tdomain%fdom%nglltot - 1
                    if (dof_corrected_fdom(idx)) fdom_corr(idx) = fdom_corr(idx) + ck*fdom_wout(idx)
                end do
            end if

            ! Advance w for the next iteration, but ONLY at DOFs actually
            ! recomputed this iteration (dof_touched_*(k)) -- everywhere
            ! else must keep its previous value untouched (required by the
            ! shrinking-halo invariant: iteration k+1 needs w correct at
            ! exactly the region iteration k just computed, nothing wider).
            associate (idxs => dof_touched_sdom(k)%idx)
                do idx = 1, size(idxs)
                    sdom_w(idxs(idx),:) = sdom_wout(idxs(idx),:)
                end do
            end associate
            associate (idxs => dof_touched_fdom(k)%idx)
                do idx = 1, size(idxs)
                    fdom_w(idxs(idx)) = fdom_wout(idxs(idx))
                end do
            end associate
        end do

        ! 4) Apply the correction (corrected DOFs only) and resync
        ! Veloc = (Depla_new - depl(n))/dt there so the next step's
        ! classic-Newmark base call sees a consistent pair again.
        do idx = 0, Tdomain%sdom%nglltot - 1
            if (.not. dof_corrected_sdom(idx)) cycle
            Tdomain%sdom%champs(0)%Depla(idx,:) = Tdomain%sdom%champs(0)%Depla(idx,:) + sdom_corr(idx,:)
            Tdomain%sdom%champs(0)%Veloc(idx,:) = (Tdomain%sdom%champs(0)%Depla(idx,:) - sdom_un(idx,:)) / dt
        end do
        do idx = 0, Tdomain%fdom%nglltot - 1
            if (.not. dof_corrected_fdom(idx)) cycle
            Tdomain%fdom%champs(0)%Phi(idx) = Tdomain%fdom%champs(0)%Phi(idx) + fdom_corr(idx)
            Tdomain%fdom%champs(0)%VelPhi(idx) = (Tdomain%fdom%champs(0)%Phi(idx) - fdom_un(idx)) / dt
        end do
    end subroutine NewmarkModified

    !> A = -Minv*K, matrix-free: given w (module-level sdom_w/fdom_w),
    !! overwrite sdom_wout/fdom_wout with A*w, restricted to
    !! block_active_sdom(jiter)/block_active_fdom(jiter) (whole blocks --
    !! calcul_forces_iso/aniso/fluid only ever operate on a full VCHUNK
    !! block, there is no per-element call in SEM3D). No MPI exchange
    !! (mono-rank only, enforced at startup by irons_dtcrit.F90). No
    !! external forcing (Compute_external_forces is NOT called) -- this is a
    !! pure linear operator application, chainable to get A^2*w, etc.
    subroutine apply_a_regional(Tdomain, jiter)
        use m_calcul_forces
        use m_calcul_forces_fluid
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: jiter
        integer :: ib, bnum, ngll, i, j, k, ee, idx, n
        real(fpp), allocatable :: Fox(:,:,:,:), Foy(:,:,:,:), Foz(:,:,:,:), Depla(:,:,:,:,:)
        real(fpp), allocatable :: FFl(:,:,:,:), Phi(:,:,:,:)

        ! Zero exactly the DOFs recomputed this iteration (everything else
        ! in sdom_wout/fdom_wout keeps last iteration's stale value, which is
        ! never read again -- see the shrinking-halo comment above).
        associate (idxs => dof_touched_sdom(jiter)%idx)
            do idx = 1, size(idxs)
                sdom_wout(idxs(idx),:) = 0._fpp
            end do
        end associate
        associate (idxs => dof_touched_fdom(jiter)%idx)
            do idx = 1, size(idxs)
                fdom_wout(idxs(idx)) = 0._fpp
            end do
        end associate

        ngll = Tdomain%sdom%ngll
        if (ngll > 0 .and. size(block_active_sdom(jiter)%idx) > 0) then
            allocate(Fox(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1))
            allocate(Foy, mold=Fox); allocate(Foz, mold=Fox)
            allocate(Depla(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2))
            do ib = 1, size(block_active_sdom(jiter)%idx)
                bnum = block_active_sdom(jiter)%idx(ib)
                do k = 0, ngll-1
                    do j = 0, ngll-1
                        do i = 0, ngll-1
                            do ee = 0, VCHUNK-1
                                n = Tdomain%sdom%Idom_(i,j,k,bnum,ee)
                                Depla(ee,i,j,k,:) = sdom_w(n,:)
                            end do
                        end do
                    end do
                end do
                if (Tdomain%sdom%aniso) then
                    call calcul_forces_aniso(Tdomain%sdom, bnum, Fox, Foy, Foz, Depla)
                else
                    call calcul_forces_iso(Tdomain%sdom, bnum, Fox, Foy, Foz, Depla)
                end if
                ! Accumulate -K*w, same sign convention as forces_int_solid
                ! (dvdt%Veloc -= Fox/Foy/Foz).
                do k = 0, ngll-1
                    do j = 0, ngll-1
                        do i = 0, ngll-1
                            do ee = 0, VCHUNK-1
                                n = Tdomain%sdom%Idom_(i,j,k,bnum,ee)
                                sdom_wout(n,0) = sdom_wout(n,0) - Fox(ee,i,j,k)
                                sdom_wout(n,1) = sdom_wout(n,1) - Foy(ee,i,j,k)
                                sdom_wout(n,2) = sdom_wout(n,2) - Foz(ee,i,j,k)
                            end do
                        end do
                    end do
                end do
            end do
            deallocate(Fox, Foy, Foz, Depla)
        end if

        ngll = Tdomain%fdom%ngll
        if (ngll > 0 .and. size(block_active_fdom(jiter)%idx) > 0) then
            allocate(FFl(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1))
            allocate(Phi, mold=FFl)
            do ib = 1, size(block_active_fdom(jiter)%idx)
                bnum = block_active_fdom(jiter)%idx(ib)
                do k = 0, ngll-1
                    do j = 0, ngll-1
                        do i = 0, ngll-1
                            do ee = 0, VCHUNK-1
                                n = Tdomain%fdom%Idom_(i,j,k,bnum,ee)
                                Phi(ee,i,j,k) = fdom_w(n)
                            end do
                        end do
                    end do
                end do
                call calcul_forces_fluid(Tdomain%fdom, ngll, bnum, FFl, Phi)
                ! Same sign convention as forces_int_fluid (ForcesFl -= Fo_Fl).
                do k = 0, ngll-1
                    do j = 0, ngll-1
                        do i = 0, ngll-1
                            do ee = 0, VCHUNK-1
                                n = Tdomain%fdom%Idom_(i,j,k,bnum,ee)
                                fdom_wout(n) = fdom_wout(n) - FFl(ee,i,j,k)
                            end do
                        end do
                    end do
                end do
            end do
            deallocate(FFl, Phi)
        end if

        ! A*w = Minv*(accumulated -K*w); MassMat is already 1/M (see
        ! define_arrays.F90). Then zero Dirichlet DOFs, mirroring
        ! newmark_corrector_solid/fluid's own Dirichlet zeroing so absorbing
        ! (Dirichlet) boundaries see the same physics here as in the base
        ! (classic Newmark) step.
        associate (idxs => dof_touched_sdom(jiter)%idx)
            do idx = 1, size(idxs)
                sdom_wout(idxs(idx),:) = sdom_wout(idxs(idx),:) * Tdomain%sdom%MassMat(idxs(idx))
            end do
        end associate
        do idx = 0, Tdomain%sdom%n_dirich - 1
            sdom_wout(Tdomain%sdom%dirich(idx),:) = 0._fpp
        end do

        associate (idxs => dof_touched_fdom(jiter)%idx)
            do idx = 1, size(idxs)
                fdom_wout(idxs(idx)) = fdom_wout(idxs(idx)) * Tdomain%fdom%MassMat(idxs(idx))
            end do
        end associate
        do idx = 0, Tdomain%fdom%n_dirich - 1
            fdom_wout(Tdomain%fdom%dirich(idx)) = 0._fpp
        end do
    end subroutine apply_a_regional

end module snewmark_modified

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
