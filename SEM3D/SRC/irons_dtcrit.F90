!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file irons_dtcrit.F90
!!\brief Rigorous critical time step bound (Irons & Treharne, 1971) via
!!  matrix-free power iteration on the per-element eigenproblem
!!  K_e*phi = omega^2*M_e*phi. Purely diagnostic: does not change the
!!  simulation's actual dt, only prints alongside it. Processes one
!!  VCHUNK-sized block at a time (SEM3D's native vectorization unit -- the
!!  Fox/Foy/Foz force kernels only operate on a whole block, never on a
!!  single element), so memory stays O(1 block) regardless of mesh size.
!!  Only non-PML, non-DG, non-attenuated-kernel CG elements (solid iso/aniso,
!!  fluid iso/aniso) are covered; PML/DG elements are skipped and counted.
!<
#include "index.h"
module m_irons_dtcrit
    use constants
    implicit none

    integer, parameter :: IRONS_MAX_ITER = 2000
    real(fpp), parameter :: IRONS_RTOL = 1e-10_fpp

contains

    !> Entry point, called from Compute_Courant (courant.f90) once at startup,
    !! after init_materials has filled Density_/Lambda_/Mu_/Jacob_/InvGrad_.
    subroutine compute_irons_dtcrit(Tdomain, rg)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer :: ierr
        real(fpp) :: dt_loc, dt_min, dt_used, ratio, dt_courant_irons
        integer :: n_not_converged, n_skipped, n_total, n_covered
        integer :: n_not_converged_g, n_skipped_g, n_total_g

        dt_loc = huge(1._fpp)
        n_not_converged = 0
        n_covered = 0

        if (Tdomain%sdom%nbelem > 0) then
            call irons_domain_solid(Tdomain%sdom, dt_loc, n_not_converged)
            n_covered = n_covered + Tdomain%sdom%nbelem
        end if
        if (Tdomain%fdom%nbelem > 0) then
            call irons_domain_fluid(Tdomain%fdom, dt_loc, n_not_converged)
            n_covered = n_covered + Tdomain%fdom%nbelem
        end if

        n_total = Tdomain%n_elem
        n_skipped = n_total - n_covered

        call MPI_AllReduce(dt_loc, dt_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)
        ! n_skipped/n_not_converged/n_total above are this rank's local
        ! elements only -- reduce (sum) across ranks so the printed
        ! diagnostics describe the whole mesh, not just rank 0's partition.
        call MPI_AllReduce(n_skipped, n_skipped_g, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur, ierr)
        call MPI_AllReduce(n_not_converged, n_not_converged_g, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur, ierr)
        call MPI_AllReduce(n_total, n_total_g, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur, ierr)

        dt_used = Tdomain%TimeD%dtmin
        dt_courant_irons = Tdomain%TimeD%courant * dt_min
        if (dt_courant_irons > 0) then
            ratio = dt_used / dt_courant_irons
        else
            ratio = -1._fpp
        end if

        if (rg == 0) then
            write (*,*) "[Irons dt_crit] eigen-based dt (theoretical bound 2/omega_max) = ", dt_min
            write (*,*) "[Irons dt_crit] eigen-based dt with courant safety factor      = ", dt_courant_irons
            write (*,*) "[Irons dt_crit] steps to cover Duration at dt_Irons: ", &
                int(Tdomain%TimeD%duration/dt_min), " (current heuristic: ", Tdomain%TimeD%ntimeMax, ")"
            write (*,*) "[Irons dt_crit] ratio dt_heuristico_usado / dt_Irons_com_courant = ", ratio
            write (*,*) "[Irons dt_crit] elements skipped (PML/DG): ", n_skipped_g, " / ", n_total_g
            write (*,*) "[Irons dt_crit] elements not converged in", IRONS_MAX_ITER, "iters: ", n_not_converged_g
        end if
    end subroutine compute_irons_dtcrit

    !> Solid CG domain (iso/aniso, unrelaxed elastic operator -- attenuation
    !! ignored on purpose: it only makes the medium look stiffer, i.e. dt
    !! smaller, the safe direction). Reuses calcul_forces_iso/calcul_forces_aniso
    !! (SEM3D/SRC/Solid/calcul_forces.F90) unmodified as the K_e*w operator.
    subroutine irons_domain_solid(dom, dt_loc, n_not_converged)
        use champs_solid
        use m_calcul_forces
        implicit none
        type(domain_solid), intent(inout) :: dom
        real(fpp), intent(inout) :: dt_loc
        integer, intent(inout) :: n_not_converged
        !
        integer :: bnum, ee, i, j, k, c, iter, ngll, n_valid
        real(fpp), allocatable :: minv_sqrt(:,:,:,:)
        real(fpp), allocatable :: u(:,:,:,:,:), y(:,:,:,:,:), Depla(:,:,:,:,:)
        real(fpp), allocatable :: Fox(:,:,:,:), Foy(:,:,:,:), Foz(:,:,:,:)
        real(fpp), dimension(0:VCHUNK-1) :: raw_lambda, best_lambda, lam1, lam2, aitken_prev, norm_y
        logical, dimension(0:VCHUNK-1) :: done
        real(fpp) :: omega, aitken, denom

        ngll = dom%ngll
        allocate(minv_sqrt(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1))
        allocate(u(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        allocate(y, mold=u)
        allocate(Depla, mold=u)
        allocate(Fox(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1))
        allocate(Foy, mold=Fox)
        allocate(Foz, mold=Fox)

        do bnum = 0, dom%nblocks - 1
            n_valid = min(VCHUNK, dom%nbelem - bnum*VCHUNK)

            do ee = 0, VCHUNK-1
                if (ee >= n_valid) then
                    minv_sqrt(ee,:,:,:) = 1._fpp
                    u(ee,:,:,:,:) = 0._fpp
                    cycle
                end if
                do k = 0, ngll-1
                    do j = 0, ngll-1
                        do i = 0, ngll-1
                            minv_sqrt(ee,i,j,k) = 1._fpp / sqrt(dom%gllw(i)*dom%gllw(j)*dom%gllw(k)* &
                                dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee))
                            ! deterministic, non-symmetric start vector: cheap, avoids
                            ! touching the global RNG state (used elsewhere for randomField).
                            u(ee,i,j,k,0) = sin(1.3_fpp*i + 0.7_fpp*j + 1.9_fpp*k + 1.0_fpp)
                            u(ee,i,j,k,1) = cos(0.9_fpp*i + 1.1_fpp*j + 0.6_fpp*k + 0.5_fpp)
                            u(ee,i,j,k,2) = sin(0.5_fpp*i + 1.7_fpp*j + 1.1_fpp*k + 2.2_fpp)
                        end do
                    end do
                end do
                u(ee,:,:,:,:) = u(ee,:,:,:,:) / sqrt(sum(u(ee,:,:,:,:)**2))
            end do

            lam1 = 0._fpp; lam2 = 0._fpp; aitken_prev = 0._fpp; best_lambda = 0._fpp
            done = .false.
            do iter = 1, IRONS_MAX_ITER
                do c = 0, 2
                    Depla(:,:,:,:,c) = minv_sqrt * u(:,:,:,:,c)
                end do
                if (dom%aniso) then
                    call calcul_forces_aniso(dom, bnum, Fox, Foy, Foz, Depla)
                else
                    call calcul_forces_iso(dom, bnum, Fox, Foy, Foz, Depla)
                end if
                y(:,:,:,:,0) = minv_sqrt * Fox
                y(:,:,:,:,1) = minv_sqrt * Foy
                y(:,:,:,:,2) = minv_sqrt * Foz

                do ee = 0, n_valid-1
                    if (done(ee)) cycle
                    raw_lambda(ee) = sum(u(ee,:,:,:,:) * y(ee,:,:,:,:))
                    norm_y(ee) = sqrt(sum(y(ee,:,:,:,:)**2))
                    if (norm_y(ee) > tiny(1._fpp)) u(ee,:,:,:,:) = y(ee,:,:,:,:) / norm_y(ee)

                    ! Aitken delta-squared extrapolation of the raw Rayleigh-quotient
                    ! sequence -- lam1/lam2 always hold RAW lambdas (never the
                    ! extrapolated value) so the recurrence stays well-defined.
                    if (iter >= 3) then
                        denom = raw_lambda(ee) - 2._fpp*lam1(ee) + lam2(ee)
                        if (abs(denom) > IRONS_RTOL*max(abs(raw_lambda(ee)),tiny(1._fpp))) then
                            aitken = raw_lambda(ee) - (raw_lambda(ee)-lam1(ee))**2 / denom
                        else
                            aitken = raw_lambda(ee) ! denominator degenerate: converged to roundoff
                        end if
                        best_lambda(ee) = aitken
                        if (iter >= 4 .and. abs(aitken-aitken_prev(ee)) < IRONS_RTOL*max(abs(aitken),tiny(1._fpp))) &
                            done(ee) = .true.
                        aitken_prev(ee) = aitken
                    else
                        best_lambda(ee) = raw_lambda(ee)
                    end if
                    lam2(ee) = lam1(ee)
                    lam1(ee) = raw_lambda(ee)
                end do

                if (all(done(0:n_valid-1))) exit
            end do

            do ee = 0, n_valid-1
                if (.not. done(ee)) n_not_converged = n_not_converged + 1
                omega = sqrt(max(best_lambda(ee), tiny(1._fpp)))
                dt_loc = min(dt_loc, 2._fpp/omega)
            end do
        end do
    end subroutine irons_domain_solid

    !> Fluid CG domain (velocity-potential formulation, iso/aniso density).
    !! Reuses calcul_forces_fluid (dispatches iso/aniso internally) as the
    !! K_e*w operator; local mass = Whei*Jacob_ divided by the fluid bulk
    !! modulus, stored in the Lambda_ field (same field name reused across
    !! domain types).
    subroutine irons_domain_fluid(dom, dt_loc, n_not_converged)
        use champs_fluid
        use m_calcul_forces_fluid
        implicit none
        type(domain_fluid), intent(inout) :: dom
        real(fpp), intent(inout) :: dt_loc
        integer, intent(inout) :: n_not_converged
        !
        integer :: bnum, ee, i, j, k, iter, ngll, n_valid
        real(fpp), allocatable :: minv_sqrt(:,:,:,:), u(:,:,:,:), y(:,:,:,:), Phi(:,:,:,:), FFl(:,:,:,:)
        real(fpp), dimension(0:VCHUNK-1) :: raw_lambda, best_lambda, lam1, lam2, aitken_prev, norm_y
        logical, dimension(0:VCHUNK-1) :: done
        real(fpp) :: omega, aitken, denom

        ngll = dom%ngll
        allocate(minv_sqrt(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1))
        allocate(u, mold=minv_sqrt)
        allocate(y, mold=minv_sqrt)
        allocate(Phi, mold=minv_sqrt)
        allocate(FFl, mold=minv_sqrt)

        do bnum = 0, dom%nblocks - 1
            n_valid = min(VCHUNK, dom%nbelem - bnum*VCHUNK)

            do ee = 0, VCHUNK-1
                if (ee >= n_valid) then
                    minv_sqrt(ee,:,:,:) = 1._fpp
                    u(ee,:,:,:) = 0._fpp
                    cycle
                end if
                do k = 0, ngll-1
                    do j = 0, ngll-1
                        do i = 0, ngll-1
                            minv_sqrt(ee,i,j,k) = 1._fpp / sqrt(dom%gllw(i)*dom%gllw(j)*dom%gllw(k)* &
                                dom%Jacob_(i,j,k,bnum,ee) / dom%Lambda_(i,j,k,bnum,ee))
                            u(ee,i,j,k) = sin(1.3_fpp*i + 0.7_fpp*j + 1.9_fpp*k + 1.0_fpp)
                        end do
                    end do
                end do
                u(ee,:,:,:) = u(ee,:,:,:) / sqrt(sum(u(ee,:,:,:)**2))
            end do

            lam1 = 0._fpp; lam2 = 0._fpp; aitken_prev = 0._fpp; best_lambda = 0._fpp
            done = .false.
            do iter = 1, IRONS_MAX_ITER
                Phi = minv_sqrt * u
                call calcul_forces_fluid(dom, ngll, bnum, FFl, Phi)
                y = minv_sqrt * FFl

                do ee = 0, n_valid-1
                    if (done(ee)) cycle
                    raw_lambda(ee) = sum(u(ee,:,:,:) * y(ee,:,:,:))
                    norm_y(ee) = sqrt(sum(y(ee,:,:,:)**2))
                    if (norm_y(ee) > tiny(1._fpp)) u(ee,:,:,:) = y(ee,:,:,:) / norm_y(ee)

                    ! Aitken delta-squared extrapolation -- see irons_domain_solid.
                    if (iter >= 3) then
                        denom = raw_lambda(ee) - 2._fpp*lam1(ee) + lam2(ee)
                        if (abs(denom) > IRONS_RTOL*max(abs(raw_lambda(ee)),tiny(1._fpp))) then
                            aitken = raw_lambda(ee) - (raw_lambda(ee)-lam1(ee))**2 / denom
                        else
                            aitken = raw_lambda(ee)
                        end if
                        best_lambda(ee) = aitken
                        if (iter >= 4 .and. abs(aitken-aitken_prev(ee)) < IRONS_RTOL*max(abs(aitken),tiny(1._fpp))) &
                            done(ee) = .true.
                        aitken_prev(ee) = aitken
                    else
                        best_lambda(ee) = raw_lambda(ee)
                    end if
                    lam2(ee) = lam1(ee)
                    lam1(ee) = raw_lambda(ee)
                end do

                if (all(done(0:n_valid-1))) exit
            end do

            do ee = 0, n_valid-1
                if (.not. done(ee)) n_not_converged = n_not_converged + 1
                omega = sqrt(max(best_lambda(ee), tiny(1._fpp)))
                dt_loc = min(dt_loc, 2._fpp/omega)
            end do
        end do
    end subroutine irons_domain_fluid

end module m_irons_dtcrit

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
