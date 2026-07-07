!! ============================================================================
!! solid <-> fluid-aniso coupling for SEM2D.
!!
!! STATUS: EXPERIMENTAL -- vertex DOFs at S-F corners are not coupled (small
!! but nonzero contribution). S-F faces split across MPI ranks are not handled
!! (ensure mesh partitioning places the full interface on one rank). Verified
!! for single-rank CG non-PML runs only.
!!
!! DOF separation is done at the face level: per-side (_sol/_flu) split arrays on
!! every is_sf_iface face carry the solid displacement and the fluid potential
!! separately, so the CG assembly never SUMS Ux(solid)+phi(fluid) at the interface.
!!
!! PHYSICS (mirror of SEM3D solid_fluid_coupling.f90, anisotropy-agnostic at the
!! interface -- rho^-1 lives only in the volume kernel):
!!   StoF: Forces_flu(idxF) += sum_j BtN(j) * Veloc_sol(idxS,j)
!!   FtoS: Forces_sol(idxS,j) -= BtN(j) * VelPhi_flu(idxF)        [VelPhi = -pressure]
!! BtN = outward unit normal (solid->fluid) * line-Jacobian * GLL weight on the face.
!! ============================================================================
module solid_fluid_coupling_2d

    use constants
    use sdomain
    use selement
    implicit none

    type sf_iface_2d
        integer :: nface = -1            ! the changing_media face
        integer :: e_sol = -1, e_flu = -1   ! solid-side, fluid-side element
        integer :: wf_sol = -1, wf_flu = -1 ! which local face on each element
        integer :: ngll = 0
        ! local element (i,j) indices of the face GLL points, per side:
        integer, allocatable :: is_sol(:), js_sol(:)
        integer, allocatable :: is_flu(:), js_flu(:)
        real(fpp), allocatable :: BtN(:,:) ! (0:1,0:ngll-1) outward normal * lineJac * weight
    end type sf_iface_2d

    type(sf_iface_2d), allocatable :: sfi(:)
    integer :: n_sfi = 0

contains

    !-----------------------------------------------------------------------
    !> Enumerate solid<->fluid-aniso interface faces and build BtN + the per-side
    !! face-GLL -> element-(i,j) maps. Call once at init (after define_arrays).
    subroutine build_sf_interface_2d(Tdomain)
        use shape_lin, only: compute_Jacobian_1D
        type(domain), intent(inout) :: Tdomain
        integer :: nf, e0, e1, esol, eflu, nfound, k
        integer :: wf_sol, wf_flu, mat_sol
        integer :: ngll, ngllx_sol, ngllz_sol, ngllx_flu, ngllz_flu
        integer :: p
        logical :: a0, a1, logic_sol, logic_flu
        real(fpp) :: Jac1D, sign_n

        ! ---- count S-F interface faces ----
        nfound = 0
        do nf = 0, Tdomain%n_face-1
            if (.not. Tdomain%sFace(nf)%changing_media) cycle
            e0 = Tdomain%sFace(nf)%Near_Element(0)
            e1 = Tdomain%sFace(nf)%Near_Element(1)
            if (e0 < 0 .or. e1 < 0) cycle
            a0 = allocated(Tdomain%specel(e0)%AcoeffFl)
            a1 = allocated(Tdomain%specel(e1)%AcoeffFl)
            if (a0 .eqv. a1) cycle
            nfound = nfound + 1
        end do
        n_sfi = nfound
        if (n_sfi == 0) return
        allocate(sfi(n_sfi))

        ! ---- fill sf_iface_2d entries ----
        k = 0
        do nf = 0, Tdomain%n_face-1
            if (.not. Tdomain%sFace(nf)%changing_media) cycle
            e0 = Tdomain%sFace(nf)%Near_Element(0)
            e1 = Tdomain%sFace(nf)%Near_Element(1)
            if (e0 < 0 .or. e1 < 0) cycle
            a0 = allocated(Tdomain%specel(e0)%AcoeffFl)
            a1 = allocated(Tdomain%specel(e1)%AcoeffFl)
            if (a0 .eqv. a1) cycle
            k = k + 1
            sfi(k)%nface = nf
            if (a0) then
                sfi(k)%e_flu = e0;  sfi(k)%e_sol = e1
            else
                sfi(k)%e_flu = e1;  sfi(k)%e_sol = e0
            end if
            ngll = Tdomain%sFace(nf)%ngll
            sfi(k)%ngll = ngll

            ! ---- wf_sol / wf_flu ----
            ! Which_face(0) is the face index from Near_Element(0)'s perspective.
            ! Match e_sol/e_flu to Near_Element(0/1).
            if (sfi(k)%e_sol == e0) then
                wf_sol = Tdomain%sFace(nf)%Which_face(0)
                wf_flu = Tdomain%sFace(nf)%Which_face(1)
            else
                wf_sol = Tdomain%sFace(nf)%Which_face(1)
                wf_flu = Tdomain%sFace(nf)%Which_face(0)
            end if
            sfi(k)%wf_sol = wf_sol
            sfi(k)%wf_flu = wf_flu

            ! ---- logic flags for (i,j) mapping ----
            ! Near_Element(0) -> logic=.true.; Near_Element(1) -> logic=coherency
            if (sfi(k)%e_sol == e0) then
                logic_sol = .true.
            else
                logic_sol = Tdomain%sFace(nf)%coherency
            end if
            if (sfi(k)%e_flu == e0) then
                logic_flu = .true.
            else
                logic_flu = Tdomain%sFace(nf)%coherency
            end if

            esol = sfi(k)%e_sol;  eflu = sfi(k)%e_flu
            ngllx_sol = Tdomain%specel(esol)%ngllx
            ngllz_sol = Tdomain%specel(esol)%ngllz
            ngllx_flu = Tdomain%specel(eflu)%ngllx
            ngllz_flu = Tdomain%specel(eflu)%ngllz

            allocate(sfi(k)%is_sol(0:ngll-1), sfi(k)%js_sol(0:ngll-1))
            allocate(sfi(k)%is_flu(0:ngll-1), sfi(k)%js_flu(0:ngll-1))
            allocate(sfi(k)%BtN(0:1, 0:ngll-1))
            sfi(k)%BtN = 0._fpp

            ! ---- (i,j) maps for p=0:ngll-1 ----
            ! wf=0->j=0; wf=1->i=ngllx-1; wf=2->j=ngllz-1; wf=3->i=0
            ! logic=.true.  -> p maps directly (no reversal)
            ! logic=.false. -> p is reversed (ngll-1-p)
            do p = 0, ngll-1
                ! Solid
                select case (wf_sol)
                case (0)
                    if (logic_sol) then
                        sfi(k)%is_sol(p) = p;           sfi(k)%js_sol(p) = 0
                    else
                        sfi(k)%is_sol(p) = ngll-1-p;    sfi(k)%js_sol(p) = 0
                    end if
                case (1)
                    if (logic_sol) then
                        sfi(k)%is_sol(p) = ngllx_sol-1; sfi(k)%js_sol(p) = p
                    else
                        sfi(k)%is_sol(p) = ngllx_sol-1; sfi(k)%js_sol(p) = ngll-1-p
                    end if
                case (2)
                    if (logic_sol) then
                        sfi(k)%is_sol(p) = p;            sfi(k)%js_sol(p) = ngllz_sol-1
                    else
                        sfi(k)%is_sol(p) = ngll-1-p;    sfi(k)%js_sol(p) = ngllz_sol-1
                    end if
                case (3)
                    if (logic_sol) then
                        sfi(k)%is_sol(p) = 0;            sfi(k)%js_sol(p) = p
                    else
                        sfi(k)%is_sol(p) = 0;            sfi(k)%js_sol(p) = ngll-1-p
                    end if
                end select
                ! Fluid
                select case (wf_flu)
                case (0)
                    if (logic_flu) then
                        sfi(k)%is_flu(p) = p;            sfi(k)%js_flu(p) = 0
                    else
                        sfi(k)%is_flu(p) = ngll-1-p;    sfi(k)%js_flu(p) = 0
                    end if
                case (1)
                    if (logic_flu) then
                        sfi(k)%is_flu(p) = ngllx_flu-1; sfi(k)%js_flu(p) = p
                    else
                        sfi(k)%is_flu(p) = ngllx_flu-1; sfi(k)%js_flu(p) = ngll-1-p
                    end if
                case (2)
                    if (logic_flu) then
                        sfi(k)%is_flu(p) = p;            sfi(k)%js_flu(p) = ngllz_flu-1
                    else
                        sfi(k)%is_flu(p) = ngll-1-p;    sfi(k)%js_flu(p) = ngllz_flu-1
                    end if
                case (3)
                    if (logic_flu) then
                        sfi(k)%is_flu(p) = 0;            sfi(k)%js_flu(p) = p
                    else
                        sfi(k)%is_flu(p) = 0;            sfi(k)%js_flu(p) = ngll-1-p
                    end if
                end select
            end do

            ! ---- BtN for interior nodes p=1:ngll-2 ----
            ! BtN(c,p) = sign_n * Normal(c) * Jac1D * GLLw1d(p)
            ! where Normal is the unit outward normal from Near_Element(0).
            ! sign_n = +1 if e_sol==Near_Element(0) (normal points solid->fluid), -1 otherwise.
            ! Jac1D = face_length/2 = arc-length Jacobian for [-1,1] reference segment.
            ! GLLw1d: for wf=0,2 (face along x), use GLLwx; for wf=1,3, use GLLwz.
            call compute_Jacobian_1D(Tdomain, nf, Jac1D)
            if (sfi(k)%e_sol == e0) then
                sign_n = 1._fpp
            else
                sign_n = -1._fpp
            end if
            mat_sol = Tdomain%specel(esol)%mat_index
            ! p=0 and p=ngll-1 are vertex DOFs; skip (BtN stays 0 for vertices).
            do p = 1, ngll-2
                if (wf_sol == 0 .or. wf_sol == 2) then
                    ! face varies in i (x) direction
                    sfi(k)%BtN(0,p) = sign_n * Tdomain%sFace(nf)%Normal(0) * Jac1D * &
                                       Tdomain%sSubdomain(mat_sol)%GLLwx(p)
                    sfi(k)%BtN(1,p) = sign_n * Tdomain%sFace(nf)%Normal(1) * Jac1D * &
                                       Tdomain%sSubdomain(mat_sol)%GLLwx(p)
                else
                    ! wf=1 or 3: face varies in j (z) direction
                    sfi(k)%BtN(0,p) = sign_n * Tdomain%sFace(nf)%Normal(0) * Jac1D * &
                                       Tdomain%sSubdomain(mat_sol)%GLLwz(p)
                    sfi(k)%BtN(1,p) = sign_n * Tdomain%sFace(nf)%Normal(1) * Jac1D * &
                                       Tdomain%sSubdomain(mat_sol)%GLLwz(p)
                end if
            end do

            ! ---- Allocate and initialise split face DOF arrays ----
            Tdomain%sFace(nf)%is_sf_iface = .true.
            allocate(Tdomain%sFace(nf)%Veloc_sol(1:ngll-2, 0:1))
            allocate(Tdomain%sFace(nf)%V0_sol   (1:ngll-2, 0:1))
            allocate(Tdomain%sFace(nf)%Forces_sol(1:ngll-2, 0:1))
            allocate(Tdomain%sFace(nf)%MassMat_sol(1:ngll-2))
            allocate(Tdomain%sFace(nf)%Veloc_flu(1:ngll-2, 0:1))
            allocate(Tdomain%sFace(nf)%V0_flu   (1:ngll-2, 0:1))
            allocate(Tdomain%sFace(nf)%Forces_flu(1:ngll-2, 0:1))
            allocate(Tdomain%sFace(nf)%MassMat_flu(1:ngll-2))
            Tdomain%sFace(nf)%Veloc_sol  = 0._fpp
            Tdomain%sFace(nf)%V0_sol     = 0._fpp
            Tdomain%sFace(nf)%Forces_sol = 0._fpp
            Tdomain%sFace(nf)%Veloc_flu  = 0._fpp
            Tdomain%sFace(nf)%V0_flu     = 0._fpp
            Tdomain%sFace(nf)%Forces_flu = 0._fpp

            ! ---- MassMat per side (from element%MassMat at face boundary nodes) ----
            ! Both solid CG and fluid-aniso CG store the inverse mass in element%MassMat
            ! (used by the GALERKIN_CONT corrector, Element.F90). Accumulate from the
            ! solid element -> MassMat_sol, and from the fluid element -> MassMat_flu.
            do p = 1, ngll-2
                Tdomain%sFace(nf)%MassMat_sol(p) = &
                    Tdomain%specel(esol)%MassMat(sfi(k)%is_sol(p), sfi(k)%js_sol(p))
                Tdomain%sFace(nf)%MassMat_flu(p) = &
                    Tdomain%specel(eflu)%MassMat(sfi(k)%is_flu(p), sfi(k)%js_flu(p))
            end do

        end do

        if (Tdomain%Mpi_var%my_rank == 0) &
            write(*,'(a,i0,a)') ' [aniso] solid-fluid coupling: ', n_sfi, ' interface face(s) found.'

    end subroutine build_sf_interface_2d

    !-----------------------------------------------------------------------
    !> Apply the interface coupling each Newmark step (after internal forces, before
    !! the corrector). Mirror of SEM3D StoF + FtoS, on the split face DOFs.
    subroutine apply_sf_coupling_2d(Tdomain)
        type(domain), intent(inout) :: Tdomain
        integer :: m, p, nf
        real(fpp) :: vn, pf

        do m = 1, n_sfi
            nf = sfi(m)%nface
            do p = 1, sfi(m)%ngll-2    ! interior face nodes only; vertex DOFs skipped
                ! StoF: solid normal velocity -> fluid phi RHS
                vn = sfi(m)%BtN(0,p) * Tdomain%sFace(nf)%Veloc_sol(p,0) &
                   + sfi(m)%BtN(1,p) * Tdomain%sFace(nf)%Veloc_sol(p,1)
                Tdomain%sFace(nf)%Forces_flu(p,0) = Tdomain%sFace(nf)%Forces_flu(p,0) + vn

                ! FtoS: fluid pressure = -VelPhi -> solid force (outward normal convention)
                pf = Tdomain%sFace(nf)%Veloc_flu(p,0)   ! VelPhi = -pressure
                Tdomain%sFace(nf)%Forces_sol(p,0) = Tdomain%sFace(nf)%Forces_sol(p,0) &
                                                   - sfi(m)%BtN(0,p) * pf
                Tdomain%sFace(nf)%Forces_sol(p,1) = Tdomain%sFace(nf)%Forces_sol(p,1) &
                                                   - sfi(m)%BtN(1,p) * pf
            end do
        end do

    end subroutine apply_sf_coupling_2d

end module solid_fluid_coupling_2d
