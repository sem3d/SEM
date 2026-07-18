!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file irons_dtcrit.F90
!!\brief Rigorous critical time step bound (Irons & Treharne, 1971) via
!!  matrix-free power iteration on the per-element eigenproblem
!!  K_e*phi = omega^2*M_e*phi. Purely diagnostic: does not change the
!!  simulation's actual dt, only prints alongside it.
!<
module m_irons_dtcrit
    use constants
    implicit none

    integer, parameter :: IRONS_MAX_ITER = 2000
    real(fpp), parameter :: IRONS_RTOL = 1e-10_fpp

contains

    !> Entry point: loops over all local (this-rank) CG, non-PML elements,
    !! reduces the global critical dt across ranks, and prints the report.
    !! Called once at startup, after define_arrays (main.F90), so that
    !! Elem%Jacob/Density/invKappa2d/Acoeff/AcoeffFl are already built.
    subroutine compute_irons_dtcrit(Tdomain)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: n, mat, ierr, rg
        real(fpp) :: dt_loc, dt_elem, dt_min, dt_used, ratio, dt_courant_irons
        integer :: n_skipped, n_not_converged, n_total, n_skipped_g, n_not_converged_g, n_total_g
        logical :: converged

        rg = Tdomain%Mpi_Var%my_rank
        n_total = Tdomain%n_elem
        n_skipped = 0
        n_not_converged = 0
        dt_loc = huge(1._fpp)

        do n = 0, Tdomain%n_elem - 1
            if (Tdomain%specel(n)%PML .or. Tdomain%specel(n)%type_DG /= GALERKIN_CONT) then
                n_skipped = n_skipped + 1
                cycle
            end if
            mat = Tdomain%specel(n)%mat_index
            call irons_element_dt(Tdomain, n, mat, dt_elem, converged)
            if (.not. converged) n_not_converged = n_not_converged + 1
            dt_loc = min(dt_loc, dt_elem)
        end do

        call MPI_AllReduce(dt_loc, dt_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)
        ! n_skipped/n_not_converged/n_total above are this rank's local elements
        ! only -- reduce (sum) across ranks so the printed diagnostics describe
        ! the whole mesh, not just rank 0's local partition.
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

    !> Per-element critical dt: dispatches to the solid or acoustic
    !! matrix-free power iteration depending on Elem%acoustic.
    subroutine irons_element_dt(Tdomain, n, mat, dt_elem, converged)
        use sdomain
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n, mat
        real(fpp), intent(out) :: dt_elem
        logical, intent(out) :: converged
        !
        real(fpp) :: omega

        if (Tdomain%specel(n)%acoustic) then
            call irons_power_iteration_acoustic(Tdomain, n, mat, omega, converged)
        else
            call irons_power_iteration_solid(Tdomain, n, mat, omega, converged)
        end if
        dt_elem = 2._fpp / omega
    end subroutine irons_element_dt

    !> Solid (2-component displacement) matrix-free power iteration.
    !! Reuses compute_InternalForces_Elem as the K_e*w operator; the
    !! element's local, UNASSEMBLED mass is rebuilt from Whei*Density*Jacob
    !! (Elem%MassMat is assembled+inverted+interior-only by this point in the
    !! init sequence, so it cannot be reused for the local eigenproblem).
    subroutine irons_power_iteration_solid(Tdomain, n, mat, omega, converged)
        use sdomain
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n, mat
        real(fpp), intent(out) :: omega
        logical, intent(out) :: converged
        !
        integer :: ngllx, ngllz, i, j, iter
        real(fpp), dimension(0:Tdomain%specel(n)%ngllx-1,0:Tdomain%specel(n)%ngllz-1) :: minv_sqrt
        real(fpp), dimension(0:Tdomain%specel(n)%ngllx-1,0:Tdomain%specel(n)%ngllz-1,0:1) :: u, y, saved_forces
        real(fpp) :: raw_lambda, best_lambda, lam1, lam2, aitken, aitken_prev, denom, norm_y

        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz

        do j = 0, ngllz-1
            do i = 0, ngllx-1
                minv_sqrt(i,j) = 1._fpp / sqrt( &
                    Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwz(j) * &
                    Tdomain%specel(n)%Density(i,j) * Tdomain%specel(n)%Jacob(i,j))
                ! deterministic, non-symmetric start vector: cheap and avoids
                ! any dependence on / interference with the global RNG state
                ! used elsewhere (e.g. randomField material generation).
                u(i,j,0) = sin(1.3_fpp*i + 0.7_fpp*j + 1.0_fpp)
                u(i,j,1) = cos(0.9_fpp*i + 1.1_fpp*j + 0.5_fpp)
            end do
        end do
        u = u / sqrt(sum(u*u))

        saved_forces = Tdomain%specel(n)%Forces(:,:,0:1)

        lam1 = 0._fpp; lam2 = 0._fpp; aitken_prev = 0._fpp; best_lambda = 0._fpp
        converged = .false.
        do iter = 1, IRONS_MAX_ITER
            Tdomain%specel(n)%Forces(:,:,0) = minv_sqrt * u(:,:,0)
            Tdomain%specel(n)%Forces(:,:,1) = minv_sqrt * u(:,:,1)
            call compute_InternalForces_Elem(Tdomain%specel(n), &
                Tdomain%sSubDomain(mat)%hprimex, Tdomain%sSubDomain(mat)%hTprimex, &
                Tdomain%sSubDomain(mat)%hprimez, Tdomain%sSubDomain(mat)%hTprimez)
            ! compute_InternalForces_Elem returns -K_e*w (Correction_Elem_Veloc
            ! later does Veloc = V0 + dt*MassMatInv*Forces, i.e. Forces already
            ! carries the F_ext - K*u convention) -- negate to recover +K_e*w.
            y(:,:,0) = -minv_sqrt * Tdomain%specel(n)%Forces(:,:,0)
            y(:,:,1) = -minv_sqrt * Tdomain%specel(n)%Forces(:,:,1)

            raw_lambda = sum(u*y)
            norm_y = sqrt(sum(y*y))
            if (norm_y < tiny(1._fpp)) then
                best_lambda = raw_lambda
                exit
            end if
            u = y / norm_y

            ! Aitken delta-squared extrapolation of the (linearly convergent)
            ! raw Rayleigh-quotient sequence -- the extrapolated value, not
            ! the raw one, is what has to stabilize to IRONS_RTOL. lam1/lam2
            ! always hold RAW lambdas so the recurrence stays well-defined.
            if (iter >= 3) then
                denom = raw_lambda - 2._fpp*lam1 + lam2
                if (abs(denom) > IRONS_RTOL*max(abs(raw_lambda),tiny(1._fpp))) then
                    aitken = raw_lambda - (raw_lambda-lam1)**2 / denom
                else
                    aitken = raw_lambda ! denominator degenerate: already converged to roundoff
                end if
                best_lambda = aitken
                if (iter >= 4 .and. abs(aitken-aitken_prev) < IRONS_RTOL*max(abs(aitken),tiny(1._fpp))) then
                    converged = .true.
                    exit
                end if
                aitken_prev = aitken
            else
                best_lambda = raw_lambda
            end if
            lam2 = lam1
            lam1 = raw_lambda
        end do

        Tdomain%specel(n)%Forces(:,:,0:1) = saved_forces
        omega = sqrt(max(best_lambda,tiny(1._fpp)))
    end subroutine irons_power_iteration_solid

    !> Acoustic (scalar velocity-potential) matrix-free power iteration.
    !! Reuses compute_InternalForcesFl_Elem; local mass = Whei*invKappa2d*Jacob.
    subroutine irons_power_iteration_acoustic(Tdomain, n, mat, omega, converged)
        use sdomain
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: n, mat
        real(fpp), intent(out) :: omega
        logical, intent(out) :: converged
        !
        integer :: ngllx, ngllz, i, j, iter
        real(fpp), dimension(0:Tdomain%specel(n)%ngllx-1,0:Tdomain%specel(n)%ngllz-1) :: minv_sqrt, u, y
        real(fpp), dimension(0:Tdomain%specel(n)%ngllx-1,0:Tdomain%specel(n)%ngllz-1,0:1) :: saved_forces
        real(fpp) :: raw_lambda, best_lambda, lam1, lam2, aitken, aitken_prev, denom, norm_y

        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz

        do j = 0, ngllz-1
            do i = 0, ngllx-1
                minv_sqrt(i,j) = 1._fpp / sqrt( &
                    Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwz(j) * &
                    Tdomain%specel(n)%invKappa2d(i,j) * Tdomain%specel(n)%Jacob(i,j))
                u(i,j) = sin(1.3_fpp*i + 0.7_fpp*j + 1.0_fpp)
            end do
        end do
        u = u / sqrt(sum(u*u))

        saved_forces = Tdomain%specel(n)%Forces(:,:,0:1)

        lam1 = 0._fpp; lam2 = 0._fpp; aitken_prev = 0._fpp; best_lambda = 0._fpp
        converged = .false.
        do iter = 1, IRONS_MAX_ITER
            Tdomain%specel(n)%Forces(:,:,0) = minv_sqrt * u
            call compute_InternalForcesFl_Elem(Tdomain%specel(n), &
                Tdomain%sSubDomain(mat)%hprimex, Tdomain%sSubDomain(mat)%hTprimex, &
                Tdomain%sSubDomain(mat)%hprimez, Tdomain%sSubDomain(mat)%hTprimez)
            ! same -K_e*w convention as compute_InternalForces_Elem -- negate.
            y = -minv_sqrt * Tdomain%specel(n)%Forces(:,:,0)

            raw_lambda = sum(u*y)
            norm_y = sqrt(sum(y*y))
            if (norm_y < tiny(1._fpp)) then
                best_lambda = raw_lambda
                exit
            end if
            u = y / norm_y

            ! Aitken delta-squared extrapolation -- see irons_power_iteration_solid.
            if (iter >= 3) then
                denom = raw_lambda - 2._fpp*lam1 + lam2
                if (abs(denom) > IRONS_RTOL*max(abs(raw_lambda),tiny(1._fpp))) then
                    aitken = raw_lambda - (raw_lambda-lam1)**2 / denom
                else
                    aitken = raw_lambda
                end if
                best_lambda = aitken
                if (iter >= 4 .and. abs(aitken-aitken_prev) < IRONS_RTOL*max(abs(aitken),tiny(1._fpp))) then
                    converged = .true.
                    exit
                end if
                aitken_prev = aitken
            else
                best_lambda = raw_lambda
            end if
            lam2 = lam1
            lam1 = raw_lambda
        end do

        Tdomain%specel(n)%Forces(:,:,0:1) = saved_forces
        omega = sqrt(max(best_lambda,tiny(1._fpp)))
    end subroutine irons_power_iteration_acoustic

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
