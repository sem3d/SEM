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
    use m_mod_eq_zcrit
    use m_cost_optimizer_common
    use m_modified_newmark_logger
    implicit none

    integer, parameter :: IRONS_MAX_ITER = 2000
    real(fpp), parameter :: IRONS_RTOL = 1e-10_fpp

    ! Regional modified-equation Newmark: fraction of (local, per-rank)
    ! elements with the smallest critical dt that seed the "core" needing the
    ! order-m correction, and the geometric buffer radius around each core
    ! element (multiple of that element's dist_max) that also gets corrected
    ! for a smooth order-1/order-m transition. ponytail: hardcoded rather than
    ! wired into input.spec -- promote to a config key if these ever need
    ! per-run tuning instead of the values agreed on 2026-07-22.
    real(fpp), parameter :: MODIFIED_FRACTION = 0.05_fpp
    real(fpp), parameter :: MODIFIED_BUFFER_MULT = 3.0_fpp

contains

    !> Entry point: loops over all local (this-rank) CG, non-PML elements,
    !! reduces the global critical dt across ranks, and prints the report.
    !! Called once at startup, after define_arrays (main.F90), so that
    !! Elem%Jacob/Density/invKappa2d/Acoeff/AcoeffFl are already built.
    !!
    !! If Tdomain%TimeD%modified (newmark_modified=true in input.spec), the
    !! order-m bound (see ModifiedEquationZCrit) DRIVES the real dtmin used
    !! by the simulation, overriding the courant.F90 heuristic set earlier
    !! in the startup sequence - order 1 (classic Newmark's own eigenvalue
    !! bound) leaves dtmin untouched otherwise, matching prior behaviour.
    subroutine compute_irons_dtcrit(Tdomain)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: n, mat, ierr, rg, mm
        real(fpp) :: dt_loc, dt_elem, dt_min, dt_used, ratio, dt_courant_irons
        real(fpp) :: zc, dt_target_g, dt_courant_irons_order
        integer :: n_skipped, n_not_converged, n_total, n_skipped_g, n_not_converged_g, n_total_g
        logical :: converged
        real(fpp), dimension(:), allocatable :: dt_elem_loc
        integer :: n_included
        real(fpp) :: dt_sum_loc, dt_max_loc

        rg = Tdomain%Mpi_Var%my_rank
        n_total = Tdomain%n_elem
        n_skipped = 0
        n_not_converged = 0
        dt_loc = huge(1._fpp)
        dt_sum_loc = 0._fpp
        dt_max_loc = 0._fpp
        n_included = 0
        allocate(dt_elem_loc(0:Tdomain%n_elem-1))

        do n = 0, Tdomain%n_elem - 1
            if (Tdomain%specel(n)%PML .or. Tdomain%specel(n)%type_DG /= GALERKIN_CONT) then
                n_skipped = n_skipped + 1
                dt_elem_loc(n) = -1._fpp ! marks "excluded" for report_local_order_requirements
                cycle
            end if
            mat = Tdomain%specel(n)%mat_index
            call irons_element_dt(Tdomain, n, mat, dt_elem, converged)
            if (.not. converged) n_not_converged = n_not_converged + 1
            dt_loc = min(dt_loc, dt_elem)
            dt_sum_loc = dt_sum_loc + dt_elem
            dt_max_loc = max(dt_max_loc, dt_elem)
            n_included = n_included + 1
            dt_elem_loc(n) = dt_elem
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

        ! Diagnostic only (does not change the simulation): if the mesh has
        ! a small minority of elements dragging dt_min down for everyone, this
        ! reports how many elements would need which modified-equation order
        ! to tolerate a larger, more "typical" dt - and an IDEALIZED (best
        ! case, ignoring the neighbourhood-growth overhead a truncated
        ! regional implementation adds) speedup estimate.
        block
            real(fpp) :: dt_sum_g, dt_max_g, dt_mean
            integer :: n_included_g
            call MPI_AllReduce(dt_sum_loc, dt_sum_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, Tdomain%communicateur, ierr)
            call MPI_AllReduce(dt_max_loc, dt_max_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, Tdomain%communicateur, ierr)
            call MPI_AllReduce(n_included, n_included_g, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur, ierr)
            if (n_included_g > 0) then
                dt_mean = dt_sum_g / real(n_included_g,fpp)
                call report_local_order_requirements(Tdomain, dt_elem_loc, "mean", dt_mean, dt_min, 5)
                call report_local_order_requirements(Tdomain, dt_elem_loc, "max (best case)", dt_max_g, dt_min, 5)
            end if
        end block

        if (Tdomain%TimeD%modified) then
            block
                use snewmark_modified, only : build_regional_halo_with_orders
                integer, dimension(:), allocatable :: elem_order_opt
                if (rg == 0) call log_header_modified_newmark("2D")
                allocate(elem_order_opt(0:Tdomain%n_elem-1))
                call optimize_cost_and_orders(Tdomain, dt_elem_loc, dt_target_g, elem_order_opt)
                call build_regional_halo_with_orders(Tdomain, elem_order_opt)
                deallocate(elem_order_opt)
            end block

            zc = get_zcrit(Tdomain%TimeD%modified_order)
            dt_courant_irons_order = Tdomain%TimeD%courant * dt_target_g

            if (rg == 0) then
                write (*,*) "[Irons dt_crit] optimal max_modified_order = ", Tdomain%TimeD%modified_order, "  z_crit = ", zc
                write (*,*) "[Irons dt_crit] optimal target dt (drives the simulation) = ", dt_target_g
                write (*,*) "[Irons dt_crit] optimal target dt with courant safety factor = ", dt_courant_irons_order
            end if

            do mat = 0, Tdomain%n_mat - 1
                Tdomain%sSubdomain(mat)%Dt = dt_courant_irons_order
            enddo
            Tdomain%TimeD%dtmin = dt_courant_irons_order
            if (Tdomain%TimeD%dtmin > 0) then
                Tdomain%TimeD%ntimeMax = int(Tdomain%TimeD%duration/Tdomain%TimeD%dtmin)
            else
                write (*,*) "Your dt min is zero : verify it"
                stop
            endif
        end if
        deallocate(dt_elem_loc)
    end subroutine compute_irons_dtcrit

    !> Diagnostic only: for a candidate dt_target, classify every local
    !! (order-1-basis) element by the minimum modified-equation order m in
    !! [1,max_order] such that z_crit(m) >= 4*(dt_target/dt_elem)^2, i.e.
    !! the smallest order that would let that element tolerate dt_target
    !! (elements with dt_elem >= dt_target already need no help: order 1).
    !! Reports the histogram (reduced across ranks) and an IDEALIZED
    !! speedup estimate cost_baseline/cost_regional = (dt_target/dt_min) *
    !! n_total / sum(m_e) - a best-case bound that ignores the cost of the
    !! neighbourhood growth a truncated regional implementation needs to
    !! pay at the boundary of the "needs boost" region.
    subroutine report_local_order_requirements(Tdomain, dt_elem_loc, label, dt_target, dt_min, max_order)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(in) :: Tdomain
        real(fpp), dimension(0:Tdomain%n_elem-1), intent(in) :: dt_elem_loc
        character(len=*), intent(in) :: label
        real(fpp), intent(in) :: dt_target, dt_min
        integer, intent(in) :: max_order

        integer :: n, m, rg, ierr
        integer, dimension(0:10) :: hist_loc, hist_g ! index 0 = "impossible within max_order"
        integer :: m_needed, sum_m_loc, sum_m_g, n_total_loc, n_total_g
        real(fpp) :: need_ratio, zc

        rg = Tdomain%Mpi_Var%my_rank
        hist_loc = 0
        sum_m_loc = 0
        n_total_loc = 0

        do n = 0, Tdomain%n_elem - 1
            if (dt_elem_loc(n) < 0._fpp) cycle ! PML/DG-skipped, not part of this diagnostic
            n_total_loc = n_total_loc + 1
            if (dt_elem_loc(n) >= dt_target) then
                m_needed = 1
            else
                need_ratio = 4._fpp * (dt_target/dt_elem_loc(n))**2
                m_needed = 0 ! 0 = impossible within max_order
                do m = 1, max_order
                    zc = ModifiedEquationZCrit(m)
                    if (zc >= need_ratio) then
                        m_needed = m
                        exit
                    end if
                end do
            end if
            if (m_needed == 0) then
                hist_loc(0) = hist_loc(0) + 1
                sum_m_loc = sum_m_loc + max_order ! conservative: would need > max_order
            else
                if (m_needed <= 10) hist_loc(m_needed) = hist_loc(m_needed) + 1
                sum_m_loc = sum_m_loc + m_needed
            end if
        end do

        call MPI_AllReduce(hist_loc, hist_g, 11, MPI_INTEGER, MPI_SUM, Tdomain%communicateur, ierr)
        call MPI_AllReduce(sum_m_loc, sum_m_g, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur, ierr)
        call MPI_AllReduce(n_total_loc, n_total_g, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur, ierr)

        if (rg == 0 .and. n_total_g > 0) then
            write (*,*) "[Irons order-diag] target dt = ", trim(label), " = ", dt_target
            write (*,*) "[Irons order-diag]   order 1 (no boost needed): ", hist_g(1), " / ", n_total_g
            do m = 2, min(max_order,10)
                if (hist_g(m) > 0) write (*,*) "[Irons order-diag]   order ", m, " needed: ", hist_g(m)
            end do
            if (hist_g(0) > 0) write (*,*) "[Irons order-diag]   IMPOSSIBLE within order ", max_order, ": ", hist_g(0)
            write (*,*) "[Irons order-diag]   idealized best-case speedup (no truncation overhead) = ", &
                (real(n_total_g,fpp)/real(sum_m_g,fpp)) * (dt_target/dt_min)
        end if
    end subroutine report_local_order_requirements

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

    !> S(z) = sum_{k=1}^m 2*(-z)^k/(2k)!, the modified-equation amplification
    !> Speedup & Cost Optimization: Evaluates candidate target time steps dt_target
    !! (global percentiles of dt_elem_loc), determines minimum required element orders m_i,
    !! resolves halos via Distributed BFS, and selects the (dt_target, {m_i}) configuration
    !! that minimizes total computational cost C = (1 / dt_target) * sum(m_i).
    subroutine optimize_cost_and_orders(Tdomain, dt_elem_loc, dt_target_opt, elem_order_opt)
        use sdomain
        use mpi
        use snewmark_modified, only : resolver_halos_from_orders
        implicit none

        type(domain), intent(inout) :: Tdomain
        real(fpp), dimension(0:Tdomain%n_elem-1), intent(in) :: dt_elem_loc
        real(fpp), intent(out) :: dt_target_opt
        integer, dimension(0:Tdomain%n_elem-1), intent(out) :: elem_order_opt

        integer :: n, i, c, ierr, rg, n_procs, n_local, n_global
        integer :: cand_idx, m_candidate, feasible_loc, feasible_g
        real(fpp) :: r_req, local_cost, global_cost, total_cost, min_cost
        real(fpp) :: dt_target_cand
        logical :: cand_feasible

        integer, dimension(:), allocatable :: recv_counts, displs
        real(fpp), dimension(:), allocatable :: dt_global, dt_candidates
        integer, dimension(:), allocatable :: m_base_local, m_work_local

        integer, parameter :: NUM_PERCENTILES = 10
        real(fpp), dimension(NUM_PERCENTILES), parameter :: PERCENTILES = &
            (/ 0.00_fpp, 0.02_fpp, 0.05_fpp, 0.10_fpp, 0.20_fpp, 0.30_fpp, 0.40_fpp, 0.50_fpp, 0.70_fpp, 0.90_fpp /)

        rg = Tdomain%Mpi_Var%my_rank
        n_procs = Tdomain%Mpi_Var%n_proc
        n_local = Tdomain%n_elem

        allocate(m_base_local(0:n_local-1))
        allocate(m_work_local(0:n_local-1))

        allocate(recv_counts(0:n_procs-1))
        allocate(displs(0:n_procs-1))
        call MPI_Allgather(n_local, 1, MPI_INTEGER, recv_counts, 1, MPI_INTEGER, Tdomain%communicateur, ierr)

        displs(0) = 0
        do i = 1, n_procs - 1
            displs(i) = displs(i-1) + recv_counts(i-1)
        end do
        n_global = sum(recv_counts)

        allocate(dt_global(0:n_global-1))
        call MPI_Allgatherv(dt_elem_loc, n_local, MPI_DOUBLE_PRECISION, &
                            dt_global, recv_counts, displs, MPI_DOUBLE_PRECISION, &
                            Tdomain%communicateur, ierr)

        if (rg == 0) then
            call quicksort_real(dt_global, 0, n_global - 1)
        end if
        call MPI_Bcast(dt_global, n_global, MPI_DOUBLE_PRECISION, 0, Tdomain%communicateur, ierr)

        allocate(dt_candidates(1:NUM_PERCENTILES))
        do c = 1, NUM_PERCENTILES
            cand_idx = min(n_global - 1, max(0, int(PERCENTILES(c) * real(n_global, fpp))))
            dt_candidates(c) = dt_global(cand_idx)
        end do

        deallocate(dt_global, recv_counts, displs)

        min_cost = huge(1._fpp)
        dt_target_opt = dt_candidates(1)
        elem_order_opt = 1

        do c = 1, NUM_PERCENTILES
            dt_target_cand = dt_candidates(c)

            call compute_element_base_orders(dt_elem_loc, dt_target_cand, m_base_local, cand_feasible)

            ! Check global feasibility across all MPI ranks
            feasible_loc = merge(1, 0, cand_feasible)
            call MPI_Allreduce(feasible_loc, feasible_g, 1, MPI_INTEGER, MPI_MIN, Tdomain%communicateur, ierr)
            if (feasible_g == 0) then
                if (rg == 0) call log_candidate_eval(c, dt_target_cand, 0._fpp, 0._fpp, .false., .false.)
                cycle
            end if

            m_work_local = m_base_local
            call resolver_halos_from_orders(Tdomain, m_work_local)

            local_cost = 0._fpp
            do n = 0, n_local - 1
                local_cost = local_cost + real(m_work_local(n), fpp)
            end do

            call MPI_Allreduce(local_cost, global_cost, 1, MPI_DOUBLE_PRECISION, MPI_SUM, Tdomain%communicateur, ierr)

            total_cost = (1._fpp / dt_target_cand) * global_cost

            if (total_cost < min_cost) then
                min_cost = total_cost
                dt_target_opt = dt_target_cand
                elem_order_opt = m_work_local
                if (rg == 0) call log_candidate_eval(c, dt_target_cand, global_cost, total_cost, .true., .true.)
            else
                if (rg == 0) call log_candidate_eval(c, dt_target_cand, global_cost, total_cost, .false., .true.)
            end if
        end do

        deallocate(dt_candidates, m_base_local, m_work_local)

        if (rg == 0) then
            call log_optimal_selection(dt_target_opt, min_cost, 0._fpp, 1.0_fpp)
        end if
    end subroutine optimize_cost_and_orders

    function fact2k(k) result(f)
        implicit none
        integer, intent(in) :: k
        real(fpp) :: f
        integer :: i
        f = 1._fpp
        do i = 2, 2*k
            f = f * real(i,fpp)
        end do
    end function fact2k

    !> Pick which elements need the order-m correction (the local
    !! worst MODIFIED_FRACTION by critical dt) plus a MODIFIED_BUFFER_MULT*
    !! dist_max geometric ring around them, so NewmarkModified can add the
    !! correction only there while the rest of the mesh runs plain order-1
    !! Newmark at a LARGER dt (dt_target_g) instead of everyone paying for
    !! the single worst element (the old dt_min-driven formula this
    !! replaces). Single-rank only (enforced by check_modified_newmark) --
    !! no shared rank-boundary face/vertex can end up disagreeing on
    !! %modified between two ranks.
    subroutine select_modified_region(Tdomain, dt_elem_loc, n_included, mm, dt_target_g)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        real(fpp), dimension(0:Tdomain%n_elem-1), intent(in) :: dt_elem_loc
        integer, intent(in) :: n_included, mm
        real(fpp), intent(out) :: dt_target_g

        real(fpp), dimension(:), allocatable :: sorted_dt
        real(fpp), dimension(:,:), allocatable :: centroid
        logical, dimension(:), allocatable :: is_core
        integer :: n, k, i, n_worst, ierr, n_valid, n_core, n_modified
        real(fpp) :: dt_target_loc, zc, need_ratio, dist2, buf2
        logical :: any_unstable, any_unstable_g

        ! 1) Local target dt: sort the valid (non-PML/DG) critical dts
        ! ascending; the boundary just past the worst MODIFIED_FRACTION is
        ! the dt every "typical" (non-core) element already tolerates at
        ! order 1.
        allocate(sorted_dt(0:max(n_included,1)-1))
        n_valid = 0
        do n = 0, Tdomain%n_elem - 1
            if (dt_elem_loc(n) < 0._fpp) cycle
            sorted_dt(n_valid) = dt_elem_loc(n)
            n_valid = n_valid + 1
        end do

        if (n_valid == 0) then
            dt_target_loc = huge(1._fpp)
        else
            call quicksort_real(sorted_dt, 0, n_valid - 1)
            if (n_valid == 1) then
                dt_target_loc = sorted_dt(0)
            else
                n_worst = ceiling(MODIFIED_FRACTION * real(n_valid,fpp))
                n_worst = max(1, min(n_worst, n_valid - 1))
                dt_target_loc = sorted_dt(n_worst)
            end if
        end if
        deallocate(sorted_dt)

        call MPI_AllReduce(dt_target_loc, dt_target_g, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)

        ! 2) Flag "core": elements that cannot survive dt_target_g at order 1.
        allocate(is_core(0:Tdomain%n_elem-1))
        is_core = .false.
        n_core = 0
        do n = 0, Tdomain%n_elem - 1
            Tdomain%specel(n)%modified = .false.
            if (dt_elem_loc(n) < 0._fpp) cycle
            if (dt_elem_loc(n) < dt_target_g) then
                is_core(n) = .true.
                Tdomain%specel(n)%modified = .true.
                n_core = n_core + 1
            end if
        end do

        ! 3) Refuse to run silently unstable: verify the configured order
        ! actually stabilizes every core element at dt_target_g.
        zc = ModifiedEquationZCrit(mm)
        any_unstable = .false.
        do n = 0, Tdomain%n_elem - 1
            if (.not. is_core(n)) cycle
            need_ratio = 4._fpp * (dt_target_g/dt_elem_loc(n))**2
            if (zc < need_ratio) then
                any_unstable = .true.
                write(*,*) "[Irons dt_crit] WARNING: element ", n, " dt_elem=", dt_elem_loc(n), &
                    " cannot be stabilized at dt_target=", dt_target_g, " by modified_order=", mm
            end if
        end do
        call MPI_AllReduce(any_unstable, any_unstable_g, 1, MPI_LOGICAL, MPI_LOR, Tdomain%communicateur, ierr)
        if (any_unstable_g) then
            call MPI_Barrier(Tdomain%communicateur, ierr)
            STOP "ERROR : newmark_modified regional selection -- modified_order too low for the elements it selected as core (see WARNING lines above). Raise newmark_modified_order or lower MODIFIED_FRACTION in irons_dtcrit.F90."
        end if

        ! 4) Geometric buffer: MODIFIED_BUFFER_MULT*dist_max(core) ring
        ! around each core element also gets flagged, for a smooth
        ! order-1/order-m transition instead of an abrupt one at the core
        ! boundary.
        call dist_max_elem(Tdomain)
        allocate(centroid(0:1, 0:Tdomain%n_elem-1))
        do n = 0, Tdomain%n_elem - 1
            centroid(:,n) = element_centroid(Tdomain, n)
        end do

        do n = 0, Tdomain%n_elem - 1
            if (Tdomain%specel(n)%modified) cycle  ! already core
            do i = 0, Tdomain%n_elem - 1
                if (.not. is_core(i)) cycle
                buf2 = (MODIFIED_BUFFER_MULT * Tdomain%specel(i)%dist_max)**2
                dist2 = (centroid(0,n)-centroid(0,i))**2 + (centroid(1,n)-centroid(1,i))**2
                if (dist2 <= buf2) then
                    Tdomain%specel(n)%modified = .true.
                    exit
                end if
            end do
        end do
        deallocate(centroid, is_core)

        ! 5) Derive face/vertex flags: modified only if EVERY touching
        ! element is modified -- a shared DOF cannot be half-corrected.
        do n = 0, Tdomain%n_face - 1
            Tdomain%sFace(n)%modified = .true.
        end do
        do n = 0, Tdomain%n_vertex - 1
            Tdomain%sVertex(n)%modified = .true.
        end do
        do n = 0, Tdomain%n_elem - 1
            if (Tdomain%specel(n)%modified) cycle
            do k = 0, 3
                Tdomain%sFace(Tdomain%specel(n)%Near_Face(k))%modified = .false.
                Tdomain%sVertex(Tdomain%specel(n)%Near_Vertex(k))%modified = .false.
            end do
        end do

        n_modified = count(Tdomain%specel(:)%modified)
        if (Tdomain%Mpi_var%my_rank == 0) then
            write(*,*) "[Irons dt_crit] regional selection: ", n_core, " core element(s) (dt < target), ", &
                n_modified, " total modified (core+buffer) / ", Tdomain%n_elem
        end if
    end subroutine select_modified_region

    !> Element centroid (average of its Control_Nodes physical coordinates) --
    !! cheap geometric proxy, consistent with dist_max_elem's own use of
    !! Control_Nodes/Coord_nodes (no dependency on GLL/Jacobian being built).
    function element_centroid(Tdomain, n) result(c)
        use sdomain
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: n
        real(fpp), dimension(0:1) :: c
        integer :: i, ipoint

        c = 0._fpp
        do i = 0, Tdomain%n_nodes - 1
            ipoint = Tdomain%specel(n)%Control_Nodes(i)
            c = c + Tdomain%Coord_nodes(0:1, ipoint)
        end do
        c = c / real(Tdomain%n_nodes, fpp)
    end function element_centroid

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
