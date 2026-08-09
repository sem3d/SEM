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
    use m_mod_eq_zcrit
    use m_cost_optimizer_common
    use m_modified_newmark_logger
    implicit none

    integer, parameter :: IRONS_MAX_ITER = 2000
    real(fpp), parameter :: IRONS_RTOL = 1e-10_fpp

    ! Regional modified-equation Newmark: BFS graph-distance (face-adjacency
    ! hops) from each element to the nearest %modified element, capped at
    ! modified_order-1, plus the per-element size/centroid proxies used to
    ! grow the geometric buffer ring in select_modified_region. See
    ! SEM2D/SRC/irons_dtcrit.F90's module-level comment on elem_hop for why
    ! this must exist (correctness of the shrinking-halo matrix-free
    ! iteration in NewmarkModified.f90's apply_a_regional).
    integer, dimension(:), allocatable, save :: elem_hop
    real(fpp), dimension(:), allocatable, save :: elem_dist_max
    real(fpp), dimension(:,:), allocatable, save :: elem_centroid3d_arr

    ! Regional modified-equation Newmark: fraction of (local, single-rank)
    ! elements with the smallest critical dt that seed the "core" needing the
    ! order-m correction, and the geometric buffer radius around each core
    ! element (multiple of that element's dist_max proxy) that also gets
    ! corrected for a smooth order-1/order-m transition. ponytail: hardcoded
    ! rather than wired into input.spec, same as SEM2D's irons_dtcrit.F90 --
    ! promote to a config key if these ever need per-run tuning.
    real(fpp), parameter :: MODIFIED_FRACTION = 0.05_fpp
    real(fpp), parameter :: MODIFIED_BUFFER_MULT = 3.0_fpp

    ! Per-domain block/DOF activation derived by select_modified_region:
    ! block_active_*(j)%idx = block indices touched by apply_a_regional's
    ! j-th (of mm) matrix-free iteration (>=1 element in that block has
    ! hop<=mm-j); dof_corrected_* = DOFs receiving the correction (mm-th
    ! iteration only, i.e. the final core+buffer set), false for any DOF
    ! shared with a non-modified element or a solid-fluid interface.
    type :: idx_list_t
        integer, dimension(:), allocatable :: idx
    end type idx_list_t
    type(idx_list_t), dimension(:), allocatable, save :: block_active_sdom, block_active_fdom
    ! dof_touched_*(j)%idx = every DOF belonging to any element in a
    ! block_active_*(j) block (i.e. the whole block, including any "bonus"
    ! lane sharing the block with a truly hop<=dmax element -- harmless
    ! over-computation, see NewmarkModified.f90's apply_a_regional). Used to
    ! zero/mass-multiply exactly the right subset each iteration without
    ! disturbing untouched DOFs, which must keep their previous iteration's
    ! value (required by the shrinking-halo invariant).
    type(idx_list_t), dimension(:), allocatable, save :: dof_touched_sdom, dof_touched_fdom
    logical, dimension(:), allocatable, save :: dof_corrected_sdom, dof_corrected_fdom

contains

    !> Entry point, called from Compute_Courant (courant.f90) once at startup,
    !! after init_materials has filled Density_/Lambda_/Mu_/Jacob_/InvGrad_.
    subroutine compute_irons_dtcrit(Tdomain, rg)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer :: ierr, n, dom, mm
        real(fpp) :: dt_loc, dt_min, dt_used, ratio, dt_courant_irons
        real(fpp) :: dt_target_g, dt_courant_irons_order
        integer :: n_not_converged, n_skipped, n_total, n_covered
        integer :: n_not_converged_g, n_skipped_g, n_total_g
        real(fpp), dimension(:), allocatable :: dt_elem_loc
        real(fpp), dimension(:), allocatable :: dt_elem_sdom, dt_elem_fdom
        real(fpp) :: dt_sum_loc, dt_max_loc

        dt_loc = huge(1._fpp)
        n_not_converged = 0
        n_covered = 0
        dt_sum_loc = 0._fpp
        dt_max_loc = 0._fpp
        allocate(dt_elem_loc(0:Tdomain%n_elem-1))
        dt_elem_loc = -1._fpp

        if (Tdomain%sdom%nbelem > 0) then
            allocate(dt_elem_sdom(0:Tdomain%sdom%nbelem-1))
            call irons_domain_solid(Tdomain%sdom, dt_loc, n_not_converged, dt_elem_sdom)
            n_covered = n_covered + Tdomain%sdom%nbelem
        end if
        if (Tdomain%fdom%nbelem > 0) then
            allocate(dt_elem_fdom(0:Tdomain%fdom%nbelem-1))
            call irons_domain_fluid(Tdomain%fdom, dt_loc, n_not_converged, dt_elem_fdom)
            n_covered = n_covered + Tdomain%fdom%nbelem
        end if

        ! Scatter per-domain-local dt back into the global-element-indexed
        ! array via specel(n)%domain/%lnum (lnum is 0-based within that
        ! element's domain) -- needed by select_modified_region later.
        do n = 0, Tdomain%n_elem - 1
            dom = Tdomain%specel(n)%domain
            if (dom == DM_SOLID_CG) then
                dt_elem_loc(n) = dt_elem_sdom(Tdomain%specel(n)%lnum)
            else if (dom == DM_FLUID_CG) then
                dt_elem_loc(n) = dt_elem_fdom(Tdomain%specel(n)%lnum)
            end if
            if (dt_elem_loc(n) >= 0._fpp) then
                dt_sum_loc = dt_sum_loc + dt_elem_loc(n)
                dt_max_loc = max(dt_max_loc, dt_elem_loc(n))
            end if
        end do

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
            call MPI_AllReduce(n_covered, n_included_g, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur, ierr)
            if (n_included_g > 0) then
                dt_mean = dt_sum_g / real(n_included_g,fpp)
                call report_local_order_requirements(Tdomain, dt_elem_loc, "mean", dt_mean, dt_min, 5)
                call report_local_order_requirements(Tdomain, dt_elem_loc, "max (best case)", dt_max_g, dt_min, 5)
            end if
        end block

        if (Tdomain%TimeD%modified) then
            if (Tdomain%nb_procs > 1) then
                STOP "ERROR : newmark_modified=true does not support MPI (n_proc>1) yet -- run with a single process."
            end if
            ! No beta/gamma parametrization check here (unlike SEM2D's
            ! NewmarkModified.F90): SEM3D's newmark_corrector_solid/fluid
            ! hardcode the explicit leapfrog update (Veloc += dt*Forces;
            ! Depla += dt*Veloc) unconditionally -- alpha/beta/gamma from
            ! time_scheme are silently overridden to 0.5/0.5/1 right after
            ! being read (read_input.f90) and never actually parametrize the
            ! corrector, so the base (classic) Newmark call NewmarkModified
            ! reuses for its k=1 term is always exactly the leapfrog step
            ! this recursion assumes, regardless of the input.spec values.

            mm = Tdomain%TimeD%modified_order
            call select_modified_region(Tdomain, dt_elem_loc, n_covered, mm, dt_target_g)
            dt_courant_irons_order = Tdomain%TimeD%courant * dt_target_g

            if (rg == 0) then
                write (*,*) "[Irons dt_crit] newmark_modified_order = ", mm
                write (*,*) "[Irons dt_crit] regional target dt (order-1, drives the simulation) = ", dt_target_g
                write (*,*) "[Irons dt_crit] regional target dt with courant safety factor = ", dt_courant_irons_order
            end if

            if (dt_courant_irons_order <= 0._fpp) then
                write (*,*) "Your regional target dt is zero : verify it"
                stop
            end if
            Tdomain%TimeD%dtmin = dt_courant_irons_order
            Tdomain%TimeD%ntimeMax = int(Tdomain%TimeD%duration/Tdomain%TimeD%dtmin)
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
    !! pay at the boundary of the "needs boost" region. Ported verbatim
    !! from SEM2D/SRC/irons_dtcrit.F90.
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

        rg = Tdomain%rank
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

    !> Solid CG domain (iso/aniso, unrelaxed elastic operator -- attenuation
    !! ignored on purpose: it only makes the medium look stiffer, i.e. dt
    !! smaller, the safe direction). Reuses calcul_forces_iso/calcul_forces_aniso
    !! (SEM3D/SRC/Solid/calcul_forces.F90) unmodified as the K_e*w operator.
    subroutine irons_domain_solid(dom, dt_loc, n_not_converged, dt_elem)
        use champs_solid
        use m_calcul_forces
        implicit none
        type(domain_solid), intent(inout) :: dom
        real(fpp), intent(inout) :: dt_loc
        integer, intent(inout) :: n_not_converged
        real(fpp), dimension(0:dom%nbelem-1), intent(out) :: dt_elem
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
                dt_elem(bnum*VCHUNK + ee) = 2._fpp/omega
            end do
        end do
    end subroutine irons_domain_solid

    !> Fluid CG domain (velocity-potential formulation, iso/aniso density).
    !! Reuses calcul_forces_fluid (dispatches iso/aniso internally) as the
    !! K_e*w operator; local mass = Whei*Jacob_ divided by the fluid bulk
    !! modulus, stored in the Lambda_ field (same field name reused across
    !! domain types).
    subroutine irons_domain_fluid(dom, dt_loc, n_not_converged, dt_elem)
        use champs_fluid
        use m_calcul_forces_fluid
        implicit none
        type(domain_fluid), intent(inout) :: dom
        real(fpp), intent(inout) :: dt_loc
        integer, intent(inout) :: n_not_converged
        real(fpp), dimension(0:dom%nbelem-1), intent(out) :: dt_elem
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
                dt_elem(bnum*VCHUNK + ee) = 2._fpp/omega
            end do
        end do
    end subroutine irons_domain_fluid

    !> S(z) = sum_{k=1}^m 2*(-z)^k/(2k)!, the modified-equation amplification
    !! term for NewmarkModified.f90's order-m recursion (single mode
    !! A*v=-omega^2*v, z=(omega*dt)^2). Stability (bounded, oscillatory
    !! roots) holds iff S(z) in [-4,0]. Pure math, ported verbatim from
    !! SEM2D/SRC/irons_dtcrit.F90.
    function ModifiedEquationS(z, m) result(s)
        implicit none
        real(fpp), intent(in) :: z
        integer, intent(in) :: m
        real(fpp) :: s
        integer :: k
        s = 0._fpp
        do k = 1, m
            s = s + 2._fpp*(-z)**k/fact2k(k)
        end do
    end function ModifiedEquationS

    !> Smallest positive z where S(z) exits [-4,0]: coarse forward scan
    !! (step 0.01) to bracket the crossing, then bisection. z_crit(1)=4
    !! recovers the plain-leapfrog CFL bound dt=2/omega_max. Order does NOT
    !! increase z_crit monotonically: [4, 12, 7.57, 21.48, 9.53] for m=1..5
    !! (odd orders are less stable than their even neighbours) - verified
    !! against direct scalar time-stepping in the labcorrea/FEM MATLAB
    !! prototype (StabilityClass.ModifiedEquationZCrit).
    function ModifiedEquationZCrit(m) result(zc)
        implicit none
        integer, intent(in) :: m
        real(fpp) :: zc
        real(fpp) :: dz, z, zlo, zhi, zmid
        integer :: it

        dz = 0.01_fpp
        z = 0._fpp
        do while (mod_eq_viol(z,m) <= 0._fpp)
            z = z + dz
            if (z > 1000._fpp) then
                write(*,*) "ModifiedEquationZCrit: no instability boundary found up to z=1000 for order m=", m
                stop
            endif
        end do

        zlo = z - dz; zhi = z
        do it = 1, 100
            zmid = 0.5_fpp*(zlo+zhi)
            if (mod_eq_viol(zmid,m) <= 0._fpp) then
                zlo = zmid
            else
                zhi = zmid
            endif
        end do
        zc = 0.5_fpp*(zlo+zhi)
    end function ModifiedEquationZCrit

    function mod_eq_viol(z, m) result(v)
        implicit none
        real(fpp), intent(in) :: z
        integer, intent(in) :: m
        real(fpp) :: v, s
        s = ModifiedEquationS(z,m)
        v = max(s, -4._fpp-s)
    end function mod_eq_viol

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

    !> Element centroid (average of its Control_nodes physical coordinates).
    function element_centroid3d(Tdomain, n) result(c)
        use sdomain
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: n
        real(fpp), dimension(0:2) :: c
        integer :: i, ipoint, nnodes
        nnodes = size(Tdomain%specel(n)%Control_nodes)
        c = 0._fpp
        do i = 0, nnodes - 1
            ipoint = Tdomain%specel(n)%Control_nodes(i)
            c = c + Tdomain%Coord_nodes(0:2, ipoint)
        end do
        c = c / real(nnodes, fpp)
    end function element_centroid3d

    !> Element size proxy: max pairwise distance between Control_nodes
    !! (cheap O(n_nodes^2) with n_nodes = 8 or 27, computed once at startup).
    !! No SEM3D equivalent of SEM2D's dist_max field exists, so this is
    !! computed directly instead of reusing a stored value.
    function element_dist_max3d(Tdomain, n) result(dmax)
        use sdomain
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: n
        real(fpp) :: dmax
        integer :: i, j, nnodes, ip, jp
        real(fpp) :: d
        nnodes = size(Tdomain%specel(n)%Control_nodes)
        dmax = 0._fpp
        do i = 0, nnodes - 2
            ip = Tdomain%specel(n)%Control_nodes(i)
            do j = i+1, nnodes - 1
                jp = Tdomain%specel(n)%Control_nodes(j)
                d = sqrt(sum((Tdomain%Coord_nodes(0:2,ip) - Tdomain%Coord_nodes(0:2,jp))**2))
                if (d > dmax) dmax = d
            end do
        end do
    end function element_dist_max3d

    !> Precompute elem_hop(n) = graph-distance (face-adjacency hops) from
    !! element n to the nearest %modified element, capped at mm-1, -1 if
    !! never reached. Mirrors SEM2D's build_regional_halo but walks up to 6
    !! faces per (hex) element instead of 4, via sFace%elem_0/elem_1 (SEM3D
    !! faces store the two neighbouring LOCAL element numbers directly,
    !! unlike SEM2D's sFace%Near_Element which SEM3D has no equivalent of).
    subroutine build_regional_halo3d(Tdomain, mm)
        use sdomain
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: mm
        integer, dimension(:), allocatable :: queue
        integer :: qhead, qtail, n, k, nf, e0, e1, other, d

        if (allocated(elem_hop)) deallocate(elem_hop)
        allocate(elem_hop(0:Tdomain%n_elem-1))
        elem_hop = -1

        allocate(queue(0:Tdomain%n_elem-1))
        qhead = 0; qtail = 0
        do n = 0, Tdomain%n_elem - 1
            if (Tdomain%specel(n)%modified) then
                elem_hop(n) = 0
                queue(qtail) = n
                qtail = qtail + 1
            end if
        end do

        do while (qhead < qtail)
            n = queue(qhead)
            qhead = qhead + 1
            d = elem_hop(n)
            if (d >= mm - 1) cycle  ! do not expand past the widest ever-needed hop
            do k = 0, 5
                nf = Tdomain%specel(n)%Near_Faces(k)
                e0 = Tdomain%sFace(nf)%elem_0
                e1 = Tdomain%sFace(nf)%elem_1
                if (e0 == n) then
                    other = e1
                else
                    other = e0
                end if
                if (other < 0) cycle  ! domain/mesh boundary, no neighbour on this side
                if (elem_hop(other) < 0) then
                    elem_hop(other) = d + 1
                    queue(qtail) = other; qtail = qtail + 1
                end if
            end do
        end do
        deallocate(queue)
    end subroutine build_regional_halo3d

    !> Pick which sdom/fdom elements need the order-m correction (worst
    !! MODIFIED_FRACTION by critical dt, pooled across both domains) plus a
    !! MODIFIED_BUFFER_MULT*dist_max geometric ring around them. Single-rank
    !! only (enforced by the caller) -- dt_target_g is simply the local
    !! value, no AllReduce needed (unlike SEM2D, which supports MPI here).
    subroutine select_modified_region(Tdomain, dt_elem_loc, n_included, mm, dt_target_g)
        use sdomain
        implicit none
        type(domain), intent(inout) :: Tdomain
        real(fpp), dimension(0:Tdomain%n_elem-1), intent(in) :: dt_elem_loc
        integer, intent(in) :: n_included, mm
        real(fpp), intent(out) :: dt_target_g

        real(fpp), dimension(:), allocatable :: sorted_dt
        logical, dimension(:), allocatable :: is_core
        integer :: n, i, n_worst, n_valid, n_core, n_modified
        real(fpp) :: dt_target, zc, need_ratio, dist2, buf2
        logical :: any_unstable

        ! 1) Target dt: sort valid critical dts ascending; the boundary just
        ! past the worst MODIFIED_FRACTION is the dt every "typical"
        ! (non-core) element already tolerates at order 1.
        allocate(sorted_dt(0:max(n_included,1)-1))
        n_valid = 0
        do n = 0, Tdomain%n_elem - 1
            if (dt_elem_loc(n) < 0._fpp) cycle
            sorted_dt(n_valid) = dt_elem_loc(n)
            n_valid = n_valid + 1
        end do
        if (n_valid == 0) then
            dt_target = huge(1._fpp)
        else
            call quicksort_real(sorted_dt, 0, n_valid - 1)
            if (n_valid == 1) then
                dt_target = sorted_dt(0)
            else
                n_worst = ceiling(MODIFIED_FRACTION * real(n_valid,fpp))
                n_worst = max(1, min(n_worst, n_valid - 1))
                dt_target = sorted_dt(n_worst)
            end if
        end if
        deallocate(sorted_dt)
        dt_target_g = dt_target  ! mono-rank: local == global, no AllReduce needed

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
        if (any_unstable) then
            STOP "ERROR : newmark_modified regional selection -- modified_order too low for the elements it selected as core (see WARNING lines above). Raise newmark_modified_order or lower MODIFIED_FRACTION in irons_dtcrit.F90."
        end if

        ! 4) Geometric buffer: MODIFIED_BUFFER_MULT*dist_max(core) ring
        ! around each core element also gets flagged, for a smooth
        ! order-1/order-m transition instead of an abrupt one at the core
        ! boundary.
        if (allocated(elem_dist_max)) deallocate(elem_dist_max)
        if (allocated(elem_centroid3d_arr)) deallocate(elem_centroid3d_arr)
        allocate(elem_dist_max(0:Tdomain%n_elem-1))
        allocate(elem_centroid3d_arr(0:2, 0:Tdomain%n_elem-1))
        elem_dist_max = 0._fpp
        do n = 0, Tdomain%n_elem - 1
            elem_centroid3d_arr(:,n) = element_centroid3d(Tdomain, n)
            if (is_core(n)) elem_dist_max(n) = element_dist_max3d(Tdomain, n)
        end do
        do n = 0, Tdomain%n_elem - 1
            if (Tdomain%specel(n)%modified) cycle  ! already core
            do i = 0, Tdomain%n_elem - 1
                if (.not. is_core(i)) cycle
                buf2 = (MODIFIED_BUFFER_MULT * elem_dist_max(i))**2
                dist2 = sum((elem_centroid3d_arr(:,n) - elem_centroid3d_arr(:,i))**2)
                if (dist2 <= buf2) then
                    Tdomain%specel(n)%modified = .true.
                    exit
                end if
            end do
        end do
        deallocate(is_core)

        ! 5) Build the BFS halo now that %modified is final.
        call build_regional_halo3d(Tdomain, mm)

        ! 6) Derive block_active(j) per domain and dof_corrected (mm only).
        call derive_block_and_dof_activation(Tdomain, mm)

        n_modified = count(Tdomain%specel(:)%modified)
        write(*,*) "[Irons dt_crit] regional selection: ", n_core, " core element(s) (dt < target), ", &
            n_modified, " total modified (core+buffer) / ", Tdomain%n_elem
    end subroutine select_modified_region

    !> For each domain (sdom, fdom) and each iteration j=1..mm, build the
    !! list of block indices touched (>=1 element with hop<=mm-j maps into
    !! that block) -- the force kernels (calcul_forces_iso/aniso,
    !! calcul_forces_fluid) only ever run on a whole block, so this is the
    !! coarsest granularity "regional" activation can use for them. Also
    !! builds dof_corrected (mm-th iteration only): a DOF is corrected only
    !! if every element touching it (via Idom_) is %modified, and it is not
    !! part of a solid-fluid interface (those coexist untouched, same as
    !! PML/CPML/DG which are never %modified in the first place since
    !! compute_irons_dtcrit only fills dt_elem_loc for sdom/fdom).
    subroutine derive_block_and_dof_activation(Tdomain, mm)
        use sdomain
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: mm
        integer :: n, j, dom, bnum, ee, i, jy, k, idx, dmax
        logical, dimension(:), allocatable :: block_touched
        logical, dimension(:), allocatable :: dof_touched_mask

        if (allocated(block_active_sdom)) deallocate(block_active_sdom)
        if (allocated(block_active_fdom)) deallocate(block_active_fdom)
        if (allocated(dof_touched_sdom)) deallocate(dof_touched_sdom)
        if (allocated(dof_touched_fdom)) deallocate(dof_touched_fdom)
        allocate(block_active_sdom(mm), block_active_fdom(mm))
        allocate(dof_touched_sdom(mm), dof_touched_fdom(mm))

        do j = 1, mm
            dmax = mm - j
            if (Tdomain%sdom%nblocks > 0) then
                allocate(block_touched(0:Tdomain%sdom%nblocks-1))
                block_touched = .false.
                do n = 0, Tdomain%n_elem - 1
                    if (Tdomain%specel(n)%domain /= DM_SOLID_CG) cycle
                    if (elem_hop(n) < 0 .or. elem_hop(n) > dmax) cycle
                    bnum = Tdomain%specel(n)%lnum / VCHUNK
                    block_touched(bnum) = .true.
                end do
                call pack_true_indices(block_touched, block_active_sdom(j)%idx)

                ! dof_touched_sdom(j): every DOF of every element in a
                ! touched block (whole block, see module-level comment).
                allocate(dof_touched_mask(0:Tdomain%sdom%nglltot-1))
                dof_touched_mask = .false.
                do n = 0, Tdomain%n_elem - 1
                    if (Tdomain%specel(n)%domain /= DM_SOLID_CG) cycle
                    bnum = Tdomain%specel(n)%lnum / VCHUNK
                    if (.not. block_touched(bnum)) cycle
                    ee = mod(Tdomain%specel(n)%lnum, VCHUNK)
                    do k = 0, Tdomain%sdom%ngll - 1
                        do jy = 0, Tdomain%sdom%ngll - 1
                            do i = 0, Tdomain%sdom%ngll - 1
                                dof_touched_mask(Tdomain%sdom%Idom_(i,jy,k,bnum,ee)) = .true.
                            end do
                        end do
                    end do
                end do
                call pack_true_indices(dof_touched_mask, dof_touched_sdom(j)%idx)
                deallocate(dof_touched_mask, block_touched)
            else
                allocate(block_active_sdom(j)%idx(0))
                allocate(dof_touched_sdom(j)%idx(0))
            end if

            if (Tdomain%fdom%nblocks > 0) then
                allocate(block_touched(0:Tdomain%fdom%nblocks-1))
                block_touched = .false.
                do n = 0, Tdomain%n_elem - 1
                    if (Tdomain%specel(n)%domain /= DM_FLUID_CG) cycle
                    if (elem_hop(n) < 0 .or. elem_hop(n) > dmax) cycle
                    bnum = Tdomain%specel(n)%lnum / VCHUNK
                    block_touched(bnum) = .true.
                end do
                call pack_true_indices(block_touched, block_active_fdom(j)%idx)

                allocate(dof_touched_mask(0:Tdomain%fdom%nglltot-1))
                dof_touched_mask = .false.
                do n = 0, Tdomain%n_elem - 1
                    if (Tdomain%specel(n)%domain /= DM_FLUID_CG) cycle
                    bnum = Tdomain%specel(n)%lnum / VCHUNK
                    if (.not. block_touched(bnum)) cycle
                    ee = mod(Tdomain%specel(n)%lnum, VCHUNK)
                    do k = 0, Tdomain%fdom%ngll - 1
                        do jy = 0, Tdomain%fdom%ngll - 1
                            do i = 0, Tdomain%fdom%ngll - 1
                                dof_touched_mask(Tdomain%fdom%Idom_(i,jy,k,bnum,ee)) = .true.
                            end do
                        end do
                    end do
                end do
                call pack_true_indices(dof_touched_mask, dof_touched_fdom(j)%idx)
                deallocate(dof_touched_mask, block_touched)
            else
                allocate(block_active_fdom(j)%idx(0))
                allocate(dof_touched_fdom(j)%idx(0))
            end if
        end do

        ! dof_corrected: start from "every owning element modified", derived
        ! by clearing every DOF touched by a non-modified sdom/fdom element.
        if (allocated(dof_corrected_sdom)) deallocate(dof_corrected_sdom)
        if (allocated(dof_corrected_fdom)) deallocate(dof_corrected_fdom)
        allocate(dof_corrected_sdom(0:Tdomain%sdom%nglltot-1))
        allocate(dof_corrected_fdom(0:Tdomain%fdom%nglltot-1))
        dof_corrected_sdom = .true.
        dof_corrected_fdom = .true.

        do n = 0, Tdomain%n_elem - 1
            if (Tdomain%specel(n)%modified) cycle
            dom = Tdomain%specel(n)%domain
            if (dom /= DM_SOLID_CG .and. dom /= DM_FLUID_CG) cycle
            bnum = Tdomain%specel(n)%lnum / VCHUNK
            ee = mod(Tdomain%specel(n)%lnum, VCHUNK)
            if (dom == DM_SOLID_CG) then
                do k = 0, Tdomain%sdom%ngll - 1
                    do jy = 0, Tdomain%sdom%ngll - 1
                        do i = 0, Tdomain%sdom%ngll - 1
                            dof_corrected_sdom(Tdomain%sdom%Idom_(i,jy,k,bnum,ee)) = .false.
                        end do
                    end do
                end do
            else
                do k = 0, Tdomain%fdom%ngll - 1
                    do jy = 0, Tdomain%fdom%ngll - 1
                        do i = 0, Tdomain%fdom%ngll - 1
                            dof_corrected_fdom(Tdomain%fdom%Idom_(i,jy,k,bnum,ee)) = .false.
                        end do
                    end do
                end do
            end if
        end do

        ! Exclude solid-fluid interface DOFs entirely (they coexist
        ! untouched, same policy as PML/CPML/DG). surf0%map(i) and
        ! surf1%map(i) refer to the same physical interface point on the
        ! solid/fluid sides respectively, so surf0%nbtot == surf1%nbtot.
        do idx = 0, Tdomain%SF%intSolFlu%surf0%nbtot - 1
            dof_corrected_sdom(Tdomain%SF%intSolFlu%surf0%map(idx)) = .false.
        end do
        do idx = 0, Tdomain%SF%intSolFlu%surf1%nbtot - 1
            dof_corrected_fdom(Tdomain%SF%intSolFlu%surf1%map(idx)) = .false.
        end do
    end subroutine derive_block_and_dof_activation

    !> Indices n with mask(n) true.
    subroutine pack_true_indices(mask, idx)
        implicit none
        logical, dimension(0:), intent(in) :: mask
        integer, dimension(:), allocatable, intent(out) :: idx
        integer :: n, cnt
        cnt = count(mask)
        allocate(idx(cnt))
        cnt = 0
        do n = 0, size(mask) - 1
            if (mask(n)) then
                cnt = cnt + 1
                idx(cnt) = n
            end if
        end do
    end subroutine pack_true_indices

    !> In-place ascending quicksort (Hoare partition) -- used to find the
    !! MODIFIED_FRACTION percentile boundary in select_modified_region.
    recursive subroutine quicksort_real(a, lo, hi)
        implicit none
        real(fpp), dimension(0:), intent(inout) :: a
        integer, intent(in) :: lo, hi
        integer :: i, j
        real(fpp) :: pivot, tmp
        if (lo >= hi) return
        pivot = a((lo+hi)/2)
        i = lo; j = hi
        do
            do while (a(i) < pivot)
                i = i + 1
            end do
            do while (a(j) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                tmp = a(i); a(i) = a(j); a(j) = tmp
                i = i + 1; j = j - 1
            end if
            if (i > j) exit
        end do
        if (lo < j) call quicksort_real(a, lo, j)
        if (i < hi) call quicksort_real(a, i, hi)
    end subroutine quicksort_real

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
