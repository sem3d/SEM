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
    use m_modified_newmark_logger
    use m_cost_optimizer_common
    implicit none

    integer, parameter :: IRONS_MAX_ITER = 2000
    real(fpp), parameter :: IRONS_RTOL = 1e-10_fpp

    ! Regional modified-equation Newmark: elem_hop(n) = mm - elem_order(n),
    ! the per-element correction order chosen by optimize_cost_and_orders
    ! (already grown to a consistent halo, including across rank boundaries,
    ! by resolver_halos_from_orders_3d), re-expressed as a hop-from-deepest
    ! distance so derive_block_and_dof_activation's elem_hop(n)<=dmax check
    ! is algebraically identical to elem_order(n)>=j. See
    ! SEM2D/SRC/irons_dtcrit.F90's module-level comment on elem_hop for the
    ! shrinking-halo correctness argument this mirrors.
    integer, dimension(:), allocatable, save :: elem_hop

    ! Per-domain block/DOF activation derived by derive_block_and_dof_activation:
    ! block_active_*(j)%idx = block indices touched by apply_a_regional's
    ! j-th (of mm) matrix-free iteration (>=1 element in that block has
    ! hop<=mm-j); dof_corrected_* (declared further below) is each DOF's own
    ! correction order, capped at 1 (excluded) for any DOF shared with an
    ! order-1 element or a solid-fluid interface.
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
    ! dof_corrected_*(idx) = correction order for this DOF: MIN(elem_order)
    ! over every sdom/fdom element touching it (0 would mean "excluded",
    ! but the floor is 1 -- a DOF shared with any order-1 element gets no
    ! k=2.. correction, same as before this was graded, see
    ! derive_block_and_dof_activation). NewmarkModified.f90's correction
    ! loop gates each k=2..m term by k <= dof_corrected_*(idx).
    integer, dimension(:), allocatable, save :: dof_corrected_sdom, dof_corrected_fdom

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
        ! element's domain) -- needed by optimize_cost_and_orders later.
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
            ! No beta/gamma parametrization check here (unlike SEM2D's
            ! NewmarkModified.F90): SEM3D's newmark_corrector_solid/fluid
            ! hardcode the explicit leapfrog update (Veloc += dt*Forces;
            ! Depla += dt*Veloc) unconditionally -- alpha/beta/gamma from
            ! time_scheme are silently overridden to 0.5/0.5/1 right after
            ! being read (read_input.f90) and never actually parametrize the
            ! corrector, so the base (classic) Newmark call NewmarkModified
            ! reuses for its k=1 term is always exactly the leapfrog step
            ! this recursion assumes, regardless of the input.spec values.

            block
                integer, dimension(:), allocatable :: elem_order_opt
                integer :: mm_loc
                if (rg == 0) call log_header_modified_newmark("3D")
                allocate(elem_order_opt(0:Tdomain%n_elem-1))
                call optimize_cost_and_orders(Tdomain, dt_elem_loc, dt_target_g, elem_order_opt)
                ! maxval(elem_order_opt) is only this rank's own local
                ! elements -- every rank must agree on the SAME mm (it sizes
                ! block_active_*(mm)/dof_touched_*(mm) and drives
                ! NewmarkModified's "do k=1,mm" loop; a per-rank-different mm
                ! desyncs apply_a_regional's cross-rank exchange, which then
                ! deadlocks the first rank still iterating once its faster
                ! neighbour has moved past NewmarkModified entirely).
                mm_loc = maxval(elem_order_opt)
                call MPI_AllReduce(mm_loc, mm, 1, MPI_INTEGER, MPI_MAX, Tdomain%communicateur, ierr)
                Tdomain%TimeD%modified_order = mm
                do n = 0, Tdomain%n_elem - 1
                    Tdomain%specel(n)%modified = (elem_order_opt(n) > 1)
                end do
                ! elem_hop(n) = mm - elem_order_opt(n), directly from the
                ! optimizer's already cross-rank-consistent per-element
                ! orders (resolver_halos_from_orders_3d already grew them
                ! correctly) -- elem_hop(n)<=dmax in
                ! derive_block_and_dof_activation is algebraically identical
                ! to elem_order(n)>=j for dmax=mm-j, so no other change is
                ! needed there.
                if (allocated(elem_hop)) deallocate(elem_hop)
                allocate(elem_hop(0:Tdomain%n_elem-1))
                do n = 0, Tdomain%n_elem - 1
                    elem_hop(n) = mm - elem_order_opt(n)
                end do
                call derive_block_and_dof_activation(Tdomain, mm, elem_order_opt)
                deallocate(elem_order_opt)
            end block
            dt_courant_irons_order = Tdomain%TimeD%courant * dt_target_g

            if (rg == 0) then
                write (*,*) "[Irons dt_crit] optimal newmark_modified_order = ", mm
                write (*,*) "[Irons dt_crit] optimal target dt (drives the simulation) = ", dt_target_g
                write (*,*) "[Irons dt_crit] optimal target dt with courant safety factor = ", dt_courant_irons_order
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

    !> (2k)! -- used by NewmarkModified.f90's order-m coefficient recursion.
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

    !> Point (Control_nodes/Coord_nodes global index) -> touching-elements
    !! reverse map, used by resolver_halos_from_orders_3d's local BFS
    !! (see its own docstring for why corner-node adjacency, not
    !! Near_Faces, is the right graph to walk).
    subroutine build_point_to_elems(Tdomain, point_to_elems)
        use sdomain
        implicit none
        type(domain), intent(in) :: Tdomain
        type(idx_list_t), dimension(:), allocatable, intent(out) :: point_to_elems
        integer, dimension(:), allocatable :: point_count
        integer :: n, i, p, nnodes, ipoint

        allocate(point_count(0:Tdomain%n_glob_nodes-1))
        point_count = 0
        do n = 0, Tdomain%n_elem - 1
            nnodes = size(Tdomain%specel(n)%Control_nodes)
            do i = 0, nnodes - 1
                ipoint = Tdomain%specel(n)%Control_nodes(i)
                point_count(ipoint) = point_count(ipoint) + 1
            end do
        end do
        allocate(point_to_elems(0:Tdomain%n_glob_nodes-1))
        do p = 0, Tdomain%n_glob_nodes - 1
            allocate(point_to_elems(p)%idx(point_count(p)))
        end do
        point_count = 0
        do n = 0, Tdomain%n_elem - 1
            nnodes = size(Tdomain%specel(n)%Control_nodes)
            do i = 0, nnodes - 1
                ipoint = Tdomain%specel(n)%Control_nodes(i)
                point_count(ipoint) = point_count(ipoint) + 1
                point_to_elems(ipoint)%idx(point_count(ipoint)) = n
            end do
        end do
        deallocate(point_count)
    end subroutine build_point_to_elems

    subroutine free_point_to_elems(point_to_elems, n_glob_nodes)
        implicit none
        type(idx_list_t), dimension(:), allocatable, intent(inout) :: point_to_elems
        integer, intent(in) :: n_glob_nodes
        integer :: p
        do p = 0, n_glob_nodes - 1
            deallocate(point_to_elems(p)%idx)
        end do
        deallocate(point_to_elems)
    end subroutine free_point_to_elems

    !> Growth of a per-element order array: propagates target_m =
    !! elem_order(n)-1 to every element sharing a corner control-node with n
    !! (via build_point_to_elems) within this rank, so a
    !! shrinking-halo iteration is never handed a neighbour more than one
    !! order below what it needs. Also propagates across rank boundaries:
    !! each partition-boundary sdom/fdom DOF's order requirement (max over
    !! the local elements touching it) is exchanged via the same
    !! Comm_data%IGiveS/IGiveF index lists comm_forces (Newmark.f90) uses
    !! for Veloc/ForcesFl, merged by MAX (comm_take_data_max_1) instead of
    !! SUM since every rank owning that DOF must agree on the largest order
    !! asked for, not their total. Alternates local BFS and boundary
    !! exchange (mirrors SEM2D's resolver_halos_from_orders) until neither
    !! finds anything left to raise, globally (MPI_Allreduce(LOR)).
    subroutine resolver_halos_from_orders_3d(Tdomain, elem_order)
        use sdomain
        use mpi
        use scomm
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, dimension(0:Tdomain%n_elem-1), intent(inout) :: elem_order

        type(idx_list_t), dimension(:), allocatable :: point_to_elems
        integer, dimension(:), allocatable :: queue
        integer :: qhead, qtail, qsize, n, i, j, ipoint, other, m_elem, target_m, nnodes
        real(fpp), dimension(:), allocatable :: dof_order_sdom, dof_order_fdom
        integer :: bnum, ee, ix, iy, iz, ngll, idx, k, ierr
        logical :: local_work, global_work

        call build_point_to_elems(Tdomain, point_to_elems)

        qsize = max(Tdomain%n_elem * 10, 100) ! generous size for reactivation
        allocate(queue(0:qsize-1))
        qhead = 0; qtail = 0
        do n = 0, Tdomain%n_elem - 1
            if (elem_order(n) > 1) then
                queue(qtail) = n
                qtail = qtail + 1
            end if
        end do

        if (Tdomain%sdom%nglltot > 0) allocate(dof_order_sdom(0:Tdomain%sdom%nglltot-1))
        if (Tdomain%fdom%nglltot > 0) allocate(dof_order_fdom(0:Tdomain%fdom%nglltot-1))

        global_work = .true.
        do while (global_work)
            ! --- Local BFS over shared corner control-nodes ---
            do while (qhead < qtail)
                n = queue(qhead)
                qhead = qhead + 1
                m_elem = elem_order(n)
                target_m = m_elem - 1
                if (target_m <= 1) cycle
                nnodes = size(Tdomain%specel(n)%Control_nodes)
                do i = 0, nnodes - 1
                    ipoint = Tdomain%specel(n)%Control_nodes(i)
                    do j = 1, size(point_to_elems(ipoint)%idx)
                        other = point_to_elems(ipoint)%idx(j)
                        if (elem_order(other) < target_m) then
                            elem_order(other) = target_m
                            if (qtail < qsize) then
                                queue(qtail) = other
                                qtail = qtail + 1
                            end if
                        end if
                    end do
                end do
            end do
            qhead = 0; qtail = 0
            local_work = .false.

            ! --- MPI boundary exchange & reactivation ---
            if (Tdomain%Comm_data%ncomm > 0) then
                if (allocated(dof_order_sdom)) then
                    dof_order_sdom = 1._fpp
                    ngll = Tdomain%sdom%ngll
                    do n = 0, Tdomain%n_elem - 1
                        if (Tdomain%specel(n)%domain /= DM_SOLID_CG) cycle
                        bnum = Tdomain%specel(n)%lnum / VCHUNK
                        ee = mod(Tdomain%specel(n)%lnum, VCHUNK)
                        do iz = 0, ngll-1
                            do iy = 0, ngll-1
                                do ix = 0, ngll-1
                                    idx = Tdomain%sdom%Idom_(ix,iy,iz,bnum,ee)
                                    dof_order_sdom(idx) = max(dof_order_sdom(idx), real(elem_order(n), fpp))
                                end do
                            end do
                        end do
                    end do
                end if
                if (allocated(dof_order_fdom)) then
                    dof_order_fdom = 1._fpp
                    ngll = Tdomain%fdom%ngll
                    do n = 0, Tdomain%n_elem - 1
                        if (Tdomain%specel(n)%domain /= DM_FLUID_CG) cycle
                        bnum = Tdomain%specel(n)%lnum / VCHUNK
                        ee = mod(Tdomain%specel(n)%lnum, VCHUNK)
                        do iz = 0, ngll-1
                            do iy = 0, ngll-1
                                do ix = 0, ngll-1
                                    idx = Tdomain%fdom%Idom_(ix,iy,iz,bnum,ee)
                                    dof_order_fdom(idx) = max(dof_order_fdom(idx), real(elem_order(n), fpp))
                                end do
                            end do
                        end do
                    end do
                end if

                do n = 0, Tdomain%Comm_data%ncomm - 1
                    k = 0
                    if (allocated(dof_order_sdom)) &
                        call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                            Tdomain%Comm_data%Data(n)%IGiveS, dof_order_sdom, k)
                    if (allocated(dof_order_fdom)) &
                        call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                            Tdomain%Comm_data%Data(n)%IGiveF, dof_order_fdom, k)
                    Tdomain%Comm_data%Data(n)%nsend = k
                end do

                call exchange_sem_var(Tdomain, 870, Tdomain%Comm_data)

                do n = 0, Tdomain%Comm_data%ncomm - 1
                    k = 0
                    if (allocated(dof_order_sdom)) &
                        call comm_take_data_max_1(Tdomain%Comm_data%Data(n)%Take, &
                            Tdomain%Comm_data%Data(n)%IGiveS, dof_order_sdom, k)
                    if (allocated(dof_order_fdom)) &
                        call comm_take_data_max_1(Tdomain%Comm_data%Data(n)%Take, &
                            Tdomain%Comm_data%Data(n)%IGiveF, dof_order_fdom, k)
                end do

                ! Raise local elements touching a boundary DOF whose merged
                ! order requires one more hop out from it than they have.
                if (allocated(dof_order_sdom)) then
                    ngll = Tdomain%sdom%ngll
                    do n = 0, Tdomain%n_elem - 1
                        if (Tdomain%specel(n)%domain /= DM_SOLID_CG) cycle
                        bnum = Tdomain%specel(n)%lnum / VCHUNK
                        ee = mod(Tdomain%specel(n)%lnum, VCHUNK)
                        target_m = 0
                        do iz = 0, ngll-1
                            do iy = 0, ngll-1
                                do ix = 0, ngll-1
                                    idx = Tdomain%sdom%Idom_(ix,iy,iz,bnum,ee)
                                    target_m = max(target_m, int(dof_order_sdom(idx)) - 1)
                                end do
                            end do
                        end do
                        if (target_m > elem_order(n)) then
                            elem_order(n) = target_m
                            local_work = .true.
                            if (qtail < qsize) then
                                queue(qtail) = n
                                qtail = qtail + 1
                            end if
                        end if
                    end do
                end if
                if (allocated(dof_order_fdom)) then
                    ngll = Tdomain%fdom%ngll
                    do n = 0, Tdomain%n_elem - 1
                        if (Tdomain%specel(n)%domain /= DM_FLUID_CG) cycle
                        bnum = Tdomain%specel(n)%lnum / VCHUNK
                        ee = mod(Tdomain%specel(n)%lnum, VCHUNK)
                        target_m = 0
                        do iz = 0, ngll-1
                            do iy = 0, ngll-1
                                do ix = 0, ngll-1
                                    idx = Tdomain%fdom%Idom_(ix,iy,iz,bnum,ee)
                                    target_m = max(target_m, int(dof_order_fdom(idx)) - 1)
                                end do
                            end do
                        end do
                        if (target_m > elem_order(n)) then
                            elem_order(n) = target_m
                            local_work = .true.
                            if (qtail < qsize) then
                                queue(qtail) = n
                                qtail = qtail + 1
                            end if
                        end if
                    end do
                end if
            end if

            local_work = local_work .or. (qhead < qtail)
            call MPI_Allreduce(local_work, global_work, 1, MPI_LOGICAL, MPI_LOR, Tdomain%communicateur, ierr)
        end do

        if (allocated(dof_order_sdom)) deallocate(dof_order_sdom)
        if (allocated(dof_order_fdom)) deallocate(dof_order_fdom)
        deallocate(queue)
        call free_point_to_elems(point_to_elems, Tdomain%n_glob_nodes)
    end subroutine resolver_halos_from_orders_3d

    !> Speedup & Cost Optimization: evaluates candidate target time steps
    !! dt_target (global percentiles of dt_elem_loc), determines minimum
    !! required element orders m_i, resolves halos via
    !! resolver_halos_from_orders_3d, and selects the (dt_target, {m_i})
    !! configuration that minimizes total computational cost
    !! C = (1 / dt_target) * sum(m_i). Ported from SEM2D/SRC/irons_dtcrit.F90.
    subroutine optimize_cost_and_orders(Tdomain, dt_elem_loc, dt_target_opt, elem_order_opt)
        use sdomain
        use mpi
        implicit none

        type(domain), intent(inout) :: Tdomain
        real(fpp), dimension(0:Tdomain%n_elem-1), intent(in) :: dt_elem_loc
        real(fpp), intent(out) :: dt_target_opt
        integer, dimension(0:Tdomain%n_elem-1), intent(out) :: elem_order_opt

        integer :: n, i, c, ierr, rg, n_procs, n_local, n_global
        integer :: cand_idx, feasible_loc, feasible_g
        real(fpp) :: local_cost, global_cost, total_cost, min_cost
        real(fpp) :: dt_target_cand
        logical :: cand_feasible

        integer, dimension(:), allocatable :: recv_counts, displs
        real(fpp), dimension(:), allocatable :: dt_global, dt_candidates
        integer, dimension(:), allocatable :: m_base_local, m_work_local

        integer, parameter :: NUM_PERCENTILES = 10
        real(fpp), dimension(NUM_PERCENTILES), parameter :: PERCENTILES = &
            (/ 0.00_fpp, 0.02_fpp, 0.05_fpp, 0.10_fpp, 0.20_fpp, 0.30_fpp, 0.40_fpp, 0.50_fpp, 0.70_fpp, 0.90_fpp /)

        rg = Tdomain%rank
        n_procs = Tdomain%nb_procs
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
            call resolver_halos_from_orders_3d(Tdomain, m_work_local)

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

    !> For each domain (sdom, fdom) and each iteration j=1..mm, build the
    !! list of block indices touched (>=1 element with hop<=mm-j maps into
    !! that block) -- the force kernels (calcul_forces_iso/aniso,
    !! calcul_forces_fluid) only ever run on a whole block, so this is the
    !! coarsest granularity "regional" activation can use for them. Also
    !! builds dof_corrected: each DOF's own correction order, MIN(elem_order)
    !! over every element touching it (via Idom_) -- 1 (excluded) for any DOF
    !! shared with an order-1 element or a solid-fluid interface (those
    !! coexist untouched, same as PML/CPML/DG which are never %modified in
    !! the first place since compute_irons_dtcrit only fills dt_elem_loc for
    !! sdom/fdom).
    subroutine derive_block_and_dof_activation(Tdomain, mm, elem_order)
        use sdomain
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: mm
        integer, dimension(0:Tdomain%n_elem-1), intent(in) :: elem_order
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

        ! dof_corrected: per-DOF correction order = MIN(elem_order) over
        ! every sdom/fdom element touching that DOF -- generalizes the old
        ! true/false gate ("every owning element modified, or none") to a
        ! per-DOF order ceiling. A DOF shared with any order-1 element is
        ! capped at 1 (no k=2.. term applied there, same as the old
        ! "false"); a DOF whose every owner shares the same order m reduces
        ! to the old "true" (gets every k=2..m term).
        if (allocated(dof_corrected_sdom)) deallocate(dof_corrected_sdom)
        if (allocated(dof_corrected_fdom)) deallocate(dof_corrected_fdom)
        allocate(dof_corrected_sdom(0:Tdomain%sdom%nglltot-1))
        allocate(dof_corrected_fdom(0:Tdomain%fdom%nglltot-1))
        dof_corrected_sdom = mm
        dof_corrected_fdom = mm

        do n = 0, Tdomain%n_elem - 1
            dom = Tdomain%specel(n)%domain
            if (dom /= DM_SOLID_CG .and. dom /= DM_FLUID_CG) cycle
            bnum = Tdomain%specel(n)%lnum / VCHUNK
            ee = mod(Tdomain%specel(n)%lnum, VCHUNK)
            if (dom == DM_SOLID_CG) then
                do k = 0, Tdomain%sdom%ngll - 1
                    do jy = 0, Tdomain%sdom%ngll - 1
                        do i = 0, Tdomain%sdom%ngll - 1
                            idx = Tdomain%sdom%Idom_(i,jy,k,bnum,ee)
                            dof_corrected_sdom(idx) = min(dof_corrected_sdom(idx), elem_order(n))
                        end do
                    end do
                end do
            else
                do k = 0, Tdomain%fdom%ngll - 1
                    do jy = 0, Tdomain%fdom%ngll - 1
                        do i = 0, Tdomain%fdom%ngll - 1
                            idx = Tdomain%fdom%Idom_(i,jy,k,bnum,ee)
                            dof_corrected_fdom(idx) = min(dof_corrected_fdom(idx), elem_order(n))
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
            dof_corrected_sdom(Tdomain%SF%intSolFlu%surf0%map(idx)) = 1
        end do
        do idx = 0, Tdomain%SF%intSolFlu%surf1%nbtot - 1
            dof_corrected_fdom(Tdomain%SF%intSolFlu%surf1%map(idx)) = 1
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
