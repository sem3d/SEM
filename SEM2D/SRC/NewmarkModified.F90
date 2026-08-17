!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file NewmarkModified.F90
!!\brief Matrix-free "modified equation" Newmark (Meddeb, "Seismic design of
!! dams", 2025, eq. 1.28-1.30), a NEW time-integration routine - Newmark.F90
!! is untouched. Only supports meshes with no PML/CPML, no fault, and no
!! solid-fluid interface (enforced at startup by check_modified_newmark in
!! check_inputs_and_mesh.F90).
!!
!! Derivation (see also App/@WaveMatrixFreeClass and App/@StabilityClass in
!! the labcorrea/FEM MATLAB prototype, ModifiedEquationZCrit in
!! irons_dtcrit.F90 for the matching critical-dt bound):
!!   depl(n+1) = 2*depl(n) - depl(n-1) + dt^2*Minv*(Fext-K*depl(n))
!!             + sum_{k=2}^m c_k * A^k * depl(n),   A = -Minv*K,
!!             c_k = 2*dt^(2k)/(2k)!
!! The first line is EXACTLY what the existing velocity-form Newmark(beta=0)
!! already computes in one call (see WaveMatrixFreeClass.leapfrog / the m=1
!! reduction documented in SolverClass.NewmarkModified in the MATLAB repo).
!! So this routine: (1) snapshots depl(n) for every element/face/vertex,
!! (2) calls the UNMODIFIED classic Newmark for the base step, (3) adds the
!! k=2..m correction terms computed via k-1 extra matrix-free applications
!! of A (reusing compute_InternalForces_Elem/Fl, getInternalF_el2f and the
!! existing MPI exchange pattern - never assembling a global K), and
!! (4) resyncs Veloc = (Displ_new - depl(n))/dt so the NEXT step's classic
!! Newmark call sees a consistent (Displ,Veloc) pair again.
!<

module snewmark_modified
    use sdomain
    use mpi
    use snewmark, only : Newmark
    implicit none

    type :: elem_scratch_t
        real(fpp), dimension(:,:,:), allocatable :: u
    end type elem_scratch_t
    type :: face_scratch_t
        real(fpp), dimension(:,:), allocatable :: u
    end type face_scratch_t
    type :: idx_list_t
        integer, dimension(:), allocatable :: idx
    end type idx_list_t

    ! Snapshot of depl(n) (elem_un/face_un/vert_un) and the accumulated
    ! k=2..m correction (elem_corr/face_corr/vert_corr), one entry per
    ! element/face; vertices are fixed-size (0:1) so a plain array suffices.
    type(elem_scratch_t), dimension(:), allocatable, save :: elem_un, elem_corr
    type(face_scratch_t), dimension(:), allocatable, save :: face_un, face_corr
    real(fpp), dimension(:,:), allocatable, save :: vert_un, vert_corr
    logical, save :: scratch_ready = .false.

    ! Regional halo: elem_hop(n) = graph-distance (face-adjacency
    ! hops) from element n to the nearest %modified element, -1 if beyond
    ! modified_order-1 (never touched by any apply_a_global iteration).
    ! elem_active(j)/face_active(j)/vert_active(j), j=1..modified_order, are
    ! the index lists apply_a_global visits at iteration j: exactly
    ! {hop <= modified_order-j}. This SHRINKS as j grows -- required for
    ! correctness, not just an optimization: computing A*w correctly at
    ! hop<=d needs w correct at hop<=d+1 (one application of A only mixes
    ! face/vertex neighbours), so iteration j must recompute exactly one hop
    ! narrower than what iteration j-1 left correctly computed. A fixed-width
    ! halo reused at every iteration would re-touch the outer layer using a
    ! neighbour that was never (correctly) updated, silently corrupting the
    ! value the next iteration needs. elem_active(modified_order) is exactly
    ! the %modified set (hop=0) -- also used to restrict the
    ! correction-accumulation and final-commit loops in NewmarkModified.
    integer, dimension(:), allocatable, save :: elem_hop
    type(idx_list_t), dimension(:), allocatable, save :: elem_active, face_active, vert_active

contains

    subroutine NewmarkModified(Tdomain)
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: n, k, m, idx
        real(fpp) :: dt, ck

        if (.not. Tdomain%TimeD%velocity_scheme) return

        call ensure_scratch(Tdomain)

        m = Tdomain%TimeD%modified_order
        dt = Tdomain%TimeD%dtmin

        ! Snapshot depl(n) before the base step overwrites Displ/Veloc.
        do n = 0, Tdomain%n_elem - 1
            elem_un(n)%u = Tdomain%specel(n)%Displ
        enddo
        do n = 0, Tdomain%n_face - 1
            face_un(n)%u = Tdomain%sFace(n)%Displ
        enddo
        do n = 0, Tdomain%n_vertex - 1
            vert_un(n,0:1) = Tdomain%sVertex(n)%Displ
        enddo

        ! Base order-1 (plain leapfrog) step: reuse the existing, unmodified
        ! classic Newmark - this is exactly the k=1 term of the recursion.
        call Newmark(Tdomain)

        if (m < 2) return

        ! k=2..m correction: seed w:=depl(n) into the Forces scratch, then
        ! chain A*w -> A^2*w -> ... via apply_a_global, accumulating
        ! sum_{k=2}^m c_k*A^k*depl(n).
        call seed_forces_from_un(Tdomain)

        ! Only elem_active(2)/face_active(2)/vert_active(2) (elements with
        ! order >= 2) ever get read back out of elem_corr/face_corr/vert_corr
        ! in the commit step below, so zero those.
        if (allocated(elem_active) .and. size(elem_active) >= 2) then
            do idx = 1, size(elem_active(2)%idx)
                elem_corr(elem_active(2)%idx(idx))%u = 0._fpp
            enddo
            do idx = 1, size(face_active(2)%idx)
                face_corr(face_active(2)%idx(idx))%u = 0._fpp
            enddo
            do idx = 1, size(vert_active(2)%idx)
                vert_corr(vert_active(2)%idx(idx),0:1) = 0._fpp
            enddo
        endif

        do k = 1, m
            call apply_a_global(Tdomain, k, m)
            if (k >= 2) then
                ck = 2._fpp*dt**(2*k)/fact2k(k)
                do idx = 1, size(elem_active(k)%idx)
                    n = elem_active(k)%idx(idx)
                    elem_corr(n)%u = elem_corr(n)%u + ck*elem_forces_interior(Tdomain,n)
                enddo
                do idx = 1, size(face_active(k)%idx)
                    n = face_active(k)%idx(idx)
                    face_corr(n)%u = face_corr(n)%u + ck*Tdomain%sFace(n)%Forces
                enddo
                do idx = 1, size(vert_active(k)%idx)
                    n = vert_active(k)%idx(idx)
                    vert_corr(n,0:1) = vert_corr(n,0:1) + ck*Tdomain%sVertex(n)%Forces(0:1)
                enddo
            endif
        enddo

        ! Apply the correction to elements/faces/vertices with order >= 2
        ! and resync Veloc = (Displ_new - depl(n))/dt EVERYWHERE so the next step's
        ! classic-Newmark base call sees a consistent pair.
        if (allocated(elem_active) .and. size(elem_active) >= 2) then
            do idx = 1, size(elem_active(2)%idx)
                n = elem_active(2)%idx(idx)
                Tdomain%specel(n)%Displ = Tdomain%specel(n)%Displ + elem_corr(n)%u
            enddo
        endif
        do n = 0, Tdomain%n_elem - 1
            Tdomain%specel(n)%Veloc = (Tdomain%specel(n)%Displ - elem_un(n)%u)/dt
        enddo
        if (allocated(face_active) .and. size(face_active) >= 2) then
            do idx = 1, size(face_active(2)%idx)
                n = face_active(2)%idx(idx)
                Tdomain%sFace(n)%Displ = Tdomain%sFace(n)%Displ + face_corr(n)%u
            enddo
        endif
        do n = 0, Tdomain%n_face - 1
            Tdomain%sFace(n)%Veloc = (Tdomain%sFace(n)%Displ - face_un(n)%u)/dt
        enddo
        if (allocated(vert_active) .and. size(vert_active) >= 2) then
            do idx = 1, size(vert_active(2)%idx)
                n = vert_active(2)%idx(idx)
                Tdomain%sVertex(n)%Displ = Tdomain%sVertex(n)%Displ + vert_corr(n,0:1)
            enddo
        endif
        do n = 0, Tdomain%n_vertex - 1
            Tdomain%sVertex(n)%Veloc = (Tdomain%sVertex(n)%Displ - vert_un(n,0:1))/dt
        enddo
    end subroutine NewmarkModified

    !> Allocate the module-level scratch on first use (after allocate_domain
    !! has sized Displ/Veloc/Forces, so mold= copies the right shapes).
    subroutine ensure_scratch(Tdomain)
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: n

        if (scratch_ready) return

        allocate(elem_un(0:Tdomain%n_elem-1))
        allocate(elem_corr(0:Tdomain%n_elem-1))
        do n = 0, Tdomain%n_elem - 1
            allocate(elem_un(n)%u, mold=Tdomain%specel(n)%Displ)
            allocate(elem_corr(n)%u, mold=Tdomain%specel(n)%Displ)
        enddo

        allocate(face_un(0:Tdomain%n_face-1))
        allocate(face_corr(0:Tdomain%n_face-1))
        do n = 0, Tdomain%n_face - 1
            allocate(face_un(n)%u, mold=Tdomain%sFace(n)%Displ)
            allocate(face_corr(n)%u, mold=Tdomain%sFace(n)%Displ)
        enddo

        allocate(vert_un(0:Tdomain%n_vertex-1,0:1))
        allocate(vert_corr(0:Tdomain%n_vertex-1,0:1))

        ! elem_active/face_active/vert_active/elem_hop (graded per-element
        ! orders) are already built by compute_irons_dtcrit's call to
        ! build_regional_halo_with_orders (irons_dtcrit.F90), which always
        ! runs before the time loop starts whenever Tdomain%TimeD%modified is
        ! true -- the only condition under which NewmarkModified (and hence
        ! this routine) is ever called (main.F90). Re-deriving it here from
        ! Tdomain%specel(n)%modified would flatten every modified element to
        ! the same uniform max order, discarding the cost-optimizer's
        ! per-element grading -- do not rebuild it.

        scratch_ready = .true.
    end subroutine ensure_scratch

    !> Precompute elem_order, elem_hop and per-iteration index lists using pre-assigned
    !! per-element orders, executing Distributed Multi-Source BFS across MPI partition boundaries.
    subroutine build_regional_halo_with_orders(Tdomain, elem_order_in)
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, dimension(0:Tdomain%n_elem-1), intent(in) :: elem_order_in

        integer, dimension(:), allocatable :: elem_order, face_order, vert_order
        type(idx_list_t), dimension(:), allocatable :: vert_to_elems
        integer :: n, k, nf, nv, e0, e1, j, req_m
        integer :: max_m_local, max_m_global, ierr

        allocate(elem_order(0:Tdomain%n_elem-1))
        elem_order = elem_order_in

        ! Resolve halos via Distributed BFS
        call resolver_halos_from_orders(Tdomain, elem_order)

        ! Determine global maximum order
        max_m_local = maxval(elem_order)
        if (Tdomain%Mpi_var%n_proc > 1) then
            call MPI_ALLREDUCE(max_m_local, max_m_global, 1, MPI_INTEGER, MPI_MAX, Tdomain%communicateur, ierr)
        else
            max_m_global = max_m_local
        end if

        ! Update Tdomain%TimeD%modified_order
        Tdomain%TimeD%modified_order = max_m_global

        ! Flag elements as modified if order > 1
        do n = 0, Tdomain%n_elem - 1
            Tdomain%specel(n)%modified = (elem_order(n) > 1)
        end do

        ! Build vert_to_elems mapping for derived face/vert orders
        allocate(vert_to_elems(0:Tdomain%n_vertex-1))
        allocate(vert_order(0:Tdomain%n_vertex-1))
        vert_order = 0

        do n = 0, Tdomain%n_elem - 1
            do k = 0, 3
                nv = Tdomain%specel(n)%Near_Vertex(k)
                vert_order(nv) = vert_order(nv) + 1
            end do
        end do
        do nv = 0, Tdomain%n_vertex - 1
            allocate(vert_to_elems(nv)%idx(vert_order(nv)))
        end do
        vert_order = 0
        do n = 0, Tdomain%n_elem - 1
            do k = 0, 3
                nv = Tdomain%specel(n)%Near_Vertex(k)
                vert_order(nv) = vert_order(nv) + 1
                vert_to_elems(nv)%idx(vert_order(nv)) = n
            end do
        end do

        ! Derived face_order and vert_order
        allocate(face_order(0:Tdomain%n_face-1))
        face_order = 1
        vert_order = 1

        do nf = 0, Tdomain%n_face - 1
            e0 = Tdomain%sFace(nf)%Near_Element(0)
            e1 = Tdomain%sFace(nf)%Near_Element(1)
            req_m = elem_order(e0)
            if (e1 > -1) req_m = max(req_m, elem_order(e1))
            face_order(nf) = req_m
        end do

        do nv = 0, Tdomain%n_vertex - 1
            req_m = 1
            do j = 1, size(vert_to_elems(nv)%idx)
                e0 = vert_to_elems(nv)%idx(j)
                req_m = max(req_m, elem_order(e0))
            end do
            vert_order(nv) = req_m
        end do

        ! Store elem_hop for backward compatibility with apply_a_global
        if (allocated(elem_hop)) deallocate(elem_hop)
        allocate(elem_hop(0:Tdomain%n_elem-1))
        do n = 0, Tdomain%n_elem - 1
            elem_hop(n) = max_m_global - elem_order(n)
        end do

        ! Pack active lists for each iteration j = 1 .. max_m_global
        if (allocated(elem_active)) deallocate(elem_active)
        if (allocated(face_active)) deallocate(face_active)
        if (allocated(vert_active)) deallocate(vert_active)
        allocate(elem_active(max_m_global), face_active(max_m_global), vert_active(max_m_global))
        do j = 1, max_m_global
            call pack_indices_order(elem_order, j, elem_active(j)%idx)
            call pack_indices_order(face_order, j, face_active(j)%idx)
            call pack_indices_order(vert_order, j, vert_active(j)%idx)
        end do

        deallocate(face_order, vert_order, elem_order)
        do nv = 0, Tdomain%n_vertex - 1
            deallocate(vert_to_elems(nv)%idx)
        end do
        deallocate(vert_to_elems)
    end subroutine build_regional_halo_with_orders

    !> Distributed Multi-Source BFS: Resolves order halos (m_nbr >= m - 1)
    !! across local mesh and MPI partition boundaries.
    subroutine resolver_halos_from_orders(Tdomain, elem_order)
        use mpi
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, dimension(0:Tdomain%n_elem-1), intent(inout) :: elem_order

        integer, dimension(:), allocatable :: queue
        type(idx_list_t), dimension(:), allocatable :: vert_to_elems
        integer, dimension(:), allocatable :: vert_order
        integer :: qhead, qtail, qsize, n, k, nf, nv, e0, e1, j, m_elem, target_m
        integer :: ierr, i_proc, i_send, req_m
        integer :: nv_aus, n_face_pointed, n_faces, n_verts
        logical :: local_work, global_work, updated

        ! Shared buffers for MPI exchange of face and vertex orders
        integer, dimension(:), allocatable :: send_buf, recv_buf
        integer :: tag_send, tag_recv
        integer, dimension(MPI_STATUS_SIZE) :: status

        ! Build vert_to_elems mapping for quick node-neighbor queries
        allocate(vert_to_elems(0:Tdomain%n_vertex-1))
        allocate(vert_order(0:Tdomain%n_vertex-1))
        vert_order = 0

        ! Count elements per vertex
        do n = 0, Tdomain%n_elem - 1
            do k = 0, 3
                nv = Tdomain%specel(n)%Near_Vertex(k)
                vert_order(nv) = vert_order(nv) + 1
            end do
        end do
        do nv = 0, Tdomain%n_vertex - 1
            allocate(vert_to_elems(nv)%idx(vert_order(nv)))
        end do
        vert_order = 0
        do n = 0, Tdomain%n_elem - 1
            do k = 0, 3
                nv = Tdomain%specel(n)%Near_Vertex(k)
                vert_order(nv) = vert_order(nv) + 1
                vert_to_elems(nv)%idx(vert_order(nv)) = n
            end do
        end do
        deallocate(vert_order)

        ! 2. Initialize BFS Queue
        qsize = max(Tdomain%n_elem * 10, 100) ! generous size for reactivation
        allocate(queue(0:qsize-1))
        qhead = 0; qtail = 0

        do n = 0, Tdomain%n_elem - 1
            if (elem_order(n) > 1) then
                queue(qtail) = n
                qtail = qtail + 1
            end if
        end do

        ! 3. Outer Loop: Alternates Local BFS Queue and MPI Boundary Exchange
        global_work = .true.
        do while (global_work)
            updated = .false.

            ! --- Local BFS Queue Processing ---
            do while (qhead < qtail)
                n = queue(qhead)
                qhead = qhead + 1
                m_elem = elem_order(n)
                target_m = m_elem - 1
                if (target_m <= 1) cycle

                ! Propagate to all elements sharing any vertex with element n
                do k = 0, 3
                    nv = Tdomain%specel(n)%Near_Vertex(k)
                    do j = 1, size(vert_to_elems(nv)%idx)
                        e0 = vert_to_elems(nv)%idx(j)
                        if (elem_order(e0) < target_m) then
                            elem_order(e0) = target_m
                            updated = .true.
                            if (qtail < qsize) then
                                queue(qtail) = e0
                                qtail = qtail + 1
                            end if
                        end if
                    end do
                end do
            end do

            ! Reset queue pointers once drained
            qhead = 0
            qtail = 0

            ! --- MPI Boundary Exchange & Reactivation ---
            if (Tdomain%Mpi_var%n_proc > 1) then
                do i_proc = 0, Tdomain%n_communications - 1
                    i_send = Tdomain%Communication_list(i_proc)
                    n_faces = Tdomain%sWall(i_proc)%n_faces
                    n_verts = Tdomain%sWall(i_proc)%n_vertices

                    if (n_faces + n_verts == 0) cycle

                    allocate(send_buf(n_faces + n_verts))
                    allocate(recv_buf(n_faces + n_verts))

                    ! Pack face max orders
                    do nf = 0, n_faces - 1
                        n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
                        e0 = Tdomain%sFace(n_face_pointed)%Near_Element(0)
                        e1 = Tdomain%sFace(n_face_pointed)%Near_Element(1)
                        req_m = elem_order(e0)
                        if (e1 > -1) req_m = max(req_m, elem_order(e1))
                        send_buf(nf + 1) = req_m
                    end do

                    ! Pack vertex max orders
                    do nv = 0, n_verts - 1
                        nv_aus = Tdomain%sWall(i_proc)%Vertex_List(nv)
                        req_m = 1
                        do j = 1, size(vert_to_elems(nv_aus)%idx)
                            e0 = vert_to_elems(nv_aus)%idx(j)
                            req_m = max(req_m, elem_order(e0))
                        end do
                        send_buf(n_faces + nv + 1) = req_m
                    end do

                    tag_send = i_send * Tdomain%MPI_var%n_proc + Tdomain%MPI_var%my_rank + 850
                    tag_recv = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 850

                    call MPI_SENDRECV(send_buf, n_faces + n_verts, MPI_INTEGER, i_send, tag_send, &
                                      recv_buf, n_faces + n_verts, MPI_INTEGER, i_send, tag_recv, &
                                      Tdomain%communicateur, status, ierr)

                    ! Unpack face max orders & reactivate
                    do nf = 0, n_faces - 1
                        n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
                        target_m = recv_buf(nf + 1) - 1
                        if (target_m <= 1) cycle
                        e0 = Tdomain%sFace(n_face_pointed)%Near_Element(0)
                        if (elem_order(e0) < target_m) then
                            elem_order(e0) = target_m
                            updated = .true.
                            if (qtail < qsize) then
                                queue(qtail) = e0
                                qtail = qtail + 1
                            end if
                        end if
                        e1 = Tdomain%sFace(n_face_pointed)%Near_Element(1)
                        if (e1 > -1) then
                            if (elem_order(e1) < target_m) then
                                elem_order(e1) = target_m
                                updated = .true.
                                if (qtail < qsize) then
                                    queue(qtail) = e1
                                    qtail = qtail + 1
                                end if
                            end if
                        end if
                    end do

                    ! Unpack vertex max orders & reactivate
                    do nv = 0, n_verts - 1
                        nv_aus = Tdomain%sWall(i_proc)%Vertex_List(nv)
                        target_m = recv_buf(n_faces + nv + 1) - 1
                        if (target_m <= 1) cycle
                        do j = 1, size(vert_to_elems(nv_aus)%idx)
                            e0 = vert_to_elems(nv_aus)%idx(j)
                            if (elem_order(e0) < target_m) then
                                elem_order(e0) = target_m
                                updated = .true.
                                if (qtail < qsize) then
                                    queue(qtail) = e0
                                    qtail = qtail + 1
                                end if
                            end if
                        end do
                    end do

                    deallocate(send_buf, recv_buf)
                end do
            end if

            ! --- Check Global Consensus ---
            local_work = updated .or. (qhead < qtail)
            if (Tdomain%Mpi_var%n_proc > 1) then
                call MPI_ALLREDUCE(local_work, global_work, 1, MPI_LOGICAL, MPI_LOR, Tdomain%communicateur, ierr)
            else
                global_work = local_work
            end if
        end do

        deallocate(queue)
        do nv = 0, Tdomain%n_vertex - 1
            deallocate(vert_to_elems(nv)%idx)
        end do
        deallocate(vert_to_elems)
    end subroutine resolver_halos_from_orders

    !> Indices n (0-based) with order(n) >= target_order.
    subroutine pack_indices_order(order_arr, target_order, idx)
        implicit none
        integer, dimension(0:), intent(in) :: order_arr
        integer, intent(in) :: target_order
        integer, dimension(:), allocatable, intent(out) :: idx
        integer :: n, cnt

        cnt = 0
        do n = 0, size(order_arr) - 1
            if (order_arr(n) >= target_order) cnt = cnt + 1
        end do
        allocate(idx(cnt))
        cnt = 0
        do n = 0, size(order_arr) - 1
            if (order_arr(n) >= target_order) then
                cnt = cnt + 1
                idx(cnt) = n
            end if
        end do
    end subroutine pack_indices_order

    !> Copy the snapshot depl(n) into Forces (element interior + face +
    !! vertex) - mirrors Prediction_Elem/Face/Vertex_Veloc's "Forces:=Displ"
    !! line, but from the saved snapshot instead of the (already advanced by
    !! the base Newmark step) live Displ.
    subroutine seed_forces_from_un(Tdomain)
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: n, ngllx, ngllz

        do n = 0, Tdomain%n_elem - 1
            ngllx = Tdomain%specel(n)%ngllx; ngllz = Tdomain%specel(n)%ngllz
            Tdomain%specel(n)%Forces(1:ngllx-2,1:ngllz-2,0:1) = elem_un(n)%u
        enddo
        do n = 0, Tdomain%n_face - 1
            Tdomain%sFace(n)%Forces = face_un(n)%u
        enddo
        do n = 0, Tdomain%n_vertex - 1
            Tdomain%sVertex(n)%Forces(0:1) = vert_un(n,0:1)
        enddo
    end subroutine seed_forces_from_un

    !> Element Forces' interior slice (helper so the sizes documented above
    !! line up with elem_corr(n)%u, which is mold=Displ i.e. interior-only).
    function elem_forces_interior(Tdomain, n) result(f)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: n
        real(fpp), dimension(1:Tdomain%specel(n)%ngllx-2,1:Tdomain%specel(n)%ngllz-2,0:1) :: f
        f = Tdomain%specel(n)%Forces(1:Tdomain%specel(n)%ngllx-2,1:Tdomain%specel(n)%ngllz-2,0:1)
    end function elem_forces_interior

    !> A = -Minv*K, matrix-free: given w already sitting in Forces (element
    !! interior + face + vertex), overwrite Forces with A*w. No external
    !! forcing is added (Compute_external_forces is NOT called) - this is a
    !! pure linear operator application, chainable (call again to get A^2*w,
    !! etc.). Assumes no PML/CPML/fault/solid-fluid-interface elements
    !! exist in the mesh (guaranteed by check_modified_newmark).
    !!
    !! Regional: iteration jiter (1-based, out of mm total) only
    !! touches elem_active(jiter)/face_active(jiter)/vert_active(jiter) ==
    !! {hop <= mm-jiter} -- NOT the whole mesh. See the module-level comment
    !! on elem_hop for why this must shrink with jiter, not stay fixed.
    subroutine apply_a_global(Tdomain, jiter, mm)
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: jiter, mm

        integer :: n, nf, nv, nv_aus, nelem, w_face, ngllx, ngllz, ngll, i, j, mat, idx, dmax
        integer :: n_face_pointed, tag_send, tag_receive, i_send, i_stock, i_proc, ierr
        integer, dimension(MPI_STATUS_SIZE) :: status

        dmax = mm - jiter

        ! 1) Gather face/vertex w into element boundary rows (get_Displ_fv2el
        ! reads Face/Vertex%Forces, which already hold w from the caller's
        ! seed or the previous apply_a_global call), then apply the local
        ! element stiffness in place (Forces := -K_e*w, per-element).
        ! Restricted to elem_active(jiter): every one of its face/vertex
        ! neighbours has hop<=dmax+1, guaranteed correctly set by iteration
        ! jiter-1 (or by seed_forces_from_un, exact everywhere, for jiter=1).
        do idx = 1, size(elem_active(jiter)%idx)
            n = elem_active(jiter)%idx(idx)
            mat = Tdomain%specel(n)%mat_index
            call get_Displ_fv2el(Tdomain, n)
            if (Tdomain%specel(n)%acoustic) then
                call compute_InternalForcesFl_Elem(Tdomain%specel(n), &
                    Tdomain%sSubDomain(mat)%hprimex, Tdomain%sSubDomain(mat)%hTprimex, &
                    Tdomain%sSubDomain(mat)%hprimez, Tdomain%sSubDomain(mat)%hTprimez)
            else
                call compute_InternalForces_Elem(Tdomain%specel(n), &
                    Tdomain%sSubDomain(mat)%hprimex, Tdomain%sSubDomain(mat)%hTprimex, &
                    Tdomain%sSubDomain(mat)%hprimez, Tdomain%sSubDomain(mat)%hTprimez)
            endif
        enddo

        ! 2) Scatter/sum element boundary -K_e*w into faces and vertices
        ! (verbatim pattern from Newmark.F90's "Communication of Forces"),
        ! restricted to face_active(jiter). A face at the outer edge of this
        ! iteration's halo can have ONE neighbour in elem_active(jiter) and
        ! the other outside it (never recomputed this iteration, still
        ! holding a stale/wrong-type value) -- only scatter from a side
        ! whose element was actually recomputed here. The resulting
        ! "one-sided" face value is never read by iteration jiter+1 (that
        ! face's hop is one step beyond what iteration jiter+1 needs), so
        ! this is exact, not an approximation.
        do idx = 1, size(face_active(jiter)%idx)
            nf = face_active(jiter)%idx(idx)
            Tdomain%sFace(nf)%Forces = 0
            nelem = Tdomain%sFace(nf)%Near_element(0)
            if (elem_hop(nelem) >= 0 .and. elem_hop(nelem) <= dmax) then
                w_face = Tdomain%sFace(nf)%Which_face(0)
                call getInternalF_el2f(Tdomain, nelem, nf, w_face, .true.)
            endif
            nelem = Tdomain%sFace(nf)%Near_element(1)
            if (nelem > -1) then
                if (elem_hop(nelem) >= 0 .and. elem_hop(nelem) <= dmax) then
                    w_face = Tdomain%sFace(nf)%Which_face(1)
                    call getInternalF_el2f(Tdomain, nelem, nf, w_face, Tdomain%sFace(nf)%coherency)
                endif
            endif
        enddo

        do idx = 1, size(vert_active(jiter)%idx)
            Tdomain%sVertex(vert_active(jiter)%idx(idx))%Forces = 0
        enddo

        do idx = 1, size(elem_active(jiter)%idx)
            n = elem_active(jiter)%idx(idx)
            ngllx = Tdomain%specel(n)%ngllx; ngllz = Tdomain%specel(n)%ngllz
            nv = Tdomain%specel(n)%Near_Vertex(0)
            Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%Forces(0:1) + Tdomain%specel(n)%Forces(0,0,0:1)
            nv = Tdomain%specel(n)%Near_Vertex(1)
            Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%Forces(0:1) + Tdomain%specel(n)%Forces(ngllx-1,0,0:1)
            nv = Tdomain%specel(n)%Near_Vertex(2)
            Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%Forces(0:1) + Tdomain%specel(n)%Forces(ngllx-1,ngllz-1,0:1)
            nv = Tdomain%specel(n)%Near_Vertex(3)
            Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%Forces(0:1) + Tdomain%specel(n)%Forces(0,ngllz-1,0:1)
        enddo

        ! 3) MPI exchange: sum face/vertex Forces across ranks. Verbatim
        ! pattern from Newmark.F90 (including the PML-face/PML-vertex
        ! sections of the Send/Receive buffer, which are 0-iteration
        ! no-ops here since check_modified_newmark guarantees no PML
        ! elements/faces/vertices exist).
        do n = 0, Tdomain%n_communications - 1
            do nv = 0, Tdomain%sWall(n)%n_vertices-1
                nv_aus = Tdomain%sWall(n)%Vertex_List(nv)
                Tdomain%sVertex(nv_aus)%Double_Value(0:1) = Tdomain%sVertex(nv_aus)%Forces(0:1)
            enddo
        enddo

        do i_proc = 0, Tdomain%n_communications - 1
            i_send = Tdomain%Communication_list (i_proc)
            i_stock = 0

            do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
                n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
                ngll = Tdomain%sFace(n_face_pointed)%ngll
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    Tdomain%sWall(i_proc)%Send_data_2(i_stock:i_stock+ngll-3,0:1) = Tdomain%sFace(n_face_pointed)%Forces (1:ngll-2,0:1)
                else
                    do j = 1, ngll-2
                        Tdomain%sWall(i_proc)%Send_data_2(i_stock+j-1,0:1) = Tdomain%sFace(n_face_pointed)%Forces ( ngll-1-j,0:1)
                    enddo
                endif
                i_stock = i_stock + ngll - 2
            enddo

            do nv = 0, Tdomain%sWall(i_proc)%n_vertices - 1
                nv_aus =  Tdomain%sWall(i_proc)%Vertex_List(nv)
                Tdomain%sWall(i_proc)%Send_data_2 (i_stock,0:1)  = Tdomain%sVertex(nv_aus)%Double_Value(0:1)
                i_stock = i_stock + 1
            enddo

            tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank + 950
            tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 950

            ! MPI_SENDRECV (not separate SEND+RECV): deadlock-safe regardless of
            ! message size, matching the fix in Newmark.F90 (commit 3126710e,
            ! "[BUGFIX] deadlocks in MPI comm sem2D") - a separate blocking
            ! SEND+RECV pair deadlocks once Send_data_2 exceeds the MPI
            ! eager-message threshold (two mutual partners both send first).
            call MPI_SENDRECV (Tdomain%sWall(i_proc)%Send_data_2, 2*Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION, i_send, tag_send, &
                Tdomain%sWall(i_proc)%Receive_data_2, 2*Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION, i_send, tag_receive, &
                Tdomain%communicateur, status, ierr )

            i_stock = 0
            do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
                n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
                ngll = Tdomain%sFace(n_face_pointed)%ngll
                if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                    Tdomain%sFace(n_face_pointed)%Forces (1:ngll-2,0:1) =  Tdomain%sFace(n_face_pointed)%Forces ( 1:ngll-2,0:1)  + &
                        Tdomain%sWall(i_proc)%Receive_data_2(i_stock:i_stock+ngll-3,0:1)
                else
                    do j = 1, ngll-2
                        Tdomain%sFace(n_face_pointed)%Forces (ngll-1-j,0:1) =  Tdomain%sFace(n_face_pointed)%Forces (ngll-1-j,0:1)+ &
                            Tdomain%sWall(i_proc)%Receive_data_2(i_stock+j-1,0:1)
                    enddo
                endif
                i_stock = i_stock + ngll - 2
            enddo

            do nv = 0, Tdomain%sWall(i_proc)%n_vertices - 1
                nv_aus =  Tdomain%sWall(i_proc)%Vertex_List(nv)
                Tdomain%sVertex(nv_aus)%Forces(0:1) = Tdomain%sVertex(nv_aus)%Forces(0:1) + &
                    Tdomain%sWall(i_proc)%Receive_data_2 (i_stock,0:1)
                i_stock = i_stock + 1
            enddo
        enddo

        ! 4) A*w = Minv*(Forces), Forces so far == -K*w assembled (see
        ! irons_dtcrit.F90's note on compute_InternalForces_Elem's sign
        ! convention); MassMat is already 1/M (see define_arr.F90). Then
        ! reproduce the Abs/reflex zeroing that Correction_Face/Vertex_Veloc
        ! apply, so absorbing/reflecting boundaries see the same physics
        ! here as in the base (classic Newmark) step. Restricted to
        ! elem_active(jiter)/face_active(jiter)/vert_active(jiter) -- same
        ! halo as steps 1-2, nothing else was touched this iteration.
        do idx = 1, size(elem_active(jiter)%idx)
            n = elem_active(jiter)%idx(idx)
            ngllx = Tdomain%specel(n)%ngllx; ngllz = Tdomain%specel(n)%ngllz
            do i = 0, 1
                Tdomain%specel(n)%Forces(1:ngllx-2,1:ngllz-2,i) = &
                    Tdomain%specel(n)%MassMat(:,:) * Tdomain%specel(n)%Forces(1:ngllx-2,1:ngllz-2,i)
            enddo
        enddo

        do idx = 1, size(face_active(jiter)%idx)
            nf = face_active(jiter)%idx(idx)
            ngll = Tdomain%sFace(nf)%ngll
            do i = 0, 1
                Tdomain%sFace(nf)%Forces(:,i) = Tdomain%sFace(nf)%MassMat(:) * Tdomain%sFace(nf)%Forces(:,i)
            enddo
            if (Tdomain%sFace(nf)%Abs .or. Tdomain%sFace(nf)%reflex) Tdomain%sFace(nf)%Forces = 0
        enddo

        do idx = 1, size(vert_active(jiter)%idx)
            nv = vert_active(jiter)%idx(idx)
            do i = 0, 1
                Tdomain%sVertex(nv)%Forces(i) = Tdomain%sVertex(nv)%MassMat * Tdomain%sVertex(nv)%Forces(i)
            enddo
            if (Tdomain%sVertex(nv)%Abs .or. Tdomain%sVertex(nv)%reflex) Tdomain%sVertex(nv)%Forces = 0
        enddo
    end subroutine apply_a_global

    function fact2k(k) result(f)
        implicit none
        integer, intent(in) :: k
        real(fpp) :: f
        integer :: i
        f = 1._fpp
        do i = 2, 2*k
            f = f * real(i,fpp)
        enddo
    end function fact2k

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
