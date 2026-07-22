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

        ! Only elem_active(m)/face_active(m)/vert_active(m) (== the
        ! %modified set) ever get read back out of elem_corr/face_corr/
        ! vert_corr (in the commit step below), so only those need zeroing.
        do idx = 1, size(elem_active(m)%idx)
            elem_corr(elem_active(m)%idx(idx))%u = 0._fpp
        enddo
        do idx = 1, size(face_active(m)%idx)
            face_corr(face_active(m)%idx(idx))%u = 0._fpp
        enddo
        do idx = 1, size(vert_active(m)%idx)
            vert_corr(vert_active(m)%idx(idx),0:1) = 0._fpp
        enddo

        do k = 1, m
            call apply_a_global(Tdomain, k, m)
            if (k >= 2) then
                ck = 2._fpp*dt**(2*k)/fact2k(k)
                do idx = 1, size(elem_active(m)%idx)
                    n = elem_active(m)%idx(idx)
                    elem_corr(n)%u = elem_corr(n)%u + ck*elem_forces_interior(Tdomain,n)
                enddo
                do idx = 1, size(face_active(m)%idx)
                    n = face_active(m)%idx(idx)
                    face_corr(n)%u = face_corr(n)%u + ck*Tdomain%sFace(n)%Forces
                enddo
                do idx = 1, size(vert_active(m)%idx)
                    n = vert_active(m)%idx(idx)
                    vert_corr(n,0:1) = vert_corr(n,0:1) + ck*Tdomain%sVertex(n)%Forces(0:1)
                enddo
            endif
        enddo

        ! Apply the correction (regional: only elem_active(m)/face_active(m)/
        ! vert_active(m) == %modified -- everywhere else keeps the plain
        ! base-Newmark result computed above untouched) and resync
        ! Veloc = (Displ_new - depl(n))/dt EVERYWHERE so the next step's
        ! classic-Newmark base call sees a consistent pair (a no-op resync
        ! where corr=0, since base Newmark(beta=0,gamma=0.5) already
        ! satisfies this identity on its own).
        do idx = 1, size(elem_active(m)%idx)
            n = elem_active(m)%idx(idx)
            Tdomain%specel(n)%Displ = Tdomain%specel(n)%Displ + elem_corr(n)%u
        enddo
        do n = 0, Tdomain%n_elem - 1
            Tdomain%specel(n)%Veloc = (Tdomain%specel(n)%Displ - elem_un(n)%u)/dt
        enddo
        do idx = 1, size(face_active(m)%idx)
            n = face_active(m)%idx(idx)
            Tdomain%sFace(n)%Displ = Tdomain%sFace(n)%Displ + face_corr(n)%u
        enddo
        do n = 0, Tdomain%n_face - 1
            Tdomain%sFace(n)%Veloc = (Tdomain%sFace(n)%Displ - face_un(n)%u)/dt
        enddo
        do idx = 1, size(vert_active(m)%idx)
            n = vert_active(m)%idx(idx)
            Tdomain%sVertex(n)%Displ = Tdomain%sVertex(n)%Displ + vert_corr(n,0:1)
        enddo
        do n = 0, Tdomain%n_vertex - 1
            Tdomain%sVertex(n)%Veloc = (Tdomain%sVertex(n)%Displ - vert_un(n,0:1))/dt
        enddo
    end subroutine NewmarkModified

    !> Allocate the module-level scratch on first use (after allocate_domain
    !! has sized Displ/Veloc/Forces, so mold= copies the right shapes).
    subroutine ensure_scratch(Tdomain)
        implicit none
        type(domain), intent(in) :: Tdomain
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

        ! m<2 never reaches the correction loop (early return in
        ! NewmarkModified) so the halo is never read -- skip building it.
        if (Tdomain%TimeD%modified_order >= 2) &
            call build_regional_halo(Tdomain, Tdomain%TimeD%modified_order)

        scratch_ready = .true.
    end subroutine ensure_scratch

    !> Precompute, once, elem_hop (BFS graph-distance from the
    !! %modified set, capped at mm-1) and the per-iteration index lists
    !! elem_active(j)/face_active(j)/vert_active(j), j=1..mm -- see the
    !! module-level comment above for why these must SHRINK with j rather
    !! than use one fixed halo for every iteration.
    subroutine build_regional_halo(Tdomain, mm)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: mm

        integer, dimension(:), allocatable :: queue, face_hop, vert_hop
        integer :: qhead, qtail, n, k, nf, nv, e0, e1, j, d, hop_e0, hop_e1

        allocate(elem_hop(0:Tdomain%n_elem-1))
        allocate(face_hop(0:Tdomain%n_face-1))
        allocate(vert_hop(0:Tdomain%n_vertex-1))
        elem_hop = -1
        face_hop = -1
        vert_hop = -1

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
            do k = 0, 3
                nf = Tdomain%specel(n)%Near_Face(k)
                e0 = Tdomain%sFace(nf)%Near_Element(0)
                e1 = Tdomain%sFace(nf)%Near_Element(1)
                if (e0 == n) then
                    if (e1 > -1) then
                        if (elem_hop(e1) < 0) then
                            elem_hop(e1) = d + 1
                            queue(qtail) = e1; qtail = qtail + 1
                        end if
                    end if
                else
                    if (elem_hop(e0) < 0) then
                        elem_hop(e0) = d + 1
                        queue(qtail) = e0; qtail = qtail + 1
                    end if
                end if
            end do
        end do
        deallocate(queue)

        ! Derived face hop = min over its (1 or 2) neighbouring elements'
        ! hop, treating "-1 / never reached" as +infinity unless both sides
        ! are -1 (then the face itself is never reached either).
        do nf = 0, Tdomain%n_face - 1
            e0 = Tdomain%sFace(nf)%Near_Element(0)
            e1 = Tdomain%sFace(nf)%Near_Element(1)
            hop_e0 = elem_hop(e0)
            if (e1 > -1) then
                hop_e1 = elem_hop(e1)
            else
                hop_e1 = -1
            end if
            if (hop_e0 < 0) then
                face_hop(nf) = hop_e1
            else if (hop_e1 < 0) then
                face_hop(nf) = hop_e0
            else
                face_hop(nf) = min(hop_e0, hop_e1)
            end if
        end do

        ! Derived vertex hop = min over every element touching that vertex.
        do n = 0, Tdomain%n_elem - 1
            d = elem_hop(n)
            if (d < 0) cycle
            do k = 0, 3
                nv = Tdomain%specel(n)%Near_Vertex(k)
                if (vert_hop(nv) < 0 .or. d < vert_hop(nv)) vert_hop(nv) = d
            end do
        end do

        allocate(elem_active(mm), face_active(mm), vert_active(mm))
        do j = 1, mm
            call pack_indices(elem_hop, mm - j, elem_active(j)%idx)
            call pack_indices(face_hop, mm - j, face_active(j)%idx)
            call pack_indices(vert_hop, mm - j, vert_active(j)%idx)
        end do
        deallocate(face_hop, vert_hop)
    end subroutine build_regional_halo

    !> Indices n (0-based) with 0 <= hop(n) <= dmax.
    subroutine pack_indices(hop, dmax, idx)
        implicit none
        integer, dimension(0:), intent(in) :: hop
        integer, intent(in) :: dmax
        integer, dimension(:), allocatable, intent(out) :: idx
        integer :: n, cnt

        cnt = 0
        do n = 0, size(hop) - 1
            if (hop(n) >= 0 .and. hop(n) <= dmax) cnt = cnt + 1
        end do
        allocate(idx(cnt))
        cnt = 0
        do n = 0, size(hop) - 1
            if (hop(n) >= 0 .and. hop(n) <= dmax) then
                cnt = cnt + 1
                idx(cnt) = n
            end if
        end do
    end subroutine pack_indices

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
