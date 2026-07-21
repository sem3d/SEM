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

    ! Snapshot of depl(n) (elem_un/face_un/vert_un) and the accumulated
    ! k=2..m correction (elem_corr/face_corr/vert_corr), one entry per
    ! element/face; vertices are fixed-size (0:1) so a plain array suffices.
    type(elem_scratch_t), dimension(:), allocatable, save :: elem_un, elem_corr
    type(face_scratch_t), dimension(:), allocatable, save :: face_un, face_corr
    real(fpp), dimension(:,:), allocatable, save :: vert_un, vert_corr
    logical, save :: scratch_ready = .false.

contains

    subroutine NewmarkModified(Tdomain)
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: n, k, m
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

        do n = 0, Tdomain%n_elem - 1
            elem_corr(n)%u = 0._fpp
        enddo
        do n = 0, Tdomain%n_face - 1
            face_corr(n)%u = 0._fpp
        enddo
        vert_corr = 0._fpp

        do k = 1, m
            call apply_a_global(Tdomain)
            if (k >= 2) then
                ck = 2._fpp*dt**(2*k)/fact2k(k)
                do n = 0, Tdomain%n_elem - 1
                    elem_corr(n)%u = elem_corr(n)%u + ck*elem_forces_interior(Tdomain,n)
                enddo
                do n = 0, Tdomain%n_face - 1
                    face_corr(n)%u = face_corr(n)%u + ck*Tdomain%sFace(n)%Forces
                enddo
                do n = 0, Tdomain%n_vertex - 1
                    vert_corr(n,0:1) = vert_corr(n,0:1) + ck*Tdomain%sVertex(n)%Forces(0:1)
                enddo
            endif
        enddo

        ! Apply the correction and resync Veloc = (Displ_new - depl(n))/dt so
        ! the next step's classic-Newmark base call sees a consistent pair.
        do n = 0, Tdomain%n_elem - 1
            Tdomain%specel(n)%Displ = Tdomain%specel(n)%Displ + elem_corr(n)%u
            Tdomain%specel(n)%Veloc = (Tdomain%specel(n)%Displ - elem_un(n)%u)/dt
        enddo
        do n = 0, Tdomain%n_face - 1
            Tdomain%sFace(n)%Displ = Tdomain%sFace(n)%Displ + face_corr(n)%u
            Tdomain%sFace(n)%Veloc = (Tdomain%sFace(n)%Displ - face_un(n)%u)/dt
        enddo
        do n = 0, Tdomain%n_vertex - 1
            Tdomain%sVertex(n)%Displ = Tdomain%sVertex(n)%Displ + vert_corr(n,0:1)
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

        scratch_ready = .true.
    end subroutine ensure_scratch

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
    subroutine apply_a_global(Tdomain)
        implicit none
        type(domain), intent(inout) :: Tdomain

        integer :: n, nf, nv, nv_aus, nelem, w_face, ngllx, ngllz, ngll, i, j, mat
        integer :: n_face_pointed, tag_send, tag_receive, i_send, i_stock, i_proc, ierr
        integer, dimension(MPI_STATUS_SIZE) :: status

        ! 1) Gather face/vertex w into element boundary rows (get_Displ_fv2el
        ! reads Face/Vertex%Forces, which already hold w from the caller's
        ! seed or the previous apply_a_global call), then apply the local
        ! element stiffness in place (Forces := -K_e*w, per-element).
        do n = 0, Tdomain%n_elem - 1
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
        ! (verbatim pattern from Newmark.F90's "Communication of Forces").
        do nf = 0, Tdomain%n_face - 1
            nelem = Tdomain%sFace(nf)%Near_element(0)
            w_face = Tdomain%sFace(nf)%Which_face(0)
            Tdomain%sFace(nf)%Forces = 0
            call getInternalF_el2f(Tdomain, nelem, nf, w_face, .true.)
            nelem = Tdomain%sFace(nf)%Near_element(1)
            if (nelem > -1) then
                w_face = Tdomain%sFace(nf)%Which_face(1)
                call getInternalF_el2f(Tdomain, nelem, nf, w_face, Tdomain%sFace(nf)%coherency)
            endif
        enddo

        do nv = 0, Tdomain%n_vertex - 1
            Tdomain%sVertex(nv)%Forces = 0
        enddo

        do n = 0, Tdomain%n_elem - 1
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
        ! here as in the base (classic Newmark) step.
        do n = 0, Tdomain%n_elem - 1
            ngllx = Tdomain%specel(n)%ngllx; ngllz = Tdomain%specel(n)%ngllz
            do i = 0, 1
                Tdomain%specel(n)%Forces(1:ngllx-2,1:ngllz-2,i) = &
                    Tdomain%specel(n)%MassMat(:,:) * Tdomain%specel(n)%Forces(1:ngllx-2,1:ngllz-2,i)
            enddo
        enddo

        do nf = 0, Tdomain%n_face - 1
            ngll = Tdomain%sFace(nf)%ngll
            do i = 0, 1
                Tdomain%sFace(nf)%Forces(:,i) = Tdomain%sFace(nf)%MassMat(:) * Tdomain%sFace(nf)%Forces(:,i)
            enddo
            if (Tdomain%sFace(nf)%Abs .or. Tdomain%sFace(nf)%reflex) Tdomain%sFace(nf)%Forces = 0
        enddo

        do nv = 0, Tdomain%n_vertex - 1
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
