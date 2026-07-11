!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file PML_def.F90
!!\brief La routine PML_definition applique les conditions initiales declares dans le fichier d'entree.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type (Domain), intent (INOUT) Tdomain
!<


subroutine PML_definition (Tdomain)

    ! Modified 01/06/2005 Gaetano Festa

    use sdomain
    use mpi
    implicit none

    type (Domain), intent (INOUT) :: Tdomain

    ! local variables
    integer :: n, mat,n_el0, n_el1, nv, nf, n_pml_faces, i, j
    integer, dimension (:), allocatable :: FacePML_List
    logical, dimension (:), allocatable :: Logical_PML_Vertices,  FacePML_Coherency
    logical, dimension (:), allocatable :: Logical_CPML_Vertices, Logical_ADEPML_Vertices
    logical, dimension (:), allocatable :: is_wall_face

    do n = 0, Tdomain%n_elem -1
        mat = Tdomain%specel(n)%mat_index
        Tdomain%specel(n)%PML = .false.
        Tdomain%specel(n)%CPML = .false.
        Tdomain%specel(n)%ADEPML = .false.
        if (Tdomain%sSubDomain(mat)%material_type == "P" .or. &
            Tdomain%sSubDomain(mat)%material_type == "L") then
            Tdomain%specel(n)%PML = .true.
            select case (Tdomain%pml_type)
            case (1)
                Tdomain%specel(n)%CPML = .false.
                Tdomain%specel(n)%ADEPML = .false.
            case (2)
                STOP "Wrong choice for PML types : FPML has been disabled"
                Tdomain%specel(n)%CPML = .false.
                Tdomain%specel(n)%ADEPML = .false.
            case (3)
                Tdomain%specel(n)%CPML = .true.
                Tdomain%specel(n)%ADEPML = .false.
            case (4)
                Tdomain%specel(n)%CPML = .false.
                Tdomain%specel(n)%ADEPML = .true.
            case default
                STOP "Wrong choice for PML types : it should be 1, 2, or 3"
            end select
        endif
    enddo

    ! Partition-wall faces also have Near_Element(1) == -1 locally, but they are INTERNAL
    ! faces of the global mesh: they must not receive boundary flags (abs/freesurf/reflex)
    ! and their PML status is decided by the cross-rank handshake below, mirroring the
    ! serial both-elements rule. In serial (n_communications = 0) the mask is all false.
    allocate (is_wall_face(0:Tdomain%n_face-1))
    is_wall_face = .false.
    do n = 0, Tdomain%n_communications-1
        do i = 0, Tdomain%sWall(n)%n_faces-1
            is_wall_face(Tdomain%sWall(n)%Face_List(i)) = .true.
        enddo
    enddo

    ! ASSIGNING Boundary Condition Type
    ! Flags on absorbing/free surface Faces
    do n = 0, Tdomain%n_face-1
        n_el0 = Tdomain%sFace(n)%Near_Element(0)
        n_el1 = Tdomain%sFace(n)%Near_Element(1)
        if (n_el1 == -1 .and. .not. Tdomain%sFace(n)%is_sf_iface .and. .not. is_wall_face(n)) then
            if (Tdomain%type_bc == DG_BC_ABS) then
                Tdomain%sFace(n)%abs      = .true.
                Tdomain%sFace(n)%freesurf = .false.
                Tdomain%sFace(n)%reflex   = .false.
            elseif (Tdomain%type_bc == DG_BC_FREE) then
                Tdomain%sFace(n)%abs      = .false.
                Tdomain%sFace(n)%freesurf = .true.
                Tdomain%sFace(n)%reflex   = .false.
            elseif (Tdomain%type_bc == DG_BC_REFL) then
                Tdomain%sFace(n)%abs      = .false.
                Tdomain%sFace(n)%freesurf = .false.
                Tdomain%sFace(n)%reflex   = .true.
            endif
        else
            Tdomain%sFace(n)%abs = .false.
            Tdomain%sFace(n)%freesurf = .false.
            Tdomain%sFace(n)%reflex   = .false.
        endif
        i = Tdomain%sFace(n)%Near_Vertex(0)
        j = Tdomain%sFace(n)%Near_Vertex(1)
        Tdomain%sVertex(i)%abs    = Tdomain%sFace(n)%abs
        Tdomain%sVertex(j)%abs    = Tdomain%sFace(n)%abs
        Tdomain%sVertex(i)%reflex = Tdomain%sFace(n)%reflex
        Tdomain%sVertex(j)%reflex = Tdomain%sFace(n)%reflex
    enddo

    ! Special Treatment for PML
    do n = 0, Tdomain%n_face-1
        Tdomain%sFace(n)%PML = .false.
        Tdomain%sFace(n)%CPML = .false.
        Tdomain%sFace(n)%ADEPML = .false.
        !Tdomain%sFace(n)%Abs = .false.
        n_el0 = Tdomain%sFace(n)%Near_Element(0)
        n_el1 = Tdomain%sFace(n)%Near_Element(1)
        ! Dealing With PML Faces at the border of the domain.
        if (n_el1 > -1) then
            if (Tdomain%specel(n_el0)%PML .and. Tdomain%specel(n_el1)%PML) Tdomain%sFace(n)%PML = .true.
            if (Tdomain%specel(n_el0)%CPML .and. Tdomain%specel(n_el1)%CPML) Tdomain%sFace(n)%CPML = .true.
            if (Tdomain%specel(n_el0)%ADEPML .and. Tdomain%specel(n_el1)%ADEPML) Tdomain%sFace(n)%ADEPML = .true.
        else if (.not. is_wall_face(n)) then
            ! True domain-boundary face: rigid termination on the outer PML edge.
            ! Wall faces must NOT enter here (they are internal; a one-sided Reflex
            ! zeroes Forces on one rank only and desynchronizes the shared DOFs).
            if (Tdomain%specel(n_el0)%PML) then
                mat = Tdomain%specel(n_el0)%mat_index
                if (Tdomain%sFace(n)%which_face(0) == 0 .and. Tdomain%sSubdomain(mat)%Pz &
                    .and. Tdomain%sSubdomain(mat)%Down) Tdomain%sFace(n)%Reflex = .true.
                if (Tdomain%sFace(n)%which_face(0) == 1 .and. Tdomain%sSubdomain(mat)%Px  &
                    .and.(.not.  Tdomain%sSubdomain(mat)%Left)) Tdomain%sFace(n)%Reflex = .true.
                if (Tdomain%sFace(n)%which_face(0) == 2 .and. Tdomain%sSubdomain(mat)%Pz &
                    .and. (.not. Tdomain%sSubdomain(mat)%Down)) Tdomain%sFace(n)%Reflex = .true.
                if (Tdomain%sFace(n)%which_face(0) == 3 .and. Tdomain%sSubdomain(mat)%Px  &
                    .and. Tdomain%sSubdomain(mat)%Left) Tdomain%sFace(n)%Reflex = .true.
                if (Tdomain%sFace(n)%Reflex) then
                    Tdomain%sFace(n)%Abs = .false. ; Tdomain%sFace(n)%Freesurf = .false.
                    i=Tdomain%sFace(n)%Near_Vertex(0) ; j=Tdomain%sFace(n)%Near_Vertex(1)
                    Tdomain%sVertex(i)%Reflex = .true. ; Tdomain%sVertex(j)%Reflex = .true.
                endif
            endif
        endif
    enddo

    ! Partition-wall faces: exchange the LOCAL element flags with the neighbour rank and
    ! apply the same rule as serial (line above): PML iff BOTH adjacent elements are PML
    ! (idem CPML/ADEPML). Face_List order is matched across paired walls (same guarantee
    ! the Forces exchange relies on). Replaces a directional which_face/Px/Pz heuristic
    ! that missed faces parallel to the damping direction and all corner-block (Px.and.Pz)
    ! faces, leaving them -- and their vertices -- running the regular corrector in MPI.
    do n = 0, Tdomain%n_communications-1
        block
            integer :: k, nfw, partner, ierr
            integer :: mpistat(MPI_STATUS_SIZE)
            integer, dimension(:), allocatable :: loc_flag, rem_flag
            nfw = Tdomain%sWall(n)%n_faces
            if (nfw > 0) then
                allocate (loc_flag(0:nfw-1), rem_flag(0:nfw-1))
                loc_flag = 0; rem_flag = 0
                do k = 0, nfw-1
                    nf = Tdomain%sWall(n)%Face_List(k)
                    n_el0 = Tdomain%sFace(nf)%Near_Element(0)
                    if (Tdomain%specel(n_el0)%PML)    loc_flag(k) = loc_flag(k) + 1
                    if (Tdomain%specel(n_el0)%CPML)   loc_flag(k) = loc_flag(k) + 2
                    if (Tdomain%specel(n_el0)%ADEPML) loc_flag(k) = loc_flag(k) + 4
                enddo
                partner = Tdomain%Communication_list(n)
                call MPI_Sendrecv (loc_flag, nfw, MPI_INTEGER, partner, 810, &
                                   rem_flag, nfw, MPI_INTEGER, partner, 810, &
                                   Tdomain%communicateur, mpistat, ierr)
                do k = 0, nfw-1
                    nf = Tdomain%sWall(n)%Face_List(k)
                    if (iand(loc_flag(k),1) > 0 .and. iand(rem_flag(k),1) > 0) Tdomain%sFace(nf)%PML = .true.
                    if (iand(loc_flag(k),2) > 0 .and. iand(rem_flag(k),2) > 0) Tdomain%sFace(nf)%CPML = .true.
                    if (iand(loc_flag(k),4) > 0 .and. iand(rem_flag(k),4) > 0) Tdomain%sFace(nf)%ADEPML = .true.
                enddo
                deallocate (loc_flag, rem_flag)
            endif
        end block
    enddo

    ! Define PML faces that need to communicate
    do n= 0, Tdomain%n_communications-1
        n_pml_faces = 0
        allocate (FacePML_List (0:Tdomain%SWall(n)%n_faces-1))
        allocate (FacePML_Coherency (0:Tdomain%SWall(n)%n_faces-1))
        do i = 0, Tdomain%sWall(n)%n_faces-1
            nf = Tdomain%sWall(n)%Face_List(i)
            if (Tdomain%sFace(nf)%PML) then
                FacePML_List(n_pml_faces) = nf
                FacePML_Coherency(n_pml_faces) =Tdomain%sWall(n)%Face_Coherency(i)
                n_pml_faces = n_pml_faces + 1

            endif
        enddo
        Tdomain%sWall(n)%n_pml_faces = n_pml_faces
        allocate (Tdomain%sWall(n)%FacePML_List(0:n_pml_faces-1))
        allocate (Tdomain%sWall(n)%FacePML_Coherency(0:n_pml_faces-1))
        Tdomain%sWall(n)%FacePML_List(0:n_pml_faces-1) =FacePML_List(0:n_pml_faces-1)
        Tdomain%sWall(n)%FacePML_Coherency(0:n_pml_faces-1) =FacePML_Coherency(0:n_pml_faces-1)
        deallocate (FacePML_List)
        deallocate (FacePML_Coherency)
    enddo

    allocate (Logical_PML_vertices(0:Tdomain%n_vertex-1))
    allocate (Logical_CPML_vertices(0:Tdomain%n_vertex-1))
    allocate (Logical_ADEPML_vertices(0:Tdomain%n_vertex-1))
    Logical_PML_vertices = .true.
    Logical_CPML_vertices = .true.
    Logical_ADEPML_vertices = .true.

    do n = 0, Tdomain%n_face-1
        if (.not. Tdomain%sFace(n)%PML) then
            nv = Tdomain%sFace(n)%Near_Vertex(0)
            if (nv<0 .or. nv >=Tdomain%n_vertex) stop "ERROR - PML_def : invalid near_vertex"
            Logical_PML_Vertices(nv) = .false.
            nv = Tdomain%sFace(n)%Near_Vertex(1)
            if (nv<0 .or. nv >=Tdomain%n_vertex) stop "ERROR - PML_def : invalid near_vertex"
            Logical_PML_Vertices(nv) = .false.
        endif
        if (.not. Tdomain%sFace(n)%CPML) then
            nv = Tdomain%sFace(n)%Near_Vertex(0)
            if (nv<0 .or. nv >=Tdomain%n_vertex) stop "ERROR - PML_def : invalid near_vertex"
            Logical_CPML_Vertices(nv) = .false.
            nv = Tdomain%sFace(n)%Near_Vertex(1)
            if (nv<0 .or. nv >=Tdomain%n_vertex) stop "ERROR - PML_def : invalid near_vertex"
            Logical_CPML_Vertices(nv) = .false.
        endif
        if (.not. Tdomain%sFace(n)%ADEPML) then
            nv = Tdomain%sFace(n)%Near_Vertex(0)
            Logical_ADEPML_Vertices(nv) = .false.
            nv = Tdomain%sFace(n)%Near_Vertex(1)
            Logical_ADEPML_Vertices(nv) = .false.
        endif

    enddo

    ! Reconcile vertex flags across ranks. The cascade above only sees LOCAL faces: a wall
    ! vertex whose non-PML (or PML) faces live on the neighbour rank ends up with a flag
    ! different from serial, and possibly different between the ranks sharing it. Exchange
    ! a per-vertex bitmask (Vertex_List order is matched across paired walls) and combine:
    ! AND for PML/CPML/ADEPML (serial rule: all incident faces), OR for reflex (outer PML
    ! termination seen by only one of the owning ranks). Snapshots are sent (not the live
    ! flags) so corner vertices shared by 3+ ranks combine order-independently.
    if (Tdomain%n_communications > 0) then
        block
            integer :: k, nvw, partner, ierr
            integer :: mpistat(MPI_STATUS_SIZE)
            integer, dimension(:), allocatable :: snap, loc_flag, rem_flag
            allocate (snap(0:Tdomain%n_vertex-1))
            do k = 0, Tdomain%n_vertex-1
                snap(k) = 0
                if (Logical_PML_Vertices(k))    snap(k) = snap(k) + 1
                if (Logical_CPML_Vertices(k))   snap(k) = snap(k) + 2
                if (Logical_ADEPML_Vertices(k)) snap(k) = snap(k) + 4
                if (Tdomain%sVertex(k)%reflex)  snap(k) = snap(k) + 8
            enddo
            do n = 0, Tdomain%n_communications-1
                nvw = Tdomain%sWall(n)%n_vertices
                if (nvw > 0) then
                    allocate (loc_flag(0:nvw-1), rem_flag(0:nvw-1))
                    do k = 0, nvw-1
                        loc_flag(k) = snap(Tdomain%sWall(n)%Vertex_List(k))
                    enddo
                    rem_flag = 0
                    partner = Tdomain%Communication_list(n)
                    call MPI_Sendrecv (loc_flag, nvw, MPI_INTEGER, partner, 820, &
                                       rem_flag, nvw, MPI_INTEGER, partner, 820, &
                                       Tdomain%communicateur, mpistat, ierr)
                    do k = 0, nvw-1
                        nv = Tdomain%sWall(n)%Vertex_List(k)
                        if (iand(rem_flag(k),1) == 0) Logical_PML_Vertices(nv)    = .false.
                        if (iand(rem_flag(k),2) == 0) Logical_CPML_Vertices(nv)   = .false.
                        if (iand(rem_flag(k),4) == 0) Logical_ADEPML_Vertices(nv) = .false.
                        if (iand(rem_flag(k),8) >  0) Tdomain%sVertex(nv)%reflex  = .true.
                    enddo
                    deallocate (loc_flag, rem_flag)
                endif
            enddo
            deallocate (snap)
        end block
    endif

    do n = 0, Tdomain%n_vertex-1
        Tdomain%sVertex(n)%PML = Logical_PML_Vertices (n)
        ! Be careful to the following line which is designed to avoid unusefull
        ! computations in define_array.F90
        Tdomain%sVertex(n)%CPML = Logical_CPML_Vertices (n) .or. Logical_ADEPML_Vertices (n)
        Tdomain%sVertex(n)%ADEPML = Logical_ADEPML_Vertices (n)
    enddo

    deallocate (Logical_PML_Vertices)
    deallocate (Logical_CPML_Vertices)
    deallocate (Logical_ADEPML_Vertices)

    ! Define PML vertices that need to communicate (split-field DumpMass + Forces1/2). Done
    ! AFTER the flags are reconciled above, so the ranks sharing a vertex agree on its PML
    ! status and paired walls build lists of matching size. A PML vertex shared across ranks
    ! must sum its split DumpMass (setup) and Forces1/Forces2 (each step) with the neighbour
    ! -- otherwise its corrector (Correction_Vertex_PML_Veloc) uses incomplete values.
    do n = 0, Tdomain%n_communications-1
        block
            integer :: n_pml_vertices, k
            integer, dimension(:), allocatable :: VertexPML_List
            allocate (VertexPML_List(0:max(Tdomain%sWall(n)%n_vertices-1,0)))
            n_pml_vertices = 0
            do k = 0, Tdomain%sWall(n)%n_vertices-1
                nv = Tdomain%sWall(n)%Vertex_List(k)
                if (Tdomain%sVertex(nv)%PML .and. (.not. Tdomain%sVertex(nv)%CPML) &
                                            .and. (.not. Tdomain%sVertex(nv)%ADEPML)) then
                    VertexPML_List(n_pml_vertices) = nv
                    n_pml_vertices = n_pml_vertices + 1
                endif
            enddo
            Tdomain%sWall(n)%n_pml_vertices = n_pml_vertices
            allocate (Tdomain%sWall(n)%VertexPML_List(0:max(n_pml_vertices-1,0)))
            if (n_pml_vertices > 0) &
                Tdomain%sWall(n)%VertexPML_List(0:n_pml_vertices-1) = VertexPML_List(0:n_pml_vertices-1)
            deallocate (VertexPML_List)
        end block
    enddo

    deallocate (is_wall_face)

    return
end subroutine PML_definition

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
