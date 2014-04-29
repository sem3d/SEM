!>
!!\file PML_def.F90
!!\brief La routine PML_definition applique les conditions initiales déclarées dans le fichier d'entrée.
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
    implicit none

    type (Domain), intent (INOUT) :: Tdomain

    ! local variables
    integer :: n, mat,n_el0, n_el1, nv, nf, n_pml_faces, i
    integer, dimension (:), allocatable :: FacePML_List
    logical, dimension (:), allocatable :: Logical_PML_Vertices, Logical_Abs_Vertices,FacePML_Coherency
    logical, dimension (:), allocatable :: Logical_FPML_Vertices, Logical_CPML_Vertices

    do n = 0, Tdomain%n_elem -1
        mat = Tdomain%specel(n)%mat_index
        Tdomain%specel(n)%PML = .false.
        Tdomain%specel(n)%FPML = .false.
        if (Tdomain%sSubDomain(mat)%material_type == "P" ) then
            Tdomain%specel(n)%PML = .true.
            select case (Tdomain%pml_type)
            case (1)
                Tdomain%specel(n)%FPML = .false.
                Tdomain%specel(n)%CPML = .false.
            case (2)
                Tdomain%specel(n)%FPML = .true.
                Tdomain%specel(n)%CPML = .false.
            case (3)
                Tdomain%specel(n)%FPML = .false.
                Tdomain%specel(n)%CPML = .true.
            case default
                STOP "Wrong choice for PML types : it should be 1, 2, or 3"
            end select
            !if (Tdomain%specel(n)%PML .and. Tdomain%sSubDomain(mat)%Filtering ) Tdomain%specel(n)%FPML =.true.
        endif
    enddo

    do n = 0, Tdomain%n_face-1
        Tdomain%sFace(n)%PML = .false.
        Tdomain%sFace(n)%FPML = .false.
        Tdomain%sFace(n)%CPML = .false.
        Tdomain%sFace(n)%Abs = .false.
        n_el0 = Tdomain%sFace(n)%Near_Element(0)
        n_el1 = Tdomain%sFace(n)%Near_Element(1)
        if (n_el1 > -1) then
            if (Tdomain%specel(n_el0)%PML .and. Tdomain%specel(n_el1)%PML) Tdomain%sFace(n)%PML = .true.
            if (Tdomain%specel(n_el0)%FPML .and. Tdomain%specel(n_el1)%FPML) Tdomain%sFace(n)%FPML = .true.
            if (Tdomain%specel(n_el0)%CPML .and. Tdomain%specel(n_el1)%CPML) Tdomain%sFace(n)%CPML = .true.
        else
            if (Tdomain%specel(n_el0)%PML) then
                mat = Tdomain%specel(n_el0)%mat_index
                if (Tdomain%sFace(n)%which_face(0) == 0 .and. Tdomain%sSubdomain(mat)%Pz &
                    .and. Tdomain%sSubdomain(mat)%Down) Tdomain%sFace(n)%Abs = .true.
                if (Tdomain%sFace(n)%which_face(0) == 1 .and. Tdomain%sSubdomain(mat)%Px  &
                    .and.(.not.  Tdomain%sSubdomain(mat)%Left)) Tdomain%sFace(n)%Abs = .true.
                if (Tdomain%sFace(n)%which_face(0) == 2 .and. Tdomain%sSubdomain(mat)%Pz &
                    .and. (.not. Tdomain%sSubdomain(mat)%Down)) Tdomain%sFace(n)%Abs = .true.
                if (Tdomain%sFace(n)%which_face(0) == 3 .and. Tdomain%sSubdomain(mat)%Px  &
                    .and. Tdomain%sSubdomain(mat)%Left) Tdomain%sFace(n)%Abs = .true.
            endif


        endif
    enddo

    ! Define PML faces that need to communicate
    do n= 0, Tdomain%n_communications-1
        n_pml_faces = 0
        allocate (FacePML_List (0:Tdomain%SWall(n)%n_faces-1))
        allocate (FacePML_Coherency (0:Tdomain%SWall(n)%n_faces-1))
        do i = 0, Tdomain%sWall(n)%n_faces-1
            nf = Tdomain%sWall(n)%Face_List(i)
            n_el0 = Tdomain%sFace(nf)%Near_Element(0)
            if (Tdomain%specel(n_el0)%PML) then
                mat = Tdomain%specel(n_el0)%mat_index
                if (Tdomain%sFace(nf)%which_face(0) == 0 .and. Tdomain%sSubdomain(mat)%Px &
                    .and. (.not. Tdomain%sSubdomain(mat)%Pz)) Tdomain%sFace(nf)%PML = .true.
                if (Tdomain%sFace(nf)%which_face(0) == 2 .and. Tdomain%sSubdomain(mat)%Px &
                    .and. (.not. Tdomain%sSubdomain(mat)%Pz)) Tdomain%sFace(nf)%PML = .true.
                if (Tdomain%sFace(nf)%which_face(0) == 1 .and. Tdomain%sSubdomain(mat)%Pz  &
                    .and. (.not. Tdomain%sSubdomain(mat)%Px)) Tdomain%sFace(nf)%PML = .true.
                if (Tdomain%sFace(nf)%which_face(0) ==3 .and. Tdomain%sSubdomain(mat)%Pz &
                    .and. (.not. Tdomain%sSubdomain(mat)%Px)) Tdomain%sFace(nf)%PML = .true.
            endif
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
    allocate (Logical_FPML_vertices(0:Tdomain%n_vertex-1))
    allocate (Logical_CPML_vertices(0:Tdomain%n_vertex-1))
    allocate (Logical_Abs_vertices(0:Tdomain%n_vertex-1))
    Logical_PML_vertices = .true.
    Logical_FPML_vertices = .true.
    Logical_CPML_vertices = .true.
    Logical_Abs_vertices = .false.

    do n = 0, Tdomain%n_face-1
        if (.not. Tdomain%sFace(n)%PML) then
            nv = Tdomain%sFace(n)%Near_Vertex(0)
            Logical_PML_Vertices(nv) = .false.
            nv = Tdomain%sFace(n)%Near_Vertex(1)
            Logical_PML_Vertices(nv) = .false.
        endif
        if (.not. Tdomain%sFace(n)%FPML) then
            nv = Tdomain%sFace(n)%Near_Vertex(0)
            Logical_FPML_Vertices(nv) = .false.
            nv = Tdomain%sFace(n)%Near_Vertex(1)
            Logical_FPML_Vertices(nv) = .false.
        endif
        if (.not. Tdomain%sFace(n)%CPML) then
            nv = Tdomain%sFace(n)%Near_Vertex(0)
            Logical_CPML_Vertices(nv) = .false.
            nv = Tdomain%sFace(n)%Near_Vertex(1)
            Logical_CPML_Vertices(nv) = .false.
        endif
        if (Tdomain%sFace(n)%Abs) then
            nv = Tdomain%sFace(n)%Near_Vertex(0)
            Logical_Abs_Vertices(nv) = .true.
            nv = Tdomain%sFace(n)%Near_Vertex(1)
            Logical_Abs_Vertices(nv) = .true.
        endif

    enddo

    do n = 0, Tdomain%n_vertex-1
        Tdomain%sVertex(n)%PML = Logical_PML_Vertices (n)
        Tdomain%sVertex(n)%FPML = Logical_FPML_Vertices (n)
        Tdomain%sVertex(n)%CPML = Logical_CPML_Vertices (n)
        Tdomain%sVertex(n)%Abs = Logical_Abs_Vertices (n)
    enddo

    deallocate (Logical_PML_Vertices)
    deallocate (Logical_FPML_Vertices)
    deallocate (Logical_CPML_Vertices)
    deallocate (Logical_Abs_Vertices )
    return
end subroutine PML_definition
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
