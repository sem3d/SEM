!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file read_mesh.F90
!!\brief Contient la subroutine read_mesh().
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Assure la lecture du maillage (mesh_file) et lecture des caractéristiques des matériaux.
!!
!! Ecriture de data/sem/mesh_echo
!! La racine du nom du maillage doit etre saisie dans Parametrage/sem/input.spec
!!
!! \param type (domain), intent (INOUT) Tdomain
!<


subroutine read_mesh(tDomain)

    use sdomain
    use semdatafiles
    use treceivers, only : read_receiver_file

    implicit none
    character (len=MAX_FILE_SIZE) :: fnamef

    integer :: i,j,j_aus,n_aus,k_aus,w_face
    type (domain), intent (INOUT) :: Tdomain

    integer , dimension (:), allocatable :: Ipointer
    integer :: n_dime

    ! Read Mesh properties
    ! The reading is local to the grid

    call semname_read_mesh_rank(Tdomain%mesh_file,Tdomain%Mpi_var%my_rank,fnamef)

    write(*,*)"ouverture fichier:",trim(adjustl(fnamef))
    open (12,file=fnamef, status="old", form="formatted")
    read (12,*) n_dime
    if (n_dime /= 2) then
        write (*,*) "A dimension different from 2 is not yet taken into account"
        stop
    endif
    read (12,*) Tdomain%n_glob_nodes
    allocate (Tdomain%Coord_nodes(0:1,0:Tdomain%n_glob_nodes-1))
    do i = 0,Tdomain%n_glob_nodes-1
        read (12,*) (Tdomain%Coord_nodes(j,i), j=0,1)
    enddo

    read (12,*) Tdomain%n_mat
    read (12,*) Tdomain%n_elem
    read (12,*) Tdomain%n_nodes
    read (12,*)
    allocate (Tdomain%specel(0:Tdomain%n_elem-1))

    j = Tdomain%n_nodes + 9
    allocate (Ipointer (0:j-1))
    do i = 0,  Tdomain%n_elem - 1
        allocate (Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1))
        read (12,*) (Ipointer (j_aus),j_aus=0,j-1)
        Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1) = Ipointer (0:Tdomain%n_nodes-1)
        j_aus = Tdomain%n_nodes
        Tdomain%specel(i)%Near_Face(0:3) = Ipointer(j_aus:j_aus+3)
        Tdomain%specel(i)%Near_Vertex(0:3) = Ipointer(j_aus+4:j_aus+7)
        Tdomain%specel(i)%mat_index = Ipointer(j_aus+8)
        Tdomain%specel(i)%OUTPUT = .true.
    enddo
    deallocate (Ipointer)
    read (12,*)
    read (12,*) Tdomain%n_face
    read (12,*)
    allocate (Tdomain%sFace(0:Tdomain%n_face-1))


    do i = 0, Tdomain%n_face-1
        read(12,*) Tdomain%sface(i)%Near_element(0:1), Tdomain%sface(i)%Which_face(0:1), Tdomain%sface(i)%Near_Vertex(0:1)
    enddo
    read (12,*)
    read (12,*) Tdomain%n_vertex
    read (12,*)
    allocate (Tdomain%svertex(0:Tdomain%n_vertex-1))
    do i = 0, Tdomain%n_vertex-1
        read(12,*) Tdomain%svertex(i)%Glob_numbering
    enddo

    ! Information about super-object
    if (Tdomain%logicD%super_object) then
        read(12,*)
        n_aus = 0
        do i = 0, Tdomain%n_super_object-1
            read(12,*) j_aus   ! the index is used to say if it is a fault (for mortar act here)
            if (j_aus == 0 ) then
                Tdomain%logicD%super_object_local_present = .true.
                read(12,*) Tdomain%sFault(n_aus)%n_face
                allocate (Tdomain%sFault(n_aus)%fFace(0:Tdomain%sFault(n_aus)%n_face-1))
                do j = 0, Tdomain%sFault(n_aus)%n_face -1
                    read(12,*)  Tdomain%sFault(n_aus)%fFace(j)%Face_UP, Tdomain%sFault(n_aus)%fFace(j)%Face_DOWN, &
                        Tdomain%sFault(n_aus)%fFace(j)%Face_to_Vertex(0:1), k_aus
                    Tdomain%sFault(n_aus)%fFace(j)%Coherency = .true.
                    if (k_aus == 0 ) Tdomain%sFault(n_aus)%fFace(j)%Coherency = .false.
                enddo
                read (12,*)
                read(12,*) Tdomain%sFault(n_aus)%n_vertex
                allocate (Tdomain%sFault(n_aus)%fVertex(0:Tdomain%sFault(n_aus)%n_vertex-1))
                do j = 0, Tdomain%sFault(n_aus)%n_vertex -1
                    read(12,*)  Tdomain%sFault(n_aus)%fVertex(j)%Vertex_UP,Tdomain%sFault(n_aus)%fVertex(j)%Vertex_DOWN
                enddo
            else if (j_aus == -1) then
                Tdomain%logicD%super_object_local_present = .false.
                Tdomain%sFault(n_aus)%n_face = 0; Tdomain%sFault(n_aus)%n_vertex = 0
            endif
            n_aus = n_aus + 1
        enddo
    endif

    read (12,*)
    read (12,*)
    read (12,*) n_aus
    if (n_aus /= Tdomain%MPi_var%n_proc) then
        write (*,*) "The number of selected processors is not the same as chosen for the mesh partition"
        write (*,*) " For processor " , Tdomain%Mpi_var%my_rank, " the two values are : " ,  Tdomain%Mpi_var%n_proc, n_aus
        stop
    endif
    read (12,*) Tdomain%n_communications
    allocate (Tdomain%Communication_List(0:Tdomain%n_communications-1))
    allocate (Tdomain%sWall(0:Tdomain%n_communications-1))
    do j = 0, Tdomain%n_communications-1
        read (12,*) Tdomain%Communication_List(j)
        read (12,*)
        read (12,*)
        read (12,*) Tdomain%sWall(j)%n_faces
        allocate (Tdomain%sWall(j)%Face_List(0:Tdomain%sWall(j)%n_faces-1))
        allocate (Tdomain%sWall(j)%Face_Coherency(0:Tdomain%sWall(j)%n_faces-1))
        do i = 0, Tdomain%sWall(j)%n_faces-1
            read (12,*) Tdomain%sWall(j)%Face_List(i), Tdomain%sWall(j)%Face_Coherency(i)
        enddo
        read (12,*)
        read (12,*)
        read (12,*) Tdomain%sWall(j)%n_vertices
        allocate (Tdomain%sWall(j)%Vertex_List(0:Tdomain%sWall(j)%n_vertices-1))
        do i = 0, Tdomain%sWall(j)%n_vertices-1
            read (12,*) Tdomain%sWall(j)%Vertex_List(i)
        enddo
        if (Tdomain%logicD%super_object) then
            read (12,*)
            read (12,*)
            read (12,*) Tdomain%sWall(j)%n_vertex_superobject
            if (Tdomain%sWall(j)%n_vertex_superobject > 0) &
                allocate (Tdomain%sWall(j)%Vertex_SuperObject_List(0:Tdomain%sWall(j)%n_vertex_superobject-1))
            !	                  allocate (Tdomain%sWall(j)%Vertex_SuperObject_List(0:Tdomain%sWall(j)%n_vertices-1))
            do i = 0, Tdomain%sWall(j)%n_vertex_superobject - 1
                read (12,*) Tdomain%sWall(j)%Vertex_SuperObject_List(i)
            enddo
        endif
        read(12,*)
    enddo
    close (12)

    ! Define a coherency for the face ordering and numbering
    do i =0, Tdomain%n_face-1
        j_aus = Tdomain%sFace(i)%Near_Vertex(0)
        n_aus = Tdomain%sFace(i)%Near_Element(0)
        w_face = Tdomain%sFace(i)%Which_face(0)

        if (w_face == 0 .or. w_face == 3 ) then
            if (Tdomain%specel(n_aus)%Control_Nodes(0) /= Tdomain%sVertex(j_aus)%Glob_numbering)  then
                j_aus = Tdomain%sFace(i)%Near_vertex(0)
                Tdomain%sFace(i)%Near_vertex(0) = Tdomain%sFace(i)%Near_vertex(1)
                Tdomain%sFace(i)%Near_vertex(1) = j_aus
            endif
        else if (w_face == 1 ) then
            if (Tdomain%specel(n_aus)%Control_Nodes(1) /= Tdomain%sVertex(j_aus)%Glob_numbering)  then
                j_aus = Tdomain%sFace(i)%Near_vertex(0)
                Tdomain%sFace(i)%Near_vertex(0) = Tdomain%sFace(i)%Near_vertex(1)
                Tdomain%sFace(i)%Near_vertex(1) = j_aus
            endif
        else
            if (Tdomain%specel(n_aus)%Control_Nodes(3) /= Tdomain%sVertex(j_aus)%Glob_numbering)  then
                j_aus = Tdomain%sFace(i)%Near_vertex(0)
                Tdomain%sFace(i)%Near_vertex(0) = Tdomain%sFace(i)%Near_vertex(1)
                Tdomain%sFace(i)%Near_vertex(1) = j_aus
            endif
        endif
        j_aus = Tdomain%sFace(i)%Near_Vertex(0)
        n_aus = Tdomain%sFace(i)%Near_Element(1)
        Tdomain%sFace(i)%coherency = .false.
        if (n_aus > -1 ) then
            w_face = Tdomain%sFace(i)%Which_face(1)
            if (w_face == 0 .or. w_face == 3 ) then
                if (Tdomain%specel(n_aus)%Control_Nodes(0) == Tdomain%sVertex(j_aus)%Glob_numbering)  &
                    Tdomain%sFace(i)%coherency = .true.
            else if  (w_face == 1 ) then
                if (Tdomain%specel(n_aus)%Control_Nodes(1) == Tdomain%sVertex(j_aus)%Glob_numbering)   &
                    Tdomain%sFace(i)%coherency = .true.
            else
                if (Tdomain%specel(n_aus)%Control_Nodes(3) == Tdomain%sVertex(j_aus)%Glob_numbering) &
                    Tdomain%sFace(i)%coherency = .true.
            endif
        elseif (Tdomain%sFace(i)%Near_Element(1) == -1) then
            ! Set Coherency for faces on the Boundaries
            Tdomain%sFace(i)%coherency = .true.
        endif
    enddo

    if (Tdomain%Implicitness == TIME_INTEG_SEMI_IMPLICIT) &
        call set_Vertex_Valence (Tdomain)

    call read_material_file(Tdomain)

    call read_receiver_file(Tdomain)

end subroutine read_mesh




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
