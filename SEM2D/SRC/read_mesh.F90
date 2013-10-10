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
    ! Modified by Gaetano Festa 01/06/05
    ! Introducing parallel features 12/10/05
    !

    implicit none
    character (len=MAX_FILE_SIZE) :: fnamef
    character (len=10) :: auxiliary_name

    integer :: i,j,j_aus,n_aus,k_aus,w_face,npml,mat
    type (domain), intent (INOUT) :: Tdomain

    integer , dimension (:), allocatable :: Ipointer
    real :: dtmin, Qp, Qs
    integer :: n_dime

    ! Read Mesh properties
    ! The reading is local to the grid

    call semname_read_mesh_rank(Tdomain%mesh_file,Tdomain%Mpi_var%my_rank,fnamef)

    write(*,*)"ouverture fichier:",fnamef
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


!    ! lecture obsolete car valeurs pas utilisees
!    read (12,*) Tdomain%n_line
!    read (12,*)
!    allocate (Tdomain%Name_Line(0:Tdomain%n_line-1) )
!    allocate (Tdomain%Line_index(0:Tdomain%n_line-1) )
!    do i = 0, Tdomain%n_line -1
!        read (12,*) j, auxiliary_name
!        !Tdomain%Name_line(i) = auxiliary_name (1:1)
!        !Tdomain%Line_index(i) = j
!    enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


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

    if (Tdomain%logicD%run_echo) then

        call semname_read_mesh_echo(Tdomain%Mpi_var%my_rank,fnamef)

        open(92,file=fnamef, form="formatted", status="unknown")
        write (92,*) Tdomain%n_glob_nodes, "Number of global nodes"
        do i = 0,Tdomain%n_glob_nodes-1
            write (92,*) (Tdomain%Coord_nodes(j,i), j=0,1)
        enddo
        write (92,*) Tdomain%n_mat, " Number of subdomains"
        write (92,*) Tdomain%n_line, "  Number of lines"
        do i = 0, Tdomain%n_line -1
            write (92,*) Tdomain%Name_line(i)
        enddo
        write (92,*) Tdomain%n_elem, "Number of elements"
        write (92,*) Tdomain%n_nodes, "Number of nodes"

        j = Tdomain%n_nodes + 5
        do i = 0,  Tdomain%n_elem - 1
            write (92,*) Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1),  &
                Tdomain%specel(i)%Near_Face(0:3),  Tdomain%specel(i)%Near_Vertex(0:3), Tdomain%specel(i)%mat_index
        enddo
        write (92,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write (92,*) Tdomain%n_face, "Number of faces"
        do i = 0, Tdomain%n_face-1
            write(92,*) Tdomain%sface(i)%Near_element(0:1), Tdomain%sface(i)%Which_face(0:1), Tdomain%sface(i)%Near_Vertex(0:1)
        enddo
        write (92,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write (92,*) Tdomain%n_vertex, "Number of Vertices"
        do i = 0, Tdomain%n_vertex-1
            write(92,*) Tdomain%svertex(i)%Glob_numbering
        enddo

    endif

    if (Tdomain%logicD%super_object) then
        write(92,*) " !!!!!!!!!! Here a super object"
        do i = 0, Tdomain%n_fault-1
            write (92,*) "Here we have a fault"
            if (Tdomain%logicD%super_object_local_present) then
                write (92,*) Tdomain%sFault(i)%n_face, "    # Number of faces"
                do j = 0, Tdomain%sFault(i)%n_face -1
                    write(92,*)  Tdomain%sFault(i)%fFace(j)%Face_UP, Tdomain%sFault(i)%fFace(j)%Face_DOWN, &
                        Tdomain%sFault(i)%fFace(j)%Face_to_Vertex(0:1), Tdomain%sFault(i)%fFace(j)%Coherency
                enddo
                write (92,*) Tdomain%sFault(i)%n_vertex, "    # Number of vertices"
                do j = 0, Tdomain%sFault(i)%n_vertex -1
                    write(92,*)  Tdomain%sFault(i)%fVertex(j)%Vertex_UP,Tdomain%sFault(i)%fVertex(j)%Vertex_DOWN
                enddo
            else
                write (92,*) "But not in this processor"
            endif

        enddo
    endif

    write (92,*)
    write (92,*)
    write (92,*) Tdomain%MPi_var%n_proc,  "Number of total processors"
    write (92,*) Tdomain%n_communications,  " Number  of effective communications"
    do j = 0, Tdomain%n_communications-1
        write (92,*) Tdomain%Communication_List(j), "Number of processor"
        write (92,*) " !!!!!!!!!!"
        write (92,*)  " Define faces"
        write (92,*) Tdomain%sWall(j)%n_faces,  "number of faces"
        do i = 0, Tdomain%sWall(j)%n_faces-1
            write (92,*) Tdomain%sWall(j)%Face_List(i)
        enddo
        write (92,*)" !!!!!!!!!!"
        write (92,*)  " Define vertices"
        write (92,*) Tdomain%sWall(j)%n_vertices
        do i = 0, Tdomain%sWall(j)%n_vertices-1
            write (92,*) Tdomain%sWall(j)%Vertex_List(i)
        enddo
        write (92,*) " !!!!!!!!!!!!!!"
        if (Tdomain%logicD%super_object) then
            write (92,*) " Number of super object vertices to communicate "
            write (92,*) Tdomain%sWall(j)%n_vertex_superobject, "Number of super Vertices"
            do i = 0, Tdomain%sWall(j)%n_vertex_superobject - 1
                write (92,*) Tdomain%sWall(j)%Vertex_SuperObject_List(i)
            enddo
            write (92,*) " !!!!!!!!!!!!!!"
        endif
    enddo
    close (92)


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
        endif
    enddo

    ! Read material properties
    npml = 0
    allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))

    call semname_read_inputmesh_parametrage(Tdomain%material_file,fnamef)
    open (13,file=fnamef, status="old", form="formatted")

    read (13,*) n_aus  ! Number of material in file
    if (n_aus<Tdomain%n_mat) then
        write(*,*) "Material file doesn't contain enough definitions"
        stop 1
    end if
    do i = 0, Tdomain%n_mat-1
        read (13,*) Tdomain%sSubDomain(i)%material_type, Tdomain%sSubDomain(i)%Pspeed, &
            Tdomain%sSubDomain(i)%Sspeed, Tdomain%sSubDomain(i)%dDensity, &
            Tdomain%sSubDomain(i)%NGLLx, n_aus, Tdomain%sSubDomain(i)%NGLLz, Tdomain%sSubDomain(i)%Dt, &
            Qp, Qs
        Tdomain%sSubDomain(i)%n_loc_dim = 2
        Tdomain%sSubdomain(i)%wpml = -1
        if ( Tdomain%sSubDomain(i)%NGLLx == Tdomain%sSubDomain(i)%NGLLz)  Tdomain%sSubDomain(i)%n_loc_dim = 1
        if (Tdomain%sSubDomain(i)%material_type == "P" )  then
            Tdomain%sSubDomain(i)%wpml = npml
            npml = npml + 1
        endif
    enddo
    Tdomain%any_PML  = .false.
    if (npml > 0 ) then
        Tdomain%any_PML = .true.
        read(13,*); read(13,*)
        do i = 0,Tdomain%n_mat-1
            if (Tdomain%sSubdomain(i)%material_type == "P" ) then
                read (13,*) Tdomain%sSubdomain(i)%Filtering,  Tdomain%sSubdomain(i)%npow, Tdomain%sSubdomain(i)%Apow, &
                    Tdomain%sSubdomain(i)%Px, Tdomain%sSubdomain(i)%Left, Tdomain%sSubdomain(i)%Pz,  &
                    Tdomain%sSubdomain(i)%Down, Tdomain%sSubdomain(i)%freq, Tdomain%sSubdomain(i)%k
            endif
        enddo
    endif

    close (13)

    do i = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(i)%mat_index
        Tdomain%specel(i)%ngllx = Tdomain%sSubDomain(mat)%NGLLx
        Tdomain%specel(i)%ngllz = Tdomain%sSubDomain(mat)%NGLLz
    enddo

    do i = 0, Tdomain%n_face-1
        n_aus = Tdomain%sFace(i)%Near_Element(0)
        w_face = Tdomain%sFace(i)%Which_face(0)
        mat = Tdomain%specel(n_aus)%mat_index
        if (w_face == 0 .or. w_face==2) then
            Tdomain%sFace(i)%ngll = Tdomain%sSubDomain(mat)%ngllx
        else
            Tdomain%sFace(i)%ngll = Tdomain%sSubDomain(mat)%ngllz
        endif
        Tdomain%sFace(i)%mat_index = mat
        n_aus = Tdomain%sFace(i)%Near_Vertex(0)
        Tdomain%sVertex(n_aus)%mat_index = mat
        n_aus = Tdomain%sFace(i)%Near_Vertex(1)
        Tdomain%sVertex(n_aus)%mat_index = mat
    enddo

    if (Tdomain%logicD%super_object_local_present) then
        do i = 0, Tdomain%n_fault-1
            do j = 0, Tdomain%sFault(i)%n_face-1
                n_aus = Tdomain%sFault(i)%fFace(j)%Face_UP
                Tdomain%sFault(i)%fFace(j)%ngll = Tdomain%sFace(n_aus)%ngll
                w_face = Tdomain%sFace(n_aus)%Which_face(0)
                if (w_face == 0 .or. w_face ==2) then
                    Tdomain%sFault(i)%fFace(j)%mat_dir = 1
                else
                    Tdomain%sFault(i)%fFace(j)%mat_dir = 0
                endif
                Tdomain%sFault(i)%fFace(j)%mat_index = Tdomain%sFace(n_aus)%mat_index
                n_aus = Tdomain%sFault(i)%fFace(j)%Face_to_Vertex(0)
                Tdomain%sFault(i)%fVertex(n_aus)%mat_index =  Tdomain%sFault(i)%fFace(j)%mat_index
                n_aus = Tdomain%sFault(i)%fFace(j)%Face_to_Vertex(1)
                Tdomain%sFault(i)%fVertex(n_aus)%mat_index =  Tdomain%sFault(i)%fFace(j)%mat_index
            enddo
            do j = 0, Tdomain%sFault(i)%n_vertex-1
                Tdomain%sFault(i)%fVertex(j)%Termination = .false.
                if (Tdomain%sFault(i)%fVertex(j)%Vertex_UP == -3) &
                    Tdomain%sFault(i)%fVertex(j)%Termination = .true.
            enddo
        enddo
    endif

    if (Tdomain%logicD%run_echo .and. Tdomain%MPI_var%my_rank ==0) then

        call semname_read_mesh_material_echo(fnamef)
        open(93,file=fnamef, form="formatted", status="unknown")

        write (93,*) "Material properties"
        write (93,*) Tdomain%n_mat, "   Number of material"
        write (93,*) "Type, Pspeed, Sspeed, Density, Dt, NGLLx, NGLLz"
        do i = 0, Tdomain%n_mat-1
            write (93,*) Tdomain%sSubDomain(i)%material_type, Tdomain%sSubDomain(i)%Pspeed, &
                Tdomain%sSubDomain(i)%Sspeed, Tdomain%sSubDomain(i)%DDensity, &
                Tdomain%sSubDomain(i)%Dt, Tdomain%sSubDomain(i)%NGLLx, Tdomain%sSubDomain(i)%NGLLz
        enddo
        if (Tdomain%any_PML ) then
            write (93,*) "Definition of some PML conditions"
            write (93,*) " Filtering,  np-power,A-power, x-direction, left increase, z-direction, down-increase, cut-off frequency"
            do i = 0,Tdomain%n_mat-1
                if (Tdomain%sSubdomain(i)%material_type == "P" ) then
                    write (93,*)  Tdomain%sSubdomain(i)%Filtering,  Tdomain%sSubdomain(i)%npow, Tdomain%sSubdomain(i)%Apow, &
                        Tdomain%sSubdomain(i)%Px, Tdomain%sSubdomain(i)%Left, Tdomain%sSubdomain(i)%Pz,  &
                        Tdomain%sSubdomain(i)%Down, Tdomain%sSubdomain(i)%freq,Tdomain%sSubdomain(i)%k
                endif
            enddo
            close (93)
        endif
    endif

    do i = 0, Tdomain%n_mat-1
        if (Tdomain%sSubdomain(i)%material_type == "P" .and. Tdomain%sSubdomain(i)%Filtering ) &
            Tdomain%sSubdomain(i)%freq = exp (-Tdomain%sSubdomain(i)%freq*Tdomain%sSubdomain(i)%dt/2)
    enddo

    dtmin =1e20
    do i = 0,Tdomain%n_mat-1
        if (Tdomain%sSubDomain(i)%Dt < dtmin ) dtmin = Tdomain%sSubDomain(i)%Dt
    enddo
    Tdomain%TimeD%dtmin = dtmin
    if (dtmin > 0) then
        Tdomain%TimeD%ntimeMax = int (Tdomain%TimeD%Duration/dtmin)
    else
        write (*,*) "Your dt min is zero : verify it"
        stop
    endif

    do i = 0, Tdomain%n_mat-1
        call Lame_coefficients (Tdomain%sSubDomain(i))
    enddo

    if (Tdomain%logicD%save_trace) then

        call semname_read_inputmesh_parametrage(Tdomain%station_file,fnamef)
        open(14,file=fnamef, status="old")

        read (14,*) Tdomain%n_receivers
        read (14,*)
        allocate (Tdomain%sReceiver(0:Tdomain%n_receivers-1))
        do i = 0, Tdomain%n_receivers-1
            read(14,*) Tdomain%sReceiver(i)%Xrec, Tdomain%sReceiver(i)%Zrec
        enddo
        close (14)

        if (Tdomain%logicD%run_echo .and. Tdomain%Mpi_var%my_rank ==0) then

            call semname_read_mesh_station_echo(fnamef)
            open(94,file=fnamef, status="unknown")
            write (94,*) Tdomain%n_receivers
            write (94,*)  "For any receivers, listed x and z coordinates"
            do i = 0, Tdomain%n_receivers-1
                write (94,*) Tdomain%sReceiver(i)%Xrec, Tdomain%sReceiver(i)%Zrec
            enddo
            close (94)
        endif
    endif




end subroutine read_mesh
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
