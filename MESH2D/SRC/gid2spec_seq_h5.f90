subroutine gid2spec_seq(fnamout,log_created)
    use hdf5
    use sem_hdf5
    implicit none

    logical :: logic_1, logic_2, super_object_present
    logical, intent(in)  :: log_created
    integer, parameter :: N_MAX_POINTS = 10000
    integer :: ndim_space, nglob_ctlpt, n_mat, n_elem, n_nods, n_line, n_vertex_on_interface, idummy
    integer :: i, j, n, ndim_stor, ndim_out
    integer :: i_ex,iunit,icount
    integer :: k, j1, j2, k1, k2, nn, n_vertex, n_faces
    integer :: n_MAX_faces, n_MAX_vertex
    integer, external :: max_value
    integer, dimension (:), allocatable ::  Vertex_to_glob, N_Valid_Vertex
    integer, dimension (:), allocatable ::  ind_mat,n_elem_mat
    integer, dimension (:,:), allocatable :: Ipointer,Build_Faces, El_to_Face, El_to_Vertex, Face_to_El
    integer, dimension (:,:), allocatable :: Face_to_Vertex, Face_to_El_What

    real, dimension (:,:), allocatable :: Gcoord
    character (len=10) :: answer
    character (len=20) :: fnamef
    character (len=50) :: name_cubit
    character (len=50), intent(in) :: fnamout

    ! Super Object variables
    logical , dimension (:), allocatable :: So_vertex
    logical, dimension (:,:), allocatable :: Super_Object

    integer :: low_material, up_material, n_so_vertices, n_faces_super_object
    integer, dimension (1:N_MAX_POINTS) :: Aus_int
    integer, dimension (:), allocatable :: Super_Object_Coherency
    integer, dimension (:,:), allocatable :: Super_Object_to_Face,Super_Object_Vertex, Super_Object_Face_to_Vertex
    integer, dimension (:,:), allocatable :: Nodes_on_boundary, Super_Object_Which_Face
    integer, dimension (:,:), allocatable :: Super_Object_UP_Face_to_Vertex, Super_Object_DOWN_Face_to_Vertex

    real, parameter :: tolerance = 1e-3
    real :: dist
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! This Program reads a GID output file and translates
    ! it into an input file for 2DSPEC for sequential computations
    ! It allows for super object to be analyzed
    ! For this aim the initial super-object is just a fault
    !
    ! Gaetano Festa, last modification 06/10/2005
    ! The parameter variable N_MAX_POINTS is connected to the maximum length of the
    ! input boundary file
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    print*,"    -> File creation"

    super_object_present = .false.
    if (super_object_present) then
        write (*,*) " Introduce the name of the Boundary name"
        read (*,*) fnamef
        open (12,file=fnamef,status="old")
        n_vertex_on_interface = 0
        do i = 1,N_MAX_POINTS
            read (12,*, end=80) Aus_int (i)
            n_vertex_on_interface  = n_vertex_on_interface  + 1
        enddo
80      close(12)

        write (*,*) "Number of points on the super_object interface", n_vertex_on_interface
        allocate (Nodes_on_Boundary(0:n_vertex_on_interface -1,0:1))
        Nodes_on_Boundary (0:n_vertex_on_interface -1,0) = Aus_int(1:n_vertex_on_interface) - 1
    else
        n_vertex_on_interface = 0
    endif

    if(.not. log_created)then
        fnamef="temp.txt"

        open (11,file=fnamef,status="old")
        read (11,*)                 ! Name of the mesh
        read (11,*)    ndim_space   ! Dimension
        read (11,*)                 ! Element Type
        read (11,*)    n_nods       ! Number of control points for each element
        read (11,*)    nglob_ctlpt  ! Number of points
        read (11,*)    n_elem       ! Number of elements
        read (11,*)    n_mat        ! Number of materials

        write(*,*) '     Nb of Points, Elements, Material:',nglob_ctlpt, n_elem, n_mat
        n_line = 1


        read (11,*)
        read (11,*)  ! Begin to read Coordinates
        allocate (Gcoord (0:ndim_space-1,0:nglob_ctlpt-1+n_vertex_on_interface))

        do i = 0,nglob_ctlpt-1
            read (11,*)  idummy,(Gcoord(j,i),j=0,ndim_space-1)
        enddo
        read (11,*)  ! end coordinates
        read (11,*)

        read (11,*)  ! Begin to read Elements
        ndim_stor = n_nods+1

        if ( ndim_space ==2 ) then
            ndim_out = ndim_stor + 8
        else
            write(*,*) "Warning 3D"
        endif
        allocate (Ipointer(0:ndim_out-1,0:n_elem-1))

        do n = 0,n_elem-1
            read (11,*)  idummy,(Ipointer(j,n),j=0,ndim_stor-1)
        enddo
        read (11,*)  ! endElements

        close (11)

    else   ! Cubit 2D file
        print*,"Name of the Cubit file? "
        read*, name_cubit
        iunit = 10
        open(iunit,file=name_cubit,status="old",iostat=i_ex)
        if(i_ex > 0) stop " Cubit file does not exist"
        write(*,*)
        write(*,*) "  --> Number of control nodes per element: "
        read(*,*) n_nods
        if(n_nods /= 4) stop "In mesh_properties:    &
            & Only elements with 4 nodes for the time being."

        call lec1_cubit(iunit,nglob_ctlpt,n_mat)
        allocate(ind_mat(n_mat),n_elem_mat(n_mat))
        call lec2_cubit(iunit,n_mat,n_elem,n_elem_mat,ind_mat)
        read(10,*)   ! Header
        read(10,*)   ! Name of the mesh
        read(10,*)   ! *NODE
        write(*,*) "  --> Number of control points, elements, materials: "
        write(*,*) "    ",nglob_ctlpt, n_elem, n_mat
        write(*,*)
        !- co-ordinates of control points
        allocate(Gcoord(0:1,0:nglob_ctlpt-1))
        do i = 0,nglob_ctlpt-1
            read(10,*) idummy,Gcoord(0,i),Gcoord(1,i)
        enddo

        !- general index for each control point of an element
        ndim_stor = n_nods+1
        ndim_out = ndim_stor + 8
        ndim_space = 2

        allocate(Ipointer(0:ndim_out-1,0:n_elem-1))

        do i = 0,n_mat-1
            read(10,*)    ! *ELEMENTS
            do n = 0,n_elem_mat(i+1)-1
                read(10,*) icount,(Ipointer(j,icount-1),j=0,n_nods-1)
                Ipointer(n_nods,icount-1) = ind_mat(i+1)
            end do
        end do
        !- end of reading Cubit file
        close(10)
        deallocate(ind_mat,n_elem_mat)
        write(*,*)
        write(*,*) "  - END of Cubit-style file READING -"
        write(*,*) "****************************************"
        write(*,*)

    end if


    ! If a super-object is present the number of points should be incremented
    idummy = nglob_ctlpt
    do i = 0, n_vertex_on_interface-1
        Nodes_on_Boundary(i,1) = idummy
        Gcoord(:,nglob_ctlpt+i) = Gcoord(:,Nodes_on_boundary(i,0))
        idummy = idummy + 1
    enddo
    nglob_ctlpt = idummy
    if (super_object_present) write (*,*), "Because of doubling, total number of points changed to   ",nglob_ctlpt

    ! Change the numbering of elements

    do n = 0,n_elem-1
        if (n_nods == 4) then
            do j = 0,3
                Ipointer(j,n) = Ipointer (j,n) - 1
            enddo
        else if (n_nods == 8) then
            do j = 0,7
                Ipointer(j,n) = Ipointer (j,n) - 1
            enddo
        endif
    enddo

    ! Rotate elements
    call rotate (Ipointer,Gcoord,ndim_out,n_elem,nglob_ctlpt,n_nods,4)

    if (super_object_present) then
        allocate (Super_Object(0:1,0:n_elem-1)); allocate (Super_Object_Which_Face(0:1,0:n_elem-1))
        Super_Object = .false.;   Super_Object_Which_Face = -1

        write (*,*) "introduce low_material and up_material indices"
        read(*,*) low_material, up_material
        do n = 0, n_elem -1
            if (Ipointer ( ndim_stor-1,n)== low_material) then
                do j = 0,3
                    j1 = Ipointer(j,n)
                    j2 = Ipointer(j+1,n)
                    if (j == 3) j2 = Ipointer(0,n)
                    logic_1 = .false.
                    do_nn_back01 : do nn = 0,n_elem-1
                        if (Ipointer(ndim_stor-1,nn) == up_material ) then
                            do k = 0,3
                                k1 = Ipointer(k,nn)
                                k2 = Ipointer(k+1,nn)
                                if (k == 3) k2 = Ipointer(0,nn)
                                if (k1 == j1 .and. k2 == j2 ) logic_1 = .true.
                                if (k1 == j2 .and. k2 == j1 ) logic_1 = .true.
                                if (logic_1 ) then
                                    Super_Object(1,n) = .true.
                                    Super_Object(0,nn) = .true.
                                    Super_Object_Which_Face (0,n) = j
                                    Super_Object_Which_Face (0,nn) = k
                                    Super_Object_Which_Face (1,nn) = n
                                    exit do_nn_back01
                                endif
                            enddo
                        endif
                    enddo do_nn_back01
                enddo
            endif
        enddo
        ! Change indices on the lower interface
        do n = 0, n_elem-1
            if (Ipointer ( ndim_stor-1,n)== low_material) then
                do k = 0,3
                    do_nn_interior : do  nn = 0, n_vertex_on_interface -1
                        if (Ipointer (k,n) == Nodes_on_boundary(nn,0)) then
                            Ipointer (k,n) = Nodes_on_boundary(nn,1)
                            exit do_nn_interior
                        endif
                    enddo do_nn_interior
                enddo
            endif
        enddo
    endif

    ! Compute local files

    ! Now the structure of the interior of any domain should be performed in the same way as previously
    ! Chose the effective vertices among the control nodes
    n_MAX_Vertex = 4 * n_elem
    allocate (Vertex_to_Glob(0:n_MAX_Vertex-1))
    allocate (N_valid_Vertex(0:nglob_ctlpt-1))
    N_valid_Vertex = -1

    allocate (El_to_Vertex(0:3,0:n_elem-1))

    n_vertex = 0
    do n = 0,n_elem - 1
        do j = 0,3
            idummy = Ipointer(j,n)
            if (N_valid_Vertex(idummy) < 0) then
                N_valid_Vertex(idummy) = n_vertex
                El_to_Vertex(j,n) = n_vertex
                Vertex_to_Glob(n_vertex) = idummy
                n_vertex = n_vertex + 1
            else
                El_to_Vertex(j,n) = N_Valid_Vertex(idummy)
            endif
        enddo
    enddo

    ! Define faces properties
    n_MAX_faces = 4 * n_elem
    allocate (El_to_face(0:3,0:n_elem-1))
    allocate (Face_to_El(0:1,0:n_MAX_faces-1))
    allocate (Face_to_El_What(0:1,0:n_MAX_faces-1))
    allocate (Face_to_Vertex(0:1,0:n_MAX_faces-1))

    Face_to_El = -1; Face_to_El_What = -1;
    n_faces = 0
    do n = 0, n_elem-1
        do j = 0,3
            j1 = Ipointer(j,n)
            j2 = Ipointer(j+1,n)
            if (j == 3) j2 = Ipointer(0,n)

            logic_1 = .false.
            do_nn_back1 : do nn = 0,n-1
                do k = 0,3
                    k1 = Ipointer(k,nn)
                    k2 = Ipointer(k+1,nn)
                    if (k == 3) k2 = Ipointer(0,nn)
                    if (k1 == j1 .and. k2 == j2 ) logic_1 = .true.
                    if (k1 == j2 .and. k2 == j1 ) logic_1 = .true.
                    if (logic_1 ) then
                        El_to_face (j,n) = El_to_face (k,nn)
                        Face_to_el (1,El_to_face (k,nn) ) = n
                        Face_to_el_What (1,El_to_face (k,nn) ) = j
                        exit do_nn_back1
                    endif

                enddo
            enddo do_nn_back1
            if (.not. logic_1) then
                El_to_face (j,n) = n_faces
                Face_to_el (0,n_faces ) = n
                Face_to_el_What (0,n_faces ) = j
                Face_to_Vertex (0,n_faces) = N_valid_Vertex(j1)
                Face_to_Vertex (1,n_faces) = N_valid_Vertex(j2)
                n_faces = n_faces + 1
            endif
        enddo
    enddo

    do n = 0, n_elem-1
        if (n_nods == 4 ) then
            Ipointer (12,n) = Ipointer (ndim_stor-1,n) - 1
            Ipointer (4:7,n) = El_to_face (0:3,n)
            Ipointer (8:11,n) = El_to_Vertex(0:3,n)
        else if (n_nods == 8 ) then
            Ipointer (16,n) = Ipointer (ndim_stor-1,n) -1
            Ipointer (8:11,n) = El_to_face (0:3,n)
            Ipointer(12:15,n) = El_to_Vertex(0:3,n)
        endif
    enddo

    write(*,*) "     Nb of elmts, faces and vertices: ", n_elem, n_faces, n_vertex
    allocate (Build_faces(0:5,0:n_faces-1))
    do n = 0, n_faces - 1
        Build_faces (0:1,n) = Face_to_el (0:1, n)
        Build_faces (2:3,n) = Face_to_el_What (0:1, n)
        Build_faces (4:5,n) = Face_to_Vertex (0:1, n)
    enddo

    ! Define the ordering of faces
    do n = 0, n_faces - 1
        nn = Build_Faces(0,n)
        k = Build_Faces(2,n)
        if (k == 0 .and. Vertex_to_Glob(Build_Faces(4,n)) == Ipointer(1,nn) ) then
            k1 = Build_Faces(4,n); Build_Faces(4,n) =  Build_Faces(5,n); Build_Faces(5,n)= k1
        else if (k == 1 .and. Vertex_to_Glob(Build_Faces(4,n)) == Ipointer(2,nn) ) then
            k1 = Build_Faces(4,n); Build_Faces(4,n) =  Build_Faces(5,n); Build_Faces(5,n)= k1
        else if (k == 2 .and. Vertex_to_Glob(Build_Faces(4,n)) == Ipointer(2,nn) ) then
            k1 = Build_Faces(4,n); Build_Faces(4,n) =  Build_Faces(5,n); Build_Faces(5,n)= k1
        else if (k == 3 .and. Vertex_to_Glob(Build_Faces(4,n)) == Ipointer(3,nn) ) then
            k1 = Build_Faces(4,n); Build_Faces(4,n) =  Build_Faces(5,n); Build_Faces(5,n)= k1
        endif
    enddo

    ! Define super Object properties
    if (super_object_present) then
        n_faces_super_object = 0
        do n = 0, n_elem -1
            if (Super_Object (0,n)) n_faces_super_object = n_faces_super_object + 1
        enddo
        write (*,*) "Number of faces in the super object are:  ", n_faces_super_object
        idummy = 0
        allocate (Super_Object_to_Face (0:1, 0:n_faces_super_object-1))
        allocate (Super_Object_Coherency(0:n_faces_super_object-1))
        Super_Object_to_Face = -1
        do n = 0, n_elem-1
            if (Super_Object(0,n)) then
                k = Super_Object_Which_Face(0,n)
                Super_Object_to_Face(0,idummy) = El_to_Face(k,n)
                nn = Super_Object_Which_Face(1,n)
                k = Super_Object_Which_Face(0,nn)
                Super_Object_to_Face(1,idummy) = El_to_Face(k,nn)
                Super_Object(1,nn) = .false.
                idummy = idummy + 1
            endif
        enddo
        if (idummy /= n_faces_super_object) then
            write (*,*) "There is a problem in counting super object faces"
            stop
        endif
        if (n_faces_super_object > 0) then
            allocate (Super_Object_UP_Face_to_Vertex (0:1,0:n_faces_super_object-1))
            allocate (Super_Object_DOWN_Face_to_Vertex (0:1,0:n_faces_super_object-1))
            allocate (Super_Object_Face_to_Vertex(0:1,0:n_faces_super_object-1))
        endif
        do n = 0, n_faces_super_object-1
            k = Super_Object_to_Face(0,n)
            Super_Object_UP_Face_to_Vertex (0:1,n) = Build_Faces(4:5,k)
            k = Super_Object_to_Face(1,n)
            Super_Object_DOWN_Face_to_Vertex (0:1,n) = Build_Faces(4:5,k)

            k1 = Super_Object_UP_Face_to_Vertex (0,n)
            k2 = Super_Object_UP_Face_to_Vertex (1,n)
            j1 = Super_Object_DOWN_Face_to_Vertex (0,n)
            j2 = Super_Object_DOWN_Face_to_Vertex (1,n)
            logic_1 = .false.; logic_2 = .false.
            if (k1 == j1 .or. k1 == j2) logic_1 = .true.
            if (k2 == j2 .or.  k2 == j1) logic_2 = .true.
            if (logic_1) then    ! Termination on k1
                if (k1 == j1) then
                    Super_Object_Coherency(n) = 1
                else
                    Super_Object_Coherency(n) = 0
                endif
            else if (logic_2) then  ! Termination on k2
                if (k2 == j2) then
                    Super_Object_Coherency(n) = 1
                else
                    Super_Object_Coherency(n) = 0
                endif
            else
                Super_Object_Coherency(n) = -1
                j = Vertex_to_Glob(j1)
                k = Vertex_to_Glob(k1)
                dist = (Gcoord (0,j) - Gcoord (0,k))**2 + (Gcoord (1,j) - Gcoord (1,k))**2
                if (dist < tolerance ) Super_Object_Coherency (n) = 1
                j = Vertex_to_Glob (j2)
                dist = (Gcoord (0,j) - Gcoord (0,k))**2 + (Gcoord (1,j) - Gcoord (1,k))**2
                if (dist < tolerance ) Super_Object_Coherency (n) = 0
                if (super_object_coherency(n) == -1) then
                    write (*,*) "There is a problem in defining the coherency"
                    stop
                endif
            endif
        enddo

        allocate (So_vertex(0:n_vertex-1))
        So_Vertex = .false.
        n_so_vertices = 0
        do n = 0, n_faces_super_object - 1
            k1 = Super_Object_UP_Face_to_Vertex (0,n)
            if (.not. So_Vertex (k1)) then
                So_Vertex(k1) = .true.
                Super_Object_Face_to_Vertex (0,n) = n_so_vertices
                n_so_vertices = n_so_vertices + 1
            else
                do_interior1 : do nn = 0, n-1
                    j1 = Super_Object_UP_Face_to_Vertex (0,nn)
                    j2 = Super_Object_UP_Face_to_Vertex (1,nn)
                    if (j1 == k1) then
                        Super_Object_Face_to_Vertex (0,n) = Super_Object_Face_to_Vertex(0,nn)
                        exit do_interior1
                    endif
                    if (j2 == k1) then
                        Super_Object_Face_to_Vertex (0,n) = Super_Object_Face_to_Vertex(1,nn)
                        exit do_interior1
                    endif
                enddo do_interior1
            endif
            k1 = Super_Object_UP_Face_to_Vertex (1,n)
            if (.not. So_Vertex (k1)) then
                So_Vertex(k1) = .true.
                Super_Object_Face_to_Vertex (1,n) = n_so_vertices
                n_so_vertices = n_so_vertices + 1
            else
                do_interior2 : do nn = 0, n-1
                    j1 = Super_Object_UP_Face_to_Vertex (0,nn)
                    j2 = Super_Object_UP_Face_to_Vertex (1,nn)
                    if (j1 == k1) then
                        Super_Object_Face_to_Vertex (1,n) = Super_Object_Face_to_Vertex(0,nn)
                        exit do_interior2
                    endif
                    if (j2 == k1) then
                        Super_Object_Face_to_Vertex (1,n) = Super_Object_Face_to_Vertex(1,nn)
                        exit do_interior2
                    endif
                enddo do_interior2
            endif
        enddo

        if (n_so_vertices > 0) then
            allocate (Super_Object_Vertex(0:1, 0:n_so_vertices -1))
            Super_Object_Vertex = -1
            So_Vertex = .false.
            idummy = 0.
            do n = 0, n_faces_super_object - 1
                k1 = Super_Object_UP_Face_to_Vertex (0,n)
                j1 = Super_Object_DOWN_Face_to_Vertex (0,n)
                if (.not. So_Vertex (k1)) then
                    So_Vertex(k1) = .true.
                    if (k1 == j1 ) then
                        Super_Object_Vertex (0:1, idummy) = -3
                    else
                        Super_Object_Vertex (0,idummy) = k1
                        Super_Object_Vertex (1,idummy) = j1
                    endif
                    idummy = idummy + 1
                endif
                k1 = Super_Object_UP_Face_to_Vertex (1,n)
                j1 = Super_Object_DOWN_Face_to_Vertex (1,n)
                if (.not. So_Vertex (k1)) then
                    So_Vertex(k1) = .true.
                    if (k1 == j1 ) then
                        Super_Object_Vertex (0:1, idummy) = -3
                    else
                        Super_Object_Vertex (0,idummy) = k1
                        Super_Object_Vertex (1,idummy) = j1
                    endif
                    idummy = idummy + 1
                endif
            enddo
        endif
    endif

    ! Write Local Output
    open (12,file=trim(fnamout),status="unknown",form="formatted")

    write (12,202)    ndim_space, " ! Space dimension "
    write (12,202)    nglob_ctlpt, " ! Number of local nodes "

    do i = 0,nglob_ctlpt-1
        write (12,201)  (Gcoord(j,i),j=0,ndim_space-1)
    enddo

    write (12,202)  n_mat, " ! Number of materials "

    ! Artificially created

    write (12,202)  1, " ! Number of line "
    write (12,*) "Write the line number and name"
    write (12,202) 1,"L1"

    write (12,202)   n_elem, "! number of local elements"
    write (12,202)   n_nods, "! number of nodes"
    write (12,203)  "! For any element, listed the reference to the global list of elements, the faces, the material"
    do n = 0,n_elem - 1
        if (n_nods == 4 ) then
            write (12,204) (Ipointer(k,n),k=0,ndim_out-1)
        else if  (n_nods == 8 ) then
            write (12,205) (Ipointer(k,n),k=0,ndim_out-1)
        endif
    enddo
    write (12,203)  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write (12,202)   n_faces, "! number of faces"
    write (12,203)  "! For any face, listed to the elements, what face of the element and the vertices"
    do n = 0,n_faces-1
        write (12,206) (Build_faces(k,n),k=0,5)
    enddo
    write (12,203)  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write (12,202)   n_vertex, "! number of vertices"
    write (12,203)  "! For any vertex, address to the global numbering"
    do n = 0,n_vertex-1
        write (12,*) Vertex_to_Glob(n)
    enddo
    if (super_object_present) then
        write (12,203)  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write (12,202)   0, "! Type of super-object, 0 fault, 1 Mortar"
        write (12,202)   n_faces_super_object , "Number of faces contained in the super-object"
        do n = 0, n_faces_super_object-1
            write (12,207)   Super_Object_to_Face(0:1,n), Super_Object_Face_to_Vertex(0:1,n), Super_Object_Coherency(n)
        enddo
        write (12,203)  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        write (12,202)   n_so_vertices, "Number of vertices contained in the super-object"
        do n = 0, n_so_vertices -1
            write (12,207)   Super_Object_Vertex(0:1,n)
        enddo
    endif
    write (12,203)  "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    write (12,203) " Have a good run"
    close (12)

    deallocate (Vertex_to_Glob)
    deallocate (N_Valid_Vertex)
    deallocate (El_to_Face)
    deallocate (Face_to_El)
    deallocate (Face_to_El_What)
    deallocate (Face_to_Vertex)
    deallocate (Build_faces)
    deallocate (Gcoord)
    deallocate (Ipointer)
    if (super_object_present) then
        deallocate (Super_Object_to_Face)
        deallocate (Super_Object_Coherency)
        if (n_faces_super_object > 0) then
            deallocate (Super_Object_UP_Face_to_Vertex)
            deallocate (Super_Object_DOWN_Face_to_Vertex)
            deallocate (Super_Object_Face_to_Vertex)
        endif
        if (n_so_vertices > 0) deallocate (Super_Object_Vertex)
        deallocate (Super_Object)
    endif


201 format (3e15.7)
202 format (i8,10x,a)
203 format (a)
204 format (13i8)
205 format (17i8)
206 format (6i8)
207 format (5i8)
208 format (2i8)
209 format (i8)

end subroutine gid2spec_seq
! #########################################################
subroutine rotate (Ipointer,Gcoord,ndim_stor,n_elem,nglob_ctlpt,n_nods,nv)

    implicit none

    integer, intent (IN) :: ndim_stor,n_elem,nglob_ctlpt, n_nods,nv
    integer, dimension (0:ndim_stor-1,0:n_elem-1), intent (INOUT) :: Ipointer
    real, dimension (0:1,0:nglob_ctlpt-1), intent (IN) :: Gcoord

    !local variables
    integer :: k,kp,kf,kff,n,i
    integer, dimension (0:3) :: Ivect
    real, dimension (0:3) :: XVect, Zvect
    real :: sign
    real, external :: compute_area
    logical :: logic1

    ! The condition for the coordinates are
    ! for the first point x < x (following); z< z (previous)
    ! for the second z <z (following)
    ! for the third x < x (following)
    ! This should be satisfy the convexity of the element

    do n = 0,n_elem-1
        do k= 0,3
            Xvect(k) = Gcoord (0,Ipointer(k,n))
            Zvect(k) = Gcoord (1,Ipointer(k,n))
        enddo

        sign =  compute_area(Xvect, Zvect)

        if (sign < 0) then
            if (n_nods == 4 ) then
                do i = 0,3
                    Ivect (i) = Ipointer (3-i,n)
                enddo
                Ipointer (0:3,n) = Ivect
            else
                do i = 0,3
                    Ivect (i) = Ipointer (3-i,n)
                enddo
                Ipointer (0:3,n) = Ivect
                do i = 0,3
                    Ivect (i) = Ipointer (7-i,n)
                enddo
                Ipointer (4:7,n) = Ivect
            endif
            do k= 0,3
                Xvect(k) = Gcoord (0,Ipointer(k,n))
                Zvect(k) = Gcoord (1,Ipointer(k,n))
            enddo
        endif

        k_do : do k  = 0,3
            kf = k + 1
            kp = k - 1
            kff = k + 2
            if (k == 3) then
                kf = 0; kff = 1;
            endif
            if (k == 0 ) kp = 3
            if (k == 2) kff = 0
            logic1 = Xvect(k) < Xvect(kf)     .and. Zvect(k) < Zvect(kp)
            logic1 = logic1 .and. Zvect(kf) < Zvect(kff)
            logic1 = logic1 .and. Xvect(kff) > Xvect(kp)
            if (logic1) exit k_do
        enddo k_do

        if (k /= 0 ) then

            if (n_nods == 4 ) then
                Ivect(0:3) = Ipointer (0:3,n)
                Ivect = cshift (Ivect,k)
                Ipointer (0:3,n) = Ivect (0:3)
            else if (n_nods == 8 ) then
                Ivect(0:3) = Ipointer (0:3,n)
                Ivect = cshift (Ivect,k)
                Ipointer (0:3,n) = Ivect (0:3)
                Ivect(0:3) = Ipointer (4:7,n)
                Ivect = cshift (Ivect,k)
                Ipointer (4:7,n) = Ivect (0:3)
                Ivect(0:3) = Ipointer(n_nods+nv:n_nods+nv+3,n)
                Ivect = cshift (Ivect,k)
                Ipointer(n_nods+nv:n_nods+nv+3,n) = Ivect (0:3)
            endif
        endif
    enddo
    return
end subroutine rotate

!#####################################################
real function compute_area (X,z)
    implicit none
    real, dimension (0:3) :: x,z

    compute_area = x(0)*z(1) + x(1)*z(2) + x(2)*z(0) - z(1)*x(2) - z(2)*x(0) - z(0) * x(1)
    return
end function compute_area
! #######################################################
integer function max_value (array,n)
    implicit none

    integer, dimension (0:n-1) ::array
    integer :: n,i, r_max

    r_max = array(0)
    do i  = 1, n-1
        if (r_max < array(i)) r_max = array(i)
    enddo
    max_value = r_max
    return
end function max_value
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
