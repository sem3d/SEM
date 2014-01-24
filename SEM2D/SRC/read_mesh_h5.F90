!>
!!\file read_mesh_h5.F90
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


subroutine read_mesh_h5(tDomain)
    use sdomain
    use semdatafiles
    use HDF5
    use sem_hdf5

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    !
    character (len=MAX_FILE_SIZE) :: fnamef
    integer :: i,j,j_aus,n_aus,k_aus,w_face
    integer , dimension (:), allocatable :: Ipointer
    integer :: n_dime
    !
    integer, allocatable, dimension(:,:) :: itemp2, itemp2b
    integer, allocatable, dimension(:)   :: itemp, itempb
    real,    allocatable, dimension(:,:) :: rtemp2
    integer(HID_T) :: fid, proc_id
    integer :: hdferr, ierr, n_proc_mesh
    ! Read Mesh properties
    ! The reading is local to the grid

    call semname_read_mesh_rank(Tdomain%mesh_file,Tdomain%Mpi_var%my_rank,fnamef)
    write(*,*)"ouverture fichier:",trim(adjustl(fnamef))
    call init_hdf5()
    !
    call h5fopen_f(trim(adjustl(fnamef))//".h5", H5F_ACC_RDONLY_F, fid, hdferr)


    call read_attr_int(fid, "ndim", n_dime)
    call read_attr_int(fid, "n_processors", n_proc_mesh)
    call read_attr_int(fid, "n_materials", Tdomain%n_mat)
    call read_attr_int(fid, "n_elements",  Tdomain%n_elem)
    call read_attr_int(fid, "n_edges",     Tdomain%n_face)
    call read_attr_int(fid, "n_vertices",  Tdomain%n_vertex)


    if (n_dime /= 2) then
        write (*,*) "A dimension different from 2 is not yet taken into account"
        stop
    endif
    if (n_proc_mesh /= Tdomain%MPi_var%n_proc) then
        write (*,*) "The number of selected processors is not the same as chosen for the mesh partition"
        write (*,*) " For processor " , Tdomain%Mpi_var%my_rank, " the two values are : " ,  Tdomain%Mpi_var%n_proc, n_aus
        stop
    endif
    ! Global nodes' coordinates for each proc.
    !
    call read_dataset(fid, "nodes", rtemp2)
    Tdomain%n_glob_nodes = size(rtemp2,2)
    allocate (Tdomain%Coord_nodes(0:1,0:Tdomain%n_glob_nodes-1))
    Tdomain%Coord_nodes = rtemp2
    deallocate(rtemp2)

    ! Elements (material and solid or fluid) and pml direction/flag
    !
    call read_dataset(fid, "material", itemp2)
    if (Tdomain%n_elem /= size(itemp2,2)) then
        write(*,*) "N_elem:", Tdomain%n_elem
        write(*,*) "itemp:", size(itemp2,1), size(itemp2,2)
        stop "Incoherent number of elements"
    end if

    allocate (Tdomain%specel(0:Tdomain%n_elem-1))
    do i=0,Tdomain%n_elem-1
        call init_element(Tdomain%specel(i))
        Tdomain%specel(i)%mat_index = itemp2(1,i+1)
        Tdomain%specel(i)%OUTPUT = .true.
    enddo

    ! Read element's definitions
    ! n_nodes : number of control nodes (4 or 8)
    call read_dataset(fid, "elements", itemp2)
    Tdomain%n_nodes = size(itemp2,1)
    do i = 0, Tdomain%n_elem-1
        allocate(Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1))
        Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1) = itemp2(:,i+1)
    enddo
    deallocate(itemp2)

    ! Near face
    call read_dataset(fid, "edges", itemp2)
    do i = 0,  Tdomain%n_elem - 1
        Tdomain%specel(i)%Near_Face(0:3) = itemp2(:,i+1)
        ! TODO: compute %sFace(i)%near_element / which_face... near_vertex
    enddo
    deallocate(itemp2)

    ! Near Vertex
    call read_dataset(fid, "vertices", itemp2)
    do i = 0,  Tdomain%n_elem - 1
        Tdomain%specel(i)%Near_Vertex(0:3) = itemp2(:,i+1)
    enddo
    deallocate(itemp2)


    allocate (Tdomain%sFace(0:Tdomain%n_face-1))
    call read_dataset(fid, "faces_elem", itemp2)
    do i = 0,  Tdomain%n_face - 1
        Tdomain%sface(i)%Near_element(0:1) = itemp2(:,i+1)
    enddo
    deallocate(itemp2)

    call read_dataset(fid, "faces_which", itemp2)
    do i = 0,  Tdomain%n_face - 1
        Tdomain%sface(i)%Which_face(0:1) = itemp2(:,i+1)
    enddo
    deallocate(itemp2)

    call read_dataset(fid, "faces_vertex", itemp2)
    do i = 0,  Tdomain%n_face - 1
        Tdomain%sface(i)%Near_Vertex(0:1) = itemp2(:,i+1)
    enddo
    deallocate(itemp2)


    allocate (Tdomain%svertex(0:Tdomain%n_vertex-1))
    do i = 0, Tdomain%n_vertex-1
        !!! TODO shape8
        Tdomain%svertex(i)%Glob_numbering = i
    enddo

    !! TODO COMMUNICATIONS

    call coherency_mesh_h5(Tdomain)

    call read_material_file(Tdomain)

    Tdomain%logicD%save_trace = .false.

end subroutine read_mesh_h5


subroutine coherency_mesh_h5(Tdomain)
    use sdomain
    implicit none
    type (domain), intent (INOUT) :: Tdomain

    integer :: i, v0, v1, e0, e1, w_face0, w_face1

    ! Define a coherency for the face ordering and numbering
    do i =0, Tdomain%n_face-1
        v0 = Tdomain%sFace(i)%Near_Vertex(0)
        v1 = Tdomain%sFace(i)%Near_Vertex(1)
        e0 = Tdomain%sFace(i)%Near_Element(0)
        w_face0 = Tdomain%sFace(i)%Which_face(0)

        if (w_face0 == 0 .or. w_face0 == 3 ) then
            if (Tdomain%specel(e0)%Control_Nodes(0) /= Tdomain%sVertex(v0)%Glob_numbering)  then
                Tdomain%sFace(i)%Near_vertex(0) = v1
                Tdomain%sFace(i)%Near_vertex(1) = v0
            endif
        else if (w_face0 == 1 ) then
            if (Tdomain%specel(e0)%Control_Nodes(1) /= Tdomain%sVertex(v0)%Glob_numbering)  then
                Tdomain%sFace(i)%Near_vertex(0) = v1
                Tdomain%sFace(i)%Near_vertex(1) = v0
            endif
        else
            if (Tdomain%specel(e0)%Control_Nodes(3) /= Tdomain%sVertex(v0)%Glob_numbering)  then
                Tdomain%sFace(i)%Near_vertex(0) = v1
                Tdomain%sFace(i)%Near_vertex(1) = v0
            endif
        endif

        v0 = Tdomain%sFace(i)%Near_Vertex(0)
        v1 = Tdomain%sFace(i)%Near_Vertex(1)
        e1 = Tdomain%sFace(i)%Near_Element(1)
        Tdomain%sFace(i)%coherency = .false.
        if (e1 > -1 ) then
            w_face1 = Tdomain%sFace(i)%Which_face(1)
            if (w_face1 == 0 .or. w_face1 == 3 ) then
                if (Tdomain%specel(e1)%Control_Nodes(0) == Tdomain%sVertex(v0)%Glob_numbering)  &
                    Tdomain%sFace(i)%coherency = .true.
            else if  (w_face1 == 1 ) then
                if (Tdomain%specel(e1)%Control_Nodes(1) == Tdomain%sVertex(v0)%Glob_numbering)   &
                    Tdomain%sFace(i)%coherency = .true.
            else
                if (Tdomain%specel(e1)%Control_Nodes(3) == Tdomain%sVertex(v0)%Glob_numbering) &
                    Tdomain%sFace(i)%coherency = .true.
            endif
        endif
    enddo
end subroutine coherency_mesh_h5

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
