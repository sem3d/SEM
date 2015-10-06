!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module mesh3d
    use sdomain
    implicit none
contains

    subroutine  read_interface(inter, fid, pfx)
        use hdf5
        use sem_hdf5
        implicit none
        type(inter_num), intent(inout) :: inter
        integer(HID_T), intent(in) :: fid
        character(len=*), intent(in) :: pfx
        !
        integer, dimension(:,:), allocatable :: itemp2

        call read_attr_int(fid, "n_"//trim(pfx)//"_faces",    inter%surf0%n_faces)
        call read_attr_int(fid, "n_"//trim(pfx)//"_edges",    inter%surf0%n_edges)
        call read_attr_int(fid, "n_"//trim(pfx)//"_vertices", inter%surf0%n_vertices)

        inter%surf1%n_faces = inter%surf0%n_faces
        inter%surf1%n_edges = inter%surf0%n_edges
        inter%surf1%n_vertices = inter%surf0%n_vertices

        call allocate_interface(inter)

        if (inter%surf0%n_faces/=0) then
            call read_dataset(fid, trim(pfx)//"_faces", itemp2)
            write(*,*) pfx//" Faces:", size(itemp2,1), size(itemp2,2)
            inter%surf0%if_faces = itemp2(1,:)
            inter%surf1%if_faces = itemp2(2,:)
            deallocate(itemp2)
            call read_dataset(fid, trim(pfx)//"_orient", itemp2)
            inter%surf0%if_norm = itemp2(1,:)
            inter%surf1%if_norm = itemp2(2,:)
            deallocate(itemp2)
        end if
        if (inter%surf0%n_edges/=0) then
            call read_dataset(fid, trim(pfx)//"_edges", itemp2)
            write(*,*) pfx//" Edges:", size(itemp2,1), size(itemp2,2)
            inter%surf0%if_edges = itemp2(1,:)
            inter%surf1%if_edges = itemp2(2,:)
            deallocate(itemp2)
        end if
        if (inter%surf0%n_vertices/=0) then
            call read_dataset(fid, trim(pfx)//"_vertices", itemp2)
            write(*,*) pfx//" Vertices:", size(itemp2,1), size(itemp2,2)
            inter%surf0%if_vertices = itemp2(1,:)
            inter%surf1%if_vertices = itemp2(2,:)
            deallocate(itemp2)
        end if
    end subroutine read_interface

subroutine read_mesh_file_h5(Tdomain)
    use sdomain
    use mpi
    use hdf5
    use sem_hdf5
    use sem_c_bindings
    use displayCarvalhol
    use semdatafiles, only : MAX_FILE_SIZE
    use constants
    implicit none
    !
    type(domain), intent(inout) :: Tdomain
    integer :: rg, num_processors
    integer :: i, j, k, mat, nod
    logical :: neumann_log
    !
    integer(HID_T) :: fid, proc_id
    integer :: hdferr, ierr
    integer :: code
    integer, allocatable, dimension(:,:) :: itemp2, itemp2b
    integer, allocatable, dimension(:)   :: itemp, itempb
    real,    allocatable, dimension(:,:) :: rtemp2
    !real,    allocatable, dimension(:)   :: rtemp
    character(len=10) :: proc_grp
    integer, allocatable, dimension(:)   :: nb_elems_per_proc
    character(Len=MAX_FILE_SIZE) :: fname
    double precision, dimension(:), allocatable :: tempGlobMin, tempGlobMax
    !
    rg = Tdomain%rank
    !
    call init_hdf5()
    !
    fname = trim(adjustl(Tdomain%mesh_file))//".h5"
    if (sem_check_file(trim(fname))==0) then
        write (*,*) "Process ",rg, " can't open his mesh_file ", trim(fname)
        stop
    endif

    call h5fopen_f(trim(fname), H5F_ACC_RDONLY_F, fid, hdferr)
    !
    !-- Reading mesh properties
    call read_attr_int (fid, "ndims", Tdomain%n_dime)
    Tdomain%logicD%solid_fluid = .false.
    Tdomain%logicD%all_fluid = .false.
    neumann_log = .false.
    Tdomain%logicD%Neumann_local_present = .false.
    call read_attr_int(fid, "n_processors", num_processors)
    call read_attr_int(fid, "n_materials", Tdomain%n_mat)
    call read_attr_int(fid, "n_elements",  Tdomain%n_elem)
    call read_attr_int(fid, "n_faces",     Tdomain%n_face)
    call read_attr_int(fid, "n_edges",     Tdomain%n_edge)
    call read_attr_int(fid, "n_vertices",  Tdomain%n_vertex)
    !
    ! solid - pml interface
    call init_interface(Tdomain%intSolPml)
    ! fluid - pml interface
    call init_interface(Tdomain%intFluPml)
    ! solid - fluid interface
    call init_interface(Tdomain%SF%intSolFlu)
    call init_interface(Tdomain%SF%intSolFluPml)
    !
    if(Tdomain%n_dime /=3)   &
        stop "No general code for the time being: only 3D propagation"
    !
    if(Tdomain%nb_procs/=num_processors) then
        write(*,*) "Mesh file for rank",rg,"has incoherent number of processors"
        write(*,*) "you should check you meshfiles in sem/*"
        stop 1
    endif
    !
    if(neumann_log .neqv. Tdomain%logicD%Neumann) then
        write(*,*) rg,"neumann (input.spec)=",Tdomain%logicD%Neumann
        write(*,*) rg,"neumann (mesh file)=",neumann_log
        stop "Introduction of Neumann B.C.: mesh and input files not in coincidence."
    endif

    !Subdomains allocation and Min/Max bounds initialization
    allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))
    do i = 0, Tdomain%n_mat-1
        allocate(Tdomain%sSubDomain(i)%MinBound(0:Tdomain%n_dime-1))
        allocate(Tdomain%sSubDomain(i)%MaxBound(0:Tdomain%n_dime-1))

        Tdomain%sSubDomain(i)%MinBound = MAX_DOUBLE
        Tdomain%sSubDomain(i)%MaxBound = -MAX_DOUBLE
    end do

    ! Global nodes' coordinates for each proc.
    call read_dataset(fid, "local_nodes", rtemp2)
    Tdomain%n_glob_nodes = size(rtemp2,2)
    allocate (Tdomain%Coord_nodes(0:Tdomain%n_dime-1,0:Tdomain%n_glob_nodes-1))
    Tdomain%Coord_nodes = rtemp2
    deallocate(rtemp2)

    ! Elements (material and solid or fluid and if fluid: Dirichlet boundary?)
    call read_dataset(fid, "material", itemp)

    if (Tdomain%n_elem /= size(itemp)) then
        write(*,*) "N_elem:", Tdomain%n_elem
        write(*,*) "itemp:", size(itemp)
        stop "Incoherent number of elements"
    end if

    allocate(Tdomain%specel(0:Tdomain%n_elem-1))
    allocate (Tdomain%subD_exist(0:Tdomain%n_mat-1))
    Tdomain%subD_exist(:) = .false.

    do i=0,Tdomain%n_elem-1
        call init_element(Tdomain%specel(i))
        Tdomain%specel(i)%mat_index = itemp(i+1)
        Tdomain%specel(i)%OUTPUT = .true.
        Tdomain%subD_exist(itemp(i+1)) = .true.
    end do

    deallocate(itemp)

    ! Read elements definitions
    ! n_nodes : number of control nodes (8 or 27)
    call read_dataset(fid, "elements", itemp2)
    Tdomain%n_nodes = size(itemp2,1)
    do i = 0, Tdomain%n_elem-1
        allocate(Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1))
        Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1) = itemp2(:,i+1)
        mat = Tdomain%specel(i)%mat_index

        do j = 0, Tdomain%n_nodes-1
            nod = Tdomain%specel(i)%Control_Nodes(j)
            do k = 0, Tdomain%n_dime-1
                !Max
                if(Tdomain%Coord_nodes(k, nod) > Tdomain%sSubDomain(mat)%MaxBound(k)) &
                Tdomain%sSubDomain(mat)%MaxBound(k) = Tdomain%Coord_nodes(k, nod)
                !Min
                if(Tdomain%Coord_nodes(k, nod) < Tdomain%sSubDomain(mat)%MinBound(k)) &
                Tdomain%sSubDomain(mat)%MinBound(k) = Tdomain%Coord_nodes(k, nod)
            end do
        end do
    end do

    allocate(tempGlobMin(0:Tdomain%n_dime -1))
    allocate(tempGlobMax(0:Tdomain%n_dime -1))
    do mat = 0, Tdomain%n_mat - 1
        do k = 0, Tdomain%n_dime -1
            call MPI_ALLREDUCE(Tdomain%sSubDomain(mat)%MinBound(k),tempGlobMin(k),1, &
                MPI_DOUBLE_PRECISION,MPI_MIN,Tdomain%communicateur,code)
            call MPI_ALLREDUCE(Tdomain%sSubDomain(mat)%MaxBound(k),tempGlobMax(k),1, &
                MPI_DOUBLE_PRECISION,MPI_MAX,Tdomain%communicateur,code)
            Tdomain%sSubDomain(mat)%MinBound(k) = tempGlobMin(k)
            Tdomain%sSubDomain(mat)%MaxBound(k) = tempGlobMax(k)
        end do
    end do
    deallocate(tempGlobMin)
    deallocate(tempGlobMax)

    deallocate(itemp2)

    ! Faces and elements properties related to faces
    call read_dataset(fid, "faces", itemp2)
    call read_dataset(fid, "faces_def", itemp2b)

    allocate(Tdomain%sFace(0:Tdomain%n_face-1))
    do i=0,Tdomain%n_face-1
        call init_face(Tdomain%sFace(i))
        Tdomain%sFace(i)%inodes = itemp2b(:,i+1)
    enddo
    do i = 0, Tdomain%n_elem-1
        Tdomain%specel(i)%Near_Faces(0:5) = itemp2(:,i+1)
    enddo
    deallocate(itemp2, itemp2b)

    ! Edges
    call read_dataset(fid, "edges", itemp2)
    call read_dataset(fid, "edges_def", itemp2b)

    allocate(Tdomain%sEdge(0:Tdomain%n_edge-1))
    do i=0,Tdomain%n_edge-1
        call init_edge(Tdomain%sEdge(i))
        Tdomain%sEdge(i)%inodes = itemp2b(:,i+1)
    enddo
    do i = 0, Tdomain%n_elem-1
        Tdomain%specel(i)%Near_Edges(0:11) = itemp2(:,i+1)
    enddo
    deallocate(itemp2, itemp2b)

    ! Vertices
    call read_dataset(fid, "vertices", itemp2)

    allocate(Tdomain%sVertex(0:Tdomain%n_vertex-1))
    do i=0,Tdomain%n_vertex-1
        call init_vertex(Tdomain%sVertex(i))
    enddo
    do i = 0,Tdomain%n_elem-1
        Tdomain%specel(i)%Near_Vertices(0:7) = itemp2(:,i+1)
    enddo
    deallocate(itemp2)

    call read_interface(Tdomain%intSolPml, fid, "spml")
    call read_interface(Tdomain%intFluPml, fid, "fpml")
    ! Solid-fluid properties, eventually
    call read_interface(Tdomain%SF%intSolFlu, fid, "sf")
    call read_interface(Tdomain%SF%intSolFluPml, fid, "sfpml")
    Tdomain%logicD%SF_local_present = ( &
        Tdomain%SF%intSolFlu%surf0%n_faces + &
        Tdomain%SF%intSolFlu%surf0%n_edges + &
        Tdomain%SF%intSolFlu%surf0%n_vertices + &
        Tdomain%SF%intSolFluPml%surf0%n_faces + &
        Tdomain%SF%intSolFluPml%surf0%n_edges + &
        Tdomain%SF%intSolFluPml%surf0%n_vertices) /= 0

    ! Neumann B.C. properties, eventually
    if(Tdomain%logicD%Neumann_local_present)then
        ! Neumann properties
        ! Neumann faces
        call read_dataset(fid, "neu_face_near_edges", itemp2)
        ! Edges and their orientation, for each Neumann face
        call read_dataset(fid, "neu_face_near_edges_orient", itemp2b)
        Tdomain%Neumann%Neu_n_faces = size(itemp2,2)
        allocate(Tdomain%Neumann%Neu_face(0:Tdomain%Neumann%Neu_n_faces-1))
        do i = 0, Tdomain%Neumann%Neu_n_faces-1
            Tdomain%Neumann%Neu_face(i)%Near_Edges(0:3) = itemp2(:,i+1)
            Tdomain%Neumann%Neu_face(i)%Near_Edges_Orient(0:3) = itemp2b(:,i+1)
        end do
        deallocate(itemp2, itemp2b)
        ! Vertices for each Neumann face
        call read_dataset(fid, "neu_face_near_vertices", itemp2)
        ! associated face
        call read_dataset(fid, "neu_face_glob_interface", itemp)
        do i = 0, Tdomain%Neumann%Neu_n_faces-1
            Tdomain%Neumann%Neu_face(i)%Near_Vertices(0:3) = itemp2(:,i+1)
            Tdomain%Neumann%Neu_face(i)%Face = itemp(i+1)
        end do
        deallocate(itemp, itemp2)
        ! Neumann edges
        call read_dataset(fid, "neu_edge_glob_interface", itemp)
        Tdomain%Neumann%Neu_n_edges = size(itemp,1)
        allocate(Tdomain%Neumann%Neu_edge(0:Tdomain%Neumann%Neu_n_edges-1))
        do i = 0, Tdomain%Neumann%Neu_n_edges-1
            Tdomain%Neumann%Neu_edge(i)%Edge = itemp(i+1)
        end do
        deallocate(itemp)
        ! Neumann vertices
        call read_dataset(fid,"neu_vertex_glob_interface", itemp)
        Tdomain%Neumann%Neu_n_vertices = size(itemp,1)
        allocate(Tdomain%Neumann%Neu_vertex(0:Tdomain%Neumann%Neu_n_vertices-1))
        do i = 0, Tdomain%Neumann%Neu_n_vertices-1
            Tdomain%Neumann%Neu_vertex(i)%Vertex = itemp(i+1)
        end do
        deallocate(itemp)
    end if

    ! Interproc communications
    Tdomain%tot_comm_proc = 0
    call read_attr_int(fid, "tot_comm_proc", Tdomain%tot_comm_proc)
    allocate (Tdomain%sComm(0:Tdomain%nb_procs-1))
    do i = 0,Tdomain%tot_comm_proc-1
        write(proc_grp,"(a,I4.4)") "Proc", i
        call h5gopen_f(fid, trim(adjustl(proc_grp)), proc_id, hdferr)
        call read_attr_int(proc_id, "proc_dest", Tdomain%sComm(i)%dest)
        call read_attr_int(proc_id, "n_faces", Tdomain%sComm(i)%nb_faces)
        call read_attr_int(proc_id, "n_edges", Tdomain%sComm(i)%nb_edges)
        call read_attr_int(proc_id, "n_vertices", Tdomain%sComm(i)%nb_vertices)
        if(Tdomain%logicD%SF_local_present)then
            call read_attr_int(proc_id, "n_sf_faces", Tdomain%sComm(i)%SF_nf_shared)
            call read_attr_int(proc_id, "n_sf_edges", Tdomain%sComm(i)%SF_ne_shared)
            call read_attr_int(proc_id, "n_sf_vertices", Tdomain%sComm(i)%SF_nv_shared)
        else
            Tdomain%sComm(i)%SF_nf_shared = 0
            Tdomain%sComm(i)%SF_ne_shared = 0
            Tdomain%sComm(i)%SF_nv_shared = 0
        end if
        if(Tdomain%logicD%Neumann_local_present)then
            call read_attr_int(proc_id, "n_neu_edges", Tdomain%sComm(i)%Neu_ne_shared)
            call read_attr_int(proc_id, "n_neu_vertices", Tdomain%sComm(i)%Neu_nv_shared)
        else
            Tdomain%sComm(i)%Neu_ne_shared = 0
            Tdomain%sComm(i)%Neu_nv_shared = 0
        end if
        if(Tdomain%sComm(i)%nb_faces > 0)then
            allocate(Tdomain%sComm(i)%faces(0:Tdomain%sComm(i)%nb_faces-1))
            call read_dataset(proc_id, "faces", itemp)
            do j = 0,Tdomain%sComm(i)%nb_faces-1
                Tdomain%sComm(i)%faces(j)        = itemp(j+1)
            enddo
            deallocate(itemp)
        else
            nullify(Tdomain%sComm(i)%faces)
        endif
        if(Tdomain%sComm(i)%nb_edges > 0)then
            call read_dataset(proc_id, "edges", itemp)
            allocate(Tdomain%sComm(i)%edges(0:Tdomain%sComm(i)%nb_edges-1))
            do j = 0,Tdomain%sComm(i)%nb_edges-1
                Tdomain%sComm(i)%edges(j) = itemp(j+1)
            enddo
            deallocate(itemp)
        else
            nullify(Tdomain%sComm(i)%edges)
            nullify(Tdomain%sComm(i)%orient_edges)
        endif
        if(Tdomain%sComm(i)%nb_vertices > 0)then
            call read_dataset(proc_id, "vertices", itemp)
            allocate(Tdomain%sComm(i)%vertices(0:Tdomain%sComm(i)%nb_vertices-1))
            do j = 0,Tdomain%sComm(i)%nb_vertices-1
                Tdomain%sComm(i)%vertices(j) = itemp(j+1)
            enddo
            deallocate(itemp)
        endif
        if(Tdomain%logicD%SF_local_present)then
            if(Tdomain%sComm(i)%SF_nf_shared > 0)then
                call read_dataset(proc_id, "sf_faces", itemp)
                allocate(Tdomain%sComm(i)%SF_faces_shared(0:Tdomain%sComm(i)%SF_nf_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_nf_shared-1
                    Tdomain%sComm(i)%SF_faces_shared(j) = itemp(j+1)
                enddo
                deallocate(itemp)
            else
                nullify(Tdomain%sComm(i)%SF_faces_shared)
            endif
            if(Tdomain%sComm(i)%SF_ne_shared > 0)then
                call read_dataset(proc_id,"sf_edges", itemp)
                call read_dataset(proc_id,"sf_edges_map", itempb)
                allocate(Tdomain%sComm(i)%SF_edges_shared(0:Tdomain%sComm(i)%SF_ne_shared-1))
                allocate(Tdomain%sComm(i)%SF_mapping_edges_shared(0:Tdomain%sComm(i)%SF_ne_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_ne_shared-1
                    Tdomain%sComm(i)%SF_edges_shared(j) = itemp(j+1)
                    Tdomain%sComm(i)%SF_mapping_edges_shared(j) = itempb(j+1)
                enddo
                deallocate(itemp, itempb)
            else
                nullify(Tdomain%sComm(i)%SF_edges_shared)
                nullify(Tdomain%sComm(i)%SF_mapping_edges_shared)
            endif
            if(Tdomain%sComm(i)%SF_nv_shared > 0)then
                call read_dataset(proc_id, "sf_vertices", itemp)
                allocate(Tdomain%sComm(i)%SF_vertices_shared(0:Tdomain%sComm(i)%SF_nv_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_nv_shared-1
                    Tdomain%sComm(i)%SF_vertices_shared(j) = itemp(j+1)
                enddo
                deallocate(itemp)
            else
                nullify(Tdomain%sComm(i)%SF_vertices_shared)
            endif
        end if
        if(Tdomain%logicD%Neumann_local_present)then
            if(Tdomain%sComm(i)%Neu_ne_shared > 0)then
                call read_dataset(proc_id,"neu_edges", itemp)
                call read_dataset(proc_id,"neu_edges_map", itempb)
                allocate(Tdomain%sComm(i)%Neu_edges_shared(0:Tdomain%sComm(i)%Neu_ne_shared-1))
                allocate(Tdomain%sComm(i)%Neu_mapping_edges_shared(0:Tdomain%sComm(i)%Neu_ne_shared-1))
                do j = 0,Tdomain%sComm(i)%Neu_ne_shared-1
                    Tdomain%sComm(i)%Neu_edges_shared(j) = itemp(j+1)
                    Tdomain%sComm(i)%Neu_mapping_edges_shared(j) = itempb(j+1)
                enddo
                deallocate(itemp, itempb)
            else
                nullify(Tdomain%sComm(i)%Neu_edges_shared)
                nullify(Tdomain%sComm(i)%Neu_mapping_edges_shared)
            endif
            if(Tdomain%sComm(i)%Neu_nv_shared > 0)then
                call read_dataset(proc_id,"neu_vertices", itemp)
                allocate(Tdomain%sComm(i)%Neu_vertices_shared(0:Tdomain%sComm(i)%Neu_nv_shared-1))
                do j = 0,Tdomain%sComm(i)%Neu_nv_shared-1
                    Tdomain%sComm(i)%Neu_vertices_shared(j) = itemp(j+1)
                enddo
                deallocate(itemp)
            else
                nullify(Tdomain%sComm(i)%Neu_vertices_shared)
            endif
        else
            nullify(Tdomain%sComm(i)%Neu_edges_shared)
            nullify(Tdomain%sComm(i)%Neu_mapping_edges_shared)
            nullify(Tdomain%sComm(i)%Neu_vertices_shared)
        end if

    end do

    !write(*,*) "Mesh read correctly for proc #", rg
    if(rg == 0)then
        if(Tdomain%logicD%solid_fluid)then
            write(*,*) "  --> Propagation in solid-fluid media."
        else if(Tdomain%logicD%all_fluid)then
            write(*,*) "  --> Propagation in fluid media."
        else
            write(*,*) "  --> Propagation in solid media."
        end if
    end if
    allocate(nb_elems_per_proc(0:8*((Tdomain%nb_procs-1)/8+1)))
    nb_elems_per_proc = 0
    call MPI_Gather(Tdomain%n_elem, 1, MPI_INTEGER, nb_elems_per_proc, 1, MPI_INTEGER, 0, Tdomain%communicateur, ierr)
    if (rg==0) then
        write(*,*) "Mesh read correctly, elements per proc:"
        do i=0,Tdomain%nb_procs-1,8
            write(*,'(I5.5,a,8I6)') i, ":", nb_elems_per_proc(i:i+7)
        end do
    end if
    deallocate(nb_elems_per_proc)
    !write(*,*) rg, "NFACES=", Tdomain%n_face

end subroutine read_mesh_file_h5


end module mesh3d

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
