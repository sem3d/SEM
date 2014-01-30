module msnapshots
    use sdomain
    use hdf5
    use sem_hdf5
    use semdatafiles
    use mpi
    use constants
    !use mfields
    use orientation
    implicit none
contains

    subroutine compute_saved_elements(Tdomain, irenum, nnodes)
        type (domain), intent (INOUT):: Tdomain
        integer, allocatable, dimension(:), intent(out) :: irenum ! maps Iglobnum to file node number
        integer, intent(out) :: nnodes
        integer :: n, i, k, ngllx, ngllz, ig, gn, ne

        allocate(irenum(0:Tdomain%n_glob_points-1))

        irenum = -1
        ig = 0
        ne = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            ne = ne + 1
            do k = 0,ngllz - 1
                do i = 0,ngllx - 1
                    gn = Tdomain%specel(n)%Iglobnum(i,k)
                    if (irenum(gn) == -1) then
                        irenum(gn) = ig
                        ig = ig + 1
                    end if
                end do
            end do
        end do
        nnodes = ig
    end subroutine compute_saved_elements

    subroutine create_dir_sorties(Tdomain, rg, isort)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer,intent(in) :: isort, rg
        character(Len=MAX_FILE_SIZE) :: temp
        character(len=MAX_FILE_SIZE+10) :: creer_dir
        integer :: code

        if (rg==0) then
            call semname_snap_result_dir(isort, temp)
            creer_dir = "mkdir -p "//temp
            call system(creer_dir)
        end if
        call mpi_barrier(Tdomain%communicateur, code)
    end subroutine create_dir_sorties

    !>
    !! Ecrit la geometrie pour les sorties dans un fichier HDF5
    !! Le format est destine a etre relu par paraview / xdmf
    !<
    subroutine write_snapshot_geom(Tdomain, rg)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: rg
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        integer :: hdferr, code
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes

        call init_hdf5()
        if (rg==0) then
            call system("mkdir -p " // path_results)
        end if
        call mpi_barrier(Tdomain%communicateur, code)
        call semname_snap_geom_file(rg, fnamef)

        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

        call compute_saved_elements(Tdomain, irenum, nnodes)

        call write_global_nodes(Tdomain, fid, irenum, nnodes)
        
        call write_elem_connectivity(Tdomain, fid, irenum)

        call write_constant_fields(Tdomain, fid, irenum, nnodes)

        call h5fclose_f(fid, hdferr)

        if (rg==0) call write_master_xdmf(Tdomain)
    end subroutine write_snapshot_geom


    subroutine write_global_nodes(Tdomain, fid, irenum, nnodes)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, intent(in) :: nnodes
        integer, dimension(:), intent(in), allocatable :: irenum
        !
        real, dimension(:,:), allocatable :: nodes
        integer(HSIZE_T), dimension(2) :: dims
        integer :: n
        integer(HID_T) :: nodes_id
        integer :: hdferr
        
        allocate(nodes(2,0:nnodes-1))
        do n = 0, Tdomain%n_glob_points-1
            if (irenum(n)>=0) then
                nodes(:,irenum(n)) = Tdomain%GlobCoord(:,n)
            end if
        end do

        dims(1) = 2
        dims(2) = nnodes
        call create_dset_2d(fid, "Nodes", H5T_IEEE_F64LE, dims(1), dims(2), nodes_id)
        call h5dwrite_f(nodes_id, H5T_NATIVE_DOUBLE, nodes, dims, hdferr)
        call h5dclose_f(nodes_id, hdferr)
        deallocate(nodes)
    end subroutine write_global_nodes

    subroutine write_elem_connectivity(Tdomain, fid, irenum)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, dimension(:), intent(in), allocatable :: irenum
        !
        integer(HID_T) :: elem_id, mat_id, ngll_id, globnum_id, elem_num_id
        integer :: ngllx, ngllz
        integer(HSIZE_T), dimension(2) :: dims
        integer, dimension(:,:), allocatable :: data
        integer, dimension(:), allocatable :: mat, elem_num, iglobnum
        integer, dimension(2,0:Tdomain%n_elem-1) :: ngll
        integer :: count, ig, nglobnum
        integer :: i, k, n, nb_elem
        integer :: hdferr

        ! First we count the number of hexaedrons
        count = 0
        nglobnum = 0
        k = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll(1,k) = Tdomain%specel(n)%ngllx
            ngll(2,k) = Tdomain%specel(n)%ngllz
            ! Max number of global points (can count elem vertices twice)
            nglobnum = nglobnum + ngll(1,k)*ngll(2,k)
            ! Number of subelements
            count = count+(ngll(1,k)-1)*(ngll(2,k)-1)
            k = k + 1
        enddo
        nb_elem = k
        !! Nombre de points de gauss par element
        call create_dset_2d(fid, "NGLL", H5T_STD_I16LE, 2, nb_elem, ngll_id)
        dims(1) = 2
        dims(2) = nb_elem
        call h5dwrite_f(ngll_id, H5T_NATIVE_INTEGER, ngll, dims, hdferr)
        call h5dclose_f(ngll_id, hdferr)

        Tdomain%n_quad = count
        allocate( data(1:4,0:count-1))
        allocate( mat(0:count-1))
        allocate( elem_num(0:count-1))

        call create_dset_2d(fid, "Elements", H5T_STD_I32LE, 4, count, elem_id)
        call create_dset(fid, "Material", H5T_STD_I32LE, count, mat_id)
        call create_dset(fid, "ElemID", H5T_STD_I32LE, count, elem_num_id)
        !call create_dset(fid, "ElemID", H5T_STD_I32LE, Tdomain%n_elem, elem_num_id)
        call create_dset(fid, "Iglobnum", H5T_STD_I32LE, nglobnum, globnum_id)

        allocate (iglobnum(nglobnum))
        dims(1) = 4
        dims(2) = count
        count = 0
        ig = 1
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            do k = 0,ngllz-2
                do i = 0,ngllx-2
                    data(1,count) = irenum(Tdomain%specel(n)%Iglobnum(i+0,k+0))
                    data(2,count) = irenum(Tdomain%specel(n)%Iglobnum(i+1,k+0))
                    data(3,count) = irenum(Tdomain%specel(n)%Iglobnum(i+1,k+1))
                    data(4,count) = irenum(Tdomain%specel(n)%Iglobnum(i+0,k+1))
                    mat(count) = Tdomain%specel(n)%mat_index
                    elem_num(count) = n
                    count=count+1
                end do
            end do
            do k = 0,ngllz - 1
                do i = 0,ngllx - 1
                    iglobnum(ig) = Tdomain%specel(n)%Iglobnum(i,k)
                    ig = ig + 1
                end do
            end do
        end do
        dims(1) = nglobnum
        dims(2) = 0
        call h5dwrite_f(globnum_id, H5T_NATIVE_INTEGER, iglobnum, dims, hdferr)
        deallocate(iglobnum)
        call h5dclose_f(globnum_id, hdferr)
        call h5dwrite_f(elem_id, H5T_NATIVE_INTEGER, data, dims, hdferr)
        call h5dclose_f(elem_id, hdferr)
        dims(1) = count
        dims(2) = 0
        call h5dwrite_f(mat_id, H5T_NATIVE_INTEGER, mat, dims, hdferr)
        call h5dclose_f(mat_id, hdferr)
        call h5dwrite_f(elem_num_id, H5T_NATIVE_INTEGER, elem_num, dims, hdferr)
        call h5dclose_f(elem_num_id, hdferr)
        deallocate(data)
        !deallocate(elem_num)
    end subroutine write_elem_connectivity

    subroutine save_field_h5(Tdomain, rg, isort)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: rg, isort
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid, displ_id, veloc_id, press_id
        integer(HSIZE_T), dimension(2) :: dims
        real, dimension(:,:),allocatable :: displ, veloc
        real, dimension(:), allocatable :: press
        real, dimension(:,:,:),allocatable :: field_displ, field_veloc
        integer, dimension(:), allocatable :: valence
        integer :: hdferr
        integer :: ngllx, ngllz, idx
        integer :: i, k, n
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes
        

        call create_dir_sorties(Tdomain, rg, isort)
        call semname_snap_result_file(rg, isort, fnamef)

        call compute_saved_elements(Tdomain, irenum, nnodes)


        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

        dims(1) = 3
        dims(2) = nnodes
        call create_dset_2d(fid, "displ", H5T_IEEE_F64LE, 3, nnodes, displ_id)
        call create_dset_2d(fid, "veloc", H5T_IEEE_F64LE, 3, nnodes, veloc_id)
        call create_dset(fid, "pressure", H5T_IEEE_F64LE, nnodes, press_id)

        allocate(displ(0:2,0:nnodes-1))
        allocate(veloc(0:2,0:nnodes-1))
        allocate(press(0:nnodes-1))
        allocate(valence(0:nnodes-1))

        displ(2,:) = 0
        veloc(2,:) = 0
        ngllx = 0
        ngllz = 0
        valence(:) = 0
        veloc(:,:) = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            if (ngllx /= Tdomain%specel(n)%ngllx .or. &
                ngllz /= Tdomain%specel(n)%ngllz) then
                ngllx = Tdomain%specel(n)%ngllx
                ngllz = Tdomain%specel(n)%ngllz
                if (allocated(field_displ)) deallocate(field_displ)
                if (allocated(field_veloc)) deallocate(field_veloc)
                allocate(field_displ(0:ngllx-1,0:ngllz-1,2))
                allocate(field_veloc(0:ngllx-1,0:ngllz-1,2))
            endif
            call gather_elem_displ(Tdomain, n, field_displ)
            call gather_elem_veloc(Tdomain, n, field_veloc)

            do k = 0,ngllz-1
                do i = 0,ngllx-1
                    idx = irenum(Tdomain%specel(n)%Iglobnum(i,k))
                    valence(idx) = valence(idx)+1
                    displ(0:1,idx) = field_displ(i,k,:)
                    veloc(0:1,idx) = veloc(:,idx)+field_veloc(i,k,:)
                end do
            end do
        end do
        ! normalization
        do i = 0,nnodes-1
            if (valence(i)/=0) then
                veloc(:,i) = veloc(:,i)/valence(i)
            else
                write(*,*) "Elem",i," non traite"
            end if
        end do

        call h5dwrite_f(displ_id, H5T_NATIVE_DOUBLE, displ, dims, hdferr)
        call h5dwrite_f(veloc_id, H5T_NATIVE_DOUBLE, veloc, dims, hdferr)

        call h5dclose_f(displ_id, hdferr)
        call h5dclose_f(veloc_id, hdferr)
        call h5dclose_f(press_id, hdferr)
        call h5fclose_f(fid, hdferr)
        deallocate(displ,veloc,valence)

        call write_xdmf(Tdomain, rg, isort, nnodes)
    end subroutine save_field_h5

    subroutine write_master_xdmf(Tdomain)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer :: n_procs, nelem
        character (len=MAX_FILE_SIZE) :: fnamef
        integer :: rg
        n_procs = Tdomain%mpi_var%n_proc
        nelem = Tdomain%n_elem
        call semname_xdmf_master(fnamef)

        open (61,file=fnamef,status="unknown",form="formatted")
        write(61,"(a)") '<?xml version="1.0" ?>'
        write(61,"(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write(61,"(a)") '<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
        write(61,"(a)") '<Domain>'
        write(61,"(a)") '<Grid CollectionType="Spatial" GridType="Collection">'
        !!! XXX: recuperer le nom par semname_*
        do rg=0,n_procs-1
            write(61,"(a,I4.4,a)") '<xi:include href="mesh.',rg,'.xmf"/>'
        end do
        write(61,"(a)") '</Grid>'
        write(61,"(a)") '</Domain>'
        write(61,"(a)") '</Xdmf>'
        close(61)

    end subroutine write_master_xdmf

    subroutine write_xdmf(Tdomain, rg, isort, nnodes)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer, intent(in) :: rg, isort, nnodes
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer :: i, nn, ne
        real :: time
        call semname_xdmf(rg, fnamef)

        nn = nnodes
        ne = Tdomain%n_quad
        open (61,file=fnamef,status="unknown",form="formatted")
        write(61,"(a)") '<?xml version="1.0" ?>'
        write(61,"(a)") '<Grid CollectionType="Temporal" GridType="Collection">'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Mat" Format="HDF" Datatype="Int"  Dimensions="',ne, &
            '">geometry',rg,'.h5:/Material</DataItem>'

        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Mass" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',rg,'.h5:/Mass</DataItem>'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="ID" Format="HDF" Datatype="Int"  Dimensions="',ne, &
            '">geometry',rg,'.h5:/ElemID</DataItem>'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Jac" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',rg,'.h5:/Jac</DataItem>'
        time = 0
        do i=1,isort
            write(61,"(a,I4.4,a,I4.4,a)") '<Grid Name="mesh.',i,'.',rg,'">'
            write(61,"(a,F20.10,a)") '<Time Value="', time,'"/>'
            write(61,"(a,I8,a)") '<Topology Type="Quadrilateral" NumberOfElements="',ne,'">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Int" Dimensions="',ne,' 4">'
            write(61,"(a,I4.4,a)") 'geometry',rg,'.h5:/Elements'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Topology>'
            write(61,"(a)") '<Geometry Type="XY">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 2">'
            write(61,"(a,I4.4,a)") 'geometry',rg,'.h5:/Nodes'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Geometry>'
            write(61,"(a,I4.4,a)") '<Attribute Name="Displ" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/displ'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I4.4,a)") '<Attribute Name="Veloc" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/veloc'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '<Attribute Name="Domain" Center="Grid" AttributeType="Scalar" Dimensions="1">'
            write(61,"(a,I4,a)") '<DataItem Format="XML" Datatype="Int"  Dimensions="1">',rg,'</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I4.4,a)") '<Attribute Name="Mat" Center="Cell" AttributeType="Scalar" Dimensions="',ne,'">'
            !write(61,"(a,I8,a,I4.4,a)") '<DataItem Format="HDF" Datatype="Int"  Dimensions="',ne, &
            !    '">geometry',rg,'.h5:/Material</DataItem>'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="Mat"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I4.4,a)") '<Attribute Name="Mass" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="Mass"]</DataItem>'
!            write(61,"(a,I8,a,I4.4,a)") '<DataItem Format="HDF" Datatype="Int"  Dimensions="',nn, &
!                '">geometry',rg,'.h5:/Mass</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I4.4,a)") '<Attribute Name="ID" Center="Cell" AttributeType="Scalar" Dimensions="',ne,'">'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="ID"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I4.4,a)") '<Attribute Name="Jac" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="Jac"]</DataItem>'
!            write(61,"(a,I8,a,I4.4,a)") '<DataItem Format="HDF" Datatype="Int"  Dimensions="',nn, &
!                '">geometry',rg,'.h5:/Jac</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '</Grid>'
            ! XXX inexact pour l'instant
            time = time+Tdomain%TimeD%time_snapshots
        end do
        write(61,"(a)") '</Grid>'
        close(61)
    end subroutine write_xdmf

    subroutine write_constant_fields(Tdomain, fid, irenum, nnodes)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, dimension(:), intent(in), allocatable :: irenum
        integer, intent(in) :: nnodes
        !
        integer(HID_T) :: mass_id, jac_id
        integer(HSIZE_T), dimension(1) :: dims
        real, dimension(:),allocatable :: mass, jac
        integer :: hdferr
        integer :: ngllx, ngllz, idx
        integer :: i, k, n
        
        call create_dset(fid, "Mass", H5T_IEEE_F64LE, nnodes, mass_id)
        call create_dset(fid, "Jac",  H5T_IEEE_F64LE, nnodes, jac_id)

        dims(1) = Tdomain%n_glob_points
        allocate(mass(0:nnodes-1))
        allocate(jac(0:nnodes-1))
        ! mass
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            do k = 1,ngllz-2
                do i = 1,ngllx-2
                    idx = irenum(Tdomain%specel(n)%Iglobnum(i,k))
                    mass(idx) = Tdomain%specel(n)%MassMat(i,k)
                end do
            end do
        end do
        ! jac
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            do k = 0,ngllz-1
                do i = 0,ngllx-1
                    idx = irenum(Tdomain%specel(n)%Iglobnum(i,k))
                    jac(idx) = Tdomain%specel(n)%Jacob(i,k)
                end do
            end do
        end do


!        do n = 0,Tdomain%n_face-1
!            ngllx = Tdomain%sface(n)%ngll
!            do i = 1,ngllx-2
!                idx = irenum(Tdomain%sface(n)%Iglobnum_Face(i,j))
!                if (idx>=0) mass(idx) = Tdomain%sface(n)%MassMat(i,j)
!            end do
!        end do
!
!        do n = 0,Tdomain%n_vertex-1
!            idx = irenum(Tdomain%svertex(n)%Iglobnum_Vertex)
!            if (idx>=0) mass(idx) = Tdomain%svertex(n)%MassMat
!        end do

        call h5dwrite_f(mass_id, H5T_NATIVE_DOUBLE, mass, dims, hdferr)
        call h5dclose_f(mass_id, hdferr)
        call h5dwrite_f(jac_id, H5T_NATIVE_DOUBLE, jac, dims, hdferr)
        call h5dclose_f(jac_id, hdferr)
        deallocate(mass,jac)

    end subroutine write_constant_fields

end module msnapshots
