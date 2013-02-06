
module msnapshots
    use sdomain
    use hdf5
    use sem_hdf5
    use semdatafiles
    use mpi
    implicit none
contains

    subroutine create_dir_sorties(Tdomain, rg, isort)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer,intent(in) :: isort, rg
        character(Len=MAX_FILE_SIZE) :: temp
        character(len=MAX_FILE_SIZE+10) :: creer_dir
        integer :: code

        call semname_snap_result_dir(isort, temp)
        creer_dir = "mkdir -p "//temp

        if (rg==0) call system(creer_dir)
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
        integer(HID_T) :: fid, nodes_id
        integer(HSIZE_T), dimension(2) :: dims
        integer :: hdferr, code

        if (rg==0) then
            call system("mkdir -p " // path_results)
        end if
        call mpi_barrier(Tdomain%communicateur, code)
        call semname_snap_geom_file(rg, fnamef)

        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

        dims(1) = 3
        dims(2) = Tdomain%n_glob_points
        call create_dset_2d(fid, "Nodes", H5T_IEEE_F64LE, 3, Tdomain%n_glob_points, nodes_id)
        call h5dwrite_f(nodes_id, H5T_NATIVE_DOUBLE, Tdomain%GlobCoord, dims, hdferr)
        call h5dclose_f(nodes_id, hdferr)

        call write_elem_connectivity(Tdomain, fid)

        call write_constant_fields(Tdomain, fid)

        call h5fclose_f(fid, hdferr)

        if (rg==0) call write_master_xdmf(Tdomain)
    end subroutine write_snapshot_geom


    subroutine write_elem_connectivity(Tdomain, fid)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer(HID_T) :: elem_id, mat_id
        integer :: ngllx, nglly, ngllz
        integer(HSIZE_T), dimension(2) :: dims
        integer, dimension(:,:), allocatable :: data
        integer, dimension(:), allocatable :: mat
        integer :: count, i, j, k, n
        integer :: hdferr

        ! First we count the number of hexaedrons
        count = 0
        do n = 0,Tdomain%n_elem-1
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz

            count = count+(ngllx-1)*(nglly-1)*(ngllz-1)
        enddo

        Tdomain%n_hexa = count
        allocate( data(1:8,0:count-1))
        allocate( mat(0:count-1))

        call create_dset_2d(fid, "Elements", H5T_STD_I32LE, 8, count, elem_id)
        call create_dset(fid, "Material", H5T_STD_I32LE, count, mat_id)

        dims(1) = 8
        dims(2) = count
        count = 0
        do n = 0,Tdomain%n_elem-1
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            do k = 0,ngllz-2
                do j = 0,nglly-2
                    do i = 0,ngllx-2
                        data(1,count) = Tdomain%specel(n)%Iglobnum(i+0,j+0,k+0)
                        data(2,count) = Tdomain%specel(n)%Iglobnum(i+1,j+0,k+0)
                        data(3,count) = Tdomain%specel(n)%Iglobnum(i+1,j+1,k+0)
                        data(4,count) = Tdomain%specel(n)%Iglobnum(i+0,j+1,k+0)
                        data(5,count) = Tdomain%specel(n)%Iglobnum(i+0,j+0,k+1)
                        data(6,count) = Tdomain%specel(n)%Iglobnum(i+1,j+0,k+1)
                        data(7,count) = Tdomain%specel(n)%Iglobnum(i+1,j+1,k+1)
                        data(8,count) = Tdomain%specel(n)%Iglobnum(i+0,j+1,k+1)
                        mat(count) = Tdomain%specel(n)%mat_index
                        count=count+1
                    end do
                end do
            end do
        end do
        call h5dwrite_f(elem_id, H5T_NATIVE_INTEGER, data, dims, hdferr)
        call h5dclose_f(elem_id, hdferr)
        dims(1)=count
        dims(2)=1
        call h5dwrite_f(mat_id, H5T_NATIVE_INTEGER, mat, dims, hdferr)
        call h5dclose_f(mat_id, hdferr)
        deallocate(data)
    end subroutine write_elem_connectivity

    subroutine save_field_h5(Tdomain, rg, isort)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: rg, isort
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid, displ_id, veloc_id
        integer(HSIZE_T), dimension(2) :: dims
        real, dimension(:,:),allocatable :: displ, veloc
        integer :: hdferr
        integer :: ngllx, nglly, ngllz, idx
        integer :: i, j, k, n
        
        call semname_snap_result_file(rg, isort, fnamef)

        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

        dims(1) = 3
        dims(2) = Tdomain%n_glob_points
        call create_dset_2d(fid, "displ", H5T_IEEE_F64LE, 3, Tdomain%n_glob_points, displ_id)
        call create_dset_2d(fid, "veloc", H5T_IEEE_F64LE, 3, Tdomain%n_glob_points, veloc_id)

        allocate(displ(0:2,0:Tdomain%n_glob_points-1))
        allocate(veloc(0:2,0:Tdomain%n_glob_points-1))

        do n = 0,Tdomain%n_elem-1
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        idx = Tdomain%specel(n)%Iglobnum(i,j,k)
                        displ(:,idx) = Tdomain%specel(n)%Displ(i,j,k,:)
                        veloc(:,idx) = Tdomain%specel(n)%Veloc(i,j,k,:)
                    end do
                end do
            end do
        end do

        do n = 0,Tdomain%n_face-1
            ngllx = Tdomain%sface(n)%ngll1
            nglly = Tdomain%sface(n)%ngll2
            do j = 1,nglly-2
                do i = 1,ngllx-2
                    idx = Tdomain%sface(n)%Iglobnum_Face(i,j)
                    displ(:,idx) = Tdomain%sface(n)%Displ(i,j,:)
                    veloc(:,idx) = Tdomain%sface(n)%Veloc(i,j,:)
                end do
            end do
        end do

        do n = 0,Tdomain%n_edge-1
            ngllx = Tdomain%sedge(n)%ngll
            do i = 1,ngllx-2
                idx = Tdomain%sedge(n)%Iglobnum_Edge(i)
                displ(:,idx) = Tdomain%sedge(n)%Displ(i,:)
                veloc(:,idx) = Tdomain%sedge(n)%Veloc(i,:)
            end do
        end do

        do n = 0,Tdomain%n_vertex-1
            idx = Tdomain%svertex(n)%Iglobnum_Vertex
            displ(:,idx) = Tdomain%svertex(n)%Displ(:)
            veloc(:,idx) = Tdomain%svertex(n)%Veloc(:)
        end do

        call h5dwrite_f(displ_id, H5T_NATIVE_DOUBLE, displ, dims, hdferr)
        call h5dwrite_f(veloc_id, H5T_NATIVE_DOUBLE, veloc, dims, hdferr)

        call h5dclose_f(displ_id, hdferr)
        call h5dclose_f(veloc_id, hdferr)
        call h5fclose_f(fid, hdferr)
        deallocate(displ,veloc)

        call write_xdmf(Tdomain, rg, isort)
    end subroutine save_field_h5

    subroutine write_master_xdmf(Tdomain)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer :: n_procs, nelem,rg
        character (len=MAX_FILE_SIZE) :: fnamef
        integer :: n
        n_procs = Tdomain%n_proc
        nelem = Tdomain%n_elem
        call semname_xdmf_master(fnamef)

        open (61,file=fnamef,status="unknown",form="formatted")
        write(61,"(a)") '<?xml version="1.0" ?>'
        write(61,"(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write(61,"(a)") '<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
        write(61,"(a)") '<Domain>'
        write(61,"(a)") '<Grid CollectionType="Spatial" GridType="Collection">'
        !!! XXX: recuperer le nom par semname_*
        do n=0,n_procs-1
            write(61,"(a,I4.4,a)") '<xi:include href="mesh.',n,'.xmf"/>'
        end do
        write(61,"(a)") '</Grid>'
        write(61,"(a)") '</Domain>'
        write(61,"(a)") '</Xdmf>'
        close(61)

    end subroutine write_master_xdmf

    subroutine write_xdmf(Tdomain, rg, isort)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer, intent(in) :: rg, isort
        character (len=MAX_FILE_SIZE) :: fnamef
        integer :: i, nn, ne
        real :: time
        call semname_xdmf(rg, fnamef)

        nn = Tdomain%n_glob_points
        ne = Tdomain%n_hexa
        open (61,file=fnamef,status="unknown",form="formatted")
        write(61,"(a)") '<?xml version="1.0" ?>'
        write(61,"(a)") '<Grid CollectionType="Temporal" GridType="Collection">'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Mat" Format="HDF" Datatype="Int"  Dimensions="',ne, &
            '">geometry',rg,'.h5:/Material</DataItem>'

        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Mass" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',rg,'.h5:/Mass</DataItem>'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Jac" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',rg,'.h5:/Jac</DataItem>'
        time = 0
        do i=1,isort
            write(61,"(a,I4.4,a)") '<Grid Name="mesh.',rg,'">'
            write(61,"(a,F20.10,a)") '<Time Value="', time,'"/>'
            write(61,"(a,I8,a)") '<Topology Type="Hexahedron" NumberOfElements="',ne,'">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Int" Dimensions="',ne,' 8">'
            write(61,"(a,I4.4,a)") 'geometry',rg,'.h5:/Elements'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Topology>'
            write(61,"(a)") '<Geometry Type="XYZ">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a)") 'geometry',rg,'.h5:/Nodes'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Geometry>'
            write(61,"(a)") '<Attribute Name="Displ" Center="Node" AttributeType="Vector">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/displ'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '<Attribute Name="Veloc" Center="Node" AttributeType="Vector">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/veloc'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '<Attribute Name="Domain" Center="Grid" AttributeType="Scalar">'
            write(61,"(a,I4,a)") '<DataItem Format="XML" Datatype="Int"  Dimensions="1">',rg,'</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '<Attribute Name="Mat" Center="Cell" AttributeType="Scalar">'
            !write(61,"(a,I8,a,I4.4,a)") '<DataItem Format="HDF" Datatype="Int"  Dimensions="',ne, &
            !    '">geometry',rg,'.h5:/Material</DataItem>'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="Mat"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '<Attribute Name="Mass" Center="Node" AttributeType="Scalar">'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="Mass"]</DataItem>'
!            write(61,"(a,I8,a,I4.4,a)") '<DataItem Format="HDF" Datatype="Int"  Dimensions="',nn, &
!                '">geometry',rg,'.h5:/Mass</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a)") '<Attribute Name="Jac" Center="Node" AttributeType="Scalar">'
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

    subroutine write_constant_fields(Tdomain, fid)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer(HID_T) :: mass_id, jac_id
        integer(HSIZE_T), dimension(1) :: dims
        real, dimension(:),allocatable :: mass, jac
        integer :: hdferr
        integer :: ngllx, nglly, ngllz, idx
        integer :: i, j, k, n
        
        call create_dset(fid, "Mass", H5T_IEEE_F64LE, Tdomain%n_glob_points, mass_id)
        call create_dset(fid, "Jac",  H5T_IEEE_F64LE, Tdomain%n_glob_points, jac_id)

        dims(1) = Tdomain%n_glob_points
        allocate(mass(0:Tdomain%n_glob_points-1))
        allocate(jac(0:Tdomain%n_glob_points-1))
        ! mass
        do n = 0,Tdomain%n_elem-1
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        idx = Tdomain%specel(n)%Iglobnum(i,j,k)
                        mass(idx) = Tdomain%specel(n)%MassMat(i,j,k)
                    end do
                end do
            end do
        end do
        ! jac
        do n = 0,Tdomain%n_elem-1
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        idx = Tdomain%specel(n)%Iglobnum(i,j,k)
                        jac(idx) = Tdomain%specel(n)%Jacob(i,j,k)
                    end do
                end do
            end do
        end do


        do n = 0,Tdomain%n_face-1
            ngllx = Tdomain%sface(n)%ngll1
            nglly = Tdomain%sface(n)%ngll2
            do j = 1,nglly-2
                do i = 1,ngllx-2
                    idx = Tdomain%sface(n)%Iglobnum_Face(i,j)
                    mass(idx) = Tdomain%sface(n)%MassMat(i,j)
                end do
            end do
        end do

        do n = 0,Tdomain%n_edge-1
            ngllx = Tdomain%sedge(n)%ngll
            do i = 1,ngllx-2
                idx = Tdomain%sedge(n)%Iglobnum_Edge(i)
                mass(idx) = Tdomain%sedge(n)%MassMat(i)
            end do
        end do

        do n = 0,Tdomain%n_vertex-1
            idx = Tdomain%svertex(n)%Iglobnum_Vertex
            mass(idx) = Tdomain%svertex(n)%MassMat
        end do

        call h5dwrite_f(mass_id, H5T_NATIVE_DOUBLE, mass, dims, hdferr)
        call h5dclose_f(mass_id, hdferr)
        call h5dwrite_f(jac_id, H5T_NATIVE_DOUBLE, jac, dims, hdferr)
        call h5dclose_f(jac_id, hdferr)
        deallocate(mass,jac)

    end subroutine write_constant_fields

end module msnapshots
