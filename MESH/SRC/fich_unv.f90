module fich_unv
    use point_list
    implicit none

    type(PointTable) :: unv_points
contains
    subroutine unv_find_bloc(unit, bloc_id, found)
        integer, intent(in) :: unit
        integer, intent(in) :: bloc_id
        logical, intent(out) :: found
        !
        character(len=132)  :: str
        integer :: bloc
        ! Fast forward until a specific UNV bloc id is found
        rewind(unit)
        do
            read(10,*,end=100) str
            if (trim(adjustl(str))/="-1") cycle
            read(10,*) bloc
            if (bloc==bloc_id) exit
            ! read until next -1 closing the block d
            do
                read(10,*,end=100) str
                if (trim(adjustl(str))=="-1") exit
            end do
        end do
100     if (bloc==bloc_id) then
            found = .true.
        else
            found = .false.
        end if
    end subroutine unv_find_bloc

    subroutine unv_2411_info(unit, count, min_label, max_label)
        integer, intent(in) :: unit
        integer, intent(out) :: count, min_label, max_label
        integer :: label, exp_coor, disp_coor, color
        character(len=132)  :: str
        real :: x,y,z

        min_label = 2000000000
        max_label = -1
        count = 0
        do
            read(10,'(A)',end=100) str
            if (trim(adjustl(str))=="-1") exit
            read(str,*) label, exp_coor, disp_coor, color
            !write(*,*) "XXX", label, exp_coor, disp_coor, color
            read(10,'(A)',end=100) str
            read(str,*)  x,y,z
            min_label = min(label, min_label)
            max_label = max(label, max_label)
            count = count + 1
        end do
100     return
    end subroutine unv_2411_info

    subroutine unv_2412_info(unit, count)
        integer, intent(in) :: unit
        integer, intent(out) :: count
        integer :: ielem, etype, phys_prop, mat_prop, color, nnodes
        character(len=132)  :: str
        real :: x,y,z

        count = 0
        do
            read(10,'(A)',end=100) str
            if (trim(adjustl(str))=="-1") exit
            read(str,*) ielem, etype, phys_prop, mat_prop, color, nnodes
            read(10,'(A)',end=100) str
            count = count + 1
        end do
100     return
    end subroutine unv_2412_info
    
    
    subroutine unv_2477_info(unit, el_group, count)
        integer, intent(in) :: unit
        integer, intent(out) :: count
        integer, intent(inout), dimension(0:) :: el_group
        integer group_number, group_cons, group_res, group_load, group_dof, group_temp, group_cont, group_elem
        character(len=132)  :: str
        integer :: n
        
        read(10,'(A)',end=100)
        read(10,'(A)',end=100)
        read(10,'(A)',end=100)
        count = 0

        do
            read(10,'(A)',end=100)str
            if (trim(adjustl(str))=="-1") exit
            read(str,*) group_number, group_cons, group_res, group_load, group_dof, group_temp, group_cont, group_elem
            read(10,'(A)',end=100)str

            el_group(count) = group_elem

            if (mod(group_elem,2)==0) then
               do n = 0, group_elem/2-1
                  read(10,'(A)',end=100)str
               end do
            else
               do n = 0,(group_elem+1)/2-1
                  read(10,'(A)',end=100)str
               end do
            end if

            count = count+1
        end do

100     return
    end subroutine unv_2477_info

    !-----------------------------------------
    subroutine unv_2411_read(unit, points, index_map)
        integer, intent(in) :: unit
        integer, intent(inout), dimension(:) :: index_map
        type(PointTable),intent(inout) :: points
        integer :: label, exp_coor, disp_coor, color, id
        character(len=132)  :: str
        real :: x,y,z
        do
            read(10,'(A)',end=100) str
            if (trim(adjustl(str))=="-1") exit
            read(str,*) label, exp_coor, disp_coor, color
            read(10,'(A)',end=100) str
            read(str,*)  x,y,z
            call add_point(points, x, y, z, id)
            index_map(label) = id 
        end do
100     return
    end subroutine unv_2411_read

    subroutine unv_2412_read(unit, eidx, index_map, Ipointer, Epointer)
        integer, intent(in) :: unit
        integer, intent(inout) :: eidx
        integer, intent(in), dimension(:) :: index_map
        integer, dimension(0:,0:), intent(inout) :: Ipointer
        integer, dimension(0:), intent(inout) :: Epointer
        integer :: ielem, etype, mat_prop, phys_prop, color, nnodes
        integer, dimension(0:7) :: nodes
        character(len=132)  :: str
        integer :: i

        do
            read(10,'(A)',end=100) str
            if (trim(adjustl(str))=="-1") exit
            read(str,*) ielem, etype,phys_prop, mat_prop, color, nnodes
            read(10,'(A)',end=100) str

            read(str,*) nodes
            do i=0,7
                Ipointer(i, eidx) = index_map(nodes(i))
            end do
            Epointer(eidx) = ielem
            eidx = eidx + 1
            
        end do
100     return
    end subroutine unv_2412_read
    
    subroutine unv_2477_read(unit, Material,matmax, el_group, Epointer)
        integer, intent(in) :: unit
        integer, intent(in), dimension(0:) ::  el_group
        integer, intent(in), dimension(0:) ::  Epointer
        integer, intent(out), dimension(0:) :: Material
        integer :: type_code1 ,tag1 ,nodeid1 ,num_ele1,type_code2 ,tag2 ,nodeid2 ,num_ele2
        integer, intent(in) :: matmax
        character(len=132)  :: str
        integer :: n , i, k

        read(10,'(A)', end=100) str
        read(10,'(A)', end=100) str
        read(10,'(A)', end=100) str

        n = 0
        do  
          read(10,'(A)', end=100) str
            if (trim(adjustl(str))=="-1") exit
            read(10,'(A)', end=100) str         
            if (mod(el_group(n),2) == 0) then
                do i = 0, (el_group(n)/2-1)
                    read(10,'(A)', end=100) str
                    read(str,*) type_code1 ,tag1 ,nodeid1 ,num_ele1,type_code2 ,tag2 ,nodeid2 ,num_ele2
                    do k = 0, size(Material,1)-1
                        if (Epointer(k) == tag1 .OR. Epointer(k) == tag2) then
                            Material(k) = n
                        end if
                   end do
                end do
            else
                do i = 0, ((el_group(n)+1)/2-1)
                    if ((el_group(n)+1)/2-1 ==  i) then
                        read(10,'(A)', end=100) str
                        read(str,*) type_code1 ,tag1 ,nodeid1 ,num_ele1
                        do k = 0, size(Material,1)-1
                            if (Epointer(k) == tag1) then
                                Material(k) = n
                            end if
                        end do
                    else
                        read(10,'(A)', end=100) str
                        read(str,*) type_code1 ,tag1 ,nodeid1 ,num_ele1,type_code2 ,tag2 ,nodeid2 ,num_ele2
                        do k = 0, size(Material,1)-1
                            if (Epointer(k) == tag1 .OR. Epointer(k) == tag2) then
                                Material(k) = n
                            end if
                        end do
                    end if
                end do
            end if
          n=n+1
        end do
100     return
    end subroutine unv_2477_read               
    
    !-----------------------------------------
    subroutine lec_init_unv(namefiles)

        implicit none
        character(len=60), dimension(0:), intent(out) :: namefiles
        integer    :: i, nfile

        nfile = size(namefiles)
        do i = 0, nfile-1
            write(*,*) "      --> File # ",i
            read*,namefiles(i)
        end do

    end subroutine lec_init_unv
    !------------------------------------------
    subroutine lec_unv(unv_files, n_points, n_elem, Material, Ipointer, xco, yco, zco, n_blocks, n_mat)
        character(len=*), dimension(0:), intent(in)  :: unv_files
        integer, intent(out) :: n_points, n_elem
        integer, intent(out), allocatable, dimension(:) :: Material
        real, dimension(:), allocatable, intent(out) :: xco,yco,zco
        integer, dimension(:,:), intent(out), allocatable :: Ipointer
        integer, intent(out) :: n_blocks
        integer, intent(in) :: n_mat
        !
        integer :: i, nfile, loc_count, matmax, count, eidx, npts, eidx0
        integer :: min_label, max_label, loc_min_label, loc_max_label
        logical :: found
        integer :: ios
        integer, dimension(:), allocatable :: ncount, loc_index, elcount, mat_block
        type(PointTable) :: points
        real, allocatable, dimension(:,:) :: coords
        integer, dimension(:), allocatable :: el_group, Epointer
        integer :: flag_2477, mat_count

        nfile = size(unv_files)
        allocate(elcount(0:nfile-1))
        allocate(ncount(0:nfile-1))
        allocate(mat_block(0:nfile-1))
        allocate(el_group(n_mat))

        n_elem = 0
        npts = 0
        n_blocks=0
        matmax = 0
        flag_2477 = 1

        do i = 0, nfile-1
            open(10,file=unv_files(i),status="old",position="rewind",action="read",iostat=ios)
            call unv_find_bloc(10, 2411, found)
            if (.not. found) then
                stop "Datablock 2411 not found"
            end if
            call unv_2411_info(10, count, loc_min_label, loc_max_label)
            write(*,*) unv_files(i), " npoints=", count, " min=", loc_min_label, " max=", loc_max_label
            ncount(i) = loc_max_label
            npts = npts + count
            write(*,*) "Find 2412..."
            call unv_find_bloc(10, 2412, found)
            if (.not. found) then
                stop "Datablock 2412 not found"
            end if
            call unv_2412_info(10, elcount(i))
            write(*,*) "NELEM=", elcount(i)
            n_elem = n_elem + elcount(i)
            call unv_find_bloc(10, 2477, found)
            if (.not. found) then
                flag_2477 = 0
                !stop "Datablock 2477 not found"
                mat_block(i) = 1
                write(*,*)"No physical blocks"
            else
                call unv_2477_info(10, el_group,mat_block(i))
            end if
            matmax = matmax+mat_block(i)
            write(*,*) "NDOMAIN=",matmax
            n_blocks = matmax
            close(10)
        end do
        ! ---------------------------------

        call init_point_table(points, npts)
        allocate(Ipointer(0:7,0:n_elem-1))
        allocate(Material(0:n_elem-1))
        allocate(Epointer(0:n_elem-1))
        eidx = 0
        mat_count = 0

        do i = 0, nfile-1
            open(10,file=unv_files(i),status="old",position="rewind",action="read",iostat=ios)
            call unv_find_bloc(10, 2411, found)
            allocate(loc_index(0:ncount(i)))
            loc_index(:) = -1
            call unv_2411_read(10, points, loc_index)
            ncount(i) = loc_max_label
            call unv_find_bloc(10, 2412, found)
            eidx0 = eidx
            call unv_2412_read(10, eidx, loc_index, Ipointer, Epointer)
            deallocate(loc_index)
            if (flag_2477==1) then
                call unv_find_bloc(10, 2477, found)
                call unv_2477_read(10, Material, matmax, el_group, Epointer)
            else
                if (i==0) then
                    Material(0:elcount(i)-1) = i
                else           
                    mat_count = elcount(i-1)+mat_count
                    Material(mat_count:mat_count+elcount(i)-1) = i
                end if
             end if
            close(10)
        end do

        call fill_table_array(points, coords)
        allocate(xco(0:points%count-1))
        allocate(yco(0:points%count-1))
        allocate(zco(0:points%count-1))
        xco(:) = coords(0,:)
        yco(:) = coords(1,:)
        zco(:) = coords(2,:)
        deallocate(coords)
        deallocate(mat_block)
        deallocate(el_group)
        deallocate(Epointer)
        n_points = points%count
    end subroutine lec_unv

    subroutine lec_hdf5(unv_files, n_points, n_elem, Material, Ipointer, xco, yco, zco, n_blocks)
        use sem_hdf5
        character(len=*), dimension(0:), intent(in)  :: unv_files
        integer, intent(out) :: n_points, n_elem
        integer, intent(out), allocatable, dimension(:) :: Material
        real, dimension(:), allocatable, intent(out) :: xco,yco,zco
        integer, dimension(:,:), intent(out), allocatable :: Ipointer
        integer, intent(out) :: n_blocks
        !
        integer(HID_T) :: fid
        integer :: hdferr
        real, allocatable, dimension(:,:) :: coords
        integer, allocatable, dimension(:,:) :: elems
        integer, allocatable, dimension(:) :: mat
        logical :: h8exist, h27exist

        call init_hdf5()
        ! Supporte un seul fichier
        call h5fopen_f(unv_files(0), H5F_ACC_RDONLY_F, fid, hdferr)
        call read_dataset(fid, "/Nodes", coords)
        n_points = size(coords,2)
        write(*,*) "Coords", n_points, " x ", size(coords,1)
        call h5lexists_f(fid, "/Sem3D/Hexa8",  h8exist,  hdferr)
        call h5lexists_f(fid, "/Sem3D/Hexa27", h27exist, hdferr)
        if (h8exist) then
            call read_dataset(fid, "/Sem3D/Mat",   mat)
            call read_dataset(fid, "/Sem3D/Hexa8", elems)
        else if (h27exist) then
            call read_dataset(fid, "/Sem3D/Mat",    mat)
            call read_dataset(fid, "/Sem3D/Hexa27", elems)
        else
            write(*,*) "ERROR : lec_hdf5 can not find elements"
        endif
        n_elem = size(elems,2)
        write(*,*) "Elem", n_elem, " x ", size(elems,1)
        write(*,*) "Indexes:", minval(elems), " -> ", maxval(elems)
        n_blocks = maxval(mat)+1

        allocate(xco(0:n_points-1))
        allocate(yco(0:n_points-1))
        allocate(zco(0:n_points-1))
        allocate(Ipointer(0:7,0:n_elem-1))
        Ipointer(0:3,:) = elems(5:8,:)
        Ipointer(4:7,:) = elems(1:4,:)
        deallocate(elems)
        allocate(Material(0:n_elem-1))
        Material(:) = mat(:)
        deallocate(mat)
        xco(:) = coords(1,:)
        yco(:) = coords(2,:)
        zco(:) = coords(3,:)
        write(*,*) "Xrange:", minval(xco), " to ", maxval(xco)
        write(*,*) "Yrange:", minval(yco), " to ", maxval(yco)
        write(*,*) "Zrange:", minval(zco), " to ", maxval(zco)
        deallocate(coords)
    end subroutine lec_hdf5

end module fich_unv
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
