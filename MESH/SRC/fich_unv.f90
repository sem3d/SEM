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
    subroutine lec_unv_struct(namefiles,nnod,n_mat_elem,n_elem)
        implicit none
        character(len=*), dimension(0:), intent(in)  :: namefiles
        integer, intent(out)                         :: nnod, n_elem
        integer, dimension(0:), intent(out)          :: n_mat_elem
        integer   :: i,j,k, nfile, nflag, nflag_prec,ios
        integer, dimension(:), allocatable           :: flag_line
        character(len=132)  :: str
        character(len=2), parameter  :: moinz1 = "-1"


        nfile = size(namefiles)

        ! number of flags in the files : should be equal
        do i = 0, nfile-1
            j = 0
            open(10,file=namefiles(i),status="old",position="rewind",action="read",iostat=ios)
            if(ios /=0) stop "Pb to open a unv file"
            do
                read(10,*,end = 100) str
                !      i = i+1   ! line index
                if(lge(trim(adjustl(str)),moinz1) .and. lle(trim(adjustl(str)),moinz1))then
                    j = j+1
                end if
            end do
100         close(10)
            nflag = j
            if(i > 0)then
                if(nflag /= nflag_prec) stop "Compatibility problem between .unv files"
            end if
            nflag_prec = nflag
        end do
        ! OK: same number of parts in .unv files

        ! Now: number of nodes
        nflag = 0 ; i = 0
        open(10,file=namefiles(0),status="old",position="rewind",action="read",iostat=ios)
        if(ios /=0) stop "Pb to open a .unv file"
        do
            read(10,*,end = 200) str
            i = i+1   ! line index
            if(lge(trim(adjustl(str)),moinz1) .and. lle(trim(adjustl(str)),moinz1))then
                nflag = nflag+1
                if(nflag == 1)then
                    read(10,*) ; i = 0
                end if
            end if
            if(nflag == 2)then
                i = i-1 ; exit
            end if
        end do
200     close(10)
        nnod = i/2

        ! now the number of elements
        n_mat_elem(0:) = 0

        do i = 0, nfile-1
            open(10,file=namefiles(i),status="old",position="rewind",action="read",iostat=ios)
            if(ios /=0) stop "Pb to open a .unv file"
            do k = 1,2*nnod+2   ! read all node lines
                read(10,*,end = 300)
            end do
            read(10,*) k ; if(k /= -1) stop "Pb with a unv file"
            read(10,*) k ; if(k /= -1) stop "Pb with a unv file"
            read(10,*) k ; if(k /= 2412) stop "Pb with a unv file"
            j = 0
            do
                read(10,*,end = 300) str
                j = j+1
                if(lge(trim(adjustl(str)),moinz1) .and. lle(trim(adjustl(str)),moinz1))then
                    j = j-1
                end if
            end do
300         close(10)
            n_mat_elem(i) = j/2
        end do

        n_elem = SUM(n_mat_elem(0:))

    end subroutine lec_unv_struct
    !----------------------------------------------------------
    subroutine lec_unv_final(namefiles,n_mat_elem,xc,yc,zc,Ipoint)

        implicit none
        character(len=*), dimension(0:), intent(in)  :: namefiles
        integer, dimension(0:), intent(in)           :: n_mat_elem
        real, dimension(0:), intent(out)             :: xc,yc,zc
        integer, dimension(0:,0:), intent(out)       :: Ipoint
        integer                                      :: nnod,i,j,jj,k,kk,nelem,ios, nfile


        !- nodes coordinates
        open(10,file=namefiles(0),status="old",position="rewind",action="read",iostat=ios)
        if(ios /= 0) stop "Pb with .unv file opening"
        read(10,*) ; read(10,*)
        nnod = size(xc)
        do i = 0,nnod-1
            read(10,*)
            read(10,*) xc(i),yc(i),zc(i)
        end do
        close(10)

        !- nodes for each element
        nfile = size(namefiles)
        nelem = size(Ipoint,2)
        j = 0    ! elem index
        do i = 0, nfile-1
            open(10,file=namefiles(i),status="old",position="rewind",action="read",iostat=ios)
            if(ios /= 0) stop "Pb with .unv file opening"
            do k = 0,nnod
                read(10,*) ; read(10,*)
            end do
            read(10,*) ; read(10,*) ; read(10,*)

            ! now the Ipointer
            do k = 0,n_mat_elem(i)-1
                read(10,*)
                read(10,*) (Ipoint(kk,j),kk = 0,7)
                do kk = 0,7
                    Ipoint(kk,j) = Ipoint(kk,j)-1
                end do
                j = j+1
 
            end do
            close(10)
        end do

        if(j /= nelem) stop "Pb with number of elements for .unv files"

    end subroutine lec_unv_final


    subroutine lec_unv_v2(unv_files, n_points, n_elem, Material, Ipointer, xco, yco, zco, n_blocks, n_mat)
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
        n_points = points%count
    end subroutine lec_unv_v2


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

        call init_hdf5()
        ! Supporte un seul fichier
        call h5fopen_f(unv_files(0), H5F_ACC_RDONLY_F, fid, hdferr)
        !open(10,file=unv_files(i),status="old",position="rewind",action="read",iostat=ios)
        call read_dataset(fid, "Nodes", coords)
        n_points = size(coords,2)
        write(*,*) "Coords", n_points, " x ", size(coords,1)
        call read_dataset(fid, "Mat", mat)
        call read_dataset(fid, "Elements", elems)
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
