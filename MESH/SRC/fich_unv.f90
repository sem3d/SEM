module fich_unv

    implicit none


contains

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
            if(ios /=0) stop "Pb to open a .unv file"
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


end module fich_unv
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
