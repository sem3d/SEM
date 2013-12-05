subroutine lec1_cubit(indunit,ln,n_mat)

    implicit none

    character(len=70)    :: ligne
    integer, intent(in)  :: indunit
    integer, intent(out) :: ln, n_mat
    integer              :: i


    i = 0

    ! counting of global nodes
    do
        i = i+1
        read(indunit,"(a)",end=100) ligne
        if(ligne(1:2) == "*N")then
            i = 0
        end if
        if(ligne(1:2) == "*E") then
            ln = i-1
            exit
        endif
    end do


    ! Number of materials
100 n_mat = 0
    rewind(indunit)
    do
        read(indunit,"(a)",end=101) ligne
        if(ligne(1:2) == "*E") n_mat = n_mat+1
    end do

101 return

end subroutine lec1_cubit
!--------------------------------------------
subroutine lec2_cubit(indunit,nm,n_elem,elem,ind_mat)

    implicit none

    integer, intent(in) :: indunit,nm
    character(len=70)   :: ligne
    character           :: c
    integer   :: l,i,n_mat,nn
    integer, intent(out) :: n_elem
    integer, intent(inout), dimension(nm) :: elem, ind_mat

    n_elem = 0

    ! Number of elements per block
    elem(:) = 0
    nn = 0
    rewind(indunit)
    do
        read(indunit,"(a)",end=102) ligne
        if(ligne(1:2) == "*E") nn = nn+1  ! material increment
        if(nn < 1) cycle
        elem(nn) = elem(nn)+1
    end do

102 elem(:) = elem(:)-1

    ! total number of elements
    n_elem = SUM(elem)


    ! Indices of materials
    ind_mat(:) = 0
    rewind(indunit)
    nn = 0
    do
        read(indunit,"(a)",end=121) ligne
        if(ligne(1:2) == "*E")then
            nn = nn+1
            l = len_trim(ligne)
            c = ligne(l-1:l-1)
            i = iachar(c)-iachar('0')
            if(i>=1 .and. i<=9) ind_mat(nn) = 10*i
            c = ligne(l:l)
            i = iachar(c)-iachar('0')
            if(i>=1 .and. i<=9) ind_mat(nn) = ind_mat(nn)+i
        end if
    end do


121 rewind(indunit)

    return

end subroutine lec2_cubit
!---------------------------------------------------------
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
