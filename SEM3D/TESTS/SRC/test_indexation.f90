program test_indexation
    implicit none
    call test_orient()
    call test_index()

end program test_indexation


subroutine test_orient()
    use mindex
    implicit none
    !
    integer, dimension(0:3) :: face1, face2
    integer :: node, dir

    ! TEST 1
    face1 = (/ 0, 1, 2, 3 /)
    face2 = (/ 3, 0, 1, 2 /)
    call rel_orient(face1, face2, node, dir)
    write(*,*) "test_orient(1): Got node=", node, "dir=", dir
    if (node/=1 .or. dir/=1) then
        write(*,*) "test_orient(1): FAIL"
        stop 1
    else
        write(*,*) "test_orient(1): OK"
    end if
    ! TEST 2
    face1 = (/ 0, 1, 2, 3 /)
    face2 = (/ 3, 2, 1, 0 /)
    call rel_orient(face1, face2, node, dir)
    write(*,*) "test_orient(2): Got node=", node, "dir=", dir
    if (node/=3 .or. dir/=-1) then
        write(*,*) "test_orient(2): FAIL"
        stop 1
    else
        write(*,*) "test_orient(2): OK"
    end if
end subroutine test_orient


subroutine test_index()
    use mindex
    implicit none
    !
    integer, dimension(0:2) :: ngll
    integer, dimension(0:3) :: face1, face2
    integer, dimension(0:2) :: i0, di, dj
    integer :: nf
    ngll = (/3,4,5/)

    face1 = (/ 0, 1, 2, 3 /)
    face2 = (/ 3, 0, 1, 2 /)
    do nf=0,5
        call ind_elem_face(ngll, nf, face1, face2, i0, di, dj)
        write(*,*) "test_index(1): nf=", nf
        write(*,*) "test_index(1): i0=", i0
        write(*,*) "test_index(1): di=", di
        write(*,*) "test_index(1): dj=", dj
        write(*,*) "--"
    end do

end subroutine test_index
!---------------------------------------------------------------------------
!-------------------------------------------------------------------------

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
