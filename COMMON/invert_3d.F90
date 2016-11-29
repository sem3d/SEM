!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file invert_3d.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine invert_3d (a,Det)
    use constants, only : fpp
    implicit none

    real(fpp), dimension (0:2,0:2), intent (INOUT) :: a
    real(fpp), intent (OUT) :: Det

    real(fpp), dimension(0:2,0:2) :: Inverse_A


    Det = a(0,0)* (a(1,1)*a(2,2) - a(2,1)*a(1,2) ) + a(1,0) * (a(0,2) * a(2,1) - a(0,1) * a(2,2) ) + a(2,0)  * (a(0,1)*a(1,2) - a(0,2)*a(1,1) )

    Inverse_A (0,0) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
    Inverse_A (1,0) = a(0,2) * a(2,1) - a(0,1) * a(2,2)
    Inverse_A (2,0) = a(0,1)*a(1,2) - a(0,2)*a(1,1)

    Inverse_A (0,1) = a(1,2)*a(2,0) - a(1,0)*a(2,2)
    Inverse_A (1,1) = a(0,0) * a(2,2) - a(0,2) * a(2,0)
    Inverse_A (2,1) = a(0,2)*a(1,0) - a(0,0)*a(1,2)

    Inverse_A (0,2) = a(1,0)*a(2,1) - a(1,1)*a(2,0)
    Inverse_A (1,2) = a(0,1) * a(2,0) - a(2,1) * a(0,0)
    Inverse_A (2,2) = a(1,1)*a(0,0) - a(0,1)*a(1,0)

    A = TRANSPOSE (Inverse_A)
    A = A/ Det

    return
end subroutine invert_3d

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
