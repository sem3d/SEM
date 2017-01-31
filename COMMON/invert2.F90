!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file invert2.F90
!!\brief Contient la subroutine invert2.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Inverse une matrice 2X2 en utilisant la methode de Sarrus.
!!
!! \param real, dimension (0:1,0:1), intent (INOUT) A
!! \param real, intent (OUT) det
!<


subroutine invert2 (A, det)

    !-------------------------------------------------------------------------------
    ! Subroutine to compute a 2x2 inverse matrix
    ! using the definition of inverse
    !
    ! Modified Gaetano Festa 01/06/2005
    ! -------------------------------------------------------------------------------
    use constants
    implicit none

    real(fpp), dimension (0:1,0:1), intent (INOUT) :: A
    real(fpp), intent (OUT) :: det

    real(fpp) :: aus

    ! compute determinant of A using Sarrus law

    det = A(0,0)*A(1,1)
    det = det - A(0,1)*A(1,0)

    aus = A (0,0)
    A (0,0) = A(1,1)
    A (1,0) = -A(1,0)
    A (0,1) = -A(0,1)
    A (1,1) = aus

    ! Compute the inverse

    A = A /det

    return
end subroutine invert2

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
