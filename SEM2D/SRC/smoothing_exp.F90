!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file smoothing_exp.F90
!!\brief Contient la subroutine smoothing_exp.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Evalue le coefficient de lissage smooth.
!!
!! \param integer, intent (IN) n1
!! \param integer, intent (IN) n2
!! \param integer, intent (IN) n_point1
!! \param real, intent (IN) sigma
!! \param real, dimension (0:n1,0:n2), intent (INOUT) smooth
!<


subroutine smoothing_exp (smooth, n1, n2, n_point1,sigma)
    ! Modified by Gaetano Festa 19/05/2005

    implicit none
    integer, intent (IN) :: n1, n2, n_point1
    real, intent (IN) :: sigma
    real, dimension (0:n1,0:n2), intent (INOUT) :: smooth

    integer :: i,j
    real :: weight, value_m, distance0, value0, dis
    real, dimension (0:n1,0:1) :: store

    Store(0:n1,0:1)  = Smooth(0:n1,0:1)
    do i = 0, n_point1-2
        distance0 = store (i,0); value0 = Store(i,1)
        value_m = value0; weight = 1.
        do j = 1, i
            dis = abs (Store(i+j,0) -distance0)
            dis = exp(-4*dis**2/sigma**2)
            value_m = value_m + Store(i+j,1) *dis
            weight = weight + dis
            dis = abs (Store(i-j,0) -distance0)
            dis = exp(-4*dis**2/sigma**2)
            value_m = value_m + Store(i-j,1) *dis
            weight = weight + dis
        enddo
        Smooth (i,1) = value_m /weight
    enddo
    do i = n1-n_point1+2, n1
        distance0 = store (i,0); value0 = Store(i,1)
        value_m = value0; weight = 1.
        do j = 1, n1-i
            dis = abs (Store(i+j,0) -distance0)
            dis = exp(-4*dis**2/sigma**2)
            value_m = value_m + Store(i+j,1) *dis
            weight = weight + dis
            dis = abs (Store(i-j,0) -distance0)
            dis = exp(-4*dis**2/sigma**2)
            value_m = value_m + Store(i-j,1) *dis
            weight = weight + dis
        enddo
        Smooth (i,1) = value_m /weight
    enddo
    if (n_point1> 0) then
        do i = n_point1-1, n1-n_point1+1
            distance0 = store (i,0); value0 = Store(i,1)
            value_m = value0; weight = 1.
            do j = 1, n_point1-1
                dis = abs (Store(i+j,0) -distance0)
                dis = exp(-4*dis**2/sigma**2)
                value_m = value_m + Store(i+j,1) *dis
                weight = weight + dis
                dis = abs (Store(i-j,0) -distance0)
                dis = exp(-4*dis**2/sigma**2)
                value_m = value_m + Store(i-j,1) *dis
                weight = weight + dis
            enddo
            Smooth (i,1) = value_m /weight
        enddo
    endif
    return

end subroutine smoothing_exp

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
