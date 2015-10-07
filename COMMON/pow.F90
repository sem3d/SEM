!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \fn function pow (x,vp,npm,dx,A,np)
!! \brief calcule le facteur d'attenuation PML
!!
!! \param real x
!! \param real vp
!! \param real dx
!! \param real A
!! \param integer npm
!! \param integer np
!<


real function pow (x,vp,npm,dx,A,np)

    implicit none

    real :: x,vp,dx,A
    integer :: npm,np

    ! Local variables

    real :: rnpm,pp1

    rnpm = float (npm)


    pp1 = dx
    pp1 = 1 /pp1
    pow = A * vp * pp1 * (x/rnpm)**np

    return
end function pow

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
