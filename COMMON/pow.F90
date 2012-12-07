!>
!!\file pow.F90
!!\brief Contient la fonction pow.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \fn function pow (x,vp,npm,dx,A,np)
!! \brief
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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
