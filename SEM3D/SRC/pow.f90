real function pow (x,vp,npm,dx,A,np)

    implicit none

    real :: x,vp,dx,A
    integer :: npm,np

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
