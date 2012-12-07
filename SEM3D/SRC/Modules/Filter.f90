!>
!! \file Filter.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

module sfilter

    ! Defining filter parameter

    type filter

       real, dimension  (0:2) :: nc
       real, dimension (0:1) :: dc

    end type filter

end module sfilter
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
