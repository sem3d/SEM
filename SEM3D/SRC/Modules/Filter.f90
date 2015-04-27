!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
