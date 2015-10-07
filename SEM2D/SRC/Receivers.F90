!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Receivers.F90
!!\brief Contains type definition for receivers
!!
!<

module sreceivers
    type :: receiver

       logical :: located_here, on_vertex
       integer :: Nr, Nv
       real :: xRec, zRec,xi,eta
       real, dimension(:,:), pointer :: Interp_coeff
       integer, dimension(:), allocatable :: near_faces
       character(Len=100) :: name
    end type receiver

end module sreceivers

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
