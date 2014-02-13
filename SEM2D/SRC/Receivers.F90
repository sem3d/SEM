!>
!!\file Receivers.F90
!!\brief Contains type definition for receivers
!!
!<

module sreceivers
    type :: receiver

       logical :: located_here
       integer :: Nr
       real :: xRec, zRec,xi,eta
       real, dimension(:,:), pointer :: Interp_coeff
       character(Len=100) :: name
    end type receiver

end module sreceivers
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
