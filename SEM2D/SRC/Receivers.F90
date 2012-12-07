!>
!!\file Receivers.F90
!!\brief
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module sreceivers

    ! Modified by gaetano Festa 02/06/05
    ! Minor revisions Gaetano Festa (MPI) 13/10/05
    type :: receiver

       logical :: located_here
       integer :: Nr
       real :: xRec, zRec,xi,eta
       real, dimension(:,:), pointer :: Interp_coeff

    end type receiver

end module sreceivers
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
