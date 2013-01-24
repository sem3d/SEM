module sreceivers

    type :: receiver
       integer :: elem, proc, flag, ndt
       integer , dimension(0:2) :: gll
       real :: xRec,yRec,zRec, xi,eta,zeta
       real, dimension(:), pointer:: StoreTrace_Fl
       real, dimension(:,:), pointer:: StoreTrace
       real, dimension(:,:,:), pointer:: pol,coeff_fl
       real, dimension(:,:,:,:), pointer:: coeff
    end type receiver
end module sreceivers
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
