!>
!!\file Receivers.f90
!!\brief
!!
!<

module sreceivers

    type :: receiver
       integer :: elem, proc, flag, ndt
       real :: refcolat,reflong,radius,realcolat,reallong
       real :: xRec,yRec,zRec, xi,eta,zeta

       real, dimension(0:2) :: gll
       real, dimension(0:2,0:2) :: Passage
       real, dimension(:), pointer :: cosgamma,singamma
       real, dimension(:,:), pointer :: StoreTrace
       real, dimension(:,:,:), pointer :: pol
       real, dimension(:,:,:,:), pointer :: coeff

    end type receiver
end module sreceivers
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
