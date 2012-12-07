!>
!!\file Receivers.f90
!!\brief
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module sreceivers

    type :: receiver

       ! Modified by Gaetano Festa 31/01/2005
       ! Modified by Paul Cupillard 05/06/2005

       !  modif mariotti fevrier 2007 cea
       ! flag en plus
       integer :: elem, proc, flag
       real :: refcolat,reflong,radius,realcolat,reallong
       !  modif mariotti fevrier 2007 cea
       real :: xRec,yRec,zRec, xi,eta,zeta

       real, dimension(0:2) :: gll
       real, dimension(0:2,0:2) :: Passage
       !! initialement allocatable, ne passe pas avec gfortran -> pointer
       !!real, dimension(:), allocatable :: cosgamma,singamma
       !!real, dimension(:,:), allocatable :: StoreTrace
       !!real, dimension(:,:,:), allocatable :: pol
       !!real, dimension(:,:,:,:), allocatable :: coeff
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
