module sreceivers

type :: receiver


integer :: elem, proc, flag, ndt
integer , dimension(0:2) :: gll
real :: xRec,yRec,zRec, xi,eta,zeta
real, dimension(:,:), pointer:: StoreTrace
real, dimension(:,:,:), pointer:: pol
real, dimension(:,:,:,:), pointer:: coeff

end type 
end module sreceivers
