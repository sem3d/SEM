subroutine get_PMLprediction_f2el (Tdomain, n, bega, dt)

use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n
real, intent(IN) :: dt, bega

integer :: nf, ngllx, nglly, ngllz, ngll1, ngll2, i, j


ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly;  ngllz = Tdomain%specel(n)%ngllz

nf = Tdomain%specel(n)%near_faces(0)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(0) == 0 ) then
   do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(i,j,0,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(i,j,0,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(i,j,0,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 1 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-i,j,0,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1-i,j,0,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-i,j,0,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 2 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(i,nglly-1-j,0,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(i,nglly-1-j,0,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(i,nglly-1-j,0,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 3 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,0,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,0,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,0,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 4 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(j,i,0,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(j,i,0,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(j,i,0,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 5 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-j,i,0,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1-j,i,0,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-j,i,0,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 6 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(j,nglly-1-i,0,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(j,nglly-1-i,0,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(j,nglly-1-i,0,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 7 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,0,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,0,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,0,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
endif



nf = Tdomain%specel(n)%near_faces(1)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(1) == 0 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(i,0,j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(i,0,j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(i,0,j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 1 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-i,0,j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(ngllx-1-i,0,j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-i,0,j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 2 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(i,0,ngllz-1-j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(i,0,ngllz-1-j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(i,0,ngllz-1-j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 3 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1-j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1-j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1-j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 4 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(j,0,i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(j,0,i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(j,0,i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2) 
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 5 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-j,0,i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(ngllx-1-j,0,i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-j,0,i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 6 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(j,0,ngllz-1-i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(j,0,ngllz-1-i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(j,0,ngllz-1-i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 7 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-j,0,ngllz-1-i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(ngllx-1-j,0,ngllz-1-i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-j,0,ngllz-1-i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
endif

nf = Tdomain%specel(n)%near_faces(2)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(2) == 0 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1,i,j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1,i,j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1,i,j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 1 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 2 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1,i,ngllz-1-j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1,i,ngllz-1-j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1,i,ngllz-1-j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 3 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1-j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1-j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1-j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 4 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1,j,i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1,j,i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1,j,i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 5 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1,1:nglly-1-j,i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1,1:nglly-1-j,i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1,1:nglly-1-j,i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 6 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1,j,ngllz-1-i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1,j,ngllz-1-i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1,j,ngllz-1-i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 7 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,ngllz-1-i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,ngllz-1-i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,ngllz-1-i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
endif


nf = Tdomain%specel(n)%near_faces(3)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(3) == 0 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(i,nglly-1,j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(i,nglly-1,j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(i,nglly-1,j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 1 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 2 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(i,nglly-1,ngllz-1-j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(i,nglly-1,ngllz-1-j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(i,nglly-1,ngllz-1-j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 3 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1-j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1-j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1-j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 4 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(j,nglly-1,i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(j,nglly-1,i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(j,nglly-1,i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 5 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 6 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(j,nglly-1,ngllz-1-i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(j,nglly-1,ngllz-1-i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(j,nglly-1,ngllz-1-i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 7 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,ngllz-1-i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,ngllz-1-i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,ngllz-1-i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
endif

nf = Tdomain%specel(n)%near_faces(4)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(4) == 0 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(0,i,j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(0,i,j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(0,i,j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 1 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(0,nglly-1-i,j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(0,nglly-1-i,j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(0,nglly-1-i,j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 2 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(0,i,ngllz-1-j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(0,i,ngllz-1-j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(0,i,ngllz-1-j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 3 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1-j,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1-j,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1-j,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 4 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(0,j,i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(0,j,i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(0,j,i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 5 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(0,nglly-1-j,i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(0,nglly-1-j,i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(0,nglly-1-j,i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 6 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(0,j,ngllz-1-i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(0,j,ngllz-1-i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(0,j,ngllz-1-i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 7 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(0,nglly-1-j,ngllz-1-i,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0) 
      Tdomain%specel(n)%Forces(0,nglly-1-j,ngllz-1-i,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1) 
      Tdomain%specel(n)%Forces(0,nglly-1-j,ngllz-1-i,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
endif


nf = Tdomain%specel(n)%near_faces(5)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(5) == 0 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(i,j,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(i,j,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(i,j,ngllz-1,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 1 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-i,j,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1-i,j,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-i,j,ngllz-1,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 2 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(i,nglly-1-j,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(i,nglly-1-j,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(i,nglly-1-j,ngllz-1,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 3 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1-j,ngllz-1,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 4 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(j,i,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(j,i,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(j,i,ngllz-1,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 5 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-j,i,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1-j,i,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-j,i,ngllz-1,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 6 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(j,nglly-1-i,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(j,nglly-1-i,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(j,nglly-1-i,ngllz-1,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 7 ) then
  do i=1,ngll1-2
    do j=1,ngll2-2
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(i,j,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,0)
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(i,j,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,1)
      Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-i,ngllz-1,2) =  Tdomain%sFace(nf)%Veloc(i,j,2) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(i,j,2)
    enddo
  enddo
endif

return
end subroutine get_PMLprediction_f2el
