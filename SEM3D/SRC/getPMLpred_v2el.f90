subroutine get_PMLprediction_v2el (Tdomain, n, bega, dt)

use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n
real, intent(IN) :: dt, bega

integer :: nv, ngllx, nglly, ngllz


ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly;  ngllz = Tdomain%specel(n)%ngllz

nv = Tdomain%specel(n)%near_vertices(0)
Tdomain%specel(n)%Forces(0,0,0,:) = Tdomain%sVertex(nv)%Veloc(:) + dt * (0.5 - bega) *  Tdomain%sVertex(nv)%Accel(:)

nv = Tdomain%specel(n)%near_vertices(1)
Tdomain%specel(n)%Forces(ngllx-1,0,0,:) = Tdomain%sVertex(nv)%Veloc(:) + dt * (0.5 - bega) *  Tdomain%sVertex(nv)%Accel(:)

nv = Tdomain%specel(n)%near_vertices(2)
Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,:) = Tdomain%sVertex(nv)%Veloc(:) + dt * (0.5 - bega) *  Tdomain%sVertex(nv)%Accel(:)

nv = Tdomain%specel(n)%near_vertices(3)
Tdomain%specel(n)%Forces(0,nglly-1,0,:) = Tdomain%sVertex(nv)%Veloc(:) + dt * (0.5 - bega) *  Tdomain%sVertex(nv)%Accel(:)

nv = Tdomain%specel(n)%near_vertices(4)
Tdomain%specel(n)%Forces(0,0,ngllz-1,:) = Tdomain%sVertex(nv)%Veloc(:) + dt * (0.5 - bega) *  Tdomain%sVertex(nv)%Accel(:)

nv = Tdomain%specel(n)%near_vertices(5)
Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,:) = Tdomain%sVertex(nv)%Veloc(:) + dt * (0.5 - bega) *  Tdomain%sVertex(nv)%Accel(:)

nv = Tdomain%specel(n)%near_vertices(6)
Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,:) = Tdomain%sVertex(nv)%Veloc(:) + dt * (0.5 - bega) *  Tdomain%sVertex(nv)%Accel(:)

nv = Tdomain%specel(n)%near_vertices(7)
Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,:) = Tdomain%sVertex(nv)%Veloc(:) + dt * (0.5 - bega) *  Tdomain%sVertex(nv)%Accel(:)


return
end subroutine get_PMLprediction_v2el
