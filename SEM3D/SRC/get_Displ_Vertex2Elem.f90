subroutine get_Displ_Vertex2Elem (Tdomain,n)

use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n

integer :: nv, ngllx,nglly,ngllz
 

ngllx = Tdomain%specel(n)%ngllx; nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz

nv = Tdomain%specel(n)%near_vertices(0)
Tdomain%specel(n)%Forces(0,0,0,0) = Tdomain%sVertex(nv)%Forces(0) 
Tdomain%specel(n)%Forces(0,0,0,1) = Tdomain%sVertex(nv)%Forces(1)
Tdomain%specel(n)%Forces(0,0,0,2) = Tdomain%sVertex(nv)%Forces(2)

nv = Tdomain%specel(n)%near_vertices(1)
Tdomain%specel(n)%Forces(ngllx-1,0,0,0) = Tdomain%sVertex(nv)%Forces(0)
Tdomain%specel(n)%Forces(ngllx-1,0,0,1) = Tdomain%sVertex(nv)%Forces(1)
Tdomain%specel(n)%Forces(ngllx-1,0,0,2) = Tdomain%sVertex(nv)%Forces(2)

nv = Tdomain%specel(n)%near_vertices(2)
Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,0) = Tdomain%sVertex(nv)%Forces(0)
Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,1) = Tdomain%sVertex(nv)%Forces(1)
Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,2) = Tdomain%sVertex(nv)%Forces(2) 

nv = Tdomain%specel(n)%near_vertices(3)
Tdomain%specel(n)%Forces(0,nglly-1,0,0) = Tdomain%sVertex(nv)%Forces(0)
Tdomain%specel(n)%Forces(0,nglly-1,0,1) = Tdomain%sVertex(nv)%Forces(1)
Tdomain%specel(n)%Forces(0,nglly-1,0,2) = Tdomain%sVertex(nv)%Forces(2)

nv = Tdomain%specel(n)%near_vertices(4)
Tdomain%specel(n)%Forces(0,0,ngllz-1,0) = Tdomain%sVertex(nv)%Forces(0)
Tdomain%specel(n)%Forces(0,0,ngllz-1,1) = Tdomain%sVertex(nv)%Forces(1)
Tdomain%specel(n)%Forces(0,0,ngllz-1,2) = Tdomain%sVertex(nv)%Forces(2)

nv = Tdomain%specel(n)%near_vertices(5)
Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,0) = Tdomain%sVertex(nv)%Forces(0)
Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,1) = Tdomain%sVertex(nv)%Forces(1)
Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,2) = Tdomain%sVertex(nv)%Forces(2)

nv = Tdomain%specel(n)%near_vertices(6)
Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,0) = Tdomain%sVertex(nv)%Forces(0)
Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,1) = Tdomain%sVertex(nv)%Forces(1)
Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,2) = Tdomain%sVertex(nv)%Forces(2)

nv = Tdomain%specel(n)%near_vertices(7)
Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,0) = Tdomain%sVertex(nv)%Forces(0)
Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,1) = Tdomain%sVertex(nv)%Forces(1)
Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,2) = Tdomain%sVertex(nv)%Forces(2)


return
end subroutine get_Displ_Vertex2Elem
