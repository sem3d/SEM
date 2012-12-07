subroutine get_Displ_Edge2Elem (Tdomain,n)

use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n

integer :: i, ne, ngllx,nglly,ngllz,ngll
 

ngllx = Tdomain%specel(n)%ngllx; nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz

ne = Tdomain%specel(n)%near_edges(0)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(0) == 0 ) then
  Tdomain%specel(n)%Forces(1:ngllx-2,0,0,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0) 
  Tdomain%specel(n)%Forces(1:ngllx-2,0,0,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(1:ngllx-2,0,0,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
     Tdomain%specel(n)%Forces(ngllx-1-i,0,0,0) = Tdomain%sEdge(ne)%Forces(i,0) 
     Tdomain%specel(n)%Forces(ngllx-1-i,0,0,1) = Tdomain%sEdge(ne)%Forces(i,1) 
     Tdomain%specel(n)%Forces(ngllx-1-i,0,0,2) = Tdomain%sEdge(ne)%Forces(i,2) 
  enddo
endif

ne = Tdomain%specel(n)%near_edges(1)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(1) == 0 ) then
  Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,0,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,0,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,0,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,0,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,0,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,0,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

ne = Tdomain%specel(n)%near_edges(2)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(2) == 0 ) then
  Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,0,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,0,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,0,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2) 
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,0,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,0,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,0,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

ne = Tdomain%specel(n)%near_edges(3)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(3) == 0 ) then
  Tdomain%specel(n)%Forces(0,1:nglly-2,0,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(0,1:nglly-2,0,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(0,1:nglly-2,0,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(0,nglly-1-i,0,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(0,nglly-1-i,0,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(0,nglly-1-i,0,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

ne = Tdomain%specel(n)%near_edges(4)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(4) == 0 ) then
  Tdomain%specel(n)%Forces(ngllx-1,0,1:ngllz-2,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(ngllx-1,0,1:ngllz-2,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(ngllx-1,0,1:ngllz-2,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1-i,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1-i,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1-i,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

ne = Tdomain%specel(n)%near_edges(5)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(5) == 0 ) then
  Tdomain%specel(n)%Forces(1:ngllx-2,0,ngllz-1,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(1:ngllx-2,0,ngllz-1,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(1:ngllx-2,0,ngllz-1,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

ne = Tdomain%specel(n)%near_edges(6)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(6) == 0 ) then
  Tdomain%specel(n)%Forces(0,0,1:ngllz-2,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(0,0,1:ngllz-2,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(0,0,1:ngllz-2,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(0,0,ngllz-1-i,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(0,0,ngllz-1-i,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(0,0,ngllz-1-i,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

ne = Tdomain%specel(n)%near_edges(7)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(7) == 0 ) then
  Tdomain%specel(n)%Forces(ngllx-1,nglly-1,1:ngllz-2,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(ngllx-1,nglly-1,1:ngllz-2,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(ngllx-1,nglly-1,1:ngllz-2,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1-i,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1-i,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1-i,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

ne = Tdomain%specel(n)%near_edges(8)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(8) == 0 ) then
  Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,ngllz-1,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,ngllz-1,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,ngllz-1,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1,0) = Tdomain%sEdge(ne)%Forces(i,0) 
    Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1,1) = Tdomain%sEdge(ne)%Forces(i,1) 
    Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1,2) = Tdomain%sEdge(ne)%Forces(i,2) 
  enddo
endif

ne = Tdomain%specel(n)%near_edges(9)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(9) == 0 ) then
  Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,ngllz-1,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,ngllz-1,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,ngllz-1,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

ne = Tdomain%specel(n)%near_edges(10)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(10) == 0 ) then
  Tdomain%specel(n)%Forces(0,nglly-1,1:ngllz-2,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(0,nglly-1,1:ngllz-2,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(0,nglly-1,1:ngllz-2,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1-i,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1-i,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1-i,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

ne = Tdomain%specel(n)%near_edges(11)
ngll = Tdomain%sEdge(ne)%ngll
if ( Tdomain%specel(n)%Orient_Edges(11) == 0 ) then
  Tdomain%specel(n)%Forces(0,1:nglly-2,ngllz-1,0) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0)
  Tdomain%specel(n)%Forces(0,1:nglly-2,ngllz-1,1) = Tdomain%sEdge(ne)%Forces(1:ngll-2,1)
  Tdomain%specel(n)%Forces(0,1:nglly-2,ngllz-1,2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,2)
else
  do i=1,ngll-2
    Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1,0) = Tdomain%sEdge(ne)%Forces(i,0)
    Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1,1) = Tdomain%sEdge(ne)%Forces(i,1)
    Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1,2) = Tdomain%sEdge(ne)%Forces(i,2)
  enddo
endif

return
end subroutine get_Displ_Edge2Elem
