subroutine Comm_Forces_Vertex (Tdomain,n,ngll,ngllPML)

! Modified by Elise Delavaud 08/02/2006


use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n
integer, intent (INOUT) :: ngll,ngllPML

integer :: i,nv
 

do i = 0,Tdomain%sComm(n)%nb_vertices-1
  nv =  Tdomain%sComm(n)%vertices(i)
  Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
  ngll = ngll + 1
  if (Tdomain%sVertex(nv)%PML) then
      Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
      Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
      Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
      ngllPML = ngllPML + 1
  endif
enddo


return
end subroutine Comm_Forces_Vertex
