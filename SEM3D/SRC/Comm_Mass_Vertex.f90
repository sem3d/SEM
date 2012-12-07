subroutine Comm_Mass_Vertex (Tdomain,n,ngll,ngllPML)

! Modified by Elise Delavaud 08/02/2006


use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n
integer, intent (INOUT) :: ngll,ngllPML

integer :: i,nv
 

do i = 0,Tdomain%sComm(n)%nb_vertices-1
    nv =  Tdomain%sComm(n)%vertices(i)
    Tdomain%svertex(nv)%MassMat = Tdomain%svertex(nv)%MassMat + Tdomain%sComm(n)%Take(ngll)
    ngll = ngll + 1
    if (Tdomain%sVertex(nv)%PML) then
        Tdomain%svertex(nv)%DumpMass(0:2) = Tdomain%svertex(nv)%DumpMass(0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
         if (Tdomain%any_FPML) then
                         Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                         Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                         Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) + Tdomain%sComm(n)%TakePML(ngllPML,5)
         endif
        ngllPML = ngllPML + 1
    endif
enddo

return
end subroutine Comm_Mass_Vertex
