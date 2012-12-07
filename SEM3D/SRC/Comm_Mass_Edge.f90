subroutine Comm_Mass_Edge (Tdomain,n,ngll,ngllPML)

! Modified by Elise Delavaud 08/02/2006


use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n
integer, intent (INOUT) :: ngll,ngllPML

integer :: ngll1,i,j,ne
 

do i = 0,Tdomain%sComm(n)%nb_edges-1
  ne = Tdomain%sComm(n)%edges(i) 
  ngll1 = Tdomain%sEdge(ne)%ngll

  if ( Tdomain%sComm(n)%orient_edges(i) == 0 ) then 
     do j = 1,Tdomain%sEdge(ne)%ngll-2
        Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + Tdomain%sComm(n)%Take(ngll)
        ngll = ngll + 1
     enddo
     if (Tdomain%sEdge(ne)%PML) then
         do j = 1,Tdomain%sEdge(ne)%ngll-2
             Tdomain%sEdge(ne)%DumpMass(j,0:2) = Tdomain%sEdge(ne)%DumpMass(j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
             if (Tdomain%any_FPML) then
                         Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                         Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                         Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
             endif
             ngllPML = ngllPML + 1
         enddo
     endif

  else if ( Tdomain%sComm(n)%orient_edges(i) == 1 ) then 
     do j = 1,Tdomain%sEdge(ne)%ngll-2
        Tdomain%sEdge(ne)%MassMat(ngll1-1-j) = Tdomain%sEdge(ne)%MassMat(ngll1-1-j) + Tdomain%sComm(n)%Take(ngll)
        ngll = ngll + 1
     enddo
     if (Tdomain%sEdge(ne)%PML) then
         do j = 1,Tdomain%sEdge(ne)%ngll-2
             Tdomain%sEdge(ne)%DumpMass(ngll1-1-j,0:2) = Tdomain%sEdge(ne)%DumpMass(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
            if (Tdomain%any_FPML) then
                         Tdomain%sEdge(ne)%Ivx(ngll1-1-j) = Tdomain%sEdge(ne)%Ivx(ngll1-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                         Tdomain%sEdge(ne)%Ivy(ngll1-1-j) = Tdomain%sEdge(ne)%Ivy(ngll1-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                         Tdomain%sEdge(ne)%Ivz(ngll1-1-j) = Tdomain%sEdge(ne)%Ivz(ngll1-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
             endif
             ngllPML = ngllPML + 1
         enddo
     endif

  else 
     print*,'Pb with coherency number for edge' 

  endif 
   
enddo






return
end subroutine Comm_Mass_Edge
