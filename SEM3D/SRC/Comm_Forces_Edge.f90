subroutine Comm_Forces_Edge (Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)

! Modified by Elise Delavaud 08/02/2006


use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n
integer, intent (INOUT) :: ngll,ngllPML,ngll_F,ngllPML_F


integer :: ngll1,i,j,ne
 

do i = 0,Tdomain%sComm(n)%nb_edges-1
  ne = Tdomain%sComm(n)%edges(i)
  ngll1 = Tdomain%sEdge(ne)%ngll

  if ( Tdomain%sComm(n)%orient_edges(i) == 0 ) then 
   if(Tdomain%sEdge(ne)%solid)then  !solid part
    do j = 1,Tdomain%sEdge(ne)%ngll-2
      Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
      ngll = ngll + 1
    enddo
    if (Tdomain%sEdge(ne)%PML) then
      do j = 1,Tdomain%sEdge(ne)%ngll-2
          Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
          Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
          Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
          ngllPML = ngllPML + 1
      enddo
    endif
   else   ! fluid part
    do j = 1,Tdomain%sEdge(ne)%ngll-2
      Tdomain%sEdge(ne)%ForcesFl(j) = Tdomain%sEdge(ne)%ForcesFl(j) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
      ngll_F = ngll_F + 1
    enddo
    if (Tdomain%sEdge(ne)%PML) then
      do j = 1,Tdomain%sEdge(ne)%ngll-2
          Tdomain%sEdge(ne)%ForcesFl1(j) = Tdomain%sEdge(ne)%ForcesFl1(j) + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
          Tdomain%sEdge(ne)%ForcesFl2(j) = Tdomain%sEdge(ne)%ForcesFl2(j) + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
          Tdomain%sEdge(ne)%ForcesFl3(j) = Tdomain%sEdge(ne)%ForcesFl3(j) + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
          ngllPML_F = ngllPML_F + 1
      enddo
    endif

   end if

  else  
   if(Tdomain%sEdge(ne)%solid)then  !solid part
    do j = 1,Tdomain%sEdge(ne)%ngll-2
      Tdomain%sEdge(ne)%Forces(ngll1-1-j,0:2) = Tdomain%sEdge(ne)%Forces(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
      ngll = ngll + 1
    enddo
    if (Tdomain%sEdge(ne)%PML) then
      do j = 1,Tdomain%sEdge(ne)%ngll-2
          Tdomain%sEdge(ne)%Forces1(ngll1-1-j,0:2) = Tdomain%sEdge(ne)%Forces1(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
          Tdomain%sEdge(ne)%Forces2(ngll1-1-j,0:2) = Tdomain%sEdge(ne)%Forces2(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
          Tdomain%sEdge(ne)%Forces3(ngll1-1-j,0:2) = Tdomain%sEdge(ne)%Forces3(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
          ngllPML = ngllPML + 1
      enddo
    endif
   else   ! fluid part
    do j = 1,Tdomain%sEdge(ne)%ngll-2
      Tdomain%sEdge(ne)%ForcesFl(ngll1-1-j) = Tdomain%sEdge(ne)%ForcesFl(ngll1-1-j) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
      ngll_F = ngll_F + 1
    enddo
    if (Tdomain%sEdge(ne)%PML) then
      do j = 1,Tdomain%sEdge(ne)%ngll-2
          Tdomain%sEdge(ne)%ForcesFl1(ngll1-1-j) = Tdomain%sEdge(ne)%ForcesFl1(ngll1-1-j) + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
          Tdomain%sEdge(ne)%ForcesFl2(ngll1-1-j) = Tdomain%sEdge(ne)%ForcesFl2(ngll1-1-j) + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
          Tdomain%sEdge(ne)%ForcesFl3(ngll1-1-j) = Tdomain%sEdge(ne)%ForcesFl3(ngll1-1-j) + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
          ngllPML_F = ngllPML_F + 1
      enddo
    endif

   end if

  endif

enddo

return
end subroutine Comm_Forces_Edge
