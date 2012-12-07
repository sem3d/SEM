subroutine deallocate_domain (Tdomain)

use sdomain

implicit none

type(domain), intent (INOUT):: Tdomain

integer :: n


deallocate (Tdomain%GlobCoord)
deallocate (Tdomain%Coord_Nodes)

do n = 0,Tdomain%n_elem-1
  deallocate (Tdomain%specel(n)%Density)
  deallocate (Tdomain%specel(n)%MassMat)
  deallocate (Tdomain%specel(n)%IglobNum)
  deallocate (Tdomain%specel(n)%Control_Nodes)
  deallocate (Tdomain%specel(n)%Jacob)
  if (Tdomain%TimeD%velocity_scheme) then
      deallocate (Tdomain%specel(n)%Veloc )
      deallocate (Tdomain%specel(n)%Accel)
      deallocate (Tdomain%specel(n)%V0 )  
      deallocate (Tdomain%specel(n)%Forces)
      if (Tdomain%specel(n)%PML) then
          deallocate (Tdomain%specel(n)%Acoeff)
          deallocate (Tdomain%specel(n)%Diagonal_Stress)
          deallocate (Tdomain%specel(n)%Diagonal_Stress1)
          deallocate (Tdomain%specel(n)%Diagonal_Stress2)
          deallocate (Tdomain%specel(n)%Diagonal_Stress3)
          deallocate (Tdomain%specel(n)%Residual_Stress)
          deallocate (Tdomain%specel(n)%Residual_Stress1)
          deallocate (Tdomain%specel(n)%Residual_Stress2)
          deallocate (Tdomain%specel(n)%Veloc1)
          deallocate (Tdomain%specel(n)%Veloc2)
          deallocate (Tdomain%specel(n)%Veloc3)
          deallocate (Tdomain%specel(n)%Forces1)
          deallocate (Tdomain%specel(n)%Forces2)
          deallocate (Tdomain%specel(n)%Forces3)
          deallocate (Tdomain%specel(n)%DumpSx)
          deallocate (Tdomain%specel(n)%DumpSy)
          deallocate (Tdomain%specel(n)%DumpSz)
          deallocate (Tdomain%specel(n)%DumpVx)
          deallocate (Tdomain%specel(n)%DumpVy)
          deallocate (Tdomain%specel(n)%DumpVz)
       else
          deallocate (Tdomain%specel(n)%Acoeff)
          deallocate (Tdomain%specel(n)%Displ )  
       endif
   endif
enddo

do n = 0, Tdomain%n_face-1
    deallocate (Tdomain%sFace(n)%MassMat)
    deallocate (Tdomain%sFace(n)%Veloc)
    deallocate (Tdomain%sFace(n)%Forces)
    deallocate (Tdomain%sFace(n)%Accel)
    deallocate (Tdomain%sFace(n)%V0)  
    if (Tdomain%sFace(n)%PML) then 
        deallocate (Tdomain%sFace(n)%Forces1) 
        deallocate (Tdomain%sFace(n)%Forces2) 
        deallocate (Tdomain%sFace(n)%Forces3) 
        deallocate (Tdomain%sFace(n)%Veloc1) 
        deallocate (Tdomain%sFace(n)%Veloc2)
        deallocate (Tdomain%sFace(n)%Veloc3)  
        deallocate (Tdomain%sFace(n)%DumpVx) 
        deallocate (Tdomain%sFace(n)%DumpVy) 
        deallocate (Tdomain%sFace(n)%DumpVz) 
    else
        deallocate (Tdomain%sFace(n)%Displ)
    endif
enddo

do n = 0,Tdomain%n_edge-1
    deallocate (Tdomain%sEdge(n)%MassMat)
    deallocate (Tdomain%sEdge(n)%Veloc)
    deallocate (Tdomain%sEdge(n)%Forces)
    deallocate (Tdomain%sEdge(n)%Accel)
    deallocate (Tdomain%sEdge(n)%V0)
    if (Tdomain%sEdge(n)%PML) then
        deallocate (Tdomain%sEdge(n)%Forces1)
        deallocate (Tdomain%sEdge(n)%Forces2)
        deallocate (Tdomain%sEdge(n)%Forces3)
        deallocate (Tdomain%sEdge(n)%Veloc1)
        deallocate (Tdomain%sEdge(n)%Veloc2)
        deallocate (Tdomain%sEdge(n)%Veloc3)
        deallocate (Tdomain%sEdge(n)%DumpVx)
        deallocate (Tdomain%sEdge(n)%DumpVy)
        deallocate (Tdomain%sEdge(n)%DumpVz)
    else
        deallocate (Tdomain%sEdge(n)%Displ)
   endif
enddo

do n = 0,Tdomain%n_proc-1
    if (Tdomain%sComm(n)%ngll>0) then
        deallocate (Tdomain%sComm(n)%GiveForces)
        deallocate (Tdomain%sComm(n)%TakeForces)
    endif
    if (Tdomain%sComm(n)%ngllPML>0) then
        deallocate (Tdomain%sComm(n)%GiveForcesPML)
        deallocate (Tdomain%sComm(n)%TakeForcesPML)
    endif
enddo

do n = 0, Tdomain%n_mat-1 
    if (associated(Tdomain%sSubdomain(n)%GLLcz, Tdomain%sSubdomain(n)%GLLcx) .or. &
        associated(Tdomain%sSubdomain(n)%GLLcz, Tdomain%sSubdomain(n)%GLLcy )) then
        nullify (Tdomain%sSubdomain(n)%GLLcz)
        nullify(Tdomain%sSubdomain(n)%GLLwz)
        nullify(Tdomain%sSubdomain(n)%hprimez)
        nullify(Tdomain%sSubdomain(n)%hTprimez)
    else
        deallocate (Tdomain%sSubdomain(n)%GLLcz)
        deallocate (Tdomain%sSubdomain(n)%GLLwz)
        deallocate (Tdomain%sSubdomain(n)%hprimez)
        deallocate (Tdomain%sSubdomain(n)%hTprimez)
    endif
    if (associated(Tdomain%sSubdomain(n)%GLLcy, Tdomain%sSubdomain(n)%GLLcx) ) then
        nullify (Tdomain%sSubdomain(n)%GLLcy)
        nullify(Tdomain%sSubdomain(n)%GLLwy)
        nullify(Tdomain%sSubdomain(n)%hprimey)
        nullify(Tdomain%sSubdomain(n)%hTprimey)
    else
        deallocate (Tdomain%sSubdomain(n)%GLLcy)
        deallocate (Tdomain%sSubdomain(n)%GLLwy)
        deallocate (Tdomain%sSubdomain(n)%hprimey)
        deallocate (Tdomain%sSubdomain(n)%hTprimey)
    endif
    deallocate (Tdomain%sSubdomain(n)%GLLcx)
    deallocate (Tdomain%sSubdomain(n)%GLLwx)
    deallocate (Tdomain%sSubdomain(n)%hprimex)
    deallocate (Tdomain%sSubdomain(n)%hTprimex)
enddo

return
end subroutine deallocate_domain
