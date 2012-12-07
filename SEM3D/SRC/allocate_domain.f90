subroutine allocate_domain (Tdomain)

use sdomain

implicit none

type(domain), intent (INOUT):: Tdomain

integer :: n, nf,ne,nv, i,j,k, ngllx,nglly,ngllz, ngll1,ngll2, ngll, ngllPML, ngllSO


do n=0,Tdomain%n_elem-1

  ngllx = Tdomain%specel(n)%ngllx
  nglly = Tdomain%specel(n)%nglly
  ngllz = Tdomain%specel(n)%ngllz

  allocate (Tdomain%specel(n)%Density (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
  allocate (Tdomain%specel(n)%Lambda (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
  allocate (Tdomain%specel(n)%Mu (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
  allocate (Tdomain%specel(n)%MassMat (0:ngllx-1, 0:nglly-1, 0:ngllz-1))

  if (Tdomain%TimeD%velocity_scheme) then
      allocate (Tdomain%specel(n)%Veloc (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
      allocate (Tdomain%specel(n)%Accel (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
      allocate (Tdomain%specel(n)%V0 (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))  
      allocate (Tdomain%specel(n)%Forces( 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
      Tdomain%specel(n)%Veloc = 0   
      Tdomain%specel(n)%Accel = 0  
      Tdomain%specel(n)%V0 = 0 
      Tdomain%specel(n)%Forces = 0
      if (Tdomain%specel(n)%PML) then
          allocate (Tdomain%specel(n)%Acoeff (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:35))
          allocate (Tdomain%specel(n)%Diagonal_Stress (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%Diagonal_Stress1 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%Diagonal_Stress2 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%Diagonal_Stress3 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%Residual_Stress (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%Residual_Stress1 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%Residual_Stress2 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%Veloc1 (1:ngllx-2, 1:nglly-2,1:ngllz-2, 0:2))
          allocate (Tdomain%specel(n)%Veloc2 (1:ngllx-2, 1:nglly-2,1:ngllz-2, 0:2))
          allocate (Tdomain%specel(n)%Veloc3 (1:ngllx-2, 1:nglly-2,1:ngllz-2, 0:2))
          allocate (Tdomain%specel(n)%Forces1 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%Forces2 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%Forces3 (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%DumpSx (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:1))
          allocate (Tdomain%specel(n)%DumpSy (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:1))
          allocate (Tdomain%specel(n)%DumpSz (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:1))
          allocate (Tdomain%specel(n)%DumpMass (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:2))
          allocate (Tdomain%specel(n)%DumpVx (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:1))
          allocate (Tdomain%specel(n)%DumpVy (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:1))
          allocate (Tdomain%specel(n)%DumpVz (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:1))
          Tdomain%specel(n)%Diagonal_Stress = 0
          Tdomain%specel(n)%Diagonal_Stress1 = 0
          Tdomain%specel(n)%Diagonal_Stress2 = 0
          Tdomain%specel(n)%Diagonal_Stress3 = 0
          Tdomain%specel(n)%Residual_Stress = 0
          Tdomain%specel(n)%Residual_Stress1 = 0
          Tdomain%specel(n)%Residual_Stress2 = 0
          Tdomain%specel(n)%Veloc1 = 0
          Tdomain%specel(n)%Veloc2 = 0
          Tdomain%specel(n)%Veloc3 = 0

	  if (Tdomain%specel(n)%FPML) then
	       allocate (Tdomain%specel(n)%Isx(0:ngllx-1,0:nglly-1,0:ngllz-1))
	       allocate (Tdomain%specel(n)%Isy(0:ngllx-1,0:nglly-1,0:ngllz-1))
	       allocate (Tdomain%specel(n)%Isz(0:ngllx-1,0:nglly-1,0:ngllz-1))
	       allocate (Tdomain%specel(n)%Ivx(0:ngllx-1,0:nglly-1,0:ngllz-1))
	       allocate (Tdomain%specel(n)%Ivy(0:ngllx-1,0:nglly-1,0:ngllz-1))
	       allocate (Tdomain%specel(n)%Ivz(0:ngllx-1,0:nglly-1,0:ngllz-1))
	       allocate (Tdomain%specel(n)%I_Diagonal_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
	       allocate (Tdomain%specel(n)%I_Diagonal_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
	       allocate (Tdomain%specel(n)%I_Diagonal_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
	       allocate (Tdomain%specel(n)%I_Residual_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
	       allocate (Tdomain%specel(n)%I_Residual_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
	       allocate (Tdomain%specel(n)%IVeloc1(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
	       allocate (Tdomain%specel(n)%IVeloc2(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
	       allocate (Tdomain%specel(n)%IVeloc3(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
	       Tdomain%specel(n)%I_Diagonal_Stress1 = 0.
	       Tdomain%specel(n)%I_Diagonal_Stress2 = 0.
	       Tdomain%specel(n)%I_Diagonal_Stress3 = 0.
	       Tdomain%specel(n)%I_Residual_Stress1 = 0.
	       Tdomain%specel(n)%I_Residual_Stress2 = 0.
	       Tdomain%specel(n)%IVeloc1 = 0.
	       Tdomain%specel(n)%IVeloc2 = 0.
	       Tdomain%specel(n)%IVeloc3 = 0.
	  endif
       else
          allocate (Tdomain%specel(n)%Acoeff (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:44))
          allocate (Tdomain%specel(n)%Displ (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))  
          Tdomain%specel(n)%Displ = 0    
       endif
   endif
enddo

do n = 0, Tdomain%n_face-1
  ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
  allocate (Tdomain%sFace(n)%MassMat (1:ngll1-2, 1:ngll2-2))
  allocate (Tdomain%sFace(n)%Veloc (1:ngll1-2, 1:ngll2-2, 0:2))
  allocate (Tdomain%sFace(n)%Forces (1:ngll1-2, 1:ngll2-2, 0:2))
  allocate (Tdomain%sFace(n)%Accel (1:ngll1-2, 1:ngll2-2, 0:2))
  allocate (Tdomain%sFace(n)%V0 (1:ngll1-2, 1:ngll2-2, 0:2))  
  Tdomain%sFace(n)%MassMat = 0
  Tdomain%sFace(n)%Veloc = 0   
  Tdomain%sFace(n)%Accel = 0  
  Tdomain%sFace(n)%V0 = 0  
  Tdomain%sFace(n)%Forces = 0  
  if (Tdomain%sFace(n)%PML ) then 
        allocate (Tdomain%sFace(n)%Forces1 (1:ngll1-2, 1:ngll2-2, 0:2)) 
        allocate (Tdomain%sFace(n)%Forces2 (1:ngll1-2, 1:ngll2-2, 0:2)) 
        allocate (Tdomain%sFace(n)%Forces3 (1:ngll1-2, 1:ngll2-2, 0:2)) 
        allocate (Tdomain%sFace(n)%DumpMass (1:ngll1-2, 1:ngll2-2, 0:2))
        allocate (Tdomain%sFace(n)%Veloc1 (1:ngll1-2, 1:ngll2-2,0:2)) 
        allocate (Tdomain%sFace(n)%Veloc2 (1:ngll1-2, 1:ngll2-2,0:2))
        allocate (Tdomain%sFace(n)%Veloc3 (1:ngll1-2, 1:ngll2-2,0:2))  
        allocate (Tdomain%sFace(n)%DumpVx (1:ngll1-2, 1:ngll2-2, 0:1)) 
        allocate (Tdomain%sFace(n)%DumpVy (1:ngll1-2, 1:ngll2-2, 0:1)) 
        allocate (Tdomain%sFace(n)%DumpVz (1:ngll1-2, 1:ngll2-2, 0:1)) 
        Tdomain%sFace(n)%DumpMass = 0.
        Tdomain%sFace(n)%Veloc1 = 0.
        Tdomain%sFace(n)%Veloc2 = 0.
        Tdomain%sFace(n)%Veloc3 = 0.
          if (Tdomain%sFace(n)%FPML ) then
	     allocate (Tdomain%sFace(n)%Ivx(1:ngll1-2,1:ngll2-1))
	     allocate (Tdomain%sFace(n)%Ivy(1:ngll1-2,1:ngll2-1))
	     allocate (Tdomain%sFace(n)%Ivz(1:ngll1-2,1:ngll2-1))
	     allocate (Tdomain%sFace(n)%IVeloc1(1:ngll1-2,1:ngll2-1,0:2))
	     allocate (Tdomain%sFace(n)%IVeloc2(1:ngll1-2,1:ngll2-1,0:2))
	     allocate (Tdomain%sFace(n)%IVeloc3(1:ngll1-2,1:ngll2-1,0:2))
	     Tdomain%sFace(n)%IVeloc1 = 0.
	     Tdomain%sFace(n)%IVeloc2 = 0.
	     Tdomain%sFace(n)%IVeloc3 = 0.
	     Tdomain%sFace(n)%Ivx = 0.
	     Tdomain%sFace(n)%Ivy = 0.
	     Tdomain%sFace(n)%Ivz = 0.
       endif
   else
       allocate (Tdomain%sFace(n)%Displ (1:ngll1-2, 1:ngll2-2, 0:2))
       Tdomain%sFace(n)%Displ = 0    
   endif
enddo

do n = 0,Tdomain%n_edge-1
    ngll = Tdomain%sEdge(n)%ngll
    allocate (Tdomain%sEdge(n)%MassMat (1:ngll-2))
    allocate (Tdomain%sEdge(n)%Veloc (1:ngll-2, 0:2))
    allocate (Tdomain%sEdge(n)%Forces (1:ngll-2, 0:2))
    allocate (Tdomain%sEdge(n)%Accel (1:ngll-2, 0:2))
    allocate (Tdomain%sEdge(n)%V0 (1:ngll-2, 0:2))
    Tdomain%sEdge(n)%MassMat = 0
    Tdomain%sEdge(n)%Veloc = 0
    Tdomain%sEdge(n)%Accel = 0
    Tdomain%sEdge(n)%V0 = 0
    Tdomain%sEdge(n)%Forces = 0
    if (Tdomain%sEdge(n)%PML) then
        allocate (Tdomain%sEdge(n)%Forces1 (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%Forces2 (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%Forces3 (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%DumpMass (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%Veloc1 (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%Veloc2 (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%Veloc3 (1:ngll-2, 0:2))
        allocate (Tdomain%sEdge(n)%DumpVx (1:ngll-2, 0:1))
        allocate (Tdomain%sEdge(n)%DumpVy (1:ngll-2, 0:1))
        allocate (Tdomain%sEdge(n)%DumpVz (1:ngll-2, 0:1))
        Tdomain%sEdge(n)%DumpMass = 0
        Tdomain%sEdge(n)%Veloc1 = 0
        Tdomain%sEdge(n)%Veloc2 = 0
        Tdomain%sEdge(n)%Veloc3 = 0
        if (Tdomain%sEdge(n)%FPML ) then
	     allocate (Tdomain%sEdge(n)%Ivx(1:ngll-2))
	     allocate (Tdomain%sEdge(n)%Ivy(1:ngll-2))
	     allocate (Tdomain%sEdge(n)%Ivz(1:ngll-2))
	     allocate (Tdomain%sEdge(n)%IVeloc1(1:ngll-2,0:2))
	     allocate (Tdomain%sEdge(n)%IVeloc2(1:ngll-2,0:2))
	     allocate (Tdomain%sEdge(n)%IVeloc3(1:ngll-2,0:2))
	     Tdomain%sEdge(n)%IVeloc1 = 0.
	     Tdomain%sEdge(n)%IVeloc2 = 0.
	     Tdomain%sEdge(n)%IVeloc3 = 0.
	     Tdomain%sEdge(n)%Ivx = 0.
	     Tdomain%sEdge(n)%Ivy = 0.
	     Tdomain%sEdge(n)%Ivz = 0.
        endif
   else
        allocate (Tdomain%sEdge(n)%Displ (1:ngll-2, 0:2))
        Tdomain%sEdge(n)%Displ = 0
    endif
enddo

do n = 0,Tdomain%n_vertex-1
    Tdomain%sVertex(n)%MassMat = 0
    Tdomain%sVertex(n)%Veloc = 0
    Tdomain%sVertex(n)%Accel = 0
    Tdomain%sVertex(n)%V0 = 0
    Tdomain%sVertex(n)%Forces = 0
   if (Tdomain%sVertex(n)%PML) then
        allocate (Tdomain%sVertex(n)%Veloc1(0:2))
        allocate (Tdomain%sVertex(n)%Veloc2(0:2))
        allocate (Tdomain%sVertex(n)%Veloc3(0:2))
        allocate (Tdomain%sVertex(n)%DumpVx(0:2))
        allocate (Tdomain%sVertex(n)%DumpVy(0:2))
        allocate (Tdomain%sVertex(n)%DumpVz(0:2))
        allocate (Tdomain%sVertex(n)%DumpMass(0:2))
        Tdomain%sVertex(n)%Veloc1 = 0
        Tdomain%sVertex(n)%Veloc2 = 0
        Tdomain%sVertex(n)%Veloc3 = 0
        Tdomain%sVertex(n)%DumpMass = 0
        if (Tdomain%sVertex(n)%FPML) then
             allocate (Tdomain%sVertex(n)%IVeloc1(0:2))
             allocate (Tdomain%sVertex(n)%IVeloc2(0:2))
             allocate (Tdomain%sVertex(n)%IVeloc3(0:2))
             allocate (Tdomain%sVertex(n)%Ivx(0:0))
             allocate (Tdomain%sVertex(n)%Ivy(0:0))
             allocate (Tdomain%sVertex(n)%Ivz(0:0))
             Tdomain%sVertex(n)%IVeloc1 = 0
             Tdomain%sVertex(n)%IVeloc2 = 0
             Tdomain%sVertex(n)%IVeloc3 = 0
             Tdomain%sVertex(n)%Ivx= 0
             Tdomain%sVertex(n)%Ivy= 0
             Tdomain%sVertex(n)%Ivz= 0
        endif
    else
        Tdomain%sVertex(n)%Displ = 0
    endif
enddo

if ( Tdomain%n_proc > 1 ) then
do n = 0,Tdomain%n_proc-1

    ngll = 0
    ngllPML = 0
    ngllSO = 0
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        ngll1 = Tdomain%sFace(nf)%ngll1; ngll2 = Tdomain%sFace(nf)%ngll2 
        ngll = ngll + (ngll1-2)*(ngll2-2)
        if (Tdomain%sFace(nf)%PML) then
          ngllPML = ngllPML + (ngll1-2)*(ngll2-2)
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        ngll = ngll + Tdomain%sEdge(ne)%ngll-2
        if (Tdomain%sEdge(ne)%PML)  then
           ngllPML = ngllPML + Tdomain%sEdge(ne)%ngll-2
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv = Tdomain%sComm(n)%vertices(i) 
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%PML) then
          ngllPML = ngllPML + 1
        endif
    enddo  
    do i = 0,Tdomain%sComm(n)%nb_edges_so-1
        ne = Tdomain%sComm(n)%edges_SO(i)
        ngllSO = ngllSO + Tdomain%sPlaneW%pEdge(ne)%ngll-2
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices_so-1
        ngllSO = ngllSO + 1
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges_neu-1
        ne = Tdomain%sComm(n)%edges_Neu(i)
        ngllSO = ngllSO + Tdomain%sNeu%nEdge(ne)%ngll-2
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices_neu-1
        ngllSO = ngllSO + 1
    enddo
    if (ngll>0) then
        allocate (Tdomain%sComm(n)%Give (0:ngll-1))
        allocate (Tdomain%sComm(n)%Take (0:ngll-1))
        allocate (Tdomain%sComm(n)%GiveForces (0:ngll-1, 0:2))
        allocate (Tdomain%sComm(n)%TakeForces (0:ngll-1, 0:2))
    endif
    if (ngllPML>0) then
        if (Tdomain%any_FPML ) then
           allocate (Tdomain%sComm(n)%GivePML (0:ngllPML-1, 0:5))
           allocate (Tdomain%sComm(n)%TakePML (0:ngllPML-1, 0:5))
        else
           allocate (Tdomain%sComm(n)%GivePML (0:ngllPML-1, 0:2))
           allocate (Tdomain%sComm(n)%TakePML (0:ngllPML-1, 0:2))
        endif
        allocate (Tdomain%sComm(n)%GiveForcesPML (0:ngllPML-1, 1:3, 0:2))
        allocate (Tdomain%sComm(n)%TakeForcesPML (0:ngllPML-1, 1:3, 0:2))
    endif
    if (ngllSO>0) then
        allocate (Tdomain%sComm(n)%GiveSO (0:ngllSO-1,0:2))
        allocate (Tdomain%sComm(n)%TakeSO (0:ngllSO-1,0:2))
    endif
    Tdomain%sComm(n)%ngll = ngll
    Tdomain%sComm(n)%ngllPML = ngllPML
    Tdomain%sComm(n)%ngllSO = ngllSO
enddo
endif


return
end subroutine allocate_domain
