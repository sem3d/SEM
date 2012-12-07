subroutine allocate_domain(Tdomain,rank)

    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)  :: rank
    integer :: n,nf,ne,nv,i,j,k,ngllx,nglly,ngllz,ngll1,ngll2,   &
                ngll,ngllPML,ngllSO,ngllNeu,ngllSF,ngll_F,ngllPML_F


do n = 0,Tdomain%n_elem-1

 ngllx = Tdomain%specel(n)%ngllx
 nglly = Tdomain%specel(n)%nglly
 ngllz = Tdomain%specel(n)%ngllz

 allocate(Tdomain%specel(n)%Density(0:ngllx-1,0:nglly-1,0:ngllz-1))
 allocate(Tdomain%specel(n)%Lambda(0:ngllx-1,0:nglly-1,0:ngllz-1))
 allocate(Tdomain%specel(n)%Mu(0:ngllx-1,0:nglly-1, 0:ngllz-1))
 allocate(Tdomain%specel(n)%MassMat(0:ngllx-1,0:nglly-1,0:ngllz-1))


 if(Tdomain%specel(n)%solid)then  ! SOLID PART
  allocate(Tdomain%specel(n)%Density(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(Tdomain%specel(n)%Lambda(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(Tdomain%specel(n)%Mu(0:ngllx-1,0:nglly-1, 0:ngllz-1))
  allocate(Tdomain%specel(n)%MassMat(0:ngllx-1,0:nglly-1,0:ngllz-1))

  if(Tdomain%TimeD%velocity_scheme)then
      allocate(Tdomain%specel(n)%Veloc(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
      allocate(Tdomain%specel(n)%Accel(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
      allocate(Tdomain%specel(n)%V0(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))  
      allocate(Tdomain%specel(n)%Forces(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
      Tdomain%specel(n)%Veloc = 0d0   
      Tdomain%specel(n)%Accel = 0d0
      Tdomain%specel(n)%V0 = 0d0
      Tdomain%specel(n)%Forces = 0d0
      if(Tdomain%specel(n)%PML)then
          allocate(Tdomain%specel(n)%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:35))
          allocate(Tdomain%specel(n)%Diagonal_Stress(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%Diagonal_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%Diagonal_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%Diagonal_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%Residual_Stress(0:ngllx-1, 0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%Residual_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%Residual_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%Veloc1(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
          allocate(Tdomain%specel(n)%Veloc2(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
          allocate(Tdomain%specel(n)%Veloc3(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
          allocate(Tdomain%specel(n)%Forces1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%Forces2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%Forces3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%DumpSx(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
          allocate(Tdomain%specel(n)%DumpSy(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
          allocate(Tdomain%specel(n)%DumpSz(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
          allocate(Tdomain%specel(n)%DumpMass(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%DumpVx(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1))
          allocate(Tdomain%specel(n)%DumpVy(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1))
          allocate(Tdomain%specel(n)%DumpVz(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1))
          Tdomain%specel(n)%Diagonal_Stress = 0d0
          Tdomain%specel(n)%Diagonal_Stress1 = 0d0
          Tdomain%specel(n)%Diagonal_Stress2 = 0d0
          Tdomain%specel(n)%Diagonal_Stress3 = 0d0
          Tdomain%specel(n)%Residual_Stress = 0d0
          Tdomain%specel(n)%Residual_Stress1 = 0d0
          Tdomain%specel(n)%Residual_Stress2 = 0d0
          Tdomain%specel(n)%Veloc1 = 0d0
          Tdomain%specel(n)%Veloc2 = 0d0
          Tdomain%specel(n)%Veloc3 = 0d0

          if(Tdomain%specel(n)%FPML)then
              allocate(Tdomain%specel(n)%Isx(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate(Tdomain%specel(n)%Isy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate(Tdomain%specel(n)%Isz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate(Tdomain%specel(n)%Ivx(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate(Tdomain%specel(n)%Ivy(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate(Tdomain%specel(n)%Ivz(0:ngllx-1,0:nglly-1,0:ngllz-1))
              allocate(Tdomain%specel(n)%I_Diagonal_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
              allocate(Tdomain%specel(n)%I_Diagonal_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
              allocate(Tdomain%specel(n)%I_Diagonal_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
              allocate(Tdomain%specel(n)%I_Residual_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
              allocate(Tdomain%specel(n)%I_Residual_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
              allocate(Tdomain%specel(n)%IVeloc1(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
              allocate(Tdomain%specel(n)%IVeloc2(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
              allocate(Tdomain%specel(n)%IVeloc3(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
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
          allocate(Tdomain%specel(n)%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:44))
          allocate(Tdomain%specel(n)%Displ(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))  
          Tdomain%specel(n)%Displ = 0    
      endif
  endif
 else   ! FLUID PART
  if(Tdomain%TimeD%velocity_scheme)then
      allocate(Tdomain%specel(n)%VelPhi(1:ngllx-2,1:nglly-2,1:ngllz-2))
      allocate(Tdomain%specel(n)%AccelPhi(1:ngllx-2,1:nglly-2,1:ngllz-2))
      allocate(Tdomain%specel(n)%VelPhi0(1:ngllx-2,1:nglly-2,1:ngllz-2))  
      allocate(Tdomain%specel(n)%ForcesFl(0:ngllx-1,0:nglly-1,0:ngllz-1))
      Tdomain%specel(n)%VelPhi = 0d0   
      Tdomain%specel(n)%AccelPhi = 0d0
      Tdomain%specel(n)%VelPhi0 = 0d0
      Tdomain%specel(n)%ForcesFl = 0d0
      if(Tdomain%specel(n)%PML)then
          allocate(Tdomain%specel(n)%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:17))
          allocate(Tdomain%specel(n)%Veloc1(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
          allocate(Tdomain%specel(n)%Veloc2(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
          allocate(Tdomain%specel(n)%Veloc3(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
          allocate(Tdomain%specel(n)%VelPhi1(1:ngllx-2,1:nglly-2,1:ngllz-2))
          allocate(Tdomain%specel(n)%VelPhi2(1:ngllx-2,1:nglly-2,1:ngllz-2))
          allocate(Tdomain%specel(n)%VelPhi3(1:ngllx-2,1:nglly-2,1:ngllz-2))
          allocate(Tdomain%specel(n)%ForcesFl1(0:ngllx-1,0:nglly-1,0:ngllz-1))
          allocate(Tdomain%specel(n)%ForcesFl2(0:ngllx-1,0:nglly-1,0:ngllz-1))
          allocate(Tdomain%specel(n)%ForcesFl3(0:ngllx-1,0:nglly-1,0:ngllz-1))
          allocate(Tdomain%specel(n)%DumpSx(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
          allocate(Tdomain%specel(n)%DumpSy(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
          allocate(Tdomain%specel(n)%DumpSz(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
          allocate(Tdomain%specel(n)%DumpMass(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
          allocate(Tdomain%specel(n)%DumpVx(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1))
          allocate(Tdomain%specel(n)%DumpVy(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1))
          allocate(Tdomain%specel(n)%DumpVz(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1))
          Tdomain%specel(n)%Veloc1 = 0d0
          Tdomain%specel(n)%Veloc2 = 0d0
          Tdomain%specel(n)%Veloc3 = 0d0
          Tdomain%specel(n)%VelPhi1 = 0d0
          Tdomain%specel(n)%VelPhi2 = 0d0
          Tdomain%specel(n)%VelPhi3 = 0d0
      else
          allocate(Tdomain%specel(n)%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:5))
          allocate(Tdomain%specel(n)%Phi(1:ngllx-2,1:nglly-2,1:ngllz-2))  
          Tdomain%specel(n)%Phi = 0    
      endif
  endif

 end if
enddo


do n = 0, Tdomain%n_face-1
 ngll1 = Tdomain%sFace(n)%ngll1 ; ngll2 = Tdomain%sFace(n)%ngll2
 if(Tdomain%sFace(n)%solid)then    ! SOLID PART
  allocate(Tdomain%sFace(n)%MassMat(1:ngll1-2,1:ngll2-2))
  allocate(Tdomain%sFace(n)%Veloc(1:ngll1-2,1:ngll2-2,0:2))
  allocate(Tdomain%sFace(n)%Forces(1:ngll1-2,1:ngll2-2,0:2))
  allocate(Tdomain%sFace(n)%Accel(1:ngll1-2,1:ngll2-2,0:2))
  allocate(Tdomain%sFace(n)%V0(1:ngll1-2,1:ngll2-2,0:2))  
  Tdomain%sFace(n)%MassMat = 0d0
  Tdomain%sFace(n)%Veloc = 0d0
  Tdomain%sFace(n)%Accel = 0d0
  Tdomain%sFace(n)%V0 = 0d0
  Tdomain%sFace(n)%Forces = 0d0  
  if(Tdomain%sFace(n)%PML)then 
        allocate(Tdomain%sFace(n)%Forces1(1:ngll1-2,1:ngll2-2,0:2)) 
        allocate(Tdomain%sFace(n)%Forces2(1:ngll1-2,1:ngll2-2,0:2)) 
        allocate(Tdomain%sFace(n)%Forces3(1:ngll1-2,1:ngll2-2,0:2)) 
        allocate(Tdomain%sFace(n)%DumpMass(1:ngll1-2,1:ngll2-2,0:2))
        allocate(Tdomain%sFace(n)%Veloc1(1:ngll1-2,1:ngll2-2,0:2)) 
        allocate(Tdomain%sFace(n)%Veloc2(1:ngll1-2,1:ngll2-2,0:2))
        allocate(Tdomain%sFace(n)%Veloc3(1:ngll1-2,1:ngll2-2,0:2))  
        allocate(Tdomain%sFace(n)%DumpVx(1:ngll1-2,1:ngll2-2,0:1)) 
        allocate(Tdomain%sFace(n)%DumpVy(1:ngll1-2,1:ngll2-2,0:1)) 
        allocate(Tdomain%sFace(n)%DumpVz(1:ngll1-2,1:ngll2-2,0:1)) 
        Tdomain%sFace(n)%DumpMass = 0.
        Tdomain%sFace(n)%Veloc1 = 0.
        Tdomain%sFace(n)%Veloc2 = 0.
        Tdomain%sFace(n)%Veloc3 = 0.
        if(Tdomain%sFace(n)%FPML)then
            allocate(Tdomain%sFace(n)%Ivx(1:ngll1-2,1:ngll2-1))
            allocate(Tdomain%sFace(n)%Ivy(1:ngll1-2,1:ngll2-1))
            allocate(Tdomain%sFace(n)%Ivz(1:ngll1-2,1:ngll2-1))
            allocate(Tdomain%sFace(n)%IVeloc1(1:ngll1-2,1:ngll2-1,0:2))
            allocate(Tdomain%sFace(n)%IVeloc2(1:ngll1-2,1:ngll2-1,0:2))
            allocate(Tdomain%sFace(n)%IVeloc3(1:ngll1-2,1:ngll2-1,0:2))
            Tdomain%sFace(n)%IVeloc1 = 0.
            Tdomain%sFace(n)%IVeloc2 = 0.
            Tdomain%sFace(n)%IVeloc3 = 0.
            Tdomain%sFace(n)%Ivx = 0.
            Tdomain%sFace(n)%Ivy = 0.
            Tdomain%sFace(n)%Ivz = 0.
       endif
   else
       allocate(Tdomain%sFace(n)%Displ(1:ngll1-2,1:ngll2-2,0:2))
       Tdomain%sFace(n)%Displ = 0d0  
   endif
 else   ! FLUID PART
  allocate(Tdomain%sFace(n)%MassMat(1:ngll1-2,1:ngll2-2))
  allocate(Tdomain%sFace(n)%VelPhi(1:ngll1-2,1:ngll2-2))
  allocate(Tdomain%sFace(n)%ForcesFl(1:ngll1-2,1:ngll2-2))
  allocate(Tdomain%sFace(n)%AccelPhi(1:ngll1-2,1:ngll2-2))
  allocate(Tdomain%sFace(n)%VelPhi0(1:ngll1-2,1:ngll2-2))  
  Tdomain%sFace(n)%MassMat = 0d0
  Tdomain%sFace(n)%VelPhi = 0d0
  Tdomain%sFace(n)%AccelPhi = 0d0
  Tdomain%sFace(n)%VelPhi0 = 0d0
  Tdomain%sFace(n)%ForcesFl = 0d0  
  if(Tdomain%sFace(n)%PML)then 
        allocate(Tdomain%sFace(n)%ForcesFl1(1:ngll1-2,1:ngll2-2)) 
        allocate(Tdomain%sFace(n)%ForcesFl2(1:ngll1-2,1:ngll2-2)) 
        allocate(Tdomain%sFace(n)%ForcesFl3(1:ngll1-2,1:ngll2-2)) 
        allocate(Tdomain%sFace(n)%DumpMass(1:ngll1-2,1:ngll2-2,0:2))
        allocate(Tdomain%sFace(n)%VelPhi1(1:ngll1-2,1:ngll2-2)) 
        allocate(Tdomain%sFace(n)%VelPhi2(1:ngll1-2,1:ngll2-2))
        allocate(Tdomain%sFace(n)%VelPhi3(1:ngll1-2,1:ngll2-2))  
        allocate(Tdomain%sFace(n)%DumpVx(1:ngll1-2,1:ngll2-2,0:1)) 
        allocate(Tdomain%sFace(n)%DumpVy(1:ngll1-2,1:ngll2-2,0:1)) 
        allocate(Tdomain%sFace(n)%DumpVz(1:ngll1-2,1:ngll2-2,0:1)) 
        Tdomain%sFace(n)%DumpMass = 0.
        Tdomain%sFace(n)%VelPhi1 = 0.
        Tdomain%sFace(n)%VelPhi2 = 0.
        Tdomain%sFace(n)%VelPhi3 = 0.
   else
       allocate(Tdomain%sFace(n)%Phi(1:ngll1-2,1:ngll2-2))
       Tdomain%sFace(n)%Phi = 0d0  
   endif

 end if
enddo


do n = 0,Tdomain%n_edge-1
  ngll = Tdomain%sEdge(n)%ngll
  if(Tdomain%sEdge(n)%solid)then   ! SOLID PART
    allocate(Tdomain%sEdge(n)%MassMat(1:ngll-2))
    allocate(Tdomain%sEdge(n)%Veloc(1:ngll-2,0:2))
    allocate(Tdomain%sEdge(n)%Forces(1:ngll-2,0:2))
    allocate(Tdomain%sEdge(n)%Accel(1:ngll-2,0:2))
    allocate(Tdomain%sEdge(n)%V0(1:ngll-2,0:2))
    Tdomain%sEdge(n)%MassMat = 0d0
    Tdomain%sEdge(n)%Veloc = 0d0
    Tdomain%sEdge(n)%Accel = 0d0
    Tdomain%sEdge(n)%V0 = 0d0
    Tdomain%sEdge(n)%Forces = 0d0
    if(Tdomain%sEdge(n)%PML)then
        allocate(Tdomain%sEdge(n)%Forces1(1:ngll-2,0:2))
        allocate(Tdomain%sEdge(n)%Forces2(1:ngll-2,0:2))
        allocate(Tdomain%sEdge(n)%Forces3(1:ngll-2,0:2))
        allocate(Tdomain%sEdge(n)%DumpMass(1:ngll-2, 0:2))
        allocate(Tdomain%sEdge(n)%Veloc1(1:ngll-2,0:2))
        allocate(Tdomain%sEdge(n)%Veloc2(1:ngll-2,0:2))
        allocate(Tdomain%sEdge(n)%Veloc3(1:ngll-2,0:2))
        allocate(Tdomain%sEdge(n)%DumpVx(1:ngll-2,0:1))
        allocate(Tdomain%sEdge(n)%DumpVy(1:ngll-2,0:1))
        allocate(Tdomain%sEdge(n)%DumpVz(1:ngll-2,0:1))
        Tdomain%sEdge(n)%DumpMass = 0d0
        Tdomain%sEdge(n)%Veloc1 = 0d0
        Tdomain%sEdge(n)%Veloc2 = 0d0
        Tdomain%sEdge(n)%Veloc3 = 0d0
        if(Tdomain%sEdge(n)%FPML)then
            allocate(Tdomain%sEdge(n)%Ivx(1:ngll-2))
            allocate(Tdomain%sEdge(n)%Ivy(1:ngll-2))
            allocate(Tdomain%sEdge(n)%Ivz(1:ngll-2))
            allocate(Tdomain%sEdge(n)%IVeloc1(1:ngll-2,0:2))
            allocate(Tdomain%sEdge(n)%IVeloc2(1:ngll-2,0:2))
            allocate(Tdomain%sEdge(n)%IVeloc3(1:ngll-2,0:2))
            Tdomain%sEdge(n)%IVeloc1 = 0.
            Tdomain%sEdge(n)%IVeloc2 = 0.
            Tdomain%sEdge(n)%IVeloc3 = 0.
            Tdomain%sEdge(n)%Ivx = 0.
            Tdomain%sEdge(n)%Ivy = 0.
            Tdomain%sEdge(n)%Ivz = 0.
        endif
   else
        allocate(Tdomain%sEdge(n)%Displ(1:ngll-2,0:2))
        Tdomain%sEdge(n)%Displ = 0d0
    endif
 else   ! FLUID PART
    allocate(Tdomain%sEdge(n)%MassMat(1:ngll-2))
    allocate(Tdomain%sEdge(n)%VelPhi(1:ngll-2))
    allocate(Tdomain%sEdge(n)%ForcesFl(1:ngll-2))
    allocate(Tdomain%sEdge(n)%AccelPhi(1:ngll-2))
    allocate(Tdomain%sEdge(n)%VelPhi0(1:ngll-2))
    Tdomain%sEdge(n)%MassMat = 0d0
    Tdomain%sEdge(n)%VelPhi = 0d0
    Tdomain%sEdge(n)%AccelPhi = 0d0
    Tdomain%sEdge(n)%VelPhi0 = 0d0
    Tdomain%sEdge(n)%ForcesFl = 0d0
    if(Tdomain%sEdge(n)%PML)then
        allocate(Tdomain%sEdge(n)%ForcesFl1(1:ngll-2))
        allocate(Tdomain%sEdge(n)%ForcesFl2(1:ngll-2))
        allocate(Tdomain%sEdge(n)%ForcesFl3(1:ngll-2))
        allocate(Tdomain%sEdge(n)%DumpMass(1:ngll-2,0:2))
        allocate(Tdomain%sEdge(n)%VelPhi1(1:ngll-2))
        allocate(Tdomain%sEdge(n)%VelPhi2(1:ngll-2))
        allocate(Tdomain%sEdge(n)%VelPhi3(1:ngll-2))
        allocate(Tdomain%sEdge(n)%DumpVx(1:ngll-2,0:1))
        allocate(Tdomain%sEdge(n)%DumpVy(1:ngll-2,0:1))
        allocate(Tdomain%sEdge(n)%DumpVz(1:ngll-2,0:1))
        Tdomain%sEdge(n)%DumpMass = 0d0
        Tdomain%sEdge(n)%VelPhi1 = 0d0
        Tdomain%sEdge(n)%VelPhi2 = 0d0
        Tdomain%sEdge(n)%VelPhi3 = 0d0
   else
        allocate(Tdomain%sEdge(n)%Phi(1:ngll-2))
        Tdomain%sEdge(n)%Phi = 0d0
    endif

 end if
enddo


do n = 0,Tdomain%n_vertex-1
  if(Tdomain%sVertex(n)%solid)then  ! SOLID PART
    Tdomain%sVertex(n)%MassMat = 0d0
    Tdomain%sVertex(n)%Veloc = 0d0
    Tdomain%sVertex(n)%Accel = 0d0
    Tdomain%sVertex(n)%V0 = 0d0
    Tdomain%sVertex(n)%Forces = 0d0
    if(Tdomain%sVertex(n)%PML)then
        allocate(Tdomain%sVertex(n)%Veloc1(0:2))
        allocate(Tdomain%sVertex(n)%Veloc2(0:2))
        allocate(Tdomain%sVertex(n)%Veloc3(0:2))
        allocate(Tdomain%sVertex(n)%DumpVx(0:2))
        allocate(Tdomain%sVertex(n)%DumpVy(0:2))
        allocate(Tdomain%sVertex(n)%DumpVz(0:2))
        allocate(Tdomain%sVertex(n)%DumpMass(0:2))
        Tdomain%sVertex(n)%Veloc1 = 0d0
        Tdomain%sVertex(n)%Veloc2 = 0d0
        Tdomain%sVertex(n)%Veloc3 = 0d0
        Tdomain%sVertex(n)%DumpMass = 0d0
        if(Tdomain%sVertex(n)%FPML)then
             allocate(Tdomain%sVertex(n)%IVeloc1(0:2))
             allocate(Tdomain%sVertex(n)%IVeloc2(0:2))
             allocate(Tdomain%sVertex(n)%IVeloc3(0:2))
             allocate(Tdomain%sVertex(n)%Ivx(0:0))
             allocate(Tdomain%sVertex(n)%Ivy(0:0))
             allocate(Tdomain%sVertex(n)%Ivz(0:0))
             Tdomain%sVertex(n)%IVeloc1 = 0d0
             Tdomain%sVertex(n)%IVeloc2 = 0d0
             Tdomain%sVertex(n)%IVeloc3 = 0d0
             Tdomain%sVertex(n)%Ivx = 0d0
             Tdomain%sVertex(n)%Ivy = 0d0
             Tdomain%sVertex(n)%Ivz = 0d0
        endif
    else
        Tdomain%sVertex(n)%Displ = 0d0
    endif
  else   ! FLUID PART
    Tdomain%sVertex(n)%MassMat = 0d0
    Tdomain%sVertex(n)%VelPhi = 0d0
    Tdomain%sVertex(n)%AccelPhi = 0d0
    Tdomain%sVertex(n)%VelPhi0 = 0d0
    Tdomain%sVertex(n)%ForcesFl = 0d0
    if(Tdomain%sVertex(n)%PML)then
        allocate(Tdomain%sVertex(n)%DumpVx(0:2))
        allocate(Tdomain%sVertex(n)%DumpVy(0:2))
        allocate(Tdomain%sVertex(n)%DumpVz(0:2))
        allocate(Tdomain%sVertex(n)%DumpMass(0:2))
        Tdomain%sVertex(n)%VelPhi1 = 0d0
        Tdomain%sVertex(n)%VelPhi2 = 0d0
        Tdomain%sVertex(n)%VelPhi3 = 0d0
        Tdomain%sVertex(n)%DumpMass = 0d0
    else
        Tdomain%sVertex(n)%Phi = 0d0
    endif
  end if
enddo

! solid-fluid normal terms
if(Tdomain%logicD%SF_local_present)then
    do nf = 0,Tdomain%SF%SF_n_faces-1
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        allocate(Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,0:2))
        Tdomain%SF%SF_Face(nf)%BtN = 0d0
    enddo

    do ne = 0,Tdomain%SF%SF_n_edges-1
        ngll = Tdomain%SF%SF_Edge(ne)%ngll
        allocate(Tdomain%SF%SF_Edge(ne)%BtN(1:ngll-2,0:2))
        Tdomain%SF%SF_Edge(ne)%BtN = 0d0
    enddo

    do nv = 0,Tdomain%SF%SF_n_vertices-1
        Tdomain%SF%SF_Vertex(nv)%BtN = 0d0
    enddo

end if


!------------------------
! Inter-proc communications
!------------------------
if(Tdomain%n_proc > 1)then
    do n = 0,Tdomain%n_proc-1
        ngll = 0
        ngll_F = 0
        ngllPML = 0
        ngllPML_F = 0
        ngllNeu = 0
        ngllSF = 0
        do i = 0,Tdomain%sComm(n)%nb_faces-1
            nf = Tdomain%sComm(n)%faces(i)
            ngll1 = Tdomain%sFace(nf)%ngll1 ; ngll2 = Tdomain%sFace(nf)%ngll2 
            if(Tdomain%sFace(nf)%solid)then  ! solid face
                ngll = ngll + (ngll1-2)*(ngll2-2)
                if(Tdomain%sFace(nf)%PML)then
                    ngllPML = ngllPML + (ngll1-2)*(ngll2-2)
                endif
            else   ! fluid face
                ngll_F = ngll_F + (ngll1-2)*(ngll2-2)
                if(Tdomain%sFace(nf)%PML)then
                    ngllPML_F = ngllPML_F + (ngll1-2)*(ngll2-2)
                endif
            end if
        enddo
        do i = 0,Tdomain%sComm(n)%nb_edges-1
            ne = Tdomain%sComm(n)%edges(i)
            if(Tdomain%sEdge(ne)%solid)then  ! solid edge
                ngll = ngll + Tdomain%sEdge(ne)%ngll-2
                if(Tdomain%sEdge(ne)%PML) then
                    ngllPML = ngllPML + Tdomain%sEdge(ne)%ngll-2
                endif
            else   ! fluid edge
                ngll_F = ngll_F + Tdomain%sEdge(ne)%ngll-2
                if(Tdomain%sEdge(ne)%PML) then
                    ngllPML_F = ngllPML_F + Tdomain%sEdge(ne)%ngll-2
                endif
            end if
        enddo
        do i = 0,Tdomain%sComm(n)%nb_vertices-1
            nv = Tdomain%sComm(n)%vertices(i) 
            if(Tdomain%sVertex(nv)%solid)then  ! solid vertex
                ngll = ngll + 1
                if(Tdomain%sVertex(nv)%PML) then
                    ngllPML = ngllPML + 1
                endif
            else    ! fluid vertex
                ngll_F = ngll_F + 1
                if(Tdomain%sVertex(nv)%PML) then
                    ngllPML_F = ngllPML_F + 1
                endif
            end if
        enddo
  ! Neumann
        do i = 0,Tdomain%sComm(n)%Neu_ne_shared-1
            ne = Tdomain%sComm(n)%Neu_edges_shared(i)
            ngllNeu = ngllNeu + Tdomain%Neumann%Neu_Edge(ne)%ngll-2
        enddo
        do i = 0,Tdomain%sComm(n)%Neu_nv_shared-1
            ngllNeu = ngllNeu + 1
        enddo
  ! Solid/fluid
        do i = 0,Tdomain%sComm(n)%SF_nf_shared-1
            nf = Tdomain%sComm(n)%SF_faces_shared(i)
            ngll1 = Tdomain%SF%SF_Face(nf)%ngll1; ngll2 = Tdomain%SF%SF_Face(nf)%ngll2 
            ngllSF = ngllSF + (ngll1-2)*(ngll2-2)
        enddo
        do i = 0,Tdomain%sComm(n)%SF_ne_shared-1
            ne = Tdomain%sComm(n)%SF_edges_shared(i)
            ngllSF = ngllSF + Tdomain%SF%SF_Edge(ne)%ngll-2
        enddo
        do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
            ngllSF = ngllSF + 1
        enddo

 ! allocations for inter-proc communications
  ! for mass communications
        if(ngll+ngll_F > 0)then
            allocate(Tdomain%sComm(n)%Give(0:ngll+ngll_F-1))
            allocate(Tdomain%sComm(n)%Take(0:ngll+ngll_F-1))
        end if
        if(ngllPML+ngllPML_F > 0)then
            if(Tdomain%any_FPML)then
                allocate(Tdomain%sComm(n)%GivePML(0:ngllPML+ngllPML_F-1,0:5))
                allocate(Tdomain%sComm(n)%TakePML(0:ngllPML+ngllPML_F-1,0:5))
            else
                allocate(Tdomain%sComm(n)%GivePML(0:ngllPML+ngllPML_F-1,0:2))
                allocate(Tdomain%sComm(n)%TakePML(0:ngllPML+ngllPML_F-1,0:2))
            end if
        end if
  ! solid
        if(ngll > 0)then
            allocate(Tdomain%sComm(n)%GiveForces(0:ngll-1,0:2))
            allocate(Tdomain%sComm(n)%TakeForces(0:ngll-1,0:2))
        endif
        if(ngllPML > 0)then
            allocate(Tdomain%sComm(n)%GiveForcesPML(0:ngllPML-1,1:3,0:2))
            allocate(Tdomain%sComm(n)%TakeForcesPML(0:ngllPML-1,1:3,0:2))
        endif
  ! fluid
        if(ngll_F > 0)then
            allocate(Tdomain%sComm(n)%GiveForcesFl(0:ngll_F-1))
            allocate(Tdomain%sComm(n)%TakeForcesFl(0:ngll_F-1))
        endif
        if(ngllPML_F > 0)then
            allocate(Tdomain%sComm(n)%GiveForcesPMLFl(0:ngllPML_F-1,1:3))
            allocate(Tdomain%sComm(n)%TakeForcesPMLFl(0:ngllPML_F-1,1:3))
        endif
  ! Neumann
        if(ngllNeu > 0)then
            allocate(Tdomain%sComm(n)%GiveNeu(0:ngllNeu-1,0:2))
            allocate(Tdomain%sComm(n)%TakeNeu(0:ngllNeu-1,0:2))
        endif
  ! solid-fluid interface
        if(ngllSF > 0)then
            allocate(Tdomain%sComm(n)%GiveForcesSF_FtoS(0:ngllSF-1,0:2))
            allocate(Tdomain%sComm(n)%TakeForcesSF_FtoS(0:ngllSF-1,0:2))
            allocate(Tdomain%sComm(n)%GiveForcesSF_StoF(0:ngllSF-1))
            allocate(Tdomain%sComm(n)%TakeForcesSF_StoF(0:ngllSF-1))
        endif

        Tdomain%sComm(n)%ngll = ngll
        Tdomain%sComm(n)%ngll_F = ngll_F
        Tdomain%sComm(n)%ngll_tot = ngll+ngll_F
        Tdomain%sComm(n)%ngllPML = ngllPML
        Tdomain%sComm(n)%ngllPML_F = ngllPML_F
        Tdomain%sComm(n)%ngllPML_tot = ngllPML+ngllPML_F
        Tdomain%sComm(n)%ngllNeu = ngllNeu
        Tdomain%sComm(n)%ngllSF = ngllSF


  enddo
endif

return
end subroutine allocate_domain
