subroutine Define_Arrays(Tdomain, rg)

use sdomain

implicit none

include 'mpif.h'

type (domain), intent (INOUT), target :: Tdomain
integer, intent(IN) :: rg

integer :: n, mat, ngllx,nglly,ngllz, ngll1,ngll2, ngll,ngllPML,ngllSO, i,j,k, n_elem, nf,ne,nv,nv_aus, idef, code
integer :: shift, I_give_to, I_take_from, n_rings
integer, parameter :: etiquette=100
integer, dimension(mpi_status_size) :: statut
real :: vp, ri,rj,rk, dx
real, external :: pow
real, dimension (:,:,:), allocatable :: xix,xiy,xiz, etax,etay,etaz, zetax,zetay,zetaz, Jac
real, dimension (:,:,:), allocatable :: Rlam,Rmu,RKmod, Whei, LocMassMat, wx,wy,wz, Id, Store_Btn
logical :: sortie


do n = 0,Tdomain%n_elem-1
  
   mat = Tdomain%specel(n)%mat_index
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz

   Tdomain%specel(n)%Density = Tdomain%sSubDomain(mat)%Ddensity
   Tdomain%specel(n)%Lambda = Tdomain%sSubDomain(mat)%DLambda
   Tdomain%specel(n)%Mu = Tdomain%sSubDomain(mat)%DMu

   allocate (Jac(0:ngllx-1,0:nglly-1,0:ngllz-1))
      
   allocate (xix (0:ngllx-1,0:nglly-1,0:ngllz-1)) 
   allocate (xiy (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (xiz (0:ngllx-1,0:nglly-1,0:ngllz-1))    
  
   allocate (etax (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (etay (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (etaz (0:ngllx-1,0:nglly-1,0:ngllz-1))

   allocate (zetax (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (zetay (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (zetaz (0:ngllx-1,0:nglly-1,0:ngllz-1))

   allocate (Whei (0:ngllx-1,0:nglly-1,0:ngllz-1))

   allocate (RKmod (0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (Rlam(0:ngllx-1,0:nglly-1,0:ngllz-1))
   allocate (Rmu(0:ngllx-1,0:nglly-1,0:ngllz-1))
    
   do k = 0, ngllz -1 
       do j = 0,nglly-1 
           do i = 0,ngllx-1
                Whei (i,j,k) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwy(j) &
                               * Tdomain%sSubdomain(mat)%GLLwz(k)
           enddo
       enddo
   enddo
      
   xix = Tdomain%specel(n)%InvGrad(:,:,:,0,0)
   xiy = Tdomain%specel(n)%InvGrad(:,:,:,1,0)
   xiz = Tdomain%specel(n)%InvGrad(:,:,:,2,0)
      
   etax = Tdomain%specel(n)%InvGrad(:,:,:,0,1)
   etay = Tdomain%specel(n)%InvGrad(:,:,:,1,1)
   etaz = Tdomain%specel(n)%InvGrad(:,:,:,2,1)

   zetax = Tdomain%specel(n)%InvGrad(:,:,:,0,2)
   zetay = Tdomain%specel(n)%InvGrad(:,:,:,1,2)
   zetaz = Tdomain%specel(n)%InvGrad(:,:,:,2,2)
   
   Jac  = Tdomain%specel(n)%Jacob
   Rlam = Tdomain%specel(n)%Lambda
   Rmu  = Tdomain%specel(n)%Mu
   RKmod = Rlam + 2. * Rmu

   Tdomain%specel(n)%MassMat = Whei*Tdomain%specel(n)%Density*Jac
 
   ! PML implementation

   if (.not. Tdomain%specel(n)%PML) then

     Tdomain%specel(n)%Acoeff(:,:,:,0) = -Whei*(RKmod*xix**2+Rmu*(xiy**2+xiz**2))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,1) = -Whei*(RKmod*xix*etax+Rmu*(xiy*etay+xiz*etaz))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,2) = -Whei*(rKmod*xix*zetax+Rmu*(xiy*zetay+xiz*zetaz))*Jac 
     Tdomain%specel(n)%Acoeff(:,:,:,3) = -Whei*(Rlam+Rmu)*xix*xiy*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,4) = -Whei*(Rlam*xix*etay+Rmu*xiy*etax)*Jac  
     Tdomain%specel(n)%Acoeff(:,:,:,5) = -Whei*(Rlam*xix*zetay+Rmu*xiy*zetax)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,6) = -Whei*(Rlam+Rmu)*xix*xiz*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,7) = -Whei*(Rlam*xix*etaz+Rmu*xiz*etax)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,8) = -Whei*(Rlam*xix*zetaz+rmu*xiz*zetax)*Jac

     Tdomain%specel(n)%Acoeff(:,:,:,9) = -Whei*(RKmod*etax**2+Rmu* (etay**2+etaz**2))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,10) = -Whei*(RKmod*etax*zetax+Rmu* (etay*zetay+etaz*zetaz))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,11) = -Whei*(Rlam*etax*xiy+Rmu*etay*xix)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,12) = -Whei*(Rlam+Rmu)*etay*etax*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,13) = -Whei*(Rlam*etax*zetay+Rmu*etay*zetax)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,14) = -Whei*(Rlam*etax*xiz+Rmu*etaz*xix)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,15) = -Whei*(Rlam+Rmu)*etaz*etax*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,16) = -Whei*(Rlam*etax*zetaz+Rmu*etaz*zetax)*Jac

     Tdomain%specel(n)%Acoeff(:,:,:,17) = -Whei*(RKmod*zetax**2+Rmu*  (zetay**2+zetaz**2))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,18) = -Whei*(Rlam*zetax*xiy+Rmu*zetay*xix)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,19) = -Whei*(Rlam*zetax*etay+Rmu*zetay*etax)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,20) = -Whei*(Rlam+Rmu)*zetax*zetay*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,21) = -Whei*(Rlam*zetax*xiz+Rmu*zetaz*xix)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,22) = -Whei*(Rlam*zetax*etaz+Rmu*zetaz*etax)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,23) = -Whei*(Rlam+Rmu)*zetax*zetaz*Jac

     Tdomain%specel(n)%Acoeff(:,:,:,24) = -Whei*(RKmod*xiy**2+Rmu* (xix**2+xiz**2))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,25) = -Whei*(RKmod*xiy*etay+Rmu* (xix*etax+xiz*etaz))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,26) = -Whei*(RKmod*xiy*zetay+Rmu*  (xix*zetax+xiz*zetaz))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,27) = -Whei*(Rlam+Rmu)*xiy*xiz*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,28) = -Whei*(Rlam*etaz*xiy+Rmu*etay*xiz)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,29) = -Whei*(Rlam*zetaz*xiy+Rmu*zetay*xiz)*Jac

     Tdomain%specel(n)%Acoeff(:,:,:,30) = -Whei*(RKmod*etay**2+Rmu* (etax**2+etaz**2))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,31) = -Whei*(RKmod*zetay*etay+Rmu* (zetax*etax+zetaz*etaz))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,32) = -Whei*(Rlam*etay*xiz+Rmu*etaz*xiy)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,33) = -Whei*(Rlam+Rmu)*etay*etaz*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,34) = -Whei*(Rlam*zetaz*etay+Rmu*zetay*etaz)*Jac

     Tdomain%specel(n)%Acoeff(:,:,:,35) = -Whei*(RKmod*zetay**2+Rmu* (zetax**2+zetaz**2))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,36) = -Whei*(Rlam*xiz*zetay+Rmu*xiy*zetaz)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,37) = -Whei*(Rlam*zetay*etaz+Rmu*zetaz*etay)*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,38) = -Whei*(Rlam+Rmu)*zetay*zetaz*Jac

     Tdomain%specel(n)%Acoeff(:,:,:,39) = -Whei*(RKmod*xiz**2+Rmu*  (xix**2+xiy**2))*Jac 
     Tdomain%specel(n)%Acoeff(:,:,:,40) = -Whei*(RKmod*xiz*etaz+Rmu* (xix*etax+xiy*etay))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,41) = -Whei*(RKmod*xiz*zetaz+Rmu* (xix*zetax+xiy*zetay))*Jac
                                       
     Tdomain%specel(n)%Acoeff(:,:,:,42) = -Whei*(RKmod*etaz**2+Rmu* (etax**2+etay**2))*Jac
     Tdomain%specel(n)%Acoeff(:,:,:,43) = -Whei*(RKmod*zetaz*etaz+Rmu* (zetax*etax+zetay*etay))*Jac
                                       
     Tdomain%specel(n)%Acoeff(:,:,:,44) = -Whei*(RKmod*zetaz**2+Rmu* (zetax**2+zetay**2))*Jac
                                       
   else

     Tdomain%specel(n)%Acoeff(:,:,:,0) = RKmod *xix
     Tdomain%specel(n)%Acoeff(:,:,:,1) = RKmod *etax
     Tdomain%specel(n)%Acoeff(:,:,:,2) = RKmod *zetax

     Tdomain%specel(n)%Acoeff(:,:,:,3) = RLam *xiy
     Tdomain%specel(n)%Acoeff(:,:,:,4) = RLam *etay
     Tdomain%specel(n)%Acoeff(:,:,:,5) = RLam *zetay

     Tdomain%specel(n)%Acoeff(:,:,:,6) = RLam *xiz
     Tdomain%specel(n)%Acoeff(:,:,:,7) = RLam *etaz
     Tdomain%specel(n)%Acoeff(:,:,:,8) = RLam *zetaz

     Tdomain%specel(n)%Acoeff(:,:,:,9) = RLam *xix
     Tdomain%specel(n)%Acoeff(:,:,:,10) = RLam *etax
     Tdomain%specel(n)%Acoeff(:,:,:,11) = RLam *zetax

     Tdomain%specel(n)%Acoeff(:,:,:,12) = RKmod *xiy
     Tdomain%specel(n)%Acoeff(:,:,:,13) = RKmod *etay
     Tdomain%specel(n)%Acoeff(:,:,:,14) = RKmod *zetay

     Tdomain%specel(n)%Acoeff(:,:,:,15) = RKmod *xiz
     Tdomain%specel(n)%Acoeff(:,:,:,16) = RKmod *etaz
     Tdomain%specel(n)%Acoeff(:,:,:,17) = RKmod *zetaz

     Tdomain%specel(n)%Acoeff(:,:,:,18) = RMu *xix
     Tdomain%specel(n)%Acoeff(:,:,:,19) = RMu *etax
     Tdomain%specel(n)%Acoeff(:,:,:,20) = RMu *zetax

     Tdomain%specel(n)%Acoeff(:,:,:,21) = RMu *xiy
     Tdomain%specel(n)%Acoeff(:,:,:,22) = RMu *etay
     Tdomain%specel(n)%Acoeff(:,:,:,23) = RMu *zetay

     Tdomain%specel(n)%Acoeff(:,:,:,24) = RMu *xiz
     Tdomain%specel(n)%Acoeff(:,:,:,25) = RMu *etaz
     Tdomain%specel(n)%Acoeff(:,:,:,26) = RMu *zetaz

     Tdomain%specel(n)%Acoeff(:,:,:,27) = -Whei * xix * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,28) = -Whei * xiy * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,29) = -Whei * xiz * Jac

     Tdomain%specel(n)%Acoeff(:,:,:,30) = -Whei * etax * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,31) = -Whei * etay * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,32) = -Whei * etaz * Jac

     Tdomain%specel(n)%Acoeff(:,:,:,33) = -Whei * zetax * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,34) = -Whei * zetay * Jac
     Tdomain%specel(n)%Acoeff(:,:,:,35) = -Whei * zetaz * Jac

     allocate (wx (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (wy (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (wz (0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate (Id (0:ngllx-1,0:nglly-1,0:ngllz-1))

     if (Tdomain%sSubDomain(mat)%Px) then
         idef = Tdomain%specel(n)%Iglobnum(0,0,0)
         dx = Tdomain%GlobCoord(0,idef)
         idef = Tdomain%specel(n)%Iglobnum(ngllx-1,0,0)
         dx = abs(Tdomain%GlobCoord(0,idef) - dx) 
         if (Tdomain%sSubDomain(mat)%Left) then
            do i = 0,ngllx-1
               ri = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcx(ngllx-1-i)) * float(ngllx-1)
               vp = Rkmod(i,0,0) / Tdomain%specel(n)%Density(i,0,0)
               vp = sqrt(vp)
               wx(i,0:nglly-1,0:ngllz-1) = pow(ri, vp, ngllx-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                               Tdomain%sSubdomain(mat)%npow)
            enddo
         else 
            do i = 0,ngllx-1
               ri = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcx(i)) * float(ngllx-1)
               vp = Rkmod(i,0,0) / Tdomain%specel(n)%Density(i,0,0)
               vp = sqrt(vp)
               wx(i,0:nglly-1,0:ngllz-1) = pow(ri, vp, ngllx-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                               Tdomain%sSubdomain(mat)%npow)
            enddo
         endif
     else 
         wx = 0.
     endif  

     if (Tdomain%sSubDomain(mat)%Py) then
         idef = Tdomain%specel(n)%Iglobnum(0,0,0)
         dx = Tdomain%GlobCoord(1,idef)
         idef = Tdomain%specel(n)%Iglobnum(0,nglly-1,0)
         dx = abs(Tdomain%GlobCoord (1,idef) - dx) 
         if (Tdomain%sSubDomain(mat)%Forward) then
            do j = 0,nglly-1
               rj = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcy(nglly-1-j)) * float(nglly-1)
               vp = Rkmod(0,j,0) / Tdomain%specel(n)%Density(0,j,0)
               vp = sqrt(vp)
               wy(0:ngllx-1,j,0:ngllz-1) = pow(rj, vp, nglly-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                                Tdomain%sSubdomain(mat)%npow)
            enddo
         else 
            do j = 0,nglly-1
               rj = 0.5 * (1 + Tdomain%sSubDomain(mat)%GLLcy(j)) * float(nglly-1)
               vp = Rkmod(0,j,0) / Tdomain%specel(n)%Density(0,j,0)
               vp = sqrt(vp)
               wy(0:ngllx-1,j,0:ngllz-1) = pow(rj, vp, nglly-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                                Tdomain%sSubdomain(mat)%npow)
            enddo
         endif
     else 
         wy = 0.
     endif  

     if (Tdomain%sSubDomain(mat)%Pz) then
         idef = Tdomain%specel(n)%Iglobnum(0,0,0)
         dx = Tdomain%GlobCoord(2,idef)
         idef = Tdomain%specel(n)%Iglobnum(0,0,ngllz-1)
         dx = abs(Tdomain%GlobCoord(2,idef) - dx)
         if (Tdomain%sSubDomain(mat)%Down) then
            do k = 0,ngllz-1
               rk = 0.5 * (1 + Tdomain%sSubdomain(mat)%GLLcz(ngllz-1-k)) * float(ngllz-1)
               vp = Rkmod(0,0,k) / Tdomain%specel(n)%Density(0,0,k)
               vp = sqrt(vp)
               wz(0:ngllx-1,0:nglly-1,k) = pow(rk, vp, ngllz-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                               Tdomain%sSubdomain(mat)%npow)
            enddo
         else 
            do k = 0,ngllz-1
               rk = 0.5 * (1 + Tdomain%sSubdomain(mat)%GLLcz(k)) * float(ngllz-1)
               vp = Rkmod(0,0,k) / Tdomain%specel(n)%Density(0,0,k)
               vp = sqrt(vp)
               wz(0:ngllx-1,0:nglly-1,k) = pow(rk, vp, ngllz-1, dx, Tdomain%sSubdomain(mat)%Apow, &
                                              Tdomain%sSubdomain(mat)%npow)
            enddo
         endif
     else 
         wz = 0.
     endif

     Id = 1.

! Strong formulation for stresses
       if (Tdomain%specel(n)%FPML) then

         Tdomain%specel(n)%DumpSx(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wx * Tdomain%sSubdomain(mat)%freq 
         Tdomain%specel(n)%DumpSx (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSx (:,:,:,1)
         Tdomain%specel(n)%DumpSx (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wx * Tdomain%sSubdomain(mat)%freq) *  Tdomain%specel(n)%DumpSx(:,:,:,1)

         Tdomain%specel(n)%DumpSy(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wy * Tdomain%sSubdomain(mat)%freq 
         Tdomain%specel(n)%DumpSy (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSy (:,:,:,1)
         Tdomain%specel(n)%DumpSy (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wy * Tdomain%sSubdomain(mat)%freq) *  Tdomain%specel(n)%DumpSy(:,:,:,1)

         Tdomain%specel(n)%DumpSz(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wz * Tdomain%sSubdomain(mat)%freq 
         Tdomain%specel(n)%DumpSz (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSz (:,:,:,1)
         Tdomain%specel(n)%DumpSz (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wz * Tdomain%sSubdomain(mat)%freq) *  Tdomain%specel(n)%DumpSz(:,:,:,1)

	 Tdomain%specel(n)%DumpMass(:,:,:,0) =  Tdomain%specel(n)%Density * Whei * Jac * wx * & 
	     Tdomain%sSubdomain(mat)%Dt * 0.5 * Tdomain%sSubdomain(mat)%freq
	 Tdomain%specel(n)%DumpMass(:,:,:,1) =  Tdomain%specel(n)%Density * Whei * Jac * wy * & 
	     Tdomain%sSubdomain(mat)%Dt * 0.5 * Tdomain%sSubdomain(mat)%freq
	 Tdomain%specel(n)%DumpMass(:,:,:,2) =  Tdomain%specel(n)%Density * Whei * Jac * wz * & 
	     Tdomain%sSubdomain(mat)%Dt * 0.5 * Tdomain%sSubdomain(mat)%freq
	    
         Tdomain%specel(n)%Isx = Tdomain%sSubdomain(mat)%Dt * wx * Tdomain%sSubdomain(mat)%freq * Tdomain%specel(n)%DumpSx(:,:,:,1)     
         Tdomain%specel(n)%Isy = Tdomain%sSubdomain(mat)%Dt * wy * Tdomain%sSubdomain(mat)%freq * Tdomain%specel(n)%DumpSy(:,:,:,1)     
         Tdomain%specel(n)%Isz = Tdomain%sSubdomain(mat)%Dt * wz * Tdomain%sSubdomain(mat)%freq * Tdomain%specel(n)%DumpSz(:,:,:,1)     
	 	     
  	 Tdomain%specel(n)%Ivx = Tdomain%specel(n)%Density * Whei * Tdomain%sSubdomain(mat)%Dt * wx * Jac * Tdomain%sSubdomain(mat)%freq 
  	 Tdomain%specel(n)%Ivy = Tdomain%specel(n)%Density * Whei * Tdomain%sSubdomain(mat)%Dt * wy * Jac * Tdomain%sSubdomain(mat)%freq 
	 Tdomain%specel(n)%Ivz = Tdomain%specel(n)%Density * Whei * Tdomain%sSubdomain(mat)%Dt * wz * Jac * Tdomain%sSubdomain(mat)%freq 
       else

         Tdomain%specel(n)%DumpSx(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wx
         Tdomain%specel(n)%DumpSx (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSx (:,:,:,1)
         Tdomain%specel(n)%DumpSx (:,:,:,0) = (Id- Tdomain%sSubdomain(mat)%Dt * 0.5 * wx) * Tdomain%specel(n)%DumpSx(:,:,:,1)     

         Tdomain%specel(n)%DumpSy(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wy
         Tdomain%specel(n)%DumpSy (:,:,:,1) = 1./ Tdomain%specel(n)%DumpSy (:,:,:,1) 
         Tdomain%specel(n)%DumpSy (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wy) * Tdomain%specel(n)%DumpSy(:,:,:,1)

         Tdomain%specel(n)%DumpSz(:,:,:,1) = Id + 0.5 * Tdomain%sSubdomain(mat)%Dt * wz
         Tdomain%specel(n)%DumpSz(:,:,:,1)  = 1./ Tdomain%specel(n)%DumpSz (:,:,:,1) 
         Tdomain%specel(n)%DumpSz (:,:,:,0) = (Id - Tdomain%sSubdomain(mat)%Dt * 0.5 * wz) * Tdomain%specel(n)%DumpSz(:,:,:,1)
	 
         Tdomain%specel(n)%DumpMass(:,:,:,0) = 0.5 * Tdomain%specel(n)%Density * Whei * &
                                              Tdomain%sSubdomain(mat)%Dt * wx * Jac
         Tdomain%specel(n)%DumpMass(:,:,:,1) = 0.5 * Tdomain%specel(n)%Density * Whei * &
                                              Tdomain%sSubdomain(mat)%Dt * wy * Jac
         Tdomain%specel(n)%DumpMass(:,:,:,2) = 0.5 * Tdomain%specel(n)%Density * Whei * &
                                             Tdomain%sSubdomain(mat)%Dt * wz * Jac
      endif
   deallocate (wx,wy,wz,Id )

   endif

   deallocate (Jac, xix, xiy, xiz, etax, etay, etaz, zetax, zetay, zetaz, Whei, RKmod, Rmu, Rlam)  

enddo


! Mass and DumpMass Communications inside Processors
do n = 0,Tdomain%n_elem-1
    call get_Mass_Elem2Face(Tdomain,n)
    call get_Mass_Elem2Edge(Tdomain,n)
    call get_Mass_Elem2Vertex(Tdomain,n)
enddo

! Invert Mass Matrix expression
do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz
   allocate (LocMassMat(1:ngllx-2,1:nglly-2,1:ngllz-2))
   LocMassMat(:,:,:) = Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2)

   if (Tdomain%specel(n)%PML) then
      Tdomain%specel(n)%DumpVx (:,:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,0)
      Tdomain%specel(n)%DumpVx (:,:,:,1) = 1./ Tdomain%specel(n)%DumpVx (:,:,:,1) 
      Tdomain%specel(n)%DumpVx (:,:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,0)
      Tdomain%specel(n)%DumpVx (:,:,:,0) = Tdomain%specel(n)%DumpVx (:,:,:,0) * Tdomain%specel(n)%DumpVx (:,:,:,1)

      Tdomain%specel(n)%DumpVy (:,:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,1)
      Tdomain%specel(n)%DumpVy (:,:,:,1) = 1./ Tdomain%specel(n)%DumpVy (:,:,:,1) 
      Tdomain%specel(n)%DumpVy (:,:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,1)
      Tdomain%specel(n)%DumpVy (:,:,:,0) = Tdomain%specel(n)%DumpVy (:,:,:,0) * Tdomain%specel(n)%DumpVy (:,:,:,1)

      Tdomain%specel(n)%DumpVz (:,:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,2)
      Tdomain%specel(n)%DumpVz (:,:,:,1) = 1./ Tdomain%specel(n)%DumpVz (:,:,:,1) 
      Tdomain%specel(n)%DumpVz (:,:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,1:ngllz-2,2)
      Tdomain%specel(n)%DumpVz (:,:,:,0) = Tdomain%specel(n)%DumpVz (:,:,:,0) * Tdomain%specel(n)%DumpVz (:,:,:,1)

      if (Tdomain%specel(n)%FPML) then
        LocMassMat = Tdomain%specel(n)%Ivx(1:ngllx-2,1:nglly-2,1:ngllz-2)
        deallocate (Tdomain%specel(n)%Ivx)
	allocate (Tdomain%specel(n)%Ivx(1:ngllx-2,1:nglly-2,1:ngllz-2) )
	Tdomain%specel(n)%Ivx = LocMassMat *  Tdomain%specel(n)%DumpVx (:,:,:,1)
        LocMassMat = Tdomain%specel(n)%Ivy(1:ngllx-2,1:nglly-2,1:ngllz-2)
        deallocate (Tdomain%specel(n)%Ivy)
	allocate (Tdomain%specel(n)%Ivy(1:ngllx-2,1:nglly-2,1:ngllz-2) )
	Tdomain%specel(n)%Ivy = LocMassMat *  Tdomain%specel(n)%DumpVy (:,:,:,1)
        LocMassMat = Tdomain%specel(n)%Ivz(1:ngllx-2,1:nglly-2,1:ngllz-2)
        deallocate (Tdomain%specel(n)%Ivz)
	allocate (Tdomain%specel(n)%Ivz(1:ngllx-2,1:nglly-2,1:ngllz-2) )
	Tdomain%specel(n)%Ivz = LocMassMat *  Tdomain%specel(n)%DumpVz (:,:,:,1)
      endif
      deallocate (Tdomain%specel(n)%DumpMass)
   endif

   LocMassMat(:,:,:) = Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2)

   LocMassmat = 1./ LocMassMat
   deallocate (Tdomain%specel(n)%MassMat) 
   allocate (Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,1:ngllz-2) )
   Tdomain%specel(n)%MassMat = LocMassMat
   deallocate (LocMassMat)
   deallocate (Tdomain%specel(n)%Lambda)
   deallocate (Tdomain%specel(n)%Mu)  
   deallocate (Tdomain%specel(n)%InvGrad)

enddo


! Define super objects properties (Btn)
if (Tdomain%logicD%super_object_local_present) then
 if (Tdomain%super_object_type == "P") then

   do nf = 0, Tdomain%sPlaneW%n_faces-1
      ngll1 = Tdomain%sPlaneW%pFace(nf)%ngll1
      ngll2 = Tdomain%sPlaneW%pFace(nf)%ngll2
      mat = Tdomain%sPlaneW%pFace(nf)%mat_index

      if (Tdomain%sPlaneW%pFace(nf)%dir == 0 .or. Tdomain%sPlaneW%pFace(nf)%dir == 5 ) then
          do j = 0,ngll2-1
            do i = 0,ngll1-1
              Tdomain%sPlaneW%pFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwy(j)* &
                                                       Tdomain%sPlaneW%pFace(nf)%normal(i,j,0:2)
            enddo
          enddo
      else if (Tdomain%sPlaneW%pFace(nf)%dir ==1 .or. Tdomain%sPlaneW%pFace(nf)%dir ==3 ) then
          do j = 0,ngll2-1
            do i = 0,ngll1-1
              Tdomain%sPlaneW%pFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwz(j)* &
                                                       Tdomain%sPlaneW%pFace(nf)%normal(i,j,0:2)
            enddo
          enddo
      else
          do j = 0,ngll2-1
            do i = 0,ngll1-1
              Tdomain%sPlaneW%pFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwy(i) * Tdomain%sSubdomain(mat)%GLLwz(j)* &
                                                       Tdomain%sPlaneW%pFace(nf)%normal(i,j,0:2)
            enddo
          enddo
      endif
      
      ! Internal communication of Btn
      do i = 0,3   
        ne = Tdomain%sPlaneW%pFace(nf)%Near_Edges(i)
        ngll = Tdomain%sPlaneW%pEdge(ne)%ngll
        if ( Tdomain%sPlaneW%pFace(nf)%Orient_Edges(i) == 0 ) then
          select case (i)
            case (0)          
              Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) + &
                                                            Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,0,0:2)
            case (1)          
              Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) + &
                                                            Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1,1:ngll2-2,0:2)
            case (2)          
              Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) + &
                                                            Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,ngll2-1,0:2)
            case (3)          
              Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2) + &
                                                            Tdomain%sPlaneW%pFace(nf)%Btn(0,1:ngll2-2,0:2)
          end select
        else  
          select case (i)
            case (0)         
             do j=1,ngll-2
              Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) + &
                                                            Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1-j,0,0:2)
             enddo
            case (1) 
             do j=1,ngll-2
              Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) + &
                                                            Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1,ngll2-1-j,0:2)
             enddo
            case (2)
             do j=1,ngll-2          
              Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) + &
                                                            Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1-j,ngll2-1,0:2)
             enddo
            case (3)
             do j=1,ngll-2          
              Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) + &
                                                            Tdomain%sPlaneW%pFace(nf)%Btn(0,ngll2-1-j,0:2)
             enddo
          end select
        endif

        nv = Tdomain%sPlaneW%pFace(nf)%Near_Vertices(i)
        select case (i)
            case (0)          
              Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) + &
                                                     Tdomain%sPlaneW%pFace(nf)%Btn(0,0,0:2)
            case (1)          
              Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) + &
                                                     Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1,0,0:2)
            case (2)          
              Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) + &
                                                     Tdomain%sPlaneW%pFace(nf)%Btn(ngll1-1,ngll2-1,0:2)
            case (3)          
              Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) + &
                                                     Tdomain%sPlaneW%pFace(nf)%Btn(0,ngll2-1,0:2)
        end select
      enddo

      allocate (Store_Btn(1:ngll1-2,1:ngll2-2,0:2))
      Store_Btn (1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2)
      deallocate(Tdomain%sPlaneW%pFace(nf)%Btn)
      allocate(Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2))
      Tdomain%sPlaneW%pFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2) = Store_Btn(1:ngll1-2,1:ngll2-2,0:2)
      deallocate (Store_Btn)
      deallocate (Tdomain%sPlaneW%pFace(nf)%normal)

   enddo

 endif
endif


! Define Neumann properties (Btn)
if (Tdomain%logicD%neumann_local_present) then

   do nf = 0, Tdomain%sNeu%n_faces-1

      ngll1 = Tdomain%sNeu%nFace(nf)%ngll1
      ngll2 = Tdomain%sNeu%nFace(nf)%ngll2
      mat = Tdomain%sNeu%nFace(nf)%mat_index

      if (Tdomain%sNeu%nFace(nf)%dir == 0 .or. Tdomain%sNeu%nFace(nf)%dir == 5 ) then
          do j = 0,ngll2-1
            do i = 0,ngll1-1
              Tdomain%sNeu%nFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwy(j)* &
                                                       Tdomain%sNeu%nFace(nf)%normal(i,j,0:2)
            enddo
          enddo
      else if (Tdomain%sNeu%nFace(nf)%dir ==1 .or. Tdomain%sNeu%nFace(nf)%dir ==3 ) then
          do j = 0,ngll2-1
            do i = 0,ngll1-1
              Tdomain%sNeu%nFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwx(i) * Tdomain%sSubdomain(mat)%GLLwz(j)* &
                                                       Tdomain%sNeu%nFace(nf)%normal(i,j,0:2)
            enddo
          enddo
      else
          do j = 0,ngll2-1
            do i = 0,ngll1-1
              Tdomain%sNeu%nFace(nf)%Btn(i,j,0:2) = Tdomain%sSubdomain(mat)%GLLwy(i) * Tdomain%sSubdomain(mat)%GLLwz(j)* &
                                                       Tdomain%sNeu%nFace(nf)%normal(i,j,0:2)
            enddo
          enddo
      endif
      
      ! Internal communication of Btn
      do i = 0,3   
        ne = Tdomain%sNeu%nFace(nf)%Near_Edges(i)
        ngll = Tdomain%sNeu%nEdge(ne)%ngll
        if ( Tdomain%sNeu%nFace(nf)%Orient_Edges(i) == 0 ) then
          select case (i)
            case (0)          
              Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) + &
                                                            Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,0,0:2)
            case (1)          
              Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) + &
                                                            Tdomain%sNeu%nFace(nf)%Btn(ngll1-1,1:ngll2-2,0:2)
            case (2)          
              Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) + &
                                                            Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,ngll2-1,0:2)
            case (3)          
              Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2) + &
                                                            Tdomain%sNeu%nFace(nf)%Btn(0,1:ngll2-2,0:2)
          end select
        else  
          select case (i)
            case (0)         
             do j=1,ngll-2
              Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + &
                                                            Tdomain%sNeu%nFace(nf)%Btn(ngll1-1-j,0,0:2)
             enddo
            case (1) 
             do j=1,ngll-2
              Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + &
                                                            Tdomain%sNeu%nFace(nf)%Btn(ngll1-1,ngll2-1-j,0:2)
             enddo
            case (2)
             do j=1,ngll-2          
              Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + &
                                                            Tdomain%sNeu%nFace(nf)%Btn(ngll1-1-j,ngll2-1,0:2)
             enddo
            case (3)
             do j=1,ngll-2          
              Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + &
                                                            Tdomain%sNeu%nFace(nf)%Btn(0,ngll2-1-j,0:2)
             enddo
          end select
        endif

        nv = Tdomain%sNeu%nFace(nf)%Near_Vertices(i)
        select case (i)
            case (0)          
              Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + &
                                                     Tdomain%sNeu%nFace(nf)%Btn(0,0,0:2)
            case (1)          
              Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + &
                                                     Tdomain%sNeu%nFace(nf)%Btn(ngll1-1,0,0:2)
            case (2)          
              Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + &
                                                     Tdomain%sNeu%nFace(nf)%Btn(ngll1-1,ngll2-1,0:2)
            case (3)          
              Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + &
                                                     Tdomain%sNeu%nFace(nf)%Btn(0,ngll2-1,0:2)
        end select
      enddo

      allocate (Store_Btn(1:ngll1-2,1:ngll2-2,0:2))
      Store_Btn (1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2)
      deallocate(Tdomain%sNeu%nFace(nf)%Btn)
      allocate(Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2))
      Tdomain%sNeu%nFace(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2) = Store_Btn(1:ngll1-2,1:ngll2-2,0:2)
      deallocate (Store_Btn)
      deallocate (Tdomain%sNeu%nFace(nf)%normal)

   enddo
endif

! MPI communications
if ( Tdomain%n_proc > 1 ) then
do n = 0,Tdomain%n_proc-1
    ngll = 0
    ngllPML = 0
    ngllSO = 0
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        do j = 1,Tdomain%sFace(nf)%ngll2-2
            do k = 1,Tdomain%sFace(nf)%ngll1-2
                Tdomain%sComm(n)%Give(ngll) = Tdomain%sFace(nf)%MassMat(k,j)
                ngll = ngll + 1
            enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sFace(nf)%DumpMass(k,j,0:2)
                    if (Tdomain%any_FPML) then
                          Tdomain%sComm(n)%GivePML(ngllPML,3) =  Tdomain%sFace(nf)%Ivx(k,j)
                          Tdomain%sComm(n)%GivePML(ngllPML,4) =  Tdomain%sFace(nf)%Ivy(k,j)
                          Tdomain%sComm(n)%GivePML(ngllPML,5) =  Tdomain%sFace(nf)%Ivz(k,j)
                    endif
                    ngllPML = ngllPML + 1
                enddo
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        do j = 1,Tdomain%sEdge(ne)%ngll-2
            Tdomain%sComm(n)%Give(ngll) = Tdomain%sEdge(ne)%MassMat(j)
            ngll = ngll + 1
        enddo
        if (Tdomain%sEdge(ne)%PML) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sEdge(ne)%DumpMass(j,0:2)
                if (Tdomain%any_FPML) then
                          Tdomain%sComm(n)%GivePML(ngllPML,3) =  Tdomain%sEdge(ne)%Ivx(j)
                          Tdomain%sComm(n)%GivePML(ngllPML,4) =  Tdomain%sEdge(ne)%Ivy(j)
                          Tdomain%sComm(n)%GivePML(ngllPML,5) =  Tdomain%sEdge(ne)%Ivz(j)
                endif
                ngllPML = ngllPML + 1
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        Tdomain%sComm(n)%Give(ngll) = Tdomain%svertex(nv)%MassMat
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%PML) then
            Tdomain%sComm(n)%GivePML(ngllPML,0:2) = Tdomain%sVertex(nv)%DumpMass(0:2)
                if (Tdomain%any_FPML) then
                          Tdomain%sComm(n)%GivePML(ngllPML,3) =  Tdomain%sVertex(nv)%Ivx(0)
                          Tdomain%sComm(n)%GivePML(ngllPML,4) =  Tdomain%sVertex(nv)%Ivy(0)
                          Tdomain%sComm(n)%GivePML(ngllPML,5) =  Tdomain%sVertex(nv)%Ivz(0)
                endif
            ngllPML = ngllPML + 1
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges_so-1
        ne = Tdomain%sComm(n)%edges_SO(i)
        do j = 1,Tdomain%sPlaneW%pEdge(ne)%ngll-2
            Tdomain%sComm(n)%GiveSO(ngllSO,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) 
            ngllSO = ngllSO + 1
        enddo
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices_so-1
       nv = Tdomain%sComm(n)%vertices_SO(i)
       Tdomain%sComm(n)%GiveSO(ngllSO,0:2) = Tdomain%sPlaneW%pVertex(nv)%Btn(0:2)  
       ngllSO = ngllSO + 1    
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges_neu-1
        ne = Tdomain%sComm(n)%edges_Neu(i)
        do j = 1,Tdomain%sNeu%nEdge(ne)%ngll-2
            Tdomain%sComm(n)%GiveSO(ngllSO,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2)
            ngllSO = ngllSO + 1
        enddo
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices_neu-1
       nv = Tdomain%sComm(n)%vertices_Neu(i)
       Tdomain%sComm(n)%GiveSO(ngllSO,0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2)
       ngllSO = ngllSO + 1    
    enddo
enddo

n = Tdomain%n_proc
do shift = 1,n-1
    I_give_to = rg + shift
    if (I_give_to > n-1)   I_give_to = I_give_to - n
    I_take_from = rg - shift
    if (I_take_from < 0)   I_take_from = I_take_from + n
    if (mod(n,shift)==0 .and. shift/=1) then
        n_rings = shift
    else if (mod(n,n-shift)==0 .and. shift/=n-1) then
        n_rings = n-shift
    else if (mod(n,2)==0 .and. mod(shift,2)==0) then
        n_rings = 2
    else
        n_rings = 1
    endif
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngll>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%Give, Tdomain%sComm(I_give_to)%ngll, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngll>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%Take, Tdomain%sComm(I_take_from)%ngll, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngll>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%Take, Tdomain%sComm(I_take_from)%ngll, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngll>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%Give, Tdomain%sComm(I_give_to)%ngll, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
    call MPI_BARRIER (MPI_COMM_WORLD, code)
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngllPML>0) then
                if (Tdomain%any_FPML) then
                    call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 6*Tdomain%sComm(I_give_to)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                else
                    call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 3*Tdomain%sComm(I_give_to)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                endif
            endif
            if (Tdomain%sComm(I_take_from)%ngllPML>0) then
                if (Tdomain%any_FPML) then
                       call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 6*Tdomain%sComm(I_take_from)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                else
                       call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 3*Tdomain%sComm(I_take_from)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                endif
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngllPML>0) then
                         if (Tdomain%any_FPML) then
                              call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 6*Tdomain%sComm(I_take_from)%ngllPML, &
                                   MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                        else
                              call MPI_RECV (Tdomain%sComm(I_take_from)%TakePML, 3*Tdomain%sComm(I_take_from)%ngllPML, &
                                   MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                        endif
                    endif
                    if (Tdomain%sComm(I_give_to)%ngllPML>0) then
                      if (Tdomain%any_FPML) then
                         call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 6*Tdomain%sComm(I_give_to)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                      else
                         call MPI_SEND (Tdomain%sComm(I_give_to)%GivePML, 3*Tdomain%sComm(I_give_to)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                      endif
                    endif
                endif
            enddo
        endif
    enddo
    call MPI_BARRIER (MPI_COMM_WORLD, code)
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngllSO>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%GiveSO, 3*Tdomain%sComm(I_give_to)%ngllSO, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngllSO>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%TakeSO, 3*Tdomain%sComm(I_take_from)%ngllSO, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngllSO>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%TakeSO, 3*Tdomain%sComm(I_take_from)%ngllSO, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngllSO>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%GiveSO, 3*Tdomain%sComm(I_give_to)%ngllSO, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
    call MPI_BARRIER (MPI_COMM_WORLD, code)
enddo

do n = 0,Tdomain%n_proc-1
  ngll = 0
  ngllPML = 0
  ngllSO = 0
  call Comm_Mass_Face(Tdomain,n,ngll,ngllPML)
  call Comm_Mass_Edge(Tdomain,n,ngll,ngllPML)
  call Comm_Mass_Vertex(Tdomain,n,ngll,ngllPML)

  ! Super Object  (PlaneW or Fault) 
  do i = 0,Tdomain%sComm(n)%nb_edges_so-1
    ne = Tdomain%sComm(n)%edges_SO(i) 
    ngll1 = Tdomain%sPlaneW%pEdge(ne)%ngll
    if ( Tdomain%sComm(n)%orient_edges_SO(i) == 0 ) then 
      do j = 1,Tdomain%sPlaneW%pEdge(ne)%ngll-2
        Tdomain%sPlaneW%pEdge(ne)%Btn(j,0:2) = Tdomain%splaneW%pEdge(ne)%Btn(j,0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
        ngllSO = ngllSO + 1
      enddo
    else if ( Tdomain%sComm(n)%orient_edges_SO(i) == 1 ) then 
      do j = 1,Tdomain%splaneW%pEdge(ne)%ngll-2
        Tdomain%sPlaneW%pEdge(ne)%Btn(ngll1-1-j,0:2) = Tdomain%sPlaneW%pEdge(ne)%Btn(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
        ngllSO = ngllSO + 1
      enddo
    else 
      print*,'Pb with coherency number for edge' 
    endif    
  enddo
  do i = 0,Tdomain%sComm(n)%nb_vertices_so-1
    nv = Tdomain%sComm(n)%vertices_SO(i) 
    Tdomain%sPlaneW%pVertex(nv)%Btn(0:2) = Tdomain%splaneW%pVertex(nv)%Btn(0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
    ngllSO = ngllSO + 1
  enddo

  ! Neumann
  do i = 0,Tdomain%sComm(n)%nb_edges_neu-1
    ne = Tdomain%sComm(n)%edges_Neu(i) 
    ngll1 = Tdomain%sNeu%nEdge(ne)%ngll
    if ( Tdomain%sComm(n)%orient_edges_Neu(i) == 0 ) then 
      do j = 1,Tdomain%sNeu%nEdge(ne)%ngll-2
        Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(j,0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
        ngllSO = ngllSO + 1
      enddo
    else if ( Tdomain%sComm(n)%orient_edges_Neu(i) == 1 ) then 
      do j = 1,Tdomain%sNeu%nEdge(ne)%ngll-2
        Tdomain%sNeu%nEdge(ne)%Btn(ngll1-1-j,0:2) = Tdomain%sNeu%nEdge(ne)%Btn(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
        ngllSO = ngllSO + 1
      enddo
    else 
      print*,'Pb with coherency number for edge' 
    endif 
  enddo
  do i = 0,Tdomain%sComm(n)%nb_vertices_neu-1
    nv = Tdomain%sComm(n)%vertices_Neu(i) 
    Tdomain%sNeu%nVertex(nv)%Btn(0:2) = Tdomain%sNeu%nVertex(nv)%Btn(0:2) + Tdomain%sComm(n)%TakeSO(ngllSO,0:2)
    ngllSO = ngllSO + 1
  enddo

  if (Tdomain%sComm(n)%ngll>0) then
      deallocate (Tdomain%sComm(n)%Give)
      deallocate (Tdomain%sComm(n)%Take)
  endif
  if (Tdomain%sComm(n)%ngllPML>0) then
      deallocate (Tdomain%sComm(n)%GivePML)
      deallocate (Tdomain%sComm(n)%TakePML)
  endif
  if (Tdomain%sComm(n)%ngllSO>0) then
      deallocate (Tdomain%sComm(n)%GiveSO)
      deallocate (Tdomain%sComm(n)%TakeSO)
  endif

enddo
endif

do nf = 0,Tdomain%n_face-1
    if (Tdomain%sFace(nf)%PML) then
        Tdomain%sFace(nf)%DumpVx(:,:,1) = Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,:,0)
        Tdomain%sFace(nf)%DumpVx(:,:,1) = 1./Tdomain%sFace(nf)%DumpVx(:,:,1)
       Tdomain%sFace(nf)%DumpVx(:,:,0) = Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,:,0)
        Tdomain%sFace(nf)%DumpVx(:,:,0) = Tdomain%sFace(nf)%DumpVx(:,:,0) * Tdomain%sFace(nf)%DumpVx(:,:,1)     

        Tdomain%sFace(nf)%DumpVy(:,:,1) = Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,:,1)
        Tdomain%sFace(nf)%DumpVy(:,:,1) = 1./Tdomain%sFace(nf)%DumpVy (:,:,1)
        Tdomain%sFace(nf)%DumpVy(:,:,0) = Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,:,1)
        Tdomain%sFace(nf)%DumpVy(:,:,0) = Tdomain%sFace(nf)%DumpVy(:,:,0) * Tdomain%sFace(nf)%DumpVy(:,:,1)         

        Tdomain%sFace(nf)%DumpVz(:,:,1) = Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,:,2)
        Tdomain%sFace(nf)%DumpVz(:,:,1) = 1./Tdomain%sFace(nf)%DumpVz(:,:,1)
        Tdomain%sFace(nf)%DumpVz(:,:,0) = Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,:,2)
        Tdomain%sFace(nf)%DumpVz(:,:,0) = Tdomain%sFace(nf)%DumpVz(:,:,0) * Tdomain%sFace(nf)%DumpVz(:,:,1)      

        if (Tdomain%sFace(nf)%FPML) then
              Tdomain%sFace(nf)%Ivx(:,:) = Tdomain%sFace(nf)%Ivx(:,:) *  Tdomain%sFace(nf)%DumpVx(:,:,1)
              Tdomain%sFace(nf)%Ivy(:,:) = Tdomain%sFace(nf)%Ivy(:,:) *  Tdomain%sFace(nf)%DumpVy(:,:,1)
              Tdomain%sFace(nf)%Ivz(:,:) = Tdomain%sFace(nf)%Ivz(:,:) *  Tdomain%sFace(nf)%DumpVz(:,:,1)
       endif
         deallocate (Tdomain%sFace(nf)%DumpMass)
    endif
    Tdomain%sFace(nf)%MassMat = 1./ Tdomain%sFace(nf)%MassMat
enddo

do ne = 0,Tdomain%n_edge-1
    if (Tdomain%sEdge(ne)%PML) then
        Tdomain%sEdge(ne)%DumpVx(:,1) = Tdomain%sEdge(ne)%MassMat + Tdomain%sEdge(ne)%DumpMass(:,0)
        Tdomain%sEdge(ne)%DumpVx(:,1) = 1./Tdomain%sEdge(ne)%DumpVx(:,1)
        Tdomain%sEdge(ne)%DumpVx(:,0) = Tdomain%sEdge(ne)%MassMat - Tdomain%sEdge(ne)%DumpMass(:,0)
        Tdomain%sEdge(ne)%DumpVx(:,0) = Tdomain%sEdge(ne)%DumpVx(:,0) * Tdomain%sEdge(ne)%DumpVx(:,1)

        Tdomain%sEdge(ne)%DumpVy(:,1) = Tdomain%sEdge(ne)%MassMat + Tdomain%sEdge(ne)%DumpMass(:,1)
        Tdomain%sEdge(ne)%DumpVy(:,1) = 1./Tdomain%sEdge(ne)%DumpVy(:,1)
        Tdomain%sEdge(ne)%DumpVy(:,0) = Tdomain%sEdge(ne)%MassMat - Tdomain%sEdge(ne)%DumpMass(:,1)
        Tdomain%sEdge(ne)%DumpVy(:,0) = Tdomain%sEdge(ne)%DumpVy(:,0) * Tdomain%sEdge(ne)%DumpVy(:,1)

        Tdomain%sEdge(ne)%DumpVz(:,1) = Tdomain%sEdge(ne)%MassMat + Tdomain%sEdge(ne)%DumpMass(:,2)
        Tdomain%sEdge(ne)%DumpVz(:,1) = 1./Tdomain%sEdge(ne)%DumpVz(:,1)
        Tdomain%sEdge(ne)%DumpVz(:,0) = Tdomain%sEdge(ne)%MassMat - Tdomain%sEdge(ne)%DumpMass(:,2)
        Tdomain%sEdge(ne)%DumpVz(:,0) = Tdomain%sEdge(ne)%DumpVz(:,0) * Tdomain%sEdge(ne)%DumpVz(:,1)

        if (Tdomain%sEdge(ne)%FPML) then
              Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) *  Tdomain%sEdge(ne)%DumpVx(:,1)
              Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) *  Tdomain%sEdge(ne)%DumpVy(:,1)
              Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) *  Tdomain%sEdge(ne)%DumpVz(:,1)
       endif

        deallocate (Tdomain%sEdge(ne)%DumpMass) 
    endif
    Tdomain%sEdge(ne)%MassMat = 1./ Tdomain%sEdge(ne)%MassMat
enddo

do nv = 0,Tdomain%n_vertex-1
    if (Tdomain%sVertex(nv)%PML) then
        Tdomain%sVertex(nv)%DumpVx(1) = Tdomain%sVertex(nv)%MassMat + Tdomain%sVertex(nv)%DumpMass(0)
        Tdomain%sVertex(nv)%DumpVx(1) = 1./Tdomain%sVertex(nv)%DumpVx(1)
        Tdomain%sVertex(nv)%DumpVx(0) = Tdomain%sVertex(nv)%MassMat - Tdomain%sVertex(nv)%DumpMass(0)
        Tdomain%sVertex(nv)%DumpVx(0) = Tdomain%sVertex(nv)%DumpVx(0) * Tdomain%sVertex(nv)%DumpVx(1)

        Tdomain%sVertex(nv)%DumpVy(1) = Tdomain%sVertex(nv)%MassMat + Tdomain%sVertex(nv)%DumpMass(1)
        Tdomain%sVertex(nv)%DumpVy(1) = 1./Tdomain%sVertex(nv)%DumpVy(1)
        Tdomain%sVertex(nv)%DumpVy(0) = Tdomain%sVertex(nv)%MassMat - Tdomain%sVertex(nv)%DumpMass(1)
        Tdomain%sVertex(nv)%DumpVy(0) = Tdomain%sVertex(nv)%DumpVy(0) * Tdomain%sVertex(nv)%DumpVy(1)

        Tdomain%sVertex(nv)%DumpVz(1) = Tdomain%sVertex(nv)%MassMat + Tdomain%sVertex(nv)%DumpMass(2)
        Tdomain%sVertex(nv)%DumpVz(1) = 1./Tdomain%sVertex(nv)%DumpVz(1)
        Tdomain%sVertex(nv)%DumpVz(0) = Tdomain%sVertex(nv)%MassMat - Tdomain%sVertex(nv)%DumpMass(2)
        Tdomain%sVertex(nv)%DumpVz(0) = Tdomain%sVertex(nv)%DumpVz(0) * Tdomain%sVertex(nv)%DumpVz(1)
        if (Tdomain%sVertex(nv)%FPML) then
              Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) *  Tdomain%sVertex(nv)%DumpVx(1)
              Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) *  Tdomain%sVertex(nv)%DumpVy(1)
              Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) *  Tdomain%sVertex(nv)%DumpVz(1)
        endif
    endif
    Tdomain%sVertex(nv)%MassMat = 1./ Tdomain%sVertex(nv)%MassMat
enddo

! Mass Mat fo SO
if (Tdomain%logicD%super_object_local_present) then
 if (Tdomain%super_object_type == "P") then

   do nf = 0, Tdomain%sPlaneW%n_faces-1     
      nv_aus = Tdomain%sPlaneW%pFace(nf)%Face_UP
      Tdomain%sPlaneW%pFace(nf)%MassMat_Up = Tdomain%sFace(nv_aus)%MassMat
      nv_aus = Tdomain%sPlaneW%pFace(nf)%Face_DOWN
      Tdomain%sPlaneW%pFace(nf)%MassMat_Down =  Tdomain%sFace(nv_aus)%MassMat
   enddo
   do ne = 0, Tdomain%sPlaneW%n_edges-1  
      nv_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_UP
      Tdomain%sPlaneW%pEdge(ne)%MassMat_Up = Tdomain%sEdge(nv_aus)%MassMat
      nv_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_DOWN
      Tdomain%sPlaneW%pEdge(ne)%MassMat_Down =  Tdomain%sEdge(nv_aus)%MassMat
   enddo
   do nv = 0, Tdomain%sPlaneW%n_vertices-1
      nv_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP
      Tdomain%sPlaneW%pVertex(nv)%MassMat_Up = Tdomain%sVertex(nv_aus)%MassMat
      nv_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_DOWN
      Tdomain%sPlaneW%pVertex(nv)%MassMat_Down =  Tdomain%sVertex(nv_aus)%MassMat
   enddo

 endif
endif

return
end subroutine

