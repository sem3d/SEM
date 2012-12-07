module selement

    implicit none

    type :: element
        integer :: mat_index, ngllx, nglly, ngllz
        integer, dimension (:), pointer :: Control_nodes
        integer, dimension (0:5) :: Near_Faces, Orient_Faces
        integer, dimension (0:11) :: Near_Edges, Orient_Edges
        integer, dimension (0:7) :: Near_Vertices
        integer, dimension (:,:,:), pointer :: Iglobnum,Num
        real, dimension (:,:,:), pointer :: Jacob, Density, Lambda, Mu, MassMat
        real, dimension(:,:,:,:), pointer :: ACoeff, Forces,Veloc,Displ,Accel,V0, Diagonal_Stress, Residual_Stress
        real, dimension(:,:,:,:,:), pointer :: InvGrad
      ! fluid part
        real, dimension(:,:,:), pointer:: Phi,VelPhi0,VelPhi,AccelPhi
     ! PML allocation 
        logical :: PML, FPML
        real, dimension(:,:,:,:), pointer :: Diagonal_Stress1, Diagonal_Stress2, Diagonal_Stress3 
        real, dimension(:,:,:,:), pointer :: Residual_Stress1, Residual_Stress2
        real, dimension(:,:,:,:), pointer :: DumpSx,DumpSy,DumpSz
        real, dimension(:,:,:,:), pointer :: Forces1,Forces2,Forces3, Veloc1,Veloc2,Veloc3
        real, dimension(:,:,:,:), pointer :: DumpVx,DumpVy,DumpVz, DumpMass
        real, dimension(:,:,:,:), pointer :: I_Diagonal_Stress1, I_Diagonal_Stress2, I_Diagonal_Stress3
        real, dimension(:,:,:,:), pointer :: I_Residual_Stress1, I_Residual_Stress2
        real, dimension(:,:,:,:), pointer :: Iveloc1, Iveloc2, Iveloc3
        real, dimension(:,:,:), pointer :: Isx, Isy, Isz, Ivx, Ivy, Ivz
      ! fluid part
        real, dimension(:,:,:), pointer:: ForcesFl,ForcesFl1,ForcesFl2,ForcesFl3,VelPhi1,VelPhi2,VelPhi3
        
     ! solid-fluid
        logical  :: solid
    end type element


contains 

!--------------------------------------------------------------
subroutine Prediction_Elem_Veloc(Elem,alpha,bega,gam1,dt)

    implicit none
    type(Element), intent(inout) :: Elem
    real, intent(in) :: alpha,bega,gam1,dt
    integer :: ngllx, nglly, ngllz

   
ngllx = Elem%ngllx ; nglly = Elem%nglly ; ngllz = Elem%ngllz

Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = Elem%Displ + dt * Elem%Veloc +    &
                                                 dt**2*(0.5d0-bega)*Elem%Accel
Elem%V0 = Elem%Veloc
Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2) = alpha*Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2)+   &
                                                 (1d0-alpha) * Elem%Displ 

return
end subroutine Prediction_Elem_Veloc
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine Prediction_Elem_VelPhi(Elem,alpha,bega,gam1,dt)

    implicit none
    type(Element), intent(inout) :: Elem
    real, intent(in) :: alpha,bega,gam1,dt
    integer :: ngllx, nglly, ngllz

   
ngllx = Elem%ngllx ; nglly = Elem%nglly ; ngllz = Elem%ngllz

Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2) = Elem%Phi+dt*Elem%VelPhi+dt**2*(0.5d0-bega)*Elem%AccelPhi
Elem%VelPhi0 = Elem%VelPhi
Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2) = alpha*Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2) +   &
                                               (1d0-alpha) * Elem%Phi

return
end subroutine Prediction_Elem_VelPhi
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine Correction_Elem_Veloc(Elem,bega,gam1,dt)

    implicit none
    type(Element), intent(inout) :: Elem
    real, intent(in) :: bega, gam1, dt
    integer :: ngllx, nglly, ngllz, i


ngllx = Elem%ngllx ;  nglly = Elem%nglly ; ngllz = Elem%ngllz
do i = 0,2
    Elem%Forces(1:ngllx-2,1:nglly-2, 1:ngllz-2,i) = Elem%MassMat(:,:,:)*    &
                                  Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
enddo
Elem%Veloc(:,:,:,:)  = Elem%v0(:,:,:,:)+ dt * Elem%Forces(1:ngllx-2,1:nglly-2,1:ngllz-2,:)
Elem%Accel  = Elem%Accel + gam1 /dt * (Elem%Veloc-Elem%V0)
Elem%Displ = Elem%Displ + bega *dt * (Elem%Veloc+Elem%V0)

return
end subroutine Correction_Elem_Veloc
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine Correction_Elem_VelPhi(Elem,bega,gam1,dt)

    implicit none
    type(Element), intent(inout) :: Elem
    real, intent(in) :: bega, gam1, dt
    integer :: ngllx, nglly, ngllz, i


ngllx = Elem%ngllx ;  nglly = Elem%nglly ; ngllz = Elem%ngllz

Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2) = Elem%MassMat(:,:,:)*    &
                                  Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2)
Elem%VelPhi(:,:,:) = Elem%VelPhi0(:,:,:)+ dt * Elem%ForcesFl(1:ngllx-2,1:nglly-2,1:ngllz-2)
Elem%AccelPhi  = Elem%AccelPhi + gam1/dt*(Elem%VelPhi-Elem%VelPhi0)
Elem%Phi = Elem%Phi + bega*dt*(Elem%VelPhi+Elem%VelPhi0)

return
end subroutine Correction_Elem_VelPhi
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine compute_InternalForces_Elem(Elem,hprimex,htprimex,hprimey,htprimey,hprimez,htprimez)

    implicit none

    type(Element), intent(inout) :: Elem
    real, dimension(0:Elem%ngllx-1,0:Elem%ngllx-1), intent(in) :: hprimex,htprimex
    real, dimension(0:Elem%nglly-1,0:Elem%nglly-1), intent(in) :: hprimey,htprimey
    real, dimension(0:Elem%ngllz-1,0:Elem%ngllz-1), intent(in) :: hprimez,hTprimez

    integer :: n_z, m1, m2, m3
    real, dimension(0:Elem%ngllx-1,0:Elem%nglly-1) :: s0z
    real, dimension(0:Elem%nglly-1,0:Elem%ngllz-1) :: s0x
    real, dimension(0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1) :: dUx_dxi, dUx_deta, dUx_dzeta, &
                                                                     dUy_dxi, dUy_deta, dUy_dzeta, &
                                                                     dUz_dxi, dUz_deta, dUz_dzeta, &
                                                                     t1, s0, Uxloc, Uyloc, Uzloc

m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz


!- gradient at GLL points
! dUx_(dxi,deta,dzeta)
call DGEMM('N','N',m1,m2*m3,m1,1.,htprimex,m1,Elem%Forces(0,0,0,0),m1,0.,dUx_dxi,m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM('N','N',m1,m2,m2,1.,Elem%Forces(0,0,n_z,0),m1,hprimey,m2,0.,dUx_deta(0,0,n_z),m1)
enddo
call DGEMM('N','N',m1*m2,m3,m3,1.,Elem%Forces(0,0,0,0),m1*m2,hprimez,m3,0.,dUx_dzeta,m1*m2)
! dUy_(dxi,deta,dzeta)
call DGEMM('N','N',m1,m2*m3,m1,1.,htprimex,m1,Elem%Forces(0,0,0,1),m1,0.,dUy_dxi,m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM('N','N',m1,m2,m2,1.,Elem%Forces(0,0,n_z,1),m1,hprimey,m2,0.,dUy_deta(0,0,n_z),m1)
enddo
call DGEMM('N','N',m1*m2,m3,m3,1.,Elem%Forces(0,0,0,1),m1*m2,hprimez,m3,0.,dUy_dzeta,m1*m2)
! dUz_(dxi,deta,dzeta)
call DGEMM('N','N',m1,m2*m3,m1,1.,htprimex,m1,Elem%Forces(0,0,0,2),m1,0.,dUz_dxi,m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM('N','N',m1,m2,m2,1.,Elem%Forces(0,0,n_z,2),m1,hprimey,m2,0.,dUz_deta(0,0,n_z),m1)
enddo
call DGEMM('N','N',m1*m2,m3,m3,1.,Elem%Forces(0,0,0,2),m1*m2,hprimez,m3,0.,dUz_dzeta,m1*m2)


!- Internal forces
t1 = Elem%Acoeff(:,:,:,0)*dUx_dxi + Elem%Acoeff(:,:,:,1)*dUx_deta + Elem%Acoeff(:,:,:,2)*dUx_dzeta +  &
     Elem%Acoeff(:,:,:,3)*dUy_dxi + Elem%Acoeff(:,:,:,4)*dUy_deta + Elem%Acoeff(:,:,:,5)*dUy_dzeta +  &
     Elem%Acoeff(:,:,:,6)*dUz_dxi + Elem%Acoeff(:,:,:,7)*dUz_deta + Elem%Acoeff(:,:,:,8)*dUz_dzeta 

call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Uxloc,m1)

t1 = Elem%Acoeff(:,:,:,1)*dUx_dxi + Elem%Acoeff(:,:,:,9)*dUx_deta + Elem%Acoeff(:,:,:,10)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,11)*dUy_dxi + Elem%Acoeff(:,:,:,12)*dUy_deta + Elem%Acoeff(:,:,:,13)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,14)*dUz_dxi + Elem%Acoeff(:,:,:,15)*dUz_deta + Elem%Acoeff(:,:,:,16)*dUz_dzeta

do n_z = 0,Elem%ngllz-1
     call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
enddo
Uxloc = s0 + Uxloc

t1 = Elem%Acoeff(:,:,:,2)*dUx_dxi + Elem%Acoeff(:,:,:,10)*dUx_deta + Elem%Acoeff(:,:,:,17)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,18)*dUy_dxi + Elem%Acoeff(:,:,:,19)*dUy_deta + Elem%Acoeff(:,:,:,20)*dUy_dzeta +&
     Elem%Acoeff(:,:,:,21)*dUz_dxi + Elem%ACoeff(:,:,:,22)*dUz_deta + Elem%Acoeff(:,:,:,23)*dUz_dzeta

call DGEMM('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
Uxloc = s0 + Uxloc

t1 = Elem%Acoeff(:,:,:,3)*dUx_dxi + Elem%Acoeff(:,:,:,11)*dUx_deta + Elem%Acoeff(:,:,:,18)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,24)*dUy_dxi + Elem%Acoeff(:,:,:,25)*dUy_deta + Elem%Acoeff(:,:,:,26)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,27)*dUz_dxi + Elem%Acoeff(:,:,:,28)*dUz_deta + Elem%Acoeff(:,:,:,29)*dUz_dzeta

call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Uyloc,m1)

t1 = Elem%Acoeff(:,:,:,4)*dUx_dxi + Elem%Acoeff(:,:,:,12)*dUx_deta + Elem%Acoeff(:,:,:,19)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,25)*dUy_dxi + Elem%Acoeff(:,:,:,30)*dUy_deta + Elem%Acoeff(:,:,:,31)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,32)*dUz_dxi + Elem%Acoeff(:,:,:,33)*dUz_deta + Elem%Acoeff(:,:,:,34)*dUz_dzeta 

do n_z = 0,Elem%ngllz-1
     call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
enddo
Uyloc = s0 + Uyloc

t1 = Elem%Acoeff(:,:,:,5)*dUx_dxi + Elem%Acoeff(:,:,:,13)*dUx_deta + Elem%Acoeff(:,:,:,20)* dUx_dzeta +&
     Elem%Acoeff(:,:,:,26)*dUy_dxi + Elem%Acoeff(:,:,:,31)*dUy_deta + Elem%Acoeff(:,:,:,35)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,36)*dUz_dxi + Elem%Acoeff(:,:,:,37)*dUz_deta + Elem%Acoeff(:,:,:,38)*dUz_dzeta

call DGEMM('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
Uyloc = s0 + Uyloc

t1 = Elem%Acoeff(:,:,:,6)*dUx_dxi + Elem%Acoeff(:,:,:,14)*dUx_deta + Elem%Acoeff(:,:,:,21)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,27)*dUy_dxi + Elem%Acoeff(:,:,:,32)*dUy_deta + Elem%Acoeff(:,:,:,36)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,39)*dUz_dxi + Elem%Acoeff(:,:,:,40)*dUz_deta + Elem%Acoeff(:,:,:,41)*dUz_dzeta

call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Uzloc,m1)

t1 = Elem%Acoeff(:,:,:,7)*dUx_dxi + Elem%Acoeff(:,:,:,15)*dUx_deta + Elem%Acoeff(:,:,:,22)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,28)*dUy_dxi + Elem%Acoeff(:,:,:,33)*dUy_deta + Elem%Acoeff(:,:,:,37)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,40)*dUz_dxi + Elem%Acoeff(:,:,:,42)*dUz_deta + Elem%Acoeff(:,:,:,43)*dUz_dzeta
      
do n_z = 0,Elem%ngllz-1
     call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
enddo
Uzloc = s0 + Uzloc
     
t1 = Elem%Acoeff(:,:,:,8)*dUx_dxi + Elem%Acoeff(:,:,:,16)*dUx_deta + Elem%Acoeff(:,:,:,23)*dUx_dzeta + &
     Elem%Acoeff(:,:,:,29)*dUy_dxi + Elem%Acoeff(:,:,:,34)*dUy_deta + Elem%Acoeff(:,:,:,38)*dUy_dzeta + &
     Elem%Acoeff(:,:,:,41)*dUz_dxi + Elem%Acoeff(:,:,:,43)*dUz_deta + Elem%Acoeff(:,:,:,44)*dUz_dzeta        

call DGEMM ('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
Uzloc = Uzloc + s0

Elem%Forces(:,:,:,0) = Uxloc
Elem%Forces(:,:,:,1) = Uyloc
Elem%Forces(:,:,:,2) = Uzloc

return
end subroutine compute_InternalForces_Elem
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine compute_InternalForces_Elem_Fluid(Elem,hprimex,htprimex,hprimey,htprimey,hprimez,htprimez)

    implicit none

    type(Element), intent(inout) :: Elem
    real, dimension(0:Elem%ngllx-1,0:Elem%ngllx-1), intent(in) :: hprimex,htprimex
    real, dimension(0:Elem%nglly-1,0:Elem%nglly-1), intent(in) :: hprimey,htprimey
    real, dimension(0:Elem%ngllz-1,0:Elem%ngllz-1), intent(in) :: hprimez,hTprimez
    integer :: n_z,m1,m2,m3
    real, dimension(0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1) ::    &
                                      dPhi_dxi,dPhi_deta,dPhi_dzeta,t1,s0,Floc

m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

!- Modification: potential -> density*potentiel
Elem%ForcesFl(:,:,:) = Elem%density(:,:,:)*Elem%ForcesFl(:,:,:)

!- gradients at GLLs points
 ! d(rho*Phi)_dxi
   call DGEMM('N','N',m1,m2*m3,m1,1d0,htprimex,m1,Elem%ForcesFl(0,0,0),m1,0d0,dPhi_dxi,m1)
 ! d(rho*Phi)_deta
   do n_z = 0,Elem%ngllz-1
     call DGEMM('N','N',m1,m2,m2,1d0,Elem%ForcesFl(0,0,n_z),m1,hprimey,m2,0d0,dPhi_deta(0,0,n_z),m1)
   enddo
 ! d(rho*Phi)_dzeta
   call DGEMM('N','N',m1*m2,m3,m3,1d0,Elem%ForcesFl(0,0,0),m1*m2,hprimez,m3,0d0,dPhi_dzeta,m1*m2)


!- Internal Forces
t1 = Elem%Acoeff(:,:,:,0)*dPhi_dxi + Elem%Acoeff(:,:,:,1)*dPhi_deta + Elem%Acoeff(:,:,:,2)*dPhi_dzeta
call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,t1(0,0,0),m1,0.,Floc,m1)

t1 = Elem%Acoeff(:,:,:,1)*dPhi_dxi + Elem%Acoeff(:,:,:,3)*dPhi_deta + Elem%Acoeff(:,:,:,4)*dPhi_dzeta
do n_z = 0,Elem%ngllz-1
     call DGEMM('N','N',m1,m2,m2,1.,t1(0,0,n_z),m1,htprimey,m2,0.,s0(0,0,n_z),m1)
enddo
Floc = s0 + Floc

t1 = Elem%Acoeff(:,:,:,2)*dPhi_dxi + Elem%Acoeff(:,:,:,4)*dPhi_deta + Elem%Acoeff(:,:,:,5)*dPhi_dzeta
call DGEMM('N','N',m1*m2,m3,m3,1.,t1(0,0,0),m1*m2,htprimez,m3,0.,s0,m1*m2)
Floc = s0 + Floc

!-
Elem%ForcesFl(:,:,:) = Floc

return
end subroutine compute_InternalForces_Elem_Fluid
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine Prediction_Elem_PML_Veloc(Elem,bega,dt,hTprimex,Hprimey,Hprimez,rg,n)

implicit none

type (Element), intent (INOUT) :: Elem
real, intent (IN) :: bega, dt
real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprimex
real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
integer, intent (IN) :: rg, n

real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1) :: s0z
real, dimension (0:Elem%nglly-1, 0:Elem%ngllz-1) :: s0x
real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: Vxloc,Vyloc,Vzloc, dVx_dxi,dVx_deta,dVx_dzeta, & 
                                                                    dVy_dxi,dVy_deta,dVy_dzeta, dVz_dxi,dVz_deta,dVz_dzeta

integer :: m1, m2,m3 ,n_z,n_x,i,j,k


m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz

Elem%Forces(1:m1-2,1:m2-2, 1:m3-2, 0:2)  = Elem%Veloc(:,:,:,:) + dt *(0.5-bega) *Elem%Accel(:,:,:,:)


call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,0) ,m1, 0., dVx_dxi, m1 )
do n_z = 0,Elem%ngllz-1
   call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,0), m1, hprimey ,m2, 0., dVx_deta(0,0,n_z),m1 )
enddo
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,0), m1*m2, hprimez ,m3, 0., dVx_dzeta, m1*m2 )

call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,1) ,m1, 0., dVy_dxi, m1 )
do n_z = 0,Elem%ngllz-1
   call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,1), m1, hprimey ,m2, 0., dVy_deta(0,0,n_z),m1 )
enddo
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,1), m1*m2, hprimez ,m3, 0., dVy_dzeta, m1*m2 )

call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,2) ,m1, 0., dVz_dxi, m1 )
do n_z = 0,Elem%ngllz-1
   call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,2), m1, hprimey ,m2, 0., dVz_deta(0,0,n_z),m1)
enddo
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,2), m1*m2, hprimez ,m3, 0., dVz_dzeta, m1*m2 )


Elem%Diagonal_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,0) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta + Elem%Acoeff(:,:,:,2) * dVx_dzeta)
Elem%Diagonal_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,0) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta)
Elem%Diagonal_Stress3 (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,0) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta)

Elem%Diagonal_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,1) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta)
Elem%Diagonal_Stress2 (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,1) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta + Elem%Acoeff(:,:,:,14) * dVy_dzeta)
Elem%Diagonal_Stress3 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,1) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta)

Elem%Diagonal_Stress1 (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,2) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta)
Elem%Diagonal_Stress2 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,2) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta)
Elem%Diagonal_Stress3 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,2) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,15) * dVz_dxi + Elem%Acoeff(:,:,:,16) * dVz_deta + Elem%Acoeff(:,:,:,17) * dVz_dzeta)

Elem%Diagonal_Stress = Elem%Diagonal_Stress1 + Elem%Diagonal_Stress2 + Elem%Diagonal_Stress3

Elem%Residual_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,0) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVy_dxi + Elem%Acoeff(:,:,:,19) * dVy_deta + Elem%Acoeff(:,:,:,20) * dVy_dzeta)
Elem%Residual_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,0) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVx_dxi + Elem%Acoeff(:,:,:,22) * dVx_deta + Elem%Acoeff(:,:,:,23) * dVx_dzeta)

Elem%Residual_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,1) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta + Elem%Acoeff(:,:,:,20) * dVz_dzeta)
Elem%Residual_Stress2 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,1) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVx_dxi + Elem%Acoeff(:,:,:,25) * dVx_deta + Elem%Acoeff(:,:,:,26) * dVx_dzeta)

Elem%Residual_Stress1 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,2) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVz_dxi + Elem%Acoeff(:,:,:,22) * dVz_deta + Elem%Acoeff(:,:,:,23) * dVz_dzeta)
Elem%Residual_Stress2 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,2) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVy_dxi + Elem%Acoeff(:,:,:,25) * dVy_deta + Elem%Acoeff(:,:,:,26) * dVy_dzeta)

Elem%Residual_Stress = Elem%Residual_Stress1 + Elem%Residual_Stress2 

return
end subroutine Prediction_Elem_PML_Veloc
!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine Prediction_Elem_PML_VelPhi(Elem,bega,dt,hTprimex,Hprimey,Hprimez)
  ! same as previously, but for fluid part
    implicit none

    type(Element), intent(inout) :: Elem
    real, intent(in) :: bega, dt
    real, dimension(0:Elem%ngllx-1,0:Elem%ngllx-1), intent(in) :: hTprimex
    real, dimension(0:Elem%nglly-1,0:Elem%nglly-1), intent(in) :: hprimey
    real, dimension(0:Elem%ngllz-1,0:Elem%ngllz-1), intent(in) :: hprimez

    real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVelPhi_dxi,dVelPhi_deta,dVelPhi_dzeta 
    integer :: m1,m2,m3,n_z,n_x,i,j,k


m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

! prediction in the element
Elem%ForcesFl(1:m1-2,1:m2-2,1:m3-2) = Elem%VelPhi(:,:,:) +dt*(0.5-bega)*Elem%AccelPhi(:,:,:)
! potential -> -pressure
Elem%ForcesFl(:,:,:) = Elem%Density(:,:,:)*Elem%VelPhi(:,:,:)

! d(rho*Phi)_d(xi,eta,zeta)
call DGEMM('N','N',m1,m2*m3,m1,1.,htprimex,m1,Elem%ForcesFl(0,0,0),m1,0.,dVelPhi_dxi,m1)
do n_z = 0,Elem%ngllz-1
   call DGEMM('N','N',m1,m2,m2,1.,Elem%ForcesFl(0,0,n_z),m1,hprimey,m2,0.,dVelPhi_deta(0,0,n_z),m1)
enddo
call DGEMM('N','N',m1*m2,m3,m3,1.,Elem%ForcesFl(0,0,0),m1*m2,hprimez,m3,0.,dVelPhi_dzeta,m1*m2)

! prediction for (physical) velocity (which is the equivalent of a stress, here)
! V_x^x
Elem%Veloc1(:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Veloc1(:,:,:,0) + Elem%DumpSx(:,:,:,1) * Dt *  &
        (Elem%Acoeff(:,:,:,0) * dVelPhi_dxi + Elem%Acoeff(:,:,:,1) * dVelPhi_deta + Elem%Acoeff(:,:,:,2) * dVelPhi_dzeta)
! V_x^y
Elem%Veloc2(:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Veloc2(:,:,:,0)
! V_x^z
Elem%Veloc3(:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Veloc3(:,:,:,0)
! V_y^x
Elem%Veloc1(:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Veloc1(:,:,:,1)
! V_y^y
Elem%Veloc2(:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Veloc2(:,:,:,1) + Elem%DumpSy(:,:,:,1) * Dt *  &
        (Elem%Acoeff(:,:,:,3) * dVelPhi_dxi + Elem%Acoeff(:,:,:,4) * dVelPhi_deta + Elem%Acoeff(:,:,:,5) * dVelPhi_dzeta)
! V_y^z
Elem%Veloc3(:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Veloc3(:,:,:,1)
! V_z^x
Elem%Veloc1(:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Veloc1(:,:,:,2)
! V_z^y
Elem%Veloc2(:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Veloc2(:,:,:,2)
! V_z^z
Elem%Veloc3(:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Veloc3(:,:,:,2) + Elem%DumpSz(:,:,:,1) * Dt *  &
        (Elem%Acoeff(:,:,:,6) * dVelPhi_dxi + Elem%Acoeff(:,:,:,7) * dVelPhi_deta + Elem%Acoeff(:,:,:,8) * dVelPhi_dzeta)

! total velocity vector after dumping = sum of splitted parts
Elem%Veloc(:,:,:,:) = Elem%Veloc1(:,:,:,:) + Elem%Veloc2(:,:,:,:) + Elem%Veloc3(:,:,:,:)

return
end subroutine Prediction_Elem_PML_VelPhi
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
subroutine Prediction_Elem_FPML_Veloc(Elem,bega,dt,hTprimex,Hprimey,Hprimez,rg,n,fil)

implicit none

type (Element), intent (INOUT) :: Elem
real, intent (IN) :: bega, dt, fil
real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprimex
real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey
real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
integer, intent (IN) :: rg, n

real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1) :: s0z
real, dimension (0:Elem%nglly-1, 0:Elem%ngllz-1) :: s0x
real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: Vxloc,Vyloc,Vzloc, dVx_dxi,dVx_deta,dVx_dzeta, & 
                                                                    dVy_dxi,dVy_deta,dVy_dzeta, dVz_dxi,dVz_deta,dVz_dzeta, Stress_ausiliar

integer :: m1, m2,m3 ,n_z,n_x,i,j,k
real :: fil2


m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz
fil2 = fil**2

Elem%Forces(1:m1-2,1:m2-2, 1:m3-2, 0:2)  = Elem%Veloc(:,:,:,:) + dt *(0.5-bega) *Elem%Accel(:,:,:,:)


call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,0) ,m1, 0., dVx_dxi, m1 )
do n_z = 0,Elem%ngllz-1
   call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,0), m1, hprimey ,m2, 0., dVx_deta(0,0,n_z),m1 )
enddo
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,0), m1*m2, hprimez ,m3, 0., dVx_dzeta, m1*m2 )

call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,1) ,m1, 0., dVy_dxi, m1 )
do n_z = 0,Elem%ngllz-1
   call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,1), m1, hprimey ,m2, 0., dVy_deta(0,0,n_z),m1 )
enddo
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,1), m1*m2, hprimez ,m3, 0., dVy_dzeta, m1*m2 )

call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., htprimex, m1,Elem%Forces(0,0,0,2) ,m1, 0., dVz_dxi, m1 )
do n_z = 0,Elem%ngllz-1
   call DGEMM ( 'N', 'N', m1, m2, m2, 1., Elem%Forces(0,0,n_z,2), m1, hprimey ,m2, 0., dVz_deta(0,0,n_z),m1)
enddo
call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., Elem%Forces(0,0,0,2), m1*m2, hprimez ,m3, 0., dVz_dzeta, m1*m2 )


Stress_Ausiliar = Elem%Diagonal_Stress1 (:,:,:,0)
Elem%Diagonal_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,0) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,0) * dVx_dxi + Elem%Acoeff(:,:,:,1) * dVx_deta + &
                              Elem%Acoeff(:,:,:,2) * dVx_dzeta) + Elem%Isx * Elem%I_Diagonal_stress1 (:,:,:,0)
Elem%I_Diagonal_Stress1  (:,:,:,0)= Fil2* Elem%I_Diagonal_Stress1(:,:,:,0) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Diagonal_Stress1 (:,:,:,0) )                  

Stress_Ausiliar = Elem%Diagonal_Stress2 (:,:,:,0)
Elem%Diagonal_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,0) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + & 
                           Elem%Acoeff(:,:,:,5) * dVy_dzeta) + Elem%Isy * Elem%I_Diagonal_stress2 (:,:,:,0)
Elem%I_Diagonal_Stress2  (:,:,:,0)= Fil2* Elem%I_Diagonal_Stress2(:,:,:,0) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Diagonal_Stress2 (:,:,:,0) )                  

Stress_Ausiliar = Elem%Diagonal_Stress3 (:,:,:,0)
Elem%Diagonal_Stress3 (:,:,:,0) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,0) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + &
                              Elem%Acoeff(:,:,:,8) * dVz_dzeta) + Elem%Isz * Elem%I_Diagonal_stress3 (:,:,:,0)
Elem%I_Diagonal_Stress3  (:,:,:,0)= Fil2* Elem%I_Diagonal_Stress3(:,:,:,0) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Diagonal_Stress3 (:,:,:,0) )    

Stress_Ausiliar = Elem%Diagonal_Stress1 (:,:,:,1)
Elem%Diagonal_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,1) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta) + Elem%Isx * Elem%I_Diagonal_stress1 (:,:,:,1)
Elem%I_Diagonal_Stress1  (:,:,:,1)= Fil2* Elem%I_Diagonal_Stress1(:,:,:,1) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Diagonal_Stress1 (:,:,:,1) )                  

Stress_Ausiliar = Elem%Diagonal_Stress2 (:,:,:,1)
Elem%Diagonal_Stress2 (:,:,:,1) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,1) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,12) * dVy_dxi + Elem%Acoeff(:,:,:,13) * dVy_deta + Elem%Acoeff(:,:,:,14) * dVy_dzeta) + Elem%Isy * Elem%I_Diagonal_stress2 (:,:,:,1)
Elem%I_Diagonal_Stress2  (:,:,:,1)= Fil2* Elem%I_Diagonal_Stress2(:,:,:,1) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Diagonal_Stress2 (:,:,:,1) )                  

Stress_Ausiliar = Elem%Diagonal_Stress3 (:,:,:,1)
Elem%Diagonal_Stress3 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,1) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,6) * dVz_dxi + Elem%Acoeff(:,:,:,7) * dVz_deta + Elem%Acoeff(:,:,:,8) * dVz_dzeta) + Elem%Isz * Elem%I_Diagonal_stress3 (:,:,:,1)
Elem%I_Diagonal_Stress3  (:,:,:,1)= Fil2* Elem%I_Diagonal_Stress3(:,:,:,1) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Diagonal_Stress3 (:,:,:,1) )   

Stress_Ausiliar = Elem%Diagonal_Stress1 (:,:,:,2)
Elem%Diagonal_Stress1 (:,:,:,2) = Elem%DumpSx(:,:,:,0) * Elem%Diagonal_Stress1 (:,:,:,2) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,9) * dVx_dxi + Elem%Acoeff(:,:,:,10) * dVx_deta + Elem%Acoeff(:,:,:,11) * dVx_dzeta) + Elem%Isx * Elem%I_Diagonal_stress1 (:,:,:,2)
Elem%I_Diagonal_Stress1  (:,:,:,2)= Fil2* Elem%I_Diagonal_Stress1(:,:,:,2) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Diagonal_Stress1 (:,:,:,2) )                  

Stress_Ausiliar = Elem%Diagonal_Stress2 (:,:,:,2)
Elem%Diagonal_Stress2 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Diagonal_Stress2 (:,:,:,2) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,3) * dVy_dxi + Elem%Acoeff(:,:,:,4) * dVy_deta + Elem%Acoeff(:,:,:,5) * dVy_dzeta)+ Elem%Isy * Elem%I_Diagonal_stress2 (:,:,:,2)
Elem%I_Diagonal_Stress2  (:,:,:,2)= Fil2* Elem%I_Diagonal_Stress2(:,:,:,2) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Diagonal_Stress2 (:,:,:,2) )

Stress_Ausiliar = Elem%Diagonal_Stress3 (:,:,:,2)
Elem%Diagonal_Stress3 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Diagonal_Stress3 (:,:,:,2) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,15) * dVz_dxi + Elem%Acoeff(:,:,:,16) * dVz_deta + Elem%Acoeff(:,:,:,17) * dVz_dzeta) + Elem%Isz * Elem%I_Diagonal_stress3 (:,:,:,2)
Elem%I_Diagonal_Stress3  (:,:,:,2)= Fil2* Elem%I_Diagonal_Stress3(:,:,:,2) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Diagonal_Stress3 (:,:,:,2) )   


Elem%Diagonal_Stress = Elem%Diagonal_Stress1 + Elem%Diagonal_Stress2 + Elem%Diagonal_Stress3


Stress_Ausiliar = Elem%Residual_Stress1 (:,:,:,0)
Elem%Residual_Stress1 (:,:,:,0) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,0) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVy_dxi + Elem%Acoeff(:,:,:,19) * dVy_deta + Elem%Acoeff(:,:,:,20) * dVy_dzeta) + Elem%Isx * Elem%I_Residual_stress1(:,:,:,0)
Elem%I_Residual_Stress1  (:,:,:,0)= Fil2* Elem%I_Residual_Stress1(:,:,:,0) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Residual_Stress1 (:,:,:,0) )             

Stress_Ausiliar = Elem%Residual_Stress2 (:,:,:,0)
Elem%Residual_Stress2 (:,:,:,0) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,0) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVx_dxi + Elem%Acoeff(:,:,:,22) * dVx_deta + Elem%Acoeff(:,:,:,23) * dVx_dzeta) + Elem%Isy * Elem%I_Residual_stress2(:,:,:,0)
Elem%I_Residual_Stress2  (:,:,:,0) = Fil2* Elem%I_Residual_Stress2(:,:,:,0) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Residual_Stress2 (:,:,:,0) ) 

Stress_Ausiliar = Elem%Residual_Stress1 (:,:,:,1)
Elem%Residual_Stress1 (:,:,:,1) = Elem%DumpSx(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,1) + &
                    Elem%DumpSx(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,18) * dVz_dxi + Elem%Acoeff(:,:,:,19) * dVz_deta + Elem%Acoeff(:,:,:,20) * dVz_dzeta) + Elem%Isx * Elem%I_Residual_stress1(:,:,:,1)
Elem%I_Residual_Stress1  (:,:,:,1)= Fil2* Elem%I_Residual_Stress1(:,:,:,1) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Residual_Stress1 (:,:,:,1) )   

Stress_Ausiliar = Elem%Residual_Stress2 (:,:,:,1)
Elem%Residual_Stress2 (:,:,:,1) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,1) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVx_dxi + Elem%Acoeff(:,:,:,25) * dVx_deta + Elem%Acoeff(:,:,:,26) * dVx_dzeta) + Elem%Isz * Elem%I_Residual_stress2(:,:,:,1)
Elem%I_Residual_Stress2  (:,:,:,1) = Fil2* Elem%I_Residual_Stress2(:,:,:,1) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Residual_Stress2 (:,:,:,1) ) 

Stress_Ausiliar = Elem%Residual_Stress1 (:,:,:,2)
Elem%Residual_Stress1 (:,:,:,2) = Elem%DumpSy(:,:,:,0) * Elem%Residual_Stress1 (:,:,:,2) + &
                    Elem%DumpSy(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,21) * dVz_dxi + Elem%Acoeff(:,:,:,22) * dVz_deta + Elem%Acoeff(:,:,:,23) * dVz_dzeta) + Elem%Isy * Elem%I_Residual_stress1(:,:,:,2)
Elem%I_Residual_Stress1  (:,:,:,2)= Fil2* Elem%I_Residual_Stress1(:,:,:,2) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Residual_Stress1 (:,:,:,2) ) 

Stress_Ausiliar = Elem%Residual_Stress2 (:,:,:,2)
Elem%Residual_Stress2 (:,:,:,2) = Elem%DumpSz(:,:,:,0) * Elem%Residual_Stress2 (:,:,:,2) + &
                    Elem%DumpSz(:,:,:,1) * Dt * (Elem%Acoeff(:,:,:,24) * dVy_dxi + Elem%Acoeff(:,:,:,25) * dVy_deta + Elem%Acoeff(:,:,:,26) * dVy_dzeta)+ Elem%Isz * Elem%I_Residual_stress2(:,:,:,2)
Elem%I_Residual_Stress2  (:,:,:,2) = Fil2* Elem%I_Residual_Stress2(:,:,:,2) + 0.5 * (1.-Fil2) * & 
              (Stress_Ausiliar +Elem%Residual_Stress2 (:,:,:,2) ) 

Elem%Residual_Stress = Elem%Residual_Stress1 + Elem%Residual_Stress2 

return
end subroutine Prediction_Elem_FPML_Veloc
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine Correction_Elem_PML_Veloc(Elem,dt)

    implicit none

    type(Element), intent(inout) :: Elem
    real, intent(in) :: dt
    integer :: ngllx, nglly, ngllz, i 


ngllx = Elem%ngllx ; nglly = Elem%nglly ; ngllz = Elem%ngllz

do i = 0,2
    Elem%Veloc1(:,:,:,i) = Elem%DumpVx(:,:,:,0) * Elem%Veloc1(:,:,:,i) +    &
             dt * Elem%DumpVx(:,:,:,1)*Elem%Forces1(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
    Elem%Veloc2(:,:,:,i) = Elem%DumpVy(:,:,:,0) * Elem%Veloc2(:,:,:,i) +    &
             dt * Elem%DumpVy(:,:,:,1)*Elem%Forces2(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
    Elem%Veloc3(:,:,:,i) = Elem%DumpVz(:,:,:,0) * Elem%Veloc3(:,:,:,i) +    &
             dt * Elem%DumpVz(:,:,:,1)*Elem%Forces3(1:ngllx-2,1:nglly-2,1:ngllz-2,i)
enddo

Elem%Veloc = Elem%Veloc1 + Elem%Veloc2 + Elem%Veloc3

return
end subroutine Correction_Elem_PML_Veloc
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine Correction_Elem_PML_VelPhi(Elem,dt)

    implicit none

    type(Element), intent(inout) :: Elem
    real, intent(in) :: dt
    integer :: ngllx, nglly, ngllz, i 


ngllx = Elem%ngllx ; nglly = Elem%nglly ; ngllz = Elem%ngllz

Elem%VelPhi1(:,:,:) = Elem%DumpVx(:,:,:,0) * Elem%VelPhi1(:,:,:) +    &
          dt * Elem%DumpVx(:,:,:,1)*Elem%ForcesFl1(1:ngllx-2,1:nglly-2,1:ngllz-2)
Elem%VelPhi2(:,:,:) = Elem%DumpVy(:,:,:,0) * Elem%VelPhi2(:,:,:) +    &
          dt * Elem%DumpVy(:,:,:,1)*Elem%ForcesFl2(1:ngllx-2,1:nglly-2,1:ngllz-2)
Elem%VelPhi3(:,:,:) = Elem%DumpVz(:,:,:,0) * Elem%VelPhi3(:,:,:) +    &
          dt * Elem%DumpVz(:,:,:,1)*Elem%ForcesFl3(1:ngllx-2,1:nglly-2,1:ngllz-2)

Elem%VelPhi = Elem%VelPhi1 + Elem%VelPhi2 + Elem%VelPhi3

return
end subroutine Correction_Elem_PML_VelPhi
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine Correction_Elem_FPML_Veloc (Elem, dt, fil)

implicit none

type (Element), intent (INOUT) :: Elem
real, intent (IN) :: dt, fil

integer :: ngllx, nglly, ngllz, i 
real :: fil2
real, dimension (1:Elem%ngllx-2,1:Elem%nglly-2,1:Elem%ngllz-2) :: Ausiliar_velocity

ngllx = Elem%ngllx; nglly = Elem%nglly; ngllz = Elem%ngllz
fil2 = fil**2

do i = 0,2
   Ausiliar_Velocity =  Elem%Veloc1(:,:,:,i)
    Elem%Veloc1(:,:,:,i) = Elem%DumpVx(:,:,:,0) * Elem%Veloc1(:,:,:,i) + dt * Elem%DumpVx(:,:,:,1)*Elem%Forces1(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%Ivx * Elem%Iveloc1 (:,:,:,i) 
    Elem%Iveloc1 (:,:,:,i)  = Fil2*Elem%Iveloc1 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%Veloc1(:,:,:,i))

    Ausiliar_Velocity =  Elem%Veloc2(:,:,:,i)
    Elem%Veloc2(:,:,:,i) = Elem%DumpVy(:,:,:,0) * Elem%Veloc2(:,:,:,i) + dt * Elem%DumpVy(:,:,:,1)*Elem%Forces2(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%Ivy * Elem%Iveloc2 (:,:,:,i) 
    Elem%Iveloc2 (:,:,:,i)  = Fil2*Elem%Iveloc2 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%Veloc2(:,:,:,i))

    Ausiliar_Velocity =  Elem%Veloc3(:,:,:,i)
    Elem%Veloc3(:,:,:,i) = Elem%DumpVz(:,:,:,0) * Elem%Veloc3(:,:,:,i) + dt * Elem%DumpVz(:,:,:,1)*Elem%Forces3(1:ngllx-2,1:nglly-2,1:ngllz-2,i) + Elem%Ivz * Elem%Iveloc3 (:,:,:,i) 
    Elem%Iveloc3 (:,:,:,i)  = Fil2*Elem%Iveloc3 (:,:,:,i)  + 0.5 * (1-Fil2) * (Ausiliar_Velocity +  Elem%Veloc3(:,:,:,i))
enddo

Elem%Veloc = Elem%Veloc1 + Elem%Veloc2 + Elem%Veloc3

return
end subroutine Correction_Elem_FPML_Veloc

!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine  compute_InternalForces_PML_Elem (Elem, hprimex, hTprimey,htprimez)

implicit none

type (Element), intent (INOUT) :: Elem
real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprimex
real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hTprimey
real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez

integer :: m1, m2, m3, n_z
real, dimension ( 0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1)  :: s0,s1


m1 = Elem%ngllx; m2 = Elem%nglly;  m3= Elem%ngllz

s0 = Elem%Acoeff(:,:,:,27) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,28) * Elem%Residual_Stress(:,:,:,0) + & 
     Elem%Acoeff(:,:,:,29) * Elem%Residual_Stress(:,:,:,1)
call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
Elem%Forces1(:,:,:,0) = s1

s0 = Elem%Acoeff(:,:,:,30) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,31) * Elem%Residual_Stress(:,:,:,0) + & 
     Elem%Acoeff(:,:,:,32) * Elem%Residual_Stress(:,:,:,1)
do n_z = 0,m3-1
     call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
enddo
Elem%Forces2(:,:,:,0) = s1

s0 = Elem%Acoeff(:,:,:,33) * Elem%Diagonal_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,34) * Elem%Residual_Stress(:,:,:,0) + & 
     Elem%Acoeff(:,:,:,35) * Elem%Residual_Stress(:,:,:,1)

call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
Elem%Forces3(:,:,:,0) = s1

s0 = Elem%Acoeff(:,:,:,27) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,28) * Elem%Diagonal_Stress(:,:,:,1) + & 
     Elem%Acoeff(:,:,:,29) * Elem%Residual_Stress(:,:,:,2)
call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
Elem%Forces1(:,:,:,1) = s1

s0 = Elem%Acoeff(:,:,:,30) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,31) * Elem%Diagonal_Stress(:,:,:,1) + & 
     Elem%Acoeff(:,:,:,32) * Elem%Residual_Stress(:,:,:,2)
do n_z = 0,m3-1
     call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
enddo
Elem%Forces2(:,:,:,1) = s1

s0 = Elem%Acoeff(:,:,:,33) * Elem%Residual_Stress(:,:,:,0) + Elem%Acoeff(:,:,:,34) * Elem%Diagonal_Stress(:,:,:,1) + &
     Elem%Acoeff(:,:,:,35) * Elem%Residual_Stress(:,:,:,2)

call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
Elem%Forces3(:,:,:,1) = s1


s0 = Elem%Acoeff(:,:,:,27) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,28) * Elem%Residual_Stress(:,:,:,2) + & 
           Elem%Acoeff(:,:,:,29) * Elem%Diagonal_Stress(:,:,:,2)
call DGEMM ( 'N', 'N', m1, m2*m3, m1, 1., hprimex, m1,s0(0,0,0) ,m1, 0., s1, m1 )
Elem%Forces1(:,:,:,2) = s1

s0 = Elem%Acoeff(:,:,:,30) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,31) * Elem%Residual_Stress(:,:,:,2) + & 
     Elem%Acoeff(:,:,:,32) * Elem%Diagonal_Stress(:,:,:,2)
do n_z = 0,m3-1
    call DGEMM ( 'N', 'N', m1, m2, m2, 1.,s0(0,0,n_z), m1, htprimey ,m2, 0., s1(0,0,n_z),m1 )
enddo
Elem%Forces2(:,:,:,2) = s1

s0 = Elem%Acoeff(:,:,:,33) * Elem%Residual_Stress(:,:,:,1) + Elem%Acoeff(:,:,:,34) * Elem%Residual_Stress(:,:,:,2) + & 
     Elem%Acoeff(:,:,:,35) * Elem%Diagonal_Stress(:,:,:,2)

call DGEMM ( 'N', 'N', m1*m2, m3, m3, 1., s0(0,0,0), m1*m2, htprimez ,m3, 0., s1, m1*m2 )
Elem%Forces3(:,:,:,2) = s1

Elem%Forces = Elem%Forces1 + Elem%Forces2 + Elem%Forces3

return
end subroutine compute_InternalForces_PML_Elem
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------
subroutine compute_InternalForces_PML_Elem_Fl(Elem,hprimex,hTprimey,htprimez)

    implicit none

    type(Element), intent(inout) :: Elem
    real, dimension(0:Elem%ngllx-1,0:Elem%ngllx-1), intent(in) :: hprimex
    real, dimension(0:Elem%nglly-1,0:Elem%nglly-1), intent(in) :: hTprimey
    real, dimension(0:Elem%ngllz-1,0:Elem%ngllz-1), intent(in) :: hTprimez

    integer :: m1, m2, m3, n_z
    real, dimension(0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1)  :: s0,s1


m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

! forces associated to V_x
s0 = Elem%Acoeff(:,:,:,9) * Elem%Veloc(:,:,:,0)
call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,s0(0,0,0),m1,0.,s1,m1)
Elem%ForcesFl1(:,:,:) = s1

s0 = Elem%Acoeff(:,:,:,10) * Elem%Veloc(:,:,:,0)
do n_z = 0,m3-1
     call DGEMM('N','N',m1,m2,m2,1.,s0(0,0,n_z),m1, htprimey,m2,0.,s1(0,0,n_z),m1)
enddo
Elem%ForcesFl1(:,:,:) = s1+Elem%ForcesFl1(:,:,:)

s0 = Elem%Acoeff(:,:,:,11) * Elem%Veloc(:,:,:,0)
call DGEMM('N','N',m1*m2,m3,m3,1.,s0(0,0,0),m1*m2,htprimez,m3,0.,s1,m1*m2)
Elem%ForcesFl1(:,:,:) = s1+Elem%ForcesFl1(:,:,:)

! forces associated to V_y
s0 = Elem%Acoeff(:,:,:,12) * Elem%Veloc(:,:,:,1)
call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,s0(0,0,0),m1,0.,s1,m1)
Elem%ForcesFl2(:,:,:) = s1

s0 = Elem%Acoeff(:,:,:,13) * Elem%Veloc(:,:,:,1)
do n_z = 0,m3-1
     call DGEMM('N','N',m1,m2,m2,1.,s0(0,0,n_z),m1, htprimey,m2,0.,s1(0,0,n_z),m1)
enddo
Elem%ForcesFl2(:,:,:) = s1+Elem%ForcesFl2(:,:,:)

s0 = Elem%Acoeff(:,:,:,14) * Elem%Veloc(:,:,:,1)
call DGEMM('N','N',m1*m2,m3,m3,1.,s0(0,0,0),m1*m2,htprimez,m3,0.,s1,m1*m2)
Elem%ForcesFl2(:,:,:) = s1+Elem%ForcesFl2(:,:,:)

! forces associated to V_z
s0 = Elem%Acoeff(:,:,:,15) * Elem%Veloc(:,:,:,2)
call DGEMM('N','N',m1,m2*m3,m1,1.,hprimex,m1,s0(0,0,0),m1,0.,s1,m1)
Elem%ForcesFl3(:,:,:) = s1

s0 = Elem%Acoeff(:,:,:,16) * Elem%Veloc(:,:,:,2)
do n_z = 0,m3-1
     call DGEMM('N','N',m1,m2,m2,1.,s0(0,0,n_z),m1, htprimey,m2,0.,s1(0,0,n_z),m1)
enddo
Elem%ForcesFl3(:,:,:) = s1+Elem%ForcesFl3(:,:,:)

s0 = Elem%Acoeff(:,:,:,17) * Elem%Veloc(:,:,:,2)
call DGEMM('N','N',m1*m2,m3,m3,1.,s0(0,0,0),m1*m2,htprimez,m3,0.,s1,m1*m2)
Elem%ForcesFl3(:,:,:) = s1+Elem%ForcesFl3(:,:,:)


Elem%ForcesFl(:,:,:) = Elem%ForcesFl1(:,:,:) + Elem%ForcesFl2(:,:,:) + Elem%ForcesFl3(:,:,:)

return
end subroutine compute_InternalForces_PML_Elem_Fl
!------------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------------

end module selement

