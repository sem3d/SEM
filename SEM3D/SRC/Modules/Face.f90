module sfaces


type :: face

logical :: PML, Abs, FPML

integer :: ngll1, ngll2, dir, Which_Elem
integer, dimension (:), pointer :: FaceNum
integer, dimension (:,:), pointer :: Iglobnum_Face

real, dimension (:,:), pointer  :: MassMat
real, dimension (:,:,:), pointer :: Forces, Displ, Veloc, Accel, V0
real, dimension (:,:,:), pointer :: Forces1, Forces2, Forces3, Veloc1, Veloc2, Veloc3
real, dimension (:,:,:), pointer :: DumpVx, DumpVy, DumpVz, DumpMass
real, dimension (:,:,:), pointer :: IVeloc1, IVeloc2, IVeloc3
real, dimension (:,:), pointer :: Ivx, Ivy, Ivz

end type
 
contains

! ############################################################
subroutine Prediction_Face_Veloc (F, alpha, bega, gam1, dt)

implicit none

type (Face), intent (INOUT) :: F
real, intent (IN) :: alpha, bega, gam1, dt
  
integer :: ngll1, ngll2


ngll1 = F%ngll1 ; ngll2 =F%ngll2
 
F%Forces(:,:,:) = F%Displ(:,:,:) + dt * F%Veloc(:,:,:) + dt**2 * (0.5 - bega) * F%Accel(:,:,:)
F%V0(:,:,:) = F%Veloc(:,:,:)
F%Forces(:,:,:) = alpha * F%Forces(:,:,:) + (1-alpha) * F%Displ(:,:,:) 

return
end subroutine Prediction_Face_Veloc

! ###########################################################
subroutine Correction_Face_Veloc (F, bega, gam1, dt)

implicit none

type (Face), intent (INOUT) :: F
real, intent (IN) :: bega, gam1, dt

integer :: i, ngll1, ngll2


ngll1 = F%ngll1 ; ngll2 = F%ngll2

do i = 0,2
    F%Forces(:,:,i) = F%MassMat(:,:) * F%Forces(:,:,i)
enddo
F%Veloc(:,:,:) = F%v0(:,:,:) + dt * F%Forces(:,:,:)
F%Accel(:,:,:) = F%Accel(:,:,:) + gam1 / dt * (F%Veloc(:,:,:)-F%V0(:,:,:))
F%Displ(:,:,:) = F%Displ(:,:,:) + bega * dt * (F%Veloc(:,:,:)+F%V0(:,:,:))

return
end subroutine Correction_Face_Veloc

! ##########################################################
subroutine Correction_Face_PML_Veloc (F, dt)

implicit none

type (Face), intent (INOUT) :: F
real, intent (IN) :: dt

integer :: i


do i = 0,2
    F%Veloc1(:,:,i) = F%DumpVx(:,:,0) * F%Veloc1(:,:,i) + dt * F%DumpVx(:,:,1) * F%Forces1(:,:,i)
    F%Veloc2(:,:,i) = F%DumpVy(:,:,0) * F%Veloc2(:,:,i) + dt * F%DumpVy(:,:,1) * F%Forces2(:,:,i)
    F%Veloc3(:,:,i) = F%DumpVz(:,:,0) * F%Veloc3(:,:,i) + dt * F%DumpVz(:,:,1) * F%Forces3(:,:,i)
enddo

F%Veloc = F%Veloc1 + F%Veloc2 + F%Veloc3

if (F%Abs) then
    F%Veloc = 0
endif

return
end subroutine Correction_Face_PML_Veloc

! ############################################################
subroutine Correction_Face_FPML_Veloc (F, dt, fil)

implicit none

type (Face), intent (INOUT) :: F
real, intent (IN) :: dt,fil

integer :: i
real :: fil2
real, dimension (1:F%ngll1-2,1:F%ngll2-2) :: Ausiliar_velocity


do i = 0,2
    Ausiliar_Velocity = F%Veloc1(:,:,i)
    F%Veloc1(:,:,i) = F%DumpVx(:,:,0) * F%Veloc1(:,:,i) + dt * F%DumpVx(:,:,1) * F%Forces1(:,:,i) + F%Ivx * F%Iveloc1(:,:,i)
    F%Iveloc1(:,:,i) = Fil2 * F%Iveloc1(:,:,i) + 0.5 * (1-Fil2) * (Ausiliar_Velocity + F%Veloc1(:,:,i))

    Ausiliar_Velocity = F%Veloc2(:,:,i)
    F%Veloc2(:,:,i) = F%DumpVy(:,:,0) * F%Veloc2(:,:,i) + dt * F%DumpVy(:,:,1) * F%Forces2(:,:,i) + F%Ivy * F%Iveloc2(:,:,i)
    F%Iveloc2(:,:,i) = Fil2 * F%Iveloc2(:,:,i) + 0.5 * (1-Fil2) * (Ausiliar_Velocity + F%Veloc2(:,:,i))

    Ausiliar_Velocity = F%Veloc3(:,:,i)
    F%Veloc3(:,:,i) = F%DumpVz(:,:,0) * F%Veloc3(:,:,i) + dt * F%DumpVz(:,:,1) * F%Forces3(:,:,i) + F%Ivz * F%Iveloc3(:,:,i)
    F%Iveloc3(:,:,i) = Fil2 * F%Iveloc3(:,:,i) + 0.5 * (1-Fil2) * (Ausiliar_Velocity + F%Veloc3(:,:,i))
enddo

F%Veloc = F%Veloc1 + F%Veloc2 + F%Veloc3

if (F%Abs) then
    F%Veloc = 0
endif

return
end subroutine Correction_Face_FPML_Veloc

! ############################################################



subroutine get_vel_face (F, Vfree, ngll1, ngll2, dt, logic, orient, nf)
implicit none

integer, intent (IN) :: ngll1, ngll2, orient, nf
real, intent (IN) :: dt
type (Face), intent (IN) :: F
real, dimension (1:ngll1-2,1:ngll2-2,0:2), intent (INOUT) :: Vfree
logical, intent (IN) :: logic 

integer :: i,j

if (.not. F%PML) then

 if (logic) then
   select case (orient)
     case (0)
       do i = 0,2
          Vfree(1:ngll1-2,1:ngll2-2,i) = Vfree(1:ngll1-2,1:ngll2-2,i) - ( F%V0(1:ngll1-2,1:ngll2-2,i) + &
                                         dt*F%MassMat(1:ngll1-2,1:ngll2-2)*F%Forces(1:ngll1-2,1:ngll2-2,i) )
       enddo
      case (1)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(ngll1-1-i,j,0:2) + dt*F%MassMat(ngll1-1-i,j)*F%Forces(ngll1-1-i,j,0:2) )
         enddo
       enddo
      case (2)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(i,ngll2-1-j,0:2) + dt*F%MassMat(i,ngll2-1-j)*F%Forces(i,ngll2-1-j,0:2) )
         enddo
       enddo
      case (3)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(ngll1-1-i,ngll2-1-j,0:2) + dt*F%MassMat(ngll1-1-i,ngll2-1-j)*F%Forces(ngll1-1-i,ngll2-1-j,0:2) )
         enddo
       enddo
      case (4)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(j,i,0:2) + dt*F%MassMat(j,i)*F%Forces(j,i,0:2) )
         enddo
       enddo
      case (5)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(ngll1-1-j,i,0:2) + dt*F%MassMat(ngll1-1-j,i)*F%Forces(ngll1-1-j,i,0:2) )
         enddo
       enddo
      case (6)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(j,ngll2-1-i,0:2) + dt*F%MassMat(j,ngll2-1-i)*F%Forces(j,ngll2-1-i,0:2) )
         enddo
       enddo
      case (7)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%V0(ngll1-1-j,ngll2-1-i,0:2) + dt*F%MassMat(ngll1-1-j,ngll2-1-i)*F%Forces(ngll1-1-j,ngll2-1-i,0:2) )
         enddo
       enddo
     end select
 else
   do i = 0,2
        Vfree(1:ngll1-2,1:ngll2-2,i) =  F%V0(1:ngll1-2,1:ngll2-2,i) + dt*F%MassMat(1:ngll1-2,1:ngll2-2)*F%Forces(1:ngll1-2,1:ngll2-2,i)
   enddo
 endif

else

 if (logic) then
   select case (orient)
     case (0)
       do i = 0,2
          Vfree(1:ngll1-2,1:ngll2-2,i) = Vfree(1:ngll1-2,1:ngll2-2,i) - ( F%DumpVz(1:ngll1-2,1:ngll2-2,0) * F%Veloc3(1:ngll1-2,1:ngll2-2,i) + &
                                                                  dt * F%DumpVz(1:ngll1-2,1:ngll2-2,1) * F%Forces3(1:ngll1-2,1:ngll2-2,i) )
       enddo
      case (1)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(ngll1-1-i,j,0) * F%Veloc3(ngll1-1-i,j,0:2) + &
                                               dt * F%DumpVz(ngll1-1-i,j,1) * F%Forces3(ngll1-1-i,j,0:2) )
         enddo
       enddo
      case (2)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(i,ngll2-1-j,0) * F%Veloc3(i,ngll2-1-j,0:2) + &
                                               dt * F%DumpVz(i,ngll2-1-j,1) * F%Forces3(i,ngll2-1-j,0:2) )
         enddo
       enddo
      case (3)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(ngll1-1-i,ngll2-1-j,0) * F%Veloc3(ngll1-1-i,ngll2-1-j,0:2) + &
                                               dt * F%DumpVz(ngll1-1-i,ngll2-1-j,1) * F%Forces3(ngll1-1-i,ngll2-1-j,0:2) )
         enddo
       enddo
      case (4)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(j,i,0) * F%Veloc3(j,i,0:2) + &
                                               dt * F%DumpVz(j,i,1) * F%Forces3(j,i,0:2) )
         enddo
       enddo
      case (5)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(ngll2-1-j,i,0) * F%Veloc3(ngll2-1-j,i,0:2) + &
                                               dt * F%DumpVz(ngll2-1-j,i,1) * F%Forces3(ngll2-1-j,i,0:2) )
         enddo
       enddo
      case (6)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(j,ngll1-1-i,0) * F%Veloc3(j,ngll1-1-i,0:2) + &
                                               dt * F%DumpVz(j,ngll1-1-i,1) * F%Forces3(j,ngll1-1-i,0:2) )
         enddo
       enddo
      case (7)
       do j = 1, ngll2-2
         do i = 1, ngll1-2
	   Vfree(i,j,0:2) = Vfree(i,j,0:2) - ( F%DumpVz(ngll2-1-j,ngll1-1-i,0) * F%Veloc3(ngll2-1-j,ngll1-1-i,0:2) + &
                                               dt * F%DumpVz(ngll2-1-j,ngll1-1-i,1) * F%Forces3(ngll2-1-j,ngll1-1-i,0:2) )
         enddo
       enddo
     end select
 else
   do i = 0,2
        Vfree(1:ngll1-2,1:ngll2-2,i) =  F%DumpVz(1:ngll1-2,1:ngll2-2,0) * F%Veloc3(1:ngll1-2,1:ngll2-2,i) + &
                                                                  dt * F%DumpVz(1:ngll1-2,1:ngll2-2,1) * F%Forces3(1:ngll1-2,1:ngll2-2,i)
   enddo
 endif

endif

return
end subroutine get_vel_face

! ############################################################


end module sfaces
