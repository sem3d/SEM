subroutine fluid_velocity(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,InvGrad,    &
                          density,phi,veloc)
  ! gives the physical particle velocity in the fluid = 1/dens grad(dens.Phi)
    implicit none
    integer, intent(in)  :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
    real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
    real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: density,phi
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(out) :: Veloc
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: dphi_dx,dphi_dy,dphi_dz

    ! physical gradient
    call physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
                        density*phi,dphi_dx,dphi_dy,dphi_dz)
    
    !
    Veloc(:,:,:,0) = dphi_dx(:,:,:)/density(:,:,:) 
    Veloc(:,:,:,1) = dphi_dy(:,:,:)/density(:,:,:) 
    Veloc(:,:,:,2) = dphi_dz(:,:,:)/density(:,:,:) 
   
end subroutine fluid_velocity

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

subroutine grad_displ_solid(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,InvGrad,displ,    &
           dUx_dx,dUx_dy,dUx_dz,dUy_dx,dUy_dy,dUy_dz,dUz_dx,dUz_dy,dUz_dz)
  ! gives the displacement gradient in the solid part
    implicit none
    integer, intent(in)  :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
    real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
    real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(in) :: displ 
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dUx_dx,dUx_dy,dUx_dz,  &
                                           dUy_dx,dUy_dy,dUy_dz,dUz_dx,dUz_dy,dUz_dz

    call physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
                        displ(:,:,:,0),dUx_dx,dUx_dy,dUx_dz)
    call physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
                        displ(:,:,:,1),dUy_dx,dUy_dy,dUy_dz)
    call physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
                        displ(:,:,:,2),dUz_dx,dUz_dy,dUz_dz)


end subroutine grad_displ_solid

!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------

subroutine pressure_solid(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,InvGrad,displ,    &
           lambda,mu,pressure)
  ! gives the isotropic term of the stress tensor = pressure = -1/3 trace(stress tensor)
  !                         = -kappa * div(u)
    implicit none
    integer, intent(in)  :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
    real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
    real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(in) :: displ 
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: mu,lambda
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: pressure
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: dUx_dx,dUx_dy,dUx_dz,  &
                                           dUy_dx,dUy_dy,dUy_dz,dUz_dx,dUz_dy,dUz_dz


    call grad_displ_solid(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,InvGrad,displ,    &
           dUx_dx,dUx_dy,dUx_dz,dUy_dx,dUy_dy,dUy_dz,dUz_dx,dUz_dy,dUz_dz)

    pressure(:,:,:) = -(lambda(:,:,:)+2d0/3d0*mu(:,:,:))*(dUx_dx(:,:,:)+dUy_dy(:,:,:)+dUz_dz(:,:,:))

end subroutine pressure_solid
