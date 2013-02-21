    subroutine elem_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,Scalp,    &
                      dScalp_dxi,dScalp_deta,dScalp_dzeta)
    !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dScalp_dxi,dScalp_deta,dScalp_dzeta
        integer   :: n_z

      ! d(Scalp)_dxi
        call DGEMM('N','N',ngllx,nglly*ngllz,ngllx,1d0,htprimex,ngllx,Scalp(:,:,:),ngllx,0d0,dScalp_dxi,ngllx)
      ! d(Scalp)_deta
        do n_z = 0,ngllz-1
            call DGEMM('N','N',ngllx,nglly,nglly,1d0,Scalp(:,:,n_z),ngllx,hprimey,nglly,0d0,dScalp_deta(:,:,n_z),ngllx)
        enddo
      ! d(Scalp)_dzeta
        call DGEMM('N','N',ngllx*nglly,ngllz,ngllz,1d0,Scalp(:,:,:),ngllx*nglly,hprimez,ngllz,0d0,dScalp_dzeta,ngllx*nglly)


    end subroutine elem_part_deriv

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
    subroutine physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
                        Scalp,dScalp_dx,dScalp_dy,dScalp_dz)
    !- partial derivatives in the (x,y,z) space, of the scalar Scalp
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dScalp_dx,dScalp_dy,dScalp_dz
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: dScalp_dxi,dScalp_deta,dScalp_dzeta
  
    !- in the reference element
        call elem_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,Scalp,    &
                             dScalp_dxi,dScalp_deta,dScalp_dzeta)

    !- in the physical domain
    ! dScalp_dX
        dScalp_dx(:,:,:) = dScalp_dxi(:,:,:)*InvGrad(:,:,:,0,0) + dScalp_deta(:,:,:)*InvGrad(:,:,:,0,1) +  &
                  dScalp_dzeta(:,:,:)*InvGrad(:,:,:,0,2)
    ! dScalp_dY
        dScalp_dy(:,:,:) = dScalp_dxi(:,:,:)*InvGrad(:,:,:,1,0) + dScalp_deta(:,:,:)*InvGrad(:,:,:,1,1) +  &
                  dScalp_dzeta(:,:,:)*InvGrad(:,:,:,1,2)
    ! dScalp_dZ
        dScalp_dz(:,:,:) = dScalp_dxi(:,:,:)*InvGrad(:,:,:,2,0) + dScalp_deta(:,:,:)*InvGrad(:,:,:,2,1) +  &
                  dScalp_dzeta(:,:,:)*InvGrad(:,:,:,2,2)



    end subroutine physical_part_deriv
