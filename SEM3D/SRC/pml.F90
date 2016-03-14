module pml
    use constants
    implicit none

    contains

    subroutine define_alpha_PML(Coord, dir, ngll, vp, pml_width, pml_pos, Apow, npow, alpha)
        !- routine determines attenuation profile in an PML layer (see Festa & Vilotte)
        !   dir = attenuation's direction, ldir_attenu = the logical giving the orientation
        integer, intent(in) :: dir, ngll, npow
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1, 0:2), intent(in) :: Coord
        real(fpp), intent(in)  :: Apow
        real(fpp), dimension(0:2), intent(in) :: pml_pos, pml_width
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: vp
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: alpha
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1)  :: ri
        real(fpp) :: invdh, coef
        if (pml_width(dir)==0d0) then
            alpha(:,:,:) = 0d0
            return
        end if

        ! Until we can do better, the plane equation is only parallel to one axis
        ! when more coefficient are != 0 it means we have a corner
        invdh = 1d0/abs(pml_width(dir))
        coef = 1/pml_width(dir)
        ri = coef*(Coord(:,:,:,dir)-pml_pos(dir))
        alpha = Apow * Vp * invdh *  (ri)**npow
    end subroutine define_alpha_PML

    subroutine define_PML_DumpInit(ngll,dt,alpha,&
        MassMat,DumpS,DumpMass)
        !- defining parameters related to stresses and mass matrix elements, in the case of
        !    a PML, along a given splitted direction:
        integer, intent(in)  :: ngll
        real, intent(in) :: dt
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: alpha
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: MassMat
        real, dimension(0:1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: DumpS
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: DumpMass

        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1)  :: Id

        Id = 1d0

        DumpS(1,:,:,:) = Id + 0.5d0*dt*alpha
        DumpS(1,:,:,:) = 1d0/DumpS(1,:,:,:)
        DumpS(0,:,:,:) = (Id - 0.5d0*dt*alpha)*DumpS(1,:,:,:)
        DumpMass(:,:,:) = 0.5d0*MassMat(:,:,:)*alpha(:,:,:)*dt

        return
    end subroutine define_PML_DumpInit

end module pml
