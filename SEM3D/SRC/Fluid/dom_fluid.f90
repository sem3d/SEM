!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_fluid
    use sdomain
    use constants
    use champs_fluid
    use selement
    use ssubdomains
    implicit none

contains

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
            phi,dphi_dx,dphi_dy,dphi_dz)

        !
        Veloc(:,:,:,0) = dphi_dx(:,:,:)/density(:,:,:)
        Veloc(:,:,:,1) = dphi_dy(:,:,:)/density(:,:,:)
        Veloc(:,:,:,2) = dphi_dz(:,:,:)/density(:,:,:)

    end subroutine fluid_velocity

    subroutine get_fluid_dom_var(Tdomain, el, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain)                               :: TDomain
        integer, dimension(0:8)                    :: out_variables
        type(element)                              :: el
        real(fpp), dimension(:,:,:), allocatable   :: phi
        real(fpp), dimension(:,:,:), allocatable   :: vphi
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp)                                  :: P_energy, S_energy, eps_vol
        real(fpp), dimension(0:5)                  :: eps_dev
        real(fpp), dimension(0:5)                  :: sig_dev
        !
        logical :: flag_gradU
        integer :: nx, ny, nz, i, j, k, ind, mat

        flag_gradU = (out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS) + &
            out_variables(OUT_EPS_VOL) + &
            out_variables(OUT_EPS_DEV) + &
            out_variables(OUT_STRESS_DEV)) /= 0

        nx = el%ngllx
        ny = el%nglly
        nz = el%ngllz

        do k=0,nz-1
            do j=0,ny-1
                do i=0,nx-1
                    ind = el%Idom(i,j,k)

                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        if(.not. allocated(fieldU)) allocate(fieldU(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldU(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_VITESSE) == 1) then
                        if(.not. allocated(fieldV)) allocate(fieldV(0:nx-1,0:ny-1,0:nz-1,0:2))
                        if(.not. allocated(phi))    allocate(phi(0:nx-1,0:ny-1,0:nz-1))
                        phi(i,j,k) = Tdomain%fdom%champs0%Phi(ind)
                        mat = el%mat_index
                        call fluid_velocity(nx,ny,nz,Tdomain%sSubdomain(mat)%htprimex,              &
                            Tdomain%sSubdomain(mat)%hprimey,Tdomain%sSubdomain(mat)%hprimez, &
                            el%InvGrad,el%density,phi,fieldV)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:nx-1,0:ny-1,0:nz-1,0:2))
                        if(.not. allocated(vphi))   allocate(vphi(0:nx-1,0:ny-1,0:nz-1))
                        vphi(i,j,k) = Tdomain%fdom%champs0%VelPhi(ind)
                        mat = el%mat_index
                        call fluid_velocity(nx,ny,nz,Tdomain%sSubdomain(mat)%htprimex,              &
                            Tdomain%sSubdomain(mat)%hprimey,Tdomain%sSubdomain(mat)%hprimez, &
                            el%InvGrad,el%density,vphi,fieldA)
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:nx-1,0:ny-1,0:nz-1))
                        fieldP(i,j,k) = -Tdomain%fdom%champs0%VelPhi(ind)
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        eps_vol = 0.
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        P_energy = 0.
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        S_energy = 0.
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        eps_dev = 0.
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        sig_dev = 0.
                    end if
                enddo
            enddo
        enddo
        if(allocated(phi))  deallocate(phi)
        if(allocated(vphi)) deallocate(vphi)
    end subroutine get_fluid_dom_var

    subroutine forces_int_fluid(Elem, mat, htprimex, hprimey, htprimey, hprimez, htprimez,  &
        champs1)

        type (Element), intent (INOUT) :: Elem
        type (subdomain), intent(IN) :: mat
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: htprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        type(champsfluid), intent(inout) :: champs1

        integer :: m1,m2,m3, i,j,k
        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dPhiX,dPhiY,dPhiZ,Fo_Fl

        real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: Phi


        m1 = Elem%ngllx;   m2 = Elem%nglly;   m3 = Elem%ngllz

        ! d(rho*Phi)_dX
        ! d(rho*Phi)_dY
        ! d(rho*Phi)_dZ
        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    Phi(i,j,k) = champs1%Phi(Elem%Idom(i,j,k))
                enddo
            enddo
        enddo
        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad, Phi, &
            dPhiX, dPhiY, dPhiZ)

        ! internal forces
        call calcul_forces_fluid(Fo_Fl,                &
            Elem%Invgrad, &
            htprimex,htprimey,htprimez, &
            Elem%Jacob,mat%GLLwx,mat%GLLwy,mat%GLLwz, &
            dPhiX,dPhiY,dPhiZ,       &
            Elem%Density,            &
            m1,m2,m3)

        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    champs1%ForcesFl(Elem%Idom(i,j,k)) = champs1%ForcesFl(Elem%Idom(i,j,k))-Fo_Fl(i,j,k)
                enddo
            enddo
        enddo


        return
    end subroutine forces_int_fluid

end module dom_fluid

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
