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
#include "index.h"

contains

  subroutine allocate_dom_fluid (Tdomain, dom)
        implicit none
        type(domain) :: TDomain
        type(domain_fluid), intent (INOUT) :: dom
        !
        integer nbelem, ngllx, nglly, ngllz
        !

        dom%ngllx = Tdomain%specel(0)%ngllx ! Temporaire: ngll* doit passer sur le domaine a terme
        dom%nglly = Tdomain%specel(0)%nglly ! Temporaire: ngll* doit passer sur le domaine a terme
        dom%ngllz = Tdomain%specel(0)%ngllz ! Temporaire: ngll* doit passer sur le domaine a terme

        nbelem  = dom%nbelem
        if(nbelem == 0) return ! Do not allocate if not needed (save allocation/RAM)
        ngllx   = dom%ngllx
        nglly   = dom%nglly
        ngllz   = dom%ngllz

        allocate(dom%Density_(0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))
        allocate(dom%Lambda_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))
        allocate(dom%Mu_     (0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))
        allocate(dom%Kappa_  (0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))

        allocate (dom%Jacob_  (        0:ngllx-1,0:nglly-1,0:ngllz-1,0:nbelem-1))
        allocate (dom%InvGrad_(0:2,0:2,0:ngllx-1,0:nglly-1,0:ngllz-1,0:nbelem-1))

        ! Allocation et initialisation de champs0 et champs1 pour les fluides
        if (dom%ngll /= 0) then
            allocate(dom%champs0%ForcesFl(0:dom%ngll-1))
            allocate(dom%champs0%Phi     (0:dom%ngll-1))
            allocate(dom%champs0%VelPhi  (0:dom%ngll-1))
            allocate(dom%champs1%ForcesFl(0:dom%ngll-1))
            allocate(dom%champs1%Phi     (0:dom%ngll-1))
            allocate(dom%champs1%VelPhi  (0:dom%ngll-1))

            dom%champs0%ForcesFl = 0d0
            dom%champs0%Phi = 0d0
            dom%champs0%VelPhi = 0d0

            ! Allocation de MassMat pour les fluides
            allocate(dom%MassMat(0:dom%ngll-1))
            dom%MassMat = 0d0
        endif
    end subroutine allocate_dom_fluid

    subroutine deallocate_dom_fluid (dom)
        implicit none
        type(domain_fluid), intent (INOUT) :: dom

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_Mu     )) deallocate(dom%m_Mu     )
        if(allocated(dom%m_Kappa  )) deallocate(dom%m_Kappa  )

        if(allocated(dom%m_Jacob  )) deallocate(dom%m_Jacob  )
        if(allocated(dom%m_InvGrad)) deallocate(dom%m_InvGrad)

        if(allocated(dom%champs0%ForcesFl)) deallocate(dom%champs0%ForcesFl)
        if(allocated(dom%champs0%Phi     )) deallocate(dom%champs0%Phi     )
        if(allocated(dom%champs0%VelPhi  )) deallocate(dom%champs0%VelPhi  )
        if(allocated(dom%champs1%ForcesFl)) deallocate(dom%champs1%ForcesFl)
        if(allocated(dom%champs1%Phi     )) deallocate(dom%champs1%Phi     )
        if(allocated(dom%champs1%VelPhi  )) deallocate(dom%champs1%VelPhi  )

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)
    end subroutine deallocate_dom_fluid

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
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
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
                        call fluid_velocity(nx,ny,nz,Tdomain%sSubdomain(mat)%htprimex,                  &
                            Tdomain%sSubdomain(mat)%hprimey,Tdomain%sSubdomain(mat)%hprimez,            &
                            Tdomain%fdom%InvGrad_(:,:,:,:,:,el%lnum),Tdomain%fdom%Density_(:,:,:,el%lnum),&
                            phi,fieldV)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:nx-1,0:ny-1,0:nz-1,0:2))
                        if(.not. allocated(vphi))   allocate(vphi(0:nx-1,0:ny-1,0:nz-1))
                        vphi(i,j,k) = Tdomain%fdom%champs0%VelPhi(ind)
                        mat = el%mat_index
                        call fluid_velocity(nx,ny,nz,Tdomain%sSubdomain(mat)%htprimex,                  &
                            Tdomain%sSubdomain(mat)%hprimey,Tdomain%sSubdomain(mat)%hprimez,            &
                            Tdomain%fdom%InvGrad_(:,:,:,:,:,el%lnum),Tdomain%fdom%Density_(:,:,:,el%lnum),&
                            vphi,fieldA)
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:nx-1,0:ny-1,0:nz-1))
                        fieldP(i,j,k) = -Tdomain%fdom%champs0%VelPhi(ind)
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        if(.not. allocated(eps_vol)) allocate(eps_vol(0:nx-1,0:ny-1,0:nz-1))
                        eps_vol(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        if(.not. allocated(P_energy)) allocate(P_energy(0:nx-1,0:ny-1,0:nz-1))
                        P_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        if(.not. allocated(S_energy)) allocate(S_energy(0:nx-1,0:ny-1,0:nz-1))
                        S_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        if(.not. allocated(eps_dev)) allocate(eps_dev(0:nx-1,0:ny-1,0:nz-1,0:5))
                        eps_dev(i,j,k,:) = 0.
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        if(.not. allocated(sig_dev)) allocate(sig_dev(0:nx-1,0:ny-1,0:nz-1,0:5))
                        sig_dev(i,j,k,:) = 0.
                    end if
                enddo
            enddo
        enddo
        if(allocated(phi))  deallocate(phi)
        if(allocated(vphi)) deallocate(vphi)
    end subroutine get_fluid_dom_var

    subroutine forces_int_fluid(dom, mat, htprimex, hprimey, htprimey, hprimez, htprimez,  &
        champs1, Elem, lnum)

        type(domain_fluid), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: Elem
        type (subdomain), intent(IN) :: mat
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: htprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        type(champsfluid), intent(inout) :: champs1
        integer :: lnum

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
        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez, &
                                 dom%InvGrad_(:,:,:,:,:,lnum), Phi, &
                                 dPhiX, dPhiY, dPhiZ)

        ! internal forces
        call calcul_forces_fluid(dom,lnum,Fo_Fl,htprimex,htprimey,htprimez,&
             mat%GLLwx,mat%GLLwy,mat%GLLwz,dPhiX,dPhiY,dPhiZ,m1,m2,m3)
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
