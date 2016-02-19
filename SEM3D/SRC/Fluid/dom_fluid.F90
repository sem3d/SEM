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

    subroutine allocate_dom_fluid (dom)
        implicit none
        type(domain_fluid), intent (INOUT) :: dom
        !
        integer nbelem, ngll
        !

        nbelem = dom%nbelem
        if(nbelem == 0) return ! Do not allocate if not needed (save allocation/RAM)
        ngll   = dom%ngll

        allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nbelem-1))
        allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nbelem-1))

        allocate (dom%Jacob_  (        0:ngll-1,0:ngll-1,0:ngll-1,0:nbelem-1))
        allocate (dom%InvGrad_(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nbelem-1))

        allocate(dom%Idom_(0:ngll-1,0:ngll-1,0:ngll-1,0:nbelem-1))

        ! Allocation et initialisation de champs0 et champs1 pour les fluides
        if (dom%nglltot /= 0) then
            allocate(dom%champs0%ForcesFl(0:dom%nglltot-1))
            allocate(dom%champs0%Phi     (0:dom%nglltot-1))
            allocate(dom%champs0%VelPhi  (0:dom%nglltot-1))
            allocate(dom%champs1%ForcesFl(0:dom%nglltot-1))
            allocate(dom%champs1%Phi     (0:dom%nglltot-1))
            allocate(dom%champs1%VelPhi  (0:dom%nglltot-1))

            dom%champs0%ForcesFl = 0d0
            dom%champs0%Phi = 0d0
            dom%champs0%VelPhi = 0d0

            ! Allocation de MassMat pour les fluides
            allocate(dom%MassMat(0:dom%nglltot-1))
            dom%MassMat = 0d0
        endif
    end subroutine allocate_dom_fluid

    subroutine deallocate_dom_fluid (dom)
        implicit none
        type(domain_fluid), intent (INOUT) :: dom

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )

        if(allocated(dom%m_Jacob  )) deallocate(dom%m_Jacob  )
        if(allocated(dom%m_InvGrad)) deallocate(dom%m_InvGrad)

        if(allocated(dom%m_Idom)) deallocate(dom%m_Idom)

        if(allocated(dom%champs0%ForcesFl)) deallocate(dom%champs0%ForcesFl)
        if(allocated(dom%champs0%Phi     )) deallocate(dom%champs0%Phi     )
        if(allocated(dom%champs0%VelPhi  )) deallocate(dom%champs0%VelPhi  )
        if(allocated(dom%champs1%ForcesFl)) deallocate(dom%champs1%ForcesFl)
        if(allocated(dom%champs1%Phi     )) deallocate(dom%champs1%Phi     )
        if(allocated(dom%champs1%VelPhi  )) deallocate(dom%champs1%VelPhi  )

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)
    end subroutine deallocate_dom_fluid

    subroutine fluid_velocity(ngll,htprime,InvGrad,density,phi,veloc)
        ! gives the physical particle velocity in the fluid = 1/dens grad(dens.Phi)
        implicit none
        integer, intent(in)  :: ngll
        real, dimension(0:ngll-1,0:ngll-1), intent(in) :: htprime
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: density,phi
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(out) :: Veloc
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: dphi_dx,dphi_dy,dphi_dz

        ! physical gradient
        call physical_part_deriv(ngll,htprime,InvGrad,phi,dphi_dx,dphi_dy,dphi_dz)
        Veloc(:,:,:,0) = dphi_dx(:,:,:)/density(:,:,:)
        Veloc(:,:,:,1) = dphi_dy(:,:,:)/density(:,:,:)
        Veloc(:,:,:,2) = dphi_dz(:,:,:)/density(:,:,:)
    end subroutine fluid_velocity

    subroutine get_fluid_dom_var(Tdomain, dom, el, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain)                               :: TDomain
        type(domain_fluid), intent(inout)          :: dom
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
        integer :: ngll, i, j, k, ind, mat

        flag_gradU = (out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS) + &
            out_variables(OUT_EPS_VOL) + &
            out_variables(OUT_EPS_DEV) + &
            out_variables(OUT_STRESS_DEV)) /= 0

        ngll = dom%ngll

        if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
            if(.not. allocated(fieldU)) allocate(fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
            fieldU(:,:,:,:) = 0d0
        end if

        if (out_variables(OUT_VITESSE) == 1) then
            if(.not. allocated(fieldV)) allocate(fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
            if(.not. allocated(phi))    allocate(phi(0:ngll-1,0:ngll-1,0:ngll-1))
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,el%lnum)
                        phi(i,j,k) = dom%champs0%Phi(ind)
                    enddo
                enddo
            enddo
            mat = el%mat_index
            call fluid_velocity(ngll,Tdomain%sSubdomain(mat)%htprime,                  &
                dom%InvGrad_(:,:,:,:,:,el%lnum),dom%Density_(:,:,:,el%lnum),phi,fieldV)
        end if

        if (out_variables(OUT_ACCEL) == 1) then
            if(.not. allocated(fieldA)) allocate(fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
            if(.not. allocated(vphi))   allocate(vphi(0:ngll-1,0:ngll-1,0:ngll-1))
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,el%lnum)
                        vphi(i,j,k) = dom%champs0%VelPhi(ind)
                    enddo
                enddo
            enddo
            mat = el%mat_index
            call fluid_velocity(ngll,Tdomain%sSubdomain(mat)%htprime,               &
                dom%InvGrad_(:,:,:,:,:,el%lnum),dom%Density_(:,:,:,el%lnum),vphi,fieldA)
        end if

        if (out_variables(OUT_PRESSION) == 1) then
            if(.not. allocated(fieldP)) allocate(fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,el%lnum)
                        fieldP(i,j,k) = -dom%champs0%VelPhi(ind)
                    enddo
                enddo
            enddo
        end if

        if (out_variables(OUT_EPS_VOL) == 1) then
            if(.not. allocated(eps_vol)) allocate(eps_vol(0:ngll-1,0:ngll-1,0:ngll-1))
            eps_vol(:,:,:) = 0.
        end if

        if (out_variables(OUT_ENERGYP) == 1) then
            if(.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
            P_energy(:,:,:) = 0.
        end if

        if (out_variables(OUT_ENERGYS) == 1) then
            if(.not. allocated(S_energy)) allocate(S_energy(0:ngll-1,0:ngll-1,0:ngll-1))
            S_energy(:,:,:) = 0.
        end if

        if (out_variables(OUT_EPS_DEV) == 1) then
            if(.not. allocated(eps_dev)) allocate(eps_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
            eps_dev(:,:,:,:) = 0.
        end if

        if (out_variables(OUT_STRESS_DEV) == 1) then
            if(.not. allocated(sig_dev)) allocate(sig_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
            sig_dev(:,:,:,:) = 0.
        end if
        if(allocated(phi))  deallocate(phi)
        if(allocated(vphi)) deallocate(vphi)
    end subroutine get_fluid_dom_var

    subroutine init_material_properties_fluid(dom, lnum, i, j, k, density, lambda)
        type(domain_fluid), intent(inout) :: dom
        integer, intent(in) :: lnum
        integer, intent(in) :: i, j, k ! -1 means :
        real(fpp), intent(in) :: density
        real(fpp), intent(in) :: lambda

        if (i==-1 .and. j==-1 .and. k==-1) then
            dom%Density_(:,:,:,lnum) = density
            dom%Lambda_ (:,:,:,lnum) = lambda
        else
            dom%Density_(i,j,k,lnum) = density
            dom%Lambda_ (i,j,k,lnum) = lambda
        end if
    end subroutine init_material_properties_fluid

    subroutine init_local_mass_fluid(dom,specel,i,j,k,ind,Whei)
        type(domain_fluid), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real Whei

        ! Fluid : inertial term ponderation by the inverse of the bulk modulus

        specel%MassMat(i,j,k) = Whei*dom%Jacob_(i,j,k,specel%lnum)/dom%Lambda_(i,j,k,specel%lnum)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_fluid

    subroutine forces_int_fluid(dom, mat, champs1, lnum)
        use m_calcul_forces_fluid
        type(domain_fluid), intent (INOUT) :: dom
        type (subdomain), intent(IN) :: mat
        type(champsfluid), intent(inout) :: champs1
        integer :: lnum

        integer :: ngll,i,j,k
        real, dimension(0:dom%ngll-1, 0:dom%ngll-1, 0:dom%ngll-1) :: Fo_Fl,Phi

        ngll = dom%ngll

        ! d(rho*Phi)_dX
        ! d(rho*Phi)_dY
        ! d(rho*Phi)_dZ
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    Phi(i,j,k) = champs1%Phi(dom%Idom_(i,j,k,lnum))
                enddo
            enddo
        enddo

        ! internal forces
        call calcul_forces_fluid(dom,mat,lnum,Fo_Fl,mat%GLLw,Phi)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    champs1%ForcesFl(dom%Idom_(i,j,k,lnum)) = champs1%ForcesFl(dom%Idom_(i,j,k,lnum))-Fo_Fl(i,j,k)
                enddo
            enddo
        enddo
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
