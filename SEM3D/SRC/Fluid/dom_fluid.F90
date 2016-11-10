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
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_fluid), intent (INOUT) :: dom
        !
        integer :: nbelem, ngll, nblocks
        !

        ngll   = dom%ngll
        nbelem = dom%nbelem
        write(*,*) "DOM_FLUID ngll   = ", ngll
        write(*,*) "DOM_FLUID nbelem = ", nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call init_dombase(dom)

        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)
            nblocks = dom%nblocks
            allocate(dom%IDensity_(0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
        end if
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
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - fluid domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_fluid

    subroutine deallocate_dom_fluid (dom)
        implicit none
        type(domain_fluid), intent (INOUT) :: dom

        if(allocated(dom%m_IDensity)) deallocate(dom%m_IDensity)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )

        if(allocated(dom%champs0%ForcesFl)) deallocate(dom%champs0%ForcesFl)
        if(allocated(dom%champs0%Phi     )) deallocate(dom%champs0%Phi     )
        if(allocated(dom%champs0%VelPhi  )) deallocate(dom%champs0%VelPhi  )
        if(allocated(dom%champs1%ForcesFl)) deallocate(dom%champs1%ForcesFl)
        if(allocated(dom%champs1%Phi     )) deallocate(dom%champs1%Phi     )
        if(allocated(dom%champs1%VelPhi  )) deallocate(dom%champs1%VelPhi  )

        call deallocate_dombase(dom)
    end subroutine deallocate_dom_fluid

    subroutine fluid_velocity(ngll,hprime,InvGrad,idensity,phi,veloc)
        ! gives the physical particle velocity in the fluid = 1/dens grad(dens.Phi)
        implicit none
        integer, intent(in)  :: ngll
        real, dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: idensity,phi
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(out) :: Veloc
        real, dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: dphi_dx,dphi_dy,dphi_dz

        ! physical gradient
        call physical_part_deriv(ngll,hprime,InvGrad,phi,dphi_dx,dphi_dy,dphi_dz)
        Veloc(:,:,:,0) = dphi_dx(:,:,:) * idensity(:,:,:)
        Veloc(:,:,:,1) = dphi_dy(:,:,:) * idensity(:,:,:)
        Veloc(:,:,:,2) = dphi_dz(:,:,:) * idensity(:,:,:)
    end subroutine fluid_velocity

    subroutine get_fluid_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain)                               :: TDomain
        type(domain_fluid), intent(inout)          :: dom
        integer, dimension(0:), intent(in)         :: out_variables
        integer, intent(in)                        :: lnum
        real(fpp), dimension(:,:,:), allocatable   :: phi
        real(fpp), dimension(:,:,:), allocatable   :: vphi
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        !
        logical :: flag_gradU
        integer :: ngll, i, j, k, ind
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

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
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        phi(i,j,k) = dom%champs0%Phi(ind)
                    enddo
                enddo
            enddo
            call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                 dom%IDensity_(:,:,:,bnum,ee),phi,fieldV)
        end if

        if (out_variables(OUT_ACCEL) == 1) then
            if(.not. allocated(fieldA)) allocate(fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
            if(.not. allocated(vphi))   allocate(vphi(0:ngll-1,0:ngll-1,0:ngll-1))
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        vphi(i,j,k) = dom%champs0%VelPhi(ind)
                    enddo
                enddo
            enddo
            call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                 dom%IDensity_(:,:,:,bnum,ee),vphi,fieldA)
        end if

        if (out_variables(OUT_PRESSION) == 1) then
            if(.not. allocated(fieldP)) allocate(fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
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
            P_energy(:,:,:) = 0. !TODO
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


    subroutine get_fluid_dom_elem_energy(dom, lnum, P_energy, S_energy)
        use deriv3d
        implicit none
        !
        type(domain_fluid), intent(inout)          :: dom
        integer, intent(in)                        :: lnum
        real(fpp), dimension(:,:,:), allocatable, intent(inout) :: P_energy, S_energy

        integer                  :: ngll, i, j, k, ind
        real(fpp)                :: DXX, DXY, DXZ
        real(fpp)                :: DYX, DYY, DYZ
        real(fpp)                :: DZX, DZY, DZZ
        real(fpp)                :: xmu, xlambda, xkappa
        real(fpp)                :: onemSbeta, onemPbeta
        real(fpp)                :: elem_energy_P, elem_energy_S
        real(fpp)                :: xeps_vol
        real, dimension(0:2,0:2) :: invgrad_ijk
        !
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)


        ngll = dom%ngll

        if(.not. allocated(S_energy)) allocate(S_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if(.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        S_energy = 0.0d0
        P_energy = 0.0d0 !TODO

    end subroutine get_fluid_dom_elem_energy

    subroutine init_material_properties_fluid(dom, lnum, mat, density, lambda)
        type(domain_fluid), intent(inout) :: dom
        integer, intent(in) :: lnum
        type (subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: density
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: lambda
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        dom%IDensity_(:,:,:,bnum,ee) = 1d0/density
        dom%Lambda_ (:,:,:,bnum,ee) = lambda
    end subroutine init_material_properties_fluid

    subroutine init_local_mass_fluid(dom,specel,i,j,k,ind,Whei)
        type(domain_fluid), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer, intent(in) :: i,j,k,ind
        real(kind=fpp), intent(in) :: Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ! Fluid : inertial term ponderation by the inverse of the bulk modulus

        specel%MassMat(i,j,k) = Whei*dom%Jacob_(i,j,k,bnum,ee)/dom%Lambda_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_fluid

    subroutine forces_int_fluid(dom, champs1, bnum)
        use m_calcul_forces_fluid
        type(domain_fluid), intent (INOUT) :: dom
        type(champsfluid), intent(inout) :: champs1
        integer, intent(in) :: bnum
        !
        integer :: ngll,i,j,k,e,ee,idx
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1, 0:dom%ngll-1, 0:dom%ngll-1) :: Fo_Fl,Phi
        real(fpp) :: val

        ngll = dom%ngll

        ! d(rho*Phi)_dX
        ! d(rho*Phi)_dY
        ! d(rho*Phi)_dZ
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        Phi(ee,i,j,k) = champs1%Phi(idx)
                        Fo_Fl(ee,i,j,k) = 0d0
                    enddo
                enddo
            enddo
        enddo

        ! internal forces
        call calcul_forces_fluid(dom,dom%ngll,bnum,Fo_Fl,Phi)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        e = bnum*VCHUNK+ee
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        val = champs1%ForcesFl(idx)
                        val = val - Fo_Fl(ee,i,j,k)
                        champs1%ForcesFl(idx) = val
                    enddo
                enddo
            enddo
        enddo
    end subroutine forces_int_fluid

    subroutine newmark_predictor_fluid(dom)
        type(domain_fluid), intent (INOUT) :: dom
        !
        dom%champs1%VelPhi = dom%champs0%VelPhi
        dom%champs1%Phi    = dom%champs0%Phi
        dom%champs1%ForcesFl = 0d0
    end subroutine newmark_predictor_fluid

    subroutine newmark_corrector_fluid(dom, dt)
        type(domain_fluid), intent (INOUT) :: dom
        double precision :: dt
        !
        integer  :: n,  indpml

        dom%champs0%ForcesFl = dom%champs1%ForcesFl * dom%MassMat
        dom%champs0%VelPhi = (dom%champs0%VelPhi + dt * dom%champs0%ForcesFl)
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs0%VelPhi(indpml) = 0.
        enddo
        dom%champs0%Phi = dom%champs0%Phi + dt * dom%champs0%VelPhi
    end subroutine newmark_corrector_fluid

    function fluid_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        type(domain_fluid), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        Pspeed = sqrt(dom%Lambda_(i,j,k,bnum,ee)*dom%IDensity_(i,j,k,bnum,ee))
    end function fluid_Pspeed
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
