!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!! Also carries anisotropic-density materials (Capdeville & Cance 2015):
!!     (1/kappa) phi_tt = d_i( rho^{-1}_ij d_j phi ) ,   v_i = rho^{-1}_ij d_j phi
!! selected per-domain by dom%aniso (mirrors dom_solid's Cij_/aniso split), generalising
!! the isotropic scalar IDensity_=1/rho to the symmetric tensor m_IDensTensor. See
!! calcul_forces_fluid_aniso.inc for the kernel and DOC/source/anisotropic_fluid.rst for
!! the physics writeup.

module dom_fluid
    use sdomain
    use constants
    use champs_fluid
    use selement
    use ssubdomains
    implicit none
#include "index.h"

contains

    subroutine allocate_champs_fluid(dom, i)
        type(domain_fluid), intent (INOUT) :: dom
        integer, intent(in) :: i
        allocate(dom%champs(i)%ForcesFl(0:dom%nglltot))
        allocate(dom%champs(i)%Phi     (0:dom%nglltot))
        allocate(dom%champs(i)%VelPhi  (0:dom%nglltot))

        dom%champs(i)%ForcesFl = 0d0
        dom%champs(i)%Phi = 0d0
        dom%champs(i)%VelPhi = 0d0
    end subroutine allocate_champs_fluid

    subroutine allocate_dom_fluid (Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_fluid), intent (INOUT) :: dom
        !
        integer :: nbelem, ngll, nblocks, i
        !

        ngll   = dom%ngll
        nbelem = dom%nbelem
        dom%aniso = Tdomain%aniso
!        write(*,*) "DOM_FLUID ngll   = ", ngll
!        write(*,*) "DOM_FLUID nbelem = ", nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call init_dombase(dom)

        ! Mirror
        !!! GB dom%use_mirror = Tdomain%use_mirror
        dom%mirror_type = Tdomain%mirror_type

        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)
            nblocks = dom%nblocks
            allocate(dom%IDensity_(0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            if (dom%aniso) then
                ! Anisotropic-density materials (Capdeville & Cance 2015): rho^{-1}_ij
                ! tensor, mirrors dom_solid%Cij_ -- allocated model-wide whenever any
                ! material (solid or fluid) is aniso. Isotropic fluid elements in the
                ! same run still fill it with the diagonal special case (see
                ! init_material_properties_fluid).
                allocate(dom%m_IDensTensor(IND_DIJKE(0:5,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1,0:VCHUNK-1)))
            end if
        end if
        ! Allocation et initialisation de champs0 et champs1 pour les fluides
        if (dom%nglltot /= 0) then
            do i=0,Tdomain%TimeD%nsubsteps
                call allocate_champs_fluid(dom, i)
            end do
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - fluid domain (rank 0, local) : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_fluid

    subroutine deallocate_dom_fluid (dom)
        implicit none
        type(domain_fluid), intent (INOUT) :: dom
        !
        integer :: i
        if(allocated(dom%m_IDensity)) deallocate(dom%m_IDensity)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_IDensTensor)) deallocate(dom%m_IDensTensor)

        do i=0,1
            if(allocated(dom%champs(i)%ForcesFl)) deallocate(dom%champs(i)%ForcesFl)
            if(allocated(dom%champs(i)%Phi     )) deallocate(dom%champs(i)%Phi     )
            if(allocated(dom%champs(i)%VelPhi  )) deallocate(dom%champs(i)%VelPhi  )
        end do
        call deallocate_dombase(dom)
    end subroutine deallocate_dom_fluid

    subroutine fluid_velocity(ngll,hprime,InvGrad,idensity,phi,veloc)
        ! gives the physical particle velocity in the fluid = 1/dens grad(dens.Phi)
        implicit none
        integer, intent(in)  :: ngll
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:2), intent(in) :: InvGrad
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: idensity,phi
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(out) :: Veloc
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: dphi_dx,dphi_dy,dphi_dz

        ! physical gradient
        call physical_part_deriv(ngll,hprime,InvGrad,phi,dphi_dx,dphi_dy,dphi_dz)
        Veloc(:,:,:,0) = dphi_dx(:,:,:) * idensity(:,:,:)
        Veloc(:,:,:,1) = dphi_dy(:,:,:) * idensity(:,:,:)
        Veloc(:,:,:,2) = dphi_dz(:,:,:) * idensity(:,:,:)
    end subroutine fluid_velocity

    subroutine fluid_aniso_velocity(ngll,hprime,InvGrad,IDensTensor,phi,veloc)
        ! Physical particle velocity in the anisotropic-density fluid: v_i = rho^{-1}_ij d_j phi
        implicit none
        integer, intent(in) :: ngll
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:2), intent(in) :: InvGrad
        real(fpp), dimension(0:5,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: IDensTensor
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: phi
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(out) :: Veloc
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: dphi_dx,dphi_dy,dphi_dz
        integer :: i,j,k
        real(fpp) :: B11,B22,B33,B12,B13,B23,gx,gy,gz

        call physical_part_deriv(ngll,hprime,InvGrad,phi,dphi_dx,dphi_dy,dphi_dz)
        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    B11 = IDensTensor(0,i,j,k); B22 = IDensTensor(1,i,j,k); B33 = IDensTensor(2,i,j,k)
                    B12 = IDensTensor(3,i,j,k); B13 = IDensTensor(4,i,j,k); B23 = IDensTensor(5,i,j,k)
                    gx = dphi_dx(i,j,k); gy = dphi_dy(i,j,k); gz = dphi_dz(i,j,k)
                    Veloc(i,j,k,0) = B11*gx + B12*gy + B13*gz
                    Veloc(i,j,k,1) = B12*gx + B22*gy + B23*gz
                    Veloc(i,j,k,2) = B13*gx + B23*gy + B33*gz
                enddo
            enddo
        enddo
    end subroutine fluid_aniso_velocity

    subroutine get_fluid_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, K_energy, D_energy, eps_vol, eps_dev, sig_dev, dUdX)
        use deriv3d
        implicit none
        !
        type(domain_fluid), intent(inout)          :: dom
        integer, dimension(0:), intent(in)         :: out_variables
        integer, intent(in)                        :: lnum
        real(fpp), dimension(:,:,:), allocatable   :: phi
        real(fpp), dimension(:,:,:), allocatable   :: vphi
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldU
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldV
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldA
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:8) :: dUdX
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: fieldP
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: P_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: K_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: D_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: sig_dev
        real(fpp) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,divU!
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: gx,gy,gz
        logical :: flag_gradU
        integer :: ngll, i, j, k, ind
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        flag_gradU = (out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYK) + &
            out_variables(OUT_DUDX) + &
            out_variables(OUT_EPS_VOL) + &
            out_variables(OUT_EPS_DEV) + &
            out_variables(OUT_STRESS_DEV) + &
            out_variables(OUT_ENERGYD)) /= 0

        ngll = dom%ngll
        allocate(phi(0:ngll-1,0:ngll-1,0:ngll-1))
        allocate(vphi(0:ngll-1,0:ngll-1,0:ngll-1))

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    phi(i,j,k) = dom%champs(0)%Phi(ind)
                enddo
            enddo
        enddo
        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    vphi(i,j,k) = dom%champs(0)%VelPhi(ind)
                enddo
            enddo
        enddo

        if (out_variables(OUT_PRESSION) == 1) then
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
#ifdef CPML
                        fieldP(i,j,k) = -dom%champs(0)%ForcesFl(ind)
#else
                        fieldP(i,j,k) = -dom%champs(0)%VelPhi(ind)
#endif
                    enddo
                enddo
            enddo
        end if

        if (out_variables(OUT_DEPLA) == 1) then
#ifdef CPML
            call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                 dom%IDensity_(:,:,:,bnum,ee),phi,fieldU)
#else
            fieldU(:,:,:,:) = 0.
#endif
        end if

        if (out_variables(OUT_VITESSE) == 1.or.flag_gradU) then
            if (dom%aniso) then
#ifdef CPML
                call fluid_aniso_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                     dom%m_IDensTensor(IND_DIJKE(0:5,:,:,:,bnum,ee)),vphi,fieldV)
#else
                call fluid_aniso_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                     dom%m_IDensTensor(IND_DIJKE(0:5,:,:,:,bnum,ee)),phi,fieldV)
#endif
            else
#ifdef CPML
                call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                     dom%IDensity_(:,:,:,bnum,ee),vphi,fieldV)
#else
                call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                     dom%IDensity_(:,:,:,bnum,ee),phi,fieldV)
#endif
            end if
        end if

        if (out_variables(OUT_ACCEL) == 1) then
            if (dom%aniso) then
#ifdef CPML
                call fluid_aniso_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                     dom%m_IDensTensor(IND_DIJKE(0:5,:,:,:,bnum,ee)),-fieldP,fieldA)
#else
                call fluid_aniso_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                     dom%m_IDensTensor(IND_DIJKE(0:5,:,:,:,bnum,ee)),vphi,fieldA)
#endif
            else
#ifdef CPML
                call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                     dom%IDensity_(:,:,:,bnum,ee),-fieldP,fieldA)
#else
                call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                     dom%IDensity_(:,:,:,bnum,ee),vphi,fieldA)
#endif
            end if
        end if

        if (flag_gradU) then
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        invgrad_ijk = dom%InvGrad_(:,:,i,j,k,bnum,ee) ! cache for performance
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                            invgrad_ijk,fieldV(:,:,:,0),DXX,DYX,DZX)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                            invgrad_ijk,fieldV(:,:,:,1),DXY,DYY,DZY)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                            invgrad_ijk,fieldV(:,:,:,2),DXZ,DYZ,DZZ)
                        divU = DXX+DYY+DZZ

                        if (out_variables(OUT_DUDX) == 1) then
                            dUdX(i,j,k,0) = DXX
                            dUdX(i,j,k,1) = DXY
                            dUdX(i,j,k,2) = DXZ
                            dUdX(i,j,k,3) = DYX
                            dUdX(i,j,k,4) = DYY
                            dUdX(i,j,k,5) = DYZ
                            dUdX(i,j,k,6) = DZX
                            dUdX(i,j,k,7) = DZY
                            dUdX(i,j,k,8) = DZZ
                        end if

                        if (out_variables(OUT_EPS_VOL) == 1) then
                            eps_vol(:,:,:) = divU
                        end if
                    end do
                end do
            end do
        end if
        if (out_variables(OUT_ENERGYP) == 1) then
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
#ifdef CPML
                        fieldP(i,j,k) = -dom%champs(0)%ForcesFl(ind)
                        P_energy(i,j,k) = 0.
#else
                        fieldP(i,j,k) = -dom%champs(0)%VelPhi(ind)
                        P_energy(i,j,k) = 0.5*fieldP(i,j,k)*fieldP(i,j,k)/dom%Lambda_(i,j,k,bnum,ee)
#endif
                   enddo
               enddo
            enddo
        end if

        if (out_variables(OUT_ENERGYK) == 1) then
            if (dom%aniso) then
                ! K_energy = 1/2 (grad phi).(rho^{-1} grad phi) = 1/2 v.(grad phi)
                call fluid_aniso_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                                 dom%m_IDensTensor(IND_DIJKE(0:5,:,:,:,bnum,ee)),phi,fieldV)
                call physical_part_deriv(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),phi,gx,gy,gz)
                do k=0,ngll-1
                    do j=0,ngll-1
                        do i=0,ngll-1
                            K_energy(i,j,k) = 0.5d0*(fieldV(i,j,k,0)*gx(i,j,k) + &
                                                     fieldV(i,j,k,1)*gy(i,j,k) + &
                                                     fieldV(i,j,k,2)*gz(i,j,k))
                        enddo
                    enddo
                enddo
            else
                call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                                 dom%IDensity_(:,:,:,bnum,ee),phi,fieldV)
                do k=0,ngll-1
                    do j=0,ngll-1
                        do i=0,ngll-1
                            K_energy(i,j,k) = 0.5*(fieldV(i,j,k,0)**2+fieldV(i,j,k,1)**2+fieldV(i,j,k,2)**2)/dom%IDensity_(i,j,k,bnum,ee)
                        enddo
                    enddo
                enddo
            end if
        end if

        if (out_variables(OUT_ENERGYD) == 1) then
            D_energy(:,:,:,:) = 0
        end if

        if (out_variables(OUT_EPS_DEV) == 1) then
            eps_dev(:,:,:,:) = 0.
        end if

        if (out_variables(OUT_STRESS_DEV) == 1) then
            sig_dev(:,:,:,:) = 0.
        end if
        deallocate(phi)        
        deallocate(vphi)
    end subroutine get_fluid_dom_var


    subroutine get_fluid_dom_elem_energy(dom, lnum, P_energy, K_energy, D_energy)

        use deriv3d
        implicit none
        !
        type(domain_fluid), intent(inout)          :: dom
        integer, intent(in)                        :: lnum
        real(fpp), dimension(:,:,:), allocatable, intent(inout) :: P_energy, K_energy
        real(fpp), dimension(:,:,:,:), allocatable, intent(inout) :: D_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldV
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: fieldP
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: gx,gy,gz
        real(fpp), dimension(:,:,:), allocatable   :: phi
        integer                  :: ngll, i, j, k, ind
        !
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)


        ngll = dom%ngll
        allocate(phi(0:ngll-1,0:ngll-1,0:ngll-1))

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    phi(i,j,k) = dom%champs(0)%Phi(ind)
#ifdef CPML
                    fieldP(i,j,k) = -dom%champs(0)%ForcesFl(ind)
#else
                    fieldP(i,j,k) = -dom%champs(0)%VelPhi(ind)
#endif                                                                                                                                                                       
                enddo
            enddo
        enddo

        !Allocation
        if(.not. allocated(K_energy)) allocate(K_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if(.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if(.not. allocated(P_energy)) allocate(D_energy(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        K_energy = -1       
        P_energy = -1
        D_energy = 0

        if (dom%aniso) then
            call fluid_aniso_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                                dom%m_IDensTensor(IND_DIJKE(0:5,:,:,:,bnum,ee)),phi,fieldV)
            call physical_part_deriv(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),phi,gx,gy,gz)
            ! Then, get the energies.
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        P_energy(i,j,k) = 0.5*fieldP(i,j,k)*fieldP(i,j,k)/dom%Lambda_(i,j,k,bnum,ee)
                        ! K_energy = 1/2 (grad phi).(rho^{-1} grad phi) = 1/2 v.(grad phi)
                        K_energy(i,j,k) = 0.5d0*(fieldV(i,j,k,0)*gx(i,j,k) + &
                                                 fieldV(i,j,k,1)*gy(i,j,k) + &
                                                 fieldV(i,j,k,2)*gz(i,j,k))
                    enddo
                enddo
            enddo
        else
            call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                                dom%IDensity_(:,:,:,bnum,ee),phi,fieldV)
            ! Then, get the energies.
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        P_energy(i,j,k) = 0.5*fieldP(i,j,k)*fieldP(i,j,k)/dom%Lambda_(i,j,k,bnum,ee)
                        K_energy(i,j,k) = 0.5*(fieldV(i,j,k,0)**2+fieldV(i,j,k,1)**2+fieldV(i,j,k,2)**2)/dom%IDensity_(i,j,k,bnum,ee)
                    enddo
                enddo
            enddo
        end if

        deallocate(phi)
    end subroutine get_fluid_dom_elem_energy

    subroutine init_domain_fluid(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_fluid), intent(inout) :: dom

        dom%dt = Tdomain%TimeD%dtmin
    end subroutine init_domain_fluid

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

    subroutine init_material_properties_fluid_aniso(dom, lnum, mat, IDens, lambda)
        ! IDens(0:5) = inverse-density tensor rho^{-1}_ij at each GLL point (11,22,33,12,13,23)
        ! lambda     = scalar bulk modulus kappa
        type(domain_fluid), intent(inout) :: dom
        integer, intent(in) :: lnum
        type(subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:5,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: IDens
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: lambda
        !
        integer :: bnum, ee, m
        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)

        do m = 0, 5
            dom%m_IDensTensor(IND_DIJKE(m,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,bnum,ee)) = IDens(m,:,:,:)
        end do
        dom%Lambda_(:,:,:,bnum,ee) = lambda
    end subroutine init_material_properties_fluid_aniso

    subroutine init_local_mass_fluid(dom,specel,i,j,k,ind,Whei)
        type(domain_fluid), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer, intent(in) :: i,j,k,ind
        real(fpp), intent(in) :: Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ! Fluid : inertial term ponderation by the inverse of the bulk modulus

        specel%MassMat(i,j,k) = Whei*dom%Jacob_(i,j,k,bnum,ee)/dom%Lambda_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_fluid

    subroutine forces_int_fluid(dom, field, bnum)
        use m_calcul_forces_fluid
        type(domain_fluid), intent (INOUT) :: dom
        type(champsfluid), intent(inout) :: field
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
                        Phi(ee,i,j,k) = field%Phi(idx)
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
                        val = field%ForcesFl(idx)
                        val = val - Fo_Fl(ee,i,j,k)
                        field%ForcesFl(idx) = val
                    enddo
                enddo
            enddo
        enddo
    end subroutine forces_int_fluid

    subroutine forces_int_fluid_mirror_dump(dom,field,bnum)
        use m_calcul_forces_fluid

        type(domain_fluid),intent(inout) :: dom
        type(champsfluid),intent(inout) :: field
        integer,intent(in) :: bnum
        integer :: lnum,ngll,i,j,k,ee,idx,idx_m
        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fo_Fl
        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Phi
        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: VelPhi

        ngll = dom%ngll
        lnum = bnum*VCHUNK

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        idx_m = dom%mirror_fl%map(lnum+ee,i,j,k)
                        Phi(ee,i,j,k) = field%Phi(idx)
                        VelPhi(ee,i,j,k) = field%VelPhi(idx)
                        if (idx_m>0) then
                            dom%mirror_fl%fields(1,idx_m) = Phi(ee,i,j,k)
                            dom%mirror_fl%fields(3,idx_m) = VelPhi(ee,i,j,k)
                        end if
                    enddo
                enddo
            enddo
        enddo

        Fo_Fl = 0.0_fpp
        call calcul_forces_fluid(dom,dom%ngll,bnum,Fo_Fl,Phi)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        idx_m = dom%mirror_fl%map(lnum+ee,i,j,k)
                        if (idx_m>0) dom%mirror_fl%fields(2,idx_m) = Fo_Fl(ee,i,j,k)
                        field%ForcesFl(idx) = field%ForcesFl(idx)-Fo_Fl(ee,i,j,k)
                    enddo
                enddo
            enddo
        enddo

    end subroutine forces_int_fluid_mirror_dump

    subroutine forces_int_fluid_mirror_load(dom,field,bnum)
        use m_calcul_forces_fluid

        type(domain_fluid),intent(inout) :: dom
        type(champsfluid),intent(inout) :: field
        integer,intent(in) :: bnum
        integer :: lnum,ngll,i,j,k,ee,idx,idx_m
        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Phi
        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fo_Fl

        ngll = dom%ngll
        lnum = bnum*VCHUNK

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        idx_m = dom%mirror_fl%map(lnum+ee,i,j,k)
                        Phi(ee,i,j,k) = field%Phi(idx)
                        if (idx_m>0) Phi(ee,i,j,k) = Phi(ee,i,j,k)+dom%mirror_fl%fields(1,idx_m) &
                            *dom%mirror_fl%winfunc(idx_m)
                    enddo
                enddo
            enddo
        enddo

        Fo_Fl = 0.0_fpp
        call calcul_forces_fluid(dom,dom%ngll,bnum,Fo_Fl,Phi)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        idx_m = dom%mirror_fl%map(lnum+ee,i,j,k)
                        if (idx_m>0) Fo_Fl(ee,i,j,k) = Fo_Fl(ee,i,j,k)-dom%mirror_fl%fields(2,idx_m) &
                            *dom%mirror_fl%winfunc(idx_m)
                        field%ForcesFl(idx) = field%ForcesFl(idx)-Fo_Fl(ee,i,j,k)
                    enddo
                enddo
            enddo
        enddo

    end subroutine forces_int_fluid_mirror_load

    subroutine newmark_predictor_fluid(dom, f0, f1)
        type(domain_fluid), intent (INOUT) :: dom
        integer, intent(in) :: f0, f1
        !
        dom%champs(f1)%VelPhi   = dom%champs(f0)%VelPhi
        dom%champs(f1)%Phi      = dom%champs(f0)%Phi
        dom%champs(f1)%ForcesFl = 0d0
    end subroutine newmark_predictor_fluid

    subroutine newmark_corrector_fluid(dom, dt, f0, f1)
        type(domain_fluid), intent (INOUT) :: dom
        real(fpp), intent(in) :: dt
        integer, intent(in) :: f0, f1
        !
        integer  :: n,  indpml

        dom%champs(f0)%ForcesFl = dom%champs(f1)%ForcesFl * dom%MassMat
        dom%champs(f0)%VelPhi = (dom%champs(f0)%VelPhi + dt * dom%champs(f0)%ForcesFl)
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs(f0)%VelPhi(indpml) = 0.
        enddo
        dom%champs(f0)%Phi = dom%champs(f0)%Phi + dt * dom%champs(f0)%VelPhi
    end subroutine newmark_corrector_fluid

    function fluid_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        type(domain_fluid), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed
        real(fpp) :: B11,B22,B33,B12,B13,B23,maxB
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        if (dom%aniso) then
            ! Upper bound on phase velocity: sqrt(kappa * lambda_max(rho^{-1})), with
            ! lambda_max(rho^{-1}) bounded above by the Gershgorin radius of the tensor.
            B11 = dom%m_IDensTensor(IND_DIJKE(0,i,j,k,bnum,ee))
            B22 = dom%m_IDensTensor(IND_DIJKE(1,i,j,k,bnum,ee))
            B33 = dom%m_IDensTensor(IND_DIJKE(2,i,j,k,bnum,ee))
            B12 = dom%m_IDensTensor(IND_DIJKE(3,i,j,k,bnum,ee))
            B13 = dom%m_IDensTensor(IND_DIJKE(4,i,j,k,bnum,ee))
            B23 = dom%m_IDensTensor(IND_DIJKE(5,i,j,k,bnum,ee))
            maxB = max(abs(B11)+abs(B12)+abs(B13), &
                       abs(B12)+abs(B22)+abs(B23), &
                       abs(B13)+abs(B23)+abs(B33))
            Pspeed = sqrt(maxB * dom%Lambda_(i,j,k,bnum,ee))
        else
            Pspeed = sqrt(dom%Lambda_(i,j,k,bnum,ee)*dom%IDensity_(i,j,k,bnum,ee))
        end if
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
