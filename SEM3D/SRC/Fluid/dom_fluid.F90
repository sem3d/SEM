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
#include "gllopt.h"

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
        integer :: nbelem, ngll, nblocks, i, ntemps
        !

        ngll   = dom%ngll
        nbelem = dom%nbelem
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
#if defined(OPENACC) || TEST_FLUID_ACC==1
            ntemps = nblocks
#else
            ntemps = 1
#endif
            allocate(dom%ForcesFl (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1,0:ntemps-1))
            allocate(dom%Phi      (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1,0:ntemps-1))
        end if
        ! Allocation et initialisation de champs0 et champs1 pour les fluides
        if (dom%nglltot /= 0) then
            do i=0,Tdomain%TimeD%nsubsteps
                call allocate_champs_fluid(dom, i)
            end do
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - fluid domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_fluid

    subroutine deallocate_dom_fluid (dom)
        implicit none
        type(domain_fluid), intent (INOUT) :: dom
        !
        integer :: i
        if(allocated(dom%m_IDensity)) deallocate(dom%m_IDensity)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )

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

    subroutine get_fluid_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev, dUdX)
        !$acc routine worker
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
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: S_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: sig_dev
        real(fpp) :: DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,divU!
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        logical :: flag_gradU
        integer :: ngll, i, j, k, ind
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        flag_gradU = (out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS) + &
            out_variables(OUT_DUDX) + &
            out_variables(OUT_EPS_VOL) + &
            out_variables(OUT_EPS_DEV) + &
            out_variables(OUT_STRESS_DEV)) /= 0

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
#ifdef CPML
            call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                 dom%IDensity_(:,:,:,bnum,ee),vphi,fieldV)
#else
            call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                 dom%IDensity_(:,:,:,bnum,ee),phi,fieldV)
#endif
        end if

        if (out_variables(OUT_ACCEL) == 1) then
#ifdef CPML
            call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                 dom%IDensity_(:,:,:,bnum,ee),-fieldP,fieldA)
#else
            call fluid_velocity(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),&
                 dom%IDensity_(:,:,:,bnum,ee),vphi,fieldA)
#endif
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
            P_energy(:,:,:) = 0. !TODO
        end if

        if (out_variables(OUT_ENERGYS) == 1) then
            S_energy(:,:,:) = 0.
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


    subroutine get_fluid_dom_elem_energy(dom, lnum, P_energy, S_energy)
        !$acc routine worker
        use deriv3d
        implicit none
        !
        type(domain_fluid), intent(inout)          :: dom
        integer, intent(in)                        :: lnum
        real(fpp), dimension(:,:,:), allocatable, intent(inout) :: P_energy, S_energy

        integer                  :: ngll
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

    subroutine init_domain_fluid(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_fluid), intent(inout) :: dom

        dom%dt = Tdomain%TimeD%dtmin
    end subroutine init_domain_fluid

    subroutine start_domain_fluid(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_fluid), intent(inout) :: dom
        !
        integer :: i,ns

        !$acc  enter data copyin(dom, dom%champs) &
        !$acc&            copyin(dom%MassMat, dom%dirich) &
        !$acc&            copyin(dom%Phi, dom%ForcesFl)  &
        !$acc&            copyin(dom%m_IDensity, dom%gllw, dom%hprime)&
        !$acc&            copyin(dom%m_Idom, dom%m_Jacob, dom%m_InvGrad) &
        !$acc&
        do i = 0,1
            !$acc enter data  copyin(dom%champs(i)%Phi, dom%champs(i)%VelPhi, dom%champs(i)%ForcesFl)
        end do

        do ns = 0, Tdomain%n_source-1
            if(Tdomain%rank == Tdomain%sSource(ns)%proc)then
                if(Tdomain%sSource(ns)%i_type_source == 3) then
                    !$acc enter data      copyin(Tdomain%sSource(ns)%ExtForce)
                end if
            end if
        end do
    end subroutine start_domain_fluid

    subroutine stop_domain_fluid(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_fluid), intent(inout) :: dom
        !
        integer :: i, ns

        !$acc  exit data delete(dom, dom%champs) &
        !$acc&           delete(dom%MassMat, dom%dirich) &
        !$acc&           delete(dom%Phi, dom%ForcesFl)  &
        !$acc&           delete(dom%m_IDensity, dom%gllw, dom%hprime)&
        !$acc&           delete(dom%m_Idom, dom%m_Jacob, dom%m_InvGrad) &
        !$acc&
        do i = 0,1
            !$acc exit data delete(dom%champs(i)%Phi, dom%champs(i)%VelPhi, dom%champs(i)%ForcesFl)
        end do
        do ns = 0, Tdomain%n_source-1
            if(Tdomain%rank == Tdomain%sSource(ns)%proc)then
                if(Tdomain%sSource(ns)%i_type_source == 3) then
                    !$acc exit  data      delete(Tdomain%sSource(ns)%ExtForce)
                end if
            end if
        end do

    end subroutine stop_domain_fluid

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
        real(fpp), intent(in) :: Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ! Fluid : inertial term ponderation by the inverse of the bulk modulus

        specel%MassMat(i,j,k) = Whei*dom%Jacob_(i,j,k,bnum,ee)/dom%Lambda_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_fluid


    subroutine forces_int_fluid_mainloop(dom, i0, i1, m_dump, m_load, m_expl, m_recalc)
        type(domain_fluid), intent (INOUT) :: dom
        integer, intent(IN) :: i0, i1
        logical, intent(IN) :: m_dump, m_load, m_expl, m_recalc
        !
        if (m_dump .or. m_load) then
            call dispatch_forces_int_mirror_fl(dom, dom%champs(i0), dom%champs(i1), m_dump, m_expl, m_recalc)
        else
            call dispatch_forces_int_fl(dom, dom%champs(i0), dom%champs(i1))
        endif
    end subroutine forces_int_fluid_mainloop

    subroutine dispatch_forces_int_fl(dom, var, dvdt)
        use m_calcul_forces_fluid
        type(domain_fluid), intent (INOUT) :: dom
        type(champsfluid), intent(in) :: var
        type(champsfluid), intent(inout) :: dvdt
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_fl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_5(calcul_forces_fl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_6(calcul_forces_fl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_7(calcul_forces_fl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_8(calcul_forces_fl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_9(calcul_forces_fl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_N(calcul_forces_fl,,(dom,dom%ngll,var,dvdt))
        end select
    end subroutine dispatch_forces_int_fl

    subroutine dispatch_forces_int_mirror_fl(dom, var, dvdt, m_dump, m_expl, m_recalc)
        use m_calcul_forces_fluid
        type(domain_fluid), intent (INOUT) :: dom
        type(champsfluid), intent(in) :: var
        type(champsfluid), intent(inout) :: dvdt
        logical, intent(IN) :: m_dump, m_expl, m_recalc
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_fl,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_5(calcul_forces_fl,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_6(calcul_forces_fl,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_7(calcul_forces_fl,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_8(calcul_forces_fl,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_9(calcul_forces_fl,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_N(calcul_forces_fl,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
        end select
    end subroutine dispatch_forces_int_mirror_fl

!    subroutine forces_int_fluid(dom, field, bnum)
!        use m_calcul_forces_fluid
!        type(domain_fluid), intent (INOUT) :: dom
!        type(champsfluid), intent(inout) :: field
!        integer, intent(in) :: bnum
!        !
!        integer :: ngll,i,j,k,e,ee,idx
!        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1, 0:dom%ngll-1, 0:dom%ngll-1) :: Fo_Fl,Phi
!        real(fpp) :: val
!
!        ngll = dom%ngll
!
!        ! d(rho*Phi)_dX
!        ! d(rho*Phi)_dY
!        ! d(rho*Phi)_dZ
!        do k = 0,ngll-1
!            do j = 0,ngll-1
!                do i = 0,ngll-1
!                    do ee = 0, VCHUNK-1
!                        idx = dom%Idom_(i,j,k,bnum,ee)
!                        Phi(ee,i,j,k) = field%Phi(idx)
!                        Fo_Fl(ee,i,j,k) = 0d0
!                    enddo
!                enddo
!            enddo
!        enddo
!
!        ! internal forces
!        call calcul_forces_fluid(dom,dom%ngll,bnum,Fo_Fl,Phi)
!
!        do k = 0,ngll-1
!            do j = 0,ngll-1
!                do i = 0,ngll-1
!                    do ee = 0, VCHUNK-1
!                        e = bnum*VCHUNK+ee
!                        idx = dom%Idom_(i,j,k,bnum,ee)
!                        val = field%ForcesFl(idx)
!                        val = val - Fo_Fl(ee,i,j,k)
!                        field%ForcesFl(idx) = val
!                    enddo
!                enddo
!            enddo
!        enddo
!    end subroutine forces_int_fluid

!    subroutine forces_int_fluid_mirror_dump(dom,field,bnum)
!        use m_calcul_forces_fluid
!
!        type(domain_fluid),intent(inout) :: dom
!        type(champsfluid),intent(inout) :: field
!        integer,intent(in) :: bnum
!        integer :: lnum,ngll,i,j,k,ee,idx,idx_m
!        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fo_Fl
!        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Phi
!        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: VelPhi
!
!        ngll = dom%ngll
!        lnum = bnum*VCHUNK
!
!        do k = 0,ngll-1
!            do j = 0,ngll-1
!                do i = 0,ngll-1
!                    do ee = 0, VCHUNK-1
!                        idx = dom%Idom_(i,j,k,bnum,ee)
!                        idx_m = dom%mirror_fl%map(lnum+ee,i,j,k)
!                        Phi(ee,i,j,k) = field%Phi(idx)
!                        VelPhi(ee,i,j,k) = field%VelPhi(idx)
!                        if (idx_m>0) then
!                            dom%mirror_fl%fields(1,idx_m) = Phi(ee,i,j,k)
!                            dom%mirror_fl%fields(3,idx_m) = VelPhi(ee,i,j,k)
!                        end if
!                    enddo
!                enddo
!            enddo
!        enddo
!
!        Fo_Fl = 0.0_fpp
!        call calcul_forces_fluid(dom,dom%ngll,bnum,Fo_Fl,Phi)
!
!        do k = 0,ngll-1
!            do j = 0,ngll-1
!                do i = 0,ngll-1
!                    do ee = 0, VCHUNK-1
!                        idx = dom%Idom_(i,j,k,bnum,ee)
!                        idx_m = dom%mirror_fl%map(lnum+ee,i,j,k)
!                        if (idx_m>0) dom%mirror_fl%fields(2,idx_m) = Fo_Fl(ee,i,j,k)
!                        field%ForcesFl(idx) = field%ForcesFl(idx)-Fo_Fl(ee,i,j,k)
!                    enddo
!                enddo
!            enddo
!        enddo
!
!    end subroutine forces_int_fluid_mirror_dump

!    subroutine forces_int_fluid_mirror_load(dom,field,bnum)
!        use m_calcul_forces_fluid
!
!        type(domain_fluid),intent(inout) :: dom
!        type(champsfluid),intent(inout) :: field
!        integer,intent(in) :: bnum
!        integer :: lnum,ngll,i,j,k,ee,idx,idx_m
!        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Phi
!        real(fpp),dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fo_Fl
!
!        ngll = dom%ngll
!        lnum = bnum*VCHUNK
!
!        do k = 0,ngll-1
!            do j = 0,ngll-1
!                do i = 0,ngll-1
!                    do ee = 0, VCHUNK-1
!                        idx = dom%Idom_(i,j,k,bnum,ee)
!                        idx_m = dom%mirror_fl%map(lnum+ee,i,j,k)
!                        Phi(ee,i,j,k) = field%Phi(idx)
!                        if (idx_m>0) Phi(ee,i,j,k) = Phi(ee,i,j,k)+dom%mirror_fl%fields(1,idx_m) &
!                            *dom%mirror_fl%winfunc(idx_m)
!                    enddo
!                enddo
!            enddo
!        enddo
!
!        Fo_Fl = 0.0_fpp
!        call calcul_forces_fluid(dom,dom%ngll,bnum,Fo_Fl,Phi)
!
!        do k = 0,ngll-1
!            do j = 0,ngll-1
!                do i = 0,ngll-1
!                    do ee = 0, VCHUNK-1
!                        idx = dom%Idom_(i,j,k,bnum,ee)
!                        idx_m = dom%mirror_fl%map(lnum+ee,i,j,k)
!                        if (idx_m>0) Fo_Fl(ee,i,j,k) = Fo_Fl(ee,i,j,k)-dom%mirror_fl%fields(2,idx_m) &
!                            *dom%mirror_fl%winfunc(idx_m)
!                        field%ForcesFl(idx) = field%ForcesFl(idx)-Fo_Fl(ee,i,j,k)
!                    enddo
!                enddo
!            enddo
!        enddo
!
!    end subroutine forces_int_fluid_mirror_load

    subroutine newmark_predictor_fluid(dom, f0, f1)
        type(domain_fluid), intent (INOUT) :: dom
        integer, intent(in) :: f0, f1
        !
        integer :: i
!!         !$acc parallel loop async(1) present(dom,dom%champs) &
!!         !$acc&   present(dom%champs(f0)%Phi,dom%champs(f0)%VelPhi) &
!!         !$acc&   present(dom%champs(f1)%Phi,dom%champs(f1)%VelPhi,dom%champs(f1)%ForcesFl)
!!         do i=0,dom%nglltot-1
!!             !dom%champs(f1)%VelPhi(i)   = dom%champs(f0)%VelPhi(i)
!!             !dom%champs(f1)%Phi(i)      = dom%champs(f0)%Phi(i)
!!             dom%champs(f1)%ForcesFl(i) = 0d0
!!         end do
!!         !$acc end parallel
        !$acc kernels async(1)
        dom%champs(f1)%ForcesFl = 0d0
        !$acc end kernels

        !! BUG!!
!!        !$acc kernels async(1)
!!        dom%champs(f1)%VelPhi   = dom%champs(f0)%VelPhi
!!        dom%champs(f1)%Phi      = dom%champs(f0)%Phi
!!        dom%champs(f1)%ForcesFl = 0d0
!!        !$acc end kernels

        !! OK
!!        !$acc kernels async(1)
!!        dom%champs(f1)%VelPhi   = dom%champs(f0)%VelPhi
!!        !$acc end kernels
!!        !$acc kernels async(1)
!!        dom%champs(f1)%Phi      = dom%champs(f0)%Phi
!!        !$acc end kernels
!!        !$acc kernels async(1)
!!        dom%champs(f1)%ForcesFl = 0d0
!!        !$acc end kernels

        !! OK
!!        !$acc kernels async(1) present(dom%champs(f1)%ForcesFl)
!!        dom%champs(f1)%VelPhi(0:dom%nglltot-1)   = dom%champs(f0)%VelPhi(0:dom%nglltot-1)
!!        dom%champs(f1)%Phi(0:dom%nglltot-1)      = dom%champs(f0)%Phi(0:dom%nglltot-1)
!!        dom%champs(f1)%ForcesFl(0:dom%nglltot-1) = 0d0
!!        !$acc end kernels
    end subroutine newmark_predictor_fluid

    subroutine newmark_corrector_fluid(dom, dt, f0, f1)
        type(domain_fluid), intent (INOUT) :: dom
        real(fpp), intent(in) :: dt
        integer, intent(in) :: f0, f1
        !
        integer  :: n,  indpml, count, nglltot
        !
        real(fpp) :: val
        count = dom%n_dirich
        nglltot = dom%nglltot
        !$acc kernels async(1)  present(dom,dom%champs, dom%massmat) &
        !$acc&     present(dom%champs(f1)%ForcesFl) &
        !$acc&     present(dom%champs(f0)%VelPhi,dom%champs(f1)%Phi)
        dom%champs(f0)%VelPhi = dom%champs(f0)%VelPhi + dt * dom%champs(f1)%ForcesFl * dom%MassMat
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs(f0)%VelPhi(indpml) = 0.
        enddo
        dom%champs(f0)%Phi = dom%champs(f0)%Phi + dt * dom%champs(f0)%VelPhi
        !$acc end kernels

!        !$acc update host(dom%champs(f1)%ForcesFl,dom%champs(f0)%ForcesFl) wait(1)
!        write(*,*) "nglltot:", nglltot
!        write(*,"(A,E16.9,E16.9,E16.9)") "src:", dom%MassMat(35850),dom%champs(f0)%ForcesFl(35850),dom%champs(f1)%ForcesFl(35850)

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

    subroutine apply_source_fluid(src, dom, i1, ft, lnum)
        use sdomain
        implicit none
        type(domain_fluid),intent(inout) :: dom
        type(Source),intent(inout) :: src 
        real(fpp), intent(in) :: ft
        integer, intent(in) :: i1
        integer, intent(in) :: lnum
        integer :: bnum, ee
        integer :: ngll, i, j, k, idx
        real(kind=fpp) :: val

        ngll = dom%ngll
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        !$acc parallel loop gang vector collapse(3) async(1)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    idx = dom%Idom_(i,j,k,bnum,ee)
                    val = dom%champs(i1)%ForcesFl(idx)
                    val = val + ft*src%ExtForce(i,j,k,0)
                    dom%champs(i1)%ForcesFl(idx) = val
                enddo
            enddo
        enddo
!!        write(*,*) "source:"
!!        !$acc update host(dom%champs(i1)%ForcesFl) wait(1)
!!        i=2
!!        j=2
!!        k=2
!!        idx = dom%Idom_(i,j,k,bnum,ee)
!!        val = dom%champs(i1)%ForcesFl(idx)
!!        !val = src%ExtForce(i,j,k,0)
!!        write(*,"(A,I2,I2,I2,A,I5,A,E16.9)") "src:",i,j,k,"idx:",idx,":",val

    end subroutine apply_source_fluid

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
