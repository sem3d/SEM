!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_fluidpml
    use sdomain
    use constants
    use champs_fluidpml
    use selement
    use sdomain
    use ssubdomains
    implicit none
#include "index.h"

contains

    subroutine allocate_dom_fluidpml (Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_fluidpml) :: dom
        !
        integer nbelem, ngll
        !

        nbelem = dom%nbelem
        ngll   = dom%ngll
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call compute_gll_data(ngll, dom%gllc, dom%gllw, dom%hprime, dom%htprime)

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)

            nbelem = CHUNK*((nbelem+CHUNK-1)/CHUNK)
            dom%nbelem_alloc = nbelem

            allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nbelem-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nbelem-1))

            allocate (dom%Jacob_  (        0:ngll-1,0:ngll-1,0:ngll-1,0:nbelem-1))
            allocate (dom%InvGrad_(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nbelem-1))

            allocate(dom%Idom_(0:ngll-1,0:ngll-1,0:ngll-1,0:nbelem-1))
            dom%m_Idom = 0

            if(Tdomain%TimeD%velocity_scheme)then
                allocate(dom%PMLVeloc_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nbelem-1))
                dom%PMLVeloc_(:,:,:,:,:) = 0d0
                allocate(dom%PMLDumpSx_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nbelem-1))
                allocate(dom%PMLDumpSy_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nbelem-1))
                allocate(dom%PMLDumpSz_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nbelem-1))
                dom%PMLDumpSx_(:,:,:,:,:) = 0d0
                dom%PMLDumpSy_(:,:,:,:,:) = 0d0
                dom%PMLDumpSz_(:,:,:,:,:) = 0d0
            endif
        end if

        ! Allocation et initialisation de champs0 pour les PML fluides
        if (dom%nglltot /= 0) then
            allocate(dom%champs1%fpml_Forces(0:dom%nglltot-1,0:2))
            allocate(dom%champs0%fpml_VelPhi(0:dom%nglltot-1,0:2))
            allocate(dom%champs0%fpml_Phi   (0:dom%nglltot-1,0:2))
            allocate(dom%champs1%fpml_VelPhi(0:dom%nglltot-1,0:2))
            allocate(dom%champs0%fpml_DumpV (0:dom%nglltot-1,0:1,0:2))
            dom%champs1%fpml_Forces = 0d0
            dom%champs0%fpml_VelPhi = 0d0
            dom%champs0%fpml_Phi = 0d0
            dom%champs0%fpml_DumpV = 0d0

            ! Allocation de MassMat pour les PML fluides
            allocate(dom%MassMat(0:dom%nglltot-1))
            dom%MassMat = 0d0

            allocate(dom%DumpMass(0:dom%nglltot-1,0:2))
            dom%DumpMass = 0d0
        endif
    end subroutine allocate_dom_fluidpml

    subroutine deallocate_dom_fluidpml (dom)
        implicit none
        type(domain_fluidpml) :: dom

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )

        if(allocated(dom%m_Jacob  )) deallocate(dom%m_Jacob  )
        if(allocated(dom%m_InvGrad)) deallocate(dom%m_InvGrad)

        if(allocated(dom%m_Idom)) deallocate(dom%m_Idom)

        if(allocated(dom%gllc))    deallocate(dom%gllc)
        if(allocated(dom%gllw))    deallocate(dom%gllw)
        if(allocated(dom%hprime))  deallocate(dom%hprime)
        if(allocated(dom%htprime)) deallocate(dom%htprime)

        if(allocated(dom%m_PMLVeloc )) deallocate(dom%m_PMLVeloc )
        if(allocated(dom%m_PMLDumpSx)) deallocate(dom%m_PMLDumpSx)
        if(allocated(dom%m_PMLDumpSy)) deallocate(dom%m_PMLDumpSy)
        if(allocated(dom%m_PMLDumpSz)) deallocate(dom%m_PMLDumpSz)

        if(allocated(dom%champs1%fpml_Forces)) deallocate(dom%champs1%fpml_Forces)
        if(allocated(dom%champs0%fpml_VelPhi)) deallocate(dom%champs0%fpml_VelPhi)
        if(allocated(dom%champs0%fpml_Phi   )) deallocate(dom%champs0%fpml_Phi   )
        if(allocated(dom%champs1%fpml_VelPhi)) deallocate(dom%champs1%fpml_VelPhi)
        if(allocated(dom%champs0%fpml_DumpV )) deallocate(dom%champs0%fpml_DumpV )

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)

        if(allocated(dom%DumpMass)) deallocate(dom%DumpMass)
    end subroutine deallocate_dom_fluidpml

    subroutine get_fluidpml_dom_var(dom, el, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain_fluidpml), intent(inout)       :: dom
        integer, dimension(0:8)                    :: out_variables
        type(element)                              :: el
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        !
        logical :: flag_gradU
        integer :: ngll, i, j, k, ind

        flag_gradU = (out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS) + &
            out_variables(OUT_EPS_VOL) + &
            out_variables(OUT_EPS_DEV) + &
            out_variables(OUT_STRESS_DEV)) /= 0

        ngll = dom%ngll

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,el%lnum)

                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        if(.not. allocated(fieldU)) allocate(fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldU(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_VITESSE) == 1) then
                        if(.not. allocated(fieldV)) allocate(fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldV(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldA(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
                        fieldP(i,j,k) = 0d0
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        if(.not. allocated(eps_vol)) allocate(eps_vol(0:ngll-1,0:ngll-1,0:ngll-1))
                        eps_vol(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        if(.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                        P_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        if(.not. allocated(S_energy)) allocate(S_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                        S_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        if(.not. allocated(eps_dev)) allocate(eps_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
                        eps_dev(i,j,k,:) = 0.
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        if(.not. allocated(sig_dev)) allocate(sig_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
                        sig_dev(i,j,k,:) = 0.
                    end if
                enddo
            enddo
        enddo
    end subroutine get_fluidpml_dom_var

    subroutine init_material_properties_fluidpml(dom, lnum, i, j, k, density, lambda)
        type(domain_fluidpml), intent(inout) :: dom
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
    end subroutine init_material_properties_fluidpml

    subroutine init_local_mass_fluidpml(dom,specel,i,j,k,ind,Whei)
        type(domain_fluidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real Whei

        ! Fluid : inertial term ponderation by the inverse of the bulk modulus

        specel%MassMat(i,j,k) = Whei*dom%Jacob_(i,j,k,specel%lnum)/dom%Lambda_(i,j,k,specel%lnum)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_fluidpml
#if 0
    subroutine forces_int_flu_pml(dom, champs1, lnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        type(champsfluidpml), intent(inout) :: champs1
        integer :: lnum
        !
        integer :: ngll
        integer :: i, j, k, l, ind
        real :: sum_vx, sum_vy, sum_vz, acoeff
        real, dimension(0:2,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: ForcesFl

        ngll = dom%ngll

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
                    sum_vx = 0d0
                    sum_vy = 0d0
                    sum_vz = 0d0
                    do l = 0,ngll-1
                        acoeff = - dom%hprime(i,l)*dom%GLLw(l)*dom%GLLw(j)*dom%GLLw(k)*dom%Jacob_(l,j,k,lnum)
                        sum_vx = sum_vx + acoeff*dom%InvGrad_(0,0,l,j,k,lnum)*dom%PMLVeloc_(l,j,k,0,lnum)
                        sum_vy = sum_vy + acoeff*dom%InvGrad_(1,0,l,j,k,lnum)*dom%PMLVeloc_(l,j,k,1,lnum)
                        sum_vz = sum_vz + acoeff*dom%InvGrad_(2,0,l,j,k,lnum)*dom%PMLVeloc_(l,j,k,2,lnum)
                    end do
                    ForcesFl(0,i,j,k) = sum_vx
                    ForcesFl(1,i,j,k) = sum_vy
                    ForcesFl(2,i,j,k) = sum_vz
                end do
            end do
        end do
        do k = 0,ngll-1
            do l = 0,ngll-1
                do j = 0,ngll-1
                    do i=0,ngll-1
                        acoeff = - dom%hprime(j,l)*dom%GLLw(i)*dom%GLLw(l)*dom%GLLw(k)*dom%Jacob_(i,l,k,lnum)
                        sum_vx = acoeff*dom%InvGrad_(0,1,i,l,k,lnum)*dom%PMLVeloc_(i,l,k,0,lnum)
                        sum_vy = acoeff*dom%InvGrad_(1,1,i,l,k,lnum)*dom%PMLVeloc_(i,l,k,1,lnum)
                        sum_vz = acoeff*dom%InvGrad_(2,1,i,l,k,lnum)*dom%PMLVeloc_(i,l,k,2,lnum)
                        ForcesFl(0,i,j,k) = ForcesFl(0,i,j,k) + sum_vx
                        ForcesFl(1,i,j,k) = ForcesFl(1,i,j,k) + sum_vy
                        ForcesFl(2,i,j,k) = ForcesFl(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do
        ! TODO reorder loops ?
        do l = 0,ngll-1
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i=0,ngll-1
                        acoeff = - dom%hprime(k,l)*dom%GLLw(i)*dom%GLLw(j)*dom%GLLw(l)*dom%Jacob_(i,j,l,lnum)
                        sum_vx = acoeff*dom%InvGrad_(0,2,i,j,l,lnum)*dom%PMLVeloc_(i,j,l,0,lnum)
                        sum_vy = acoeff*dom%InvGrad_(1,2,i,j,l,lnum)*dom%PMLVeloc_(i,j,l,1,lnum)
                        sum_vz = acoeff*dom%InvGrad_(2,2,i,j,l,lnum)*dom%PMLVeloc_(i,j,l,2,lnum)
                        ForcesFl(0,i,j,k) = ForcesFl(0,i,j,k) + sum_vx
                        ForcesFl(1,i,j,k) = ForcesFl(1,i,j,k) + sum_vy
                        ForcesFl(2,i,j,k) = ForcesFl(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do

        ! Assemblage
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    ind = dom%Idom_(i,j,k,lnum)
                    ! We should have atomic adds with openmp // here
                    champs1%fpml_Forces(ind,0) = champs1%fpml_Forces(ind,0) + ForcesFl(0,i,j,k)
                    champs1%fpml_Forces(ind,1) = champs1%fpml_Forces(ind,1) + ForcesFl(1,i,j,k)
                    champs1%fpml_Forces(ind,2) = champs1%fpml_Forces(ind,2) + ForcesFl(2,i,j,k)
                enddo
            enddo
        enddo
    end subroutine forces_int_flu_pml
#else
        subroutine forces_int_flu_pml(dom, champs1, lnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        type(champsfluidpml), intent(inout) :: champs1
        integer :: lnum
        !
        integer :: ngll
        integer :: i, j, k, l, ind,e,ee
        real(fpp) :: acoeff
        real(fpp), dimension(0:CHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: FFl
        ngll = dom%ngll

        FFl(:,:,:,:,0) = 0d0
        ! Fx
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
#ifdef SEM_VEC
!$omp simd linear(ee) safelen(CHUNK) private(acoeff)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,lnum)
                    acoeff = -dom%GLLw(i)*dom%GLLw(j)*dom%GLLw(k)*dom%Jacob_(i,j,k,e)*dom%PMLVeloc_(i,j,k,0,e)
                    do l = 0,ngll-1
                        FFl(ee,l,j,k,0) = FFl(ee,l,j,k,0) + dom%hprime(l,i)*acoeff*dom%InvGrad_(0,0,i,j,k,e)
                        FFl(ee,i,l,k,0) = FFl(ee,i,l,k,0) + dom%hprime(l,j)*acoeff*dom%InvGrad_(0,1,i,j,k,e)
                        FFl(ee,i,j,l,0) = FFl(ee,i,j,l,0) + dom%hprime(l,k)*acoeff*dom%InvGrad_(0,2,i,j,k,e)
                    end do
                    END_SUBELEM_LOOP()
                end do
            end do
        end do
        FFl(:,:,:,:,1) = 0d0
        ! Fy
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
#ifdef SEM_VEC
!$omp simd linear(ee) safelen(CHUNK) private(acoeff)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,lnum)
                    acoeff = -dom%GLLw(i)*dom%GLLw(j)*dom%GLLw(k)*dom%Jacob_(i,j,k,e)*dom%PMLVeloc_(i,j,k,1,e)
                    do l = 0,ngll-1
                        FFl(ee,l,j,k,1) = FFl(ee,l,j,k,1) + dom%hprime(l,i)*acoeff*dom%InvGrad_(1,0,i,j,k,e)
                        FFl(ee,i,l,k,1) = FFl(ee,i,l,k,1) + dom%hprime(l,j)*acoeff*dom%InvGrad_(1,1,i,j,k,e)
                        FFl(ee,i,j,l,1) = FFl(ee,i,j,l,1) + dom%hprime(l,k)*acoeff*dom%InvGrad_(1,2,i,j,k,e)
                    end do
                    END_SUBELEM_LOOP()
                end do
            end do
        end do
        FFl(:,:,:,:,2) = 0d0
        ! Fz
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
#ifdef SEM_VEC
!$omp simd linear(ee) safelen(CHUNK) private(acoeff)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,lnum)
                    acoeff = -dom%GLLw(i)*dom%GLLw(j)*dom%GLLw(k)*dom%Jacob_(i,j,k,e)*dom%PMLVeloc_(i,j,k,2,e)
                    do l = 0,ngll-1
                        FFl(ee,l,j,k,2) = FFl(ee,l,j,k,2) + dom%hprime(l,i)*acoeff*dom%InvGrad_(2,0,i,j,k,e)
                        FFl(ee,i,l,k,2) = FFl(ee,i,l,k,2) + dom%hprime(l,j)*acoeff*dom%InvGrad_(2,1,i,j,k,e)
                        FFl(ee,i,j,l,2) = FFl(ee,i,j,l,2) + dom%hprime(l,k)*acoeff*dom%InvGrad_(2,2,i,j,k,e)
                    end do
                    END_SUBELEM_LOOP()
                end do
            end do
        end do

        ! Assemblage
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP(e,ee,lnum)
                    if (e>=dom%nbelem) exit
                    ind = dom%Idom_(i,j,k,e)
                    ! We should have atomic adds with openmp // here
                    champs1%fpml_Forces(ind,0) = champs1%fpml_Forces(ind,0) + FFl(ee,i,j,k,0)
                    champs1%fpml_Forces(ind,1) = champs1%fpml_Forces(ind,1) + FFl(ee,i,j,k,1)
                    champs1%fpml_Forces(ind,2) = champs1%fpml_Forces(ind,2) + FFl(ee,i,j,k,2)
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    end subroutine forces_int_flu_pml
#endif

    subroutine pred_flu_pml(dom, dt, champs1, lnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        real(fpp), intent(in) :: dt
        type(champsfluidpml), intent(inout) :: champs1
        integer :: lnum
        !
        real(fpp) :: dVelPhi_dx, dVelPhi_dy, dVelPhi_dz
        integer :: ngll
        integer :: i, j, k, l, ind, e, ee
        real, dimension(0:CHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: VelPhi
        real(fpp) :: dVPhi_dxi,dVPhi_deta,dVPhi_dzeta
        real(fpp) :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3

        ngll = dom%ngll

        ! prediction in the element
        ! We do the sum V1+V2+V3 and F1+F2+F3 here
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP(e,ee,lnum)
                    ind = dom%Idom_(i,j,k,e)
                    VelPhi(ee,i,j,k) = champs1%fpml_VelPhi(ind,0) + &
                        champs1%fpml_VelPhi(ind,1) + &
                        champs1%fpml_VelPhi(ind,2)
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
#ifdef SEM_VEC
!$omp simd linear(ee) safelen(CHUNK) private(dVPhi_dxi,dVPhi_deta,dVPhi_dzeta,dVelPhi_dx, dVelPhi_dy, dVelPhi_dz)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,lnum)
                    ! d(rho*Phi)_d(xi,eta,zeta)
                    dVPhi_dxi   = 0D0
                    dVPhi_deta  = 0D0
                    dVPhi_dzeta = 0D0

                    DO L = 0, ngll-1
                        dVPhi_dxi   = dVPhi_dxi  +VelPhi(ee,L,J,K)*dom%hprime(L,I)
                        dVPhi_deta  = dVPhi_deta +VelPhi(ee,I,L,K)*dom%hprime(L,J)
                        dVPhi_dzeta = dVPhi_dzeta+VelPhi(ee,I,J,L)*dom%hprime(L,K)
                    END DO

                    xi1 = dom%InvGrad_(0,0,i,j,k,e)
                    xi2 = dom%InvGrad_(1,0,i,j,k,e)
                    xi3 = dom%InvGrad_(2,0,i,j,k,e)
                    et1 = dom%InvGrad_(0,1,i,j,k,e)
                    et2 = dom%InvGrad_(1,1,i,j,k,e)
                    et3 = dom%InvGrad_(2,1,i,j,k,e)
                    ga1 = dom%InvGrad_(0,2,i,j,k,e)
                    ga2 = dom%InvGrad_(1,2,i,j,k,e)
                    ga3 = dom%InvGrad_(2,2,i,j,k,e)
                    !- in the physical domain
                    dVelPhi_dx = dVPhi_dxi*xi1 + dVPhi_deta*et1 + dVPhi_dzeta*ga1
                    dVelPhi_dy = dVPhi_dxi*xi2 + dVPhi_deta*et2 + dVPhi_dzeta*ga2
                    dVelPhi_dz = dVPhi_dxi*xi3 + dVPhi_deta*et3 + dVPhi_dzeta*ga3

                    ! prediction for (physical) velocity (which is the equivalent of a stress, here)
                    ! V_x^x
                    dom%PMLVeloc_(i,j,k,0,e) = dom%PMLDumpSx_(i,j,k,0,e) * dom%PMLVeloc_(i,j,k,0,e) + &
                        dom%PMLDumpSx_(i,j,k,1,e) * Dt * dVelPhi_dx
                    ! V_x^y = 0
                    ! V_x^z = 0
                    ! V_y^x = 0
                    ! V_y^y
                    dom%PMLVeloc_(i,j,k,1,e) = dom%PMLDumpSy_(i,j,k,0,e) * dom%PMLVeloc_(i,j,k,1,e) + &
                        dom%PMLDumpSy_(i,j,k,1,e) * Dt * dVelPhi_dy
                    ! V_y^z = 0
                    ! V_z^x = 0
                    ! V_z^y = 0
                    ! V_z^z
                    dom%PMLVeloc_(i,j,k,2,e) = dom%PMLDumpSz_(i,j,k,0,e) * dom%PMLVeloc_(i,j,k,2,e) + &
                        dom%PMLDumpSz_(i,j,k,1,e) * Dt * dVelPhi_dz
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    end subroutine Pred_Flu_Pml
end module dom_fluidpml

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
