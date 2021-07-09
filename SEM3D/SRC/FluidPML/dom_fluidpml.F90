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
    use pml
    implicit none
#include "index.h"

contains

    subroutine allocate_champs_fluidpml(dom, i)
        type(domain_fluidpml), intent(inout) :: dom
        integer, intent(in) :: i

        allocate(dom%champs(i)%fpml_VelPhi(0:dom%nglltot,0:2))
        allocate(dom%champs(i)%fpml_Phi   (0:dom%nglltot,0:2))
        allocate(dom%champs(i)%fpml_Forces(0:dom%nglltot,0:2))
        dom%champs(i)%fpml_VelPhi = 0d0
        dom%champs(i)%fpml_Phi = 0d0
        dom%champs(i)%fpml_Forces = 0d0
    end subroutine allocate_champs_fluidpml

    subroutine allocate_dom_fluidpml (Tdomain, dom)
        use gll3d
        implicit none
        type(domain), intent(inout) :: TDomain
        type(domain_fluidpml), intent(inout) :: dom
        !
        integer :: nbelem, ngll, nblocks
        !

        nbelem = dom%nbelem
        ngll   = dom%ngll
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call init_dombase(dom)

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)

            nblocks = dom%nblocks

            allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1,0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1,0:VCHUNK-1))

            if(Tdomain%TimeD%velocity_scheme)then
                allocate(dom%PMLVeloc_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nblocks-1,0:VCHUNK-1))
                dom%PMLVeloc_(:,:,:,:,:,:) = 0d0
                allocate(dom%PMLDumpSx_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nblocks-1,0:VCHUNK-1))
                allocate(dom%PMLDumpSy_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nblocks-1,0:VCHUNK-1))
                allocate(dom%PMLDumpSz_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nblocks-1,0:VCHUNK-1))
                dom%PMLDumpSx_(:,:,:,:,:,:) = 0d0
                dom%PMLDumpSy_(:,:,:,:,:,:) = 0d0
                dom%PMLDumpSz_(:,:,:,:,:,:) = 0d0
            endif
        end if

        ! Allocation et initialisation de champs0 pour les PML fluides
        if (dom%nglltot /= 0) then
            call allocate_champs_fluidpml(dom, 0)
            call allocate_champs_fluidpml(dom, 1)
            allocate(dom%DumpV (0:dom%nglltot,0:1,0:2))
            allocate(dom%DumpMass(0:dom%nglltot,0:2))
            dom%DumpV = 0d0
            dom%DumpMass = 0d0
            dom%DumpV(dom%nglltot,1,:) = 1d0
            dom%DumpMass(dom%nglltot,:) = 1d0
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - fluid pml domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_fluidpml

    subroutine deallocate_dom_fluidpml (dom)
        implicit none
        type(domain_fluidpml) :: dom
        !
        integer :: i

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )

        if(allocated(dom%m_PMLVeloc )) deallocate(dom%m_PMLVeloc )
        if(allocated(dom%m_PMLDumpSx)) deallocate(dom%m_PMLDumpSx)
        if(allocated(dom%m_PMLDumpSy)) deallocate(dom%m_PMLDumpSy)
        if(allocated(dom%m_PMLDumpSz)) deallocate(dom%m_PMLDumpSz)

        do i=0,1
            if(allocated(dom%champs(i)%fpml_VelPhi)) deallocate(dom%champs(i)%fpml_VelPhi)
            if(allocated(dom%champs(i)%fpml_Phi   )) deallocate(dom%champs(i)%fpml_Phi   )
        end do

        if(allocated(dom%DumpMass)) deallocate(dom%DumpMass)
        if(allocated(dom%DumpV )) deallocate(dom%DumpV )

        call deallocate_dombase(dom)
    end subroutine deallocate_dom_fluidpml

    subroutine init_domain_fluidpml(Tdomain, dom)
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_fluidpml), intent(inout) :: dom
        !
        ! TODO : useless, kill this method, needed for build compatibility SolidPML / SolidCPML
    end subroutine init_domain_fluidpml

    subroutine get_fluidpml_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain_fluidpml), intent(inout)       :: dom
        integer, dimension(0:), intent(in)         :: out_variables
        integer                                    :: lnum
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldU
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldV
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldA
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: fieldP
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: P_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: S_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:8) :: dUdX
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: sig_dev
        !
        logical :: flag_gradU
        integer :: ngll, i, j, k, ind
        !
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

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)

                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        fieldU(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_VITESSE) == 1) then
                        fieldV(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        fieldA(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        fieldP(i,j,k) = 0d0
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        eps_vol(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        P_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        S_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        eps_dev(i,j,k,:) = 0.
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        sig_dev(i,j,k,:) = 0.
                    end if

                    if (out_variables(OUT_DUDX) == 1) then
                       dUdX(i,j,k,:) = 0.
                    end if
                enddo
            enddo
        enddo
    end subroutine get_fluidpml_dom_var

    subroutine init_material_properties_fluidpml(dom, lnum, mat, density, lambda)
        type(domain_fluidpml), intent(inout) :: dom
        integer, intent(in) :: lnum
        type (subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: density
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: lambda
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        dom%Density_(:,:,:,bnum,ee) = density
        dom%Lambda_ (:,:,:,bnum,ee) = lambda
    end subroutine init_material_properties_fluidpml

    subroutine init_local_mass_fluidpml(dom,specel,i,j,k,ind,Whei)
        type(domain_fluidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real(fpp) :: Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ! Fluid : inertial term ponderation by the inverse of the bulk modulus

        specel%MassMat(i,j,k) = Whei*dom%Jacob_(i,j,k,bnum,ee)/dom%Lambda_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_fluidpml

    subroutine forces_int_flu_pml_dim(dom, champs1, bnum, d)
        type (domain_fluidpml), intent (INOUT) :: dom
        type(champsfluidpml), intent(inout) :: champs1
        integer :: bnum, d
        !
        integer :: ngll
        integer :: i, j, k, l, ind,e,ee
        real(fpp) :: acoeff
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: FFl
        ngll = dom%ngll

        FFl(:,:,:,:) = 0d0
        ! F(dim=d)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
#if VCHUNK>1
!$omp simd linear(ee) safelen(VCHUNK) private(acoeff)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    acoeff = -dom%GLLw(i)*dom%GLLw(j)*dom%GLLw(k)*dom%Jacob_(i,j,k,bnum,ee) &
                        * dom%PMLVeloc_(i,j,k,d,bnum,ee)
                    do l = 0,ngll-1
                        FFl(ee,l,j,k) = FFl(ee,l,j,k) + dom%hprime(l,i)*acoeff*dom%InvGrad_(d,0,i,j,k,bnum,ee)
                        FFl(ee,i,l,k) = FFl(ee,i,l,k) + dom%hprime(l,j)*acoeff*dom%InvGrad_(d,1,i,j,k,bnum,ee)
                        FFl(ee,i,j,l) = FFl(ee,i,j,l) + dom%hprime(l,k)*acoeff*dom%InvGrad_(d,2,i,j,k,bnum,ee)
                    end do
                    END_SUBELEM_LOOP()
                end do
            end do
        end do

        ! Assemblage
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    ! We should have atomic adds with openmp // here
                    champs1%fpml_Forces(ind,d) = champs1%fpml_Forces(ind,d) + FFl(ee,i,j,k)
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    end subroutine forces_int_flu_pml_dim

    subroutine forces_int_flu_pml(dom, champs1, bnum, Tdomain)
        type (domain_fluidpml), intent (INOUT) :: dom
        type(champsfluidpml), intent(inout) :: champs1
        type (domain), intent (INOUT), target :: Tdomain
        integer :: bnum
        !
        call forces_int_flu_pml_dim(dom, champs1, bnum, 0)
        call forces_int_flu_pml_dim(dom, champs1, bnum, 1)
        call forces_int_flu_pml_dim(dom, champs1, bnum, 2)
    end subroutine forces_int_flu_pml

    subroutine pred_flu_pml(dom, dt, champs1, bnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        real(fpp), intent(in) :: dt
        type(champsfluidpml), intent(inout) :: champs1
        integer :: bnum
        !
        real(fpp) :: dVelPhi_dx, dVelPhi_dy, dVelPhi_dz
        integer :: ngll
        integer :: i, j, k, l, ind, e, ee
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: VelPhi
        real(fpp) :: dVPhi_dxi,dVPhi_deta,dVPhi_dzeta
        real(fpp) :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3

        ngll = dom%ngll

        ! prediction in the element
        ! We do the sum V1+V2+V3 and F1+F2+F3 here
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
#if VCHUNK>1
!$omp simd linear(ee) safelen(VCHUNK) private(ind)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    ind = dom%Idom_(i,j,k,bnum,ee)
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
#if VCHUNK>1
!$omp simd linear(ee) safelen(VCHUNK) private(dVPhi_dxi,dVPhi_deta,dVPhi_dzeta,dVelPhi_dx, dVelPhi_dy, dVelPhi_dz)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    ! d(rho*Phi)_d(xi,eta,zeta)
                    dVPhi_dxi   = 0D0
                    dVPhi_deta  = 0D0
                    dVPhi_dzeta = 0D0

                    DO L = 0, ngll-1
                        dVPhi_dxi   = dVPhi_dxi  +VelPhi(ee,L,J,K)*dom%hprime(L,I)
                        dVPhi_deta  = dVPhi_deta +VelPhi(ee,I,L,K)*dom%hprime(L,J)
                        dVPhi_dzeta = dVPhi_dzeta+VelPhi(ee,I,J,L)*dom%hprime(L,K)
                    END DO

                    xi1 = dom%InvGrad_(0,0,i,j,k,bnum,ee)
                    xi2 = dom%InvGrad_(1,0,i,j,k,bnum,ee)
                    xi3 = dom%InvGrad_(2,0,i,j,k,bnum,ee)
                    et1 = dom%InvGrad_(0,1,i,j,k,bnum,ee)
                    et2 = dom%InvGrad_(1,1,i,j,k,bnum,ee)
                    et3 = dom%InvGrad_(2,1,i,j,k,bnum,ee)
                    ga1 = dom%InvGrad_(0,2,i,j,k,bnum,ee)
                    ga2 = dom%InvGrad_(1,2,i,j,k,bnum,ee)
                    ga3 = dom%InvGrad_(2,2,i,j,k,bnum,ee)
                    !- in the physical domain
                    dVelPhi_dx = dVPhi_dxi*xi1 + dVPhi_deta*et1 + dVPhi_dzeta*ga1
                    dVelPhi_dy = dVPhi_dxi*xi2 + dVPhi_deta*et2 + dVPhi_dzeta*ga2
                    dVelPhi_dz = dVPhi_dxi*xi3 + dVPhi_deta*et3 + dVPhi_dzeta*ga3

                    ! prediction for (physical) velocity (which is the equivalent of a stress, here)
                    ! V_x^x
                    dom%PMLVeloc_(i,j,k,0,bnum,ee) = dom%PMLDumpSx_(i,j,k,0,bnum,ee) * &
                        dom%PMLVeloc_(i,j,k,0,bnum,ee) + &
                        dom%PMLDumpSx_(i,j,k,1,bnum,ee) * Dt * dVelPhi_dx
                    ! V_x^y = 0
                    ! V_x^z = 0
                    ! V_y^x = 0
                    ! V_y^y
                    dom%PMLVeloc_(i,j,k,1,bnum,ee) = dom%PMLDumpSy_(i,j,k,0,bnum,ee) * &
                        dom%PMLVeloc_(i,j,k,1,bnum,ee) + &
                        dom%PMLDumpSy_(i,j,k,1,bnum,ee) * Dt * dVelPhi_dy
                    ! V_y^z = 0
                    ! V_z^x = 0
                    ! V_z^y = 0
                    ! V_z^z
                    dom%PMLVeloc_(i,j,k,2,bnum,ee) = dom%PMLDumpSz_(i,j,k,0,bnum,ee) * &
                        dom%PMLVeloc_(i,j,k,2,bnum,ee) + &
                        dom%PMLDumpSz_(i,j,k,1,bnum,ee) * Dt * dVelPhi_dz
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    end subroutine Pred_Flu_Pml

    subroutine init_fluidpml_properties(Tdomain,specel,mat)
        type (domain), intent (INOUT), target :: Tdomain
        type (element), intent(inout) :: specel
        type (subdomain), intent(in) :: mat
        !
        integer :: ngll, lnum
        real(fpp), dimension(:,:,:), allocatable :: temp_PMLx,temp_PMLy
        real(fpp), dimension(:,:,:), allocatable :: wx,wy,wz
        real(fpp) :: dt
        real(fpp), dimension(:,:,:,:), allocatable :: PMLDumpMass
        integer :: i,j,k,idx,m,ind
        real(fpp), dimension(:,:,:), allocatable   :: Vp
        real(fpp), dimension(:,:,:,:), allocatable :: coords
        integer :: bnum, ee

        dt = Tdomain%TimeD%dtmin
        lnum = specel%lnum
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        ngll = domain_ngll(Tdomain, specel%domain)

        allocate(Vp(0:ngll-1,0:ngll-1,0:ngll-1))
        Vp = sqrt(Tdomain%fpmldom%Lambda_(:,:,:,bnum,ee)/Tdomain%fpmldom%Density_(:,:,:,bnum,ee))

        ! PML case: valid for solid and fluid parts

        !- definition of the attenuation coefficient in PMLs (alpha in the literature)
        allocate(wx(0:ngll-1,0:ngll-1,0:ngll-1))
        allocate(wy(0:ngll-1,0:ngll-1,0:ngll-1))
        allocate(wz(0:ngll-1,0:ngll-1,0:ngll-1))

        allocate(coords(0:ngll-1,0:ngll-1,0:ngll-1,0:2))

        DO K=0,ngll-1
            DO J=0,ngll-1
                DO I=0,ngll-1
                    idx = specel%Iglobnum(I,J,K)
                    coords(I,J,K,:) = Tdomain%GlobCoord(:,idx)
                END DO
            END DO
        END DO
        call define_alpha_PML(coords, 0, ngll, Vp, mat%pml_width, mat%pml_pos, mat%Apow, mat%npow, wx)
        call define_alpha_PML(coords, 1, ngll, Vp, mat%pml_width, mat%pml_pos, mat%Apow, mat%npow, wy)
        call define_alpha_PML(coords, 2, ngll, Vp, mat%pml_width, mat%pml_pos, mat%Apow, mat%npow, wz)

        !- M-PMLs
        if(Tdomain%logicD%MPML)then
            allocate(temp_PMLx(0:ngll-1,0:ngll-1,0:ngll-1))
            allocate(temp_PMLy(0:ngll-1,0:ngll-1,0:ngll-1))
            temp_PMLx(:,:,:) = wx(:,:,:)
            temp_PMLy(:,:,:) = wy(:,:,:)
            wx(:,:,:) = wx(:,:,:)+Tdomain%MPML_coeff*(wy(:,:,:)+wz(:,:,:))
            wy(:,:,:) = wy(:,:,:)+Tdomain%MPML_coeff*(temp_PMLx(:,:,:)+wz(:,:,:))
            wz(:,:,:) = wz(:,:,:)+Tdomain%MPML_coeff*(temp_PMLx(:,:,:)+temp_PMLy(:,:,:))
            deallocate(temp_PMLx,temp_PMLy)
        end if

        allocate(PMLDumpMass(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        PMLDumpMass = 0d0

        !- strong formulation for stresses. Dumped mass elements, convolutional terms.
        ! Compute DumpS(x,y,z) and DumpMass(0,1,2)
        call define_PML_DumpInit(ngll,dt,wx,specel%MassMat, &
            Tdomain%fpmldom%PMLDumpSx_(:,:,:,:,bnum,ee),PMLDumpMass(:,:,:,0))
        call define_PML_DumpInit(ngll,dt,wy,specel%MassMat, &
            Tdomain%fpmldom%PMLDumpSy_(:,:,:,:,bnum,ee),PMLDumpMass(:,:,:,1))
        call define_PML_DumpInit(ngll,dt,wz,specel%MassMat, &
            Tdomain%fpmldom%PMLDumpSz_(:,:,:,:,bnum,ee),PMLDumpMass(:,:,:,2))
        deallocate(wx,wy,wz)

        ! Assemble dump mass
        do m = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        ind = specel%Idom(i,j,k)
                        Tdomain%fpmldom%DumpMass(ind,m) =   Tdomain%fpmldom%DumpMass(ind,m) &
                                                          + PMLDumpMass(i,j,k,m)
                    enddo
                enddo
            enddo
        enddo
        if(allocated(PMLDumpMass)) deallocate(PMLDumpMass)

        !! XXX
        Tdomain%fpmldom%PMLDumpSx_(:,:,:,1,bnum,ee) = Tdomain%fpmldom%PMLDumpSx_(:,:,:,1,bnum,ee) / &
                                                   Tdomain%fpmldom%Density_(:,:,:  ,bnum,ee)
        Tdomain%fpmldom%PMLDumpSy_(:,:,:,1,bnum,ee) = Tdomain%fpmldom%PMLDumpSy_(:,:,:,1,bnum,ee) / &
                                                   Tdomain%fpmldom%Density_(:,:,:  ,bnum,ee)
        Tdomain%fpmldom%PMLDumpSz_(:,:,:,1,bnum,ee) = Tdomain%fpmldom%PMLDumpSz_(:,:,:,1,bnum,ee) / &
                                                   Tdomain%fpmldom%Density_(:,:,:  ,bnum,ee)

        deallocate(Vp)
    end subroutine init_fluidpml_properties

    subroutine finalize_fluidpml_properties(Tdomain,dom)
      type (domain), intent (INOUT), target :: Tdomain
      type (domain_fluidpml), intent (INOUT), target :: dom
      !
      call define_PML_DumpEnd(dom%nglltot, dom%MassMat, dom%DumpMass, dom%DumpV)
    end subroutine finalize_fluidpml_properties

    subroutine newmark_predictor_fluidpml(dom, Tdomain, i0, i1)
        type(domain_fluidpml), intent (INOUT) :: dom
        type(domain), intent(inout)   :: Tdomain
        integer, intent(in) :: i0, i1
        !
        integer :: n, indpml, indflu
        real(fpp) :: bega, dt

        bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
        dt = Tdomain%TimeD%dtmin

        dom%champs(i1)%fpml_Forces = 0.
        do n = 0,Tdomain%intFluPml%surf0%nbtot-1
            ! Couplage à l'interface fluide / PML
            indflu = Tdomain%intFluPml%surf0%map(n)
            indpml = Tdomain%intFluPml%surf1%map(n)
            dom%champs(i0)%fpml_VelPhi(indpml,0) = Tdomain%fdom%champs(i0)%VelPhi(indflu)
            dom%champs(i0)%fpml_VelPhi(indpml,1) = 0.
            dom%champs(i0)%fpml_VelPhi(indpml,2) = 0.
        enddo
        ! Prediction
        dom%champs(i1)%fpml_Velphi = dom%champs(i0)%fpml_VelPhi + dt*(0.5-bega)*dom%champs(i1)%fpml_Forces
    end subroutine newmark_predictor_fluidpml

    subroutine newmark_corrector_fluidpml(dom, dt, i0, i1)
        type(domain_fluidpml), intent (INOUT) :: dom
        real(fpp), intent(in) :: dt
        integer, intent(in) :: i0, i1
        !
        integer  :: n,  indpml

        dom%champs(i0)%fpml_VelPhi(:,:) = dom%DumpV(:,0,:) * dom%champs(i0)%fpml_VelPhi(:,:) + &
                                       dt * dom%DumpV(:,1,:) * dom%champs(i1)%fpml_Forces(:,:)
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs(i0)%fpml_VelPhi(indpml,0) = 0.
            dom%champs(i0)%fpml_VelPhi(indpml,1) = 0.
            dom%champs(i0)%fpml_VelPhi(indpml,2) = 0.
        enddo
        dom%champs(i0)%fpml_Phi = dom%champs(i0)%fpml_Phi + dt*dom%champs(i0)%fpml_VelPhi
    end subroutine newmark_corrector_fluidpml

    function fluidpml_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        type(domain_fluidpml), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        Pspeed = sqrt(dom%Lambda_(i,j,k,bnum,ee)/dom%Density_(i,j,k,bnum,ee))
    end function fluidpml_Pspeed
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
