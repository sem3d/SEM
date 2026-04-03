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
#include "optims.h"
#include "loops.h"
#include "gllopt.h"

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
        integer :: nbelem, ngll, nblocks, i
        !

        nbelem = dom%nbelem
        ngll   = dom%ngll
        if (ngll == 0) return ! Domain doesn''t exist anywhere
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
            do i=0,Tdomain%TimeD%nsubsteps
                call allocate_champs_fluidpml(dom, i)
            end do
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
        dom%dt = Tdomain%TimeD%dtmin
    end subroutine init_domain_fluidpml

    subroutine start_domain_fluidpml(Tdomain, fpmldom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_fluidpml), intent(inout) :: fpmldom
        !
        integer :: i

        !$acc  enter data copyin(fpmldom, fpmldom%champs) &
        !$acc  copyin(fpmldom%DumpMass, fpmldom%DumpV, fpmldom%m_Lambda, fpmldom%m_Density) &
        !$acc  copyin(fpmldom%m_PMLVeloc) &
        !$acc  copyin(fpmldom%m_PMLDumpSx) &
        !$acc  copyin(fpmldom%m_PMLDumpSy) &
        !$acc  copyin(fpmldom%m_PMLDumpSz) &
        !$acc&
        do i = 0,1
            !$acc enter data  &
            !$acc& copyin(fpmldom%champs(i)%fpml_Phi)    &
            !$acc& copyin(fpmldom%champs(i)%fpml_VelPhi) &
            !$acc& copyin(fpmldom%champs(i)%fpml_Forces)
        end do

    end subroutine start_domain_fluidpml

    subroutine stop_domain_fluidpml(Tdomain, fpmldom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_fluidpml), intent(inout) :: fpmldom
        !
        integer :: i

        !$acc  exit data delete(fpmldom, fpmldom%champs) &
        !$acc  delete(fpmldom%DumpMass, fpmldom%DumpV, fpmldom%m_Lambda, fpmldom%m_Density) &
        !$acc  delete(fpmldom%m_PMLVeloc) &
        !$acc  delete(fpmldom%m_PMLDumpSx) &
        !$acc  delete(fpmldom%m_PMLDumpSy) &
        !$acc  delete(fpmldom%m_PMLDumpSz) &
        !$acc&
        do i = 0,1
            !$acc exit data  delete(fpmldom%champs(i)%fpml_Phi, fpmldom%champs(i)%fpml_VelPhi, fpmldom%champs(i)%fpml_Forces)
        end do

    end subroutine stop_domain_fluidpml

    subroutine get_fluidpml_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        !$acc routine worker
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

    subroutine forces_int_fluid_pml_mainloop(dom, i0, i1)
        use m_calcul_forces_fluid_pml
        type (domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i0, i1

        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_fpml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_5(calcul_forces_fpml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_6(calcul_forces_fpml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_7(calcul_forces_fpml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_8(calcul_forces_fpml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_9(calcul_forces_fpml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_N(calcul_forces_fpml,,(dom,dom%ngll,dom%champs(i1)))
        end select

    end subroutine forces_int_fluid_pml_mainloop


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

    subroutine newmark_predictor_fluidpml(dom, Tdomain, f0, f1)
        type(domain_fluidpml), intent (INOUT) :: dom
        type (domain), intent (INOUT) :: Tdomain
        integer, intent(in) :: f0, f1
        !
        call newmark_predictor_fluidpml_sub(dom, Tdomain, dom%nglltot, &
            Tdomain%fdom%nglltot, Tdomain%fdom%champs(f0)%VelPhi, &
            dom%champs(f0)%fpml_VelPhi, dom%champs(f1)%fpml_VelPhi, &
            dom%champs(f1)%fpml_Forces, &
            Tdomain%intFluPml%surf0%nbtot, Tdomain%intFluPml%surf0%map, Tdomain%intFluPml%surf1%map )
    end subroutine newmark_predictor_fluidpml

    subroutine newmark_predictor_fluidpml_sub(dom, Tdomain, nglltot, ngllfluid, &
        fluid_VelPhi, fpml_VelPhi0, fpml_VelPhi1, fpml_Forces, &
        niface, map_flu, map_pml)
        type(domain_fluidpml), intent (INOUT) :: dom
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: nglltot, ngllfluid, niface
        real(fpp), intent(in), dimension(0:ngllfluid) :: fluid_VelPhi
        real(fpp), intent(inout), dimension(0:dom%nglltot,0:2) :: fpml_VelPhi0, fpml_VelPhi1, fpml_Forces
        integer, intent(in), dimension(0:niface-1) :: map_flu, map_pml
        !
        integer :: n, i, indpml, indflu

        fpml_Forces = 0.
        !$acc parallel loop async(1) &
        !$acc&  copyin(map_pml, map_flu, fluid_VelPhi) &
        !$acc&  present(fpml_VelPhi0) &
        !$acc&  private(indpml, indflu)
        do n = 0,niface-1
            ! Couplage avec l''interface fluide / PML
            indflu = map_flu(n)
            indpml = map_pml(n)
            fpml_VelPhi0(indpml,0) = fluid_VelPhi(indflu)
            fpml_VelPhi0(indpml,1) = 0.
            fpml_VelPhi0(indpml,2) = 0.
        enddo
        !$acc end parallel loop
        ! Prediction

        !$acc parallel loop async(1) collapse(2) &
        !$acc&  copyin(fpml_VelPhi0) &
        !$acc&  copyout(fpml_VelPhi1)
        do i = 0,2
            do n = 0, nglltot
                fpml_Velphi1(n,i) = fpml_VelPhi0(n,i)
            end do
        end do
        !$acc end parallel loop
    end subroutine newmark_predictor_fluidpml_sub

    subroutine newmark_corrector_fluidpml(dom, dt, i0, i1)
        type(domain_fluidpml), intent (INOUT) :: dom
        real(fpp), intent(in) :: dt
        integer, intent(in) :: i0, i1
        !

        call newmark_corrector_fluidpml_sub(dom, dt, dom%nglltot, &
            dom%champs(i0)%fpml_Phi, dom%champs(i0)%fpml_VelPhi, &
            dom%champs(i1)%fpml_Forces, dom%DumpV, &
            dom%n_dirich, dom%dirich)


!        integer  :: n,  indpml
!        dom%champs(i0)%fpml_VelPhi(:,:) = dom%DumpV(:,0,:) * dom%champs(i0)%fpml_VelPhi(:,:) + &
!                                       dt * dom%DumpV(:,1,:) * dom%champs(i1)%fpml_Forces(:,:)
!        do n = 0, dom%n_dirich-1
!            indpml = dom%dirich(n)
!            dom%champs(i0)%fpml_VelPhi(indpml,0) = 0.
!            dom%champs(i0)%fpml_VelPhi(indpml,1) = 0.
!            dom%champs(i0)%fpml_VelPhi(indpml,2) = 0.
!        enddo
!        dom%champs(i0)%fpml_Phi = dom%champs(i0)%fpml_Phi + dt*dom%champs(i0)%fpml_VelPhi
    end subroutine newmark_corrector_fluidpml

    subroutine newmark_corrector_fluidpml_sub(dom, dt, nglltot, &
        fpml_Phi, fpml_VelPhi, fpml_Forces, DumpV, &
        n_dirich, dirich )
        type(domain_fluidpml), intent (INOUT) :: dom
        real(fpp), intent(in) :: dt
        integer, intent(in) :: nglltot, n_dirich
        integer, intent(in), dimension(0:n_dirich-1) :: dirich
        real(fpp), intent(in), dimension(0:nglltot,0:1,0:2) :: DumpV
        real(fpp), intent(in), dimension(0:nglltot,0:2) :: fpml_Forces
        real(fpp), intent(inout), dimension(0:nglltot,0:2) :: fpml_Phi, fpml_VelPhi
        !
        integer  :: n,  indpml, i

        !! XXX tester integrer les condition dirichlet dans DumpV(dirich(n)) = 0
        !!

        !$acc  parallel loop &
        !$acc& async(1) &
        !$acc& collapse(2) &
        !$acc& copyin(DumpV,fpml_Forces) &
        !$acc& present(fpml_VelPhi) &
        !$acc& firstprivate(nglltot,dt)
        do i=0,2
            do n=0,nglltot
                fpml_VelPhi(n,i) = DumpV(n,0,i) * fpml_VelPhi(n,i) + &
                    dt * DumpV(n,1,i) * fpml_Forces(n,i)
            end do
        end do
        !$acc end parallel loop

        !$acc  parallel loop &
        !$acc& async(1) &
        !$acc& copyin(dirich) &
        !$acc& present(fpml_VelPhi) &
        !$acc& firstprivate(n_dirich)
        do n = 0, n_dirich-1
            indpml = dirich(n)
            fpml_VelPhi(indpml,0) = 0.
            fpml_VelPhi(indpml,1) = 0.
            fpml_VelPhi(indpml,2) = 0.
        end do
        !$acc end parallel loop

        !$acc  parallel loop &
        !$acc& async(1) &
        !$acc& collapse(2) &
        !$acc& copyin(fpml_VelPhi) &
        !$acc& present(fpml_Phi) &
        !$acc& firstprivate(nglltot,dt)
        do i=0,2
            do n=0,nglltot
                fpml_Phi(n,i) = fpml_Phi(n,i) + dt* fpml_VelPhi(n,i)
            end do
        end do
        !$acc end parallel loop
    end subroutine newmark_corrector_fluidpml_sub

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
