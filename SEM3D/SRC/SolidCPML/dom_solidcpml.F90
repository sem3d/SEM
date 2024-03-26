!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

!! Convolutional Perfectly Matched Layers implemented according 2 references:
!! Ref1. Improved forward wave propagation and adjoint-based sensitivity kernel calculations using a numerically stable finite-element PML
!!       Zhinan Xie, Dimitri Komatitsch, Roland Martin, Rene Matzen
!!       Geophysical Journal International, 2014, 198, 1714-1747
!! Ref2. An efficient finite element time-domain formulation for the elastic second-order wave equation: a non-split complex frequency shifted convolutional PML
!!       Rene Matzen
!!       International Journal For Numerical Methods In Engineering, 2011, 88, 951-973

module dom_solidpml
    use constants
    use sdomain
    use champs_solidpml
    use ssubdomains
    implicit none
#include "index.h"

contains

    subroutine allocate_champs_solidcpml(dom, i)
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: i

        ! Allocate ONE more gll than necessary to use as a dummy target for
        ! indirections for fake elements.
        allocate(dom%champs(i)%Forces(0:dom%nglltot,0:2))
        allocate(dom%champs(i)%Depla (0:dom%nglltot,0:2))
        allocate(dom%champs(i)%Veloc (0:dom%nglltot,0:2))
        dom%champs(i)%Depla  = 0d0
        dom%champs(i)%Veloc  = 0d0
        dom%champs(i)%Forces = 0d0
    end subroutine allocate_champs_solidcpml

    subroutine allocate_dom_solidpml (Tdomain, dom)
        use pml
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_solidpml), intent (INOUT) :: dom
        !
        integer nbelem, ngll, nblocks, dir1_count, dir2_count, nbtot_SF, i
        !

        ngll   = dom%ngll
        if (ngll == 0) return ! Domain doesn't exist anywhere
        call cpml_reorder_elements(Tdomain, dom, DM_SOLID_CG_PML)
        nbelem = dom%nbelem
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call init_dombase(dom)

        call cpml_allocate_multi_dir(Tdomain, dom)

        dir1_count = dom%dir1_count
        dir2_count = dom%dir2_count
        if (dir1_count>0) then
            allocate(dom%R1_1(0:2, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:dir1_count-1))
            allocate(dom%R2_1(0:8, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:dir1_count-1))
            dom%R1_1 = 0d0
            dom%R2_1 = 0d0
        end if
        if (dir2_count>0) then
            allocate(dom%R1_2(0:2, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:dir2_count-1))
            allocate(dom%R2_2(0:8, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:dir2_count-1))
            dom%R1_2 = 0d0
            dom%R2_2 = 0d0
        end if

        ! We reset the count here, since we will use them to renumber the specel that
        ! have more than one direction of attenuation
        dom%dir1_count = 0
        dom%dir2_count = 0

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)

            nblocks = dom%nblocks
            allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Mu_     (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))

            allocate(dom%I1      (0:VCHUNK-1, 0:nblocks-1))
            allocate(dom%I2      (0:VCHUNK-1, 0:nblocks-1))
            allocate(dom%D0      (0:VCHUNK-1, 0:nblocks-1))
            allocate(dom%D1      (0:VCHUNK-1, 0:nblocks-1))
            allocate(dom%Alpha_0 (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            allocate(dom%dxi_k_0 (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            allocate(dom%Kappa_0 (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            ! Allocation des Ri pour les PML solides (i = 0...5)
            allocate(dom%R1_0(0:VCHUNK-1,0:2, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            allocate(dom%R2_0(0:VCHUNK-1,0:8, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            dom%R1_0 = 0d0
            dom%R2_0 = 0d0
            dom%I1(:,:) = -1
            dom%I2(:,:) = -1
            dom%D0(:,:) = 0
            dom%D1(:,:) = 0
            dom%Kappa_0 = 1.
        end if

        ! Allocation et initialisation de champs0 pour les PML solides
        if (dom%nglltot /= 0) then
            do i=0,Tdomain%TimeD%nsubsteps
                call allocate_champs_solidcpml(dom, i)
            end do

            ! Allocation de DumpMat pour les PML solides
            allocate(dom%DumpMat(0:dom%nglltot))
            dom%DumpMat = 0d0

            ! Allocation de MasUMat pour les PML solides
            allocate(dom%MasUMat(0:dom%nglltot))
            dom%MasUMat = 0d0

        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid cpml domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"

        if (Tdomain%logicD%save_snapshots) then
            allocate(dom%FDump(0:dom%nglltot, 0:2))
            allocate(dom%FMasU(0:dom%nglltot, 0:2))
            allocate(dom%Fint (0:dom%nglltot, 0:2))
            dom%FDump = 0.
            dom%FMasU = 0.
            dom%Fint  = 0.
        end if

        ! CPML parameters initialisation: for the very first implementation, parameters are hard-coded.

        dom%cpml_c = 1.0
        dom%cpml_n = Tdomain%config%cpml_n
        dom%cpml_rc = Tdomain%config%cpml_rc
        dom%cpml_kappa_0 = Tdomain%config%cpml_kappa0
        dom%cpml_kappa_1 = Tdomain%config%cpml_kappa1
        dom%cpml_integ = Tdomain%config%cpml_integ_type
        dom%cpml_one_root = Tdomain%config%cpml_one_root
        dom%alphamax = 0.
        if(Tdomain%rank==0) then
            write(*,*) "INFO - solid cpml domain : kappa0 ", dom%cpml_kappa_0, " kappa1 ", dom%cpml_kappa_1
        endif

        nbtot_SF = Tdomain%SF%intSolFluPml%surf0%nbtot
        call allocate_dombase_cpml(dom, nbtot_SF)
        if (nbtot_SF >= 0) allocate(dom%R_0_SF(0:1, 0:nbtot_SF-1))
        if (nbtot_SF >= 0) allocate(dom%R_1_SF(0:1, 0:nbtot_SF-1))
        dom%R_0_SF = 0.
        dom%R_1_SF = 0.
    end subroutine allocate_dom_solidpml

    subroutine deallocate_dom_solidpml (dom)
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        !
        integer :: i

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_Mu     )) deallocate(dom%m_Mu     )

        call deallocate_dombase_cpml(dom)
        if (allocated(dom%R_0_SF)) deallocate(dom%R_0_SF)
        if (allocated(dom%R_1_SF)) deallocate(dom%R_1_SF)

        do i=0,1
            if(allocated(dom%champs(i)%Depla )) deallocate(dom%champs(i)%Depla )
            if(allocated(dom%champs(i)%Veloc )) deallocate(dom%champs(i)%Veloc )
            if(allocated(dom%champs(i)%Forces)) deallocate(dom%champs(i)%Forces )
        end do

        if(allocated(dom%DumpMat)) deallocate(dom%DumpMat)
        if(allocated(dom%MasUMat)) deallocate(dom%MasUMat)

        if(allocated(dom%R1_0)) deallocate(dom%R1_0)
        if(allocated(dom%R2_0)) deallocate(dom%R2_0)

        if(allocated(dom%FDump)) deallocate(dom%FDump)
        if(allocated(dom%FMasU)) deallocate(dom%FMasU)
        if(allocated(dom%Fint))  deallocate(dom%Fint )
    end subroutine deallocate_dom_solidpml

    subroutine get_solidpml_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        !$acc routine worker
        implicit none
        !
        type(domain_solidpml)              :: dom
        integer, dimension(0:), intent(in) :: out_variables
        integer                            :: lnum
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldU, fieldV, fieldA
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: fieldP
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: P_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: S_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: sig_dev
        !
        logical                  :: flag_gradU
        integer                  :: ngll, i, j, k, ind
        real(fpp)                :: DXX, DXY, DXZ
        real(fpp)                :: DYX, DYY, DYZ
        real(fpp)                :: DZX, DZY, DZZ
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        flag_gradU = (out_variables(OUT_PRESSION)    + &
                      out_variables(OUT_ENERGYP)     + &
                      out_variables(OUT_ENERGYS)     + &
                      out_variables(OUT_EPS_VOL)     + &
                      out_variables(OUT_EPS_DEV)     + &
                      out_variables(OUT_STRESS_DEV)) /= 0

        ngll = dom%ngll

        ! First, get displacement.

        if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        fieldU(i,j,k,:) = dom%champs(0)%Depla(ind,:)
                    enddo
                enddo
            enddo
        end if

        ! Then, get other variables.

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ! Compute gradU with displacement if needed.

                    if (flag_gradU) then
                        invgrad_ijk = dom%InvGrad_(:,:,i,j,k,bnum,ee) ! cache for performance

                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,0),DXX,DYX,DZX)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,1),DXY,DYY,DZY)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,2),DXZ,DYZ,DZZ)
                    end if

                    ! Get other variables.

                    ind = dom%Idom_(i,j,k,bnum,ee)

                    if (out_variables(OUT_VITESSE) == 1) then
                        fieldV(i,j,k,:) = dom%champs(0)%Veloc(ind,:)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        fieldA(i,j,k,:) = dom%Massmat(ind) * dom%champs(1)%Forces(ind,:)
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        fieldP(i,j,k) = 0d0
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1 .or. &
                        out_variables(OUT_ENERGYP) == 1 .or. &
                        out_variables(OUT_EPS_DEV) == 1 .or. &
                        out_variables(OUT_STRESS_DEV) == 1 ) then
                        eps_vol(i,j,k) = DXX + DYY + DZZ
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        P_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        S_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        eps_dev(i,j,k,0) = DXX - eps_vol(i,j,k) / 3
                        eps_dev(i,j,k,1) = DYY - eps_vol(i,j,k) / 3
                        eps_dev(i,j,k,2) = DZZ - eps_vol(i,j,k) / 3
                        eps_dev(i,j,k,3) = 0.5 * (DXY + DYX)
                        eps_dev(i,j,k,4) = 0.5 * (DZX + DXZ)
                        eps_dev(i,j,k,5) = 0.5 * (DZY + DYZ)
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        sig_dev(i,j,k,:) = 0.
                    end if
                enddo
            enddo
        enddo
    end subroutine get_solidpml_dom_var

    subroutine get_solidpml_rfields(Tdomain, dom, n, lnum, outputs)
        !$acc routine worker
        use msnapdata, only : output_var_t
        use m_calcul_forces_solidpml
        type(domain), intent(inout) :: TDomain
        type(domain_solidpml), intent(inout) :: dom
        type(output_var_t), intent(inout) :: outputs
        integer, intent(in) :: lnum, n
        !
        integer :: bnum, ee, i1, i2
        integer :: i, j, k, idx, ind, dir0, dir1
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        !
        do k=0,dom%ngll-1
            do j=0,dom%ngll-1
                do i=0,dom%ngll-1
                    idx = outputs%irenum(Tdomain%specel(n)%Iglobnum(i,j,k))
                    dir0 = dom%D0(ee,bnum)
                    dir1 = dom%D1(ee,bnum)
                    select case(dir0)
                    case (0)
                        outputs%R1_0(:,idx) = dom%R1_0(ee, :, i, j, k, bnum)
                        outputs%R2_0_dX(0,idx) = dom%R2_0(ee, dxx, i, j, k, bnum)
                        outputs%R2_0_dX(1,idx) = dom%R2_0(ee, dyx, i, j, k, bnum)
                        outputs%R2_0_dX(2,idx) = dom%R2_0(ee, dzx, i, j, k, bnum)
                        outputs%R2_0_dY(0,idx) = dom%R2_0(ee, dxy, i, j, k, bnum)
                        outputs%R2_0_dY(1,idx) = dom%R2_0(ee, dyy, i, j, k, bnum)
                        outputs%R2_0_dY(2,idx) = dom%R2_0(ee, dzy, i, j, k, bnum)
                        outputs%R2_0_dZ(0,idx) = dom%R2_0(ee, dxz, i, j, k, bnum)
                        outputs%R2_0_dZ(1,idx) = dom%R2_0(ee, dyz, i, j, k, bnum)
                        outputs%R2_0_dZ(2,idx) = dom%R2_0(ee, dzz, i, j, k, bnum)
                    case (1)
                        outputs%R1_1(:,idx) = dom%R1_0(ee, :, i, j, k, bnum)
                        outputs%R2_1_dX(0,idx) = dom%R2_0(ee, dxx, i, j, k, bnum)
                        outputs%R2_1_dX(1,idx) = dom%R2_0(ee, dyx, i, j, k, bnum)
                        outputs%R2_1_dX(2,idx) = dom%R2_0(ee, dzx, i, j, k, bnum)
                        outputs%R2_1_dY(0,idx) = dom%R2_0(ee, dxy, i, j, k, bnum)
                        outputs%R2_1_dY(1,idx) = dom%R2_0(ee, dyy, i, j, k, bnum)
                        outputs%R2_1_dY(2,idx) = dom%R2_0(ee, dzy, i, j, k, bnum)
                        outputs%R2_1_dZ(0,idx) = dom%R2_0(ee, dxz, i, j, k, bnum)
                        outputs%R2_1_dZ(1,idx) = dom%R2_0(ee, dyz, i, j, k, bnum)
                        outputs%R2_1_dZ(2,idx) = dom%R2_0(ee, dzz, i, j, k, bnum)
                    case (2)
                        outputs%R1_2(:,idx) = dom%R1_0(ee, :, i, j, k, bnum)
                        outputs%R2_2_dX(0,idx) = dom%R2_0(ee, dxx, i, j, k, bnum)
                        outputs%R2_2_dX(1,idx) = dom%R2_0(ee, dyx, i, j, k, bnum)
                        outputs%R2_2_dX(2,idx) = dom%R2_0(ee, dzx, i, j, k, bnum)
                        outputs%R2_2_dY(0,idx) = dom%R2_0(ee, dxy, i, j, k, bnum)
                        outputs%R2_2_dY(1,idx) = dom%R2_0(ee, dyy, i, j, k, bnum)
                        outputs%R2_2_dY(2,idx) = dom%R2_0(ee, dzy, i, j, k, bnum)
                        outputs%R2_2_dZ(0,idx) = dom%R2_0(ee, dxz, i, j, k, bnum)
                        outputs%R2_2_dZ(1,idx) = dom%R2_0(ee, dyz, i, j, k, bnum)
                        outputs%R2_2_dZ(2,idx) = dom%R2_0(ee, dzz, i, j, k, bnum)
                    end select
                    i1 = dom%I1(ee,bnum)
                    i2 = dom%I2(ee,bnum)
                    if (i1/=-1) then
                        select case(dir1)
                        case (1)
                            outputs%R1_1(:,idx)    = dom%R1_1(  :, i, j, k, i1)
                            outputs%R2_1_dX(0,idx) = dom%R2_1(dxx, i, j, k, i1)
                            outputs%R2_1_dX(1,idx) = dom%R2_1(dyx, i, j, k, i1)
                            outputs%R2_1_dX(2,idx) = dom%R2_1(dzx, i, j, k, i1)
                            outputs%R2_1_dY(0,idx) = dom%R2_1(dxy, i, j, k, i1)
                            outputs%R2_1_dY(1,idx) = dom%R2_1(dyy, i, j, k, i1)
                            outputs%R2_1_dY(2,idx) = dom%R2_1(dzy, i, j, k, i1)
                            outputs%R2_1_dZ(0,idx) = dom%R2_1(dxz, i, j, k, i1)
                            outputs%R2_1_dZ(1,idx) = dom%R2_1(dyz, i, j, k, i1)
                            outputs%R2_1_dZ(2,idx) = dom%R2_1(dzz, i, j, k, i1)
                        case (2)
                            outputs%R1_2(:,idx)    = dom%R1_1(  :, i, j, k, i1)
                            outputs%R2_2_dX(0,idx) = dom%R2_1(dxx, i, j, k, i1)
                            outputs%R2_2_dX(1,idx) = dom%R2_1(dyx, i, j, k, i1)
                            outputs%R2_2_dX(2,idx) = dom%R2_1(dzx, i, j, k, i1)
                            outputs%R2_2_dY(0,idx) = dom%R2_1(dxy, i, j, k, i1)
                            outputs%R2_2_dY(1,idx) = dom%R2_1(dyy, i, j, k, i1)
                            outputs%R2_2_dY(2,idx) = dom%R2_1(dzy, i, j, k, i1)
                            outputs%R2_2_dZ(0,idx) = dom%R2_1(dxz, i, j, k, i1)
                            outputs%R2_2_dZ(1,idx) = dom%R2_1(dyz, i, j, k, i1)
                            outputs%R2_2_dZ(2,idx) = dom%R2_1(dzz, i, j, k, i1)
                        end select
                    end if
                    if (i2/=-1) then
                        outputs%R1_2(:,idx)    = dom%R1_2(  :, i, j, k, i2)
                        outputs%R2_2_dX(0,idx) = dom%R2_2(dxx, i, j, k, i2)
                        outputs%R2_2_dX(1,idx) = dom%R2_2(dyx, i, j, k, i2)
                        outputs%R2_2_dX(2,idx) = dom%R2_2(dzx, i, j, k, i2)
                        outputs%R2_2_dY(0,idx) = dom%R2_2(dxy, i, j, k, i2)
                        outputs%R2_2_dY(1,idx) = dom%R2_2(dyy, i, j, k, i2)
                        outputs%R2_2_dY(2,idx) = dom%R2_2(dzy, i, j, k, i2)
                        outputs%R2_2_dZ(0,idx) = dom%R2_2(dxz, i, j, k, i2)
                        outputs%R2_2_dZ(1,idx) = dom%R2_2(dyz, i, j, k, i2)
                        outputs%R2_2_dZ(2,idx) = dom%R2_2(dzz, i, j, k, i2)
                    end if
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    if(allocated(dom%FDump)) outputs%FDump(0:2, idx) = dom%FDump(ind, 0:2)
                    if(allocated(dom%FMasU)) outputs%FMasU(0:2, idx) = dom%FMasU(ind, 0:2)
                    if(allocated(dom%Fint))  outputs%Fint (0:2, idx) = dom%Fint (ind, 0:2)
                end do
            end do
        end do
    end subroutine get_solidpml_rfields

    subroutine init_domain_solidpml(Tdomain, dom)
        use pml
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solidpml), intent(inout) :: dom
        !
        real(fpp) :: fmax


        call copy_cpml_coordinates(Tdomain, dom, DM_SOLID_CG_PML)

        ! Store dt for ADE equations
        dom%dt = Tdomain%TimeD%dtmin

        ! Handle on materials (to get/access pml_pos and pml_width)
        dom%sSubDomain => Tdomain%sSubDomain

        ! Compute alphamax (from fmax)
        fmax = Tdomain%TimeD%fmax
        if (fmax < 0.) stop "SolidCPML : fmax < 0."
        dom%alphamax = M_PI * fmax
    end subroutine init_domain_solidpml

    subroutine start_domain_solidpml(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solidpml), intent(inout) :: dom
        !
        integer :: i, ns

        !$acc  enter data copyin(dom, dom%champs) &
        !$acc&
        do i = 0,1
            !$acc enter data  copyin(dom%champs(i)%Depla, dom%champs(i)%Forces)
        end do

    end subroutine start_domain_solidpml

    subroutine stop_domain_solidpml(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solidpml), intent(inout) :: dom
        !
        integer :: i, ns

        !$acc  exit data delete(dom, dom%champs) &
        !$acc&
        do i = 0,1
            !$acc exit data  delete(dom%champs(i)%Depla, dom%champs(i)%Forces)
        end do

    end subroutine stop_domain_solidpml

    subroutine init_material_properties_solidpml(dom, lnum, mat, density, lambda, mu)
        type(domain_solidpml), intent(inout) :: dom
        integer, intent(in) :: lnum
        type (subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: density
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: lambda
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: mu
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        dom%Density_(:,:,:,bnum,ee) = density
        dom%Lambda_ (:,:,:,bnum,ee) = lambda
        dom%Mu_     (:,:,:,bnum,ee) = mu
    end subroutine init_material_properties_solidpml

    subroutine compute_dxi_alpha_kappa(dom, xyz, i, j, k, bnum, ee, mi, alpha, kappa, dxi)
        type(domain_solidpml), intent(inout) :: dom
        integer, intent(in) :: xyz
        integer, intent(in) :: i, j, k
        integer, intent(in) :: bnum, ee
        integer, intent(in) :: mi
        real(fpp), intent(out) :: alpha, kappa, dxi
        !
        real(fpp) :: xi, xoverl, d0, pspeed, pml_width
        integer :: lnum

        ! d0: (75) from Ref1
        ! dxi: (74) from Ref1
        pml_width = abs(dom%sSubDomain(mi)%pml_width(xyz))
        xi = abs(dom%GlobCoord(xyz,dom%Idom_(i,j,k,bnum,ee)) - dom%sSubDomain(mi)%pml_pos(xyz))
        xoverl = xi/pml_width
        if (xoverl > 1d0) xoverl = 1d0
        if (xoverl < 0d0) xoverl = 0d0

        lnum = bnum*VCHUNK+ee
        pspeed = solidpml_pspeed(dom,lnum,i,j,k)
        d0 = -1.*(dom%cpml_n+1)*Pspeed*log(dom%cpml_rc)
        d0 = d0/(2*pml_width)

        kappa = dom%cpml_kappa_0 + dom%cpml_kappa_1 * xoverl

        dxi   = dom%cpml_c*d0*(xoverl)**dom%cpml_n / kappa
        alpha = dom%alphamax*(1. - xoverl) ! alpha*: (76) from Ref1
    end subroutine compute_dxi_alpha_kappa

    ! Compute parameters for the first direction of attenuation (maybe the only one)
    subroutine compute_dxi_alpha_kappa_dir0(dom, xyz, i, j, k, bnum, ee, mi)
        type(domain_solidpml), intent(inout) :: dom
        integer :: xyz
        integer :: i, j, k
        integer :: bnum, ee
        integer :: mi
        !
        real(fpp) :: alpha, kappa, dxi

        call compute_dxi_alpha_kappa(dom, xyz, i, j, k, bnum, ee, mi, alpha, kappa, dxi)

        dom%D0(ee,bnum) = xyz
        dom%Kappa_0(ee,i,j,k,bnum) = kappa
        dom%dxi_k_0(ee,i,j,k,bnum) = dxi
        dom%Alpha_0(ee,i,j,k,bnum) = alpha
    end subroutine compute_dxi_alpha_kappa_dir0

    ! Compute parameters for the second direction of attenuation
    subroutine compute_dxi_alpha_kappa_dir1(dom, xyz, i, j, k, bnum, ee, mi)
        type(domain_solidpml), intent(inout) :: dom
        integer :: xyz
        integer :: i, j, k
        integer :: bnum, ee
        integer :: mi
        !
        real(fpp) :: alpha, kappa, dxi
        integer :: nd

        nd = get_dir1_index(dom, ee, bnum)
        call compute_dxi_alpha_kappa(dom, xyz, i, j, k, bnum, ee, mi, alpha, kappa, dxi)

        if (kappa==0d0) stop 1
        dom%D1(ee,bnum) = xyz
        dom%Kappa_1(i,j,k,nd) = kappa
        dom%dxi_k_1(i,j,k,nd) = dxi
        dom%Alpha_1(i,j,k,nd) = alpha
    end subroutine compute_dxi_alpha_kappa_dir1

    ! Compute parameters for the third direction of attenuation
    subroutine compute_dxi_alpha_kappa_dir2(dom, xyz, i, j, k, bnum, ee, mi)
        type(domain_solidpml), intent(inout) :: dom
        integer :: xyz
        integer :: i, j, k
        integer :: bnum, ee
        integer :: mi
        !
        real(fpp) :: alpha, kappa, dxi
        integer :: nd

        nd = get_dir2_index(dom, ee, bnum)
        call compute_dxi_alpha_kappa(dom, xyz, i, j, k, bnum, ee, mi, alpha, kappa, dxi)

        if (kappa==0d0) stop 1

        dom%Kappa_2(i,j,k,nd) = kappa
        dom%dxi_k_2(i,j,k,nd) = dxi
        dom%Alpha_2(i,j,k,nd) = alpha
    end subroutine compute_dxi_alpha_kappa_dir2

    function get_dir1_index(dom, ee, bnum) result(nd)
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: ee, bnum
        !
        integer :: nd
        nd = dom%I1(ee,bnum)
        if (nd==-1) then
            nd = dom%dir1_count
            dom%I1(ee,bnum)=nd
            dom%dir1_count = nd+1
        endif
    end function get_dir1_index
    !
    function get_dir2_index(dom, ee, bnum) result(nd)
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: ee, bnum
        !
        integer :: nd
        nd = dom%I2(ee,bnum)
        if (nd==-1) then
            nd = dom%dir2_count
            dom%I2(ee,bnum)=nd
            dom%dir2_count = nd+1
        endif
    end function get_dir2_index
    ! TODO : renommer init_local_mass_solidpml... en init_global_mass_solidpml ? Vu qu'on y met a jour la masse globale !?
    !        attention, ceci impacte le build (compatibilité avec SolidPML)
    subroutine init_local_mass_solidpml(dom,specel,i,j,k,ind,Whei)
        use pml, only : cpml_only_one_root
        type(domain_solidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real(fpp) :: Whei
        !
        integer :: bnum, ee
        real(fpp) :: k0, k1, k2, a0, a1, a2, d0, d1, d2
        real(fpp) :: a0b, a1b, a2b ! solidcpml_a0b_a1b_a2b
        real(fpp) :: mass_0
        integer :: mi, ndir, i1, i2

        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ind = dom%Idom_(i,j,k,bnum,ee)
        mi = specel%mat_index
        ndir = 1
        if (.not. dom%sSubDomain(mi)%dom == DM_SOLID_CG_PML) &
            stop "init_geometric_properties_solidpml : material is not a PML material"

        ! Compute alpha, beta, kappa
        if (dom%sSubDomain(mi)%pml_width(0)/=0) then
            ! DIR X
            call compute_dxi_alpha_kappa_dir0(dom, 0, i, j, k, bnum, ee, mi)
            if (dom%sSubDomain(mi)%pml_width(1)/=0) then
                ! DIR X+Y
                call compute_dxi_alpha_kappa_dir1(dom, 1, i, j, k, bnum, ee, mi)
                if (dom%sSubDomain(mi)%pml_width(2)/=0) then
                    ! DIR X+Y+Z
                    call compute_dxi_alpha_kappa_dir2(dom, 2, i, j, k, bnum, ee, mi)
                    ndir = 3
                else
                    ndir = 2
                endif
            else if (dom%sSubDomain(mi)%pml_width(2)/=0) then
                ! DIR X+Z
                call compute_dxi_alpha_kappa_dir1(dom, 2, i, j, k, bnum, ee, mi)
                ndir = 2
            endif
        else if (dom%sSubDomain(mi)%pml_width(1)/=0) then
            ! DIR Y
            call compute_dxi_alpha_kappa_dir0(dom, 1, i, j, k, bnum, ee, mi)
            if (dom%sSubDomain(mi)%pml_width(2)/=0) then
                ! DIR Y+Z
                call compute_dxi_alpha_kappa_dir1(dom, 2, i, j, k, bnum, ee, mi)
                ndir = 2
            else
                ndir = 1
            endif
        else if (dom%sSubDomain(mi)%pml_width(2)/=0) then
            ! DIR Z
            call compute_dxi_alpha_kappa_dir0(dom, 2, i, j, k, bnum, ee, mi)
            ndir = 1
        else
            stop 1
        endif
        call cpml_only_one_root(dom, i, j, k, bnum, ee, ndir)

        ! gamma_ab  defined after (12c) in Ref1
        ! gamma_abc defined after (12c) in Ref1
        ! Delta 2d derivative term from L : (12a) or (14a) from Ref1
        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        if (k0==0d0) stop 1
        select case(ndir)
        case (1)
            a0b = k0
            a1b = k0*d0
            a2b = -k0*d0*a0
        case (2)
            i1 = dom%I1(ee,bnum)
            k1 = dom%Kappa_1(i,j,k,i1)
            a1 = dom%Alpha_1(i,j,k,i1)
            d1 = dom%dxi_k_1(i,j,k,i1)
            a0b = k0*k1
            a1b = a0b*(d0+d1)
            a2b = a0b*(d0*(d1-a0) - d1*a1)
        case (3)
            i1 = dom%I1(ee,bnum)
            k1 = dom%Kappa_1(i,j,k,i1)
            a1 = dom%Alpha_1(i,j,k,i1)
            d1 = dom%dxi_k_1(i,j,k,i1)
            i2 = dom%I2(ee,bnum)
            k2 = dom%Kappa_2(i,j,k,i2)
            a2 = dom%Alpha_2(i,j,k,i2)
            d2 = dom%dxi_k_2(i,j,k,i2)
            ! Delta term from L : (12a) or (14a) from Ref1
            a0b = k0*k1*k2
            a1b = a0b*(d0+d1+d2)
            a2b = a0b*(d0*(d1-a0) + d1*(d2-a1) + d2*(d0-a2))
        case default
            a0b = 0
            a1b = 0
            a2b = 0
        end select
        mass_0 = Whei*dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)
        ! Delta 1st derivative term from L : (12a) or (14a) from Ref1
        if (a0b==0d0.or.mass_0==0d0) then
            stop 1
        endif

        dom%MassMat(ind) = dom%MassMat(ind) + (a0b+0.5d0*dom%dt*a1b)*mass_0
        dom%DumpMat(ind) = dom%DumpMat(ind) + a1b*mass_0
        dom%MasUMat(ind) = dom%MasUMat(ind) + a2b*mass_0
    end subroutine init_local_mass_solidpml

    subroutine pred_sol_pml(dom, dt, champs1, bnum)
        implicit none

        type(domain_solidpml), intent(inout) :: dom
        type(champssolidpml), intent(inout) :: champs1
        real(fpp), intent(in) :: dt
        integer :: bnum
        !
        ! Useless, kept for compatibility with SolidPML (build), can be deleted later on. TODO : kill this method.
    end subroutine pred_sol_pml

    subroutine forces_int_sol_pml_mainloop(dom, i1)
        type(domain_solidpml), intent(inout) :: dom
        integer, intent(in) :: i1
        integer :: bnum
        !
        do bnum = 0,dom%nblocks-1
            call forces_int_sol_pml(dom, dom%champs(i1), bnum)
        end do
    end subroutine forces_int_sol_pml_mainloop

    subroutine forces_int_sol_pml(dom, champs1, bnum)
!        use sdomain
        use m_calcul_forces_solidpml
        type(domain_solidpml), intent(inout) :: dom
        type(champssolidpml), intent(inout) :: champs1
        integer :: bnum
        !
        integer :: ngll,i,j,k,i_dir,e,ee,idx
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fox,Foy,Foz
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: Depla
        real(fpp) :: wk, wjk, wijk, kijk
        real(fpp), dimension(0:2) :: R

        ngll = dom%ngll

        do i_dir = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        do ee = 0, VCHUNK-1
                            idx = dom%Idom_(i,j,k,bnum,ee)
                            Depla(ee,i,j,k,i_dir) = champs1%Depla(idx,i_dir)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        Fox = 0d0
        Foy = 0d0
        Foz = 0d0
        call calcul_forces_solidpml(dom,bnum,Fox,Foy,Foz,Depla)

        do k = 0,ngll-1
            wk = dom%gllw(k)
            do j = 0,ngll-1
                wjk = wk*dom%gllw(j)
                do i = 0,ngll-1
                    wijk = wjk*dom%gllw(i)
                    do ee = 0, VCHUNK-1
                        e = bnum*VCHUNK+ee
                        if (e>=dom%nbelem) exit
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        kijk = wijk*dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)
                        call compute_L_convolution_terms(dom, i, j, k, bnum, ee, Depla(ee,i,j,k,0:2), R)
                        champs1%Forces(idx,0) = champs1%Forces(idx,0)-Fox(ee,i,j,k)-kijk*R(0)
                        champs1%Forces(idx,1) = champs1%Forces(idx,1)-Foy(ee,i,j,k)-kijk*R(1)
                        champs1%Forces(idx,2) = champs1%Forces(idx,2)-Foz(ee,i,j,k)-kijk*R(2)
                    enddo
                enddo
            enddo
        enddo
    end subroutine forces_int_sol_pml

    subroutine init_solidpml_properties(Tdomain,specel,mat)
        type (domain), intent (INOUT), target :: Tdomain
        type (element), intent(inout) :: specel
        type (subdomain), intent(in) :: mat
        !
        ! Useless, kept for compatibility with SolidPML (build), can be deleted later on. TODO : kill this method.
    end subroutine init_solidpml_properties

    subroutine finalize_solidpml_properties(Tdomain, dom)
        type (domain), intent (INOUT), target :: Tdomain
        type (domain_solidpml), intent (INOUT), target :: dom
        !
        call init_solid_fluid_coupling(Tdomain, dom)
    end subroutine finalize_solidpml_properties

    subroutine init_solid_fluid_coupling(Tdomain, dom)
        use mrenumber, only : get_surface_numbering
        implicit none
        type (domain), intent (INOUT), target :: Tdomain
        type (domain_solidpml), intent (INOUT), target :: dom
        !
        integer, dimension(:), allocatable :: renum ! renum(idom) gives n such that surf0%map(n) = idom
        integer :: ngll, el, i, j, k, idxsf, ee, bnum

        ngll = dom%ngll
        call get_surface_numbering(Tdomain, Tdomain%SF%intSolFluPML%surf0, DM_SOLID_CG_PML, renum)
        do el=0,Tdomain%n_elem-1
            if (Tdomain%specel(el)%domain /= DM_SOLID_CG_PML) cycle
            do i=0,ngll-1
                do j=0,ngll-1
                    do k=0,ngll-1
                        idxsf = renum(Tdomain%specel(el)%Idom(i,j,k))
                        if (idxsf==-1) cycle
                        bnum = Tdomain%specel(el)%lnum/VCHUNK
                        ee = mod(Tdomain%specel(el)%lnum,VCHUNK)
                        call setup_dombase_cpml(dom, idxsf, i, j, k, ee, bnum)
                        ! Better :
                        ! compute directly here b0bar, b1bar, b2bar, alpha0, alpha1
                    end do
                end do
            end do
        end do
        deallocate(renum)
    end subroutine init_solid_fluid_coupling

    subroutine select_terms(a1,a2,a3,t,t1,t2,t3)
        real(fpp), intent(in)  :: a1, a2, a3, t
        real(fpp), intent(out) :: t1, t2, t3
        t1 = 1d0
        if (a2==a1) then
            t2 = t
        else
            t2 = 1d0
        end if
        if (a3==a1) then
            if (a3==a2) then
                t3 = t*t
            else
                t3 = t
            end if
        else
            if (a3==a2) then
                t3 = t
            else
                t3 = 1d0
            end if
        end if
    end subroutine select_terms

    subroutine newmark_predictor_solidpml(dom, Tdomain, f0, f1)
        type(domain_solidpml), intent (INOUT) :: dom
        type (domain), intent (INOUT) :: Tdomain
        integer :: f0, f1
        !
        integer :: n, indpml, indsol

        ! Reset forces
        dom%champs(f1)%Forces = 0d0

        ! Coupling at solid PML interface
        do n = 0,Tdomain%intSolPml%surf0%nbtot-1
            indsol = Tdomain%intSolPml%surf0%map(n)
            indpml = Tdomain%intSolPml%surf1%map(n)
            dom%champs(f0)%Veloc(indpml,:) = Tdomain%sdom%champs(f0)%Veloc(indsol,:)
            dom%champs(f0)%Depla(indpml,:) = Tdomain%sdom%champs(f0)%Depla(indsol,:)
        enddo

        ! The prediction will be based on the current state
        dom%champs(f1)%Depla = dom%champs(f0)%Depla
        dom%champs(f1)%Veloc = dom%champs(f0)%Veloc
    end subroutine newmark_predictor_solidpml

    subroutine newmark_corrector_solidpml(dom, dt, t, f0, f1)
        type(domain_solidpml), intent (INOUT) :: dom
        double precision :: dt, t
        integer :: f0, f1
        !
        integer :: i_dir, n, indpml
        !
        real(fpp) :: V, F ! Displacement, Velocity, Force

        ! Note: the solid domain use a leap-frop time scheme => we get D_n+3/2 and V_n+1
        !       to compute V_n+2, we need F_n+1. To compute F_n+1, we need D_n+1 and V_n+1
        !       D_n+1 is estimated with D_n+3/2 and V_n+1

        ! Update velocity: compute V_n+2
        do i_dir = 0,2
            do n = 0,dom%nglltot-1
                ! Get V = V_n+1
                V = dom%champs(f0)%Veloc(n,i_dir) ! V = V_n+1

!                ! Estimate D = D_n+1
!                D = dom%champs(f0)%Depla(n,i_dir) - V*0.5*dt

                ! Save forces contributions for snapshots
                if(allocated(dom%FDump)) dom%FDump(n, i_dir) = - dom%DumpMat(n)*V
                if(allocated(dom%FMasU)) dom%FMasU(n, i_dir) = - dom%MasUMat(n)*dom%champs(f0)%Depla(n,i_dir)
                if(allocated(dom%Fint))  dom%Fint (n, i_dir) =   dom%champs(f1)%Forces(n,i_dir)

                ! Compute F_n+1 : (61a) from Ref1 with F = -F (as we add -Fo* in forces_int_sol_pml)
                F =    &
                    dom%MassMat(n)*(- dom%DumpMat(n)*V - dom%MasUMat(n)*dom%champs(f0)%Depla(n,i_dir)  &
                    + dom%champs(f1)%Forces(n,i_dir))

                ! Compute V_n+2
                ! since dom%MassMat = 1./dom%MassMat (define_arrays inverse_mass_mat)
                dom%champs(f0)%Veloc(n,i_dir) = V + dt*F

                ! Update displacement: compute D_n+5/2
                dom%champs(f0)%Depla(n,i_dir) = dom%champs(f0)%Depla(n,i_dir) + dt * dom%champs(f0)%Veloc(n,i_dir)
            end do
        end do


        ! Note: do NOT apply (dirichlet) BC for PML
        !       if PML absorption would be turned off <=> solid domain without dirichlet BC (neumann only)
        ! yes we can
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs(f0)%Veloc(indpml,:) = 0.
            dom%champs(f0)%Depla(indpml,:) = 0.
        enddo

    end subroutine newmark_corrector_solidpml

    function solidpml_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        type(domain_solidpml), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: lbd, mu, rho, Pspeed
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        lbd = dom%Lambda_(i,j,k,bnum,ee)
        mu = dom%Mu_(i,j,k,bnum,ee)
        rho = dom%Density_(i,j,k,bnum,ee)
        Pspeed = sqrt((lbd+2.*mu)/rho)
    end function solidpml_Pspeed

    subroutine couplage_pml_solid(Tdomain, sdom, spmldom, i0, i1)
        use dom_solid
        implicit none
        type(domain), intent(inout)  :: Tdomain
        type(domain_solid), intent (inout) :: sdom
        type(domain_solidpml), intent (inout) :: spmldom
        integer, intent(in) :: i0, i1
        !
        integer :: n, indsol, indpml, i
        real(fpp) :: force
        !
        do n = 0,Tdomain%intSolPml%surf0%nbtot-1
            indsol = Tdomain%intSolPml%surf0%map(n)
            indpml = Tdomain%intSolPml%surf1%map(n)
            do i=0,2
                force = spmldom%champs(i1)%Forces(indpml,i) &
                    - spmldom%DumpMat(indpml)*spmldom%champs(i0)%Veloc(indpml,i) &
                    - spmldom%MasUMat(indpml)*spmldom%champs(i0)%Depla(indpml,i)
                sdom%champs(i1)%Veloc(indsol,i) = sdom%champs(i1)%Veloc(indsol,i) + force
            enddo
        enddo

    end subroutine couplage_pml_solid

end module dom_solidpml

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
