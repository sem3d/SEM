!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_fluidpml
    use sdomain
    use constants
    use champs_fluidpml
    use selement
    use ssubdomains
    implicit none
#include "index.h"

contains

    subroutine allocate_champs_fluidcpml(dom, i)
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: i

        ! Allocate ONE more gll than necessary to use as a dummy target for
        ! indirections for fake elements.
        allocate(dom%champs(i)%ForcesFl(0:dom%nglltot-1))
        allocate(dom%champs(i)%Phi     (0:dom%nglltot-1))
        allocate(dom%champs(i)%VelPhi  (0:dom%nglltot-1))

        dom%champs(i)%ForcesFl = 0d0
        dom%champs(i)%Phi = 0d0
        dom%champs(i)%VelPhi = 0d0
    end subroutine allocate_champs_fluidcpml

    subroutine allocate_multi_dir_pml(Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_fluidpml), intent (INOUT) :: dom
        !
        integer :: dir1_count, dir2_count, ndir, mi, n, ngll
        dir1_count = 0
        dir2_count = 0
        ngll = dom%ngll
        do n=0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%domain/=DM_SOLID_PML) cycle
            ndir = 0
            mi = Tdomain%specel(n)%mat_index
            if (Tdomain%sSubDomain(mi)%pml_width(0)/=0d0) ndir = ndir + 1
            if (Tdomain%sSubDomain(mi)%pml_width(1)/=0d0) ndir = ndir + 1
            if (Tdomain%sSubDomain(mi)%pml_width(2)/=0d0) ndir = ndir + 1
            if (ndir>=2) dir1_count = dir1_count + 1
            if (ndir>=3) dir2_count = dir2_count + 1
        end do
        if (dir1_count>0) then
            allocate(dom%Alpha_1(0:ngll-1,0:ngll-1,0:ngll-1,0:dir1_count-1))
            allocate(dom%Kappa_1(0:ngll-1,0:ngll-1,0:ngll-1,0:dir1_count-1))
            allocate(dom%dxi_k_1(0:ngll-1,0:ngll-1,0:ngll-1,0:dir1_count-1))

            allocate(dom%R1_1(0:ngll-1, 0:ngll-1, 0:ngll-1, 0:dir1_count-1))
            allocate(dom%R2_1(0:3, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:dir1_count-1))
            dom%R1_1 = 0d0
            dom%R2_1 = 0d0
        endif
        if (dir2_count>0) then
            allocate(dom%Alpha_2(0:ngll-1,0:ngll-1,0:ngll-1,0:dir2_count-1))
            allocate(dom%Kappa_2(0:ngll-1,0:ngll-1,0:ngll-1,0:dir2_count-1))
            allocate(dom%dxi_k_2(0:ngll-1,0:ngll-1,0:ngll-1,0:dir2_count-1))

            allocate(dom%R1_2(0:ngll-1, 0:ngll-1, 0:ngll-1, 0:dir2_count-1))
            allocate(dom%R2_2(0:3, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:dir2_count-1))
            dom%R1_2 = 0d0
            dom%R2_2 = 0d0
        end if
#ifdef DBG
        ! Write infos for all procs as all procs have different informations !... A bit messy output but no other way
        write(*,*) "INFO - fluid cpml domain : ", dir1_count, " elems attenuated in 2 directions on proc", Tdomain%rank
        write(*,*) "INFO - fluid cpml domain : ", dir2_count, " elems attenuated in 3 directions on proc", Tdomain%rank
#endif
    end subroutine allocate_multi_dir_pml

    subroutine allocate_dom_fluidpml (Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_fluidpml), intent (INOUT) :: dom
        !
        integer :: nbelem, ngll, nblocks
        !

        ngll   = dom%ngll
        nbelem = dom%nbelem
        write(*,*) "DOM_fluidpml ngll   = ", ngll
        write(*,*) "DOM_fluidpml nbelem = ", nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call init_dombase(dom)

        call allocate_multi_dir_pml(Tdomain, dom)
        ! We reset the count here, since we will use them to renumber the specel that
        ! have more than one direction of attenuation
        dom%dir1_count = 0
        dom%dir2_count = 0

        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)
            nblocks = dom%nblocks
            allocate(dom%IDensity_(0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))

            allocate(dom%I1      (0:VCHUNK-1, 0:nblocks-1))
            allocate(dom%I2      (0:VCHUNK-1, 0:nblocks-1))
            allocate(dom%D0      (0:VCHUNK-1, 0:nblocks-1))
            allocate(dom%D1      (0:VCHUNK-1, 0:nblocks-1))
            allocate(dom%Alpha_0 (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            allocate(dom%dxi_k_0 (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            allocate(dom%Kappa_0 (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            dom%I1(:,:) = -1
            dom%I2(:,:) = -1
            dom%D0(:,:) = 0
            dom%D1(:,:) = 0
            dom%Kappa_0 = 1.

            allocate(dom%R1_0(0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            allocate(dom%R2_0(0:VCHUNK-1, 0:3, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            dom%R1_0 = 0d0
            dom%R2_0 = 0d0

            allocate(dom%PhiOld(0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            dom%PhiOld = 0.
            allocate(dom%DPhiOld(0:VCHUNK-1, 0:2, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            dom%DPhiOld = 0.
        end if
        ! Allocation et initialisation de champs0 et champs1 pour les fluides
        if (dom%nglltot /= 0) then
            call allocate_champs_fluidcpml(dom, 0)
            call allocate_champs_fluidcpml(dom, 1)

            ! Allocation de DumpMat pour les PML solides
            allocate(dom%DumpMat(0:dom%nglltot-1))
            dom%DumpMat = 0d0

            ! Allocation de MasUMat pour les PML solides
            allocate(dom%MasUMat(0:dom%nglltot-1))
            dom%MasUMat = 0d0
        endif

        ! CPML parameters initialisation: for the very first implementation, parameters are hard-coded.

        dom%cpml_c = 1.0
        dom%cpml_n = Tdomain%config%cpml_n
        dom%cpml_rc = Tdomain%config%cpml_rc
        dom%cpml_kappa_0 = Tdomain%config%cpml_kappa0
        dom%cpml_kappa_1 = Tdomain%config%cpml_kappa1
        dom%alphamax = 0.
        if(Tdomain%rank==0) write(*,*) "INFO - fluidpml domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_fluidpml

    subroutine deallocate_dom_fluidpml (dom)
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        !
        integer :: i

        if(allocated(dom%m_IDensity)) deallocate(dom%m_IDensity)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )

        if(allocated(dom%Alpha_0)) deallocate(dom%Alpha_0)
        if(allocated(dom%dxi_k_0)) deallocate(dom%dxi_k_0)
        if(allocated(dom%Kappa_0)) deallocate(dom%Kappa_0)

        if(allocated(dom%Alpha_1)) deallocate(dom%Alpha_1)
        if(allocated(dom%dxi_k_1)) deallocate(dom%dxi_k_1)
        if(allocated(dom%Kappa_1)) deallocate(dom%Kappa_1)

        if(allocated(dom%Alpha_2)) deallocate(dom%Alpha_2)
        if(allocated(dom%dxi_k_2)) deallocate(dom%dxi_k_2)
        if(allocated(dom%Kappa_2)) deallocate(dom%Kappa_2)

        if(allocated(dom%I1)) deallocate(dom%I1)
        if(allocated(dom%I2)) deallocate(dom%I2)
        if(allocated(dom%D0)) deallocate(dom%D0)
        if(allocated(dom%D1)) deallocate(dom%D1)

        do i=0,1
            if(allocated(dom%champs(i)%ForcesFl)) deallocate(dom%champs(i)%ForcesFl)
            if(allocated(dom%champs(i)%Phi     )) deallocate(dom%champs(i)%Phi     )
            if(allocated(dom%champs(i)%VelPhi  )) deallocate(dom%champs(i)%VelPhi  )
        end do

        if(allocated(dom%DumpMat)) deallocate(dom%DumpMat)
        if(allocated(dom%MasUMat)) deallocate(dom%MasUMat)

        call deallocate_dombase(dom)
    end subroutine deallocate_dom_fluidpml

    subroutine init_domain_fluidpml(Tdomain, dom)
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_fluidpml), intent(inout) :: dom
        !
        real(fpp) :: fmax
        integer :: i, j, k, n, indL, indG

        ! Handle on node global coords : mandatory to compute distances in the PML (compute_dxi_alpha_kappa)
        ! TODO precompute usefull coeffs instead of copying coords...
        allocate(dom%GlobCoord(0:2,0:dom%nglltot-1))
        do n=0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%domain/=DM_FLUID_PML) cycle
            do k = 0,dom%ngll-1
                do j = 0,dom%ngll-1
                    do i = 0,dom%ngll-1
                        indG = Tdomain%specel(n)%Iglobnum(i,j,k)
                        indL = Tdomain%specel(n)%Idom(i,j,k)
                        dom%GlobCoord(:,indL) = Tdomain%GlobCoord(:,indG)
                    end do
                end do
            end do
        end do

        ! Store dt for ADE equations
        dom%dt = Tdomain%TimeD%dtmin

        ! Handle on materials (to get/access pml_pos and pml_width)
        dom%sSubDomain => Tdomain%sSubDomain

        ! Compute alphamax (from fmax)
        fmax = Tdomain%TimeD%fmax
        if (fmax < 0.) stop "FluidCPML : fmax < 0."
        dom%alphamax = M_PI * fmax
    end subroutine init_domain_fluidpml

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

    subroutine get_fluidpml_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain)                               :: TDomain
        type(domain_fluidpml), intent(inout)          :: dom
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
                        phi(i,j,k) = dom%champs(0)%Phi(ind)
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
                        vphi(i,j,k) = dom%champs(0)%VelPhi(ind)
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
                        fieldP(i,j,k) = -dom%champs(0)%VelPhi(ind)
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
    end subroutine get_fluidpml_dom_var


    subroutine get_fluidpml_dom_elem_energy(dom, lnum, P_energy, S_energy)
        use deriv3d
        implicit none
        !
        type(domain_fluidpml), intent(inout)          :: dom
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

    end subroutine get_fluidpml_dom_elem_energy

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

        dom%IDensity_(:,:,:,bnum,ee) = 1d0/density
        dom%Lambda_ (:,:,:,bnum,ee) = lambda
    end subroutine init_material_properties_fluidpml

    subroutine compute_dxi_alpha_kappa(dom, xyz, i, j, k, bnum, ee, mi, alpha, kappa, dxi)
        type(domain_fluidpml), intent(inout) :: dom
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
        pspeed = fluidpml_pspeed(dom,lnum,i,j,k)
        d0 = -1.*(dom%cpml_n+1)*Pspeed*log(dom%cpml_rc)
        d0 = d0/(2*pml_width)

        kappa = dom%cpml_kappa_0 + dom%cpml_kappa_1 * xoverl

        dxi   = dom%cpml_c*d0*(xoverl)**dom%cpml_n / kappa
        alpha = dom%alphamax*(1. - xoverl) ! alpha*: (76) from Ref1
    end subroutine compute_dxi_alpha_kappa

    ! Compute parameters for the first direction of attenuation (maybe the only one)
    subroutine compute_dxi_alpha_kappa_dir0(dom, xyz, i, j, k, bnum, ee, mi)
        type(domain_fluidpml), intent(inout) :: dom
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
        type(domain_fluidpml), intent(inout) :: dom
        integer :: xyz
        integer :: i, j, k
        integer :: bnum, ee
        integer :: mi
        !
        real(fpp) :: alpha, kappa, dxi
        integer :: nd

        nd = get_dir1_index(dom, ee, bnum)
        call compute_dxi_alpha_kappa(dom, xyz, i, j, k, bnum, ee, mi, alpha, kappa, dxi)
        dom%D1(ee,bnum) = xyz
        dom%Kappa_1(i,j,k,nd) = kappa
        dom%dxi_k_1(i,j,k,nd) = dxi
        dom%Alpha_1(i,j,k,nd) = alpha
    end subroutine compute_dxi_alpha_kappa_dir1

    ! Compute parameters for the third direction of attenuation
    subroutine compute_dxi_alpha_kappa_dir2(dom, xyz, i, j, k, bnum, ee, mi)
        type(domain_fluidpml), intent(inout) :: dom
        integer :: xyz
        integer :: i, j, k
        integer :: bnum, ee
        integer :: mi
        !
        real(fpp) :: alpha, kappa, dxi
        integer :: nd

        nd = get_dir2_index(dom, ee, bnum)
        call compute_dxi_alpha_kappa(dom, xyz, i, j, k, bnum, ee, mi, alpha, kappa, dxi)
        dom%Kappa_2(i,j,k,nd) = kappa
        dom%dxi_k_2(i,j,k,nd) = dxi
        dom%Alpha_2(i,j,k,nd) = alpha
    end subroutine compute_dxi_alpha_kappa_dir2

    function get_dir1_index(dom, ee, bnum) result(nd)
        type(domain_fluidpml), intent (INOUT) :: dom
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

    function get_dir2_index(dom, ee, bnum) result(nd)
        type(domain_fluidpml), intent (INOUT) :: dom
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

    subroutine init_local_mass_fluidpml(dom,specel,i,j,k,ind,Whei)
        type(domain_fluidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real(kind=fpp), intent(in) :: Whei
        !
        integer :: bnum, ee
        real(fpp) :: k0, k1, k2, a0, a1, a2, d0, d1, d2
        real(fpp) :: a0b, a1b, a2b
        real(fpp) :: mass_0
        integer :: mi, ndir, i1, i2

        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ind = dom%Idom_(i,j,k,bnum,ee)
        mi = specel%mat_index
        ndir = 1
        if (.not. dom%sSubDomain(mi)%dom == DM_FLUID_PML) &
            stop "init_local_mass_fluidpml : material is not a PML material"

        ! Compute alpha, beta, kappa
        if (dom%sSubDomain(mi)%pml_width(0)/=0) then
            call compute_dxi_alpha_kappa_dir0(dom, 0, i, j, k, bnum, ee, mi)
            if (dom%sSubDomain(mi)%pml_width(1)/=0) then
                call compute_dxi_alpha_kappa_dir1(dom, 1, i, j, k, bnum, ee, mi)
                if (dom%sSubDomain(mi)%pml_width(2)/=0) then
                    call compute_dxi_alpha_kappa_dir2(dom, 2, i, j, k, bnum, ee, mi)
                    ndir = 3
                else
                    ndir = 2
                endif
            else if (dom%sSubDomain(mi)%pml_width(1)/=0) then
                call compute_dxi_alpha_kappa_dir1(dom, 2, i, j, k, bnum, ee, mi)
                ndir = 2
            endif
        else if (dom%sSubDomain(mi)%pml_width(1)/=0) then
            call compute_dxi_alpha_kappa_dir0(dom, 1, i, j, k, bnum, ee, mi)
            if (dom%sSubDomain(mi)%pml_width(2)/=0) then
                call compute_dxi_alpha_kappa_dir1(dom, 2, i, j, k, bnum, ee, mi)
                ndir = 2
            else
                ndir = 1
            endif
        else if (dom%sSubDomain(mi)%pml_width(2)/=0) then
            call compute_dxi_alpha_kappa_dir0(dom, 2, i, j, k, bnum, ee, mi)
            ndir = 1
        else
            stop 1
        endif

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
        end select

        ! Fluid : inertial term ponderation by the inverse of the bulk modulus

        mass_0 = Whei*dom%Jacob_(i,j,k,bnum,ee)/dom%Lambda_(i,j,k,bnum,ee)
        ! Delta 1st derivative term from L : (12a) or (14a) from Ref1
        if (a0b==0d0.or.mass_0==0d0) then
            stop 1
        endif
        dom%MassMat(ind) = dom%MassMat(ind) + a0b*mass_0
        dom%DumpMat(ind) = dom%DumpMat(ind) + a1b*mass_0
        dom%MasUMat(ind) = dom%MasUMat(ind) + a2b*mass_0
    end subroutine init_local_mass_fluidpml

    subroutine forces_int_fluidpml(dom, champs1, bnum)
        use m_calcul_forces_fluidpml
        type(domain_fluidpml), intent (INOUT) :: dom
        type(champsfluid), intent(inout) :: champs1
        integer, intent(in) :: bnum
        !
        integer :: ngll,i,j,k,e,ee,idx
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1, 0:dom%ngll-1, 0:dom%ngll-1) :: Fo_Fl,Phi
        real(fpp) :: val
        real(fpp) :: wk, wjk, wijk, kijk, R

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
        call calcul_forces_fluidpml(dom,dom%ngll,bnum,Fo_Fl,Phi)
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
                        val = champs1%ForcesFl(idx)
                        val = val - Fo_Fl(ee,i,j,k)
                        champs1%ForcesFl(idx) = val
                        call compute_L_convolution_terms(dom, i, j, k, bnum, ee, Phi(ee,i,j,k), R)
                        kijk = wijk*1./dom%IDensity_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)
                        champs1%ForcesFl(idx) = champs1%ForcesFl(idx) - kijk*R
                    enddo
                enddo
            enddo
        enddo
    end subroutine forces_int_fluidpml

    subroutine pred_flu_pml(dom, dt, champs1, bnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        real(fpp), intent(in) :: dt
        type(champsfluidpml), intent(inout) :: champs1
        integer :: bnum
        !
    end subroutine pred_flu_pml

    subroutine finalize_fluidpml_properties(Tdomain,dom)
      type (domain), intent (INOUT), target :: Tdomain
      type (domain_fluidpml), intent (INOUT), target :: dom
      !
      ! TODO : kill (needed for compatibility - build)
    end subroutine finalize_fluidpml_properties

    subroutine init_fluidpml_properties(Tdomain,specel,mat)
        type (domain), intent (INOUT), target :: Tdomain
        type (element), intent(inout) :: specel
        type (subdomain), intent(in) :: mat
        !
      ! TODO : kill (needed for compatibility - build)
    end subroutine init_fluidpml_properties

    subroutine forces_int_flu_pml(dom, champs1, bnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        type(champsfluidpml), intent(inout) :: champs1
        integer :: bnum
        !
      ! TODO : kill (needed for compatibility - build)
    end subroutine forces_int_flu_pml

    subroutine newmark_predictor_fluidpml(dom, Tdomain, f0, f1)
        type(domain_fluidpml), intent (INOUT) :: dom
        type(domain)                          :: TDomain
        integer :: f0, f1
        !
        integer :: n, indpml, indflu

        ! Reset forces
        dom%champs(f1)%ForcesFl = 0d0

        ! Coupling at solid PML interface
        do n = 0,Tdomain%intFluPml%surf0%nbtot-1
            indflu = Tdomain%intFluPml%surf0%map(n)
            indpml = Tdomain%intFluPml%surf1%map(n)
            dom%champs(f1)%VelPhi(indpml) = Tdomain%fdom%champs(f0)%VelPhi(indflu)
            dom%champs(f1)%Phi(indpml) = Tdomain%fdom%champs(f0)%Phi(indflu)
        enddo

        ! The prediction will be based on the current state
        dom%champs(f1)%VelPhi = dom%champs(f0)%VelPhi
        dom%champs(f1)%Phi    = dom%champs(f0)%Phi
    end subroutine newmark_predictor_fluidpml

    subroutine newmark_corrector_fluidpml(dom, dt, f0, f1)
        type(domain_fluidpml), intent (INOUT) :: dom
        double precision :: dt
        integer :: f0, f1
        !
        integer :: n, indpml

        dom%champs(f0)%ForcesFl = (- dom%DumpMat*dom%champs(f0)%VelPhi - dom%MasUMat*dom%champs(f0)%Phi &
                                   + dom%champs(f1)%ForcesFl) * dom%MassMat
        dom%champs(f0)%VelPhi = (dom%champs(f0)%VelPhi + dt * dom%champs(f0)%ForcesFl)
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs(f0)%VelPhi(indpml) = 0. ! Apply (dirichlet) BC for PML
        enddo

        dom%champs(f0)%Phi = dom%champs(f0)%Phi + dt * dom%champs(f0)%VelPhi
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs(f0)%Phi = 0. ! Apply (dirichlet) BC for PML
        enddo
    end subroutine newmark_corrector_fluidpml

    function fluidpml_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        type(domain_fluidpml), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        Pspeed = sqrt(dom%Lambda_(i,j,k,bnum,ee)*dom%IDensity_(i,j,k,bnum,ee))
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
