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

    subroutine allocate_multi_dir_pml(Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_solidpml), intent (INOUT) :: dom
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
        endif
        if (dir2_count>0) then
            allocate(dom%Alpha_2(0:ngll-1,0:ngll-1,0:ngll-1,0:dir2_count-1))
            allocate(dom%Kappa_2(0:ngll-1,0:ngll-1,0:ngll-1,0:dir2_count-1))
            allocate(dom%dxi_k_2(0:ngll-1,0:ngll-1,0:ngll-1,0:dir2_count-1))
        end if
!        write(*,*) "DEBUG:", dir1_count,"Elements 2 dir"
!        write(*,*) "DEBUG:", dir2_count,"Elements 3 dir"
    end subroutine allocate_multi_dir_pml

    subroutine allocate_dom_solidpml (Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_solidpml), intent (INOUT) :: dom
        !
        integer nbelem, ngll, nblocks
        !

        ngll   = dom%ngll
        nbelem = dom%nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call init_dombase(dom)

        call allocate_multi_dir_pml(Tdomain, dom)
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
            allocate(dom%R2_0(0:VCHUNK-1,0:2, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1))
            dom%R1_0 = 0d0
            dom%R2_0 = 0d0
            dom%I1(:,:) = -1
            dom%I2(:,:) = -1
            dom%D0(:,:) = 0
            dom%Kappa_0 = 1.
        end if

        ! Allocation et initialisation de champs0 pour les PML solides
        if (dom%nglltot /= 0) then
            allocate(dom%champs0%Forces(0:dom%nglltot-1,0:2))
            allocate(dom%champs0%Depla (0:dom%nglltot-1,0:2))
            allocate(dom%champs0%Veloc (0:dom%nglltot-1,0:2))
            allocate(dom%champs1%Forces(0:dom%nglltot-1,0:2))
            allocate(dom%champs1%Depla (0:dom%nglltot-1,0:2))
            allocate(dom%champs1%Veloc (0:dom%nglltot-1,0:2))
            dom%champs0%Depla  = 0d0
            dom%champs0%Veloc  = 0d0
            dom%champs0%Forces = 0d0
            dom%champs1%Depla  = 0d0
            dom%champs1%Veloc  = 0d0
            dom%champs1%Forces = 0d0

            ! Allocation de DumpMat pour les PML solides
            allocate(dom%DumpMat(0:dom%nglltot-1))
            dom%DumpMat = 0d0

            ! Allocation de MasUMat pour les PML solides
            allocate(dom%MasUMat(0:dom%nglltot-1))
            dom%MasUMat = 0d0

        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid cpml domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"

        ! CPML parameters initialisation: for the very first implementation, parameters are hard-coded.
        ! TODO : read parameters (kappa_* ?) from input.spec ?

        dom%c(:) = 1.
        dom%n(:) = Tdomain%config%cpml_n
        dom%rc = Tdomain%config%cpml_rc
        dom%cpml_kappa_0 = Tdomain%config%cpml_kappa0
        dom%cpml_kappa_1 = Tdomain%config%cpml_kappa1
        dom%alphamax = 0.
        if(Tdomain%rank==0) then
            write(*,*) "INFO - solid cpml domain : kappa0 ", dom%cpml_kappa_0, " kappa1 ", dom%cpml_kappa_1
        endif
    end subroutine allocate_dom_solidpml

    subroutine deallocate_dom_solidpml (dom)
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_Mu     )) deallocate(dom%m_Mu     )

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

        if(allocated(dom%champs0%Depla )) deallocate(dom%champs0%Depla )
        if(allocated(dom%champs0%Veloc )) deallocate(dom%champs0%Veloc )
        if(allocated(dom%champs0%Forces )) deallocate(dom%champs0%Forces )
        if(allocated(dom%champs1%Depla )) deallocate(dom%champs1%Depla )
        if(allocated(dom%champs1%Veloc )) deallocate(dom%champs1%Veloc )
        if(allocated(dom%champs1%Forces )) deallocate(dom%champs1%Forces )

        if(allocated(dom%DumpMat)) deallocate(dom%DumpMat)
        if(allocated(dom%MasUMat)) deallocate(dom%MasUMat)

        if(allocated(dom%R1_0)) deallocate(dom%R1_0)
        if(allocated(dom%R2_0)) deallocate(dom%R2_0)

        if(allocated(dom%GlobCoord)) deallocate(dom%GlobCoord)
    end subroutine deallocate_dom_solidpml

    subroutine get_solidpml_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain_solidpml)                      :: dom
        integer, dimension(0:), intent(in)         :: out_variables
        integer                                    :: lnum
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        !
        logical                  :: flag_gradU
        integer                  :: ngll, i, j, k, ind
        real(fpp)                :: DXX, DXY, DXZ
        real(fpp)                :: DYX, DYY, DYZ
        real(fpp)                :: DZX, DZY, DZZ
        real, dimension(0:2,0:2) :: invgrad_ijk
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
            if(.not. allocated(fieldU)) allocate(fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        fieldU(i,j,k,:) = dom%champs0%Depla(ind,:)
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
                        if(.not. allocated(fieldV)) allocate(fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldV(i,j,k,:) = dom%champs0%Veloc(ind,:)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldA(i,j,k,:) = dom%Massmat(ind) * dom%champs1%Forces(ind,:)
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
                        fieldP(i,j,k) = 0d0
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1 .or. &
                        out_variables(OUT_ENERGYP) == 1 .or. &
                        out_variables(OUT_EPS_DEV) == 1 .or. &
                        out_variables(OUT_STRESS_DEV) == 1 ) then
                        if(.not. allocated(eps_vol)) allocate(eps_vol(0:ngll-1,0:ngll-1,0:ngll-1))
                        eps_vol(i,j,k) = DXX + DYY + DZZ
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
                        eps_dev(i,j,k,0) = DXX - eps_vol(i,j,k) / 3
                        eps_dev(i,j,k,1) = DYY - eps_vol(i,j,k) / 3
                        eps_dev(i,j,k,2) = DZZ - eps_vol(i,j,k) / 3
                        eps_dev(i,j,k,3) = 0.5 * (DXY + DYX)
                        eps_dev(i,j,k,4) = 0.5 * (DZX + DXZ)
                        eps_dev(i,j,k,5) = 0.5 * (DZY + DYZ)
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        if(.not. allocated(sig_dev)) allocate(sig_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
                        sig_dev(i,j,k,:) = 0.
                    end if
                enddo
            enddo
        enddo
    end subroutine get_solidpml_dom_var

    subroutine init_domain_solidpml(Tdomain, dom)
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solidpml), intent(inout) :: dom
        !
        integer n, bnum, ee
        real(fpp) :: fmax
        integer :: i,j,k, indL, indG
        ! Handle on node global coords : mandatory to compute distances in the PML (compute_dxi_alpha_kappa)
        ! TODO precompute usefull coeffs instead of copying coords...
        allocate(dom%GlobCoord(0:2,0:dom%nglltot-1))
        do n=0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%domain/=DM_SOLID_PML) cycle
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

        dom%dt = Tdomain%TimeD%dtmin
        ! Handle on materials (to get/access pml_pos and pml_width)
        dom%sSubDomain => Tdomain%sSubDomain

        ! Compute alphamax (from fmax)

        ! TODO : compute max of all fmax of all procs: need to speak about that with Ludovic: where ? how ?
        ! My understanding is that all procs process different elements (PML or not), so all procs are NOT
        ! here (in init_solidpml_properties) at the same time : broadcast fmax of each proc to all procs and get
        ! the max could NOT work ?!... Right ? Wrong ? Don't know !...
        ! TODO : replace all this by fmax from read_input.c
        if (Tdomain%nb_procs /= 1) stop "ERROR : SolidCPML is limited to monoproc for now"

        fmax = Tdomain%TimeD%fmax
        if (fmax < 0.) stop "SolidCPML : fmax < 0."
        dom%alphamax = M_PI * fmax
    end subroutine init_domain_solidpml

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
        real(fpp) :: xi, xoverl, d0, pspeed
        integer :: lnum

        ! d0: (75) from Ref1
        ! dxi: (74) from Ref1
        xi = abs(dom%GlobCoord(xyz,dom%Idom_(i,j,k,bnum,ee)) - dom%sSubDomain(mi)%pml_pos(xyz))
        xoverl = xi/abs(dom%sSubDomain(mi)%pml_width(xyz))
        if (xoverl > 1) xoverl = 1d0

        lnum = bnum*VCHUNK+ee
        pspeed = solidpml_pspeed(dom,lnum,i,j,k)
        d0 = -1.*(dom%n(xyz)+1)*Pspeed*log(dom%rc)/log(10d0)
        d0 = d0/abs(2*dom%sSubDomain(mi)%pml_width(xyz))

        kappa = dom%cpml_kappa_0 + dom%cpml_kappa_1 * xoverl

        dxi   = dom%c(xyz)*d0*(xoverl)**dom%n(xyz) / kappa
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
        type(domain_solidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real Whei
        !
        integer :: bnum, ee
        real(fpp) :: k0, k1, k2, a0, a1, a2, d0, d1, d2
        real(fpp) :: a0b, a1b, a2b ! solidcpml_a0b_a1b_a2b
        real(fpp) :: mass_0
        integer :: mi, nd, ndir, i1, i2

        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ind = dom%Idom_(i,j,k,bnum,ee)
        mi = specel%mat_index
        ndir = 1
        if (.not. dom%sSubDomain(mi)%dom == DM_SOLID_PML) &
            stop "init_geometric_properties_solidpml : material is not a PML material"

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
            a2b = k0*d0*a0
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
            a0b = k0*k1*k2
            a1b = a0b*(d0+d1+d2)
            ! Delta term from L : (12a) or (14a) from Ref1
            a2b = a0b*(d0*(d1-a0) + d1*(d2-a1) + d2*(d0-a2))
        end select
        mass_0 = Whei*dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)
        ! Delta 1st derivative term from L : (12a) or (14a) from Ref1
        if (a0b==0d0.or.mass_0==0d0) then
            stop 1
        endif
        dom%MassMat(ind) = dom%MassMat(ind) + a0b*mass_0
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

    subroutine forces_int_sol_pml(dom, champs1, bnum, Tdomain)
        use sdomain
        use m_calcul_forces_solidpml
        type(domain_solidpml), intent(inout) :: dom
        type(champssolidpml), intent(inout) :: champs1
        integer :: bnum
        type (domain), intent (INOUT), target :: Tdomain
        !
        integer :: ngll,i,j,k,i_dir,e,ee,idx
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fox,Foy,Foz
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: Depla
        real(fpp) :: Rx, Ry, Rz, wk, wjk, wijk

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
                        call compute_L_convolution_terms(dom, i, j, k, bnum, ee, Depla(ee,i,j,k,0:2), &
                            Rx, Ry, Rz)
                        champs1%Forces(idx,0) = champs1%Forces(idx,0)-Fox(ee,i,j,k)-wijk*Rx
                        champs1%Forces(idx,1) = champs1%Forces(idx,1)-Foy(ee,i,j,k)-wijk*Ry
                        champs1%Forces(idx,2) = champs1%Forces(idx,2)-Foz(ee,i,j,k)-wijk*Rz
                        !R(idx,0) += wijk*(R1+R2+R3)
                    enddo
                enddo
            enddo
        enddo
        ! Update convolution terms
        call update_convolution_terms(dom, champs1, bnum, Tdomain)
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
        ! Useless, kept for compatibility with SolidPML (build), can be deleted later on. TODO : kill this method.
    end subroutine finalize_solidpml_properties

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

    subroutine update_convolution_terms(dom, champs1, bnum, Tdomain)
        type(domain_solidpml), intent (INOUT) :: dom
        type(champssolidpml), intent(inout) :: champs1
        integer :: bnum
        type (domain), intent (INOUT), target :: Tdomain
        !
        real(fpp) :: a1, a2, a3, ui, t1, t2, t3, t, dt
        integer :: i_dir, i, j, k, ee, idx

!        t  = Tdomain%timeD%rtime
!        dt = Tdomain%timeD%dtmin
!        do i_dir = 0,2
!            do k = 0,dom%ngll-1
!                do j = 0,dom%ngll-1
!                    do i = 0,dom%ngll-1
!                        do ee = 0, VCHUNK-1
!                            idx = dom%Idom_(i,j,k,bnum,ee)
!                            a1 = dom%Alpha(idx,0)
!                            a2 = dom%Alpha(idx,1)
!                            a3 = dom%Alpha(idx,2)
!                            ui = champs1%Depla(idx,i_dir)
!                            call select_terms(a1,a2,a3,t,t1,t2,t3)
!                            dom%R4(ee,i_dir,i,j,k,bnum) = dom%R4(ee,i_dir,i,j,k,bnum)*(1-dt*b1)+dt*ui
!                            dom%R5(ee,i_dir,i,j,k,bnum) = dom%R5(ee,i_dir,i,j,k,bnum)*(1-dt*b2)+dt*ui
!                            dom%R6(ee,i_dir,i,j,k,bnum) = dom%R6(ee,i_dir,i,j,k,bnum)*(1-dt*b3)+dt*ui
!                        end do
!                    end do
!                end do
!            end do
!        end do
    end subroutine update_convolution_terms

    subroutine newmark_predictor_solidpml(dom, Tdomain)
        type(domain_solidpml), intent (INOUT) :: dom
        type (domain), intent (INOUT) :: Tdomain
        !
        integer :: n, indpml, indsol

        ! Reset forces
        dom%champs1%Forces = 0d0

        ! Coupling at solid PML interface
        do n = 0,Tdomain%intSolPml%surf0%nbtot-1
            indsol = Tdomain%intSolPml%surf0%map(n)
            indpml = Tdomain%intSolPml%surf1%map(n)
            dom%champs0%Veloc(indpml,:) = Tdomain%sdom%champs0%Veloc(indsol,:)
            dom%champs0%Depla(indpml,:) = Tdomain%sdom%champs0%Depla(indsol,:)
        enddo

        ! The prediction will be based on the current state
        dom%champs1%Depla = dom%champs0%Depla
        dom%champs1%Veloc = dom%champs0%Veloc
    end subroutine newmark_predictor_solidpml

    subroutine newmark_corrector_solidpml(dom, dt, t)
        type(domain_solidpml), intent (INOUT) :: dom
        double precision :: dt, t
        !
        integer :: i_dir, n
        !
        real(fpp) :: D, V, F ! Displacement, Velocity, Force

        ! Note: the solid domain use a leap-frop time scheme => we get D_n+3/2 and V_n+1
        !       to compute V_n+2, we need F_n+1. To compute F_n+1, we need D_n+1 and V_n+1
        !       D_n+1 is estimated with D_n+3/2 and V_n+1

        ! Update velocity: compute V_n+2
        do i_dir = 0,2
            do n = 0,dom%nglltot-1
                ! Get V = V_n+1
                V = dom%champs0%Veloc(n,i_dir) ! V = V_n+1

                ! Estimate D = D_n+1
                D = dom%champs0%Depla(n,i_dir) - V*0.5*dt

                ! Compute F_n+1 : (61a) from Ref1 with F = -F (as we add -Fo* in forces_int_sol_pml)
                F =    &
                    dom%MassMat(n)*(- dom%DumpMat(n)*V - dom%MasUMat(n)*dom%champs0%Depla(n,i_dir)  &
                    + dom%champs1%Forces(n,i_dir))

                ! Compute V_n+2
                ! since dom%MassMat = 1./dom%MassMat (define_arrays inverse_mass_mat)
                dom%champs0%Veloc(n,i_dir) = V + dt*F

                ! Update displacement: compute D_n+5/2
                dom%champs0%Depla(n,i_dir) = dom%champs0%Depla(n,i_dir) + dt * dom%champs0%Veloc(n,i_dir)
            end do
        end do


        ! Note: do NOT apply (dirichlet) BC for PML
        !       if PML absorption would be turned off <=> solid domain without dirichlet BC (neumann only)
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
