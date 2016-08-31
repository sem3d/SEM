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

#include "dom_solidcpml_macro.F90"

module dom_solidpml
    use constants
    use sdomain
    use champs_solidpml
    use ssubdomains
    implicit none
#include "index.h"

contains

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
        call compute_gll_data(ngll, dom%gllc, dom%gllw, dom%hprime, dom%htprime)

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)

            nblocks = ((nbelem+VCHUNK-1)/VCHUNK)
            dom%nblocks = nblocks

            allocate (dom%Jacob_  (        0:ngll-1,0:ngll-1,0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            allocate (dom%InvGrad_(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))

            allocate(dom%Idom_(0:ngll-1,0:ngll-1,0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            dom%m_Idom = 0
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

            ! Allocation de MassMat pour les PML solides
            allocate(dom%MassMat(0:dom%nglltot-1))
            dom%MassMat = 0d0

            ! Allocation de DumpMat pour les PML solides
            allocate(dom%DumpMat(0:dom%nglltot-1))
            dom%DumpMat = 0d0

            ! Allocation de MasUMat pour les PML solides
            allocate(dom%MasUMat(0:dom%nglltot-1))
            dom%MasUMat = 0d0

            ! Allocation des Ri pour les PML solides
            allocate(dom%R1(0:dom%nglltot-1,0:2))
            dom%R1 = 0d0
            allocate(dom%R2(0:dom%nglltot-1,0:2))
            dom%R2 = 0d0
            allocate(dom%R3(0:dom%nglltot-1,0:2))
            dom%R3 = 0d0
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid cpml domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"

        ! CPML parameters initialisation: for the very first implementation, parameters are hard-coded.
        ! TODO : read parameters (kappa_* ?) from input.spec ?

        dom%c = 1.
        dom%n = 2
        dom%r_c = 0.001
        dom%kappa_0 = 1; dom%kappa_1 = 0;
        dom%alphamax = 0.
    end subroutine allocate_dom_solidpml

    subroutine deallocate_dom_solidpml (dom)
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom

        if(allocated(dom%m_Jacob  )) deallocate(dom%m_Jacob  )
        if(allocated(dom%m_InvGrad)) deallocate(dom%m_InvGrad)

        if(allocated(dom%m_Idom)) deallocate(dom%m_Idom)

        if(allocated(dom%gllc))    deallocate(dom%gllc)
        if(allocated(dom%gllw))    deallocate(dom%gllw)
        if(allocated(dom%hprime))  deallocate(dom%hprime)
        if(allocated(dom%htprime)) deallocate(dom%htprime)

        if(allocated(dom%champs0%Depla )) deallocate(dom%champs0%Depla )
        if(allocated(dom%champs0%Veloc )) deallocate(dom%champs0%Veloc )
        if(allocated(dom%champs0%Forces )) deallocate(dom%champs0%Forces )
        if(allocated(dom%champs1%Depla )) deallocate(dom%champs1%Depla )
        if(allocated(dom%champs1%Veloc )) deallocate(dom%champs1%Veloc )
        if(allocated(dom%champs1%Forces )) deallocate(dom%champs1%Forces )

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)
        if(allocated(dom%DumpMat)) deallocate(dom%DumpMat)
        if(allocated(dom%MasUMat)) deallocate(dom%MasUMat)

        if(allocated(dom%R1)) deallocate(dom%R1)
        if(allocated(dom%R2)) deallocate(dom%R2)
        if(allocated(dom%R3)) deallocate(dom%R3)
    end subroutine deallocate_dom_solidpml

    subroutine get_solidpml_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain_solidpml)                      :: dom
        integer, dimension(0:8)                    :: out_variables
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
        ! Handle on node global coords : mandatory to compute distances in the PML (solidcpml_abk)
        dom%GlobCoord => Tdomain%GlobCoord ! Pointer to coord (avoid allocate + copy, just point to it)
        ! Handle on materials
        dom%sSubDomain => Tdomain%sSubDomain
    end subroutine init_domain_solidpml

    subroutine init_material_properties_solidpml(dom, lnum, i, j, k, density, lambda, mu)
        type(domain_solidpml), intent(inout) :: dom
        integer, intent(in) :: lnum
        integer, intent(in) :: i, j, k ! -1 means :
        real(fpp), intent(in) :: density
        real(fpp), intent(in) :: lambda
        real(fpp), intent(in) :: mu
        !
        ! Useless, kept for compatibility with SolidPML (build), can be deleted later on. TODO : kill this method.
    end subroutine init_material_properties_solidpml

    ! TODO : renommer init_local_mass_solidpml... en init_global_mass_solidpml ? Vu qu'on y met a jour la masse globale !?
    !        attention, ceci impacte le build (compatibilité avec SolidPML)
    subroutine init_local_mass_solidpml(dom,specel,i,j,k,ind,Whei)
        type(domain_solidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real Whei
        !
        integer :: bnum, ee
        real(fpp) :: xi, xoverl, dxi, d0, alpha(0:2), beta(0:2), kappa(0:2) ! solidcpml_abk
        real(fpp) :: g0, g1, g2 ! solidcpml_gamma_ab
        real(fpp) :: g101, g212, g002 ! solidcpml_gamma_abc
        real(fpp) :: a0b, a1b, a2b
        real(fpp) :: density
        integer :: mi

        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ind = dom%Idom_(i,j,k,bnum,ee)
        if (.not. dom%sSubDomain(specel%mat_index)%material_type == "P") &
            stop "init_geometric_properties_solidpml : material is not a PML material"
        density = dom%sSubDomain(specel%mat_index)%Ddensity

        ! Compute alpha, beta, kappa
        mi = specel%mat_index
        solidcpml_abk(0,i,j,k,bnum,ee,mi)
        solidcpml_abk(1,i,j,k,bnum,ee,mi)
        solidcpml_abk(2,i,j,k,bnum,ee,mi)

        ! Delta 2d derivative term from L : (12a) or (14a) from Ref1
        a0b = kappa(0)*kappa(1)*kappa(2)
        dom%MassMat(ind) = dom%MassMat(ind) + density*a0b*dom%Jacob_(i,j,k,bnum,ee)*Whei
        if (abs(dom%MassMat(ind)) < solidcpml_eps) stop "ERROR : MassMat is null" ! Check

        ! Delta 1st derivative term from L : (12a) or (14a) from Ref1
        solidcpml_gamma_ab(g0,beta,0,alpha,0)
        solidcpml_gamma_ab(g1,beta,1,alpha,1)
        solidcpml_gamma_ab(g2,beta,2,alpha,2)
        a1b = a0b*(g0+g1+g2)
        dom%DumpMat(ind) = dom%DumpMat(ind) + density*a1b*dom%Jacob_(i,j,k,bnum,ee)*Whei

        ! Delta term from L : (12a) or (14a) from Ref1
        solidcpml_gamma_abc(g101,beta,1,alpha,0,alpha,1)
        solidcpml_gamma_abc(g212,beta,2,alpha,1,alpha,2)
        solidcpml_gamma_abc(g002,beta,0,alpha,0,alpha,2)
        a2b = a0b*(g0*g101+g1*g212+g2*g002)
        dom%MasUMat(ind) = dom%MasUMat(ind) + density*a2b*dom%Jacob_(i,j,k,bnum,ee)*Whei
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
        call calcul_forces_solidpml(dom,bnum,Fox,Foy,Foz,Depla,Tdomain)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        e = bnum*VCHUNK+ee
                        if (e>=dom%nbelem) exit
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        champs1%Forces(idx,0) = champs1%Forces(idx,0)-Fox(ee,i,j,k)
                        champs1%Forces(idx,1) = champs1%Forces(idx,1)-Foy(ee,i,j,k)
                        champs1%Forces(idx,2) = champs1%Forces(idx,2)-Foz(ee,i,j,k)
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
        real(fpp) :: fmax
        integer :: nsrc

        ! Compute alphamax (from fmax)

        ! TODO : compute max of all fmax of all procs: need to speak about that with Ludovic: where ? how ?
        ! My understanding is that all procs process different elements (PML or not), so all procs are NOT
        ! here (in init_solidpml_properties) at the same time : broadcast fmax of each proc to all procs and get
        ! the max could NOT work ?!... Right ? Wrong ? Don't know !...
        ! TODO : replace all this by fmax from read_input.c
        if (Tdomain%nb_procs /= 1) stop "ERROR : SolidCPML is limited to monoproc for now"

        fmax = -1.
        do nsrc = 0, Tdomain%n_source -1
            if (fmax < Tdomain%sSource(nsrc)%cutoff_freq) then
                fmax = Tdomain%sSource(nsrc)%cutoff_freq
            end if
        end do
        if (fmax < 0.) stop "SolidCPML : fmax < 0."
        dom%alphamax = M_PI * fmax
    end subroutine finalize_solidpml_properties

    subroutine update_convolution_terms(dom)
        type(domain_solidpml), intent (INOUT) :: dom
        ! TODO : compute / update dom%m_R1, dom%m_R2, dom%m_R3
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

        ! Update convolution terms
        call update_convolution_terms(dom)
    end subroutine newmark_predictor_solidpml

    subroutine newmark_corrector_solidpml(dom, dt)
        type(domain_solidpml), intent (INOUT) :: dom
        double precision :: dt
        !
        integer :: i_dir, n, indpml
        real(fpp) :: V(0:dom%nglltot-1), A(0:dom%nglltot-1) ! Velocity, Acceleration

        ! Update champs1 velocity from champs0
        do i_dir = 0,2
            ! Estimate V = V_n+3/2
            V(:) = dom%champs0%Veloc(:,i_dir) ! V = V_n+1
            V(:) = 0.5 * ( V(:) + (dom%champs0%Depla(:,i_dir)-dom%DeplaPrev(:,i_dir))/dt ) ! V is corrected with U_n+1/2 to estimate V_n+3/2

            ! Compute acceleration
            A(:) =   dom%R1(:,i_dir)     + dom%R2(:,i_dir)                           + dom%R3(:,i_dir)           &
                   - dom%DumpMat(:)*V(:) - dom%MasUMat(:)*dom%champs0%Depla(:,i_dir) - dom%champs1%Forces(:,i_dir) ! (61a) from Ref1

            ! Compute V_n+2
            dom%champs0%Veloc(:,i_dir) = dom%champs0%Veloc(:,i_dir) + &
                                         dt * dom%MassMat(:) * A(:) ! dom%MassMat = 1./dom%MassMat (define_arrays inverse_mass_mat)
        enddo

        ! Update current state
        dom%champs0%Depla = dom%champs0%Depla + dt * dom%champs0%Veloc

        ! Apply BC for PML (dirichlet)
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs0%Veloc(indpml,:) = 0.
            dom%champs0%Depla(indpml,:) = 0.
        enddo
    end subroutine newmark_corrector_solidpml
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
