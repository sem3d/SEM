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
    implicit none
#include "index.h"

    !! CPML parameters: for the very first implementation, parameters are hard-coded. TODO : read parameters (kappa_* ?) from input.spec ?
    real(fpp), private, parameter :: c_x = 1., c_y = 1., c_z = 1.
    integer,   private, parameter :: n_x = 2,  n_y = 2,  n_z = 2
    real(fpp), private, parameter :: r_c = 0.001
    integer,   private, parameter :: kappa_0 = 1, kappa_1 = 0
    real(fpp), private, parameter :: L_x = -1., L_y = -1., L_z = -1.

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

            allocate(dom%Density_(      0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Cij_    (0:20, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))

            allocate (dom%Jacob_  (        0:ngll-1,0:ngll-1,0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            allocate (dom%InvGrad_(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))

            allocate(dom%Idom_(0:ngll-1,0:ngll-1,0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            dom%m_Idom = 0

            allocate(dom%R1_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nblocks-1,0:VCHUNK-1))
            allocate(dom%R2_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nblocks-1,0:VCHUNK-1))
            allocate(dom%R3_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nblocks-1,0:VCHUNK-1))
        end if

        ! Allocation et initialisation de champs0 pour les PML solides
        if (dom%nglltot /= 0) then
            allocate(dom%champs0%Depla (0:dom%nglltot-1,0:2))
            allocate(dom%champs0%Veloc (0:dom%nglltot-1,0:2))
            allocate(dom%champs1%Depla (0:dom%nglltot-1,0:2))
            allocate(dom%champs1%Veloc (0:dom%nglltot-1,0:2))
            dom%champs0%Depla  = 0d0
            dom%champs0%Veloc  = 0d0
            dom%champs1%Depla  = 0d0
            dom%champs1%Veloc  = 0d0

            allocate(dom%Forces(0:dom%nglltot-1,0:2))
            dom%Forces = 0d0

            ! Allocation de MassMat pour les PML solides
            allocate(dom%MassMat(0:dom%nglltot-1))
            dom%MassMat = 0d0
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid cpml domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_solidpml

    subroutine deallocate_dom_solidpml (dom)
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Cij    )) deallocate(dom%m_Cij    )

        if(allocated(dom%m_Jacob  )) deallocate(dom%m_Jacob  )
        if(allocated(dom%m_InvGrad)) deallocate(dom%m_InvGrad)

        if(allocated(dom%m_Idom)) deallocate(dom%m_Idom)

        if(allocated(dom%m_R1)) deallocate(dom%m_R1)
        if(allocated(dom%m_R2)) deallocate(dom%m_R2)
        if(allocated(dom%m_R3)) deallocate(dom%m_R3)

        if(allocated(dom%gllc))    deallocate(dom%gllc)
        if(allocated(dom%gllw))    deallocate(dom%gllw)
        if(allocated(dom%hprime))  deallocate(dom%hprime)
        if(allocated(dom%htprime)) deallocate(dom%htprime)

        if(allocated(dom%champs0%Depla )) deallocate(dom%champs0%Depla )
        if(allocated(dom%champs0%Veloc )) deallocate(dom%champs0%Veloc )
        if(allocated(dom%champs1%Depla )) deallocate(dom%champs1%Depla )
        if(allocated(dom%champs1%Veloc )) deallocate(dom%champs1%Veloc )

        if(allocated(dom%Forces )) deallocate(dom%Forces )
        if(allocated(dom%MassMat)) deallocate(dom%MassMat)
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
                        fieldV(i,j,k,:) = dom%champs0%Veloc(ind,0) + &
                                          dom%champs0%Veloc(ind,1) + &
                                          dom%champs0%Veloc(ind,2)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldA(i,j,k,:) = dom%Massmat(ind) * ( dom%Forces(ind,0) + &
                                                               dom%Forces(ind,1) + &
                                                               dom%Forces(ind,2) )
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

    subroutine init_material_properties_solidpml(dom, lnum, i, j, k, density, lambda, mu)
        type(domain_solidpml), intent(inout) :: dom
        integer, intent(in) :: lnum
        integer, intent(in) :: i, j, k ! -1 means :
        real(fpp), intent(in) :: density
        real(fpp), intent(in) :: lambda
        real(fpp), intent(in) :: mu
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        if (i==-1 .and. j==-1 .and. k==-1) then
            dom%Density_(:,:,:,bnum,ee) = density
        else
            dom%Density_(i,j,k,bnum,ee) = density
        end if
        ! TODO : compute dom%Cij
    end subroutine init_material_properties_solidpml

    subroutine init_local_mass_solidpml(dom,specel,i,j,k,ind,Whei)
        type(domain_solidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real Whei
        !
        integer :: bnum, ee
        real(fpp) :: ab2

        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ! Delta term from L : (12a) or (14a) from Ref1

        ab2 = 1. ! TODO : compute ab2 !...
        dom%MassMat(ind) = dom%MassMat(ind) + ab2*dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)*Whei

        ! TODO: compute dom%MassMat
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

    subroutine forces_int_sol_pml(dom, champs1, bnum)
        use m_calcul_forces_solidpml
        type(domain_solidpml), intent(inout) :: dom
        type(champssolidpml), intent(inout) :: champs1
        integer :: bnum
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
        call calcul_forces_solidpml(dom,bnum,Fox,Foy,Foz,Depla)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        e = bnum*VCHUNK+ee
                        if (e>=dom%nbelem) exit
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        dom%Forces(idx,0) = dom%Forces(idx,0)-Fox(ee,i,j,k)
                        dom%Forces(idx,1) = dom%Forces(idx,1)-Foy(ee,i,j,k)
                        dom%Forces(idx,2) = dom%Forces(idx,2)-Foz(ee,i,j,k)
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

    subroutine finalize_solidpml_properties(dom)
        type (domain_solidpml), intent (INOUT), target :: dom
        !
        ! Useless, kept for compatibility with SolidPML (build), can be deleted later on. TODO : kill this method.
    end subroutine finalize_solidpml_properties

    subroutine update_convolution_terms(dom)
        type(domain_solidpml), intent (INOUT) :: dom
        ! TODO : compute / update dom%m_R1, dom%m_R2, dom%m_R3
    end subroutine

    subroutine newmark_predictor_solidpml(dom, Tdomain)
        type(domain_solidpml), intent (INOUT) :: dom
        type (domain), intent (INOUT) :: Tdomain
        !
        dom%champs1%Depla = dom%champs0%Depla
        dom%champs1%Veloc = dom%champs0%Veloc
        dom%Forces = 0d0
        call update_convolution_terms(dom)
    end subroutine newmark_predictor_solidpml

    subroutine newmark_corrector_solidpml(dom, dt)
        type(domain_solidpml), intent (INOUT) :: dom
        double precision :: dt
        !
        integer :: i_dir, n, indpml
        ! Update velocity in champs0 (Note: dom%MassMat = 1./dom%MassMat in define_arrays::inverse_mass_mat)
        do i_dir = 0,2
            dom%champs0%Veloc(:,i_dir) = dom%champs0%Veloc(:,i_dir) + &
                                         dt * ( dom%Forces(:,i_dir) * dom%MassMat(:) ) ! dt * acceleration
        enddo
        ! Apply BC (dirichlet)
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs0%Veloc(indpml,:) = 0.
        enddo
        ! Update displacement in champs0
        dom%champs0%Depla = dom%champs0%Depla + dt * dom%champs0%Veloc
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
