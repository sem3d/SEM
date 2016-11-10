!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_solid
    use constants
    use champs_solid
    implicit none
#include "index.h"

contains

    subroutine allocate_dom_solid (Tdomain, dom)
        use sdomain
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_solid), intent (INOUT) :: dom
        !
        integer :: nbelem, nblocks, ngll, n_solid
        logical :: aniso,nl_flag
        !

        ngll    = dom%ngll
        nbelem  = dom%nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call init_dombase(dom)

        aniso       = Tdomain%aniso
        n_solid     = Tdomain%n_sls
        dom%n_sls   = n_solid
        dom%aniso   = Tdomain%aniso
        nl_flag     = Tdomain%nl_flag

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! Do not allocate if not needed (save allocation/RAM)
            ! Wexo can have glls without elements
            nblocks = dom%nblocks

            allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Mu_     (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Kappa_  (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))

            if (nl_flag) then
                write(*,*) "nonlinear activated!",allocated(dom%nl_param)
                ! nonlinear parameters
                allocate(dom%nl_param)
                write(*,*) allocated(dom%nl_param)
                allocate(dom%nl_param%LMC)
                allocate(dom%nl_param%LMC%syld_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                allocate(dom%nl_param%LMC%ckin_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                allocate(dom%nl_param%LMC%kkin_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                allocate(dom%nl_param%LMC%rinf_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                allocate(dom%nl_param%LMC%biso_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                ! internal variables
                allocate(dom%radius_      (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                allocate(dom%stress_  (0:5,0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                allocate(dom%center_  (0:5,0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                allocate(dom%strain_  (0:5,0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                allocate(dom%plstrain_(0:5,0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                dom%strain_  (:,:,:,:,:,:) = zero
                dom%plstrain_(:,:,:,:,:,:) = zero
                dom%stress_  (:,:,:,:,:,:) = zero
                dom%center_  (:,:,:,:,:,:) = zero
                dom%radius_  (:,:,:,:,:)   = zero
            end if

            if (aniso) then
                allocate (dom%Cij_ (0:20, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            endif
            if (n_solid>0) then
                if (aniso) then
                    allocate (dom%Q_(0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                else
                    allocate (dom%Qs_         (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                    allocate (dom%Qp_         (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                    allocate (dom%onemPbeta_  (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                    allocate (dom%epsilonvol_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                    dom%epsilonvol_(:,:,:,:,:) = 0
                    allocate (dom%R_vol_       (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                    dom%R_vol_(:,:,:,:,:,:) = 0
                endif
                allocate (dom%onemSbeta_     (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%epsilondev_xx_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%epsilondev_yy_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%epsilondev_xy_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%epsilondev_xz_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%epsilondev_yz_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                dom%epsilondev_xx_(:,:,:,:,:) = 0
                dom%epsilondev_yy_(:,:,:,:,:) = 0
                dom%epsilondev_xy_(:,:,:,:,:) = 0
                dom%epsilondev_xz_(:,:,:,:,:) = 0
                dom%epsilondev_yz_(:,:,:,:,:) = 0
                allocate (dom%R_xx_           (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%R_yy_           (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%R_xy_           (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%R_xz_           (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%R_yz_           (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                !
                allocate (dom%omega_tau_s     (0:n_solid-1))
                allocate (dom%agamma_mu_      (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%agamma_kappa_   (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                !
                dom%R_xx_(:,:,:,:,:,:) = 0
                dom%R_yy_(:,:,:,:,:,:) = 0
                dom%R_xy_(:,:,:,:,:,:) = 0
                dom%R_xz_(:,:,:,:,:,:) = 0
                dom%R_yz_(:,:,:,:,:,:) = 0
            endif ! n_solid
        end if
        ! Allocation et initialisation de champs0 et champs1 pour les solides
        if (dom%nglltot /= 0) then
            ! Allocate ONE more gll than necessary to use as a dummy target for
            ! indirections for fake elements.
            allocate(dom%champs0%Forces(0:dom%nglltot,0:2))
            allocate(dom%champs0%Depla (0:dom%nglltot,0:2))
            allocate(dom%champs0%Veloc (0:dom%nglltot,0:2))
            allocate(dom%champs1%Forces(0:dom%nglltot,0:2))
            allocate(dom%champs1%Depla (0:dom%nglltot,0:2))
            allocate(dom%champs1%Veloc (0:dom%nglltot,0:2))

            dom%champs0%Forces = 0d0
            dom%champs0%Depla = 0d0
            dom%champs0%Veloc = 0d0
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_solid

    subroutine deallocate_dom_solid (dom,nl_flag)
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        logical, intent(in) :: nl_flag

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_Mu     )) deallocate(dom%m_Mu     )
        if(allocated(dom%m_Kappa  )) deallocate(dom%m_Kappa  )

        if(allocated(dom%m_Cij            )) deallocate (dom%m_Cij            )
        if(allocated(dom%m_Q              )) deallocate (dom%m_Q              )
        if(allocated(dom%m_Qs             )) deallocate (dom%m_Qs             )
        if(allocated(dom%m_Qp             )) deallocate (dom%m_Qp             )
        if(allocated(dom%m_onemPbeta      )) deallocate (dom%m_onemPbeta      )
        if(allocated(dom%m_epsilonvol     )) deallocate (dom%m_epsilonvol     )
        if(allocated(dom%m_R_vol          )) deallocate (dom%m_R_vol          )
        if(allocated(dom%m_onemSbeta      )) deallocate (dom%m_onemSbeta      )
        if(allocated(dom%m_epsilondev_xx  )) deallocate (dom%m_epsilondev_xx  )
        if(allocated(dom%m_epsilondev_yy  )) deallocate (dom%m_epsilondev_yy  )
        if(allocated(dom%m_epsilondev_xy  )) deallocate (dom%m_epsilondev_xy  )
        if(allocated(dom%m_epsilondev_xz  )) deallocate (dom%m_epsilondev_xz  )
        if(allocated(dom%m_epsilondev_yz  )) deallocate (dom%m_epsilondev_yz  )
        if(allocated(dom%m_R_xx           )) deallocate (dom%m_R_xx           )
        if(allocated(dom%m_R_yy           )) deallocate (dom%m_R_yy           )
        if(allocated(dom%m_R_xy           )) deallocate (dom%m_R_xy           )
        if(allocated(dom%m_R_xz           )) deallocate (dom%m_R_xz           )
        if(allocated(dom%m_R_yz           )) deallocate (dom%m_R_yz           )

        if(allocated(dom%champs0%Forces)) deallocate(dom%champs0%Forces)
        if(allocated(dom%champs0%Depla )) deallocate(dom%champs0%Depla )
        if(allocated(dom%champs0%Veloc )) deallocate(dom%champs0%Veloc )
        if(allocated(dom%champs1%Forces)) deallocate(dom%champs1%Forces)
        if(allocated(dom%champs1%Depla )) deallocate(dom%champs1%Depla )
        if(allocated(dom%champs1%Veloc )) deallocate(dom%champs1%Veloc )

        ! nonlinear parameters TODO
        if (nl_flag) then
            if(allocated(dom%nl_param%LMC%m_syld))  deallocate(dom%nl_param%LMC%m_syld)
            if(allocated(dom%nl_param%LMC%m_biso))  deallocate(dom%nl_param%LMC%m_biso)
            if(allocated(dom%nl_param%LMC%m_rinf))  deallocate(dom%nl_param%LMC%m_rinf)
            if(allocated(dom%nl_param%LMC%m_Ckin))  deallocate(dom%nl_param%LMC%m_Ckin)
            if(allocated(dom%nl_param%LMC%m_kkin))  deallocate(dom%nl_param%LMC%m_kkin)
            if(allocated(dom%nl_param%LMC       ))  deallocate(dom%nl_param%LMC)
            if(allocated(dom%nl_param           ))  deallocate(dom%nl_param)
            ! nonlinear internal variables
            if(allocated(dom%m_radius)) deallocate(dom%m_radius)
            if(allocated(dom%m_stress)) deallocate(dom%m_stress)
            if(allocated(dom%m_center)) deallocate(dom%m_center)
            if(allocated(dom%m_strain)) deallocate(dom%m_strain)
            if(allocated(dom%m_plstrain)) deallocate(dom%m_plstrain)
        endif

        call deallocate_dombase(dom)

        end subroutine deallocate_dom_solid

    subroutine get_solid_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy,&
        eps_vol, eps_dev, sig_dev, nl_flag, eps_dev_pl)
        use deriv3d
        implicit none
        !
        type(domain_solid), intent(inout)          :: dom
        integer, intent(in), dimension(0:)         :: out_variables
        integer, intent(in)                        :: lnum
        logical, intent(in)                        :: nl_flag
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev_pl
        !
        logical                  :: flag_gradU
        integer                  :: ngll, i, j, k, ind, nblocks
        real(fpp)                :: EXY, EXZ, EYZ
        real(fpp)                :: DXX, DXY, DXZ
        real(fpp)                :: DYX, DYY, DYZ
        real(fpp)                :: DZX, DZY, DZZ
        real(fpp)                :: divU
        real(fpp)                :: xmu, xlambda, xkappa
        real(fpp)                :: x2mu, xlambda2mu
        real(fpp)                :: onemSbeta, onemPbeta
        real(fpp), dimension(0:20) :: CC
        real, dimension(0:2,0:2) :: invgrad_ijk
        !
        integer :: bnum, ee, flag_gradUint

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        flag_gradUint = 0
        flag_gradUint = out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS)

        if (.not.nl_flag) then
            flag_gradUint = flag_gradUint     + &
                out_variables(OUT_PRESSION)   + &
                out_variables(OUT_EPS_VOL)    + &
                out_variables(OUT_EPS_DEV)    + &
                out_variables(OUT_STRESS_DEV)
        endif

        flag_gradU = flag_gradUint /= 0
        ngll = dom%ngll

        if (out_variables(OUT_VITESSE) == 1) then
            if(allocated(fieldV)) then
                if(any(shape(fieldV(:,:,:,0)) /= ngll)) deallocate(fieldV)
            end if
            if(.not. allocated(fieldV)) allocate(fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:OUT_VAR_DIMS_3D(OUT_VITESSE)-1))
        else
            if(.not. allocated(fieldV)) allocate(fieldV(0:0,0:0,0:0,0:0))
        end if
        !
        if (out_variables(OUT_ACCEL) == 1) then
            if(allocated(fieldA)) then
                if(any(shape(fieldA(:,:,:,0)) /= ngll)) deallocate(fieldA)
            end if
            if(.not. allocated(fieldA)) allocate(fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:OUT_VAR_DIMS_3D(OUT_ACCEL)-1))
        else
            if(.not. allocated(fieldA)) allocate(fieldA(0:0,0:0,0:0,0:0))
        end if
        !
        if (out_variables(OUT_DEPLA) == 1) then
            if(allocated(fieldU)) then
                if(any(shape(fieldU(:,:,:,0)) /= ngll)) deallocate(fieldU)
            end if
            if(.not. allocated(fieldU)) allocate(fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:OUT_VAR_DIMS_3D(OUT_DEPLA)-1))
        else
            if(.not. allocated(fieldU)) allocate(fieldU(0:0,0:0,0:0,0:0))
        end if
        !
        if (out_variables(OUT_PRESSION) == 1) then
            if(allocated(fieldP)) then
                if(any(shape(fieldP) /= ngll)) deallocate(fieldP)
            end if
            if(.not. allocated(fieldP)) allocate(fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
        else
            if(.not. allocated(fieldP)) allocate(fieldP(0:0,0:0,0:0))
        end if
        !
        if (out_variables(OUT_EPS_VOL) == 1) then
            if(allocated(eps_vol)) then
                if(any(shape(eps_vol) /= ngll)) deallocate(eps_vol)
            end if
            if(.not. allocated(eps_vol)) allocate(eps_vol(0:ngll-1,0:ngll-1,0:ngll-1))
        else
            if(.not. allocated(eps_vol)) allocate(eps_vol(0:0,0:0,0:0))
        end if
        !
        if (out_variables(OUT_ENERGYP) == 1) then
            if(allocated(P_energy)) then
                if(any(shape(P_energy) /= ngll)) deallocate(P_energy)
            end if
            if(.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        else
            if(.not. allocated(P_energy)) allocate(P_energy(0:0,0:0,0:0))
        end if
        !
        if (out_variables(OUT_ENERGYS) == 1) then
            if(allocated(S_energy)) then
                if(any(shape(S_energy) /= ngll)) deallocate(S_energy)
            end if
            if(.not. allocated(S_energy)) allocate(S_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        else
            if(.not. allocated(S_energy)) allocate(S_energy(0:0,0:0,0:0))
        end if
        !
        if (out_variables(OUT_EPS_DEV) == 1) then
            if(allocated(eps_dev)) then
                if(any(shape(eps_dev(:,:,:,0)) /= ngll)) deallocate(eps_dev)
            end if
            if(.not. allocated(eps_dev)) allocate(eps_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:OUT_VAR_DIMS_3D(OUT_EPS_DEV)-1))
        else
            if(.not. allocated(eps_dev)) allocate(eps_dev(0:0,0:0,0:0,0:0))
        end if
        !
        if (out_variables(OUT_EPS_DEV_PL) == 1) then
            if(allocated(eps_dev_pl)) then
                if(any(shape(eps_dev_pl(:,:,:,0)) /= ngll)) deallocate(eps_dev_pl)
            end if
            if(.not. allocated(eps_dev_pl)) allocate(eps_dev_pl(0:ngll-1,0:ngll-1,0:ngll-1,0:OUT_VAR_DIMS_3D(OUT_EPS_DEV_PL)-1))
        else
            if(.not. allocated(eps_dev_pl)) allocate(eps_dev_pl(0:0,0:0,0:0,0:0))
        end if
        !
        if (out_variables(OUT_STRESS_DEV) == 1) then
            if(allocated(sig_dev)) then
                if(any(shape(sig_dev(:,:,:,0)) /= ngll)) deallocate(sig_dev)
            end if
            if(.not. allocated(sig_dev)) allocate(sig_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:OUT_VAR_DIMS_3D(OUT_STRESS_DEV)-1))
        else
            if(.not. allocated(sig_dev)) allocate(sig_dev(0:0,0:0,0:0,0:0))
        end if

        !print*, "BEFORE GET VELOCITY"
        ! GET VELOCITY (no gradient required)
        if (out_variables(OUT_VITESSE) == 1) then
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        fieldV(i,j,k,:) = dom%champs0%Veloc(ind,:)
                    enddo
                enddo
            enddo
        endif

        !print*, "BEFORE GET ACC"
        ! GET ACCLELERATION (no gradient required)
        if (out_variables(OUT_ACCEL) == 1) then
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        fieldA(i,j,k,:) = dom%champs0%Forces(ind,:)
                    enddo
                enddo
            enddo
        endif

        !print*, "BEFORE GET DISP"
        ! GET DISPLACEMENT (direct output or for gradient calculation)
        if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        fieldU(i,j,k,:) = dom%champs0%Depla(ind,:)
                    enddo
                enddo
            enddo
        end if


        !print*, "BEFORE GET OTHER"
        ! GET OTHER VARIABLES (gradient required)
        if (flag_gradU) then
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1

                        !print*, "IND "
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        ! GRADIENT

                        !print*, "BEFORE GRADIENT"
                        invgrad_ijk = dom%InvGrad_(:,:,i,j,k,bnum,ee) ! cache for performance
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,0),DXX,DYX,DZX)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,1),DXY,DYY,DZY)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,2),DXZ,DYZ,DZZ)
                        divU = DXX+DYY+DZZ
                        ! PRESSION
                        !print*, "BEFORE PRESSION"
                        if (out_variables(OUT_PRESSION) == 1) then
                            fieldP(i,j,k) = -(dom%Lambda_(i,j,k,bnum,ee)+two*M_1_3*dom%Mu_(i,j,k,bnum,ee))*divU
                        endif
                        ! EPSVOL
                        !print*, "BEFORE EPSVOL"
                        if (out_variables(OUT_EPS_VOL) == 1) then
                            eps_vol(i,j,k) = divU
                        end if
                        ! ELASTIC-VISCOELASTIC MODULA
                        !print*, "BEFORE ELASTIC-VISCO MODULA"
                        if (out_variables(OUT_ENERGYP) == 1 .or. &
                            out_variables(OUT_ENERGYS) == 1 .or. &
                            out_variables(OUT_STRESS_DEV) == 1) then
                            if (dom%aniso) then
                                CC = dom%Cij_(:,i,j,k,bnum,ee)
                            else
                                xmu     = dom%Mu_    (i,j,k,bnum,ee)
                                xlambda = dom%Lambda_(i,j,k,bnum,ee)
                                xkappa  = dom%Kappa_ (i,j,k,bnum,ee)
                                if (dom%n_sls>0) then
                                    onemSbeta = dom%onemSbeta_(i,j,k,bnum,ee)
                                    onemPbeta = dom%onemPbeta_(i,j,k,bnum,ee)
                                    xmu    = xmu * onemSbeta
                                    xkappa = xkappa * onemPbeta
                                endif
                                x2mu       = two * xmu
                                xlambda2mu = xlambda + x2mu
                            end if
                        endif
                        ![TODOLUCIANO]
                        ! P-ENERGY
                        !print*, "BEFORE P_EN"
                        if (out_variables(OUT_ENERGYP) == 1) then
                            P_energy(i,j,k) = zero
                            !half * xlambda2mu * eps_vol(i,j,k)**2
                        end if
                        ! S-ENERGY
                        !print*, "BEFORE S_EN"
                        if (out_variables(OUT_ENERGYS) == 1) then
                            S_energy(i,j,k) = zero
!                            if (.not. dom%aniso) then
!                                S_energy(i,j,k) =  zero
!                                ! half*xmu * (       DXY**2 + DYX**2 &
!                                                              +     DXZ**2 + DZX**2 &
!                                                              +     DYZ**2 + DZY**2 &
!                                                              - 2 * DXY * DYX     &
!                                                              - 2 * DXZ * DZX     &
!                                                              - 2 * DYZ * DZY )
                        end if
                        ! DEVIATORIC STRAIN

                        !print*, "BEFORE DEV STRAIN"
                        if (out_variables(OUT_EPS_DEV) == 1) then
                            eps_dev(i,j,k,0:5) = zero
                            eps_dev(i,j,k,0) = DXX - M_1_3 * divU
                            eps_dev(i,j,k,1) = DYY - M_1_3 * divU
                            eps_dev(i,j,k,2) = DZZ - M_1_3 * divU
                            eps_dev(i,j,k,3) = (DXY + DYX)
                            eps_dev(i,j,k,4) = (DZX + DXZ)
                            eps_dev(i,j,k,5) = (DZY + DYZ)
                        endif
                        ! DEVIATORIC STRESS
                        !print*, "BEFORE DEV STRESS"
                        if (out_variables(OUT_STRESS_DEV) == 1) then
                            sig_dev(i,j,k,0:5) = zero
                            if (dom%aniso) then
                                EXY = (DXY+DYX)*M_SQRT1_2
                                EYZ = (DYZ+DZY)*M_SQRT1_2
                                EXZ = (DXZ+DZX)*M_SQRT1_2
                                sig_dev(i,j,k,0) = DXX*CC( 0) + DYY*CC( 1) + DZZ*CC( 2) + EYZ*CC( 3) + EXZ*CC( 4) + EXY*CC( 5)
                                sig_dev(i,j,k,1) = DXX*CC( 1) + DYY*CC( 6) + DZZ*CC( 7) + EYZ*CC( 8) + EXZ*CC( 9) + EXY*CC(10)
                                sig_dev(i,j,k,2) = DXX*CC( 2) + DYY*CC( 7) + DZZ*CC(11) + EYZ*CC(12) + EXZ*CC(13) + EXY*CC(14)
                                sig_dev(i,j,k,3) = DXX*CC( 5) + DYY*CC(10) + DZZ*CC(14) + EYZ*CC(17) + EXZ*CC(19) + EXY*CC(20)
                                sig_dev(i,j,k,4) = DXX*CC( 4) + DYY*CC( 9) + DZZ*CC(13) + EYZ*CC(16) + EXZ*CC(18) + EXY*CC(19)
                                sig_dev(i,j,k,5) = DXX*CC( 3) + DYY*CC( 8) + DZZ*CC(12) + EYZ*CC(15) + EXZ*CC(16) + EXY*CC(17)
                            else
                                sig_dev(i,j,k,0) = x2mu * (DXX - M_1_3 * divU)
                                sig_dev(i,j,k,1) = x2mu * (DYY - M_1_3 * divU)
                                sig_dev(i,j,k,2) = x2mu * (DZZ - M_1_3 * divU)
                                sig_dev(i,j,k,3) =  xmu * (DXY + DYX)
                                sig_dev(i,j,k,4) =  xmu * (DXZ + DZX)
                                sig_dev(i,j,k,5) =  xmu * (DYZ + DZY)
                            endif
                        endif
                    enddo
                enddo
            enddo

            !print*, "AFTER BOUCLE LINEAR "
        else

            !print*, "NON LINEAR -----------------------------------"
            do k=0,ngll-1
                do j=0,ngll-1
                    do i=0,ngll-1
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        ! PRESSION
                        if (out_variables(OUT_PRESSION) == 1) then
                            fieldP(i,j,k) = -sum(dom%stress_(0:2,i,j,k,bnum,ee))*M_1_3
                        endif
                        ! EPSVOL
                        if (out_variables(OUT_EPS_VOL) == 1) then
                            eps_vol(i,j,k) = -sum(dom%strain_(0:2,i,j,k,bnum,ee))
                        end if
                        ![TODOLUCIANO]
                        ! P-ENERGY
                        if (out_variables(OUT_ENERGYP) == 1) then
                            P_energy(i,j,k) = zero
                        end if
                        ! S-ENERGY
                        if (out_variables(OUT_ENERGYS) == 1) then
                            S_energy(i,j,k) = zero
!                            if (.not. dom%aniso) then
!                                S_energy(i,j,k) =  zero
!                                ! half*xmu * (       DXY**2 + DYX**2 &
!                                                              +     DXZ**2 + DZX**2 &
!                                                              +     DYZ**2 + DZY**2 &
!                                                              - 2 * DXY * DYX     &
!                                                              - 2 * DXZ * DZX     &
!                                                              - 2 * DYZ * DZY )
                        end if
                        ! DEVIATORIC TOTAL STRAIN
                        if (out_variables(OUT_EPS_DEV) == 1) then
                            eps_dev(i,j,k,:)   = zero
                            eps_dev(i,j,k,:)   = dom%strain_(:,i,j,k,bnum,ee)
                            eps_dev(i,j,k,0:2) = eps_dev(i,j,k,0:2)-sum(eps_dev(i,j,k,0:2))*M_1_3
                        endif
                        ! DEVIATORIC PLASTIC STRAIN
                        if (out_variables(OUT_EPS_DEV_PL) == 1) then
                            eps_dev_pl(i,j,k,:)   = zero
                            eps_dev_pl(i,j,k,:)   = dom%plstrain_(:,i,j,k,bnum,ee)
                            eps_dev_pl(i,j,k,0:2) = eps_dev_pl(i,j,k,0:2)-sum(eps_dev_pl(i,j,k,0:2))*M_1_3
                        endif
                        ! DEVIATORIC STRESS
                        if (out_variables(OUT_STRESS_DEV) == 1) then
                            sig_dev(i,j,k,:)   = zero
                            sig_dev(i,j,k,:)   = dom%stress_(:,i,j,k,bnum,ee)
                            sig_dev(i,j,k,0:2) = sig_dev(i,j,k,0:2) - sum(sig_dev(i,j,k,0:2))*M_1_3
                        endif
                    enddo
                enddo
            enddo
        end if
        !
    end subroutine get_solid_dom_var


    subroutine get_solid_dom_elem_energy(dom, lnum, P_energy, S_energy, R_energy, C_energy)
        use deriv3d
        implicit none
        !
        type(domain_solid), intent(inout)          :: dom
        integer, intent(in)                        :: lnum
        real(fpp), dimension(:,:,:), allocatable, intent(inout) :: P_energy, S_energy, R_energy !R_energy = Residual energy (tend to zero as propagation takes place)
        real(fpp), dimension(:,:,:), allocatable, intent(inout) :: C_energy !Cinetic energy
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV

        integer                  :: ngll, i, j, k, ind
        real(fpp)                :: xmu, xlambda, xkappa, xdensity
        real(fpp)                :: onemSbeta, onemPbeta
        real(fpp)                :: xeps_vol
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        real(fpp), dimension(0:2) ::xvel
        real :: dUx_dx,dUx_dy,dUx_dz,  &
                dUy_dx,dUy_dy,dUy_dz,  &
                dUz_dx,dUz_dy,dUz_dz
        !
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)


        ngll = dom%ngll

        !Dellocation
        if(allocated(S_energy)) then
            if(size(S_energy) /= ngll*ngll*ngll) deallocate(S_energy)
        end if

        if(allocated(P_energy)) then
            if(size(P_energy) /= ngll*ngll*ngll) deallocate(P_energy)
        end if

        if(allocated(R_energy)) then
            if(size(R_energy) /= ngll*ngll*ngll) deallocate(R_energy)
        end if

        if(allocated(C_energy)) then
            if(size(C_energy) /= ngll*ngll*ngll) deallocate(C_energy)
        end if

        if(allocated(fieldU)) then
            if(size(fieldU) /= ngll*ngll*ngll) deallocate(fieldU)
        end if

        if(allocated(fieldV)) then
            if(size(fieldV) /= ngll*ngll*ngll) deallocate(fieldV)
        end if

        !Allocation
        if(.not. allocated(S_energy)) allocate(S_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if(.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if(.not. allocated(R_energy)) allocate(R_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if(.not. allocated(C_energy)) allocate(C_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        S_energy = -1
        P_energy = -1
        if (dom%aniso) return


        if(.not. allocated(fieldU)) allocate(fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        if(.not. allocated(fieldV)) allocate(fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:2))

        ! First, get displacement.
        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    fieldU(i,j,k,:) = dom%champs0%Depla(ind,:)
                    fieldV(i,j,k,:) = dom%champs0%Veloc(ind,:)
                enddo
            enddo
        enddo

        ! Then, get the energies.
        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ! Compute gradU with displacement.
                    invgrad_ijk = dom%InvGrad_(:,:,i,j,k,bnum,ee) ! cache for performance

                    call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                         invgrad_ijk,fieldU(:,:,:,0),dUx_dx,dUx_dy,dUx_dz)
                    call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                         invgrad_ijk,fieldU(:,:,:,1),dUy_dx,dUy_dy,dUy_dz)
                    call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                         invgrad_ijk,fieldU(:,:,:,2),dUz_dx,dUz_dy,dUz_dz)


                    ! Get other variables.
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    xeps_vol = dUx_dx + dUy_dy + dUz_dz

                    xmu     = dom%Mu_    (i,j,k,bnum,ee)
                    xlambda = dom%Lambda_(i,j,k,bnum,ee)
                    xkappa  = dom%Kappa_ (i,j,k,bnum,ee)
                    xdensity = dom%Density_ (i,j,k,bnum,ee)
                    xvel     = fieldV(i,j,k,:)

                    if (dom%n_sls>0) then
                        onemSbeta = dom%onemSbeta_(i,j,k,bnum,ee)
                        onemPbeta = dom%onemPbeta_(i,j,k,bnum,ee)
                        xmu    = xmu * onemSbeta
                        xkappa = xkappa * onemPbeta
                    endif

                    P_energy(i,j,k) = ((0.5d0*xlambda) + xmu) * xeps_vol**2d0
                    S_energy(i,j,k) = xmu/2.0d0 * (                       &
                                                    (dUz_dy - dUy_dz)**2d0  &
                                                  + (dUx_dz - dUz_dx)**2d0  &
                                                  + (dUy_dx - dUx_dy)**2d0  &
                                                  )
                    R_energy(i,j,k) = 2.0d0*xmu*(dUx_dy*dUy_dx + dUx_dz*dUz_dx + dUy_dz*dUz_dy) &
                                     -2.0d0*xmu*(dUx_dx*dUy_dy + dUx_dx*dUz_dz + dUy_dy*dUz_dz)

                    C_energy(i,j,k) = 0.5d0*xdensity*(xvel(0)**2.0d0 + xvel(1)**2.0d0 + xvel(2)**2.0d0)

                    !PAPER: The Energy Partitioning and the Diffusive Character of the Seismic Coda, Shapiro et al, 2000

                enddo
            enddo
        enddo

        if(allocated(fieldU)) deallocate(fieldU)
        if(allocated(fieldV)) deallocate(fieldV)

    end subroutine get_solid_dom_elem_energy


    subroutine init_material_properties_solid(dom, lnum, mat, density, lambda, mu)
        use ssubdomains
        type(domain_solid), intent(inout) :: dom
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
        if (mat%deftype==MATDEF_MU_SYLD_RHO) then
            dom%nl_param%LMC%syld_(:,:,:,bnum,ee) = mat%DSyld
            dom%nl_param%LMC%ckin_(:,:,:,bnum,ee) = mat%DCkin
            dom%nl_param%LMC%kkin_(:,:,:,bnum,ee) = mat%DKkin
            dom%nl_param%LMC%rinf_(:,:,:,bnum,ee) = mat%DRinf
            dom%nl_param%LMC%biso_(:,:,:,bnum,ee) = mat%DBiso
        endif

        if (dom%n_sls>0)  then
            dom%Kappa_  (:,:,:,bnum,ee) = lambda + 2d0*mu/3d0
            if (dom%aniso) then
                dom%Q_(:,:,:,bnum,ee) = mat%Qmu
            else
                dom%Qs_(:,:,:,bnum,ee) = mat%Qmu
                dom%Qp_(:,:,:,bnum,ee) = mat%Qpression
            endif
        endif

    end subroutine init_material_properties_solid

    subroutine init_material_tensor_solid(dom, lnum, mat, density, Cij)
        use ssubdomains
        type(domain_solid), intent(inout) :: dom
        integer, intent(in) :: lnum
        type (subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: density
        real(fpp), dimension(1:6,1:6,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(in) :: Cij

        integer :: idef, ii, jj
        integer :: i, j, k
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        do i=0,dom%ngll-1
            do j=0,dom%ngll-1
                do k=0,dom%ngll-1
                    idef = 0
                    do ii = 1,6
                        do jj = ii,6
                            dom%Cij_(idef,i,j,k,bnum,ee) = Cij(ii,jj,i,j,k)
                            idef = idef + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo
        dom%Density_(:,:,:,bnum,ee) = density
    end subroutine init_material_tensor_solid

    subroutine init_local_mass_solid(dom,specel,i,j,k,ind,Whei)
        use selement
        type(domain_solid), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ! Solid.

        specel%MassMat(i,j,k) = Whei*dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_solid

    subroutine forces_int_solid(dom, champs1, bnum, nl_flag)
        use m_calcul_forces
        use m_calcul_forces_atn
        use m_calcul_forces_nl

        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(inout) :: champs1
        integer, intent(in) :: bnum
        !
        logical,    intent(in) :: nl_flag
        !
        integer :: ngll,i,j,k,i_dir,e,ee,idx
        integer :: n_solid
        logical :: aniso
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fox,Foy,Foz
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: Depla

        n_solid = dom%n_sls
        aniso   = dom%aniso
        ngll    = dom%ngll

        if (nl_flag) then
            ! PLASTIC CASE: VELOCITY PREDICTION
            do i_dir = 0,2
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            do ee = 0, VCHUNK-1
                                idx = dom%Idom_(i,j,k,bnum,ee)
                                Depla(ee,i,j,k,i_dir) = champs1%Veloc(idx,i_dir)
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        else
            ! ELASTIC CASE: DISPLACEMENT PREDICTION
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
        endif
        Fox = 0d0
        Foy = 0d0
        Foz = 0d0

        if (aniso) then
            if (n_solid>0) then
                call calcul_forces_aniso_atn(dom,bnum,Fox,Foy,Foz,Depla)
            else
                call calcul_forces_aniso(dom,bnum,Fox,Foy,Foz,Depla)
            end if
        else
            if (n_solid>0) then
                call calcul_forces_iso_atn(dom,bnum,Fox,Foy,Foz,Depla)
            else
                if (nl_flag) then
                    call calcul_forces_nl(dom,bnum,Fox,Foy,Foz,Depla)
                else
                    call calcul_forces_iso(dom,bnum,Fox,Foy,Foz,Depla)
                endif
            end if
        end if

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        e = bnum*VCHUNK+ee
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        champs1%Forces(idx,0) = champs1%Forces(idx,0)-Fox(ee,i,j,k)
                        champs1%Forces(idx,1) = champs1%Forces(idx,1)-Foy(ee,i,j,k)
                        champs1%Forces(idx,2) = champs1%Forces(idx,2)-Foz(ee,i,j,k)
                    enddo
                enddo
            enddo
        enddo
    end subroutine forces_int_solid

    subroutine compute_planeW_Exafield(lnum,ctime,Tdomain)
        use sdomain
        use Surface_prbl_type

        implicit none
        type(domain),               intent(inout) :: Tdomain
        real(kind=8),               intent(in  )  :: ctime
        integer,                    intent(in   ) :: lnum
        real(kind=8), dimension(0:2):: coord, displ, veloc, accel
        real(kind=8)                :: PWspeed
        character(len=20)           :: char
        integer                     :: ns, im, i, j, k, dom, ipw, ngll
        integer                     :: bnum, ee, ind, ss

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        do ns=1,size(Tdomain%list_PWBC)
           ipw     = Tdomain%list_PWBC(ns)
           im      = Tdomain%nsurfsource(ipw)%mat_index
           dom     = Tdomain%sSubDomain(im)%dom
           ngll     = Tdomain%sSubDomain(im)%NGLL
           write(char,*) Tdomain%nsurfsource(ipw)%index(1)

           do ss=0,size(Tdomain%sSurfaces)-1
              if (Tdomain%sSurfaces(ss)%name=="surface"//adjustl(char(:len_trim(char)))) &
                  exit
           end do
           PWspeed = Tdomain%sSurfaces(ss)%Elastic%PWspeed
           select case (dom)
                  case (DM_SOLID)
                       ngll     = Tdomain%sdom%ngll
                       do k=0,ngll-1
                          do j=0,ngll-1
                             do i=0,ngll-1
                                ind = Tdomain%sdom%Idom_(i,j,k,bnum,ee)
                                coord = Tdomain%GlobCoord(:,ind)-Tdomain%nsurfsource(ipw)%scoord(:)
                                call PlaneWavedispl(Tdomain%nsurfsource(ipw),coord,ctime,PWspeed, displ,veloc,accel)
                               if (allocated(Tdomain%sdom%champs0%Depla))  Tdomain%sdom%champs0%Depla(ind,:)  = Tdomain%sdom%champs0%Depla(ind,:) + displ
                               if (allocated(Tdomain%sdom%champs0%Veloc))  Tdomain%sdom%champs0%Veloc(ind,:)  = Tdomain%sdom%champs0%Veloc(ind,:) + veloc
                               if (allocated(Tdomain%sdom%champs0%Forces)) Tdomain%sdom%champs0%Forces(ind,:) = Tdomain%sdom%champs0%Forces(ind,:) + accel
                             enddo
                          enddo
                       enddo
                   case(DM_FLUID)
                      ! pas encore implémenté
           end select
        enddo

    end subroutine compute_planeW_Exafield

    subroutine newmark_predictor_solid(dom)
        type(domain_solid), intent (INOUT) :: dom
        !
        dom%champs1%Depla = dom%champs0%Depla
        dom%champs1%Veloc = dom%champs0%Veloc
        dom%champs1%Forces = 0d0
    end subroutine newmark_predictor_solid

    subroutine newmark_corrector_solid(dom, dt)
        type(domain_solid), intent (INOUT) :: dom
        double precision :: dt
        !
        integer :: i_dir, n, indpml
        do i_dir = 0,2
            do n = 0,dom%nglltot-1
                dom%champs0%Forces(n,i_dir) = dom%champs1%Forces(n,i_dir) * dom%MassMat(n)
                dom%champs0%Veloc(n,i_dir) = dom%champs0%Veloc(n,i_dir) + dt * dom%champs0%Forces(n,i_dir)
            end do
        enddo
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs0%Veloc(indpml,:) = 0.
        enddo
        dom%champs0%Depla = dom%champs0%Depla + dt * dom%champs0%Veloc
    end subroutine newmark_corrector_solid

    function solid_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        type(domain_solid), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed, M
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        M = dom%Lambda_(i,j,k,bnum,ee) + 2.*dom%Mu_(i,j,k,bnum,ee)
        Pspeed = sqrt(M/dom%Density_(i,j,k,bnum,ee))
    end function solid_Pspeed
end module dom_solid

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
