!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_solid
    use constants
    use champs_solid
    implicit none
#include "index.h"
#include "gllopt.h"
contains

    subroutine allocate_champs_solid(dom, i)
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: i

        ! Allocate ONE more gll than necessary to use as a dummy target for
        ! indirections for fake elements.
        allocate(dom%champs(i)%Depla (0:dom%nglltot,0:2))
        allocate(dom%champs(i)%Veloc (0:dom%nglltot,0:2))

        dom%champs(i)%Depla = 0d0
        dom%champs(i)%Veloc = 0d0
    end subroutine allocate_champs_solid
    subroutine allocate_dom_solid (Tdomain, dom)
        use sdomain
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_solid), intent (INOUT) :: dom
        !
        integer :: nbelem, nblocks, ngll, n_solid, i
        integer :: ntemps, nprops
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
        ! Mirror
        !!! GB dom%use_mirror = Tdomain%use_mirror
        dom%mirror_type = Tdomain%mirror_type

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! Do not allocate if not needed (save allocation/RAM)
            ! Wexo can have glls without elements
            nblocks = dom%nblocks

            ! number of material properties
            ! indexed as 0:nprops with 0 being rho
            if (nl_flag) then
                ! rho, lam, mu, kappa, qp, qs + 5 (Syld,CKin,KKin,Rinf,Biso)
                nprops = 10
            else
                if (aniso) then
                    ! rho, lam, mu, kappa, qp, qs + 21
                    nprops = cAniso+20
                else
                    if (n_solid>0) then
                        ! rho, lam, mu, kappa, qp, qs
                        nprops = 5
                    else
                        ! rho, lam, mu
                        nprops = 2
                    endif
                endif
            endif
            dom%nprops = nprops
            allocate (dom%props (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nprops, 0:nblocks-1))
#if defined(OPENACC) || TEST_FORCE==1
            ntemps = nblocks
#else
            ntemps = 1
#endif
            allocate(dom%Forces (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1,0:2,0:ntemps-1))
            allocate(dom%Depla  (0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1,0:2,0:ntemps-1))
            ! TODO: only for mirrors
            allocate(dom%Veloc(0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1,0:2,0:ntemps-1))
            allocate(dom%Sigma(0:VCHUNK-1, 0:ngll-1, 0:ngll-1, 0:ngll-1,0:5,0:ntemps-1))
            dom%props = 0. ! May not be initialized if run without attenuation

            if (nl_flag) then
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

            if (n_solid>0) then
                !!! <GB> temporary modif
                !!!
                !!! if (aniso) then
                !!!     allocate (dom%Q_(0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                !!! else
                !!!     allocate (dom%Qs_         (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                !!!     allocate (dom%Qp_         (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                !!!     allocate (dom%onemPbeta_  (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                !!!     allocate (dom%epsilonvol_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                !!!     dom%epsilonvol_(:,:,:,:,:) = 0
                !!!     allocate (dom%R_vol_       (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                !!!     dom%R_vol_(:,:,:,:,:,:) = 0
                !!! endif
                allocate (dom%onemPbeta_  (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%epsilonvol_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                dom%epsilonvol_(:,:,:,:,:) = 0
                allocate (dom%R_vol_       (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
                dom%R_vol_(:,:,:,:,:,:) = 0
                !!! </GB>
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
                allocate (dom%R_xx_ (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%R_yy_ (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%R_xy_ (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%R_xz_ (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
                allocate (dom%R_yz_ (0:n_solid-1, 0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
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
            do i=0,Tdomain%TimeD%nsubsteps
                call allocate_champs_solid(dom, i)
            end do
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_solid

    subroutine deallocate_dom_solid (dom,nl_flag)
        implicit none
        type(domain_solid), intent (INOUT) :: dom
        logical, intent(in) :: nl_flag
        integer :: i

        if(allocated(dom%props)) deallocate(dom%props)

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

        do i=0,size(dom%champs)-1
            if(allocated(dom%champs(i)%Depla )) deallocate(dom%champs(i)%Depla )
            if(allocated(dom%champs(i)%Veloc )) deallocate(dom%champs(i)%Veloc )
        end do
        ! nonlinear parameters TODO
        if (nl_flag) then
            ! nonlinear internal variables
            if(allocated(dom%m_radius)) deallocate(dom%m_radius)
            if(allocated(dom%m_stress)) deallocate(dom%m_stress)
            if(allocated(dom%m_center)) deallocate(dom%m_center)
            if(allocated(dom%m_strain)) deallocate(dom%m_strain)
            if(allocated(dom%m_plstrain)) deallocate(dom%m_plstrain)
        endif

        call deallocate_dombase(dom)

        end subroutine deallocate_dom_solid

    subroutine get_solid_dom_grad_lambda(dom, lnum, grad_La)
        !$acc routine worker
        use deriv3d
        implicit none
        !
        type(domain_solid), intent(inout)     :: dom
        integer, intent(in)                   :: lnum
        real(fpp), intent(inout), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: grad_La 
        integer                  :: ngll, i, j, k
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: xlambda
        real(fpp)                :: DLambdaDX,DLambdaDY,DLambdaDZ
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        ngll = dom%ngll

        grad_La = 0d0
        xlambda = dom%props(ee,:,:,:,cLambda,bnum) 
        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    invgrad_ijk = dom%InvGrad_(:,:,i,j,k,bnum,ee) ! cache for performance
                    call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                        invgrad_ijk,xlambda,DLambdaDX,DLambdaDY,DLambdaDZ)
                    grad_La(i,j,k,0)=DLambdaDX
                    grad_La(i,j,k,1)=DLambdaDY
                    grad_La(i,j,k,2)=DLambdaDZ
                end do
            end do
        end do
        return
    end subroutine get_solid_dom_grad_lambda

    subroutine get_solid_dom_grad_mu(dom, lnum, grad_Mu)
        !$acc routine worker
        use deriv3d
        implicit none
        !
        type(domain_solid), intent(inout)     :: dom
        integer, intent(in)                   :: lnum
        real(fpp), intent(inout), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: grad_Mu
        integer                  :: ngll, i, j, k
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: xmu
        real(fpp)                :: DMuDX,DMuDY,DMuDZ
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        ngll = dom%ngll

        grad_Mu = 0d0
        xmu = dom%props(ee,:,:,:,cMu,bnum)
        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    invgrad_ijk = dom%InvGrad_(:,:,i,j,k,bnum,ee) ! cache for performance
                    call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                        invgrad_ijk,xmu,DMuDX,DMuDY,DMuDZ)
                    grad_Mu(i,j,k,0)=DMuDX
                    grad_Mu(i,j,k,1)=DMuDY
                    grad_Mu(i,j,k,2)=DMuDZ
                end do
            end do
        end do
        return
    end subroutine get_solid_dom_grad_mu

    subroutine get_solid_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy,&
        eps_vol, eps_dev, sig_dev, dUdX, nl_flag, eps_dev_pl)
        !$acc routine worker
        use deriv3d
        implicit none
        !
        type(domain_solid), intent(inout)     :: dom
        integer, intent(in), dimension(0:)    :: out_variables
        integer, intent(in)                   :: lnum
        logical, intent(in)                   :: nl_flag
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldU, fieldV, fieldA
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: fieldP
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: P_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: S_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:8) :: dUdX
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:6) :: eps_dev_pl
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: sig_dev
        !
        logical                  :: flag_gradU
        integer                  :: ngll, i, j, k, ind
        real(fpp)                :: EXY, EXZ, EYZ
        real(fpp)                :: DXX, DXY, DXZ
        real(fpp)                :: DYX, DYY, DYZ
        real(fpp)                :: DZX, DZY, DZZ
        real(fpp)                :: divU
        real(fpp)                :: xmu, xlambda, xkappa
        real(fpp)                :: x2mu, xlambda2mu
        real(fpp)                :: onemSbeta, onemPbeta
        real(fpp), dimension(0:20) :: CC
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        !
        integer :: bnum, ee, flag_gradUint

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        xlambda = 0d0
        xmu = 0d0
        x2mu = 0d0
        flag_gradUint = 0
        flag_gradUint = out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS)

        if (.not.nl_flag) then
            flag_gradUint = flag_gradUint     + &
                out_variables(OUT_PRESSION)   + &
                out_variables(OUT_DUDX)       + &
                out_variables(OUT_EPS_VOL)    + &
                out_variables(OUT_EPS_DEV)    + &
                out_variables(OUT_STRESS_DEV)
        endif

        flag_gradU = flag_gradUint /= 0
        ngll = dom%ngll

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        fieldU(i,j,k,:) = dom%champs(0)%Depla(ind,:)
                    end if
                end do
            end do
        end do

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    if (flag_gradU .or. (out_variables(OUT_DUDX) == 1)) then
                        invgrad_ijk = dom%InvGrad_(:,:,i,j,k,bnum,ee) ! cache for performance
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,0),DXX,DYX,DZX)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,1),DXY,DYY,DZY)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,2),DXZ,DYZ,DZZ)
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
                    end if
                    if (out_variables(OUT_VITESSE) == 1) fieldV(i,j,k,:) = dom%champs(0)%Veloc(ind,:)
                    if (out_variables(OUT_ACCEL) == 1) fieldA(i,j,k,:) = dom%champs(1)%Veloc(ind,:)
                    if (flag_gradU .and. .not. nl_flag) then
                        ! PRESSION
                        if (out_variables(OUT_PRESSION) == 1) then
                            xkappa = dom%props(ee,i,j,k,cLambda,bnum) + &
                                two*M_1_3*dom%props(ee,i,j,k,cMu,bnum)
                            fieldP(i,j,k) = -xkappa*divU
                        endif
                        ! EPSVOL
                        if (out_variables(OUT_EPS_VOL) == 1) then
                            eps_vol(i,j,k) = divU
                        end if
                        ! ELASTIC-VISCOELASTIC MODULA
                        if (out_variables(OUT_ENERGYP) == 1 .or. &
                            out_variables(OUT_ENERGYS) == 1 .or. &
                            out_variables(OUT_STRESS_DEV) == 1) then
                            if (dom%aniso) then
                                CC = dom%props(ee,i,j,k,1:21,bnum)
                            else
                                xmu     = dom%props(ee,i,j,k,cMu,bnum)
                                xlambda = dom%props(ee,i,j,k,cLambda,bnum)
                                xkappa  = dom%props(ee,i,j,k,cKappa,bnum)
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
                        ! P-ENERGY
                        if (out_variables(OUT_ENERGYP) == 1) then !XXX: aniso
                            P_energy(i,j,k) = ((0.5d0*xlambda) + xmu) * divU**2d0
                        end if
                        ! S-ENERGY
                        if (out_variables(OUT_ENERGYS) == 1) then
                            S_energy(i,j,k) = xmu/2.0d0 * (           &
                                                    (DZY - DYZ)**2d0  &
                                                  + (DXZ - DZX)**2d0  &
                                                  + (DYX - DXY)**2d0  &
                                                  )
                        end if
                        ! DEVIATORIC STRAIN
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
                    end if
                    if (nl_flag) then
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
                        end if
                        ! TOTAL STRAIN
                        if (out_variables(OUT_EPS_DEV) == 1) then
                            dUdX(i,j,k,0) = dom%strain_(0,i,j,k,bnum,ee)  !DXX
                            dUdX(i,j,k,1) = dom%strain_(3,i,j,k,bnum,ee)  !DXY
                            dUdX(i,j,k,2) = dom%strain_(4,i,j,k,bnum,ee)  !DXZ
                            dUdX(i,j,k,3) = dom%strain_(3,i,j,k,bnum,ee)  !DYX
                            dUdX(i,j,k,4) = dom%strain_(1,i,j,k,bnum,ee)  !DYY
                            dUdX(i,j,k,5) = dom%strain_(5,i,j,k,bnum,ee)  !DYZ
                            dUdX(i,j,k,6) = dom%strain_(4,i,j,k,bnum,ee)  !DZX
                            dUdX(i,j,k,7) = dom%strain_(5,i,j,k,bnum,ee)  !DZY
                            dUdX(i,j,k,8) = dom%strain_(2,i,j,k,bnum,ee)  !DZZ
                        endif
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
                    end if
                enddo
            enddo
        enddo
        !
    end subroutine get_solid_dom_var


    subroutine get_solid_dom_elem_energy(dom, lnum, P_energy, S_energy, R_energy, C_energy)
        !$acc routine worker
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
        real(fpp) :: dUx_dx,dUx_dy,dUx_dz
        real(fpp) :: dUy_dx,dUy_dy,dUy_dz
        real(fpp) :: dUz_dx,dUz_dy,dUz_dz
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
                    fieldU(i,j,k,:) = dom%champs(0)%Depla(ind,:)
                    fieldV(i,j,k,:) = dom%champs(0)%Veloc(ind,:)
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

                    xmu      = dom%props(ee,i,j,k,cMu,bnum)
                    xlambda  = dom%props(ee,i,j,k,cLambda,bnum)
                    xkappa   = dom%props(ee,i,j,k,cKappa,bnum)
                    xdensity = dom%props(ee,i,j,k,cRho,bnum)
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

    subroutine init_domain_solid(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solid), intent(inout) :: dom

        dom%dt = Tdomain%TimeD%dtmin
    end subroutine init_domain_solid

    subroutine start_domain_solid(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solid), intent(inout) :: dom
        !
        integer :: i, ns

        !$acc  enter data copyin(dom, dom%champs) &
        !$acc&            copyin(dom%MassMat, dom%dirich) &
        !$acc&      copyin(dom%props, dom%gllw, dom%hprime)&
        !$acc&      copyin(dom%m_Idom, dom%m_Jacob, dom%m_InvGrad) &
        !$acc&      create(dom%Forces, dom%Depla, dom%Sigma)
        do i = 0,1
            !$acc enter data  copyin(dom%champs(i)%Depla, dom%champs(i)%veloc)
        end do
        

        do ns = 0, Tdomain%n_source-1
            if(Tdomain%rank == Tdomain%sSource(ns)%proc)then
                if(Tdomain%sSource(ns)%i_type_source == 1 .or. Tdomain%sSource(ns)%i_type_source == 2) then
                    !$acc enter data      copyin(Tdomain%sSource(ns)%ExtForce)
                end if
            end if
        end do

    end subroutine start_domain_solid

    subroutine stop_domain_solid(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solid), intent(inout) :: dom
        !
        integer :: i, ns

        do ns = 0, Tdomain%n_source-1
            if(Tdomain%rank == Tdomain%sSource(ns)%proc)then
                if(Tdomain%sSource(ns)%i_type_source == 1 .or. Tdomain%sSource(ns)%i_type_source == 2) then
                    !$acc exit data      delete(Tdomain%sSource(ns)%ExtForce)
                end if
            end if
        end do

        do i = 0,1
            !$acc exit data  delete(dom%champs(i)%Depla, Tdomain%sdom%champs(i)%veloc)
        end do

        !$acc  exit data  delete(dom%MassMat, dom%dirich) &
        !$acc&            delete(dom%props, dom%gllw, dom%hprime)&
        !$acc&            delete(dom%m_Idom, dom%m_Jacob, dom%m_InvGrad) &
        !$acc&            delete(dom%Forces, dom%Depla, dom%Sigma) &
        !$acc&            delete(dom%champs, dom)
    end subroutine stop_domain_solid

    
    subroutine init_material_properties_solid(dom, lnum, mat, density, lambda, mu, &
        Qkappa, Qmu, nlkp, nl_flag)
        use ssubdomains
        type(domain_solid), intent(inout) :: dom
        integer, intent(in) :: lnum
        type (subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: density
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: lambda
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: mu
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Qkappa
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Qmu
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: nlkp
        logical, intent(in) :: nl_flag
        real(fpp),parameter :: gamma_el = real(1.0d-5,fpp)
        real(fpp),parameter :: gamma_pl = real(1.0d-4,fpp)
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        dom%props(ee,:,:,:,cRho,bnum) = density
        dom%props(ee,:,:,:,cLambda,bnum) = lambda
        dom%props(ee,:,:,:,cMu,bnum) = mu
        
        if (nl_flag) then
            if (mat%deftype.eq.MATDEF_NLKP_VS_RHO) then
                dom%props(ee,:,:,:,cSyld,bnum) = gamma_el*mu*sqrt(3.0d0)
                dom%props(ee,:,:,:,cCKin,bnum) = mu
                dom%props(ee,:,:,:,cKKin,bnum) = 1.0d0/(sqrt(3.0d0)*(nlkp-gamma_el))
                dom%props(ee,:,:,:,cRinf,bnum) = mat%DRinf 
                dom%props(ee,:,:,:,cBiso,bnum) = mat%DBiso
            else
                dom%props(ee,:,:,:,cSyld,bnum) = gamma_el*nlkp*sqrt(3.0d0)
                dom%props(ee,:,:,:,cCKin,bnum) = nlkp
                dom%props(ee,:,:,:,CKKin,bnum) = 1.0d0/(sqrt(3.0d0)*(nlkp-gamma_el))
                dom%props(ee,:,:,:,CRinf,bnum) = 0.0D0
                dom%props(ee,:,:,:,CBiso,bnum) = 0.0D0
            endif
        endif

        if (dom%n_sls>0)  then
            dom%props(ee,:,:,:,cKappa,bnum) = lambda + 2d0*mu/3d0
            dom%props(ee,:,:,:,cQp,bnum) = Qkappa
            dom%props(ee,:,:,:,cQs,bnum) = Qmu
        endif

    end subroutine init_material_properties_solid

    subroutine init_material_tensor_solid(dom, lnum, mat, density, lambda, mu, Qkappa, Qmu, Cij)
        use ssubdomains
        type(domain_solid), intent(inout) :: dom
        integer, intent(in) :: lnum
        type (subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: density
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: lambda
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: mu
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Qkappa, Qmu
        real(fpp), dimension(1:6,1:6,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1), intent(in) :: Cij

        integer :: idef, ii, jj
        integer :: i, j, k
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        dom%props(ee,:,:,:,cRho   ,bnum) = density
        dom%props(ee,:,:,:,cLambda,bnum) = lambda
        dom%props(ee,:,:,:,cMu    ,bnum) = mu

        do i=0,dom%ngll-1
            do j=0,dom%ngll-1
                do k=0,dom%ngll-1
                    idef = 0
                    do ii = 1,6
                        do jj = ii,6
                            dom%props(ee,i,j,k,cAniso+idef,bnum) = Cij(ii,jj,i,j,k)
                            idef = idef + 1
                        enddo
                    enddo
                enddo
            enddo
        enddo

        if (dom%n_sls>0)  then
            dom%props(ee,:,:,:,cKappa,bnum) = lambda + 2d0*mu/3d0
            !!! <GB> uniform cover of att
            dom%props(ee,:,:,:,cQs,bnum) = Qmu
            dom%props(ee,:,:,:,cQp,bnum) = Qkappa
            !!!if (dom%aniso) then
            !!!    dom%Q_(:,:,:,bnum,ee) = Qmu
            !!!else
            !!!    dom%Qs_(:,:,:,bnum,ee) = Qmu
            !!!    dom%Qp_(:,:,:,bnum,ee) = Qkappa
            !!!endif
            !!! <GB>
        endif

    end subroutine init_material_tensor_solid

    subroutine init_local_mass_solid(dom,specel,i,j,k,ind,Whei)
        use selement
        type(domain_solid), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real(fpp) :: Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ! Solid.

        specel%MassMat(i,j,k) = Whei*dom%props(ee,i,j,k,cRho,bnum)*dom%Jacob_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_solid
      
    subroutine forces_int_solid_mainloop(dom, i0, i1, nlflag, m_dump, m_load, m_expl, m_recalc)
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(IN) :: i0, i1
        logical, intent(IN) :: nlflag
        logical, intent(IN) :: m_dump, m_load, m_expl, m_recalc
        !
        integer :: n_solid
        logical :: aniso
        n_solid = dom%n_sls
        aniso = dom%aniso

        if (nlflag) then
            if (n_solid>0) then
                call dispatch_forces_int_aniso_nl(dom, dom%champs(i0), dom%champs(i1))
            else
                call dispatch_forces_int_iso_nl(dom, dom%champs(i0), dom%champs(i1))
            endif
        else if (m_dump .or. m_load) then
            !!! WITH MIRROR
            if (aniso) then
                if (n_solid>0) then
                    call dispatch_forces_int_mirror_aniso_atn(dom, dom%champs(i0), dom%champs(i1), m_dump, m_expl, m_recalc)
                else
                    call dispatch_forces_int_mirror_aniso(dom, dom%champs(i0), dom%champs(i1), m_dump, m_expl, m_recalc)
                end if
            else
                if (n_solid>0) then
                    call dispatch_forces_int_mirror_iso_atn(dom, dom%champs(i0), dom%champs(i1), m_dump, m_expl, m_recalc)
                else
                    call dispatch_forces_int_mirror_iso(dom, dom%champs(i0), dom%champs(i1), m_dump, m_expl, m_recalc)
                endif
            end if            
        else
            if (aniso) then
                if (n_solid>0) then
                    call dispatch_forces_int_aniso_atn(dom, dom%champs(i0), dom%champs(i1))
                else
                    call dispatch_forces_int_aniso(dom, dom%champs(i0), dom%champs(i1))
                end if
            else
                if (n_solid>0) then
                    call dispatch_forces_int_iso_atn(dom, dom%champs(i0), dom%champs(i1))
                else
                    call dispatch_forces_int_iso(dom, dom%champs(i0), dom%champs(i1))
                endif
            end if

        endif
    end subroutine forces_int_solid_mainloop

    subroutine dispatch_forces_int_iso(dom, var, dvdt)
        use m_calcul_forces_iso
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_iso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_5(calcul_forces_iso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_6(calcul_forces_iso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_7(calcul_forces_iso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_8(calcul_forces_iso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_9(calcul_forces_iso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_N(calcul_forces_iso,,(dom,dom%ngll,var,dvdt))
        end select
    end subroutine dispatch_forces_int_iso

    subroutine dispatch_forces_int_iso_atn(dom, var, dvdt)
        use m_calcul_forces_iso_atn
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_iso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_5(calcul_forces_iso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_6(calcul_forces_iso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_7(calcul_forces_iso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_8(calcul_forces_iso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_9(calcul_forces_iso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_N(calcul_forces_iso_atn,,(dom,dom%ngll,var,dvdt))
        end select
    end subroutine dispatch_forces_int_iso_atn
    
    subroutine dispatch_forces_int_aniso(dom, var, dvdt)
        use m_calcul_forces_aniso
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_aniso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_5(calcul_forces_aniso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_6(calcul_forces_aniso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_7(calcul_forces_aniso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_8(calcul_forces_aniso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_9(calcul_forces_aniso,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_N(calcul_forces_aniso,,(dom,dom%ngll,var,dvdt))
        end select
    end subroutine dispatch_forces_int_aniso

     subroutine dispatch_forces_int_aniso_atn(dom, var, dvdt)
        use m_calcul_forces_aniso_atn
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_aniso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_5(calcul_forces_aniso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_6(calcul_forces_aniso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_7(calcul_forces_aniso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_8(calcul_forces_aniso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_9(calcul_forces_aniso_atn,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_N(calcul_forces_aniso_atn,,(dom,dom%ngll,var,dvdt))
        end select
    end subroutine dispatch_forces_int_aniso_atn

    subroutine dispatch_forces_int_iso_nl(dom, var, dvdt)
        use m_calcul_forces_iso_nl
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_iso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_5(calcul_forces_iso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_6(calcul_forces_iso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_7(calcul_forces_iso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_8(calcul_forces_iso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_9(calcul_forces_iso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_N(calcul_forces_iso_nl,,(dom,dom%ngll,var,dvdt))
        end select
    end subroutine dispatch_forces_int_iso_nl

    subroutine dispatch_forces_int_aniso_nl(dom, var, dvdt)
        use m_calcul_forces_aniso_nl
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_aniso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_5(calcul_forces_aniso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_6(calcul_forces_aniso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_7(calcul_forces_aniso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_8(calcul_forces_aniso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_9(calcul_forces_aniso_nl,,(dom,dom%ngll,var,dvdt))
            NGLLDISPATCHCALL_N(calcul_forces_aniso_nl,,(dom,dom%ngll,var,dvdt))
        end select
    end subroutine dispatch_forces_int_aniso_nl

    subroutine dispatch_forces_int_mirror_iso(dom, var, dvdt, m_dump, m_expl, m_recalc)
        use m_calcul_forces_iso
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        logical, intent(IN) :: m_dump, m_expl, m_recalc
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_iso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_5(calcul_forces_iso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_6(calcul_forces_iso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_7(calcul_forces_iso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_8(calcul_forces_iso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_9(calcul_forces_iso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_N(calcul_forces_iso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
        end select
    end subroutine dispatch_forces_int_mirror_iso

    subroutine dispatch_forces_int_mirror_aniso(dom, var, dvdt, m_dump, m_expl, m_recalc)
        use m_calcul_forces_aniso
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        logical, intent(IN) :: m_dump, m_expl, m_recalc
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_aniso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_5(calcul_forces_aniso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_6(calcul_forces_aniso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_7(calcul_forces_aniso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_8(calcul_forces_aniso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_9(calcul_forces_aniso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_N(calcul_forces_aniso,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
        end select
    end subroutine dispatch_forces_int_mirror_aniso
    
    subroutine dispatch_forces_int_mirror_iso_atn(dom, var, dvdt, m_dump, m_expl, m_recalc)
        use m_calcul_forces_iso_atn
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        logical, intent(IN) :: m_dump, m_expl, m_recalc
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_iso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_5(calcul_forces_iso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_6(calcul_forces_iso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_7(calcul_forces_iso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_8(calcul_forces_iso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_9(calcul_forces_iso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_N(calcul_forces_iso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
        end select
    end subroutine dispatch_forces_int_mirror_iso_atn

    subroutine dispatch_forces_int_mirror_aniso_atn(dom, var, dvdt, m_dump, m_expl, m_recalc)
        use m_calcul_forces_aniso_atn
        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(in) :: var
        type(champssolid), intent(inout) :: dvdt
        logical, intent(IN) :: m_dump, m_expl, m_recalc
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_aniso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_5(calcul_forces_aniso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_6(calcul_forces_aniso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_7(calcul_forces_aniso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_8(calcul_forces_aniso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_9(calcul_forces_aniso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
            NGLLDISPATCHCALL_N(calcul_forces_aniso_atn,_mirror,(dom,dom%ngll,var,dvdt,m_dump,m_expl,m_recalc))
        end select
    end subroutine dispatch_forces_int_mirror_aniso_atn

    


!!    ! XXX WTF?
!!    subroutine compute_planeW_Exafield(lnum,ctime,Tdomain,f)
!!        use sdomain
!!        use Surface_prbl_type
!!
!!        implicit none
!!        type(domain),            intent(inout) :: Tdomain
!!        real(fpp),               intent(in  )  :: ctime
!!        integer,                 intent(in   ) :: lnum
!!        integer,                 intent(in)    :: f
!!        !
!!        real(fpp), dimension(0:2):: coord, displ, veloc, accel
!!        real(fpp)                :: PWspeed
!!        character(len=20)           :: char
!!        integer                     :: ns, im, i, j, k, dom, ipw, ngll
!!        integer                     :: bnum, ee, ind, ss
!!
!!        bnum = lnum/VCHUNK
!!        ee = mod(lnum,VCHUNK)
!!
!!        do ns=1,size(Tdomain%list_PWBC)
!!            ipw     = Tdomain%list_PWBC(ns)
!!            im      = Tdomain%nsurfsource(ipw)%mat_index
!!            dom     = Tdomain%sSubDomain(im)%dom
!!            ngll     = Tdomain%sSubDomain(im)%NGLL
!!            write(char,*) Tdomain%nsurfsource(ipw)%index(1)
!!
!!            do ss=0,size(Tdomain%sSurfaces)-1
!!                if (Tdomain%sSurfaces(ss)%name=="surface"//adjustl(char(:len_trim(char)))) &
!!                    exit
!!            end do
!!            PWspeed = Tdomain%sSurfaces(ss)%Elastic%PWspeed
!!            select case (dom)
!!            case (DM_SOLID_CG)
!!                ngll     = Tdomain%sdom%ngll
!!                do k=0,ngll-1
!!                    do j=0,ngll-1
!!                        do i=0,ngll-1
!!                            ind = Tdomain%sdom%Idom_(i,j,k,bnum,ee)
!!                            coord = Tdomain%GlobCoord(:,ind)-Tdomain%nsurfsource(ipw)%scoord(:)
!!                            call PlaneWavedispl(Tdomain%nsurfsource(ipw),coord,ctime,PWspeed, displ,veloc,accel)
!!                            Tdomain%sdom%champs(f)%Depla(ind,:)  = Tdomain%sdom%champs(f)%Depla(ind,:) + displ
!!                            Tdomain%sdom%champs(f)%Veloc(ind,:)  = Tdomain%sdom%champs(f)%Veloc(ind,:) + veloc
!!                            Tdomain%sdom%champs(f)%Forces(ind,:) = Tdomain%sdom%champs(f)%Forces(ind,:) + accel
!!                        enddo
!!                    enddo
!!                enddo
!!            case(DM_FLUID_CG)
!!                ! pas encore implémenté
!!            end select
!!        enddo
!!
!!    end subroutine compute_planeW_Exafield

    subroutine newmark_predictor_solid(dom, f0, f1)
        type(domain_solid), intent (INOUT) :: dom
        integer :: f0, f1
        !$acc kernels async(1)
        dom%champs(f1)%Veloc = 0d0
        !$acc end kernels
    end subroutine newmark_predictor_solid


    ! predictor (1.veloc)
    ! forces    (wait 1.veloc)  -> host(1.veloc)
    ! source    1.veloc sur host
    ! corrector (device 1.veloc)

    
    subroutine newmark_corrector_solid(dom, dt, f0, f1)
        type(domain_solid), intent (INOUT) :: dom
        real(fpp):: dt
        integer, intent(in) :: f0, f1
        !
        integer :: n, idx, count

        count = dom%n_dirich
        !!$acc update device(dom%champs(f1)%Veloc) async(1)
        !$acc parallel loop  async(1)
        do n = 0, count-1
            idx = dom%dirich(n)
            dom%champs(f1)%Veloc(idx,:) = 0.
            !dom%champs(f1)%Depla(idx,:) = 0.
        enddo
        !$acc end parallel

        count = dom%nglltot
        call newmark_corrector_time_scheme(count,dt,dom%MassMat,dom%champs(f0)%Depla, &
            dom%champs(f0)%Veloc,dom%champs(f1)%Veloc)
    end subroutine newmark_corrector_solid

    subroutine newmark_corrector_time_scheme(count,dt,invmass,depla,veloc,accel)
        integer, intent(in) :: count
        real(fpp), intent(in) :: dt
        real(fpp), dimension(0:count,0:2), intent(inout) :: depla, veloc, accel
        real(fpp), dimension(0:count), intent(in) :: invmass
        !
        integer :: n, i_dir
        real(fpp) :: acc, vel
        !$acc parallel loop collapse(2)  async(1)
        do i_dir = 0,2
            !$omp simd linear(n)
            do n = 0,count-1
                acc = accel(n,i_dir) * invmass(n)
                vel = veloc(n,i_dir)
                accel(n,i_dir) = acc
                vel = vel + dt * acc
                veloc(n,i_dir) = vel
                depla(n,i_dir) = depla(n,i_dir) + dt * vel
            end do
        enddo
        !$acc end parallel
   
    end subroutine newmark_corrector_time_scheme

    function solid_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        type(domain_solid), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed, M
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        M = dom%props(ee,i,j,k,cLambda,bnum) + 2.*dom%props(ee,i,j,k,cMu,bnum)
        Pspeed = sqrt(M/dom%props(ee,i,j,k,cRho,bnum))
    end function solid_Pspeed


    subroutine lddrk_init_solid(dom, f0)
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: f0

        if (dom%nglltot==0) return
        dom%champs(f0)%Veloc = 0d0 ! XXX
    end subroutine lddrk_init_solid

    subroutine lddrk_update_solid(dom, f0, f1, f2, dt, cb, cg)
        
        type(domain_solid), intent (INOUT) :: dom
        integer, intent(in) :: f0, f1, f2
        real(fpp), intent(in) :: cb, cg, dt
        integer :: i, n
        
        ! f2  contains forces computation  ie f2 = dU/dt
        ! f1 : w(n+1) = cb*w(n) + dt*dU/dt
        ! f0 : U(n+1) = U(n) + cg*w(n+1)
        ! Only solid for starters
        do i = 0,2
            do n = 0,dom%nglltot-1
                dom%champs(f2)%Veloc(n,i) =    dom%champs(f2)%Veloc(n,i) * dom%MassMat(n)
                dom%champs(f1)%Veloc(n,i) = cb*dom%champs(f1)%Veloc(n,i) + dt*dom%champs(f2)%Veloc(n,i)
                dom%champs(f1)%Depla(n,i) = cb*dom%champs(f1)%Depla(n,i) + dt*dom%champs(f0)%Veloc(n,i)
                dom%champs(f0)%Depla(n,i) =    dom%champs(f0)%Depla(n,i) + cg*dom%champs(f1)%Depla(n,i)
                dom%champs(f0)%Veloc(n,i) =    dom%champs(f0)%Veloc(n,i) + cg*dom%champs(f1)%Veloc(n,i)
            end do
        end do
        
    end subroutine lddrk_update_solid

    subroutine apply_source_solid(src, dom, i1, ft, lnum)
        use sdomain
        implicit none
        type(domain_solid),intent(inout) :: dom
        type(Source),intent(inout) :: src 
        real(fpp), intent(in) :: ft
        integer, intent(in) :: i1
        integer, intent(in) :: lnum
        integer :: bnum, ee
        integer :: ngll, i, j, k, i_dir, idx
        real(kind=fpp) :: force

        ngll = dom%ngll
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        !$acc parallel loop gang vector collapse(4) async(1)
        do i_dir = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx   = dom%Idom_(i,j,k,bnum,ee)
                        force = dom%champs(i1)%Veloc(idx, i_dir) + ft*src%ExtForce(i,j,k,i_dir)
                        dom%champs(i1)%Veloc(idx, i_dir) = force
                    enddo
                enddo
            enddo
        enddo
    end subroutine apply_source_solid
end module dom_solid
