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
        logical :: aniso
        !

        ngll    = dom%ngll
        nbelem  = dom%nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call compute_gll_data(ngll, dom%gllc, dom%gllw, dom%hprime, dom%htprime)

        aniso   = Tdomain%aniso
        n_solid = Tdomain%n_sls
        dom%n_sls = n_solid
        dom%aniso = Tdomain%aniso

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! Do not allocate if not needed (save allocation/RAM)
            ! Wexo can have glls without elements
            nblocks = ((nbelem+VCHUNK-1)/VCHUNK)
            dom%nblocks = nblocks

            allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Mu_     (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Kappa_  (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))

            allocate (dom%Jacob_  (        0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate (dom%InvGrad_(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1, 0:VCHUNK-1))

            allocate(dom%Idom_(0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            dom%m_Idom = 0

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
            allocate(dom%champs0%Forces(0:dom%nglltot-1,0:2))
            allocate(dom%champs0%Depla (0:dom%nglltot-1,0:2))
            allocate(dom%champs0%Veloc (0:dom%nglltot-1,0:2))
            allocate(dom%champs1%Forces(0:dom%nglltot-1,0:2))
            allocate(dom%champs1%Depla (0:dom%nglltot-1,0:2))
            allocate(dom%champs1%Veloc (0:dom%nglltot-1,0:2))

            dom%champs0%Forces = 0d0
            dom%champs0%Depla = 0d0
            dom%champs0%Veloc = 0d0

            ! Allocation de MassMat pour les solides
            allocate(dom%MassMat(0:dom%nglltot-1))
            dom%MassMat = 0d0
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_solid

    subroutine deallocate_dom_solid (dom)
        implicit none
        type(domain_solid), intent (INOUT) :: dom

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_Mu     )) deallocate(dom%m_Mu     )
        if(allocated(dom%m_Kappa  )) deallocate(dom%m_Kappa  )

        if(allocated(dom%m_Jacob  )) deallocate(dom%m_Jacob  )
        if(allocated(dom%m_InvGrad)) deallocate(dom%m_InvGrad)

        if(allocated(dom%m_Idom)) deallocate(dom%m_Idom)

        if(allocated(dom%gllc))    deallocate(dom%gllc)
        if(allocated(dom%gllw))    deallocate(dom%gllw)
        if(allocated(dom%hprime))  deallocate(dom%hprime)
        if(allocated(dom%htprime)) deallocate(dom%htprime)

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

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)
    end subroutine deallocate_dom_solid

    subroutine get_solid_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        use deriv3d
        implicit none
        !
        type(domain_solid), intent(inout)          :: dom
        integer, intent(in), dimension(0:8)        :: out_variables
        integer, intent(in)                        :: lnum
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
        real(fpp)                :: xmu, xlambda, xkappa
        real(fpp)                :: x2mu, xlambda2mu
        real(fpp)                :: onemSbeta, onemPbeta
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
                        fieldA(i,j,k,:) = dom%champs0%Forces(ind,:)
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
                        fieldP(i,j,k) = -(dom%Lambda_(i,j,k,bnum,ee)&
                                          +2d0/3d0*dom%Mu_(i,j,k,bnum,ee))&
                                        *(DXX+DYY+DZZ)
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1 .or. &
                        out_variables(OUT_ENERGYP) == 1 .or. &
                        out_variables(OUT_EPS_DEV) == 1 .or. &
                        out_variables(OUT_STRESS_DEV) == 1 ) then
                        if(.not. allocated(eps_vol)) allocate(eps_vol(0:ngll-1,0:ngll-1,0:ngll-1))
                        eps_vol(i,j,k) = DXX + DYY + DZZ
                    end if

                    if (.not. dom%aniso) then
                        xmu     = dom%Mu_    (i,j,k,bnum,ee)
                        xlambda = dom%Lambda_(i,j,k,bnum,ee)
                        xkappa  = dom%Kappa_ (i,j,k,bnum,ee)

                        if (dom%n_sls>0) then
                            onemSbeta = dom%onemSbeta_(i,j,k,bnum,ee)
                            onemPbeta = dom%onemPbeta_(i,j,k,bnum,ee)
                            xmu    = xmu * onemSbeta
                            xkappa = xkappa * onemPbeta
                        endif
                        x2mu       = 2. * xmu
                        xlambda2mu = xlambda + x2mu
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        if(.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                        P_energy(i,j,k) = 0.
                        if (.not. dom%aniso) then
                            P_energy(i,j,k) = .5 * xlambda2mu * eps_vol(i,j,k)**2
                        end if
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        if(.not. allocated(S_energy)) allocate(S_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                        S_energy(i,j,k) = 0.
                        if (.not. dom%aniso) then
                            S_energy(i,j,k) =   xmu/2 * (       DXY**2 + DYX**2 &
                                                          +     DXZ**2 + DZX**2 &
                                                          +     DYZ**2 + DZY**2 &
                                                          - 2 * DXY * DYX     &
                                                          - 2 * DXZ * DZX     &
                                                          - 2 * DYZ * DZY )
                        end if
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
                        sig_dev(i,j,k,0:5) = 0.
                        if (.not. dom%aniso) then
                            sig_dev(i,j,k,0) = x2mu * (DXX - eps_vol(i,j,k) * M_1_3)
                            sig_dev(i,j,k,1) = x2mu * (DYY - eps_vol(i,j,k) * M_1_3)
                            sig_dev(i,j,k,2) = x2mu * (DZZ - eps_vol(i,j,k) * M_1_3)
                            sig_dev(i,j,k,3) = xmu * (DXY + DYX)
                            sig_dev(i,j,k,4) = xmu * (DXZ + DZX)
                            sig_dev(i,j,k,5) = xmu * (DYZ + DZY)
                        end if
                    end if
                enddo
            enddo
        enddo
    end subroutine get_solid_dom_var

    subroutine init_material_properties_solid(dom, lnum, i, j, k, density, lambda, mu, kappa, mat)
        use ssubdomains
        type(domain_solid), intent(inout) :: dom
        integer, intent(in) :: lnum
        integer, intent(in) :: i, j, k ! -1 means :
        real(fpp), intent(in) :: density
        real(fpp), intent(in) :: lambda
        real(fpp), intent(in) :: mu
        real(fpp), intent(in) :: kappa
        type (subdomain), intent(in), optional :: mat
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        if (i==-1 .and. j==-1 .and. k==-1) then
            dom%Density_(:,:,:,bnum,ee) = density
            dom%Lambda_ (:,:,:,bnum,ee) = lambda
            dom%Kappa_  (:,:,:,bnum,ee) = kappa
            dom%Mu_     (:,:,:,bnum,ee) = mu
        else
            dom%Density_(i,j,k,bnum,ee) = density
            dom%Lambda_ (i,j,k,bnum,ee) = lambda
            dom%Kappa_  (i,j,k,bnum,ee) = kappa
            dom%Mu_     (i,j,k,bnum,ee) = mu
        end if

        if (present(mat)) then
            if (dom%n_sls>0)  then
                if (dom%aniso) then
                    dom%Q_(:,:,:,bnum,ee) = mat%Qmu
                else
                    dom%Qs_(:,:,:,bnum,ee) = mat%Qmu
                    dom%Qp_(:,:,:,bnum,ee) = mat%Qpression
                endif
            endif
        endif
    end subroutine init_material_properties_solid

    subroutine init_material_tensor_solid(dom, lnum, i, j, k, density, Cij)
        type(domain_solid), intent(inout) :: dom
        integer, intent(in) :: lnum
        integer, intent(in) :: i, j, k
        real(fpp), intent(in) :: density
        real(fpp), dimension(1:6,1:6), intent(in) :: Cij

        integer :: idef, ii, jj
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        idef = 0
        do ii = 1,6
            do jj = ii,6
                dom%Cij_(idef,i,j,k,bnum,ee) = Cij(ii,jj)
                idef = idef + 1
            enddo
        enddo
        dom%Density_(i,j,k,bnum,ee) = density
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

    subroutine forces_int_solid(dom, champs1, bnum)
        use m_calcul_forces
        use m_calcul_forces_atn

        type(domain_solid), intent (INOUT) :: dom
        type(champssolid), intent(inout) :: champs1
        integer, intent(in) :: bnum
        !
        integer :: ngll,i,j,k,i_dir,e,ee,idx
        integer :: n_solid
        logical :: aniso
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fox,Foy,Foz
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: Depla

        n_solid = dom%n_sls
        aniso = dom%aniso
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
                call calcul_forces_iso(dom,bnum,Fox,Foy,Foz,Depla)
            end if
        end if

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
    end subroutine forces_int_solid
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
