!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_solid
    use sdomain
    use constants
    use champs_solid
    use selement
    use ssubdomains
    implicit none
#include "index.h"

contains

    subroutine allocate_dom_solid (Tdomain, dom)
        implicit none
        type(domain) :: TDomain
        type(domain_solid), intent (INOUT) :: dom
        !
        integer nbelem, ngllx, nglly, ngllz, n_solid
        logical aniso
        !
        dom%ngllx = Tdomain%specel(0)%ngllx ! Temporaire: ngll* doit passer sur le domaine a terme
        dom%nglly = Tdomain%specel(0)%nglly ! Temporaire: ngll* doit passer sur le domaine a terme
        dom%ngllz = Tdomain%specel(0)%ngllz ! Temporaire: ngll* doit passer sur le domaine a terme

        nbelem  = dom%nbelem
        if(nbelem == 0) return ! Do not allocate if not needed (save allocation/RAM)
        ngllx   = dom%ngllx
        nglly   = dom%nglly
        ngllz   = dom%ngllz
        aniso   = Tdomain%aniso
        n_solid = Tdomain%n_sls

        allocate(dom%Density_(0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))
        allocate(dom%Lambda_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))
        allocate(dom%Mu_     (0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))
        allocate(dom%Kappa_  (0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))

        allocate (dom%Jacob_  (        0:ngllx-1,0:nglly-1,0:ngllz-1,0:nbelem-1))
        allocate (dom%InvGrad_(0:2,0:2,0:ngllx-1,0:nglly-1,0:ngllz-1,0:nbelem-1))

        allocate(dom%Idom_(0:ngllx-1,0:nglly-1,0:ngllz-1,0:nbelem-1))

        if (aniso) then
            allocate (dom%Cij_ (0:20, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
        endif
        if (n_solid>0) then
            if (aniso) then
                allocate (dom%Q (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            else
                allocate (dom%Qs (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
                allocate (dom%Qp (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
                allocate (dom%onemPbeta (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
                allocate (dom%epsilonvol_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
                dom%epsilonvol_ = 0
                allocate (dom%factor_common_P (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
                allocate (dom%alphaval_P (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
                allocate (dom%betaval_P (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
                allocate (dom%gammaval_P (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
                allocate (dom%R_vol_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
                dom%R_vol_ = 0
            endif
            allocate (dom%onemSbeta (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%epsilondev_xx_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%epsilondev_yy_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%epsilondev_xy_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%epsilondev_xz_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%epsilondev_yz_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            dom%epsilondev_xx_ = 0
            dom%epsilondev_yy_ = 0
            dom%epsilondev_xy_ = 0
            dom%epsilondev_xz_ = 0
            dom%epsilondev_yz_ = 0
            allocate (dom%factor_common_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%alphaval_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%betaval_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%gammaval_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%R_xx_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%R_yy_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%R_xy_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%R_xz_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            allocate (dom%R_yz_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:nbelem-1))
            dom%R_xx_ = 0
            dom%R_yy_ = 0
            dom%R_xy_ = 0
            dom%R_xz_ = 0
            dom%R_yz_ = 0
        endif ! n_solid

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

        if(allocated(dom%m_Cij           )) deallocate (dom%m_Cij           )
        if(allocated(dom%Q               )) deallocate (dom%Q               )
        if(allocated(dom%Qs              )) deallocate (dom%Qs              )
        if(allocated(dom%Qp              )) deallocate (dom%Qp              )
        if(allocated(dom%onemPbeta       )) deallocate (dom%onemPbeta       )
        if(allocated(dom%epsilonvol_     )) deallocate (dom%epsilonvol_     )
        if(allocated(dom%factor_common_P )) deallocate (dom%factor_common_P )
        if(allocated(dom%alphaval_P      )) deallocate (dom%alphaval_P      )
        if(allocated(dom%betaval_P       )) deallocate (dom%betaval_P       )
        if(allocated(dom%gammaval_P      )) deallocate (dom%gammaval_P      )
        if(allocated(dom%R_vol_          )) deallocate (dom%R_vol_          )
        if(allocated(dom%onemSbeta       )) deallocate (dom%onemSbeta       )
        if(allocated(dom%epsilondev_xx_  )) deallocate (dom%epsilondev_xx_  )
        if(allocated(dom%epsilondev_yy_  )) deallocate (dom%epsilondev_yy_  )
        if(allocated(dom%epsilondev_xy_  )) deallocate (dom%epsilondev_xy_  )
        if(allocated(dom%epsilondev_xz_  )) deallocate (dom%epsilondev_xz_  )
        if(allocated(dom%epsilondev_yz_  )) deallocate (dom%epsilondev_yz_  )
        if(allocated(dom%factor_common_3 )) deallocate (dom%factor_common_3 )
        if(allocated(dom%alphaval_3      )) deallocate (dom%alphaval_3      )
        if(allocated(dom%betaval_3       )) deallocate (dom%betaval_3       )
        if(allocated(dom%gammaval_3      )) deallocate (dom%gammaval_3      )
        if(allocated(dom%R_xx_           )) deallocate (dom%R_xx_           )
        if(allocated(dom%R_yy_           )) deallocate (dom%R_yy_           )
        if(allocated(dom%R_xy_           )) deallocate (dom%R_xy_           )
        if(allocated(dom%R_xz_           )) deallocate (dom%R_xz_           )
        if(allocated(dom%R_yz_           )) deallocate (dom%R_yz_           )

        if(allocated(dom%champs0%Forces)) deallocate(dom%champs0%Forces)
        if(allocated(dom%champs0%Depla )) deallocate(dom%champs0%Depla )
        if(allocated(dom%champs0%Veloc )) deallocate(dom%champs0%Veloc )
        if(allocated(dom%champs1%Forces)) deallocate(dom%champs1%Forces)
        if(allocated(dom%champs1%Depla )) deallocate(dom%champs1%Depla )
        if(allocated(dom%champs1%Veloc )) deallocate(dom%champs1%Veloc )

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)
    end subroutine deallocate_dom_solid

    subroutine get_solid_dom_var(Tdomain, el, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain)                               :: TDomain
        integer, dimension(0:8)                    :: out_variables
        type(element)                              :: el
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        !
        logical                                  :: flag_gradU
        integer                                  :: nx, ny, nz, i, j, k, ind, mat
        real(fpp), dimension(:,:,:), allocatable :: DXX, DXY, DXZ
        real(fpp), dimension(:,:,:), allocatable :: DYX, DYY, DYZ
        real(fpp), dimension(:,:,:), allocatable :: DZX, DZY, DZZ
        real(fpp), dimension (:,:), allocatable  :: htprimex, hprimey, hprimez
        real(fpp)                                :: xmu, xlambda, xkappa
        real(fpp)                                :: x2mu, xlambda2mu
        real(fpp)                                :: onemSbeta, onemPbeta

        flag_gradU = (out_variables(OUT_PRESSION)    + &
                      out_variables(OUT_ENERGYP)     + &
                      out_variables(OUT_ENERGYS)     + &
                      out_variables(OUT_EPS_VOL)     + &
                      out_variables(OUT_EPS_DEV)     + &
                      out_variables(OUT_STRESS_DEV)) /= 0

        nx = el%ngllx
        ny = el%nglly
        nz = el%ngllz
        mat = el%mat_index

        if (flag_gradU) then
            allocate(DXX(0:nx-1,0:ny-1,0:nz-1))
            allocate(DXY(0:nx-1,0:ny-1,0:nz-1))
            allocate(DXZ(0:nx-1,0:ny-1,0:nz-1))
            allocate(DYX(0:nx-1,0:ny-1,0:nz-1))
            allocate(DYY(0:nx-1,0:ny-1,0:nz-1))
            allocate(DYZ(0:nx-1,0:ny-1,0:nz-1))
            allocate(DZX(0:nx-1,0:ny-1,0:nz-1))
            allocate(DZY(0:nx-1,0:ny-1,0:nz-1))
            allocate(DZZ(0:nx-1,0:ny-1,0:nz-1))
            allocate(hTprimex(0:nx-1,0:nx-1))
            allocate(hprimey(0:ny-1,0:ny-1))
            allocate(hprimez(0:nz-1,0:nz-1))
            hTprimex=Tdomain%sSubDomain(mat)%hTprimex
            hprimey=Tdomain%sSubDomain(mat)%hprimey
            hprimez=Tdomain%sSubDomain(mat)%hprimez
        end if

        ! First, get displacement.

        do k=0,nz-1
            do j=0,ny-1
                do i=0,nx-1
                    ind = Tdomain%sdom%Idom_(i,j,k,el%lnum)

                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        if(.not. allocated(fieldU)) allocate(fieldU(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldU(i,j,k,:) = Tdomain%sdom%champs0%Depla(ind,:)
                    end if

                enddo
            enddo
        enddo

        ! Then, compute gradU with displacement if needed.

        if (flag_gradU) then
            call physical_part_deriv(nx,ny,nz,htprimex,hprimey,hprimez,&
                 Tdomain%sdom%InvGrad_(:,:,:,:,:,el%lnum),fieldU(:,:,:,0),DXX,DYX,DZX)
            call physical_part_deriv(nx,ny,nz,htprimex,hprimey,hprimez,&
                 Tdomain%sdom%InvGrad_(:,:,:,:,:,el%lnum),fieldU(:,:,:,1),DXY,DYY,DZY)
            call physical_part_deriv(nx,ny,nz,htprimex,hprimey,hprimez,&
                 Tdomain%sdom%InvGrad_(:,:,:,:,:,el%lnum),fieldU(:,:,:,2),DXZ,DYZ,DZZ)
        end if

        ! Then, get other variables.

        do k=0,nz-1
            do j=0,ny-1
                do i=0,nx-1
                    ind = Tdomain%sdom%Idom_(i,j,k,el%lnum)

                    if (out_variables(OUT_VITESSE) == 1) then
                        if(.not. allocated(fieldV)) allocate(fieldV(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldV(i,j,k,:) = Tdomain%sdom%champs0%Veloc(ind,:)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldA(i,j,k,:) = Tdomain%sdom%champs0%Forces(ind,:)
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:nx-1,0:ny-1,0:nz-1))
                        fieldP(i,j,k) = -(Tdomain%sdom%Lambda_(i,j,k,el%lnum)&
                                          +2d0/3d0*Tdomain%sdom%Mu_(i,j,k,el%lnum))&
                                        *(DXX(i,j,k)+DYY(i,j,k)+DZZ(i,j,k))
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        if(.not. allocated(eps_vol)) allocate(eps_vol(0:nx-1,0:ny-1,0:nz-1))
                        eps_vol(i,j,k) = DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k)
                    end if

                    if (.not. Tdomain%aniso) then
                        xmu     = Tdomain%sdom%Mu_    (i,j,k,el%lnum)
                        xlambda = Tdomain%sdom%Lambda_(i,j,k,el%lnum)
                        xkappa  = Tdomain%sdom%Kappa_ (i,j,k,el%lnum)

                        if (Tdomain%n_sls>0) then
                            onemSbeta = Tdomain%sdom%onemSbeta(i,j,k,el%lnum)
                            onemPbeta = Tdomain%sdom%onemPbeta(i,j,k,el%lnum)
                            xmu    = xmu * onemSbeta
                            xkappa = xkappa * onemPbeta
                        endif
                        x2mu       = 2. * xmu
                        xlambda2mu = xlambda + x2mu
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        if(.not. allocated(P_energy)) allocate(P_energy(0:nx-1,0:ny-1,0:nz-1))
                        P_energy(i,j,k) = 0.
                        if (.not. Tdomain%aniso) then
                            P_energy(i,j,k) = .5 * xlambda2mu * eps_vol(i,j,k)**2
                        end if
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        if(.not. allocated(S_energy)) allocate(S_energy(0:nx-1,0:ny-1,0:nz-1))
                        S_energy(i,j,k) = 0.
                        if (.not. Tdomain%aniso) then
                            S_energy(i,j,k) =   xmu/2 * ( DXY(i,j,k)**2 + DYX(i,j,k)**2 &
                                              +   DXZ(i,j,k)**2 + DZX(i,j,k)**2 &
                                              +   DYZ(i,j,k)**2 + DZY(i,j,k)**2 &
                                              - 2 * DXY(i,j,k) * DYX(i,j,k)     &
                                              - 2 * DXZ(i,j,k) * DZX(i,j,k)     &
                                              - 2 * DYZ(i,j,k) * DZY(i,j,k))
                        end if
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        if(.not. allocated(eps_dev)) allocate(eps_dev(0:nx-1,0:ny-1,0:nz-1,0:5))
                        eps_dev(i,j,k,0) = DXX(i,j,k) - eps_vol(i,j,k) / 3
                        eps_dev(i,j,k,1) = DYY(i,j,k) - eps_vol(i,j,k) / 3
                        eps_dev(i,j,k,2) = DZZ(i,j,k) - eps_vol(i,j,k) / 3
                        eps_dev(i,j,k,3) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
                        eps_dev(i,j,k,4) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
                        eps_dev(i,j,k,5) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        if(.not. allocated(sig_dev)) allocate(sig_dev(0:nx-1,0:ny-1,0:nz-1,0:5))
                        sig_dev(i,j,k,0:5) = 0.
                        if (.not. Tdomain%aniso) then
                            sig_dev(i,j,k,0) = x2mu * (DXX(i,j,k) - eps_vol(i,j,k) * M_1_3)
                            sig_dev(i,j,k,1) = x2mu * (DYY(i,j,k) - eps_vol(i,j,k) * M_1_3)
                            sig_dev(i,j,k,2) = x2mu * (DZZ(i,j,k) - eps_vol(i,j,k) * M_1_3)
                            sig_dev(i,j,k,3) = xmu * (DXY(i,j,k) + DYX(i,j,k))
                            sig_dev(i,j,k,4) = xmu * (DXZ(i,j,k) + DZX(i,j,k))
                            sig_dev(i,j,k,5) = xmu * (DYZ(i,j,k) + DZY(i,j,k))
                        end if
                    end if
                enddo
            enddo
        enddo

        if (allocated(DXX))      deallocate(DXX)
        if (allocated(DXY))      deallocate(DXY)
        if (allocated(DXZ))      deallocate(DXZ)
        if (allocated(DYX))      deallocate(DYX)
        if (allocated(DYY))      deallocate(DYY)
        if (allocated(DYZ))      deallocate(DYZ)
        if (allocated(DZX))      deallocate(DZX)
        if (allocated(DZY))      deallocate(DZY)
        if (allocated(DZZ))      deallocate(DZZ)
        if (allocated(hTprimex)) deallocate(hTprimex)
        if (allocated(hprimey))  deallocate(hprimey)
        if (allocated(hprimez))  deallocate(hprimez)
    end subroutine get_solid_dom_var

    subroutine init_material_properties_solid(dom, lnum, i, j, k, density, lambda, mu, kappa, Tdomain, mat)
        type(domain_solid), intent(inout) :: dom
        integer, intent(in) :: lnum
        integer, intent(in) :: i, j, k ! -1 means :
        real(fpp), intent(in) :: density
        real(fpp), intent(in) :: lambda
        real(fpp), intent(in) :: mu
        real(fpp), intent(in) :: kappa
        type(domain), intent(in), optional :: Tdomain
        type (subdomain), intent(in), optional :: mat

        if (i==-1 .and. j==-1 .and. k==-1) then
            dom%Density_(:,:,:,lnum) = density
            dom%Lambda_ (:,:,:,lnum) = lambda
            dom%Kappa_  (:,:,:,lnum) = kappa
            dom%Mu_     (:,:,:,lnum) = mu
        else
            dom%Density_(i,j,k,lnum) = density
            dom%Lambda_ (i,j,k,lnum) = lambda
            dom%Kappa_  (i,j,k,lnum) = kappa
            dom%Mu_     (i,j,k,lnum) = mu
        end if

        if (present(Tdomain) .and. present(mat)) then
            if (Tdomain%n_sls>0)  then
                if (Tdomain%aniso) then
                    dom%Q = mat%Qmu
                else
                    dom%Qs = mat%Qmu
                    dom%Qp = mat%Qpression
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

        idef = 0
        do ii = 1,6
            do jj = ii,6
                dom%Cij_(idef,i,j,k,lnum) = Cij(ii,jj)
                idef = idef + 1
            enddo
        enddo
        dom%Density_(i,j,k,lnum) = density
    end subroutine init_material_tensor_solid

    subroutine init_local_mass_solid(dom,specel,i,j,k,ind,Whei)
        type(domain_solid), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real Whei

        ! Solid.

        specel%MassMat(i,j,k) = Whei*dom%Density_(i,j,k,specel%lnum)*dom%Jacob_(i,j,k,specel%lnum)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_solid

    subroutine forces_int_solid(dom, mat, htprimex, hprimey, htprimey, hprimez, htprimez,  &
               n_solid, aniso, champs1, lnum)

        use attenuation_solid

        type(domain_solid), intent (INOUT) :: dom
        type (subdomain), intent(IN) :: mat
        real, dimension (0:dom%ngllx-1, 0:dom%ngllx-1), intent (IN) :: htprimex
        real, dimension (0:dom%nglly-1, 0:dom%nglly-1), intent (IN) :: hprimey, htprimey
        real, dimension (0:dom%ngllz-1, 0:dom%ngllz-1), intent (IN) :: hprimez, hTprimez
        integer, intent(IN) :: n_solid
        logical, intent(IN) :: aniso
        type(champssolid), intent(inout) :: champs1
        integer :: lnum

        integer :: m1,m2,m3, i,j,k,i_dir
        real :: epsilon_trace_over_3
        real, dimension (0:dom%ngllx-1, 0:dom%nglly-1, 0:dom%ngllz-1) ::  DXX,DXY,DXZ, &
            DYX,DYY,DYZ,DZX,DZY,DZZ,Fox,Foy,Foz

        real, dimension(:,:,:), allocatable :: epsilondev_xx_loc, epsilondev_yy_loc, &
            epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
        real, dimension(:,:,:), allocatable :: epsilonvol_loc
        real, dimension(0:dom%ngllx-1, 0:dom%nglly-1, 0:dom%ngllz-1,0:2) :: Depla

        m1 = dom%ngllx;   m2 = dom%nglly;   m3 = dom%ngllz

        do i_dir = 0,2
            do k = 0,m3-1
                do j = 0,m2-1
                    do i = 0,m1-1
                        Depla(i,j,k,i_dir) = champs1%Depla(dom%Idom_(i,j,k,lnum),i_dir)
                    enddo
                enddo
            enddo
        enddo

        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,&
             dom%InvGrad_(:,:,:,:,:,lnum),Depla(:,:,:,0),dxx,dyx,dzx)
        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,&
             dom%InvGrad_(:,:,:,:,:,lnum),Depla(:,:,:,1),dxy,dyy,dzy)
        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,&
             dom%InvGrad_(:,:,:,:,:,lnum),Depla(:,:,:,2),dxz,dyz,dzz)

        if (n_solid>0) then
            if (aniso) then
            else
                allocate (epsilonvol_loc(0:m1-1,0:m2-1,0:m3-1))
            endif
            allocate (epsilondev_xx_loc(0:m1-1,0:m2-1,0:m3-1))
            allocate (epsilondev_yy_loc(0:m1-1,0:m2-1,0:m3-1))
            allocate (epsilondev_xy_loc(0:m1-1,0:m2-1,0:m3-1))
            allocate (epsilondev_xz_loc(0:m1-1,0:m2-1,0:m3-1))
            allocate (epsilondev_yz_loc(0:m1-1,0:m2-1,0:m3-1))
            do i = 0,m1-1
                do j = 0,m2-1
                    do k = 0,m3-1
                        epsilon_trace_over_3 = 0.333333333333333333333333333333d0 * (DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k))
                        if (aniso) then
                        else
                            epsilonvol_loc(i,j,k) = DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k)
                        endif
                        epsilondev_xx_loc(i,j,k) = DXX(i,j,k) - epsilon_trace_over_3
                        epsilondev_yy_loc(i,j,k) = DYY(i,j,k) - epsilon_trace_over_3
                        epsilondev_xy_loc(i,j,k) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
                        epsilondev_xz_loc(i,j,k) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
                        epsilondev_yz_loc(i,j,k) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
                    enddo
                enddo
            enddo
        endif

        if (aniso) then
            if (n_solid>0) then
                call calcul_forces_aniso_att(dom,lnum,Fox,Foy,Foz,htprimex,htprimey,htprimez,&
                     mat%GLLwx,mat%GLLwy,mat%GLLwz,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,m1,m2,m3,n_solid)
                call attenuation_aniso_update(dom,lnum,epsilondev_xx_loc,epsilondev_yy_loc, &
                     epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc,m1,m2,m3,n_solid)
                deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
            else
                call calcul_forces_aniso(dom,lnum,Fox,Foy,Foz,htprimex,htprimey,htprimez, &
                     mat%GLLwx, mat%GLLwy, mat%GLLwz,DXX,DXY,DXZ,DYX,DYY,DYZ, DZX,DZY,DZZ,m1,m2,m3)
            endif
        else
            if (n_solid>0) then
                call calcul_forces_att(dom,lnum,Fox,Foy,Foz,htprimex,htprimey,htprimez, &
                     mat%GLLwx, mat%GLLwy, mat%GLLwz,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,m1,m2,m3, n_solid)
                call attenuation_update(dom,lnum,epsilondev_xx_loc,epsilondev_yy_loc, &
                     epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc,m1,m2,m3,n_solid,epsilonvol_loc)
                deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
                deallocate(epsilonvol_loc)
            else
                call calcul_forces(dom,lnum,Fox,Foy,Foz,htprimex,htprimey,htprimez, &
                     mat%GLLwx,mat%GLLwy,mat%GLLwz,DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ,m1,m2,m3)
            endif
        endif

        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    champs1%Forces(dom%Idom_(i,j,k,lnum),0) = champs1%Forces(dom%Idom_(i,j,k,lnum),0)-Fox(i,j,k)
                    champs1%Forces(dom%Idom_(i,j,k,lnum),1) = champs1%Forces(dom%Idom_(i,j,k,lnum),1)-Foy(i,j,k)
                    champs1%Forces(dom%Idom_(i,j,k,lnum),2) = champs1%Forces(dom%Idom_(i,j,k,lnum),2)-Foz(i,j,k)
                enddo
            enddo
        enddo

        return
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
