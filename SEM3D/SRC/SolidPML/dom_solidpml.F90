!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_solidpml
    use sdomain
    use constants
    use champs_solidpml
    use selement
    use ssubdomains
    use sdomain
    implicit none
#include "index.h"

contains

    subroutine allocate_dom_solidpml (Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_solidpml), intent (INOUT) :: dom
        !
        integer nbelem, ngll
        !

        nbelem = dom%nbelem
        if(nbelem == 0) return ! Do not allocate if not needed (save allocation/RAM)
        ngll   = dom%ngll

        allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nbelem-1))
        allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nbelem-1))
        allocate(dom%Mu_     (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nbelem-1))

        allocate (dom%Jacob_  (        0:ngll-1,0:ngll-1,0:ngll-1,0:nbelem-1))
        allocate (dom%InvGrad_(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nbelem-1))

        allocate(dom%Idom_(0:ngll-1,0:ngll-1,0:ngll-1,0:nbelem-1))
        dom%m_Idom = 0
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call compute_gll_data(ngll, dom%gllc, dom%gllw, dom%hprime, dom%htprime)

        if(Tdomain%TimeD%velocity_scheme)then
            allocate(dom%Diagonal_Stress_ (0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nbelem-1))
            allocate(dom%Diagonal_Stress1_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nbelem-1))
            allocate(dom%Diagonal_Stress2_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nbelem-1))
            allocate(dom%Diagonal_Stress3_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nbelem-1))
            allocate(dom%Residual_Stress_ (0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nbelem-1))
            allocate(dom%Residual_Stress1_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nbelem-1))
            allocate(dom%Residual_Stress2_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nbelem-1))
            allocate(dom%Residual_Stress3_(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nbelem-1))
            dom%Diagonal_Stress_(:,:,:,:,:)  = 0d0
            dom%Diagonal_Stress1_(:,:,:,:,:) = 0d0
            dom%Diagonal_Stress2_(:,:,:,:,:) = 0d0
            dom%Diagonal_Stress3_(:,:,:,:,:) = 0d0
            dom%Residual_Stress_(:,:,:,:,:)  = 0d0
            dom%Residual_Stress1_(:,:,:,:,:) = 0d0
            dom%Residual_Stress2_(:,:,:,:,:) = 0d0
            dom%Residual_Stress3_(:,:,:,:,:) = 0d0
            allocate(dom%PMLDumpSx_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nbelem-1))
            allocate(dom%PMLDumpSy_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nbelem-1))
            allocate(dom%PMLDumpSz_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nbelem-1))
            dom%PMLDumpSx_(:,:,:,:,:) = 0d0
            dom%PMLDumpSy_(:,:,:,:,:) = 0d0
            dom%PMLDumpSz_(:,:,:,:,:) = 0d0
        endif

        ! Allocation et initialisation de champs0 pour les PML solides
        if (dom%nglltot /= 0) then
            allocate(dom%champs1%ForcesPML(0:dom%nglltot-1,0:2,0:2))
            allocate(dom%champs0%VelocPML (0:dom%nglltot-1,0:2,0:2))
            allocate(dom%champs1%VelocPML (0:dom%nglltot-1,0:2,0:2))
            allocate(dom%champs0%DumpV    (0:dom%nglltot-1,0:1,0:2))
            dom%champs1%ForcesPML = 0d0
            dom%champs0%VelocPML = 0d0
            dom%champs0%DumpV = 0d0

            ! Allocation de MassMat pour les PML solides
            allocate(dom%MassMat(0:dom%nglltot-1))
            dom%MassMat = 0d0

            allocate(dom%DumpMass(0:dom%nglltot-1,0:2))
            dom%DumpMass = 0d0
        endif
    end subroutine allocate_dom_solidpml

    subroutine deallocate_dom_solidpml (dom)
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_Mu     )) deallocate(dom%m_Mu     )

        if(allocated(dom%m_Jacob  )) deallocate(dom%m_Jacob  )
        if(allocated(dom%m_InvGrad)) deallocate(dom%m_InvGrad)

        if(allocated(dom%m_Idom)) deallocate(dom%m_Idom)

        if(allocated(dom%m_Diagonal_Stress )) deallocate(dom%m_Diagonal_Stress )
        if(allocated(dom%m_Diagonal_Stress1)) deallocate(dom%m_Diagonal_Stress1)
        if(allocated(dom%m_Diagonal_Stress2)) deallocate(dom%m_Diagonal_Stress2)
        if(allocated(dom%m_Diagonal_Stress3)) deallocate(dom%m_Diagonal_Stress3)
        if(allocated(dom%m_Residual_Stress )) deallocate(dom%m_Residual_Stress )
        if(allocated(dom%m_Residual_Stress1)) deallocate(dom%m_Residual_Stress1)
        if(allocated(dom%m_Residual_Stress2)) deallocate(dom%m_Residual_Stress2)
        if(allocated(dom%m_Residual_Stress3)) deallocate(dom%m_Residual_Stress3)
        if(allocated(dom%m_PMLDumpSx  ))      deallocate(dom%m_PMLDumpSx  )
        if(allocated(dom%m_PMLDumpSy  ))      deallocate(dom%m_PMLDumpSy  )
        if(allocated(dom%m_PMLDumpSz  ))      deallocate(dom%m_PMLDumpSz  )

        if(allocated(dom%champs1%ForcesPML)) deallocate(dom%champs1%ForcesPML)
        if(allocated(dom%champs0%VelocPML )) deallocate(dom%champs0%VelocPML )
        if(allocated(dom%champs1%VelocPML )) deallocate(dom%champs1%VelocPML )
        if(allocated(dom%champs0%DumpV    )) deallocate(dom%champs0%DumpV    )

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)

        if(allocated(dom%DumpMass)) deallocate(dom%DumpMass)
    end subroutine deallocate_dom_solidpml

    subroutine get_solidpml_dom_var(dom, el, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain_solidpml)                      :: dom
        integer, dimension(0:8)                    :: out_variables
        type(element)                              :: el
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        !
        logical :: flag_gradU
        integer :: ngll, i, j, k, ind

        flag_gradU = (out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS) + &
            out_variables(OUT_EPS_VOL) + &
            out_variables(OUT_EPS_DEV) + &
            out_variables(OUT_STRESS_DEV)) /= 0

        ngll = dom%ngll

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,el%lnum)

                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        if(.not. allocated(fieldU)) allocate(fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldU(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_VITESSE) == 1) then
                        if(.not. allocated(fieldV)) allocate(fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldV(i,j,k,:) = dom%champs0%VelocPml(ind,:,0) + &
                                          dom%champs0%VelocPml(ind,:,1) + &
                                          dom%champs0%VelocPml(ind,:,2)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldA(i,j,k,:) = dom%Massmat(ind) * &
                            ( dom%champs1%ForcesPml(ind,:,0) + &
                            dom%champs1%ForcesPml(ind,:,1) + &
                            dom%champs1%ForcesPml(ind,:,2) )
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
                        fieldP(i,j,k) = 0d0
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        if(.not. allocated(eps_vol)) allocate(eps_vol(0:ngll-1,0:ngll-1,0:ngll-1))
                        eps_vol(i,j,k) = 0.
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
                        eps_dev(i,j,k,:) = 0.
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

        if (i==-1 .and. j==-1 .and. k==-1) then
            dom%Density_(:,:,:,lnum) = density
            dom%Lambda_ (:,:,:,lnum) = lambda
            dom%Mu_     (:,:,:,lnum) = mu
        else
            dom%Density_(i,j,k,lnum) = density
            dom%Lambda_ (i,j,k,lnum) = lambda
            dom%Mu_     (i,j,k,lnum) = mu
        end if
    end subroutine init_material_properties_solidpml

    subroutine init_local_mass_solidpml(dom,specel,i,j,k,ind,Whei)
        type(domain_solidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real Whei

        ! Solid.

        specel%MassMat(i,j,k) = Whei*dom%Density_(i,j,k,specel%lnum)*dom%Jacob_(i,j,k,specel%lnum)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_solidpml

    subroutine forces_int_sol_pml(dom, mat, champs1, lnum)
        type(domain_solidpml), intent(inout) :: dom
        type (subdomain), intent(IN) :: mat
        type(champssolidpml), intent(inout) :: champs1
        integer :: lnum
        !
        integer :: ngll
        integer :: i, j, k, l, ind
        real :: sum_vx, sum_vy, sum_vz, acoeff
        real, dimension(0:2,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)  :: Forces1, Forces2, Forces3

        ngll = dom%ngll

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
                    sum_vx = 0d0
                    sum_vy = 0d0
                    sum_vz = 0d0
                    do l = 0,ngll-1
                        acoeff = - dom%hprime(i,l)*dom%GLLw(l)*dom%GLLw(j)*dom%GLLw(k)*dom%Jacob_(l,j,k,lnum)
                        sum_vx = sum_vx + acoeff*dom%InvGrad_(0,0,l,j,k,lnum)*dom%Diagonal_Stress_(l,j,k,0,lnum)
                        sum_vx = sum_vx + acoeff*dom%InvGrad_(1,0,l,j,k,lnum)*dom%Residual_Stress_(l,j,k,0,lnum)
                        sum_vx = sum_vx + acoeff*dom%InvGrad_(2,0,l,j,k,lnum)*dom%Residual_Stress_(l,j,k,1,lnum)

                        sum_vy = sum_vy + acoeff*dom%InvGrad_(0,0,l,j,k,lnum)*dom%Residual_Stress_(l,j,k,0,lnum)
                        sum_vy = sum_vy + acoeff*dom%InvGrad_(1,0,l,j,k,lnum)*dom%Diagonal_Stress_(l,j,k,1,lnum)
                        sum_vy = sum_vy + acoeff*dom%InvGrad_(2,0,l,j,k,lnum)*dom%Residual_Stress_(l,j,k,2,lnum)

                        sum_vz = sum_vz + acoeff*dom%InvGrad_(0,0,l,j,k,lnum)*dom%Residual_Stress_(l,j,k,1,lnum)
                        sum_vz = sum_vz + acoeff*dom%InvGrad_(1,0,l,j,k,lnum)*dom%Residual_Stress_(l,j,k,2,lnum)
                        sum_vz = sum_vz + acoeff*dom%InvGrad_(2,0,l,j,k,lnum)*dom%Diagonal_Stress_(l,j,k,2,lnum)
                    end do
                    Forces1(0,i,j,k) = sum_vx
                    Forces1(1,i,j,k) = sum_vy
                    Forces1(2,i,j,k) = sum_vz
                end do
            end do
        end do

        do k = 0,ngll-1
            Forces2(:,:,:,k) = 0d0
            do l = 0,ngll-1
                do j = 0,ngll-1
                    do i=0,ngll-1
                        acoeff = - dom%hprime(j,l)*dom%GLLw(i)*dom%GLLw(l)*dom%GLLw(k)*dom%Jacob_(i,l,k,lnum)
                        sum_vx = acoeff*(dom%InvGrad_(0,1,i,l,k,lnum)*dom%Diagonal_Stress_(i,l,k,0,lnum) + &
                                         dom%InvGrad_(1,1,i,l,k,lnum)*dom%Residual_Stress_(i,l,k,0,lnum) + &
                                         dom%InvGrad_(2,1,i,l,k,lnum)*dom%Residual_Stress_(i,l,k,1,lnum))

                        sum_vy = acoeff*(dom%InvGrad_(0,1,i,l,k,lnum)*dom%Residual_Stress_(i,l,k,0,lnum) + &
                                         dom%InvGrad_(1,1,i,l,k,lnum)*dom%Diagonal_Stress_(i,l,k,1,lnum) + &
                                         dom%InvGrad_(2,1,i,l,k,lnum)*dom%Residual_Stress_(i,l,k,2,lnum))

                        sum_vz = acoeff*(dom%InvGrad_(0,1,i,l,k,lnum)*dom%Residual_Stress_(i,l,k,1,lnum) + &
                                         dom%InvGrad_(1,1,i,l,k,lnum)*dom%Residual_Stress_(i,l,k,2,lnum) + &
                                         dom%InvGrad_(2,1,i,l,k,lnum)*dom%Diagonal_Stress_(i,l,k,2,lnum))
                        Forces2(0,i,j,k) = Forces2(0,i,j,k) + sum_vx
                        Forces2(1,i,j,k) = Forces2(1,i,j,k) + sum_vy
                        Forces2(2,i,j,k) = Forces2(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do
        ! TODO reorder loops ?
        Forces3 = 0
        do l = 0,ngll-1
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i=0,ngll-1
                        acoeff = - dom%hprime(k,l)*dom%GLLw(i)*dom%GLLw(j)*dom%GLLw(l)*dom%Jacob_(i,j,l,lnum)
                        sum_vx = acoeff*(dom%InvGrad_(0,2,i,j,l,lnum)*dom%Diagonal_Stress_(i,j,l,0,lnum) + &
                                         dom%InvGrad_(1,2,i,j,l,lnum)*dom%Residual_Stress_(i,j,l,0,lnum) + &
                                         dom%InvGrad_(2,2,i,j,l,lnum)*dom%Residual_Stress_(i,j,l,1,lnum))

                        sum_vy = acoeff*(dom%InvGrad_(0,2,i,j,l,lnum)*dom%Residual_Stress_(i,j,l,0,lnum) + &
                                         dom%InvGrad_(1,2,i,j,l,lnum)*dom%Diagonal_Stress_(i,j,l,1,lnum) + &
                                         dom%InvGrad_(2,2,i,j,l,lnum)*dom%Residual_Stress_(i,j,l,2,lnum))

                        sum_vz = acoeff*(dom%InvGrad_(0,2,i,j,l,lnum)*dom%Residual_Stress_(i,j,l,1,lnum) + &
                                         dom%InvGrad_(1,2,i,j,l,lnum)*dom%Residual_Stress_(i,j,l,2,lnum) + &
                                         dom%InvGrad_(2,2,i,j,l,lnum)*dom%Diagonal_Stress_(i,j,l,2,lnum))
                        Forces3(0,i,j,k) = Forces3(0,i,j,k) + sum_vx
                        Forces3(1,i,j,k) = Forces3(1,i,j,k) + sum_vy
                        Forces3(2,i,j,k) = Forces3(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do

        ! Assemblage
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    ind = dom%Idom_(i,j,k,lnum)
                    champs1%ForcesPML(ind,:,0) = champs1%ForcesPML(ind,:,0) + Forces1(:,i,j,k)
                    champs1%ForcesPML(ind,:,1) = champs1%ForcesPML(ind,:,1) + Forces2(:,i,j,k)
                    champs1%ForcesPML(ind,:,2) = champs1%ForcesPML(ind,:,2) + Forces3(:,i,j,k)
                enddo
            enddo
        enddo
    end subroutine forces_int_sol_pml

    subroutine pred_sol_pml(dom, mat, dt, champs1, lnum)
        implicit none

        type(domain_solidpml), intent(inout) :: dom
        type (subdomain), intent(IN) :: mat
        type(champssolidpml), intent(inout) :: champs1
        real, intent(in) :: dt
        integer :: lnum
        !
        real, dimension(0:dom%ngll-1, 0:dom%ngll-1, 0:dom%ngll-1) :: dVx_dx, dVx_dy, dVx_dz
        real, dimension(0:dom%ngll-1, 0:dom%ngll-1, 0:dom%ngll-1) :: dVy_dx, dVy_dy, dVy_dz
        real, dimension(0:dom%ngll-1, 0:dom%ngll-1, 0:dom%ngll-1) :: dVz_dx, dVz_dy, dVz_dz
        integer :: ngll
        integer :: i, j, k, ind, i_dir
        real, dimension (:,:,:,:), allocatable :: Veloc

        ngll = dom%ngll

        allocate(Veloc(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        do i_dir = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        ind = dom%Idom_(i,j,k,lnum)
                        Veloc(i,j,k,i_dir) = champs1%VelocPML(ind,i_dir,0) + &
                            champs1%VelocPML(ind,i_dir,1) + &
                            champs1%VelocPML(ind,i_dir,2)
                    enddo
                enddo
            enddo
        enddo

        ! partial of velocity components with respect to xi,eta,zeta
        call physical_part_deriv(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,lnum),Veloc(:,:,:,0),dVx_dx,dVx_dy,dVx_dz)
        call physical_part_deriv(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,lnum),Veloc(:,:,:,1),dVy_dx,dVy_dy,dVy_dz)
        call physical_part_deriv(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,lnum),Veloc(:,:,:,2),dVz_dx,dVz_dy,dVz_dz)
        deallocate(Veloc)

        ! Stress_xx
        dom%Diagonal_Stress1_(:,:,:,0,lnum) = dom%PMLDumpSx_(:,:,:,0,lnum)*dom%Diagonal_Stress1_(:,:,:,0,lnum) + &
                                              dom%PMLDumpSx_(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum)+2*dom%Mu_(:,:,:,lnum))*dVx_dx
        dom%Diagonal_Stress2_(:,:,:,0,lnum) = dom%PMLDumpSy_(:,:,:,0,lnum)*dom%Diagonal_Stress2_(:,:,:,0,lnum) + &
                                              dom%PMLDumpSy_(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVy_dy
        dom%Diagonal_Stress3_(:,:,:,0,lnum) = dom%PMLDumpSz_(:,:,:,0,lnum)*dom%Diagonal_Stress3_(:,:,:,0,lnum) + &
                                              dom%PMLDumpSz_(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVz_dz

        ! Stress_yy
        dom%Diagonal_Stress1_(:,:,:,1,lnum) = dom%PMLDumpSx_(:,:,:,0,lnum)*dom%Diagonal_Stress1_(:,:,:,1,lnum) + &
                                              dom%PMLDumpSx_(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVx_dx
        dom%Diagonal_Stress2_(:,:,:,1,lnum) = dom%PMLDumpSy_(:,:,:,0,lnum)*dom%Diagonal_Stress2_(:,:,:,1,lnum) + &
                                              dom%PMLDumpSy_(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum)+2*dom%Mu_(:,:,:,lnum))*dVy_dy
        dom%Diagonal_Stress3_(:,:,:,1,lnum) = dom%PMLDumpSz_(:,:,:,0,lnum)*dom%Diagonal_Stress3_(:,:,:,1,lnum) + &
                                              dom%PMLDumpSz_(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVz_dz

        ! Stress_zz
        dom%Diagonal_Stress1_(:,:,:,2,lnum) = dom%PMLDumpSx_(:,:,:,0,lnum)*dom%Diagonal_Stress1_(:,:,:,2,lnum) + &
                                              dom%PMLDumpSx_(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVx_dx
        dom%Diagonal_Stress2_(:,:,:,2,lnum) = dom%PMLDumpSy_(:,:,:,0,lnum)*dom%Diagonal_Stress2_(:,:,:,2,lnum) + &
                                              dom%PMLDumpSy_(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVy_dy
        dom%Diagonal_Stress3_(:,:,:,2,lnum) = dom%PMLDumpSz_(:,:,:,0,lnum)*dom%Diagonal_Stress3_(:,:,:,2,lnum) + &
                                              dom%PMLDumpSz_(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum)+2*dom%Mu_(:,:,:,lnum))*dVz_dz

        dom%Diagonal_Stress_(:,:,:,:,lnum) = dom%Diagonal_Stress1_(:,:,:,:,lnum) + &
                                             dom%Diagonal_Stress2_(:,:,:,:,lnum) + &
                                             dom%Diagonal_Stress3_(:,:,:,:,lnum)

        ! Stress_xy
        dom%Residual_Stress1_(:,:,:,0,lnum) = dom%PMLDumpSx_(:,:,:,0,lnum)*dom%Residual_Stress1_(:,:,:,0,lnum) + &
                                              dom%PMLDumpSx_(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVy_dx
        dom%Residual_Stress2_(:,:,:,0,lnum) = dom%PMLDumpSy_(:,:,:,0,lnum)*dom%Residual_Stress2_(:,:,:,0,lnum) + &
                                              dom%PMLDumpSy_(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVx_dy
        dom%Residual_Stress3_(:,:,:,0,lnum) = dom%PMLDumpSz_(:,:,:,0,lnum)*dom%Residual_Stress3_(:,:,:,0,lnum)

        ! Stress_xz
        dom%Residual_Stress1_ (:,:,:,1,lnum) = dom%PMLDumpSx_(:,:,:,0,lnum)*dom%Residual_Stress1_(:,:,:,1,lnum) + &
                                               dom%PMLDumpSx_(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVz_dx
        dom%Residual_Stress2_ (:,:,:,1,lnum) = dom%PMLDumpSy_(:,:,:,0,lnum)*dom%Residual_Stress2_(:,:,:,1,lnum)
        dom%Residual_Stress3_ (:,:,:,1,lnum) = dom%PMLDumpSz_(:,:,:,0,lnum)*dom%Residual_Stress3_(:,:,:,1,lnum) + &
                                               dom%PMLDumpSz_(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVx_dz

        ! Stress_yz
        dom%Residual_Stress1_ (:,:,:,2,lnum) = dom%PMLDumpSx_(:,:,:,0,lnum)*dom%Residual_Stress1_(:,:,:,2,lnum)
        dom%Residual_Stress2_ (:,:,:,2,lnum) = dom%PMLDumpSy_(:,:,:,0,lnum)*dom%Residual_Stress2_(:,:,:,2,lnum) + &
                                               dom%PMLDumpSy_(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVz_dy
        dom%Residual_Stress3_ (:,:,:,2,lnum) = dom%PMLDumpSz_(:,:,:,0,lnum)*dom%Residual_Stress3_(:,:,:,2,lnum) + &
                                               dom%PMLDumpSz_(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVy_dz

        dom%Residual_Stress_(:,:,:,:,lnum) = dom%Residual_Stress1_(:,:,:,:,lnum) + &
                                             dom%Residual_Stress2_(:,:,:,:,lnum) + &
                                             dom%Residual_Stress3_(:,:,:,:,lnum)
    end subroutine pred_sol_pml

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
