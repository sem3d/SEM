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
        implicit none
        type(domain) :: TDomain
        type(domain_solidpml), intent (INOUT) :: dom
        !
        integer nbelem, ngllx, nglly, ngllz
        !

        dom%ngllx = Tdomain%specel(0)%ngllx ! Temporaire: ngll* doit passer sur le domaine a terme
        dom%nglly = Tdomain%specel(0)%nglly ! Temporaire: ngll* doit passer sur le domaine a terme
        dom%ngllz = Tdomain%specel(0)%ngllz ! Temporaire: ngll* doit passer sur le domaine a terme

        nbelem  = dom%nbelem
        if(nbelem == 0) return ! Do not allocate if not needed (save allocation/RAM)
        ngllx   = dom%ngllx
        nglly   = dom%nglly
        ngllz   = dom%ngllz

        allocate(dom%Density_(0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))
        allocate(dom%Lambda_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))
        allocate(dom%Mu_     (0:ngllx-1, 0:nglly-1, 0:ngllz-1,0:nbelem-1))

        allocate (dom%Jacob_  (        0:ngllx-1,0:nglly-1,0:ngllz-1,0:nbelem-1))
        allocate (dom%InvGrad_(0:2,0:2,0:ngllx-1,0:nglly-1,0:ngllz-1,0:nbelem-1))

        allocate(dom%Idom_(0:ngllx-1,0:nglly-1,0:ngllz-1,0:nbelem-1))

        if(Tdomain%TimeD%velocity_scheme)then
            allocate(dom%Diagonal_Stress (0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:nbelem-1))
            allocate(dom%Diagonal_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:nbelem-1))
            allocate(dom%Diagonal_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:nbelem-1))
            allocate(dom%Diagonal_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:nbelem-1))
            allocate(dom%Residual_Stress (0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:nbelem-1))
            allocate(dom%Residual_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:nbelem-1))
            allocate(dom%Residual_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:nbelem-1))
            allocate(dom%Residual_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:nbelem-1))
            dom%Diagonal_Stress  = 0d0
            dom%Diagonal_Stress1 = 0d0
            dom%Diagonal_Stress2 = 0d0
            dom%Diagonal_Stress3 = 0d0
            dom%Residual_Stress  = 0d0
            dom%Residual_Stress1 = 0d0
            dom%Residual_Stress2 = 0d0
            dom%Residual_Stress3 = 0d0
            allocate(dom%PMLDumpSx  (0:ngllx-1,0:nglly-1,0:ngllz-1,0:1,0:nbelem-1))
            allocate(dom%PMLDumpSy  (0:ngllx-1,0:nglly-1,0:ngllz-1,0:1,0:nbelem-1))
            allocate(dom%PMLDumpSz  (0:ngllx-1,0:nglly-1,0:ngllz-1,0:1,0:nbelem-1))
            dom%PMLDumpSx   = 0d0
            dom%PMLDumpSy   = 0d0
            dom%PMLDumpSz   = 0d0
        endif

        ! Allocation et initialisation de champs0 pour les PML solides
        if (dom%ngll /= 0) then
            allocate(dom%champs1%ForcesPML(0:dom%ngll-1,0:2,0:2))
            allocate(dom%champs0%VelocPML (0:dom%ngll-1,0:2,0:2))
            allocate(dom%champs1%VelocPML (0:dom%ngll-1,0:2,0:2))
            allocate(dom%champs0%DumpV    (0:dom%ngll-1,0:1,0:2))
            dom%champs1%ForcesPML = 0d0
            dom%champs0%VelocPML = 0d0
            dom%champs0%DumpV = 0d0

            ! Allocation de MassMat pour les PML solides
            allocate(dom%MassMat(0:dom%ngll-1))
            dom%MassMat = 0d0

            allocate(dom%DumpMass(0:dom%ngll-1,0:2))
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

        if(allocated(dom%Diagonal_Stress )) deallocate(dom%Diagonal_Stress )
        if(allocated(dom%Diagonal_Stress1)) deallocate(dom%Diagonal_Stress1)
        if(allocated(dom%Diagonal_Stress2)) deallocate(dom%Diagonal_Stress2)
        if(allocated(dom%Diagonal_Stress3)) deallocate(dom%Diagonal_Stress3)
        if(allocated(dom%Residual_Stress )) deallocate(dom%Residual_Stress )
        if(allocated(dom%Residual_Stress1)) deallocate(dom%Residual_Stress1)
        if(allocated(dom%Residual_Stress2)) deallocate(dom%Residual_Stress2)
        if(allocated(dom%Residual_Stress3)) deallocate(dom%Residual_Stress3)
        if(allocated(dom%PMLDumpSx  ))      deallocate(dom%PMLDumpSx  )
        if(allocated(dom%PMLDumpSy  ))      deallocate(dom%PMLDumpSy  )
        if(allocated(dom%PMLDumpSz  ))      deallocate(dom%PMLDumpSz  )

        if(allocated(dom%champs1%ForcesPML)) deallocate(dom%champs1%ForcesPML)
        if(allocated(dom%champs0%VelocPML )) deallocate(dom%champs0%VelocPML )
        if(allocated(dom%champs1%VelocPML )) deallocate(dom%champs1%VelocPML )
        if(allocated(dom%champs0%DumpV    )) deallocate(dom%champs0%DumpV    )

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)

        if(allocated(dom%DumpMass)) deallocate(dom%DumpMass)
    end subroutine deallocate_dom_solidpml

    subroutine get_solidpml_dom_var(Tdomain, el, out_variables, &
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
        logical :: flag_gradU
        integer :: nx, ny, nz, i, j, k, ind

        flag_gradU = (out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS) + &
            out_variables(OUT_EPS_VOL) + &
            out_variables(OUT_EPS_DEV) + &
            out_variables(OUT_STRESS_DEV)) /= 0

        nx = el%ngllx
        ny = el%nglly
        nz = el%ngllz

        do k=0,nz-1
            do j=0,ny-1
                do i=0,nx-1
                    ind = Tdomain%spmldom%Idom_(i,j,k,el%lnum)

                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        if(.not. allocated(fieldU)) allocate(fieldU(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldU(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_VITESSE) == 1) then
                        if(.not. allocated(fieldV)) allocate(fieldV(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldV(i,j,k,:) = Tdomain%spmldom%champs0%VelocPml(ind,:,0) + &
                                          Tdomain%spmldom%champs0%VelocPml(ind,:,1) + &
                                          Tdomain%spmldom%champs0%VelocPml(ind,:,2)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldA(i,j,k,:) = Tdomain%spmldom%Massmat(ind) * &
                            ( Tdomain%spmldom%champs1%ForcesPml(ind,:,0) + &
                            Tdomain%spmldom%champs1%ForcesPml(ind,:,1) + &
                            Tdomain%spmldom%champs1%ForcesPml(ind,:,2) )
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:nx-1,0:ny-1,0:nz-1))
                        fieldP(i,j,k) = 0d0
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        if(.not. allocated(eps_vol)) allocate(eps_vol(0:nx-1,0:ny-1,0:nz-1))
                        eps_vol(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        if(.not. allocated(P_energy)) allocate(P_energy(0:nx-1,0:ny-1,0:nz-1))
                        P_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        if(.not. allocated(S_energy)) allocate(S_energy(0:nx-1,0:ny-1,0:nz-1))
                        S_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        if(.not. allocated(eps_dev)) allocate(eps_dev(0:nx-1,0:ny-1,0:nz-1,0:5))
                        eps_dev(i,j,k,:) = 0.
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        if(.not. allocated(sig_dev)) allocate(sig_dev(0:nx-1,0:ny-1,0:nz-1,0:5))
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
        integer :: m1, m2, m3
        integer :: i, j, k, l, ind
        real :: sum_vx, sum_vy, sum_vz, acoeff
        real, dimension(0:2,0:dom%ngllx-1,0:dom%nglly-1,0:dom%ngllz-1)  :: Forces1, Forces2, Forces3

        m1 = dom%ngllx ; m2 = dom%nglly ; m3 = dom%ngllz

        do k = 0,m3-1
            do j = 0,m2-1
                do i=0,m1-1
                    sum_vx = 0d0
                    sum_vy = 0d0
                    sum_vz = 0d0
                    do l = 0,m1-1
                        acoeff = - mat%hprimex(i,l)*mat%GLLwx(l)*mat%GLLwy(j)*mat%GLLwz(k)*dom%Jacob_(l,j,k,lnum)
                        sum_vx = sum_vx + acoeff*dom%InvGrad_(0,0,l,j,k,lnum)*dom%Diagonal_Stress(l,j,k,0,lnum)
                        sum_vx = sum_vx + acoeff*dom%InvGrad_(1,0,l,j,k,lnum)*dom%Residual_Stress(l,j,k,0,lnum)
                        sum_vx = sum_vx + acoeff*dom%InvGrad_(2,0,l,j,k,lnum)*dom%Residual_Stress(l,j,k,1,lnum)

                        sum_vy = sum_vy + acoeff*dom%InvGrad_(0,0,l,j,k,lnum)*dom%Residual_Stress(l,j,k,0,lnum)
                        sum_vy = sum_vy + acoeff*dom%InvGrad_(1,0,l,j,k,lnum)*dom%Diagonal_Stress(l,j,k,1,lnum)
                        sum_vy = sum_vy + acoeff*dom%InvGrad_(2,0,l,j,k,lnum)*dom%Residual_Stress(l,j,k,2,lnum)

                        sum_vz = sum_vz + acoeff*dom%InvGrad_(0,0,l,j,k,lnum)*dom%Residual_Stress(l,j,k,1,lnum)
                        sum_vz = sum_vz + acoeff*dom%InvGrad_(1,0,l,j,k,lnum)*dom%Residual_Stress(l,j,k,2,lnum)
                        sum_vz = sum_vz + acoeff*dom%InvGrad_(2,0,l,j,k,lnum)*dom%Diagonal_Stress(l,j,k,2,lnum)
                    end do
                    Forces1(0,i,j,k) = sum_vx
                    Forces1(1,i,j,k) = sum_vy
                    Forces1(2,i,j,k) = sum_vz
                end do
            end do
        end do

        do k = 0,m3-1
            Forces2(:,:,:,k) = 0d0
            do l = 0,m2-1
                do j = 0,m2-1
                    do i=0,m1-1
                        acoeff = - mat%hprimey(j,l)*mat%GLLwx(i)*mat%GLLwy(l)*mat%GLLwz(k)*dom%Jacob_(i,l,k,lnum)
                        sum_vx = acoeff*(dom%InvGrad_(0,1,i,l,k,lnum)*dom%Diagonal_Stress(i,l,k,0,lnum) + &
                                         dom%InvGrad_(1,1,i,l,k,lnum)*dom%Residual_Stress(i,l,k,0,lnum) + &
                                         dom%InvGrad_(2,1,i,l,k,lnum)*dom%Residual_Stress(i,l,k,1,lnum))

                        sum_vy = acoeff*(dom%InvGrad_(0,1,i,l,k,lnum)*dom%Residual_Stress(i,l,k,0,lnum) + &
                                         dom%InvGrad_(1,1,i,l,k,lnum)*dom%Diagonal_Stress(i,l,k,1,lnum) + &
                                         dom%InvGrad_(2,1,i,l,k,lnum)*dom%Residual_Stress(i,l,k,2,lnum))

                        sum_vz = acoeff*(dom%InvGrad_(0,1,i,l,k,lnum)*dom%Residual_Stress(i,l,k,1,lnum) + &
                                         dom%InvGrad_(1,1,i,l,k,lnum)*dom%Residual_Stress(i,l,k,2,lnum) + &
                                         dom%InvGrad_(2,1,i,l,k,lnum)*dom%Diagonal_Stress(i,l,k,2,lnum))
                        Forces2(0,i,j,k) = Forces2(0,i,j,k) + sum_vx
                        Forces2(1,i,j,k) = Forces2(1,i,j,k) + sum_vy
                        Forces2(2,i,j,k) = Forces2(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do
        ! TODO reorder loops ?
        Forces3 = 0
        do l = 0,m3-1
            do k = 0,m3-1
                do j = 0,m2-1
                    do i=0,m1-1
                        acoeff = - mat%hprimez(k,l)*mat%GLLwx(i)*mat%GLLwy(j)*mat%GLLwz(l)*dom%Jacob_(i,j,l,lnum)
                        sum_vx = acoeff*(dom%InvGrad_(0,2,i,j,l,lnum)*dom%Diagonal_Stress(i,j,l,0,lnum) + &
                                         dom%InvGrad_(1,2,i,j,l,lnum)*dom%Residual_Stress(i,j,l,0,lnum) + &
                                         dom%InvGrad_(2,2,i,j,l,lnum)*dom%Residual_Stress(i,j,l,1,lnum))

                        sum_vy = acoeff*(dom%InvGrad_(0,2,i,j,l,lnum)*dom%Residual_Stress(i,j,l,0,lnum) + &
                                         dom%InvGrad_(1,2,i,j,l,lnum)*dom%Diagonal_Stress(i,j,l,1,lnum) + &
                                         dom%InvGrad_(2,2,i,j,l,lnum)*dom%Residual_Stress(i,j,l,2,lnum))

                        sum_vz = acoeff*(dom%InvGrad_(0,2,i,j,l,lnum)*dom%Residual_Stress(i,j,l,1,lnum) + &
                                         dom%InvGrad_(1,2,i,j,l,lnum)*dom%Residual_Stress(i,j,l,2,lnum) + &
                                         dom%InvGrad_(2,2,i,j,l,lnum)*dom%Diagonal_Stress(i,j,l,2,lnum))
                        Forces3(0,i,j,k) = Forces3(0,i,j,k) + sum_vx
                        Forces3(1,i,j,k) = Forces3(1,i,j,k) + sum_vy
                        Forces3(2,i,j,k) = Forces3(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do

        ! Assemblage
        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
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
        real, dimension(0:dom%ngllx-1, 0:dom%nglly-1, 0:dom%ngllz-1) :: dVx_dx, dVx_dy, dVx_dz
        real, dimension(0:dom%ngllx-1, 0:dom%nglly-1, 0:dom%ngllz-1) :: dVy_dx, dVy_dy, dVy_dz
        real, dimension(0:dom%ngllx-1, 0:dom%nglly-1, 0:dom%ngllz-1) :: dVz_dx, dVz_dy, dVz_dz
        integer :: m1, m2, m3
        integer :: i, j, k, ind, i_dir
        real, dimension (:,:,:,:), allocatable :: Veloc

        m1 = dom%ngllx ; m2 = dom%nglly ; m3 = dom%ngllz

        allocate(Veloc(0:m1-1,0:m2-1,0:m3-1,0:2))
        do i_dir = 0,2
            do k = 0,m3-1
                do j = 0,m2-1
                    do i = 0,m1-1
                        ind = dom%Idom_(i,j,k,lnum)
                        Veloc(i,j,k,i_dir) = champs1%VelocPML(ind,i_dir,0) + &
                            champs1%VelocPML(ind,i_dir,1) + &
                            champs1%VelocPML(ind,i_dir,2)
                    enddo
                enddo
            enddo
        enddo

        ! partial of velocity components with respect to xi,eta,zeta
        call physical_part_deriv(m1,m2,m3,mat%htprimex,mat%hprimey,mat%hprimez,&
             dom%InvGrad_(:,:,:,:,:,lnum), Veloc(:,:,:,0),dVx_dx,dVx_dy,dVx_dz)
        call physical_part_deriv(m1,m2,m3,mat%htprimex,mat%hprimey,mat%hprimez,&
             dom%InvGrad_(:,:,:,:,:,lnum), Veloc(:,:,:,1),dVy_dx,dVy_dy,dVy_dz)
        call physical_part_deriv(m1,m2,m3,mat%htprimex,mat%hprimey,mat%hprimez,&
             dom%InvGrad_(:,:,:,:,:,lnum), Veloc(:,:,:,2),dVz_dx,dVz_dy,dVz_dz)

        deallocate(Veloc)

        ! Stress_xx
        dom%Diagonal_Stress1(:,:,:,0,lnum) = dom%PMLDumpSx(:,:,:,0,lnum)*dom%Diagonal_Stress1(:,:,:,0,lnum) + &
                                             dom%PMLDumpSx(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum)+2*dom%Mu_(:,:,:,lnum))*dVx_dx
        dom%Diagonal_Stress2(:,:,:,0,lnum) = dom%PMLDumpSy(:,:,:,0,lnum)*dom%Diagonal_Stress2(:,:,:,0,lnum) + &
                                             dom%PMLDumpSy(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVy_dy
        dom%Diagonal_Stress3(:,:,:,0,lnum) = dom%PMLDumpSz(:,:,:,0,lnum)*dom%Diagonal_Stress3(:,:,:,0,lnum) + &
                                             dom%PMLDumpSz(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVz_dz

        ! Stress_yy
        dom%Diagonal_Stress1(:,:,:,1,lnum) = dom%PMLDumpSx(:,:,:,0,lnum)*dom%Diagonal_Stress1(:,:,:,1,lnum) + &
                                             dom%PMLDumpSx(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVx_dx
        dom%Diagonal_Stress2(:,:,:,1,lnum) = dom%PMLDumpSy(:,:,:,0,lnum)*dom%Diagonal_Stress2(:,:,:,1,lnum) + &
                                             dom%PMLDumpSy(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum)+2*dom%Mu_(:,:,:,lnum))*dVy_dy
        dom%Diagonal_Stress3(:,:,:,1,lnum) = dom%PMLDumpSz(:,:,:,0,lnum)*dom%Diagonal_Stress3(:,:,:,1,lnum) + &
                                             dom%PMLDumpSz(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVz_dz

        ! Stress_zz
        dom%Diagonal_Stress1(:,:,:,2,lnum) = dom%PMLDumpSx(:,:,:,0,lnum)*dom%Diagonal_Stress1(:,:,:,2,lnum) + &
                                             dom%PMLDumpSx(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVx_dx
        dom%Diagonal_Stress2(:,:,:,2,lnum) = dom%PMLDumpSy(:,:,:,0,lnum)*dom%Diagonal_Stress2(:,:,:,2,lnum) + &
                                             dom%PMLDumpSy(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum))*dVy_dy
        dom%Diagonal_Stress3(:,:,:,2,lnum) = dom%PMLDumpSz(:,:,:,0,lnum)*dom%Diagonal_Stress3(:,:,:,2,lnum) + &
                                             dom%PMLDumpSz(:,:,:,1,lnum)*Dt*(dom%Lambda_(:,:,:,lnum)+2*dom%Mu_(:,:,:,lnum))*dVz_dz

        dom%Diagonal_Stress(:,:,:,:,lnum) = dom%Diagonal_Stress1(:,:,:,:,lnum) + &
                                            dom%Diagonal_Stress2(:,:,:,:,lnum) + &
                                            dom%Diagonal_Stress3(:,:,:,:,lnum)

        ! Stress_xy
        dom%Residual_Stress1 (:,:,:,0,lnum) = dom%PMLDumpSx(:,:,:,0,lnum)*dom%Residual_Stress1 (:,:,:,0,lnum) + &
                                              dom%PMLDumpSx(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVy_dx
        dom%Residual_Stress2 (:,:,:,0,lnum) = dom%PMLDumpSy(:,:,:,0,lnum)*dom%Residual_Stress2 (:,:,:,0,lnum) + &
                                              dom%PMLDumpSy(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVx_dy
        dom%Residual_Stress3 (:,:,:,0,lnum) = dom%PMLDumpSz(:,:,:,0,lnum)*dom%Residual_Stress3 (:,:,:,0,lnum)

        ! Stress_xz
        dom%Residual_Stress1 (:,:,:,1,lnum) = dom%PMLDumpSx(:,:,:,0,lnum)*dom%Residual_Stress1 (:,:,:,1,lnum) + &
                                              dom%PMLDumpSx(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVz_dx
        dom%Residual_Stress2 (:,:,:,1,lnum) = dom%PMLDumpSy(:,:,:,0,lnum)*dom%Residual_Stress2 (:,:,:,1,lnum)
        dom%Residual_Stress3 (:,:,:,1,lnum) = dom%PMLDumpSz(:,:,:,0,lnum)*dom%Residual_Stress3 (:,:,:,1,lnum) + &
                                              dom%PMLDumpSz(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVx_dz

        ! Stress_yz
        dom%Residual_Stress1 (:,:,:,2,lnum) = dom%PMLDumpSx(:,:,:,0,lnum)*dom%Residual_Stress1 (:,:,:,2,lnum)
        dom%Residual_Stress2 (:,:,:,2,lnum) = dom%PMLDumpSy(:,:,:,0,lnum)*dom%Residual_Stress2 (:,:,:,2,lnum) + &
                                              dom%PMLDumpSy(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVz_dy
        dom%Residual_Stress3 (:,:,:,2,lnum) = dom%PMLDumpSz(:,:,:,0,lnum)*dom%Residual_Stress3 (:,:,:,2,lnum) + &
                                              dom%PMLDumpSz(:,:,:,1,lnum)*Dt*(dom%Mu_(:,:,:,lnum))*dVy_dz

        dom%Residual_Stress(:,:,:,:,lnum) = dom%Residual_Stress1(:,:,:,:,lnum) + &
                                            dom%Residual_Stress2(:,:,:,:,lnum) + &
                                            dom%Residual_Stress3(:,:,:,:,lnum)
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
