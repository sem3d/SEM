!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!! Anisotropic-density acoustic domain. Velocity-potential formulation, modelled on
!! dom_fluid.F90, with the anisotropy carried by the DENSITY (Capdeville & Cance 2015):
!!     (1/kappa) phi_tt = d_i( rho^{-1}_ij d_j phi ) ,   v_i = rho^{-1}_ij d_j phi
!! m_IDensTensor = rho^{-1}_ij (0:5 -> 11,22,33,12,13,23) ; m_Lambda = kappa (scalar).

module dom_fluid_aniso
    use sdomain
    use constants
    use champs_fluid_aniso
    use selement
    use ssubdomains
    implicit none
#include "index.h"

contains

    subroutine allocate_champs_fluid_aniso(dom, i)
        type(domain_fluid_aniso), intent(INOUT) :: dom
        integer, intent(in) :: i
        allocate(dom%champs(i)%ForcesFl(0:dom%nglltot))
        allocate(dom%champs(i)%Phi     (0:dom%nglltot))
        allocate(dom%champs(i)%VelPhi  (0:dom%nglltot))

        dom%champs(i)%ForcesFl = 0d0
        dom%champs(i)%Phi      = 0d0
        dom%champs(i)%VelPhi   = 0d0
    end subroutine allocate_champs_fluid_aniso

    subroutine allocate_dom_fluid_aniso(Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: Tdomain
        type(domain_fluid_aniso), intent(INOUT) :: dom
        !
        integer :: nbelem, ngll, nblocks, i

        ngll   = dom%ngll
        nbelem = dom%nbelem
        if (ngll == 0) return

        call init_dombase(dom)

        if (nbelem /= 0) then
            nblocks = dom%nblocks
            ! 6 independent components of rho^{-1}_ij: 11,22,33,12,13,23
            allocate(dom%m_IDensTensor(IND_DIJKE(0:5,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1,0:VCHUNK-1)))
            allocate(dom%m_Lambda(IND_IJKE(0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1,0:VCHUNK-1)))
        end if

        if (dom%nglltot /= 0) then
            do i = 0, Tdomain%TimeD%nsubsteps
                call allocate_champs_fluid_aniso(dom, i)
            end do
        endif

        if (Tdomain%rank==0) write(*,*) "INFO - fluid aniso domain (rank 0, local) : ", &
            dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_fluid_aniso

    subroutine deallocate_dom_fluid_aniso(dom)
        implicit none
        type(domain_fluid_aniso), intent(INOUT) :: dom
        !
        integer :: i
        if (allocated(dom%m_IDensTensor)) deallocate(dom%m_IDensTensor)
        if (allocated(dom%m_Lambda))      deallocate(dom%m_Lambda)

        do i = 0, 1
            if (allocated(dom%champs(i)%ForcesFl)) deallocate(dom%champs(i)%ForcesFl)
            if (allocated(dom%champs(i)%Phi     )) deallocate(dom%champs(i)%Phi     )
            if (allocated(dom%champs(i)%VelPhi  )) deallocate(dom%champs(i)%VelPhi  )
        end do
        call deallocate_dombase(dom)
    end subroutine deallocate_dom_fluid_aniso

    subroutine fluid_aniso_velocity(ngll,hprime,InvGrad,IDensTensor,phi,veloc)
        ! Physical particle velocity in the anisotropic-density fluid:
        !   v_i = rho^{-1}_ij d_j phi
        use deriv3d
        implicit none
        integer, intent(in) :: ngll
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:2), intent(in) :: InvGrad
        real(fpp), dimension(0:5,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: IDensTensor
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: phi
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1,0:2), intent(out) :: Veloc
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: dphi_dx,dphi_dy,dphi_dz
        integer :: i,j,k
        real(fpp) :: B11,B22,B33,B12,B13,B23,gx,gy,gz

        call physical_part_deriv(ngll,hprime,InvGrad,phi,dphi_dx,dphi_dy,dphi_dz)
        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    B11 = IDensTensor(0,i,j,k); B22 = IDensTensor(1,i,j,k); B33 = IDensTensor(2,i,j,k)
                    B12 = IDensTensor(3,i,j,k); B13 = IDensTensor(4,i,j,k); B23 = IDensTensor(5,i,j,k)
                    gx = dphi_dx(i,j,k); gy = dphi_dy(i,j,k); gz = dphi_dz(i,j,k)
                    Veloc(i,j,k,0) = B11*gx + B12*gy + B13*gz
                    Veloc(i,j,k,1) = B12*gx + B22*gy + B23*gz
                    Veloc(i,j,k,2) = B13*gx + B23*gy + B33*gz
                enddo
            enddo
        enddo
    end subroutine fluid_aniso_velocity

    subroutine init_domain_fluid_aniso(Tdomain, dom)
        type(domain), intent(INOUT), target :: Tdomain
        type(domain_fluid_aniso), intent(inout) :: dom
        dom%dt = Tdomain%TimeD%dtmin
    end subroutine init_domain_fluid_aniso

    subroutine init_material_properties_fluid_aniso(dom, lnum, mat, IDens, lambda)
        ! IDens(0:5) = inverse-density tensor rho^{-1}_ij at each GLL point (11,22,33,12,13,23)
        ! lambda     = scalar bulk modulus kappa
        type(domain_fluid_aniso), intent(inout) :: dom
        integer, intent(in) :: lnum
        type(subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:5,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: IDens
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: lambda
        !
        integer :: bnum, ee, m
        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)

        do m = 0, 5
            dom%m_IDensTensor(IND_DIJKE(m,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,bnum,ee)) = IDens(m,:,:,:)
        end do
        dom%m_Lambda(IND_IJKE(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,bnum,ee)) = lambda
    end subroutine init_material_properties_fluid_aniso

    subroutine init_local_mass_fluid_aniso(dom, specel, i, j, k, ind, Whei)
        ! Inertial term weighted by the inverse of the bulk modulus (= regular fluid)
        type(domain_fluid_aniso), intent(INOUT) :: dom
        type(Element), intent(INOUT) :: specel
        integer, intent(in) :: i, j, k, ind
        real(fpp), intent(in) :: Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee   = mod(specel%lnum,VCHUNK)

        specel%MassMat(i,j,k) = Whei * dom%Jacob_(i,j,k,bnum,ee) / dom%Lambda_(i,j,k,bnum,ee)
        dom%MassMat(ind) = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_fluid_aniso

    subroutine forces_int_fluid_aniso(dom, field, bnum)
        use m_calcul_forces_fluid_aniso
        type(domain_fluid_aniso), intent(INOUT) :: dom
        type(champsfluid_aniso), intent(inout) :: field
        integer, intent(in) :: bnum
        !
        integer :: ngll, i, j, k, ee, idx
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fo_Fl, Phi
        real(fpp) :: val

        ngll = dom%ngll

        do k = 0, ngll-1
            do j = 0, ngll-1
                do i = 0, ngll-1
                    do ee = 0, VCHUNK-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        Phi(ee,i,j,k)   = field%Phi(idx)
                        Fo_Fl(ee,i,j,k) = 0d0
                    enddo
                enddo
            enddo
        enddo

        call calcul_forces_fluid_aniso(dom, dom%ngll, bnum, Fo_Fl, Phi)

        do k = 0, ngll-1
            do j = 0, ngll-1
                do i = 0, ngll-1
                    do ee = 0, VCHUNK-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        val = field%ForcesFl(idx)
                        val = val - Fo_Fl(ee,i,j,k)
                        field%ForcesFl(idx) = val
                    enddo
                enddo
            enddo
        enddo
    end subroutine forces_int_fluid_aniso

    subroutine newmark_predictor_fluid_aniso(dom, f0, f1)
        type(domain_fluid_aniso), intent(INOUT) :: dom
        integer, intent(in) :: f0, f1
        dom%champs(f1)%VelPhi   = dom%champs(f0)%VelPhi
        dom%champs(f1)%Phi      = dom%champs(f0)%Phi
        dom%champs(f1)%ForcesFl = 0d0
    end subroutine newmark_predictor_fluid_aniso

    subroutine newmark_corrector_fluid_aniso(dom, dt, f0, f1)
        type(domain_fluid_aniso), intent(INOUT) :: dom
        real(fpp), intent(in) :: dt
        integer, intent(in) :: f0, f1
        !
        integer :: n, indpml

        dom%champs(f0)%ForcesFl = dom%champs(f1)%ForcesFl * dom%MassMat
        dom%champs(f0)%VelPhi   = dom%champs(f0)%VelPhi + dt * dom%champs(f0)%ForcesFl
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs(f0)%VelPhi(indpml) = 0d0
        enddo
        dom%champs(f0)%Phi = dom%champs(f0)%Phi + dt * dom%champs(f0)%VelPhi
    end subroutine newmark_corrector_fluid_aniso

    function fluid_aniso_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        ! Upper bound on phase velocity: vp = sqrt( kappa * lambda_max(rho^{-1}) ).
        ! lambda_max(rho^{-1}) bounded above by the Gershgorin radius of the symmetric tensor.
        type(domain_fluid_aniso), intent(IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed
        real(fpp) :: B11,B22,B33,B12,B13,B23, maxB, kappa
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)

        B11 = dom%m_IDensTensor(IND_DIJKE(0,i,j,k,bnum,ee))
        B22 = dom%m_IDensTensor(IND_DIJKE(1,i,j,k,bnum,ee))
        B33 = dom%m_IDensTensor(IND_DIJKE(2,i,j,k,bnum,ee))
        B12 = dom%m_IDensTensor(IND_DIJKE(3,i,j,k,bnum,ee))
        B13 = dom%m_IDensTensor(IND_DIJKE(4,i,j,k,bnum,ee))
        B23 = dom%m_IDensTensor(IND_DIJKE(5,i,j,k,bnum,ee))
        kappa = dom%Lambda_(i,j,k,bnum,ee)

        maxB = max(abs(B11)+abs(B12)+abs(B13), &
                   abs(B12)+abs(B22)+abs(B23), &
                   abs(B13)+abs(B23)+abs(B33))
        Pspeed = sqrt(maxB * kappa)
    end function fluid_aniso_Pspeed

    subroutine get_fluid_aniso_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, K_energy, D_energy, eps_vol, eps_dev, sig_dev, dUdX)
        use deriv3d
        implicit none
        type(domain_fluid_aniso), intent(inout)   :: dom
        integer, dimension(0:), intent(in)         :: out_variables
        integer, intent(in)                        :: lnum
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldU, fieldV, fieldA
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:8) :: dUdX
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: fieldP, P_energy, K_energy, eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: D_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev, sig_dev
        !
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: phi, vphi, gx, gy, gz
        integer :: ngll, i, j, k, ind, bnum, ee

        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)
        ngll = dom%ngll

        fieldU=0d0; fieldV=0d0; fieldA=0d0; fieldP=0d0
        P_energy=0d0; K_energy=0d0; D_energy=0d0
        eps_vol=0d0; eps_dev=0d0; sig_dev=0d0; dUdX=0d0

        do k=0,ngll-1; do j=0,ngll-1; do i=0,ngll-1
            ind = dom%Idom_(i,j,k,bnum,ee)
            phi(i,j,k)  = dom%champs(0)%Phi(ind)
            vphi(i,j,k) = dom%champs(0)%VelPhi(ind)
        enddo; enddo; enddo

        if (out_variables(OUT_PRESSION) == 1 .or. out_variables(OUT_ENERGYP) == 1) then
            do k=0,ngll-1; do j=0,ngll-1; do i=0,ngll-1
                fieldP(i,j,k) = -vphi(i,j,k)
            enddo; enddo; enddo
        end if

        ! velocity v = rho^{-1} grad(phi) ; acceleration a = rho^{-1} grad(VelPhi)
        if (out_variables(OUT_VITESSE) == 1 .or. out_variables(OUT_ENERGYK) == 1) then
            call fluid_aniso_velocity(ngll, dom%hprime, dom%InvGrad_(:,:,:,:,:,bnum,ee), &
                 dom%m_IDensTensor(IND_DIJKE(0:5,:,:,:,bnum,ee)), phi, fieldV)
        end if
        if (out_variables(OUT_ACCEL) == 1) then
            call fluid_aniso_velocity(ngll, dom%hprime, dom%InvGrad_(:,:,:,:,:,bnum,ee), &
                 dom%m_IDensTensor(IND_DIJKE(0:5,:,:,:,bnum,ee)), vphi, fieldA)
        end if

        if (out_variables(OUT_ENERGYP) == 1) then
            do k=0,ngll-1; do j=0,ngll-1; do i=0,ngll-1
                P_energy(i,j,k) = 0.5d0*fieldP(i,j,k)*fieldP(i,j,k)/dom%Lambda_(i,j,k,bnum,ee)
            enddo; enddo; enddo
        end if

        if (out_variables(OUT_ENERGYK) == 1) then
            ! K_energy = 1/2 (grad phi).(rho^{-1} grad phi) = 1/2 v.(grad phi)
            call physical_part_deriv(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),phi,gx,gy,gz)
            do k=0,ngll-1; do j=0,ngll-1; do i=0,ngll-1
                K_energy(i,j,k) = 0.5d0*(fieldV(i,j,k,0)*gx(i,j,k) + &
                                         fieldV(i,j,k,1)*gy(i,j,k) + &
                                         fieldV(i,j,k,2)*gz(i,j,k))
            enddo; enddo; enddo
        end if
    end subroutine get_fluid_aniso_dom_var

    subroutine get_fluid_aniso_dom_elem_energy(dom, lnum, P_energy, K_energy, D_energy)
        use deriv3d
        type(domain_fluid_aniso), intent(inout) :: dom
        integer, intent(in) :: lnum
        real(fpp), dimension(:,:,:), allocatable, intent(inout) :: P_energy, K_energy
        real(fpp), dimension(:,:,:,:), allocatable, intent(inout) :: D_energy
        !
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: phi, gx, gy, gz
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: vel
        integer :: ngll, i, j, k, ind, bnum, ee
        real(fpp) :: pval

        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)
        ngll = dom%ngll

        if (.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if (.not. allocated(K_energy)) allocate(K_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if (.not. allocated(D_energy)) allocate(D_energy(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        D_energy = 0d0

        do k=0,ngll-1; do j=0,ngll-1; do i=0,ngll-1
            ind = dom%Idom_(i,j,k,bnum,ee)
            phi(i,j,k) = dom%champs(0)%Phi(ind)
        enddo; enddo; enddo

        call fluid_aniso_velocity(ngll, dom%hprime, dom%InvGrad_(:,:,:,:,:,bnum,ee), &
             dom%m_IDensTensor(IND_DIJKE(0:5,:,:,:,bnum,ee)), phi, vel)
        call physical_part_deriv(ngll,dom%hprime,dom%InvGrad_(:,:,:,:,:,bnum,ee),phi,gx,gy,gz)

        do k=0,ngll-1; do j=0,ngll-1; do i=0,ngll-1
            ind  = dom%Idom_(i,j,k,bnum,ee)
            pval = -dom%champs(0)%VelPhi(ind)
            P_energy(i,j,k) = 0.5d0*pval*pval/dom%Lambda_(i,j,k,bnum,ee)
            K_energy(i,j,k) = 0.5d0*(vel(i,j,k,0)*gx(i,j,k) + &
                                     vel(i,j,k,1)*gy(i,j,k) + &
                                     vel(i,j,k,2)*gz(i,j,k))
        enddo; enddo; enddo
    end subroutine get_fluid_aniso_dom_elem_energy

end module dom_fluid_aniso

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
