!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

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
        allocate(dom%champs(i)%ForcesP(0:dom%nglltot))
        allocate(dom%champs(i)%P      (0:dom%nglltot))
        allocate(dom%champs(i)%VelP   (0:dom%nglltot))
        allocate(dom%champs(i)%Vel  (0:2, 0:dom%nglltot))

        dom%champs(i)%ForcesP = 0d0
        dom%champs(i)%P       = 0d0
        dom%champs(i)%VelP    = 0d0
        dom%champs(i)%Vel     = 0d0
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
            ! 6 independent components of Kij: K11,K22,K33,K12,K13,K23
            allocate(dom%m_Kij(IND_DIJKE(0:5,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1,0:VCHUNK-1)))
            allocate(dom%m_Rho(IND_IJKE(0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1,0:VCHUNK-1)))
        end if

        if (dom%nglltot /= 0) then
            do i = 0, Tdomain%TimeD%nsubsteps
                call allocate_champs_fluid_aniso(dom, i)
            end do
        endif

        if (Tdomain%rank==0) write(*,*) "INFO - fluid aniso domain : ", &
            dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_fluid_aniso

    subroutine deallocate_dom_fluid_aniso(dom)
        implicit none
        type(domain_fluid_aniso), intent(INOUT) :: dom
        !
        integer :: i
        if (allocated(dom%m_Kij)) deallocate(dom%m_Kij)
        if (allocated(dom%m_Rho)) deallocate(dom%m_Rho)

        do i = 0, 1
            if (allocated(dom%champs(i)%ForcesP)) deallocate(dom%champs(i)%ForcesP)
            if (allocated(dom%champs(i)%P      )) deallocate(dom%champs(i)%P      )
            if (allocated(dom%champs(i)%VelP   )) deallocate(dom%champs(i)%VelP   )
            if (allocated(dom%champs(i)%Vel    )) deallocate(dom%champs(i)%Vel    )
        end do
        call deallocate_dombase(dom)
    end subroutine deallocate_dom_fluid_aniso

    subroutine init_domain_fluid_aniso(Tdomain, dom)
        type(domain), intent(INOUT), target :: Tdomain
        type(domain_fluid_aniso), intent(inout) :: dom

        dom%dt = Tdomain%TimeD%dtmin
    end subroutine init_domain_fluid_aniso

    subroutine init_material_properties_fluid_aniso(dom, lnum, mat, rho, Kij)
        ! Kij(0:5): K11,K22,K33,K12,K13,K23 at each GLL point
        type(domain_fluid_aniso), intent(inout) :: dom
        integer, intent(in) :: lnum
        type(subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: rho
        real(fpp), intent(in), dimension(0:5,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Kij
        !
        integer :: bnum, ee, m
        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)

        dom%m_Rho(IND_IJKE(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,bnum,ee)) = rho
        do m = 0, 5
            dom%m_Kij(IND_DIJKE(m,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,bnum,ee)) = Kij(m,:,:,:)
        end do
    end subroutine init_material_properties_fluid_aniso

    subroutine init_local_mass_fluid_aniso(dom, specel, i, j, k, ind, Whei)
        ! Mass matrix: M = int rho dOmega  (scalar density)
        type(domain_fluid_aniso), intent(INOUT) :: dom
        type(Element), intent(INOUT) :: specel
        integer, intent(in) :: i, j, k, ind
        real(fpp), intent(in) :: Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee   = mod(specel%lnum,VCHUNK)

        specel%MassMat(i,j,k) = Whei * dom%m_Jacob(IND_IJKE(i,j,k,bnum,ee)) &
                                      * dom%m_Rho(IND_IJKE(i,j,k,bnum,ee))
        dom%MassMat(ind) = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_fluid_aniso

    subroutine forces_int_fluid_aniso(dom, field, bnum)
        use m_calcul_forces_fluid_aniso
        type(domain_fluid_aniso), intent(INOUT) :: dom
        type(champsfluid_aniso), intent(inout) :: field
        integer, intent(in) :: bnum
        !
        integer :: ngll, i, j, k, e, ee, idx
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: FFl, Plocal

        ngll = dom%ngll

        do k = 0, ngll-1
            do j = 0, ngll-1
                do i = 0, ngll-1
                    do ee = 0, VCHUNK-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        Plocal(ee,i,j,k) = field%P(idx)
                        FFl(ee,i,j,k)    = 0d0
                    enddo
                enddo
            enddo
        enddo

        call calcul_forces_fluid_aniso(dom, dom%ngll, bnum, FFl, Plocal)

        do k = 0, ngll-1
            do j = 0, ngll-1
                do i = 0, ngll-1
                    do ee = 0, VCHUNK-1
                        e   = bnum*VCHUNK + ee
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        field%ForcesP(idx) = field%ForcesP(idx) - FFl(ee,i,j,k)
                    enddo
                enddo
            enddo
        enddo
    end subroutine forces_int_fluid_aniso

    subroutine newmark_predictor_fluid_aniso(dom, f0, f1)
        type(domain_fluid_aniso), intent(INOUT) :: dom
        integer, intent(in) :: f0, f1

        dom%champs(f1)%VelP    = dom%champs(f0)%VelP
        dom%champs(f1)%P       = dom%champs(f0)%P
        dom%champs(f1)%ForcesP = 0d0
        dom%champs(f1)%Vel     = dom%champs(f0)%Vel
    end subroutine newmark_predictor_fluid_aniso

    subroutine newmark_corrector_fluid_aniso(dom, dt, f0, f1)
        type(domain_fluid_aniso), intent(INOUT) :: dom
        real(fpp), intent(in) :: dt
        integer, intent(in) :: f0, f1
        !
        integer :: n, inddir, ind
        real(fpp), allocatable :: FVel(:,:)

        ! acc = ForcesP * MassMat_inv  (MassMat stores 1/M after inversion in define_arrays)
        dom%champs(f0)%ForcesP = dom%champs(f1)%ForcesP * dom%MassMat
        dom%champs(f0)%VelP    = dom%champs(f0)%VelP + dt * dom%champs(f0)%ForcesP
        do n = 0, dom%n_dirich-1
            inddir = dom%dirich(n)
            dom%champs(f0)%VelP(inddir) = 0d0
        enddo
        dom%champs(f0)%P = dom%champs(f0)%P + dt * dom%champs(f0)%VelP

        ! Update particle velocity: rho dv/dt = -grad(p)
        allocate(FVel(0:2, 0:dom%nglltot))
        FVel = 0d0
        call compute_vel_forces_fluid_aniso(dom, f0, FVel)
        do ind = 0, dom%nglltot-1
            dom%champs(f0)%Vel(:,ind) = dom%champs(f0)%Vel(:,ind) + dt * dom%MassMat(ind) * FVel(:,ind)
        enddo
        deallocate(FVel)
    end subroutine newmark_corrector_fluid_aniso

    function fluid_aniso_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        ! Upper bound on phase velocity via Gershgorin theorem on Kij/rho
        type(domain_fluid_aniso), intent(IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed
        real(fpp) :: K11,K22,K33,K12,K13,K23, maxK, rho
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)

        K11 = dom%m_Kij(IND_DIJKE(0,i,j,k,bnum,ee))
        K22 = dom%m_Kij(IND_DIJKE(1,i,j,k,bnum,ee))
        K33 = dom%m_Kij(IND_DIJKE(2,i,j,k,bnum,ee))
        K12 = dom%m_Kij(IND_DIJKE(3,i,j,k,bnum,ee))
        K13 = dom%m_Kij(IND_DIJKE(4,i,j,k,bnum,ee))
        K23 = dom%m_Kij(IND_DIJKE(5,i,j,k,bnum,ee))
        rho = dom%m_Rho(IND_IJKE(i,j,k,bnum,ee))

        maxK = max(abs(K11)+abs(K12)+abs(K13), &
                   abs(K12)+abs(K22)+abs(K23), &
                   abs(K13)+abs(K23)+abs(K33))
        Pspeed = sqrt(maxK / rho)
    end function fluid_aniso_Pspeed

    subroutine compute_vel_forces_fluid_aniso(dom, fidx, FVel)
        ! Assembles the weak gradient: FVel(d,ind) -= int w*J * dP/dXd dOmega
        ! Called once per timestep from the corrector to update particle velocity.
        type(domain_fluid_aniso), intent(in) :: dom
        integer, intent(in) :: fidx
        real(fpp), dimension(0:2, 0:dom%nglltot), intent(inout) :: FVel
        !
        integer :: bnum, ee, i, j, k, l, ind, ngll
        real(fpp) :: dP_dxi, dP_deta, dP_dzeta, dP_dX, dP_dY, dP_dZ
        real(fpp) :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3, wJ

        ngll = dom%ngll
        do bnum = 0, dom%nblocks-1
            do ee = 0, VCHUNK-1
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            dP_dxi=0d0; dP_deta=0d0; dP_dzeta=0d0
                            do l = 0,ngll-1
                                dP_dxi   = dP_dxi   + dom%hprime(l,i)*dom%champs(fidx)%P(dom%Idom_(l,j,k,bnum,ee))
                                dP_deta  = dP_deta  + dom%hprime(l,j)*dom%champs(fidx)%P(dom%Idom_(i,l,k,bnum,ee))
                                dP_dzeta = dP_dzeta + dom%hprime(l,k)*dom%champs(fidx)%P(dom%Idom_(i,j,l,bnum,ee))
                            enddo
                            xi1=dom%InvGrad_(0,0,i,j,k,bnum,ee); xi2=dom%InvGrad_(1,0,i,j,k,bnum,ee); xi3=dom%InvGrad_(2,0,i,j,k,bnum,ee)
                            et1=dom%InvGrad_(0,1,i,j,k,bnum,ee); et2=dom%InvGrad_(1,1,i,j,k,bnum,ee); et3=dom%InvGrad_(2,1,i,j,k,bnum,ee)
                            ga1=dom%InvGrad_(0,2,i,j,k,bnum,ee); ga2=dom%InvGrad_(1,2,i,j,k,bnum,ee); ga3=dom%InvGrad_(2,2,i,j,k,bnum,ee)
                            dP_dX = dP_dxi*xi1 + dP_deta*et1 + dP_dzeta*ga1
                            dP_dY = dP_dxi*xi2 + dP_deta*et2 + dP_dzeta*ga2
                            dP_dZ = dP_dxi*xi3 + dP_deta*et3 + dP_dzeta*ga3
                            wJ = dom%Jacob_(i,j,k,bnum,ee)*dom%gllw(i)*dom%gllw(j)*dom%gllw(k)
                            ind = dom%Idom_(i,j,k,bnum,ee)
                            FVel(0,ind) = FVel(0,ind) - wJ*dP_dX
                            FVel(1,ind) = FVel(1,ind) - wJ*dP_dY
                            FVel(2,ind) = FVel(2,ind) - wJ*dP_dZ
                        enddo
                    enddo
                enddo
            enddo
        enddo
    end subroutine compute_vel_forces_fluid_aniso

    subroutine get_fluid_aniso_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, K_energy, D_energy, eps_vol, eps_dev, sig_dev, dUdX)
        type(domain_fluid_aniso), intent(inout) :: dom
        integer, intent(in) :: lnum
        integer, dimension(0:), intent(in) :: out_variables
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldU, fieldV, fieldA
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:8) :: dUdX
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: fieldP
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: P_energy, K_energy, eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: D_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev, sig_dev
        !
        integer :: ngll, i, j, k, ind, bnum, ee
        real(fpp) :: p_val, K_eff, rho_val

        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)
        ngll = dom%ngll

        fieldU=0d0; fieldV=0d0; fieldA=0d0; fieldP=0d0
        P_energy=0d0; K_energy=0d0; D_energy=0d0
        eps_vol=0d0; eps_dev=0d0; sig_dev=0d0; dUdX=0d0

        do k=0,ngll-1; do j=0,ngll-1; do i=0,ngll-1
            ind = dom%Idom_(i,j,k,bnum,ee)
            p_val = dom%champs(0)%P(ind)
            fieldP(i,j,k)   = p_val
            fieldV(i,j,k,0) = dom%champs(0)%Vel(0,ind)
            fieldV(i,j,k,1) = dom%champs(0)%Vel(1,ind)
            fieldV(i,j,k,2) = dom%champs(0)%Vel(2,ind)
            K_eff = (dom%m_Kij(IND_DIJKE(0,i,j,k,bnum,ee)) + &
                     dom%m_Kij(IND_DIJKE(1,i,j,k,bnum,ee)) + &
                     dom%m_Kij(IND_DIJKE(2,i,j,k,bnum,ee))) / 3d0
            rho_val = dom%m_Rho(IND_IJKE(i,j,k,bnum,ee))
            P_energy(i,j,k) = 0.5d0*p_val*p_val/K_eff
            K_energy(i,j,k) = 0.5d0*rho_val*(fieldV(i,j,k,0)**2 + fieldV(i,j,k,1)**2 + fieldV(i,j,k,2)**2)
        enddo; enddo; enddo
    end subroutine get_fluid_aniso_dom_var

    subroutine get_fluid_aniso_dom_elem_energy(dom, lnum, P_energy, K_energy, D_energy)
        type(domain_fluid_aniso), intent(inout) :: dom
        integer, intent(in) :: lnum
        real(fpp), dimension(:,:,:), allocatable, intent(inout) :: P_energy, K_energy
        real(fpp), dimension(:,:,:,:), allocatable, intent(inout) :: D_energy
        !
        integer :: ngll, i, j, k, ind, bnum, ee
        real(fpp) :: p_val, K_eff, rho_val, vx, vy, vz

        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)
        ngll = dom%ngll

        if (.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if (.not. allocated(K_energy)) allocate(K_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if (.not. allocated(D_energy)) allocate(D_energy(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        D_energy = 0d0

        do k=0,ngll-1; do j=0,ngll-1; do i=0,ngll-1
            ind = dom%Idom_(i,j,k,bnum,ee)
            p_val   = dom%champs(0)%P(ind)
            vx      = dom%champs(0)%Vel(0,ind)
            vy      = dom%champs(0)%Vel(1,ind)
            vz      = dom%champs(0)%Vel(2,ind)
            K_eff   = (dom%m_Kij(IND_DIJKE(0,i,j,k,bnum,ee)) + &
                       dom%m_Kij(IND_DIJKE(1,i,j,k,bnum,ee)) + &
                       dom%m_Kij(IND_DIJKE(2,i,j,k,bnum,ee))) / 3d0
            rho_val = dom%m_Rho(IND_IJKE(i,j,k,bnum,ee))
            P_energy(i,j,k) = 0.5d0*p_val*p_val/K_eff
            K_energy(i,j,k) = 0.5d0*rho_val*(vx*vx + vy*vy + vz*vz)
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
