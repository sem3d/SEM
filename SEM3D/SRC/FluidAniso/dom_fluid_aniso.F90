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

        dom%champs(i)%ForcesP = 0d0
        dom%champs(i)%P       = 0d0
        dom%champs(i)%VelP    = 0d0
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
    end subroutine newmark_predictor_fluid_aniso

    subroutine newmark_corrector_fluid_aniso(dom, dt, f0, f1)
        type(domain_fluid_aniso), intent(INOUT) :: dom
        real(fpp), intent(in) :: dt
        integer, intent(in) :: f0, f1
        !
        integer :: n, inddir

        ! acc = ForcesP * MassMat_inv  (MassMat stores 1/M after inversion in define_arrays)
        dom%champs(f0)%ForcesP = dom%champs(f1)%ForcesP * dom%MassMat
        dom%champs(f0)%VelP    = dom%champs(f0)%VelP + dt * dom%champs(f0)%ForcesP
        do n = 0, dom%n_dirich-1
            inddir = dom%dirich(n)
            dom%champs(f0)%VelP(inddir) = 0d0
        enddo
        dom%champs(f0)%P = dom%champs(f0)%P + dt * dom%champs(f0)%VelP
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
