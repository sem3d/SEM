!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!----------------------------------------------------------------
!----------------------------------------------------------------

! For CPML, the reference used is :
! Ref3. Application to ocean-acoustics simulations with solid ocean bottoms, Journal of the Acoustical Society of America, vol. 140(1), p.165-175
!       Zhinan Xie, René Matzen, Paul Cristini, Dimitri Komatitsch and Roland Martin, A perfectly matched layer for fluid-solid problems:

module sf_coupling
    use constants, only : fpp
    implicit none

#ifdef CPML

#endif

contains

#ifdef CPML

    ! M_ij = delta_ij . F^-1[s0 s1 s2 / si]
    subroutine compute_convolution_StoF(dom, f0, idxS, idxSF, mu)
        use champs_solidpml
        use m_calcul_forces_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: f0, idxS, idxSF
        real(fpp), intent(out) :: mu(0:2)
        !
        real(fpp) :: k0, a0, d0, k1, a1, d1, cf0, cf1, cf2
        real(fpp), dimension(0:2) :: b0bar, b1bar, b2bar
        integer :: r
        integer :: dir0, dir1
        integer :: Ix, Iy, Iz

        mu(0:2) = 0 ! No convolution

        dir0 = dom%D0_SF(idxSF)
        dir1 = dom%D1_SF(idxSF)

        k0 = dom%Kappa_SF(0, idxSF)
        a0 = dom%Alpha_SF(0, idxSF)
        d0 = dom%dxi_k_SF(0, idxSF)

        k1 = dom%Kappa_SF(1, idxSF)
        a1 = dom%Alpha_SF(1, idxSF)
        d1 = dom%dxi_k_SF(1, idxSF)

        ! Update convolution terms for R0
        do r=0,2
            call cpml_compute_coefs(CPML_INTEG, a0, dom%dt, cf0, cf1, cf2)
            dom%R_0_SF(r, idxSF) = cf0*dom%R_0_SF(r,idxSF)+cf1*dom%champs(f0)%Depla(idxS, r)&
                                  +cf2*dom%champs(f0)%Depla(idxS, r)

            call cpml_compute_coefs(CPML_INTEG, a1, dom%dt, cf0, cf1, cf2)
            if(.not. isclose(a0, a1)) then
                dom%R_1_SF(r, idxSF) = cf0*dom%R_1_SF(r,idxSF)+cf1*dom%champs(f0)%Depla(idxS, r)&
                                      +cf2*dom%champs(f0)%Depla(idxS, r)
            else
                dom%R_1_SF(r, idxSF) = cf0*dom%R_1_SF(r,idxSF)+cf1*dom%R_0_SF(r, idxSF)&
                                      +cf2*dom%R_0_SF(r, idxSF)
            end if
        end do

        if (dom%I1_SF(idxSF) == -1) then
            select case(dir0)
            case(0)
                Ix = 0
                Iy = 1
                Iz = 2
            case(1)
                Ix = 1
                Iy = 0
                Iz = 2
            case(2)
                Ix = 2
                Iy = 1
                Iz = 0
            case default
                stop 1
            end select

            b0bar(Ix) = 1.
            b0bar(Iy) = k0
            b0bar(Iz) = k0
            b1bar(Ix) = 0.
            b1bar(Iy) = k0*d0
            b1bar(Iz) = k0*d0
            b2bar(Ix) = 0.
            b2bar(Iy) = 0.
            b2bar(Iz) = 0.
        else
            select case(dir0)
            case(0)
                Ix = 0
                if(dir1==1) then
                    Iy = 1
                    Iz = 2
                else
                    Iy = 2
                    Iz = 1
                end if
            case(1)
                Ix = 1
                Iy = 2
                Iz = 0
            case default
                stop 1
            end select

            if(.not. isclose(a0, a1)) then
                b0bar(Ix) = k1    ! F-1(s1)
                b0bar(Iy) = k0    ! F-1(s0)
                b0bar(Iz) = k0*k1 ! F-1(s0s1)
                b1bar(Ix) = 0.
                b1bar(Iy) = k0*d0
                b1bar(Iz) = k0*k1*d0*(a0-a1-d1)/(a0-a1)
                b2bar(Ix) = k1*d1
                b2bar(Iy) = 0.
                b2bar(Iz) = k0*k1*d1*(a0-a1+d0)/(a0-a1)
            else
                b0bar(Ix) = k1    ! F-1(s1)
                b0bar(Iy) = k0    ! F-1(s0)
                b0bar(Iz) = k0*k1 ! F-1(s0s1)
                b1bar(Ix) = k1*d1
                b1bar(Iy) = k0*d0
                b1bar(Iz) = k0*k1*(d0+d1)
                b2bar(Ix) = 0.
                b2bar(Iy) = 0.
                b2bar(Iz) = k0*k1*d0*d1
            end if
        end if

        mu(:) =         b0bar(:) * dom%champs(f0)%Depla(idxS, :)
        mu(:) = mu(:) + b1bar(:) * dom%R_0_SF(:, idxSF)
        mu(:) = mu(:) + b2bar(:) * dom%R_1_SF(:, idxSF)
    end subroutine compute_convolution_StoF

    ! N_ij = delta_ij . F^-1[-w^2 s0 s1 s2 / si]
    subroutine compute_convolution_FtoS(dom, f0, idxF, idxSF, nphi)
        use champs_fluidpml
        use m_calcul_forces_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: f0, idxF, idxSF
        real(fpp), intent(out) :: nphi(0:2)
        !
        real(fpp) :: k0, a0, d0, k1, a1, d1, cf0, cf1, cf2
        real(fpp), dimension(0:2) :: a0bar, a1bar, a2bar, a3bar, a4bar
        real(fpp) :: AccPhi
        integer :: dir0, dir1
        integer :: Is0, Is1, Is0s1

        nphi(0:2) = 0 ! No convolution

        ! Convolute the first direction to attenuate (A.5*) from Ref1.

        dir0 = dom%D0_SF(idxSF)
        dir1 = dom%D1_SF(idxSF)

        k0 = dom%Kappa_SF(0, idxSF)
        a0 = dom%Alpha_SF(0, idxSF)
        d0 = dom%dxi_k_SF(0, idxSF)

        a0bar = k0
        a1bar = a0bar * d0
        a2bar = a0bar * d0 * (-a0)
        a3bar = a0bar * a0**2 * d0
        a4bar = 0.

        a0bar(dir0) = 0.
        a1bar(dir0) = 0.
        a2bar(dir0) = 0.
        a3bar(dir0) = 0.

        ! Update convolution terms for R0
        call cpml_compute_coefs(CPML_INTEG, a0, dom%dt, cf0, cf1, cf2)
        dom%R_0_SF(idxSF) = cf0*dom%R_0_SF(idxSF)+cf1*dom%champs(f0)%Phi(idxF)&
                           +cf2*dom%champs(f0)%Phi(idxF)

        ! Convolute the second direction to attenuate (A.5*) from Ref1.

        if(dom%I1_SF(idxSF) /= -1) then
            k1 = dom%Kappa_SF(1, idxSF)
            a1 = dom%Alpha_SF(1, idxSF)
            d1 = dom%dxi_k_SF(1, idxSF)

            if (dir0==0.and.dir1==1) then
                ! For X=0,Y=1
                Is0   = 1
                Is1   = 0
                Is0s1 = 2
            else if (dir0==0.and.dir1==2) then
                ! For X=0,Z=1
                Is0   = 2
                Is1   = 0
                Is0s1 = 1
            else if (dir0==1.and.dir1==2) then
                ! For Y=0,Z=1
                Is0   = 2
                Is1   = 1
                Is0s1 = 0
            else
                stop 1
            end if

            a0bar(Is0)   = k0
            a0bar(Is1)   = k1
            a0bar(Is0s1) = k0 * k1
            a1bar(Is0)   = a0bar(Is0)   * d0
            a1bar(Is1)   = a0bar(Is1)   * d1
            a1bar(Is0s1) = a0bar(Is0s1) * (d0 + d1)
            a2bar(Is0)   = a0bar(Is0)   * d0 * (-a0)
            a2bar(Is1)   = a0bar(Is1)   * d1 * (-a1)
            a2bar(Is0s1) = a0bar(Is0s1) * (d0 * d1 - a0 * d0 - a1 * d1)
            a3bar(Is0)   = a0bar(Is0)   * a0 ** 2 * d0
            a3bar(Is1)   = a0bar(Is1)   * a1 ** 2 * d1
            a3bar(Is0s1) = a0bar(Is0s1) * a0 ** 2 * d0 * (a1 - a0 + d1) / (a1 - a0)
            a4bar(Is0)   = a0bar(Is0)   * a1 ** 2 * d1
            a4bar(Is1)   = a0bar(Is1)   * a0 ** 2 * d0
            a4bar(Is0s1) = a0bar(Is0s1) * a1 ** 2 * d1 * (a0 - a1 + d0) / (a0 - a1)

            ! Update convolution terms for R1
            call cpml_compute_coefs(CPML_INTEG, a1, dom%dt, cf0, cf1, cf2)
            dom%R_1_SF(idxSF) = cf0*dom%R_1_SF(idxSF)+cf1*dom%champs(f0)%Phi(idxF)&
                               +cf2*dom%champs(f0)%Phi(idxF)

        end if

        AccPhi = dom%champs(f0)%ForcesFl(idxF) * dom%MassMat(idxF)
        nphi(:) =           a0bar * AccPhi
        nphi(:) = nphi(:) + a1bar * dom%champs(f0)%VelPhi(idxF)
        nphi(:) = nphi(:) + a2bar * dom%champs(f0)%Phi(idxF)
        nphi(:) = nphi(:) + a3bar * dom%R_0_SF(idxSF)
        nphi(:) = nphi(:) + a4bar * dom%R_1_SF(idxSF)
    end subroutine compute_convolution_FtoS

#endif

subroutine StoF_coupling(Tdomain, f0, f1)
    ! from solid to fluid: velocity (dot) normal
    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: f0, f1
    !
    integer:: ngll_sf, ngll_sf_pml
    integer  :: i,j
    integer :: idxS, idxF
    real(fpp) :: vn
    real(fpp), dimension(0:2) :: btn
#ifdef CPML
    real(fpp) :: mu(0:2)
#else
    real(fpp) :: vn1, vn2, vn3
#endif

    ngll_sf = Tdomain%SF%intSolFlu%surf0%nbtot
    ngll_sf_pml = Tdomain%SF%intSolFluPml%surf0%nbtot

    do i = 0,ngll_sf-1
        idxS = Tdomain%SF%intSolFlu%surf0%map(i)
        idxF = Tdomain%SF%intSolFlu%surf1%map(i)
        BtN = Tdomain%SF%SF_Btn(:,i)
        vn = 0.
        do j = 0,2
#ifdef CPML
            Tdomain%fdom%champs(f1)%ForcesFl(idxF) = Tdomain%fdom%champs(f1)%ForcesFl(idxF) &
                                                     + (BtN(j) * Tdomain%sdom%champs(f0)%Depla(idxS,j))
#else
            Tdomain%fdom%champs(f1)%ForcesFl(idxF) = Tdomain%fdom%champs(f1)%ForcesFl(idxF) &
                                                     + (BtN(j) * Tdomain%sdom%champs(f0)%Veloc(idxS,j))
#endif
        enddo
    enddo

    do i = 0,ngll_sf_pml-1
        idxS = Tdomain%SF%intSolFluPml%surf0%map(i)
        idxF = Tdomain%SF%intSolFluPml%surf1%map(i)
        BtN = Tdomain%SF%SFPml_Btn(:,i)
#ifdef CPML
        ! u_f.n = [M(t)*u_s].n (32) from Ref3. with u_f = grad(phi)
        call compute_convolution_StoF(Tdomain%spmldom, f0, idxS, i, mu)
        vn = (BtN(0) * mu(0)) + (BtN(1) * mu(1)) + (BtN(2) * mu(2))
        Tdomain%fpmldom%champs(f1)%ForcesFl(idxF) = Tdomain%fpmldom%champs(f1)%ForcesFl(idxF) + vn
#else
        vn1 = 0.
        vn2 = 0.
        vn3 = 0.
        do j = 0,2
            vn1 = vn1 + (BtN(j) * Tdomain%spmldom%champs(f0)%VelocPML(idxS,j,0))
            vn2 = vn2 + (BtN(j) * Tdomain%spmldom%champs(f0)%VelocPML(idxS,j,1))
            vn3 = vn3 + (BtN(j) * Tdomain%spmldom%champs(f0)%VelocPML(idxS,j,2))
        enddo
        Tdomain%fpmldom%champs(f1)%fpml_Forces(idxF,0) = Tdomain%fpmldom%champs(f1)%fpml_Forces(idxF,0) + vn1
        Tdomain%fpmldom%champs(f1)%fpml_Forces(idxF,1) = Tdomain%fpmldom%champs(f1)%fpml_Forces(idxF,1) + vn2
        Tdomain%fpmldom%champs(f1)%fpml_Forces(idxF,2) = Tdomain%fpmldom%champs(f1)%fpml_Forces(idxF,2) + vn3
#endif
    enddo
end subroutine StoF_coupling
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine FtoS_coupling(Tdomain, f0, f1)
    ! from fluid to solid: normal times pressure (= -rho . VelPhi)
    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: f0, f1
    !
    real(fpp), dimension(0:2) :: BtN
    integer :: ngll_sf, ngll_sf_pml
    integer :: i,j
    integer :: idxS, idxF
#ifdef CPML
    real(fpp) :: nphi(0:2)
#endif

    ngll_sf = Tdomain%SF%intSolFlu%surf0%nbtot
    ngll_sf_pml = Tdomain%SF%intSolFluPml%surf0%nbtot

    do i = 0,ngll_sf-1
        idxS = Tdomain%SF%intSolFlu%surf0%map(i)
        idxF = Tdomain%SF%intSolFlu%surf1%map(i)
        BtN = Tdomain%SF%SF_Btn(:,i)
        do j = 0,2
#ifdef CPML
            ! Le potentiel Phi est tel que d2Phi/dt2 = -p
            Tdomain%sdom%champs(f1)%Forces(idxS,j) = Tdomain%sdom%champs(f1)%Forces(idxS,j) &
                                                     - (BtN(j) * Tdomain%fdom%champs(f0)%ForcesFl(idxF))
#else
            ! Le potentiel Phi est tel que dPhi/dt = -p
            Tdomain%sdom%champs(f1)%Forces(idxS,j) = Tdomain%sdom%champs(f1)%Forces(idxS,j) &
                                                     - (BtN(j) * Tdomain%fdom%champs(f0)%VelPhi(idxF))
#endif
        enddo
    enddo

    do i = 0,ngll_sf_pml-1
        idxS = Tdomain%SF%intSolFluPml%surf0%map(i)
        idxF = Tdomain%SF%intSolFluPml%surf1%map(i)
        BtN = Tdomain%SF%SFPml_Btn(:,i)
#ifdef CPML
        ! sigma_s.n = [N(t)*Phi].n (32) from Ref3.
        call compute_convolution_FtoS(Tdomain%fpmldom, f0, idxF, i, nphi)
        Tdomain%spmldom%champs(f1)%Forces(idxS,0) = Tdomain%spmldom%champs(f1)%Forces(idxS,0) - (BtN(0)*nphi(0))
        Tdomain%spmldom%champs(f1)%Forces(idxS,1) = Tdomain%spmldom%champs(f1)%Forces(idxS,1) - (BtN(1)*nphi(1))
        Tdomain%spmldom%champs(f1)%Forces(idxS,2) = Tdomain%spmldom%champs(f1)%Forces(idxS,2) - (BtN(2)*nphi(2))
#else
        do j = 0,2
            Tdomain%spmldom%champs(f1)%ForcesPML(idxS,j,0) = Tdomain%spmldom%champs(f1)%ForcesPML(idxS,j,0) &
                                                          - (BtN(j)*Tdomain%fpmldom%champs(f0)%fpml_VelPhi(idxF,0))
            Tdomain%spmldom%champs(f1)%ForcesPML(idxS,j,1) = Tdomain%spmldom%champs(f1)%ForcesPML(idxS,j,1) &
                                                          - (BtN(j)*Tdomain%fpmldom%champs(f0)%fpml_VelPhi(idxF,1))
            Tdomain%spmldom%champs(f1)%ForcesPML(idxS,j,2) = Tdomain%spmldom%champs(f1)%ForcesPML(idxS,j,2) &
                                                          - (BtN(j)*Tdomain%fpmldom%champs(f0)%fpml_VelPhi(idxF,2))
        enddo
#endif
    enddo
end subroutine FtoS_coupling
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
end module sf_coupling
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
