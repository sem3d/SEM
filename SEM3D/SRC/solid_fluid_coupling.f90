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

    subroutine calc_directions(dir0, dir1, flag0, Ix, Iy, Iz)
        integer, intent(in) :: dir0, dir1, flag0
        integer, intent(out) :: Ix, Iy, Iz

        Ix = dir0
        if (flag0 == -1) then
            ! only one direction of attenuation
            ! y / z should have a symetric role
            ! keep the frame orientation direct
            select case(dir0)
            case(0)
                Iy = 1
                Iz = 2
            case(1)
                Iy = 0
                Iz = 2
            case(2)
                Iy = 1
                Iz = 0
            case default
                stop 1
            end select
        else
            ! two directions
            Iy = dir1
            select case(dir0)
            case(0)
                if(dir1==1) then
                    Iz = 2
                else
                    Iz = 1
                end if
            case(1)
                Iz = 0
            case default
                stop 1
            end select
        end if
    end subroutine calc_directions

    ! M_ij = delta_ij . F^-1[s0 s1 s2 / si]
    subroutine compute_convolution_StoF(dom, f0, idxS, idxSF, mu)
        use champs_solidpml
        use m_calcul_forces_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: f0, idxS, idxSF
        real(fpp), intent(out) :: mu(0:2)
        !
        real(fpp) :: k0, a0, d0, k1, a1, d1, cf0, cf1
        real(fpp) :: R0y, R0z, dR0y, dR0z, R0yn, R0zn
        real(fpp) :: R1x, R1z, dR1x, dR1z, R1xn, R1zn
        real(fpp) :: Ux, Uy, Uz
        real(fpp) :: b1bar, b2bar
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
        call calc_directions(dir0, dir1, dom%I1_SF(idxSF), Ix, Iy, Iz)

        Ux = dom%champs(f0)%Depla(idxS, Ix)
        Uy = dom%champs(f0)%Depla(idxS, Iy)
        Uz = dom%champs(f0)%Depla(idxS, Iz)

        ! Update convolution terms for R0
        ! R0(0) = exp(-a0.t)*Uy ; R0(1) = exp(-a0.t)*Uz ;
        call cpml_coefs_midpoint2(a0, dom%dt, cf0, cf1)
        R0y = dom%R_0_SF(0, idxSF)
        R0z = dom%R_0_SF(1, idxSF)
        dR0y = cf0*R0y + cf1*Uy
        dR0z = cf0*R0z + cf1*Uz
        dom%R_0_SF(0, idxSF) = R0y + dR0y
        dom%R_0_SF(1, idxSF) = R0z + dR0z
        R0yn = R0y + 0.5*dR0y
        R0zn = R0z + 0.5*dR0z


        if (dom%I1_SF(idxSF) /= -1) then
            ! Update convolution terms for R1
            ! R1(0) = exp(-a1.t)*Ux ; R1(1) = exp(-a1.t)*Uz ;
            ! if a0==a1 : R1(1) = t.exp(-a1.t)*Uz ;
            ! if no attenuation we should have a1=0 hence cf0=0
            !write(*,*) idxSF, "S->F:0", k0, a0, d0
            !write(*,*) idxSF, "S->F:1", k1, a1, d1
            call cpml_coefs_midpoint2(a1, dom%dt, cf0, cf1)
            R1x = dom%R_1_SF(0, idxSF)
            R1z = dom%R_1_SF(1, idxSF)
            dR1x = cf0*R1x + cf1*Ux
            if(.not. isclose(a0, a1)) then
                dR1z = cf0*R1z + cf1*Uz
            else
                dR1z = cf0*R1z + cf1*R0zn
            end if
            dom%R_1_SF(0, idxSF) = R1x + dR1x
            dom%R_1_SF(1, idxSF) = R1z + dR1z
            R1xn = R1x + 0.5*dR1x
            R1zn = R1z + 0.5*dR1z
        else
            R1xn = 0
            R1zn = 0
        end if

        if (dom%I1_SF(idxSF) == -1) then
            ! From here, attenuation is along Ix.
            mu(Ix) = 1.*Ux
            mu(Iy) = k0*Uy + k0*d0*R0yn
            mu(Iz) = k0*Uz + k0*d0*R0zn
        else
            ! From here, attenuation is along Ix, Iy.
            if(.not. isclose(a0, a1)) then
                b1bar = k0*k1*d0*(a0-a1-d1)/(a0-a1)
                b2bar = k0*k1*d1*(a0-a1+d0)/(a0-a1)
            else
                b1bar = k0*k1*(d0+d1)
                b2bar = k0*k1*d0*d1
            end if
            mu(Ix) =    k1*Ux + k1*d1*R1xn
            mu(Iy) =    k0*Uy + k0*d0*R0yn
            mu(Iz) = k0*k1*Uz + b1bar*R0zn + b2bar*R1zn
            mu(0:2) = 0.
        end if
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
        real(fpp) :: k0, a0, d0, k1, a1, d1, cf0, cf1, R0, R0n, dR0, R1, R1n, dR1
        real(fpp) :: AccPhi, VelPhi, Phi
        integer :: dir0, dir1
        integer :: Ix, Iy, Iz

        nphi(0:2) = 0 ! No convolution

        ! Convolute the first direction to attenuate (A.5*) from Ref1.

        dir0 = dom%D0_SF(idxSF)
        dir1 = dom%D1_SF(idxSF)

        k0 = dom%Kappa_SF(0, idxSF)
        a0 = dom%Alpha_SF(0, idxSF)
        d0 = dom%dxi_k_SF(0, idxSF)

        k1 = dom%Kappa_SF(1, idxSF)
        a1 = dom%Alpha_SF(1, idxSF)
        d1 = dom%dxi_k_SF(1, idxSF)
        call calc_directions(dir0, dir1, dom%I1_SF(idxSF), Ix, Iy, Iz)

        AccPhi = dom%champs(f0)%ForcesFl(idxF)
        VelPhi = dom%champs(f0)%VelPhi(idxF)
        Phi    = dom%champs(f0)%Phi(idxF)
        ! Update convolution terms for R0
        ! R0(0) = exp(-a0.t)*X
        call cpml_coefs_midpoint2(a0, dom%dt, cf0, cf1)
        R0 = dom%R_0_SF(idxSF)
        dR0 = cf0*R0 + cf1*Phi
        dom%R_0_SF(idxSF) = R0 + dR0
        R0n = R0 + 0.5*dR0

        ! Convolute the second direction to attenuate (A.5*) from Ref1.
        if(dom%I1_SF(idxSF) /= -1) then
            !write(*,*) idxSF, "F->S:0", k0, a0, d0
            !write(*,*) idxSF, "F->S:1", k1, a1, d1
            ! Update convolution terms for R1
            ! R1(0) = exp(-a1.t)*X R1(1) = t.exp(-a1.t)*X
            ! R1(1) used only if a0==a1
            call cpml_coefs_midpoint2(a1, dom%dt, cf0, cf1)
            R1 = dom%R_1_SF(idxSF)
            if(.not. isclose(a0, a1)) then
                dR1 = cf0*R1 + cf1*Phi
            else
                dR1 = cf0*R1 + cf1*R0n
            end if
            dom%R_1_SF(idxSF) = R1 + dR1
            R1n = R1 + 0.5*dR1

            if(.not. isclose(a0, a1)) then
                nphi(Ix) = k1*(AccPhi + d1*(VelPhi - a1*(Phi + a1*R1n)))
                nphi(Iy) = k0*(AccPhi + d0*(VelPhi - a0*(Phi + a0*R0n)))
                nphi(Iz) = k0*k1*(AccPhi + (d0+d1)*VelPhi + (d0*d1-a0*d0-a1*d1)*Phi)
                nphi(Iz) = nphi(Iz) + k0*k1*(a0*a0*d0*(a0-a1-d1)*R0n+a1*a1*d1*(a0-a1+d0)*R1n )/(a0-a1)
            else
                nphi(Ix) = k1*(AccPhi + d1*(VelPhi - a1*(Phi + a1*R0n))) ! R0n is right since R1n is /=
                nphi(Iy) = k0*(AccPhi + d0*(VelPhi - a0*(Phi + a0*R0n)))
                nphi(Iz) = k0*k1*(AccPhi + (d0+d1)*VelPhi + (d0*d1-a0*d0-a1*d1)*Phi)
                nphi(Iz) = nphi(Iz) + k0*k1*a0*( (a0*(d0+d1)-2*d0*d1)*R0n + a0*d0*d1*R1n)
            end if
            nphi(0:2) = 0.
        else
            nphi(Ix) =     AccPhi
            nphi(Iy) = k0*(AccPhi + d0*(VelPhi - a0*(Phi + a0*R0n)))
            nphi(Iz) = nphi(Iy)
        end if
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
