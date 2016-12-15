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

    integer, parameter :: dXX=0, dXY=1, dXZ=2
    integer, parameter :: dYX=3, dYY=4, dYZ=5
    integer, parameter :: dZX=6, dZY=7, dZZ=8

#endif

contains

#ifdef CPML

    ! M_ij = delta_ij . F^-1[s0 s1 s2 / si]
    subroutine compute_convolution_StoF(dom, f0, idxS, idxSF, ee, bnum, mu)
        use champs_solidpml
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: f0, idxS, idxSF, ee, bnum
        real(fpp), intent(out) :: mu(0:2)
        !
        real(fpp) :: k0, a0, d0, b0, k1, a1, d1, b1
        real(fpp) :: b0bar, b1bar, b2bar
        integer :: i1, i2, dir

        mu(0:2) = 1 ! No convolution

        ! Convolute the first direction to attenuate (A.20*) from Ref1. (s2 = 1, s3 cancels)

        i1 = dom%I1(ee,bnum)
        if (i1 /= -1) then
            dir = -1
            if (i1 == 0) dir = dXX
            if (i1 == 1) dir = dYY
            if (i1 == 2) dir = dZZ
            if (dir == -1) stop "compute_convolution_StoF - invalid i1"

            k0 = dom%Kappa_SF(idxSF, i1)
            a0 = dom%Alpha_SF(idxSF, i1)
            d0 = dom%dxi_k_SF(idxSF, i1)
            b0 = a0 + d0

            b0bar = k0
            b1bar = -b0bar * (-d0)

            mu(dir) =           b0bar * dom%champs(f0)%Depla(idxS, dir)
            mu(dir) = mu(dir) + b1bar * dom%champs(f0)%Depla(idxS, dir) * dom%R2_SF(i1, idxSF)
        end if

        ! Convolute the second direction to attenuate (A.20*) from Ref1. (s2 != 1, s3 cancels)

        i2 = dom%I2(ee,bnum)
        if (i2 /= -1) then
            dir = -1
            if (i2 == 0) dir = dXX
            if (i2 == 1) dir = dYY
            if (i2 == 2) dir = dZZ
            if (dir == -1) stop "compute_convolution_StoF - invalid i2"
            if (i1 == -1) stop "compute_convolution_StoF - invalid i1 for i2"

            k1 = dom%Kappa_SF(idxSF, i2)
            a1 = dom%Alpha_SF(idxSF, i2)
            d1 = dom%dxi_k_SF(idxSF, i2)
            b1 = a1 + d1

            b0bar = k0 * k1
            b1bar = -b0bar * (-d0) * (a0 - b1) / (a0 - a1)
            b2bar = -b0bar * (-d1) * (a1 - b0) / (a1 - a0)

            mu(dir) =           b0bar * dom%champs(f0)%Depla(idxS, dir)
            mu(dir) = mu(dir) + b1bar * dom%champs(f0)%Depla(idxS, dir) * dom%R2_SF(i1, idxSF)
            mu(dir) = mu(dir) + b2bar * dom%champs(f0)%Depla(idxS, dir) * dom%R2_SF(i2, idxSF)
        end if
    end subroutine compute_convolution_StoF

    ! N_ij = delta_ij . F^-1[-w^2 s0 s1 s2 / si]
    subroutine compute_convolution_FtoS(dom, f0, idxF, idxSF, ee, bnum, nphi)
        use champs_fluidpml
        implicit none
        type(domain_fluidpml), intent (INOUT) :: dom
        integer, intent(in) :: f0, idxF, idxSF, ee, bnum
        real(fpp), intent(out) :: nphi(0:2)
        !
        real(fpp) :: k0, a0, d0, b0, k1, a1, d1, b1
        real(fpp) :: a0bar, a1bar, a2bar, a3bar, a4bar
        real(fpp) :: AccPhi
        integer :: i1, i2, dir

        nphi(0:2) = 1 ! No convolution

        ! Convolute the first direction to attenuate (A.5*) from Ref1. (s2 = 1, s3 cancels)

        i1 = dom%I1(ee,bnum)
        if (i1 /= -1) then
            dir = -1
            if (i1 == 0) dir = dXX
            if (i1 == 1) dir = dYY
            if (i1 == 2) dir = dZZ
            if (dir == -1) stop "compute_convolution_FtoS - invalid i1"

            k0 = dom%Kappa_SF(idxSF, i1)
            a0 = dom%Alpha_SF(idxSF, i1)
            d0 = dom%dxi_k_SF(idxSF, i1)
            b0 = a0 + d0

            a0bar = k0
            a1bar = a0bar * d0
            a2bar = a0bar * d0 * (-a0)
            a3bar = a0bar * a0**2 * d0

            AccPhi = dom%champs(f0)%ForcesFl(idxF) * dom%MassMat(idxF)
            nphi(dir) =             a0bar * AccPhi
            nphi(dir) = nphi(dir) + a1bar * dom%champs(f0)%VelPhi(idxF)
            nphi(dir) = nphi(dir) + a2bar * dom%champs(f0)%Phi(idxF)
            nphi(dir) = nphi(dir) + a3bar * dom%R1_SF(i1, idxSF)
        end if

        ! Convolute the second direction to attenuate (A.5*) from Ref1. (s2 != 1, s3 cancels)

        i2 = dom%i2(ee,bnum)
        if (i2 /= -1) then
            dir = -1
            if (i2 == 0) dir = dXX
            if (i2 == 1) dir = dYY
            if (i2 == 2) dir = dZZ
            if (dir == -1) stop "compute_convolution_FtoS - invalid i2"
            if (i1 == -1) stop "compute_convolution_FtoS - invalid i1 for i2"

            k1 = dom%Kappa_SF(idxSF, i2)
            a1 = dom%Alpha_SF(idxSF, i2)
            d1 = dom%dxi_k_SF(idxSF, i2)
            b1 = a1 + d1

            a0bar = k0 * k1
            a1bar = a0bar * d0 + a0bar * d1
            a2bar = a0bar * d0 * (d1 - a0) + a0bar * d1 * (-a1)
            a3bar = a0bar * a0**2 * d0 * (b1 - a0) / (a1 - a0)
            a4bar = a0bar * a1**2 * d1 * (b0 - a1) / (a0 - a1)

            AccPhi = dom%champs(f0)%ForcesFl(idxF) * dom%MassMat(idxF)
            nphi(dir) =             a0bar * AccPhi
            nphi(dir) = nphi(dir) + a1bar * dom%champs(f0)%VelPhi(idxF)
            nphi(dir) = nphi(dir) + a2bar * dom%champs(f0)%Phi(idxF)
            nphi(dir) = nphi(dir) + a3bar * dom%R1_SF(i1, idxSF)
            nphi(dir) = nphi(dir) + a4bar * dom%R1_SF(i2, idxSF)
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
                                                     + (BtN(j) * Tdomain%sdom%champs(f0)%Forces(idxS,j))
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
        !call compute_convolution_StoF(Tdomain%spmldom, f0, idxS, ??, ??, ??, ??, ??, mu)
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
        !call compute_convolution_FtoS(Tdomain%fpmldom, f0, idxF, ??, ??, ??, ??, ??, nphi)
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
