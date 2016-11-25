!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!----------------------------------------------------------------
!----------------------------------------------------------------
module sf_coupling
    implicit none

contains
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
    real :: vn, vn1, vn2, vn3
    real, dimension(0:2) :: btn

    ngll_sf = Tdomain%SF%intSolFlu%surf0%nbtot
    ngll_sf_pml = Tdomain%SF%intSolFluPml%surf0%nbtot

    do i = 0,ngll_sf-1
        idxS = Tdomain%SF%intSolFlu%surf0%map(i)
        idxF = Tdomain%SF%intSolFlu%surf1%map(i)
        BtN = Tdomain%SF%SF_Btn(:,i)
        vn = 0.
        do j = 0,2
            Tdomain%fdom%champs(f1)%ForcesFl(idxF) = Tdomain%fdom%champs(f1)%ForcesFl(idxF) &
                                                  + (BtN(j) * Tdomain%sdom%champs(f0)%Veloc(idxS,j))
        enddo
    enddo

    do i = 0,ngll_sf_pml-1
        idxS = Tdomain%SF%intSolFluPml%surf0%map(i)
        idxF = Tdomain%SF%intSolFluPml%surf1%map(i)
        BtN = Tdomain%SF%SFPml_Btn(:,i)
        vn1 = 0.
        vn2 = 0.
        vn3 = 0.
#ifdef CPML
        vn1 = vn1 + (BtN(0) * Tdomain%spmldom%champs(f0)%Veloc(idxS,0))
        vn2 = vn2 + (BtN(1) * Tdomain%spmldom%champs(f0)%Veloc(idxS,1))
        vn3 = vn3 + (BtN(2) * Tdomain%spmldom%champs(f0)%Veloc(idxS,2))
#else
        do j = 0,2
            vn1 = vn1 + (BtN(j) * Tdomain%spmldom%champs(f0)%VelocPML(idxS,j,0))
            vn2 = vn2 + (BtN(j) * Tdomain%spmldom%champs(f0)%VelocPML(idxS,j,1))
            vn3 = vn3 + (BtN(j) * Tdomain%spmldom%champs(f0)%VelocPML(idxS,j,2))
        enddo
#endif
#ifdef CPML
        !TODO
        !Tdomain%fpmldom%champs(f1)%ForcesFl(idxF,0) = Tdomain%fpmldom%champs(f1)%ForcesFl(idxF,0) + vn1
        !Tdomain%fpmldom%champs(f1)%ForcesFl(idxF,1) = Tdomain%fpmldom%champs(f1)%ForcesFl(idxF,1) + vn2
        !Tdomain%fpmldom%champs(f1)%ForcesFl(idxF,2) = Tdomain%fpmldom%champs(f1)%ForcesFl(idxF,2) + vn3
#else
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
    real, dimension(0:2) :: BtN
    integer :: ngll_sf, ngll_sf_pml
    integer :: i,j
    integer :: idxS, idxF

    ngll_sf = Tdomain%SF%intSolFlu%surf0%nbtot
    ngll_sf_pml = Tdomain%SF%intSolFluPml%surf0%nbtot

    do i = 0,ngll_sf-1
        idxS = Tdomain%SF%intSolFlu%surf0%map(i)
        idxF = Tdomain%SF%intSolFlu%surf1%map(i)
        BtN = Tdomain%SF%SF_Btn(:,i)
        do j = 0,2
            Tdomain%sdom%champs(f1)%Forces(idxS,j) = Tdomain%sdom%champs(f1)%Forces(idxS,j) &
                                                  - (BtN(j) * Tdomain%fdom%champs(f0)%VelPhi(idxF))
        enddo
    enddo

    do i = 0,ngll_sf_pml-1
        idxS = Tdomain%SF%intSolFluPml%surf0%map(i)
        idxF = Tdomain%SF%intSolFluPml%surf1%map(i)
        BtN = Tdomain%SF%SFPml_Btn(:,i)
#ifdef CPML
        !Tdomain%spmldom%Forces(idxS,0) = Tdomain%spmldom%Forces(idxS,0) &
        !                                 - (BtN(0)*Tdomain%fpmldom%champs(f0)%fpml_VelPhi(idxF,0))
        !Tdomain%spmldom%Forces(idxS,1) = Tdomain%spmldom%Forces(idxS,1) &
        !                                 - (BtN(1)*Tdomain%fpmldom%champs(f0)%fpml_VelPhi(idxF,1))
        !Tdomain%spmldom%Forces(idxS,2) = Tdomain%spmldom%Forces(idxS,2) &
        !                                 - (BtN(2)*Tdomain%fpmldom%champs(f0)%fpml_VelPhi(idxF,2))
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
