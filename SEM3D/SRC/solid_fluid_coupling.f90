!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine StoF_coupling(Tdomain, BtN, champs0, champs1)
    ! from solid to fluid: velocity (dot) normal
    use sdomain
    use scomm
    use schamps
    implicit none

    type(domain), intent(inout) :: Tdomain
    real, intent(in), dimension(0:Tdomain%SF%ngll-1,0:2) :: BtN
    type(champs), intent(in) :: champs0
    type(champs), intent(inout) :: champs1
    !
    integer:: ngll_sf, ngll_sf_pml
    integer  :: i,j,k,n,idx
    real, dimension(0:Tdomain%SF%ngll-1) :: vn
    real, dimension(0:Tdomain%SF%ngll_pml-1) :: vn1, vn2, vn3

    ngll_sf = Tdomain%SF%ngll
    ngll_sf_pml = Tdomain%SF%ngll_pml

    vn(:) = 0d0
    do i = 0,ngll_sf-1
        idx = Tdomain%sf%SF_IGlobSol(i)
        if (idx >= 0) then
            do j = 0,2
                vn(i) = vn(i) + (BtN(i,j) * champs0%Veloc(idx,j))
            enddo
        end if
    enddo
    vn1 = 0.
    vn2 = 0.
    vn3 = 0.
    do i = 0,ngll_sf_pml-1
        idx = Tdomain%sf%SF_IGlobSol_pml(i)
        if (idx >= 0) then
            !! sum veloc1+veloc2+veloc3
            do j = 0,2
                vn1(i) = vn1(i) + (BtN(i,j) * champs0%VelocPML(idx+0,j))
                vn2(i) = vn2(i) + (BtN(i,j) * champs0%VelocPML(idx+1,j))
                vn3(i) = vn3(i) + (BtN(i,j) * champs0%VelocPML(idx+2,j))
            enddo
        end if
    enddo

    if (Tdomain%Comm_SolFlu%ncomm > 0)then
        ! Give
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            k = 0
            do i=0,Tdomain%Comm_SolFlu%Data(n)%ndata-1
                idx = Tdomain%Comm_SolFlu%Data(n)%IGiveS(i)
                Tdomain%Comm_SolFlu%Data(n)%Give(k) = vn(idx)
                k=k+1
            end do
            Tdomain%Comm_SolFlu%Data(n)%nsend = k
        enddo

        ! Exchange interproc of vn
        call exchange_sem_var(Tdomain, 301, Tdomain%Comm_SolFlu)

        ! Take
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            k = 0
            do i=0,Tdomain%Comm_SolFlu%Data(n)%ndata-1
                idx = Tdomain%Comm_SolFlu%Data(n)%ITakeS(i)
                vn(idx) = vn(idx) + Tdomain%Comm_SolFlu%Data(n)%Take(k)
                k=k+1
            enddo
        enddo
    endif

    do i = 0,ngll_sf-1
        idx = TDomain%SF%SF_IGlobFlu(i)
        if (idx >= 0) then
            champs1%ForcesFl(idx) = champs1%ForcesFl(idx) + vn(i)
        end if
    enddo

    do i = 0,ngll_sf_pml-1
        idx = TDomain%SF%SF_IGlobFlu_pml(i)
        if (idx >= 0) then
            champs1%fpml_Forces(idx+0) = champs1%fpml_Forces(idx+0) + vn1(i)
            champs1%fpml_Forces(idx+1) = champs1%fpml_Forces(idx+1) + vn2(i)
            champs1%fpml_Forces(idx+2) = champs1%fpml_Forces(idx+2) + vn3(i)
        end if
    enddo

end subroutine StoF_coupling
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine FtoS_coupling(Tdomain, BtN, champs0, champs1)
    ! from fluid to solid: normal times pressure (= -rho . VelPhi)
    use sdomain
    use scomm
    implicit none

    type(domain), intent(inout) :: Tdomain
    real, intent(in), dimension(0:Tdomain%SF%ngll-1,0:2) :: BtN
    type(champs), intent(in) :: champs0
    type(champs), intent(inout) :: champs1
    !
    integer :: ngll_sf, ngll_sf_pml
    integer  :: i,j,k,n,idx
    real, dimension(0:2, 0:Tdomain%SF%ngll-1) :: pn, pn1, pn2, pn3

    ngll_sf = Tdomain%SF%ngll
    ngll_sf_pml = Tdomain%SF%ngll_pml
    pn = 0d0
    do i = 0,ngll_sf-1
        idx = Tdomain%SF%SF_IGlobFlu(i)
        if (idx >= 0) then
            do j = 0,2
                pn(j,i) = - (BtN(i,j) * champs0%VelPhi(idx))
            enddo
        end if
    enddo
    pn1 = 0d0
    pn2 = 0d0
    pn3 = 0d0
    do i = 0,ngll_sf_pml-1
        idx = Tdomain%SF%SF_IGlobFlu_pml(i)
        if (idx >= 0) then
            do j = 0,2
                pn1(j,i) =  - (BtN(i,j) * champs0%fpml_VelPhi(idx+0))
                pn2(j,i) =  - (BtN(i,j) * champs0%fpml_VelPhi(idx+1))
                pn3(j,i) =  - (BtN(i,j) * champs0%fpml_VelPhi(idx+2))
            enddo
        end if
    enddo

    if (Tdomain%Comm_SolFlu%ncomm > 0)then
        ! Give
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            k = 0
            do i=0,Tdomain%Comm_SolFlu%Data(n)%ndata-1
                idx = Tdomain%Comm_SolFlu%Data(n)%IGiveS(i)
                Tdomain%Comm_SolFlu%Data(n)%Give(k+0) = pn(0,idx)
                Tdomain%Comm_SolFlu%Data(n)%Give(k+1) = pn(1,idx)
                Tdomain%Comm_SolFlu%Data(n)%Give(k+2) = pn(2,idx)
                k=k+3
            end do
            Tdomain%Comm_SolFlu%Data(n)%nsend = k
        enddo

        ! Exchange interproc of vn
        call exchange_sem_var(Tdomain, 302, Tdomain%Comm_SolFlu)

        ! Take
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            k = 0
            do i=0,Tdomain%Comm_SolFlu%Data(n)%ndata-1
                idx = Tdomain%Comm_SolFlu%Data(n)%ITakeS(i)
                pn(0,idx) = pn(0,idx) + Tdomain%Comm_SolFlu%Data(n)%Take(k+0)
                pn(1,idx) = pn(1,idx) + Tdomain%Comm_SolFlu%Data(n)%Take(k+1)
                pn(2,idx) = pn(2,idx) + Tdomain%Comm_SolFlu%Data(n)%Take(k+2)
                k=k+3
            enddo
        enddo
    endif

    do i = 0,ngll_sf-1
        idx = Tdomain%SF%SF_IGlobSol(i)
        if (idx >= 0) then
            do j = 0,2
                champs1%Forces(idx,j) = champs1%Forces(idx,j) + pn(j,i)
            enddo
        end if
    enddo
    do i = 0,ngll_sf_pml-1
        idx = Tdomain%SF%SF_IGlobSol_pml(i)
        if (idx >= 0) then
            do j = 0,2
                champs1%ForcesPML(idx+0,j) = champs1%ForcesPML(idx+0,j) + pn1(j,i)
                champs1%ForcesPML(idx+1,j) = champs1%ForcesPML(idx+1,j) + pn2(j,i)
                champs1%ForcesPML(idx+2,j) = champs1%ForcesPML(idx+2,j) + pn3(j,i)
            enddo
        end if
    enddo

end subroutine FtoS_coupling
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

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
