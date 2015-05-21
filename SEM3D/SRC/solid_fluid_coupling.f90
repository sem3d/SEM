!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine StoF_coupling(Tdomain, Veloc, BtN, ForcesFl)
    ! from solid to fluid: velocity (dot) normal
    use sdomain
    use scomm
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer:: ngll_sf
    real, intent(in), dimension(0:Tdomain%ngll_s-1,0:2) :: Veloc
    real, intent(in), dimension(0:Tdomain%SF%ngll-1,0:2) :: BtN
    real, intent(inout), dimension(0:Tdomain%ngll_f-1) :: ForcesFl
    integer  :: i,j,k,n,idx
    real, dimension(0:Tdomain%SF%ngll-1) :: vn

    ngll_sf = Tdomain%SF%ngll

    vn(:) = 0d0
    do i = 0,ngll_sf-1
        idx = Tdomain%sf%SF_IGlobSol(i)
        if (idx >= 0) then
            do j = 0,2
                vn(i) = vn(i) + (BtN(i,j) * Veloc(idx,j))
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
            ForcesFl(idx) = ForcesFl(idx) + vn(i)
        end if
    enddo

end subroutine StoF_coupling
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine FtoS_coupling(Tdomain, BtN, VelPhi, Forces)
    ! from fluid to solid: normal times pressure (= -rho . VelPhi)
    use sdomain
    use scomm
    implicit none

    type(domain), intent(inout) :: Tdomain
    real, intent(in), dimension(0:Tdomain%SF%ngll-1,0:2) :: BtN
    real, intent(in), dimension(0:Tdomain%ngll_f-1) :: VelPhi
    real, intent(inout), dimension(0:Tdomain%ngll_f-1,0:2) :: Forces
    !
    integer :: ngll_sf
    integer  :: i,j,k,n,idx
    real, dimension(0:2, 0:Tdomain%SF%ngll-1) :: pn

    ngll_sf = Tdomain%SF%ngll
    pn = 0d0
    do i = 0,ngll_sf-1
        idx = Tdomain%SF%SF_IGlobFlu(i)
        do j = 0,2
            pn(j,i) = - (BtN(i,j) * VelPhi(idx))
        enddo
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
        do j = 0,2
            Forces(idx,j) = Forces(idx,j) + pn(j,i)
        enddo
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
