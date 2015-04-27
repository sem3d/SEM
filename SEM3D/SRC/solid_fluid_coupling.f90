!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine StoF_coupling(Tdomain,ngll_sol, ngll_flu, SF_ngll, SF_IGlobSol, SF_IGlobFlu, Veloc, BtN, ForcesFl)
    ! from solid to fluid: velocity (dot) normal
    use sdomain
    use scomm
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: ngll_sol, ngll_flu, SF_ngll
    integer, intent(in), dimension(0:SF_ngll-1) :: SF_IGlobSol, SF_IGlobFlu
    real, intent(in), dimension(0:ngll_sol-1,0:2) :: Veloc
    real, intent(in), dimension(0:SF_ngll-1,0:2) :: BtN
    real, intent(inout), dimension(0:ngll_flu-1) :: ForcesFl
    integer  :: i,j,k,n,idx
    real, dimension(0:SF_ngll-1) :: vn

    vn(:) = 0d0
    do i = 0,SF_ngll-1
        if (SF_IGlobSol(i) < 0) cycle ! solid face not on this proc
        do j = 0,2
            vn(i) = vn(i) + (BtN(i,j) * Veloc(SF_IGlobSol(i),j))
        enddo
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

    do i = 0,SF_ngll-1
        if (SF_IGlobFlu(i) < 0) cycle ! fluid face not on this proc
        ForcesFl(SF_IGlobFlu(i)) = ForcesFl(SF_IGlobFlu(i)) + vn(i)
    enddo

end subroutine StoF_coupling
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine FtoS_coupling(Tdomain,ngll_sol, ngll_flu, SF_ngll, SF_IGlobSol, SF_IGlobFlu, BtN, Save_forces, Save_depla, VelPhi, Forces, Depla)
    ! from fluid to solid: normal times pressure (= -rho . VelPhi)
    use sdomain
    use scomm
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: ngll_sol, ngll_flu, SF_ngll
    integer, intent(in), dimension(0:SF_ngll-1) :: SF_IGlobSol, SF_IGlobFlu
    real, intent(in), dimension(0:SF_ngll-1,0:2) :: Save_forces, Save_depla, BtN
    real, intent(in), dimension(0:ngll_flu-1) :: VelPhi
    real, intent(inout), dimension(0:ngll_sol-1,0:2) :: Forces, Depla

    integer  :: i,j,k,n,idx
    real, dimension(0:SF_ngll-1) :: pn

    pn(:) = 0d0
    do i = 0,SF_ngll-1
        do j = 0,2
            pn(i) = pn(i) + (BtN(i,j) * VelPhi(SF_IGlobFlu(i)))
        enddo
    enddo

    if (Tdomain%Comm_SolFlu%ncomm > 0)then
        ! Give
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            k = 0
            do i=0,Tdomain%Comm_SolFlu%Data(n)%ndata-1
                idx = Tdomain%Comm_SolFlu%Data(n)%IGiveS(i)
                Tdomain%Comm_SolFlu%Data(n)%Give(k) = pn(idx)
                k=k+1
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
                pn(idx) = pn(idx) + Tdomain%Comm_SolFlu%Data(n)%Take(k)
                k=k+1
            enddo
        enddo
    endif

    do j = 0,2
        do i = 0,SF_ngll-1
            Forces(SF_IGlobSol(i),j) = Save_forces(i,j) + pn(i)
            Depla(SF_IGlobSol(i),j) = Save_depla(i,j)
        enddo
    enddo

end subroutine FtoS_coupling
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine SF_solid_values_saving(ngll_sol, SF_ngll, SF_IGlobSol, Forces, Depla, Save_forces, Save_depla)
    ! saves values of fields Forces, Displ on the solid side of a SF object -
    !    for they are flushed in the 1st correction phase
    implicit none

    integer, intent(in) :: ngll_sol, SF_ngll
    integer, dimension(0:SF_ngll-1), intent(in) :: SF_IGlobSol
    real, dimension(0:ngll_sol-1,0:2), intent(in) :: Forces, Depla
    real, dimension(0:SF_ngll-1,0:2), intent(inout) :: Save_forces, Save_depla
    integer  :: i,idx
    
    do i = 0, SF_ngll-1
        idx = SF_IGlobSol(i)
        if (idx /= -1) then
            Save_forces(i,:) = Forces(SF_IGlobSol(i),:)
            Save_depla(i,:) = Depla(SF_IGlobSol(i),:)
        else
            Save_forces(i,:) = 0.
            Save_depla(i,:) = 0.
        endif
    enddo

    return
end subroutine SF_solid_values_saving
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine Newmark_recorrect_solid(ngll_sol, SF_ngll, dt, SF_IGlobSol, MassMatSol, Forces, Veloc, Depla)
    implicit none

    integer, intent(in) :: ngll_sol, SF_ngll
    real, intent(in) :: dt
    integer, intent(in), dimension(0:SF_ngll-1) :: SF_IGlobSol
    real, intent(in), dimension(0:ngll_sol-1) :: MassMatSol
    real, intent(inout), dimension(0:ngll_sol-1,0:2) :: Forces, Veloc, Depla
    integer  :: i,j,n
    
    do j = 0,2
        do n = 0,SF_ngll-1
            i = SF_IGlobSol(n)
            Forces(i,j) = Forces(i,j) * MassMatSol(i)
            Veloc(i,j) = Veloc(i,j) + dt * Forces(i,j)
            Depla(i,j) = Depla(i,j) + dt * Veloc(i,j)
        enddo
    enddo

end subroutine Newmark_recorrect_solid

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
