!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file courant.f90
!!\brief Contient la subroutine Compute_Courant.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \fn subroutine Compute_Courant (Tdomain)
!! \brief Vérifie que le pas de temps est adapté aux déplacements évalués
!!
!! \param type (Domain), intent (IN) Tdomain
!<

subroutine update_minmax(Tdomain, idef0, idef1, dxmin, dxmax)
    use sdomain
    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent(in) :: idef0, idef1
    real, intent(inout) :: dxmin, dxmax
    !
    real :: dx

    if (idef0==idef1) then
        write(*,*) "Numbering problem"
        stop 1
    end if
    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
    if(dx < dxmin) dxmin = dx
    if(dx > dxmax) dxmax = dx
end subroutine update_minmax

subroutine Compute_Courant (Tdomain,rg)

    use sdomain
    use mpi
    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent(IN) :: rg

    integer :: i,j,k,n,ngll,idef0,idef1, mat, ierr
    real :: dxmin,courant,courant_max
    real :: dt_min, dt, dt_loc
    real :: floc_max, dxmax, f_max

    courant = Tdomain%TimeD%courant
    courant_max = 0
    dt_loc = huge(1.)
    floc_max  = huge(1.)
    do n = 0, Tdomain%n_elem -1
        dxmin = 1e10
        dxmax = 0

        ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
        do k = 0, ngll-2
            do j = 0, ngll - 2
                do i = 0, ngll -2

                    idef0 = Tdomain%specel(n)%Iglobnum(i,j,k)
                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j,k)
                    call update_minmax(Tdomain, idef0, idef1, dxmin, dxmax)
                    idef1 = Tdomain%specel(n)%Iglobnum(i,j+1,k)
                    call update_minmax(Tdomain, idef0, idef1, dxmin, dxmax)
                    idef1 = Tdomain%specel(n)%Iglobnum(i,j,k+1)
                    call update_minmax(Tdomain, idef0, idef1, dxmin, dxmax)
                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j+1,k)
                    call update_minmax(Tdomain, idef0, idef1, dxmin, dxmax)
                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j,k+1)
                    call update_minmax(Tdomain, idef0, idef1, dxmin, dxmax)
                    idef1 = Tdomain%specel(n)%Iglobnum(i,j+1,k+1)
                    call update_minmax(Tdomain, idef0, idef1, dxmin, dxmax)
                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j+1,k+1)
                    call update_minmax(Tdomain, idef0, idef1, dxmin, dxmax)
                enddo
            enddo
        enddo

        dxmin = sqrt(dxmin)
        mat = Tdomain%specel(n)%mat_index
        dt_loc = min(dt_loc, dxmin/Tdomain%sSubdomain(mat)%Pspeed)
        dxmax = sqrt(dxmax)
        floc_max  = min(floc_max,  Tdomain%sSubdomain(mat)%Sspeed/30/dxmax)
    enddo

    if(rg==0) write(*,*) 'dt_loc',dt_loc
    call MPI_AllReduce (dt_loc, dt_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)

    call MPI_AllReduce (floc_max, f_max, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)
    if(rg==0) write(*,*) 'f_max',f_max

    dt = courant * dt_min

    if (rg == 0 ) write (*,*) "For Courant =",courant," the time step is now ",dt_min,dt

    !! Affectation a tous les materiaux
    !!
    Tdomain%TimeD%dtmin = dt
    if (Tdomain%TimeD%dtmin > 0) then
        Tdomain%TimeD%ntimeMax = int (Tdomain%TimeD%Duration/Tdomain%TimeD%dtmin)
    else
        write (*,*) "Your dt min is zero : verify it"
        stop
    endif

    return
end subroutine Compute_Courant

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
