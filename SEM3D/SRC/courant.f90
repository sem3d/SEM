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
!! \brief V�rifie que le pas de temps est adapt� aux d�placements �valu�s
!!
!! \param type (Domain), intent (IN) Tdomain
!<
subroutine Compute_Courant (Tdomain,rg)

    use sdomain
    use mpi
    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent(IN) :: rg

    integer :: i,j,k, n, ngllx,nglly,ngllz, idef0,idef1, mat, ierr
    real :: dx,dxmin,courant,courant_max
    real :: dt_min, dt, dt_loc
    ! start modifs
    real :: floc_max, dx_max, f_max
    !end modifs

    courant = Tdomain%TimeD%courant
    courant_max = 0
    dt_loc = huge(1.)
    floc_max  = huge(1.)
    do n = 0, Tdomain%n_elem -1
        dxmin = 1e10
        ! start modifs
        dx_max= 0
        ! end modifs

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 0, ngllz-2
            do j = 0, nglly - 2
                do i = 0, ngllx -2

                    idef0 = Tdomain%specel(n)%Iglobnum(i,j,k)
                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j,k)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    ! start modifs
                    if(dx > dx_max) dx_max = dx
                    ! end modifs
                    idef1 = Tdomain%specel(n)%Iglobnum(i,j+1,k)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    ! start modifs
                    if(dx > dx_max) dx_max = dx
                    ! end modifs
                    idef1 = Tdomain%specel(n)%Iglobnum(i,j,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    ! start modifs
                    if(dx > dx_max) dx_max = dx
                    ! end modifs
                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j+1,k)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    ! start modifs
                    if(dx > dx_max) dx_max = dx
                    ! end modifs
                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    ! start modifs
                    if(dx > dx_max) dx_max = dx
                    ! end modifs
                    idef1 = Tdomain%specel(n)%Iglobnum(i,j+1,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    ! start modifs
                    if(dx > dx_max) dx_max = dx
                    ! end modifs
                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j+1,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    ! start modifs
                    if(dx > dx_max) dx_max = dx
                    ! end modifs
                enddo
            enddo
        enddo

        dxmin = sqrt(dxmin)
        mat = Tdomain%specel(n)%mat_index
        dt_loc = min(dt_loc, dxmin/Tdomain%sSubdomain(mat)%Pspeed)
        !!    courant = (Tdomain%sSubdomain(mat)%Pspeed * Tdomain%sSubdomain(mat)%Dt)/ dxmin
        !!    if (courant > courant_max) courant_max = courant

        ! start modifs
        dx_max = sqrt(dx_max)
        floc_max  = min(floc_max,  Tdomain%sSubdomain(mat)%Sspeed/30/dx_max)
        ! end modifs

    enddo

    if(rg==0) write(*,*) 'dt_loc',dt_loc
    call MPI_AllReduce (dt_loc, dt_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)

    ! start modifs
    call MPI_AllReduce (floc_max, f_max, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)
    if(rg==0) write(*,*) 'f_max',f_max
    ! end modifs

    dt = courant * dt_min

    if (rg == 0 ) write (*,*) "For Courant =",courant," the time step is now ",dt_min,dt

!!!!write (*,*) "The Courant number in your simulation is ", courant_max  , dxmin
!!!!if (courant_max  > 0.5 .and. courant_max  < 1) then
!!!!   write (*,*) " Warning: you have a Courant larger than 0.5 "
!!!!else if (courant_max  > 1 ) then
!!!!   write (*,*) " You have a Courant larger than 1, you are unstable "
!!!!   stop
!!!!endif

    !! Affectation a tous les materiaux
    !!
    do n=0,Tdomain%n_elem - 1
        mat = Tdomain%specel(n)%mat_index
        Tdomain%sSubdomain(mat)%Dt = dt
    enddo
    Tdomain%TimeD%dtmin = dt
    if (Tdomain%TimeD%dtmin > 0) then
        Tdomain%TimeD%ntimeMax = int (Tdomain%TimeD%Duration/Tdomain%TimeD%dtmin) !!oubli
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
