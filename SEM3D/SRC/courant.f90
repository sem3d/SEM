!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \fn subroutine Compute_Courant (Tdomain)
!! \brief Calcul le pas de temps adapte au maillage et Vs/Vp
!! \brief Vérifie que le pas de temps est adapté aux déplacements évalués
!!
!! \param type (Domain), intent (IN) Tdomain
!<
module mCourant

    implicit none
    public :: Compute_Courant
contains
    subroutine update_minmax(Tdomain, idef0, idef1, dxmin, dxmax)
        use sdomain
        implicit none
        type (Domain), intent (INOUT) :: Tdomain
        integer, intent(in) :: idef0, idef1
        real(fpp), intent(inout) :: dxmin, dxmax
        !
        real(fpp) :: dx

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
        use constants
        use dom_solid
        use dom_fluid
        use dom_solidpml
        use dom_fluidpml
        use mpi
        implicit none
        type (Domain), intent (INOUT) :: Tdomain
        integer, intent(IN) :: rg

        integer :: i,j,k,n,ngll,idef0,idef1, mat, ierr
        integer :: lnum
        real(fpp) :: dxmin,courant,courant_max
        real(fpp) :: dt_min, dt, dt_loc
        real(fpp) :: floc_max, dxmax, f_max
        real(fpp) :: pspeed, maxPspeed, sumPspeed
        real(fpp), dimension(0:Tdomain%n_mat-1) :: avgPspeed_mat, dxmin_mat, dt_loc_mat
        integer, dimension(0:Tdomain%n_mat-1) :: ngll_mat
        logical :: use_average


        courant = Tdomain%TimeD%courant
        courant_max = 0
        dt_loc = huge(1.)
        floc_max  = huge(1.)
        dxmin_mat = huge(1.)
        avgPspeed_mat = 0d0
        ngll_mat = 0
        use_average = Tdomain%use_avg
        !print*, "Tdomain%use_avg = ", Tdomain%use_avg

        do n = 0, Tdomain%n_elem -1
            dxmin = 1e10
            dxmax = 0

            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            do k = 0, ngll-2
                do j = 0, ngll-2
                    do i = 0, ngll-2

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
            ! Compute max Pspeed on the cell
            maxPspeed = 0d0
            sumPspeed = 0d0
            lnum = Tdomain%specel(n)%lnum
            do k = 0, ngll-1
                do j = 0, ngll-1
                    do i = 0, ngll-1
                        Pspeed = 0d0
                        select case(Tdomain%specel(n)%domain)
                        case (DM_SOLID)
                            Pspeed = solid_pspeed(Tdomain%sdom, lnum, i, j, k)
                        case (DM_SOLID_PML)
                            Pspeed = solidpml_pspeed(Tdomain%spmldom, lnum, i, j, k)
                        case (DM_FLUID)
                            Pspeed = fluid_pspeed(Tdomain%fdom, lnum, i, j, k)
                        case (DM_FLUID_PML)
                            Pspeed = fluidpml_pspeed(Tdomain%fpmldom, lnum, i, j, k)
                        end select
                        if (Pspeed>maxPspeed) maxPspeed = Pspeed
                        sumPspeed = sumPspeed + Pspeed
                    end do
                end do
            end do

            dxmin = sqrt(dxmin)
            mat = Tdomain%specel(n)%mat_index
            dt_loc = min(dt_loc, dxmin/maxPspeed)
            dxmax = sqrt(dxmax)
            floc_max  = min(floc_max,  Tdomain%sSubdomain(mat)%Sspeed/30/dxmax)
            ngll_mat(mat) = ngll_mat(mat) + (ngll*ngll*ngll)
            avgPspeed_mat(mat) = avgPspeed_mat(mat) + sumPspeed
            dxmin_mat(mat) = min(dxmin, dxmin_mat(mat))
        enddo

        if(use_average) then
            if(rg==0) write(*,*) 'WARNING!!using average'
            dt_loc_mat = huge(1.)
            where(ngll_mat > 0) 
                    avgPspeed_mat = avgPspeed_mat/dble(ngll_mat)
                    dt_loc_mat = dxmin_mat/avgPspeed_mat
            end where
            dt_loc = minval(dt_loc_mat)
        end if

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
end module mCourant
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
