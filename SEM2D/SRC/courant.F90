!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file courant.F90
!!\brief Contient la subroutine Compute_Courant().
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Vérifie que le pas de temps est adapté
!!
!! Calcul du nombre de Courant
!! \param type (Domain), intent (IN) Tdomain
!<


subroutine Compute_Courant (Tdomain)
    use sdomain
    use mpi
    implicit none

    type (Domain), intent (INOUT) :: Tdomain

    integer :: i, j, ki, kj, n,ngllx, ngllz, idef0, idef1,mat, ierr
    real :: dx,dxmin, courant, courant_max
    real :: dt_min, dt, dt_loc

    courant = Tdomain%TimeD%courant

    courant_max = 0
    dt_loc = huge(1.)
    do n = 0, Tdomain%n_elem - 1
        dxmin = 1e10
        ngllx =  Tdomain%specel(n)%ngllx
        ngllz =  Tdomain%specel(n)%ngllz
        do j = 0, ngllz -1
            do i = 0, ngllx-1
                idef0 = Tdomain%specel(n)%Iglobnum(i,j)
                if (i < ngllx-1) then
                    do ki = i+1, ngllx-1
                        idef1 = Tdomain%specel(n)%Iglobnum(ki,j)
                        dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                            (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2
                        if (dx < dxmin ) then
                            dxmin = dx
                        endif
                    enddo
                endif
                if (j < ngllz-1) then
                    do kj  = j+1, ngllz-1
                        do ki  = 0, ngllx-1
                            idef1 = Tdomain%specel(n)%Iglobnum(ki,kj)
                            dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                                (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2
                            if (dx < dxmin ) then
                                dxmin = dx
                            endif
                            if(abs(dx)<1.e-5) then
                                print *,'idef',idef0,idef1,Tdomain%Globcoord(0:1,idef0)
                            endif
                        enddo
                    enddo
                endif
            enddo
        enddo
        dxmin = sqrt (dxmin)
        mat = Tdomain%specel(n)%mat_index
        dt_loc = min(dt_loc, dxmin/Tdomain%sSubdomain(mat)%Pspeed)   !! ajout calcul pas de temps local
    enddo


    print*,'dt_loc',dt_loc,dxmin
    call MPI_AllReduce (dt_loc, dt_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)

    dt = courant * dt_min

    if (Tdomain%Mpi_Var%my_rank == 0 ) write (*,*) "For Courant =",courant," the time step is now ",dt

    !! Affectation a tous les materiaux             !!! ajout 01/10
    !!
    do mat=0,Tdomain%n_mat - 1
        Tdomain%sSubdomain(mat)%Dt = dt
    enddo
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
