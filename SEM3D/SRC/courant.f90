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
subroutine Compute_Courant (Tdomain,rg)

    use sdomain
    use mpi

    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent(IN) :: rg

    integer :: i,j,k, n, ngllx,nglly,ngllz, idef0,idef1, mat, ierr
    real :: dx,dxmin, courant,courant_max
    real :: dt_min, dt, dt_loc

    !  courant = 0.05 !0.1 !0.4 !0.2 !0.5       !! 0.4 semblait ok mais uniformisation avec Sem2d
    courant = 0.2
    !    courant = 0.1
    !    test pour can1mka
    ! si le pas de temps de mka est trop petit par rapport a celui de sem
    ! il peut se produire une divergence lors d un calcul couple
    ! cas deja produit pour un rapport de 50 entre les deux pas de temps
    !           courant = 0.005
    ! pour volvi2b et volvi2c
    !       courant = 0.05
    courant_max = 0
    dt_loc = huge(1.)
    do n = 0, Tdomain%n_elem -1
        dxmin = 1e10

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 0, ngllz-2
            do j = 0, nglly - 2
                do i = 0, ngllx -2

                    idef0 = Tdomain%specel(n)%Iglobnum (i,j,k)
                    idef1 = Tdomain%specel(n)%Iglobnum (i+1,j,k)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if (dx < dxmin ) dxmin = dx
                    idef1 = Tdomain%specel(n)%Iglobnum (i,j+1,k)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if (dx < dxmin ) dxmin = dx
                    idef1 = Tdomain%specel(n)%Iglobnum (i,j,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if (dx < dxmin ) dxmin = dx

                    idef1 = Tdomain%specel(n)%Iglobnum (i+1,j+1,k)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if (dx < dxmin ) dxmin = dx
                    idef1 = Tdomain%specel(n)%Iglobnum (i+1,j,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if (dx < dxmin ) dxmin = dx
                    idef1 = Tdomain%specel(n)%Iglobnum (i,j+1,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if (dx < dxmin ) dxmin = dx

                    idef1 = Tdomain%specel(n)%Iglobnum (i+1,j+1,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if (dx < dxmin ) dxmin = dx
                enddo
            enddo
        enddo

        dxmin = sqrt (dxmin)
        mat = Tdomain%specel(n)%mat_index
        dt_loc = min(dt_loc, dxmin/Tdomain%sSubdomain(mat)%Pspeed)
        !!    courant = (Tdomain%sSubdomain(mat)%Pspeed * Tdomain%sSubdomain(mat)%Dt)/ dxmin
        !!    if (courant > courant_max) courant_max = courant
    enddo

    if(rg==0) &
        write(*,*) 'dt_loc',dt_loc
    call MPI_AllReduce (dt_loc, dt_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)

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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
