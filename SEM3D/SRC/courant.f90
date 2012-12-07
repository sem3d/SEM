subroutine Compute_Courant(Tdomain,rank)

    use sdomain
    implicit none

    type(domain), intent(in) :: Tdomain
    integer, intent(in)  :: rank
    integer :: i,j,k,n,ngllx,nglly,ngllz,idef0,idef1,mat
    real :: dx,dxmin,courant,courant_max


    courant_max = 0
    do n = 0, Tdomain%n_elem -1
        dxmin = 1e10

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
                    idef1 = Tdomain%specel(n)%Iglobnum(i,j+1,k)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    idef1 = Tdomain%specel(n)%Iglobnum(i,j,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx

                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j+1,k)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                    idef1 = Tdomain%specel(n)%Iglobnum(i,j+1,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx

                    idef1 = Tdomain%specel(n)%Iglobnum(i+1,j+1,k+1)
                    dx = (Tdomain%Globcoord(0,idef0) - Tdomain%Globcoord(0,idef1))**2 + &
                        (Tdomain%Globcoord(1,idef0) - Tdomain%Globcoord(1,idef1))**2 + &
                        (Tdomain%Globcoord(2,idef0) - Tdomain%Globcoord(2,idef1))**2
                    if(dx < dxmin) dxmin = dx
                enddo
            enddo
        enddo

        dxmin = sqrt(dxmin)
        mat = Tdomain%specel(n)%mat_index
        courant = (Tdomain%sSubdomain(mat)%Pspeed*Tdomain%sSubdomain(mat)%Dt)/dxmin
        if(courant > courant_max) courant_max = courant
    enddo

    write(*,*) "Proc #",rank
    write(*,*) "The Courant number of the simulation is: ",courant_max
    if(courant_max > 0.5 .and. courant_max < 1) then
        write(*,*) " Warning: you have a Courant larger than 0.5 "
    else if(courant_max > 1)then
        write(*,*) " You have a Courant larger than 1, instability"
        stop
    endif

    return
end subroutine Compute_Courant
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
