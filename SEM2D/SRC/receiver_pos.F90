!>
!!\file receiver_pos.F90
!!\brief Contient la subroutine ReceiverPosition.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

subroutine ReceiverPosition(Tdomain)

    use sdomain
    use shape_lin
    use shape_quad
    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables

    integer :: i,j,n,nel,nnelem, nrec,idef, nind,ngllx, ngllz, nmin, mat, ipoint
    integer, parameter :: nimax=50,njmax=50
    integer, dimension (0:5) :: nreceiv

    real :: Dmin,Dist
    real :: eta1,xi1
    real :: outx,outz,dximax, detamax
    real, dimension (0:7) :: xc,zc

    logical :: inosol,inner

    ! ################################################
    ! Modified by Gaetano Festa 15/9/2004
    ! Some minor changes for compatibility with MPI 13/10/2005

    ! Compute the nearest point to the real location of the receivers in GLL scheme

    nel = Tdomain%n_elem
    do nrec = 0, Tdomain%n_receivers -1

        Dmin = 1e20
        nmin = -1
        do n = 0,Tdomain%n_glob_points-1

            Dist = (Tdomain%GlobCoord(0,n)-Tdomain%sReceiver(nrec)%Xrec)**2  +   &
                (Tdomain%GlobCoord(1,n)- Tdomain%sReceiver(nrec)%zrec)**2

            if ( Dist < Dmin ) then
                nmin = n
                Dmin = Dist
            endif
        enddo
        if (nmin==-1) then
            write(*,*) "Incorrect receiver positions"
            stop
        end if

        ! Search for the element
        nind = 0
        do n = 0,nel-1
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            do j = 0,ngllz-1
                do i = 0,ngllx-1
                    idef = Tdomain%specel(n)%Iglobnum(i,j)
                    if (idef == nmin ) then
                        nreceiv(nind) = n
                        nind = nind + 1
                        if (nind > 6 ) then
                            write (*,*) "Maximum allowed connection is 6! Some elements have a larger connection"
                            write (*,*) " Return to continue, and Ctrl C to quit"
                            read  (*,*)
                        endif
                    endif
                enddo
            enddo
        enddo

        if (Tdomain%n_nodes == 4) then
            inner = .false.
            do n = 0,nind-1
                nnelem = nreceiv(n)
                do i=0, 3
                    ipoint = Tdomain%specel(nnelem)%Control_Nodes(i)
                    xc(i)= Tdomain%Coord_nodes(0,ipoint)
                    zc(i)= Tdomain%Coord_nodes(1,ipoint)
                end do
                call verify_in_quad (xc,zc, Tdomain%sReceiver(nrec)%Xrec,Tdomain%sReceiver(nrec)%Zrec, inner)
                if (inner) exit
            enddo
            if (inner) then
                nreceiv(0) = nreceiv(n)
                call shape4_local_coord(xc, zc, &
                    Tdomain%sReceiver(nrec)%Xrec, Tdomain%sReceiver(nrec)%Zrec, &
                    xi1, eta1, inosol)

                if (inosol) then
                    !              write (*,*)  "Receiver",nrec, " is not in the processor" , Tdomain%Mpi_var%my_rank
                    Tdomain%sReceiver(nrec)%located_here = .false.
                else
                    write (*,*)  "Receiver",nrec, " is in the processor" , Tdomain%Mpi_var%my_rank
                    write (*,*) "At the location (Element, xi,eta) : "
                    write (*,*) nreceiv (0), xi1, eta1
                    Tdomain%sReceiver(nrec)%located_here = .true.
                    Tdomain%sReceiver(nrec)%nr = nreceiv (0)
                    Tdomain%sReceiver(nrec)%xi = xi1
                    Tdomain%sReceiver(nrec)%eta = eta1
                endif
            end if
        else if (Tdomain%n_nodes == 8 ) then
            dximax = 2./nimax; detamax = 2./njmax
            do11_n : do n = 0,nind-1
                nnelem = nreceiv(n)
                do i=0, 7
                    ipoint = Tdomain%specel(nnelem)%Control_Nodes(i)
                    xc(i)= Tdomain%Coord_nodes(0,ipoint)
                    zc(i)= Tdomain%Coord_nodes(1,ipoint)
                end do
                call shape8_local_coord(xc, zc, &
                    Tdomain%sReceiver(nrec)%Xrec, Tdomain%sReceiver(nrec)%Zrec, &
                    xi1, eta1, inosol)
                if (.not. inosol) exit do11_n
            enddo do11_n
            if (inosol) then
                Tdomain%sReceiver(nrec)%located_here = .false.
                !	  write (*,*)  "Receiver",nrec, " is not in the processor" , Tdomain%Mpi_var%my_rank

            else
                write (*,*)  "Receiver",nrec, " is in the processor" , Tdomain%Mpi_var%my_rank
                Tdomain%sReceiver(nrec)%located_here = .true.
                Tdomain%sReceiver(nrec)%nr = nreceiv (n)
                Tdomain%sReceiver(nrec)%xi = xi1
                Tdomain%sReceiver(nrec)%eta = eta1
                write (*,*) "At the location (Element, xi,eta) : ",xi1, eta1
            endif
        endif
        if (Tdomain%sReceiver(nrec)%located_here) then
            nnelem = Tdomain%sReceiver(nrec)%nr
            mat = Tdomain%specel(nnelem)%mat_index
            ngllx = Tdomain%specel(nnelem)%ngllx
            ngllz =Tdomain%specel(nnelem)%ngllz
            allocate (Tdomain%sReceiver(nrec)%Interp_Coeff(0:ngllx-1,0:ngllz-1))

            do j = 0,ngllz-1
                call  pol_lagrange (ngllz,Tdomain%sSubdomain(mat)%GLLcz,j,Tdomain%sReceiver(nrec)%eta,outz)
                do i = 0,ngllx -1
                    call  pol_lagrange (ngllx,Tdomain%sSubdomain(mat)%GLLcx,i,Tdomain%sReceiver(nrec)%xi,outx)
                    Tdomain%sReceiver(nrec)%Interp_Coeff(i,j) = outx*outz
                enddo
            enddo
        endif
    enddo

    i = 0
    do nrec = 0, Tdomain%n_receivers-1
        if (Tdomain%sReceiver(nrec)%located_here) i= i + 1
    enddo
    if (i > 0) then
        allocate (Tdomain%Store_Trace(0:1,0:i-1,0:Tdomain%TimeD%ntimeMax-1))
    endif
    return
end subroutine ReceiverPosition
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
