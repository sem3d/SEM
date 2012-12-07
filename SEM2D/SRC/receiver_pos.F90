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

    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables

    integer :: i,j,n,nel,nnelem, nrec,idef, nind,ngllx, ngllz, nmin, mat, ipoint
    integer, parameter :: nimax=50,njmax=50
    integer, dimension (0:5) :: nreceiv

    real :: Dmin,Dist,x0,x1,x2,x3,z0,z1,z2,z3,a1,a2,b1,b2,c1,c2,d1,d2
    real :: alpha,beta,gamm,delta,eta1,eta2,xi1,outx,outz,dximax, detamax,xi2
    real :: x4,x5,x6,x7,z4,z5,z6,z7
    real, dimension (0:3) :: xc,zc

    logical :: inosol,inner

    ! ################################################
    ! Modified by Gaetano Festa 15/9/2004
    ! Some minor changes for compatibility with MPI 13/10/2005

    ! Compute the nearest point to the real location of the receivers in GLL scheme

    nel = Tdomain%n_elem
    do nrec = 0, Tdomain%n_receivers -1

        Dmin = 1e10

        do n = 0,Tdomain%n_glob_points-1

            Dist = (Tdomain%GlobCoord(0,n)-Tdomain%sReceiver(nrec)%Xrec)**2  +   &
                (Tdomain%GlobCoord(1,n)- Tdomain%sReceiver(nrec)%zrec)**2

            if ( Dist < Dmin ) then
                nmin = n
                Dmin = Dist
            endif
        enddo

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
            do1 : do n = 0,nind-1
                nnelem = nreceiv(n)
                ipoint = Tdomain%specel(nnelem)%Control_Nodes(0)
                xc(0)= Tdomain%Coord_nodes(0,ipoint); zc(0)= Tdomain%Coord_nodes(1,ipoint)
                ipoint = Tdomain%specel(nnelem)%Control_Nodes(1)
                xc(1)= Tdomain%Coord_nodes(0,ipoint); zc(1)= Tdomain%Coord_nodes(1,ipoint)
                ipoint = Tdomain%specel(nnelem)%Control_Nodes(2)
                xc(2)= Tdomain%Coord_nodes(0,ipoint); zc(2)= Tdomain%Coord_nodes(1,ipoint)
                ipoint = Tdomain%specel(nnelem)%Control_Nodes(3)
                xc(3)= Tdomain%Coord_nodes(0,ipoint); zc(3)= Tdomain%Coord_nodes(1,ipoint)
                call verify_in_quad (xc,zc, Tdomain%sReceiver(nrec)%Xrec,Tdomain%sReceiver(nrec)%Zrec, inner)
                if (inner ) exit do1
            enddo do1
            if (inner) then
                nreceiv(0) = nreceiv(n)
                x0= xc(0); x1= xc(1); x2= xc(2); x3= xc(3)
                z0= zc(0); z1= zc(1); z2= zc(2); z3= zc(3)

                a1 = 4 * Tdomain%sReceiver(nrec)%Xrec - x0 -x1- x2 - x3
                b1 = x0 - x1 + x3 - x2
                c1 = x0 + x1 - x2 - x3
                d1 = -x0 + x1 + x3 - x2
                a2 = 4 * Tdomain%sReceiver(nrec)%Zrec - z0 -z1- z2 - z3
                b2 = z0 - z1 + z3 - z2
                c2 = z0 + z1 - z2 - z3
                d2 = -z0 + z1 + z3 - z2
                alpha = c1*d2 - d1*c2  ; beta = a1*d2 - b1*c2 + c1*b2 - d1*a2; gamm  =a1*b2 - a2*b1
                if (alpha ==0 ) then
                    eta1 = -gamm/beta
                    if (d2 == 0 .and. b2==0) then
                        xi1 = -(a1 + c1*eta1)/(b1+d1*eta1)
                    else
                        xi1 = -(a2 + c2*eta1)/(b2+d2*eta1)
                    endif

                    inosol = xi1 <=1 .and. xi1>=-1 .and. eta1>=-1 .and. eta1<=1
                    inosol =.not. inosol
                else
                    delta = beta**2 - 4* alpha*gamm
                    if (delta < 0) then
                        inosol = .true.
                    else
                        eta1 = 0.5 * (- beta + sqrt (delta) )/ alpha
                        inosol = .true.
                        if (eta1 <= 1 .and. eta1 >=-1) then
                            xi1 = -(a2 + c2*eta1)/(b2+d2*eta1)
                            if (xi1 <=1 .and. xi1 >= -1) inosol = .false.
                        endif
                        eta2 =  0.5 * (- beta - sqrt (delta) )/ alpha
                        if (eta2 <= 1 .and. eta2 >=-1 .and. inosol) then
                            xi1 = -(a2 + c2*eta2)/(b2+d2*eta2)
                            if (xi1 <=1 .and. xi1 >= -1) inosol = .false.
                            eta1 = eta2
                        endif
                    endif
                endif
            else
                inosol = .true.
            endif
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

        else if (Tdomain%n_nodes == 8 ) then
            dximax = 2./nimax; detamax = 2./njmax
            do11_n : do n = 0,nind-1
                nnelem = nreceiv(n)
                inner = .false.
                do j = 0,njmax-1
                    do i = 0, nimax-1
                        xi1 = i*dximax -1; xi2 = xi1 + dximax
                        eta1 = j * detamax-1 ; eta2 = eta1 + detamax
                        xi1 = xi1 -dximax/nimax; xi2 = xi2 + dximax/nimax;
                        eta1 = eta1 - detamax/njmax; eta2 = eta2 + detamax/njmax

                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(0)
                        x0 = Tdomain%Coord_Nodes(0,ipoint);  z0 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(1)
                        x1 = Tdomain%Coord_Nodes(0,ipoint);  z1 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(2)
                        x2 = Tdomain%Coord_Nodes(0,ipoint);  z2 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(3)
                        x3 = Tdomain%Coord_Nodes(0,ipoint);  z3 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(4)
                        x4 = Tdomain%Coord_Nodes(0,ipoint);  z4 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(5)
                        x5 = Tdomain%Coord_Nodes(0,ipoint);  z5 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(6)
                        x6 = Tdomain%Coord_Nodes(0,ipoint);  z6 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(7)
                        x7 = Tdomain%Coord_Nodes(0,ipoint);  z7 = Tdomain%Coord_Nodes(1,ipoint)

                        xc(0) = 0.25 * ( -x0 * (1.-xi1)*(1.-eta1)*(1+xi1+eta1) - x1 * (1.+xi1)*(1.-eta1)*(1-xi1+eta1) - x2 * &
                            (1.+xi1)*(1.+eta1)*(1-xi1-eta1) - x3 * (1.-xi1)*(1.+eta1)*(1+xi1-eta1) ) + 0.5 * ( x4 * (1.-xi1**2)*(1.-eta1) + &
                            x5 * (1.+xi1)*(1.-eta1**2) +x6 * (1.-xi1**2)*(1.+eta1) + x7 * (1.-xi1)*(1.-eta1**2) )
                        zc(0) = 0.25 * ( -z0 * (1.-xi1)*(1.-eta1)*(1+xi1+eta1) - z1 * (1.+xi1)*(1.-eta1)*(1-xi1+eta1) - z2 * &
                            (1.+xi1)*(1.+eta1)*(1-xi1-eta1) -  z3 * (1.-xi1)*(1.+eta1)*(1+xi1-eta1) ) + 0.5 * ( z4 * (1.-xi1**2)*(1.-eta1) + &
                            z5 * (1.+xi1)*(1.-eta1**2) +z6 * (1.-xi1**2)*(1.+eta1) + z7 * (1.-xi1)*(1.-eta1**2) )
                        xc(1) = 0.25 * ( -x0 * (1.-xi2)*(1.-eta1)*(1+xi2+eta1) - x1 * (1.+xi2)*(1.-eta1)*(1-xi2+eta1) - x2 * &
                            (1.+xi2)*(1.+eta1)*(1-xi2-eta1) -  x3 * (1.-xi2)*(1.+eta1)*(1+xi2-eta1) ) + 0.5 * ( x4 * (1.-xi2**2)*(1.-eta1) + &
                            x5 * (1.+xi2)*(1.-eta1**2) +x6 * (1.-xi2**2)*(1.+eta1) + x7 * (1.-xi2)*(1.-eta1**2) )
                        zc(1) = 0.25 * ( -z0 * (1.-xi2)*(1.-eta1)*(1+xi2+eta1) - z1 * (1.+xi2)*(1.-eta1)*(1-xi2+eta1) - z2 * &
                            (1.+xi2)*(1.+eta1)*(1-xi2-eta1) -  z3 * (1.-xi2)*(1.+eta1)*(1+xi2-eta1) ) + 0.5 * ( z4 * (1.-xi2**2)*(1.-eta1) + &
                            z5 * (1.+xi2)*(1.-eta1**2) +z6 * (1.-xi2**2)*(1.+eta1) + z7 * (1.-xi2)*(1.-eta1**2) )
                        xc(2) = 0.25 * ( -x0 * (1.-xi2)*(1.-eta2)*(1+xi2+eta2) - x1 * (1.+xi2)*(1.-eta2)*(1-xi2+eta2) - x2 * &
                            (1.+xi2)*(1.+eta2)*(1-xi2-eta2) -  x3 * (1.-xi2)*(1.+eta2)*(1+xi2-eta2) ) + 0.5 * ( x4 * (1.-xi2**2)*(1.-eta2) + &
                            x5 * (1.+xi2)*(1.-eta2**2) +x6 * (1.-xi2**2)*(1.+eta2) + x7 * (1.-xi2)*(1.-eta2**2) )
                        zc(2) = 0.25 * ( -z0 * (1.-xi2)*(1.-eta2)*(1+xi2+eta2) - z1 * (1.+xi2)*(1.-eta2)*(1-xi2+eta2) - z2 * &
                            (1.+xi2)*(1.+eta2)*(1-xi2-eta2) - z3 * (1.-xi2)*(1.+eta2)*(1+xi2-eta2) ) + 0.5 * ( z4 * (1.-xi2**2)*(1.-eta2) + &
                            z5 * (1.+xi2)*(1.-eta2**2) +z6 * (1.-xi2**2)*(1.+eta2) + z7 * (1.-xi2)*(1.-eta2**2) )
                        xc(3) = 0.25 * ( -x0 * (1.-xi1)*(1.-eta2)*(1+xi1+eta2) - x1 * (1.+xi1)*(1.-eta2)*(1-xi1+eta2) - x2 * &
                            (1.+xi1)*(1.+eta2)*(1-xi1-eta2) -  x3 * (1.-xi1)*(1.+eta2)*(1+xi1-eta2) ) + 0.5 * ( x4 * (1.-xi1**2)*(1.-eta2) +&
                            x5 * (1.+xi1)*(1.-eta2**2) +x6 * (1.-xi1**2)*(1.+eta2) + x7 * (1.-xi1)*(1.-eta2**2) )
                        zc(3) = 0.25 * ( -z0 * (1.-xi1)*(1.-eta2)*(1+xi1+eta2) - z1 * (1.+xi1)*(1.-eta2)*(1-xi1+eta2) - z2 *  &
                            (1.+xi1)*(1.+eta2)*(1-xi1-eta2) -  z3 * (1.-xi1)*(1.+eta2)*(1+xi1-eta2) ) +  &
                            0.5 * ( z4 * (1.-xi1**2)*(1.-eta2) + z5 * (1.+xi1)*(1.-eta2**2) +z6 * (1.-xi1**2)*(1.+eta2) +  &
                            z7 * (1.-xi1)*(1.-eta2**2) )

                        call verify_in_quad (Xc,Zc,Tdomain%sReceiver(nrec)%Xrec,Tdomain%sReceiver(nrec)%Zrec, inner)
                        if (inner) exit do11_n
                    enddo
                enddo
            enddo do11_n
            if (.not. inner ) then
                Tdomain%sReceiver(nrec)%located_here = .false.
                !	  write (*,*)  "Receiver",nrec, " is not in the processor" , Tdomain%Mpi_var%my_rank

            else
                write (*,*)  "Receiver",nrec, " is in the processor" , Tdomain%Mpi_var%my_rank
                Tdomain%sReceiver(nrec)%located_here = .true.
                Tdomain%sReceiver(nrec)%nr = nreceiv (n)
                xi1 = 0.5 * (xi1+xi2); eta1 = 0.5 *( eta1+eta2)
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
