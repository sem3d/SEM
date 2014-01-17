!>
!!\file source_pos.F90
!!\brief Contient la subroutine SourcePosition.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

subroutine SourcePosition(Tdomain)

    use sdomain

    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables

    integer :: nel,n, nsour,i,j,nmin,nind, mind, mat, idef, nnelem,ipoint
    integer :: ngllx,ngllz
    integer, parameter :: nimax = 50, njmax = 50
    integer, dimension (0:5) :: nsource
    real :: Dmin, Dist,a1,b1,c1,d1,a2,b2,c2,d2,alpha,beta,gamm,delta
    real :: eta1,eta2,xi1,xi2,xi,eta,wxi,weta,dwdxi,dwdeta,Jac,dximax,detamax
    real :: x0,x1,x2,x3,x4,x5,x6,x7,z0,z1,z2,z3,z4,z5,z6,z7
    real, dimension (0:3) :: xc,zc
    real, dimension (0:1,0:1) :: LocInvGrad
    real, dimension (0:5) :: xis, etas

    logical :: inner,inosol


    ! ########################################
    ! Compute the nearest point to the real location of the source in GLL scheme
    ! Search for the element
    ! Verify inside which element the source is located
    !
    ! major modification  by Gaetano Festa 06/01/2005
    ! changes for the MPI compatibility 13/10/2005
    ! ##########################################

    nel = Tdomain%n_elem

    do nsour = 0, Tdomain%n_source -1
        Dmin = 1e10
        ! Find the nearest point (gauss node) to the source
        do n = 0,Tdomain%n_glob_points-1

            Dist = (Tdomain%GlobCoord(0,n)-Tdomain%Ssource(nsour)%Xsource)**2  +   &
                (Tdomain%GlobCoord(1,n)-Tdomain%Ssource(nsour)%Zsource)**2

            if ( Dist < Dmin ) then
                nmin = n
                Dmin = Dist
            endif
        enddo

        ! Search for the elements containing this point (allow a maximum of 6)
        nind = 0
        do n = 0,nel-1
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            do j = 0,ngllz-1
                do i = 0,ngllx-1
                    idef = Tdomain%specel(n)%Iglobnum(i,j)
                    if (idef == nmin ) then
                        nsource(nind) = n
                        nind = nind + 1
                        if (nind > 6 ) then
                            write (*,*) "Maximum allowed connection is 6! Some elements have a larger connection"
                        endif
                    endif
                enddo
            enddo
        enddo

        ! Verify inside which element the source is
        if (Tdomain%n_nodes == 4) then
            mind = 0
            do n = 0, nind-1
                nnelem = nsource(n)
                ipoint = Tdomain%specel(nnelem)%Control_Nodes(0); xc(0)= Tdomain%Coord_nodes(0,ipoint); zc(0)= Tdomain%Coord_nodes(1,ipoint)
                ipoint = Tdomain%specel(nnelem)%Control_Nodes(1); xc(1)= Tdomain%Coord_nodes(0,ipoint); zc(1)= Tdomain%Coord_nodes(1,ipoint)
                ipoint = Tdomain%specel(nnelem)%Control_Nodes(2); xc(2)= Tdomain%Coord_nodes(0,ipoint); zc(2)= Tdomain%Coord_nodes(1,ipoint)
                ipoint = Tdomain%specel(nnelem)%Control_Nodes(3); xc(3)= Tdomain%Coord_nodes(0,ipoint); zc(3)= Tdomain%Coord_nodes(1,ipoint)
                call verify_in_quad (xc,zc, Tdomain%sSource(nsour)%Xsource,  Tdomain%sSource(nsour)%Zsource, inner)
                if (inner) then
                    nsource(mind) = nsource (n)
                    mind = mind + 1
                endif
            enddo
            Tdomain%sSource(nsour)%ine = mind

            if (mind /=0) then

                ! Determine xi-eta coordinates

                allocate (Tdomain%sSource(nsour)%Elem (0:mind-1)  )
                do n = 0,mind-1
                    nnelem = nsource(n)
                    ipoint = Tdomain%specel(nnelem)%Control_Nodes(0); x0= Tdomain%Coord_nodes(0,ipoint); z0= Tdomain%Coord_nodes(1,ipoint)
                    ipoint = Tdomain%specel(nnelem)%Control_Nodes(1); x1= Tdomain%Coord_nodes(0,ipoint); z1= Tdomain%Coord_nodes(1,ipoint)
                    ipoint = Tdomain%specel(nnelem)%Control_Nodes(2); x2= Tdomain%Coord_nodes(0,ipoint); z2= Tdomain%Coord_nodes(1,ipoint)
                    ipoint = Tdomain%specel(nnelem)%Control_Nodes(3); x3= Tdomain%Coord_nodes(0,ipoint); z3= Tdomain%Coord_nodes(1,ipoint)
                    a1 = 4 * Tdomain%sSource(nsour)%Xsource - x0 -x1- x2 - x3
                    b1 = x0 - x1 + x3 - x2
                    c1 = x0 + x1 - x2 - x3
                    d1 = -x0 + x1 + x3 - x2
                    a2 = 4 * Tdomain%sSource(nsour)%Zsource - z0 -z1- z2 - z3
                    b2 = z0 - z1 + z3 - z2
                    c2 = z0 + z1 - z2 - z3
                    d2 = -z0 + z1 + z3 - z2
                    alpha = c1*d2 - d1*c2  ; beta = a1*d2 - b1*c2 + c1*b2 - d1*a2; gamm  =a1*b2 - a2*b1
                    if (abs(alpha)<1e-7 ) then
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
                            write (*,*)  "No solution for the location"
                            write (*,*) " Return to continue, and Ctrl C to quit"
                            read  (*,*)
                        endif
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
                    if (inosol) then
                        write (*,*)  "No solution found for the source    ",nsour
                        stop
                    endif
                    Tdomain%sSource(nsour)%Elem(n)%nr = nsource (n)
                    Tdomain%sSource(nsour)%Elem(n)%xi = xi1
                    Tdomain%sSource(nsour)%Elem(n)%eta = eta1

                enddo


                if (Tdomain%sSource(nsour)%i_type_source == 1) then  ! Pulse directional force
                    do n = 0, Tdomain%sSource(nsour)%ine-1
                        nnelem = Tdomain%sSource(nsour)%Elem(n)%nr
                        mat = Tdomain%specel(nnelem)%mat_index
                        ngllx = Tdomain%specel(nnelem)%ngllx
                        ngllz = Tdomain%specel(nnelem)%ngllz
                        allocate  (Tdomain%sSource(nsour)%Elem(n)%ExtForce(0:ngllx-1,0:ngllz-1))
                        do j = 0,ngllz-1
                            call pol_lagrange (ngllz, Tdomain%sSubdomain(mat)%GLLcz, j, Tdomain%sSource(nsour)%Elem(n)%eta,weta)
                            do i = 0,ngllx-1
                                call pol_lagrange (ngllx, Tdomain%sSubdomain(mat)%GLLcx, i, Tdomain%sSource(nsour)%Elem(n)%xi, wxi )
                                Tdomain%sSource(nsour)%Elem(n)%ExtForce (i,j) = wxi*weta
                            enddo
                        enddo
                    enddo

                else if (Tdomain%sSource(nsour)%i_type_source  == 2) then
                    do n = 0, Tdomain%sSource(nsour)%ine-1
                        nnelem = Tdomain%sSource(nsour)%Elem(n)%nr
                        mat = Tdomain%specel(nnelem)%mat_index
                        ngllx = Tdomain%sSubdomain(mat)%ngllx
                        ngllz = Tdomain%sSubdomain(mat)%ngllz
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(0); x0= Tdomain%Coord_nodes(0,ipoint); z0= Tdomain%Coord_nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(1); x1= Tdomain%Coord_nodes(0,ipoint); z1= Tdomain%Coord_nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(2); x2= Tdomain%Coord_nodes(0,ipoint); z2= Tdomain%Coord_nodes(1,ipoint)
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(3); x3= Tdomain%Coord_nodes(0,ipoint); z3= Tdomain%Coord_nodes(1,ipoint)
                        eta = Tdomain%sSource(nsour)%Elem(n)%eta
                        xi = Tdomain%sSource(nsour)%Elem(n)%xi

                        LocInvGrad(0,0) = 0.25 * ( (x1-x0) * (1-eta) + (x2-x3) * (1+eta) )
                        LocInvGrad(1,0) = 0.25 * ( (x3-x0) * (1-xi) + (x2-x1) * (1+xi) )
                        LocInvGrad(0,1) = 0.25 * ( (z1-z0) * (1-eta) + (z2-z3) * (1+eta) )
                        LocInvGrad(1,1) = 0.25 * ( (z3-z0) * (1-xi) + (z2-z1) * (1+xi) )

                        call invert2 (LocInvGrad, Jac )
                        Tdomain%sSource(nsour)%Elem(n)%Scoeff  = LocInvGrad

                        allocate  (Tdomain%sSource(nsour)%Elem(n)%Explosion(0:ngllx-1,0:ngllz-1,0:1))
                        do j = 0,ngllz-1
                            call pol_lagrange (ngllz, Tdomain%sSubdomain(mat)%GLLcz, j, eta,weta)
                            call DERIVLAG (Tdomain%sSubdomain(mat)%GLLcz, ngllz, j, eta, dwdeta)
                            do i = 0,ngllx-1
                                call pol_lagrange (ngllx, Tdomain%sSubdomain(mat)%GLLcx, i, xi, wxi )
                                call DERIVLAG ( Tdomain%sSubdomain(mat)%GLLcx, ngllx, i, xi, dwdxi)
                                Tdomain%sSource(nsour)%Elem(n)%Explosion (i,j,0) = Tdomain%sSource(nsour)%Elem(n)%Scoeff(0,0) * dwdxi * weta + &
                                    Tdomain%sSource(nsour)%Elem(n)%Scoeff(0,1) * dwdeta * wxi
                                Tdomain%sSource(nsour)%Elem(n)%Explosion (i,j,1) =Tdomain%sSource(nsour)%Elem(n)%Scoeff(1,0) * dwdxi * weta + &
                                    Tdomain%sSource(nsour)%Elem(n)%Scoeff(1,1) * dwdeta * wxi
                            enddo
                        enddo
                    enddo
                endif

            endif

        else if (Tdomain%n_nodes == 8 ) then
            mind = 0
            dximax = 2./nimax; detamax = 2./njmax
            do n = 0, nind-1
                nnelem = nsource(n)
                inner = .false.
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
                do12_jmax : do j = 0,njmax-1
                    do i = 0, nimax-1
                        xi1 = i*dximax -1; xi2 = xi1 + dximax
                        eta1 = j * detamax-1 ; eta2 = eta1 + detamax
                        xi1 = xi1 -dximax/nimax; xi2 = xi2 + dximax/nimax;
                        eta1 = eta1 - detamax/njmax; eta2 = eta2 + detamax/njmax

                        xc(0) = 0.25 * ( -x0 * (1.-xi1)*(1.-eta1)*(1+xi1+eta1) - x1 * (1.+xi1)*(1.-eta1)*(1-xi1+eta1) - x2 *  &
                            (1.+xi1)*(1.+eta1)*(1-xi1-eta1) -  x3 * (1.-xi1)*(1.+eta1)*(1+xi1-eta1) ) + 0.5 * ( x4 * (1.-xi1**2)*(1.-eta1) + &
                            x5 * (1.+xi1)*(1.-eta1**2) +x6 * (1.-xi1**2)*(1.+eta1) + x7 * (1.-xi1)*(1.-eta1**2) )
                        zc(0) = 0.25 * ( -z0 * (1.-xi1)*(1.-eta1)*(1+xi1+eta1) - z1 * (1.+xi1)*(1.-eta1)*(1-xi1+eta1) - z2 *   &
                            (1.+xi1)*(1.+eta1)*(1-xi1-eta1) -  z3 * (1.-xi1)*(1.+eta1)*(1+xi1-eta1) ) + 0.5 * ( z4 * (1.-xi1**2)*(1.-eta1) +  &
                            z5 * (1.+xi1)*(1.-eta1**2) +z6 * (1.-xi1**2)*(1.+eta1) + z7 * (1.-xi1)*(1.-eta1**2) )
                        xc(1) = 0.25 * ( -x0 * (1.-xi2)*(1.-eta1)*(1+xi2+eta1) - x1 * (1.+xi2)*(1.-eta1)*(1-xi2+eta1) - x2 * &
                            (1.+xi2)*(1.+eta1)*(1-xi2-eta1) -  x3 * (1.-xi2)*(1.+eta1)*(1+xi2-eta1) ) + 0.5 * ( x4 * (1.-xi2**2)*(1.-eta1) + &
                            x5 * (1.+xi2)*(1.-eta1**2) +x6 * (1.-xi2**2)*(1.+eta1) + x7 * (1.-xi2)*(1.-eta1**2) )
                        zc(1) = 0.25 * ( -z0 * (1.-xi2)*(1.-eta1)*(1+xi2+eta1) - z1 * (1.+xi2)*(1.-eta1)*(1-xi2+eta1) - z2 *  &
                            (1.+xi2)*(1.+eta1)*(1-xi2-eta1) -  z3 * (1.-xi2)*(1.+eta1)*(1+xi2-eta1) ) + 0.5 * ( z4 * (1.-xi2**2)*(1.-eta1) + &
                            z5 * (1.+xi2)*(1.-eta1**2) +z6 * (1.-xi2**2)*(1.+eta1) + z7 * (1.-xi2)*(1.-eta1**2) )
                        xc(2) = 0.25 * ( -x0 * (1.-xi2)*(1.-eta2)*(1+xi2+eta2) - x1 * (1.+xi2)*(1.-eta2)*(1-xi2+eta2) - x2 * &
                            (1.+xi2)*(1.+eta2)*(1-xi2-eta2) -  x3 * (1.-xi2)*(1.+eta2)*(1+xi2-eta2) ) + 0.5 * ( x4 * (1.-xi2**2)*(1.-eta2) +&
                            x5 * (1.+xi2)*(1.-eta2**2) +x6 * (1.-xi2**2)*(1.+eta2) + x7 * (1.-xi2)*(1.-eta2**2) )
                        zc(2) = 0.25 * ( -z0 * (1.-xi2)*(1.-eta2)*(1+xi2+eta2) - z1 * (1.+xi2)*(1.-eta2)*(1-xi2+eta2) - z2 * &
                            (1.+xi2)*(1.+eta2)*(1-xi2-eta2) -  z3 * (1.-xi2)*(1.+eta2)*(1+xi2-eta2) ) + 0.5 * ( z4 * (1.-xi2**2)*(1.-eta2) +&
                            z5 * (1.+xi2)*(1.-eta2**2) +z6 * (1.-xi2**2)*(1.+eta2) + z7 * (1.-xi2)*(1.-eta2**2) )
                        xc(3) = 0.25 * ( -x0 * (1.-xi1)*(1.-eta2)*(1+xi1+eta2) - x1 * (1.+xi1)*(1.-eta2)*(1-xi1+eta2) - x2 * &
                            (1.+xi1)*(1.+eta2)*(1-xi1-eta2) -  x3 * (1.-xi1)*(1.+eta2)*(1+xi1-eta2) ) + 0.5 * ( x4 * (1.-xi1**2)*(1.-eta2) + &
                            x5 * (1.+xi1)*(1.-eta2**2) +x6 * (1.-xi1**2)*(1.+eta2) + x7 * (1.-xi1)*(1.-eta2**2) )
                        zc(3) = 0.25 * ( -z0 * (1.-xi1)*(1.-eta2)*(1+xi1+eta2) - z1 * (1.+xi1)*(1.-eta2)*(1-xi1+eta2) - z2 * &
                            (1.+xi1)*(1.+eta2)*(1-xi1-eta2) -  z3 * (1.-xi1)*(1.+eta2)*(1+xi1-eta2) ) + 0.5 * ( z4 * (1.-xi1**2)*(1.-eta2) + &
                            z5 * (1.+xi1)*(1.-eta2**2) +z6 * (1.-xi1**2)*(1.+eta2) + z7 * (1.-xi1)*(1.-eta2**2) )

                        call verify_in_quad (Xc,Zc,Tdomain%sSource(nsour)%Xsource,Tdomain%sSource(nsour)%Zsource, inner)
                        !print *, xc, zc, Tdomain%MPI_var%My_rank
                        if (inner) exit do12_jmax
                    enddo
                enddo do12_jmax
                if (inner) then
                    nsource (mind) = nsource (n)
                    xis(mind) = 0.5 * (xi1+xi2); etas(mind) = 0.5 * (eta1+eta2)
                    mind = mind +1
                endif
            enddo
            Tdomain%sSource(nsour)%ine = mind;
            if ( mind /= 0 ) then

                allocate ( Tdomain%sSource(nsour)%Elem(0:mind-1) )
                Tdomain%sSource(nsour)%Elem(0:mind-1)%nr = nsource(0:mind-1)
                Tdomain%sSource(nsour)%Elem(0:mind-1)%xi = xis(0:mind-1)
                Tdomain%sSource(nsour)%Elem(0:mind-1)%eta = etas(0:mind-1)

                if (Tdomain%sSource(nsour)%i_type_source == 1) then  ! Pulse directional force
                    do n = 0, Tdomain%sSource(nsour)%ine-1
                        nnelem = Tdomain%sSource(nsour)%Elem(n)%nr
                        mat = Tdomain%specel(nnelem)%mat_index
                        ngllx = Tdomain%specel(nnelem)%ngllx
                        ngllz = Tdomain%specel(nnelem)%ngllz
                        allocate  (Tdomain%sSource(nsour)%Elem(n)%ExtForce(0:ngllx-1,0:ngllz-1))
                        do j = 0,ngllz-1
                            call pol_lagrange (ngllz, Tdomain%sSubdomain(mat)%GLLcz, j, Tdomain%sSource(nsour)%Elem(n)%eta,weta)
                            do i = 0,ngllx-1
                                call pol_lagrange (ngllx, Tdomain%sSubdomain(mat)%GLLcx, i, Tdomain%sSource(nsour)%Elem(n)%xi, wxi )
                                Tdomain%sSource(nsour)%Elem(n)%ExtForce (i,j) = wxi*weta
                            enddo
                        enddo
                    enddo



                else if ( Tdomain%sSource(nsour)%i_type_source == 2) then   ! Explosive source diagonal moment considered
                    do n = 0, Tdomain%sSource(nsour)%ine-1
                        nnelem = Tdomain%sSource(nsour)%Elem(n)%nr
                        mat = Tdomain%specel(nnelem)%mat_index
                        ngllx = Tdomain%specel(nnelem)%ngllx
                        ngllz = Tdomain%specel(nnelem)%ngllz

                        ipoint = Tdomain%specel(n)%Control_Nodes(0)
                        x0 = Tdomain%Coord_Nodes(0,ipoint);  z0 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(n)%Control_Nodes(1)
                        x1 = Tdomain%Coord_Nodes(0,ipoint);  z1 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(n)%Control_Nodes(2)
                        x2 = Tdomain%Coord_Nodes(0,ipoint);  z2 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(n)%Control_Nodes(3)
                        x3 = Tdomain%Coord_Nodes(0,ipoint);  z3 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(n)%Control_Nodes(4)
                        x4 = Tdomain%Coord_Nodes(0,ipoint);  z4 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(n)%Control_Nodes(5)
                        x5 = Tdomain%Coord_Nodes(0,ipoint);  z5 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(n)%Control_Nodes(6)
                        x6 = Tdomain%Coord_Nodes(0,ipoint);  z6 = Tdomain%Coord_Nodes(1,ipoint)
                        ipoint = Tdomain%specel(n)%Control_Nodes(7)
                        x7 = Tdomain%Coord_Nodes(0,ipoint);  z7 = Tdomain%Coord_Nodes(1,ipoint)

                        eta = Tdomain%sSource(nsour)%Elem(n)%eta
                        xi = Tdomain%sSource(nsour)%Elem(n)%xi

                        LocInvGrad(0,0) = 0.25 * (x0 *(1-eta)*(2*xi+eta) + x1 *(1-eta)*(2*xi-eta) + x2 *(1+eta)*(2*xi+eta)+ &
                            x3 *(1+eta)*(2*xi-eta)) - x4 * xi*(1-eta) - x6 *xi *(1+eta) + 0.5* (x5-x7)* (1-eta**2)
                        LocInvGrad(1,0) = 0.25 * (x0 *(1-xi)*(2*eta+xi) - x1 *(1+xi)*(xi-2*eta) + x2 *(1+xi)*(2*eta+xi)-  &
                            x3 *(1-xi)*(xi-2*eta)) - x5 *eta*(1+xi) - x7 *eta *(1-xi) + 0.5* (x6-x4)* (1-xi**2)
                        LocInvGrad(0,1) = 0.25 * (z0 *(1-eta)*(2*xi+eta) + z1 *(1-eta)*(2*xi-eta) + z2 *(1+eta)*(2*xi+eta)+ &
                            z3 *(1+eta)*(2*xi-eta)) - z4 * xi*(1-eta) - z6 *xi *(1+eta) + 0.5* (z5-z7)* (1-eta**2)
                        LocInvGrad(1,1) = 0.25 * (z0 *(1-xi)*(2*eta+xi) - z1 *(1+xi)*(xi-2*eta) + z2 *(1+xi)*(2*eta+xi)-  &
                            z3 *(1-xi)*(xi-2*eta)) - z5 *eta*(1+xi) - z7 *eta *(1-xi) + 0.5* (z6-z4)* (1-xi**2)
                        call invert2 (LocInvGrad, Jac)
                        Tdomain%sSource(nsour)%Elem(n)%Scoeff  = LocInvGrad
                        allocate  (Tdomain%sSource(nsour)%Elem(n)%Explosion(0:ngllx-1,0:ngllz-1,0:1))
                        do j = 0,ngllz-1
                            call pol_lagrange (ngllz, Tdomain%sSubdomain(mat)%GLLcz, j, eta,weta)
                            call DERIVLAG (Tdomain%sSubdomain(mat)%GLLcz, ngllz, j, eta, dwdeta)
                            do i = 0,ngllx-1
                                call pol_lagrange (ngllx, Tdomain%sSubdomain(mat)%GLLcx, i, xi, wxi )
                                call DERIVLAG ( Tdomain%sSubdomain(mat)%GLLcx, ngllx, i, xi, dwdxi)
                                Tdomain%sSource(nsour)%Elem(n)%Explosion (i,j,0) = Tdomain%sSource(nsour)%Elem(n)%Scoeff(0,0) * dwdxi * weta + &
                                    Tdomain%sSource(nsour)%Elem(n)%Scoeff(0,1) * dwdeta * wxi
                                Tdomain%sSource(nsour)%Elem(n)%Explosion (i,j,1) =Tdomain%sSource(nsour)%Elem(n)%Scoeff(1,0) * dwdxi * weta + &
                                    Tdomain%sSource(nsour)%Elem(n)%Scoeff(1,1) * dwdeta * wxi
                            enddo
                        enddo
                    enddo
                endif
            endif
        endif   ! if on nnods

        ! if (Tdomain%sSource(nsour)%ine == 0) then
        !       write (*,*) "No source found in the processor " ,Tdomain%Mpi_var%my_rank
        !else
        if (Tdomain%sSource(nsour)%ine /= 0) then
            write (*,*) "Source ", nsour , " is in the processor ", Tdomain%Mpi_var%my_rank
            write (*,*) " Location of the source is (Element, xi,eta) : "
            do n = 0,Tdomain%sSource(nsour)%ine -1
                write (*,*) Tdomain%sSource(nsour)%Elem(n)%nr, Tdomain%sSource(nsour)%Elem(n)%xi, Tdomain%sSource(nsour)%Elem(n)%eta
            enddo
        endif
    enddo ! number of sources
    return
end subroutine SourcePosition
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
