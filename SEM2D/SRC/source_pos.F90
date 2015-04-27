!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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
    use shape_lin
    use shape_quad
    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables

    integer :: nel,n, nsour,i,j,nmin,nind, mind, idef, nnelem,ipoint
    integer :: ngllx,ngllz
    integer, dimension (0:5) :: nsource
    real :: Dmin, Dist, eta1, xi1
    real, dimension (0:7) :: xc,zc
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
        nmin = 0
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
                do i=0, 3
                    ipoint = Tdomain%specel(nnelem)%Control_Nodes(i)
                    xc(i)= Tdomain%Coord_nodes(0,ipoint)
                    zc(i)= Tdomain%Coord_nodes(1,ipoint)
                end do
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
                    do i=0, 3
                        ipoint = Tdomain%specel(nnelem)%Control_Nodes(i)
                        xc(i)= Tdomain%Coord_nodes(0,ipoint)
                        zc(i)= Tdomain%Coord_nodes(1,ipoint)
                    end do
                    call shape4_local_coord(xc, zc, Tdomain%sSource(nsour)%Xsource, &
                        Tdomain%sSource(nsour)%Zsource, xi1, eta1, inosol)

                    if (inosol) then
                        write(*,*) "Source not found"
                        stop
                    else
                        Tdomain%sSource(nsour)%Elem(n)%nr = nsource (n)
                        Tdomain%sSource(nsour)%Elem(n)%xi = xi1
                        Tdomain%sSource(nsour)%Elem(n)%eta = eta1
                    end if
                enddo


                if (Tdomain%sSource(nsour)%i_type_source == 1) then  ! Pulse directional force
                    call source_excit_pulse(Tdomain, Tdomain%sSource(nsour))
                else if (Tdomain%sSource(nsour)%i_type_source  == 2) then
                    call calc_shape4_coeffs(Tdomain, Tdomain%sSource(nsour))
                    call source_excit_moment(Tdomain, Tdomain%sSource(nsour))
                endif
            endif

        else if (Tdomain%n_nodes == 8 ) then
            mind = 0
            do n = 0, nind-1
                nnelem = nsource(n)
                inner = .false.
                do i=0, 7
                    ipoint = Tdomain%specel(nnelem)%Control_Nodes(i)
                    xc(i)= Tdomain%Coord_nodes(0,ipoint)
                    zc(i)= Tdomain%Coord_nodes(1,ipoint)
                end do
                call shape8_local_coord(xc, zc, Tdomain%sSource(nsour)%Xsource, &
                    Tdomain%sSource(nsour)%Zsource, xi1, eta1, inosol)
                if (inosol) then
                    write(*,*) "Source not found"
                    stop
                else
                    nsource (mind) = nsource (n)
                    xis(mind) = xi1
                    etas(mind) = eta1
                    mind = mind +1
                endif
            enddo
            Tdomain%sSource(nsour)%ine = mind;
            if ( mind /= 0 ) then

                allocate ( Tdomain%sSource(nsour)%Elem(0:mind-1) )
                Tdomain%sSource(nsour)%Elem(0:mind-1)%nr = nsource(0:mind-1)
                Tdomain%sSource(nsour)%Elem(0:mind-1)%xi = xis(0:mind-1)
                Tdomain%sSource(nsour)%Elem(0:mind-1)%eta = etas(0:mind-1)

                if (Tdomain%sSource(nsour)%i_type_source == 1) then
                    ! Pulse directional force
                    call source_excit_pulse(Tdomain, Tdomain%sSource(nsour))
                else if ( Tdomain%sSource(nsour)%i_type_source == 2) then
                    ! Explosive source diagonal moment considered
                    call calc_shape8_coeffs(Tdomain, Tdomain%sSource(nsour))
                    call source_excit_moment(Tdomain, Tdomain%sSource(nsour))
                endif
            endif
        endif   ! if on nnods

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

subroutine source_excit_pulse(Tdomain, src)
    use sdomain
    use ssources
    use ssubdomains
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(Subdomain), pointer :: mat
    type(Source), intent(inout) :: src
    integer :: n, i, j, ngllx, ngllz, nmat, nnelem
    real :: weta, wxi

    do n = 0, src%ine-1
        nnelem = src%Elem(n)%nr
        nmat = Tdomain%specel(nnelem)%mat_index
        mat => Tdomain%sSubdomain(nmat)
        ngllx = mat%ngllx
        ngllz = mat%ngllz
        allocate  (src%Elem(n)%ExtForce(0:ngllx-1,0:ngllz-1,0:1))
        do j = 0,ngllz-1
            call pol_lagrange (ngllz, mat%GLLcz, j, src%Elem(n)%eta,weta)
            do i = 0,ngllx-1
                call pol_lagrange (ngllx, mat%GLLcx, i, src%Elem(n)%xi, wxi )
                src%Elem(n)%ExtForce (i,j,0) = wxi*weta*src%dir(1)
                src%Elem(n)%ExtForce (i,j,1) = wxi*weta*src%dir(2)
            enddo
        enddo
    enddo

end subroutine source_excit_pulse

subroutine source_excit_moment(Tdomain, src)
    use sdomain
    use ssources
    use ssubdomains
    use selement
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(Subdomain), pointer :: mat
    type(element), pointer :: elem
    type(Source), intent(inout) :: src
    integer :: n, i, j, ngllx, ngllz, nmat, nnelem
    real, dimension(0:1, 0:1) :: InvGrad, M
    real :: xi, eta, wxi, weta, dwdxi, dwdeta

    M = src%moment
    do n = 0, src%ine-1
        nnelem = src%Elem(n)%nr
        elem => Tdomain%specel(nnelem)
        nmat = elem%mat_index
        mat => Tdomain%sSubdomain(nmat)
        ngllx = mat%ngllx
        ngllz = mat%ngllz

        InvGrad = src%Elem(n)%Scoeff
        eta = src%Elem(n)%eta
        xi = src%Elem(n)%xi

        allocate  (src%Elem(n)%ExtForce(0:ngllx-1,0:ngllz-1,0:1))
        do j = 0,ngllz-1
            call pol_lagrange (ngllz, mat%GLLcz, j, eta,weta)
            call DERIVLAG (mat%GLLcz, ngllz, j, eta, dwdeta)
            do i = 0,ngllx-1
                call pol_lagrange (ngllx, mat%GLLcx, i, xi, wxi )
                call DERIVLAG ( mat%GLLcx, ngllx, i, xi, dwdxi)
                src%Elem(n)%ExtForce (i,j,0) = &
                    (InvGrad(0,0)*dwdxi*weta + InvGrad(0,1)*dwdeta*wxi)*M(0,0) + &
                    (InvGrad(1,0)*dwdxi*weta + InvGrad(1,1)*dwdeta*wxi)*M(0,1)
                src%Elem(n)%ExtForce (i,j,1) = &
                    (InvGrad(0,0)*dwdxi*weta + InvGrad(0,1)*dwdeta*wxi)*M(1,0) + &
                    (InvGrad(1,0)*dwdxi*weta + InvGrad(1,1)*dwdeta*wxi)*M(1,1)
            enddo
        enddo
    enddo
end subroutine source_excit_moment

subroutine calc_shape4_coeffs(Tdomain, src)
    use sdomain
    use ssources
    use ssubdomains
    use selement
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(element), pointer :: elem
    type(Source), intent(inout) :: src
    real, dimension(0:1, 0:1) :: InvGrad
    integer :: n, ipoint, nnelem
    real :: xi, eta, jac
    real :: x0, x1, x2, x3
    real :: z0, z1, z2, z3

    do n = 0, src%ine-1
        nnelem = src%Elem(n)%nr
        elem => Tdomain%specel(nnelem)

        ipoint = elem%Control_Nodes(0); x0= Tdomain%Coord_nodes(0,ipoint); z0= Tdomain%Coord_nodes(1,ipoint)
        ipoint = elem%Control_Nodes(1); x1= Tdomain%Coord_nodes(0,ipoint); z1= Tdomain%Coord_nodes(1,ipoint)
        ipoint = elem%Control_Nodes(2); x2= Tdomain%Coord_nodes(0,ipoint); z2= Tdomain%Coord_nodes(1,ipoint)
        ipoint = elem%Control_Nodes(3); x3= Tdomain%Coord_nodes(0,ipoint); z3= Tdomain%Coord_nodes(1,ipoint)
        eta = src%Elem(n)%eta
        xi = src%Elem(n)%xi

        InvGrad(0,0) = 0.25 * ( (x1-x0) * (1-eta) + (x2-x3) * (1+eta) )
        InvGrad(1,0) = 0.25 * ( (x3-x0) * (1-xi) + (x2-x1) * (1+xi) )
        InvGrad(0,1) = 0.25 * ( (z1-z0) * (1-eta) + (z2-z3) * (1+eta) )
        InvGrad(1,1) = 0.25 * ( (z3-z0) * (1-xi) + (z2-z1) * (1+xi) )

        call invert2 (InvGrad, Jac )
        src%Elem(n)%Scoeff = InvGrad
    enddo
end subroutine calc_shape4_coeffs

subroutine calc_shape8_coeffs(Tdomain, src)
    use sdomain
    use ssources
    use ssubdomains
    use selement
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(element), pointer :: elem
    type(Source), intent(inout) :: src
    real, dimension(0:1, 0:1) :: InvGrad
    integer :: n, ipoint, nnelem
    real :: xi, eta, jac
    real :: x0, x1, x2, x3, x4, x5, x6, x7
    real :: z0, z1, z2, z3, z4, z5, z6, z7

    do n = 0, src%ine-1
        nnelem = src%Elem(n)%nr
        elem => Tdomain%specel(nnelem)

        ipoint = elem%Control_Nodes(0)
        x0 = Tdomain%Coord_Nodes(0,ipoint);  z0 = Tdomain%Coord_Nodes(1,ipoint)
        ipoint = elem%Control_Nodes(1)
        x1 = Tdomain%Coord_Nodes(0,ipoint);  z1 = Tdomain%Coord_Nodes(1,ipoint)
        ipoint = elem%Control_Nodes(2)
        x2 = Tdomain%Coord_Nodes(0,ipoint);  z2 = Tdomain%Coord_Nodes(1,ipoint)
        ipoint = elem%Control_Nodes(3)
        x3 = Tdomain%Coord_Nodes(0,ipoint);  z3 = Tdomain%Coord_Nodes(1,ipoint)
        ipoint = elem%Control_Nodes(4)
        x4 = Tdomain%Coord_Nodes(0,ipoint);  z4 = Tdomain%Coord_Nodes(1,ipoint)
        ipoint = elem%Control_Nodes(5)
        x5 = Tdomain%Coord_Nodes(0,ipoint);  z5 = Tdomain%Coord_Nodes(1,ipoint)
        ipoint = elem%Control_Nodes(6)
        x6 = Tdomain%Coord_Nodes(0,ipoint);  z6 = Tdomain%Coord_Nodes(1,ipoint)
        ipoint = elem%Control_Nodes(7)
        x7 = Tdomain%Coord_Nodes(0,ipoint);  z7 = Tdomain%Coord_Nodes(1,ipoint)

        eta = src%Elem(n)%eta
        xi = src%Elem(n)%xi

        InvGrad(0,0) = 0.25 * (x0 *(1-eta)*(2*xi+eta) + x1 *(1-eta)*(2*xi-eta) + x2 *(1+eta)*(2*xi+eta)+ &
            x3 *(1+eta)*(2*xi-eta)) - x4 * xi*(1-eta) - x6 *xi *(1+eta) + 0.5* (x5-x7)* (1-eta**2)
        InvGrad(1,0) = 0.25 * (x0 *(1-xi)*(2*eta+xi) - x1 *(1+xi)*(xi-2*eta) + x2 *(1+xi)*(2*eta+xi)-  &
            x3 *(1-xi)*(xi-2*eta)) - x5 *eta*(1+xi) - x7 *eta *(1-xi) + 0.5* (x6-x4)* (1-xi**2)
        InvGrad(0,1) = 0.25 * (z0 *(1-eta)*(2*xi+eta) + z1 *(1-eta)*(2*xi-eta) + z2 *(1+eta)*(2*xi+eta)+ &
            z3 *(1+eta)*(2*xi-eta)) - z4 * xi*(1-eta) - z6 *xi *(1+eta) + 0.5* (z5-z7)* (1-eta**2)
        InvGrad(1,1) = 0.25 * (z0 *(1-xi)*(2*eta+xi) - z1 *(1+xi)*(xi-2*eta) + z2 *(1+xi)*(2*eta+xi)-  &
            z3 *(1-xi)*(xi-2*eta)) - z5 *eta*(1+xi) - z7 *eta *(1-xi) + 0.5* (z6-z4)* (1-xi**2)
        call invert2 (InvGrad, Jac)

        src%Elem(n)%Scoeff = InvGrad
    enddo
end subroutine calc_shape8_coeffs




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
