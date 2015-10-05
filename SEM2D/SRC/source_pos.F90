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
    use mlocations2d
    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables

    integer, parameter :: NMAXEL=20
    integer, dimension(NMAXEL) :: elems
    double precision, dimension(0:1,NMAXEL) :: coordloc
    double precision, parameter :: EPS = 1D-13
    integer :: nsour, nmax, i, n_el
    double precision :: xc, zc, xi, eta
    logical :: inside

    ! Find the nearest point to the real location of the sources in GLL scheme

    do nsour = 0, Tdomain%n_source -1
        ! Find source location

        xc = Tdomain%sSource(nsour)%xsource
        zc = Tdomain%sSource(nsour)%zsource
        nmax = NMAXEL
        call find_location(Tdomain, xc, zc, nmax, elems, coordloc)

        Tdomain%sSource(nsour)%located_here = .false.
        do i=1,nmax! When the source is in the mesh
            inside = .true.
            n_el = elems(i)
            xi   = coordloc(0,i)
            eta  = coordloc(1,i)
            if (xi<(-1-EPS) .or. eta<(-1-EPS)) inside = .false.
            if (xi>( 1+EPS) .or. eta>( 1+EPS)) inside = .false.
            if (inside) then
                write (*,'(a,i5,a,i4,a,f10.5,a,f10.5,a,i10,a,f20.17,a,f20.17,a)') " Source ", nsour, " : found on proc. ", Tdomain%Mpi_var%my_rank,     &
                                                                                  ", (xc, zc) : (", xc, ", ", zc, ") <=> (element, xi, eta) : (", n_el, &
                                                                                  ", ", xi, ", ", eta, ")"
                Tdomain%sSource(nsour)%located_here = .true.
                Tdomain%sSource(nsour)%nr           = n_el
                Tdomain%sSource(nsour)%xi           = xi
                Tdomain%sSource(nsour)%eta          = eta

                ! Customizing according to source type

                if (Tdomain%sSource(nsour)%i_type_source == 1) then       ! Pulse directional force
                    call source_excit_pulse(Tdomain, Tdomain%sSource(nsour))
                else if (Tdomain%sSource(nsour)%i_type_source  == 2) then ! Explosive source diagonal moment considered
                    if (Tdomain%n_nodes == 4) call calc_shape4_coeffs(Tdomain, Tdomain%sSource(nsour))
                    if (Tdomain%n_nodes == 8) call calc_shape8_coeffs(Tdomain, Tdomain%sSource(nsour))
                    call source_excit_moment(Tdomain, Tdomain%sSource(nsour))
                endif
            end if
        end do
    enddo
end subroutine SourcePosition

subroutine source_excit_pulse(Tdomain, src)
    use sdomain
    use ssources
    use ssubdomains
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(Subdomain), pointer :: mat
    type(Source), intent(inout) :: src
    integer :: i, j, ngllx, ngllz, nmat, nnelem
    real :: weta, wxi

    nnelem = src%nr
    nmat = Tdomain%specel(nnelem)%mat_index
    mat => Tdomain%sSubdomain(nmat)
    ngllx = mat%ngllx
    ngllz = mat%ngllz
    allocate(src%ExtForce(0:ngllx-1,0:ngllz-1,0:1))
    do j = 0,ngllz-1
        call pol_lagrange (ngllz, mat%GLLcz, j, src%eta,weta)
        do i = 0,ngllx-1
            call pol_lagrange (ngllx, mat%GLLcx, i, src%xi, wxi )
            src%ExtForce (i,j,0) = wxi*weta*src%dir(1)
            src%ExtForce (i,j,1) = wxi*weta*src%dir(2)
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
    integer :: i, j, ngllx, ngllz, nmat, nnelem
    real, dimension(0:1, 0:1) :: InvGrad, M
    real :: xi, eta, wxi, weta, dwdxi, dwdeta

    M = src%moment
    nnelem = src%nr
    elem => Tdomain%specel(nnelem)
    nmat = elem%mat_index
    mat => Tdomain%sSubdomain(nmat)
    ngllx = mat%ngllx
    ngllz = mat%ngllz

    InvGrad = src%Scoeff
    eta = src%eta
    xi = src%xi

    allocate(src%ExtForce(0:ngllx-1,0:ngllz-1,0:1))
    do j = 0,ngllz-1
        call pol_lagrange (ngllz, mat%GLLcz, j, eta,weta)
        call DERIVLAG (mat%GLLcz, ngllz, j, eta, dwdeta)
        do i = 0,ngllx-1
            call pol_lagrange (ngllx, mat%GLLcx, i, xi, wxi )
            call DERIVLAG ( mat%GLLcx, ngllx, i, xi, dwdxi)
            src%ExtForce (i,j,0) = &
                (InvGrad(0,0)*dwdxi*weta + InvGrad(0,1)*dwdeta*wxi)*M(0,0) + &
                (InvGrad(1,0)*dwdxi*weta + InvGrad(1,1)*dwdeta*wxi)*M(0,1)
            src%ExtForce (i,j,1) = &
                (InvGrad(0,0)*dwdxi*weta + InvGrad(0,1)*dwdeta*wxi)*M(1,0) + &
                (InvGrad(1,0)*dwdxi*weta + InvGrad(1,1)*dwdeta*wxi)*M(1,1)
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
    integer :: ipoint, nnelem
    real :: xi, eta, jac
    real :: x0, x1, x2, x3
    real :: z0, z1, z2, z3

    nnelem = src%nr
    elem => Tdomain%specel(nnelem)

    ipoint = elem%Control_Nodes(0); x0= Tdomain%Coord_nodes(0,ipoint); z0= Tdomain%Coord_nodes(1,ipoint)
    ipoint = elem%Control_Nodes(1); x1= Tdomain%Coord_nodes(0,ipoint); z1= Tdomain%Coord_nodes(1,ipoint)
    ipoint = elem%Control_Nodes(2); x2= Tdomain%Coord_nodes(0,ipoint); z2= Tdomain%Coord_nodes(1,ipoint)
    ipoint = elem%Control_Nodes(3); x3= Tdomain%Coord_nodes(0,ipoint); z3= Tdomain%Coord_nodes(1,ipoint)
    eta = src%eta
    xi = src%xi

    InvGrad(0,0) = 0.25 * ( (x1-x0) * (1-eta) + (x2-x3) * (1+eta) )
    InvGrad(1,0) = 0.25 * ( (x3-x0) * (1-xi) + (x2-x1) * (1+xi) )
    InvGrad(0,1) = 0.25 * ( (z1-z0) * (1-eta) + (z2-z3) * (1+eta) )
    InvGrad(1,1) = 0.25 * ( (z3-z0) * (1-xi) + (z2-z1) * (1+xi) )

    call invert2 (InvGrad, Jac )
    src%Scoeff = InvGrad
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
    integer :: ipoint, nnelem
    real :: xi, eta, jac
    real :: x0, x1, x2, x3, x4, x5, x6, x7
    real :: z0, z1, z2, z3, z4, z5, z6, z7

    nnelem = src%nr
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

    eta = src%eta
    xi = src%xi

    InvGrad(0,0) = 0.25 * (x0 *(1-eta)*(2*xi+eta) + x1 *(1-eta)*(2*xi-eta) + x2 *(1+eta)*(2*xi+eta)+ &
        x3 *(1+eta)*(2*xi-eta)) - x4 * xi*(1-eta) - x6 *xi *(1+eta) + 0.5* (x5-x7)* (1-eta**2)
    InvGrad(1,0) = 0.25 * (x0 *(1-xi)*(2*eta+xi) - x1 *(1+xi)*(xi-2*eta) + x2 *(1+xi)*(2*eta+xi)-  &
        x3 *(1-xi)*(xi-2*eta)) - x5 *eta*(1+xi) - x7 *eta *(1-xi) + 0.5* (x6-x4)* (1-xi**2)
    InvGrad(0,1) = 0.25 * (z0 *(1-eta)*(2*xi+eta) + z1 *(1-eta)*(2*xi-eta) + z2 *(1+eta)*(2*xi+eta)+ &
        z3 *(1+eta)*(2*xi-eta)) - z4 * xi*(1-eta) - z6 *xi *(1+eta) + 0.5* (z5-z7)* (1-eta**2)
    InvGrad(1,1) = 0.25 * (z0 *(1-xi)*(2*eta+xi) - z1 *(1+xi)*(xi-2*eta) + z2 *(1+xi)*(2*eta+xi)-  &
        z3 *(1-xi)*(xi-2*eta)) - z5 *eta*(1+xi) - z7 *eta *(1-xi) + 0.5* (z6-z4)* (1-xi**2)
    call invert2 (InvGrad, Jac)

    src%Scoeff = InvGrad
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
