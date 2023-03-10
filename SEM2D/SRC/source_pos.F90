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
    use constants
    use mlocations2d
    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables

    integer, parameter :: NMAXEL=20
    integer, dimension(NMAXEL) :: elems
    real(fpp), dimension(0:1,NMAXEL) :: coordloc
    real(fpp), parameter :: EPS = 1D-13
    integer   :: nsour, nmax, i, n_el
    real(fpp) :: xc, zc, xi, eta
    logical :: inside

    ! Find the nearest point to the real location of the sources in GLL scheme

    do nsour = 0, Tdomain%n_source -1
        ! Find source location
        xc = Tdomain%sSource(nsour)%xsource
        zc = Tdomain%sSource(nsour)%zsource
        nmax = NMAXEL
        call find_location(Tdomain, xc, zc, nmax, elems, coordloc)

        Tdomain%sSource(nsour)%located_here = .false.
        Tdomain%sSource(nsour)%ine = 0
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
                Tdomain%sSource(nsour)%ine          = 1 ! Source spread over one element by default (spread over more elements later on if necessary)
                allocate(Tdomain%sSource(nsour)%Elem(0:Tdomain%sSource(nsour)%ine-1))
                Tdomain%sSource(nsour)%Elem(0)%nr   = n_el
                Tdomain%sSource(nsour)%Elem(0)%xi   = xi
                Tdomain%sSource(nsour)%Elem(0)%eta  = eta

                exit ! Choose the first one to be consistent with previous behavior
            end if
        end do

        ! Customizing according to source type

        if (Tdomain%sSource(nsour)%located_here) then
            if (Tdomain%sSource(nsour)%i_type_source == 1) then       ! Pulse directional force
                call source_excit_pulse(Tdomain, Tdomain%sSource(nsour))
            else if (Tdomain%sSource(nsour)%i_type_source  == 2) then ! Explosive source diagonal moment considered
                if (Tdomain%n_nodes == 4) call calc_shape4_coeffs(Tdomain, Tdomain%sSource(nsour))
                if (Tdomain%n_nodes == 8) call calc_shape8_coeffs(Tdomain, Tdomain%sSource(nsour))
                call source_excit_moment(Tdomain, Tdomain%sSource(nsour))
            else if (Tdomain%sSource(nsour)%i_type_source  == 4) then ! Dirac projected
                call source_dirac_projected(Tdomain, Tdomain%sSource(nsour))
            else if (Tdomain%sSource(nsour)%i_type_source  == 6) then ! Source on strain equation
                call source_excit_strain(Tdomain, Tdomain%sSource(nsour))
            else if  (Tdomain%sSource(nsour)%i_type_source == 5) then
                call source_space_gaussian(Tdomain, Tdomain%sSource(nsour))
            endif
        end if
    enddo
end subroutine SourcePosition

subroutine source_excit_pulse(Tdomain, src)
    use sdomain
    use constants
    use ssources
    use ssubdomains
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(Subdomain), pointer :: mat
    type(Source), intent(inout) :: src
    integer :: n, i, j, ngllx, ngllz, nmat, nnelem
    real(fpp) :: weta, wxi

    do n = 0, src%ine-1
        nnelem = src%Elem(n)%nr
        nmat = Tdomain%specel(nnelem)%mat_index
        mat => Tdomain%sSubdomain(nmat)
        ngllx = mat%ngllx
        ngllz = mat%ngllz

        allocate(src%Elem(n)%ExtForce(0:ngllx-1,0:ngllz-1,0:1))
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


! Following subroutine added for testing JP proposition
subroutine source_excit_strain(Tdomain, src)
    use sdomain
    use constants
    use ssources
    use ssubdomains
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(Subdomain), pointer :: mat
    type(Source), intent(inout) :: src
    real(fpp), dimension(0:1, 0:1)   :: M
    integer :: n, i, j, ngllx, ngllz, nmat, nnelem
    real(fpp) :: weta, wxi, invE, nu

    do n = 0, src%ine-1
        nnelem = src%Elem(n)%nr
        nmat = Tdomain%specel(nnelem)%mat_index
        mat => Tdomain%sSubdomain(nmat)
        ngllx = mat%ngllx
        ngllz = mat%ngllz
        M = src%moment
        ! Computing parameters for the compliance matrix :
        !invE = (mat%Dlambda + mat%Dmu) / (mat%Dmu*(3.*mat%Dlambda + 2.*mat%Dmu))
        !nu   =  mat%Dlambda / (2.*mat%Dmu + 2.*mat%Dlambda)
        invE = 1. ! A SUPPRIMER !!!!!!
        nu   = 0. ! A SUPPRIMER !!!!!!
        src%Elem(n)%invE = invE
        src%Elem(n)%nu   = nu
        allocate  (src%Elem(n)%ExtForce(0:ngllx-1,0:ngllz-1,0:2))
        do j = 0,ngllz-1
            call pol_lagrange (ngllz, mat%GLLcz, j, src%Elem(n)%eta,weta)
            do i = 0,ngllx-1
                call pol_lagrange (ngllx, mat%GLLcx, i, src%Elem(n)%xi, wxi )
                src%Elem(n)%ExtForce (i,j,0) = wxi*weta*invE*(M(0,0)-nu*M(1,1))
                src%Elem(n)%ExtForce (i,j,1) = wxi*weta*invE*(M(1,1)-nu*M(0,0))
                src%Elem(n)%ExtForce (i,j,2) = wxi*weta*invE*(1.+nu)*M(0,1)
            enddo
        enddo
    enddo

end subroutine source_excit_strain


subroutine source_excit_moment(Tdomain, src)
    use sdomain
    use constants
    use ssources
    use ssubdomains
    use selement
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(Subdomain), pointer :: mat
    type(element), pointer :: elem
    type(Source), intent(inout) :: src
    integer :: n, i, j, ngllx, ngllz, nmat, nnelem
    real(fpp), dimension(0:1, 0:1) :: InvGrad, M
    real(fpp) :: xi, eta, wxi, weta, dwdxi, dwdeta

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

        allocate(src%Elem(n)%ExtForce(0:ngllx-1,0:ngllz-1,0:1))
        do j = 0,ngllz-1
            call pol_lagrange (ngllz, mat%GLLcz, j, eta,weta)
            call der_lagrange (ngllz, mat%GLLcz, j, eta, dwdeta)
            do i = 0,ngllx-1
                call pol_lagrange (ngllx, mat%GLLcx, i, xi, wxi )
                call der_lagrange (ngllx, mat%GLLcx, i, xi, dwdxi)
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

! ############################################################
!>
!! \brief This subroutine computes the distribution of the forces on an element
!! approximating a Dirac function by a finite sum of Legendre Polynomials.
!!
!! \param type (Element), intent (INOUT) Tdomain
!! \param type (source),  intent (INOUT) src
!<
subroutine source_dirac_projected(Tdomain, src)
    use sdomain
    use constants
    use ssources
    use ssubdomains
    use selement
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(Source), intent(inout) :: src
    type(Subdomain), pointer :: mat
    type(element), pointer :: elem
    real(fpp), dimension(0:1, 0:1) :: M
    real(fpp), dimension(:,:), allocatable :: Dirac_projected, whei, xix, xiz, etax, etaz, aux1, aux2
    real(fpp), dimension(:),   allocatable :: Fsurf1, Fsurf2
    integer  :: n, i, j, ngllx, ngllz, nelem, nmat, imin, imax
    real(fpp):: xi, eta

    M = src%moment
    do n = 0, src%ine-1 ! loop on all elements
        nelem = src%Elem(n)%nr
        nmat = Tdomain%specel(nelem)%mat_index
        mat => Tdomain%sSubdomain(nmat)
        elem => Tdomain%specel(nelem)
        eta = src%Elem(n)%eta
        xi = src%Elem(n)%xi
        ngllx = elem%ngllx
        ngllz = elem%ngllz
        allocate (xix  (0:ngllx-1,0:ngllz-1))
        allocate (xiz  (0:ngllx-1,0:ngllz-1))
        allocate (etax (0:ngllx-1,0:ngllz-1))
        allocate (etaz (0:ngllx-1,0:ngllz-1))
        allocate (whei (0:ngllx-1,0:ngllz-1))
        allocate (aux1 (0:ngllx-1,0:ngllz-1))
        allocate (aux2 (0:ngllx-1,0:ngllz-1))
        allocate (Dirac_projected(0:ngllx-1,0:ngllz-1))
        xix  = elem%InvGrad(:,:,0,0)
        xiz  = elem%InvGrad(:,:,1,0)
        etax = elem%InvGrad(:,:,0,1)
        etaz = elem%InvGrad(:,:,1,1)
        do j = 0,ngllz-1
            do i = 0,ngllx-1
                Whei (i,j) = mat%GLLwx(i)* mat%GLLwz(j)
            enddo
        enddo
        call project_dirac_on_Legendre(Dirac_projected,mat%gllcx,mat%gllcz,xi,eta,ngllx,ngllz)
        allocate  (src%Elem(n)%ExtForce(0:ngllx-1,0:ngllz-1,0:1))
        aux1 = (M(0,0) * xix + M(0,1) * xiz) * whei(:,:) * Dirac_projected(:,:) !* elem%Jacob(:,:)
        aux2 = (M(0,0) *etax + M(0,1) *etaz) * whei(:,:) * Dirac_projected(:,:) !* elem%Jacob(:,:)
        src%Elem(n)%ExtForce(:,:,0) = MATMUL(mat%hprimex,aux1) + MATMUL(aux2,mat%hTprimez)
        aux1 = (M(1,0) * xix + M(1,1) * xiz) * whei(:,:) * Dirac_projected(:,:) !* elem%Jacob(:,:)
        aux2 = (M(1,0) *etax + M(1,1) *etaz) * whei(:,:) * Dirac_projected(:,:) !* elem%Jacob(:,:)
        src%Elem(n)%ExtForce(:,:,1) = MATMUL(mat%hprimex,aux1) + MATMUL(aux2,mat%hTprimez)

        ! contribution integrale surfacique :
        allocate(Fsurf1(0:2*(ngllx+ngllz)-1))
        allocate(Fsurf2(0:2*(ngllx+ngllz)-1))
        Fsurf1 = (M(0,0)*Elem%Normal_Nodes(:,0) + M(0,1)*Elem%Normal_Nodes(:,1)) &
                       * Elem%Coeff_integr_Faces(:) * 0.
        Fsurf2 = (M(1,0)*Elem%Normal_Nodes(:,0) + M(1,1)*Elem%Normal_Nodes(:,1)) &
                       * Elem%Coeff_integr_Faces(:) * 0.
        ! For the Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        src%Elem(n)%ExtForce(0:ngllx-1,0,0) = src%Elem(n)%ExtForce(0:ngllx-1,0,0) &
                  - Dirac_projected(0:ngllx-1,0) * Fsurf1(imin:imax)
        src%Elem(n)%ExtForce(0:ngllx-1,0,1) = src%Elem(n)%ExtForce(0:ngllx-1,0,1) &
                  - Dirac_projected(0:ngllx-1,0) * Fsurf2(imin:imax)
        ! For the right Face :
        call get_iminimax(Elem,1,imin,imax)
        src%Elem(n)%ExtForce(ngllx-1,0:ngllz-1,0) = src%Elem(n)%ExtForce(ngllx-1,0:ngllz-1,0) &
                  - Dirac_projected(ngllx-1,0:ngllz-1) * Fsurf1(imin:imax)
        src%Elem(n)%ExtForce(ngllx-1,0:ngllz-1,1) = src%Elem(n)%ExtForce(ngllx-1,0:ngllz-1,1) &
                  - Dirac_projected(ngllx-1,0:ngllz-1) * Fsurf2(imin:imax)
        ! For the Top Face :
        call get_iminimax(Elem,2,imin,imax)
        src%Elem(n)%ExtForce(0:ngllx-1,ngllz-1,0) = src%Elem(n)%ExtForce(0:ngllx-1,ngllz-1,0) &
                  - Dirac_projected(0:ngllx-1,ngllz-1) * Fsurf1(imin:imax)
        src%Elem(n)%ExtForce(0:ngllx-1,ngllz-1,1) = src%Elem(n)%ExtForce(0:ngllx-1,ngllz-1,1) &
                  - Dirac_projected(0:ngllx-1,ngllz-1) * Fsurf2(imin:imax)
        ! For the Left Face :
        call get_iminimax(Elem,3,imin,imax)
        src%Elem(n)%ExtForce(0,0:ngllz-1,0) = src%Elem(n)%ExtForce(0,0:ngllz-1,0) &
                  - Dirac_projected(ngllx-1,0:ngllz-1) * Fsurf1(imin:imax)
        src%Elem(n)%ExtForce(0,0:ngllz-1,1) = src%Elem(n)%ExtForce(0,0:ngllz-1,1) &
                  - Dirac_projected(0,0:ngllz-1) * Fsurf2(imin:imax)
        deallocate(Dirac_projected, whei, xix, xiz, etax, etaz, aux1, aux2, Fsurf1, Fsurf2)
    enddo

end subroutine source_dirac_projected


! ############################################################
!>
!! \brief This subroutine uses the decomposition of the dirac
!! function Delta(x-xs,y-ys) on the Legendre basis on the source
!! element and then the subroutine computes the values of this
!! decomposition on all the gll nodes of the element.
!! Then everything is stored in the matrix : Dirac_projected
!!
!! \param type (Element), intent (INOUT) Tdomain
!! \param type (source),  intent (INOUT) src
!<
subroutine project_dirac_on_Legendre(Dirac_projected,GLLcx,GLLcz,xi,eta,ngllx,ngllz)

    use constants
    use splib, only : valepo

    real(fpp), dimension (0:ngllx-1, 0:ngllz-1), intent(INOUT) :: Dirac_projected
    real(fpp), dimension (0:ngllx-1), intent(IN) :: GLLcx
    real(fpp), dimension (0:ngllz-1), intent(IN) :: GLLcz
    real(fpp),    intent(IN) :: xi, eta
    integer, intent(IN) :: ngllx, ngllz
    integer  :: i, j, n, m
    real(fpp):: interp_at_gll, d1, d2, resx, resy, xi_i, eta_j
    real(fpp), dimension (0:ngllx-1) :: Pn_xs
    real(fpp), dimension (0:ngllz-1) :: Pm_ys

    ! Compute values of the Legendre polynomials Pn on the sources coords (eta,xi)
    do n=0,ngllx-1
        call VALEPO(n,xi, resx,d1,d2)
        Pn_xs(n) = resx
    enddo
    do n=0,ngllz-1
        call VALEPO(n,eta,resy,d1,d2)
        Pm_ys(n) = resy
    enddo
    do i=0,ngllx-1
        do j=0,ngllz-1
            ! For each gll (eta_i,xi_j), polynomial projection is computed
            interp_at_gll = 0.
            xi_i = GLLcx(i)
            eta_j = GLLcz(j)
            do m=0,ngllx-1
                do n=0,ngllz-1
                    call VALEPO(m,xi_i, resx,d1,d2)
                    call VALEPO(n,eta_j,resy,d1,d2)
                    interp_at_gll = interp_at_gll + (m+0.5)*(n+0.5)*Pn_xs(m)*Pm_ys(n)*resx*resy
                enddo
            enddo
            Dirac_projected(i,j) = interp_at_gll
        enddo
    enddo
end subroutine project_dirac_on_Legendre


! ############################################################
!>
!! \brief This subroutine deals with a source which in spatially
!! mollified by a Gaussian function.
!! In a first step, the subroutine finds the elements that are in the scope
!! of the Gaussian, and therefore these elements are condidered as sources.
!! In a second step, the matrices ExtForce(:,:,:) are computed for these
!! elements using a analytic computation of the derivatives of the Gaussian.
!!
!! \param type (Element), intent (INOUT) Tdomain
!! \param type (source),  intent (INOUT) src
!<
subroutine source_space_gaussian(Tdomain, src)
    use sdomain
    use constants
    use ssources
    use ssubdomains
    use selement
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(Source), intent(inout) :: src
    type(Subdomain), pointer    :: mat
    real(fpp), dimension(0:1, 0:1)   :: M
    real(fpp) :: Tol, Dx, Dz, Dist, res, dGauss_dx, dGauss_dz, WheiJac
    integer   :: n, ng, ns, nel, nmat, mind, ngllx, ngllz, i, j

    M = src%moment
    nel  = Tdomain%n_elem
    mind = 0

    ! Fist step : finding elements concerned by the gaussian Source
    do n = 0,nel-1
       Tol = 1.E-7
       ngllx = Tdomain%specel(n)%ngllx
       ngllz = Tdomain%specel(n)%ngllz
       Tdomain%specel(n)%is_source = .false.
       do j = 0,ngllz-1
          do i = 0,ngllx-1
             ng = Tdomain%specel(n)%Iglobnum(i,j)
             Dist = (Tdomain%GlobCoord(0,ng)-src%Xsource)**2  +  &
                    (Tdomain%GlobCoord(1,ng)-src%Zsource)**2
             call Gaussian2d(Dist,src%sigma,res)
             if ((res > Tol) .and. (.not. Tdomain%specel(n)%is_source)) then
                Tdomain%specel(n)%is_source = .true.
                mind = mind +1
             endif
          enddo
       enddo
    enddo

    write (*,*) "Gaussian Source, elements concerned by the space source Gaussian : ", mind
    src%ine = mind
    deallocate(src%Elem)
    allocate (src%Elem(0:mind-1))
    ns = 0
    ! Second step : computing coefficient source matrix ExtForces for source elements
    do n = 0,nel-1
       if (Tdomain%specel(n)%is_source) then
          src%Elem(ns)%nr = n
          ngllx = Tdomain%specel(n)%ngllx
          ngllz = Tdomain%specel(n)%ngllz
          nmat  = Tdomain%specel(n)%mat_index
          mat  => Tdomain%sSubdomain(nmat)
          allocate(src%Elem(ns)%ExtForce(0:ngllx-1,0:ngllz-1,0:1))
          do j = 0,ngllz-1
             do i = 0,ngllx-1
                ng = Tdomain%specel(n)%Iglobnum(i,j)
                WheiJac = Tdomain%specel(n)%Jacob(i,j) * mat%GLLwx(i) * mat%GLLwz(j)
                Dx = Tdomain%GlobCoord(0,ng)-src%Xsource
                Dz = Tdomain%GlobCoord(1,ng)-src%Zsource
                Dist = Dx**2 + Dz**2
                call Gaussian2d(Dist,src%sigma,res)
                dGauss_dx = -Dx / (src%sigma**2) * res
                dGauss_dz = -Dz / (src%sigma**2) * res
                src%Elem(ns)%ExtForce(i,j,0) = -wheiJac * (dGauss_dx*M(0,0) + dGauss_dz*M(0,1))
                src%Elem(ns)%ExtForce(i,j,1) = -wheiJac * (dGauss_dx*M(1,0) + dGauss_dz*M(1,1))
             enddo
          enddo
          ns = ns+1
       endif
    enddo
    return
end subroutine source_space_gaussian

! ############################################################

subroutine calc_shape4_coeffs(Tdomain, src)
    use sdomain
    use constants
    use ssources
    use ssubdomains
    use selement
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(element), pointer :: elem
    type(Source), intent(inout) :: src
    real(fpp), dimension(0:1, 0:1) :: InvGrad
    integer :: n, ipoint, nnelem
    real(fpp) :: xi, eta, jac
    real(fpp) :: x0, x1, x2, x3
    real(fpp) :: z0, z1, z2, z3

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
    use constants
    use ssources
    use ssubdomains
    use selement
    implicit none
    type(Domain), intent(inout) :: Tdomain
    type(element), pointer :: elem
    type(Source), intent(inout) :: src
    real(fpp), dimension(0:1, 0:1) :: InvGrad
    integer :: n, ipoint, nnelem
    real(fpp) :: xi, eta, jac
    real(fpp) :: x0, x1, x2, x3, x4, x5, x6, x7
    real(fpp) :: z0, z1, z2, z3, z4, z5, z6, z7

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



subroutine  Gaussian2d(Dist,sigma,res)

  use constants
  implicit none

  real(fpp), intent(in)  :: Dist, sigma
  real(fpp), intent(out) :: res
  real(fpp)              :: PI

  PI = Acos(-1.)
  res = 1./(sigma * sqrt(2.*PI)) * exp(-Dist/(2.*sigma**2))

end subroutine Gaussian2d


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
