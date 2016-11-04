!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file shape8.F90
!!\brief Contient la subroutine shape8.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Evalue les fonctions de forme.
!!
!! \param type(domain),target, intent (INOUT) Tdomain
!<

module shape_quad
    use constants
    implicit none
contains

    subroutine shape8(Tdomain)

        use sdomain
        use semdatafiles

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain


        ! local variables

        integer :: i_aus,n, mat,ngllx,ngllz,i,j,ipoint
        real, dimension(0:1,0:7) :: coord
        real :: xi,eta,xp,zp, Jac
        real, dimension (0:1,0:1) :: LocInvGrad


        ! Modified by Gaetano Festa, 27/05/05
        !---------------------------------------------------------------------------
        ! Shape functions are derived from "Finite Elements and Approximations"
        ! by Zienkiewicz, O.C. and Morgan, K.
        ! John Wiley and Sons, 1983
        ! --------------------------------------------------------------------------

        allocate (Tdomain%GlobCoord(0:1,0:Tdomain%n_glob_points-1))

        do n = 0,Tdomain%n_elem - 1
            do i=0,7
                i_aus = Tdomain%specel(n)%Control_Nodes(i)
                coord(0,i) = Tdomain%Coord_Nodes(0,i_aus)
                coord(1,i) = Tdomain%Coord_Nodes(1,i_aus)
            end do

            mat = Tdomain%specel(n)%mat_index
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz

            allocate (Tdomain%specel(n)%Jacob(0:ngllx-1,0:ngllz-1) )
            allocate (Tdomain%specel(n)%InvGrad(0:ngllx-1,0:ngllz-1,0:1,0:1) )

            do j = 0,ngllz - 1
                eta =   Tdomain%sSubdomain(mat)%GLLcz (j)
                do i = 0,ngllx - 1
                    xi = Tdomain%sSubdomain(mat)%GLLcx (i)

                    ! Computation of global coordinates

                    ipoint = Tdomain%specel(n)%Iglobnum(i,j)
                    call shape8_local2global(coord, xi, eta, xp, zp)
                    Tdomain%GlobCoord(0,ipoint) = xp; Tdomain%GlobCoord(1,ipoint) = zp

                    ! Computation of the derivative matrix

                    call shape8_local2jacob(coord, xi, eta, LocInvGrad)
                    call invert2(LocInvGrad, Jac)
                    Tdomain%specel(n)%InvGrad (i,j,0:1,0:1) = LocInvGrad(0:1,0:1)
                    Tdomain%specel(n)%Jacob (i,j) = Jac
                enddo
            enddo
            call buildNormal(Tdomain,n)
        enddo

        if (Tdomain%logicD%super_object_local_present) then
            do n = 0, Tdomain%n_fault-1
                call shape8_manage_super_object(Tdomain, n)
            end do
        endif
        return
    end subroutine shape8

    subroutine shape8_manage_super_object(Tdomain, n)
        use sdomain
        implicit none
        type(domain),target, intent (INOUT) :: Tdomain
        integer, intent(IN) :: n
        integer :: i, j, i_aus
        real, dimension (:), pointer :: GLLc_face
        real, dimension (:,:), allocatable :: Store_normal
        real :: ds_local, normal_0, normal_1, normalization
        real :: x0, x1, x2, z0, z1, z2, xp, zp
        integer :: Face_up, mat, ngll, nv, nv2
        real :: zin,xin,loc_distance


        do  j = 0, Tdomain%sFault(n)%n_face-1
            ngll = Tdomain%sFault(n)%fFace(j)%ngll
            Face_up = Tdomain%sFault(n)%fFace(j)%Face_UP
            nv2 = Tdomain%sFace(Face_up)%Near_element(0)
            mat = Tdomain%specel(nv2)%mat_index
            if (Tdomain%specel(nv2)%near_face(0) == Face_up .or. Tdomain%specel(nv2)%near_face(2) == face_up) then
                GLLc_face => Tdomain%sSubdomain(mat)%GLLcx
            else
                GLLc_face => Tdomain%sSubdomain(mat)%GLLcz
            endif

            if (Tdomain%specel(nv2)%near_face(0) == Face_up) then
                i_aus = Tdomain%specel(nv2)%Control_Nodes(0);
                x0 = Tdomain%Coord_Nodes(0,i_aus);  z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(1);
                x1 = Tdomain%Coord_Nodes(0,i_aus);  z1 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(4);
                x2 = Tdomain%Coord_Nodes(0,i_aus);  z2 = Tdomain%Coord_Nodes(1,i_aus)
            else if (Tdomain%specel(nv2)%near_face(1) == Face_up) then
                i_aus = Tdomain%specel(nv2)%Control_Nodes(1);
                x0 = Tdomain%Coord_Nodes(0,i_aus);  z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(2);
                x1 = Tdomain%Coord_Nodes(0,i_aus);  z1 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(5);
                x2 = Tdomain%Coord_Nodes(0,i_aus);  z2 = Tdomain%Coord_Nodes(1,i_aus)
            else if (Tdomain%specel(nv2)%near_face(2) == Face_up) then
                i_aus = Tdomain%specel(nv2)%Control_Nodes(2);
                x0 = Tdomain%Coord_Nodes(0,i_aus);  z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(3);
                x1 = Tdomain%Coord_Nodes(0,i_aus);  z1 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(6);
                x2 = Tdomain%Coord_Nodes(0,i_aus);  z2 = Tdomain%Coord_Nodes(1,i_aus)
            else
                i_aus = Tdomain%specel(nv2)%Control_Nodes(3);
                x0 = Tdomain%Coord_Nodes(0,i_aus);  z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(0);
                x1 = Tdomain%Coord_Nodes(0,i_aus);  z1 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(7);
                x2 = Tdomain%Coord_Nodes(0,i_aus);  z2 = Tdomain%Coord_Nodes(1,i_aus)
            endif

            allocate (Tdomain%sFault(n)%fFace(j)%ds(0:ngll-1))
            allocate (Tdomain%sFault(n)%fFace(j)%normal(0:ngll-1,0:1))
            allocate (Tdomain%sFault(n)%fFace(j)%distance(0:ngll-1))

            loc_distance = 0.
            xin = x0; zin = z0

            do i = 0, ngll-1
                xp = x1 * (GLLc_face(i)+0.5) + x0 * (GLLc_face(i)-0.5) - 2*x2* GLLc_face(i)
                zp  = z1 * (GLLc_face(i)+0.5) + z0 * (GLLc_face(i)-0.5) - 2*z2* GLLc_face(i)

                ds_local = xp**2 + zp**2
                ds_local = sqrt(ds_local)
                normal_0 = zp/ds_local; normal_1 = -xp/ds_local
                Tdomain%sFault(n)%fFace(j)%ds(i) = ds_local
                Tdomain%sFault(n)%fFace(j)%normal(i,0) = normal_0; Tdomain%sFault(n)%fFace(j)%normal(i,1) = normal_1

                Tdomain%sFault(n)%fFace(j)%distance(i) = loc_distance
                if (i < ngll-1) then
                    xp = 0.5 * GLLc_face(i+1)*(GLLc_face(i+1)+1) *x1 +  0.5 * GLLc_face(i+1)*(-GLLc_face(i+1)+1) *x0 + &
                        x2 * (1 -  GLLc_face(i+1)**2)
                    zp = 0.5 * GLLc_face(i+1)*(GLLc_face(i+1)+1) *z1 +  0.5 * GLLc_face(i+1)*(-GLLc_face(i+1)+1) *z0 + &
                        z2 * (1 -  GLLc_face(i+1)**2)
                    loc_distance = loc_distance + sqrt((xp-xin)**2 + (zp-zin)**2)
                    xin = xp; zin = zp
                endif
            enddo
            Tdomain%sFault(n)%fFace(j)%X_Vertex(0) = x0; Tdomain%sFault(n)%fFace(j)%X_Vertex(1) = x1
            Tdomain%sFault(n)%fFace(j)%Z_Vertex(0) = z0; Tdomain%sFault(n)%fFace(j)%Z_Vertex(1) = z1
        enddo

        do j = 0, Tdomain%sFault(n)%n_vertex-1
            Tdomain%sFault(n)%fVertex(j)%normal(0:1) = 0
        enddo

        do j = 0, Tdomain%sFault(n)%n_face-1
            ngll = Tdomain%sFault(n)%fFace(j)%ngll
            nv = Tdomain%sFault(n)%fFace(j)%Face_To_Vertex(0)
            Tdomain%sFault(n)%fVertex(nv)%Normal(0:1) = Tdomain%sFault(n)%fVertex(nv)%Normal(0:1) + &
                Tdomain%sFault(n)%fFace(j)%normal(0,0:1)
            nv = Tdomain%sFault(n)%fFace(j)%Face_To_Vertex(1)
            Tdomain%sFault(n)%fVertex(nv)%Normal(0:1) = Tdomain%sFault(n)%fVertex(nv)%Normal(0:1) + &
                Tdomain%sFault(n)%fFace(j)%normal(ngll-1,0:1)
        enddo

        do j = 0, Tdomain%sFault(n)%n_vertex-1
            normalization = Tdomain%sFault(n)%fVertex(j)%normal(0)**2+Tdomain%sFault(n)%fVertex(j)%normal(1)**2
            normalization = sqrt(normalization)
            Tdomain%sFault(n)%fVertex(j)%normal(:) = Tdomain%sFault(n)%fVertex(j)%normal(:)/normalization
        enddo

        do j = 0, Tdomain%sFault(n)%n_face-1
            ngll = Tdomain%sFault(n)%fFace(j)%ngll
            allocate (Store_normal(1:ngll-2,0:1))
            Store_normal (1:ngll-2,0:1) = Tdomain%sFault(n)%fFace(j)%normal(1:ngll-2,0:1)
            deallocate(Tdomain%sFault(n)%fFace(j)%normal)
            allocate(Tdomain%sFault(n)%fFace(j)%normal(1:ngll-2,0:1))
            Tdomain%sFault(n)%fFace(j)%normal(1:ngll-2,0:1) = Store_normal (1:ngll-2,0:1)
            deallocate (Store_normal)
        enddo

        ! Valid only if points are ordered
        ds_local = 0
        do j =  0, Tdomain%sFault(n)%n_face-1
            ngll = Tdomain%sFault(n)%fFace(j)%ngll
            do i= 0, ngll-1
                Tdomain%sFault(n)%fFace(j)%distance(i) = ds_local +  Tdomain%sFault(n)%fFace(j)%distance(i)
            enddo
            ds_local = Tdomain%sFault(n)%fFace(j)%distance(ngll-1)
            nv = Tdomain%sFault(n)%fFace(j)%Face_to_vertex(0)
            Tdomain%sFault(n)%fvertex(nv)%distance = Tdomain%sFault(n)%fFace(j)%distance(0)
            nv = Tdomain%sFault(n)%fFace(j)%Face_to_vertex(1)
            Tdomain%sFault(n)%fvertex(nv)%distance = Tdomain%sFault(n)%fFace(j)%distance(ngll-1)
            allocate (Store_normal(0:ngll-1,0:0))
            Store_normal(:,0) =  Tdomain%sFault(n)%fFace(j)%distance
            allocate ( Tdomain%sFault(n)%fFace(j)%distance(1:ngll-2))
            Tdomain%sFault(n)%fFace(j)%distance(1:ngll-2)= Store_normal (1:ngll-2,0)
            deallocate (Store_normal)
        enddo

    end subroutine shape8_manage_super_object

    subroutine buildNormal(Tdomain,n_elem)

        use sdomain

        implicit none

        ! subroutine arguments
        type(domain),target, intent (INOUT) :: Tdomain
        integer, intent (IN) :: n_elem

        ! local variables
        real, dimension (0:1,0:1) :: LocInvGrad
        real, dimension(0:1,0:7) :: coord
        integer :: i, j, i_aus, ngllx, ngllz, npg, k, mat, nf
        real    :: eta, xi, tx, tz, nx, nz, n_norm, w

        ! allocate space to store normals of each gauss point of each face of the element
        if (n_elem .lt. 0 .or. n_elem .ge. Tdomain%n_elem) stop "shape8 - buildNormal : invalid element"
        ngllx = Tdomain%specel(n_elem)%ngllx
        ngllz = Tdomain%specel(n_elem)%ngllz
        npg = 2*(ngllx+ngllz) ! number of gauss points over all element faces
        if (.not. allocated(Tdomain%specel(n_elem)%Normal_Nodes))       allocate(Tdomain%specel(n_elem)%Normal_Nodes(0:npg-1,0:1))
        if (.not. allocated(Tdomain%specel(n_elem)%Coeff_Integr_Faces)) allocate(Tdomain%specel(n_elem)%Coeff_Integr_Faces(0:npg-1))

        ! initialization
        mat = Tdomain%specel(n_elem)%mat_index
        k = 0 ! k scans gauss points face by face (with respect to each face order)
        do i=0,7
            i_aus = Tdomain%specel(n_elem)%Control_Nodes(i)
            coord(0,i) = Tdomain%Coord_Nodes(0,i_aus)
            coord(1,i) = Tdomain%Coord_Nodes(1,i_aus)
        end do

        ! face 0
        j = 0
        do i = 0, ngllx - 1
            ! Tangent from derivatives along xi
            xi  = Tdomain%sSubdomain(mat)%GLLcx(i)
            eta = Tdomain%sSubdomain(mat)%GLLcz(j)
            call shape8_local2jacob(coord, xi, eta, LocInvGrad)
            tx = LocInvGrad(0, 0) ! dx/dxi
            tz = LocInvGrad(0, 1) ! dz/dxi
            ! Normal build from tangent
            nx = tz ; nz = -tx
            n_norm = sqrt(nx**2 + nz**2)
            if (abs(n_norm) .le. 1.e-12) stop "buildNormal KO Face 0"
            ! Store normal at gauss point
            Tdomain%specel(n_elem)%Normal_Nodes(k, 0) = nx / n_norm
            Tdomain%specel(n_elem)%Normal_Nodes(k, 1) = nz / n_norm
            ! Impact surfacic integral (normal estimation)
            w = Tdomain%sSubdomain(mat)%GLLwx(i)
            Tdomain%specel(n_elem)%Coeff_Integr_Faces(k) = w * n_norm
            nf = Tdomain%specel(n_elem)%Near_Face(0)
            if (i ==       0) Tdomain%sface(nf)%Coeff_Integr_Ends(0) = Tdomain%specel(n_elem)%Coeff_Integr_Faces(k)
            if (i == ngllx-1) Tdomain%sface(nf)%Coeff_Integr_Ends(1) = Tdomain%specel(n_elem)%Coeff_Integr_Faces(k)
            ! Go to next gauss point
            k = k+1
        end do

        ! face 1
        i = ngllx - 1
        do j = 0, ngllz - 1
            ! Tangent from derivatives along eta
            xi  = Tdomain%sSubdomain(mat)%GLLcx(i)
            eta = Tdomain%sSubdomain(mat)%GLLcz(j)
            call shape8_local2jacob(coord, xi, eta, LocInvGrad)
            tx = LocInvGrad(1, 0) ! dx/deta
            tz = LocInvGrad(1, 1) ! dz/deta
            ! Normal build from tangent
            nx = tz ; nz = -tx
            n_norm = sqrt(nx**2 + nz**2)
            if (abs(n_norm) .le. 1.e-12) stop "buildNormal KO Face 1"
            ! Store normal at gauss point
            Tdomain%specel(n_elem)%Normal_Nodes(k, 0) = nx / n_norm
            Tdomain%specel(n_elem)%Normal_Nodes(k, 1) = nz / n_norm
            ! Impact surfacic integral (normal estimation)
            w = Tdomain%sSubdomain(mat)%GLLwz(j)
            Tdomain%specel(n_elem)%Coeff_Integr_Faces(k) = w * n_norm
            nf = Tdomain%specel(n_elem)%Near_Face(1)
            if (j ==       0) Tdomain%sface(nf)%Coeff_Integr_Ends(0) = Tdomain%specel(n_elem)%Coeff_Integr_Faces(k)
            if (j == ngllz-1) Tdomain%sface(nf)%Coeff_Integr_Ends(1) = Tdomain%specel(n_elem)%Coeff_Integr_Faces(k)
            ! Go to next gauss point
            k = k+1
        end do

        ! face 2
        j = ngllz - 1
        do i = 0, ngllx - 1
            ! Tangent from derivatives along eta
            xi  = Tdomain%sSubdomain(mat)%GLLcx(i)
            eta = Tdomain%sSubdomain(mat)%GLLcz(j)
            call shape8_local2jacob(coord, xi, eta, LocInvGrad)
            tx = LocInvGrad(0, 0) ! dx/dxi
            tz = LocInvGrad(0, 1) ! dz/dxi
            ! Normal build from tangent
            nx = -tz ; nz = tx
            n_norm = sqrt(nx**2 + nz**2)
            if (abs(n_norm) .le. 1.e-12) stop "buildNormal KO Face 2"
            ! Store normal at gauss point
            Tdomain%specel(n_elem)%Normal_Nodes(k, 0) = nx / n_norm
            Tdomain%specel(n_elem)%Normal_Nodes(k, 1) = nz / n_norm
            ! Impact surfacic integral (normal estimation)
            w = Tdomain%sSubdomain(mat)%GLLwx(i)
            Tdomain%specel(n_elem)%Coeff_Integr_Faces(k) = w * n_norm
            nf = Tdomain%specel(n_elem)%Near_Face(2)
            if (i ==       0) Tdomain%sface(nf)%Coeff_Integr_Ends(0) = Tdomain%specel(n_elem)%Coeff_Integr_Faces(k)
            if (i == ngllx-1) Tdomain%sface(nf)%Coeff_Integr_Ends(1) = Tdomain%specel(n_elem)%Coeff_Integr_Faces(k)
            ! Go to next gauss point
            k = k+1
        end do

        ! face 3
        i = 0
        do j = 0, ngllz - 1
            ! Tangent from derivatives along eta
            xi  = Tdomain%sSubdomain(mat)%GLLcx(i)
            eta = Tdomain%sSubdomain(mat)%GLLcz(j)
            call shape8_local2jacob(coord, xi, eta, LocInvGrad)
            tx = LocInvGrad(1, 0) ! dx/deta
            tz = LocInvGrad(1, 1) ! dz/deta
            ! Normal build from tangent
            nx = -tz ; nz = tx
            n_norm = sqrt(nx**2 + nz**2)
            if (abs(n_norm) .le. 1.e-12) stop "buildNormal KO Face 3"
            ! Store normal at gauss point
            Tdomain%specel(n_elem)%Normal_Nodes(k, 0) = nx / n_norm
            Tdomain%specel(n_elem)%Normal_Nodes(k, 1) = nz / n_norm
            ! Impact surfacic integral (normal estimation)
            w = Tdomain%sSubdomain(mat)%GLLwz(j)
            Tdomain%specel(n_elem)%Coeff_Integr_Faces(k) = w * n_norm
            nf = Tdomain%specel(n_elem)%Near_Face(3)
            if (j ==       0) Tdomain%sface(nf)%Coeff_Integr_Ends(0) = Tdomain%specel(n_elem)%Coeff_Integr_Faces(k)
            if (j == ngllz-1) Tdomain%sface(nf)%Coeff_Integr_Ends(1) = Tdomain%specel(n_elem)%Coeff_Integr_Faces(k)
            ! Go to next gauss point
            k = k+1
        end do
    end subroutine buildNormal


    subroutine shape8_global2local(coord, xa, za, xi, eta, ok)
        double precision, dimension(0:1,0:7), intent(in)  :: coord
        double precision, intent(in) :: xa, za
        double precision, intent(out) :: xi, eta
        logical, intent(out) :: ok
        !
        integer :: niter
        double precision, dimension(0:1) :: xin, xout, xref
        ok = .true.
        xin(0) = 0.
        xin(1) = 0.
        xref(0) = xa
        xref(1) = za
        call shape8_simple_newton(coord, xref, xin, xout, niter)
        if (niter==1000 .or. niter<0) ok=.false.
        xi  = xout(0)
        eta = xout(1)
    end subroutine shape8_global2local

    subroutine shape8_simple_newton(nodes, xref, xin, xout, nit)
        double precision, dimension(0:1,0:7), intent(in) :: nodes
        double precision, dimension(0:1), intent(in) :: xref, xin
        double precision, dimension(0:1), intent(out) :: xout
        integer, intent(out) :: nit
        !
        double precision, dimension(0:1,0:1) :: jac
        double precision, dimension(0:1) :: x
        double precision :: xa, za, err, Det
        integer, parameter :: niter=1000
        integer :: i
        xout = xin
        do i=1,niter
            call shape8_local2global(nodes, xout(0), xout(1), xa, za)
            call shape8_local2jacob (nodes, xout(0), xout(1), Jac)
            call invert2 (Jac, Det)
            xa = xref(0)-xa
            za = xref(1)-za
            x(0) = Jac(0,0)*xa + Jac(1,0)*za
            x(1) = Jac(0,1)*xa + Jac(1,1)*za
            err = x(0)**2+x(1)**2
            if (err<1e-24) exit ! eps^2 with eps=1e-12 : for consistency with shape4
            xout = xout + x
        end do
        nit = i
    end subroutine shape8_simple_newton

    subroutine shape8_local2global(coord, xi, eta, xa, za)
        double precision, dimension(0:1,0:7), intent(in) :: coord
        double precision, intent(in) :: xi, eta
        double precision, intent(out) :: xa, za
        !
        real, dimension(0:7) :: xc, zc ! Used for clarity (-O2 should be able to optimize this)
        xc = coord(0,:)
        zc = coord(1,:)
        !- computation of the global coordinates from the local coordinates
        xa = 0.25*(-xc(0)*(1.-xi)*(1.-eta)*(1+xi+eta)  &
                   -xc(1)*(1.+xi)*(1.-eta)*(1-xi+eta)  &
                   -xc(2)*(1.+xi)*(1.+eta)*(1-xi-eta)  &
                   -xc(3)*(1.-xi)*(1.+eta)*(1+xi-eta)) &
            + 0.5*( xc(4)*(1.-xi**2)*(1.-eta) +        &
                    xc(5)*(1.+xi)*(1.-eta**2) +        &
                    xc(6)*(1.-xi**2)*(1.+eta) +        &
                    xc(7)*(1.-xi)*(1.-eta**2) )
        za = 0.25*(-zc(0)*(1.-xi)*(1.-eta)*(1+xi+eta)  &
                   -zc(1)*(1.+xi)*(1.-eta)*(1-xi+eta)  &
                   -zc(2)*(1.+xi)*(1.+eta)*(1-xi-eta)  &
                   -zc(3)*(1.-xi)*(1.+eta)*(1+xi-eta)) &
            + 0.5*( zc(4)*(1.-xi**2)*(1.-eta) +        &
                    zc(5)*(1.+xi)*(1.-eta**2) +        &
                    zc(6)*(1.-xi**2)*(1.+eta) +        &
                    zc(7)*(1.-xi)*(1.-eta**2) )
    end subroutine shape8_local2global

    subroutine shape8_local2jacob(coord, xi, eta, jac)
        double precision, dimension(0:1,0:7), intent(in) :: coord
        double precision, intent(in) :: xi, eta
        double precision, dimension(0:1,0:1), intent(out) :: jac
        !
        real, dimension(0:7) :: xc, zc ! Used for clarity (-O2 should be able to optimize this)
        xc = coord(0,:)
        zc = coord(1,:)
        !- computation of the derivative matrix, dx_(jj)/dxi_(ii)
        jac(0,0) = 0.25*( xc(0)*(1-eta)*(2*xi+eta)    &
                         +xc(1)*(1-eta)*(2*xi-eta)    &
                         +xc(2)*(1+eta)*(2*xi+eta)    &
                         +xc(3)*(1+eta)*(2*xi-eta))   &
                   -xc(4)*xi*(1-eta)-xc(6)*xi*(1+eta) &
                   +0.5*(xc(5)-xc(7))*(1-eta**2)
        jac(1,0) = 0.25*( xc(0)*(1-xi)*(2*eta+xi)     &
                         -xc(1)*(1+xi)*(xi-2*eta)     &
                         +xc(2)*(1+xi)*(2*eta+xi)     &
                         -xc(3)*(1-xi)*(xi-2*eta))    &
                   -xc(5)*eta*(1+xi)-xc(7)*eta*(1-xi) &
                   +0.5*(xc(6)-xc(4))*(1-xi**2)
        jac(0,1) = 0.25*( zc(0)*(1-eta)*(2*xi+eta)    &
                         +zc(1)*(1-eta)*(2*xi-eta)    &
                         +zc(2)*(1+eta)*(2*xi+eta)    &
                         +zc(3)*(1+eta)*(2*xi-eta))   &
                   -zc(4)*xi*(1-eta)-zc(6)*xi*(1+eta) &
                   +0.5*(zc(5)-zc(7))*(1-eta**2)
        jac(1,1) = 0.25*( zc(0)*(1-xi)*(2*eta+xi)     &
                         -zc(1)*(1+xi)*(xi-2*eta)     &
                         +zc(2)*(1+xi)*(2*eta+xi)     &
                         -zc(3)*(1-xi)*(xi-2*eta))    &
                   -zc(5)*eta*(1+xi)-zc(7)*eta*(1-xi) &
                   +0.5*(zc(6)-zc(4))*(1-xi**2)
    end subroutine shape8_local2jacob

end module shape_quad

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
