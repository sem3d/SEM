!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file shape4.F90
!!\brief Contient la subroutine shape4.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module shape_lin

    use constants

contains

    ! #########################################
    subroutine compute_normals(Tdomain,nf)

        use sdomain

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain
        integer, intent (IN) :: nf

        ! local variables
        real(fpp), dimension(0:1) :: normal_aux
        integer :: i, n_elem, ngll, nv0, nv1

        n_elem = Tdomain%sFace(nf)%Near_Element(0)
        ngll   = Tdomain%sFace(nf)%ngll
        nv0    = Tdomain%sFace(nf)%Near_Vertex(0)
        nv1    = Tdomain%sFace(nf)%Near_Vertex(1)
        allocate(Tdomain%sFace(nf)%Normal(0:1))
        allocate(Tdomain%sFace(nf)%Normal_Nodes(0:ngll-1,0:1))

        ! Computation of an unique Face normal
        call n_from_vertices(Tdomain, normal_aux, nv0, nv1)
        Tdomain%sFace(nf)%Normal = normal_aux
        ! Check if normal is outward...
        if (Tdomain%sFace(nf)%Which_Face(0) .GE. 2) then
            Tdomain%sFace(nf)%Normal(:) = - Tdomain%sFace(nf)%Normal(:)
        endif

        ! Computation of a field of normals for each Gauss points
        ! For shape4, elements are linears, and the normal for the
        ! Gauss points are the same than the Face Normal
        do i = 0,ngll-1
            Tdomain%sFace(nf)%Normal_Nodes(i,0) = Tdomain%sFace(nf)%Normal(0)
            Tdomain%sFace(nf)%Normal_Nodes(i,1) = Tdomain%sFace(nf)%Normal(1)
        end do

    end subroutine compute_normals


! #########################################
    subroutine compute_normals_elem(Tdomain, nelem)

        use sdomain

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain
        integer, intent (IN) :: nelem

        ! local variables
        integer              :: i, j, k, ngllx, ngllz, ngll, nglltot, nv0, nv1
        real(fpp), dimension(0:1) :: normal_aux

        ngllx  = Tdomain%specel(nelem)%ngllx
        ngllz  = Tdomain%specel(nelem)%ngllz
        nglltot= 2*(ngllx+ngllz)
        k = 0

        allocate(Tdomain%specel(nelem)%Normal_nodes(0:nglltot-1,0:1))

        do i=0,3
            nv0  = Tdomain%specel(nelem)%Near_Vertex(i)
            nv1  = Tdomain%specel(nelem)%Near_Vertex(mod(i+1,4))
            call n_from_vertices(Tdomain,normal_aux,nv0,nv1)
            if (i==0 .OR. i==2) then
                ngll = ngllx
            else
                ngll = ngllz
            endif
            do j=0,ngll-1
                Tdomain%specel(nelem)%Normal_Nodes(k,0) = normal_aux(0)
                Tdomain%specel(nelem)%Normal_Nodes(k,1) = normal_aux(1)
                k = k+1
            enddo
        enddo

    end subroutine compute_normals_elem


    ! #########################################
    !>
    !!\file shape4.F90
    !!\brief
    !!\version 1.0
    !!\date 01/04/2014
    !! This function computes a normal from 2 vertices
    !! Be carefull : this normal may not be outward !
    !<
    subroutine n_from_vertices(Tdomain,n_vec,nv0,nv1)

        use sDomain
        implicit none

        type (domain),target, intent (IN)    :: Tdomain
        real(fpp), dimension(0:1), intent (INOUT) :: n_vec
        integer, intent(IN)        :: nv0, nv1
        integer                    :: n0, n1
        real(fpp)                  :: dx, dy, nx, ny, norm

        n0 = Tdomain%sVertex(nv0)%Glob_numbering
        n1 = Tdomain%sVertex(nv1)%Glob_numbering
        dx = Tdomain%coord_nodes(0,n1) - Tdomain%coord_nodes(0,n0)
        dy = Tdomain%coord_nodes(1,n1) - Tdomain%coord_nodes(1,n0)
        nx  = dy ; ny = -dx
        norm = sqrt(nx*nx + ny*ny)
        n_vec(0) = nx / norm
        n_vec(1) = ny / norm

    end subroutine n_from_vertices


    ! #########################################
    !>
    !!\file shape4.F90
    !!\brief
    !!\version 1.0
    !!\date 17/01/2014
    !! This subroutine computes the Jacobian of the transformation
    !! from a reference segment [-1,1] to the actual face nf. This
    !! tranformation is, in the linear case, just homothetic, therefore
    !! the Jacobian is the same for all the nodes of the face, and its
    !! value is equal to elongation coefficient.
    !<
    subroutine compute_Jacobian_1D(Tdomain,nf,Jac1D)

        use sdomain

        implicit none

        type(domain), intent (IN) :: Tdomain
        integer, intent (IN) :: nf
        real(fpp), intent (INOUT) :: Jac1D

        ! local variables
        integer :: nv0, nv1, n0, n1

        nv0= Tdomain%sFace(nf)%Near_Vertex(0)
        nv1= Tdomain%sFace(nf)%Near_Vertex(1)
        n0 = Tdomain%sVertex(nv0)%Glob_numbering
        n1 = Tdomain%sVertex(nv1)%Glob_numbering

        Jac1D = 0.5 * sqrt((Tdomain%Coord_Nodes(0,n0)-Tdomain%Coord_Nodes(0,n1))**2 &
            +(Tdomain%Coord_Nodes(1,n0)-Tdomain%Coord_Nodes(1,n1))**2 )

    end subroutine compute_Jacobian_1D

    ! #########################################

    subroutine shape4(Tdomain)

        use sdomain
        use semdatafiles

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain

        ! local variables

        integer :: i_aus, n, mat, ngllx, ngllz, i, j, ipoint, imin, imax
        real(fpp), dimension(0:1,0:3) :: coord
        real(fpp) :: xi, eta, xp, zp, Jac
        real(fpp), dimension (0:1,0:1) :: LocInvGrad
        integer   :: nf
        real(fpp) :: Jac1D

        ! Modified by Gaetano Festa, 26/05/05
        !---------------------------------------------------------------------------------------------------------------
        ! Shape functions are derived from "Finite Elements and Approximations"
        ! by Zienkiewicz, O.C. and Morgan, K.
        ! John Wiley and Sons, 1983
        ! --------------------------------------------------------------------------------------------------------------

        allocate (Tdomain%GlobCoord(0:1,0:Tdomain%n_glob_points-1))

        do n = 0,Tdomain%n_elem - 1
            i_aus = Tdomain%specel(n)%Control_Nodes(0); coord(0,0) = Tdomain%Coord_Nodes(0,i_aus); coord(1,0) = Tdomain%Coord_Nodes(1,i_aus)
            i_aus = Tdomain%specel(n)%Control_Nodes(1); coord(0,1) = Tdomain%Coord_Nodes(0,i_aus); coord(1,1) = Tdomain%Coord_Nodes(1,i_aus)
            i_aus = Tdomain%specel(n)%Control_Nodes(2); coord(0,2) = Tdomain%Coord_Nodes(0,i_aus); coord(1,2) = Tdomain%Coord_Nodes(1,i_aus)
            i_aus = Tdomain%specel(n)%Control_Nodes(3); coord(0,3) = Tdomain%Coord_Nodes(0,i_aus); coord(1,3) = Tdomain%Coord_Nodes(1,i_aus)

            mat = Tdomain%specel(n)%mat_index
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz

            allocate (Tdomain%specel(n)%Jacob(0:ngllx-1,0:ngllz-1) )
            allocate (Tdomain%specel(n)%InvGrad(0:ngllx-1,0:ngllz-1,0:1,0:1) )

            do j = 0,ngllz - 1
                eta = Tdomain%sSubdomain(mat)%GLLcz(j)
                do i = 0,ngllx - 1
                    xi = Tdomain%sSubdomain(mat)%GLLcx(i)

                    ! Computation of global coordinates

                    call shape4_local2global(coord, xi, eta, xp, zp)
                    ipoint = Tdomain%specel(n)%Iglobnum(i,j)
                    Tdomain%GlobCoord(0,ipoint) = xp; Tdomain%GlobCoord(1,ipoint) = zp

                    ! Computation of the derivative matrix

                    call shape4_local2jacob(coord, xi, eta, LocInvGrad)
                    call invert2(LocInvGrad, Jac)
                    Tdomain%specel(n)%InvGrad(i,j,0:1,0:1) = LocInvGrad(0:1,0:1)
                    Tdomain%specel(n)%Jacob(i,j) = Jac
                enddo
            enddo

            ! Computation of coefficients of integration for element boundaries
            allocate(Tdomain%specel(n)%Coeff_Integr_Faces(0:2*(ngllx+ngllz)-1))

            ! Bottom boundary :
            imin = 0 ; imax = ngllx-1
            nf = Tdomain%specel(n)%Near_Face(0)
            call compute_Jacobian_1D(Tdomain,nf,Jac1D)
            Tdomain%specel(n)%Coeff_Integr_Faces(imin:imax) = Tdomain%sSubdomain(mat)%GLLwx(:) * Jac1D
            if (Tdomain%sFace(nf)%Near_Element(0) == n) then
                Tdomain%sFace(nf)%Coeff_Integr_Ends(0) = Tdomain%specel(n)%Coeff_Integr_Faces(imin)
                Tdomain%sFace(nf)%Coeff_Integr_Ends(1) = Tdomain%specel(n)%Coeff_Integr_Faces(imax)
            endif

            ! Right boundary :
            imin = imax+1 ; imax = imax + ngllz
            nf = Tdomain%specel(n)%Near_Face(1)
            call compute_Jacobian_1D(Tdomain,nf,Jac1D)
            Tdomain%specel(n)%Coeff_Integr_Faces(imin:imax) = Tdomain%sSubdomain(mat)%GLLwz(:) * Jac1D
            if (Tdomain%sFace(nf)%Near_Element(0) == n) then
                Tdomain%sFace(nf)%Coeff_Integr_Ends(0) = Tdomain%specel(n)%Coeff_Integr_Faces(imin)
                Tdomain%sFace(nf)%Coeff_Integr_Ends(1) = Tdomain%specel(n)%Coeff_Integr_Faces(imax)
            endif

            ! Top boundary :
            imin = imax+1 ; imax = imax + ngllx
            nf = Tdomain%specel(n)%Near_Face(2)
            call compute_Jacobian_1D(Tdomain,nf,Jac1D)
            Tdomain%specel(n)%Coeff_Integr_Faces(imin:imax) = Tdomain%sSubdomain(mat)%GLLwx(:) * Jac1D
            if (Tdomain%sFace(nf)%Near_Element(0) == n) then
                Tdomain%sFace(nf)%Coeff_Integr_Ends(0) = Tdomain%specel(n)%Coeff_Integr_Faces(imin)
                Tdomain%sFace(nf)%Coeff_Integr_Ends(1) = Tdomain%specel(n)%Coeff_Integr_Faces(imax)
            endif

            ! Left boundary :
            imin = imax+1 ; imax = imax + ngllz
            nf = Tdomain%specel(n)%Near_Face(3)
            call compute_Jacobian_1D(Tdomain,nf,Jac1D)
            Tdomain%specel(n)%Coeff_Integr_Faces(imin:imax) = Tdomain%sSubdomain(mat)%GLLwz(:) * Jac1D
            if (Tdomain%sFace(nf)%Near_Element(0) == n) then
                Tdomain%sFace(nf)%Coeff_Integr_Ends(0) = Tdomain%specel(n)%Coeff_Integr_Faces(imin)
                Tdomain%sFace(nf)%Coeff_Integr_Ends(1) = Tdomain%specel(n)%Coeff_Integr_Faces(imax)
            endif
        enddo

        ! Compute Normals for Faces with DG :
        do nf=0,Tdomain%n_face-1
            if (Tdomain%sFace(nf)%Type_flux .NE. FLUX_NONE) then
                call compute_normals(Tdomain,nf)
            end if
        end do


        ! Compute Normals for Elements HDG :
        do n=0,Tdomain%n_elem-1
            !if (Tdomain%specel(n)%Type_DG .EQ. GALERKIN_HDG_RP) then
                call compute_normals_elem(Tdomain,n)
            !end if
        end do

        if (Tdomain%logicD%super_object_local_present) then
            do n = 0, Tdomain%n_fault-1
                call shape4_manage_super_object(Tdomain, n)
            enddo
        endif
        return
    end subroutine shape4

    subroutine shape4_manage_super_object(Tdomain, n)
        use sdomain
        implicit none
        type(domain),target, intent (INOUT) :: Tdomain
        integer, intent(IN) :: n
        integer :: i, j, i_aus
        real(fpp), dimension (:), pointer :: GLLc_face
        real(fpp), dimension (:,:), allocatable :: Store_normal
        real(fpp) :: ds_local, normal_0, normal_1, normalization
        real(fpp) :: x0, x1, z0, z1
        integer :: Face_up, mat, ngll, nv, nv2

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
                i_aus = Tdomain%specel(nv2)%Control_Nodes(0)
                x0 = Tdomain%Coord_Nodes(0,i_aus)
                z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(1)
                x1 = Tdomain%Coord_Nodes(0,i_aus)
                z1 = Tdomain%Coord_Nodes(1,i_aus)
            else if (Tdomain%specel(nv2)%near_face(1) == Face_up) then
                i_aus = Tdomain%specel(nv2)%Control_Nodes(1)
                x0 = Tdomain%Coord_Nodes(0,i_aus)
                z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(2)
                x1 = Tdomain%Coord_Nodes(0,i_aus)
                z1 = Tdomain%Coord_Nodes(1,i_aus)
            else if (Tdomain%specel(nv2)%near_face(2) == Face_up) then
                i_aus = Tdomain%specel(nv2)%Control_Nodes(2)
                x0 = Tdomain%Coord_Nodes(0,i_aus)
                z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(3)
                x1 = Tdomain%Coord_Nodes(0,i_aus)
                z1 = Tdomain%Coord_Nodes(1,i_aus)
            else if (Tdomain%specel(nv2)%near_face(3) == Face_up) then
                i_aus = Tdomain%specel(nv2)%Control_Nodes(3)
                x0 = Tdomain%Coord_Nodes(0,i_aus)
                z0 = Tdomain%Coord_Nodes(1,i_aus)
                i_aus = Tdomain%specel(nv2)%Control_Nodes(0)
                x1 = Tdomain%Coord_Nodes(0,i_aus)
                z1 = Tdomain%Coord_Nodes(1,i_aus)
            else
                stop "Inconsistency with Face_up"
            endif

            allocate (Tdomain%sFault(n)%fFace(j)%ds(0:ngll-1))
            allocate (Tdomain%sFault(n)%fFace(j)%normal(0:ngll-1,0:1))
            allocate (Tdomain%sFault(n)%fFace(j)%distance(0:ngll-1))

            ds_local =  (x1-x0)**2 + (z1-z0)**2
            ds_local =  Sqrt(ds_local)
            normal_0 = (z1 - z0)/ds_local; normal_1 = (x0-x1)/ds_local

            do  i = 0, ngll-1
                Tdomain%sFault(n)%fFace(j)%distance(i) = ds_local * (1+GLLc_face(i))
            enddo

            Tdomain%sFault(n)%fFace(j)%X_Vertex(0) = x0; Tdomain%sFault(n)%fFace(j)%X_Vertex(1) = x1
            Tdomain%sFault(n)%fFace(j)%Z_Vertex(0) = z0; Tdomain%sFault(n)%fFace(j)%Z_Vertex(1) = z1
            Tdomain%sFault(n)%fFace(j)%ds(:) = 0.5 * ds_local
            Tdomain%sFault(n)%fFace(j)%normal(0:ngll-1,0) = normal_0; Tdomain%sFault(n)%fFace(j)%normal(0:ngll-1,1) = normal_1
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
        enddo;

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
            ds_local = ds_local + (Tdomain%sFault(n)%fFace(j)%distance(ngll-1) - Tdomain%sFault(n)%fFace(j)%distance(0))
            nv = Tdomain%sFault(n)%fFace(j)%Face_to_vertex(0)
            Tdomain%sFault(n)%fvertex(nv)%distance = Tdomain%sFault(n)%fFace(j)%distance(0)
            nv = Tdomain%sFault(n)%fFace(j)%Face_to_vertex(1)
            Tdomain%sFault(n)%fvertex(nv)%distance = Tdomain%sFault(n)%fFace(j)%distance(ngll-1)
            allocate (Store_normal(0:ngll-1,0:0))
            Store_normal(:,0) =  Tdomain%sFault(n)%fFace(j)%distance
            deallocate ( Tdomain%sFault(n)%fFace(j)%distance)
            allocate ( Tdomain%sFault(n)%fFace(j)%distance(1:ngll-2))
            Tdomain%sFault(n)%fFace(j)%distance(1:ngll-2)= Store_normal (1:ngll-2,0)
            deallocate (Store_normal)
        enddo

    end subroutine shape4_manage_super_object

    subroutine shape4_global2local(coord, xa, za, xi, eta, ok)
        real(fpp), dimension(0:1,0:3), intent(in)  :: coord
        real(fpp), intent(in) :: xa, za
        real(fpp), intent(out) :: xi, eta
        logical, intent(out) :: ok
        !
        integer :: niter
        real(fpp), dimension(0:1) :: xin, xout, xref
        ok = .true.
        xin(0) = 0.
        xin(1) = 0.
        xref(0) = xa
        xref(1) = za
        call shape4_simple_newton(coord, xref, xin, xout, niter)
        if (niter==1000 .or. niter<0) ok=.false.
        xi  = xout(0)
        eta = xout(1)
    end subroutine shape4_global2local

    subroutine shape4_simple_newton(nodes, xref, xin, xout, nit)
        real(fpp), dimension(0:1,0:3), intent(in) :: nodes
        real(fpp), dimension(0:1), intent(in) :: xref, xin
        real(fpp), dimension(0:1), intent(out) :: xout
        integer, intent(out) :: nit
        !
        real(fpp), dimension(0:1,0:1) :: jac
        real(fpp), dimension(0:1) :: x
        real(fpp) :: xa, za, err, Det
        integer, parameter :: niter=1000
        integer :: i
        xout = xin
        do i=1,niter
            call shape4_local2global(nodes, xout(0), xout(1), xa, za)
            call shape4_local2jacob (nodes, xout(0), xout(1), Jac)
            call invert2 (Jac, Det)
            xa = xref(0)-xa
            za = xref(1)-za
            x(0) = Jac(0,0)*xa + Jac(1,0)*za
            x(1) = Jac(0,1)*xa + Jac(1,1)*za
            err = x(0)**2+x(1)**2
            if (err<1e-24) exit ! eps^2 with eps=1e-12 : need high precision to recover previous results (basics tests)
            xout = xout + x
        end do
        nit = i
    end subroutine shape4_simple_newton

    subroutine shape4_local2global(coord, xi, eta, xa, za)
        real(fpp), dimension(0:1,0:3), intent(in) :: coord
        real(fpp), intent(in) :: xi, eta
        real(fpp), intent(out) :: xa, za
        !
        real(fpp) x0, x1, x2, x3, z0, z1, z2, z3 ! Used for clarity (-O2 should be able to optimize this)
        x0 = coord(0,0); x1 = coord (0,1); x2 = coord (0,2); x3 = coord (0,3);
        z0 = coord(1,0); z1 = coord (1,1); z2 = coord (1,2); z3 = coord (1,3);
        !- computation of the global coordinates from the local coordinates
        xa = 0.25*(x0 * (1.-xi)*(1.-eta) + x1 * (1.+xi)*(1.-eta) + x2 * (1.+xi)*(1.+eta) + x3 * (1.-xi)*(1.+eta))
        za = 0.25*(z0 * (1.-xi)*(1.-eta) + z1 * (1.+xi)*(1.-eta) + z2 * (1.+xi)*(1.+eta) + z3 * (1.-xi)*(1.+eta))
    end subroutine shape4_local2global

    subroutine shape4_local2jacob(coord, xi, eta, jac)
        real(fpp), dimension(0:1,0:3), intent(in)  :: coord
        real(fpp), intent(in) :: xi, eta
        real(fpp), dimension(0:1,0:1), intent(out) :: jac
        !
        real(fpp) x0, x1, x2, x3, z0, z1, z2, z3 ! Used for clarity (-O2 should be able to optimize this)
        x0 = coord(0,0); x1 = coord (0,1); x2 = coord (0,2); x3 = coord (0,3);
        z0 = coord(1,0); z1 = coord (1,1); z2 = coord (1,2); z3 = coord (1,3);
        !- computation of the derivative matrix, dx_(jj)/dxi_(ii)
        jac(0,0) = 0.25*((x1-x0) * (1-eta) + (x2-x3) * (1+eta))
        jac(1,0) = 0.25*((x3-x0) * (1-xi)  + (x2-x1) * (1+xi) )
        jac(0,1) = 0.25*((z1-z0) * (1-eta) + (z2-z3) * (1+eta))
        jac(1,1) = 0.25*((z3-z0) * (1-xi)  + (z2-z1) * (1+xi) )
    end subroutine shape4_local2jacob

end module shape_lin

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
