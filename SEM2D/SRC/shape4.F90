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


contains

    ! #########################################

    subroutine shape4(Tdomain)

        use sdomain
        use semdatafiles

        implicit none

        type(domain),target, intent (INOUT) :: Tdomain

        ! local variables

        integer :: i_aus, n, mat, ngllx, ngllz, i, j, ipoint
        real, dimension(0:1,0:3) :: coord
        real :: xi, eta, xp, zp, Jac
        real, dimension (0:1,0:1) :: LocInvGrad

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
        enddo

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
        real, dimension (:), pointer :: GLLc_face
        real, dimension (:,:), allocatable :: Store_normal
        real :: ds_local, normal_0, normal_1, normalization
        real :: x0, x1, z0, z1
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
        double precision, dimension(0:1,0:3), intent(in)  :: coord
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
        call shape4_simple_newton(coord, xref, xin, xout, niter)
        if (niter==1000 .or. niter<0) ok=.false.
        xi  = xout(0)
        eta = xout(1)
    end subroutine shape4_global2local

    subroutine shape4_simple_newton(nodes, xref, xin, xout, nit)
        double precision, dimension(0:1,0:3), intent(in) :: nodes
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
        double precision, dimension(0:1,0:3), intent(in) :: coord
        double precision, intent(in) :: xi, eta
        double precision, intent(out) :: xa, za
        !
        double precision x0, x1, x2, x3, z0, z1, z2, z3 ! Used for clarity (-O2 should be able to optimize this)
        x0 = coord(0,0); x1 = coord (0,1); x2 = coord (0,2); x3 = coord (0,3);
        z0 = coord(1,0); z1 = coord (1,1); z2 = coord (1,2); z3 = coord (1,3);
        !- computation of the global coordinates from the local coordinates
        xa = 0.25*(x0 * (1.-xi)*(1.-eta) + x1 * (1.+xi)*(1.-eta) + x2 * (1.+xi)*(1.+eta) + x3 * (1.-xi)*(1.+eta))
        za = 0.25*(z0 * (1.-xi)*(1.-eta) + z1 * (1.+xi)*(1.-eta) + z2 * (1.+xi)*(1.+eta) + z3 * (1.-xi)*(1.+eta))
    end subroutine shape4_local2global

    subroutine shape4_local2jacob(coord, xi, eta, jac)
        double precision, dimension(0:1,0:3), intent(in)  :: coord
        double precision, intent(in) :: xi, eta
        double precision, dimension(0:1,0:1), intent(out) :: jac
        !
        double precision x0, x1, x2, x3, z0, z1, z2, z3 ! Used for clarity (-O2 should be able to optimize this)
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
