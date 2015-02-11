!>
!! \file shape27.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<
module mshape27
    implicit none
contains
    !>
    !! shape27: alloue et calcule la jacobienne et l'inverse du gradient de ??
    !<
    subroutine shape27_init(Tdomain)

        ! Modified by Paul Cupillard 26/06/2006


        use sdomain

        implicit none

        type(domain), target, intent (INOUT) :: Tdomain

        integer :: n, ngllx,nglly,ngllz, mat, i,j,k, ipoint
        real :: xi,eta,zeta, xp,yp,zp, Jac
        real, dimension(0:2,0:26) :: coord
        real, dimension(0:2,0:2) :: LocInvGrad


        allocate (Tdomain%GlobCoord(0:2,0:Tdomain%n_glob_points-1))

        do n = 0,Tdomain%n_elem - 1

            do i = 0,26
                j = Tdomain%specel(n)%Control_Nodes(i)
                coord(0:2,i) = Tdomain%Coord_Nodes(0:2,j)
            enddo

            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            mat   = Tdomain%specel(n)%mat_index

            allocate (Tdomain%specel(n)%Jacob(0:ngllx-1,0:nglly-1,0:ngllz-1) )
            allocate (Tdomain%specel(n)%InvGrad(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2) )

            do k = 0,ngllz - 1
                zeta =  Tdomain%sSubdomain(mat)%GLLcz (k)
                do j = 0,nglly - 1
                    eta =  Tdomain%sSubdomain(mat)%GLLcy (j)
                    do i = 0,ngllx - 1
                        xi = Tdomain%sSubdomain(mat)%GLLcx (i)
                        call shape27_local2global(coord, xi, eta, zeta, xp, yp, zp)
                        ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)
                        Tdomain%GlobCoord (0,ipoint) = xp
                        Tdomain%GlobCoord (1,ipoint) = yp
                        Tdomain%GlobCoord (2,ipoint) = zp

                        call shape27_local2jacob(coord, xi, eta, zeta, LocInvGrad)
                        call invert_3d (LocInvGrad, Jac)
                        Tdomain%specel(n)%Jacob(i,j,k) = Jac
                        Tdomain%specel(n)%InvGrad (i,j,k,0:2,0:2) = LocInvGrad(0:2,0:2)
                    enddo
                enddo
            enddo
        enddo
    end subroutine shape27_init
    !---------------------------------------------------------------------------
    real function shape27_func(which_nod,xi,eta,zeta)

        integer :: which_nod
        real :: xi,eta,zeta

        shape27_func = 0.0
        select case (which_nod)
        case (0)
            shape27_func = 0.125 * xi*(xi-1) * eta*(eta-1) * zeta*(zeta-1)
        case (1)
            shape27_func = 0.125 * xi*(xi+1) * eta*(eta-1) * zeta*(zeta-1)
        case (2)
            shape27_func = 0.125 * xi*(xi+1) * eta*(eta+1) * zeta*(zeta-1)
        case (3)
            shape27_func = 0.125 * xi*(xi-1) * eta*(eta+1) * zeta*(zeta-1)
        case (4)
            shape27_func = 0.125 * xi*(xi-1) * eta*(eta-1) * zeta*(zeta+1)
        case (5)
            shape27_func = 0.125 * xi*(xi+1) * eta*(eta-1) * zeta*(zeta+1)
        case (6)
            shape27_func = 0.125 * xi*(xi+1) * eta*(eta+1) * zeta*(zeta+1)
        case (7)
            shape27_func = 0.125 * xi*(xi-1) * eta*(eta+1) * zeta*(zeta+1)
        case (8)
            shape27_func = 0.25 * (1-xi**2) * eta*(eta-1) * zeta*(zeta-1)
        case (9)
            shape27_func = 0.25 * xi*(xi+1) * (1-eta**2) * zeta*(zeta-1)
        case (10)
            shape27_func = 0.25 * (1-xi**2) * eta*(eta+1) * zeta*(zeta-1)
        case (11)
            shape27_func = 0.25 * xi*(xi-1) * (1-eta**2) * zeta*(zeta-1)
        case (12)
            shape27_func = 0.25 * xi*(xi-1) * eta*(eta-1) * (1-zeta**2)
        case (13)
            shape27_func = 0.25 * xi*(xi+1) * eta*(eta-1) * (1-zeta**2)
        case (14)
            shape27_func = 0.25 * xi*(xi+1) * eta*(eta+1) * (1-zeta**2)
        case (15)
            shape27_func = 0.25 * xi*(xi-1) * eta*(eta+1) * (1-zeta**2)
        case (16)
            shape27_func = 0.25 * (1-xi**2) * eta*(eta-1) * zeta*(zeta+1)
        case (17)
            shape27_func = 0.25 * xi*(xi+1) * (1-eta**2) * zeta*(zeta+1)
        case (18)
            shape27_func = 0.25 * (1-xi**2) * eta*(eta+1) * zeta*(zeta+1)
        case (19)
            shape27_func = 0.25 * xi*(xi-1) * (1-eta**2) * zeta*(zeta+1)
        case(20)
            shape27_func = 0.5 * (1-xi**2) * (1-eta**2) * zeta*(zeta-1)
        case(21)
            shape27_func = 0.5 * (1-xi**2) * eta*(eta-1) * (1-zeta**2)
        case(22)
            shape27_func = 0.5 * xi*(xi+1) * (1-eta**2) * (1-zeta**2)
        case(23)
            shape27_func = 0.5 * (1-xi**2) * eta*(eta+1) * (1-zeta**2)
        case(24)
            shape27_func = 0.5 * xi*(xi-1) * (1-eta**2) * (1-zeta**2)
        case(25)
            shape27_func = 0.5 * (1-xi**2) * (1-eta**2) * zeta*(zeta+1)
        case(26)
            shape27_func = (1-xi**2) * (1-eta**2) * (1-zeta**2)
        end select

        return
    end function shape27_func
    !---------------------------------------------------------------------------
    real function shape27_derivfunc(which_nod,xi,eta,zeta,compo)

        integer :: which_nod, compo
        real :: xi,eta,zeta

        real, dimension(0:2) :: df

        select case (which_nod)
        case (0)
            df(0) = 0.125 * (2*xi-1) * eta*(eta-1) * zeta*(zeta-1)
            df(1) = 0.125 * xi*(xi-1) * (2*eta-1) * zeta*(zeta-1)
            df(2) = 0.125 * xi*(xi-1) * eta*(eta-1) * (2*zeta-1)
        case (1)
            df(0) = 0.125 * (2*xi+1) * eta*(eta-1) * zeta*(zeta-1)
            df(1) = 0.125 * xi*(xi+1) * (2*eta-1) * zeta*(zeta-1)
            df(2) = 0.125 * xi*(xi+1) * eta*(eta-1) * (2*zeta-1)
        case (2)
            df(0) = 0.125 * (2*xi+1) * eta*(eta+1) * zeta*(zeta-1)
            df(1) = 0.125 * xi*(xi+1) * (2*eta+1) * zeta*(zeta-1)
            df(2) = 0.125 * xi*(xi+1) * eta*(eta+1) * (2*zeta-1)
        case (3)
            df(0) = 0.125 * (2*xi-1) * eta*(eta+1) * zeta*(zeta-1)
            df(1) = 0.125 * xi*(xi-1) * (2*eta+1) * zeta*(zeta-1)
            df(2) = 0.125 * xi*(xi-1) * eta*(eta+1) * (2*zeta-1)
        case (4)
            df(0) = 0.125 * (2*xi-1) * eta*(eta-1) * zeta*(zeta+1)
            df(1) = 0.125 * xi*(xi-1) * (2*eta-1) * zeta*(zeta+1)
            df(2) = 0.125 * xi*(xi-1) * eta*(eta-1) * (2*zeta+1)
        case (5)
            df(0) = 0.125 * (2*xi+1) * eta*(eta-1) * zeta*(zeta+1)
            df(1) = 0.125 * xi*(xi+1) * (2*eta-1) * zeta*(zeta+1)
            df(2) = 0.125 * xi*(xi+1) * eta*(eta-1) * (2*zeta+1)
        case (6)
            df(0) = 0.125 * (2*xi+1) * eta*(eta+1) * zeta*(zeta+1)
            df(1) = 0.125 * xi*(xi+1) * (2*eta+1) * zeta*(zeta+1)
            df(2) = 0.125 * xi*(xi+1) * eta*(eta+1) * (2*zeta+1)
        case (7)
            df(0) = 0.125 * (2*xi-1) * eta*(eta+1) * zeta*(zeta+1)
            df(1) = 0.125 * xi*(xi-1) * (2*eta+1) * zeta*(zeta+1)
            df(2) = 0.125 * xi*(xi-1) * eta*(eta+1) * (2*zeta+1)
        case (8)
            df(0) = 0.25 * (-2*xi) * eta*(eta-1) * zeta*(zeta-1)
            df(1) = 0.25 * (1-xi**2) * (2*eta-1) * zeta*(zeta-1)
            df(2) = 0.25 * (1-xi**2) * eta*(eta-1) * (2*zeta-1)
        case (9)
            df(0) = 0.25 * (2*xi+1) * (1-eta**2) * zeta*(zeta-1)
            df(1) = 0.25 * xi*(xi+1) * (-2*eta) * zeta*(zeta-1)
            df(2) = 0.25 * xi*(xi+1) * (1-eta**2) * (2*zeta-1)
        case (10)
            df(0) = 0.25 * (-2*xi) * eta*(eta+1) * zeta*(zeta-1)
            df(1) = 0.25 * (1-xi**2) * (2*eta+1) * zeta*(zeta-1)
            df(2) = 0.25 * (1-xi**2) * eta*(eta+1) * (2*zeta-1)
        case (11)
            df(0) = 0.25 * (2*xi-1) * (1-eta**2) * zeta*(zeta-1)
            df(1) = 0.25 * xi*(xi-1) * (-2*eta) * zeta*(zeta-1)
            df(2) = 0.25 * xi*(xi-1) * (1-eta**2) * (2*zeta-1)
        case (12)
            df(0) = 0.25 * (2*xi-1) * eta*(eta-1) * (1-zeta**2)
            df(1) = 0.25 * xi*(xi-1) * (2*eta-1) * (1-zeta**2)
            df(2) = 0.25 * xi*(xi-1) * eta*(eta-1) * (-2*zeta)
        case (13)
            df(0) = 0.25 * (2*xi+1) * eta*(eta-1) * (1-zeta**2)
            df(1) = 0.25 * xi*(xi+1) * (2*eta-1) * (1-zeta**2)
            df(2) = 0.25 * xi*(xi+1) * eta*(eta-1) * (-2*zeta)
        case (14)
            df(0) = 0.25 * (2*xi+1) * eta*(eta+1) * (1-zeta**2)
            df(1) = 0.25 * xi*(xi+1) * (2*eta+1) * (1-zeta**2)
            df(2) = 0.25 * xi*(xi+1) * eta*(eta+1) * (-2*zeta)
        case (15)
            df(0) = 0.25 * (2*xi-1) * eta*(eta+1) * (1-zeta**2)
            df(1) = 0.25 * xi*(xi-1) * (2*eta+1) * (1-zeta**2)
            df(2) = 0.25 * xi*(xi-1) * eta*(eta+1) * (-2*zeta)
        case (16)
            df(0) = 0.25 * (-2*xi) * eta*(eta-1) * zeta*(zeta+1)
            df(1) = 0.25 * (1-xi**2) * (2*eta-1) * zeta*(zeta+1)
            df(2) = 0.25 * (1-xi**2) * eta*(eta-1) * (2*zeta+1)
        case (17)
            df(0) = 0.25 * (2*xi+1) * (1-eta**2) * zeta*(zeta+1)
            df(1) = 0.25 * xi*(xi+1) * (-2*eta) * zeta*(zeta+1)
            df(2) = 0.25 * xi*(xi+1) * (1-eta**2) * (2*zeta+1)
        case (18)
            df(0) = 0.25 * (-2*xi) * eta*(eta+1) * zeta*(zeta+1)
            df(1) = 0.25 * (1-xi**2) * (2*eta+1) * zeta*(zeta+1)
            df(2) = 0.25 * (1-xi**2) * eta*(eta+1) * (2*zeta+1)
        case (19)
            df(0) = 0.25 * (2*xi-1) * (1-eta**2) * zeta*(zeta+1)
            df(1) = 0.25 * xi*(xi-1) * (-2*eta) * zeta*(zeta+1)
            df(2) = 0.25 * xi*(xi-1) * (1-eta**2) * (2*zeta+1)
        case(20)
            df(0) = 0.5 * (-2*xi) * (1-eta**2) * zeta*(zeta-1)
            df(1) = 0.5 * (1-xi**2) * (-2*eta) * zeta*(zeta-1)
            df(2) = 0.5 * (1-xi**2) * (1-eta**2) * (2*zeta-1)
        case(21)
            df(0) = 0.5 * (-2*xi) * eta*(eta-1) * (1-zeta**2)
            df(1) = 0.5 * (1-xi**2) * (2*eta-1) * (1-zeta**2)
            df(2) = 0.5 * (1-xi**2) * eta*(eta-1) * (-2*zeta)
        case(22)
            df(0) = 0.5 * (2*xi+1) * (1-eta**2) * (1-zeta**2)
            df(1) = 0.5 * xi*(xi+1) * (-2*eta) * (1-zeta**2)
            df(2) = 0.5 * xi*(xi+1) * (1-eta**2) * (-2*zeta)
        case(23)
            df(0) = 0.5 * (-2*xi) * eta*(eta+1) * (1-zeta**2)
            df(1) = 0.5 * (1-xi**2) * (2*eta+1) * (1-zeta**2)
            df(2) = 0.5 * (1-xi**2) * eta*(eta+1) * (-2*zeta)
        case(24)
            df(0) = 0.5 * (2*xi-1) * (1-eta**2) * (1-zeta**2)
            df(1) = 0.5 * xi*(xi-1) * (-2*eta) * (1-zeta**2)
            df(2) = 0.5 * xi*(xi-1) * (1-eta**2) * (-2*zeta)
        case(25)
            df(0) = 0.5 * (-2*xi) * (1-eta**2) * zeta*(zeta+1)
            df(1) = 0.5 * (1-xi**2) * (-2*eta) * zeta*(zeta+1)
            df(2) = 0.5 * (1-xi**2) * (1-eta**2) * (2*zeta+1)
        case(26)
            df(0) = (-2*xi) * (1-eta**2) * (1-zeta**2)
            df(1) = (1-xi**2) * (-2*eta) * (1-zeta**2)
            df(2) = (1-xi**2) * (1-eta**2) * (-2*zeta)
        end select

        shape27_derivfunc = df(compo)

        return
    end function shape27_derivfunc
    !---------------------------------------------------------------------------
    subroutine shape27_local2global(coord, xi, eta, zeta, xa, ya, za)
        real, dimension(0:2,0:26), intent(in) :: coord
        real, intent(in) :: xi, eta, zeta
        real, intent(out) :: xa, ya, za
        !
        integer :: i
        double precision :: f
        !
        xa = 0;   ya = 0;   za = 0;
        do i = 0,26
            f = shape27_func(i,xi,eta,zeta)
            xa = xa + coord(0,i)*f
            ya = ya + coord(1,i)*f
            za = za + coord(2,i)*f
        enddo
    end subroutine shape27_local2global
    !---------------------------------------------------------------------------
    subroutine shape27_local2jacob(coord, xi, eta, zeta, jac)
        double precision, dimension(0:2,0:7), intent(in)  :: coord
        double precision, intent(in) :: xi, eta, zeta
        double precision, dimension(0:2,0:2), intent(out) :: jac
        !
        double precision :: f
        integer :: i, j
        !- computation of the derivative matrix, J(i,j)= dx_(jj)/dxi_(ii)
        jac = 0.
        do i = 0,26
            do j = 0,2
                f = shape27_derivfunc(i,xi,eta,zeta,j)
                jac(j,0) = jac(j,0) + coord(0,i)*f
                jac(j,1) = jac(j,1) + coord(1,i)*f
                jac(j,2) = jac(j,2) + coord(2,i)*f
            enddo
        enddo
    end subroutine shape27_local2jacob
    !---------------------------------------------------------------------------
    double precision function shape27_min(dim, nn, x, nodes, xref)
        integer, intent(in) :: dim, nn
        double precision, dimension(0:dim-1), intent(in) :: x, xref
        double precision, dimension(0:dim-1,0:nn-1), intent(in) :: nodes
        double precision :: xa, ya, za, xrx, xry, xrz
        call shape27_local2global(nodes, x(0), x(1), x(2), xa, ya, za)
        xrx = 1d0/max(abs(xref(0)),1d0)
        xry = 1d0/max(abs(xref(1)),1d0)
        xrz = 1d0/max(abs(xref(2)),1d0)
        shape27_min =  (xrx*(xa-xref(0)))**2 + (xry*(ya-xref(1)))**2 + (xrz*(za-xref(2)))**2
 !       write(*,*) "minfun=", shape27_min, " at ", x(0), x(1), x(2)
    end function shape27_min
    !---------------------------------------------------------------------------
    subroutine shape27_mingrad(dim, nn, x, nodes, xref, grad)
        integer, intent(in) :: dim, nn
        double precision, dimension(0:dim-1), intent(in) :: x, xref
        double precision, dimension(0:dim-1,0:nn-1), intent(in) :: nodes
        double precision, dimension(0:dim-1), intent(out) :: grad
        double precision, dimension(0:dim-1,0:dim-1) :: jac
        double precision :: xa, ya, za, xrx, xry, xrz, fval
        call shape27_local2global(nodes, x(0), x(1), x(2), xa, ya, za)
        call shape27_local2jacob(nodes, x(0), x(1), x(2), jac)
        xrx = 1d0/max(xref(0)**2,1d0)
        xry = 1d0/max(xref(1)**2,1d0)
        xrz = 1d0/max(xref(2)**2,1d0)
        grad(0) = 2*(jac(0,0)*xrx*(xa-xref(0))+jac(0,1)*xry*(ya-xref(1))+jac(0,2)*xrz*(za-xref(2)))
        grad(1) = 2*(jac(1,0)*xrx*(xa-xref(0))+jac(1,1)*xry*(ya-xref(1))+jac(1,2)*xrz*(za-xref(2)))
        grad(2) = 2*(jac(2,0)*xrx*(xa-xref(0))+jac(2,1)*xry*(ya-xref(1))+jac(2,2)*xrz*(za-xref(2)))
!        write(*,*) "dminfun/dx=", grad(0)
!        write(*,*) "dminfun/dy=", grad(1)
!        write(*,*) "dminfun/dz=", grad(2)
    end subroutine shape27_mingrad
    !---------------------------------------------------------------------------
    subroutine simple_newton_27(nodes, xref, xin, xout, nit)
        double precision, dimension(0:2), intent(in) :: xref, xin
        double precision, dimension(0:2), intent(out) :: xout
        integer, intent(out) :: nit
        double precision, dimension(0:2,0:26), intent(in) :: nodes
        double precision, dimension(0:2,0:2) :: jac
        double precision, dimension(0:2) :: x
        double precision :: xa, ya, za, err, Det
        integer, parameter :: niter=1000
        integer :: i
        xout = xin
        do i=1,niter
            call shape27_local2global(nodes, xout(0), xout(1), xout(2), xa, ya, za)
            call shape27_local2jacob(nodes, xout(0), xout(1), xout(2), Jac)
            call invert_3d (Jac, Det)
            xa = xref(0)-xa
            ya = xref(1)-ya
            za = xref(2)-za
            x(0) = Jac(0,0)*xa + Jac(1,0)*ya + Jac(2,0)*za
            x(1) = Jac(0,1)*xa + Jac(1,1)*ya + Jac(2,1)*za
            x(2) = Jac(0,2)*xa + Jac(1,2)*ya + Jac(2,2)*za
            err = x(0)**2+x(1)**2+x(2)**2
            if (err<1e-12) exit
            xout = xout + x
        end do
        nit = i
    end subroutine simple_newton_27
    !---------------------------------------------------------------------------
    subroutine shape27_global2local(coord, xa, ya, za, xi, eta, zeta, ok)
        use mleastsq
        double precision, dimension(0:2,0:26), intent(in)  :: coord
        double precision, intent(in) :: xa, ya, za
        double precision, intent(out) :: xi, eta, zeta
        logical, intent(out) :: ok
        !
        integer :: niter
        double precision, dimension(0:2) :: xin, xout, xref
        ok = .true.
        xin(0) = 0.
        xin(1) = 0.
        xin(2) = 0.
        xref(0) = xa
        xref(1) = ya
        xref(2) = za
        !call minimize_cg(3, 27, xin, coord, xref, shape27_min, shape27_mingrad, 0.1D0, xout, niter)
        call simple_newton_27(coord, xref, xin, xout, niter)
        if (niter==1000 .or. niter<0) ok=.false.
        xi = xout(0)
        eta = xout(1)
        zeta = xout(2)
    end subroutine shape27_global2local

end module mshape27
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
