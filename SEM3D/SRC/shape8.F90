! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file shape8.f90
!!\brief Contient la subroutine shape8.
!!
!<
module mshape8
    use constants, only : fpp
    use sdomain
    implicit none
#include "index.h"

contains
    subroutine shape8_init(Tdomain)
        use scomm
        use mrenumber
        use surface_load, only : get_surf_gll_coord
        use DimensionalShape, only : surfaceSource
        implicit none

        type(domain), intent(inout)          :: Tdomain
        integer                              :: n,ngll,i,j,k,ipoint
        real(fpp)                            :: Jac, xi, eta, zeta
        real(fpp), dimension(0:2,0:2)        :: LocInvGrad
        real(fpp), dimension(0:2,0:7)        :: coord  ! coordinates of nodes
        real(fpp), dimension(:), allocatable :: GLLc
        real(fpp)                            :: xp, yp, zp
        real(fpp), parameter                 :: XEPS=1e-10
        real(fpp)                            :: dxp, dyp, dzp
        integer                              :: bnum, ee
        integer                              :: i_surf

        ! Tdomain%n_glob_points is the number of degrees of fredom
        allocate(Tdomain%GlobCoord(0:2,0:Tdomain%n_glob_points-1))
        Tdomain%GlobCoord = 0d0

        do n = 0,Tdomain%n_elem-1
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee = mod(Tdomain%specel(n)%lnum,VCHUNK)

            ! coordinates of control nodes (which are vertices also)
            call nodes_coord_8(Tdomain%specel(n)%Control_Nodes(0:),Tdomain%n_glob_nodes,    &
                Tdomain%Coord_Nodes,coord)

            ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
            call domain_gllc(Tdomain, Tdomain%specel(n)%domain, GLLc)

            ! coordinates of GLL points, and values of Jacobian and dX_dxi at each GLL point.
            do k = 0,ngll-1
                zeta = GLLc(k)
                do j = 0,ngll-1
                    eta = GLLc(j)
                    do i = 0,ngll-1
                        xi = GLLc(i)

                        ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)

                        call shape8_local2global(coord, xi, eta, zeta, xp, yp, zp)
                        if (Tdomain%GlobCoord(0,ipoint)/=0d0 .or. &
                            Tdomain%GlobCoord(1,ipoint)/=0d0 .or. &
                            Tdomain%GlobCoord(2,ipoint)/=0d0) then
                            dxp = abs(xp-Tdomain%GlobCoord(0,ipoint))
                            dyp = abs(yp-Tdomain%GlobCoord(1,ipoint))
                            dzp = abs(zp-Tdomain%GlobCoord(2,ipoint))
                            if (dxp>XEPS .or. dyp>XEPS .or. dzp>XEPS) then
                                !write(*,*) "DIFF", n, ipoint, Tdomain%GlobCoord(:,ipoint), ":", dxp, dyp, dzp
                            end if
                        end if
                        !
                        Tdomain%GlobCoord(0,ipoint) = xp
                        Tdomain%GlobCoord(1,ipoint) = yp
                        Tdomain%GlobCoord(2,ipoint) = zp

                        call shape8_local2jacob(coord, xi, eta, zeta, LocInvGrad)
                        call invert_3d(LocInvGrad,Jac)
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID)
                                Tdomain%sdom%Jacob_     (        i,j,k,bnum,ee) = Jac
                                Tdomain%sdom%InvGrad_   (0:2,0:2,i,j,k,bnum,ee) = LocInvGrad(0:2,0:2)
                            case (DM_FLUID)
                                Tdomain%fdom%Jacob_     (        i,j,k,bnum,ee) = Jac
                                Tdomain%fdom%InvGrad_   (0:2,0:2,i,j,k,bnum,ee) = LocInvGrad(0:2,0:2)
                            case (DM_SOLID_PML)
                                Tdomain%spmldom%Jacob_  (        i,j,k,bnum,ee) = Jac
                                Tdomain%spmldom%InvGrad_(0:2,0:2,i,j,k,bnum,ee) = LocInvGrad(0:2,0:2)
                            case (DM_FLUID_PML)
                                Tdomain%fpmldom%Jacob_  (        i,j,k,bnum,ee) = Jac
                                Tdomain%fpmldom%InvGrad_(0:2,0:2,i,j,k,bnum,ee) = LocInvGrad(0:2,0:2)
                            case default
                                stop "unknown domain"
                        end select
                    enddo
                enddo
            enddo ! end of loops onto GLL points inside an element
            deallocate(GLLc)
        enddo ! end of loop onto elements

        ! Surface Conditions : normal vectors
        if (Tdomain%logicD%surfBC) then
            do i_surf=0,size(Tdomain%sSurfaces)-1
               select case (Tdomain%sSurfaces(i_surf)%domain)
                  case (DM_SOLID)
                     call compute_normals(Tdomain, Tdomain%sSurfaces(i_surf)%surf_sl, DM_SOLID, Tdomain%sSurfaces(i_surf)%Surf_BtN)
                     call get_surf_gll_coord(Tdomain%sSurfaces(i_surf)%surf_sl,Tdomain, Tdomain%sSurfaces(i_surf)%coord)
                     call surfaceSource(Tdomain%sSurfaces(i_surf)%surf_sl%nbtot,Tdomain,Tdomain%sSurfaces(i_surf))
                  case (DM_FLUID)
                     call compute_normals(Tdomain, Tdomain%sSurfaces(i_surf)%surf_fl, DM_FLUID, Tdomain%sSurfaces(i_surf)%Surf_BtN)
                     call get_surf_gll_coord(Tdomain%sSurfaces(i_surf)%surf_fl,Tdomain, Tdomain%sSurfaces(i_surf)%coord)
                     call surfaceSource(Tdomain%sSurfaces(i_surf)%surf_fl%nbtot,Tdomain,Tdomain%sSurfaces(i_surf))
                  case (DM_SOLID_PML)
                     call compute_normals(Tdomain, Tdomain%sSurfaces(i_surf)%surf_spml, DM_SOLID_PML,Tdomain%sSurfaces(i_surf)%Surf_BtN)
                     call get_surf_gll_coord(Tdomain%sSurfaces(i_surf)%surf_spml,Tdomain, Tdomain%sSurfaces(i_surf)%coord)
                     call surfaceSource(Tdomain%sSurfaces(i_surf)%surf_spml%nbtot,Tdomain,Tdomain%sSurfaces(i_surf))
                  case (DM_FLUID_PML)
                     call compute_normals(Tdomain, Tdomain%sSurfaces(i_surf)%surf_fpml, DM_FLUID_PML,Tdomain%sSurfaces(i_surf)%Surf_BtN)
                     call get_surf_gll_coord(Tdomain%sSurfaces(i_surf)%surf_fpml,Tdomain, Tdomain%sSurfaces(i_surf)%coord)
                     call surfaceSource(Tdomain%sSurfaces(i_surf)%surf_fpml%nbtot,Tdomain,Tdomain%sSurfaces(i_surf))
                end select
              enddo
         endif

        ! Solid-Fluid interfaces : normal vectors
        if(Tdomain%logicD%SF_local_present)then
            call compute_normals(Tdomain, Tdomain%SF%intSolFlu%surf1, DM_FLUID, Tdomain%SF%SF_BtN)
            call compute_normals(Tdomain, Tdomain%SF%intSolFluPml%surf1, DM_FLUID_PML, Tdomain%SF%SFpml_BtN)
!            call dump_sf_btn(Tdomain,"BEFORE  ")
            call exchange_sf_normals(Tdomain)
!            call dump_sf_btn(Tdomain,"AFTER   ")
        endif ! Solid-Fluid interface


        if (Tdomain%sdom%nbelem>0) Tdomain%sdom%m_Jacob    = abs(Tdomain%sdom%m_Jacob   )
        if (Tdomain%fdom%nbelem>0) Tdomain%fdom%m_Jacob    = abs(Tdomain%fdom%m_Jacob   )
        if (Tdomain%spmldom%nbelem>0) Tdomain%spmldom%m_Jacob = abs(Tdomain%spmldom%m_Jacob)
        if (Tdomain%fpmldom%nbelem>0) Tdomain%fpmldom%m_Jacob = abs(Tdomain%fpmldom%m_Jacob)

    end subroutine shape8_init
    !-------------------------------------------------------------------------
    subroutine dump_sf_btn(Tdomain,s)
        use mrenumber
        type(domain), intent(inout) :: Tdomain
        integer, allocatable, dimension(:) :: renum
        character(len=8),intent(in) :: s
        !
        integer :: nf, i, j, idx, idom, iglob, ngll
        real(kind=fpp) :: nx, ny, nz, px, py, pz
        call get_surface_numbering(Tdomain,Tdomain%SF%intSolFlu%surf1, DM_FLUID, renum)

        do nf=0,Tdomain%n_face-1
            ngll = Tdomain%sFace(nf)%ngll
            if (Tdomain%sFace(nf)%domain/=DM_FLUID) cycle
            do j=0,ngll-1
                do i=0,ngll-1
                    idom = Tdomain%sFace(nf)%Idom(i,j)
                    iglob = Tdomain%sFace(nf)%Iglobnum_Face(i,j)
                    idx = renum(idom)
                    if (idx==-1) cycle
                    px = Tdomain%GlobCoord(0,iglob)
                    py = Tdomain%GlobCoord(1,iglob)
                    pz = Tdomain%GlobCoord(2,iglob)
                    nx = Tdomain%SF%SF_BtN(0,idx)
                    ny = Tdomain%SF%SF_BtN(1,idx)
                    nz = Tdomain%SF%SF_BtN(2,idx)
                    write(*,"(A8,I4,A,I4,A,I4,A,F10.3,F10.3,F10.3,A,F10.3)") s,Tdomain%rank, "/", i, ":", idx, &
                        "(",px, py, pz,")->", nz
                end do
            end do
        end do
    end subroutine dump_sf_btn
    !-------------------------------------------------------------------------
    subroutine compute_normals(Tdomain, surf, dom, BtN)
        use mrenumber
        use splib
        implicit none
        type(domain), intent(inout) :: Tdomain
        type(surf_num), intent(inout) :: surf
        integer, intent(in) :: dom
        real(fpp), allocatable, dimension(:,:), intent(out) :: BtN
        !
        integer, allocatable, dimension(:) :: renum
        integer :: nf, nfs
        integer :: ngll
        integer :: i, j, idx
        real(fpp), allocatable, dimension(:) :: gllc
        real(fpp), allocatable, dimension(:) :: gllw
        real(fpp), dimension(0:2, 0:3) :: nodes
        real(fpp), dimension(0:2) :: normal
        real(fpp) :: orient

        if (surf%nbtot==0) return
        allocate(BtN(0:2, 0:surf%nbtot-1))

        Btn(:,:) = 0.
        call get_surface_numbering(Tdomain, surf, dom, renum)
        ngll = domain_ngll(Tdomain, dom)
        call domain_gllc(Tdomain, dom, gllc)
        call domain_gllw(Tdomain, dom, gllw)
        ! SF faces: the normal is taken outward from the fluid elements
        do nf = 0, surf%n_faces-1
            nfs = surf%if_faces(nf)
            orient = -1d0*surf%if_norm(nf)
            ! We don't treat orphan faces, they will be handled by communications
            if (Tdomain%sFace(nfs)%orphan) cycle
            do i = 0,3
                nodes(:,i) = Tdomain%Coord_nodes(:,Tdomain%sFace(nfs)%inodes(i))
            end do
            do j=0,ngll-1
                do i=0,ngll-1
                    idx = Tdomain%sFace(nfs)%Idom(i,j)
                    if (idx==-1) then
                        write(*,*) "ERROR, getting face gll:", dom, Tdomain%sFace(nfs)%domain
                        write(*,*) "IND:", Tdomain%sFace(nfs)%Idom
                        stop 1
                    end if
                    idx = renum(idx)
                    call normal_face(nodes, gllc(i), gllc(j), normal)
                    !
                    if (idx==-1) then
                        write(*,*) "ERROR, getting face gll:", dom, Tdomain%sFace(nfs)%domain
                        write(*,*) "IND:", Tdomain%sFace(nfs)%Idom
                        stop 1
                    end if
                    BtN(:,idx) = BtN(:,idx) + orient*normal*gllw(i)*gllw(j)
                end do
            end do
        end do
        if (allocated(gllc)) deallocate(gllc)
        if (allocated(gllw)) deallocate(gllw)
        deallocate(renum)
    end subroutine compute_normals
    !---------------------------------------------------------------------------
    subroutine normal_face(nodes, xi, eta, normal)
        use shape_geom_3d
        real(fpp), dimension(0:2, 0:3), intent(in) :: nodes
        real(fpp), intent(in) :: xi, eta
        real(fpp), dimension(0:2), intent(out) :: normal
        !
        real(fpp), dimension(0:2) :: d_xi
        real(fpp), dimension(0:2) :: d_eta
        integer :: i
        do i=0,2
            d_xi(i) = 0.25*( &
                - nodes(i,0)*(1-eta) &
                + nodes(i,1)*(1-eta) &
                + nodes(i,2)*(1+eta) &
                - nodes(i,3)*(1+eta))
            d_eta(i) = 0.25*( &
                - nodes(i,0)*(1-xi) &
                - nodes(i,1)*(1+xi) &
                + nodes(i,2)*(1+xi) &
                + nodes(i,3)*(1-xi) )
        end do
        call cross_prod(d_xi, d_eta, normal)
        ! We don't normalize, we want normal = n.dS
    end subroutine normal_face
    !---------------------------------------------------------------------------
    subroutine nodes_coord_8(Control_Nodes,n_glob_nodes,Coord_Nodes,coord)
        ! gives the coordinates of the 8 nodes, for a given element
        implicit none
        integer, dimension(0:7), intent(in)   :: Control_Nodes
        integer, intent(in)  :: n_glob_nodes
        real(fpp), dimension(0:2,0:n_glob_nodes-1), intent(in)  :: coord_nodes
        real(fpp), dimension(0:2,0:7), intent(out)  :: coord
        integer   :: n,i_n

        do n = 0,7
            i_n = Control_Nodes(n)
            coord(0,n) = Coord_Nodes(0,i_n)
            coord(1,n) = Coord_Nodes(1,i_n)
            coord(2,n) = Coord_Nodes(2,i_n)
        end do

    end subroutine nodes_coord_8
    !---------------------------------------------------------------------------
    subroutine shape8_local2global(coord, xi, eta, zeta, xa, ya, za)
        real(fpp), dimension(0:2,0:7), intent(in)  :: coord
        real(fpp), intent(in) :: xi, eta, zeta
        real(fpp), intent(out) :: xa, ya, za
        xa = 0.125 * (coord(0,0)*(1-xi)*(1-eta)*(1-zeta) + coord(0,1)*(1+xi)*(1-eta)*(1-zeta) + &
            coord(0,2)*(1+xi)*(1+eta)*(1-zeta) + coord(0,3)*(1-xi)*(1+eta)*(1-zeta) + &
            coord(0,4)*(1-xi)*(1-eta)*(1+zeta) + coord(0,5)*(1+xi)*(1-eta)*(1+zeta) + &
            coord(0,6)*(1+xi)*(1+eta)*(1+zeta) + coord(0,7)*(1-xi)*(1+eta)*(1+zeta))
        ya = 0.125 * (coord(1,0)*(1-xi)*(1-eta)*(1-zeta) + coord(1,1)*(1+xi)*(1-eta)*(1-zeta) + &
            coord(1,2)*(1+xi)*(1+eta)*(1-zeta) + coord(1,3)*(1-xi)*(1+eta)*(1-zeta) + &
            coord(1,4)*(1-xi)*(1-eta)*(1+zeta) + coord(1,5)*(1+xi)*(1-eta)*(1+zeta) + &
            coord(1,6)*(1+xi)*(1+eta)*(1+zeta) + coord(1,7)*(1-xi)*(1+eta)*(1+zeta))
        za = 0.125 * (coord(2,0)*(1-xi)*(1-eta)*(1-zeta) + coord(2,1)*(1+xi)*(1-eta)*(1-zeta) + &
            coord(2,2)*(1+xi)*(1+eta)*(1-zeta) + coord(2,3)*(1-xi)*(1+eta)*(1-zeta) + &
            coord(2,4)*(1-xi)*(1-eta)*(1+zeta) + coord(2,5)*(1+xi)*(1-eta)*(1+zeta) + &
            coord(2,6)*(1+xi)*(1+eta)*(1+zeta) + coord(2,7)*(1-xi)*(1+eta)*(1+zeta))
    end subroutine shape8_local2global
    !---------------------------------------------------------------------------
    subroutine shape8_local2jacob(coord, xi, eta, zeta, jac)
        real(fpp), dimension(0:2,0:7), intent(in)  :: coord
        real(fpp), intent(in) :: xi, eta, zeta
        real(fpp), dimension(0:2,0:2), intent(out) :: jac
        !- computation of the derivative matrix, dx_(jj)/dxi_(ii)
        ! dx/dxi
        jac(0,0) = 0.125 * ((coord(0,1)-coord(0,0))*(1-eta)*(1-zeta) + &
            (coord(0,2)-coord(0,3))*(1+eta)*(1-zeta) + &
            (coord(0,5)-coord(0,4))*(1-eta)*(1+zeta) + &
            (coord(0,6)-coord(0,7))*(1+eta)*(1+zeta))
        ! dx/deta
        jac(1,0) = 0.125 * ((coord(0,3)-coord(0,0))*(1-xi)*(1-zeta) + &
            (coord(0,2)-coord(0,1))*(1+xi)*(1-zeta) + &
            (coord(0,7)-coord(0,4))*(1-xi)*(1+zeta) + &
            (coord(0,6)-coord(0,5))*(1+xi)*(1+zeta))
        ! dx/dzeta
        jac(2,0) = 0.125 * ((coord(0,4)-coord(0,0))*(1-xi)*(1-eta) + &
            (coord(0,5)-coord(0,1))*(1+xi)*(1-eta) + &
            (coord(0,7)-coord(0,3))*(1-xi)*(1+eta) + &
            (coord(0,6)-coord(0,2))*(1+xi)*(1+eta))
        ! dy/dxi
        jac(0,1) = 0.125 * ((coord(1,1)-coord(1,0))*(1-eta)*(1-zeta) + &
            (coord(1,2)-coord(1,3))*(1+eta)*(1-zeta) + &
            (coord(1,5)-coord(1,4))*(1-eta)*(1+zeta) + &
            (coord(1,6)-coord(1,7))*(1+eta)*(1+zeta))
        ! dy/deta
        jac(1,1) = 0.125 * ((coord(1,3)-coord(1,0))*(1-xi)*(1-zeta) + &
            (coord(1,2)-coord(1,1))*(1+xi)*(1-zeta) + &
            (coord(1,7)-coord(1,4))*(1-xi)*(1+zeta) + &
            (coord(1,6)-coord(1,5))*(1+xi)*(1+zeta))
        ! dy/dzeta
        jac(2,1) = 0.125 * ((coord(1,4)-coord(1,0))*(1-xi)*(1-eta) + &
            (coord(1,5)-coord(1,1))*(1+xi)*(1-eta) + &
            (coord(1,7)-coord(1,3))*(1-xi)*(1+eta) + &
            (coord(1,6)-coord(1,2))*(1+xi)*(1+eta))
        ! dz/dxi
        jac(0,2) = 0.125 * ((coord(2,1)-coord(2,0))*(1-eta)*(1-zeta) + &
            (coord(2,2)-coord(2,3))*(1+eta)*(1-zeta) + &
            (coord(2,5)-coord(2,4))*(1-eta)*(1+zeta) + &
            (coord(2,6)-coord(2,7))*(1+eta)*(1+zeta))
        ! dz/deta
        jac(1,2) = 0.125 * ((coord(2,3)-coord(2,0))*(1-xi)*(1-zeta) + &
            (coord(2,2)-coord(2,1))*(1+xi)*(1-zeta) + &
            (coord(2,7)-coord(2,4))*(1-xi)*(1+zeta) + &
            (coord(2,6)-coord(2,5))*(1+xi)*(1+zeta))
        ! dz/dzeta
        jac(2,2) = 0.125 * ((coord(2,4)-coord(2,0))*(1-xi)*(1-eta) + &
            (coord(2,5)-coord(2,1))*(1+xi)*(1-eta) + &
            (coord(2,7)-coord(2,3))*(1-xi)*(1+eta) + &
            (coord(2,6)-coord(2,2))*(1+xi)*(1+eta))
    end subroutine shape8_local2jacob
    !---------------------------------------------------------------------------
    real(fpp) function shape8_min(dim, nn, x, nodes, xref)
        integer, intent(in) :: dim, nn
        real(fpp), dimension(0:dim-1), intent(in) :: x, xref
        real(fpp), dimension(0:dim-1,0:nn-1), intent(in) :: nodes
        real(fpp) :: xa, ya, za
        call shape8_local2global(nodes, x(0), x(1), x(2), xa, ya, za)
        shape8_min = (xa-xref(0))**2 + (ya-xref(1))**2 + (za-xref(2))**2
    end function shape8_min
    !---------------------------------------------------------------------------
    subroutine shape8_mingrad(dim, nn, x, nodes, xref, grad)
        integer, intent(in) :: dim, nn
        real(fpp), dimension(0:dim-1), intent(in) :: x, xref
        real(fpp), dimension(0:dim-1,0:nn-1), intent(in) :: nodes
        real(fpp), dimension(0:dim-1), intent(out) :: grad
        real(fpp), dimension(0:dim-1,0:dim-1) :: jac
        real(fpp) :: xa, ya, za
        call shape8_local2global(nodes, x(0), x(1), x(2), xa, ya, za)
        call shape8_local2jacob(nodes, x(0), x(1), x(2), jac)
        grad(0) = 2*(jac(0,0)*(xa-xref(0))+jac(0,1)*(ya-xref(1))+jac(0,2)*(za-xref(2)))
        grad(1) = 2*(jac(1,0)*(xa-xref(0))+jac(1,1)*(ya-xref(1))+jac(1,2)*(za-xref(2)))
        grad(2) = 2*(jac(2,0)*(xa-xref(0))+jac(2,1)*(ya-xref(1))+jac(2,2)*(za-xref(2)))
    end subroutine shape8_mingrad
    !---------------------------------------------------------------------------
    subroutine simple_newton_8(nodes, xref, xin, xout, nit)
        real(fpp), dimension(0:2), intent(in) :: xref, xin
        real(fpp), dimension(0:2), intent(out) :: xout
        integer, intent(out) :: nit
        real(fpp), dimension(0:2,0:7), intent(in) :: nodes
        real(fpp), dimension(0:2,0:2) :: jac
        real(fpp), dimension(0:2) :: x
        real(fpp) :: xa, ya, za, err, Det
        integer, parameter :: niter=1000
        integer :: i
        xout = xin
        do i=1,niter
            call shape8_local2global(nodes, xout(0), xout(1), xout(2), xa, ya, za)
            call shape8_local2jacob(nodes, xout(0), xout(1), xout(2), Jac)
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
    end subroutine simple_newton_8
    !---------------------------------------------------------------------------
    subroutine shape8_global2local(coord, xa, ya, za, xi, eta, zeta, ok)
        use mleastsq
        real(fpp), dimension(0:2,0:7), intent(in)  :: coord
        real(fpp), intent(in) :: xa, ya, za
        real(fpp), intent(out) :: xi, eta, zeta
        logical, intent(out) :: ok
        !
        integer :: niter
        real(fpp), dimension(0:2) :: xin, xout, xref
        ok = .true.
        xin(0) = 0.
        xin(1) = 0.
        xin(2) = 0.
        xref(0) = xa
        xref(1) = ya
        xref(2) = za
        call simple_newton_8(coord, xref, xin, xout, niter)
        if (niter==1000 .or. niter<0) ok=.false.
        xi = xout(0)
        eta = xout(1)
        zeta = xout(2)
    end subroutine shape8_global2local
end module mshape8

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
