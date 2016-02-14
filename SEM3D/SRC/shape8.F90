!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file shape8.f90
!!\brief Contient la subroutine shape8.
!!
!<
module mshape8
    use sdomain
    implicit none
#include "index.h"

contains
    subroutine shape8_init(Tdomain)
        use shape_geom_3d
        use mpi
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer :: n,  ngllx,nglly,ngllz, mat, i,j,k ,ipoint
!            nf,ne,nv
!            mat_index,which_elem,which_face,which_edge,which_vertex,
!        integer, dimension(0:3)  :: loc_vertices,loc_nodes
        real(FPP) :: Jac, xi, eta, zeta
        real(FPP), dimension(0:2,0:2) :: LocInvGrad
        real(FPP), dimension(0:2,0:7) :: coord  ! coordinates of nodes
        real(FPP) :: xp, yp, zp
        real(FPP), parameter :: XEPS=1e-10
        real(FPP) :: dxp, dyp, dzp

        ! Tdomain%n_glob_points is the number of degrees of fredom
        allocate(Tdomain%GlobCoord(0:2,0:Tdomain%n_glob_points-1))
        Tdomain%GlobCoord = 0d0

        do n = 0,Tdomain%n_elem-1

            ! coordinates of control nodes (which are vertices also)
            call nodes_coord_8(Tdomain%specel(n)%Control_Nodes(0:),Tdomain%n_glob_nodes,    &
                Tdomain%Coord_Nodes(0:,0:),coord)

            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            mat = Tdomain%specel(n)%mat_index

            ! coordinates of GLL points, and values of Jacobian and dX_dxi at each GLL point.
            do k = 0,ngllz-1
                zeta = Tdomain%sSubdomain(mat)%GLLc(k)
                do j = 0,nglly-1
                    eta = Tdomain%sSubdomain(mat)%GLLc(j)
                    do i = 0,ngllx-1
                        xi = Tdomain%sSubdomain(mat)%GLLc(i)

                        ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)

                        call shape8_local2global(coord, xi, eta, zeta, xp, yp, zp)
                        ! DEBUG CHECK
                        if (Tdomain%GlobCoord(0,ipoint)/=0d0 .or. &
                            Tdomain%GlobCoord(1,ipoint)/=0d0 .or. &
                            Tdomain%GlobCoord(2,ipoint)/=0d0) then
                            dxp = abs(xp-Tdomain%GlobCoord(0,ipoint))
                            dyp = abs(yp-Tdomain%GlobCoord(1,ipoint))
                            dzp = abs(zp-Tdomain%GlobCoord(2,ipoint))
                            if (dxp>XEPS .or. dyp>XEPS .or. dzp>XEPS) then
                                write(*,*) "DIFF", n, ipoint, Tdomain%GlobCoord(:,ipoint), ":", dxp, dyp, dzp
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
                                Tdomain%sdom%Jacob_     (        i,j,k,Tdomain%specel(n)%lnum) = Jac
                                Tdomain%sdom%InvGrad_   (0:2,0:2,i,j,k,Tdomain%specel(n)%lnum) = LocInvGrad(0:2,0:2)
                            case (DM_FLUID)
                                Tdomain%fdom%Jacob_     (        i,j,k,Tdomain%specel(n)%lnum) = Jac
                                Tdomain%fdom%InvGrad_   (0:2,0:2,i,j,k,Tdomain%specel(n)%lnum) = LocInvGrad(0:2,0:2)
                            case (DM_SOLID_PML)
                                Tdomain%spmldom%Jacob_  (        i,j,k,Tdomain%specel(n)%lnum) = Jac
                                Tdomain%spmldom%InvGrad_(0:2,0:2,i,j,k,Tdomain%specel(n)%lnum) = LocInvGrad(0:2,0:2)
                            case (DM_FLUID_PML)
                                Tdomain%fpmldom%Jacob_  (        i,j,k,Tdomain%specel(n)%lnum) = Jac
                                Tdomain%fpmldom%InvGrad_(0:2,0:2,i,j,k,Tdomain%specel(n)%lnum) = LocInvGrad(0:2,0:2)
                        end select
                    enddo
                enddo
            enddo  ! end of loops onto GLL points inside an element
        enddo    ! end of loop onto elements


        ! Neumann Boundary Conditions : normal vectors
        if(Tdomain%logicD%neumann_local_present)then
            ! call compute_normals(Tdomain, Tdomain%neu%surf, Tdomain%neu%btn)
!            ! Neumann faces
!            do nf = 0, Tdomain%Neumann%Neu_n_faces-1
!                which_face = Tdomain%Neumann%Neu_Face(nf)%Face
!                which_elem = Tdomain%sFace(which_face)%which_elem
!                ngll1 = Tdomain%Neumann%Neu_Face(nf)%ngll1
!                ngll2 = Tdomain%Neumann%Neu_Face(nf)%ngll2
!                ngllx = Tdomain%specel(which_elem)%ngllx
!                nglly = Tdomain%specel(which_elem)%nglly
!                ngllz = Tdomain%specel(which_elem)%ngllz
!                f_dir = Tdomain%Neumann%Neu_Face(nf)%dir
!                mat_index = Tdomain%specel(which_elem)%mat_index
!
!                allocate(Tdomain%Neumann%Neu_Face(nf)%normal(0:ngll1-1,0:ngll2-1,0:2))
!                do i = 0,3
!                    loc_vertices(i) =        &
!                        Tdomain%Neumann%Neu_Vertex(Tdomain%Neumann%Neu_Face(nf)%Near_Vertices(i))%vertex
!                    loc_nodes(i) = Tdomain%sVertex(loc_vertices(i))%Iglobnum_Vertex
!                end do
!
!                call coord_nodes_face(f_dir,coord,loc_nodes,Tdomain%n_glob_nodes,    &
!                    Tdomain%Coord_Nodes(0:,0:))
!                call normal_face(f_dir,ngllx,nglly,ngllz,ngll1,ngll2,coord,                     &
!                    Tdomain%sSubdomain(mat_index)%GLLcx,Tdomain%sSubdomain(mat_index)%GLLcy,   &
!                    Tdomain%sSubdomain(mat_index)%GLLcz,Tdomain%Neumann%Neu_Face(nf)%normal)
!
!                ! eventual switch of the normal direction
!                call inversion_normal(f_dir,Tdomain%specel(which_elem),ngll1,ngll2,    &
!                    Tdomain%Neumann%Neu_Face(nf)%normal)
!
!            end do
!
!            ! co-ordinates of Neumann GLL points: necessary to impose a boundary condition
!            do nf = 0,Tdomain%Neumann%Neu_n_faces-1
!                ngll1 = Tdomain%Neumann%Neu_Face(nf)%ngll1
!                ngll2 = Tdomain%Neumann%Neu_Face(nf)%ngll2
!                which_face = Tdomain%Neumann%Neu_Face(nf)%Face
!                allocate(Tdomain%Neumann%Neu_Face(nf)%Coord(1:ngll1-2,1:ngll2-2,0:2))
!                do j = 1,ngll2-2
!                    do i = 1,ngll1-2
!                        Tdomain%Neumann%Neu_Face(nf)%Coord(i,j,0:2) =    &
!                            Tdomain%GlobCoord(0:2,Tdomain%sFace(which_face)%Iglobnum_Face(i,j))
!                    enddo
!                enddo
!            end do
!            do ne = 0,Tdomain%Neumann%Neu_n_edges-1
!                ngll = Tdomain%Neumann%Neu_Edge(ne)%ngll
!                which_edge = Tdomain%Neumann%Neu_Edge(ne)%Edge
!                allocate(Tdomain%Neumann%Neu_Edge(ne)%Coord(1:ngll-2,0:2))
!                do i = 1,ngll-2
!                    Tdomain%Neumann%Neu_Edge(ne)%Coord(i,0:2) =    &
!                        Tdomain%GlobCoord(0:2,Tdomain%sEdge(which_edge)%Iglobnum_Edge(i))
!                enddo
!            end do
!            do nv = 0, Tdomain%Neumann%Neu_n_vertices-1
!                which_vertex = Tdomain%Neumann%Neu_Vertex(nv)%Vertex
!                Tdomain%Neumann%Neu_Vertex(nv)%Coord(0:2) =     &
!                    Tdomain%GlobCoord(0:2,Tdomain%sVertex(which_vertex)%Iglobnum_Vertex)
!            enddo
!
        endif ! Neumann

        ! Solid-Fluid interfaces : normal vectors
        if(Tdomain%logicD%SF_local_present)then
            call compute_normals(Tdomain, Tdomain%SF%intSolFlu%surf1, DM_FLUID, Tdomain%SF%SF_BtN)
            call compute_normals(Tdomain, Tdomain%SF%intSolFluPml%surf1, DM_FLUID_PML, Tdomain%SF%SFpml_BtN)
        endif ! Solid-Fluid interface


        if (Tdomain%sdom%nbelem>0) Tdomain%sdom%m_Jacob    = abs(Tdomain%sdom%m_Jacob   )
        if (Tdomain%fdom%nbelem>0) Tdomain%fdom%m_Jacob    = abs(Tdomain%fdom%m_Jacob   )
        if (Tdomain%spmldom%nbelem>0) Tdomain%spmldom%m_Jacob = abs(Tdomain%spmldom%m_Jacob)
        if (Tdomain%fpmldom%nbelem>0) Tdomain%fpmldom%m_Jacob = abs(Tdomain%fpmldom%m_Jacob)

    end subroutine shape8_init
    !-------------------------------------------------------------------------
    subroutine compute_normals(Tdomain, surf, dom, BtN)
        use mrenumber
        use splib
        implicit none
        type(domain), intent(inout) :: Tdomain
        type(surf_num), intent(inout) :: surf
        integer, intent(in) :: dom
        real(kind=FPP), allocatable, dimension(:,:), intent(out) :: BtN
        !
        integer, allocatable, dimension(:) :: renum
        integer :: nf, nfs
        integer :: ngll1, ngll2
        integer :: i, j, idx
        real(kind=FPP), allocatable, dimension(:) :: gllc1, gllc2
        real(kind=FPP), allocatable, dimension(:) :: pol1, pol2
        real(kind=FPP), allocatable, dimension(:) :: gllw1, gllw2
        real(kind=FPP), dimension(0:2, 0:3) :: nodes
        real(kind=FPP), dimension(0:2) :: normal
        real(kind=FPP) :: orient
        allocate(BtN(0:2, 0:surf%nbtot-1))
        ngll1 = 0
        ngll2 = 0

        Btn(:,:) = 0.
        call get_surface_numbering(Tdomain, surf, dom, renum)
        ! SF faces: the normal is taken outward from the fluid elements
        do nf = 0, surf%n_faces-1
            nfs = surf%if_faces(nf)
            orient = -1d0*surf%if_norm(nf)
            if (ngll1/=Tdomain%sFace(nfs)%ngll1) then
                if (allocated(gllc1)) deallocate(gllc1, pol1, gllw1)
                ngll1 = Tdomain%sFace(nfs)%ngll1
                allocate(pol1(0:ngll1-1))
                allocate(gllc1(0:ngll1-1))
                allocate(gllw1(0:ngll1-1))
                call zelegl(ngll1-1, gllc1, pol1)
                call welegl(ngll1-1, gllc1, pol1, gllw1)
            end if
            if (ngll2/=Tdomain%sFace(nfs)%ngll2) then
                if (allocated(gllc2)) deallocate(gllc2, pol2, gllw2)
                ngll2 = Tdomain%sFace(nfs)%ngll2
                allocate(pol2(0:ngll2-1))
                allocate(gllc2(0:ngll2-1))
                allocate(gllw2(0:ngll2-1))
                call zelegl(ngll2-1, gllc2, pol2)
                call welegl(ngll2-1, gllc2, pol2, gllw2)
            end if
            do i = 0,3
                nodes(:,i) = Tdomain%Coord_nodes(:,Tdomain%sFace(nfs)%inodes(i))
            end do
            do j=0,ngll2-1
                do i=0,ngll1-1
                    idx = Tdomain%sFace(nfs)%Idom(i,j)
                    if (idx==-1) then
                        write(*,*) "ERROR, getting face gll:", dom, Tdomain%sFace(nfs)%domain
                        write(*,*) "IND:", Tdomain%sFace(nfs)%Idom
                        stop 1
                    end if
                    idx = renum(idx)
                    call normal_face(nodes, gllc1(i), gllc2(j), normal)
                    !
                    if (idx==-1) then
                        write(*,*) "ERROR, getting face gll:", dom, Tdomain%sFace(nfs)%domain
                        write(*,*) "IND:", Tdomain%sFace(nfs)%Idom
                        stop 1
                    end if
                    BtN(:,idx) = BtN(:,idx) + orient*normal*gllw1(i)*gllw2(j)
                end do
            end do
        end do
        if (allocated(gllc1)) deallocate(gllc1)
        if (allocated(pol1)) deallocate(pol1)
        if (allocated(gllc2)) deallocate(gllc2)
        if (allocated(pol2)) deallocate(pol2)
        deallocate(renum)
    end subroutine compute_normals
    !---------------------------------------------------------------------------
    subroutine normal_face(nodes, xi, eta, normal)
        use shape_geom_3d
        real(kind=FPP), dimension(0:2, 0:3), intent(in) :: nodes
        real(kind=FPP), intent(in) :: xi, eta
        real(kind=FPP), dimension(0:2), intent(out) :: normal
        !
        real(kind=FPP), dimension(0:2) :: d_xi
        real(kind=FPP), dimension(0:2) :: d_eta
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
        real, dimension(0:2,0:n_glob_nodes-1), intent(in)  :: coord_nodes
        real, dimension(0:2,0:7), intent(out)  :: coord
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
        real, dimension(0:2,0:7), intent(in)  :: coord
        real, intent(in) :: xi, eta, zeta
        real, intent(out) :: xa, ya, za
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
        double precision, dimension(0:2,0:7), intent(in)  :: coord
        double precision, intent(in) :: xi, eta, zeta
        double precision, dimension(0:2,0:2), intent(out) :: jac
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
    double precision function shape8_min(dim, nn, x, nodes, xref)
        integer, intent(in) :: dim, nn
        double precision, dimension(0:dim-1), intent(in) :: x, xref
        double precision, dimension(0:dim-1,0:nn-1), intent(in) :: nodes
        double precision :: xa, ya, za
        call shape8_local2global(nodes, x(0), x(1), x(2), xa, ya, za)
        shape8_min = (xa-xref(0))**2 + (ya-xref(1))**2 + (za-xref(2))**2
    end function shape8_min
    !---------------------------------------------------------------------------
    subroutine shape8_mingrad(dim, nn, x, nodes, xref, grad)
        integer, intent(in) :: dim, nn
        double precision, dimension(0:dim-1), intent(in) :: x, xref
        double precision, dimension(0:dim-1,0:nn-1), intent(in) :: nodes
        double precision, dimension(0:dim-1), intent(out) :: grad
        double precision, dimension(0:dim-1,0:dim-1) :: jac
        double precision :: xa, ya, za
        call shape8_local2global(nodes, x(0), x(1), x(2), xa, ya, za)
        call shape8_local2jacob(nodes, x(0), x(1), x(2), jac)
        grad(0) = 2*(jac(0,0)*(xa-xref(0))+jac(0,1)*(ya-xref(1))+jac(0,2)*(za-xref(2)))
        grad(1) = 2*(jac(1,0)*(xa-xref(0))+jac(1,1)*(ya-xref(1))+jac(1,2)*(za-xref(2)))
        grad(2) = 2*(jac(2,0)*(xa-xref(0))+jac(2,1)*(ya-xref(1))+jac(2,2)*(za-xref(2)))
    end subroutine shape8_mingrad
    !---------------------------------------------------------------------------
    subroutine simple_newton_8(nodes, xref, xin, xout, nit)
        double precision, dimension(0:2), intent(in) :: xref, xin
        double precision, dimension(0:2), intent(out) :: xout
        integer, intent(out) :: nit
        double precision, dimension(0:2,0:7), intent(in) :: nodes
        double precision, dimension(0:2,0:2) :: jac
        double precision, dimension(0:2) :: x
        double precision :: xa, ya, za, err, Det
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
        double precision, dimension(0:2,0:7), intent(in)  :: coord
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
        call minimize_cg(3, 8, xin, coord, xref, shape8_min, shape8_mingrad, 0.001D0, xout, niter)
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
