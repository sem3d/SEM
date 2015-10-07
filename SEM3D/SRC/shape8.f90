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
    implicit none
contains
    subroutine shape8_init(Tdomain)
        use shape_geom_3d
        use sdomain
        use mpi
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer :: n,  ngll,ngllx,nglly,ngllz,ngll1,ngll2, mat, i,j,k ,ipoint,   &
            nf,ne,nv,               &
            mat_index,which_elem,which_face,which_edge,which_vertex,f_dir,indv
        integer, dimension(0:3)  :: loc_vertices,loc_nodes
        real :: xi,eta,zeta, Jac, fact_norm
        real, dimension(0:2,0:2) :: LocInvGrad
        real, dimension(0:2,0:7) :: coord  ! coordinates of nodes
        real :: xp, yp, zp

        ! Tdomain%n_glob_points is the number of degrees of fredom
        allocate(Tdomain%GlobCoord(0:2,0:Tdomain%n_glob_points-1))

        do n = 0,Tdomain%n_elem-1

            ! coordinates of control nodes (which are vertices also)
            call nodes_coord_8(Tdomain%specel(n)%Control_Nodes(0:),Tdomain%n_glob_nodes,    &
                Tdomain%Coord_Nodes(0:,0:),coord)

            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            mat = Tdomain%specel(n)%mat_index

            allocate(Tdomain%specel(n)%Jacob(0:ngllx-1,0:nglly-1,0:ngllz-1))
            allocate(Tdomain%specel(n)%InvGrad(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2))

            ! coordinates of GLL points, and values of Jacobian and dX_dxi at each GLL point.
            do k = 0,ngllz-1
                zeta = Tdomain%sSubdomain(mat)%GLLcz(k)
                do j = 0,nglly-1
                    eta = Tdomain%sSubdomain(mat)%GLLcy(j)
                    do i = 0,ngllx-1
                        xi = Tdomain%sSubdomain(mat)%GLLcx(i)

                        ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)

                        call shape8_local2global(coord, xi, eta, zeta, xp, yp, zp)
                        Tdomain%GlobCoord(0,ipoint) = xp
                        Tdomain%GlobCoord(1,ipoint) = yp
                        Tdomain%GlobCoord(2,ipoint) = zp

                        call shape8_local2jacob(coord, xi, eta, zeta, LocInvGrad)

                        call invert_3d(LocInvGrad,Jac)

                        Tdomain%specel(n)%Jacob(i,j,k) = Jac
                        Tdomain%specel(n)%InvGrad(i,j,k,0:2,0:2) = LocInvGrad(0:2,0:2)

                    enddo
                enddo
            enddo  ! end of loops onto GLL points inside an element
        enddo    ! end of loop onto elements


        ! Neumann Boundary Conditions : normal vectors
        if(Tdomain%logicD%neumann_local_present)then
            ! Neumann faces
            do nf = 0, Tdomain%Neumann%Neu_n_faces-1
                which_face = Tdomain%Neumann%Neu_Face(nf)%Face
                which_elem = Tdomain%sFace(which_face)%which_elem
                ngll1 = Tdomain%Neumann%Neu_Face(nf)%ngll1
                ngll2 = Tdomain%Neumann%Neu_Face(nf)%ngll2
                ngllx = Tdomain%specel(which_elem)%ngllx
                nglly = Tdomain%specel(which_elem)%nglly
                ngllz = Tdomain%specel(which_elem)%ngllz
                f_dir = Tdomain%Neumann%Neu_Face(nf)%dir
                mat_index = Tdomain%specel(which_elem)%mat_index

                allocate(Tdomain%Neumann%Neu_Face(nf)%normal(0:ngll1-1,0:ngll2-1,0:2))
                do i = 0,3
                    loc_vertices(i) =        &
                        Tdomain%Neumann%Neu_Vertex(Tdomain%Neumann%Neu_Face(nf)%Near_Vertices(i))%vertex
                    loc_nodes(i) = Tdomain%sVertex(loc_vertices(i))%global_numbering
                end do

                call coord_nodes_face(f_dir,coord,loc_nodes,Tdomain%n_glob_nodes,    &
                    Tdomain%Coord_Nodes(0:,0:))
                call normal_face(f_dir,ngllx,nglly,ngllz,ngll1,ngll2,coord,                     &
                    Tdomain%sSubdomain(mat_index)%GLLcx,Tdomain%sSubdomain(mat_index)%GLLcy,   &
                    Tdomain%sSubdomain(mat_index)%GLLcz,Tdomain%Neumann%Neu_Face(nf)%normal)

                ! eventual switch of the normal direction
                call inversion_normal(f_dir,Tdomain%specel(which_elem),ngll1,ngll2,    &
                    Tdomain%Neumann%Neu_Face(nf)%normal)

            end do

            ! co-ordinates of Neumann GLL points: necessary to impose a boundary condition
            do nf = 0,Tdomain%Neumann%Neu_n_faces-1
                ngll1 = Tdomain%Neumann%Neu_Face(nf)%ngll1
                ngll2 = Tdomain%Neumann%Neu_Face(nf)%ngll2
                which_face = Tdomain%Neumann%Neu_Face(nf)%Face
                allocate(Tdomain%Neumann%Neu_Face(nf)%Coord(1:ngll1-2,1:ngll2-2,0:2))
                do j = 1,ngll2-2
                    do i = 1,ngll1-2
                        Tdomain%Neumann%Neu_Face(nf)%Coord(i,j,0:2) =    &
                            Tdomain%GlobCoord(0:2,Tdomain%sFace(which_face)%Iglobnum_Face(i,j))
                    enddo
                enddo
            end do
            do ne = 0,Tdomain%Neumann%Neu_n_edges-1
                ngll = Tdomain%Neumann%Neu_Edge(ne)%ngll
                which_edge = Tdomain%Neumann%Neu_Edge(ne)%Edge
                allocate(Tdomain%Neumann%Neu_Edge(ne)%Coord(1:ngll-2,0:2))
                do i = 1,ngll-2
                    Tdomain%Neumann%Neu_Edge(ne)%Coord(i,0:2) =    &
                        Tdomain%GlobCoord(0:2,Tdomain%sEdge(which_edge)%Iglobnum_Edge(i))
                enddo
            end do
            do nv = 0, Tdomain%Neumann%Neu_n_vertices-1
                which_vertex = Tdomain%Neumann%Neu_Vertex(nv)%Vertex
                Tdomain%Neumann%Neu_Vertex(nv)%Coord(0:2) =     &
                    Tdomain%GlobCoord(0:2,Tdomain%sVertex(which_vertex)%Iglobnum_Vertex)
            enddo

        endif ! Neumann

        ! Solid-Fluid interfaces : normal vectors
        if(Tdomain%logicD%SF_local_present)then
            ! SF faces: the normal is taken outward from the fluid element !!
            do nf = 0, Tdomain%SF%SF_n_faces-1
                which_face = Tdomain%SF%SF_Face(nf)%Face(0)
                fact_norm = 1d0 ; indv = 0
                if(which_face < 0)then  ! this SF face has no fluid side on this proc
                    which_face = Tdomain%SF%SF_Face(nf)%Face(1)
                    fact_norm = -1d0 ; indv = 1
                end if
                which_elem = Tdomain%sFace(which_face)%which_elem
                ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
                ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
                ngllx = Tdomain%specel(which_elem)%ngllx
                nglly = Tdomain%specel(which_elem)%nglly
                ngllz = Tdomain%specel(which_elem)%ngllz
                f_dir = Tdomain%SF%SF_Face(nf)%dir
                mat_index = Tdomain%specel(which_elem)%mat_index

                allocate(Tdomain%SF%SF_Face(nf)%normal(0:ngll1-1,0:ngll2-1,0:2))
                do i = 0,3
                    loc_vertices(i) =        &
                        Tdomain%SF%SF_Vertex(Tdomain%SF%SF_Face(nf)%Near_Vertices(i))%Vertex(indv)
                    loc_nodes(i) = Tdomain%sVertex(loc_vertices(i))%global_numbering
                end do

                call coord_nodes_face(f_dir,coord,loc_nodes,Tdomain%n_glob_nodes,    &
                    Tdomain%Coord_Nodes(0:,0:))
                call normal_face(f_dir,ngllx,nglly,ngllz,ngll1,ngll2,coord, &
                    Tdomain%sSubdomain(mat_index)%GLLcx,Tdomain%sSubdomain(mat_index)%GLLcy,   &
                    Tdomain%sSubdomain(mat_index)%GLLcz,Tdomain%SF%SF_Face(nf)%normal)

                ! eventual switch of the normal direction
                call inversion_normal(f_dir,Tdomain%specel(which_elem),ngll1,ngll2,    &
                    Tdomain%SF%SF_Face(nf)%normal)

                ! correct direction for the normal (outwards from fluid)
                Tdomain%SF%SF_Face(nf)%normal(:,:,:) = fact_norm*Tdomain%SF%SF_Face(nf)%normal(:,:,:)
            end do


        endif ! Solid-Fluid interface


        ! Obtention of a positive Jacobian.
        do n = 0,Tdomain%n_elem - 1
            ngllx = Tdomain%specel(n)%ngllx;
            nglly = Tdomain%specel(n)%nglly;
            ngllz = Tdomain%specel(n)%ngllz;
            do k = 0,ngllz - 1
                do j = 0,nglly - 1
                    do i = 0,ngllx - 1
                        Tdomain%specel(n)%Jacob(i,j,k) = abs(Tdomain%specel(n)%Jacob(i,j,k))
                    enddo
                enddo
            enddo
            !OBS: could be rewriten
            !Tdomain%specel(n)%Jacob(:,:,:) = abs(Tdomain%specel(n)%Jacob(:,:,:))
            !Keeping only the loop over the elements (n)
        enddo


        return
    end subroutine shape8_init
    !-------------------------------------------------------------------------
    subroutine coord_nodes_face(dir,coord,local_nod,n_glob_nodes,coord_nodes)
        ! returns the coordinates for all vertices on a face of interest.
        use svertices
        implicit none
        integer, intent(in)   :: dir
        real, dimension(0:2,0:7), intent(out)  :: coord
        integer, dimension(0:3), intent(in)  :: local_nod
        integer, intent(in)  :: n_glob_nodes
        real, dimension(0:2,0:n_glob_nodes-1), intent(in)  :: coord_nodes

        coord(:,:) = 0d0
        select case(dir)
        case(0)
            coord(:,0) = coord_nodes(:,local_nod(0))
            coord(:,1) = coord_nodes(:,local_nod(1))
            coord(:,2) = coord_nodes(:,local_nod(2))
            coord(:,3) = coord_nodes(:,local_nod(3))
        case(1)
            coord(:,0) = coord_nodes(:,local_nod(0))
            coord(:,1) = coord_nodes(:,local_nod(1))
            coord(:,5) = coord_nodes(:,local_nod(2))
            coord(:,4) = coord_nodes(:,local_nod(3))
        case(2)
            coord(:,1) = coord_nodes(:,local_nod(0))
            coord(:,2) = coord_nodes(:,local_nod(1))
            coord(:,6) = coord_nodes(:,local_nod(2))
            coord(:,5) = coord_nodes(:,local_nod(3))
        case(3)
            coord(:,3) = coord_nodes(:,local_nod(0))
            coord(:,2) = coord_nodes(:,local_nod(1))
            coord(:,6) = coord_nodes(:,local_nod(2))
            coord(:,7) = coord_nodes(:,local_nod(3))
        case(4)
            coord(:,0) = coord_nodes(:,local_nod(0))
            coord(:,3) = coord_nodes(:,local_nod(1))
            coord(:,7) = coord_nodes(:,local_nod(2))
            coord(:,4) = coord_nodes(:,local_nod(3))
        case(5)
            coord(:,4) = coord_nodes(:,local_nod(0))
            coord(:,5) = coord_nodes(:,local_nod(1))
            coord(:,6) = coord_nodes(:,local_nod(2))
            coord(:,7) = coord_nodes(:,local_nod(3))
        case default
            stop 1
        end select
    end subroutine coord_nodes_face
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    subroutine inversion_normal(dir,elem,ngll1,ngll2,normal)
        use selement
        implicit none
        integer, intent(in)   :: dir,ngll1,ngll2
        type(element), intent(in)  :: elem
        real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(inout) :: normal
        integer  :: i,ngllx,nglly,ngllz

        ngllx = elem%ngllx ; nglly = elem%nglly ; ngllz = elem%ngllz

        select case(dir)
        case(0)
            do i = 0,2
                where(elem%Jacob(:,:,0) > 0)
                    normal(:,:,i) = -normal(:,:,i)
                end where
            end do
        case(1)
            do i = 0,2
                where(elem%Jacob(:,0,:) < 0)
                    normal(:,:,i) = -normal(:,:,i)
                end where
            end do
        case(2)
            do i = 0,2
                where(elem%Jacob(ngllx-1,:,:) < 0)
                    normal(:,:,i) = -normal(:,:,i)
                end where
            end do
        case(3)
            do i = 0,2
                where(elem%Jacob(:,nglly-1,:) > 0)
                    normal(:,:,i) = -normal(:,:,i)
                end where
            end do
        case(4)
            do i = 0,2
                where(elem%Jacob(0,:,:) > 0)
                    normal(:,:,i) = -normal(:,:,i)
                end where
            end do
        case(5)
            do i = 0,2
                where(elem%Jacob(:,:,ngllz-1) < 0)
                    normal(:,:,i) = -normal(:,:,i)
                end where
            end do
        end select

    end subroutine inversion_normal
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
