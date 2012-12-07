subroutine shape8(Tdomain,rank)

    use sdomain
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: rank
    integer :: n, i_aus, ngll,ngllx,nglly,ngllz,ngll1,ngll2, mat, i,j,k ,ipoint,   &
        n_face, nf,ne,nv,n_aus,nf_up,nf_down,Orient,n_up,n_down,               &
        mat_index,which_elem,which_face,which_edge,which_vertex,f_dir,code,indv
    integer, dimension(0:3)  :: loc_vertices,loc_nodes
    real :: x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7
    real :: xi,eta,zeta, xp,yp,zp, Jac,norm, fact_norm
    real, dimension(0:2,0:1) :: LocInvGradS
    real, dimension(0:2,0:2) :: LocInvGrad
    real, dimension(0:1,0:1) :: LocInvGrad2
    real, dimension(0:7) :: xco,yco,zco  ! coordinates of nodes
    real, external  :: f_p,f_xp,f_yp,f_zp,der_dx_dxi,der_dx_deta,der_dx_dzeta,     &
        der_dy_dxi,der_dy_deta,der_dy_dzeta,der_dz_dxi,der_dz_deta,der_dz_dzeta
    character(len=100) :: fnamef



    allocate(Tdomain%GlobCoord(0:2,0:Tdomain%n_glob_points-1))

    do n = 0,Tdomain%n_elem-1

        ! coordinates of control nodes (which are vertices also)
        call nodes_coord_8(Tdomain%specel(n)%Control_Nodes(0:),Tdomain%n_glob_nodes,    &
            Tdomain%Coord_Nodes(0:,0:),xco,yco,zco)

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

                    Tdomain%GlobCoord(0,ipoint) = f_p(xco,xi,eta,zeta)
                    Tdomain%GlobCoord(1,ipoint) = f_p(yco,xi,eta,zeta)
                    Tdomain%GlobCoord(2,ipoint) = f_p(zco,xi,eta,zeta)

                    !- computation of the derivative matrix, dx_(jj)/dxi_(ii)
                    LocInvGrad(0,0) = der_dx_dxi(xco,eta,zeta)
                    LocInvGrad(1,0) = der_dx_deta(xco,xi,zeta)
                    LocInvGrad(2,0) = der_dx_dzeta(xco,xi,eta)
                    LocInvGrad(0,1) = der_dy_dxi(yco,eta,zeta)
                    LocInvGrad(1,1) = der_dy_deta(yco,xi,zeta)
                    LocInvGrad(2,1) = der_dy_dzeta(yco,xi,eta)
                    LocInvGrad(0,2) = der_dz_dxi(zco,eta,zeta)
                    LocInvGrad(1,2) = der_dz_deta(zco,xi,zeta)
                    LocInvGrad(2,2) = der_dz_dzeta(zco,xi,eta)

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

            call coord_nodes_face(f_dir,xco,yco,zco,loc_nodes,Tdomain%n_glob_nodes,    &
                Tdomain%Coord_Nodes(0:,0:))
            call normal_face(f_dir,ngllx,nglly,ngllz,ngll1,ngll2,xco,yco,zco,                     &
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

            call coord_nodes_face(f_dir,xco,yco,zco,loc_nodes,Tdomain%n_glob_nodes,    &
                Tdomain%Coord_Nodes(0:,0:))
            call normal_face(f_dir,ngllx,nglly,ngllz,ngll1,ngll2,xco,yco,zco,                     &
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
        ngllx = Tdomain%specel(n)%ngllx ; nglly = Tdomain%specel(n)%nglly ; ngllz = Tdomain%specel(n)%ngllz
        do k = 0,ngllz - 1
            do j = 0,nglly - 1
                do i = 0,ngllx - 1
                    Tdomain%specel(n)%Jacob(i,j,k) = abs(Tdomain%specel(n)%Jacob(i,j,k))
                enddo
            enddo
        enddo
    enddo


    return
end subroutine shape8
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine nodes_coord_8(Control_Nodes,n_glob_nodes,Coord_Nodes,xco,yco,zco)
    ! gives the coordinates of the 8 nodes, for a given element
    implicit none
    integer, dimension(0:7), intent(in)   :: Control_Nodes
    integer, intent(in)  :: n_glob_nodes
    real, dimension(0:2,0:n_glob_nodes-1), intent(in)  :: coord_nodes
    real, dimension(0:7), intent(out)  :: xco,yco,zco
    integer   :: n,i_n

    do n = 0,7
        i_n = Control_Nodes(n)
        xco(n) = Coord_Nodes(0,i_n)
        yco(n) = Coord_Nodes(1,i_n)
        zco(n) = Coord_Nodes(2,i_n)
    end do

end subroutine nodes_coord_8
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
subroutine coord_nodes_face(dir,x,y,z,local_nod,n_glob_nodes,coord_nodes)
    ! returns the coordinates for all vertices on a face of interest.
    use svertices
    implicit none
    integer, intent(in)   :: dir
    real, dimension(0:7), intent(out)  :: x,y,z
    integer, dimension(0:3), intent(in)  :: local_nod
    integer, intent(in)  :: n_glob_nodes
    real, dimension(0:2,0:n_glob_nodes-1), intent(in)  :: coord_nodes
    integer   :: i

    x(:) = 0d0 ; y(:) = 0d0 ; z(:) = 0d0

    select case(dir)
    case(0)
        x(0) = coord_nodes(0,local_nod(0)) ; y(0) = coord_nodes(1,local_nod(0)) ; z(0) = coord_nodes(2,local_nod(0))
        x(1) = coord_nodes(0,local_nod(1)) ; y(1) = coord_nodes(1,local_nod(1)) ; z(1) = coord_nodes(2,local_nod(1))
        x(2) = coord_nodes(0,local_nod(2)) ; y(2) = coord_nodes(1,local_nod(2)) ; z(2) = coord_nodes(2,local_nod(2))
        x(3) = coord_nodes(0,local_nod(3)) ; y(3) = coord_nodes(1,local_nod(3)) ; z(3) = coord_nodes(2,local_nod(3))
    case(1)
        x(0) = coord_nodes(0,local_nod(0)) ; y(0) = coord_nodes(1,local_nod(0)) ; z(0) = coord_nodes(2,local_nod(0))
        x(1) = coord_nodes(0,local_nod(1)) ; y(1) = coord_nodes(1,local_nod(1)) ; z(1) = coord_nodes(2,local_nod(1))
        x(5) = coord_nodes(0,local_nod(2)) ; y(5) = coord_nodes(1,local_nod(2)) ; z(5) = coord_nodes(2,local_nod(2))
        x(4) = coord_nodes(0,local_nod(3)) ; y(4) = coord_nodes(1,local_nod(3)) ; z(4) = coord_nodes(2,local_nod(3))
    case(2)
        x(1) = coord_nodes(0,local_nod(0)) ; y(1) = coord_nodes(1,local_nod(0)) ; z(1) = coord_nodes(2,local_nod(0))
        x(2) = coord_nodes(0,local_nod(1)) ; y(2) = coord_nodes(1,local_nod(1)) ; z(2) = coord_nodes(2,local_nod(1))
        x(6) = coord_nodes(0,local_nod(2)) ; y(6) = coord_nodes(1,local_nod(2)) ; z(6) = coord_nodes(2,local_nod(2))
        x(5) = coord_nodes(0,local_nod(3)) ; y(5) = coord_nodes(1,local_nod(3)) ; z(5) = coord_nodes(2,local_nod(3))
    case(3)
        x(3) = coord_nodes(0,local_nod(0)) ; y(3) = coord_nodes(1,local_nod(0)) ; z(3) = coord_nodes(2,local_nod(0))
        x(2) = coord_nodes(0,local_nod(1)) ; y(2) = coord_nodes(1,local_nod(1)) ; z(2) = coord_nodes(2,local_nod(1))
        x(6) = coord_nodes(0,local_nod(2)) ; y(6) = coord_nodes(1,local_nod(2)) ; z(6) = coord_nodes(2,local_nod(2))
        x(7) = coord_nodes(0,local_nod(3)) ; y(7) = coord_nodes(1,local_nod(3)) ; z(7) = coord_nodes(2,local_nod(3))
    case(4)
        x(0) = coord_nodes(0,local_nod(0)) ; y(0) = coord_nodes(1,local_nod(0)) ; z(0) = coord_nodes(2,local_nod(0))
        x(3) = coord_nodes(0,local_nod(1)) ; y(3) = coord_nodes(1,local_nod(1)) ; z(3) = coord_nodes(2,local_nod(1))
        x(7) = coord_nodes(0,local_nod(2)) ; y(7) = coord_nodes(1,local_nod(2)) ; z(7) = coord_nodes(2,local_nod(2))
        x(4) = coord_nodes(0,local_nod(3)) ; y(4) = coord_nodes(1,local_nod(3)) ; z(4) = coord_nodes(2,local_nod(3))
    case(5)
        x(4) = coord_nodes(0,local_nod(0)) ; y(4) = coord_nodes(1,local_nod(0)) ; z(4) = coord_nodes(2,local_nod(0))
        x(5) = coord_nodes(0,local_nod(1)) ; y(5) = coord_nodes(1,local_nod(1)) ; z(5) = coord_nodes(2,local_nod(1))
        x(6) = coord_nodes(0,local_nod(2)) ; y(6) = coord_nodes(1,local_nod(2)) ; z(6) = coord_nodes(2,local_nod(2))
        x(7) = coord_nodes(0,local_nod(3)) ; y(7) = coord_nodes(1,local_nod(3)) ; z(7) = coord_nodes(2,local_nod(3))
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
            where(elem%Jacob(:,:,0) > 0) normal(:,:,i) = -normal(:,:,i)
            end do
        case(1)
            do i = 0,2
                where(elem%Jacob(:,0,:) < 0) normal(:,:,i) = -normal(:,:,i)
                end do
            case(2)
                do i = 0,2
                    where(elem%Jacob(ngllx-1,:,:) < 0) normal(:,:,i) = -normal(:,:,i)
                    end do
                case(3)
                    do i = 0,2
                        where(elem%Jacob(:,nglly-1,:) > 0) normal(:,:,i) = -normal(:,:,i)
                        end do
                    case(4)
                        do i = 0,2
                            where(elem%Jacob(0,:,:) > 0) normal(:,:,i) = -normal(:,:,i)
                            end do
                        case(5)
                            do i = 0,2
                                where(elem%Jacob(:,:,ngllz-1) < 0) normal(:,:,i) = -normal(:,:,i)
                                end do
                            end select

                        end subroutine inversion_normal
                        !------------------------------------------------------------------
                        !------------------------------------------------------------------
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
