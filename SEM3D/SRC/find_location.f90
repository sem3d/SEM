!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file find_location.f90
!!\brief Utilitaire permettant de localiser un point dans le maillage
!!
!!
!<
module mlocations3d
    use sdomain
    use mshape8
    use mshape27
    use selement
    use mindex
    implicit none
contains

    ! Finds the element(s) containing the point (x0,y0,z0)
    ! returns xi,eta,zeta coordinates within each elements
    ! nmax is the maximum number of elements returned, (can be 1) and
    ! is changed by the subroutine to reflect the number of elements found
    ! multiple elements are returned only if the point lies on elements boudaries
    ! The local coords are always within the element except when the point is outside
    ! the mesh in which case nmax is -1 on return and only one candidate is selected
    subroutine find_location(Tdomain, x0,y0,z0, nmax, elems, localcoord, nsrc)
        type(domain), intent(in) :: Tdomain
        real(fpp), intent(in) :: x0, y0, z0
        integer, intent(inout) :: nmax
        integer, dimension(nmax) :: elems
        real(fpp), dimension(0:2,nmax) :: localcoord
        integer, intent(in), optional :: nsrc
        !
        integer :: n, nnodes, i, j, k
        real(fpp) :: mindist, dist
        integer :: best_node
        real(fpp), dimension(0:2) :: p
        real(fpp) :: xi, eta, zeta
        real(fpp), allocatable, dimension(:,:) :: coord
        logical :: ok
        !! Find the closest node
        nnodes = Tdomain%n_nodes
        allocate(coord(0:2, 0:nnodes-1))
        best_node = 0
        mindist = MAX_DOUBLE
        do n = 0, Tdomain%n_glob_nodes-1
            p = Tdomain%Coord_Nodes(:, n)
            dist = (p(0)-x0)**2 + (p(1)-y0)**2 + (p(2)-z0)**2
            if (dist<mindist) then
                best_node = n
                mindist = dist
            end if
        end do
        !! Now find all the elements in contact with the best_node
        k = 0
        do n = 0, Tdomain%n_elem-1
            do i = 0,nnodes-1
                if (Tdomain%specel(n)%Control_Nodes(i)/=best_node) cycle
                !!
                ! Compute local coordinates for the node
                do j = 0, nnodes-1
                    coord(0:2, j) = Tdomain%Coord_Nodes(0:2, Tdomain%specel(n)%Control_Nodes(j))
                enddo
                if (nnodes==8) then
                    call shape8_global2local(coord, x0, y0, z0, xi, eta, zeta, ok)
                else
                    call shape27_global2local(coord, x0, y0, z0, xi, eta, zeta, ok)
                end if
                if (ok) then
                    k = k+1
                    elems(k) = n
                    localcoord(0,k) = xi
                    localcoord(1,k) = eta
                    localcoord(2,k) = zeta
                else
                    write(*,*) "find_location failed:", xi, eta, zeta, ok
                    write(*,*) "find_location for source at :", x0,y0, z0
                endif
                if (k>=nmax) exit
            end do
            if (k>=nmax) exit
        end do
        !!
        if (k<=nmax) nmax = k
        deallocate(coord)
    end subroutine find_location


    subroutine find_location_centroid(Tdomain, x0,y0,z0, nmax, elems, localcoord, nsrc)
        type(domain), intent(in) :: Tdomain
        real(fpp), intent(in) :: x0, y0, z0
        integer, intent(inout) :: nmax
        integer, dimension(nmax) :: elems
        real(fpp), dimension(0:2,nmax) :: localcoord
        integer, intent(in), optional :: nsrc
        integer :: n, nnodes, i, j, k
        real(fpp) :: mindist, distance,dist_elem,dist_faces,dist_edges,dist_nodes
        real(fpp), dimension(0:2) :: p
        real(fpp) :: xi, eta, zeta
        real(fpp), allocatable, dimension(:,:) :: coord
        logical :: ok
        integer :: n_elems
        real(fpp) :: HUGE_NUMBER
        integer, dimension(nmax+1) :: i_elems
        real(fpp), dimension(nmax+1) :: d_elems

        HUGE_NUMBER=1.0e10_fpp

        nnodes = Tdomain%n_nodes

        d_elems=HUGE_NUMBER
        i_elems=-1
        n_elems=0
        do n=0,Tdomain%n_elem-1
            !!! compute distance between source and control nodes
            dist_nodes=HUGE_NUMBER
            do j=0,nnodes-1
                p(0:2)=Tdomain%Coord_Nodes(0:2,Tdomain%specel(n)%Control_Nodes(j))
                distance=sqrt((p(0)-x0)**2+(p(1)-y0)**2+(p(2)-z0)**2)
                if (distance<dist_nodes) dist_nodes=distance
            enddo
            !!! compute distance between source and faces centroids
            dist_edges=HUGE_NUMBER
            do j=0,11
                p=0.0_fpp
                do k=0,1
                    p(0:2)=p(0:2)+Tdomain%Coord_Nodes(0:2,Tdomain%specel(n)%Control_Nodes(edge_def(k,j)))
                enddo
                p=p/2
                distance=sqrt((p(0)-x0)**2+(p(1)-y0)**2+(p(2)-z0)**2)
                if (distance<dist_edges) dist_edges=distance
            enddo
            !!! compute distance between source and faces centroids
            dist_faces=HUGE_NUMBER
            do j=0,5
                p=0.0_fpp
                do k=0,3
                    p(0:2)=p(0:2)+Tdomain%Coord_Nodes(0:2,Tdomain%specel(n)%Control_Nodes(face_def(k,j)))
                enddo
                p=p/4
                distance=sqrt((p(0)-x0)**2+(p(1)-y0)**2+(p(2)-z0)**2)
                if (distance<dist_faces) dist_faces=distance
            enddo
            !!! compute distance between source and element centroid
            p=0.0_fpp
            do j=0,nnodes-1
                p(0:2)=p(0:2)+Tdomain%Coord_Nodes(0:2,Tdomain%specel(n)%Control_Nodes(j))
            enddo
            p=p/nnodes
            dist_elem=sqrt((p(0)-x0)**2+(p(1)-y0)**2+(p(2)-z0)**2)
            !!! minimum distance
            distance=min(dist_elem,dist_faces,dist_edges,dist_nodes)
            !!! insert distance in nearest elems(nmax)
            j=nmax
            do while (j>=1)
                if (d_elems(j)<=distance) exit
                d_elems(j+1)=d_elems(j)
                i_elems(j+1)=i_elems(j)
                j=j-1
            enddo
            d_elems(j+1)=distance
            i_elems(j+1)=n
            if (j<nmax) n_elems=n_elems+1
        enddo
        if (n_elems>nmax) n_elems=nmax

        !!! remove (except one) elements through dxmax criterion
        !!!do i=2,nmax
        !!!    if (d_elems(i)>Tdomain%dxmax) then
        !!!        n_elems=i
        !!!        exit
        !!!    endif
        !!!enddo

        !!! compute local coordinates for source in nearest elems
        allocate(coord(0:2,0:nnodes-1))
        k=0
        do i=1,n_elems
            n=i_elems(i)
            do j=0,nnodes-1
                coord(0:2,j)=Tdomain%Coord_Nodes(0:2,Tdomain%specel(n)%Control_Nodes(j))
            enddo
            if (nnodes==8) then
                call shape8_global2local(coord,x0,y0,z0,xi,eta,zeta,ok)
            else
                call shape27_global2local(coord,x0,y0,z0,xi,eta,zeta,ok)
            endif
            if (ok) then
                k=k+1
                elems(k)=n
                localcoord(0,k)=xi
                localcoord(1,k)=eta
                localcoord(2,k)=zeta
            else
                write(*,*) "find_location failed:",xi,eta,zeta,ok
                write(*,*) "find_location for source at :",x0,y0,z0
            endif
            if (k>=nmax) exit
        enddo
        if (k<=nmax) nmax=k
        deallocate(coord)

    end subroutine find_location_centroid


    subroutine find_coords_elem(Tdomain, NELEMENT)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: NELEMENT
        integer :: nnodes, i
        type(element) :: Elem

        nnodes = Tdomain%n_nodes
        Elem=Tdomain%specel(NELEMENT)
        write(*,*) "control_nodes"
        do i=0,nnodes-1
            write(*,*) "node=",Elem%Control_Nodes(i)
            write(*,*) "    x=",Tdomain%Coord_Nodes(0, Elem%Control_Nodes(i))
            write(*,*) "    y=",Tdomain%Coord_Nodes(1, Elem%Control_Nodes(i))
            write(*,*) "    z=",Tdomain%Coord_Nodes(2, Elem%Control_Nodes(i))
        end do
    end subroutine find_coords_elem
end module mlocations3d

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
