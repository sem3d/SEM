!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file SourcePosition.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine SourcePosition (Tdomain)
    use sdomain
    use mpi
    use mshape8
    use mshape27
    use mlocations3d
    implicit none
    type (domain), intent(inout) :: Tdomain
    !
    integer :: rg, ierr, nnodes
    integer :: n_src, i,j
    double precision :: xs,ys,zs,xi,eta,zeta
    double precision, dimension(0:2,0:2) :: LocInvGrad
    double precision, dimension(:,:), allocatable :: coord
    double precision :: R
    integer :: src_proc, nmax, n_el
    integer, parameter :: NMAXEL=20
    integer, dimension(NMAXEL) :: elems
    double precision, dimension(0:2,NMAXEL) :: coordloc
    double precision, parameter :: EPS = 1D-13
    logical :: inside

    rg = Tdomain%rank
    nnodes = Tdomain%n_nodes
    allocate(coord(0:2, 0:nnodes-1))

    do n_src = 0, Tdomain%n_source-1
        xs = Tdomain%Ssource(n_src)%Xsource
        ys = Tdomain%Ssource(n_src)%Ysource
        zs = Tdomain%Ssource(n_src)%Zsource

        ! On trouve les elements autour de la source
        nmax = NMAXEL
        call find_location(Tdomain, xs, ys, zs, nmax, elems, coordloc)
        n_el = -1
        do i=1,nmax
            inside = .true.
            xi   = coordloc(0,i)
            eta  = coordloc(1,i)
            zeta = coordloc(2,i)
            if (xi<(-1-EPS) .or. eta<(-1-EPS) .or. zeta<(-1-EPS)) inside = .false.
            if (xi>(1+EPS) .or. eta>(1+EPS) .or. zeta>(1+EPS)) inside = .false.
            if (inside) then
                n_el = elems(i)
                exit
            end if
        end do
        ! On ignore une source fluide dans le domaine solide
        if (n_el/=-1) then
            if(Tdomain%sSource(n_src)%i_type_source == 3 .and. Tdomain%specel(elems(i))%domain/=DM_FLUID) n_el = -1
            if(Tdomain%sSource(n_src)%i_type_source /= 3 .and. Tdomain%specel(elems(i))%domain/=DM_SOLID) n_el = -1
        endif

        Tdomain%Ssource(n_src)%elem = n_el
        if (n_el/=-1) then
            Tdomain%sSource(n_src)%proc = rg
        else
            Tdomain%sSource(n_src)%proc = -1
        end if        
        ! On trouve le processeur qui contient la src
        call MPI_AllReduce(Tdomain%sSource(n_src)%proc, src_proc, 1, MPI_INTEGER, &
                           MPI_MAX, Tdomain%communicateur, ierr)

        if (src_proc==-1) then
            if (rg==0) then
                write(*,*) "The source at (", xs, ys, zs, ") doesn't appear to be on any processor"
                write(*,*) "Please verify that the source location is within the computation domain"
                write(*,*) "And the source type (solid, or fluid) is in accordance with its location"
            end if
            stop 1
        endif

        if (rg==src_proc) then
            n_el = Tdomain%sSource(n_src)%elem
            do j = 0, nnodes-1
                coord(0:2, j) = Tdomain%Coord_Nodes(0:2, Tdomain%specel(n_el)%Control_Nodes(j))
            enddo
            if (nnodes==27) then
                call shape27_local2jacob(coord, xi, eta, zeta, LocInvGrad)
            else if (nnodes==8) then
                call shape8_local2jacob(coord, xi, eta, zeta, LocInvGrad)
            endif
            Tdomain%sSource(n_src)%refcoord(0) = xi
            Tdomain%sSource(n_src)%refcoord(1) = eta
            Tdomain%sSource(n_src)%refcoord(2) = zeta
            !write(*,*) "Found source", n_src, "on proc", src_proc
            !write(*,*) "Source local coordinates:", Tdomain%sSource(n_src)%refcoord
            !write(*,*) "From position", xs, ys, zs
            call invert_3d (LocInvGrad, R)
            Tdomain%sSource(n_src)%InvGrad(0:2,0:2) = LocInvGrad(0:2,0:2)
        endif
    enddo
end subroutine SourcePosition

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

