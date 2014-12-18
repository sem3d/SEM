!>
!! \file SourcePosition.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine SourcePosition (Tdomain,rg)

    use sdomain
    use constants, only : M_PI
    use mpi
    use mshape8
    use mshape27
    implicit none

    type (domain), intent(inout) :: Tdomain
    integer, intent(in) :: rg

    integer :: n_src, near_node, i,j, n_around_elem, x,y,z, ngllx,nglly,ngllz, code, mat, n
    integer , dimension (0:20) :: el_around_node
    real :: xs,ys,zs,d,dmin,xi,eta,zeta
    !  modif mariotti fevrier 2007 cea
    real :: Xtemp,Ytemp,R
    real, dimension(:), allocatable :: distance
    real, dimension(0:2,0:2) :: LocInvGrad
    real, dimension(:,:), allocatable :: coord
    real, parameter :: DISTMAX=10000000
    integer :: src_proc


    allocate(distance(0:Tdomain%nb_procs-1))
    allocate(coord(0:2, 0:Tdomain%n_nodes-1))


    do n_src = 0, Tdomain%n_source-1


        ! On calcule les coordonnees des src dans le chunk de reference
        if (Tdomain%curve) then
            !  modif mariotti fevrier 2007 cea
            ! Passage de latitude a colatitude dans le chunk reel
            !    Tdomain%Ssource(n_src)%realcolat = 90 - Tdomain%Ssource(n_src)%realcolat
            !    if (Tdomain%Ssource(n_src)%reallong<0) Tdomain%Ssource(n_src)%reallong = 360 + Tdomain%Ssource(n_src)%reallong
            !    colat = Tdomain%Ssource(n_src)%realcolat
            !    long = Tdomain%Ssource(n_src)%reallong
            ! Passage de spherique a cartesien dans le chunk reel
            !    R = Tdomain%Ssource(n_src)%radius
            !    ct = dcos(pi*colat/180)
            !    st = dsin(pi*colat/180)
            !    cp = dcos(pi*long/180)
            !    sp = dsin(pi*long/180)
            !    xa = R * st * cp
            !    ya = R * st * sp
            !    za = R * ct
            ! Passage dans le chunk de reference
            !    tRot = transpose(Tdomain%rot)
            !    xs = tRot(0,0)*xa + tRot(0,1)*ya + tRot(0,2)*za
            !    ys = tRot(1,0)*xa + tRot(1,1)*ya + tRot(1,2)*za
            !    zs = tRot(2,0)*xa + tRot(2,1)*ya + tRot(2,2)*za
            ! Passage de cartesien a spherique dans le chunk de reference
            !    call cart2sph(xs,ys,zs,R,colat,long)
            !    Tdomain%Ssource(n_src)%refcolat = colat
            !    Tdomain%Ssource(n_src)%reflong = long
            R = Tdomain%Ssource(n_src)%Zsource
            Xtemp = tan(M_PI*Tdomain%Ssource(n_src)%Xsource/180)
            Ytemp = tan(M_PI*Tdomain%Ssource(n_src)%Ysource/180)
            D = sqrt(1 + Ytemp**2 + Xtemp**2)
            xs = R*Xtemp/D
            ys = R*Ytemp/D
            zs = R/D

        else
            xs = Tdomain%Ssource(n_src)%Xsource
            ys = Tdomain%Ssource(n_src)%Ysource
            zs = Tdomain%Ssource(n_src)%Zsource
            !if(rg==0) &
            !    print *,'Position de la source xs,ys,zs',xs,ys,zs !Gsa
        endif

        ! On trouve le noeud le plus proche de la src
        dmin = DISTMAX
        do i = 0,Tdomain%n_glob_nodes-1
            d = sqrt((Tdomain%Coord_nodes(0,i)-xs)**2 + (Tdomain%Coord_nodes(1,i)-ys)**2 + (Tdomain%Coord_nodes(2,i)-zs)**2)
            if (d <= dmin) then
                dmin = d
                near_node = i
            endif
        enddo

        ! On trouve l'element dans lequel est la src ainsi que le GLL interne le plus proche
        n_around_elem = 0
        do i = 0,Tdomain%n_elem-1
            if(Tdomain%sSource(n_src)%i_type_source == 3 .and.    &
                  Tdomain%specel(i)%solid)   cycle
            if(Tdomain%sSource(n_src)%i_type_source /= 3 .and.    &
                  .not. Tdomain%specel(i)%solid)   cycle
            do j = 0,Tdomain%n_nodes-1
                if(Tdomain%specel(i)%Control_nodes(j) == near_node) then
                    el_around_node(n_around_elem) = i
                    n_around_elem = n_around_elem + 1
                endif
            enddo
        enddo

        dmin = DISTMAX
        do i = 0,n_around_elem-1
            ngllx = Tdomain%specel(el_around_node(i))%ngllx
            nglly = Tdomain%specel(el_around_node(i))%nglly
            ngllz = Tdomain%specel(el_around_node(i))%ngllz
            if ( Tdomain%Ssource(n_src)%i_type_source==1) then
                do x = 0,ngllx-1
                    do y = 0,nglly-1
                        do z = 0,ngllz-1
                            j = Tdomain%specel(el_around_node(i))%Iglobnum(x,y,z)
                            d = sqrt((Tdomain%Globcoord(0,j)-xs)**2+(Tdomain%Globcoord(1,j)-ys)**2 + (Tdomain%Globcoord(2,j)-zs)**2)
                            if (d <= dmin) then
                                dmin = d
                                Tdomain%Ssource(n_src)%elem = el_around_node(i)
                                Tdomain%Ssource(n_src)%gll(0) = x
                                Tdomain%Ssource(n_src)%gll(1) = y
                                Tdomain%Ssource(n_src)%gll(2) = z
                            endif
                        enddo
                    enddo
                enddo
            else
                do x = 1,ngllx-2
                    do y = 1,nglly-2
                        do z = 1,ngllz-2
                            j = Tdomain%specel(el_around_node(i))%Iglobnum(x,y,z)
                            d = sqrt ((Tdomain%Globcoord(0,j)-xs)**2 + (Tdomain%Globcoord(1,j)-ys)**2 + (Tdomain%Globcoord(2,j)-zs)**2)
                            if (d <= dmin) then
                                dmin = d
                                Tdomain%Ssource(n_src)%elem = el_around_node(i)
                                Tdomain%Ssource(n_src)%gll(0) = x
                                Tdomain%Ssource(n_src)%gll(1) = y
                                Tdomain%Ssource(n_src)%gll(2) = z
                            endif
                        enddo
                    enddo
                enddo
            endif
        enddo
        ! On trouve le processeur qui contient la src
        call mpi_allgather(dmin,1,mpi_double_precision,distance,1,mpi_double_precision,Tdomain%communicateur,code)
        dmin = DISTMAX
        src_proc = -1
        do i = 0,Tdomain%nb_procs-1
            if (distance(i) < dmin) then
                dmin = distance(i)
                src_proc = i
            endif
        enddo
        write(*,*) "Found source on proc", src_proc
        if (src_proc==-1) then
            if (rg==0) then
                write(*,*) "A source doesn't appear to be on any processor"
                write(*,*) "Please verify that the source location is within the computation domain"
                write(*,*) "And the source type (solid, or fluid) is in accordance with its location"
            end if
            stop 1
        endif
        Tdomain%sSource(n_src)%proc = src_proc

        if (rg==Tdomain%sSource(n_src)%proc) then
            n = Tdomain%sSource(n_src)%elem
            do i = 0,Tdomain%n_nodes-1
                j = Tdomain%specel(n)%Control_Nodes(i)
                coord(0:2,i) = Tdomain%Coord_Nodes(0:2,j)
            enddo


            ! Dichotomie pour passer de S(x,y,z) a S(xi,eta,zeta)
            mat = Tdomain%specel(n)%mat_index

            if (Tdomain%n_nodes==27) then
                call shape27_global2local(coord, xs, ys, zs, xi, eta, zeta)
                call shape27_local2jacob(coord, xi, eta, zeta, LocInvGrad)
            else if (Tdomain%n_nodes==8) then
                call shape8_global2local(coord, xs, ys, zs, xi, eta, zeta)
                call shape8_local2jacob(coord, xi, eta, zeta, LocInvGrad)
            endif
            Tdomain%sSource(n_src)%refcoord(0) = xi
            Tdomain%sSource(n_src)%refcoord(1) = eta
            Tdomain%sSource(n_src)%refcoord(2) = zeta
            write(*,*) "Source local coordinates:", Tdomain%sSource(n_src)%refcoord
            write(*,*) "From position", xs, ys, zs
            call invert_3d (LocInvGrad, R)
            Tdomain%sSource(n_src)%InvGrad(0:2,0:2) = LocInvGrad(0:2,0:2)
        endif
    enddo

    deallocate (coord)
    deallocate (distance)

    return
end subroutine SourcePosition
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
