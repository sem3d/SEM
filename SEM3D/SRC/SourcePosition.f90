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

    implicit none

    type (domain), intent(inout) :: Tdomain
    integer, intent(in) :: rg

    integer :: n_src, near_node, i,j,k, n_around_elem, x,y,z, ngllx,nglly,ngllz, code, mat, num, n, i_count
    integer , dimension (0:20) :: el_around_node
    real :: xs,ys,zs,d,dmin,epsil,xa,ya,za,dist_xi,dist_eta,      &
        dist_zeta,xi,eta,zeta,f,ximin,ximax,etamin,etamax,zetamin,zetamax
    !  modif mariotti fevrier 2007 cea
    real :: Xtemp,Ytemp,R
    real, dimension(0:2) :: xi_search,eta_search,zeta_search, centre
    real, dimension(:), allocatable :: distance
    real, dimension(0:2,0:2) :: LocInvGrad
    real, dimension(:,:), allocatable :: coord


    allocate (distance(0:Tdomain%n_proc-1))


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
            if(rg==0) &
                print *,'Position de la source xs,ys,zs',xs,ys,zs !Gsa
        endif

        ! On trouve le noeud le plus proche de la src
        dmin = 1000000
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
            do j = 0,Tdomain%n_nodes-1
                if(Tdomain%specel(i)%Control_nodes(j) == near_node) then
                    el_around_node(n_around_elem) = i
                    n_around_elem = n_around_elem + 1
                endif
            enddo
        enddo

        dmin = 1000000
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
        if(rg==0) &
            print *,'GLL le plus proche de la source ', Tdomain%Ssource(n_src)%elem, Tdomain%Ssource(n_src)%gll(0:2) !Gsa
        ! On trouve le processeur qui contient la src
        call mpi_allgather(dmin,1,mpi_double_precision,distance,1,mpi_double_precision,Tdomain%communicateur,code)
        dmin = 10000000
        do i = 0,Tdomain%n_proc-1
            if (distance(i) < dmin) then
                dmin = distance(i)
                Tdomain%sSource(n_src)%proc = i
            endif
        enddo


        if (rg==Tdomain%sSource(n_src)%proc) then

            ! Dichotomie pour passer de S(x,y,z) a S(xi,eta,zeta)
            n = Tdomain%sSource(n_src)%elem
            mat = Tdomain%specel(n)%mat_index
            x = Tdomain%sSource(n_src)%gll(0)
            y = Tdomain%sSource(n_src)%gll(1)
            z = Tdomain%sSource(n_src)%gll(2)
            ximin = merge(Tdomain%sSubdomain(mat)%GLLcx(x),Tdomain%sSubdomain(mat)%GLLcx(x-1),x == 0)
            etamin = merge(Tdomain%sSubdomain(mat)%GLLcy(y),Tdomain%sSubdomain(mat)%GLLcy(y-1),y == 0)
            zetamin = merge(Tdomain%sSubdomain(mat)%GLLcz(z),Tdomain%sSubdomain(mat)%GLLcz(z-1),z == 0)
            ximax = merge(Tdomain%sSubdomain(mat)%GLLcx(x),Tdomain%sSubdomain(mat)%GLLcx(x+1),   &
                x == Tdomain%sSubdomain(mat)%ngllx-1)
            etamax = merge(Tdomain%sSubdomain(mat)%GLLcy(y),Tdomain%sSubdomain(mat)%GLLcy(y+1),   &
                y == Tdomain%sSubdomain(mat)%nglly-1)
            zetamax = merge(Tdomain%sSubdomain(mat)%GLLcz(z),Tdomain%sSubdomain(mat)%GLLcz(z+1),   &
                z == Tdomain%sSubdomain(mat)%ngllz-1)
            dist_xi = ximax-ximin
            dist_eta = etamax-etamin
            dist_zeta = zetamax-zetamin
            do i = 0,2
                xi_search(i) = ximin + (i+1)*dist_xi/4
                eta_search(i) = etamin + (i+1)*dist_eta/4
                zeta_search(i) = zetamin + (i+1)*dist_zeta/4
            enddo
            allocate(coord(0:Tdomain%n_nodes-1,0:2))
            do i = 0,Tdomain%n_nodes-1
                j = Tdomain%specel(n)%Control_Nodes(i)
                coord(i,0:2) = Tdomain%Coord_Nodes(0:2,j)
            enddo
            dmin = max(sqrt((coord(0,0)-coord(6,0))**2 + (coord(0,1)-coord(6,1))**2 + (coord(0,2)-coord(6,2))**2), &
                sqrt((coord(1,0)-coord(7,0))**2 + (coord(1,1)-coord(7,1))**2 + (coord(1,2)-coord(7,2))**2), &
                sqrt((coord(2,0)-coord(4,0))**2 + (coord(2,1)-coord(4,1))**2 + (coord(2,2)-coord(4,2))**2), &
                sqrt((coord(3,0)-coord(5,0))**2 + (coord(3,1)-coord(5,1))**2 + (coord(3,2)-coord(5,2))**2))
            epsil = dmin/10000.
            num = 0
            dicho : do while(dmin > epsil)
                do i = 0,2
                    do j = 0,2
                        do k = 0,2
                            ! calcul des coords dans l'espace physique
                            xi = xi_search(i);   eta = eta_search(j);   zeta = zeta_search(k)
                            if (Tdomain%n_nodes==27) then
                                xa = 0;   ya = 0;   za = 0;
                                do i_count = 0,Tdomain%n_nodes-1
                                    f = Comp_shapefunc(i_count,xi,eta,zeta)
                                    xa = xa + coord(i_count,0)*f
                                    ya = ya + coord(i_count,1)*f
                                    za = za + coord(i_count,2)*f
                                enddo
                            else if (Tdomain%n_nodes==8) then
                                xa = 0.125 * (coord(0,0)*(1-xi)*(1-eta)*(1-zeta) + coord(1,0)*(1+xi)*(1-eta)*(1-zeta) + &
                                    coord(2,0)*(1+xi)*(1+eta)*(1-zeta) + coord(3,0)*(1-xi)*(1+eta)*(1-zeta) + &
                                    coord(4,0)*(1-xi)*(1-eta)*(1+zeta) + coord(5,0)*(1+xi)*(1-eta)*(1+zeta) + &
                                    coord(6,0)*(1+xi)*(1+eta)*(1+zeta) + coord(7,0)*(1-xi)*(1+eta)*(1+zeta))
                                ya = 0.125 * (coord(0,1)*(1-xi)*(1-eta)*(1-zeta) + coord(1,1)*(1+xi)*(1-eta)*(1-zeta) + &
                                    coord(2,1)*(1+xi)*(1+eta)*(1-zeta) + coord(3,1)*(1-xi)*(1+eta)*(1-zeta) + &
                                    coord(4,1)*(1-xi)*(1-eta)*(1+zeta) + coord(5,1)*(1+xi)*(1-eta)*(1+zeta) + &
                                    coord(6,1)*(1+xi)*(1+eta)*(1+zeta) + coord(7,1)*(1-xi)*(1+eta)*(1+zeta))
                                za = 0.125 * (coord(0,2)*(1-xi)*(1-eta)*(1-zeta) + coord(1,2)*(1+xi)*(1-eta)*(1-zeta) + &
                                    coord(2,2)*(1+xi)*(1+eta)*(1-zeta) + coord(3,2)*(1-xi)*(1+eta)*(1-zeta) + &
                                    coord(4,2)*(1-xi)*(1-eta)*(1+zeta) + coord(5,2)*(1+xi)*(1-eta)*(1+zeta) + &
                                    coord(6,2)*(1+xi)*(1+eta)*(1+zeta) + coord(7,2)*(1-xi)*(1+eta)*(1+zeta))
                            endif
                            ! calcul de la distance a la source
                            d = sqrt((xa-xs)**2+(ya-ys)**2+(za-zs)**2)
                            if(d < dmin)then
                                dmin = d
                                centre(0) = xi ; centre(1) = eta ; centre(2) = zeta
                            endif
                        enddo
                    enddo
                enddo
                dist_xi = dist_xi/2.
                dist_eta = dist_eta/2.
                dist_zeta = dist_zeta/2.
                do i = 0,2
                    xi_search(i) = centre(0) + (i-1)*dist_xi/4.
                    eta_search(i) = centre(1) + (i-1)*dist_eta/4.
                    zeta_search(i) = centre(2) + (i-1)*dist_zeta/4.
                enddo
                num = num + 1
                if (num>20) exit dicho
            enddo dicho
            Tdomain%sSource(n_src)%refcoord = centre

            print*,'Finalement la position dans l''elt de ref est',centre(0:2) !Gsa
            ! Calcul de InvGrad en S(xi,eta,zeta)
            xi = centre(0);   eta = centre(1);   zeta = centre(2)
            if (Tdomain%n_nodes==27) then
                LocInvGrad = 0.
                do i = 0,Tdomain%n_nodes-1
                    do j = 0,2
                        f = Comp_derivshapefunc(i,xi,eta,zeta,j)
                        LocInvGrad(j,0) = LocInvGrad(j,0) + coord(i,0)*f
                        LocInvGrad(j,1) = LocInvGrad(j,1) + coord(i,1)*f
                        LocInvGrad(j,2) = LocInvGrad(j,2) + coord(i,2)*f
                    enddo
                enddo
            else if (Tdomain%n_nodes==8) then
                LocInvGrad(0,0) = 0.125 * ((coord(1,0)-coord(0,0))*(1-eta)*(1-zeta) + &
                                           (coord(2,0)-coord(3,0))*(1+eta)*(1-zeta) + &
                                           (coord(5,0)-coord(4,0))*(1-eta)*(1+zeta) + &
                                           (coord(6,0)-coord(7,0))*(1+eta)*(1+zeta))
                LocInvGrad(1,0) = 0.125 * ((coord(3,0)-coord(0,0))*(1-xi)*(1-zeta) + (coord(2,0)-coord(1,0))*(1+xi)*(1-zeta) + &
                    (coord(7,0)-coord(4,0))*(1-xi)*(1+zeta) + (coord(6,0)-coord(5,0))*(1+xi)*(1+zeta))
                LocInvGrad(2,0) = 0.125 * ((coord(4,0)-coord(0,0))*(1-xi)*(1-eta) + (coord(5,0)-coord(1,0))*(1+xi)*(1-eta) + &
                    (coord(7,0)-coord(3,0))*(1-xi)*(1+eta) + (coord(6,0)-coord(2,0))*(1+xi)*(1+eta))
                LocInvGrad(0,1) = 0.125 * ((coord(1,1)-coord(0,1))*(1-eta)*(1-zeta) + (coord(2,1)-coord(3,1))*(1+eta)*(1-zeta) + &
                    (coord(5,1)-coord(4,1))*(1-eta)*(1+zeta) + (coord(6,1)-coord(7,1))*(1+eta)*(1+zeta))
                LocInvGrad(1,1) = 0.125 * ((coord(3,1)-coord(0,1))*(1-xi)*(1-zeta) + (coord(2,1)-coord(1,1))*(1+xi)*(1-zeta) + &
                    (coord(7,1)-coord(4,1))*(1-xi)*(1+zeta) + (coord(6,1)-coord(5,1))*(1+xi)*(1+zeta))
                LocInvGrad(2,1) = 0.125 * ((coord(4,1)-coord(0,1))*(1-xi)*(1-eta) + (coord(5,1)-coord(1,1))*(1+xi)*(1-eta) + &
                    (coord(7,1)-coord(3,1))*(1-xi)*(1+eta) + (coord(6,1)-coord(2,1))*(1+xi)*(1+eta))
                LocInvGrad(0,2) = 0.125 * ((coord(1,2)-coord(0,2))*(1-eta)*(1-zeta) + (coord(2,2)-coord(3,2))*(1+eta)*(1-zeta) + &
                    (coord(5,2)-coord(4,2))*(1-eta)*(1+zeta) + (coord(6,2)-coord(7,2))*(1+eta)*(1+zeta))
                LocInvGrad(1,2) = 0.125 * ((coord(3,2)-coord(0,2))*(1-xi)*(1-zeta) + (coord(2,2)-coord(1,2))*(1+xi)*(1-zeta) + &
                    (coord(7,2)-coord(4,2))*(1-xi)*(1+zeta) + (coord(6,2)-coord(5,2))*(1+xi)*(1+zeta))
                LocInvGrad(2,2) = 0.125 * ((coord(4,2)-coord(0,2))*(1-xi)*(1-eta) + (coord(5,2)-coord(1,2))*(1+xi)*(1-eta) + &
                    (coord(7,2)-coord(3,2))*(1-xi)*(1+eta) + (coord(6,2)-coord(2,2))*(1+xi)*(1+eta))
            endif
            call invert_3d (LocInvGrad, R)
            Tdomain%sSource(n_src)%InvGrad(0:2,0:2) = LocInvGrad(0:2,0:2)
            deallocate (coord)

        endif


    enddo

    deallocate (distance)

    return
end subroutine SourcePosition
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
