!>
!! \file ReceiverPosition.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine ReceiverPosition (Tdomain,rg)

    ! Written by Paul Cupillard 06/06/2005

    use sdomain
    use constants, only : M_PI
    use mpi

    implicit none

    type (domain), intent (inout) :: Tdomain
    integer, intent(in) :: rg

    integer :: n_rcp, near_node, n_around_elem, i,j,k, x,y,z, ngllx,nglly,ngllz, code, mat, num, n, a, i_count
    integer , dimension (0:20) :: el_around_node
    real :: xs,ys,zs, d,dmin, R,epsil,xa,ya,za, &
        dist_xi,dist_eta,dist_zeta, xi,eta,zeta, f
    !  real :: colat,long, ct,st,cp,sp, coeff, cosdelta,sindelta
    !  integer :: n_src
    !  modif mariotti fevrier 2007 cea
    real :: Xtemp,Ytemp
    real, dimension(0:2) :: xi_search,eta_search,zeta_search, centre
    real, dimension(:), allocatable :: distance
    !real, dimension(0:2,0:2) :: tRot
    real, dimension(:,:), allocatable :: coord


    allocate (distance(0:Tdomain%n_proc-1))
    allocate (coord(0:Tdomain%n_nodes-1,0:2))

    do n_rcp = 0, Tdomain%n_receivers-1


        ! On calcule les coordonnees du rcp dans le chunk de reference
        if (Tdomain%curve) then
            !  modif mariotti fevrier 2007 cea
            ! Passage de latitude a colatitude dans le chunk reel
            !    Tdomain%sReceiver(n_rcp)%realcolat = 90 - Tdomain%sReceiver(n_rcp)%realcolat
            !    if (Tdomain%sReceiver(n_rcp)%reallong < 0.) Tdomain%sReceiver(n_rcp)%reallong = 360 + Tdomain%sReceiver(n_rcp)%reallong
            !    colat = Tdomain%sReceiver(n_rcp)%realcolat
            !    long = Tdomain%sReceiver(n_rcp)%reallong
            ! Passage de spherique a cartesien dans le chunk reel
            !    R = Tdomain%sReceiver(n_rcp)%radius
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
            !    Tdomain%sReceiver(n_rcp)%refcolat = colat
            !    Tdomain%sReceiver(n_rcp)%reflong = long
            ! Calcul de la matrice de passage du systeme cartesien au systeme Z/N/E
            !    ct = dcos(colat)
            !    st = dsin(colat)
            !    cp = dcos(long)
            !    sp = dsin(long)
            !    Tdomain%sReceiver(n_rcp)%Passage(0,0) = st*cp; Tdomain%sReceiver(n_rcp)%Passage(0,1) = st*sp; Tdomain%sReceiver(n_rcp)%Passage(0,2) = ct
            !    Tdomain%sReceiver(n_rcp)%Passage(1,0) =-ct*cp; Tdomain%sReceiver(n_rcp)%Passage(1,1) =-ct*sp; Tdomain%sReceiver(n_rcp)%Passage(1,2) = st
            !    Tdomain%sReceiver(n_rcp)%Passage(2,0) = -sp  ; Tdomain%sReceiver(n_rcp)%Passage(2,1) = cp   ; Tdomain%sReceiver(n_rcp)%Passage(2,2) = 0
            ! Calcul de l'angle de passage du systeme Z/N/E au systeme Z/R/T
            !    allocate (Tdomain%sReceiver(n_rcp)%cosgamma(0:Tdomain%n_source-1))
            !    allocate (Tdomain%sReceiver(n_rcp)%singamma(0:Tdomain%n_source-1))
            !    do n_src = 0, Tdomain%n_source-1
            !       call arc_gd_cercle(long,Tdomain%sSource(n_src)%reflong,colat,Tdomain%sSource(n_src)%refcolat,cosdelta,sindelta)
            !       Tdomain%sReceiver(n_rcp)%cosgamma(n_src) = (dcos(Tdomain%sSource(n_src)%refcolat) - cosdelta*ct) / (sindelta*st)
            !       Tdomain%sReceiver(n_rcp)%singamma(n_src) = dsqrt(1-Tdomain%sReceiver(n_rcp)%cosgamma(n_src)**2)
            !    enddo
            R = Tdomain%sReceiver(n_rcp)%zRec
            Xtemp = tan(M_PI*Tdomain%sReceiver(n_rcp)%xRec/180)
            Ytemp = tan(M_PI*Tdomain%sReceiver(n_rcp)%yRec/180)
            D = sqrt(1 + Ytemp**2 + Xtemp**2)
            xs = R*Xtemp/D
            ys = R*Ytemp/D
            zs = R/D

        else
            xs = Tdomain%sReceiver(n_rcp)%xRec
            ys = Tdomain%sReceiver(n_rcp)%yRec
            zs = Tdomain%sReceiver(n_rcp)%zRec
        endif

        ! On trouve le noeud le plus proche du rcp
        dmin = 1000000
        do i = 0,Tdomain%n_glob_nodes-1
            d = sqrt ((Tdomain%Coord_nodes(0,i)-xs)**2 + (Tdomain%Coord_nodes(1,i)-ys)**2 + (Tdomain%Coord_nodes(2,i)-zs)**2)
            if (d <= dmin) then
                dmin = d
                near_node = i
            endif
        enddo

        ! On trouve l'element dans lequel est le rcp ainsi que le GLL interne le plus proche
        n_around_elem = 0
        do i = 0,Tdomain%n_elem-1
            do j = 0,Tdomain%n_nodes-1
                if (Tdomain%specel(i)%Control_nodes(j)==near_node) then
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
            do x = 1,ngllx-2
                do y = 1,nglly-2
                    do z = 1,ngllz-2
                        j = Tdomain%specel(el_around_node(i))%Iglobnum(x,y,z)
                        d = sqrt ((Tdomain%Globcoord(0,j)-xs)**2 + (Tdomain%Globcoord(1,j)-ys)**2 + (Tdomain%Globcoord(2,j)-zs)**2)
                        if (d <= dmin) then
                            dmin = d
                            Tdomain%sReceiver(n_rcp)%elem = el_around_node(i)
                            Tdomain%sReceiver(n_rcp)%gll(0) = x
                            Tdomain%sReceiver(n_rcp)%gll(1) = y
                            Tdomain%sReceiver(n_rcp)%gll(2) = z
                        endif
                    enddo
                enddo
            enddo
        enddo

        ! On trouve le processeur qui contient le rcp
        call mpi_allgather(dmin,1,mpi_double_precision,distance,1,mpi_double_precision,Tdomain%communicateur,code)

        do i = 0,Tdomain%n_proc-1
            if (distance(i) <= dmin) then
                dmin = distance(i)
                Tdomain%sReceiver(n_rcp)%proc = i
            endif
        enddo


        if (rg==Tdomain%sReceiver(n_rcp)%proc) then

            ! Dichotomie pour passer de S(x,y,z) a S(xi,eta,zeta)
            n = Tdomain%sReceiver(n_rcp)%elem
            mat = Tdomain%specel(n)%mat_index
            x = Tdomain%sReceiver(n_rcp)%gll(0)
            y = Tdomain%sReceiver(n_rcp)%gll(1)
            z = Tdomain%sReceiver(n_rcp)%gll(2)
            dist_xi = Tdomain%sSubdomain(mat)%GLLcx(x+1) - Tdomain%sSubdomain(mat)%GLLcx(x-1)
            dist_eta = Tdomain%sSubdomain(mat)%GLLcy(y+1) - Tdomain%sSubdomain(mat)%GLLcy(y-1)
            dist_zeta = Tdomain%sSubdomain(mat)%GLLcz(z+1) - Tdomain%sSubdomain(mat)%GLLcz(z-1)
            do i = 0,2
                xi_search(i) = Tdomain%sSubdomain(mat)%GLLcx(x-1) + (i+1)*dist_xi/4
                eta_search(i) = Tdomain%sSubdomain(mat)%GLLcy(y-1) + (i+1)*dist_eta/4
                zeta_search(i) = Tdomain%sSubdomain(mat)%GLLcz(z-1) + (i+1)*dist_zeta/4
            enddo

            do i = 0,Tdomain%n_nodes-1
                j = Tdomain%specel(n)%Control_Nodes(i)
                coord(i,0:2) = Tdomain%Coord_Nodes(0:2,j)
            enddo
            dmin = max (sqrt((coord(0,0)-coord(6,0))**2 + (coord(0,1)-coord(6,1))**2 + (coord(0,2)-coord(6,2))**2), &
                sqrt((coord(1,0)-coord(7,0))**2 + (coord(1,1)-coord(7,1))**2 + (coord(1,2)-coord(7,2))**2), &
                sqrt((coord(2,0)-coord(4,0))**2 + (coord(2,1)-coord(4,1))**2 + (coord(2,2)-coord(4,2))**2), &
                sqrt((coord(3,0)-coord(5,0))**2 + (coord(3,1)-coord(5,1))**2 + (coord(3,2)-coord(5,2))**2))
            epsil = dmin/10000.
            num = 0
            dicho : do while (dmin > epsil)
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
                            d = sqrt ((xa-xs)**2 + (ya-ys)**2 + (za-zs)**2)
                            if (d < dmin) then
                                dmin = d
                                centre(0) = xi;   centre(1) = eta;   centre(2) = zeta
                            endif
                        enddo
                    enddo
                enddo
                dist_xi = dist_xi/2.
                dist_eta = dist_eta/2.
                dist_zeta = dist_zeta/2.
                do a = 0,2
                    xi_search(a) = centre(0) + (a-1)*dist_xi/4.
                    eta_search(a) = centre(1) + (a-1)*dist_eta/4.
                    zeta_search(a) = centre(2) + (a-1)*dist_zeta/4.
                enddo
                num = num + 1
                if (num>20) exit dicho
            enddo dicho

            ! Calcul en S(xi,eta,zeta) des polynomes de Lagrange associes aux GLLs
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            allocate (Tdomain%sReceiver(n_rcp)%coeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
            allocate (Tdomain%sReceiver(n_rcp)%pol(0:ngllx-1,0:nglly-1,0:ngllz-1))
            Tdomain%sReceiver(n_rcp)%pol = 1

            do x = 0,ngllx-1
                do y = 0,nglly-1
                    do z = 0,ngllz-1
                        xi = Tdomain%sSubdomain(mat)%GLLcx(x)
                        eta = Tdomain%sSubdomain(mat)%GLLcy(y)
                        zeta = Tdomain%sSubdomain(mat)%GLLcz(z)
                        do a = 0,ngllx-1
                            if (a/=x) then
                                xa = Tdomain%sSubdomain(mat)%GLLcx(a)
                                Tdomain%sReceiver(n_rcp)%pol(x,y,z) = Tdomain%sReceiver(n_rcp)%pol(x,y,z) * (centre(0)-xa)/(xi-xa)
                            endif
                        enddo
                        do a = 0,nglly-1
                            if (a/=y) then
                                ya = Tdomain%sSubdomain(mat)%GLLcy(a)
                                Tdomain%sReceiver(n_rcp)%pol(x,y,z) = Tdomain%sReceiver(n_rcp)%pol(x,y,z) * (centre(1)-ya)/(eta-ya)
                            endif
                        enddo
                        do a = 0,ngllz-1
                            if (a/=z) then
                                za = Tdomain%sSubdomain(mat)%GLLcz(a)
                                Tdomain%sReceiver(n_rcp)%pol(x,y,z) = Tdomain%sReceiver(n_rcp)%pol(x,y,z) * (centre(2)-za)/(zeta-za)
                            endif
                        enddo
                    enddo
                enddo
            enddo

        endif


    enddo

    deallocate (distance,coord)

    return
end subroutine ReceiverPosition
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
