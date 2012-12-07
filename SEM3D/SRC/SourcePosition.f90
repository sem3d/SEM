subroutine SourcePosition(Tdomain,nb_procs,rg)

    use sdomain
    use pig
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: nb_procs,rg
    integer :: n_src,i,j,k,n_around_elem,x,y,z,ngllx,nglly,ngllz,code,  &
               mat,num,n,i_count,near_node
    integer, dimension(0:20) :: el_around_node
    real :: xs,ys,zs,d,dmin,trash,epsil,xa,ya,za,dist_xi,dist_eta,dist_zeta,xi,eta,zeta,f
    real, dimension(0:nb_procs-1) :: distance
    real, dimension(0:2,0:2) :: LocInvGrad
    real, dimension(0:2) :: xi_search,eta_search,zeta_search, centre
    real, dimension(:,:), allocatable :: coord



do n_src = 0, Tdomain%n_source-1
  ! source co-ordinates
    xs = Tdomain%Ssource(n_src)%Xsource
    ys = Tdomain%Ssource(n_src)%Ysource
    zs = Tdomain%Ssource(n_src)%Zsource

 dmin = 1000d0
 do i = 0,Tdomain%n_glob_nodes-1
    d = sqrt((Tdomain%Coord_nodes(0,i)-xs)**2 + (Tdomain%Coord_nodes(1,i)-ys)**2 + (Tdomain%Coord_nodes(2,i)-zs)**2)
    if (d <= dmin) then
       dmin = d
       near_node = i
    endif
 enddo
 n_around_elem = 0
 do i = 0,Tdomain%n_elem-1
    do j = 0,Tdomain%n_nodes-1
       if(Tdomain%specel(i)%Control_nodes(j) == near_node) then
          el_around_node(n_around_elem) = i
          n_around_elem = n_around_elem + 1
       endif
    enddo
 enddo
 do i = 0,n_around_elem-1
    ngllx = Tdomain%specel(el_around_node(i))%ngllx
    nglly = Tdomain%specel(el_around_node(i))%nglly
    ngllz = Tdomain%specel(el_around_node(i))%ngllz
    do x = 0,ngllx-1
     do y = 0,nglly-1
      do z = 0,ngllz-1
         j = Tdomain%specel(el_around_node(i))%Iglobnum(x,y,z)
         d = sqrt((Tdomain%Globcoord(0,j)-xs)**2+(Tdomain%Globcoord(1,j)-ys)**2+(Tdomain%Globcoord(2,j)-zs)**2)
         if(d <= dmin)then
            dmin = d
            Tdomain%Ssource(n_src)%elem = el_around_node(i)
            Tdomain%Ssource(n_src)%gll(0) = x
            Tdomain%Ssource(n_src)%gll(1) = y
            Tdomain%Ssource(n_src)%gll(2) = z
         endif
      enddo
     enddo
    enddo
 enddo

 call MPI_ALLGATHER(dmin,1,MPI_DOUBLE_PRECISION,distance,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,code)
 do i = 0,nb_procs-1
    if(distance(i) <= dmin)then
       dmin = distance(i)
       Tdomain%sSource(n_src)%proc = i
    endif
 enddo

! from the closest GLL to coordinates in the source element (dichotomy: to be changed?)
    if(rg == Tdomain%sSource(n_src)%proc)then
        n = Tdomain%sSource(n_src)%elem
        mat = Tdomain%specel(n)%mat_index
        x = Tdomain%sSource(n_src)%gll(0)
        y = Tdomain%sSource(n_src)%gll(1)
        z = Tdomain%sSource(n_src)%gll(2)
        dist_xi = Tdomain%sSubdomain(mat)%GLLcx(x+1) - Tdomain%sSubdomain(mat)%GLLcx(x-1)
        dist_eta = Tdomain%sSubdomain(mat)%GLLcy(y+1) - Tdomain%sSubdomain(mat)%GLLcy(y-1)
        dist_zeta = Tdomain%sSubdomain(mat)%GLLcz(z+1) - Tdomain%sSubdomain(mat)%GLLcz(z-1)
        do i = 0,2
           xi_search(i) = Tdomain%sSubdomain(mat)%GLLcx(x-1) + (i+1)*dist_xi/4
           eta_search(i) = Tdomain%sSubdomain(mat)%GLLcy(y-1) + (i+1)*dist_eta/4
           zeta_search(i) = Tdomain%sSubdomain(mat)%GLLcz(z-1) + (i+1)*dist_zeta/4
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
                    if(Tdomain%n_nodes == 27)then
                       xa = 0;   ya = 0;   za = 0;
                       do i_count = 0,Tdomain%n_nodes-1
                          f = Comp_shapefunc(i_count,xi,eta,zeta)
                          xa = xa + coord(i_count,0)*f
                          ya = ya + coord(i_count,1)*f
                          za = za + coord(i_count,2)*f
                       enddo
                    else if(Tdomain%n_nodes == 8)then
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
           if(num > 20) exit dicho
        enddo dicho
        Tdomain%sSource(n_src)%RefCoord(:) = centre(:)

        ! InvGrad for the moment tensor source
        if(Tdomain%sSource(n_src)%i_type_source == 2)then
            xi = centre(0) ; eta = centre(1) ; zeta = centre(2)
            LocInvGrad(0,0) = 0.125*((coord(1,0)-coord(0,0))*(1-eta)*(1-zeta)+(coord(2,0)-coord(3,0))*(1+eta)*(1-zeta) + &
                (coord(5,0)-coord(4,0))*(1-eta)*(1+zeta) + (coord(6,0)-coord(7,0))*(1+eta)*(1+zeta))
            LocInvGrad(1,0) = 0.125*((coord(3,0)-coord(0,0))*(1-xi)*(1-zeta)+(coord(2,0)-coord(1,0))*(1+xi)*(1-zeta) + &
                (coord(7,0)-coord(4,0))*(1-xi)*(1+zeta)+(coord(6,0)-coord(5,0))*(1+xi)*(1+zeta))
            LocInvGrad(2,0) = 0.125*((coord(4,0)-coord(0,0))*(1-xi)*(1-eta) + (coord(5,0)-coord(1,0))*(1+xi)*(1-eta) + &
                (coord(7,0)-coord(3,0))*(1-xi)*(1+eta) + (coord(6,0)-coord(2,0))*(1+xi)*(1+eta))
            LocInvGrad(0,1) = 0.125*((coord(1,1)-coord(0,1))*(1-eta)*(1-zeta) + (coord(2,1)-coord(3,1))*(1+eta)*(1-zeta) + &
                (coord(5,1)-coord(4,1))*(1-eta)*(1+zeta) + (coord(6,1)-coord(7,1))*(1+eta)*(1+zeta))
            LocInvGrad(1,1) = 0.125*((coord(3,1)-coord(0,1))*(1-xi)*(1-zeta) + (coord(2,1)-coord(1,1))*(1+xi)*(1-zeta) + &
                (coord(7,1)-coord(4,1))*(1-xi)*(1+zeta) + (coord(6,1)-coord(5,1))*(1+xi)*(1+zeta))
            LocInvGrad(2,1) = 0.125*((coord(4,1)-coord(0,1))*(1-xi)*(1-eta) + (coord(5,1)-coord(1,1))*(1+xi)*(1-eta) + &
                (coord(7,1)-coord(3,1))*(1-xi)*(1+eta) + (coord(6,1)-coord(2,1))*(1+xi)*(1+eta))
            LocInvGrad(0,2) = 0.125*((coord(1,2)-coord(0,2))*(1-eta)*(1-zeta) + (coord(2,2)-coord(3,2))*(1+eta)*(1-zeta) + &
                (coord(5,2)-coord(4,2))*(1-eta)*(1+zeta) + (coord(6,2)-coord(7,2))*(1+eta)*(1+zeta))
            LocInvGrad(1,2) = 0.125*((coord(3,2)-coord(0,2))*(1-xi)*(1-zeta) + (coord(2,2)-coord(1,2))*(1+xi)*(1-zeta) + &
                (coord(7,2)-coord(4,2))*(1-xi)*(1+zeta) + (coord(6,2)-coord(5,2))*(1+xi)*(1+zeta))
            LocInvGrad(2,2) = 0.125*((coord(4,2)-coord(0,2))*(1-xi)*(1-eta) + (coord(5,2)-coord(1,2))*(1+xi)*(1-eta) + &
                (coord(7,2)-coord(3,2))*(1-xi)*(1+eta) + (coord(6,2)-coord(2,2))*(1+xi)*(1+eta))
       
           call invert_3d(LocInvGrad,trash)
           Tdomain%sSource(n_src)%InvGrad(0:2,0:2) = LocInvGrad(0:2,0:2)
           deallocate(coord)
        end if
    end if
enddo


return
end subroutine SourcePosition
