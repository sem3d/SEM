subroutine SourcePosition(Tdomain,nb_procs,rg)

use sdomain
use pig

implicit none

include 'mpif.h'

type (domain), intent(inout) :: Tdomain
integer, intent(in) :: nb_procs,rg

integer :: n_src, i,j, n_around_elem, x,y,z, ngllx,nglly,ngllz, code
integer , dimension (0:20) :: el_around_node
integer , dimension (:), allocatable :: near_node
real :: xs,ys,zs, d,dmin, R,Xtemp,Ytemp
real, dimension(0:nb_procs-1) :: distance


allocate (near_node(0:Tdomain%n_source-1))

do n_src = 0, Tdomain%n_source-1
 if (Tdomain%curve) then
    R = Tdomain%Ssource(n_src)%Zsource
    Xtemp = tan(pi*Tdomain%Ssource(n_src)%Xsource/180)
    Ytemp = tan(pi*Tdomain%Ssource(n_src)%Ysource/180)
    D = sqrt(1 + Ytemp**2 + Xtemp**2)
    xs = R*Xtemp/D
    ys = R*Ytemp/D
    zs = R/D
 else
    xs = Tdomain%Ssource(n_src)%Xsource
    ys = Tdomain%Ssource(n_src)%Ysource
    zs = Tdomain%Ssource(n_src)%Zsource
 endif
 dmin = 1000000
 do i = 0,Tdomain%n_glob_nodes-1
    d = sqrt ((Tdomain%Coord_nodes(0,i)-xs)**2 + (Tdomain%Coord_nodes(1,i)-ys)**2 + (Tdomain%Coord_nodes(2,i)-zs)**2)
    if (d <= dmin) then
       dmin = d
       near_node(n_src) = i
    endif
 enddo
 n_around_elem = 0
 do i = 0,Tdomain%n_elem-1
    do j = 0,Tdomain%n_nodes-1
       if (Tdomain%specel(i)%Control_nodes(j)==near_node(n_src)) then
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
 enddo

 call mpi_allgather(dmin,1,mpi_double_precision,distance,1,mpi_double_precision,mpi_comm_world,code)
 do i = 0,nb_procs-1
    if (distance(i) <= dmin) then
       dmin = distance(i)
       Tdomain%sSource(n_src)%proc = i
    endif
 enddo

if (rg==Tdomain%sSource(n_src)%proc) then
print*,'Source',n_src,Tdomain%Globcoord(0,Tdomain%specel(Tdomain%sSource(n_src)%elem)%Iglobnum(Tdomain%sSource(n_src)%gll(0), &
                                                     Tdomain%sSource(n_src)%gll(1),Tdomain%sSource(n_src)%gll(2))),&
             Tdomain%Globcoord(1,Tdomain%specel(Tdomain%sSource(n_src)%elem)%Iglobnum(Tdomain%sSource(n_src)%gll(0), &
                                                    Tdomain%sSource(n_src)%gll(1),Tdomain%sSource(n_src)%gll(2))),&
             Tdomain%Globcoord(2,Tdomain%specel(Tdomain%sSource(n_src)%elem)%Iglobnum(Tdomain%sSource(n_src)%gll(0), &
                                                     Tdomain%sSource(n_src)%gll(1),Tdomain%sSource(n_src)%gll(2)))
endif
enddo

deallocate (near_node)

return
end subroutine SourcePosition
