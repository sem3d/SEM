subroutine savefield (Tdomain,it,rg,icount)

use sdomain

implicit none

type (domain), intent (IN):: Tdomain
integer, intent (IN) :: it,rg,icount

! local variables
integer :: i,n,nv,nbvert
integer, dimension(:), allocatable:: count
character (len=100) :: fnamef



  write (fnamef,"(a,I2.2,a,I3.3)") "SField/Proc",rg,"fieldx",icount
  open (61, file=fnamef, status="unknown", form="formatted")
  write (fnamef,"(a,I2.2,a,I3.3)") "SField/Proc",rg,"fieldy",icount
  open (62, file=fnamef, status="unknown", form="formatted")
  write (fnamef,"(a,I2.2,a,I3.3)") "SField/Proc",rg,"fieldz",icount
  open (63, file=fnamef, status="unknown", form="formatted")


  if ( Tdomain%logicD%Save_Surface .and. Tdomain%logicD%Neumann ) then
     nbvert = Tdomain%sSurf%n_vertices+Tdomain%sNeu%n_vertices
  else if (Tdomain%logicD%Save_Surface) then
     nbvert = Tdomain%sSurf%n_vertices
  else ! Save all the vertices
    allocate (count(0:Tdomain%n_vertex-1))
    count = -1
    nbvert = 0
    do n = 0, Tdomain%n_elem - 1
     do i = 0,7
       nv = Tdomain%specel(n)%Near_Vertices(i)
       if ( count(nv) < 0 ) then
         nbvert = nbvert+1
         count(nv) = 1
       endif
     enddo
    enddo
  endif      

  write (61,*) nbvert,Tdomain%TimeD%time_snapshots,Tdomain%TimeD%rtime
  write (62,*) nbvert,Tdomain%TimeD%time_snapshots,Tdomain%TimeD%rtime
  write (63,*) nbvert,Tdomain%TimeD%time_snapshots,Tdomain%TimeD%rtime

  if (Tdomain%logicD%Save_Surface) then
    do nv = 0,Tdomain%sSurf%n_vertices-1
       n = Tdomain%sSurf%nVertex(nv)%Vertex
       write (61,*) n+1, Tdomain%svertex(n)%Veloc(0)
       write (62,*) n+1, Tdomain%svertex(n)%Veloc(1)
       write (63,*) n+1, Tdomain%svertex(n)%Veloc(2)
    enddo
    if (Tdomain%logicD%Neumann) then
      do nv = 0,Tdomain%sNeu%n_vertices-1
         n = Tdomain%sNeu%nVertex(nv)%Vertex
         write (61,*) n+1, Tdomain%svertex(n)%Veloc(0)
         write (62,*) n+1, Tdomain%svertex(n)%Veloc(1)
         write (63,*) n+1, Tdomain%svertex(n)%Veloc(2)
      enddo
    endif
  else
    count = -1
    do n = 0, Tdomain%n_elem - 1
     do i = 0,7
       nv = Tdomain%specel(n)%Near_Vertices(i)
       if ( count(nv) < 0 ) then
         write (61,*) nv+1, Tdomain%svertex(nv)%Veloc(0)
         write (62,*) nv+1, Tdomain%svertex(nv)%Veloc(1)
         write (63,*) nv+1, Tdomain%svertex(nv)%Veloc(2)
         count(nv) = 1
       endif
     enddo
    enddo
  endif

  close(61)
  close(62)
  close(63)

deallocate(count)


return

end subroutine savefield
