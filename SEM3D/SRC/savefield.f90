subroutine savefield(Tdomain,it,rg,icount)

    use sdomain
    implicit none

    type(domain), intent(in):: Tdomain
    integer, intent(in) :: it,rg,icount
    integer :: i,n,nv,nbvert,nn
    integer, dimension(:), allocatable:: count
    character (len=100) :: fnamef
    real  :: x,y,z



    write(fnamef,"(a,I2.2,a,I5.5)") "Proc",rg,"field",icount
    open(61,file=fnamef,status="unknown",form="formatted")

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
        do n = 0, Tdomain%n_elem-1
            do i = 0,7
                nv = Tdomain%specel(n)%Near_Vertices(i)
                nn = Tdomain%sVertex(nv)%global_numbering
                x = Tdomain%Coord_Nodes(0,nn)
                y = Tdomain%Coord_Nodes(1,nn)
                z = Tdomain%Coord_Nodes(2,nn)
                if ( count(nv) < 0 .and. abs(y)<0.005d0) then
                    write (61,*) x,z, Tdomain%svertex(nv)%Veloc(0)**2+   &
                        Tdomain%svertex(nv)%Veloc(1)**2+Tdomain%svertex(nv)%Veloc(2)**2
                    count(nv) = 1
                endif
            enddo
        enddo
    endif

    close(61)

    deallocate(count)

    return

end subroutine savefield
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
