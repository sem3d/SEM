subroutine define_Neumann_properties(Tdomain,rank)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)  :: rank
    integer :: ngll,ngll1,ngll2,mat_index,nf,ne,nv,nf_aus

    ! allocations for faces, edges and vertices
    do nf = 0,Tdomain%Neumann%Neu_n_faces-1
        ngll1 = Tdomain%Neumann%Neu_Face(nf)%ngll1
        ngll2 = Tdomain%Neumann%Neu_Face(nf)%ngll2
        nf_aus = Tdomain%Neumann%Neu_Face(nf)%Face
        !if(Tdomain%Neumann%Neu_Param%what_bc == "F")then
        !    Tdomain%Neumann%Neu_Face(nf)%mat_index =       &
        !        Tdomain%specel(Tdomain%sFace(nf_aus)%Which_Elem)%mat_index
        !else
        !    Tdomain%Neumann%Neu_Face(nf)%mat_index = Tdomain%Neumann%Neu_Param%mat_index
        !end if
        allocate(Tdomain%Neumann%Neu_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,0:2))
        allocate(Tdomain%Neumann%Neu_Face(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2))
        Tdomain%Neumann%Neu_Face(nf)%BtN = 0d0
    enddo

    do ne = 0,Tdomain%Neumann%Neu_n_edges-1
        ngll = Tdomain%Neumann%Neu_Edge(ne)%ngll
        allocate(Tdomain%Neumann%Neu_Edge(ne)%BtN(1:ngll-2,0:2))
        allocate(Tdomain%Neumann%Neu_Edge(ne)%Forces(1:ngll-2,0:2))
        Tdomain%Neumann%Neu_Edge(ne)%BtN = 0d0
    enddo

    do nv = 0,Tdomain%Neumann%Neu_n_vertices-1
        Tdomain%Neumann%Neu_Vertex(nv)%BtN = 0d0
    enddo

    ! elastic properties
    Tdomain%Neumann%Neu_Param%lambda = Tdomain%sSubDomain(Tdomain%Neumann%Neu_Param%mat_index)%DLambda
    Tdomain%Neumann%Neu_Param%Mu = Tdomain%sSubDomain(Tdomain%Neumann%Neu_Param%mat_index)%DMu

    if(Tdomain%Neumann%Neu_Param%wtype == 'S') then
        Tdomain%Neumann%Neu_Param%speed = Tdomain%sSubDomain(Tdomain%Neumann%Neu_Param%mat_index)%Sspeed
    else
        Tdomain%Neumann%Neu_Param%speed = Tdomain%sSubDomain(Tdomain%Neumann%Neu_Param%mat_index)%Pspeed
    endif


    return
end subroutine define_Neumann_properties
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
