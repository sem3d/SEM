!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file PML_def.f90
!!\brief La routine PML_definition applique les conditions initiales déclarées dans le fichier d'entrée.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<
subroutine PML_definition (Tdomain)

    use sdomain

    implicit none

    type (Domain), intent (INOUT) :: Tdomain

    integer :: n, nf, ne, nv, mat
    integer :: domtype

    do n = 0,Tdomain%n_face-1
        Tdomain%sFace(n)%Abs = .false.
    enddo
    do n = 0,Tdomain%n_edge-1
        Tdomain%sEdge(n)%Abs = .false.
    enddo
    do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%Abs = .false.
    enddo

    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        domtype = Tdomain%specel(n)%domain
        if (domtype/=DM_SOLID_PML .or. domtype/=DM_FLUID_PML) cycle

        if (Tdomain%sSubdomain(mat)%Px) then
            if (Tdomain%sSubdomain(mat)%Left) then
                nf = Tdomain%specel(n)%Near_Faces(4);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%Abs = .true.
            else
                nf = Tdomain%specel(n)%Near_Faces(2);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%Abs = .true.
            endif
        endif
        if (Tdomain%sSubdomain(mat)%Py) then
            if (Tdomain%sSubdomain(mat)%Forward) then
                nf = Tdomain%specel(n)%Near_Faces(1);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%Abs = .true.
            else
                nf = Tdomain%specel(n)%Near_Faces(3);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%Abs = .true.
            endif
        endif
        if (Tdomain%sSubdomain(mat)%Pz) then
            if (Tdomain%sSubdomain(mat)%Down) then
                nf = Tdomain%specel(n)%Near_Faces(0);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%Abs = .true.
            else
                nf = Tdomain%specel(n)%Near_Faces(5);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%Abs = .true.
            endif
        endif
    enddo

    return
end subroutine PML_definition

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
