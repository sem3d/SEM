subroutine get_ScalarProperty_Elem2Vertex(nv,ngllx,nglly,ngllz,rank,   &
    prop_vertex,prop_elem)
    ! general routine for the assemblage procedure: Element -> vertex
    implicit none

    integer, intent(in) :: nv,ngllx,nglly,ngllz,rank
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: prop_elem
    real, intent(out)  :: prop_vertex
    integer, dimension(0:2)  :: index_elem_v

    ! search for the relevant indices
    call ind_elem_vertex(nv,ngllx,nglly,ngllz,index_elem_v,rank)
    ! assemblage
    prop_vertex = prop_vertex +   &
        prop_elem(index_elem_v(0),index_elem_v(1),index_elem_v(2))

    return

end subroutine get_ScalarProperty_Elem2Vertex
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
