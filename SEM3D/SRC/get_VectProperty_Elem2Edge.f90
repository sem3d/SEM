subroutine get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,   &
    prop_edge,prop_elem)
    ! general routine for the assemblage procedure: Element -> edge
    implicit none

    integer, intent(in) :: ne,orient_e,ngllx,nglly,ngllz,ngll,rank
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(in) :: prop_elem
    real, dimension(1:ngll-2,0:2), intent(out)  :: prop_edge
    integer, dimension(0:4)  :: index_elem_e

    ! search for the relevant indices
    call ind_elem_edge(ne,orient_e,ngllx,nglly,ngllz,index_elem_e,rank)
    ! assemblage
    select case(ne)
    case(1,3,8,11)  ! only y-coordinate does vary
        prop_edge(1:ngll-2,0:2) = prop_edge(1:ngll-2,0:2) +   &
            prop_elem(index_elem_e(0),index_elem_e(2):index_elem_e(3):index_elem_e(4),  &
            index_elem_e(1),0:2)
    case(0,2,5,9)   ! only x-coordinate does vary
        prop_edge(1:ngll-2,0:2) = prop_edge(1:ngll-2,0:2) +   &
            prop_elem(index_elem_e(2):index_elem_e(3):index_elem_e(4),index_elem_e(0),  &
            index_elem_e(1),0:2)
    case(4,6,7,10)   ! only z-coordinate does vary
        prop_edge(1:ngll-2,0:2) = prop_edge(1:ngll-2,0:2) +   &
            prop_elem(index_elem_e(0),index_elem_e(1),                 &
            index_elem_e(2):index_elem_e(3):index_elem_e(4),0:2)
    end select

    return

end subroutine get_VectProperty_Elem2Edge
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
