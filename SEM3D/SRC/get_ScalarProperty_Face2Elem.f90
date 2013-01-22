subroutine get_ScalarProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,   &
    prop_face,prop_elem)
    ! general routine for the deassemblage procedure: Face -> element
    use mindex, only : ind_elem_face
    implicit none

    integer, intent(in)  :: nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(inout) :: prop_elem
    real, dimension(1:ngll1-2,1:ngll2-2), intent(in) :: prop_face
    integer, dimension(0:6)  :: index_elem_f
    integer  :: i


    ! search for the relevant indices
    call ind_elem_face(nf,orient_f,ngllx,nglly,ngllz,index_elem_f)

    ! deassemblage
    select case(orient_f)
    case(0,1,2,3)
        if(nf == 2 .or. nf == 4)then
            prop_elem(index_elem_f(0),                                       &
                index_elem_f(1):index_elem_f(2):index_elem_f(3),       &
                index_elem_f(4):index_elem_f(5):index_elem_f(6)) =     &
                prop_face(1:ngll1-2,1:ngll2-2)
        else if(nf == 1 .or. nf == 3)then
            prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),       &
                index_elem_f(0),                                       &
                index_elem_f(4):index_elem_f(5):index_elem_f(6)) =     &
                prop_face(1:ngll1-2,1:ngll2-2)
        else if(nf == 0 .or. nf == 5)then
            prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),       &
                index_elem_f(4):index_elem_f(5):index_elem_f(6),       &
                index_elem_f(0)) =                                     &
                prop_face(1:ngll1-2,1:ngll2-2)
        end if
    case(4,5,6,7)
        if(nf == 2 .or. nf == 4)then
            prop_elem(index_elem_f(0),                                      &
                index_elem_f(1):index_elem_f(2):index_elem_f(3),      &
                index_elem_f(4):index_elem_f(5):index_elem_f(6)) =    &
                TRANSPOSE(prop_face(1:ngll1-2,1:ngll2-2))
        else if(nf == 1 .or. nf == 3)then
            prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),      &
                index_elem_f(0),                                      &
                index_elem_f(4):index_elem_f(5):index_elem_f(6)) =    &
                TRANSPOSE(prop_face(1:ngll1-2,1:ngll2-2))
        else if(nf == 0 .or. nf == 5)then
            prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),      &
                index_elem_f(4):index_elem_f(5):index_elem_f(6),      &
                index_elem_f(0)) =                                    &
                TRANSPOSE(prop_face(1:ngll1-2,1:ngll2-2))
        end if
    end select

    return

end subroutine get_ScalarProperty_Face2Elem
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
