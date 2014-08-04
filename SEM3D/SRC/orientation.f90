module orientation

contains
#if ! NEW_GLOBAL_METHOD   
    subroutine get_VectProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,   &
        prop_vertex,prop_elem)
        ! general routine for the assemblage procedure: Element -> vertex
        use mindex, only : ind_elem_vertex
        implicit none

        integer, intent(in) :: nv,ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(inout) :: prop_elem
        real, dimension(0:2), intent(in)  :: prop_vertex
        integer, dimension(0:2)  :: index_elem_v

        ! search for the relevant indices
        call ind_elem_vertex(nv,ngllx,nglly,ngllz,index_elem_v)
        ! assemblage
        prop_elem(index_elem_v(0),index_elem_v(1),index_elem_v(2),0:2) = prop_vertex(0:2)

        return

    end subroutine get_VectProperty_Vertex2Elem


    subroutine get_VectProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,   &
        prop_edge,prop_elem)
        ! general routine for the deassemblage procedure: Edge -> element
        use mindex, only : ind_elem_edge
        implicit none

        integer, intent(in) :: ne,orient_e,ngllx,nglly,ngllz,ngll
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(inout) :: prop_elem
        real, dimension(1:ngll-2,0:2), intent(in)  :: prop_edge
        integer, dimension(0:4)  :: index_elem_e

        ! search for the relevant indices
        call ind_elem_edge(ne,orient_e,ngllx,nglly,ngllz,index_elem_e)
        ! deassemblage
        select case(ne)
        case(1,3,8,11)  ! only y-coordinate does vary
            prop_elem(index_elem_e(0),index_elem_e(2):index_elem_e(3):index_elem_e(4),  &
                index_elem_e(1),0:2) = prop_edge(1:ngll-2,0:2)
        case(0,2,5,9)   ! only x-coordinate does vary
            prop_elem(index_elem_e(2):index_elem_e(3):index_elem_e(4),index_elem_e(0),  &
                index_elem_e(1),0:2) = prop_edge(1:ngll-2,0:2)
        case(4,6,7,10)   ! only z-coordinate does vary
            prop_elem(index_elem_e(0),index_elem_e(1),                 &
                index_elem_e(2):index_elem_e(3):index_elem_e(4),0:2) = prop_edge(1:ngll-2,0:2)
        end select

        return

    end subroutine get_VectProperty_Edge2Elem


    subroutine get_VectProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,   &
        prop_face,prop_elem)
        ! general routine for the deassemblage procedure: Face -> element
        use mindex, only : ind_elem_face
        implicit none

        integer, intent(in)  :: nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(inout) :: prop_elem
        real, dimension(1:ngll1-2,1:ngll2-2,0:2), intent(in) :: prop_face
        integer, dimension(0:6)  :: index_elem_f
        integer  :: i

        ! search for the relevant indices
        call ind_elem_face(nf,orient_f,ngllx,nglly,ngllz,index_elem_f)
        ! deassemblage
        select case(orient_f)
        case(0,1,2,3)
            if(nf == 2 .or. nf == 4)then
                prop_elem(index_elem_f(0),                                 &
                    index_elem_f(1):index_elem_f(2):index_elem_f(3),       &
                    index_elem_f(4):index_elem_f(5):index_elem_f(6),0:2) = &
                    prop_face(1:ngll1-2,1:ngll2-2,0:2)
            else if(nf == 1 .or. nf == 3)then
                prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3), &
                    index_elem_f(0),                                       &
                    index_elem_f(4):index_elem_f(5):index_elem_f(6),0:2) = &
                    prop_face(1:ngll1-2,1:ngll2-2,0:2)
            else if(nf == 0 .or. nf == 5)then
                prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3), &
                    index_elem_f(4):index_elem_f(5):index_elem_f(6),       &
                    index_elem_f(0),0:2) =                                 &
                    prop_face(1:ngll1-2,1:ngll2-2,0:2)
            end if
        case(4,5,6,7)
            if(nf == 2 .or. nf == 4)then
                do i = 0,2
                    prop_elem(index_elem_f(0),                                &
                        index_elem_f(1):index_elem_f(2):index_elem_f(3),      &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6),i) =  &
                        TRANSPOSE(prop_face(1:ngll1-2,1:ngll2-2,i))
                end do
            else if(nf == 1 .or. nf == 3)then
                do i = 0,2
                    prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),&
                        index_elem_f(0),                                      &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6),i) =  &
                        TRANSPOSE(prop_face(1:ngll1-2,1:ngll2-2,i))
                end do
            else if(nf == 0 .or. nf == 5)then
                do i = 0,2
                    prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),&
                        index_elem_f(4):index_elem_f(5):index_elem_f(6),      &
                        index_elem_f(0),i) =                                  &
                        TRANSPOSE(prop_face(1:ngll1-2,1:ngll2-2,i))
                end do
            end if
        end select

        return

    end subroutine get_VectProperty_Face2Elem


    subroutine get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,   &
        prop_edge,prop_elem)
        ! general routine for the assemblage procedure: Element -> edge
        use mindex, only : ind_elem_edge
        implicit none

        integer, intent(in) :: ne,orient_e,ngllx,nglly,ngllz,ngll
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(in) :: prop_elem
        real, dimension(1:ngll-2,0:2), intent(out)  :: prop_edge
        integer, dimension(0:4)  :: index_elem_e

        ! search for the relevant indices
        call ind_elem_edge(ne,orient_e,ngllx,nglly,ngllz,index_elem_e)
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

    subroutine get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,   &
        prop_face,prop_elem)
        ! general routine for the assemblage procedure: Element -> face
        use mindex, only : ind_elem_face
        implicit none

        integer, intent(in)  :: nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(in) :: prop_elem
        real, dimension(1:ngll1-2,1:ngll2-2,0:2), intent(out) :: prop_face
        integer, dimension(0:6)  :: index_elem_f
        integer  :: i


        ! search for the relevant indices
        call ind_elem_face(nf,orient_f,ngllx,nglly,ngllz,index_elem_f)
        ! assemblage
        select case(orient_f)
        case(0,1,2,3)
            if(nf == 2 .or. nf == 4)then
                prop_face(1:ngll1-2,1:ngll2-2,0:2) =                                       &
                    prop_face(1:ngll1-2,1:ngll2-2,0:2) +                             &
                    prop_elem(index_elem_f(0),                                   &
                    index_elem_f(1):index_elem_f(2):index_elem_f(3),   &
                    index_elem_f(4):index_elem_f(5):index_elem_f(6),0:2)
            else if(nf == 1 .or. nf == 3)then
                prop_face(1:ngll1-2,1:ngll2-2,0:2) =                                       &
                    prop_face(1:ngll1-2,1:ngll2-2,0:2) +                             &
                    prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),   &
                    index_elem_f(0),                                   &
                    index_elem_f(4):index_elem_f(5):index_elem_f(6),0:2)
            else if(nf == 0 .or. nf == 5)then
                prop_face(1:ngll1-2,1:ngll2-2,0:2) =                                       &
                    prop_face(1:ngll1-2,1:ngll2-2,0:2) +                             &
                    prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),   &
                    index_elem_f(4):index_elem_f(5):index_elem_f(6),   &
                    index_elem_f(0),0:2)
            end if
        case(4,5,6,7)
            if(nf == 2 .or. nf == 4)then
                do i = 0,2
                    prop_face(1:ngll1-2,1:ngll2-2,i) =                                       &
                        prop_face(1:ngll1-2,1:ngll2-2,i) +                             &
                        TRANSPOSE(prop_elem(index_elem_f(0),                                   &
                        index_elem_f(1):index_elem_f(2):index_elem_f(3),   &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6),i))
                end do
            else if(nf == 1 .or. nf == 3)then
                do i = 0,2
                    prop_face(1:ngll1-2,1:ngll2-2,i) =                                       &
                        prop_face(1:ngll1-2,1:ngll2-2,i) +                             &
                        TRANSPOSE(prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),   &
                        index_elem_f(0),                                   &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6),i))
                end do
            else if(nf == 0 .or. nf == 5)then
                do i = 0,2
                    prop_face(1:ngll1-2,1:ngll2-2,i) =                                       &
                        prop_face(1:ngll1-2,1:ngll2-2,i) +                             &
                        TRANSPOSE(prop_elem(index_elem_f(1):index_elem_f(2):index_elem_f(3),   &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6),   &
                        index_elem_f(0),i))
                end do
            end if
        end select

        return

    end subroutine get_VectProperty_Elem2face

    subroutine get_VectProperty_Elem2Vertex(nv,ngllx,nglly,ngllz,   &
        prop_vertex,prop_elem)
        ! general routine for the assemblage procedure: Element -> vertex
        use mindex, only : ind_elem_vertex
        implicit none

        integer, intent(in) :: nv,ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(in) :: prop_elem
        real, dimension(0:2), intent(out)  :: prop_vertex
        integer, dimension(0:2)  :: index_elem_v

        ! search for the relevant indices
        call ind_elem_vertex(nv,ngllx,nglly,ngllz,index_elem_v)
        ! assemblage
        prop_vertex(0:2) = prop_vertex(0:2)+   &
            prop_elem(index_elem_v(0),index_elem_v(1),index_elem_v(2),0:2)

        return

    end subroutine get_VectProperty_Elem2Vertex
#endif
end module orientation
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
