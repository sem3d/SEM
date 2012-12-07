subroutine get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,   &
                                  rank,prop_face,prop_elem)
  ! general routine for the assemblage procedure: Element -> face
    implicit none

    integer, intent(in)  :: nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(in) :: prop_elem 
    real, dimension(1:ngll1-2,1:ngll2-2,0:2), intent(out) :: prop_face
    integer, dimension(0:6)  :: index_elem_f
    integer  :: i
    

  ! search for the relevant indices
    call ind_elem_face(nf,orient_f,ngllx,nglly,ngllz,index_elem_f,rank)
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
