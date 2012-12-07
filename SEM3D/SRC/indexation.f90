subroutine ind_elem_face(ind_face,orient_face,ngllx,nglly,ngllz,index_elem)
  !- routine giving info about index array for Iglobnum, depending on the index of
  !     the face and its orientation:
  !- "index_elem" contains the whole set of extreme indices and increment, for each
  !  space dimension that vary: index_elem(1) the min index for x-dimension in Iglobnum,
  !   index_elem(2) the max index and index(3) the increment or decrement. And the same
  !   for other varying space dimensions. For the one that does stay fixed: the first slot
  !   index_elem(0), which gives the GLL value for the fixed coordinate (0 or ngll-1).
    implicit none
    integer, intent(in)  :: ind_face,orient_face,ngllx,nglly,ngllz
    integer, dimension(0:6), intent(out)   :: index_elem

    select case(ind_face)
        case(0)   ! constant z-coordinate : ngll_z=0
            index_elem(0) = 0 
          ! only x and y related indices do vary:
            call index_face(ngllx,nglly,orient_face,index_elem(1),index_elem(2),   &
                            index_elem(3),index_elem(4),index_elem(5),index_elem(6))
        case(1)   ! constant y-coordinate : ngll_y=0
            index_elem(0) = 0 
          ! only x and z related indices do vary:
            call index_face(ngllx,ngllz,orient_face,index_elem(1),index_elem(2),   &
                            index_elem(3),index_elem(4),index_elem(5),index_elem(6))
        case(2)   ! constant x-coordinate : ngll_x=ngllx-1
            index_elem(0) = ngllx-1 
          ! only y and z related indices do vary:
            call index_face(nglly,ngllz,orient_face,index_elem(1),index_elem(2),   &
                            index_elem(3),index_elem(4),index_elem(5),index_elem(6))
        case(3)   ! constant y-coordinate : ngll_y=nglly-1
            index_elem(0) = nglly-1 
          ! only x and z related indices do vary:
            call index_face(ngllx,ngllz,orient_face,index_elem(1),index_elem(2),   &
                            index_elem(3),index_elem(4),index_elem(5),index_elem(6))
        case(4)   ! constant x-coordinate : ngll_x=0
            index_elem(0) = 0 
          ! only y and z related indices do vary:
            call index_face(nglly,ngllz,orient_face,index_elem(1),index_elem(2),   &
                            index_elem(3),index_elem(4),index_elem(5),index_elem(6))
        case(5)   ! constant z-coordinate : ngll_z=ngllz-1
            index_elem(0) = ngllz-1 
          ! only x and y related indices do vary:
            call index_face(ngllx,nglly,orient_face,index_elem(1),index_elem(2),   &
                            index_elem(3),index_elem(4),index_elem(5),index_elem(6))

    end select


end subroutine ind_elem_face
!------------------------------------------------------------
!------------------------------------------------------------
subroutine index_face(ngll1,ngll2,orient_f,ind1,ind2,ind3,ind4,ind5,ind6)
  ! returns min and max bounds for indices, and the de- or in-crement, which
  !   depend on the orientation of the face.
    implicit none
    integer, intent(in)  :: ngll1,ngll2,orient_f
    integer, intent(out) :: ind1,ind2,ind3,ind4,ind5,ind6

    select case(orient_f)
        case(0)
            ind1 = 1 ; ind2 = ngll1-2 ; ind3 = 1
            ind4 = 1 ; ind5 = ngll2-2 ; ind6 = 1
        case(1)
            ind1 = ngll1-2 ; ind2 = 1 ; ind3 = -1
            ind4 = 1 ; ind5 = ngll2-2 ; ind6 = 1
        case(2)
            ind1 = 1 ; ind2 = ngll1-2 ; ind3 = 1
            ind4 = ngll2-2 ; ind5 = 1 ; ind6 = -1
        case(3)
            ind1 = ngll1-2 ; ind2 = 1 ; ind3 = -1
            ind4 = ngll2-2 ; ind5 = 1 ; ind6 = -1
        case(4)
            ind1 = 1 ; ind2 = ngll2-2 ; ind3 = 1
            ind4 = 1 ; ind5 = ngll1-2 ; ind6 = 1
        case(5)
            ind1 = ngll2-2 ; ind2 = 1 ; ind3 = -1
            ind4 = 1 ; ind5 = ngll1-2 ; ind6 = 1
        case(6)
            ind1 = 1 ; ind2 = ngll2-2 ; ind3 = 1
            ind4 = ngll1-2 ; ind5 = 1 ; ind6 = -1
        case(7)
            ind1 = ngll2-2 ; ind2 = 1 ; ind3 = -1
            ind4 = ngll1-2 ; ind5 = 1 ; ind6 = -1
    end select
   
end subroutine index_face
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine ind_elem_edge(ind_edge,orient_edge,ngllx,nglly,ngllz,index_elem,rank)
  !- routine giving info about index array for Iglobnum, depending on the index of
  !     the edge and its orientation:
  !- "index_elem" contains the whole set of extreme indices and increment, for the
  !  only space dimension that vary: index_elem(2) the min index for this dimension in Iglobnum,
  !   index_elem(3) the max index and index(4) the increment or decrement. For the ones that
  !    do stay fixed: the first two slots index_elem(0) and index_elem(1), which give the GLL
  !    values for the fixed coordinates (0 or ngll-1).
    implicit none
    integer, intent(in)  :: ind_edge,orient_edge,ngllx,nglly,ngllz,rank
    integer, dimension(0:4), intent(out)   :: index_elem

    select case(ind_edge)
        case(0) ! constant y (ngll_y = 0) and z (ngll_z = 0)
            index_elem(0) = 0
            index_elem(1) = 0
          ! only x related index does vary:
            call index_edge(ngllx,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(1) ! constant x (ngll_x = ngll-1) and z (ngll_z = 0)
            index_elem(0) = ngllx-1
            index_elem(1) = 0
          ! only y related index does vary:
            call index_edge(nglly,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(2) ! constant y (ngll_y = ngll-1) and z (ngll_z = 0)
            index_elem(0) = nglly-1
            index_elem(1) = 0
          ! only x related index does vary:
            call index_edge(ngllx,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(3) ! constant x (ngll_x = 0) and z (ngll_z = 0)
            index_elem(0) = 0
            index_elem(1) = 0
          ! only y related index does vary:
            call index_edge(nglly,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(4) ! constant x (ngll_x = ngll-1) and y (ngll_y = 0)
            index_elem(0) = ngllx-1
            index_elem(1) = 0
          ! only z related index does vary:
            call index_edge(ngllz,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(5) ! constant y (ngll_y = 0) and z (ngll_z = ngll_z)
            index_elem(0) = 0
            index_elem(1) = ngllz-1
          ! only x related index does vary:
            call index_edge(ngllx,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(6) ! constant x (ngll_x = 0) and y (ngll_y = 0)
            index_elem(0) = 0
            index_elem(1) = 0
          ! only z related index does vary:
            call index_edge(ngllz,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(7) ! constant x (ngll_x = ngll-1) and y (ngll_y = ngll-1)
            index_elem(0) = ngllx-1
            index_elem(1) = nglly-1
          ! only z related index does vary:
            call index_edge(ngllz,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(8) ! constant x (ngll_x = ngll-1) and z (ngll_z = ngll-1)
            index_elem(0) = ngllx-1
            index_elem(1) = ngllz-1
          ! only y related index does vary:
            call index_edge(nglly,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(9) ! constant y (ngll_y = ngll-1) and z (ngll_z = ngll-1)
            index_elem(0) = nglly-1
            index_elem(1) = ngllz-1
          ! only x related index does vary:
            call index_edge(ngllx,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(10) ! constant x (ngll_x = 0) and y (ngll_y = ngll-1)
            index_elem(0) = 0
            index_elem(1) = nglly-1
          ! only z related index does vary:
            call index_edge(ngllz,orient_edge,index_elem(2),index_elem(3),index_elem(4))
        case(11) ! constant x (ngll_x = 0) and z (ngll_z = ngll-1)
            index_elem(0) = 0
            index_elem(1) = ngllz-1
          ! only y related index does vary:
            call index_edge(nglly,orient_edge,index_elem(2),index_elem(3),index_elem(4))
    end select

end subroutine ind_elem_edge
!------------------------------------------------------------
!------------------------------------------------------------
subroutine index_edge(ngll,orient_e,ind1,ind2,ind3)
  ! returns min and max bounds for indices, and the de- or in-crement, which
  !   depend on the orientation of the edge.
    implicit none
    integer, intent(in)  :: ngll,orient_e
    integer, intent(out) :: ind1,ind2,ind3

    select case(orient_e)
        case(0)
            ind1 = 1 ; ind2 = ngll-2 ; ind3 = 1
        case(1)
            ind1 = ngll-2 ; ind2 = 1 ; ind3 = -1
    end select

end subroutine index_edge
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine ind_elem_vertex(ind_vertex,ngllx,nglly,ngllz,index_elem,rank)
    implicit none
    integer, intent(in)  :: ind_vertex,ngllx,nglly,ngllz,rank
    integer, dimension(0:2), intent(out)  :: index_elem

    select case(ind_vertex)
        case(0)
            index_elem(0) = 0 ; index_elem(1) = 0 ; index_elem(2) = 0
        case(1)
            index_elem(0) = ngllx-1 ; index_elem(1) = 0 ; index_elem(2) = 0
        case(2)
            index_elem(0) = ngllx-1 ; index_elem(1) = nglly-1 ; index_elem(2) = 0
        case(3)
            index_elem(0) = 0 ; index_elem(1) = nglly-1 ; index_elem(2) = 0
        case(4)
            index_elem(0) = 0 ; index_elem(1) = 0 ; index_elem(2) = ngllz-1
        case(5)
            index_elem(0) = ngllx-1 ; index_elem(1) = 0 ; index_elem(2) = ngllz-1
        case(6)
            index_elem(0) = ngllx-1 ; index_elem(1) = nglly-1 ; index_elem(2) = ngllz-1
        case(7)
            index_elem(0) = 0 ; index_elem(1) = nglly-1 ; index_elem(2) = ngllz-1
    end select

end subroutine ind_elem_vertex
!---------------------------------------------------------------------------
!-------------------------------------------------------------------------
