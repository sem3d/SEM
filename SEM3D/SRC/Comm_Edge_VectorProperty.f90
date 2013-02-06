subroutine Comm_Edge_VectorProperty(ngll1,orient_e,exch_val,prop_edge)
    ! from vector of exchanged interproc values, to edges - vectorial case
    implicit none

    integer, intent(in)  :: ngll1,orient_e
    real, intent(in), dimension(ngll1-2,0:2)  :: exch_val
    real, dimension(1:ngll1-2,0:2), intent(inout) :: prop_edge
    integer  :: j

    select case(orient_e)
    case(0)
        do j = 1,ngll1-2
            prop_edge(j,0:2) = prop_edge(j,0:2) + exch_val(j,0:2)
        end do
    case(1)
        do j = 1,ngll1-2
            prop_edge(ngll1-1-j,0:2) = prop_edge(ngll1-1-j,0:2) + exch_val(j,0:2)
        end do
    end select

end subroutine Comm_Edge_VectorProperty

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
