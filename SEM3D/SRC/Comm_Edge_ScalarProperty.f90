subroutine Comm_Edge_ScalarProperty(ngll1,orient_e,exch_val,prop_edge)
    ! from vector of exchanged interproc values, to edges - scalar case
    implicit none

    integer, intent(in)  :: ngll1,orient_e
    real, intent(in), dimension(ngll1-2)  :: exch_val
    real, dimension(1:ngll1-2), intent(inout) :: prop_edge
    integer  :: j,ngll

    ngll = 1
    select case(orient_e)
    case(0)
        do j = 1,ngll1-2
            prop_edge(j) = prop_edge(j) + exch_val(ngll)
            ngll = ngll + 1
        end do
    case(1)
        do j = 1,ngll1-2
            prop_edge(ngll1-1-j) = prop_edge(ngll1-1-j) + exch_val(ngll)
            ngll = ngll + 1
        end do
    end select

    if(ngll /= ngll1-1) &
        stop "Pb in Comm_Edge_ScalarProperty"

    return

end subroutine Comm_Edge_ScalarProperty
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
