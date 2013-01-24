subroutine Comm_Vertex_ScalarProperty(exch_val,prop_vertex)
    ! from vector of exchanged interproc values, to vertices - scalar case
    implicit none

    real, intent(in)  :: exch_val
    real, intent(inout) :: prop_vertex

    prop_vertex = prop_vertex + exch_val

    return

end subroutine Comm_Vertex_ScalarProperty
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
