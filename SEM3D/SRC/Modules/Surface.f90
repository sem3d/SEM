module ssurf

    ! #####################################################################################
    ! #####################################################################################

    type Face_Surf

       integer :: ngll1, ngll2, mat_index, Face
       integer, dimension (0:3) :: Near_Edges, Near_Vertices

    end type Face_Surf


    type Edge_Surf

       integer :: ngll, mat_index, Edge

    end type Edge_Surf


    type Vertex_Surf

       integer :: Vertex , mat_index

    end type Vertex_Surf

    type Surf

       integer :: n_faces, n_edges, n_vertices
       type(Face_Surf), dimension (:), pointer :: nFace
       type(Edge_Surf), dimension (:), pointer :: nEdge
       type(Vertex_Surf), dimension (:), pointer :: nVertex

    end type Surf


end module ssurf
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
