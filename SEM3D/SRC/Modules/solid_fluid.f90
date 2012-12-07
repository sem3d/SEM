module solid_fluid

    implicit none

! in general, when refering to both sides of faces, edges or vertices:
!       index 0 for fluid part, 1 for the solid one

! Solid-fluid face properties
    type :: face_SF
        integer, dimension(0:1) :: Face 
                     ! index 0 for fluid element, 1 for solid one
        integer, dimension(0:3) :: Near_Edges, Near_Vertices
                     ! all 3 edges and vertices for a SF face
        integer, dimension(0:3) :: Near_Edges_Orient
        integer                 :: Orient_Face
        integer                 :: ngll1,ngll2,dir
        real, allocatable, dimension(:,:,:) :: normal, BtN
    end type face_SF
    
    type :: edge_SF
        integer, dimension(0:1) :: Edge
                     ! index 0 for fluid element, 1 for solid one
        integer                 :: Orient_Edge
        integer                 :: ngll
        real, allocatable, dimension(:,:) ::  BtN
    end type edge_SF

    type :: vertex_SF
        integer, dimension(0:1) :: vertex
                     ! index 0 for fluid element, 1 for solid one
        real, dimension(0:2) ::  BtN
    end type vertex_SF


! general SF object
    type :: SF_object
        integer  :: SF_n_faces, SF_n_edges, SF_n_vertices
        type(face_SF), dimension(:), pointer  :: SF_face
        type(edge_SF), dimension(:), pointer  :: SF_edge
        type(vertex_SF), dimension(:), pointer  :: SF_vertex
    end type SF_object

end module solid_fluid
