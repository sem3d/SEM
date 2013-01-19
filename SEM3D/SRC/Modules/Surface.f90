!>
!! \file Surface.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

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

module sbassin
    implicit none

    type :: Bassin
       !       real     :: ymin,ymax
       !     n_layer nombre de couches
       !     n_colonne nombre de colonnes en x ici uniquement
       !     x_type == 0 on remet des materiaux  homogenes dans chaque bloc
       !     x_type == 1 on met des gradients pour chaque colonne en interpolant suivant z
       !     x_type == 2 on met des gradients pour chaque colonne en interpolant suivant z
       !      uniquement dans les couches entre zmin et zmax
       integer  :: n_colonne, n_layer, x_type
       real     :: zmin,zmax
       !    x_coord correspond aux abscisses des colonnes
       real, pointer, dimension(:) :: x_coord
       !      z_layer profondeur de  linterface pour chaque x de colonne
       !      on definit egalement le materiaux par rho, Cp , Cs
       real, pointer, dimension(:,:) :: z_layer, z_rho, z_Cp, z_Cs
    end type Bassin


end module sbassin
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
