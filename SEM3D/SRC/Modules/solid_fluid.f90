!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module solid_fluid
    use sinterface
    implicit none

    ! in general, when refering to both sides of faces, edges or vertices:
    !       index 0 for fluid part, 1 for the solid one

    ! Solid-fluid face properties
!    type :: face_SF
!        integer, dimension(0:1) :: Face
!        ! index 0 for fluid element, 1 for solid one
!        integer, dimension(0:3) :: Near_Edges, Near_Vertices
!        ! all 3 edges and vertices for a SF face
!        integer, dimension(0:3) :: Near_Edges_Orient
!        integer                 :: Orient_Face
!        integer                 :: ngll1,ngll2,dir
!        logical                 :: PML
!        real, allocatable, dimension(:,:,:) :: normal, BtN
!        ! Assigns a unique (assembled, shared with edge, vertices) index to solid-fluid glls
!        integer, allocatable, dimension(:,:) :: I_sf
!        integer, allocatable, dimension(:,:,:) :: IG
!    end type face_SF
!
!    type :: edge_SF
!        integer, dimension(0:1) :: Edge
!        ! index 0 for fluid element, 1 for solid one
!        integer                 :: Orient_Edge
!        integer                 :: ngll
!        logical                 :: PML
!        real, allocatable, dimension(:,:) ::  BtN,pn,save_forces,save_displ,pn1,pn2,pn3,  &
!            save_veloc1,save_veloc2,save_veloc3
!        real, allocatable, dimension(:) ::  vn,vn1,vn2,vn3
!        integer, allocatable, dimension(:) :: I_sf
!        integer, allocatable, dimension(:,:) :: IG
!    end type edge_SF
!
!    type :: vertex_SF
!        integer, dimension(0:1) :: vertex
!        ! index 0 for fluid element, 1 for solid one
!        logical                 :: PML
!        real, dimension(0:2) ::  BtN,pn,save_forces,save_displ,pn1,pn2,pn3,   &
!            save_veloc1,save_veloc2,save_veloc3
!        real  :: vn,vn1,vn2,vn3
!        integer :: I_sf
!        integer, dimension(0:2) :: IG
!    end type vertex_SF


    ! general SF object
    type :: SF_object
        type(inter_num) :: intSolFlu  ! 0: Solid 1: fluid
        type(inter_num) :: intSolFluPml ! 0: spml 1: fpml
        real, allocatable, dimension(:,:) :: SF_BtN
        real, allocatable, dimension(:,:) :: SFPml_BtN
        ! DELETE
!        integer  :: SF_n_faces, SF_n_edges, SF_n_vertices
!        integer  :: ngll, ngll_pml
!        integer, allocatable, dimension(:) :: SF_IGlobSol, SF_IGlobSol_pml
!        integer, allocatable, dimension(:) :: SF_IGlobFlu, SF_IGlobFlu_pml
!        ! SF_IGlob(idx,0) (resp. 1) contient pour le pt idx d'interf solid fluid l'index correspondant dans
!        ! le domaine fluide (resp. solide) -1 si le point n'existe pas sur le proc. Dans ce cas
!        ! SF_IGlob(idx,2) contient le numero du pt d'interface
!        integer, allocatable, dimension(:,:) :: SF_IGlob     ! 0 fluid, 1 solid, 2 index of shared gll
!        ! Idem SF_IGlob mais pour les domaines PML
!        integer, allocatable, dimension(:,:) :: SF_IGlob_pml ! 0 fluid, 1 solid, 2 index of shared gll
!        type(face_SF), dimension(:), pointer  :: SF_face
!        type(edge_SF), dimension(:), pointer  :: SF_edge
!        type(vertex_SF), dimension(:), pointer  :: SF_vertex
    end type SF_object

end module solid_fluid

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
