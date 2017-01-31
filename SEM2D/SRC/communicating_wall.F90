!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file communicating_wall.F90
!!\brief
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module communication_object

    use constants

    type :: Communicating_wall

       logical, dimension (:), pointer :: Face_Coherency , FacePML_Coherency
       integer :: n_faces, n_vertices, n_points, n_pml_faces, n_points_pml, n_vertex_superobject
       integer, dimension (:), pointer  :: Face_list, Vertex_list, FacePML_List, Vertex_SuperObject_List
       real(fpp), dimension(:), pointer :: Send_data_1, Receive_data_1
       real(fpp), dimension(:,:),pointer:: Send_data_2, Receive_data_2

    end type Communicating_wall

end module communication_object

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
