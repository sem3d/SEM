!>
!! \file Comm.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

module scomms

    ! Written by Paul Cupillard 07/11/05

    type :: comm

       integer :: nb_faces, nb_edges, nb_vertices, nb_edges_so,nb_vertices_so,nb_edges_neu,nb_vertices_neu, ngll, ngllPML, ngllSO
       integer, dimension(:), pointer :: faces, edges, vertices, edges_SO, vertices_SO,edges_Neu, vertices_Neu
       integer, dimension(:), pointer :: orient_faces, orient_edges, orient_edges_SO, orient_edges_Neu
       real, dimension(:), pointer :: Give, Take
       real, dimension(:,:), pointer :: GivePML, TakePML, GiveSO, TakeSO
       real, dimension(:,:), pointer :: GiveForces, TakeForces
       real, dimension(:,:,:), pointer :: GiveForcesPML, TakeForcesPML

    end type comm

end module scomms
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
