module scomms


    type :: comm

       integer :: nb_faces, nb_edges, nb_vertices, nb_edges_so,nb_vertices_so,    &
           nb_edges_neu,nb_vertices_neu,ngll,ngllPML,ngllSO,ngll_F,ngllPML_F,ngll_tot,ngllPML_tot
       integer, dimension(:), pointer :: faces, edges, vertices, edges_SO, vertices_SO,edges_Neu, vertices_Neu
       integer, dimension(:), pointer :: orient_faces, orient_edges, orient_edges_SO, orient_edges_Neu
       real, dimension(:), pointer :: Give, Take, GiveForcesSF_StoF, TakeForcesSF_StoF
       real, dimension(:,:), pointer :: GivePML, TakePML, GiveSO, TakeSO, GiveNeu, TakeNeu, GiveSF, TakeSF
       real, dimension(:,:), pointer :: GiveForces, TakeForces, GiveForcesSF_FtoS, TakeForcesSF_FtoS
       real, dimension(:,:,:), pointer :: GiveForcesPML, TakeForcesPML
       ! fluid communications
       real, dimension(:,:), pointer :: GiveForcesPMLFl, TakeForcesPMLFl
       real, dimension(:), pointer :: GiveForcesFl, TakeForcesFl

       ! Solid-fluid communication properties
       integer  :: SF_nf_shared, SF_ne_shared, SF_nv_shared, ngllSF
       integer, dimension(:), pointer  :: SF_faces_shared, SF_edges_shared,      &
           SF_vertices_shared,     &  ! indices
           SF_mapping_edges_shared
       ! orientations: only for edges. Face orientation is
       ! referenced using the classical number of the face.
       ! Neumann communication properties
       integer  :: Neu_ne_shared, Neu_nv_shared, ngllNeu
       integer, dimension(:), pointer  :: Neu_edges_shared, Neu_vertices_shared,     &  ! indices
           Neu_mapping_edges_shared     ! orientations
    end type comm

end module scomms
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
