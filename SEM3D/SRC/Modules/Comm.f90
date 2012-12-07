module scomms


type :: comm

integer :: nb_faces, nb_edges, nb_vertices, nb_edges_so,nb_vertices_so,nb_edges_neu,nb_vertices_neu, ngll, ngllPML, ngllSO
integer, dimension(:), pointer :: faces, edges, vertices, edges_SO, vertices_SO,edges_Neu, vertices_Neu
integer, dimension(:), pointer :: orient_faces, orient_edges, orient_edges_SO, orient_edges_Neu
real, dimension(:), pointer :: Give, Take
real, dimension(:,:), pointer :: GivePML, TakePML, GiveSO, TakeSO
real, dimension(:,:), pointer :: GiveForces, TakeForces
real, dimension(:,:,:), pointer :: GiveForcesPML, TakeForcesPML

end type

end module scomms
