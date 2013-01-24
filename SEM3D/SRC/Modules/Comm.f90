!>
!! \file Comm.f90
!! \brief
!!
!<

module scomms

    type :: comm_vector_1d
       integer, dimension(:), allocatable :: n
       real, dimension(:,:), allocatable :: Give, Take
       integer, dimension(:), allocatable :: send_reqs
       integer, dimension(:), allocatable :: recv_reqs
    end type comm_vector_1d

    type :: comm_vector_2d
       integer :: n,m
       real, dimension(:,:), pointer :: Give, Take
    end type comm_vector_2d

    type :: comm

       integer :: nb_faces, nb_edges, nb_vertices, nb_edges_so,nb_vertices_so,    &
           nb_edges_neu,nb_vertices_neu,ngllSO,ngll_F,ngllPML_F,ngllPML_tot
       integer, dimension(:), pointer :: faces, edges, vertices, edges_SO, vertices_SO,edges_Neu, vertices_Neu
       integer, dimension(:), pointer :: orient_faces, orient_edges, orient_edges_SO, orient_edges_Neu
       ! Solide
       integer :: ngll_tot
       real, dimension(:), pointer :: Give, Take

       real, dimension(:,:), pointer :: GivePML, TakePML

       real, dimension(:), pointer :: GiveForcesSF_StoF, TakeForcesSF_StoF
       real, dimension(:,:), pointer :: GiveSO, TakeSO
       real, dimension(:,:), pointer :: GiveNeu, TakeNeu
       real, dimension(:,:), pointer :: GiveSF, TakeSF

       integer :: ngll ! *0:2
       real, dimension(:,:), pointer :: GiveForces, TakeForces
       integer :: ngllPML ! *1:3*0:2
       real, dimension(:,:,:), pointer :: GiveForcesPML, TakeForcesPML

       real, dimension(:,:), pointer :: GiveForcesSF_FtoS, TakeForcesSF_FtoS
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

contains

    subroutine allocate_comm_vector_1d(vector, nprocs, n)
        type(comm_vector_1d), intent(inout) :: vector
        integer, intent(in) :: nprocs, n

        allocate(vector%n(0:nprocs-1))
        allocate(vector%recv_reqs(0:nprocs-1))
        allocate(vector%send_reqs(0:nprocs-1))
        allocate(vector%Give(0:n-1, 0:nprocs-1))
        allocate(vector%Take(0:n-1, 0:nprocs-1))

    end subroutine allocate_comm_vector_1d
end module scomms
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
