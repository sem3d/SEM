!>
!! \file Comm.f90
!! \brief
!!
!<

module scomms

    type :: exchange_vector
       integer :: ndata ! Taille allouee de Give/Take
       integer :: nsend ! Taille a echanger pour une communication (<=ndata)
       integer :: nsol, nsolpml, nflu, nflupml ! Nombre de point de gauss a echanger pour chaque domaine
       ! si on echange 4 ddl pour solpml et 2 pour sol alors ndata=2*nsol+4*nsolpml...
       integer :: src, dest
       real, dimension(:), allocatable :: Give, Take
       integer, dimension(:), allocatable :: IGiveS, IGiveSPML, IGiveF, IGiveFPML
       integer, dimension(:), allocatable :: ITakeS, ITakeSPML, ITakeF, ITakeFPML
    end type

    type :: comm_vector
       integer :: ncomm
       type(exchange_vector), dimension(:), allocatable :: Data ! dimmension : 0,ncomm-1
       integer, dimension(:), allocatable :: send_reqs
       integer, dimension(:), allocatable :: recv_reqs
    end type comm_vector 

!     type :: comm_vector_1d
!        integer, dimension(:), allocatable :: n
!        real, dimension(:,:), allocatable :: Give, Take
!        integer, dimension(:), allocatable :: send_reqs
!        integer, dimension(:), allocatable :: recv_reqs
!     end type comm_vector_1d
! 
!     type :: comm_vector_2d
!        integer :: n,m
!        real, dimension(:,:), pointer :: Give, Take
!     end type comm_vector_2d

    type :: comm

       integer :: nb_faces, nb_edges, nb_vertices, nb_edges_so,nb_vertices_so,    &
           nb_edges_neu,nb_vertices_neu,ngllSO
#if ! NEW_GLOBAL_METHOD
       integer :: ngll_F,ngllPML_F,ngllPML_tot
#endif
       integer, dimension(:), pointer :: faces, edges, vertices, edges_SO, vertices_SO,edges_Neu, vertices_Neu
       integer, dimension(:), pointer :: orient_faces, orient_edges, orient_edges_SO, orient_edges_Neu
#if ! NEW_GLOBAL_METHOD
       ! Solide
       integer :: ngll_tot
       real, dimension(:), pointer :: Give, Take
       real, dimension(:,:), pointer :: GivePML, TakePML
#endif
       real, dimension(:), pointer :: GiveForcesSF_StoF, TakeForcesSF_StoF
       real, dimension(:,:), pointer :: GiveSO, TakeSO
       real, dimension(:,:), pointer :: GiveNeu, TakeNeu
       real, dimension(:,:), pointer :: GiveSF, TakeSF

#if ! NEW_GLOBAL_METHOD
       integer :: ngll ! *0:2
       real, dimension(:,:), pointer :: GiveForces, TakeForces
       integer :: ngllPML ! *1:3*0:2
#endif
       real, dimension(:,:,:), pointer :: GiveForcesPML, TakeForcesPML,   &
                                        GiveForcesSF_FtoS_PML,TakeForcesSF_FtoS_PML
       real, dimension(:,:), pointer :: GiveForcesSF_FtoS, TakeForcesSF_FtoS,  &
                                        GiveForcesSF_StoF_PML,TakeForcesSF_StoF_PML
#if ! NEW_GLOBAL_METHOD
       ! fluid communications
       real, dimension(:,:), pointer :: GiveForcesPMLFl, TakeForcesPMLFl
       real, dimension(:), pointer :: GiveForcesFl, TakeForcesFl
#endif
       ! Solid-fluid communication properties
       integer  :: SF_nf_shared, SF_ne_shared, SF_nv_shared, ngllSF, ngllSF_PML
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

!     subroutine allocate_comm_vector_1d(vector, nprocs, n)
!         type(comm_vector_1d), intent(inout) :: vector
!         integer, intent(in) :: nprocs, n
! 
!         allocate(vector%n(0:nprocs-1))
!         allocate(vector%recv_reqs(0:nprocs-1))
!         allocate(vector%send_reqs(0:nprocs-1))
!         allocate(vector%Give(0:n-1, 0:nprocs-1))
!         allocate(vector%Take(0:n-1, 0:nprocs-1))
! 
!     end subroutine allocate_comm_vector_1d
end module scomms
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
