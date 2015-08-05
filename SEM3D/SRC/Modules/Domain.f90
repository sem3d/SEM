!>
!!\file Domain.f90
!!\brief Contient le d�finition du type domain
!!
!<

module sdomain

    use selement
    use sfaces
    use sedges
    use svertices
    use scomms
    use ssources
    use stimeparam
    use logical_input
    use ssubdomains
    use splanew
    use sneu
    use ssurf
    use sbassin
    use solid_fluid
    use semdatafiles
    use sem_c_config

    type :: domain
       integer :: communicateur !<<< Communicator including all SEM processors
       integer :: rank          !<<< Rank of this process within this communicator
       integer :: nb_procs      !<<< Total number of SEM processors
       ! Without coupling : communicateur=communicateur_global
       ! With coupling    : communicateur : includes every processes
       integer :: communicateur_global
       ! Communicator used for output grouping. Only rank 0 of this comm produces outputs
       integer :: comm_output
       ! Nombre de processeur dans le groupe de communication associe a comm_output
       integer :: nb_output_procs
       integer :: output_rank
       integer, dimension(:), allocatable :: output_nodes, output_nodes_offset, output_elems
       ! Nombre de process par sorties pour le reassemblage
       integer :: ngroup
       ! Nombre de processeur avec qui on communique (size(sComm))
       integer :: tot_comm_proc
       ! En mode couplage : Rg du superviseur dans le communicateur global
       integer :: master_superviseur

       type(time)          :: TimeD
       type(logical_array) :: logicD
       type(planew)        :: sPlaneW
       type(Neu_object)    :: Neumann
       type(surf)          :: sSurf
       type(bassin)        :: sBassin
       type(SF_object)     :: SF
       type(source)   , dimension (:), pointer :: sSource
       type(element)  , dimension (:), pointer :: specel
       type(face)     , dimension (:), pointer :: sFace
       type(comm)     , dimension (:), pointer :: sComm
       type(edge)     , dimension (:), pointer :: sEdge
       type(vertex)   , dimension (:), pointer :: sVertex
       type(subdomain), dimension (:), pointer :: sSubDomain


       logical :: any_PML, curve, any_FPML, aniso, any_Random

       integer :: n_source, n_dime, n_glob_nodes, n_mat, n_nodes, n_receivers
       integer :: n_elem, n_face, n_edge, n_vertex, n_glob_points, n_sls
       integer :: n_hexa  !< Nombre de maille hexa ~= (ngllx-1)*(nglly-1)*(ngllz-1)*nelem
       logical, dimension(:), allocatable :: not_PML_List, subD_exist
       integer, dimension(:), allocatable :: subDComm

       real :: T1_att, T2_att, T0_modele
       real, dimension (0:2,0:2) :: rot
       real, dimension (:,:), pointer :: Coord_nodes, GlobCoord

       integer :: traces_format
       character (len=MAX_FILE_SIZE) :: Title_simulation, mesh_file,station_file,material_file,   &
           Super_object_file,neumann_file,neumann_dat,check_mesh_file
       character (len=30) :: file_bassin
       character (len=1)  :: Super_object_type

       integer, dimension(0:8) :: out_variables
       integer                 :: nReqOut ! number of required outputs
       integer :: earthchunk_isInit
       character (len=MAX_FILE_SIZE) :: earthchunk_file
       real :: earthchunk_delta_lon, earthchunk_delta_lat

       real :: MPML_coeff

       !real, dimension(:, :), allocatable :: Mu, Lambda, Dens
       ! Configuration parameters as returned from read_input.c
       type(sem_config) :: config
    end type domain

contains


    subroutine dist_max_elem(Tdomain)
        implicit none
        type (Domain), intent (INOUT) :: Tdomain
        integer :: ipoint, jpoint
        integer :: i, j, n
        real :: coor_i(0:2), coor_j(0:2)
        real :: dist_max


        do n = 0,Tdomain%n_elem-1
            dist_max = 0.

            do i=0,Tdomain%n_nodes-1
                ipoint = Tdomain%specel(n)%Control_Nodes(i)
                coor_i = Tdomain%Coord_nodes(0:2,ipoint)
                do j=i+1,Tdomain%n_nodes-1
                    jpoint = Tdomain%specel(n)%Control_Nodes(j)
                    coor_j = Tdomain%Coord_nodes(0:2,jpoint)
                    dist_max = max(dist_max, sqrt((coor_i(0)-coor_j(0))**2 + (coor_i(1)-coor_j(1))**2 &
                        + (coor_i(2)-coor_j(2))**2))
                enddo
            enddo
            Tdomain%specel(n)%dist_max = dist_max
        enddo

    end subroutine dist_max_elem


end module sdomain
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
