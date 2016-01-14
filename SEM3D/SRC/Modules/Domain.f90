!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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
    use sinterface
    use champs_solid
    use champs_solidpml
    use champs_fluid
    use champs_fluidpml
    use constants
    implicit none

    type domain_solid
        ! D'abord, les données membres qui ne sont pas modifiées

        ! Nombre de gll
        integer :: ngll

        ! Nombre d'elements dans le domaine
        integer :: nbelem

        ! MassMat pour elements solide, fluide, solide pml et fluide pml
        real, dimension(:), allocatable :: MassMat

        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champssolid) :: champs0
        type(champssolid) :: champs1
    end type domain_solid

    type domain_solidpml
        ! D'abord, les données membres qui ne sont pas modifiées

        ! Nombre de gll
        integer :: ngll

        ! Nombre d'elements dans le domaine
        integer :: nbelem

        ! MassMat pour elements solide, fluide, solide pml et fluide pml
        real, dimension(:), allocatable :: MassMat
        real, dimension(:,:), allocatable :: DumpMass

        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champssolidpml) :: champs0
        type(champssolidpml) :: champs1
    end type domain_solidpml

    type domain_fluid
        ! D'abord, les données membres qui ne sont pas modifiées

        ! Nombre de gll
        integer :: ngll

        ! Nombre d'elements dans le domaine
        integer :: nbelem

        ! MassMat pour elements solide, fluide, solide pml et fluide pml
        real, dimension(:), allocatable :: MassMat

        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champsfluid) :: champs0
        type(champsfluid) :: champs1
    end type domain_fluid

    type domain_fluidpml
        ! D'abord, les données membres qui ne sont pas modifiées

        ! Nombre de gll
        integer :: ngll

        ! Nombre d'elements dans le domaine
        integer :: nbelem

        ! MassMat pour elements solide, fluide, solide pml et fluide pml
        real(fpp), dimension(:), allocatable :: MassMat
        real(fpp), dimension(:,:), allocatable :: DumpMass
        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champsfluidpml) :: champs0
        type(champsfluidpml) :: champs1
        real(fpp), dimension(:,:,:,:,:), allocatable :: Veloc
    end type domain_fluidpml

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
       type(bassin)        :: sBassin
       type(source)   , dimension (:), pointer :: sSource
       type(element)  , dimension (:), pointer :: specel
       type(face)     , dimension (:), pointer :: sFace
       type(comm)     , dimension (:), pointer :: sComm
       type(edge)     , dimension (:), pointer :: sEdge
       type(vertex)   , dimension (:), pointer :: sVertex
       type(subdomain), dimension (:), pointer :: sSubDomain
       type(SurfaceT), dimension(:), allocatable :: sSurfaces

       logical :: aniso
       logical :: any_Random, any_PropOnFile

       integer :: n_source, n_dime, n_glob_nodes, n_mat, n_nodes, n_receivers
       integer :: n_elem, n_face, n_edge, n_vertex, n_glob_points, n_sls
       integer :: n_hexa  !< Nombre de maille hexa ~= (ngllx-1)*(nglly-1)*(ngllz-1)*nelem
       logical, dimension(:), allocatable :: not_PML_List, subD_exist

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

       ! Domains
       type(domain_solid)    :: sdom
       type(domain_solidpml) :: spmldom
       type(domain_fluid)    :: fdom
       type(domain_fluidpml) :: fpmldom

       ! Interface Solide / PML
       type(inter_num) :: intSolPml
       ! Interface Fluide / PML
       type(inter_num) :: intFluPml
       ! Interface Fluide / Solide
       type(SF_object)     :: SF

        ! Communication
        type(comm_vector) :: Comm_data      ! Comm mass et forces
        type(comm_vector) :: Comm_SolFlu    ! Comm couplage solide/fluide

       ! Configuration parameters as returned from read_input.c
       type(sem_config) :: config
    end type domain

contains

    function domain_ngll(Tdomain, dom)
        integer, intent(in) :: dom
        type(domain), intent(in) :: Tdomain
        !
        integer :: domain_ngll
        domain_ngll = 0
        select case(dom)
        case(DM_SOLID)
            domain_ngll = Tdomain%sdom%ngll
        case(DM_FLUID)
            domain_ngll = Tdomain%fdom%ngll
        case(DM_SOLID_PML)
            domain_ngll = Tdomain%spmldom%ngll
        case(DM_FLUID_PML)
            domain_ngll = Tdomain%fpmldom%ngll
        end select
    end function domain_ngll

    subroutine check_interface_orient(Tdomain, inter, xeps)
        type(domain), intent(in) :: Tdomain
        type(inter_num), intent(in) :: inter
        real(FPP), intent(in) :: xeps
        !
        integer :: nif, i, j
        integer :: nf0, nf1, ip0, ip1
        integer :: ngll1, ngll2
        real(FPP), dimension(0:2) :: dp
        logical :: bad
        !
        do nif=0,inter%surf0%n_faces-1
            nf0 = inter%surf0%if_faces(nif)
            nf1 = inter%surf1%if_faces(nif)
            ngll1 = Tdomain%sFace(nf0)%ngll1
            ngll2 = Tdomain%sFace(nf0)%ngll2
            bad = .false.
            do j=0,ngll2-1
                do i=0,ngll1-1
                    ip0 = Tdomain%sFace(nf0)%Iglobnum_Face(i,j)
                    ip1 = Tdomain%sFace(nf1)%Iglobnum_Face(i,j)
                    dp = abs(Tdomain%GlobCoord(:,ip0)- Tdomain%GlobCoord(:,ip1))
                    if (dp(0)>xeps.or.dp(1)>xeps.or.dp(2)>xeps) then
                        !write(*,*) "IF ERROR:", Tdomain%GlobCoord(:,ip0), Tdomain%GlobCoord(:,ip1), dp
                        bad = .true.
                    end if
                end do
            end do
            if (bad) then
                write(*,*) "BAD INTERFACE:"
                write(*,*) "IF0:", nf0, "[", Tdomain%sFace(nf0)%inodes, "]"
                write(*,*) "IF1:", nf1, "[", Tdomain%sFace(nf1)%inodes, "]"
            end if
        end do
    end subroutine check_interface_orient

end module sdomain

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
