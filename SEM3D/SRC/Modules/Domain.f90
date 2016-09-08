!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Domain.f90
!!\brief Contient les definition du type domain
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

       type(SurfaceParam), dimension(:), allocatable :: nsurfsource
       type(source)   , dimension (:), allocatable :: sSource
       type(element)  , dimension (:), pointer     :: specel
       type(face)     , dimension (:), allocatable :: sFace
       type(edge)     , dimension (:), allocatable :: sEdge
       type(vertex)   , dimension (:), allocatable :: sVertex
       type(subdomain), dimension (:), pointer     :: sSubDomain
       type(SurfaceT) , dimension (:), allocatable :: sSurfaces
       type(comm)     , dimension (:), allocatable :: sComm

       logical :: aniso
       logical :: any_Random, any_PropOnFile
       logical :: nl_flag
       integer :: nRandom
       integer :: n_source, n_dime, n_glob_nodes, n_mat, n_nodes, n_receivers
       integer :: n_elem, n_face, n_edge, n_vertex, n_glob_points, n_sls, n_neumannfind
       integer :: n_hexa  !< Nombre de maille hexa ~= (ngllx-1)*(nglly-1)*(ngllz-1)*nelem
       integer :: n_hexa_local !< Nombre de subelements hexa dans le proc(division aux GLLs)
       logical, dimension(:), allocatable :: not_PML_List, subD_exist
       logical :: any_sdom, any_fdom, any_spml, any_fpml

       real(fpp) :: T1_att, T2_att, T0_modele
       real(fpp), dimension (0:2,0:2) :: rot
       real(fpp), dimension (:,:), allocatable:: Coord_nodes
       real(fpp), dimension (:,:), pointer:: GlobCoord ! No allocate, use pointer: enable pointing to coord (avoid allocate + copy)

       integer :: traces_format
       character (len=MAX_FILE_SIZE) :: Title_simulation, mesh_file,station_file,material_file,   &
           Super_object_file,neumann_file,neumann_dat,check_mesh_file
       character (len=1)  :: Super_object_type

       integer, dimension(0:8) :: out_variables
       integer :: out_energy
       integer                 :: nReqOut ! number of required outputs
       integer :: earthchunk_isInit
       character (len=MAX_FILE_SIZE) :: earthchunk_file
       real(fpp) :: earthchunk_delta_lon, earthchunk_delta_lat

       real(fpp) :: MPML_coeff

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
       type(sem_material_list) :: material_list

       ! Rajouter pour le traitement des surface
       integer :: n_NEBC, n_PWBC, n_FTBC, nsurface
    end type domain

contains

    function domain_from_type_char(ch)
        character, intent(in) :: ch
        integer :: domain_from_type_char
        select case(ch)
        case('R')
            domain_from_type_char = DM_SOLID
        case('S')
            domain_from_type_char = DM_SOLID
        case('P')
            domain_from_type_char = DM_SOLID_PML
        case('F')
            domain_from_type_char = DM_FLUID
        case('L')
            domain_from_type_char = DM_FLUID_PML
        case default
            stop "Unknown material type"
        end select
    end function domain_from_type_char

    function domain_nglltot(Tdomain, dom)
        integer, intent(in) :: dom
        type(domain), intent(in) :: Tdomain
        !
        integer :: domain_nglltot
        domain_nglltot = 0
        select case(dom)
        case(DM_SOLID)
            domain_nglltot = Tdomain%sdom%nglltot
        case(DM_FLUID)
            domain_nglltot = Tdomain%fdom%nglltot
        case(DM_SOLID_PML)
            domain_nglltot = Tdomain%spmldom%nglltot
        case(DM_FLUID_PML)
            domain_nglltot = Tdomain%fpmldom%nglltot
        case default
            stop "Unknown Domain"
        end select
    end function domain_nglltot

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
        case default
            stop "Unknown Domain"
        end select
    end function domain_ngll

    subroutine domain_gllc(Tdomain, dom, GLLc)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: dom
        real(fpp), dimension(:), allocatable, intent(out) :: gllc
        !
        integer :: ngll
        select case (dom)
        case (DM_SOLID)
            ngll = Tdomain%sdom%ngll
            allocate(GLLc(0:ngll-1))
            GLLc = Tdomain%sdom%GLLc
        case (DM_FLUID)
            ngll = Tdomain%fdom%ngll
            allocate(GLLc(0:ngll-1))
            GLLc = Tdomain%fdom%GLLc
        case (DM_SOLID_PML)
            ngll = Tdomain%spmldom%ngll
            allocate(GLLc(0:ngll-1))
            GLLc = Tdomain%spmldom%GLLc
        case (DM_FLUID_PML)
            ngll = Tdomain%fpmldom%ngll
            allocate(GLLc(0:ngll-1))
            GLLc = Tdomain%fpmldom%GLLc
        case default
            stop "Unknown Domain"
        end select
    end subroutine domain_gllc

    subroutine domain_gllw(Tdomain, dom, gllw)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: dom
        real(fpp), dimension(:), allocatable, intent(out) :: gllw
        !
        integer :: ngll
        select case (dom)
        case (DM_SOLID)
            ngll = Tdomain%sdom%ngll
            allocate(gllw(0:ngll-1))
            gllw = Tdomain%sdom%GLLw
        case (DM_FLUID)
            ngll = Tdomain%fdom%ngll
            allocate(gllw(0:ngll-1))
            gllw = Tdomain%fdom%GLLw
        case (DM_SOLID_PML)
            ngll = Tdomain%spmldom%ngll
            allocate(gllw(0:ngll-1))
            gllw = Tdomain%spmldom%GLLw
        case (DM_FLUID_PML)
            ngll = Tdomain%fpmldom%ngll
            allocate(gllw(0:ngll-1))
            gllw = Tdomain%fpmldom%GLLw
        case default
            stop "Unknown Domain"
        end select
    end subroutine domain_gllw


    subroutine check_interface_orient(Tdomain, inter, xeps)
        type(domain), intent(in) :: Tdomain
        type(inter_num), intent(in) :: inter
        real(FPP), intent(in) :: xeps
        !
        integer :: nif, i, j
        integer :: nf0, nf1, ip0, ip1
        integer :: ngll
        real(FPP), dimension(0:2) :: dp, x0, x1
        logical :: bad
        !
        do nif=0,inter%surf0%n_faces-1
            nf0 = inter%surf0%if_faces(nif)
            nf1 = inter%surf1%if_faces(nif)
            ngll = Tdomain%sFace(nf0)%ngll
            bad = .false.
            do j=0,ngll-1
                do i=0,ngll-1
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
                write(*,*) "IF0: INODES", nf0, "[", Tdomain%sFace(nf0)%inodes, "]"
                write(*,*) "IF1: INODES", nf1, "[", Tdomain%sFace(nf1)%inodes, "]"
                write(*,*) "IF0: IGLOBN", nf0, "[", Tdomain%sFace(nf0)%Iglobnum_Face, "]"
                write(*,*) "IF1: IGLOBN", nf1, "[", Tdomain%sFace(nf1)%Iglobnum_Face, "]"
                write(*,*) "DX:", nf0, nf1, ":"
                do j=0,ngll-1
                    do i=0,ngll-1
                        ip0 = Tdomain%sFace(nf0)%Iglobnum_Face(i,j)
                        ip1 = Tdomain%sFace(nf1)%Iglobnum_Face(i,j)
                        if (ip0>=0) then
                            x0 = Tdomain%GlobCoord(:,ip0)
                        else
                            x0 = 0
                        end if
                        if (ip1>=0) then
                            x1 = Tdomain%GlobCoord(:,ip1)
                        else
                            x1 = 0.
                        end if
                        write(*,*) "dx:", nf0, nf1, i,j, ":", x0, x1
                    end do
                end do
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
