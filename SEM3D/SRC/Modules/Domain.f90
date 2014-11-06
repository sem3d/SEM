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

    type :: domain

       ! Communicateur incluant les processeurs SEM uniquement
       integer :: communicateur
       ! Hors couplage : communicateur=communicateur_global
       ! mode couplage : communicateur incluant tous les codes
       integer :: communicateur_global
       ! Communicateur pour le reassemblage des sorties
       integer :: comm_output
       ! Nombre de processeur dans le groupe de communication associe a comm_output
       integer :: nb_output_procs
       integer :: output_rank
       integer, dimension(:), allocatable :: output_nodes, output_nodes_offset, output_elems
       ! Nombre de process par sorties pour le reassemblage
       integer :: ngroup
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

       integer :: n_source, n_dime, n_glob_nodes, n_mat, n_nodes, n_receivers, n_proc
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

       !!
       integer :: earthchunk_isInit
       character (len=MAX_FILE_SIZE) :: earthchunk_file
       real :: earthchunk_delta_lon, earthchunk_delta_lat

       real :: MPML_coeff

       !real, dimension(:, :), allocatable :: Mu, Lambda, Dens

    end type domain

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real function Comp_shapefunc(which_nod,xi,eta,zeta)

        integer :: which_nod
        real :: xi,eta,zeta

        select case (which_nod)
        case (0)
            Comp_shapefunc = 0.125 * xi*(xi-1) * eta*(eta-1) * zeta*(zeta-1)
        case (1)
            Comp_shapefunc = 0.125 * xi*(xi+1) * eta*(eta-1) * zeta*(zeta-1)
        case (2)
            Comp_shapefunc = 0.125 * xi*(xi+1) * eta*(eta+1) * zeta*(zeta-1)
        case (3)
            Comp_shapefunc = 0.125 * xi*(xi-1) * eta*(eta+1) * zeta*(zeta-1)
        case (4)
            Comp_shapefunc = 0.125 * xi*(xi-1) * eta*(eta-1) * zeta*(zeta+1)
        case (5)
            Comp_shapefunc = 0.125 * xi*(xi+1) * eta*(eta-1) * zeta*(zeta+1)
        case (6)
            Comp_shapefunc = 0.125 * xi*(xi+1) * eta*(eta+1) * zeta*(zeta+1)
        case (7)
            Comp_shapefunc = 0.125 * xi*(xi-1) * eta*(eta+1) * zeta*(zeta+1)
        case (8)
            Comp_shapefunc = 0.25 * (1-xi**2) * eta*(eta-1) * zeta*(zeta-1)
        case (9)
            Comp_shapefunc = 0.25 * xi*(xi+1) * (1-eta**2) * zeta*(zeta-1)
        case (10)
            Comp_shapefunc = 0.25 * (1-xi**2) * eta*(eta+1) * zeta*(zeta-1)
        case (11)
            Comp_shapefunc = 0.25 * xi*(xi-1) * (1-eta**2) * zeta*(zeta-1)
        case (12)
            Comp_shapefunc = 0.25 * xi*(xi-1) * eta*(eta-1) * (1-zeta**2)
        case (13)
            Comp_shapefunc = 0.25 * xi*(xi+1) * eta*(eta-1) * (1-zeta**2)
        case (14)
            Comp_shapefunc = 0.25 * xi*(xi+1) * eta*(eta+1) * (1-zeta**2)
        case (15)
            Comp_shapefunc = 0.25 * xi*(xi-1) * eta*(eta+1) * (1-zeta**2)
        case (16)
            Comp_shapefunc = 0.25 * (1-xi**2) * eta*(eta-1) * zeta*(zeta+1)
        case (17)
            Comp_shapefunc = 0.25 * xi*(xi+1) * (1-eta**2) * zeta*(zeta+1)
        case (18)
            Comp_shapefunc = 0.25 * (1-xi**2) * eta*(eta+1) * zeta*(zeta+1)
        case (19)
            Comp_shapefunc = 0.25 * xi*(xi-1) * (1-eta**2) * zeta*(zeta+1)
        case(20)
            Comp_shapefunc = 0.5 * (1-xi**2) * (1-eta**2) * zeta*(zeta-1)
        case(21)
            Comp_shapefunc = 0.5 * (1-xi**2) * eta*(eta-1) * (1-zeta**2)
        case(22)
            Comp_shapefunc = 0.5 * xi*(xi+1) * (1-eta**2) * (1-zeta**2)
        case(23)
            Comp_shapefunc = 0.5 * (1-xi**2) * eta*(eta+1) * (1-zeta**2)
        case(24)
            Comp_shapefunc = 0.5 * xi*(xi-1) * (1-eta**2) * (1-zeta**2)
        case(25)
            Comp_shapefunc = 0.5 * (1-xi**2) * (1-eta**2) * zeta*(zeta+1)
        case(26)
            Comp_shapefunc = (1-xi**2) * (1-eta**2) * (1-zeta**2)
        end select

        return
    end function Comp_shapefunc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    real function Comp_derivshapefunc(which_nod,xi,eta,zeta,compo)

        integer :: which_nod, compo
        real :: xi,eta,zeta

        real, dimension(0:2) :: df

        select case (which_nod)
        case (0)
            df(0) = 0.125 * (2*xi-1) * eta*(eta-1) * zeta*(zeta-1)
            df(1) = 0.125 * xi*(xi-1) * (2*eta-1) * zeta*(zeta-1)
            df(2) = 0.125 * xi*(xi-1) * eta*(eta-1) * (2*zeta-1)
        case (1)
            df(0) = 0.125 * (2*xi+1) * eta*(eta-1) * zeta*(zeta-1)
            df(1) = 0.125 * xi*(xi+1) * (2*eta-1) * zeta*(zeta-1)
            df(2) = 0.125 * xi*(xi+1) * eta*(eta-1) * (2*zeta-1)
        case (2)
            df(0) = 0.125 * (2*xi+1) * eta*(eta+1) * zeta*(zeta-1)
            df(1) = 0.125 * xi*(xi+1) * (2*eta+1) * zeta*(zeta-1)
            df(2) = 0.125 * xi*(xi+1) * eta*(eta+1) * (2*zeta-1)
        case (3)
            df(0) = 0.125 * (2*xi-1) * eta*(eta+1) * zeta*(zeta-1)
            df(1) = 0.125 * xi*(xi-1) * (2*eta+1) * zeta*(zeta-1)
            df(2) = 0.125 * xi*(xi-1) * eta*(eta+1) * (2*zeta-1)
        case (4)
            df(0) = 0.125 * (2*xi-1) * eta*(eta-1) * zeta*(zeta+1)
            df(1) = 0.125 * xi*(xi-1) * (2*eta-1) * zeta*(zeta+1)
            df(2) = 0.125 * xi*(xi-1) * eta*(eta-1) * (2*zeta+1)
        case (5)
            df(0) = 0.125 * (2*xi+1) * eta*(eta-1) * zeta*(zeta+1)
            df(1) = 0.125 * xi*(xi+1) * (2*eta-1) * zeta*(zeta+1)
            df(2) = 0.125 * xi*(xi+1) * eta*(eta-1) * (2*zeta+1)
        case (6)
            df(0) = 0.125 * (2*xi+1) * eta*(eta+1) * zeta*(zeta+1)
            df(1) = 0.125 * xi*(xi+1) * (2*eta+1) * zeta*(zeta+1)
            df(2) = 0.125 * xi*(xi+1) * eta*(eta+1) * (2*zeta+1)
        case (7)
            df(0) = 0.125 * (2*xi-1) * eta*(eta+1) * zeta*(zeta+1)
            df(1) = 0.125 * xi*(xi-1) * (2*eta+1) * zeta*(zeta+1)
            df(2) = 0.125 * xi*(xi-1) * eta*(eta+1) * (2*zeta+1)
        case (8)
            df(0) = 0.25 * (-2*xi) * eta*(eta-1) * zeta*(zeta-1)
            df(1) = 0.25 * (1-xi**2) * (2*eta-1) * zeta*(zeta-1)
            df(2) = 0.25 * (1-xi**2) * eta*(eta-1) * (2*zeta-1)
        case (9)
            df(0) = 0.25 * (2*xi+1) * (1-eta**2) * zeta*(zeta-1)
            df(1) = 0.25 * xi*(xi+1) * (-2*eta) * zeta*(zeta-1)
            df(2) = 0.25 * xi*(xi+1) * (1-eta**2) * (2*zeta-1)
        case (10)
            df(0) = 0.25 * (-2*xi) * eta*(eta+1) * zeta*(zeta-1)
            df(1) = 0.25 * (1-xi**2) * (2*eta+1) * zeta*(zeta-1)
            df(2) = 0.25 * (1-xi**2) * eta*(eta+1) * (2*zeta-1)
        case (11)
            df(0) = 0.25 * (2*xi-1) * (1-eta**2) * zeta*(zeta-1)
            df(1) = 0.25 * xi*(xi-1) * (-2*eta) * zeta*(zeta-1)
            df(2) = 0.25 * xi*(xi-1) * (1-eta**2) * (2*zeta-1)
        case (12)
            df(0) = 0.25 * (2*xi-1) * eta*(eta-1) * (1-zeta**2)
            df(1) = 0.25 * xi*(xi-1) * (2*eta-1) * (1-zeta**2)
            df(2) = 0.25 * xi*(xi-1) * eta*(eta-1) * (-2*zeta)
        case (13)
            df(0) = 0.25 * (2*xi+1) * eta*(eta-1) * (1-zeta**2)
            df(1) = 0.25 * xi*(xi+1) * (2*eta-1) * (1-zeta**2)
            df(2) = 0.25 * xi*(xi+1) * eta*(eta-1) * (-2*zeta)
        case (14)
            df(0) = 0.25 * (2*xi+1) * eta*(eta+1) * (1-zeta**2)
            df(1) = 0.25 * xi*(xi+1) * (2*eta+1) * (1-zeta**2)
            df(2) = 0.25 * xi*(xi+1) * eta*(eta+1) * (-2*zeta)
        case (15)
            df(0) = 0.25 * (2*xi-1) * eta*(eta+1) * (1-zeta**2)
            df(1) = 0.25 * xi*(xi-1) * (2*eta+1) * (1-zeta**2)
            df(2) = 0.25 * xi*(xi-1) * eta*(eta+1) * (-2*zeta)
        case (16)
            df(0) = 0.25 * (-2*xi) * eta*(eta-1) * zeta*(zeta+1)
            df(1) = 0.25 * (1-xi**2) * (2*eta-1) * zeta*(zeta+1)
            df(2) = 0.25 * (1-xi**2) * eta*(eta-1) * (2*zeta+1)
        case (17)
            df(0) = 0.25 * (2*xi+1) * (1-eta**2) * zeta*(zeta+1)
            df(1) = 0.25 * xi*(xi+1) * (-2*eta) * zeta*(zeta+1)
            df(2) = 0.25 * xi*(xi+1) * (1-eta**2) * (2*zeta+1)
        case (18)
            df(0) = 0.25 * (-2*xi) * eta*(eta+1) * zeta*(zeta+1)
            df(1) = 0.25 * (1-xi**2) * (2*eta+1) * zeta*(zeta+1)
            df(2) = 0.25 * (1-xi**2) * eta*(eta+1) * (2*zeta+1)
        case (19)
            df(0) = 0.25 * (2*xi-1) * (1-eta**2) * zeta*(zeta+1)
            df(1) = 0.25 * xi*(xi-1) * (-2*eta) * zeta*(zeta+1)
            df(2) = 0.25 * xi*(xi-1) * (1-eta**2) * (2*zeta+1)
        case(20)
            df(0) = 0.5 * (-2*xi) * (1-eta**2) * zeta*(zeta-1)
            df(1) = 0.5 * (1-xi**2) * (-2*eta) * zeta*(zeta-1)
            df(2) = 0.5 * (1-xi**2) * (1-eta**2) * (2*zeta-1)
        case(21)
            df(0) = 0.5 * (-2*xi) * eta*(eta-1) * (1-zeta**2)
            df(1) = 0.5 * (1-xi**2) * (2*eta-1) * (1-zeta**2)
            df(2) = 0.5 * (1-xi**2) * eta*(eta-1) * (-2*zeta)
        case(22)
            df(0) = 0.5 * (2*xi+1) * (1-eta**2) * (1-zeta**2)
            df(1) = 0.5 * xi*(xi+1) * (-2*eta) * (1-zeta**2)
            df(2) = 0.5 * xi*(xi+1) * (1-eta**2) * (-2*zeta)
        case(23)
            df(0) = 0.5 * (-2*xi) * eta*(eta+1) * (1-zeta**2)
            df(1) = 0.5 * (1-xi**2) * (2*eta+1) * (1-zeta**2)
            df(2) = 0.5 * (1-xi**2) * eta*(eta+1) * (-2*zeta)
        case(24)
            df(0) = 0.5 * (2*xi-1) * (1-eta**2) * (1-zeta**2)
            df(1) = 0.5 * xi*(xi-1) * (-2*eta) * (1-zeta**2)
            df(2) = 0.5 * xi*(xi-1) * (1-eta**2) * (-2*zeta)
        case(25)
            df(0) = 0.5 * (-2*xi) * (1-eta**2) * zeta*(zeta+1)
            df(1) = 0.5 * (1-xi**2) * (-2*eta) * zeta*(zeta+1)
            df(2) = 0.5 * (1-xi**2) * (1-eta**2) * (2*zeta+1)
        case(26)
            df(0) = (-2*xi) * (1-eta**2) * (1-zeta**2)
            df(1) = (1-xi**2) * (-2*eta) * (1-zeta**2)
            df(2) = (1-xi**2) * (1-eta**2) * (-2*zeta)
        end select

        Comp_derivshapefunc = df(compo)

        return
    end function Comp_derivshapefunc

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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
