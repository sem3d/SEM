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
    use sextended_sources
    use stimeparam
    use logical_input
    use ssubdomains
    use splanew
    use sneu
    use ssurf
    use solid_fluid
    use semdatafiles
    use sem_c_config
    use sinterface
    use champs_solid
    use champs_solid_dg
    use champs_solidpml
    use champs_fluid
    use champs_fluidpml
    use constants
    use msnapdata, only : output_var_t
    implicit none


    type :: domain
       integer :: communicateur !<<< Communicator including all SEM processors
       integer :: rank          !<<< Rank of this process within this communicator
       integer :: nb_procs      !<<< Total number of SEM processors
       ! Without coupling : communicateur=communicateur_global
       ! With coupling    : communicateur : includes every processes
       integer :: communicateur_global
       ! Nombre de processeur avec qui on communique (size(sComm))
       integer :: tot_comm_proc
       ! En mode couplage : Rg du superviseur dans le communicateur global
       integer :: master_superviseur
       !
       ! Number of procs per snapshot groups (should be the number of MPI processes per node)
       integer :: ngroup
       ! Holds temp memory used to gather outputs for each group
       type(output_var_t)  :: SnapData
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
       type(extended_source) , dimension (:), allocatable :: sExtendSource

       logical :: aniso
       logical :: has_station
       logical :: nl_flag
       logical :: use_avg
       integer :: prot_at_time
       double precision :: time_0
       logical :: first_prot=.true.
       integer :: time_now_s
       integer :: n_source, n_extsource, n_dime, n_glob_nodes, n_mat, n_nodes, n_receivers
       integer :: n_elem, n_face, n_edge, n_vertex, n_glob_points, n_sls, n_neumannfind
       integer :: n_hexa  !< Nombre de maille hexa ~= (ngllx-1)*(nglly-1)*(ngllz-1)*nelem
       integer :: n_hexa_local !< Nombre de subelements hexa dans le proc(division aux GLLs)
       integer :: ngll
       logical, dimension(:), allocatable :: not_PML_List, subD_exist
       logical :: any_sdom, any_fdom, any_spml, any_fpml, any_sdomdg

       real(fpp) :: dxmax

       real(fpp) :: T1_att, T2_att, T0_modele
       real(fpp), dimension (0:2,0:2) :: rot
       real(fpp), dimension (:,:), allocatable:: Coord_nodes
       real(fpp), dimension (:,:), pointer:: GlobCoord ! No allocate, use pointer: enable pointing to coord (avoid allocate + copy)
       ! Mirror
       logical :: use_mirror,mirror_expl,mirror_recalc
       integer :: mirror_type

       integer :: traces_format
       character (len=MAX_FILE_SIZE) :: Title_simulation, mesh_file,material_file


       integer, dimension(0:OUT_LAST) :: out_var_capt
       integer, dimension(0:OUT_LAST) :: out_var_snap
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
       type(domain_solid_dg)    :: sdomdg

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
       integer                            :: n_NEBC, n_PWBC, n_FTBC, nsurface, n_DIRIC
       integer, dimension(:), allocatable :: list_NEBC, list_PWBC, list_FTBC, list_DIRICBC
    end type domain

contains

    function domain_from_type_char(ch)
        character, intent(in) :: ch
        integer :: domain_from_type_char
        select case(ch)
        case('S')
            domain_from_type_char = DM_SOLID_CG
        case('P')
            domain_from_type_char = DM_SOLID_CG_PML
        case('F')
            domain_from_type_char = DM_FLUID_CG
        case('L')
            domain_from_type_char = DM_FLUID_CG_PML
        case('D')
            domain_from_type_char = DM_SOLID_DG
        case('E')
            domain_from_type_char = DM_FLUID_DG
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
        case(DM_SOLID_CG)
            domain_nglltot = Tdomain%sdom%nglltot
        case(DM_FLUID_CG)
            domain_nglltot = Tdomain%fdom%nglltot
        case(DM_SOLID_CG_PML)
            domain_nglltot = Tdomain%spmldom%nglltot
        case(DM_FLUID_CG_PML)
            domain_nglltot = Tdomain%fpmldom%nglltot
        case(DM_SOLID_DG)
            domain_nglltot = Tdomain%sdomdg%nglltot
        case default
            stop "Unknown Domain, nglltot"
        end select
    end function domain_nglltot

    function domain_ngll(Tdomain, dom)
        integer, intent(in) :: dom
        type(domain), intent(in) :: Tdomain
        !
        integer :: domain_ngll
        domain_ngll = 0
        select case(dom)
        case(DM_SOLID_CG)
            domain_ngll = Tdomain%sdom%ngll
        case(DM_FLUID_CG)
            domain_ngll = Tdomain%fdom%ngll
        case(DM_SOLID_CG_PML)
            domain_ngll = Tdomain%spmldom%ngll
        case(DM_FLUID_CG_PML)
            domain_ngll = Tdomain%fpmldom%ngll
        case(DM_SOLID_DG)
            domain_ngll = Tdomain%sdomdg%ngll
        case default
            stop "Unknown Domain, ngll"
        end select
    end function domain_ngll

    subroutine domain_nelems(Tdomain, dom, nelems, nglltot)
        integer, intent(in) :: dom
        type(domain), intent(in) :: Tdomain
        integer, intent(out) :: nelems, nglltot
        !
        nelems = 0
        nglltot = 0
        select case(dom)
        case(DM_SOLID_CG)
            nelems  = Tdomain%sdom%nbelem
            nglltot = Tdomain%sdom%nglltot
        case(DM_FLUID_CG)
            nelems  = Tdomain%fdom%nbelem
            nglltot = Tdomain%fdom%nglltot
        case(DM_SOLID_CG_PML)
            nelems  = Tdomain%spmldom%nbelem
            nglltot = Tdomain%spmldom%nglltot
        case(DM_FLUID_CG_PML)
            nelems  = Tdomain%fpmldom%nbelem
            nglltot = Tdomain%fpmldom%nglltot
        case(DM_SOLID_DG)
            nelems  = Tdomain%sdomdg%nbelem
            nglltot = Tdomain%sdomdg%nglltot
        case(DM_FLUID_DG)
            nelems  = 0
            nglltot = 0
        case default
            stop "Unknown Domain, nelems"
        end select
    end subroutine domain_nelems

    subroutine domain_gllc(Tdomain, dom, GLLc)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: dom
        real(fpp), dimension(:), allocatable, intent(out) :: gllc
        !
        integer :: ngll
        select case (dom)
        case (DM_SOLID_CG)
            ngll = Tdomain%sdom%ngll
            allocate(GLLc(0:ngll-1))
            GLLc = Tdomain%sdom%GLLc
        case (DM_FLUID_CG)
            ngll = Tdomain%fdom%ngll
            allocate(GLLc(0:ngll-1))
            GLLc = Tdomain%fdom%GLLc
        case (DM_SOLID_CG_PML)
            ngll = Tdomain%spmldom%ngll
            allocate(GLLc(0:ngll-1))
            GLLc = Tdomain%spmldom%GLLc
        case (DM_FLUID_CG_PML)
            ngll = Tdomain%fpmldom%ngll
            allocate(GLLc(0:ngll-1))
            GLLc = Tdomain%fpmldom%GLLc
        case (DM_SOLID_DG)
            ngll = Tdomain%sdomdg%ngll
            allocate(GLLc(0:ngll-1))
            GLLc = Tdomain%sdomdg%GLLc
        case default
            stop "Unknown Domain, gllc"
        end select
    end subroutine domain_gllc

    subroutine domain_gllw(Tdomain, dom, gllw)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: dom
        real(fpp), dimension(:), allocatable, intent(out) :: gllw
        !
        integer :: ngll
        select case (dom)
        case (DM_SOLID_CG)
            ngll = Tdomain%sdom%ngll
            allocate(gllw(0:ngll-1))
            gllw = Tdomain%sdom%GLLw
        case (DM_FLUID_CG)
            ngll = Tdomain%fdom%ngll
            allocate(gllw(0:ngll-1))
            gllw = Tdomain%fdom%GLLw
        case (DM_SOLID_CG_PML)
            ngll = Tdomain%spmldom%ngll
            allocate(gllw(0:ngll-1))
            gllw = Tdomain%spmldom%GLLw
        case (DM_FLUID_CG_PML)
            ngll = Tdomain%fpmldom%ngll
            allocate(gllw(0:ngll-1))
            gllw = Tdomain%fpmldom%GLLw
        case (DM_SOLID_DG)
            ngll = Tdomain%sdomdg%ngll
            allocate(gllw(0:ngll-1))
            gllw = Tdomain%sdomdg%GLLw
        case default
            stop "Unknown Domain, gllw"
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
            ngll = domain_ngll(Tdomain, Tdomain%sFace(nf0)%domain)
            bad = .false.
            do j=0,ngll-1
                do i=0,ngll-1
                    ip0 = Tdomain%sFace(nf0)%Iglobnum_Face(i,j)
                    ip1 = Tdomain%sFace(nf1)%Iglobnum_Face(i,j)
                    dp = abs(Tdomain%GlobCoord(:,ip0)- Tdomain%GlobCoord(:,ip1))
                    dp = dp/Tdomain%dxmax
                    if (dp(0)>xeps.or.dp(1)>xeps.or.dp(2)>xeps) then
                        !write(*,*) "IF ERROR:", Tdomain%GlobCoord(:,ip0), Tdomain%GlobCoord(:,ip1), dp
                        bad = .true.
                    end if
                end do
            end do
            if (bad) then
                write(*,*) "BAD INTERFACE:", dp, Tdomain%dxmax
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
                        write(*,*) "dx:", nf0, nf1, i,j, ":", x0, x1, x1-x0
                    end do
                end do
            end if
        end do
    end subroutine check_interface_orient

    subroutine map_surface_faces_to_elem(Tdomain, ngll, surf, renum, map)
        use mindex, only : ind_elem_face, face_def
        type(domain), intent(in) :: Tdomain
        type(surf_num), intent(in) :: surf
        integer, intent(in) :: ngll
        integer, intent(in), dimension(:) :: renum
        integer, intent(out), allocatable, dimension(:,:,:,:) :: map
        !
        integer :: idxf, idxi, idxj, idxk
        integer :: nf, nnf, i, j, k
        integer, dimension(0:2) :: i0, di, dj
        integer, dimension(0:3) :: elface
        integer :: nel

        ! Use like that :
        !call get_surface_numbering(Tdomain, Tdomain%SF%intSolFluPml%surf0, DM_FLUID_CG_PML, renum)
        !call map_surface_faces_to_elem(Tdomain, dom%ngll, Tdomain%SF%intSolFluPml%surf0, renum, dom%sf_map)

         allocate(map(0:surf%n_faces-1, 0:ngll-1,0:ngll-1, 6)) ! 6=lnum,i,j,k,idom,imap
         map = -1 ! init to -1, need to identify and skip orphan faces
         do idxf = 0, surf%n_faces-1
             nnf = surf%if_faces(idxf)
             if (Tdomain%sFace(nnf)%orphan) cycle
             ! Since it's an interface face, there's only one element that can be associated
             nel = Tdomain%sFace(nnf)%elem_0
             if (nel==-1) nel = Tdomain%sFace(nnf)%elem_1

             ! We have a face, need to know which face we are
             do nf=0,5
                 if (nnf == Tdomain%specel(nel)%Near_Faces(nf)) exit
             end do
             do k=0,3
                 elface(k) = Tdomain%specel(nel)%Control_nodes(face_def(k,nf))
             end do
             call ind_elem_face(ngll, nf, Tdomain%sFace(nnf)%inodes, elface, i0, di, dj)

             do i=0,ngll-1
                 do j=0,ngll-1
                     idxi = i0(0)+i*di(0)+j*dj(0)
                     idxj = i0(1)+i*di(1)+j*dj(1)
                     idxk = i0(2)+i*di(2)+j*dj(2)
                     map(idxf, i, j, 0) = Tdomain%specel(nel)%lnum
                     map(idxf, i, j, 1) = idxi
                     map(idxf, i, j, 2) = idxj
                     map(idxf, i, j, 3) = idxk
                     map(idxf, i, j, 4) = Tdomain%specel(nel)%Idom(i,j,k)
                     map(idxf, i, j, 5) = renum(Tdomain%specel(nel)%Idom(i,j,k))
                 end do
             end do
         end do
    end subroutine map_surface_faces_to_elem
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
