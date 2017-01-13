!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!



module sem_c_config
    use iso_c_binding
    use sem_c_bindings

    integer, parameter :: NAME_MAX=200
    ! Ce type doit correspondre au type sem_config_t du module read_input.c **a l'ordre pres**
    type, bind(c) :: sem_config
       type(C_PTR)    :: run_name

       !! Integration
       integer(C_INT) :: type_timeinteg
       integer(C_INT) :: implicitness
       integer(C_INT) :: accel_scheme
       integer(C_INT) :: veloc_scheme
       real(C_DOUBLE) :: sim_time
       real(C_DOUBLE) :: alpha
       real(C_DOUBLE) :: beta
       real(C_DOUBLE) :: gamma
       real(C_DOUBLE) :: courant
       real(C_DOUBLE) :: fmax
       integer(C_INT) :: ngll
       integer(C_INT) :: dim

       !! Modele, maillage
       type(C_PTR)    :: mesh_file
       integer(C_INT) :: model
       integer(C_INT) :: anisotropy
       type(C_PTR)    :: mat_file
       integer(C_INT) :: nsources
       type(C_PTR)    :: source
       integer(C_INT) :: nextended_sources
       type(C_PTR)    :: extended_source
       integer(C_INT) :: is_lamb_test

       !! Capteurs
       integer(C_INT) :: save_traces
       integer(C_INT) :: traces_interval
       integer(C_INT) :: traces_format
       type(C_PTR)    :: station_file
       integer(C_INT) :: capt_loc_type
       integer(C_INT) :: post_processing

       !! Snapshots
       integer(C_INT) :: save_snap
       integer(C_INT) :: n_group_outputs
       real(C_DOUBLE) :: snap_interval
       integer(C_INT) :: n_snap_cond
       type(C_PTR)    :: snapshot_selection
       integer(C_INT) :: comp_energ

       !! Output Variables
       integer(C_INT), dimension(11) :: out_variables
       integer(C_INT) :: nl_flag

       !! Protection reprise
       integer(C_INT) :: prorep
       integer(C_INT) :: prorep_iter
       integer(C_INT) :: prorep_restart_iter

       integer(C_INT) :: verbose_level
       real(C_DOUBLE) :: mpml

       !! Amortissement
       integer(C_INT) :: nsolids
       real(C_DOUBLE), dimension(2) :: atn_band
       real(C_DOUBLE) :: atn_period

       !! PML informations
       integer(C_INT) :: pml_type
       real(C_DOUBLE) :: cpml_kappa0
       real(C_DOUBLE) :: cpml_kappa1
       integer(C_INT) :: cpml_n
       real(C_DOUBLE) :: cpml_rc

       !! Type Elements (DG)
       integer(C_INT) :: type_elem
       integer(C_INT) :: type_flux
       integer(C_INT) :: type_bc

       !!Material
       integer(C_INT) :: material_present
       integer(C_INT) :: material_type
       type(C_PTR)    :: model_file
       real(C_DOUBLE) :: delta_lon
       real(C_DOUBLE) :: delta_lat

       type(C_PTR)    :: stations

       integer(C_INT) :: nsurface
       integer(C_INT) :: surface_find
       type(C_PTR)    :: surface
    end type sem_config


    ! Ce type doit correspondre au type source_t de sem_input.h **a l'ordre pres**
    type, bind(c) :: sem_source
       type(C_PTR) :: next
       real(C_DOUBLE), dimension(3) :: coords;
       integer(C_INT) :: type;
       real(C_DOUBLE), dimension(3) :: dir;
       integer(C_INT) :: func;
       real(C_DOUBLE), dimension(6) :: moments;
       real(C_DOUBLE) :: tau;
       real(C_DOUBLE) :: freq;
       real(C_DOUBLE), dimension(4) :: band;
       real(C_DOUBLE) :: ts
       real(C_DOUBLE) :: gamma
       real(C_DOUBLE) :: amplitude
       real(C_DOUBLE) :: sigma
       type(C_PTR) :: time_file
       real(C_DOUBLE) :: Q
       real(C_DOUBLE) :: X
       real(C_DOUBLE) :: Y
       real(C_DOUBLE) :: L
       real(C_DOUBLE) :: v
       real(C_DOUBLE) :: d
       real(C_DOUBLE) :: a
    end type sem_source


    ! Ce type doit correspondre au type extended_source_t de sem_input.h **a l'ordre pres**
    type, bind(c) :: sem_extended_source
       type(C_PTR) :: next
       type(C_PTR) :: kine_file
       type(C_PTR) :: slip_file
       real(C_DOUBLE) :: dip
       real(C_DOUBLE) :: strike
       real(C_DOUBLE) :: rake
    end type sem_extended_source

    ! ce type doit correspondre au type station_def_t de sem_input.h **a l'ordre pres**
    type, bind(c) :: sem_station
       type(C_PTR) :: next
       real(C_DOUBLE), dimension(3) :: coords;
       type(C_PTR)    :: name
       integer(C_INT) :: period
    end type sem_station

    type, bind(c) :: sem_snapshot_cond
       type(C_PTR) :: next
       integer(C_INT) :: type
       integer(C_INT) :: include
       real(C_DOUBLE), dimension(6) :: box
       real(C_DOUBLE), dimension(4) :: plane
       integer(C_INT) :: material
    end type sem_snapshot_cond

    ! ce type doit correspondre au type surface_t de sem_input.h **a l'ordre pres**
    type, bind(c) :: sem_surfaces
       type(C_PTR) :: next
       integer(C_INT):: surface_list(1:40)
       integer(C_INT) :: surface_present
       integer(C_INT) :: surface_type
       integer(C_INT) :: surface_mat
       integer(C_INT) :: surface_whatbc
       real(C_DOUBLE), dimension(3) :: surface_L
       real(C_DOUBLE), dimension(3) :: surface_C
       real(C_DOUBLE) :: surface_f0
       integer(C_INT) :: surface_dim
       real(C_DOUBLE) :: surface_paravalue(1:100)
       type(C_PTR) :: surface_paramname
       integer(C_INT) :: surface_nparamvar
       integer(C_INT) :: surface_paramvar
       type(C_PTR) :: surface_source
       type(C_PTR) :: surface_funcx
       type(C_PTR) :: surface_funcy
       type(C_PTR) :: surface_funcz
       type(C_PTR) :: surface_funcxy
       type(C_PTR) :: surface_funcxz
       type(C_PTR) :: surface_funcyz
       type(C_PTR) :: surface_varia
       real(C_DOUBLE) :: amplitude
       real(C_DOUBLE) :: Rtau
       integer(C_INT) :: surface_space
       real(C_DOUBLE) :: surface_size
       type(C_PTR)    :: surface_name
       integer(C_INT) :: surface_wave
       real(C_DOUBLE) :: surface_Speed;
       real(C_DOUBLE), dimension(3) :: surface_dirU;
    end type sem_surfaces

    type, bind(c) :: sem_material
       integer(C_INT):: num
       integer(C_INT):: domain
       integer(C_INT):: deftype
       integer(C_INT):: defspatial
       !
       real(C_DOUBLE) :: rho
       real(C_DOUBLE) :: Vp
       real(C_DOUBLE) :: Vs
       real(C_DOUBLE) :: E
       real(C_DOUBLE) :: nu
       real(C_DOUBLE) :: lambda
       real(C_DOUBLE) :: kappa
       real(C_DOUBLE) :: mu
       !
       type(C_PTR)    :: filename0
       type(C_PTR)    :: filename1
       type(C_PTR)    :: filename2
       !
       real(C_DOUBLE) :: syld
       real(C_DOUBLE) :: ckin
       real(C_DOUBLE) :: kkin
       real(C_DOUBLE) :: rinf
       real(C_DOUBLE) :: biso
       !
       type(C_PTR)    :: next
    end type sem_material

    type, bind(c) :: sem_material_list
        integer(C_INT):: count
        type(C_PTR):: head
    end type sem_material_list

    interface
       subroutine read_sem_config(config, spec, err) bind(c)
           use iso_c_binding
           import :: sem_config
           type(sem_config), intent(in) :: config
           character(C_CHAR), dimension(*) :: spec
           integer(C_INT), intent(out) :: err
       end subroutine read_sem_config

       subroutine read_sem_materials(materials, spec, err) bind(c)
           use iso_c_binding
           import :: sem_material_list
           type(sem_material_list), intent(in) :: materials
           character(C_CHAR), dimension(*) :: spec
           integer(C_INT), intent(out) :: err
       end subroutine read_sem_materials

       subroutine dump_config(config) bind(c)
           use iso_c_binding
           import :: sem_config
           type(sem_config), intent(in) :: config
       end subroutine dump_config

    end interface

contains

    function fromcstr(cstr)
        use iso_c_binding
        use semdatafiles
        type(c_ptr), intent(in) :: cstr
        character(Len=MAX_FILE_SIZE) :: fromcstr
        character,pointer,dimension(:) :: ctemp
        integer, dimension(1) :: clen
        integer :: i
        clen(1) = strlen(cstr)
        call c_f_pointer(cstr, ctemp, clen)
        fromcstr=''
        do i=1,clen(1)
            fromcstr(i:i) = ctemp(i)
        end do
    end function fromcstr

end module sem_c_config

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
