

module sem_c_config
    use iso_c_binding

    integer, parameter :: NAME_MAX=200
    ! Ce type doit correspondre au type sem_config_t du module read_input.c **a l'ordre pres**
    type, bind(c) :: sem_config
       type(C_PTR)    :: run_name

       !! Integration
       integer(C_INT) :: accel_scheme
       integer(C_INT) :: veloc_scheme
       real(C_DOUBLE) :: sim_time
       real(C_DOUBLE) :: alpha
       real(C_DOUBLE) :: beta
       real(C_DOUBLE) :: gamma
       real(C_DOUBLE) :: courant

       !! Modele, maillage
       type(C_PTR)    :: mesh_file
       integer(C_INT) :: model
       integer(C_INT) :: anisotropy
       type(C_PTR)    :: mat_file
       integer(C_INT) :: nsources
       type(C_PTR)    :: source

       !! Capteurs
       integer(C_INT) :: save_traces
       integer(C_INT) :: traces_interval
       integer(C_INT) :: traces_format
       type(C_PTR)    :: station_file

       !! Snapshots
       integer(C_INT) :: save_snap
       real(C_DOUBLE) :: snap_interval
       integer(C_INT) :: n_snap_cond
       type(C_PTR)    :: snapshot_selection;

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

       !! Neumann
       integer(C_INT) :: neu_present
       integer(C_INT) :: neu_type
       integer(C_INT) :: neu_mat
       real(C_DOUBLE), dimension(3) :: neu_L
       real(C_DOUBLE), dimension(3) :: neu_C
       real(C_DOUBLE) :: neu_f0

       !!Material
       integer(C_INT) :: material_present
       integer(C_INT) :: material_type
       type(C_PTR)    :: model_file
       real(C_DOUBLE) :: delta_lon
       real(C_DOUBLE) :: delta_lat

    end type sem_config


    ! Ce type doit correspondre au type source_t du module read_input.c **a l'ordre pres**
    type, bind(c) :: sem_source
       type(C_PTR) :: next
       real(C_DOUBLE), dimension(3) :: coords;
       integer(C_INT) :: type;
       integer(C_INT) :: dir;
       integer(C_INT) :: func;
       real(C_DOUBLE), dimension(6) :: moments;
       real(C_DOUBLE) :: tau;
       real(C_DOUBLE) :: freq;
       real(C_DOUBLE), dimension(4) :: band;
       real(C_DOUBLE) :: ts
       real(C_DOUBLE) :: gamma
       real(C_DOUBLE) :: amplitude
       type(C_PTR) :: time_file
    end type sem_source

    type, bind(c) :: sem_snapshot_cond
       type(C_PTR) :: next
       integer(C_INT) :: type
       integer(C_INT) :: include
       real(C_DOUBLE), dimension(6) :: box
       integer(C_INT) :: material
    end type sem_snapshot_cond

    interface
       subroutine read_sem_config(config, spec, err) bind(c)
           use iso_c_binding
           import :: sem_config
           type(sem_config), intent(in) :: config
           character(C_CHAR), dimension(*) :: spec
           integer(C_INT), intent(out) :: err
       end subroutine read_sem_config

       subroutine dump_config(config) bind(c)
           use iso_c_binding
           import :: sem_config
           type(sem_config), intent(in) :: config
       end subroutine dump_config

       function strlen(s) bind(c)
           use iso_c_binding
           type(c_ptr), intent(in), value :: s
           integer(C_INT) :: strlen
       end function strlen
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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
