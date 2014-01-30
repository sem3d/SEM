!>
!!\file read_input.F90
!!\brief Contient la subroutine read_input().
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<


subroutine create_sem2d_sources(Tdomain, config)
    use sem_c_config
    use sdomain
    implicit none
    type(domain), intent(inout)  :: Tdomain
    type(sem_config), intent(in) :: config
    type(sem_source), pointer :: src
    integer :: nsrc

    Tdomain%n_source = config%nsources
    allocate (Tdomain%Ssource(0:Tdomain%n_source-1))

    call c_f_pointer(config%source, src)
    nsrc = 0
    do while(associated(src))
        Tdomain%Ssource(nsrc)%Xsource = src%coords(1)
        !Tdomain%Ssource(nsrc)%Ysource = src%coords(2)
        Tdomain%Ssource(nsrc)%Zsource = src%coords(3)
        Tdomain%Ssource(nsrc)%i_type_source = src%type
        ! Comportement temporel
        Tdomain%Ssource(nsrc)%i_time_function = src%func
        Tdomain%Ssource(nsrc)%cutoff_freq = src%freq ! func=2,4
        Tdomain%Ssource(nsrc)%tau_b = src%tau ! func=1,2,3,4,5
        !Tdomain%Ssource(nsrc)%fh = src%band  ! func=3
        !Tdomain%Ssource(nsrc)%gamma = src%gamma ! func=4
        !Tdomain%Ssource(nsrc)%ts = src%ts   ! func=4
        Tdomain%Ssource(nsrc)%amplitude = src%amplitude
        ! Comportement Spacial
        ! i_type_source==1
        ! DIR = Z OR y -> y
        if (src%dir==2) src%dir = 1
        Tdomain%Ssource(nsrc)%i_dir = src%dir
        ! i_type_source==2
        !Tdomain%Ssource(nsrc)%Moment(0,0) = src%moments(1)
        !Tdomain%Ssource(nsrc)%Moment(1,1) = src%moments(2)
        !Tdomain%Ssource(nsrc)%Moment(2,1) = src%moments(3)
        !Tdomain%Ssource(nsrc)%Moment(0,1) = src%moments(4)
        !Tdomain%Ssource(nsrc)%Moment(1,0) = src%moments(4)
        !Tdomain%Ssource(nsrc)%Moment(0,2) = src%moments(5)
        !Tdomain%Ssource(nsrc)%Moment(2,0) = src%moments(5)
        !Tdomain%Ssource(nsrc)%Moment(1,2) = src%moments(6)
        !Tdomain%Ssource(nsrc)%Moment(2,1) = src%moments(6)

        nsrc = nsrc + 1
        Tdomain%logicD%any_source = .true.
        Tdomain%openfilescapt = .false.
        call c_f_pointer(src%next, src)
    end do

end subroutine create_sem2d_sources

!>
!! \brief Assure la lecture des fichiers de données en entrée à partir du fichier Parametrage/sem/input.spec
!!
!! \param type (domain), intent (INOUT) Tdomain
!<

subroutine read_input (Tdomain)
    use sem_c_config
    use sdomain
    use semdatafiles
    ! Modified by Gaetano Festa 01/06/05
    ! Introducing parallel features 12/10/05

    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables
    integer :: i, n_aus
    logical :: logic_scheme

    character(Len=MAX_FILE_SIZE) :: fnamef
    type(sem_config) :: config
    integer :: code
    ! #######################################

    ! It is done by any processor
    call semname_file_input_spec(fnamef)

    write(*,*) "Opening:", trim(fnamef)
    call read_sem_config(config, trim(fnamef)//C_NULL_CHAR, code)

    Tdomain%Title_simulation = fromcstr(config%run_name)
    Tdomain%Type_TimeInteg = config%type_timeinteg
    Tdomain%TimeD%acceleration_scheme = config%accel_scheme .ne. 0
    Tdomain%TimeD%velocity_scheme = config%veloc_scheme .ne. 0
    Tdomain%TimeD%duration = config%sim_time
    Tdomain%TimeD%alpha = config%alpha
    Tdomain%TimeD%beta = config%beta
    Tdomain%TimeD%gamma = config%gamma
    Tdomain%TimeD%courant = config%courant
    Tdomain%type_elem = config%type_elem
    Tdomain%type_flux = config%type_flux
    Tdomain%type_bc   = config%type_bc
    Tdomain%mesh_file = fromcstr(config%mesh_file)
    Tdomain%material_file = fromcstr(config%mat_file)
    Tdomain%logicD%save_trace = config%save_traces .ne. 0
    Tdomain%logicD%save_snapshots = config%save_snap .ne. 0
    Tdomain%logicD%run_restart = config%prorep .ne. 0
    Tdomain%TimeD%iter_reprise = config%prorep_iter
    Tdomain%TimeD%ncheck = config%prorep_iter ! frequence de sauvegarde
    Tdomain%station_file = fromcstr(config%station_file)
    !Tdomain%TimeD%ntrace = config%traces_interval ! XXX
    Tdomain%TimeD%time_snapshots = config%snap_interval
    logic_scheme = Tdomain%TimeD%acceleration_scheme .neqv. Tdomain%TimeD%velocity_scheme
    if(.not. logic_scheme) then
        stop "Both acceleration and velocity schemes: no compatibility, chose only one."
    end if


    call create_sem2d_sources(Tdomain, config)

    Tdomain%logicD%run_echo = .false.
    Tdomain%logicD%super_object = .false.
    Tdomain%logicD%super_object_local_present = .false.
    !read (11,*) Tdomain%logicD%save_deformation
    !read (11,*) Tdomain%logicD%save_energy

    !read (11,*) Tdomain%logicD%run_exec
    !read (11,*) Tdomain%logicD%run_debug
    !read (11,*) Tdomain%logicD%run_echo
    !if (Tdomain%logicD%save_trace) then
    !    read (11,*) Tdomain%station_file
    !else
    !    read( 11,*)
    !endif

    !if (Tdomain%logicD%save_snapshots .or.Tdomain%logicD%save_deformation) then
    !    read (11,*) Tdomain%TimeD%time_snapshots
    !else
    !    read( 11,*)
    !endif

    !logic_scheme = (Tdomain%TimeD%acceleration_scheme .or. Tdomain%TimeD%velocity_scheme) .and. &
    !    ((.not. Tdomain%TimeD%acceleration_scheme) .or. ( .not. Tdomain%TimeD%velocity_scheme))
    !
    !if (.not. logic_scheme) then
    !    write (*,*) "No compatible acceleration and velocity schemes"
    !    stop
    !endif

!! TODO
    Tdomain%logicD%super_object = .false.
    !read (11,*) Tdomain%logicD%super_object
    !if (Tdomain%logicD%super_object) then
    !    read (11,*) Tdomain%n_super_object
    !    allocate (Tdomain%Super_object_Type(0:Tdomain%n_super_object-1))
    !    allocate (Tdomain%Super_object_File(0:Tdomain%n_super_object-1))
    !    do i = 0, Tdomain%n_super_object-1
    !        read (11,*) Tdomain%Super_object_type(i), Tdomain%super_object_file(i)
    !    enddo
    !endif

    Tdomain%bMailUnv = .false.
    Tdomain%bCapteur = .false.
    ! conversion dun maillage unv en maillage sem
    !read (11,*) Tdomain%bMailUnv

    ! prise en compte du fichier des capteurs
    !read (11,*) Tdomain%bCapteur

    ! demarrage d un cas de reprise
    !read (11,*) Tdomain%logicD%run_restart

    ! numero de l iteration de reprise
    !read (11,*) Tdomain%TimeD%iter_reprise
    !read (11,*)

    ! creation de fichiers de reprise
    !read (11,*) Tdomain%logicD%save_restart
    !if (Tdomain%logicD%save_restart) then
    !    read (11,*) Tdomain%TimeD%ncheck ! frequence de sauvegarde
    !endif
    !
    !close (11)

    ! If echo modality write the read parameter in a file

    if (Tdomain%logicD%run_echo .and. Tdomain%Mpi_var%my_rank == 0) then

        call semname_read_input_spec(fnamef)
        open(91,file=fnamef, form="formatted", status="unknown")

        write (91,*) Tdomain%Title_simulation
        write (91,*) Tdomain%TimeD%acceleration_scheme
        write (91,*) Tdomain%TimeD%velocity_scheme
        write (91,*) Tdomain%TimeD%duration
        write (91,*) Tdomain%TimeD%alpha
        write (91,*) Tdomain%TimeD%beta
        write (91,*) Tdomain%TimeD%gamma
        write (91,*) Tdomain%mesh_file
        write (91,*) Tdomain%material_file
        write (91,*) Tdomain%logicD%save_trace
        write (91,*) Tdomain%logicD%save_snapshots
        write (91,*) Tdomain%logicD%save_energy
        write (91,*) Tdomain%logicD%plot_grid
        write (91,*) Tdomain%logicD%run_exec
        write (91,*) Tdomain%logicD%run_debug
        write (91,*) Tdomain%logicD%run_echo

        if (Tdomain%logicD%save_trace) then
            write (91,*) Tdomain%station_file
        else
            write (91,*) " No parameter need here"
        endif

        if (Tdomain%logicD%save_snapshots) then
            write (91,*) Tdomain%TimeD%time_snapshots
        else
            write (91,*) " No parameter need here"
        endif

        if (Tdomain%logicD%super_object) then
            write (91,*) Tdomain%n_super_object
            do i = 0, Tdomain%n_super_object-1
                write (91,*) Tdomain%Super_object_type(i),"   ", Tdomain%super_object_file(i)
            enddo
        else
            write (91,*) "no super objects present"
        endif

        if (Tdomain%logicD%any_source) then
            write (91,*) Tdomain%n_source
            do i = 0, Tdomain%n_source - 1
                write (91,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Zsource
                write (91,*) Tdomain%Ssource(i)%i_type_source
                if (Tdomain%Ssource(i)%i_type_source == 1 ) then
                    write (91,*) Tdomain%Ssource(i)%i_dir
                else
                    write (91,*) "No parameter need here"
                endif
                write (91,*) Tdomain%Ssource(i)%i_time_function
                write (91,*) Tdomain%Ssource(i)%tau_b
                write (91,*) Tdomain%Ssource(i)%cutoff_freq
            enddo
        else
            write (91,*)  "No available sources "
        endif
        write (91,*) "All right, runner ?"
        close (91)

    endif

    ! Define Super_Object properties: for the moment the only proper definition is the fault
    Tdomain%n_fault = 0 !Ajout Gsa  pour variables non initialisees??
    Tdomain%logicD%save_fault_trace = .false.

    if (Tdomain%logicD%super_object) then
        n_aus = 0
        do i = 0, Tdomain%n_super_object-1
            if (Tdomain%super_object_type(i) == "F") n_aus = n_aus + 1
        enddo
        Tdomain%n_fault = n_aus
        allocate (Tdomain%SFault(0:Tdomain%n_fault-1))
    endif


    return
end subroutine read_Input
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
