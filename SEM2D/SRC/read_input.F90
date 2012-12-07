!>
!!\file read_input.F90
!!\brief Contient la subroutine read_input().
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Assure la lecture des fichiers de données en entrée à partir du fichier Parametrage/sem/input.spec
!!
!! \param type (domain), intent (INOUT) Tdomain
!<


subroutine read_input (Tdomain)

    use sdomain
    use semdatafiles
    ! Modified by Gaetano Festa 01/06/05
    ! Introducing parallel features 12/10/05

    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables
    integer :: i, n_aus
    logical :: logic_scheme

    integer :: unit_src
    logical :: trouve_src
    character(Len=MAX_FILE_SIZE) :: fnamef
    ! #######################################

    ! It is done by any processor
    call semname_read_input_input(fnamef)

    open (11,file=fnamef,form="formatted",status="old")

    read (11,*) Tdomain%Title_simulation
    read (11,*) Tdomain%TimeD%acceleration_scheme
    read (11,*) Tdomain%TimeD%velocity_scheme
    read (11,*) Tdomain%TimeD%duration
    read (11,*) Tdomain%TimeD%alpha
    read (11,*) Tdomain%TimeD%beta
    read (11,*) Tdomain%TimeD%gamma
    read (11,*) Tdomain%mesh_file
    read (11,*) Tdomain%material_file
    read (11,*) Tdomain%logicD%save_trace
    read (11,*) Tdomain%logicD%save_snapshots
    read (11,*) Tdomain%logicD%save_deformation
    read (11,*) Tdomain%logicD%save_energy

    read (11,*) Tdomain%logicD%plot_grid
    read (11,*) Tdomain%logicD%run_exec
    read (11,*) Tdomain%logicD%run_debug
    read (11,*) Tdomain%logicD%run_echo
    if (Tdomain%logicD%save_trace) then
        read (11,*) Tdomain%station_file
    else
        read( 11,*)
    endif

    if (Tdomain%logicD%save_snapshots .or.Tdomain%logicD%save_deformation) then
        read (11,*) Tdomain%TimeD%time_snapshots
    else
        read( 11,*)
    endif

    logic_scheme = (Tdomain%TimeD%acceleration_scheme .or. Tdomain%TimeD%velocity_scheme) .and. &
        ((.not. Tdomain%TimeD%acceleration_scheme) .or. ( .not. Tdomain%TimeD%velocity_scheme))

    if (.not. logic_scheme) then
        write (*,*) "No compatible acceleration and velocity schemes"
        stop
    endif
    read (11,*) Tdomain%logicD%super_object
    if (Tdomain%logicD%super_object) then
        read (11,*) Tdomain%n_super_object
        allocate (Tdomain%Super_object_Type(0:Tdomain%n_super_object-1))
        allocate (Tdomain%Super_object_File(0:Tdomain%n_super_object-1))
        do i = 0, Tdomain%n_super_object-1
            read (11,*) Tdomain%Super_object_type(i), Tdomain%super_object_file(i)
        enddo
    endif

#ifdef MKA3D
    unit_src = 21
    Tdomain%logicD%any_source = .true.
    inquire(file="./Parametrage/sem/source.dat",exist=trouve_src)
    if( .not. trouve_src) then
        Tdomain%logicD%any_source=.false.
        write(*,*) 'No source for Sem'
    else
        call semname_read_input_source(fnamef)
        open (unit=unit_src,file=fnamef,form="formatted",status="old")
    endif
#else
    unit_src = 11
    read(unit_src,*) Tdomain%logicD%any_source
#endif

    !! read(11,*) Tdomain%logicD%any_source   !! version initiale
    if (Tdomain%logicD%any_source) then
        read (unit_src,*) Tdomain%n_source
        !! read (11,*) Tdomain%n_source !! version initiale
        allocate (Tdomain%Ssource(0:Tdomain%n_source-1))
        do i = 0, Tdomain%n_source - 1
            read (unit_src,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Zsource
            read (unit_src,*) Tdomain%Ssource(i)%i_type_source
            !!     read (11,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Zsource      !! version initiale
            !!     read (11,*) Tdomain%Ssource(i)%i_type_source                            !! version initiale
            if (Tdomain%Ssource(i)%i_type_source == 1 ) then
                !!        read (11,*) Tdomain%Ssource(i)%i_dir
                read (unit_src,*) Tdomain%Ssource(i)%i_dir
            else
                !!        read (11,*)
                read (unit_src,*)
            endif
            !!     read (11,*) Tdomain%Ssource(i)%i_time_function
            !!     read (11,*) Tdomain%Ssource(i)%tau_b
            !!     read (11,*) Tdomain%Ssource(i)%cutoff_freq
            !!     read (11,*) Tdomain%Ssource(i)%amplitude
            read (unit_src,*) Tdomain%Ssource(i)%i_time_function
            read (unit_src,*) Tdomain%Ssource(i)%tau_b
            read (unit_src,*) Tdomain%Ssource(i)%cutoff_freq
            read (unit_src,*) Tdomain%Ssource(i)%amplitude
        enddo
    endif

#ifdef MKA3D
    close(unit_src)
#endif


    ! conversion dun maillage unv en maillage sem
    read (11,*) Tdomain%bMailUnv

    ! prise en compte du fichier des capteurs
    read (11,*) Tdomain%bCapteur

    ! demarrage d un cas de reprise
    read (11,*) Tdomain%logicD%run_restart

    ! numero de l iteration de reprise
    read (11,*) Tdomain%TimeD%iter_reprise
    !read (11,*)

    ! creation de fichiers de reprise
    read (11,*) Tdomain%logicD%save_restart
    if (Tdomain%logicD%save_restart) then
        read (11,*) Tdomain%TimeD%ncheck ! frequence de sauvegarde
    endif

    close (11)

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
