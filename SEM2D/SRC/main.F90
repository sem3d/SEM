!>
!!\file main.F90
!!\brief Assure l'appel au code de calcul SEM2D.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<


program main
    use mpi
    use semdatafiles
    implicit none

    character(Len=MAX_FILE_SIZE),parameter :: p_param = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_traces = "./traces"
    character(Len=MAX_FILE_SIZE),parameter :: p_results = "./res"
    character(Len=MAX_FILE_SIZE),parameter :: p_data = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_prot = "./prot"
    character(Len=MAX_FILE_SIZE),parameter :: p_prop = "./prop"
    character(Len=MAX_FILE_SIZE),parameter :: p_prop_h5 = "./prop/h5"
    character(Len=MAX_FILE_SIZE),parameter :: p_mat = "./mat"
    character(Len=MAX_FILE_SIZE),parameter :: p_mirror = "./mirror"

    call init_sem_path(p_param, p_traces, p_results, p_data, p_prot, p_mat, p_mirror)

    call sem()

end program main

subroutine  sem()
    use sdomain
    use mCapteur
    use semdatafiles
    use mpi
    use msnapshots
    use sem_c_bindings
    use shape_lin
    use shape_quad
    use treceivers
    use sglobal_energy
    use snewmark
    use snewmark_modified
    use solid_fluid_coupling_2d
    use smidpoint
    use srungekutta
    use m_irons_dtcrit

    implicit none

    type (domain), target  :: Tdomain
    integer :: ntime,i_snap, ierr
    integer :: isort, nsrc
    character(len=MAX_FILE_SIZE) :: fnamef
    integer :: getpid, pid

    real(fpp), parameter :: max_time_left=900
    real(fpp) :: remaining_time

    integer :: display_iter !! Indique si on doit faire des sortie lors de cette iteration
    real(kind=4) :: t_fin, t_ini
    double precision :: t_wall_ini, t_wall_fin
    double precision :: my_cpu_time_d, my_cpu_sq
    double precision :: min_cpu, max_cpu, sum_cpu, sum_cpu_sq
    double precision :: avg_cpu, var_cpu, std_cpu
    integer :: interrupt, rg, code, protection, n_it_max

    pid = getpid()
    !write(*,*) "SEM2D[", pid, "] : Demarrage."

    display_iter = 1

    call START_SEM(Tdomain)

    call CPU_TIME(t_ini)
    t_wall_ini = MPI_Wtime()

    rg = Tdomain%Mpi_var%my_rank
    if (rg == 0) then
        write(*,*)
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "**************************                          ************************"
        write(*,*) "**************************     SEM - 2D VERSION     ************************"
        write(*,*) "**************************                          ************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*)
        
        print*
        print*
        print*, "****************************************************************************"
        print*, "********************                                     *******************"
        print*, "********************      RUN PREPARATION : INPUT DATA,  *******************"
        print*, "********************     ELEMENTAL AND GLOBAL MACHINERY  *******************"
        print*, "********************                                     *******************"
        print*, "****************************************************************************"
        print*
    end if

    if (rg == 0) call create_sem_output_directories()
    call MPI_Barrier (Tdomain%communicateur, ierr)

    !lecture du fichier de donnee
    if (rg == 0) write (*,*) "--> READING INPUT PARAMETERS AND DATA"
    call read_input (Tdomain)

    !lecture du fichier de maillage unv avec conversion en fichier sem2D
    if (rg == 0) write (*,*) "--> DEFINING MESH PROPERTIES"
    call read_mesh_h5(Tdomain)

    ! mesh deformation (for testing purposes)
    !call rotate_mesh(Tdomain)
    !call random_mesh_deformation(Tdomain)

    if (rg == 0) write (*,*) "--> CHECKING INPUTS AND MESH"
    call check_inputs_and_mesh (Tdomain)

    if (rg == 0) write (*,*) "--> COMPUTING GAUSS-LOBATTO-LEGENDRE WEIGHTS AND ZEROES"
    call compute_GLL (Tdomain)

    if (rg == 0) write (*,*) "--> SPLITTING SOLID-FLUID INTERFACE FACES (SEPARATE DOFS)"
    call split_sf_interface_faces (Tdomain)

    if (rg == 0) write (*,*) "--> DEFINING A GLOBAL NUMBERING FOR COLLOCATION POINTS"
    call global_numbering (Tdomain)

    if (rg == 0) write (*,'(a,i1,a)') "--> COMPUTING SHAPE ",Tdomain%n_nodes," FUNCTIONS AND THEIR DERIVATIVES"
    if  (Tdomain%n_nodes == 4) then
        call shape4(TDomain)   ! Linear interpolation
    else if (Tdomain%n_nodes == 8) then
        call shape8(TDomain)  ! Quadratic interpolation
    else
        if (rg == 0) write (*,*) " Bad number of nodes for hexaedral shape "
        stop
    endif

    if (rg == 0) write (*,*) "--> COMPUTING COURANT PARAMETERS"
    call compute_Courant (Tdomain)

    if (rg == 0) write (*,*) "--> DEFINING BOUNDARY CONDITIONS AND PML PROPERTIES"
    call PML_definition (Tdomain)
    call check_modified_newmark (Tdomain)

    if (Tdomain%logicD%any_source) then
        if (rg == 0) write (*,*) "--> COMPUTING POINT-SOURCE PARAMETERS AND LOCATION"
        call SourcePosition(Tdomain)
        ! source time dependence read from a file (func=file), cf. SEM3D drive_sem.f90
        do nsrc = 0, Tdomain%n_source-1
            if (Tdomain%sSource(nsrc)%i_time_function == 5) then
                if (rg == 0) write (*,*) "Reading source time file: ", &
                    trim(Tdomain%sSource(nsrc)%time_file)
                call read_source_file(Tdomain%sSource(nsrc))
            endif
        end do
    endif

    ! Legacy ASCII (.vel) receiver output DISABLED (kept for reference): traces now go through
    ! the station/capteurs h5/txt path (create_capteurs above, gated by save_traces). Flip the
    ! `.false.` here and at the save_trace call in the time loop to re-enable the .vel files.
    if (.false. .and. Tdomain%logicD%save_trace ) then
        if (rg == 0) write (*,*) "Computing receivers parameters and locations"
        call ReceiverPosition(Tdomain)
    endif

    if (rg ==0 .and. Tdomain%logicD%super_object) write(*,*) "--> DEFINING FAULT PROPERTIES"
    if (Tdomain%logicD%super_object_local_present) then
        if (Tdomain%n_fault > 0) call define_fault_properties (Tdomain)
    endif

    if (rg == 0) write (*,*) "--> ALLOCATING FIELDS"
    call allocate_domain (Tdomain)

    if (rg == 0) write (*,*) "--> COMPUTING WALL TRANSFER QUANTITIES"
    call wall_transfer (Tdomain)

    if (rg == 0) write (*,*) "--> COMPUTING MASS MATRIX AND INTERNAL FORCES COEFFICIENTS"
    call define_arrays (Tdomain)

    if (rg == 0) write (*,*) "--> COMPUTING RIGOROUS DT_CRIT (IRONS-TREHARNE EIGENVALUE BOUND)"
    call compute_irons_dtcrit (Tdomain)

    if (rg == 0) write (*,*) "--> BUILDING SOLID-FLUID INTERFACE COUPLING"
    call build_sf_coupling (Tdomain)

    ! initialisation des temps
    Tdomain%TimeD%rtime = 0
    Tdomain%TimeD%NtimeMin = 0

    ! Traces are gated by save_traces (as in SEM3D drive_sem.f90): if it is false, no
    ! receiver output is produced at all -- neither the legacy ASCII (.vel) below nor the
    ! station/capteurs h5/txt here. traces_format then selects the capteur format.
    Tdomain%has_station = .false.
    if (Tdomain%logicD%save_trace) call create_capteurs (Tdomain)
    ! Nombre d'iterations pour schemas en temps iteratifs
    if (Tdomain%type_timeInteg==TIME_INTEG_MIDPOINT) then
        n_it_max = 0
    elseif (Tdomain%type_timeInteg==TIME_INTEG_MIDPOINT_ITER) then
        n_it_max = 1
    endif

    if (rg == 0) then
        print*
        print*
        print*, "****************************************************************************"
        print*, "*****************                                          *****************"
        print*, "*****************  INITIALIZATION OF IN/OUT INTERACTIONS:  *****************"
        print*, "*****************     RESTART, SNAPSHOTS, RECEIVERS        *****************"
        print*, "*****************                                          *****************"
        print*, "****************************************************************************"
        print*
    end if

    isort = 1


    ! traitement de la reprise
    call semname_results_temps_sem(fnamef)
    if (Tdomain%logicD%run_restart ) then
        call read_restart(Tdomain,isort)
        open (78,file = fnamef,status="unknown",form="formatted",position="append")
    else ! pas de reprise
        if (rg == 0) then
            open (78,file = fnamef,status="unknown",form="formatted",position="rewind")
            ! on supprime tous les fichiers et repertoire de protection
            call system('rm -Rf ./ProRep/sem/Prot*')
        endif
        ! Ecriture de la geometrie des snapshots
        call write_snapshot_geom(Tdomain, rg)
    endif


    if (Tdomain%logicD%save_snapshots .or. Tdomain%logicD%save_deformation) then
        Tdomain%timeD%nsnap = int(Tdomain%TimeD%time_snapshots / Tdomain%TimeD%dtmin)
        if (Tdomain%timeD%nsnap == 0) Tdomain%timeD%nsnap = 1
        if (rg==0) write(*,*) "--> SNAPSHOTS RECORDED EVERY ", Tdomain%timeD%nsnap, " iterations"
    endif


    i_snap =1;


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BOUCLE DE CALCUL EN TEMPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (rg == 0) then
        print*
        print*
        print*, "****************************************************************************"
        print*, "************************                                 *******************"
        print*, "************************        TIME STEPPING            *******************"
        print*, "************************                                 *******************"
        print*, "****************************************************************************"
        print*
        print*,"--> Duration of the run: ",Tdomain%TimeD%Duration
        print*,"--> Time step size: ",Tdomain%TimeD%dtmin
        print*,"--> Number of time steps: ",Tdomain%TimeD%ntimeMax
        print*
    end if

    call CPU_TIME( t_ini )
    protection = 0
    interrupt = 0
    do ntime= Tdomain%TimeD%NtimeMin, Tdomain%TimeD%NtimeMax-1

        Tdomain%TimeD%ntime = ntime
        protection = 0
        if (interrupt>0) then
            if (rg==0) write(*,*) "Sortie sur limite de temps..."
            exit
        end if

        if (Tdomain%type_timeInteg==TIME_INTEG_NEWMARK) then
            if (Tdomain%TimeD%modified) then
                call NewmarkModified (Tdomain)
            else
                call Newmark (Tdomain)
            endif
        else if (Tdomain%type_timeInteg==TIME_INTEG_RK4) then
            call Runge_Kutta4(Tdomain, Tdomain%TimeD%dtmin)
        else if (Tdomain%type_timeInteg==TIME_INTEG_MIDPOINT .OR. &
                 Tdomain%type_timeInteg==TIME_INTEG_MIDPOINT_ITER) then
            if (Tdomain%Implicitness==TIME_INTEG_EXPLICIT) then
                call Midpoint_impl_expl(Tdomain, Tdomain%TimeD%dtmin,n_it_max)
                !call Midpoint_SEM (Tdomain)
            elseif (Tdomain%Implicitness==TIME_INTEG_SEMI_IMPLICIT) then
                !call Midpoint_test(Tdomain, Tdomain%TimeD%dtmin,n_it_max)
                call Midpoint_impl_semi_impl(Tdomain, Tdomain%TimeD%dtmin,n_it_max)
            endif
        endif

        if (rg == 0 .and. mod(ntime, 20) == 0) then
            print *, ' Iteration  =  ', ntime, '    temps  = ', Tdomain%TimeD%rtime
        end if

        if (ntime==Tdomain%TimeD%NtimeMax-1) then
            interrupt=1
        endif

        ! Checkpoint restart
        if (Tdomain%logicD%save_restart)  then
            !checkpoint a iteration ntime
            if (mod (ntime,Tdomain%TimeD%ncheck) == 0) then
                protection=1
            endif
        endif

        call tremain( remaining_time )
        !if (rg==0) write(*,*) "remain:", remaining_time
        if (remaining_time<max_time_left) then
            interrupt = 1
        end if
        call mpi_allreduce(MPI_IN_PLACE, interrupt, 1, MPI_INTEGER, MPI_SUM, &
            Tdomain%communicateur_global, code)

        if (Tdomain%logicD%save_snapshots)  then
            i_snap = mod (ntime, Tdomain%TimeD%nsnap)
        end if

        ! Ici, on a une info globale pour interrupt, protection, i_snap
        if (interrupt>0) then
            protection = 1
        end if

        if (mod(ntime,100)==0) then
            if(Tdomain%LogicD%CompEnerg) call global_energy_generalized(Tdomain)
        endif

        if (i_snap == 0) then

            if (rg==0 .and. display_iter==1) then
                write(*,'(a34,i6.6,a8,f11.5)') "--> SEM : snapshot at iteration : ", ntime, " ,time: ", Tdomain%TimeD%rtime
            endif
            call save_field_h5(Tdomain, rg, isort)

            ! Sortie Energie totale du systeme
            !if(Tdomain%LogicD%CompEnerg) call global_energy_generalized(Tdomain)

        endif

        if (i_snap==0) then
            if (rg == 0) then
                write(78,*)isort,Tdomain%TimeD%rtime
            endif
        endif


        ! sortie des  ...
        if (Tdomain%logicD%save_fault_trace.and.i_snap==0) call save_fault_trace (Tdomain, ntime)


        ! sauvegarde des vitesses -- legacy ASCII (.vel) writer DISABLED (kept for reference;
        ! see the ReceiverPosition guard above). h5/txt traces come from save_capteur below.
        if (.false. .and. Tdomain%logicD%save_trace) call save_trace(Tdomain, ntime)

        if (Tdomain%has_station) call save_capteur(Tdomain, ntime)


        if (i_snap==0) then
            isort=isort+1  ! a faire avant le save_checkpoint
        endif

        if (protection/=0) then
            call save_checkpoint(Tdomain,Tdomain%TimeD%rtime,Tdomain%TimeD%dtmin,ntime,isort)
        endif

        ! arret des calculs sur tous les procs
        if (interrupt/=0) then
            print*,"Arret de SEM, iteration=", ntime
            exit
        endif

        ! incrementation du pas de temps
        Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin

    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FIN BOUCLE DE CALCUL EN TEMPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (Tdomain%has_station) call flushAllCapteurs(Tdomain)

    call CPU_TIME(t_fin)
    t_wall_fin = MPI_Wtime()
    my_cpu_time_d = dble(t_fin - t_ini)
    my_cpu_sq = my_cpu_time_d * my_cpu_time_d

    call MPI_Reduce(my_cpu_time_d, min_cpu, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, Tdomain%communicateur, ierr)
    call MPI_Reduce(my_cpu_time_d, max_cpu, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, Tdomain%communicateur, ierr)
    call MPI_Reduce(my_cpu_time_d, sum_cpu, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, Tdomain%communicateur, ierr)
    call MPI_Reduce(my_cpu_sq, sum_cpu_sq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, Tdomain%communicateur, ierr)

    call END_SEM(Tdomain, ntime)
    call MPI_Finalize  (ierr)

    if (rg == 0) then
        avg_cpu = sum_cpu / dble(Tdomain%Mpi_var%n_proc)
        var_cpu = (sum_cpu_sq - (sum_cpu * sum_cpu) / dble(Tdomain%Mpi_var%n_proc)) / dble(Tdomain%Mpi_var%n_proc)
        if (var_cpu < 0.d0) var_cpu = 0.d0
        std_cpu = sqrt(var_cpu)

        write (*,*) "Execution completed"
        write (*,*) "Wall Time for computation : ", t_wall_fin - t_wall_ini
        write (*,*) "CPU Time stats across ", Tdomain%Mpi_var%n_proc, " process(es):"
        write (*,*) "  Min CPU Time : ", min_cpu
        write (*,*) "  Max CPU Time : ", max_cpu
        write (*,*) "  Avg CPU Time : ", avg_cpu
        write (*,*) "  Std Dev CPU  : ", std_cpu
        close(78)
    endif

end subroutine sem


subroutine START_SEM(Tdomain)
    use sdomain
    use mpi
    implicit none
    type(domain), intent(inout) :: Tdomain
    integer :: rg, ierr
    integer :: global_rank, global_nb_proc

    call MPI_Init(ierr)
    call MPI_Comm_Rank (MPI_COMM_WORLD, global_rank, ierr)
    call MPI_Comm_size (MPI_COMM_WORLD, global_nb_proc, ierr)

    Tdomain%Mpi_var%my_rank = global_rank
    Tdomain%Mpi_var%n_proc = global_nb_proc
    Tdomain%communicateur=MPI_COMM_WORLD
    Tdomain%communicateur_global=MPI_COMM_WORLD

    rg = Tdomain%Mpi_var%my_rank

    if (rg==0) then
        open (111,file = "fin_sem", status="REPLACE")
        write(111,*) -1
        close(111)
    end if

    Tdomain%TimeD%prot_m2 = -1
    Tdomain%TimeD%prot_m1 = -1
    Tdomain%TimeD%prot_m0 = -1

end subroutine START_SEM

subroutine END_SEM(Tdomain, ntime)
    use sdomain
    implicit none
    type(domain), intent(in) :: Tdomain
    integer, intent(in) :: ntime
    integer :: rg
    rg = Tdomain%Mpi_var%my_rank
    if (rg==0) then
        open (111,file = "fin_sem", status="REPLACE")
        if (ntime >= Tdomain%TimeD%NtimeMax-1) then
            write(111,*) 1
        else
            write(111,*) 0
        end if
        close(111)
    end if

end subroutine END_SEM


!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
