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
    !use mCapteur
    use semdatafiles
    use mpi
    use msnapshots
    use sem_c_bindings
    use shape_lin
    use shape_quad
    use treceivers
    use sglobal_energy
    use snewmark
    use smidpoint
    use srungekutta

    implicit none

    type (domain), target  :: Tdomain
    integer :: ntime,i_snap, ierr
    integer :: isort
    character(len=MAX_FILE_SIZE) :: fnamef
    integer :: getpid, pid

    real(fpp), parameter :: max_time_left=900
    real(fpp) :: remaining_time

    integer :: display_iter !! Indique si on doit faire des sortie lors de cette iteration
    real(kind=4) :: t_fin, t_ini
    integer :: interrupt, rg, code, protection, n_it_max

    pid = getpid()
    !write(*,*) "SEM2D[", pid, "] : Demarrage."

    display_iter = 1

    call START_SEM(Tdomain)

    rg = Tdomain%Mpi_var%my_rank
    if (rg == 0) call create_sem_output_directories()
    call MPI_Barrier (Tdomain%communicateur, ierr)

    !lecture du fichier de donnee
    if (rg == 0) write (*,*) "Read input.spec"
    call read_input (Tdomain)

    !lecture du fichier de maillage unv avec conversion en fichier sem2D
    if (rg == 0) write (*,*) "Define mesh properties"
    call read_mesh_h5(Tdomain)

    ! mesh deformation (for testing purposes)
    !call rotate_mesh(Tdomain)
    !call random_mesh_deformation(Tdomain)

    if (rg == 0) write (*,*) "Checks the inputs and the mesh"
    call check_inputs_and_mesh (Tdomain)

    if (rg == 0) write (*,*) "Compute Gauss-Lobatto-Legendre weights and zeroes"
    call compute_GLL (Tdomain)

    if (rg == 0) write (*,*) "Define a global numbering for the collocation points"
    call global_numbering (Tdomain)

    if (rg == 0) write (*,'(x,a,i1,a)') "Computing shape",Tdomain%n_nodes," functions within thier derivatives"
    if  (Tdomain%n_nodes == 4) then
        call shape4(TDomain)   ! Linear interpolation
    else if (Tdomain%n_nodes == 8) then
        call shape8(TDomain)  ! Quadratic interpolation
    else
        if (rg == 0) write (*,*) " Bad number of nodes for hexaedral shape "
        stop
    endif

    if (rg == 0) write (*,*) " Compute Courant parameter"
    call compute_Courant (Tdomain)

    if (rg == 0) write (*,*) "Attribute Boundary Conditions and PML properties"
    call PML_definition (Tdomain)

    if (Tdomain%logicD%any_source) then
        if (rg == 0) write (*,*) "Computing point-source parameters and location"
        call SourcePosition(Tdomain)
    endif

    if (Tdomain%logicD%save_trace ) then
        if (rg == 0) write (*,*) "Computing receivers parameters and locations"
        call ReceiverPosition(Tdomain)
    endif

    if (rg ==0 .and. Tdomain%logicD%super_object) write(*,*) "Define Fault properties"
    if (Tdomain%logicD%super_object_local_present) then
        if (Tdomain%n_fault > 0) call define_fault_properties (Tdomain)
    endif

    if (rg == 0) write (*,*) " Allocate fields"
    call allocate_domain (Tdomain)

    if (rg == 0) write (*,*) "Compute tranfer quantities"
    call wall_transfer (Tdomain)

    if (rg == 0) write (*,*) " Compute mass matrix and internal forces coefficients"
    call define_arrays (Tdomain)

    ! initialisation des temps
    Tdomain%TimeD%rtime = 0
    Tdomain%TimeD%NtimeMin = 0
    ! Nombre d'iterations pour schemas en temps iteratifs
    if (Tdomain%type_timeInteg==TIME_INTEG_MIDPOINT) then
        n_it_max = 0
    elseif (Tdomain%type_timeInteg==TIME_INTEG_MIDPOINT_ITER) then
        n_it_max = 1
    endif

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
        if (rg==0) write(*,*) "Snapshot every ", Tdomain%timeD%nsnap, " iterations"
    endif


    i_snap =1;


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BOUCLE DE CALCUL EN TEMPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
            call Newmark (Tdomain)
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

            if (rg==0 .and. display_iter==1) write(*,*) "Snapshot:",isort," iteration=", ntime, " tps=", Tdomain%TimeD%rtime
            call save_field_h5(Tdomain, rg, isort)

            ! Sortie Energie totale du systeme
            !if(Tdomain%LogicD%CompEnerg) call global_energy_generalized(Tdomain)

        endif

        if (i_snap==0) then
            if (rg == 0) then
                write(78,*)isort,Tdomain%TimeD%rtime
                call semname_nb_proc(isort,fnamef)
                open (79,file = fnamef,status="replace",form="formatted")
                write(79,*) Tdomain%Mpi_var%n_proc
                close(79)
            endif
        endif


        ! sortie des  ...
        if (Tdomain%logicD%save_fault_trace.and.i_snap==0) call save_fault_trace (Tdomain, ntime)


        ! sauvegarde des vitesses
        if (Tdomain%logicD%save_trace) call save_trace(Tdomain, ntime)


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

    call END_SEM(Tdomain, ntime)
    call MPI_Finalize  (ierr)

    if (rg == 0) write (*,*) "Execution completed"
    if (rg == 0) close(78)
    if (rg == 0) call CPU_TIME( t_fin )
    if (rg == 0) write (*,*) "CPU Time for computation : ", t_fin - t_ini

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
