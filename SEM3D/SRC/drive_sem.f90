!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file drive_sem.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

module drive_sem

contains

subroutine sem(master_superviseur, communicateur, communicateur_global)
    use sdomain
    use mrenumber
    use mCapteur
    use semdatafiles
    use sem_mpi
    use msnapshots
    use mdefinitions, only : define_arrays
    use semconfig !< pour config C
    use sem_c_bindings
    use stat, only : stat_starttick, stat_stoptick, stat_init, stat_finalize, STAT_START

    implicit none

    ! Hors couplage on doit avoir -1 MPI_COMM_WORLD, MPI_COMM_WORLD
    type(MPI_Comm), intent(in) :: communicateur, communicateur_global
    integer, intent(in) :: master_superviseur
    integer :: ee, ii, j, k, n, ngll
    real(fpp) :: bnum
    type(domain) :: Tdomain
    integer :: rg, nb_procs, ntime
    integer :: isort, ierr
    real(kind=8), parameter :: max_time_left=900

    Tdomain%time_0 = MPI_Wtime();

!----------------------------------------------------------------------------------------------!
!------------------------------------    COMMUNICATORS  ---------------------------------------!
!----------------------------------------------------------------------------------------------!

    Tdomain%communicateur = communicateur
    Tdomain%communicateur_global = communicateur_global
    Tdomain%master_superviseur = master_superviseur

    call MPI_Comm_Rank (Tdomain%communicateur, rg, ierr)
    call MPI_Comm_Size (Tdomain%communicateur, nb_procs, ierr)

    Tdomain%rank = rg
    Tdomain%nb_procs = nb_procs

    call stat_init(communicateur, rg, nb_procs, .false.)
    call stat_starttick(STAT_START)

 !----------------------------------------------------------------------------------------------!
 !--------------------------------       SEM 3D - RUNNING     ----------------------------------!
 !----------------------------------------------------------------------------------------------!
    call INIT_MESSAGE(rg)
    call START_SEM(rg)

 !---------------------------------------------------------------------------------------------!
 !------------------------------    RUN PREPARATION : INPUT DATA, -----------------------------!
 !----------------------------     ELEMENTAL AND GLOBAL MACHINERY  ----------------------------!
 !---------------------------------------------------------------------------------------------!

    call RUN_PREPARED(Tdomain)
    call RUN_INIT_INTERACT(Tdomain,isort)

!---------------------------------------------------------------------------------------------!
!-------------------------------    TIME STEPPING : EVOLUTION     ----------------------------!
!---------------------------------------------------------------------------------------------!
    ngll = Tdomain%sdomdg%ngll
    do n = 0,Tdomain%n_elem-1
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do ii = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,:,:) = 0.d0
                    enddo
                enddo
            enddo
        enddo
    enddo
    call stat_stoptick(STAT_START)

    call TIME_STEPPING(Tdomain,isort,ntime)

 !---------------------------------------------------------------------------------------------!
 !-------------------------------      NORMAL  END OF THE RUN      ----------------------------!
 !---------------------------------------------------------------------------------------------!

    call END_SEM(Tdomain,ntime)
    call stat_finalize(Tdomain)

end subroutine sem
 !-----------------------------------------------------------------------------------
 !-----------------------------------------------------------------------------------
subroutine INIT_MESSAGE(rg)
    implicit none
    integer, intent(in)   :: rg

    if(rg == 0)then
        write(*,*)
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "**************************                          ************************"
        write(*,*) "**************************     SEM - 3D VERSION     ************************"
        write(*,*) "**************************                          ************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*)
    end if

end subroutine INIT_MESSAGE
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine RUN_PREPARED(Tdomain)
    use sdomain
    use mrenumber
    use mdefinitions
    use sdomain_alloc
    use attenuation
    use mCapteur
    use semdatafiles
    use sem_mpi
    use msnapshots
    use semconfig !< pour config C
    use sem_c_bindings
    use mdefinitions
    use mshape8
    use mshape27
    use mdoublecouple
    use msource_excit
    use mondelette
    use surface_input

    implicit none
    type(domain), intent(inout) :: Tdomain
    integer :: rg
    integer :: code, i

    rg = Tdomain%rank
    if(rg == 0) print*
    if(rg == 0) print*
    if(rg == 0) print*, "****************************************************************************"
    if(rg == 0) print*, "********************                                     *******************"
    if(rg == 0) print*, "********************      RUN PREPARATION : INPUT DATA,  *******************"
    if(rg == 0) print*, "********************     ELEMENTAL AND GLOBAL MACHINERY  *******************"
    if(rg == 0) print*, "********************                                     *******************"
    if(rg == 0) print*, "****************************************************************************"
    if(rg == 0) print*

    if(rg == 0) print*, "--> CREATING OUTPUT DIRECTORIES"
    if(rg == 0) call create_sem_output_directories()
    ! Read_input peut avoir besoin du repertoire de logs
    call MPI_Barrier(Tdomain%communicateur, code)

 !- reading external data: run parameters (geometry, materials, time evolution,..)
    if(rg == 0) write(*,*) "--> READING INPUT PARAMETERS AND DATA"
    call read_input(Tdomain, code)

    !- eventual plane wave (transmission process: Bielak & Cristiano 1984)
    ! XXX

 !- eventual Neumann boundary conditions
    !if (Tdomain%logicD%neumann_local_present) then
    !    if (rg == 0) write(*,*) "--> DEFINING NEUMANN PROPERTIES"
    !    call define_Neumann_properties(Tdomain)
    !endif
    if (Tdomain%logicD%surfBC) then
       if (rg == 0) write(*,*) "--> DEFINING SURFACE PROPERTIES"
       call define_surface_properties(Tdomain)
    endif

    call MPI_Barrier(Tdomain%communicateur, code)

 !- from elementary to global numbering
    if (rg == 0) write (*,*) "--> DEFINING A GLOBAL NUMBERING FOR COLLOCATION POINTS"
    call global_numbering (Tdomain)
    call MPI_Barrier(Tdomain%communicateur,code)

 !- allocation of different fields' sizes
    if (rg == 0) write (*,*) "--> ALLOCATING FIELDS"
    call allocate_domain(Tdomain)
    call MPI_Barrier(Tdomain%communicateur,code)

 !- geometrical properties for integrals' calculations
    if (rg == 0) write (*,*) "--> COMPUTING SHAPE FUNCTIONS AND THEIR DERIVATIVES"
    if (Tdomain%n_nodes == 8) then
        ! Linear interpolation
        call shape8_init(Tdomain)
    else if (Tdomain%n_nodes == 27) then
        ! Quadratic interpolation
        call shape27_init(TDomain)
    else
        write (*,*) Tdomain%n_nodes, "control points not yet implemented in the code. Wait for an upgrade"
        stop
    endif
    call MPI_Barrier(Tdomain%communicateur,code)

    !- elementary properties (mass matrices, PML factors,..) geometry
    if (rg == 0) write (*,*) "--> COMPUTING MASS MATRIX AND INTERNAL FORCES COEFFICIENTS "
    call define_arrays(Tdomain)
    call MPI_Barrier(Tdomain%communicateur,code)
    ! Need to call compute_courant first
    call check_interface_orient(Tdomain, Tdomain%intSolPml, SMALLFPP)
    call check_interface_orient(Tdomain, Tdomain%intFluPml, SMALLFPP)
    call check_interface_orient(Tdomain, Tdomain%SF%intSolFlu, SMALLFPP)
    call check_interface_orient(Tdomain, Tdomain%SF%intSolFluPml, SMALLFPP)

 !- anelastic properties
    if (Tdomain%n_sls>0) then
        if (Tdomain%aniso) then
            if (rg == 0) write (*,*) "--> COMPUTING ANISOTROPIC ATTENUATION FEATURES"
            call set_attenuation_aniso_param(Tdomain)
        else
            if (rg == 0) write (*,*) "--> COMPUTING ATTENUATION FEATURES"
            call set_attenuation_param(Tdomain)
        endif
    endif
 !- eventual classical seismic point sources: their spatial and temporal properties
    if (Tdomain%logicD%any_source) then
        if (rg == 0) write (*,*) "--> COMPUTING SOURCE PARAMETERS "
        call SourcePosition (Tdomain)
        call double_couple (Tdomain, rg)
        call source_excit(Tdomain,rg)
        call def_timefunc (Tdomain, rg)
        !- source time dependence read in a file: Modules/Source.f90
        !  valid only for one point source - to be generalized for each
        do i = 0,Tdomain%n_source-1
            if(Tdomain%sSource(i)%i_time_function == 5)then
                call read_source_file(Tdomain%sSource(i))
            elseif (Tdomain%sSource(i)%i_time_function == 15) then
                call read_source_file_h5(Tdomain%sSource(i))
            endif
        end do
    endif
 !- time: initializations. Eventual changes if restarting from a checkpoint (see
 !         routine  RUN_INIT_INTERACT)
    Tdomain%TimeD%rtime = 0
    Tdomain%TimeD%NtimeMin = 0
    call MPI_Barrier(Tdomain%communicateur,code)

end subroutine RUN_PREPARED
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine RUN_INIT_INTERACT(Tdomain,isort)
    use sdomain
    use mCapteur
    use semdatafiles
    use sem_mpi
    use mloadcheckpoint
    use msnapshots
    use semconfig !< pour config C
    use sem_c_bindings
    use smirror

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer                     :: rg, nb_procs
    integer, intent(inout)      :: isort

    integer :: code,ierr
    integer :: info_capteur

    rg = Tdomain%rank
    nb_procs = Tdomain%nb_procs

    if(rg == 0) print*
    if(rg == 0) print*
    if(rg == 0) print*, "****************************************************************************"
    if(rg == 0) print*, "*****************                                          *****************"
    if(rg == 0) print*, "*****************  INITIALIZATION OF IN/OUT INTERACTIONS:  *****************"
    if(rg == 0) print*, "*****************     RESTART, SNAPSHOTS, RECEIVERS        *****************"
    if(rg == 0) print*, "*****************                                          *****************"
    if(rg == 0) print*, "****************************************************************************"
    if(rg == 0) print*

    isort = 1   ! index for snapshots outputs


!- eventual restarting from a previous run (checkpoint)
    if (Tdomain%logicD%run_restart) then
        !! Il faudra ajouter la gravite ici #ifdef COUPLAGE
        call read_restart(Tdomain, rg, isort)
        call MPI_Barrier(Tdomain%communicateur,code)
        if(rg == 0) then
            write (*,*) "--> RESTARTING ON ALL CPUs"
        end if
    else
        ! on supprime tous les fichiers et repertoire de protection
        if(rg == 0) then
            ! Sauvegarde des donnees de post-traitement
            call system('rm -Rf '//path_prot)
            ierr = sem_mkdir(path_prot)
        end if
        Tdomain%TimeD%prot_m2 = -1
        Tdomain%TimeD%prot_m1 = -1
        Tdomain%TimeD%prot_m0 = -1
    endif
    call MPI_Barrier(Tdomain%communicateur,code)

!- snapshots
    if (Tdomain%logicD%save_snapshots)  then
        call write_snapshot_geom(Tdomain, Tdomain%SnapData)
        Tdomain%timeD%nsnap = int(Tdomain%TimeD%time_snapshots / Tdomain%TimeD%dtmin)
        Tdomain%timeD%nsnap = max(1, Tdomain%timeD%nsnap)
        if(rg == 0) write (*,*) "--> SNAPSHOTS RECORDED EVERY ", Tdomain%timeD%nsnap, " iterations"
    end if

!- eventual outputs at receivers, and their properties
    info_capteur = 0
    if(Tdomain%logicD%save_trace) call create_capteurs(Tdomain)

    if(info_capteur /= 0) Tdomain%logicD%save_trace = .false.

    call MPI_Barrier(Tdomain%communicateur,code)

end subroutine RUN_INIT_INTERACT
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine TIME_STEPPING(Tdomain,isort,ntime)
    use sdomain
    use scomm
    use mCapteur
    use semdatafiles
    use sem_mpi
    use msnapshots
    use msavecheckpoint
    use mtimestep
    use semconfig !< pour config C
    use sem_c_bindings
    use stat, only : stat_starttick, stat_stoptick, STAT_TSTEP, STAT_IO, STAT_ITER
    use dom_solid, only : start_domain_solid, stop_domain_solid
    use dom_solidpml, only : start_domain_solidpml, stop_domain_solidpml
    use dom_fluid, only : start_domain_fluid, stop_domain_fluid
    use dom_fluidpml, only : start_domain_fluidpml, stop_domain_fluidpml
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(inout)      :: isort
    integer, intent(out)        :: ntime

    integer :: rg, nb_procs
    integer :: i_snap
    real(kind=8), parameter :: max_time_left = 900
    integer :: protection
    integer :: interrupt
    logical :: sortie_capteur
    double precision :: time_now
    integer :: time_now_s
    real(kind=8) :: remaining_time
    integer :: code, n

    rg = Tdomain%rank
    nb_procs = Tdomain%nb_procs

    if(rg == 0)then
        print*
        print*
        print*, "****************************************************************************"
        print*, "************************                                 *******************"
        print*, "************************        TIME STEPPING            *******************"
        print*, "************************                                 *******************"
        print*, "****************************************************************************"
        print*
    end if

    if(rg == 0)then
        print*,"--> Duration of the run: ",Tdomain%TimeD%Duration
        print*,"--> Time step size: ",Tdomain%TimeD%dtmin
        print*,"--> Number of time steps: ",Tdomain%TimeD%ntimeMax
        print*
    end if
    Tdomain%sdom%dt = Tdomain%TimeD%dtmin
!- snapshots counters
!   (isort already defined for snapshots outputting index)
    i_snap = 1

!- end of run flag
    interrupt = 0

    startup_init_done = .true.
    !---------------------------------------------------------!
    !--------------------  LOOP UPON TIME  -------------------!
    !---------------------------------------------------------!
    !$acc data &
    !$acc& copyin(Tdomain) &
    !$acc& copyin(Tdomain%sSource) &
    !$acc& copyin(Tdomain%Comm_data) &
    !$acc& copyin(Tdomain%Comm_data%Data) &
    !$acc&
    do n = 0,Tdomain%Comm_data%ncomm-1
        !$acc  enter data create(Tdomain%Comm_data%Data(n)%Give) &
        !$acc&            create(Tdomain%Comm_data%Data(n)%Take) &
        !$acc&      copyin(Tdomain%Comm_data%Data(n)) &
        !$acc&      copyin(Tdomain%Comm_data%Data(n)%IGiveF) &
        !$acc&      copyin(Tdomain%Comm_data%Data(n)%IGiveS) &
        !$acc&      copyin(Tdomain%Comm_data%Data(n)%IGiveSDG) &
        !$acc&      copyin(Tdomain%Comm_data%Data(n)%IGiveSPML) &
        !$acc&      copyin(Tdomain%Comm_data%Data(n)%IGiveFPML) &
        !$acc&
    end do
    call start_domain_solid(Tdomain, Tdomain%sdom)
    call start_domain_solidpml(Tdomain, Tdomain%spmldom)
    call start_domain_fluid(Tdomain, Tdomain%fdom)
    call start_domain_fluidpml(Tdomain, Tdomain%fpmldom)
    do ntime = Tdomain%TimeD%NtimeMin, Tdomain%TimeD%NtimeMax-1

        protection = 0
        if (interrupt > 0) then
            if(rg == 0) write(*,*) "--> RUN INTERRUPTION."
            exit
        end if

        !---------------------------------------------------------!
        !- TIME STEPPER TO CHOSE
        !---------------------------------------------------------!
        !- Newmark reduced to leap-frog
        call stat_starttick(STAT_ITER)
        select case(Tdomain%TimeD%type_timeinteg)
        case (TIME_INTEG_NEWMARK)
            call Newmark(Tdomain, ntime)
        case (TIME_INTEG_LDDRK64)
            call Timestep_LDDRK(Tdomain, ntime)
        end select
        call stat_stoptick(STAT_ITER)
        if (Tdomain%rank==0 .and. mod(ntime,20)==0) then
            print *,' Iteration  =  ',ntime,'    temps  = ',Tdomain%TimeD%rtime
        end if

        !---------------------------------------------------------!
        !- logical end of run
        !---------------------------------------------------------!
        if(ntime == Tdomain%TimeD%NtimeMax-1)then
            interrupt = 1
        endif



        !---------------------------------------------------------!
        !- OUTPUTS AND INTERACTIONS SEM -> OTHERS
        !---------------------------------------------------------!
        !- Checkpoint restart
        if(Tdomain%logicD%save_restart)then
            ! protection for a future restart?
            if(mod(ntime,Tdomain%TimeD%ncheck) == 0) then
                if (ntime/=0) protection = 1
            end if
            if(Tdomain%prot_at_time >0 .and. (Tdomain%first_prot)) then
                time_now = MPI_Wtime();
                time_now_s = int(time_now-Tdomain%time_0) 
                call MPI_ALLREDUCE(time_now_s, Tdomain%time_now_s, 1, MPI_INTEGER, MPI_MAX, Tdomain%communicateur_global, code)
                if(Tdomain%time_now_s < Tdomain%prot_at_time) then
                    protection = 0
                else
                    if(Tdomain%first_prot) then
                        protection = 1
                        Tdomain%first_prot = .false.
                        if(Tdomain%rank == 0) then
                             print*," MAKING FIRST PROT_AT_TIME:", ntime
                             print*, '(rm -r prot/*)'
                             call system("(rm -r prot/*)")
                        end if
                    end if
                end if
            end if
        endif
        !- Time remaining
        if (mod(ntime,20)==0) then
            call stat_starttick(STAT_TSTEP)
            call TREMAIN(remaining_time)
            if(remaining_time < max_time_left) interrupt = 1
            call MPI_ALLREDUCE(MPI_IN_PLACE, interrupt, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur_global, code)
            call stat_stoptick(STAT_TSTEP)
        end if

        !- snapshotting
        if(Tdomain%logicD%save_snapshots) i_snap = mod(ntime, Tdomain%TimeD%nsnap)

        !- Ici, on a une info globale pour interrupt, protection, i_snap
        if(interrupt > 0) protection = 1

        call stat_starttick(STAT_IO)

        !---------------------------------------------------------!
        !- SNAPSHOTS
        !---------------------------------------------------------!
        if(i_snap == 0 .and. Tdomain%logicD%save_snapshots) then
            !$acc update host(Tdomain%sdom%champs(0)%Depla, Tdomain%sdom%champs(0)%Veloc, Tdomain%sdom%champs(1)%Veloc) async(2) wait(1)
            !$acc wait(2)
            call OUTPUT_SNAPSHOTS(Tdomain,ntime,isort)
        end if
        !---------------------------------------------------------!
        !- RECEIVERS'OUTPUTS
        !---------------------------------------------------------!
        call evalueSortieCapteur(ntime, sortie_capteur)

        ! sortie des quantites demandees par les capteur
        if (sortie_capteur) then
!!            !$acc update host(Tdomain%sdom%champs(0)%Depla, Tdomain%sdom%champs(0)%Veloc, Tdomain%sdom%champs(1)%Veloc) async(2) wait(1)
!!            !$acc wait(2)
            call save_capteur(Tdomain, ntime)
        end if

        !---------------------------------------------------------!
        !- SAVE TO EVENTUAL RESTART
        !---------------------------------------------------------!
        if (protection /= 0 .and. Tdomain%logicD%save_restart) then
            if(Tdomain%rank == 0) print*," SAVING PROT"
            !$acc update host(Tdomain%sdom%champs(0)%Depla, Tdomain%sdom%champs(0)%Veloc, Tdomain%sdom%champs(1)%Veloc) async(2) wait(1)
            !$acc wait(2)

            call flushAllCapteurs(Tdomain)
            call save_checkpoint(Tdomain, Tdomain%TimeD%rtime, ntime, Tdomain%TimeD%dtmin, isort)
        endif
        call stat_stoptick(STAT_IO)

        !---------------------------------------------------------!
        !-  STOPPING RUN on all procs
        !---------------------------------------------------------!
        if(interrupt /= 0) then
            if(rg == 0) print*,"--> SEM INTERRUPTED, iteration= ", ntime
            exit
        endif

        !---------------------------------------------------------!
        !----      END OF LOOP: TIME INCREASING
        !---------------------------------------------------------!
        Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin
    enddo
    call stop_domain_solid(Tdomain, Tdomain%sdom)
    call stop_domain_solidpml(Tdomain, Tdomain%spmldom)
    call stop_domain_fluid(Tdomain, Tdomain%fdom)
    call stop_domain_fluidpml(Tdomain, Tdomain%fpmldom)
    !$acc end data

end subroutine TIME_STEPPING
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine OUTPUT_SNAPSHOTS(Tdomain,ntime,isort)
    use sdomain
    use semdatafiles
    use msnapshots

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: ntime
    integer, intent(inout)      :: isort
    !
    integer :: rg

    rg = Tdomain%rank
    if(rg == 0)then
        write(*,'(a34,i6.6,a8,f11.5)') "--> SEM : snapshot at iteration : ", ntime, " ,time: ", Tdomain%TimeD%rtime
    endif
    call save_field_h5(Tdomain, isort, Tdomain%SnapData)
    isort = isort + 1  ! a faire avant le save_checkpoint

end subroutine OUTPUT_SNAPSHOTS
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine END_SEM(Tdomain,ntime)
    use sdomain
    use mCapteur
    use semdatafiles
    use semconfig !< pour config C
    use sem_c_bindings
    use sem_mpi
    implicit none

    type(domain) :: Tdomain
    integer, intent(in)  :: ntime
    integer :: rg, ierr

    rg = Tdomain%rank
    if (rg==0) then
        open (111,file = "fin_sem", status="REPLACE")
        if (ntime >= Tdomain%TimeD%NtimeMax-1) then
            write(111,*) 1
        else
            write(111,*) 0
        end if
        close(111)
    end if


    call flushAllCapteurs(Tdomain)

    if (rg == 0) write (*,*) "--> DEALLOCATING DOMAIN."
    call deallocate_domain (Tdomain)

    if (Tdomain%logicD%save_snapshots)  then
        call MPI_Comm_free(Tdomain%SnapData%comm, ierr)
    end if
end subroutine END_SEM

subroutine START_SEM(rg)
    implicit none
    integer, intent(in) :: rg
    ! Ce fichier sert d'indicateur de fin de calcul
    ! Si en fin de run on trouve :
    !   -1 : il y a eu un crash/stop ou erreur avant la fin
    !    1 : le calcul s'est bien passe et est fini
    !    0 : le calcul doit repartir en reprise pour continuer
    !
    if (rg==0) then
        open (111,file = "fin_sem", status="REPLACE")
        write(111,*) -1
        close(111)
    end if

end subroutine START_SEM

end module drive_sem

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
