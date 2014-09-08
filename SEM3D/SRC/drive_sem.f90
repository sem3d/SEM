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
    use mdefinitions
    use mrenumber
    use mCapteur
    use semdatafiles
    use mpi
    use msnapshots
    use mdefinitions, only : define_arrays
    use semconfig !< pour config C
    use sem_c_bindings
#ifdef COUPLAGE
    use scouplage
#endif

    implicit none

    ! Hors couplage on doit avoir -1 MPI_COMM_WORLD, MPI_COMM_WORLD
    integer, intent(in) :: communicateur, communicateur_global, master_superviseur

    type(domain), target :: Tdomain

    integer :: code, rg, nb_procs, ntime, i_snap, n
    integer :: i, isort, ierr
    real(kind=8) :: remaining_time
    real(kind=8), parameter :: max_time_left=900
    integer :: protection
    character(Len=MAX_FILE_SIZE) :: fnamef
#ifdef COUPLAGE
    integer :: global_rank, global_nb_proc, worldgroup, intergroup
    integer :: m_localComm, comm_super_mka
    integer :: sortie
    integer :: tag
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer :: MaxNgParDir
    integer, dimension(3) :: flags_synchro ! fin/protection/sortie
    integer :: getpid, pid
    integer, dimension(3) :: tab
    integer :: min_rank_glob_sem
    integer :: master_sup, nb_procs_sup, master_mka, nb_procs_mka
#endif
    integer :: interrupt
    integer :: group, subgroup

    call MPI_Init (ierr)

!----------------------------------------------------------------------------------------------!
!------------------------------------    COMMUNICATORS  ---------------------------------------!
!----------------------------------------------------------------------------------------------!

#ifdef COUPLAGE
    pid = getpid()
    write(*,*) "SEM3D[", pid, "] : Demarrage."
    call MPI_Comm_Rank (MPI_COMM_WORLD, global_rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, global_nb_proc, ierr)
    call MPI_Comm_split(MPI_COMM_WORLD, 2, rg, m_localComm, ierr)
    call MPI_Comm_Rank (m_localComm, rg, ierr)
    call MPI_Comm_size(m_localComm, nb_procs, ierr)

    Tdomain%communicateur=m_localComm
    Tdomain%communicateur_global=MPI_COMM_WORLD
    Tdomain%master_superviseur=0

    !! Reception des infos du superviseur
    tag=8000000+global_rank
    call MPI_Recv(tab, 2, MPI_INTEGER, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, ierr)

    !! Reception des infos de mka
    tag=8100000+global_rank
    call MPI_Recv(tab, 2, MPI_INTEGER, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, ierr)

    !! On cherche le global rank minimal des procs sem
    call MPI_Reduce(global_rank, min_rank_glob_sem, 1, MPI_INTEGER, MPI_MIN, 0, m_localComm, ierr)

    !! Envoi des infos de couplage
    if (rg == 0) then
	tab(1) = global_rank
	tab(2) = nb_procs
        tab(3) = min_rank_glob_sem

	do i=1, global_nb_proc
	    tag=8200000+i
	    call MPI_Send(tab, 3, MPI_INTEGER, i-1, tag, MPI_COMM_WORLD, ierr)
	enddo
    endif

    !! Reception des infos de sem
    tag=8200000+global_rank+1
    call MPI_Recv(tab, 3, MPI_INTEGER, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, ierr)

    call MPI_Comm_group(MPI_COMM_WORLD, worldgroup, ierr)
    call MPI_Group_incl(worldgroup, 2, tab, intergroup, ierr)
    call MPI_Comm_create(MPI_COMM_WORLD, intergroup, comm_super_mka, ierr)
#else
    Tdomain%communicateur = communicateur
    Tdomain%communicateur_global = communicateur_global
    Tdomain%master_superviseur = master_superviseur

    call MPI_Comm_Rank (Tdomain%communicateur, rg, code)
    call MPI_Comm_Size (Tdomain%communicateur, nb_procs,  code)
#endif

!----------------------------------------------------------------------------------------------!
!--------------------------------       SEM 3D - RUNNING     ----------------------------------!
!----------------------------------------------------------------------------------------------!

    call INIT_MESSAGE(rg)

!---------------------------------------------------------------------------------------------!
!------------------------------    RUN PREPARATION : INPUT DATA, -----------------------------!
!----------------------------     ELEMENTAL AND GLOBAL MACHINERY  ----------------------------!
!---------------------------------------------------------------------------------------------!

    call RUN_PREPARED(Tdomain,rg)
    call RUN_INIT_INTERACT(Tdomain,rg,isort) 

!---------------------------------------------------------------------------------------------!
!-------------------------------    TIME STEPPING : EVOLUTION     ----------------------------!
!---------------------------------------------------------------------------------------------!

    call TIME_STEPPING(Tdomain,rg,isort,ntime)

!---------------------------------------------------------------------------------------------!
!-------------------------------      NORMAL  END OF THE RUN      ----------------------------!
!---------------------------------------------------------------------------------------------!

    call END_SEM(Tdomain,rg,ntime)

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
subroutine RUN_PREPARED(Tdomain,rg)

    use sdomain
    use mCapteur
    use semdatafiles
    use mpi
    use msnapshots
    use semconfig !< pour config C
    use sem_c_bindings
#ifdef COUPLAGE
    use scouplage
#endif

    implicit none
    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    integer :: code, i, ierr, group, subgroup

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
    call read_input(Tdomain, rg, code)

!- Create subdomains communicators
    group = rg/Tdomain%ngroup
    subgroup = mod(rg,Tdomain%ngroup)
    call MPI_Comm_split(MPI_COMM_WORLD, group, subgroup, Tdomain%comm_output, ierr)
    call MPI_Comm_size(Tdomain%comm_output, Tdomain%nb_output_procs,  code)
    call MPI_Comm_rank(Tdomain%comm_output, Tdomain%output_rank, code)


!- eventual plane wave (transmission process: Bielak & Cristiano 1984)
    if (Tdomain%logicD%super_object_local_present) then
        if (Tdomain%super_object_type == "P") then
            if(rg == 0) write(*,*) "--> DEFINING PLANE WAVE PROPERTIES"
            call define_planew_properties(Tdomain)
        endif
    endif

!- eventual Neumann boundary conditions
    if (Tdomain%logicD%neumann_local_present) then
        if (rg == 0) write(*,*) "--> DEFINING NEUMANN PROPERTIES"
        call define_Neumann_properties(Tdomain,rg)
    endif
    call MPI_Barrier(Tdomain%communicateur, code)

!- discretization (collocation) points' properties
    if (rg == 0) write (*,*) "--> COMPUTING GAUSS-LOBATTO-LEGENDRE PROPERTIES"
    call compute_GLL(Tdomain)
    call MPI_Barrier(Tdomain%communicateur, code)

!- absorbing layers (PMLs)
    if (Tdomain%any_PML)then
        if (rg == 0) write (*,*) "Attribute PML properties "
        call PML_definition(Tdomain)
    endif
    call MPI_BARRIER(Tdomain%communicateur,code)

!- from elementary to global numbering
    if (rg == 0) write (*,*) "--> DEFINING A GLOBAL NUMBERING FOR COLLOCATION POINTS"
    call global_numbering (Tdomain,rg)
    call MPI_Barrier(Tdomain%communicateur,code)

!- geometrical properties for integrals' calculations 
    if (rg == 0) write (*,*) "--> COMPUTING SHAPE FUNCTIONS AND THEIR DERIVATIVES"
    if (Tdomain%n_nodes == 8) then
        call shape8(Tdomain,rg)   ! Linear interpolation
    else if (Tdomain%n_nodes == 20) then
        write (*,*) "    -> 20 control points not yet implemented in the code. Wait for an upgrade"
        stop
    else
        call shape27(TDomain)   ! Quadratic interpolation
    endif
    call MPI_Barrier(Tdomain%communicateur,code)

!- timestep value - > Courant, or Courant -> timestep 
    if (rg == 0) write (*,*) "--> COMPUTING COURANT PARAMETER"
    call compute_Courant(Tdomain,rg)
    call MPI_Barrier(Tdomain%communicateur,code)

!- allocation of different fields' sizes 
    if (rg == 0) write (*,*) "--> ALLOCATING FIELDS"
    call allocate_domain(Tdomain, rg)
    call MPI_Barrier(Tdomain%communicateur,code)

!- elementary properties (mass matrices, PML factors,..)
    if (rg == 0) write (*,*) "--> COMPUTING MASS MATRIX AND INTERNAL FORCES COEFFICIENTS "
    call define_arrays(Tdomain, rg)
    call MPI_Barrier(Tdomain%communicateur,code)

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
        call SourcePosition (Tdomain, rg)
        call double_couple (Tdomain, rg)
        call source_excit(Tdomain,rg)
        call def_timefunc (Tdomain, rg)
        !- source time dependence read in a file: Modules/Source.f90
        !  valid only for one point source - to be generalized for each
        do i = 0,Tdomain%n_source-1
            if(Tdomain%sSource(i)%i_time_function == 5)then
                call read_source_file(Tdomain%sSource(i))
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
subroutine RUN_INIT_INTERACT(Tdomain,rg,isort) 
    use sdomain
    use mCapteur
    use semdatafiles
    use mpi
    use msnapshots
    use semconfig !< pour config C
    use sem_c_bindings
#ifdef COUPLAGE
    use scouplage
#endif

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    integer, intent(inout)      :: isort

    integer :: code,nb_procs, i_snap, ierr
    integer :: info_capteur
    character(Len=MAX_FILE_SIZE) :: fnamef

    
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

!- coupling with other codes: initialization (only Mka3d for the time being)
#ifdef COUPLAGE
    call INIT_COUPLING_MKA(Tdomain,nb_procs,rg)
#endif


!- eventual restarting from a previous run (checkpoint)
    call semname_results_temps_sem(fnamef)
    if (Tdomain%logicD%run_restart) then
        !! Il faudra ajouter la gravite ici #ifdef COUPLAGE
        call read_restart(Tdomain, rg, isort)
        call MPI_Barrier(Tdomain%communicateur,code)
        if(rg == 0) write (*,*) , "--> RESTARTING ON ALL CPUs"
        open(78,file=fnamef,status="unknown",position="append")
    else
        ! Sauvegarde des donnees de post-traitement
        open(78,file=fnamef,status="unknown",position="rewind")
        ! on supprime tous les fichiers et repertoire de protection
        if(rg == 0) then
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
        call write_snapshot_geom(Tdomain, rg)
        Tdomain%timeD%nsnap = int(Tdomain%TimeD%time_snapshots / Tdomain%TimeD%dtmin)
        Tdomain%timeD%nsnap = max(1, Tdomain%timeD%nsnap)
        if(rg == 0) write (*,*) , "--> SNAPSHOTS RECORDED EVERY ", Tdomain%timeD%nsnap, " iterations"
    end if

!- eventual outputs at receivers, and their properties
    info_capteur = 0
    if(Tdomain%logicD%save_trace) call read_capteur(Tdomain, rg, info_capteur)
    if(info_capteur /= 0) Tdomain%logicD%save_trace = .false.

end subroutine RUN_INIT_INTERACT
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine INIT_COUPLING_MKA(Tdomain,nb_procs,rg)

    use sdomain
    use mCapteur
    use semdatafiles
    use mpi
    use msnapshots
    use semconfig !< pour config C
    use sem_c_bindings
    use scouplage

    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)          :: nb_procs,rg
    character(Len=MAX_FILE_SIZE) :: fnamef
    integer                      :: finSem,MaxNgParDir

    if (rg == 0) write(6,'(A)') 'Methode de couplage Mka3D/SEM3D 4: Peigne de points d''interpolation, en vitesses, systeme lineaire sur la vitesse,', &
        ' contraintes discontinues pour les ddl de couplage'
    call semname_drive_sem_listing(rg, fnamef)
    open(UNIT=50, FILE=fnamef, STATUS='UNKNOWN')
#ifdef COUPLAGE
    call initialisation_couplage(Tdomain, rg, MaxNgParDir, nb_procs)
#endif
    finSem = 0

    !- new max. number of iterations 
    Tdomain%TimeD%ntimeMax = int(Tdomain%TimeD%Duration/Tdomain%TimeD%dtmin)
   ! this block placed here - is it ok, or not?
#ifdef COUPLAGE
    call reception_surface_part_mka(Tdomain, rg)
    call reception_nouveau_pdt_sem(Tdomain, rg)
#endif


end subroutine INIT_COUPLING_MKA
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine TIME_STEPPING(Tdomain,rg,isort,ntime)
    use sdomain
    use mCapteur
    use semdatafiles
    use mpi
    use msnapshots
    use semconfig !< pour config C
    use sem_c_bindings
#ifdef COUPLAGE
    use scouplage
#endif

    implicit none

    type(domain), target :: Tdomain
    integer, intent(inout)  :: rg,isort
    integer, intent(out)  :: ntime

    integer :: code,nb_procs,i_snap, n
    integer :: i
    real(kind=8) :: remaining_time
    real(kind=8), parameter :: max_time_left = 900
    integer :: protection
    character(Len=MAX_FILE_SIZE) :: fnamef
    integer :: interrupt
    logical :: sortie_capteur


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

!- snapshots counters
!   (isort already defined for snapshots outputting index)
    i_snap = 1

!- end of run flag
    interrupt = 0

!---------------------------------------------------------!
!--------------------  LOOP UPON TIME  -------------------!
!---------------------------------------------------------!
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
        call NEWMARK(Tdomain, rg, ntime)


!---------------------------------------------------------!
    !- logical end of run
!---------------------------------------------------------!
        if(ntime == Tdomain%TimeD%NtimeMax-1)then
            interrupt = 1
        endif


        
!---------------------------------------------------------!
    !- OUTPUTS AND INTERACTIONS SEM -> OTHERS
!---------------------------------------------------------!
#ifdef COUPLAGE
    !- CODE COUPLING: Mka3d (SEM -> Mka)
        call MKA_COUPLING_OUT(Tdomain,rg,ntime,interrupt,protection,i_snap)
#else
    !- Checkpoint restart
        if(Tdomain%logicD%save_restart)then
            ! protection for a future restart?
            if(mod(ntime,Tdomain%TimeD%ncheck) == 0) protection = 1
        endif
    !- Time remaining
        call TREMAIN(remaining_time)
        if(remaining_time < max_time_left) interrupt = 1
        call MPI_ALLREDUCE(MPI_IN_PLACE, interrupt, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur_global, code)

    !- snapshotting
        if(Tdomain%logicD%save_snapshots) i_snap = mod(ntime, Tdomain%TimeD%nsnap)
#endif

    !- Ici, on a une info globale pour interrupt, protection, i_snap
        if(interrupt > 0) protection = 1


!---------------------------------------------------------!
    !- SNAPSHOTS
!---------------------------------------------------------!
        if(i_snap == 0 .and. Tdomain%logicD%save_snapshots) &
            call OUTPUT_SNAPSHOTS(Tdomain,nb_procs,rg,ntime,isort)


!---------------------------------------------------------!
    !- RECEIVERS'OUTPUTS
!---------------------------------------------------------!
        call evalueSortieCapteur(ntime, Tdomain%TimeD%rtime, sortie_capteur)
        ! sortie des quantites demandees par les capteur
        if (sortie_capteur) call save_capteur(Tdomain, ntime, rg)


!---------------------------------------------------------!
    !- SAVE TO EVENTUAL RESTART
!---------------------------------------------------------!
        if(protection /= 0)then
            call flushAllCapteurs(Tdomain, rg)
            call save_checkpoint(Tdomain, Tdomain%TimeD%rtime, ntime, rg, Tdomain%TimeD%dtmin, isort)
        endif

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


end subroutine TIME_STEPPING
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine MKA_COUPLING_OUT(Tdomain,rg,ntime,interrupt,protection,i_snap)
    use sdomain
    use mCapteur
    use semdatafiles
    use mpi
    use msnapshots
    use semconfig !< pour config C
    use sem_c_bindings
    use scouplage

    implicit none

    type(domain), target :: Tdomain
    integer, intent(in)  :: rg,ntime
    integer, intent(inout) :: interrupt,protection,i_snap

    integer :: code,n
    integer, dimension(3) :: flags_synchro ! fin/protection/sortie

#ifdef COUPLAGE
    call ENVOI_VITESSE_MKA(Tdomain, ntime, rg) !! syst. lineaire vitesse
#endif

!- traitement de l arret
!  comm entre le master_sem et le Tdomain%master_superviseur : reception de la valeur de l interruption

    flags_synchro(1) = interrupt
    flags_synchro(2) = 0      ! SEM ne declenche pas les protections
    flags_synchro(3) = 0

    call MPI_ALLREDUCE(MPI_IN_PLACE, flags_synchro, 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, code)

    interrupt = flags_synchro(1)
    protection = flags_synchro(2)
    if(flags_synchro(3) /= 0) then
        i_snap = 0
    else
        i_snap = 1
    end if
!- reinitializations of coupling fields    
    ! remise a zero dans tous les cas des forces Mka sur les faces
#ifdef COUPLAGE
    do n = 0, Tdomain%n_face-1
        Tdomain%sFace(n)%ForcesMka = 0.
    enddo
    ! remise a zero dans tous les cas des forces Mka sur les vertex
    do n = 0, Tdomain%n_vertex-1
        Tdomain%sVertex(n)%ForcesMka = 0.
    enddo
    ! remise a zero dans tous les cas des forces Mka sur les Edges
    do n = 0, Tdomain%n_edge-1
        Tdomain%sEdge(n)%ForcesMka = 0.
    enddo
#endif

end subroutine MKA_COUPLING_OUT
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine OUTPUT_SNAPSHOTS(Tdomain,nb_procs,rg,ntime,isort)
    use sdomain
    use semdatafiles
    use msnapshots

    implicit none

    type(domain), target   :: Tdomain
    integer, intent(in)    :: rg,ntime,nb_procs
    integer, intent(inout) :: isort
    character(Len=MAX_FILE_SIZE) :: fnamef

    if(rg == 0)then
        write(*,'(a34,i6.6,a8,f11.5)') "--> SEM : snapshot at iteration : ", ntime, " ,time: ", Tdomain%TimeD%rtime
    endif
    call save_field_h5(Tdomain, rg, isort)

    if(rg == 0)then
        write(78,*) isort, Tdomain%TimeD%rtime
        call semname_nb_proc(isort,fnamef)
        open (79,file = fnamef,status="UNKNOWN")
        write(79,*) nb_procs
        close(79)
    endif

    isort = isort + 1  ! a faire avant le save_checkpoint

end subroutine OUTPUT_SNAPSHOTS
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
subroutine END_SEM(Tdomain,rg,ntime)
    use sdomain
    use mCapteur
    use semdatafiles
    use semconfig !< pour config C
    use sem_c_bindings
    use scouplage

    implicit none

    type(domain), target :: Tdomain
    integer, intent(in)  :: rg,ntime

    
    open (111,file = "fin_sem")
    if (ntime >= Tdomain%TimeD%NtimeMax-1) then
        write(111,*) 1
    else
        write(111,*) 0
    end if
    close(111)


    call flushAllCapteurs(Tdomain, rg)

    if (rg == 0) write (*,*) "--> DEALLOCATING DOMAIN."
    call deallocate_domain (Tdomain, rg)
#ifdef COUPLAGE
    close(50)
#endif

end subroutine END_SEM

end module drive_sem
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
