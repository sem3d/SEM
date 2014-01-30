!>
!!\file main.F90
!!\brief Assure l'appel au code de calcul SEM2D.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

#ifndef COUPLAGE

program main
    use mpi
    use semdatafiles
    implicit none

    character(Len=MAX_FILE_SIZE),parameter :: p_param = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_traces = "./traces"
    character(Len=MAX_FILE_SIZE),parameter :: p_results = "./res"
    character(Len=MAX_FILE_SIZE),parameter :: p_data = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_prot = "./prot"

    call init_sem_path(p_param, p_traces, p_results, p_data, p_prot)
    call sem(0, MPI_COMM_WORLD, MPI_COMM_WORLD)

end program main
#endif

subroutine  sem(master_superviseur,communicateur,communicateur_global)
    use sdomain
    use mCapteur
    use semdatafiles
    use mpi
    use msnapshots
    use sem_c_bindings
    use shape_lin
    use shape_quad
#ifdef COUPLAGE
    use scouplage
#endif
    use snewmark

    implicit none
    integer, intent(in) :: communicateur,communicateur_global,master_superviseur

    type (domain), target  :: Tdomain
    integer :: ntime,i_snap, ierr
    integer :: isort
    character(len=MAX_FILE_SIZE) :: fnamef
    integer :: info_capteur
    real(kind=8) :: remaining_time
    real(kind=8), parameter :: max_time_left=900

#ifdef COUPLAGE
    integer :: groupe
    integer :: sortie
    integer :: finSem,nrec
    integer :: tag, MaxNgParFace
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer, dimension(3) :: flags_synchro ! fin/protection/sortie
#endif
    integer :: display_iter !! Indique si on doit faire des sortie lors de cette iteration
    real(kind=4), dimension(2) :: tarray
    real(kind=4) :: tref
    real(kind=4), parameter :: display_iter_time = 5.
    integer :: interrupt, rg, code, protection

    display_iter = 1

#ifdef COUPLAGE
    Tdomain%communicateur=communicateur
    Tdomain%communicateur_global=communicateur_global
    Tdomain%master_superviseur=master_superviseur
    call MPI_Comm_Group(Tdomain%communicateur,groupe,ierr)
    call MPI_Group_Size(groupe,Tdomain%Mpi_var%n_proc,ierr)
    call MPI_Group_Rank(groupe,Tdomain%Mpi_var%my_rank,ierr)
#else
    ! initialisations MPI
    call MPI_Init (ierr)
    call MPI_Comm_Rank (MPI_COMM_WORLD, Tdomain%Mpi_var%my_rank, ierr)
    call MPI_Comm_Size (MPI_COMM_WORLD, Tdomain%Mpi_var%n_proc,  ierr)
    Tdomain%communicateur=MPI_COMM_WORLD
#endif
    rg = Tdomain%Mpi_var%my_rank

#ifdef COUPLAGE
    call init_mka3d_path()
#endif
    if (rg == 0) call create_sem_output_directories()

    !lecture du fichier de donnee
    if (rg == 0) write (*,*) "Read input.spec"
    call read_input (Tdomain)

    !lecture du fichier de maillage unv avec conversion en fichier sem2D
    if (rg == 0) write (*,*) "Define mesh properties"
    call read_mesh_h5(Tdomain)

    if (rg == 0) write (*,*) "Compute Gauss-Lobatto-Legendre weights and zeroes"
    call compute_GLL (Tdomain)

    if (rg == 0) write (*,*) "Define a global numbering for the collocation points"
    call global_numbering (Tdomain)

    if (rg == 0) write (*,*) "Computing shape functions within thier derivatives"
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

    if (rg == 0) write (*,*) "Attribute PML properties"
    call PML_definition (Tdomain)

    if (Tdomain%logicD%any_source) then
        if (rg == 0) write (*,*) "Computing point-source parameters and location"
        call SourcePosition(Tdomain)
    endif

    if (Tdomain%logicD%save_trace ) then
        if (rg == 0) write (*,*) "Computing receivers parameters and locations"
        call ReceiverPosition(Tdomain)
    endif

    ! Initialisations des capteurs par S Terrana Jauary 2014
    open (90,file="capteur_sourceX",status="replace",form="formatted")
    call init_capteurs_veloc (Tdomain)

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

    isort = 1

#ifdef COUPLAGE
    if (rg == 0) write(6,'(A)') 'Methode de couplage Mka3D/SEM2D 4: Peigne de points d''interpolation, en vitesses, systeme lineaire sur la vitesse,', &
        ' contraintes discontinues pour les ddl de couplage'

    call initialisation_couplage(Tdomain, MaxNgParFace)
    finSem=0

    ! recalcul du nbre d'iteration max
    Tdomain%TimeD%ntimeMax = int (Tdomain%TimeD%Duration/Tdomain%TimeD%dtmin)

    if (Tdomain%logicD%save_trace) then
        i = 0
        do nrec = 0, Tdomain%n_receivers-1
            if (Tdomain%sReceiver(nrec)%located_here) i= i + 1
        enddo
        if (i>0) then
            deallocate(Tdomain%Store_Trace)
            allocate (Tdomain%Store_Trace(0:1,0:i-1,0:Tdomain%TimeD%ntimeMax-1))
        endif
    endif

#endif



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

    ! preparation des eventuelles sorties capteur
    info_capteur = 0
    if (Tdomain%bCapteur) call read_capteur(Tdomain,info_capteur)
    if(info_capteur /= 0) then
        Tdomain%bCapteur = .FALSE.
        sortie_capteur = .FALSE.
        sortie_capteur_vitesse = .FALSE.
        sortie_capteur_depla = .FALSE.
        sortie_capteur_deformation = .FALSE.
    endif



    if (Tdomain%logicD%save_snapshots .or. Tdomain%logicD%save_deformation) then
        Tdomain%timeD%nsnap = Tdomain%TimeD%time_snapshots / Tdomain%TimeD%dtmin
        if (Tdomain%timeD%nsnap == 0) Tdomain%timeD%nsnap = 1
        write(*,*) "Snapshot every ", Tdomain%timeD%nsnap, " iterations"
    endif


    i_snap =1;

#ifdef COUPLAGE
    if(Tdomain%TimeD%iter_reprise>0) then
        ntime = Tdomain%TimeD%iter_reprise
    else
        ntime = 0
    endif

    call reception_surface_part_mka(Tdomain)
    call reception_nouveau_pdt_sem(Tdomain)
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BOUCLE DE CALCUL EN TEMPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call dtime(tarray, tref)
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
            call Runge_Kutta4(Tdomain, ntime, Tdomain%TimeD%dtmin)
        endif

        if (ntime==Tdomain%TimeD%NtimeMax-1) then
            interrupt=1
        endif

#ifdef COUPLAGE
        call envoi_vitesse_mka(Tdomain, ntime) !! syst. lineaire vitesse

        ! traitement de l arret
        ! comm entre le master_sem et le Tdomain%master_superviseur : reception de la valeur de l interruption

        flags_synchro(1) = interrupt
        flags_synchro(2) = 0      ! SEM ne declenche pas les protections
        flags_synchro(3) = 0

        call MPI_ALLREDUCE(MPI_IN_PLACE, flags_synchro, 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, code)

        interrupt = flags_synchro(1)
        protection = flags_synchro(2)
        if (flags_synchro(3)/=0) then
            i_snap = 0
        else
            i_snap = 1
        end if
#else
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
#endif

        ! Ici, on a une info globale pour interrupt, protection, i_snap
        if (interrupt>0) then
            protection = 1
        end if

        if (Tdomain%bCapteur) call evalueSortieCapteur(ntime)

        if (i_snap == 0 .or. sortie_capteur) then


#ifdef COUPLAGE
            if (rg == 0.and.i_snap == 0 .and. display_iter==1) then
                write(*,'(a35,i8.8,a8,f10.5)') "SEM : sortie resultats iteration : ",ntime, &
                    " temps : ",Tdomain%TimeD%rtime
            endif
#endif
            if (rg==0) write(*,*) "Snapshot iteration=", ntime, " tps=", Tdomain%TimeD%rtime

            call save_field_h5(Tdomain, rg, isort)
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


        ! sauvegarde des vitesses ?
        if (Tdomain%logicD%save_trace) call save_trace(Tdomain, ntime)


        if (i_snap==0) then
            isort=isort+1  ! a faire avant le save_checkpoint
        endif
#ifdef COUPLAGE

        ! remise a zero dans tous les cas des forces Mka sur les faces
        do n = 0, Tdomain%n_face-1
            Tdomain%sFace(n)%ForcesMka = 0
        enddo
        ! remise a zero dans tous les cas des forces Mka sur les vertex
        do n = 0, Tdomain%n_vertex-1
            Tdomain%sVertex(n)%ForcesMka = 0
        enddo

#endif

        ! sortie des quantites demandees par les capteur
        if (sortie_capteur) call save_capteur(Tdomain, ntime)

        if (protection/=0) then
            !call save_checkpoint(Tdomain,Tdomain%TimeD%rtime,Tdomain%TimeD%dtmin,ntime,isort)
        endif

        ! arret des calculs sur tous les procs
        if (interrupt/=0) then
            print*,"Arret de SEM, iteration=", ntime
            exit
        endif

        ! incrementation du pas de temps
        Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin

        ! Sorties des Capteurs :
        call capteurs_veloc (Tdomain,Tdomain%TimeD%rtime)

    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FIN BOUCLE DE CALCUL EN TEMPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    if (Tdomain%logicD%save_trace) call dump_trace(Tdomain)

    if (Tdomain%bCapteur) then
        deallocate(Tdomain%GrandeurDeformation)
        deallocate(Tdomain%GrandeurVitesse)
    endif

#ifdef COUPLAGE
    ! on quitte directement SEM, la fin du parallelisme est geree par le main.C de la maquette
#else
    call MPI_Finalize  (ierr)
    if (rg == 0) write (*,*) "Execution completed"
#endif
    if (rg == 0) close(78)

end subroutine sem

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
