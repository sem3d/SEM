!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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
#ifndef COUPLAGE
    logical, parameter :: couplage = .false.
#else
    logical, parameter :: couplage = .true.
#endif

    if (.not. couplage) then
        call init_sem_path(p_param, p_traces, p_results, p_data, p_prot, p_prop, p_prop_h5)
    else
        call init_mka3d_path()
    endif

    call sem(couplage)

end program main

subroutine  sem(couplage)
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
    use scouplage
    implicit none

    logical, intent(in) :: couplage
    type (domain), target  :: Tdomain
    integer :: ntime,i_snap, ierr
    integer :: isort
    character(len=MAX_FILE_SIZE) :: fnamef
    integer :: info_capteur
    real(kind=8), parameter :: max_time_left=900
    integer :: getpid, pid

#ifdef COUPLAGE
    integer :: n
    integer, dimension(3) :: flags_synchro ! fin/protection/sortie
#else
    real(kind=8) :: remaining_time
#endif
    integer :: display_iter !! Indique si on doit faire des sortie lors de cette iteration
    real(kind=4), dimension(2) :: tarray
    real(kind=4) :: tref
    real(kind=4), parameter :: display_iter_time = 5.
    integer :: interrupt, rg, code, protection

    pid = getpid()
    !write(*,*) "SEM2D[", pid, "] : Demarrage."

    Tdomain%couplage = couplage
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

    if (couplage) then
        if (rg == 0) write(6,'(A)') 'Methode de couplage Mka3D/SEM2D 4: ', &
            'Peigne de points d''interpolation, en vitesses, systeme lineaire sur la vitesse, ', &
            'contraintes discontinues pour les ddl de couplage'

        call initialisation_couplage(Tdomain)

        ! recalcul du nbre d'iteration max
        Tdomain%TimeD%ntimeMax = int (Tdomain%TimeD%Duration/Tdomain%TimeD%dtmin)
    endif



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
        sortie_capteur_deformation = .FALSE.
    endif



    if (Tdomain%logicD%save_snapshots .or. Tdomain%logicD%save_deformation) then
        Tdomain%timeD%nsnap = int(Tdomain%TimeD%time_snapshots / Tdomain%TimeD%dtmin)
        if (Tdomain%timeD%nsnap == 0) Tdomain%timeD%nsnap = 1
        if (rg==0) write(*,*) "Snapshot every ", Tdomain%timeD%nsnap, " iterations"
    endif


    i_snap =1;

    if (couplage) then
        if(Tdomain%TimeD%iter_reprise>0) then
            ntime = Tdomain%TimeD%iter_reprise
        else
            ntime = 0
        endif
        call reception_surface_part_mka(Tdomain)
        call reception_nouveau_pdt_sem(Tdomain)
    endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BOUCLE DE CALCUL EN TEMPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call dtime(tarray, tref)
    protection = 0
    interrupt = 0
    do ntime= Tdomain%TimeD%NtimeMin, Tdomain%TimeD%NtimeMax-1

        Tdomain%TimeD%ntime = ntime
        protection = 0
        if (interrupt>0) then
            if (rg==0) write(*,*) "Sortie sur limite de temps..."
            exit
        end if

        call Newmark (Tdomain)

        if (ntime==Tdomain%TimeD%NtimeMax-1) then
            interrupt=1
        endif

        ! incrementation du pas de temps
        Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin

#ifdef COUPLAGE
        call envoi_vitesse_mka(Tdomain) !! syst. lineaire vitesse

        ! traitement de l arret
        ! comm entre le master_sem et le Tdomain%master_superviseur : reception de la valeur de l interruption

        flags_synchro(1) = interrupt
        flags_synchro(2) = 0      ! SEM ne declenche pas les protections
        flags_synchro(3) = 0

        call MPI_Allreduce(MPI_IN_PLACE, flags_synchro, 3, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, code)

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
        call MPI_Allreduce(MPI_IN_PLACE, interrupt, 1, MPI_INTEGER, MPI_SUM, &
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

            if (rg==0 .and. display_iter==1) write(*,*) "Snapshot:",isort," iteration=", ntime, " tps=", Tdomain%TimeD%rtime

            call save_field_h5(Tdomain, rg, isort)

            ! Sortie Energie totale du systeme
            if(Tdomain%LogicD%CompEnerg) call global_energy(Tdomain)

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
            call save_checkpoint(Tdomain,Tdomain%TimeD%rtime,Tdomain%TimeD%dtmin,ntime,isort)
        endif

        ! arret des calculs sur tous les procs
        if (interrupt/=0) then
            print*,"Arret de SEM, iteration=", ntime
            exit
        endif

    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FIN BOUCLE DE CALCUL EN TEMPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call END_SEM(Tdomain, ntime)

    if (Tdomain%bCapteur) then
        deallocate(Tdomain%GrandeurDeformation)
        deallocate(Tdomain%GrandeurVitesse)
    endif

    if (.not. couplage) then
        call MPI_Finalize  (ierr)
        if (rg == 0) write (*,*) "Execution completed"
    endif
    if (rg == 0) close(78)

end subroutine sem


subroutine START_SEM(Tdomain)
    use sdomain
    use mpi
    implicit none
    type(domain), intent(inout) :: Tdomain
    integer :: rg, ierr
    integer :: global_rank, global_nb_proc
#ifdef COUPLAGE
    integer :: worldgroup, intergroup
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer :: m_localComm, comm_super_mka
    integer :: min_rank_glob_sem
    integer :: n, tag
    integer, dimension(3) :: tab
#endif

    call MPI_Init(ierr)
    call MPI_Comm_Rank (MPI_COMM_WORLD, global_rank, ierr)
    call MPI_Comm_size (MPI_COMM_WORLD, global_nb_proc, ierr)

#ifdef COUPLAGE
    call MPI_Comm_split(MPI_COMM_WORLD, 2, Tdomain%Mpi_var%my_rank, m_localComm, ierr)
    call MPI_Comm_Rank (m_localComm, Tdomain%Mpi_var%my_rank, ierr)
    call MPI_Comm_size (m_localComm, Tdomain%Mpi_var%n_proc, ierr)

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
    if (Tdomain%Mpi_var%my_rank == 0) then
        tab(1) = global_rank
        tab(2) = Tdomain%Mpi_var%n_proc
        tab(3) = min_rank_glob_sem
        do n=1, global_nb_proc
            tag=8200000+n
            call MPI_Send(tab, 3, MPI_INTEGER, n-1, tag, MPI_COMM_WORLD, ierr)
        enddo
    endif

    !! Reception des infos de sem
    tag=8200000+global_rank+1
    call MPI_Recv(tab, 3, MPI_INTEGER, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, status, ierr)

    call MPI_Comm_group(MPI_COMM_WORLD, worldgroup, ierr)
    call MPI_Group_incl(worldgroup, 2, tab, intergroup, ierr)
    call MPI_Comm_create(MPI_COMM_WORLD, intergroup, comm_super_mka, ierr)

#else
    Tdomain%Mpi_var%my_rank = global_rank
    Tdomain%Mpi_var%n_proc = global_nb_proc
    Tdomain%communicateur=MPI_COMM_WORLD
    Tdomain%communicateur_global=MPI_COMM_WORLD
#endif

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
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
