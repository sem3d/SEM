!>
!! \file drive_sem.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<
subroutine  sem(master_superviseur, communicateur, communicateur_global)
    use sdomain
    use Postraitement
    use mCapteur
    use semdatafiles
    use mpi
#ifdef COUPLAGE
    use scouplage
#endif

    implicit none

    ! Hors couplage on doit avoir -1 MPI_COMM_WORLD, MPI_COMM_WORLD
    integer, intent(in) :: communicateur, communicateur_global, master_superviseur

    type(domain), target :: Tdomain

    integer :: code, rg, nb_procs, ntime, i_snap, n, icount
#ifdef MKA3D
    integer :: i, isort
    character(len=24) :: nom_dir_sorties
    integer :: info_capteur
    real(kind=8) :: remaining_time
    real(kind=8), parameter :: max_time_left=900
    integer :: protection
#endif
    character(Len=MAX_FILE_SIZE) :: fnamef
#ifdef COUPLAGE
    integer :: groupe
    integer :: sortie
    integer :: finSem
    integer :: tag
    integer, dimension (MPI_STATUS_SIZE) :: status
    character*2 :: sit
    integer :: MaxNgParDir
#endif
    integer, dimension(3) :: flags_synchro ! fin/protection/sortie
    integer :: interrupt

    Tdomain%communicateur = communicateur
    Tdomain%communicateur_global = communicateur_global
    Tdomain%master_superviseur = master_superviseur
#ifdef COUPLAGE
    call MPI_Comm_Group(Tdomain%communicateur, groupe, code)
    call MPI_Group_Size(groupe, nb_procs, code)
    call MPI_Group_Rank(groupe, rg, code)
#else
    ! initialisations MPI
    call MPI_Init (code)
    call MPI_Comm_Rank (MPI_COMM_WORLD, rg, code)
    call MPI_Comm_Size (MPI_COMM_WORLD, nb_procs,  code)
#endif


#ifdef MKA3D
    if (rg == 0) then
        write(*,*)  "execution avec DIRECTIVE=MKA3D"
    endif
    call system('mkdir -p Resultats')
    call system('mkdir -p Capteurs/sem')
#endif

    ! ##############  Begin of the program  #######################
    write (*,*) "Define mesh properties ",rg
    call read_input (Tdomain, rg, code)

    write(*,*) "----         SEM - 3d version      ------"
    write(*,*) "-----------------------------------------"
    write(*,*)

    write(*,*) "  --> Reading run properties, proc. # ",rg
    call MPI_BARRIER(MPI_COMM_WORLD, code)

    !  modif mariotti fevrier 2007 cea
    if (Tdomain%logicD%super_object_local_present) then
        if (Tdomain%super_object_type == "P") then
            if (rg == 0) write(*,*) "Define Plane Wave properties",rg
            call define_planew_properties (Tdomain)
        endif
    endif
    if (Tdomain%logicD%neumann_local_present) then
        if (rg == 0) write(*,*) "Define Neumann properties",rg
        call define_Neumann_properties(Tdomain,rg)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD, code)

    if (rg == 0) write (*,*) "Compute Gauss-Lobatto-Legendre weights and zeroes ",rg
    call compute_GLL (Tdomain)
    call MPI_BARRIER(MPI_COMM_WORLD, code)

    if (rg == 0) write (*,*) "Define a global numbering for the collocation points ",rg
    call global_numbering (Tdomain,rg)
    call MPI_BARRIER(MPI_COMM_WORLD, code)

    if (rg == 0) write (*,*) "Compute shape functions within their derivatives ",rg
    if (Tdomain%n_nodes == 8) then
        call shape8(Tdomain,rg)   ! Linear interpolation
    else if (Tdomain%n_nodes == 20) then
        write (*,*) "20 control points not yet implemented in the code. Wait for an upgrade"
        stop
    else
        call shape27(TDomain)   ! Quadratic interpolation
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,code)

    call compute_Courant(Tdomain,rg)
    call MPI_BARRIER(MPI_COMM_WORLD,code)
    if (Tdomain%any_PML)   then
        if (rg == 0) write (*,*) "Attribute PML properties ",rg
        call PML_definition (Tdomain)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,code)

    if (Tdomain%logicD%any_source) then
        if (rg == 0) write (*,*) , "Compute source parameters ",rg
        call SourcePosition (Tdomain, rg)
        call double_couple (Tdomain, rg)
        call source_excit(Tdomain,rg)
        call def_timefunc (Tdomain, rg)
        !- pour entrees en temps: directement dans Modules/Source.f90
        ! ->  lecture d'un fichier-entree pour la source; valable pour une seule source
        do i = 0,Tdomain%n_source-1
            if(Tdomain%sSource(i)%i_time_function == 5)then
                call read_source(Tdomain%sSource(i))
            endif
        end do
    endif

    if (Tdomain%logicD%save_trace) then
        if (rg == 0) write (*,*) "Compute receiver locations ",rg
        call ReceiverPosition (Tdomain, rg)
    endif
    call MPI_BARRIER(MPI_COMM_WORLD,code)

    if (rg == 0) write (*,*) "Allocate fields ",rg
    call allocate_domain (Tdomain, rg)
    call MPI_BARRIER(MPI_COMM_WORLD,code)

    if (rg == 0) write (*,*) "Compute Courant parameter ",rg
    call compute_Courant (Tdomain,rg)

    if (rg == 0) write (*,*) "Compute mass matrix and internal forces coefficients ",rg
    call define_arrays (Tdomain, rg)
    call MPI_BARRIER(MPI_COMM_WORLD,code)

    if (Tdomain%n_sls>0) then
        if (Tdomain%aniso) then
            if (rg == 0) write (*,*) "Compute attenuation features aniso ", rg
            call set_attenuation_aniso_param(Tdomain)
        else
            if (rg == 0) write (*,*) "Compute attenuation features", rg
            call set_attenuation_param(Tdomain)
        endif
    endif

#ifndef MKA3D
    !ajout post-ensight
    if (rg == 0) write (*,*) " Create an Ensight gold geometry file"
    call POST_Ensight(Tdomain, rg)
#endif
    if (rg == 0) write (*,*) "Entering the time evolution ",rg
    ! initialisation des temps
    Tdomain%TimeD%rtime = 0
    Tdomain%TimeD%NtimeMin = 0
    call MPI_BARRIER(MPI_COMM_WORLD,code)

    isort = 1
#ifdef COUPLAGE
    if (rg == 0) write(6,'(A)') 'Methode de couplage Mka3D/SEM3D 4: Peigne de points d''interpolation, en vitesses, systeme lineaire sur la vitesse,', &
        ' contraintes discontinues pour les ddl de couplage'
    call semname_drive_sem_listing(rg, fnamef)
    open(UNIT=50, FILE=fnamef, STATUS='UNKNOWN')
    call initialisation_couplage(Tdomain, rg, MaxNgParDir, nb_procs)
    finSem=0

    ! recalcul du nbre d'iteration max
    Tdomain%TimeD%ntimeMax = int (Tdomain%TimeD%Duration/Tdomain%TimeD%dtmin)
#endif
    if (rg == 0) print *,'Tdomain%TimeD%ntimeMax', Tdomain%TimeD%ntimeMax,Tdomain%TimeD%dtmin,Tdomain%TimeD%Duration

#ifdef MKA3D
    ! traitement de la reprise
    !  modif mariotti fevrier 2007 cea
    call semname_drive_sem_resulttemp(fnamef)
    if (Tdomain%logicD%run_restart) then
        !! Il faudra ajouter la gravite ici #ifdef COUPLAGE
        call read_restart(Tdomain, rg, isort)
        write (*,*) "Reprise effectuee sur processeur ",rg
        open (78,file=fnamef,status="unknown",position="append")
    else
        ! Sauvegarde des donnees de post-traitement
        call save_post_data(Tdomain, rg)
        open (78,file=fnamef,status="unknown",position="rewind")
        ! on supprime tous les fichiers et repertoire de protection
        if(rg == 0) call system('rm -Rf ./ProRep/sem/Prot*')
    endif
#endif
    !stop

    ! preparation des eventuelles sorties capteur
    info_capteur = 0
    if (Tdomain%bCapteur) call read_capteur(Tdomain, rg, info_capteur)

    if(info_capteur /= 0) then
        Tdomain%bCapteur = .FALSE.
        sortie_capteur = .FALSE.
        sortie_capteur_vitesse = .FALSE.
        sortie_capteur_depla = .FALSE.
        sortie_capteur_deformation = .FALSE.
    endif


    if (Tdomain%logicD%save_snapshots) then
        Tdomain%timeD%nsnap = Tdomain%TimeD%time_snapshots / Tdomain%TimeD%dtmin
        icount = 0
    endif

    i_snap = 1

#ifdef COUPLAGE
    call reception_surface_part_mka(Tdomain, rg)
    call reception_nouveau_pdt_sem(Tdomain, rg)
#endif


    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! BOUCLE DE CALCUL EN TEMPS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    interrupt = 0
    do ntime = Tdomain%TimeD%NtimeMin, Tdomain%TimeD%NtimeMax-1

        protection = 0
        if (interrupt>0) then
            if (rg==0) write(*,*) "Sortie sur limite de temps..."
            exit
        end if
        call Newmark (Tdomain, rg, ntime)
        !doute sur la modif de rtime. A l'origine, ligne suivante, a priori meilleur en fin de boucle
        !       Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin

        if (ntime==Tdomain%TimeD%NtimeMax-1) then
            interrupt=1
        endif

#ifdef COUPLAGE
        call envoi_vitesse_mka(Tdomain, ntime, rg) !! syst. lineaire vitesse

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
#ifdef MKA3D
        call tremain( remaining_time )
        if (rg==0) write(*,*) "remain:", remaining_time
        if (remaining_time<max_time_left) then
            interrupt = 1
        end if
        call mpi_allreduce(MPI_IN_PLACE, interrupt, 1, MPI_INTEGER, MPI_SUM, Tdomain%communicateur_global, code)
#endif

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
            if (rg == 0.and.i_snap == 0) then
                write(*,'(a35,i6.6,a8,f10.5)') "SEM : sortie resultats iteration : ", ntime, " temps : ", Tdomain%TimeD%rtime
            endif

            if ( (Tdomain%logicD%save_snapshots.and.i_snap==0).or.sortie_capteur_vitesse ) then
#ifdef MKA3D
                call savefield (Tdomain, isort, sortie_capteur_vitesse, rg, i_snap, nom_dir_sorties)
#else
                if(i_snap == 0) then
                    icount = icount+1
                    call savefield (Tdomain, ntime, sortie_capteur_vitesse, rg, icount, i_snap)
                endif
#endif

            endif


#ifdef MKA3D
            if ( (Tdomain%logicD%save_snapshots.and.i_snap==0).or.sortie_capteur_depla ) then
                call savefield_disp (Tdomain, isort, sortie_capteur_depla, rg, i_snap, nom_dir_sorties)
            endif
#endif

        endif

#ifndef MKA3D
        if (Tdomain%logicD%save_trace) then
            if (mod (ntime+1,Tdomain%TimeD%ntrace) == 0) then
                call savetrace (Tdomain,rg,int(ntime/Tdomain%TimeD%ntrace))
                do n = 0, Tdomain%n_receivers-1
                endif
            endif
        endif
#endif

#ifdef MKA3D
        if(Tdomain%logicD%save_snapshots .and. i_snap == 0) then
            write (*,*) "A new snapshot is recorded for processor",rg
            if(rg==0) then
                write(78,*) isort, Tdomain%TimeD%rtime
                call semname_nb_proc(nom_dir_sorties,fnamef)
                open (79,file = fnamef,status="UNKNOWN")
                write(79,*) nb_procs
                close(79)
            endif
        endif

        if (i_snap == 0) then
            isort = isort + 1  ! a faire avant le save_checkpoint
        endif

        ! remise a zero dans tous les cas des forces Mka sur les faces
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

        ! sortie des quantites demandees par les capteur
        if (sortie_capteur) call save_capteur(Tdomain, ntime, rg)


        if (protection/=0) then
            write (*,*) " sauvegarde  ",ntime," sur processeur ",rg
            call save_checkpoint(Tdomain, Tdomain%TimeD%rtime, ntime, rg, Tdomain%TimeD%dtmin, isort)
        endif

        ! arret des calculs sur tous les procs
        if (interrupt/=0) then
            print*,"Arret de SEM, iteration=", ntime
            exit
        endif

        ! incrementation du pas de temps !on copie Sem2d
        Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! FIN BOUCLE EN TEMPS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    write (*,*) "Deallocate fields ",rg

    if (Tdomain%bCapteur) then
        deallocate(Tdomain%GrandeurVitesse)
    endif

    if (rg == 0) write (*,*) "Deallocate domain ",rg
    call deallocate_domain (Tdomain, rg)
#ifdef COUPLAGE
    close(50)
#endif

end subroutine sem
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
