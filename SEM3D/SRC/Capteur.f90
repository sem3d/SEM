!>
!!\file Capteur.f90
!!\brief Permet de manipuler les quantités associées aux capteurs.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module mCapteur

    use sdomain
    use semdatafiles
    use mpi
    use mfields
    use sem_hdf5

    implicit none

    public :: read_capteur, save_capteur, evalueSortieCapteur, flushAllCapteurs
    private :: lireCapteur, capterPointDeGauss, sortirGrandeurSousCapteur

    ! Les fichiers capteurs sont ecrits toutes les NCAPT_CACHE sorties
    integer, parameter :: NCAPT_CACHE=100

    type :: tPtgauss
       integer :: num      ! numero du point de Gauss
       integer :: nel, i,j,k ! Numero de l'element, et numerotation interne du pt de gauss
       real    :: Coord(3) ! coordonnees du point de gauss
       real    :: Poids    ! poids du point de Gauss
       type(tPtGauss),pointer :: suivant  ! point de Gauss suivant
    end type Tptgauss


    type :: tCapteur

       integer :: numero             ! numero du capteur
       real :: rayon                 ! rayon d'action du capteur
       real :: distanceMin           ! distance au pt de Gauss le plus proche
       integer :: frequence          ! frequence de captation des grandeur
       real, dimension (3) :: Coord  ! localisation du capteur
       character(LEN=3) :: operation ! type d operation a mener sur la grandeur
       character(LEN=50) :: grandeur ! grandeur physique a extraire
       character(LEN=20) :: nom      ! nom du capteur
       integer :: nPtGauss           ! nombre de points de Gauss
       type(tPtgauss),pointer :: listePtGauss ! liste des points de Gauss associe au capteur
       type(Tcapteur),pointer :: suivant ! pour passer au capteur suivant
       integer :: type_calcul        ! 0 => pas d'interpolation, recherche du point de Gauss le plus proche,
       ! 1 => interpolation aux positions du capteur
       !   (necessite le calcul des abscisses curvilignes dans l'element de reference
       !     = processus de dichotomie pour trouver (xi,eta,zeta))
       integer :: n_el ! numero de la maille dans laquelle se trouve le capteur
       ! si le capteur est partage entre plusieurs mailles, une seule suffit (type_calcul=1)
       real :: xi, eta, zeta ! abscisses curvilignes pour le capteur en cas d'interpolation (type_calcul=1)
       integer :: numproc               ! numero du proc localisant le capteur
       integer :: icache
       real, dimension(4,NCAPT_CACHE) :: valuecache
    end type tCapteur

    integer           :: dimCapteur        ! nombre total de capteurs

    type(Tcapteur), pointer :: listeCapteur

    integer,parameter :: fileIdCapteur=200  ! id fichier capteur

    logical :: traces_h5_created
contains

    function grandeur_depla(Tdomain, PtGauss)
        type(domain), intent(in) :: Tdomain
        type(tPtgauss), intent(in) :: PtGauss
        real, dimension(3) :: grandeur_depla

    end function grandeur_depla

    function grandeur_vitesse(Tdomain, PtGauss)
        type(domain), intent(in) :: Tdomain
        type(tPtgauss), intent(in) :: PtGauss
        real, dimension(3) :: grandeur_vitesse

    end function grandeur_vitesse

    !--------------------------------------------------------------

    !>
    !! \brief La routine read_capteur() permet de charger les caractéristiques des capteurs définies dans les fichiers d'entrée.
    !! Elle est appelée une seule fois, avant la boucle en temps, afin d'initialiser les capteurs.
    !! \param type (domain), target Tdomain
    !<
    subroutine read_capteur(tDomain, rg, info_capteur)

        implicit none

        type (domain), target  :: Tdomain
        type(Tcapteur),pointer :: capteur
        type(tPtgauss),pointer :: PtGauss

        integer,parameter :: fid=98
        character(len=MAX_FILE_SIZE) :: fnamef

        integer :: rg
        integer :: ierr, tag
        integer , dimension  (MPI_STATUS_SIZE) :: status
        real :: distanceMinMin
        real,dimension(3) :: val, val0
        real,dimension(:),allocatable :: sendbuf, recvbuf
        integer info_capteur, numproc_max
        real xi, eta, zeta
        integer n_el

        ! lecture des parametres des capteurs

        !lecture du fichier des capteurs capteurs.dat
        call lireCapteur(tDomain, rg, info_capteur)

        ! En reprise on ne recree pas le fichier capteur
        if (Tdomain%TimeD%NtimeMin==0) then
            traces_h5_created = .false.
        else
            traces_h5_created = .true.
        endif

        if(info_capteur /= 0) return

        ! creer pour chaque capteur la liste du ou des points de gauss a prendre en compte
        call capterPointDeGauss(Tdomain%GlobCoord, Tdomain%n_glob_points)

        ! tableau de correspondance pour les com
        allocate(sendbuf(3*dimCapteur))
        allocate(recvbuf(3*dimCapteur))

        call dist_max_elem(Tdomain) !!ajout calcul de la distance max entre deux sommets a l'interieur d'une maille
        !!permet de limiter le nombre de mailles pouvant contenir un capteur

        ! afficher position des points de Gauss pour chaque capteur ?
        capteur=>listeCapteur

        do while (associated(capteur))

            tag=capteur%numero ! une com par capteur

            call semname_capteur_pos(capteur%nom,fnamef)

            if(capteur%type_calcul==0) then

                PtGauss=>capteur%listePtGauss

                if (capteur%rayon==0) then ! un seul point de Gauss a ecrire, mais le plus proche de tous

                    if (rg.eq.0) open(UNIT=fid,FILE=trim(fnamef),STATUS='replace')  !!modif seul proc 0
                    ! tous les proc recupere la distance min du pdg le + proche
                    call MPI_AllReduce(capteur%distanceMin, distanceMinMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, Tdomain%communicateur, ierr)
                    if (ABS(distanceMinMin - capteur%distanceMin) < 1e-5) then ! on est sur le proc qui a le pt de gauss le + proche
                        val(1) = PtGauss%coord(1) ! coord en x du pdg le plus proche
                        val(2) = PtGauss%coord(2) ! coord en y du pdg le plus proche
                        val(3) = PtGauss%coord(3) ! coord en z du pdg le plus proche
                        write(6,'(A30,I8,A13,I8)') "numero du pdg le plus proche:", PtGauss%num, " sur le proc ", rg
                        if (rg.EQ.0 ) then ! pdg le plus proche est sur le proc 0
                            val0(1) = val(1)
                            val0(2) = val(2)
                            val0(3) = val(3)
                        else ! on envoie le pdg le plus proche au proc 0
                            sendbuf(3*(dimCapteur-1)+1) = val(1)
                            sendbuf(3*(dimCapteur-1)+2) = val(2)
                            sendbuf(3*(dimCapteur-1)+3) = val(3)
                            call MPI_SEND(sendbuf, 3*dimCapteur, MPI_DOUBLE_PRECISION, 0, tag,  Tdomain%communicateur, ierr)
                        endif
                    else
                        if (rg .eq.0 ) then ! proc 0 n'a pas le pdg le + proche et doit donc le recevoir

                            call mpi_recv(recvbuf,3*dimCapteur,MPI_DOUBLE_PRECISION,mpi_any_source,tag, Tdomain%communicateur,status,ierr)

                            val0(1) = recvbuf(3*(dimCapteur-1)+1)
                            val0(2) = recvbuf(3*(dimCapteur-1)+2)
                            val0(3) = recvbuf(3*(dimCapteur-1)+3)
                        endif
                    endif

                    if (rg .eq.0 ) then
                        rewind(fid)
                        write(fid,'(4(A6,1X,E12.4,1X))') "X=", val0(1), " Y=", val0(2), " Z=", val0(3), " dist=", distanceMinMin
                        close(fid)
                    endif
                else ! les proc ecrivent tour a tour les pt de gauss

                    open(UNIT=fid,FILE=trim(fnamef),STATUS='replace')  !!modif // version initiale Capteur Sem2d
                    ! on parcourt les pdg
                    do while (associated(PtGauss))
                        write(fid,*) PtGauss%coord(1), PtGauss%coord(2), PtGauss%coord(3)
                        PtGauss=>PtGauss%suivant
                    enddo
                    close(fid) !!modif // version initiale Capteur Sem2d
                endif

                !deallocate(PtGauss)

!!!    close(fid) !!version initiale Capteur Sem2d
            elseif(capteur%type_calcul==1) then
                xi = -1.
                eta = -1.
                zeta  = -1.
                n_el = -1
                call trouve_capteur(Tdomain, rg, capteur, n_el, xi, eta, zeta)

                val(1) = xi
                val(2) = eta
                val(3) = zeta

                call MPI_AllReduce(capteur%numproc, numproc_max, 1, MPI_INTEGER, MPI_MAX, Tdomain%communicateur, ierr)

                ! attention si le capteur est partage par plusieurs procs. On choisit le proc de num max
                capteur%numproc = numproc_max
                if(rg==numproc_max) then
                    capteur%n_el = n_el ! uniquement pour l'affichage dans le fichier *position
                    capteur%xi = xi !val0(1)
                    capteur%eta = eta !val0(2)
                    capteur%zeta = zeta !val0(3)

                    open(UNIT=fid,FILE=trim(fnamef),STATUS='replace')
                    rewind(fid)
                    write(fid,'(3(A6,1X,1pe15.8,1X),A22,1X,I6)') "xi=", capteur%xi, " eta=", capteur%eta, &
                        " zeta=", capteur%zeta, " numero element Sem", capteur%n_el
                    close(fid)
                end if
                !! on reattribue au proc son numero d'element initial
            endif
            capteur=>capteur%suivant

        enddo

        deallocate(sendbuf)
        deallocate(recvbuf)

    end subroutine read_capteur


    !------------------------------------------------------------------------------------------------

    !>
    !! \brief Lit les caratéristiques des capteurs.
    !!
    !! \param type (Domain), intent (IN) Tdomain
    !<
    subroutine lireCapteur(Tdomain, rg, info_capteur)

        implicit none


        type (Domain), intent (IN) :: Tdomain

        integer            :: CodeErreur,fileId,i, rg
        character(LEN=MAX_FILE_SIZE) :: ligne,fnamef
        logical            :: status

        type(Tcapteur),pointer :: capteur
        integer info_capteur
        character(LEN=6) :: char_calcul

        ! Id des fichiers de sortie des capteur
        fileId=99

        ! nombre de capteurs dans le fichier capteurs.dat
        dimCapteur=0

        ! liste initialement vide
        nullify(listeCapteur)

        !controle d'existence du fichier
        INQUIRE(File=trim(Tdomain%station_file),Exist=status)
        info_capteur=0
        if ( .not.status ) then
            if(rg==0) &
                write (*,*)"fichier des capteurs introuvable :",trim(Tdomain%station_file)
            !!      stop
            info_capteur=1
            return
        else
            if(rg==0) then
                write(*,*) "Reading file:", trim(adjustl(Tdomain%station_file))
            end if
        end if

        open(UNIT=fileIdCapteur,IOSTAT=CodeErreur,FILE=trim(Tdomain%station_file),FORM='formatted',STATUS='old',ACTION='read')
        if (CodeErreur .ne.0 ) print*,'Pb a l''ouverture du fichier Capteur  :CodeErreur=',CodeErreur

        do ! while (.not.eof(fileIdCapteur))


            ! initialisation d un nouveau capteur
            dimCapteur=dimCapteur+1   ! nbre total de capteur
            allocate(capteur)
            capteur%nom = ""         ! ses caracteristiques par defaut
            capteur%numero = dimCapteur ! numero du capteur
            capteur%frequence = 10
            capteur%rayon = 0
            capteur%operation = ""
            capteur%grandeur = ""
            capteur%coord(:) = -1.
            capteur%n_el = -1
            capteur%numproc = -1
            capteur%type_calcul = 1 !mode de calcul du capteur
            nullify(capteur%listePtGauss)
            capteur%distanceMin = huge(1.)
            capteur%nPtGauss = 0  ! initialisation du nombre de points de Gauss a prendre en compte
            capteur%icache = 0
            read(fileIdCapteur,'(A)',end=100) ligne ! titre 1
            read(fileIdCapteur,'(A)',end=100) ligne ! titre 2


            ! recuperation du nom du capteur
            read(fileIdCapteur,'(A)',end=100) ligne


            do while (trim(ligne).ne."")


                if (ligne(1:11).eq."NOM_CAPTEUR") capteur%nom = trim(ligne(12:100))
                if (ligne(1:4).eq."FREQ") read(ligne(5:100),*,ERR=201)capteur%frequence
                if (ligne(1:5).eq."RAYON") read(ligne(6:100),*,ERR=202)capteur%rayon
                if (ligne(1:9).eq."OPERATION") capteur%operation=trim(adjustl(ligne(10:100)))
                if (ligne(1:8).eq."GRANDEUR") capteur%grandeur=trim(adjustl(ligne(9:100)))
                if (ligne(1:6).eq."COORDX") read(ligne(7:100),*,ERR=203) capteur%Coord(1)
                if (ligne(1:6).eq."COORDY") read(ligne(7:100),*,ERR=204) capteur%Coord(2)
                if (ligne(1:6).eq."COORDZ") read(ligne(7:100),*,ERR=205) capteur%Coord(3)

                if (ligne(1:11).eq."TYPE_CALCUL") then

                    read(ligne(12:100),'(A6)',ERR=206) char_calcul

                    if(trim(adjustl(char_calcul)) == 'GAUSS') then
                        capteur%type_calcul = 0
                    elseif    (trim(adjustl(char_calcul)) == 'INTERP') then
                        capteur%type_calcul = 1
                    endif

                endif

                !if (ligne(1:8).eq."MATERIAU") read(ligne(9:100),'(I4)',,ERR=202)capteur%materiau


                ! recuperation du nom du capteur
                read(fileIdCapteur,'(A)',end=50) ligne
            enddo ! fin de lecture du capteur courant


50          if (rg==0) then
                write(*,*) "capteur : ",dimCapteur," ", capteur%nom
                write(*,*) "- freq = ", capteur%frequence
                write(*,*) "- pos = ", capteur%Coord
            end if

            capteur%suivant=>listeCapteur
            listeCapteur=>capteur
            ! si c'est un nouveau run, suppression de l'eventuel fichier de sortie des capteurs
            if ( .not.Tdomain%logicD%run_restart .and. rg==0) then

                call semname_capteur_type(capteur%nom,"_deformation",fnamef)

                open(fileId,file=trim(fnamef),status="replace",form="formatted")
                close(fileId)

                call semname_capteur_type(capteur%nom,"_vitesse",fnamef)

                open(fileId,file=trim(fnamef),status="replace",form="formatted")
                close(fileId)

                call semname_capteur_type(capteur%nom,"_depla",fnamef)

                open(fileId,file=trim(fnamef),status="replace",form="formatted")
                close(fileId)

            elseif (Tdomain%logicD%run_restart.and. rg==0) then ! c'est une reprise, il faut se repositionner
                ! au bon endroit dans le fichier de capteur pour le completer ensuite a chaque iteration

                if (capteur%grandeur.eq."DEFORMATION") then ! lecture du fichier de deformation

                    call semname_capteur_type(capteur%nom,"_deformation",fnamef)
                    open(fileId,file=trim(fnamef),status="old",form="formatted")
                    do i=1,Tdomain%TimeD%NtimeMin
                        read(fileId,*,END=998)
                    enddo
                    endfile(fileId)
998                 continue ! on a atteint la fin de fichier prematurement
                    close(fileId)
                endif

                if (capteur%grandeur.eq."VITESSE") then ! lecture du fichier de vitesse

                    call semname_capteur_type(capteur%nom,"_vitesse",fnamef)
                    !     print*,"capteur.f90, pt 2 - vitesse"    ,Tdomain%TimeD%NtimeMin, trim(fnamef)
                    !     print*,"capteur.f90, pt 2 - vitesse 2",fileId,fnamef,rg
                    open(fileId,file=trim(fnamef),status="old",form="formatted")
                    do i=1,Tdomain%TimeD%NtimeMin
                        read(fileId,*,END=999)
                    enddo
                    endfile(fileId)
999                 continue ! on a atteint la fin de fichier prematurement
                    close(fileId)
                endif

                if (capteur%grandeur.eq."DEPLA") then ! lecture du fichier de deplacements

                    call semname_capteur_type(capteur%nom,"_depla",fnamef)
                    !     print*,"capteur.f90, pt 2 - vitesse"    ,Tdomain%TimeD%NtimeMin, trim(fnamef)
                    !     print*,"capteur.f90, pt 2 - vitesse 2",fileId,fnamef,rg
                    open(fileId,file=trim(fnamef),status="old",form="formatted")
                    do i=1,Tdomain%TimeD%NtimeMin
                        read(fileId,*,END=997)
                    enddo
                    endfile(fileId)
997                 continue ! on a atteint la fin de fichier prematurement
                    close(fileId)
                endif

            endif ! fin du test run vs reprise

            cycle

100         exit  ! fin du fichier atteinte

        enddo

        !    print*,"avt fermeture capteur"
        close(fileIdCapteur)

        dimCapteur=dimCapteur-1


        return


201     write(*,*)"Erreur (label 201) de format frequence capteur",dimCapteur
        info_capteur = 0 !! ajout 03/10
        !!    stop  !!initial

202     write(*,*)"Erreur (label 202) de format rayon capteur",dimCapteur
        info_capteur = 0 !! ajout 03/10
        !!    stop !!initial

203     write(*,*)"Erreur (label 203) de format X capteur",dimCapteur
        info_capteur = 0 !! ajout 03/10
        !!    stop !!initial

204     write(*,*)"Erreur (label 204) de format Y capteur",dimCapteur
        info_capteur = 0 !! ajout 03/10
        !!    stop !!initial

205     write(*,*)"Erreur (label 205) de format Z capteur",dimCapteur
        info_capteur = 0 !! ajout 03/10
        !!    stop    !!initial

206     write(*,*)"Erreur (label 206) pour type de calcul du capteur",dimCapteur
        info_capteur = 0 !! ajout 03/10
        !!    stop  !!initial

    end subroutine lireCapteur

    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !>
    !! \brief Cree pour chaque capteur la liste du ou des points de Gauss à prendre en compte
    !!
    !>

    subroutine capterPointDeGauss(coord, npg)

        implicit none

        integer :: npg                  ! nombre de points de gauss
        real, dimension(3,npg) :: coord ! coordonnees des points de gauss

        integer :: i
        real    :: dist   ! distance capteur-PtGauss

        type(tCapteur),pointer :: capteur
        type(tPtgauss),pointer :: PtGauss

        ! boucle sur tous les points de Gauss
        do i=1,npg

            allocate (capteur)
            capteur=>listeCapteur

            ! boucle sur tous les capteurs
            do while (associated(capteur))


                ! distance entre le capteur et le point de Gauss
                dist = sqrt((coord(1,i) - capteur%coord(1))**2 + (coord(2,i)-capteur%coord(2))**2+ (coord(3,i)-capteur%coord(3))**2)

                if (capteur%rayon.GT.0.) then ! on teste si le point de Gauss est dans le rayon du capteur

                    if (dist.LE.capteur%rayon) then  ! le Pt Gauss est a prendre
                        allocate(PtGauss)
                        PtGauss%coord(:) = coord(:,i)
                        PtGauss%num = i-1
                        capteur%nPtGauss = capteur%nPtGauss +1
                        PtGauss%suivant=>capteur%listePtGauss
                        capteur%listePtGauss=>PtGauss
                    endif

                else ! on prend le pt de Gauss le plus proche du capteur

                    if (dist.lt.capteur%distanceMin) then ! le pt Gauss devient le pt le plus proche
                        capteur%distanceMin = dist
                        allocate(PtGauss)
                        PtGauss%coord(:) = coord(:,i)
                        PtGauss%num = i-1
                        capteur%nPtGauss = 1
                        nullify(PtGauss%suivant)
                        capteur%listePtGauss => PtGauss
                    endif

                endif

                capteur=>capteur%suivant

            enddo


        enddo


    end subroutine capterPointDeGauss

    !---------------------------------------------------------------------
    !---------------------------------------------------------------------
    !>
    !! \fn subroutine capterPointDeGauss(coord,npg)
    !! \brief
    !!
    !>
    !! \fn subroutine evalueSortieCapteur(it)
    !! \brief
    !!
    !! \param integer it
    !<


    subroutine evalueSortieCapteur(it, time, sortie_capteur)
        implicit none
        integer, intent(in) :: it
        real*8, intent(in) :: time
        logical, intent(out) :: sortie_capteur
        type(tCapteur),pointer :: capteur

        ! boucle sur les capteurs
        allocate (capteur)
        capteur=>listeCapteur

        do while (associated(capteur))
            if (mod(it,capteur%frequence)==0) then ! on fait la sortie
                sortie_capteur = .TRUE.
            endif

            capteur=>capteur%suivant
        enddo
    end subroutine evalueSortieCapteur

    !---------------------------------------------------------------------
    !---------------------------------------------------------------------

    !>
    !! \fn subroutine save_capteur(Tdomain)
    !! \brief
    !!
    !! \param type (domain) TDomain
    !<


    subroutine save_capteur(Tdomain, ntime, rg)


        implicit none

        integer :: rg, ntime
        type (domain) :: TDomain

        type(tCapteur),pointer :: capteur
        logical :: do_flush

        do_flush = .false.
        ! boucle sur les capteurs
        capteur=>listeCapteur
        do while (associated(capteur))

            if(capteur%type_calcul==0) then
                call sortirGrandeurSousCapteur(capteur,Tdomain, rg)
            else
                call sortieGrandeurCapteur_interp(capteur,Tdomain, ntime, rg)
            endif

            if (capteur%icache==NCAPT_CACHE) do_flush = .true.

            capteur=>capteur%suivant
        enddo

        if (do_flush) call flushAllCapteurs(Tdomain,rg)

    end subroutine save_capteur

    function dset_capteur_name(capteur)
        implicit none
        type(tCapteur),pointer :: capteur
        character(len=40) :: dset_capteur_name
        dset_capteur_name = trim(adjustl(capteur%nom))//"_"//trim(adjustl(capteur%grandeur))
    end function dset_capteur_name

    subroutine create_traces_h5_skel(Tdomain)
        implicit none
        type (domain) :: TDomain
        type(tCapteur),pointer :: capteur
        character (len=MAX_FILE_SIZE) :: fnamef
        character (len=40) :: dname
        integer(HID_T) :: fid, dset_id
        integer :: hdferr

        call init_hdf5()
        call semname_tracefile_h5(fnamef)

        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

        capteur=>listeCapteur
        do while (associated(capteur))
            dname = dset_capteur_name(capteur)
            write(*,*) "Create dset:", trim(adjustl(dname))
            call create_dset_2d(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, 4, H5S_UNLIMITED_F, dset_id)
            call h5dclose_f(dset_id, hdferr)
            capteur=>capteur%suivant
        enddo

        call h5fclose_f(fid, hdferr)
    end subroutine create_traces_h5_skel

    subroutine append_traces_h5(Tdomain)
        implicit none
        type (domain) :: TDomain
        type(tCapteur),pointer :: capteur
        character (len=40) :: dname
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: dset_id, fid
        integer :: hdferr

        call semname_tracefile_h5(fnamef)

        call h5fopen_f(fnamef, H5F_ACC_RDWR_F, fid, hdferr)

        capteur=>listeCapteur
        do while (associated(capteur))
            !write(*,*) "Capteur:", capteur%nom
            dname = dset_capteur_name(capteur)
            if (capteur%icache==0) then
                capteur=>capteur%suivant
                cycle
            endif
            call h5dopen_f(fid, trim(dname), dset_id, hdferr)
            call append_dataset_2d(dset_id, capteur%valuecache(:,1:capteur%icache), hdferr)
            call h5dclose_f(dset_id, hdferr)
            capteur%icache=0
            capteur=>capteur%suivant
        enddo

        call h5fclose_f(fid, hdferr)
    end subroutine append_traces_h5

    subroutine flushAllCapteurs(Tdomain, rg)
        implicit none
        integer :: rg
        type (domain) :: TDomain
        type(tCapteur),pointer :: capteur

        !write(*,*) "::::::::::::: TRACES FORMAT ::::", Tdomain%traces_format
        ! Default unspecified value is 'text'
        if (Tdomain%traces_format == 0) Tdomain%traces_format = 1

        if (Tdomain%traces_format == 1) then
            ! boucle sur les capteurs
            capteur=>listeCapteur
            do while (associated(capteur))
                call flushCapteur(capteur,Tdomain,rg)
                capteur=>capteur%suivant
            enddo
        else
            if (rg/=0) return
            ! Sauvegarde au format hdf5
            if (.not. traces_h5_created) then
                call create_traces_h5_skel(Tdomain)
                traces_h5_created = .true.
            end if
            call append_traces_h5(Tdomain)
        end if
    end subroutine flushAllCapteurs

    subroutine flushCapteur(capteur, Tdomain, rg)
        implicit none
        integer :: rg
        type (domain) :: TDomain
        type(tCapteur),pointer :: capteur
        integer :: fileId, j, imax
        character(len=MAX_FILE_SIZE) :: fnamef

        if (rg/=0) return

        if (capteur%icache==0) return

        if (trim(capteur%grandeur).eq."DEFORMATION") then
            call semname_capteur_type(capteur%nom,"_deformation",fnamef)
            imax = 1
        endif

        if (trim(capteur%grandeur).eq."VITESSE") then
            call semname_capteur_type(capteur%nom,"_vitesse",fnamef)
            imax = 3
        endif

        if (trim(capteur%grandeur).eq."DEPLA") then
            call semname_capteur_type(capteur%nom,"_depla",fnamef)
            imax = 3
        endif

        open(fileId,file=trim(fnamef),status="unknown",form="formatted",position="append")
        do j=1,capteur%icache
            if (imax==1) then
                write(fileId,'(2(1X,E16.8E3))') capteur%valuecache(1,j), capteur%valuecache(2,j)
            else
                write(fileId,'(4(1X,E16.8E3))') capteur%valuecache(1,j), capteur%valuecache(2:4,j)
            endif
        end do
        close(fileId)
        capteur%icache = 0
    end subroutine flushCapteur


    !---------------------------------------------------------------------
    !---------------------------------------------------------------------


    !>
    !! \fn subroutine sortirGrandeurSousCapteur(capteur,Tdomain)
    !! \brief
    !!
    !! \param type(tCapteur) capteur
    !! \param type (domain) TDomain
    !<


    subroutine sortirGrandeurSousCapteur(capteur,Tdomain, rg)

        implicit none

        type(tCapteur) :: capteur
        type (domain) :: TDomain

        integer :: rg
        integer :: fileId, i, ierr, nPtGaussTotal, tag, borneInf, iproc

        type(tPtGauss),pointer :: PtGauss

        integer , dimension  (MPI_STATUS_SIZE) :: status

        real, dimension(3) :: val0, val
        real, dimension(:,:), allocatable :: grandeur

        real,dimension(:),allocatable :: recvbuf
        real,dimension(4) :: sendbuf

        ! tableau de correspondance
        allocate(recvbuf(4*Tdomain%n_proc))
        sendbuf=0.
        recvbuf=0.

        ! ETAPE 0 : initialisations
        fileId=99 ! id du fichier des sorties capteur
        allocate (grandeur(3,max(1,capteur%nPtGauss))) ! allocation du tableau de la grandeur restreinte aux pts de Gauss captes
        grandeur(:,:) = 0.                               ! si maillage vide donc pas de pdg, on fait comme si il y en avait 1


        ! ETAPE 1 : Recuperation des valeurs de la grandeur pour les pts de Gauss definis par le capteur
        allocate(PtGauss)
        PtGauss=>capteur%listePtGauss
        i=0
        do while (associated(PtGauss))

            i=i+1

            if (trim(capteur%grandeur).eq."VITESSE") then
                grandeur(:,i) = grandeur_vitesse(Tdomain, PtGauss)
            endif

            if (trim(capteur%grandeur).eq."DEPLA") then
                grandeur(:,i) = grandeur_depla(Tdomain, PtGauss)
            endif

            PtGauss=>PtGauss%suivant ! passage au point de Gauss suivant

        enddo


        ! ETAPE 3 : Realisation de l operation stipulee par le capteur

        ! en parallele, la listeCapteur est propre a chaque proc; cela est Ok qd on cherche les points
        ! de gauss dans un rayon du capteur, chaque proc recupere sur son maillage les pdg compris dans
        ! le rayon, mais cela pose un pb dans le cas du pdg le plus proche. En effet, pour chaque proc
        ! il y a un pdg le plus proche du capteur. Reste donc a determiner le + proche de tous !

        if (capteur%rayon.eq.0) then ! on ne retient que la valeur du point le plus proche


            !cas du proc 0
            if (rg.eq.0) then

                ! on remplit recvbuf avec les infos du proc 0
                recvbuf(1) = capteur%distanceMin
                recvbuf(2) = grandeur(1,1)
                recvbuf(3) = grandeur(2,1)
                recvbuf(4) = grandeur(3,1)

                ! et ensuite avec les infos des autres proc
                do iproc=1,Tdomain%n_proc-1
                    borneInf = 4*iproc+1
                    tag = 100*(iproc+1)+capteur%numero

                    call mpi_recv(recvbuf(borneInf),4,MPI_DOUBLE_PRECISION,iproc,tag, Tdomain%communicateur,status,ierr)
                enddo


                ! les comm st bloquantes, ici, on a forcement le tableau recvbuf complet
                ! on recherche le pdg le + proche et ses valeurs associees dans recvbuf
                ! le resultat est stocke dans recvbuf(1:3)
                do iproc=1,Tdomain%n_proc-1
                    borneInf = 4*iproc+1
                    if (recvbuf(borneInf).lt.recvbuf(1)) then
                        recvbuf(1)=recvbuf(borneInf)
                        recvbuf(2)=recvbuf(borneInf+1)
                        recvbuf(3)=recvbuf(borneInf+2)
                        recvbuf(4)=recvbuf(borneInf+3)
                    endif
                enddo

                val0(1)=recvbuf(2)
                val0(2)=recvbuf(3)
                val0(3)=recvbuf(4)

            else ! cas des autres procs : on envoie des infos au proc 0

                sendbuf(1)=capteur%distanceMin
                sendbuf(2) = grandeur(1,1) ! grandeur en x au pdg le plus proche
                sendbuf(3) = grandeur(2,1) ! grandeur en y au pdg le plus proche
                sendbuf(4) = grandeur(3,1) ! grandeur en z au pdg le plus proche
                tag = 100*(rg + 1) + capteur%numero

                call mpi_send(sendbuf,4,MPI_DOUBLE_PRECISION,0,tag, Tdomain%communicateur,ierr)
            endif

        else ! on est dans le cas d un capteur a rayon d action

            ! CALCUL DE MOYENNE
            if (capteur%operation.eq."MOY") then
                val(1:3) = sum(grandeur,dim=2) ! on somme la grandeur sur les pdg du proc courant

                ! on fait la somme total des grandeurs de l'ensemble des procs qu'on ramene sur proc 0
                call MPI_Reduce (val,val0,3, MPI_DOUBLE_PRECISION, MPI_SUM, 0, Tdomain%communicateur, ierr)
                ! on recupere sur proc 0 le nombre total de capteurs
                call MPI_Reduce (capteur%nPtGauss,nPtGaussTotal,1, MPI_INTEGER, MPI_SUM, 0, Tdomain%communicateur, ierr)

                ! on moyenne sur le proc 0
                if(rg==0) then
                    val0(:)=val(:)/nPtGaussTotal
                endif

            endif

            ! CALCUL DU MAX
            if (capteur%operation.eq."MAX") then
                val(:) = maxval(grandeur,dim=2)
                call MPI_Reduce (val, val0,3, MPI_DOUBLE_PRECISION, MPI_MAX,0, Tdomain%communicateur, ierr)
            endif


            ! CALCUL DU MIN
            if (capteur%operation.eq."MIN") then
                val(:) = minval(grandeur,dim=2)
                call MPI_Reduce (val, val0,3, MPI_DOUBLE_PRECISION, MPI_MIN,0, Tdomain%communicateur, ierr)
            endif

        endif



        ! ETAPE 4 : Impression du resultat dans le fichier de sortie par le proc 0
        if(rg==0)then
            i = capteur%icache+1
            capteur%valuecache(1,i) = Tdomain%TimeD%rtime
            capteur%valuecache(2:4,i) = val0(1:3)
            capteur%icache = i
        endif

        deallocate(grandeur)
        deallocate(recvbuf)


    end subroutine sortirGrandeurSousCapteur



    !! effectue l'interpolation des grandeurs dans la maille dans laquelle se trouve le capteur
    !! la maille se trouve dans un seul proc
    !! seul le proc gere l'ecriture
    !!
    subroutine sortieGrandeurCapteur_interp(capteur,Tdomain, ntime, rg)

        implicit none

        type(tCapteur) :: capteur
        type (domain) :: TDomain
        integer :: rg

        integer :: fileId, i, j, k, idim, ierr, tag

        integer , dimension  (MPI_STATUS_SIZE) :: status
        integer request

        real, dimension(3) :: val0

        real, dimension(3) :: recvbuf
        real, dimension(3) :: sendbuf
        real, dimension(3) :: grandeur
        real, dimension(:,:,:,:), allocatable :: field
        integer n_el, ngllx, nglly, ngllz, mat
        real xi, eta, zeta
        real outx, outy, outz

        integer ntime
        integer numproc

        ! tableau de correspondance
        sendbuf=0.
        recvbuf=0.

        ! ETAPE 0 : initialisations
        fileId=99 ! id du fichier des sorties capteur
        grandeur(:)=0.                               ! si maillage vide donc pas de pdg, on fait comme si il y en avait 1


        ! Recuperation du numero de la maille Sem et des abscisses
        n_el = capteur%n_el
        xi = capteur%xi
        eta = capteur%eta
        zeta = capteur%zeta

        ! ETAPE 1 : interpolations

        numproc = -1
        if((n_el/=-1) .AND. (capteur%numproc==rg)) then
            ngllx = Tdomain%specel(n_el)%ngllx
            nglly = Tdomain%specel(n_el)%nglly
            ngllz = Tdomain%specel(n_el)%ngllz
            mat = Tdomain%specel(n_el)%mat_index
            allocate(field(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
            if (trim(capteur%grandeur).eq."VITESSE") then
                call gather_elem_veloc(Tdomain, n_el, field)
            else if (trim(capteur%grandeur).eq."DEPLA") then
                call gather_elem_displ(Tdomain, n_el, field)
            end if
            do i = 0,ngllx - 1
                do j = 0,nglly - 1
                    do k = 0,ngllz - 1
                        call  pol_lagrange (ngllx,Tdomain%sSubdomain(mat)%GLLcx,i,xi,outx)
                        call  pol_lagrange (nglly,Tdomain%sSubdomain(mat)%GLLcy,j,eta,outy)
                        call  pol_lagrange (ngllz,Tdomain%sSubdomain(mat)%GLLcz,k,zeta,outz)

                        grandeur(:) = grandeur(:) + outx*outy*outz*field(i,j,k,:)
                    enddo
                enddo
            enddo
            sendbuf(1) = grandeur(1) ! grandeur en x
            sendbuf(2) = grandeur(2) ! grandeur en y
            sendbuf(3) = grandeur(3) ! grandeur en z
            tag = 100*(rg + 1)+capteur%numero

            call mpi_Isend(sendbuf,3,MPI_DOUBLE_PRECISION,0,tag, Tdomain%communicateur,request, ierr)
            numproc = rg !!numero du proc courant

        endif

        !! Si le capteur est situe sur plusieurs processeurs, on choisit un proc et on envoie les resu au proc 0.

!!!!
!!!!    ! ETAPE 3 : Realisation de l operation stipulee par le capteur
!!!!
!!!!    ! en parallele, la listeCapteur est propre a chaque proc; cela est Ok qd on cherche les points
!!!!    ! de gauss dans un rayon du capteur, chaque proc recupere sur son maillage les pdg compris dans
!!!!    ! le rayon, mais cela pose un pb dans le cas du pdg le plus proche. En effet, pour chaque proc
!!!!    ! il y a un pdg le plus proche du capteur. Reste donc a determiner le + proche de tous !
!!!!
        !cas du proc 0
        if (rg.eq.0) then
            if((capteur%numproc)>-1) then

                ! on remplit recvbuf avec les infos du proc 0
                recvbuf(:)=grandeur(:)

                ! et ensuite avec les infos des autres proc
                tag = 100*(capteur%numproc+1) + capteur%numero

                call mpi_recv(recvbuf(1),3,MPI_DOUBLE_PRECISION,capteur%numproc,tag, Tdomain%communicateur,status,ierr)
                do idim=1,3
                    val0(idim)=recvbuf(idim)
                enddo
            else
                if(ntime<= 1) write(6,'(A,A,A,1X,I3)') 'Le capteur ',capteur%nom,' n''est pas pris en compte car sur aucun proc',capteur%numproc
            endif
        endif


        ! ETAPE 4 : Impression du resultat dans le fichier de sortie par le proc 0
        if((rg==0) .AND. (capteur%numproc>-1)) then
            i = capteur%icache+1
            capteur%valuecache(1,i) = Tdomain%TimeD%rtime
            capteur%valuecache(2:4,i) = val0(1:3)
            capteur%icache = i
        endif

    end subroutine sortieGrandeurCapteur_interp

    !!on identifie la maille dans laquelle se trouve le capteur. Il peut y en avoir plusieurs,
    !! alors le capteur est sur une face, arete ou sur un sommet
    !!
    subroutine trouve_capteur(Tdomain, rg, capteur, n_el, xi, eta, zeta)

        implicit none
        type (domain), INTENT(INOUT)  :: Tdomain
        integer rg
        type(Tcapteur) :: capteur
        integer ipoint
        integer i, n, n_el
        real dist
        integer P(0:7)
        real coor(0:7,0:2)
        real eps
        real xi, eta, zeta
        logical dans_maille
        integer i_sens
        eps=1.e-4

        n_el = -1
        i_sens = 0 !sens de parcours continu
        do n=0,Tdomain%n_elem-1
            ipoint = Tdomain%specel(n)%Control_Nodes(0) !!premier noeud de la maille n

            dist = (capteur%coord(1) - Tdomain%Coord_Nodes(0,ipoint))**2 + &
                (capteur%coord(2) - Tdomain%Coord_Nodes(1,ipoint))**2 + &
                (capteur%coord(3) - Tdomain%Coord_Nodes(2,ipoint))**2
            dist = sqrt(dist)

            if(dist> Tdomain%specel(n)%dist_max) then
                !          print*,'dist',dist,n,Tdomain%specel(n)%dist_max
                cycle
            else
                !maille retenue
                P(0:7) = Tdomain%specel(n)%Control_Nodes(0:7)

                do i=0,7
                    coor(i,0:2) = Tdomain%Coord_Nodes(0:2,P(i))
                enddo

                dans_maille = test_contour_capteur(coor, capteur, i_sens)

                if( dans_maille) then
                    !! on peut retenir cette maille et calculer xi, eta
                    !! si le capteur est sur une face ou un sommet, il suffit de garder une maille et de faire l'interpolation.
                    !! quel que soit l'element on obtiendra la bonne valeur
                    write(6,*) 'on retient la maille',n,capteur%coord(1:3)
                    call calc_xi_eta_zeta_capteur(capteur, coor, xi, eta, zeta)
                    n_el = n
                    capteur%numproc = rg
                    print *,'Arrivee a return'
                    return
                else
                    cycle

                endif
            endif
        enddo

        if(n_el/=-1) then
            write(6,*) 'xi,eta,zeta trouves pour le capteur ',capteur%nom,' de positions ',capteur%coord(1), capteur%coord(2), capteur%coord(3),&
                ' sont' ,  xi, eta, zeta, 'dans l element ',n_el,' pour le proc ',rg
        endif


    end subroutine trouve_capteur


    !!
    !! calcul des abscisses curvilignes du cpateur
    !! processus de dichotomie
    !!
    subroutine calc_xi_eta_zeta_capteur(capteur, coord, xitrouve_def, etatrouve_def, zetatrouve_def)

        implicit none
        type(Tcapteur) :: capteur
        real xi_min, xi_max, xi0, eta0, eta_min, eta_max, zeta0, zeta_min, zeta_max
        real P(8,0:2), dxi, deta, dzeta, coord(0:7,0:2)
        real xfact
        real xtrouve_def, ytrouve_def, ztrouve_def, xitrouve_def, etatrouve_def, zetatrouve_def
        real dist_def, dist, dist_P
        real xi, eta, zeta
        real eps
        integer i, j, k, im, il, in, idim, npts
        integer n_it
        logical interieur
        integer, parameter :: n_itmax=20
        integer i_sens
        !! attention si le point du capteur se trouve partage entre plusieurs elements (sommets, face ou aretes)
        !! identifiation

        xi_min = -1.  !bornes min et max de la zone d'etude dans le carre de reference
        xi_max = 1.   !on coupe en 2 dans chaque direction la zone d'etude
        eta_min = -1.
        eta_max = 1.
        zeta_min = -1.
        zeta_max = 1.
        dist_def = huge(1.) !distance entre la solution et le point de depart
        !!  eps = 1.e-6 !tolerance pour accepter la solution !!trop grand pour un cas
        eps = 1.e-6 !tolerance pour accepter la solution
        n_it = 0  !nombre d'iterations
        xi0 = xi_min
        eta0 = eta_min
        zeta0 = zeta_min
        i_sens = 1 !sens de parcours non continu

        !!  print *,'coord de la maille pour capteur',coord(0,:)
        !!  print *,'coord de la maille pour capteur',coord(1,:)
        !!  print *,'coord de la maille pour capteur',coord(2,:)
        !!  print *,'coord de la maille pour capteur',coord(3,:)
        !!  print *,'coord de la maille pour capteur',coord(4,:)
        !!  print *,'coord de la maille pour capteur',coord(5,:)
        !!  print *,'coord de la maille pour capteur',coord(6,:)
        !!  print *,'coord de la maille pour capteur',coord(7,:)
        !  do while(dist_def > eps) !precedemment
        do while((dist_def > eps) .AND. (n_it<n_itmax))
            n_it = n_it +1
            !on subdivise la zone en 4 sous-zones d'etude
            do in=0,1
                do im=0,1
                    do il=0,1
                        dxi = (xi_max - xi_min)/2.
                        deta = (eta_max - eta_min)/2.
                        dzeta = (zeta_max - zeta_min)/2.
                        !     on elargit un peu le domaine pour ne pas passer a cote de certains points critiques
                        xfact = 1.01
                        dxi   = xfact*dxi
                        deta  = xfact*deta
                        dzeta = xfact*dzeta
                        !
                        xi0 = xi_min + il*dxi
                        eta0 = eta_min + im*deta
                        zeta0 = zeta_min + in*dzeta
                        npts = 0
                        dist = huge(1.)
                        !Dans les boucles sur i, j et k, on definit les points P
                        do k=0,1
                            do j=0,1
                                do i=0,1
                                    xi = xi0 + i*dxi
                                    eta = eta0 + j*deta
                                    zeta = zeta0 + k*dzeta
                                    npts = npts + 1
                                    do idim=0,2
                                        P(npts,idim)=  0.125 * ( coord(0,idim)*(1.-xi)*(1.-eta)*(1.-zeta) + coord(1,idim)*(1.+xi)*(1.-eta)*(1.-zeta) + &
                                            coord(2,idim)*(1.+xi)*(1.+eta)*(1.-zeta) + coord(3,idim)*(1.-xi)*(1.+eta)*(1.-zeta)  + &
                                            coord(4,idim)*(1.-xi)*(1.-eta)*(1.+zeta) + coord(5,idim)*(1.+xi)*(1.-eta)*(1.+zeta) + &
                                            coord(6,idim)*(1.+xi)*(1.+eta)*(1.+zeta) + coord(7,idim)*(1.-xi)*(1.+eta)*(1.+zeta) )
                                    enddo

                                    dist_P = 0.
                                    do idim=0,2
                                        dist_P = dist_P + (P(npts,idim) - capteur%coord(idim+1))**2
                                    enddo
                                    dist_P = sqrt(dist_P)
                                    !                     print*,' dist_P ',dist_P

                                    dist = min(dist, dist_P)
                                    !
                                    !on teste si P(npts) correspond au capteur. Si oui on sort
                                    !
                                    if(dist_P < eps) then
                                        xtrouve_def = P(npts,0)
                                        ytrouve_def = P(npts,1)
                                        ztrouve_def = P(npts,2)
                                        !    modif complementaire a xfact
                                        if ( xi < -1. ) xi = -1.
                                        if ( xi > 1. )  xi = 1.
                                        if ( eta < -1. ) eta = -1.
                                        if ( eta > 1. )  eta = 1.
                                        if ( zeta < -1. ) zeta = -1.
                                        if ( zeta > 1. )  zeta = 1.
                                        !   fin modif
                                        xitrouve_def = xi
                                        etatrouve_def = eta
                                        zetatrouve_def = zeta
                                        print *,'final xtrouve ytrouve ztrouve',xtrouve_def,ytrouve_def,ztrouve_def, &
                                            xitrouve_def,etatrouve_def,zetatrouve_def,' iteration ',n_it,' distance ',dist_P
                                        return
                                    endif
                                enddo
                            enddo
                        enddo
                        !! on teste si le capteur est a l'interieur du contour forme des points P
                        !! si oui on conserve xi, eta et on reduit la region
                        interieur = test_contour_capteur(P, capteur, i_sens)
                        if(interieur) then
                            xi_min = xi0
                            xi_max = xi0 + dxi
                            eta_min = eta0
                            eta_max = eta0 + deta
                            zeta_min = zeta0
                            zeta_max = zeta0 + dzeta
                            dist_def = dist_P
                        endif
                    enddo
                enddo
            enddo
        enddo

    end subroutine calc_xi_eta_zeta_capteur


    logical function test_contour_capteur(P, capteur, i_sens)

!!!!  Description: on teste si le point du capteur est a l'interieur du volume defini par les points P.
!!!!  Si le capteur est sur une face, une arete ou sur un sommet, renvoie true
!!!!  * attention : le sens de parcours dans la face n'est pas continu
!!!!
!!!!  Historique: 03/10 - Gsa Ipsis - Creation
!!!! --------------------------------------------------

        implicit none
        type(Tcapteur) :: capteur
        real P(0:7,0:2)
        real coord(0:3,0:2)
        real, parameter :: eps = 1.e-4
        integer i, idim, iface
        integer, parameter :: n_faces=6
        integer, parameter :: n_nodes=8
        real, dimension(0:2) :: pg, vn, v
        integer i_sens
        integer,dimension(0:5,0:3) :: NOEUD_FACE !Gsa

        test_contour_capteur = .true.

        !tableau de correspondance num_local de face/numero des noeuds de l'element

        if (i_sens==0) then
            NOEUD_FACE(0,:) = (/ 0, 1, 2, 3 /) !face 0 - k=0
            NOEUD_FACE(1,:) = (/ 0, 1, 5, 4 /) !face 1 - j=0
            NOEUD_FACE(2,:) = (/ 1, 2, 6, 5 /) !face 2 - i=n-1
            NOEUD_FACE(3,:) = (/ 3, 2, 6, 7 /) !face 3 - j=n-1
            NOEUD_FACE(4,:) = (/ 0, 3, 7, 4 /) !face 4 - i=0
            NOEUD_FACE(5,:) = (/ 4, 5, 6, 7 /) !face 5 - k=n-1
        elseif (i_sens==1) then
            NOEUD_FACE(0,:) = (/ 0, 1, 3, 2 /) !face 0 - k=0
            NOEUD_FACE(1,:) = (/ 0, 1, 5, 4 /) !face 1 - j=0
            NOEUD_FACE(2,:) = (/ 2, 0, 4, 6 /) !face 2 - i=n-1
            NOEUD_FACE(3,:) = (/ 2, 3, 7, 6 /) !face 3 - j=n-1
            NOEUD_FACE(4,:) = (/ 3, 1, 5, 7 /) !face 4 - i=0
            NOEUD_FACE(5,:) = (/ 4, 5, 7, 6 /) !face 5 - k=n-1
        endif

        pg = 0.  ! centre de gravite de la maille
        do i=0,n_nodes-1
            do idim=0,2
                pg(idim) = pg(idim) + P(i,idim)
            enddo
        enddo
        pg = pg/n_nodes

!!!!  write(6,*) 'Barycentre' ,pg
        do iface=0,n_faces-1
            do i=0,3
                do idim=0,2
                    coord(i,idim) = P(NOEUD_FACE(iface,i),idim)
                enddo
            enddo

            !calcul de la normale de la face
            vn(0) = (coord(1,1) - coord(0,1) ) * (coord(3,2) - coord(0,2) ) &
                - (coord(1,2) - coord(0,2) ) * (coord(3,1) - coord(0,1) )
            vn(1) = (coord(1,2) - coord(0,2) ) * (coord(3,0) - coord(0,0) ) &
                - (coord(1,0) - coord(0,0) ) * (coord(3,2) - coord(0,2) )
            vn(2) = (coord(1,0) - coord(0,0) ) * (coord(3,1) - coord(0,1) ) &
                - (coord(1,1) - coord(0,1) ) * (coord(3,0) - coord(0,0) )

            !verification de l'orientation de la normale de la face
            if(dot_product(vn, pg - coord(0,:)) > 0.) vn = -vn
            !normalisation de la normale
            if(abs(dot_product(vn,vn)) < 1.e-30) then !! attention on eleve au carre - la precision est petite
                print *,'Pb normale de face nulle - Arret',vn
                write(50,*) 'Pb normale de face nulle - Arret'
                print*,' isens ',i_sens
                print*,' capteur coord ',capteur%coord(1),capteur%coord(2),capteur%coord(3)
                print*,' noeudelement0 ',P(0,0),P(0,1),P(0,2)
                print*,' noeudelement1 ',P(1,0),P(1,1),P(1,2)
                print*,' noeudelement2 ',P(2,0),P(2,1),P(2,2)
                print*,' noeudelement3 ',P(3,0),P(3,1),P(3,2)
                print*,' noeudelement4 ',P(4,0),P(4,1),P(4,2)
                print*,' noeudelement5 ',P(5,0),P(5,1),P(5,2)
                print*,' noeudelement6 ',P(6,0),P(6,1),P(6,2)
                print*,' noeudelement7 ',P(7,0),P(7,1),P(7,2)
                print*,' noeudface0 ',coord(0,0),coord(0,1),coord(0,2)
                print*,' noeudface1 ',coord(1,0),coord(1,1),coord(1,2)
                print*,' noeudface2 ',coord(2,0),coord(2,1),coord(2,2)
                print*,' noeudface3 ',coord(3,0),coord(3,1),coord(3,2)
                print*,' vn ',vn(0),vn(1),vn(2)
                stop
            else
                vn = vn / sqrt(dot_product(vn,vn))
            endif

            do idim=0,2
                v(idim) = capteur%coord(idim+1) - coord(0,idim)
            enddo

            if(dot_product(vn, v)>0.) then
                test_contour_capteur = .false.
                !!        print *,'test faux'
                return
            endif
        enddo

        !!  print *,'test vrai'
        return
    end function test_contour_capteur



end module Mcapteur
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
