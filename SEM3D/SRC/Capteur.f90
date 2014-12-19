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
    use constants, only : NCAPT_CACHE
    implicit none

    public :: read_capteur, save_capteur, evalueSortieCapteur, flushAllCapteurs
    private :: lireCapteur, capterPointDeGauss, sortirGrandeurSousCapteur, flushCapteur


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
    subroutine read_capteur(Tdomain, info_capteur)

        implicit none

        type (domain), target  :: Tdomain
        integer :: info_capteur
        !
        type(Tcapteur),pointer :: capteur
        type(tPtgauss),pointer :: PtGauss

        integer :: rg
        integer :: ierr, tag
        integer , dimension  (MPI_STATUS_SIZE) :: status
        real :: distanceMinMin
        real,dimension(3) :: val, val0
        real,dimension(:),allocatable :: sendbuf, recvbuf
        integer numproc_max
        real xi, eta, zeta
        integer n_el

        rg = Tdomain%rank
        ! lecture des parametres des capteurs

        !lecture du fichier des capteurs capteurs.dat
        call lireCapteur(Tdomain, info_capteur)

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


            if(capteur%type_calcul==0) then

                PtGauss=>capteur%listePtGauss

                if (capteur%rayon==0) then ! un seul point de Gauss a ecrire, mais le plus proche de tous
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

                endif

            elseif(capteur%type_calcul==1) then
                xi = -1.
                eta = -1.
                zeta  = -1.
                n_el = -1
                call trouve_capteur(Tdomain, capteur, n_el, xi, eta, zeta)
                call MPI_AllReduce(capteur%numproc, numproc_max, 1, MPI_INTEGER, &
                    MPI_MAX, Tdomain%communicateur, ierr)

                ! attention si le capteur est partage par plusieurs procs. On choisit le proc de num max
                capteur%numproc = numproc_max
                if(rg==numproc_max) then
                    write(*,"(A,A,A,I5,A,I6,A,F8.4,A,F8.4,A,F8.4)") "Capteur", trim(capteur%nom), &
                        " on proc ", rg, " in elem ", n_el, " at ", xi, ",", eta, ",", zeta
                    capteur%n_el = n_el ! uniquement pour l'affichage dans le fichier *position
                    capteur%xi = xi
                    capteur%eta = eta
                    capteur%zeta = zeta

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
    subroutine lireCapteur(Tdomain, info_capteur)
        implicit none
        type (Domain), intent (IN) :: Tdomain
        integer :: info_capteur
        !
        integer            :: CodeErreur,fileId,rg
        character(LEN=MAX_FILE_SIZE) :: ligne,fnamef
        logical            :: status

        type(Tcapteur),pointer :: capteur
        character(LEN=6) :: char_calcul

        rg = Tdomain%rank
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
            if(Tdomain%rank==0) &
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
                ! Tests:
                if (capteur%frequence<1) then
                    write(*,*) "ERREUR: la periode de sous-echantillonage du capteur est nulle"
                    write(*,*) "        Il faut entrer un nombre d'iteration >0 pour le parametre FREQ"
                    stop 1
                endif
            end if

            capteur%suivant=>listeCapteur
            listeCapteur=>capteur
            ! si c'est un nouveau run, suppression de l'eventuel fichier de sortie des capteurs
            if (Tdomain%traces_format == 1) then
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

                    call semname_capteur_type(capteur%nom,"_accel",fnamef)

                    open(fileId,file=trim(fnamef),status="replace",form="formatted")
                    close(fileId)

                elseif (Tdomain%logicD%run_restart.and. rg==0) then
                    ! c'est une reprise, il faut se repositionner
                endif ! fin du test run vs reprise
            endif

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
        double precision, intent(in) :: time
        logical, intent(out) :: sortie_capteur
        type(tCapteur),pointer :: capteur

        ! boucle sur les capteurs
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


    subroutine save_capteur(Tdomain, ntime)


        implicit none

        integer :: ntime
        type (domain) :: TDomain

        type(tCapteur),pointer :: capteur
        logical :: do_flush

        do_flush = .false.
        ! boucle sur les capteurs
        capteur=>listeCapteur
        do while (associated(capteur))

            if(capteur%type_calcul==0) then
                call sortirGrandeurSousCapteur(capteur,Tdomain)
            else
                call sortieGrandeurCapteur_interp(capteur,Tdomain, ntime)
            endif

            if (capteur%icache==NCAPT_CACHE) do_flush = .true.

            capteur=>capteur%suivant
        enddo

        if (do_flush) call flushAllCapteurs(Tdomain)

    end subroutine save_capteur

    function dset_capteur_name(capteur)
        implicit none
        type(tCapteur),pointer :: capteur
        character(len=40) :: dset_capteur_name
        dset_capteur_name = trim(adjustl(capteur%nom))//"_"//trim(adjustl(capteur%grandeur))
    end function dset_capteur_name

    subroutine create_traces_h5_skel()
        use HDF5
        implicit none
        !type (domain) :: TDomain
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
            call create_dset_2d(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, int(4,HSIZE_T), int(H5S_UNLIMITED_F,HSIZE_T), dset_id)
            call h5dclose_f(dset_id, hdferr)
            capteur=>capteur%suivant
        enddo

        call h5fclose_f(fid, hdferr)
    end subroutine create_traces_h5_skel

    subroutine append_traces_h5()
        implicit none
        !type (domain) :: TDomain
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

    subroutine flushAllCapteurs(Tdomain)
        implicit none
        type (domain) :: TDomain
        type(tCapteur),pointer :: capteur

        !write(*,*) "::::::::::::: TRACES FORMAT ::::", Tdomain%traces_format
        ! Default unspecified value is 'text'
        if (Tdomain%traces_format == 0) Tdomain%traces_format = 1

        if (Tdomain%traces_format == 1) then
            ! boucle sur les capteurs
            capteur=>listeCapteur
            do while (associated(capteur))
                call flushCapteur(capteur,Tdomain%rank)
                capteur=>capteur%suivant
            enddo
        else
            if (Tdomain%rank/=0) return
            ! Sauvegarde au format hdf5
            if (.not. traces_h5_created) then
                call create_traces_h5_skel()
                traces_h5_created = .true.
            end if
            call append_traces_h5()
        end if
    end subroutine flushAllCapteurs

    subroutine flushCapteur(capteur, rg)
        implicit none
        integer, intent(in) :: rg
        type(tCapteur),pointer :: capteur
        !
        integer, parameter :: fileId=123
        integer :: j, imax
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

        if (trim(capteur%grandeur).eq."ACCEL") then
            call semname_capteur_type(capteur%nom,"_accel",fnamef)
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


    subroutine sortirGrandeurSousCapteur(capteur,Tdomain)

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

        rg = Tdomain%rank
        ! tableau de correspondance
        allocate(recvbuf(4*Tdomain%nb_procs))
        sendbuf=0.
        recvbuf=0.

        ! ETAPE 0 : initialisations
        fileId=99 ! id du fichier des sorties capteur
        allocate (grandeur(3,max(1,capteur%nPtGauss))) ! allocation du tableau de la grandeur restreinte aux pts de Gauss captes
        grandeur(:,:) = 0.                               ! si maillage vide donc pas de pdg, on fait comme si il y en avait 1


        ! ETAPE 1 : Recuperation des valeurs de la grandeur pour les pts de Gauss definis par le capteur
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
                do iproc=1,Tdomain%nb_procs-1
                    borneInf = 4*iproc+1
                    tag = 100*(iproc+1)+capteur%numero

                    call mpi_recv(recvbuf(borneInf),4,MPI_DOUBLE_PRECISION,iproc,tag, Tdomain%communicateur,status,ierr)
                enddo


                ! les comm st bloquantes, ici, on a forcement le tableau recvbuf complet
                ! on recherche le pdg le + proche et ses valeurs associees dans recvbuf
                ! le resultat est stocke dans recvbuf(1:3)
                do iproc=1,Tdomain%nb_procs-1
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
    subroutine sortieGrandeurCapteur_interp(capteur,Tdomain, ntime)
        use mpi
        implicit none

        type(tCapteur) :: capteur
        type (domain) :: TDomain
        integer :: rg

        integer :: i, j, k, ierr, tag

        integer , dimension  (MPI_STATUS_SIZE) :: status
        integer request

        real, dimension(3) :: grandeur
        real, dimension(:,:,:,:), allocatable :: field

        integer n_el, ngllx, nglly, ngllz, mat
        real xi, eta, zeta
        real outx, outy, outz

        integer ntime
        integer numproc

        rg = Tdomain%rank
        ! tableau de correspondance
        request=MPI_REQUEST_NULL

        ! ETAPE 0 : initialisations
        grandeur(:)=0. ! si maillage vide donc pas de pdg, on fait comme si il y en avait 1


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
            else if (trim(capteur%grandeur).eq."ACCEL") then
                call gather_elem_accel(Tdomain, n_el, field)
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
            tag = 100*(rg + 1)+capteur%numero

            if (rg/=0) call mpi_Isend(grandeur, 3,MPI_DOUBLE_PRECISION, 0, tag, Tdomain%communicateur, request, ierr)
            numproc = rg !!numero du proc courant
            deallocate(field)
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

                ! et ensuite avec les infos des autres proc
                tag = 100*(capteur%numproc+1) + capteur%numero

                if (capteur%numproc/=0) call mpi_recv(grandeur,3,MPI_DOUBLE_PRECISION,capteur%numproc,tag, Tdomain%communicateur,status,ierr)

                ! ETAPE 4 : Impression du resultat dans le fichier de sortie par le proc 0
                i = capteur%icache+1
                capteur%valuecache(1,i) = Tdomain%TimeD%rtime
                capteur%valuecache(2:4,i) = grandeur(1:3)
                capteur%icache = i
            else
                if(ntime<= 1) then
                    write(6,'(A,A,A,1X,I3)') 'Le capteur ',capteur%nom, &
                        ' n''est pas pris en compte car sur aucun proc',capteur%numproc
                endif
            endif

        endif

        if (request/=MPI_REQUEST_NULL) call MPI_Wait(request, status, ierr)

    end subroutine sortieGrandeurCapteur_interp

    !!on identifie la maille dans laquelle se trouve le capteur. Il peut y en avoir plusieurs,
    !! alors le capteur est sur une face, arete ou sur un sommet
    !!
    subroutine trouve_capteur(Tdomain, capteur, n_el, xi, eta, zeta)
        use mshape8
        use mshape27
        use mlocations3d
        implicit none
        type (domain), INTENT(INOUT)  :: Tdomain
        type(Tcapteur), intent(inout) :: capteur
        integer, intent(out) :: n_el
        double precision, intent(out) :: xi, eta, zeta
        !
        integer :: i
        double precision :: xc, yc, zc
        logical :: inside
        integer :: nmax
        integer, parameter :: NMAXEL=20
        integer, dimension(NMAXEL) :: elems
        double precision, dimension(0:2,NMAXEL) :: coordloc

        xc = capteur%coord(1)
        yc = capteur%coord(2)
        zc = capteur%coord(3)
        nmax = NMAXEL
        call find_location(Tdomain, xc, yc, zc, nmax, elems, coordloc)
        do i=1,nmax
            inside = .true.
            xi   = coordloc(0,i)
            eta  = coordloc(1,i)
            zeta = coordloc(2,i)
            if (xi<-1 .or. eta<-1 .or. zeta<-1) inside = .false.
            if (xi>1 .or. eta>1 .or. zeta>1) inside = .false.
            if (inside) then
                n_el = elems(i)
                capteur%numproc = Tdomain%rank
                return
            end if
        end do
        n_el = -1
        return
    end subroutine trouve_capteur

end module Mcapteur
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
