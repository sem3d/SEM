!>
!!\file Capteur.F90
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

    implicit none

    public :: read_capteur,save_capteur,evalueSortieCapteur
    private :: lireCapteur,capterPointDeGauss,sortirGrandeurSousCapteur

    type :: tPtgauss
       integer :: num      ! numero du point de Gauss
       real    :: Coord(2) ! coordonnees du point de gauss
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
       integer :: type_calcul           ! 0 => pas d'interpolation, recherche du point de Gauss le plus proche,
       ! 1 => interpolation aux positions du capteur
       !   (necessite le calcul des abscisses curvilignes dans l'element de reference)
       integer :: n_el ! numero de la maille dans laquelle se trouve le capteur
       ! si le capteur est partage entre plusieurs mailles, une seule suffit (type_calcul=1)
       real :: xi, eta ! abscisses curvilignes pour le capteur en cas d'interpolation (type_calcul=1)
       integer :: numproc               ! numero du proc localisant le capteur
    end type tCapteur


    integer           :: dimCapteur        ! nombre total de capteurs

    type(Tcapteur), pointer :: listeCapteur

    integer,parameter  :: fileIdCapteur=200  ! id fichier capteur
    logical            :: sortie_capteur
    logical            :: sortie_capteur_deformation

contains


    !>
    !! \brief La routine read_capteur() permet de charger les caractéristiques des capteurs définies dans les fichiers d'entrée.
    !! Elle est appelée une seule fois, avant la boucle en temps, afin d'initialiser les capteurs.
    !! \param type (domain), target Tdomain
    !<

    subroutine read_capteur(tDomain, info_capteur)

        implicit none

        type (domain), target  :: Tdomain
        type(Tcapteur),pointer :: capteur
        type(tPtgauss),pointer :: PtGauss
        integer,parameter :: fid=98
        character(len=MAX_FILE_SIZE) :: fnamef

        integer :: ierr,tag
        integer , dimension  (MPI_STATUS_SIZE) :: status
        real :: distanceMinMin
        real,dimension(2) :: val,val0
        real,dimension(:),allocatable :: sendbuf,recvbuf
        integer info_capteur, numproc_max
        real xi, eta
        integer n_el
        ! lecture des parametres des capteurs
        sortie_capteur = .FALSE.



        !!  call lireCapteur(tDomain)  !!initial
        call lireCapteur(tDomain, info_capteur)
        if(info_capteur /= 0) return

        ! creer pour chaque capteur la liste du ou des points de gauss a prendre en compte
        call capterPointDeGauss(Tdomain%GlobCoord,Tdomain%n_glob_points)


        ! tableau de correspondance pour les com
        allocate(sendbuf(2*dimCapteur))
        allocate(recvbuf(2*dimCapteur))


        call dist_max_elem(Tdomain) !!ajout calcul de la distance max entre deux sommets a l'interieur d'une maille
        !!permet de limiter le nombre de mailles pouvant contenir un capteur


        ! afficher position des points de Gauss pour chaque capteur ?
        allocate(capteur)
        capteur=>listeCapteur
        do while (associated(capteur))

            tag=capteur%numero ! une com par capteur

            call semname_capteur_pos(capteur%nom,fnamef)


            if(capteur%type_calcul==0) then

                open(UNIT=fid,FILE=trim(fnamef),STATUS='replace')
                allocate(PtGauss)
                PtGauss=>capteur%listePtGauss

                if (capteur%rayon==0) then ! un seul point de Gauss a ecrire, mais le plus proche de tous

                    ! tous les proc recuperent la distance min du pdg le + proche
                    call MPI_AllReduce(capteur%distanceMin,distanceMinMin,1,MPI_DOUBLE_PRECISION,MPI_MIN,Tdomain%communicateur,ierr)
                    if (distanceMinMin == capteur%distanceMin) then ! on est sur le proc qui a le pt de gauss le + proche
                        val(1)=PtGauss%coord(1) ! coord en x du pdg le plus proche
                        val(2)=PtGauss%coord(2) ! coord en y du pdg le plus proche
                        write(6,'(A30,I8,A13,I8)')"numero du pdg le plus proche:",PtGauss%num," sur le proc ",Tdomain%MPI_var%my_rank
                        if (Tdomain%MPI_var%my_rank.eq.0 ) then ! pdg le plus proche est sur le proc 0
                            val0(1)=val(1)
                            val0(2)=val(2)
                        else ! on envoie le pdg le plus proche au proc 0
                            sendbuf(2*(dimCapteur-1)+1)=val(1)
                            sendbuf(2*(dimCapteur-1)+2)=val(2)
                            call MPI_SEND(sendbuf,2*dimCapteur,MPI_DOUBLE_PRECISION,0,tag, Tdomain%communicateur,ierr)
                        endif
                    else
                        if (Tdomain%MPI_var%my_rank.eq.0 ) then ! proc 0 n'a pas le pdg le + proche et doit donc le recevoir

                            call mpi_recv(recvbuf,2*dimCapteur,MPI_DOUBLE_PRECISION,mpi_any_source,tag, Tdomain%communicateur,status,ierr)

                            val0(1)=recvbuf(2*(dimCapteur-1)+1)
                            val0(2)=recvbuf(2*(dimCapteur-1)+2)
                        endif
                    endif

                    if (Tdomain%MPI_var%my_rank.eq.0 ) then
                        rewind(fid)
                        write(fid,'(A2,E12.4,A3,E12.4,A6,E12.4)') "X=",val0(1)," Y=",val0(2)," dist=",distanceMinMin
                    endif

                else ! les proc ecrivent tour a tour les pt de gauss

                    ! on parcourt les pdg
                    do while (associated(PtGauss))
                        write(fid,*) PtGauss%coord(1),PtGauss%coord(2)
                        PtGauss=>PtGauss%suivant
                    enddo

                endif

                !deallocate(PtGauss)

                close(fid)
            elseif(capteur%type_calcul==1) then
                xi = -1.
                eta = -1.
                n_el = -1
                call trouve_capteur(Tdomain, capteur, n_el, xi, eta) !!Gsa 03/10 - interpolation du capteur

                if (Tdomain%MPI_var%my_rank.eq.0) open(UNIT=fid,FILE=trim(fnamef),STATUS='unknown')  !!modif seul proc 0
                val(1) = xi
                val(2) = eta
                call MPI_AllReduce(capteur%numproc, numproc_max, 1, MPI_INTEGER, MPI_MAX, Tdomain%communicateur, ierr)

                ! attention si le capteur est partage par plusieurs procs. On choisit le proc de num max
                capteur%numproc = numproc_max
                if(Tdomain%MPI_var%my_rank==numproc_max) then
                    capteur%n_el = n_el           !!uniquement pour l'affichage dans le fichier *position
                    capteur%xi = xi !val0(1)      !!uniquement pour l'affichage dans le fichier *position
                    capteur%eta = eta !val0(2)    !!uniquement pour l'affichage dans le fichier *position
                endif

                rewind(fid)
                write(fid,'(2(A6,1X,1pe15.8,1X),A22,1X,I6)') "xi=", capteur%xi, " eta=", capteur%eta, " numero element Sem", capteur%n_el
                close(fid)
                !          !! on reattribue au proc son numero d'element initial
                !          capteur%n_el = n_el
                !          capteur%xi = xi
                !          capteur%eta = eta
            endif
            capteur=>capteur%suivant


            !!    close(fid) !!initial
        enddo


        deallocate(sendbuf)
        deallocate(recvbuf)

        ! creation des tableaux de grandeur associees aux capteurs
        allocate(Tdomain%GrandeurDeformation(0:max(Tdomain%n_glob_points-1,0)))
        allocate(Tdomain%GrandeurVitesse(0:1,0:max(Tdomain%n_glob_points-1,0)))
        allocate(Tdomain%GrandeurDepla(0:1,0:max(Tdomain%n_glob_points-1,0)))



    end subroutine read_capteur


    !>
    !! \brief Lit les caratéristiques des capteurs.
    !!
    !! \param type (Domain), intent (IN) Tdomain
    !<
    subroutine lireCapteur(Tdomain,info_capteur)
        implicit none


        type (Domain), intent (IN) :: Tdomain

        integer            :: CodeErreur,fileId,i
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

	call semname_read_capteurs(trim(Tdomain%station_file),fnamef)
        !controle d'existence du fichier
        INQUIRE(File=trim(fnamef),Exist=status)
        info_capteur=0 !!Gsa
        if ( .not.status ) then
            write (*,*)"fichier introuvable :",trim(fnamef)
            !!      stop
            info_capteur=1
            return
        endif

        open(UNIT=fileIdCapteur,IOSTAT=CodeErreur,FILE=trim(fnamef),FORM='formatted',STATUS='old',ACTION='read')
        if (CodeErreur .ne.0 ) print*,'Ouverture du fichier Capteur  :CodeErreur=',CodeErreur


	
        do ! while (.not.eof(fileIdCapteur))


            ! initialisation d un nouveau capteur
            dimCapteur=dimCapteur+1   ! nbre total de capteur
            allocate(capteur)
            capteur%nom=""         ! ses caracteristiques par defaut
            capteur%numero = dimCapteur ! numero du capteur
            capteur%frequence=0
            capteur%rayon=0
            capteur%operation=""
            capteur%grandeur=""
            capteur%Coord(:)=-1
            capteur%n_el = -1
            capteur%numproc = -1
            capteur%type_calcul = 1 !mode de calcul du capteur
            nullify(capteur%listePtGauss)
            capteur%distanceMin = huge(1.)
            capteur%nPtGauss = 0  ! initialisation du nombre de points de Gauss a prendre en compte

            read(fileIdCapteur,'(A)',end=100)ligne ! titre 1
            read(fileIdCapteur,'(A)',end=100)ligne ! titre 2




            ! recuperation du nom du capteur
            read(fileIdCapteur,'(A)',end=100)ligne


            do while (trim(ligne).ne."")


                if (ligne(1:11).eq."NOM_CAPTEUR") capteur%nom = trim(ligne(12:100))

                if (ligne(1:4).eq."FREQ") read(ligne(5:100),*,ERR=201)capteur%frequence

                if (ligne(1:5).eq."RAYON") read(ligne(6:100),*,ERR=202)capteur%rayon

                if (ligne(1:9).eq."OPERATION") capteur%operation=trim(adjustl(ligne(10:100)))

                if (ligne(1:8).eq."GRANDEUR") capteur%grandeur=trim(adjustl(ligne(9:100)))

                if (ligne(1:6).eq."COORDX") read(ligne(7:100),*,ERR=203)capteur%Coord(1)

                if (ligne(1:6).eq."COORDY") read(ligne(7:100),*,ERR=204)capteur%Coord(2)

                if (ligne(1:6).eq."COORDZ") read(ligne(7:100),*,ERR=205)capteur%Coord(3)

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
                read(fileIdCapteur,'(A)',end=50)ligne

            enddo ! fin de lecture du capteur courant


50          if (Tdomain%MPI_var%my_rank==0) write(*,'(A7,I4,A6)')"capteur",dimCapteur,"OK"
            capteur%suivant=>listeCapteur
            listeCapteur=>capteur


            ! si c'est un nouveau run, suppression de l'eventuel fichier de sortie des capteurs
            if (.not.Tdomain%logicD%run_restart.and.Tdomain%MPI_var%my_rank==0) then

                call semname_capteur_type(capteur%nom,"_deformation",fnamef)
                fileId=99
                open(fileId,file=trim(fnamef),status="replace",form="formatted")
                close(fileId)

                call semname_capteur_type(capteur%nom,"_vitesse",fnamef)
                open(fileId,file=trim(fnamef),status="replace",form="formatted")
                close(fileId)

                call semname_capteur_type(capteur%nom,"_depla",fnamef)
                open(fileId,file=trim(fnamef),status="replace",form="formatted")
                close(fileId)

            elseif (Tdomain%logicD%run_restart.and.Tdomain%MPI_var%my_rank==0) then ! c'est une reprise, il faut se repositionner
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
        !!    stop  !!initial

206     write(*,*)"Erreur (label 206) pour type de calcul du capteur",dimCapteur
        info_capteur = 0 !! ajout 03/10
        !!    stop  !!initial


    end subroutine lireCapteur

    !>
    !! \brief Cree pour chaque capteur la liste du ou des points de Gauss à prendre en compte
    !!
    !>

    subroutine capterPointDeGauss(coord,npg)

        implicit none

        integer :: npg                  ! nombre de points de gauss
        real, dimension(2,npg) :: coord ! coordonnees des points de gauss

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
                dist = sqrt((coord(1,i)-capteur%coord(1))**2 + (coord(2,i)-capteur%coord(2))**2)

                if (capteur%rayon.gt.0) then ! on teste si le point de Gauss est dans le rayon du capteur

                    if (dist.le.capteur%rayon) then  ! le Pt Gauss est a prendre
                        allocate(PtGauss)
                        PtGauss%coord(1) = coord(1,i)
                        PtGauss%coord(2) = coord(2,i)
                        PtGauss%num = i-1
                        capteur%nPtGauss = capteur%nPtGauss +1
                        PtGauss%suivant=>capteur%listePtGauss
                        capteur%listePtGauss=>PtGauss
                    endif

                else ! on prend le pt de Gauss le plus proche du capteur

                    if (dist.lt.capteur%distanceMin) then ! le pt Gauss devient le pt le plus proche
                        capteur%distanceMin = dist
                        allocate(PtGauss)
                        PtGauss%coord(1) = coord(1,i)
                        PtGauss%coord(2) = coord(2,i)
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

    !>
    !! Appelé à chaque pas de temps
    !! \brief
    !!
    !! \param integer it
    !<


    subroutine evalueSortieCapteur(it)


        implicit none
        integer :: it
        type(tCapteur),pointer :: capteur



        sortie_capteur_deformation = .FALSE.
        sortie_capteur             = .FALSE.

        ! boucle sur les capteurs
        allocate (capteur)
        capteur=>listeCapteur

        do while (associated(capteur))

            if (mod(it,capteur%frequence)==0) then ! on fait la sortie
                sortie_capteur = .TRUE.

                if (capteur%grandeur.eq."DEFORMATION") sortie_capteur_deformation = .TRUE.

            endif

            capteur=>capteur%suivant

        enddo






    end subroutine evalueSortieCapteur




    !>
    !! \brief
    !!
    !! \param type (domain) TDomain
    !<


    subroutine save_capteur(Tdomain, ntime)


        implicit none

        type (domain) :: TDomain
        integer ntime

        type(tCapteur),pointer :: capteur


        ! boucle sur les capteurs
        allocate (capteur)
        capteur=>listeCapteur

        do while (associated(capteur))

            if(capteur%type_calcul==0) then
                call sortirGrandeurSousCapteur(capteur,Tdomain)
            else
                call sortieGrandeurCapteur_interp(capteur,Tdomain, ntime)
            endif

            capteur=>capteur%suivant

        enddo


    end subroutine save_capteur






    !>
    !! \brief
    !! \return
    !!
    !! \param type(tCapteur) capteur
    !! \param type (domain) TDomain
    !<


    subroutine sortirGrandeurSousCapteur(capteur,Tdomain)

        implicit none

        type(tCapteur) :: capteur
        type (domain) :: TDomain

        integer :: fileId,i,ierr,nPtGaussTotal,tag,borneInf,iproc

        character(len=MAX_FILE_SIZE) :: fnamef
        type(tPtGauss),pointer :: PtGauss

        integer , dimension  (MPI_STATUS_SIZE) :: status

        real, dimension(2) :: val0,val
        real, dimension(:,:), allocatable :: grandeur

        real,dimension(:),allocatable :: recvbuf
        real,dimension(3) :: sendbuf


        ! tableau de correspondance
        allocate(recvbuf(3*Tdomain%Mpi_var%n_proc))
        sendbuf=0.
        recvbuf=0.

        ! ETAPE 0 : initialisations
        fileId=99 ! id du fichier des sorties capteur
        allocate (grandeur(2,max(1,capteur%nPtGauss))) ! allocation du tableau de la grandeur restreinte aux pts de Gauss captes
        grandeur(:,:)=0.                               ! si maillage vide donc pas de pdg, on fait comme si il y en avait 1


        ! ETAPE 1 : Recuperation des valeurs de la grandeur pour les pts de Gauss definis par le capteur
        allocate(PtGauss)
        PtGauss=>capteur%listePtGauss
        i=0
        do while (associated(PtGauss))

            i=i+1
            if (trim(capteur%grandeur).eq."DEFORMATION") then
                grandeur(1,i)=tdomain%GrandeurDeformation(PtGauss%num)
            endif

            if (trim(capteur%grandeur).eq."VITESSE") then
                grandeur(1,i)=tdomain%GrandeurVitesse(0,PtGauss%num)
                grandeur(2,i)=tdomain%GrandeurVitesse(1,PtGauss%num)
            endif

            if (trim(capteur%grandeur).eq."DEPLA") then
                grandeur(1,i)=tdomain%GrandeurDepla(0,PtGauss%num)
                grandeur(2,i)=tdomain%GrandeurDepla(1,PtGauss%num)
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
            if (Tdomain%MPI_var%my_rank.eq.0) then

                ! on remplit recvbuf avec les infos du proc 0
                recvbuf(1)=capteur%distanceMin
                recvbuf(2)=grandeur(1,1)
                recvbuf(3)=grandeur(2,1)

                ! et ensuite avec les infos des autres proc
                do iproc=1,Tdomain%Mpi_var%n_proc-1
                    borneInf=3*iproc+1
                    tag = 100*(iproc+1)+capteur%numero

                    call mpi_recv(recvbuf(borneInf:borneInf+2),3,MPI_DOUBLE_PRECISION,iproc,tag, Tdomain%communicateur,status,ierr)
                enddo


                ! les comm st bloquantes, ici, on a forcement le tableau recvbuf complet
                ! on recherche le pdg le + proche et ses valeurs associees dans recvbuf
                ! le resultat est stocke dans recvbuf(1:3)
                do iproc=1,Tdomain%Mpi_var%n_proc-1
                    borneInf=3*iproc+1
                    if (recvbuf(borneInf).lt.recvbuf(1)) then
                        recvbuf(1)=recvbuf(borneInf)
                        recvbuf(2)=recvbuf(borneInf+1)
                        recvbuf(3)=recvbuf(borneInf+2)
                    endif
                enddo

                val0(1)=recvbuf(2)
                val0(2)=recvbuf(3)

            else ! cas des autres procs : on envoie des infos au proc 0

                sendbuf(1)=capteur%distanceMin


                sendbuf(2)=grandeur(1,1) ! grandeur en x au pdg le plus proche
                sendbuf(3)=grandeur(2,1) ! grandeur en x au pdg le plus proche
                tag = 100*(Tdomain%MPI_var%my_rank+1)+capteur%numero

                call mpi_send(sendbuf,3,MPI_DOUBLE_PRECISION,0,tag, Tdomain%communicateur,ierr)
            endif

        else ! on est dans le cas d un capteur a rayon d action

            ! CALCUL DE MOYENNE
            if (capteur%operation.eq."MOY") then
                val(1:2) = sum(grandeur,dim=2) ! on somme la grandeur sur les pdg du proc courant

                ! on fait la somme total des grandeurs de l'ensemble des procs qu'on ramene sur proc 0
                call MPI_Reduce (val,val0,2, MPI_DOUBLE_PRECISION, MPI_SUM, 0, Tdomain%communicateur, ierr)
                ! on recupere sur proc 0 le nombre total de capteurs
                call MPI_Reduce (capteur%nPtGauss,nPtGaussTotal,1, MPI_INTEGER, MPI_SUM, 0, Tdomain%communicateur, ierr)

                ! on moyenne sur le proc 0
                if(Tdomain%Mpi_var%my_rank==0) then
                    val0(:)=val(:)/nPtGaussTotal
                endif

            endif

            ! CALCUL DU MAX
            if (capteur%operation.eq."MAX") then
                val(:) = maxval(grandeur,dim=2)
                call MPI_Reduce (val, val0,2, MPI_DOUBLE_PRECISION, MPI_MAX,0, Tdomain%communicateur, ierr)
            endif


            ! CALCUL DU MIN
            if (capteur%operation.eq."MIN") then
                val(:) = minval(grandeur,dim=2)
                call MPI_Reduce (val, val0,2, MPI_DOUBLE_PRECISION, MPI_MIN,0, Tdomain%communicateur, ierr)
            endif

        endif



        ! ETAPE 4 : Impression du resultat dans le fichier de sortie par le proc 0
        if(Tdomain%Mpi_var%my_rank==0)then

            if (trim(capteur%grandeur).eq."DEFORMATION") then

                call semname_capteur_type(capteur%nom,"_deformation",fnamef)
                open(fileId,file=trim(fnamef),status="unknown",form="formatted",position="append")
                write(fileId,'(2(1X,E16.8))') Tdomain%TimeD%rtime,val0(1)
            endif


            if (trim(capteur%grandeur).eq."VITESSE") then

                call semname_capteur_type(capteur%nom,"_vitesse",fnamef)
                open(fileId,file=trim(fnamef),status="unknown",form="formatted",position="append")
                write(fileId,'(3(1X,E16.8))') Tdomain%TimeD%rtime,val0(1),val0(2)
            endif

            call flush(fileId) !!Gsa
            close(fileId)

            if (trim(capteur%grandeur).eq."DEPLA") then

                call semname_capteur_type(capteur%nom,"_depla",fnamef)
                open(fileId,file=trim(fnamef),status="unknown",form="formatted",position="append")
                write(fileId,'(3(1X,E16.8))') Tdomain%TimeD%rtime,val0(1),val0(2)
            endif

            close(fileId)


        endif

        deallocate(grandeur)
        deallocate(recvbuf)


    end subroutine sortirGrandeurSousCapteur


    !! effectue l'interpolation des grandeurs dans la maille dans laquelle se trouve le capteur
    !! la maille se trouve dans un seul proc
    !! seul le proc gere l'ecriture
    !!
    subroutine sortieGrandeurCapteur_interp(capteur,Tdomain, ntime)

        implicit none

        type(tCapteur) :: capteur
        type (domain) :: TDomain
        integer ntime

        integer :: fileId,i,j, ierr, tag
        character(len=MAX_FILE_SIZE) :: fnamef

        integer , dimension  (MPI_STATUS_SIZE) :: status
        integer request

        real, dimension(2) :: val0
        real, dimension(:), allocatable :: grandeur

        real,dimension(2) :: recvbuf
        real,dimension(2) :: sendbuf
        integer n_el, ipoint, ngllx, ngllz, mat
        real xi, eta
        real outx, outz

        integer numproc

        ! tableau de correspondance
        sendbuf = 0.
        recvbuf = 0.

        ! ETAPE 0 : initialisations
        fileId=99 ! id du fichier des sorties capteur
        allocate (grandeur(2)) ! allocation du tableau de la grandeur restreinte aux pts de Gauss captes
        grandeur(:)=0.                               ! si maillage vide donc pas de pdg, on fait comme si il y en avait 1


        ! Recuperation du numero de la maille Sem et des abscisses
        n_el = capteur%n_el
        xi = capteur%xi
        eta = capteur%eta

        ! ETAPE 1 : interpolations

        numproc = -1
        if((n_el/=-1) .AND. (capteur%numproc==Tdomain%MPI_var%my_rank)) then
            ngllx = Tdomain%specel(n_el)%ngllx
            ngllz = Tdomain%specel(n_el)%ngllz
            mat = Tdomain%specel(n_el)%mat_index

            do i = 0,ngllx - 1
                do j = 0,ngllz - 1
                    call  pol_lagrange (ngllx,Tdomain%sSubdomain(mat)%GLLcx,i,xi,outx)
                    call  pol_lagrange (ngllz,Tdomain%sSubdomain(mat)%GLLcz,j,eta,outz)

                    ipoint = Tdomain%specel(n_el)%Iglobnum(i,j)

                    if (trim(capteur%grandeur).eq."DEFORMATION") then
                        grandeur(1) = grandeur(1) + outx*outz*Tdomain%GrandeurDeformation(ipoint)
                    elseif (trim(capteur%grandeur).eq."VITESSE") then
                        grandeur(1) = grandeur(1) + outx*outz*Tdomain%GrandeurVitesse(0,ipoint)
                        grandeur(2) = grandeur(2) + outx*outz*Tdomain%GrandeurVitesse(1,ipoint)
                    elseif (trim(capteur%grandeur).eq."DEPLA") then
                        grandeur(1) = grandeur(1) + outx*outz*Tdomain%GrandeurDepla(0,ipoint)
                        grandeur(2) = grandeur(2) + outx*outz*Tdomain%GrandeurDepla(1,ipoint)
                    endif

                enddo
            enddo

            sendbuf(1) = grandeur(1) ! grandeur en x
            sendbuf(2) = grandeur(2) ! grandeur en y
            tag = 100*(Tdomain%MPI_var%my_rank+1)+capteur%numero

            call mpi_Isend(sendbuf,2,MPI_DOUBLE_PRECISION,0,tag, Tdomain%communicateur,request, ierr)
            numproc = Tdomain%MPI_var%my_rank  !!numero du proc courant

            call mpi_wait(request,status,ierr)
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
        if (Tdomain%MPI_var%my_rank.eq.0) then
            if ((capteur%numproc)>-1) then

                ! on remplit recvbuf avec les infos du proc 0
                recvbuf(:)=grandeur(:)

                ! et ensuite avec les infos des autres proc
                tag = 100*(capteur%numproc+1) + capteur%numero

                call mpi_recv(recvbuf(1),2,MPI_DOUBLE_PRECISION,capteur%numproc,tag, Tdomain%communicateur,status,ierr)
                val0(1)=recvbuf(1)
                val0(2)=recvbuf(2)
            else
                if(ntime<=1) write(6,'(A,A,A,1X,I3)') 'Le capteur ',capteur%nom,' n''est pas pris en compte car sur aucun proc',capteur%numproc
            endif
        endif

        ! ETAPE 4 : Impression du resultat dans le fichier de sortie par le proc 0
        if((Tdomain%Mpi_var%my_rank==0 .AND. (capteur%numproc>-1)))then

            if (trim(capteur%grandeur).eq."DEFORMATION") then

                call semname_capteur_type(capteur%nom,"_deformation",fnamef)
                open(fileId,file=trim(fnamef),status="unknown",form="formatted",position="append")
                write(fileId,'(2(1X,E16.8))') Tdomain%TimeD%rtime,val0(1)
            endif


            if (trim(capteur%grandeur).eq."VITESSE") then

                call semname_capteur_type(capteur%nom,"_vitesse",fnamef)
                open(fileId,file=trim(fnamef),status="unknown",form="formatted",position="append")
                write(fileId,'(3(1X,E16.8))') Tdomain%TimeD%rtime,val0(1),val0(2)
            endif

            close(fileId)

            if (trim(capteur%grandeur).eq."DEPLA") then

                call semname_capteur_type(capteur%nom,"_depla",fnamef)
                open(fileId,file=trim(fnamef),status="unknown",form="formatted",position="append")
                write(fileId,'(3(1X,E16.8))') Tdomain%TimeD%rtime,val0(1),val0(2)
            endif

            close(fileId)


        endif

        deallocate(grandeur)


    end subroutine sortieGrandeurCapteur_interp

    !!on identifie la maille dans laquelle se trouve le capteur. Il peut y en avoir plusieurs, alors le capteur est sur un face ou sur un sommet
    !!
    subroutine trouve_capteur(Tdomain, capteur, n_el, xi, eta)

        implicit none
        type (domain), INTENT(INOUT)  :: Tdomain
        type(Tcapteur) :: capteur
        integer ipoint
        integer i, n, n_el
        real dist
        integer P(0:3)
        real coor(0:3,0:1)
        real eps
        real xi, eta
        logical dans_maille
        integer i_sens
        eps=1.e-4

        n_el = -1
        i_sens = 0 !sens de parcours continu
        do n=0,Tdomain%n_elem-1
            ipoint = Tdomain%specel(n)%Control_Nodes(0) !!premier noeud de la maille n
            dist = (capteur%coord(1) - Tdomain%Coord_Nodes(0,ipoint))**2 + (capteur%coord(2) - Tdomain%Coord_Nodes(1,ipoint))**2
            dist = sqrt(dist)
            if(dist> Tdomain%specel(n)%dist_max) then
                !        print*,'dist',dist,n,Tdomain%specel(n)%dist_max
                cycle
            else
                !maille retenue
                P(0:3) = Tdomain%specel(n)%Control_Nodes(0:3)

                do i=0,3
                    coor(i,0:1) = Tdomain%Coord_Nodes(0:1,P(i))
                enddo

                dans_maille = test_contour_capteur(coor, capteur, i_sens)

                if( dans_maille) then
                    !! on peut retenir cette maille et calculer xi, eta
                    !! si le capteur est sur une face ou un sommet, il suffit de retenir une maille et de faire l'interpolation.
                    !! Quel que soit l'element on obtiendra la bonne valeur
                    write(6,*) 'on retient la maille',n,capteur%coord(1:2), Tdomain%Coord_Nodes(0:1,Tdomain%specel(n)%Control_Nodes(0:3))
                    call calc_xi_eta_capteur(capteur, coor, xi, eta)
                    n_el = n
                    capteur%numproc = Tdomain%MPI_var%my_rank
                    return
                else
                    cycle

                endif
            endif
        enddo

        if(n_el/=-1) then
            write(6,*) 'xi,eta trouves pour le capteur ',capteur%nom,' de positions ',capteur%coord(1), capteur%coord(2),' sont' ,  &
                xi, eta, 'dans l element ',n_el,' pour le proc ',Tdomain%MPI_var%my_rank
        endif

    end subroutine trouve_capteur


    !!
    !! calcul des abscisses curvilignes du cpateur
    !! processus de dichotomie
    !!
    subroutine calc_xi_eta_capteur(capteur, coord, xitrouve_def, etatrouve_def)

        implicit none
        type(Tcapteur) :: capteur
        real xi_min, xi_max, xi0, eta0, eta_min, eta_max, xitrouve, etatrouve
        real P(4,0:1), dxi, deta, coord(0:3,0:1), xitrouve_def, etatrouve_def
        real dist_def, dist, dist_k
        real xi, eta
        real eps
        integer i, j, k, im, il, idim
        integer n_it
        logical interieur
        integer, parameter :: n_itmax=2000
        integer i_sens
        !! attention si le point du capteur se trouve partage entre
        !! plusieurs elements (sommets, face ou aretes)
        !! identifiation

        xi_min = -1.  !bornes min et max de la zone d'etude dans le carre de reference
        xi_max = 1.   !on coupe en 2 dans chaque direction la zone d'etude
        eta_min = -1.
        eta_max = 1.
        dist_def = huge(1.) !distance entre la solution et le point de depart
        eps = 1e-6 !tolerance pour accepter la solution
        n_it = 0  !nombre d'iterations
        xi0 = xi_min
        eta0 = eta_min
        i_sens = 1 !sens de parcours non continu

        do while((dist_def > eps) .AND. (n_it<n_itmax))
            n_it = n_it +1
            !on subdivise la zone en 4 sous-zones d'etude
            do im=0,1
                do il=0,1
                    dxi = (xi_max - xi_min)/2.
                    deta = (eta_max - eta_min)/2.
                    xi0 = xi_min + il*dxi
                    eta0 = eta_min + im*deta
                    k = 0
                    dist = huge(1.)
                    !Dans les boucles sur i et j, on definit les points P
                    do j=0,1
                        do i=0,1
                            xi = xi0 + i*dxi
                            eta = eta0 + j*deta
                            k = k + 1
                            do idim=0,1
                                P(k,idim)=  0.25 * (coord(0,idim)*(1-xi)*(1-eta) + coord(1,idim)*(1+xi)*(1-eta) + &
                                    coord(2,idim)*(1+xi)*(1+eta) + coord(3,idim)*(1-xi)*(1+eta))
                            enddo
                            dist_k = sqrt((P(k,0) - capteur%coord(1))**2 + (P(k,1) - capteur%coord(2))**2)
                            dist = min(dist,dist_k)
                            !
                            !on teste si P(k) correspond au capteur. Si oui on sort
                            !
                            if(dist_k < eps) then
                                xitrouve_def = xi
                                etatrouve_def = eta
                                print *,'xtrouve ytrouve',P(k,0),P(k,1),xi,eta,n_it
                                return
                            endif

                            !! on teste si P(k) realise le min de la distance au capteur

                            if(abs(dist - dist_k) < eps) then
                                xitrouve = xi
                                etatrouve = eta
                            endif

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
                        xitrouve_def = xitrouve
                        etatrouve_def = etatrouve
                        dist_def = dist
                    endif
                enddo
            enddo
        enddo
        !  print *,'xtrouve ytrouve',xtrouve_def,ytrouve_def,xitrouve_def,etatrouve_def,n_it

    end subroutine calc_xi_eta_capteur


    logical function test_contour_capteur(P, capteur, i_sens)

!!!!  Description: on teste si le point du capteur est a l'interieur du contour defini par les points P.
!!!!  Si le capteur est sur le bord ou sur un sommet, renvoie true
!!!!  * attention : le sens de parcours dans la face (habituellement 1-2-3-4) est ici 1-2-4-3 (d'ou la variable Q)
!!!!
!!!!  Historique: 03/10 - Gsa Ipsis - Creation
!!!! --------------------------------------------------

        implicit none
        type(Tcapteur) :: capteur
        real P(4,0:1), Q(4,0:1)
        real val1, val2
        real eps
        integer i, j
        integer i_sens
        eps=1.e-7

        test_contour_capteur = .true.
        !! A cause du sens de parcours des noeuds 1-2-4-3 on cree le tableau Q
        !!
        Q = P
        if(i_sens==1) then
            Q(3,0:1) = P(4,0:1)
            Q(4,0:1) = P(3,0:1)
        endif

        val2 = (Q(1,0) - capteur%coord(1))*(Q(2,1) - capteur%coord(2)) - &
            (Q(2,0) - capteur%coord(1))*(Q(1,1) - capteur%coord(2))

        do i=2,4
            j = mod(i,4) + 1
            val1 = ( (Q(i,0) - capteur%coord(1))*(Q(j,1) - capteur%coord(2)) - &
                (Q(j,0) - capteur%coord(1))*(Q(i,1) - capteur%coord(2)) )
            if(val1*val2<0.) then
                ! if(val1*val2<-eps) then !!test trop restrictif
                test_contour_capteur = .false.
                exit
            endif
            ! si val1 <> 0, on le copie dans val2, sinon on garde l'ancien val2 qui doit etre <> 0
            if(abs(val1)>0.) then
                val2 = val1
            else
                if(abs(val2)<1.e-7) then
                    print*,'Erreur test_contour_capteur',capteur%nom,capteur%coord
                endif
            endif
        enddo
    end function test_contour_capteur

end module Mcapteur
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
