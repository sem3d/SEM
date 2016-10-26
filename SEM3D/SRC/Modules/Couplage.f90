!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Couplage.f90
!!\brief Contient les routines propres au traitement du couplage.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ref_orient

    !tableau de correspondance num_local de face/numero des noeuds de l'element
    !
    !(face 0 - k=0 | face 1 - j=0 | face 2 - i=n-1 | face 3 - j=n-1 | face 4 - i=0 | face 5 - k=n-1 )
    !
    integer,dimension(0:5,0:3)  :: NOEUD_FACE
    parameter (NOEUD_FACE = reshape &
        & ((/ 0, 0, 1, 3, 0, 4, 1, 1, 2, 2, 3, 5, 2, 5, 6, 6, 7, 6, 3, 4, 5, 7, 4, 7 /),(/6,4/)))

    !tableau de correspondance num_local de face/numero des aretes de l'element
    !
    !(face 0 - k=0 | face 1 - j=0 | face 2 - i=n-1 | face 3 - j=n-1 | face 4 - i=0 | face 5 - k=n-1 )

    integer,dimension(0:5,0:3) :: ARETE_FACE
    parameter (ARETE_FACE = reshape &
        & ((/ 0, 0, 1, 2, 3, 5, 1, 4, 7, 7, 10, 8, 2, 5, 8, 9, 11, 9, 3, 6, 4, 10, 6, 11/),(/6,4/)))

end module ref_orient

module scouplage

#ifdef COUPLAGE

    use sdomain
    use semdatafiles
    use mpi

    implicit none

    type :: arete
       real, dimension(0:2) :: coord
       integer, dimension(2) :: numface
       integer num, num_def
       integer ngll, num_arete_local !num_arete_local entre 0 et 3
       integer, dimension(0:9) :: num_imat !numeros globaux des ddl de l'arete dans la matrice de couplage !attention taille en dur (10 pts de Gauss au max ici)!!
    end type arete

    type :: sommet_gll
       real, dimension(0:2) :: coord
       integer nbface
       integer, dimension(12) :: numface !au plus un noeud est connecte a 12 faces !attention taille en dur!!
       integer num, num_def, numglobal, num_local !num_sommet_local entre 0 et 3
       integer :: num_imat !numero global du sommet dans la matrice de couplage
    end type sommet_gll

    type :: type_comm_couplage

       integer :: m_dim                                   ! dimension de l espace
       integer :: m_nb_point_de_couplage                  ! nombre de particules de couplage sur le proc courant
       integer :: m_NbFace                                ! nombre de face de couplage sur le proc courant
       integer :: m_MaxNgParFace                          ! nombre max de pts de Gauss sur une face de couplage du proc courant
       integer :: m_gl_Ngauss                             ! nombre total max  de points de gauss sur tte la surface de couplage
       integer :: m_local_ngauss                          ! nombre de pt de gauss de couplage sur ce processeur
       integer,dimension(:),pointer  :: m_numFace         ! numeros des faces de couplage
       integer,dimension(:),pointer  :: m_numMaille       ! numeros des mailles sem contenant les faces de couplage
       integer,dimension(:),pointer  :: m_tabFaceCouplage ! liste des faces de couplage
       real,dimension(:),pointer     :: m_pos_proj        ! abscisse curviligne sur la face de couplage du projete de la particule de couplage
       !       integer,dimension(:,:),pointer :: m_Pk_plus_proche ! numero local du Pk le plus proche
       real,dimension(:), pointer    :: m_surf_part       ! surface des particules Mka de couplage

    end type type_comm_couplage

    type :: type_face_couplage

       integer :: face          ! numero de la face de couplage sem dans la numerotation du proc sem qui la contient
       integer :: proc          ! numero du processeur sem qui contient cette face de couplage
       integer :: ngll, ngll1, ngll2      ! nbre de point de gauss sur la face
       integer,dimension(0:3) :: noeud    ! numero des noeuds de la face
       real,dimension(0:3,0:2) :: coord   ! coordonnees des noeuds
       integer :: elem          ! numero de la maille
       integer :: numlocal      ! numero local de la face
       integer :: nbPk          ! nombre de points du peigne

    end type type_face_couplage

    real,dimension(:), pointer :: Xpdc, Ypdc, Zpdc

    integer :: nbFaceTotal,nbFace
    type(type_face_couplage), dimension(:),pointer :: face_couplage
    type(type_comm_couplage) :: comm_couplage

contains

!!!! --------------------------------------------------
!!!! --------------------------------------------------
!!!! --------------------------------------------------
    !>
    !! \brief Alloue les tableaux utilisés lors du couplage.
    !!
    !! \param integer,intent(IN) dim
    !<


    subroutine alloue_couplage(dim)

        integer,intent(IN) :: dim

        comm_couplage%m_dim=dim

        allocate(comm_couplage%m_numFace(comm_couplage%m_nb_point_de_couplage))
        allocate(comm_couplage%m_numMaille(comm_couplage%m_nb_point_de_couplage))
        allocate(comm_couplage%m_pos_proj(comm_couplage%m_dim*comm_couplage%m_nb_point_de_couplage))

        !il y a au maximum autant de face de couplage que de points de couplage
        allocate(comm_couplage%m_tabFaceCouplage(comm_couplage%m_nb_point_de_couplage))

        allocate(comm_couplage%m_surf_part(comm_couplage%m_nb_point_de_couplage))
    end subroutine alloue_couplage

!!!! --------------------------------------------------
!!!! --------------------------------------------------
!!!! --------------------------------------------------
    !>
    !! \brief Desalloue les tableaux utilisés lors du couplage
    !!
    !<
    subroutine desalloue_couplage()

        deallocate(comm_couplage%m_numFace)
        deallocate(comm_couplage%m_numMaille)
        deallocate(comm_couplage%m_pos_proj)
        deallocate(comm_couplage%m_surf_part)
        deallocate(comm_couplage%m_tabFaceCouplage)

    end subroutine desalloue_couplage

!!!! --------------------------------------------------
!!!! --------------------------------------------------
!!!! --------------------------------------------------
    !>
    !! \brief Desalloue les tableaux utilisés lors du couplage
    !!
    !<

    !>
    !! \brief Iniatilise les données propres au couplage.
    !! Effectue notamment des communications avec le superviseur
    !! afin d'envoyer le pas de temps et de recevoir des informations sur le
    !! nombre de particules de couplage, les mailles couplées,...
    !! \param type (domain), intent(INOUT) Tdomain
    !<


    subroutine initialisation_couplage(Tdomain, MaxNgParDir)

!!!!  Description: Initialisation du couplage pour Sem3d
!!!!
!!!!  Historique: 10/08 - Gsa Ipsis - Adaptation du 2D au 3D
!!!!
!!!!  * Envoie les variables necessaires au superviseur
!!!!  * Calcul des matrices elementaires permettant la construction de la matrice de couplage
!!!!
!!!!  Descriptif des variables:
!!!!
!!!! --------------------------------------------------
        type (domain), intent(INOUT)  :: Tdomain
        integer, intent(out) :: MaxNgParDir
        !
        integer :: rg, nb_procs
        integer :: i,j,l
        integer :: tag,ierr
        integer :: np, ngll  !, ngllx, ngllz,
        integer :: numFace
        integer :: iface, bufsize,decal,nbchamps
        integer :: Ngauss, MaxNgParFace, gl_MaxNgParFace, gl_MaxNgParDir !, iout
        logical :: dejaPresent

        integer, dimension (MPI_STATUS_SIZE) :: status

        integer, dimension(:),pointer:: tabNbFace
        integer,dimension(Tdomain%nb_procs) :: displs1, count1
        integer,dimension(5) :: ibuf
        real, dimension (:),pointer :: buf

        integer nb_point_de_couplage
        integer gl_Ngauss

        real,dimension (:), pointer :: dmin_couplage


        rg = Tdomain%rank
        nb_procs = Tdomain%nb_procs
        ! envoi du pas de temps sem
        if (rg == 0) then
            tag=500000

            call MPI_SEND(Tdomain%TimeD%dtmin,1,MPI_DOUBLE_PRECISION, Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )
        endif

        ! reception du nombre de faces et de mailles de couplage
        tag=1000000+rg
        call MPI_RECV(comm_couplage%m_nb_point_de_couplage,1, MPI_INTEGER, Tdomain%master_superviseur,tag, &
            Tdomain%communicateur_global, status,ierr)

        call alloue_couplage(3)

        if (comm_couplage%m_nb_point_de_couplage.gt.0) then


            ! reception du tableau des mailles de couplage
            tag=1100000+rg
            call MPI_RECV(comm_couplage%m_numMaille,comm_couplage%m_nb_point_de_couplage, MPI_INTEGER, Tdomain%master_superviseur,tag, &
                Tdomain%communicateur_global, status,ierr)

            ! reception du tableau des faces de couplage
            tag=1200000+rg
            call MPI_RECV(comm_couplage%m_numFace,comm_couplage%m_nb_point_de_couplage, MPI_INTEGER, Tdomain%master_superviseur,tag, &
                Tdomain%communicateur_global, status,ierr)

            ! reception du tableau des projetes de couplage
            tag=1300000+rg
            call MPI_RECV(comm_couplage%m_pos_proj,comm_couplage%m_dim*comm_couplage%m_nb_point_de_couplage,MPI_DOUBLE_PRECISION,Tdomain%master_superviseur,&
                tag, Tdomain%communicateur_global, status,ierr)

        endif

        Nbface=0



        ! calcul des fonctions de forme en chaque point de couplage
        nb_point_de_couplage =  comm_couplage%m_nb_point_de_couplage

        allocate (Xpdc(nb_point_de_couplage))
        allocate (Ypdc(nb_point_de_couplage))
        allocate (Zpdc(nb_point_de_couplage))

        do np=1,comm_couplage%m_nb_point_de_couplage
            ! recuperation de la face et de l element du point de couplage courant et de ses coordonnees, du numero de materiau et du nombre de pdg en x et z
            numFace = comm_couplage%m_numFace(np)

            XpdC(np) = comm_couplage%m_pos_proj(3*(np-1)+1)
            YpdC(np) = comm_couplage%m_pos_proj(3*(np-1)+2)
            ZpdC(np) = comm_couplage%m_pos_proj(3*np)
            
            dejaPresent=.false.
            do i=1,Nbface
                if (comm_couplage%m_tabFaceCouplage(i).eq.numFace) then
                    dejaPresent=.true.
                    exit
                endif
            enddo

            ! liste des faces de couplage du proc courant
            if (.not.dejaPresent) then
                Nbface=Nbface+1
                comm_couplage%m_tabFaceCouplage(Nbface)=numFace
            endif

        enddo ! fin de boucle sur les particules de couplage du proc courant

        ! on conserve le nombre de face de couplage du proc courant
        comm_couplage%m_Nbface=Nbface

        allocate(tabNbFace(nb_procs))
        tabNbFace(:)=0

        if (rg == 0) then
            count1(:)=1
            displs1(1)=0
            do i=2,nb_procs
                displs1(i)=displs1(i-1)+count1(i-1)
            enddo
        endif

        ! tabNbFace,defini uniquement sur proc 0,  contient le nbre de face de couplage de chaque proc
        call MPI_Gatherv (Nbface,1,MPI_INTEGER,tabNbFace,count1,displs1,MPI_INTEGER,0, Tdomain%communicateur, ierr)


        ! calcul du nombre de face de couplage total
        if (rg == 0) then
            NbFaceTotal=sum(tabNbFace)
            allocate(face_couplage(NbFaceTotal))
        else
            allocate(face_couplage(Nbface))
        endif

        call MPI_Bcast(NbFaceTotal,1,MPI_INTEGER,0,Tdomain%communicateur,ierr)
        write(50,*)'Le nombre de faces de couplage total est ' ,Nbfacetotal !GSa Ipsis


        ! Determination du nbre max de points de gauss par maille et du nombre total de points de gauss de la zone de couplage sur chaque proc
        MaxNgParFace=0 ! nombre max de points de gauss par bord de maille sur le proc courant
        MaxNgParDir=0 ! nombre max de points de gauss par arete de maille sur le proc courant !Gsa Ipsis
        Ngauss = 0 ! nombre max de points de gauss sur la zone de couplage du proc courant


        ! remplissage de la structure face_couplage

        call remplit_face_couplage(Tdomain, MaxNgParFace, MaxNgParDir)

        !Il faut transmettre gl_MaxNgPardir a tous les procs car mpi_scatterv l'utilise
        call MPI_AllReduce (MaxNgParDir, gl_MaxNgParDir, 1, MPI_INTEGER, MPI_MAX, Tdomain%communicateur, ierr)

        ! on conserve le nombre max de gauss par face de couplage du proc courant
        comm_couplage%m_MaxNgParFace=MaxNgParFace

        ! calcul (sur proc 0) du nombre max total de ddl sur une face de couplage
        call MPI_Reduce (MaxNgParFace,gl_MaxNgParFace,1,MPI_INTEGER,MPI_MAX,0, Tdomain%communicateur, ierr)

        allocate (dmin_couplage(Nbface))
        call dist_min_point_de_couplage(comm_couplage%m_nb_point_de_couplage, Xpdc, Ypdc, Zpdc, dmin_couplage)

        call definit_nb_pts_Pk(Tdomain, dmin_couplage)
        deallocate(dmin_couplage)

        ! on collecte sur le proc 0 un ensemble des donnees pour ttes les faces de couplage sur les differents proc
        NbChamps = 20  ! taille des donnees pour une face_couplage
        if (rg == 0) then
            ! reception par le proc maitre des faces de couplage
            decal=0
            do i=1,nb_procs-1
                decal = decal+tabNbFace(i)
                if (tabNbFace(i+1)>0) then
                    bufsize = NbChamps*tabNbFace(i+1)
                    allocate(buf(bufsize))
                    tag = 70000+i
                    call MPI_RECV(buf,bufsize, MPI_DOUBLE_PRECISION, i,tag, Tdomain%communicateur, status,ierr)

                    do iface=1,tabNbFace(i+1)
                        face_couplage(iface+decal)%face   = int(buf(1+NbChamps*(iface-1)))
                        face_couplage(iface+decal)%proc   = int(buf(2+NbChamps*(iface-1)))
                        face_couplage(iface+decal)%ngll1  = int(buf(3+NbChamps*(iface-1)))
                        face_couplage(iface+decal)%ngll2  = int(buf(4+NbChamps*(iface-1)))
                        do j=0,3
                            face_couplage(iface+decal)%noeud(j) = int(buf(5+j+NbChamps*(iface-1)))
                        enddo
                        do l=0,2
                            do j=0,3
                                face_couplage(iface+decal)%coord(j,l) = buf(9+j+4*l+NbChamps*(iface-1))
                            enddo
                        enddo
                    enddo
                    deallocate(buf)
                endif

            enddo

        else ! envoi des eventuelles faces de couplage du proc courant au proc maitre

            if (Nbface>0) then
                bufsize = NbChamps*Nbface
                allocate(buf(bufsize))
                do iface=1,Nbface
                    buf(1+NbChamps*(iface-1)) = face_couplage(iface)%face
                    buf(2+NbChamps*(iface-1)) = face_couplage(iface)%proc
                    buf(3+NbChamps*(iface-1)) = face_couplage(iface)%ngll1
                    buf(4+NbChamps*(iface-1)) = face_couplage(iface)%ngll2
                    do j=0,3
                        buf(5+j+NbChamps*(iface-1)) = face_couplage(iface)%noeud(j)
                    enddo
                    do l=0,2
                        do j=0,3
                            buf(9+j+4*l+NbChamps*(iface-1)) = face_couplage(iface)%coord(j,l)
                        enddo
                    enddo
                enddo

                tag = 70000+rg
                call MPI_SEND(buf,bufsize,MPI_DOUBLE_PRECISION, 0,tag, Tdomain%communicateur, ierr )

                deallocate(buf)
            endif

        endif ! finsi proc non maitre


        ! tri des faces de couplage et construction des tableaux de reperage position* pour l'assemblage de la matrice
        call calc_nb_gll(Ngauss)
        comm_couplage%m_local_ngauss = Ngauss

        !Le nb total de pts de gauss a ete calcule dans ordonne_face_couplage_3D(different d'avant)
        ! calcul (sur proc 0) du nombre max total de points de gauss de la surface de couplage
        call MPI_AllReduce (Ngauss,comm_couplage%m_gl_Ngauss,1,MPI_INTEGER,MPI_SUM, Tdomain%communicateur, ierr)


        gl_Ngauss = comm_couplage%m_gl_Ngauss

        ! envoi au superviseur
        if (rg == 0) then
            tag=600000
            ibuf(1) = 3  ! dim du pb
            ibuf(2) = comm_couplage%m_gl_Ngauss
            ibuf(3) = gl_MaxNgParFace
            ibuf(4) = gl_MaxNgParDir
            ibuf(5) = 0
            call MPI_SEND(ibuf,5,MPI_INTEGER, Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )
        endif

        ibuf(1) = Nbface
        ibuf(2) = MaxNgParFace
        ibuf(3) = MaxNgParDir
        ibuf(4) = Ngauss
        ! envoi au superviseur du nb de face de couplage  et de MaxNgParFace, pour chaque proc
        tag=610000+rg
        call MPI_SEND(ibuf,4, MPI_INTEGER, Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )

        ! construction de la matrice elementaire de chaque face de couplage et envoi au superviseur, pour chaque proc sem
        do iface=1,Nbface

            numFace=comm_couplage%m_tabFaceCouplage(iface)
            ngll = Tdomain%sFace(numFace)%ngll1 *Tdomain%sFace(numFace)%ngll2

            ibuf(1) = Tdomain%sFace(numFace)%ngll1
            ibuf(2) = Tdomain%sFace(numFace)%ngll2
            ibuf(3) = face_couplage(iface)%nbPk
            ibuf(4) = face_couplage(iface)%numlocal
            ibuf(5) = numFace

            tag=670000+rg+iface-1
            call MPI_SEND(ibuf, 5, MPI_INTEGER, Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )

            tag=690000+rg+iface-1
            call MPI_SEND(face_couplage(iface)%coord, 12, MPI_DOUBLE_PRECISION, Tdomain%master_superviseur, &
                tag, Tdomain%communicateur_global, ierr )

        enddo ! fin boucle sur les faces


        ! Pas de protection/reprise a l'initiative de sem. La protection est pilotee par le superviseur
        Tdomain%logicD%save_restart =.false.
        Tdomain%logicD%run_restart =.false.
        Tdomain%TimeD%iter_reprise =0
        Tdomain%TimeD%prot_m0 =0
        Tdomain%TimeD%prot_m1 =0
        Tdomain%TimeD%prot_m2 =0


        ! reception de l iteration de reprise
        tag=900000+rg
        call MPI_RECV(Tdomain%TimeD%iter_reprise,1, MPI_INTEGER, Tdomain%master_superviseur,tag, &
            Tdomain%communicateur_global, status,ierr)

        if (Tdomain%TimeD%iter_reprise > 0) then
            Tdomain%logicD%run_restart = .true.
        endif

        ! faire toutes les deallocations possibles
        deallocate(tabNbFace)
        deallocate(comm_couplage%m_pos_proj)

    end subroutine initialisation_couplage
!!!! --------------------------------------------------
!!!! --------------------------------------------------
!!!! --------------------------------------------------



    !>
    !! Numerotation des ddl de Sem situes sur la surface de couplage
    !!
    !! Calcul du nb total de ddl (=Ngauss) pour la matrice de couplage (comprend tous les pts de Gauss de la surface de couplage)
    !!
    !! La numerotation des ddl s'effectue suivant la numerotation des faces de couplage Sem.
!!!  Pour une face consideree, apres s'etre ramene au repere de reference, on numerote  les ddl de gauche a droite
    !!  et de bas en haut, sauf s'ils n'ont pas deja ete numerotes pour une face precedente
    !!  (concerne les ddl sommets ou ddl de l'interieur des aretes).
    !!  On numerote les sommets et les aretes de la surface de couplage.
    !<
    subroutine calc_nb_gll(Ngauss)

        integer i
        integer ngll1, ngll2
        integer, intent(OUT)  :: Ngauss

        !----------------------------------------------------------------
        !----------------------------------------------------------------
        !calcul du nombre de ddl par face suivant la numerotation des faces (position_face) - On part du plus petit numero de face
        Ngauss = 0
        do i=1,NbFace
            ngll1 = face_couplage(i)%ngll1
            ngll2 = face_couplage(i)%ngll2
            Ngauss = Ngauss + ngll1*ngll2
        enddo
        write(*,*) 'SEM3D: Nombre total de point de gauss:', Ngauss

    end subroutine calc_nb_gll


!!!! --------------------------------------------------
!!!! --------------------------------------------------
!!!! --------------------------------------------------
    !>
    !! \brief Evalue les forces sur les particules couplées.
    !! -Recoit les forces MKA3d appliquées aux particules de couplage,
    !! -calcul le second membre du système couplé à résoudre,
    !! -recoit la solution du système et applique les forces en conséquence.
    !!
    !! Calcul des forces aux points de Gauss apres resolution du systeme lineaire (faite par le Superviseur)
    !!
    !! \param type (domain), intent(INOUT) Tdomain
    !! \param integer, intent(IN) ntime
    !<


    subroutine calcul_couplage_force(Tdomain,ntime)

!!!!  Description: Calcul des forces aux points de Gauss apres resolution du systeme lineaire (faite par le Superviseur)
!!!!
!!!!  Historique: 10/08 - Gsa Ipsis - Reprise du 2D en 3D
!!!!              11/09 - Gsa Ipsis - Correction bug si nb_particules_couplage=0 sur le proc courant
!!!!
!!!!  * Reception depuis le superviseur des forces des particules de couplage Mka
!!!!  * Envoi au superviseur des donnees necessaires a la resolution du systeme lineaire (nombre de points de Gauss dans les deux directions de la face,
!!!!  * des numeros globaux dans la matrice de couplage de ces ddl,...)
!!!!  * Construction du second membre du systeme lineaire
!!!!  * Envoi au superviseur de ce second membre du systeme lineaire
!!!!  * Reception depuis le superviseur des "forces imposees" = forces calculees aux points de Gauss apres resolution du systeme lineaire
!!!!  * Calcul des surfaces liees aux points de Gauss en utilisant la differentielle de la fonction d'interpolation
!!!!  * (qui transforme l'element de reference Sem en element courant) - On n'utilise pas le Jacobien ( ca y ressemble)
!!!!  * Affectation des forces imposees aux variables ForcesExt des aretes, faces, sommets de la surface de couplage
!!!!
!!!!
!!!! --------------------------------------------------

        use ref_orient

        type (domain), intent(INOUT)  :: Tdomain
        integer, intent(IN) :: ntime

        integer :: i,j, numElem,numFace, mat
        integer :: ngll, ngllx, nglly, ngllz
        integer :: iface
        integer :: tag,ierr
        integer, dimension (MPI_STATUS_SIZE) :: status

        real,dimension(0:2) :: force
        integer numlocal, numEdge, numloc_edge
        integer,dimension(0:3) :: node
        integer :: ngll1, ngll2, nglltot
        !real :: outx, outz
        real :: tsurf
        integer :: iarete, ind1, ind2
        integer :: position
        real, dimension(comm_couplage%m_dim, comm_couplage%m_local_ngauss) :: force_impose

        !   flag pour sauvegarder la surface tsurf
        integer :: sauvetsurf

        sauvetsurf = 1
        if ( .not.Tdomain%logicD%run_restart) then
            if ( ntime .gt. 1 ) then
                sauvetsurf = 0
            endif
        else
            if ( ntime-1   .gt. Tdomain%TimeD%iter_reprise ) then
                sauvetsurf = 0
            endif
        endif

        ! RECEPTION DES FORCES DU SUPERVISEUR

        ! reception du vecteur force appliquees au pt de gauss de sem, solution du systeme resolu par le superviseur
        tag=1500000+10000*Tdomain%rank + ntime
        call MPI_RECV (force_impose,comm_couplage%m_dim*comm_couplage%m_local_ngauss, MPI_DOUBLE_PRECISION,&
            Tdomain%master_superviseur,tag, Tdomain%communicateur_global, status, ierr )


        ! on depouille le vecteur force imposee : en parcourant les
        ! faces de couplage et en recuperant, selon la position de chacune, sa
        ! force

        ! on recupere numface, numlocal (numero local de face par rapport a l'element de reference) et numElem (numero de l'element)
        position = 1
        do iface=1,comm_couplage%m_nbface
            numface  = face_couplage(iface)%face
            !    numElem  = Tdomain%sFace(numFace)%Near_Element(0)
            numElem = face_couplage(iface)%elem
            numlocal = face_couplage(iface)%numlocal

            ngll1 = face_couplage(iface)%ngll1
            ngll2 = face_couplage(iface)%ngll2
            nglltot = ngll1*ngll2
            ngllx = Tdomain%specel(numElem)%ngllx
            nglly = Tdomain%specel(numElem)%nglly
            ngllz = Tdomain%specel(numElem)%ngllz

            if( ((ngll1 /= ngllx) .AND. (ngll1 /= nglly)) .OR. &
                ((ngll2 /= nglly) .AND. (ngll1 /= ngllz))) then
                write(50,*) 'Probleme ngll pour la face ', numface, &
                    ' inversion des ddl '
            endif
            mat = Tdomain%specel(numElem)%mat_index

            do i=0,3
                node(i)= Tdomain%specel(numElem)%Near_Vertices(NOEUD_FACE(numlocal,i))
            enddo

            !! pour la surface rattachee a un point de Gauss, on a
            !! besoin de la differentielle de la fonction d'interpolation (=fonction
            !! qui transforme l'element de reference en element Sem courant )
            !! (fonction restreinte a la face ici)
            !! il semble qu'on ne puisse pas directement utiliser le Jacobien (il faudrait le ramener a la face)
            !!
            do j=1,ngll2-2
                do i=1,ngll1-2
                    ! calcul des poids au point de gauss
                    ! Surface correspondant a la surface de la cellule de voronoi associee
                    if ( sauvetsurf == 1 ) then
                        tsurf = surface_gll(Tdomain, mat, numElem, numlocal, i, j)
                        Tdomain%sFace(numface)%tsurfsem(i,j) = tsurf
                    else
                        tsurf = Tdomain%sFace(numface)%tsurfsem(i,j)
                    endif
                    !                    tsurf = surface_gll(Tdomain, mat, numElem, numlocal, i, j)
                    Tdomain%sFace(numface)%ForcesExt(i,j,:) = tsurf*force_impose(:, position+i + j*ngll1)
                    !                    Tdomain%sFace(numface)%FlagMka(i,j,:) = 1
                enddo
            enddo

            !aretes
            !prendre en compte les orientations (au contraire des faces)
            do iarete=0,3
                numloc_edge = ARETE_FACE(numlocal,iarete) !numero local de l'arete (entre 0 et 11)
                numEdge = Tdomain%specel(numElem)%Near_Edges(numloc_edge)
                ngll =  Tdomain%sEdge(numEdge)%ngll

                do j=1,ngll-2
                    ! attention aux indices de InterP_coeff
                    !
                    if(iarete==0) then
                        ind1 = j
                        ind2 = 0
                    elseif(iarete==1) then
                        ind1 = ngll1-1
                        ind2 = j
                    elseif(iarete==2) then
                        ind1 = j
                        ind2 = ngll2-1
                    elseif(iarete==3) then
                        ind1 = 0
                        ind2 = j
                    endif

                    if ( sauvetsurf == 1 ) then
                        tsurf = surface_gll(Tdomain, mat, numElem, numlocal, ind1, ind2)
                        Tdomain%sEdge(numEdge)%tsurfsem(j) = tsurf
                    else
                        tsurf = Tdomain%sEdge(numEdge)%tsurfsem(j)
                    endif
                    !                    tsurf = surface_gll(Tdomain, mat, numElem, numlocal, ind1, ind2)

                    force(:) = force_impose(:,position + ind1 + ind2*ngll1)

                    if(Tdomain%specel(numElem)%Orient_Edges(numloc_edge) == 0 ) then
                        Tdomain%sEdge(numEdge)%ForcesExt(j,:) = Tdomain%sEdge(numEdge)%ForcesExt(j,:) + tsurf*force
                        !                        Tdomain%sEdge(numEdge)%FlagMka(j,:) = 1
                    elseif (Tdomain%specel(numElem)%Orient_Edges(numloc_edge) == 1 ) then
                        select case (numloc_edge)
                        case (0,2,5,9)
                            Tdomain%sEdge(numEdge)%ForcesExt(ngllx-1-j,:) = Tdomain%sEdge(numEdge)%ForcesExt(ngllx-1-j,:) + tsurf*force
                            !                            Tdomain%sEdge(numEdge)%FlagMka(ngllx-1-j,:) = 1
                        case (1,3,8,11)
                            Tdomain%sEdge(numEdge)%ForcesExt(nglly-1-j,:) = Tdomain%sEdge(numEdge)%ForcesExt(nglly-1-j,:) + tsurf*force
                            !                            Tdomain%sEdge(numEdge)%FlagMka(nglly-1-j,:) = 1
                        case (4,6,7,10)
                            Tdomain%sEdge(numEdge)%ForcesExt(ngllz-1-j,:) = Tdomain%sEdge(numEdge)%ForcesExt(ngllz-1-j,:) + tsurf*force
                            !                            Tdomain%sEdge(numEdge)%FlagMka(ngllz-1-j,:) = 1
                        end select
                    else
                        write(50,*) 'Pb couplage_force - mauvaise valeur de orient_edges'
                        stop
                    endif
                enddo
            enddo

            !sommets
            do i=0,3
                node(i)= Tdomain%specel(numElem)%Near_Vertices(NOEUD_FACE(numlocal,i))

                select case (numlocal)
                case(0,5)
                    select case(i)
                    case(0)
                        ind1 = 0; ind2 = 0
                    case(1)
                        ind1 = ngllx - 1; ind2 = 0
                    case(2)
                        ind1 = ngllx - 1; ind2 = nglly - 1
                    case(3)
                        ind1 = 0; ind2 = nglly - 1
                    end select
                case  (1,3)
                    select case(i)
                    case(0)
                        ind1 = 0;         ind2 = 0
                    case(1)
                        ind1 = ngllx - 1; ind2 = 0
                    case(2)
                        ind1 = ngllx - 1; ind2 = ngllz - 1
                    case(3)
                        ind1 = 0;         ind2 = ngllz - 1
                    end select
                case(2,4)
                    select case(i)
                    case(0)
                        ind1 = 0;         ind2 = 0
                    case(1)
                        ind1 = nglly - 1; ind2 = 0
                    case(2)
                        ind1 = nglly - 1; ind2 = ngllz - 1
                    case(3)
                        ind1 = 0;         ind2 = ngllz - 1
                    end select
                end select

                if ( sauvetsurf == 1 ) then
                    tsurf = surface_gll(Tdomain, mat, numElem, numlocal, ind1, ind2)
                    Tdomain%sVertex(node(i))%tsurfsem = tsurf
                else
                    tsurf = Tdomain%sVertex(node(i))%tsurfsem
                endif
                !                tsurf = surface_gll(Tdomain, mat, numElem, numlocal, ind1, ind2)

                force(:) = force_impose(:,position + ind1 + ind2*ngll1)

                Tdomain%sVertex(node(i))%ForcesExt(:) = Tdomain%sVertex(node(i))%ForcesExt + tsurf*force
                !                Tdomain%sVertex(node(i))%FlagMka(:) = 1

            enddo !fin boucle sur les 4 noeuds de la face de couplage

            position = position + nglltot
        enddo

    end subroutine calcul_couplage_force




!!!! --------------------------------------------------
!!!! --------------------------------------------------
!!!! --------------------------------------------------
    !>
    !! \brief Calcul de la surface d'un quadrilatere quelconque correspondant a une face de couplage
    !!
    !<

    subroutine surface(Tdomain, node, tsurf)

!!!!  Description: Calcul de la surface d'un quadrilatere quelconque correspondant a une face de couplage
!!!!
!!!!  Historique: 10/08 - Gsa Ipsis - Creation
!!!!
!!!!  * Le quadrilatere est decompose en deux triangles. L'aire du quadrilatere est la somme de l'aire des deux triangles
!!!!
!!!!
!!!! --------------------------------------------------
        type (domain), intent(INOUT)  :: Tdomain
        integer node(0:3)
        real tsurf, pvec(3)
        real :: x(0:3), y(0:3), z(0:3)   !!! coordonnees des points de la surface
        integer i

        do i=0,3
            x(i)= Tdomain%Coord_nodes(0,node(i))
            y(i)= Tdomain%Coord_nodes(1,node(i))
            z(i)= Tdomain%Coord_nodes(2,node(i))
        enddo

        pvec(1) = (y(1)-y(0))*(z(3)-z(0)) - (z(1)-z(0))*(y(3)-y(0))
        pvec(2) = (z(1)-z(0))*(x(3)-x(0)) - (x(1)-x(0))*(z(3)-z(0))
        pvec(3) = (x(1)-x(0))*(y(3)-y(0)) - (y(1)-y(0))*(x(3)-x(0))

        tsurf = 0.5* sqrt( pvec(1)**2 + pvec(2)**2 + pvec(3)**2 )

        pvec(1) = (y(1)-y(2))*(z(3)-z(2)) - (z(1)-z(2))*(y(3)-y(2))
        pvec(2) = (z(1)-z(2))*(x(3)-x(2)) - (x(1)-x(2))*(z(3)-z(2))
        pvec(3) = (x(1)-x(2))*(y(3)-y(2)) - (y(1)-y(2))*(x(3)-x(2))

        tsurf = tsurf + 0.5* sqrt( pvec(1)**2 + pvec(2)**2 + pvec(3)**2 )
    end subroutine surface

    function surface_gll(Tdomain, mat, numElem, numlocal, i, j)

        type (domain), intent(IN)  :: Tdomain
        integer, intent(IN) ::  mat, numlocal, numElem
        integer, intent(IN) :: i, j
        real xi, eta, zeta
        real dsurf
        real outx, outz
        real surface_gll

        ! calcul des poids au point de gauss
        select case (numlocal)
        case(0)
            xi = Tdomain%sSubdomain(mat)%GLLcx (i)
            eta = Tdomain%sSubdomain(mat)%GLLcy (j)
            zeta = -1  !! face  0 qui verifie zeta = -1
            outx = Tdomain%sSubdomain(mat)%GLLwx(i)
            outz = Tdomain%sSubdomain(mat)%GLLwy(j)
        case(5)
            xi = Tdomain%sSubdomain(mat)%GLLcx (i)
            eta = Tdomain%sSubdomain(mat)%GLLcy (j)
            zeta = 1 !! face 5 qui verifie zeta = 1
            outx = Tdomain%sSubdomain(mat)%GLLwx(i)
            outz = Tdomain%sSubdomain(mat)%GLLwy(j)
        case  (1)
            xi = Tdomain%sSubdomain(mat)%GLLcx (i)
            eta = -1 !! face 1 qui verifie eta = -1
            zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
            outx = Tdomain%sSubdomain(mat)%GLLwx(i)
            outz = Tdomain%sSubdomain(mat)%GLLwz(j)
        case(3)
            xi = Tdomain%sSubdomain(mat)%GLLcx (i)
            eta = 1 !! face 3 qui verifie eta = 1
            zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
            outx = Tdomain%sSubdomain(mat)%GLLwx(i)
            outz = Tdomain%sSubdomain(mat)%GLLwz(j)
        case(2)
            xi = 1 !! face 2 qui verifie xi = 1
            eta = Tdomain%sSubdomain(mat)%GLLcy (i)
            zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
            outx = Tdomain%sSubdomain(mat)%GLLwy(i)
            outz = Tdomain%sSubdomain(mat)%GLLwz(j)
        case(4)
            xi = -1 !! face 4 qui verifie xi = -1
            eta = Tdomain%sSubdomain(mat)%GLLcy (i)
            zeta = Tdomain%sSubdomain(mat)%GLLcz (j)
            outx = Tdomain%sSubdomain(mat)%GLLwy(i)
            outz = Tdomain%sSubdomain(mat)%GLLwz(j)
        end select
        dsurf = calcule_surf(Tdomain, numElem, numlocal)*.25
        surface_gll = dsurf*outx*outz

    end function surface_gll


    subroutine cross_(u, v, w)
        real, dimension(3), intent(in) :: u, v
        real, dimension(3), intent(out) :: w
        w(1) =   u(2)*v(3) - u(3)*v(2)
        w(2) = - u(1)*v(3) + u(3)*v(1)
        w(3) =   u(1)*v(2) - u(2)*v(1)
    end subroutine cross_

    function calcule_surf(Tdomain, numElem, numlocal)
        type (domain), intent(IN)  :: Tdomain
        integer, intent(IN) ::  numlocal, numElem
        real,dimension(3, 0:7) :: X
        real,dimension(3) :: a, b, c, d, u, v
        real :: calcule_surf
        integer :: i, i_aux

        !! recuperation des coordonnees des sommets de la maille Sem contenant la face de couplage
        do i=0,7
            i_aux = Tdomain%specel(numElem)%Control_Nodes(i)
            X(:,i) = Tdomain%Coord_Nodes(:,i_aux)
        enddo

        ! On calcule la somme (signee) des aires des deux triangles
        select case (numlocal)
        case(0)
            ! face 0, 1, 2, 3
            a = X(:,1)-X(:,0)
            b = X(:,3)-X(:,0)
            c = X(:,3)-X(:,2)
            d = X(:,1)-X(:,2)
        case(5)
            ! face 4 5 6 7
            a = X(:,5)-X(:,4)
            b = X(:,7)-X(:,4)
            c = X(:,7)-X(:,6)
            d = X(:,5)-X(:,6)
        case  (1)
            ! face 0 1 5 4
            a = X(:,1)-X(:,0)
            b = X(:,4)-X(:,0)
            c = X(:,4)-X(:,5)
            d = X(:,1)-X(:,5)
        case(3)
            ! face 2 3 7 6
            a = X(:,3)-X(:,2)
            b = X(:,6)-X(:,2)
            c = X(:,6)-X(:,7)
            d = X(:,3)-X(:,7)
        case(2)
            ! face 1 2 6 5
            a = X(:,2)-X(:,1)
            b = X(:,5)-X(:,1)
            c = X(:,5)-X(:,6)
            d = X(:,2)-X(:,6)
        case(4)
            ! face 0 3 7 4
            a = X(:,3)-X(:,0)
            b = X(:,4)-X(:,0)
            c = X(:,4)-X(:,7)
            d = X(:,3)-X(:,7)
        end select
        call cross_(a, b, u)
        call cross_(c, d, v)
        u = (u+v)/2.
        calcule_surf = sqrt( u(1)**2+u(2)**2+u(3)**2 )
    end function calcule_surf

    !>
    !! \brief Remplissage de la structure face_couplage
    !! - Definition du numero global des faces Sem de couplage, des numeros des noeuds les formant
    !! - des positions de ces noeuds, du numero de proc Sem de la face
    !!
    !<

    subroutine remplit_face_couplage(Tdomain, MaxNgParFace, MaxNgParDir)

        use ref_orient

        type (domain), intent(IN)  :: Tdomain
        integer, intent(out) :: MaxNgParFace, MaxNgParDir
        integer :: numFace, numElem, ngll1, ngll2, numlocal
        integer :: i, iface, j, iel, k

        do iface=1,Nbface
            numFace = comm_couplage%m_tabFaceCouplage(iface)

            ngll1 = Tdomain%sFace(numFace)%ngll1
            ngll2 = Tdomain%sFace(numFace)%ngll2

            MaxNgParFace = max(MaxNgParFace,ngll1*ngll2)
            MaxNgParDir = max(max(MaxNgParDir,ngll1),ngll2)

            face_couplage(iface)%face   = numFace
            face_couplage(iface)%proc   = Tdomain%rank
            face_couplage(iface)%ngll1   = ngll1  !3D Gsa
            face_couplage(iface)%ngll2   = ngll2  !3D Gsa

            do iel=0,Tdomain%n_elem - 1
                do k=0,5
                    if(Tdomain%specel(iel)%Near_Faces(k) == numFace) then
                        numElem = iel
                        exit
                    endif
                enddo
            enddo
            face_couplage(iface)%elem = numElem

            !un seul element suffit pour gerer la face. on n'a pas besoin du 2eme element contenant la face
            !alternative
            !    face_couplage(iface)%numElem   = !3D Gsa
            do i=0,5 !3D Gsa
                if(Tdomain%specel(numElem)%Near_Faces(i) == numFace) then
                    numlocal = i
                    exit
                endif
            enddo
            face_couplage(iface)%numlocal = numlocal

            do i=0,3 !3D Gsa
                face_couplage(iface)%noeud(i) = Tdomain%specel(numElem)%Near_Vertices(NOEUD_FACE(numlocal,i))
            enddo

            do j=0,2 !3D Gsa
                do i=0,3 !3D Gsa
                    face_couplage(iface)%coord(i,j) = Tdomain%Coord_nodes(j,face_couplage(iface)%noeud(i))
                enddo
            enddo
        enddo

    end subroutine remplit_face_couplage

    !>
    !! \brief Calcul de la plus petite distance entre particules de couplage Mka pour une meme face
    !! de couplage Sem.
    !<

    subroutine dist_min_point_de_couplage(n, Xpdc, Ypdc, Zpdc, dmin)

        integer, intent(in) :: n
        real, intent(in) :: Xpdc(n), Ypdc(n), Zpdc(n)
        integer :: iface, np1, np2
        real :: dist, dmin(Nbface)

        if(Nbface .EQ. 0) return
        dmin(:) = huge(1.)

        do iface=1,Nbface
            do np1=1,n
                if(comm_couplage%m_numFace(np1) == face_couplage(iface)%face ) then
                    do np2=np1+1,n
                        if(comm_couplage%m_numFace(np2) == comm_couplage%m_numFace(np1) ) then
                            dist = (Xpdc(np1) - Xpdc(np2))**2 + (Ypdc(np1) - Ypdc(np2))**2 + (Zpdc(np1) - Zpdc(np2))**2
                            dist = sqrt(dist)
                            dmin(iface) = min(dmin(iface), dist)
                        endif
                    enddo
                endif
            enddo
        enddo
    end subroutine dist_min_point_de_couplage


    !>
    !! \brief Definition du nombre de points d'interpolation dans une direction pour une face de couplage Sem
    !!
    !! Ce nombre depend du nombre de points de Gauss et de la plus petite distance
    !! entre les particules de couplage.
    !<


    subroutine definit_nb_pts_Pk(Tdomain, dmin_couplage)

        type (domain), intent(IN)  :: Tdomain
        real, intent(in) :: dmin_couplage(Nbface)
        integer ngll1, ngll2, numFace, numElem, numlocal
        integer i, iarete, iface, ipoint, j
        integer inoeud(0:1)
        real,dimension(0:7) :: x, y, z
        !    real Lface(0:11)
        real Lface(0:3)
        real hk
        integer,dimension(0:11,0:1) :: NOEUD_ARETE !Gsa
        integer,dimension(0:5,0:3) :: ARETE_FACE !Gsa


        !tableau de correspondance num_local de face/numero des aretes de l'element
        ARETE_FACE(0,:) = (/ 0, 1, 2, 3 /) !face 0 - k=0
        ARETE_FACE(1,:) = (/ 0, 4, 5, 6 /) !face 1 - j=0
        ARETE_FACE(2,:) = (/ 1, 7, 8, 4 /) !face 2 - i=n-1
        ARETE_FACE(3,:) = (/ 2, 7, 9, 10 /) !face 3 - j=n-1
        ARETE_FACE(4,:) = (/ 3, 10, 11, 6 /) !face 4 - i=0
        ARETE_FACE(5,:) = (/ 5, 8, 9, 11 /) !face 5 - k=n-1

        !tableau de correspondance num_local de face/numero des noeuds de l'element
        NOEUD_ARETE(0,:) =(/1,2/)
        NOEUD_ARETE(1,:) =(/2,3/)
        NOEUD_ARETE(2,:) =(/4,3/)
        NOEUD_ARETE(3,:) =(/1,4/)
        NOEUD_ARETE(4,:) =(/2,6/)
        NOEUD_ARETE(5,:) =(/5,6/)
        NOEUD_ARETE(6,:) =(/1,5/)
        NOEUD_ARETE(7,:) =(/3,7/)
        NOEUD_ARETE(8,:) =(/6,7/)
        NOEUD_ARETE(9,:) =(/8,7/)
        NOEUD_ARETE(10,:)=(/4,8/)
        NOEUD_ARETE(11,:)=(/5,8/)

        do iface=1,Nbface
            numFace = comm_couplage%m_tabFaceCouplage(iface)
            ngll1 = Tdomain%sFace(numFace)%ngll1
            ngll2 = Tdomain%sFace(numFace)%ngll2
            numElem = face_couplage(iface)%elem
            numlocal = face_couplage(iface)%numlocal
            do i=0,7
                ipoint = Tdomain%specel(numElem)%Control_Nodes(i)
                x(i)= Tdomain%Coord_nodes(0,ipoint)
                y(i)= Tdomain%Coord_nodes(1,ipoint)
                z(i)= Tdomain%Coord_nodes(2,ipoint)
            enddo

            do iarete=0,3
                do j=0,1
                    inoeud(j) = NOEUD_ARETE(ARETE_FACE(numlocal,iarete),j) - 1
                enddo
                Lface(iarete) = sqrt( (x(inoeud(0)) - x(inoeud(1)) )**2 &
                    + ( y(inoeud(0)) - y(inoeud(1)) ) **2 &
                    + ( z(inoeud(0)) - z(inoeud(1)) ) **2 )
            enddo
            hk = min(dmin_couplage(iface), minval(Lface)/(2.*max(ngll1,ngll2) + 1.))
            !! minimum 5 points d'interpolation par longueur caracteristique
            face_couplage(iface)%nbPk = 3*(floor( maxval(Lface)/hk) + 1)
            !hk = minval(Lface)/tab_Pk(iface)%nb_pts !correction
        enddo
    end subroutine definit_nb_pts_Pk


    !>
    !! \brief Reception des surfaces des particules de couplage Mka.
    !!
    !<

    subroutine reception_surface_part_mka(Tdomain)
        use sdomain
        implicit none

        type (domain), intent(IN)  :: Tdomain
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer :: tag, ierr

        ! a faire une seule fois
        if(comm_couplage%m_nb_point_de_couplage > 0) then
            tag = 510000 + 1000*Tdomain%rank + 1
            call MPI_RECV (comm_couplage%m_surf_part,comm_couplage%m_nb_point_de_couplage, &
                MPI_DOUBLE_PRECISION, Tdomain%master_superviseur, &
                tag, Tdomain%communicateur_global, status, ierr )
        endif

        deallocate(Xpdc)
        deallocate(Ypdc)
        deallocate(Zpdc)


    end subroutine reception_surface_part_mka

!!!!


!!!!  Description: Definition du second membre pour la vitesse dans le systeme lineaire
!!!!
!!!! Simple recuperation des vitesses aux points de Gauss
!!!! Envoi des donnees au Superviseur
!!!! --------------------------------------------------
    subroutine envoi_vitesse_mka(Tdomain,ntime)

        use ref_orient

        type (domain), intent(INOUT)  :: Tdomain
        integer, intent(IN) :: ntime

        integer :: tag, ierr
        integer :: i,j,numElem, numFace, numlocal, mat
        integer :: ngll
        integer :: numVertex(0:3), numEdge(0:3), numloc_edge
        integer ngll1, ngll2, ngllx, nglly
        integer ind1, ind2
        integer iarete, inod, iface
        integer idim
        real tsurf

        real,dimension(:,:),pointer :: vecu, vecu_tmp
#if 0

        allocate(vecu(comm_couplage%m_MaxNgParFace,comm_couplage%m_dim))
        allocate(vecu_tmp(comm_couplage%m_MaxNgParFace,comm_couplage%m_dim))

        do iface=1,comm_couplage%m_NbFace
            numFace=comm_couplage%m_tabFaceCouplage(iface)
            numElem = face_couplage(iface)%elem
            mat = Tdomain%specel(numElem)%mat_index
            numlocal = face_couplage(iface)%numlocal
            ngll1 = Tdomain%sFace(numFace)%ngll1
            ngll2 = Tdomain%sFace(numFace)%ngll2
            ngll = ngll1 * ngll2
            ngllx = Tdomain%specel(numElem)%ngllx
            nglly = Tdomain%specel(numElem)%nglly
            do iarete=0,3
                numloc_edge = ARETE_FACE(numlocal,iarete) !numero local de l'arete (entre 0 et 11)
                numEdge(iarete) = Tdomain%specel(numElem)%Near_Edges(numloc_edge)
            enddo

            do inod=0,3
                numVertex(inod)= Tdomain%specel(numElem)%Near_Vertices(NOEUD_FACE(numlocal,inod))
            enddo

            do idim=1,comm_couplage%m_dim
                vecu(1:ngll,idim)=0.
                do j=1,ngll2-2
                    do i=1,ngll1-2
                        vecu(i+j*ngll1+1,idim) = Tdomain%sFace(numFace)%Veloc(i,j,idim-1)
                    enddo
                enddo
                do iarete=0,3
                    ngll = Tdomain%sEdge(numEdge(iarete))%ngll
                    numloc_edge = ARETE_FACE(numlocal,iarete) !numero local de l'arete (entre 0 et 11)

                    do i=1,ngll-2
                        if(iarete==0) then
                            ind1 = i
                            ind2 = 0
                        elseif(iarete==1) then
                            ind1 = ngll1-1
                            ind2 = i
                        elseif(iarete==2) then
                            ind1 = i
                            ind2 = ngll2-1
                        elseif(iarete==3) then
                            ind1 = 0
                            ind2 = i
                        endif
                        if(Tdomain%specel(numElem)%Orient_Edges(numloc_edge) == 0 ) then
                            vecu(ind1+ind2*ngll1+1,idim) = Tdomain%sEdge(numEdge(iarete))%Veloc(i,idim-1)
                        else
                            vecu(ind1+ind2*ngll1+1,idim) = Tdomain%sEdge(numEdge(iarete))%Veloc(ngll-1-i,idim-1)
                        endif
                    enddo
                enddo

                do inod=0,3

                    select case(inod)
                    case(0)
                        ind1 = 0
                        ind2 = 0
                    case(1)
                        ind1 = ngllx - 1
                        ind2 = 0
                    case(2)
                        ind1 = ngllx - 1
                        ind2 = nglly - 1
                    case(3)
                        ind1 = 0
                        ind2 = nglly - 1
                    end select
                    vecu(ind1+ind2*ngll1+1,idim) =  Tdomain%sVertex(numVertex(inod))%Veloc(idim-1)
                enddo
            enddo

            !! si CAS  M^-1 (D^t A D) u
            vecu_tmp = vecu
            do idim=1,comm_couplage%m_dim
                do j=0,ngll2-1
                    do i=0,ngll1-1
                        tsurf = surface_gll(Tdomain, mat, numElem, numlocal, i, j)
                        vecu_tmp(i+j*ngll1+1,idim) = vecu_tmp(i+j*ngll1+1,idim)*tsurf
                    enddo
                enddo
            enddo
            ngll = ngll1 * ngll2
            do idim=1,comm_couplage%m_dim

                tag = 730000+3*(iface-1)+(idim-1)
                call MPI_SEND(vecu_tmp(1,idim), ngll, MPI_DOUBLE_PRECISION, Tdomain%master_superviseur,&
                    tag, Tdomain%communicateur_global, ierr )

            enddo
            !!       enddo
        enddo

        deallocate(vecu)
        deallocate(vecu_tmp)
#endif
    end subroutine envoi_vitesse_mka


    subroutine reception_nouveau_pdt_sem(Tdomain)

        use sdomain
        implicit none
        type (Domain), intent (INOUT) :: Tdomain
        integer :: n, mat
        real dtsem
        integer :: tag,ierr
        integer, dimension (MPI_STATUS_SIZE) :: status

        tag = 4100000+10000*Tdomain%rank
        call MPI_RECV(dtsem, 1, MPI_DOUBLE_PRECISION, Tdomain%master_superviseur,tag, &
            Tdomain%communicateur_global, status,ierr)


        do n=0,Tdomain%n_elem - 1
            mat = Tdomain%specel(n)%mat_index
            Tdomain%sSubdomain(mat)%Dt = dtsem
        enddo
        Tdomain%TimeD%dtmin = dtsem
    end subroutine reception_nouveau_pdt_sem

#endif

end module scouplage


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
