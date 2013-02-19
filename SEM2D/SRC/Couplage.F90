!>
!!\file Couplage.F90
!!\brief Contient les routines propres au traitement du couplage.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module scouplage

#ifdef MKA3D

    use sdomain
    use semdatafiles
    use mpi
    implicit none


    type :: type_comm_couplage

       integer :: m_dim                                   ! dimension de l espace
       integer :: m_nb_point_de_couplage                  ! nombre de particules de couplage sur le proc courant
       integer :: m_NbFace                                ! nombre de face de couplage sur le proc courant
       integer :: m_MaxNgParFace                          ! nombre max de pt de Gauss sur une face de couplage du proc courant
       integer :: m_gl_Ngauss                             ! nombre total max  de points de gauss sur tt la surface de couplage
       integer,dimension(:),pointer  :: m_numFace         ! numeros des faces de couplage
       integer,dimension(:),pointer  :: m_numMaille       ! numeros des mailles sem contenant les faces de couplage
       integer,dimension(:),pointer  :: m_tabFaceCouplage ! liste des faces de couplage
       real,dimension(:),pointer     :: m_pos_proj        ! abscisse curviligne sur la face de couplage du projete de la particule de couplage
       integer,dimension(:),pointer  :: m_Pk_plus_proche  ! numero local du Pk le plus proche
       real,dimension(:), pointer    :: m_surf_part       ! surface des particules Mka de couplage

    end type type_comm_couplage

    type :: type_point_a_interpoler

       real,dimension(:),pointer :: Interp_Coeff

    end type type_point_a_interpoler


    type :: type_face_couplage

       integer :: face   ! numero de la face de couplage sem dans la numerotation du proc sem qui la contient
       integer :: proc   ! numero du processeur sem qui contient cette face de couplage
       integer :: ngll   ! nbre de point de gauss sur la face
       integer :: noeud0 ! numero du premier noeud de la face
       integer :: noeud1 ! numero du 2ieme noeud de la face
       real,dimension(2) :: coord0 ! coordonnees noeud0
       real,dimension(2) :: coord1 ! coordonnees noeud1

    end type type_face_couplage


    integer :: nbFaceTotal,nbFace

    type(type_face_couplage), dimension(:),pointer :: face_couplage

    type(type_comm_couplage) :: comm_couplage

    type :: sommet_gll
       real, dimension(0:2) :: coord
       integer nbface
       integer, dimension(8) :: numface !au plus un noeud est connecte a 8 faces
       integer num, numglobal
    end type sommet_gll

    type :: type_interp_Pk_face
       integer :: nb_pts
       real, dimension(:,:), pointer :: Interp_coeff
       integer, dimension(:), pointer :: part_mka_plus_proche
       real, dimension(:), pointer :: xi, eta
       integer, dimension(:), pointer :: pdg_plus_proche
       real :: surf
    end type type_interp_Pk_face

    type(type_interp_Pk_face), dimension(:),pointer :: tab_Pk
    real,dimension(:), pointer :: Xpdc, Zpdc

contains


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
        allocate(comm_couplage%m_pos_proj(comm_couplage%m_nb_point_de_couplage))

        !il y a au maximum autant de face de couplage que de points de couplage
        allocate(comm_couplage%m_tabFaceCouplage(comm_couplage%m_nb_point_de_couplage))
        allocate(comm_couplage%m_surf_part(comm_couplage%m_nb_point_de_couplage))
    end subroutine alloue_couplage


    !>
    !! \brief Desalloue les tableaux utilisés lors du couplage
    !!
    !<
    subroutine desalloue_couplage()

        deallocate(comm_couplage%m_numFace)
        deallocate(comm_couplage%m_numMaille)
        deallocate(comm_couplage%m_surf_part)
        deallocate(comm_couplage%m_tabFaceCouplage)

    end subroutine desalloue_couplage
    !>
    !! \brief Initialise les données propres au couplage.
    !! Effectue notamment des communications avec le superviseur
    !! afin d'envoyer le pas de temps et de recevoir des informations sur le
    !! nombre de particules de couplage, les mailles couplées,...
    !! \param type (domain), intent(INOUT) Tdomain
    !<
    subroutine initialisation_couplage(Tdomain, MaxNgparFace)
        type (domain), intent(INOUT)  :: Tdomain

        integer :: i,ipoint
        integer :: tag,ierr
        integer :: wf,ngll,np
        integer :: numFace,numElem,mat
        integer :: iface,bufsize,decal,nbchamps
        integer :: Ngauss,MaxNgParFace,gl_MaxNgParFace
        real :: xi,eta
        logical :: dejaPresent

        integer, dimension (MPI_STATUS_SIZE) :: status
        integer, dimension(:),pointer:: tabNbFace
        integer,dimension(:), pointer :: idiag
        integer,dimension(:),pointer :: displs,count
        integer,dimension(5)  :: ibuf
        real, dimension (:),pointer :: buf
        real,dimension (:), pointer :: acol
        real,dimension(0:3) :: x,z

        real,dimension (:), pointer :: dmin_couplage
        logical,dimension(:), pointer :: face_deja_traitee
        character(Len=MAX_FILE_SIZE) fnamef

        allocate(displs(Tdomain%Mpi_var%n_proc),count(Tdomain%Mpi_var%n_proc))

        ! envoi du pas de temps sem
        if (Tdomain%MPI_var%my_rank == 0) then
            tag=500000
            call MPI_SEND(Tdomain%TimeD%dtmin,1,MPI_DOUBLE_PRECISION, Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )
        endif

        ! reception du nombre de faces et de mailles de couplage
        tag=1000000+Tdomain%MPI_var%my_rank
        call MPI_RECV(comm_couplage%m_nb_point_de_couplage,1, MPI_INTEGER, Tdomain%master_superviseur,tag, &
            Tdomain%communicateur_global, status,ierr)

        call alloue_couplage(3)

        if (comm_couplage%m_nb_point_de_couplage.gt.0) then


            ! reception du tableau des mailles de couplage
            tag=1100000+Tdomain%MPI_var%my_rank
            call MPI_RECV(comm_couplage%m_numMaille,comm_couplage%m_nb_point_de_couplage, MPI_INTEGER, Tdomain%master_superviseur,tag, &
                Tdomain%communicateur_global, status,ierr)

            ! reception du tableau des faces de couplage
            tag=1200000+Tdomain%MPI_var%my_rank
            call MPI_RECV(comm_couplage%m_numFace,comm_couplage%m_nb_point_de_couplage, MPI_INTEGER, Tdomain%master_superviseur,tag, &
                Tdomain%communicateur_global, status,ierr)

            ! reception du tableau des projetes de couplage
            tag=1300000+Tdomain%MPI_var%my_rank
            call MPI_RECV(comm_couplage%m_pos_proj,comm_couplage%m_nb_point_de_couplage,MPI_DOUBLE_PRECISION,Tdomain%master_superviseur,&
                tag, Tdomain%communicateur_global, status,ierr)

        endif

        Nbface=0

        allocate (Xpdc(comm_couplage%m_nb_point_de_couplage))
        allocate (Zpdc(comm_couplage%m_nb_point_de_couplage))
        call calcule_coord_glob(Tdomain, comm_couplage%m_nb_point_de_couplage, Xpdc, Zpdc)


        ! calcul des fonctions de forme en chaque point de couplage
        do np=1,comm_couplage%m_nb_point_de_couplage
            ! recuperation de la face et de l element du point de couplage courant et de ses coordonnees, du numero de materiau et du nombre de pdg en x et z
            numFace = comm_couplage%m_numFace(np)
            numElem = comm_couplage%m_numMaille(np)

            mat = Tdomain%specel(numElem)%mat_index
            ngll = Tdomain%sFace(numFace)%ngll

            ! recuperation des coor des 4 sommets de l'element
            do i=0,3
                ipoint = Tdomain%specel(numElem)%Control_Nodes(i)
                x(i)= Tdomain%Coord_nodes(0,ipoint)
                z(i)= Tdomain%Coord_nodes(1,ipoint)
            enddo

            if ( Tdomain%sface(numFace)%Near_Element(0) == numElem ) then
                wf = Tdomain%sface(numFace)%Which_face(0)
            elseif ( Tdomain%sface(numFace)%Near_Element(1) == numElem ) then
                wf = Tdomain%sface(numFace)%Which_face(1)
            endif


            ! calcul des coordonnees (xi,eta) du point de couplage
            if (wf == 0) then
                eta=-1
                xi = calcul_xi_eta(Xpdc(np),Zpdc(np),x(0),z(0),x(1),z(1))
            else if  (wf == 1) then
                ! calcule eta,xi=1
                xi=1
                eta = calcul_xi_eta(Xpdc(np),Zpdc(np),x(1),z(1),x(2),z(2))
            else if  (wf == 2) then
                ! calcule xi,eta=1
                eta=1
                xi = calcul_xi_eta (Xpdc(np),Zpdc(np),x(2),z(2),x(3),z(3))
            else if  (wf == 3) then
                ! calcule eta,xi=-1
                xi=-1
                eta = calcul_xi_eta(Xpdc(np),Zpdc(np),x(3),z(3),x(0),z(0))
            endif

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


        allocate (tab_Pk(Nbface))
        allocate(tabNbFace(Tdomain%Mpi_var%n_proc))
        tabNbFace(:)=0

        !call MPI_Gather (Nbface,1,MPI_INTEGER,tabNbFace,1,MPI_INTEGER,0, Tdomain%communicateur, ierr)
        if (Tdomain%MPI_var%my_rank == 0) then
            count(:)=1
            displs(1)=0
            do i=2,Tdomain%Mpi_var%n_proc
                displs(i)=displs(i-1)+count(i-1)
            enddo
        endif

        ! tabNbFace,defini uniquement sur proc 0,  contient le nbre de face de couplage de chaque proc
        call MPI_Gatherv (Nbface,1,MPI_INTEGER,tabNbFace,count,displs,MPI_INTEGER,0, Tdomain%communicateur, ierr)


        ! calcul du nombre de face de couplage total
        if (Tdomain%MPI_var%my_rank == 0) then
            NbFaceTotal=sum(tabNbFace)
            allocate(face_couplage(NbFaceTotal))
        else
            allocate(face_couplage(Nbface))
        endif

        call MPI_Bcast(NbFaceTotal,1,MPI_INTEGER,0,Tdomain%communicateur,ierr)

        ! Determination du nbre max de points de gauss par maille et du nombre total de points de gauss de la zone de couplage sur chaque proc
        MaxNgParFace=0 ! nombre max de points de gauss par bord de maille sur le proc courant


        ! remplissage de la structure face_couplage
        call remplit_face_couplage(Tdomain, MaxNgParFace)

        ! on conserve le nombre max de gauss par face de couplage du proc courant
        comm_couplage%m_MaxNgParFace=MaxNgParFace

        allocate (dmin_couplage(Nbface))
        call dist_min_point_de_couplage(comm_couplage%m_nb_point_de_couplage, Xpdc, Zpdc, dmin_couplage)
        call definit_nb_pts_Pk(Tdomain, dmin_couplage)
        deallocate(dmin_couplage)

        !! liste des tous les points d'interpolation du peigne et de leur projete Mka

        call semname_couplage_listepts(Tdomain%MPI_var%my_rank,fnamef)
        open(unit=75,file=fnamef,status='unknown')

        allocate(face_deja_traitee(Nbface))
        face_deja_traitee = .false.
        do np=1,comm_couplage%m_nb_point_de_couplage
            ! recuperation de la face et de l element du point de couplage courant et de ses coordonnees, du numero de materiau et du nombre de pdg en x et z
            numFace = comm_couplage%m_numFace(np)
            !on cherche iface
            do iface=1,Nbface
                if(face_couplage(iface)%face == numFace) exit
            enddo
            if(face_deja_traitee(iface)) cycle

            numElem = comm_couplage%m_numMaille(np)

            mat = Tdomain%specel(numElem)%mat_index

            call traitement_Pk(Tdomain, iface, numElem, numFace, mat)
            face_deja_traitee(iface) = .true.
        enddo
        deallocate(face_deja_traitee)
        close(75)

        allocate(comm_couplage%m_Pk_plus_proche(comm_couplage%m_nb_point_de_couplage))
        !!liste de tous les projetes Mka et de leur point d'interpolation du peigne le plus proche
        call semname_couplage_listeproj(Tdomain%MPI_var%my_rank,fnamef)
        open(unit=77,file=fnamef,status='unknown')
        do np=1,comm_couplage%m_nb_point_de_couplage
            comm_couplage%m_Pk_plus_proche(np) = trouve_Pk_plus_proche(Tdomain, comm_couplage%m_nb_point_de_couplage, np, Xpdc, Zpdc)
        enddo
        close(77)

        ! calcul (sur proc 0) du nombre max total de maille sur une face de couplage
        call MPI_Reduce (MaxNgParFace,gl_MaxNgParFace,1,MPI_INTEGER,MPI_MAX,0, Tdomain%communicateur, ierr)

        ! on collecte sur le proc 0 un ensemble des donnees pour ttes les faces de couplage sur les differents proc
        NbChamps = 9  ! taille des donnees pour une face_couplage
        if (Tdomain%MPI_var%my_rank == 0) then
            ! reception par le proc maitre des faces de couplage
            decal=0
            do i=1,Tdomain%Mpi_var%n_proc-1
                decal = decal+tabNbFace(i)
                if (tabNbFace(i+1)>0) then
                    bufsize = NbChamps*tabNbFace(i+1)
                    allocate(buf(bufsize))
                    tag = 70000+i
                    call MPI_RECV(buf,bufsize, MPI_DOUBLE_PRECISION, i,tag, Tdomain%communicateur, status,ierr)

                    do iface=1,tabNbFace(i+1)
                        face_couplage(iface+decal)%face   = buf(1+NbChamps*(iface-1))
                        face_couplage(iface+decal)%proc   = buf(2+NbChamps*(iface-1))
                        face_couplage(iface+decal)%ngll   = buf(3+NbChamps*(iface-1))
                        face_couplage(iface+decal)%noeud0 = buf(4+NbChamps*(iface-1))
                        face_couplage(iface+decal)%noeud1 = buf(5+NbChamps*(iface-1))
                        face_couplage(iface+decal)%coord0(1) = buf(6+NbChamps*(iface-1))
                        face_couplage(iface+decal)%coord0(2) = buf(7+NbChamps*(iface-1))
                        face_couplage(iface+decal)%coord1(1) = buf(8+NbChamps*(iface-1))
                        face_couplage(iface+decal)%coord1(2) = buf(9+NbChamps*(iface-1))
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
                    buf(3+NbChamps*(iface-1)) = face_couplage(iface)%ngll
                    buf(4+NbChamps*(iface-1)) = face_couplage(iface)%noeud0
                    buf(5+NbChamps*(iface-1)) = face_couplage(iface)%noeud1
                    buf(6+NbChamps*(iface-1)) = face_couplage(iface)%coord0(1)
                    buf(7+NbChamps*(iface-1)) = face_couplage(iface)%coord0(2)
                    buf(8+NbChamps*(iface-1)) = face_couplage(iface)%coord1(1)
                    buf(9+NbChamps*(iface-1)) = face_couplage(iface)%coord1(2)
                enddo

                tag = 70000+Tdomain%MPI_var%my_rank
                call MPI_SEND(buf,bufsize,MPI_DOUBLE_PRECISION, 0,tag, Tdomain%communicateur, ierr )

                deallocate(buf)
            endif

        endif ! finsi proc non maitre


        Ngauss=0 ! nombre max de points de gauss sur la zone de couplage du proc courant
        do iface=1,NbFace
            Ngauss = Ngauss + face_couplage(iface)%ngll
        enddo

        call MPI_AllReduce (Ngauss,comm_couplage%m_gl_Ngauss,1,MPI_INTEGER,MPI_SUM, Tdomain%communicateur, ierr)
        write(6,*) 'NGAUSS= ',comm_couplage%m_gl_Ngauss



        ! envoi au superviseur
        if (Tdomain%MPI_var%my_rank == 0) then
            tag=600000
            ibuf(1) = 2
            ibuf(2) = comm_couplage%m_gl_Ngauss
            ibuf(3) = gl_MaxNgParFace
            ibuf(4) = 0
            ibuf(5) = 0
            call MPI_SEND(ibuf,5,MPI_INTEGER, Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )
        endif


        ibuf(1) = Nbface
        ibuf(2) = MaxNgParFace
        ibuf(3) = 0 ! MaxNgParDir
        ibuf(4) = Ngauss

        ! envoi au superviseur du nb de face de couplage  et de MaxNgParFace, pour chaque proc
        tag=610000+Tdomain%MPI_var%my_rank
        call MPI_SEND(ibuf, 4, MPI_INTEGER, Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )

        ! allocation de la matrice elementaire
        allocate(idiag(0:MaxNgParFace))
        allocate(acol(MaxNgParFace*(MaxNgParFace+1)/2))

        ! construction de la matrice elementaire de chaque face de couplage et envoi au superviseur, pour chaque proc sem
        do iface=1,Nbface

            numFace=comm_couplage%m_tabFaceCouplage(iface)
            ngll = Tdomain%sFace(numFace)%ngll

            ibuf(1) = ngll ! ngll1
            ibuf(2) = 0    ! ngll2
            ibuf(3) = tab_Pk(iface)%nb_pts
            ibuf(4) = 0

            tag=670000+Tdomain%MPI_var%my_rank+iface-1
            call MPI_SEND(ibuf, 4, MPI_INTEGER, Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )

            tag=630000+Tdomain%MPI_var%my_rank+iface-1
            call MPI_SEND(tab_Pk(iface)%Interp_coeff(1,0), tab_Pk(iface)%nb_pts*ngll, MPI_DOUBLE_PRECISION, &
                Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )

            tag=640000+Tdomain%MPI_var%my_rank+iface-1
            call MPI_SEND(tab_Pk(iface)%part_mka_plus_proche,tab_Pk(iface)%nb_pts, MPI_INTEGER, Tdomain%master_superviseur,tag, &
                Tdomain%communicateur_global, ierr )

        enddo ! fin boucle sur les faces

        comm_couplage%m_surf_part = 0.


        ! Pas de protection/reprise a l'initiative de sem. La protection est pilotee par le superviseur
        Tdomain%logicD%save_restart =.false.
        Tdomain%logicD%run_restart =.false.
        Tdomain%TimeD%iter_reprise =0
        Tdomain%TimeD%prot_m0 =0
        Tdomain%TimeD%prot_m1 =0
        Tdomain%TimeD%prot_m2 =0


        ! reception de l iteration de reprise
        tag=900000+Tdomain%MPI_var%my_rank
        call MPI_RECV(Tdomain%TimeD%iter_reprise,1, MPI_INTEGER, Tdomain%master_superviseur,tag, &
            Tdomain%communicateur_global, status,ierr)

        if (Tdomain%TimeD%iter_reprise > 0) then
            Tdomain%logicD%run_restart = .true.
        endif

        ! faire toutes les deallocations possibles
        deallocate(tabNbFace)
        deallocate(acol)
        deallocate(idiag)
        deallocate(displs)
        deallocate(count)
        deallocate(comm_couplage%m_pos_proj)

    end subroutine initialisation_couplage

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




    !>
    !! \brief Calcule les positions des projetes des particules de couplage de Mka
    !!
    !! \param type (domain), intent(INOUT) Tdomain
    !! \param integer, intent(IN) numFace
    !! \param integer, intent(IN) numElem
    !! \param real, intent(IN) pos_proj
    !! \param real, intent(INOUT) Xpdc
    !! \param real, intent(INOUT) Zpdc
    !<
    subroutine calcule_coord(Tdomain,numFace,numElem,pos_proj,Xpdc,Zpdc)

        type (domain), intent(IN)  :: Tdomain

        integer, intent(IN) :: numFace,numElem
        real, intent(IN)    :: pos_proj
        real, intent(INOUT)   :: Xpdc,Zpdc

        integer :: ipoint1,ipoint2
        real :: x1,y1,x2,y2


        if (Tdomain%specel(numElem)%near_face(0)==numFace) then
            ipoint1 = Tdomain%specel(numElem)%Control_Nodes(0)
            ipoint2 = Tdomain%specel(numElem)%Control_Nodes(1)

        elseif  (Tdomain%specel(numElem)%near_face(1)==numFace) then
            ipoint1 = Tdomain%specel(numElem)%Control_Nodes(1)
            ipoint2 = Tdomain%specel(numElem)%Control_Nodes(2)

        elseif  (Tdomain%specel(numElem)%near_face(2)==numFace) then
            ipoint1 = Tdomain%specel(numElem)%Control_Nodes(2)
            ipoint2 = Tdomain%specel(numElem)%Control_Nodes(3)

        elseif  (Tdomain%specel(numElem)%near_face(3)==numFace) then
            ipoint1 = Tdomain%specel(numElem)%Control_Nodes(3)
            ipoint2 = Tdomain%specel(numElem)%Control_Nodes(0)
        else
            print*,"calcule_coord:pb de coherence entre numElem et numFace"
            stop 1
        endif

        ! coordonnees du projete
        x1 = Tdomain%Coord_nodes(0,ipoint1);y1=Tdomain%Coord_nodes(1,ipoint1)
        x2 = Tdomain%Coord_nodes(0,ipoint2);y2=Tdomain%Coord_nodes(1,ipoint2)
        Xpdc = x1 + pos_proj * (x2 - x1)
        Zpdc = y1 + pos_proj * (y2 - y1)

    end subroutine calcule_coord


    !>
    !! \brief Evalue les forces sur les particules couplées.
    !! -Recoit les forces MKA3d appliquées aux particules de couplage,
    !! -calcul le second membre du système couplé à résoudre,
    !! -recoit la solution du système et applique les forces en conséquence.
    !!
    !!
    !! \param type (domain), intent(INOUT) Tdomain
    !! \param integer, intent(IN) ntime
    !<

    subroutine calcul_couplage_force(Tdomain,ntime)

        type (domain), intent(INOUT)  :: Tdomain
        integer, intent(IN) :: ntime


        integer :: i,j,ipoint,numElem,numFace,mat
        integer :: ngll, wf
        integer :: iface,n0,n1,position
        integer :: tag,ierr
        integer, dimension (MPI_STATUS_SIZE) :: status

        real :: out,dist
        real :: forcex,forcey
        real,dimension(0:3) :: x,z
        real, dimension(comm_couplage%m_dim, comm_couplage%m_gl_Ngauss) :: force_impose

        ! RECEPTION DES FORCES DU SUPERVISEUR

        ! reception du vecteur force appliquees au pt de gauss de sem, solution du systeme resolu par le superviseur
        tag=1500000+10000*Tdomain%MPI_var%my_rank+ntime
        call MPI_RECV (force_impose,comm_couplage%m_dim*comm_couplage%m_gl_Ngauss, MPI_DOUBLE_PRECISION,&
            Tdomain%master_superviseur,tag, Tdomain%communicateur_global, status, ierr )


        ! on depouille le vecteur force imposee : en parcourant les faces de couplage et en recuperant, selon la position de chacune, sa force
        position=1
        do iface=1,comm_couplage%m_nbface
            numface  = face_couplage(iface)%face
            numElem  = Tdomain%sFace(numFace)%Near_Element(0)
            ngll     = face_couplage(iface)%ngll

            if ( Tdomain%sface(numFace)%Near_Element(0) == numElem ) then
                wf = Tdomain%sface(numFace)%Which_face(0)
            elseif ( Tdomain%sface(numFace)%Near_Element(1) == numElem ) then
                wf = Tdomain%sface(numFace)%Which_face(1)
            endif

            mat = Tdomain%specel(numElem)%mat_index

            n0=face_couplage(iface)%noeud0
            n1=face_couplage(iface)%noeud1

            ipoint = Tdomain%specel(numElem)%Control_Nodes(0)
            x(0)= Tdomain%Coord_nodes(0,ipoint); z(0)= Tdomain%Coord_nodes(1,ipoint)
            ipoint = Tdomain%specel(numElem)%Control_Nodes(1)
            x(1)= Tdomain%Coord_nodes(0,ipoint); z(1)= Tdomain%Coord_nodes(1,ipoint)
            ipoint = Tdomain%specel(numElem)%Control_Nodes(2)
            x(2)= Tdomain%Coord_nodes(0,ipoint); z(2)= Tdomain%Coord_nodes(1,ipoint)
            ipoint = Tdomain%specel(numElem)%Control_Nodes(3)
            x(3)= Tdomain%Coord_nodes(0,ipoint); z(3)= Tdomain%Coord_nodes(1,ipoint)



            if (wf == 0) then
                dist = sqrt((x(0)-x(1))**2+(z(0)-z(1))**2)
            else if  (wf == 1) then
                dist = sqrt((x(1)-x(2))**2+(z(1)-z(2))**2)
            else if  (wf == 2) then
                dist = sqrt((x(3)-x(2))**2+(z(3)-z(2))**2)
            else if  (wf == 3) then
                dist = sqrt((x(0)-x(3))**2+(z(0)-z(3))**2)
            endif

            do i=1,ngll-2
                ! calcul des poids au point de gauss
                if (wf == 0) then
                    out = Tdomain%sSubdomain(mat)%GLLwx(i)
                else if  (wf == 1) then
                    out = Tdomain%sSubdomain(mat)%GLLwz(i)
                else if  (wf == 2) then
                    out = Tdomain%sSubdomain(mat)%GLLwx(i)
                else if  (wf == 3) then
                    out = Tdomain%sSubdomain(mat)%GLLwz(i)
                endif

                if ( wf == 2 .or. wf == 3 ) then
                    j = ngll -i - 1
                else
                    j = i
                endif
                forcex=force_impose(1,position+j)
                forcey=force_impose(2,position+j)

                Tdomain%sFace(numface)%ForcesMka(i,0) = dist*(out/2.)*forcex
                Tdomain%sFace(numface)%ForcesMka(i,1) = dist*(out/2.)*forcey
            enddo



            ! calcul du poids de gauss sur le vertex 0 (il est identique sur le vertex 1)
            ! test mariotti avec surface correspondat a la surface de la cellule de voronoi associee
            if (wf == 0 .or. wf==2) then
                out = Tdomain%sSubdomain(mat)%GLLwx(0)
            else
                out = Tdomain%sSubdomain(mat)%GLLwz(0)
            endif

            forcex=force_impose(1,position)
            forcey=force_impose(2,position)

            Tdomain%sVertex(n0)%ForcesMka(0) = Tdomain%sVertex(n0)%ForcesMka(0)+dist*(out/2.)*forcex
            Tdomain%sVertex(n0)%ForcesMka(1) = Tdomain%sVertex(n0)%ForcesMka(1)+dist*(out/2.)*forcey

            forcex=force_impose(1,position+ngll-1)
            forcey=force_impose(2,position+ngll-1)

            Tdomain%sVertex(n1)%ForcesMka(0) = Tdomain%sVertex(n1)%ForcesMka(0)+dist*(out/2.)*forcex
            Tdomain%sVertex(n1)%ForcesMka(1) = Tdomain%sVertex(n1)%ForcesMka(1)+dist*(out/2.)*forcey

            position = position+ngll
        enddo


    end subroutine calcul_couplage_force



    !>
    !! \brief Calcul des abscisses curvilignes du projete d'un point de couplage
    !! Valeurs definies dans [-1,1]
    !!
    !! \param real xpdc Abscisse du projete du point de couplage
    !! \param real zpdc Ordonnee du projete du point de couplage
    !! \param real x0 Abscisse du 1er sommet de la face
    !! \param real z0 Ordonnee du 1er sommet de la face
    !! \param real x1 Abscisse du 2eme sommet de la face
    !! \param real z1 Ordonnee du 2eme sommet de la face
    !!
    !<

    real function calcul_xi_eta(xpdc,zpdc,x0,z0,x1,z1)


        real :: xpdc,zpdc
        real :: x0,x1,z0,z1
        real :: prodscal,norm2

        !calcul_xi_eta_old = (2*xpdc - (x1+x0))/(x1-x0)
        prodscal = (xpdc-x0)*(x1-x0)+(zpdc-z0)*(z1-z0)
        norm2     = (x1-x0)*(x1-x0) + (z1-z0)*(z1-z0)

        calcul_xi_eta = -1+2*prodscal/norm2


    end function calcul_xi_eta


    subroutine calcule_coord_glob(Tdomain, n, Xpdc, Zpdc)

        type (domain), intent(IN)  :: Tdomain
        integer, intent(in) :: n
        integer :: np
        real Xpdc(n), Zpdc(n)
        integer :: numFace,numElem

        do np=1,n
            numFace = comm_couplage%m_numFace(np)
            numElem = comm_couplage%m_numMaille(np)
            call calcule_coord(Tdomain, numFace, numElem, comm_couplage%m_pos_proj(np), Xpdc(np), Zpdc(np))
        enddo


    end subroutine calcule_coord_glob

    !>
    !! \brief Remplissage de la structure face_couplage
    !! - Definition du numero global des faces Sem de couplage, des numeros des noeuds les formant
    !! - des positions de ces noeuds, du numero de proc Sem de la face
    !!
    !<


    subroutine remplit_face_couplage(Tdomain, MaxNgParFace)

        type (domain), intent(IN)  :: Tdomain
        integer, intent(out) :: MaxNgParFace
        integer :: numFace, ngll
        integer :: iface

        ! remplissage de la structure face_couplage
        do iface=1,Nbface
            numFace=comm_couplage%m_tabFaceCouplage(iface)
            ngll = Tdomain%sFace(numFace)%ngll

            MaxNgParFace = max(MaxNgParFace,ngll)
            !    Ngauss = Ngauss + ngll !!initialement Gsa

            face_couplage(iface)%face   = numFace
            face_couplage(iface)%proc   = Tdomain%MPI_var%my_rank
            face_couplage(iface)%ngll   = ngll
            if (Tdomain%sface(numFace)%Which_face(0)<2) then ! orientation des faces dans le sens trigo
                face_couplage(iface)%noeud0 = Tdomain%sface(numFace)%Near_Vertex(0)
                face_couplage(iface)%noeud1 = Tdomain%sface(numFace)%Near_Vertex(1)
            else
                face_couplage(iface)%noeud0 = Tdomain%sface(numFace)%Near_Vertex(1)
                face_couplage(iface)%noeud1 = Tdomain%sface(numFace)%Near_Vertex(0)
            endif
            face_couplage(iface)%coord0(1) = Tdomain%Coord_nodes(0,face_couplage(iface)%noeud0)
            face_couplage(iface)%coord0(2) = Tdomain%Coord_nodes(1,face_couplage(iface)%noeud0)
            face_couplage(iface)%coord1(1) = Tdomain%Coord_nodes(0,face_couplage(iface)%noeud1)
            face_couplage(iface)%coord1(2) = Tdomain%Coord_nodes(1,face_couplage(iface)%noeud1)

        enddo

    end subroutine remplit_face_couplage


    !>
    !! \brief Definition des tableaux concernant les points du peigne
    !! Abscisse curviligne de ces points d'interpolation, matrice de projection, ..
    !!
    !<


    subroutine traitement_Pk(Tdomain, iface, numElem, numFace, mat)

        type (domain), intent(IN)  :: Tdomain
        integer, intent(in) :: iface, numElem, numFace, mat
        integer ik, nb_Pk
        integer ngll
        integer i, ip, wf
        real :: Lface, outx, outz

        ngll = Tdomain%sFace(numFace)%ngll

        if ( Tdomain%sface(numFace)%Near_Element(0) == numElem ) then
            wf = Tdomain%sface(numFace)%Which_face(0)
        elseif ( Tdomain%sface(numFace)%Near_Element(1) == numElem ) then
            wf = Tdomain%sface(numFace)%Which_face(1)
        endif

        ! calcul des coordonnees (xi,eta) du point de couplage

        nb_Pk = tab_Pk(iface)%nb_pts

        allocate(tab_Pk(iface)%Interp_coeff(nb_Pk, 0:ngll-1))
        allocate(tab_Pk(iface)%part_mka_plus_proche(nb_Pk))
        allocate(tab_Pk(iface)%xi(nb_Pk))
        allocate(tab_Pk(iface)%eta(nb_Pk))
        allocate(tab_Pk(iface)%pdg_plus_proche(nb_Pk))

        Lface = sqrt( (face_couplage(iface)%coord0(1) - face_couplage(iface)%coord1(1) )**2 + &
            ( face_couplage(iface)%coord0(2) - face_couplage(iface)%coord1(2) ) **2)
        !! surface liee a un point du peigne sur la face iface
        !! elle est identique pour tous les points du peigne de la face iface
        tab_Pk(iface)%surf = Lface/(nb_Pk + 1.)

        do ik=1, nb_Pk
            if (wf == 0) then
                tab_Pk(iface)%eta(ik) =-1
                tab_Pk(iface)%xi(ik) = -1. + 2.*ik/(nb_Pk + 1.)
            else if  (wf == 1) then
                tab_Pk(iface)%xi(ik) =1
                tab_Pk(iface)%eta(ik) = -1. + 2.*ik/(nb_Pk + 1.)
            else if  (wf == 2) then
                tab_Pk(iface)%eta(ik)=1
                tab_Pk(iface)%xi(ik) =  1. - 2.*ik/(nb_Pk + 1.)
            else if  (wf == 3) then
                tab_Pk(iface)%xi(ik)=-1
                tab_Pk(iface)%eta(ik) =  1. - 2.*ik/(nb_Pk + 1.)
            endif

            call proj_plus_proche(Tdomain, comm_couplage%m_nb_point_de_couplage, iface, tab_Pk(iface)%xi(ik), &
                tab_Pk(iface)%eta(ik), Xpdc, Zpdc, ip)
            tab_Pk(iface)%part_mka_plus_proche(ik) = ip

            if (wf == 1.or. wf ==3) then
                do i=0,ngll-1
                    call  pol_lagrange (ngll,Tdomain%sSubdomain(mat)%GLLcz,i,tab_Pk(iface)%eta(ik),outz)
                    tab_Pk(iface)%Interp_coeff(ik,i) = outz
                enddo
            else if (wf == 0.or. wf ==2) then
                do i=0,ngll-1
                    call  pol_lagrange (ngll,Tdomain%sSubdomain(mat)%GLLcx,i,tab_Pk(iface)%xi(ik),outx)
                    tab_Pk(iface)%Interp_coeff(ik,i) = outx
                enddo
            endif

        enddo

    end subroutine traitement_Pk


    !>
    !! \brief Calcul de la plus petite distance entre particules de couplage Mka pour une meme face
    !! de couplage Sem. On teste aussi par rapport aux sommets de la face Sem.
    !<

    subroutine dist_min_point_de_couplage(n, Xpdc, Zpdc, dmin)

        integer, intent(in) :: n
        real, intent(in) :: Xpdc(n), Zpdc(n)
        integer :: iface, np1, np2
        real :: dist, dmin(Nbface)

        if(Nbface .EQ. 0) return
        dmin(:) = huge(1.)

        do iface=1,Nbface
            do np1=1,n
                if(comm_couplage%m_numFace(np1) == face_couplage(iface)%face ) then
                    do np2=np1+1,n
                        if(comm_couplage%m_numFace(np2) == comm_couplage%m_numFace(np1) ) then
                            dist = (Xpdc(np1) - Xpdc(np2))**2 + (Zpdc(np1) - Zpdc(np2))**2
                            dist = sqrt(dist)
                            dmin(iface) = min(dmin(iface), dist)
                        endif
                    enddo
                    ! on utilise les sommets de la face (pour bien discretiser les maillages fins)
                    ! meme si ca ajoute des points d'interpolation
                    dist = (Xpdc(np1) - face_couplage(iface)%coord0(1))**2 + (Zpdc(np1) - face_couplage(iface)%coord0(2))**2
                    dist = sqrt(dist)
                    dmin(iface) = min(dmin(iface), dist)
                    dist = (Xpdc(np1) - face_couplage(iface)%coord1(1))**2 + (Zpdc(np1) - face_couplage(iface)%coord1(2))**2
                    dist = sqrt(dist)
                    dmin(iface) = min(dmin(iface), dist)
                endif

            enddo
        enddo
    end subroutine dist_min_point_de_couplage


    !>
    !! \brief Retourne  l'indice de la particule de couplage Mka se projetant sur la face iface
    !! la plus proche du point d'interpolation d'abscisses (xi,eta)
    !<


    !  integer function proj_plus_proche(Tdomain, n, iface, ik, Xpdc, Zpdc)
    subroutine proj_plus_proche(Tdomain, n, iface, xi, eta, Xpdc, Zpdc, ip)
        type (domain), intent(IN)  :: Tdomain
        integer, intent(in) :: iface, n
        real, intent(in) :: Xpdc(n), Zpdc(n)
        real, intent(IN) :: xi, eta
        integer ip, np
        real :: dmin, dist
        integer :: numFace, numElem, ngll, wf
        real,dimension(0:3) :: x,z
        real,dimension(0:1) :: pos_Pk
        integer ipoint, inod

        dmin = huge(1.)

        do np=1,n
            ! recuperation de la face et de l element du point de couplage courant et de ses coordonnees, du numero de materiau et du nombre de pdg en x et z
            numFace = comm_couplage%m_numFace(np)

            if(face_couplage(iface)%face == numFace) then

                numElem = comm_couplage%m_numMaille(np)

                ! recuperation des coor des 4 sommets de l'element
                do inod=0,3
                    ipoint = Tdomain%specel(numElem)%Control_Nodes(inod)
                    x(inod) = Tdomain%Coord_nodes(0,ipoint)
                    z(inod) = Tdomain%Coord_nodes(1,ipoint)
                enddo

                !ngllx = Tdomain%specel(numElem)%ngllx;  ngllz = Tdomain%specel(numElem)%ngllz

                ngll = Tdomain%sFace(numFace)%ngll

                if ( Tdomain%sface(numFace)%Near_Element(0) == numElem ) then
                    wf = Tdomain%sface(numFace)%Which_face(0)
                elseif ( Tdomain%sface(numFace)%Near_Element(1) == numElem ) then
                    wf = Tdomain%sface(numFace)%Which_face(1)
                endif

                if((wf==0) .OR. (wf==1)) then
                    pos_Pk(0) = 1./4.*( (1.-xi)*(1.-eta)*x(0) + (1.+xi)*(1.-eta)*x(1) &
                        + (1.+xi)*(1.+eta)*x(2) + (1.-xi)*(1.+eta)*x(3) )

                    pos_Pk(1) = 1./4.*( (1.-xi)*(1.-eta)*z(0) + (1.+xi)*(1.-eta)*z(1) &
                        + (1.+xi)*(1.+eta)*z(2) + (1.-xi)*(1.+eta)*z(3) )

                elseif (wf==2) then
                    pos_Pk(0) = 1./4.*( (1.-xi)*(1.-eta)*x(1) + (1.+xi)*(1.-eta)*x(0) &
                        + (1.+xi)*(1.+eta)*x(3) + (1.-xi)*(1.+eta)*x(2) )

                    pos_Pk(1) = 1./4.*( (1.-xi)*(1.-eta)*z(1) + (1.+xi)*(1.-eta)*z(0) &
                        + (1.+xi)*(1.+eta)*z(3) + (1.-xi)*(1.+eta)*z(2) )

                elseif (wf==3) then
                    pos_Pk(0) = 1./4.*( (1.-xi)*(1.-eta)*x(3) + (1.+xi)*(1.-eta)*x(2) &
                        + (1.+xi)*(1.+eta)*x(1) + (1.-xi)*(1.+eta)*x(0) )

                    pos_Pk(1) = 1./4.*( (1.-xi)*(1.-eta)*z(3) + (1.+xi)*(1.-eta)*z(2) &
                        + (1.+xi)*(1.+eta)*z(1) + (1.-xi)*(1.+eta)*z(0) )

                endif


                dist = sqrt( (pos_Pk(0)- Xpdc(np))**2 + (pos_Pk(1)- Zpdc(np))**2)
                if(dist< dmin) then
                    ip = np
                    dmin = dist
                endif
            endif
        enddo

        write(75,'(6(1pe15.8,1X),1X,I6)') pos_Pk, Xpdc(ip), Zpdc(ip), xi, eta, ip

    end subroutine proj_plus_proche

    !>
    !! \brief Retourne l'indice du point d'interpolation le plus proche d'une particule
    !! de couplage Mka (d'indice ip)
    !<

    integer function trouve_Pk_plus_proche(Tdomain, n, ip, Xpdc, Zpdc)
        type (domain), intent(IN)  :: Tdomain
        integer, intent(in) ::  n
        real, intent(in) :: Xpdc(n), Zpdc(n)
        integer ip
        real :: dmin, dist
        integer :: numFace, numElem, wf
        real,dimension(0:3) :: x,z
        real,dimension(0:1) :: pos_Pk, pos_trouve
        integer iface, ik, ipoint, itrouve, nb_Pk, inod
        real :: xi, eta
        !! trouver la position du projete - on l'a Xpdc
        !! utiliser wf
        !! on choisit le plus proche
        !! -1. + 2.*ik/(tab_Pk(iface)%nb_Pk + 2.)


        dmin = huge(1.)
        numFace = comm_couplage%m_numFace(ip)

        do iface=1,Nbface
            if(face_couplage(iface)%face == numFace) exit
        enddo
        nb_Pk = tab_Pk(iface)%nb_pts

        numElem = comm_couplage%m_numMaille(ip)
        ! recuperation des coor des 4 sommets de l'element
        do inod=0,3
            ipoint = Tdomain%specel(numElem)%Control_Nodes(inod)
            x(inod) = Tdomain%Coord_nodes(0,ipoint)
            z(inod) = Tdomain%Coord_nodes(1,ipoint)
        enddo

        if ( Tdomain%sface(numFace)%Near_Element(0) == numElem ) then
            wf = Tdomain%sface(numFace)%Which_face(0)
        elseif ( Tdomain%sface(numFace)%Near_Element(1) == numElem ) then
            wf = Tdomain%sface(numFace)%Which_face(1)
        endif

        do ik=1,nb_Pk
            ! recuperation de la face et de l element du point de
            ! couplage courant et de ses coordonnees, du numero de materiau et du
            ! nombre de pdg en x et z

            xi = tab_Pk(iface)%xi(ik)
            eta = tab_Pk(iface)%eta(ik)

            if((wf==0) .OR. (wf==1)) then
                pos_Pk(0) = 1./4.*( (1.-xi)*(1.-eta)*x(0) + (1.+xi)*(1.-eta)*x(1) &
                    + (1.+xi)*(1.+eta)*x(2) + (1.-xi)*(1.+eta)*x(3) )

                pos_Pk(1) = 1./4.*( (1.-xi)*(1.-eta)*z(0) + (1.+xi)*(1.-eta)*z(1) &
                    + (1.+xi)*(1.+eta)*z(2) + (1.-xi)*(1.+eta)*z(3) )

            elseif (wf==2) then
                pos_Pk(0) = 1./4.*( (1.-xi)*(1.-eta)*x(1) + (1.+xi)*(1.-eta)*x(0) &
                    + (1.+xi)*(1.+eta)*x(3) + (1.-xi)*(1.+eta)*x(2) )

                pos_Pk(1) = 1./4.*( (1.-xi)*(1.-eta)*z(1) + (1.+xi)*(1.-eta)*z(0) &
                    + (1.+xi)*(1.+eta)*z(3) + (1.-xi)*(1.+eta)*z(2) )

            elseif (wf==3) then
                pos_Pk(0) = 1./4.*( (1.-xi)*(1.-eta)*x(3) + (1.+xi)*(1.-eta)*x(2) &
                    + (1.+xi)*(1.+eta)*x(1) + (1.-xi)*(1.+eta)*x(0) )

                pos_Pk(1) = 1./4.*( (1.-xi)*(1.-eta)*z(3) + (1.+xi)*(1.-eta)*z(2) &
                    + (1.+xi)*(1.+eta)*z(1) + (1.-xi)*(1.+eta)*z(0) )

            endif
            dist = sqrt( (pos_Pk(0)- Xpdc(ip))**2 + (pos_Pk(1)- Zpdc(ip))**2)
            if(dist< dmin) then
                itrouve = ik
                dmin = dist
                pos_trouve = pos_Pk
            endif
        enddo

        trouve_Pk_plus_proche = itrouve
        write(77,'(4(1pe15.8,1X),1X,I6,1X,I6)') Xpdc(ip), Zpdc(ip), pos_trouve(0:1), ip, itrouve

    end function trouve_Pk_plus_proche


    !>
    !! \brief Definition du nombre de points d'interpolation pour une face de couplage Sem
    !!
    !! Ce nombre depend du nombre de points de Gauss et de la plus petite distance
    !! entre les particules de couplage.
    !<

    subroutine definit_nb_pts_Pk(Tdomain, dmin_couplage)

        type (domain), intent(IN)  :: Tdomain
        real, intent(in) :: dmin_couplage(Nbface)
        integer ngll, numFace
        integer iface
        real Lface
        real hk

        do iface=1,Nbface
            numFace = comm_couplage%m_tabFaceCouplage(iface)
            ngll = Tdomain%sFace(numFace)%ngll
            Lface = sqrt( (face_couplage(iface)%coord0(1) - face_couplage(iface)%coord1(1) )**2 + &
                ( face_couplage(iface)%coord0(2) - face_couplage(iface)%coord1(2) ) **2)
            !       hk = 1./3. *min(dmin_couplage(iface), Lface/(2.*ngll + 1.)) !! pour avoir plus de points d'interpolation
            hk = min(dmin_couplage(iface), Lface/(2.*ngll + 1.))
            !       hk = min(dmin_couplage(iface), Lface/(2.*ngll + 25.))
            !       tab_Pk(iface)%nb_pts = floor( Lface/hk) + 1
            tab_Pk(iface)%nb_pts = floor( Lface/hk)
            hk = Lface/tab_Pk(iface)%nb_pts !correction
        enddo
    end subroutine definit_nb_pts_Pk

    !---------------------------------------------
    !---------------------------------------------

    !>
    !! \brief Reception des surfaces des particules de couplage Mka.
    !!
    !<


    subroutine reception_surface_part_mka(Tdomain)

        implicit none

        type (domain), intent(IN)  :: Tdomain
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer :: tag, ierr

        ! a faire une seule fois
        if(comm_couplage%m_nb_point_de_couplage > 0) then
            tag = 510000+1000*Tdomain%MPI_var%my_rank+1
            call MPI_RECV (comm_couplage%m_surf_part,comm_couplage%m_nb_point_de_couplage, &
                MPI_DOUBLE_PRECISION, Tdomain%master_superviseur, &
                tag, Tdomain%communicateur_global, status, ierr )
        endif

        deallocate(Xpdc)
        deallocate(Zpdc)

    end subroutine reception_surface_part_mka


    !---------------------------------------------
!!!!  Description: Definition du second membre pour la vitesse dans le systeme lineaire
!!!!
!!!! Simple recuperation des vitesses aux points de Gauss
!!!! Envoi des donnees au Superviseur
!!!! --------------------------------------------------

    subroutine envoi_vitesse_mka(Tdomain,ntime)

        type (domain), intent(INOUT)  :: Tdomain
        integer, intent(IN) :: ntime

        integer :: tag, ierr
        integer :: i,numFace
        integer :: ngll
        integer :: numVertex0, numVertex1
        integer iface
        integer idim
        integer :: numElem, wf
        real :: tsurf
        real,dimension(:,:),pointer :: vecu, vecu_tmp

        allocate(vecu(comm_couplage%m_MaxNgParFace,comm_couplage%m_dim))
        allocate(vecu_tmp(comm_couplage%m_MaxNgParFace,comm_couplage%m_dim))

        vecu = 0.
        do iface=1,comm_couplage%m_NbFace

            numFace=comm_couplage%m_tabFaceCouplage(iface)
            ngll = Tdomain%sFace(numFace)%ngll
            numElem  = Tdomain%sFace(numFace)%Near_Element(0)
            wf = Tdomain%sface(numFace)%Which_face(0)

            numVertex0 = Tdomain%sface(numFace)%Near_Vertex(0)
            numVertex1 = Tdomain%sface(numFace)%Near_Vertex(1)


            do idim=1,comm_couplage%m_dim-1
                vecu(1:ngll,idim)=0.
                if (wf<2) then
                    vecu(1,idim) = Tdomain%sVertex(numVertex0)%Veloc(idim-1)
                    do i=1,ngll-2
                        vecu(i+1,idim) = Tdomain%sFace(numFace)%Veloc(i,idim-1)
                    enddo
                    vecu(ngll,idim) = Tdomain%sVertex(numVertex1)%Veloc(idim-1)

                elseif (wf>1) then
                    vecu(1,idim) = Tdomain%sVertex(numVertex1)%Veloc(idim-1)
                    do i=1,ngll-2
                        vecu(i+1,idim) = Tdomain%sFace(numFace)%Veloc(ngll-1-i,idim-1)
                    enddo
                    vecu(ngll,idim) = Tdomain%sVertex(numVertex0)%Veloc(idim-1)
                endif
            enddo

            !! si CAS  M^-1 (D^t A D) u
            vecu_tmp = vecu
            do i=0,ngll-1
                tsurf = surface_gll(Tdomain, numElem, wf, i)
                vecu_tmp(i+1,1:3) = vecu_tmp(i+1,1:3)*tsurf
            enddo

            do idim=1,comm_couplage%m_dim

                tag = 730000+1000*Tdomain%MPI_var%my_rank+10*ntime+(iface-1)+100*(idim-1)
                call MPI_SEND(vecu_tmp(1,idim), ngll, MPI_DOUBLE_PRECISION, Tdomain%master_superviseur,&
                    tag, Tdomain%communicateur_global, ierr )

            enddo
        enddo

        deallocate(vecu)
        deallocate(vecu_tmp)

    end subroutine envoi_vitesse_mka


    function surface_gll(Tdomain, numElem, wf,  i)
        use sdomain
        implicit none
        type (domain), intent(IN)  :: Tdomain
        integer, intent(IN) ::  numElem
        integer, intent(IN) :: i, wf
        real out, dist
        real surface_gll
        real,dimension(0:3) :: x,z
        integer ipoint, inod, mat

        mat = Tdomain%specel(numElem)%mat_index
        ! recuperation des coor des 4 sommets de l'element
        do inod=0,3
            ipoint = Tdomain%specel(numElem)%Control_Nodes(inod)
            x(inod) = Tdomain%Coord_nodes(0,ipoint)
            z(inod) = Tdomain%Coord_nodes(1,ipoint)
        enddo

        ! calcul des poids au point de gauss
        select case (wf)
        case(0)
            out = Tdomain%sSubdomain(mat)%GLLwx(i)
            dist = sqrt((x(0)-x(1))**2+(z(0)-z(1))**2)
        case(1)
            out = Tdomain%sSubdomain(mat)%GLLwz(i)
            dist = sqrt((x(1)-x(2))**2+(z(1)-z(2))**2)
        case(2)
            out = Tdomain%sSubdomain(mat)%GLLwx(i)
            dist = sqrt((x(3)-x(2))**2+(z(3)-z(2))**2)
        case(3)
            out = Tdomain%sSubdomain(mat)%GLLwz(i)
            dist = sqrt((x(0)-x(3))**2+(z(0)-z(3))**2)
        end select

        surface_gll = dist*out/2.


    end function surface_gll

    subroutine reception_nouveau_pdt_sem(Tdomain)

        use sdomain
        implicit none

        type (Domain), intent (INOUT) :: Tdomain
        integer :: n, mat
        real dtsem
        integer :: tag,ierr
        integer, dimension (MPI_STATUS_SIZE) :: status

        tag = 4100000+10000*Tdomain%MPI_var%my_rank
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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
