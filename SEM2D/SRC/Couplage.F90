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
       real,dimension(:), pointer    :: m_surf_part       ! surface des particules Mka de couplage

    end type type_comm_couplage


    type :: type_face_couplage

       integer :: face   ! numero de la face de couplage sem dans la numerotation du proc sem qui la contient
       integer :: proc   ! numero du processeur sem qui contient cette face de couplage
       integer :: ngll   ! nbre de point de gauss sur la face
       integer :: noeud0 ! numero du premier noeud de la face
       integer :: noeud1 ! numero du 2ieme noeud de la face
       real,dimension(2) :: coord0 ! coordonnees noeud0
       real,dimension(2) :: coord1 ! coordonnees noeud1
       integer :: nbPk   ! nombre de points du peigne

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
        integer :: i
        integer :: tag,ierr
        integer :: wf,ngll,np
        integer :: numFace,numElem
        integer :: iface,bufsize,decal,nbchamps
        integer :: Ngauss,MaxNgParFace,gl_MaxNgParFace
        logical :: dejaPresent
        integer, dimension (MPI_STATUS_SIZE) :: status
        integer, dimension(:),pointer:: tabNbFace
        integer,dimension(:),pointer :: displs,count
        integer,dimension(5)  :: ibuf
        real, dimension (:),pointer :: buf
        real,dimension (:), pointer :: dmin_couplage
        real,dimension(4) :: coord


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
        ! calcul des fonctions de forme en chaque point de couplage
        do np=1,comm_couplage%m_nb_point_de_couplage
            ! recuperation de la face et de l element du point de couplage courant et de ses coordonnees, du numero de materiau et du nombre de pdg en x et z
            numFace = comm_couplage%m_numFace(np)
            numElem = comm_couplage%m_numMaille(np)

            call calcule_coord(Tdomain, numFace, numElem, comm_couplage%m_pos_proj(np), Xpdc(np), Zpdc(np))

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


        allocate(tabNbFace(Tdomain%Mpi_var%n_proc))
        tabNbFace(:)=0

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
                        face_couplage(iface+decal)%face   = int(buf(1+NbChamps*(iface-1)))
                        face_couplage(iface+decal)%proc   = int(buf(2+NbChamps*(iface-1)))
                        face_couplage(iface+decal)%ngll   = int(buf(3+NbChamps*(iface-1)))
                        face_couplage(iface+decal)%noeud0 = int(buf(4+NbChamps*(iface-1)))
                        face_couplage(iface+decal)%noeud1 = int(buf(5+NbChamps*(iface-1)))
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


        ! construction de la matrice elementaire de chaque face de couplage et envoi au superviseur, pour chaque proc sem
        do iface=1,Nbface

            numFace=comm_couplage%m_tabFaceCouplage(iface)
            ngll = Tdomain%sFace(numFace)%ngll
            numElem  = Tdomain%sFace(numFace)%Near_Element(0)

            if ( Tdomain%sface(numFace)%Near_Element(0) == numElem ) then
                wf = Tdomain%sface(numFace)%Which_face(0)
            elseif ( Tdomain%sface(numFace)%Near_Element(1) == numElem ) then
                wf = Tdomain%sface(numFace)%Which_face(1)
            endif

            ibuf(1) = ngll ! ngll1
            ibuf(2) = 0    ! ngll2
            ibuf(3) = face_couplage(iface)%nbPk
            ibuf(4) = wf
            ibuf(5) = numFace

            tag=670000+Tdomain%MPI_var%my_rank+iface-1
            call MPI_SEND(ibuf, 5, MPI_INTEGER, Tdomain%master_superviseur,tag, Tdomain%communicateur_global, ierr )

            coord(1) = face_couplage(iface)%coord0(1)
            coord(2) = face_couplage(iface)%coord1(1)
            coord(3) = face_couplage(iface)%coord0(2)
            coord(4) = face_couplage(iface)%coord1(2)

            ! On envoie les coordonnees des sommets de la face de couplage
            tag=690000+Tdomain%MPI_var%my_rank+iface-1
            call MPI_SEND(coord, 4, MPI_DOUBLE_PRECISION, Tdomain%master_superviseur, &
                tag, Tdomain%communicateur_global, ierr )

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
            else
                stop "Internal Error, incorrect mesh information"
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
            else
                stop "Internal Error, incorrect mesh information"
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
                else
                    stop "Internal Error, incorrect mesh information"
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
            face_couplage(iface)%nbPk = int(floor( Lface/hk)*10.)
            !hk = Lface/face_couplage(iface)%nbPk !correction
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

                tag = 730000+3*(iface-1)+(idim-1)
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
        dist = 1.0
        out = 1.0 ! Suppression warning
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
