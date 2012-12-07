!>
!!\file sem.f90
!!\brief Contient le module ipg2.
!!\author
!!\version 1.0
!!\date 10/03/2009
!! Gère les noeuds, les éléments du maillage à la frontière.
!<

module ipg2

    use semdatafiles

    type :: Tsem

       integer :: fileId       ! Id du fichier sem a generer
       integer :: n_dime       ! dimension de l'espace , en principe on est en 2D
       integer :: n_glob_nodes ! nombre de noeuds du maillage
       integer :: n_elem       ! nombre d'elements du maillage
       integer :: n_nodes      ! nombre de noeuds par element
       integer :: n_vertex     ! nombre de noeuds par element
       integer :: n_edge       ! nombre d'aretes par element
       integer :: n_face       ! nombre de face du maillage
       integer :: n_mat        ! nombre de materiaux
       integer :: typeElem     ! typeElem=8 pour quad8, 20 ou 27

       ! mise en comment a cause de gfortran (allocatable ->pointer)
!!!  integer, dimension (:), allocatable   :: glob_numbering ! numeros globaux des noeuds
!!!
!!!  integer, dimension (:,:), allocatable :: Near_Faces      ! numeros globaux des faces de l'element
!!!  integer, dimension (:,:), allocatable :: Near_Edges      ! numeros globaux des aretes de l'element
!!!  integer, dimension (:,:), allocatable :: Near_Vertex    ! numeros globaux des noeuds de l'element
!!!  integer, dimension (:,:), allocatable :: Control_nodes  ! numeros des noeuds (numerotation sequentielle des points de coordonnees X Y)
!!!  integer, dimension (:,:), allocatable :: Orient_Face      ! numeros d'orientation des faces de l'element
!!!  integer, dimension (:,:), allocatable :: Orient_Edges      ! numeros d'orientation des aretes de l'element
!!!  integer, dimension (:,:), allocatable :: face_element   ! numeros du ou des elements de la face
!!!  integer, dimension (:,:), allocatable :: face_vertex    ! numero global des sommets de la face
!!!  integer, dimension (:,:), allocatable :: edge_vertex    ! numero global des sommets de l'arete
!!!  real   , dimension (:,:), allocatable :: Coord_nodes ! coordonnees X,Y,Z des noeuds


       integer, dimension (:), pointer   :: glob_numbering ! numeros globaux des noeuds

       integer, dimension (:,:), pointer :: Near_Faces      ! numeros globaux des faces de l'element
       integer, dimension (:,:), pointer :: Near_Edges      ! numeros globaux des aretes de l'element
       integer, dimension (:,:), pointer :: Near_Vertex    ! numeros globaux des noeuds de l'element
       integer, dimension (:,:), pointer :: Control_nodes  ! numeros des noeuds (numerotation sequentielle des points de coordonnees X Y)
       integer, dimension (:,:), pointer :: Orient_Face      ! numeros d'orientation des faces de l'element
       integer, dimension (:,:), pointer :: Orient_Edges      ! numeros d'orientation des aretes de l'element
       integer, dimension (:,:), pointer :: face_element   ! numeros du ou des elements de la face
       integer, dimension (:,:), pointer :: face_vertex    ! numero global des sommets de la face
       integer, dimension (:,:), pointer :: edge_vertex    ! numero global des sommets de l'arete
       real   , dimension (:,:), pointer :: Coord_nodes ! coordonnees X,Y,Z des noeuds

       character (len=MAX_FILE_SIZE) :: fichier

    end type Tsem
    character (len=MAX_FILE_SIZE) :: fnamef

    type TelementDeFace
       integer :: numElem
       integer :: numLocalFace
       !integer, dimension(0:3) :: face ! tableau des numeros globaux de faces
       type(TelementDeFace), pointer :: suivant
    end type TelementDeFace

    type TelementDArete ! NOUVEAU TYPE
       integer :: numElem
       integer :: numLocalEdge
       type(TelementDArete), pointer :: suivant
    end type TelementDArete


    type Tface   ! TYPE MODIFIE
       integer :: globnum ! numero global de la face
       integer :: n1, n2, n3, n4 ! numero globaux des noeuds de la face
       type(TelementDeFace),pointer :: listeElement
       type(Tface), pointer :: suivant
    end type Tface

    type Tedge     ! NOUVEAU TYPE
       integer :: globnum ! numero global de l'arete
       integer :: n1, n2 ! numero globaux des noeuds de l'arete
       type(TelementDArete), pointer :: listeElement
       type(Tedge), pointer :: suivant
    end type Tedge


    type Tedgefrontiere ! NOUVEAU TYPE
       integer :: iedge ! numero de la face sur la frontiere
       integer :: iVertex1, iVertex2 ! numeros des sommets de iface
       integer :: iCtrlNode1, iCtrlNode2 ! numeros des pts de controle corres. a ivertex1,2
       real    :: xVertex1, yVertex1, zVertex1, xVertex2, yVertex2, zVertex2 ! coordonnees X,Y des vertex
       type(Tedgefrontiere), pointer :: suivant
    end type Tedgefrontiere

    type Tfacefrontiere ! NOUVEAU TYPE
       integer :: iface ! numero de la face sur la frontiere
       integer :: iVertex1, iVertex2, iVertex3, iVertex4 ! numeros des sommets de iface
       integer :: iCtrlNode1, iCtrlNode2, iCtrlNode3, iCtrlNode4 ! numeros des pts de controle corres. a ivertex1,2,3,4
       real    :: Vertex(3,4) ! 4 sommets a 3 composantes
       type(Tfacefrontiere),pointer :: suivant
    end type Tfacefrontiere

    type(Tsem) :: sem

    type(Tface), pointer :: listeFace
    type(TEdge), pointer :: listeEdge

    type(Tfacefrontiere), pointer :: listeFaceFrontiere
    type(Tedgefrontiere), pointer :: listeEdgeFrontiere

contains
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !

    ! MISE EN COMMENTAIRE - on ne gere pas la frontiere encore
!!!subroutine ecrireFichierSemFrontiere
!!!
!!! implicit none
!!!
!!! integer :: i,iface,nFaceFrontiere,nVertexFrontiere,CodeErreur
!!! type(Tfacefrontiere),pointer :: nouvelleFace,chercheFace
!!!
!!! integer, dimension(:,:),allocatable :: listeVertex
!!!
!!!
!!! !liste initialement vide
!!! nullify(listeFrontiere) !listeFrontiere est un pointeur sur un type Tfacefrontiere
!!!
!!! ! nbre de faces initiallement nul
!!! nFaceFrontiere = 0
!!!
!!! ! nbre de vertex initiallement nul
!!! nVertexFrontiere = 0
!!!
!!!
!!! ! boucle sur les faces
!!! do iface=1,sem%n_face
!!!    if ( sem%face_Element(1,iface) .ne. -1 ) cycle ! ce n est pas une face frontiere, on passe
!!!
!!!    allocate(nouvelleFace)
!!!
!!!    nouvelleFace%iface = iface
!!!
!!!    ! recup des vertex associes a iface
!!!    nouvelleFace%ivertex1 = sem%face_vertex(0,iface)
!!!    nouvelleFace%ivertex2 = sem%face_vertex(1,iface)
!!!
!!!    ! recup des points de controle
!!!    nouvelleFace%iCtrlNode1 = sem%glob_numbering(nouvelleFace%ivertex1)-1
!!!    nouvelleFace%iCtrlNode2 = sem%glob_numbering(nouvelleFace%ivertex2)-1
!!!
!!!    ! recup des coordonnees des vertex
!!!    nouvelleFace%xVertex1 = sem%Coord_nodes(1,nouvelleFace%iCtrlNode1)
!!!    nouvelleFace%yVertex1 = sem%Coord_nodes(2,nouvelleFace%iCtrlNode1)
!!!    nouvelleFace%zVertex1 = sem%Coord_nodes(3,nouvelleFace%iCtrlNode1)
!!!    nouvelleFace%xVertex2 = sem%Coord_nodes(1,nouvelleFace%iCtrlNode2)
!!!    nouvelleFace%yVertex2 = sem%Coord_nodes(2,nouvelleFace%iCtrlNode2)
!!!    nouvelleFace%zVertex2 = sem%Coord_nodes(3,nouvelleFace%iCtrlNode2)
!!!    nouvelleFace%xVertex3 = sem%Coord_nodes(1,nouvelleFace%iCtrlNode3)
!!!    nouvelleFace%yVertex3 = sem%Coord_nodes(2,nouvelleFace%iCtrlNode3)
!!!    nouvelleFace%zVertex3 = sem%Coord_nodes(3,nouvelleFace%iCtrlNode3)
!!!    nouvelleFace%xVertex4 = sem%Coord_nodes(1,nouvelleFace%iCtrlNode4)
!!!    nouvelleFace%yVertex4 = sem%Coord_nodes(2,nouvelleFace%iCtrlNode4)
!!!    nouvelleFace%zVertex4 = sem%Coord_nodes(3,nouvelleFace%iCtrlNode4)
!!!
!!!
!!!    ! compteur de faces
!!!    nFaceFrontiere = nFaceFrontiere + 1
!!!
!!!    ! on empile la nouvelle face en tete de liste
!!!    nouvelleFace%suivant=>listeFrontiere
!!!
!!!    ! la liste demarre sur cette nouvelle face
!!!    listeFrontiere=>nouvelleFace
!!!
!!! enddo
!!!
!!! ! ouverture du fichier frontiere
!!! open(UNIT=sem%fileId,IOSTAT=CodeErreur,FILE=trim(sem%fichier)//"frontiere",FORM='formatted',STATUS='replace',ACTION='write')
!!! if (CodeErreur .ne.0 ) print*,'Ouverture du fichier SEM frontiere  :CodeErreur=',CodeErreur
!!!
!!!
!!! ! tableau ds lequel on va stocker la liste exacte (sans doublon) des vertex frontieres
!!! allocate(listeVertex(2*nFaceFrontiere,3))
!!!
!!!
!!! ! ecriture des donnees de faces frontieres
!!!
!!! write(sem%fileId,*) nFaceFrontiere ! nombre de faces frontieres
!!!
!!! ! on parcourt la liste frontiere
!!! chercheFace=>listeFrontiere
!!!
!!! do
!!!    if (.not.associated(chercheFace)) exit ! on est en fin de liste
!!!
!!!    ! ecriture du numero de la face
!!!    write(sem%fileId,*) chercheFace%iface, &
!!!         chercheFace%ivertex1,chercheFace%ivertex2
!!!
!!!    chercheFace=>chercheFace%suivant
!!!
!!!    ! mise a jour de la liste des vertex frontieres
!!!    if (.not.VertexFrontiereDejaCompte(chercheFace%ivertex1,nVertexFrontiere,listeVertex(1:nVertexFrontiere,1))) then
!!!       nVertexFrontiere = nVertexFrontiere +1
!!!       listeVertex(nVertexFrontiere,1) = chercheFace%ivertex1
!!!       listeVertex(nVertexFrontiere,2) = chercheFace%xVertex1
!!!       listeVertex(nVertexFrontiere,3) = chercheFace%yVertex1
!!!    endif
!!!
!!!    if (.not.VertexFrontiereDejaCompte(chercheFace%ivertex2,nVertexFrontiere,listeVertex(1:nVertexFrontiere,1))) then
!!!       nVertexFrontiere = nVertexFrontiere +1
!!!       listeVertex(nVertexFrontiere,1) = chercheFace%ivertex2
!!!       listeVertex(nVertexFrontiere,2) = chercheFace%xVertex2
!!!       listeVertex(nVertexFrontiere,3) = chercheFace%yVertex2
!!!    endif
!!!
!!! enddo
!!!
!!!
!!! write(sem%fileId,*)"------------------"
!!! ! ecriture des donnees de aux faces frontieres
!!! write(sem%fileId,*) nVertexFrontiere ! nombre de faces frontieres
!!!
!!! ! on parcourt la liste des vertex
!!!
!!! do i=1,nVertexFrontiere
!!!
!!!    ! ecriture du numero de la face
!!!    write(sem%fileId,*) listeVertex(i,1), &
!!!         listeVertex(i,2),listeVertex(i,3)
!!!
!!! enddo
!!!
!!!
!!! deallocate(listeVertex)
!!!
!!! close(sem%fileId)
!!!
!!!end subroutine ecrireFichierSemFrontiere

    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !


    !>
    !! \fn function VertexFrontiereDejaCompte(ivertex,nVertexFrontiere,listeVertex)
    !! \brief
    !!
    !! \param integer ivertex,
    !! \param integer nVertexFrontiere
    !! \param integer, dimension(1:nVertexFrontiere) listeVertex
    !<


    logical function VertexFrontiereDejaCompte(ivertex,nVertexFrontiere,listeVertex)

        integer :: ivertex,nVertexFrontiere,i
        integer, dimension(1:nVertexFrontiere) :: listeVertex


        VertexFrontiereDejaCompte=.False.
        do i=1,nVertexFrontiere
            if (listeVertex(i).eq.ivertex) then
                VertexFrontiereDejaCompte=.True.
                exit
            endif
        enddo

    end function VertexFrontiereDejaCompte

    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !


    !>
    !! \fn subroutine ecrireFichierSem
    !! \brief
    !!
    !<


    subroutine ecrireFichierSem

        implicit none

        integer :: i,j,CodeErreur
        logical curve !! Tdomain%curve est un logical - Ici on ecrit dans le fichier de maillage que curve=.false.

        call semname_sem_fichier(sem%fichier,fnamef)
        open(UNIT=sem%fileId,IOSTAT=CodeErreur,FILE=fnamef,FORM='formatted',STATUS='replace',ACTION='write')
        if (CodeErreur .ne.0 ) print*,'Ouverture du fichier SEM  :CodeErreur=',CodeErreur



        ! dimension du pb (2D)
        write(sem%fileId,*) sem%n_dime

        ! nombre de noeuds
        write(sem%fileId,*) sem%n_glob_nodes

        ! presence de courbure
        curve = .false.
        write(sem%fileId,*) curve !curve est .false. par defaut - A modifier

        ! coordonnees X Y
        do j=1,sem%n_glob_nodes
            write(sem%fileId,'(3(1X,E14.7))') ((sem%Coord_nodes(i,j)),i=1,sem%n_dime)
        enddo

        ! nombre d'elements
        write(sem%fileId,*) sem%n_elem

        ! nombre de materiaux !changement de l'ordre // SEM2D
        write(sem%fileId,*) sem%n_mat

        if(sem%n_mat > 1) then ! A modifier
            write(6,*) 'Trop de materiaux. 1 seul materiau prevu dans cette version'
            stop
        else
            do i = 0, sem%n_elem - 1
                write(sem%fileId,*) 0
            enddo
        endif
        ! nombre de sommets par element
        write(sem%fileId,*) sem%n_nodes

        ! Modif GS - Pas de commentaire en 3D ! une ligne de commentaire
        ! write(sem%fileId,'(A)')"! For any element, listed the reference to the global list of elements, the faces, the material"

        ! pour chaque element, les noeuds, les faces, les sommets, le materiau
        ! attention decalage de la numerotation, elle demarre a zero

        sem%control_nodes = sem%control_nodes - 1
        sem%Near_Faces = sem%Near_Faces - 1
        sem%Near_Edges = sem%Near_Edges - 1
        sem%near_Vertex = sem%near_Vertex - 1

        if (sem%typeElem==8) then
            do j=1,sem%n_elem
                write(sem%fileId,'(8(1X,I7))') (sem%control_nodes(i,j),i=1,sem%n_nodes)
            enddo
        else
            print*,"ecrireFichierSem : type d element non traite", sem%typeElem
            stop
        endif

        ! nombre de faces
        write(sem%fileId,*) sem%n_face

        do i = 1, sem%n_elem
            write(sem%fileId,*) (sem%Near_Faces(j,i),j=0,5)
            !      write(sem%fileId,*) 0, 0, 0, 0, 0, 0 ! A modifier GS
            write(sem%fileId,*) (sem%Orient_Face(j,i),j=0,5)
        enddo

        ! nombre d'aretes
        write(sem%fileId,*) sem%n_edge ! nb d'aretes par element ???

        do i = 1, sem%n_elem
            write(sem%fileId,'(12(I6,1X))') (sem%Near_Edges(j,i),j=0,11)
            write(sem%fileId,'(12(I6,1X))') (sem%Orient_Edges(j,i),j=0,11)
            !! il faut definir sem%Orient_Edges
            !      write(sem%fileId,*) sem%Orient_Edges(0:5)
            !     write(sem%fileId,'(12(I2,1X))') 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ! A modifier GS
        enddo

        ! nombre de sommets
        write(sem%fileId,*) sem%n_vertex !nb de sommets par element ???

        do i = 1, sem%n_elem
            write(sem%fileId,*) (sem%Near_Vertex(j,i),j=1,8)
        enddo

        ! une ligne de commentaire
        write(sem%fileId,*)"!!!!!!!!!!!!!!!"

        write(sem%fileId,*)"Define super_objects for communications"

        write(sem%fileId,*) .false. ! Pas de super objets par definition - A modifier

        ! une ligne de commentaire
        write(sem%fileId,*)"!!!!!!!!!!!!!!!"

        ! une ligne de commentaire
        write(sem%fileId,*) "Neumann condition?"

        write(sem%fileId,*) .false. !Pas de condition de Neumann pour l'instant - A modifier

        ! une ligne de commentaire
        write(sem%fileId,*)"!!!!!!!!!!!!!!!"
        ! nombre de processeurs (forcement 1)
        write(sem%fileId,*)"       1          Number of processors"

        ! nombre de processeurs (forcement 1)
        write(sem%fileId,*) 0, 0, 0, 0, 0, 0, 0
        ! une ligne de commentaire
        write(sem%fileId,*)"!!!!!!!!!!!!!!!"
        ! une ligne de commentaire
        write(sem%fileId,*) "Save_Surface?"
        write(sem%fileId,*) .false. !save_surface=.false. par defaut - A modifier
        close(sem%fileId)


    end subroutine ecrireFichierSem
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !


    !>
    !! \fn subroutine lister_faces_sem ()
    !! \brief
    !!
    !<


    subroutine lister_faces_sem ()

        implicit none

        integer :: i, iel, P(8)

        ! tableau de correspondance num local num globale des faces pour chaque element
        allocate(sem%Near_faces(0:5,sem%n_elem))

        !numero de face initialement -1, la 1ere face est en 0
        sem%n_face = 0

        !liste initialement vide
        nullify(listeFace)


        do iel=1,sem%n_elem

            do i=1,8
                P(i) = sem%Near_Vertex(i,iel)
            enddo


            ! ATTENTION A L'ORIENTATION

            ! Version precedente
            !!    call updateListeFace(P(1), P(2), P(6), P(5), iel, 0)
            !!    call updateListeFace(P(1), P(2), P(3), P(4), iel, 1)
            !!    call updateListeFace(P(5), P(1), P(4), P(8), iel, 2)
            !!    call updateListeFace(P(5), P(6), P(7), P(8), iel, 3)
            !!    call updateListeFace(P(6), P(2), P(3), P(7), iel, 4)
            !!    call updateListeFace(P(4), P(3), P(7), P(8), iel, 5)

            call updateListeFace(P(1), P(2), P(3), P(4), iel, 0)
            call updateListeFace(P(1), P(2), P(6), P(5), iel, 1)
            call updateListeFace(P(2), P(3), P(7), P(6), iel, 2)
            call updateListeFace(P(4), P(3), P(7), P(8), iel, 3)
            call updateListeFace(P(1), P(4), P(8), P(5), iel, 4)
            call updateListeFace(P(5), P(6), P(7), P(8), iel, 5)


        enddo


    end subroutine lister_faces_sem

    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !


    subroutine lister_edges_sem ()

        implicit none

        integer :: i, iel, P(8)

        ! tableau de correspondance num local num globale des edges pour chaque element
        allocate(sem%Near_Edges(0:11,sem%n_elem))


        !numero de edges initialement -1, la 1ere edge est en 0
        sem%n_edge = 0

        !liste initialement vide
        nullify(listeEdge)


        do iel=1,sem%n_elem

            do i=1,8
                P(i) = sem%Near_Vertex(i,iel)
            enddo

            !   PRINT*,'Avt appel updatelisteedge', iel
            call updateListeEdge(P(1), P(2), iel, 0)
            call updateListeEdge(P(2), P(3), iel, 1)
            call updateListeEdge(P(3), P(4), iel, 2)
            call updateListeEdge(P(4), P(1), iel, 3)
            call updateListeEdge(P(2), P(6), iel, 4)
            call updateListeEdge(P(6), P(5), iel, 5)
            call updateListeEdge(P(5), P(1), iel, 6)
            call updateListeEdge(P(3), P(7), iel, 7)
            call updateListeEdge(P(7), P(6), iel, 8)
            call updateListeEdge(P(7), P(8), iel, 9)
            call updateListeEdge(P(8), P(4), iel, 10)
            call updateListeEdge(P(8), P(5), iel, 11)

        enddo


    end subroutine lister_edges_sem
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !


    subroutine updateListeEdge(noeud1, noeud2, iel, numlocalEdge)

        integer :: noeud1, noeud2 ! numero globaux de 2 noeuds composant une edge
        integer :: iel ! numero global d element
        integer :: numlocalEdge ! numero local de la edge de element iel

        logical :: trouveEdge,trouveElem
        integer :: n1, n2 ! numero global de noeud tq n1 < n2
        type(Tedge), pointer :: nouvelleEdge, chercheEdge
        type(TelementDArete), pointer :: nouvelElement, chercheElement


        ! les numeros de noeuds sont ranges en ordre croissant dans la liste des edges, pour faciliter l'algo de recherche


        !PRINT *,'Entree updatelisteedge',noeud1,noeud2,iel,numlocaledge
        n1 = MIN(noeud1, noeud2)
        n2 = MAX(noeud1, noeud2)

        !n1 = noeud1
        !n2 = noeud2

        trouveEdge = .FALSE.
        trouveElem = .FALSE.

        chercheEdge => listeEdge !listeEdge est un pointeur sur un Tedge (dont les champs sont globnum, n1, n2, listelement (pointeur sur TelementdArete), suivant (pointeur sur Tedge) donc chercheedge est comme une instance de Tedge

        ! parcourir les elements de listeEdge pour savoir si les noeuds n1 et n2 ont deja ete stockes
        do
            if (.not.associated(chercheEdge)) exit ! on est en fin de liste
            if (chercheEdge%n1==n1 .and. chercheEdge%n2==n2) then ! edge deja prise en compte
                trouveEdge = .TRUE.
                ! si elle a ete prise en compte pour iel, rien de plus a faire
                chercheElement => chercheEdge%listeElement
                do
                    if (.not.associated(chercheElement)) exit
                    if (chercheElement%numElem == iel) then
                        trouveElem =.TRUE.
                        exit
                    endif
                    chercheElement=>chercheElement%suivant
                enddo

                if (.not.trouveElem) then
                    ! elle a ete prise en compte mais pas pour iel, il faut completer listeElement
                    allocate(nouvelElement)
                    nouvelElement%numElem = iel
                    nouvelElement%numlocalEdge = numlocalEdge
                    sem%Near_edges(numlocalEdge,iel) = chercheEdge%globnum
                    nouvelElement%suivant => chercheEdge%listeElement
                    chercheEdge%listeElement => nouvelElement
                endif
                exit
            endif
            chercheEdge => chercheEdge%suivant
        enddo

        if (trouveEdge) return
        !PRINT*,'passage listeedge'

        ! nouvelle edge
        sem%n_edge = sem%n_edge + 1
        allocate(nouvelleEdge)
        nouvelleEdge%n1 = n1
        nouvelleEdge%n2 = n2
        nouvelleEdge%globnum  =  sem%n_edge
        nullify(nouvelleEdge%listeElement)
        nouvelleEdge%suivant => listeEdge
        listeEdge => nouvelleEdge

        allocate(nouvelElement)
        nouvelElement%numElem = iel
        nouvelElement%numlocalEdge = numlocalEdge
        nouvelElement%suivant => listeEdge%listeElement
        listeEdge%listeElement => nouvelElement

        sem%Near_edges(numlocalEdge,iel) = sem%n_edge
    end subroutine updateListeEdge
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !
    subroutine updateListeFace(noeud1, noeud2, noeud3, noeud4, iel, numlocalFace)

        integer :: noeud1, noeud2, noeud3, noeud4 ! numero globaux de 2 noeuds composant une face
        integer :: iel ! numero global d element
        integer :: numlocalFace ! numero local de la face de element iel

        logical :: trouveFace,trouveElem
        integer :: n1, n2, n3, n4 ! numero global de noeud tq n1 < n2
        type(Tface),pointer :: nouvelleFace,chercheFace
        type(TelementDeFace),pointer :: nouvelElement, chercheElement


        ! les numeros de noeuds sont ranges en ordre croissant dans la liste des faces, pour faciliter l'algo de recherche
        ! VOIR S'IL FAUT RENUMEROTER LES NOEUDS
        ! if (noeud1.lt.noeud2) then
        !    n1 = noeud1
        !    n2 = noeud2   !! n1 = MIN(n1, n2) ; n2 = MAX(n1,n2)
        ! else
        !    n1 = noeud2
        !    n2 = noeud1
        ! endif

        n1 = noeud1;  n2 = noeud2;  n3 = noeud3;  n4 = noeud4 !! A VERIFIER ET A MODIFIER
        n1 = MIN(noeud1, noeud2, noeud3, noeud4)
        n4 = MAX(noeud1, noeud2, noeud3, noeud4)
        if(n1==noeud1) then
            n2 = MIN(noeud2, noeud3, noeud4)
        elseif(n1==noeud2) then
            n2 = MIN(noeud1, noeud3, noeud4)
        elseif(n1==noeud3) then
            n2 = MIN(noeud1, noeud2, noeud4)
        elseif(n1==noeud4) then
            n2 = MIN(noeud1, noeud2, noeud3)
        endif

        if(n4==noeud1) then
            n3 = MAX(noeud2, noeud3, noeud4)
        elseif(n4==noeud2) then
            n3 = MAX(noeud1, noeud3, noeud4)
        elseif(n4==noeud3) then
            n3 = MAX(noeud1, noeud2, noeud4)
        elseif(n4==noeud4) then
            n3 = MAX(noeud1, noeud2, noeud3)
        endif
        IF(.NOT.((n4>n3).AND.(n3>n2).AND.(n2>n1))) THEN
            WRITE(6,*) 'Pb detection des sommets des faces'
            STOP
        ENDIF
        trouveFace = .FALSE.
        trouveElem = .FALSE.

        chercheFace => listeFace !listeFace est un pointeur sur un Tface (dont les champs sont globnum, n1, n2, listelement (pointeur sur Telementdeface), suivant (pointeur sur Tface) donc chercheface est comme une instance de Tface

        ! parcourir les elements de listeFace pour savoir si les noeuds n1 et n2 ont deja ete stockes
        do
            if (.not.associated(chercheFace)) exit ! on est en fin de liste
            if (chercheFace%n1==n1 .and. chercheFace%n2==n2 .and. chercheFace%n3==n3 .and. chercheFace%n4==n4) then ! face deja prise en compte
                trouveFace = .TRUE.
                ! si elle a ete prise en compte pour iel, rien de plus a faire
                chercheElement => chercheFace%listeElement
                do
                    if (.not.associated(chercheElement)) exit
                    if (chercheElement%numElem.eq.iel) then
                        trouveElem = .TRUE.
                        exit
                    endif
                    chercheElement=>chercheElement%suivant
                enddo
                if (.not.trouveElem) then
                    ! elle a ete prise en compte mais pas pour iel, il faut completer listeElement
                    allocate(nouvelElement)
                    nouvelElement%numElem = iel
                    nouvelElement%numlocalFace = numlocalFace
                    sem%Near_faces(numlocalFace,iel) = chercheFace%globnum
                    nouvelElement%suivant => chercheFace%listeElement
                    chercheFace%listeElement => nouvelElement
                endif
                exit
            endif
            chercheFace => chercheFace%suivant
        enddo

        if (trouveFace) return

        ! nouvelle face
        sem%n_face = sem%n_face+1
        allocate(nouvelleFace)
        nouvelleFace%n1 = n1
        nouvelleFace%n2 = n2
        nouvelleFace%n3 = n3
        nouvelleFace%n4 = n4
        nouvelleFace%globnum = sem%n_face
        nullify(nouvelleFace%listeElement)
        nouvelleFace%suivant => listeFace
        listeFace=>nouvelleFace

        allocate(nouvelElement)
        nouvelElement%numElem = iel
        nouvelElement%numlocalFace = numlocalFace
        nouvelElement%suivant => listeFace%listeElement
        listeFace%listeElement => nouvelElement

        sem%Near_faces(numlocalFace,iel) = sem%n_face

    end subroutine updateListeFace

end module ipg2
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
