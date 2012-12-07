!>
!!\file sem.F90
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
       integer :: n_face       ! nombre de face du maillage
       integer :: n_mat        ! nombre de materiaux
       integer :: typeElem     ! typeElem=4 pour quad4 ou 8 pour quad8

       !!  integer, dimension (:), allocatable   :: glob_numbering ! numeros globaux des noeuds
       !!
       !!  integer, dimension (:,:), allocatable :: Near_Face      ! numeros globaux des faces de l'element
       !!  integer, dimension (:,:), allocatable :: Near_Vertex    ! numeros globaux des noeuds de l'element
       !!  integer, dimension (:,:), allocatable :: Control_nodes  ! numeros des noeuds (numerotation sequentielle des points de coordonnees X Y)
       !!  integer, dimension (:,:), allocatable :: face_element   ! numeros du ou des elements de la face
       !!  integer, dimension (:,:), allocatable :: which_face     ! numero local du ou des elements de la face
       !!  integer, dimension (:,:), allocatable :: face_vertex    ! numero local du ou des elements de la face
       !!
       !!  real   , dimension (:,:), allocatable :: Coord_nodes ! coordonnees X,Y des noeuds
       !!

       integer, dimension (:), pointer   :: glob_numbering ! numeros globaux des noeuds

       integer, dimension (:,:), pointer :: Near_Face      ! numeros globaux des faces de l'element
       integer, dimension (:,:), pointer :: Near_Vertex    ! numeros globaux des noeuds de l'element
       integer, dimension (:,:), pointer :: Control_nodes  ! numeros des noeuds (numerotation sequentielle des points de coordonnees X Y)
       integer, dimension (:,:), pointer :: face_element   ! numeros du ou des elements de la face
       integer, dimension (:,:), pointer :: which_face     ! numero local du ou des elements de la face
       integer, dimension (:,:), pointer :: face_vertex    ! numero local du ou des elements de la face

       real   , dimension (:,:), pointer :: Coord_nodes ! coordonnees X,Y des noeuds


       character (len=110) :: fichier

    end type Tsem
    character (len=MAX_FILE_SIZE) :: fnamef

    type TelementDeFace
       integer :: numElem
       integer :: numLocalFace
       !integer, dimension(0:3) :: face ! tableau des numeros globaux de faces
       type(TelementDeFace),pointer :: suivant
    end type TelementDeFace


    type Tface
       integer :: globnum ! numero global de la face
       integer :: n1,n2 ! numero globaux des noeuds de la face
       type(TelementDeFace),pointer :: listeElement
       type(Tface), pointer :: suivant
    end type Tface


    type Tfrontiere
       integer :: iface ! numero de la face sur la frontiere
       integer :: iVertex1,iVertex2 ! numeros des sommets de iface
       integer :: iCtrlNode1,iCtrlNode2 ! numeros des pts de controle corres. a ivertex1,2
       real    :: xVertex1,yVertex1,xVertex2,yVertex2 ! coordonnees X,Y des vertex
       logical :: coherency ! a definir
       type(Tfrontiere),pointer :: suivant
    end type Tfrontiere


    type(Tsem) :: sem

    type(Tface),pointer :: listeFace

    type(Tfrontiere),pointer :: listeFrontiere

contains
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !

    !>
    !! \brief
    !!
    !<


    subroutine ecrireFichierSemFrontiere

        implicit none

        integer :: i,iface,nFaceFrontiere,nVertexFrontiere,CodeErreur
        type(Tfrontiere),pointer :: nouvelleFace,chercheFace

        integer, dimension(:,:),allocatable :: listeVertex


        !liste initialement vide
        nullify(listeFrontiere) !listeFrontiere est un pointeur sur un type Tfrontiere

        ! nbre de faces initiallement nul
        nFaceFrontiere = 0

        ! nbre de vertex initiallement nul
        nVertexFrontiere = 0


        ! boucle sur les faces
        do iface=1,sem%n_face
            if ( sem%face_Element(1,iface) .ne. -1 ) cycle ! ce n est pas une face frontiere, on passe

            allocate(nouvelleFace)

            nouvelleFace%iface = iface

            ! recup des vertex associes a iface
            nouvelleFace%ivertex1 = sem%face_vertex(0,iface)
            nouvelleFace%ivertex2 = sem%face_vertex(1,iface)

            ! recup des points de controle
            nouvelleFace%iCtrlNode1 = sem%glob_numbering(nouvelleFace%ivertex1)-1
            nouvelleFace%iCtrlNode2 = sem%glob_numbering(nouvelleFace%ivertex2)-1

            ! recup des coordonnees des vertex
            nouvelleFace%xVertex1 = sem%Coord_nodes(1,nouvelleFace%iCtrlNode1)
            nouvelleFace%yVertex1 = sem%Coord_nodes(2,nouvelleFace%iCtrlNode1)
            nouvelleFace%xVertex2 = sem%Coord_nodes(1,nouvelleFace%iCtrlNode2)
            nouvelleFace%yVertex2 = sem%Coord_nodes(2,nouvelleFace%iCtrlNode2)


            ! coherency
            nouvelleFace%coherency = .false.
            if (nouvelleFace%iCtrlNode1 .lt. nouvelleFace%iCtrlNode2) nouvelleFace%coherency=.true.

            ! compteur de faces
            nFaceFrontiere = nFaceFrontiere + 1

            ! on empile la nouvelle face en tete de liste
            nouvelleFace%suivant=>listeFrontiere

            ! la liste demarre sur cette nouvelle face
            listeFrontiere=>nouvelleFace

        enddo

        ! ouverture du fichier frontiere
        call semname_sem_frontiere(sem%fichier,fnamef)
        open(UNIT=sem%fileId,IOSTAT=CodeErreur,FILE=fnamef,FORM='formatted',STATUS='replace',ACTION='write')
        if (CodeErreur .ne.0 ) print*,'Ouverture du fichier SEM frontiere  :CodeErreur=',CodeErreur


        ! tableau ds lequel on va stocker la liste exacte (sans doublon) des vertex frontieres
        allocate(listeVertex(2*nFaceFrontiere,3))


        ! ecriture des donnees de faces frontieres

        write(sem%fileId,*) nFaceFrontiere ! nombre de faces frontieres

        ! on parcourt la liste frontiere
        chercheFace=>listeFrontiere

        do
            if (.not.associated(chercheFace)) exit ! on est en fin de liste

            ! ecriture du numero de la face
            write(sem%fileId,*) chercheFace%iface, &
                chercheFace%ivertex1,chercheFace%ivertex2, &
                chercheFace%coherency

            chercheFace=>chercheFace%suivant

            ! mise a jour de la liste des vertex frontieres
            if (.not.VertexFrontiereDejaCompte(chercheFace%ivertex1,nVertexFrontiere,listeVertex(1:nVertexFrontiere,1))) then
                nVertexFrontiere = nVertexFrontiere +1
                listeVertex(nVertexFrontiere,1) = chercheFace%ivertex1
                listeVertex(nVertexFrontiere,2) = chercheFace%xVertex1
                listeVertex(nVertexFrontiere,3) = chercheFace%yVertex1
            endif

            if (.not.VertexFrontiereDejaCompte(chercheFace%ivertex2,nVertexFrontiere,listeVertex(1:nVertexFrontiere,1))) then
                nVertexFrontiere = nVertexFrontiere +1
                listeVertex(nVertexFrontiere,1) = chercheFace%ivertex2
                listeVertex(nVertexFrontiere,2) = chercheFace%xVertex2
                listeVertex(nVertexFrontiere,3) = chercheFace%yVertex2
            endif

        enddo


        write(sem%fileId,*)"------------------"
        ! ecriture des donnees de aux faces frontieres
        write(sem%fileId,*) nVertexFrontiere ! nombre de faces frontieres

        ! on parcourt la liste des vertex

        do i=1,nVertexFrontiere

            ! ecriture du numero de la face
            write(sem%fileId,*) listeVertex(i,1), &
                listeVertex(i,2),listeVertex(i,3)

        enddo


        deallocate(listeVertex)

        close(sem%fileId)

    end subroutine ecrireFichierSemFrontiere

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
    !! \brief
    !!
    !<


    subroutine ecrireFichierSem

        implicit none

        integer :: i,j,temp,CodeErreur

        call semname_sem_fichier(sem%fichier,fnamef)
        open(UNIT=sem%fileId,IOSTAT=CodeErreur,FILE=fnamef,FORM='formatted',STATUS='replace',ACTION='write')
        if (CodeErreur .ne.0 ) print*,'Ouverture du fichier SEM  :CodeErreur=',CodeErreur



        ! dimension du pb (2D)
        write(sem%fileId,*) sem%n_dime

        ! nombre de noeuds
        write(sem%fileId,*) sem%n_glob_nodes

        ! coordonnees X Y
        do i=1,sem%n_glob_nodes
            write(sem%fileId,'(2(1X,E14.7))') sem%Coord_nodes(1,i),sem%Coord_nodes(2,i)
        enddo

        ! nombre de materiaux
        write(sem%fileId,*)sem%n_mat

        ! nombre et nom de lignes ? => supprimes car inutilises!
        write(sem%fileId,*) "1"
        write(sem%fileId,*)" "
        write(sem%fileId,*)"       1          L1"

        ! nombre d'elements
        write(sem%fileId,*) sem%n_elem

        ! nombre de sommets par element
        write(sem%fileId,*) sem%n_nodes

        ! une ligne de commentaire
        write(sem%fileId,'(A)')"! For any element, listed the reference to the global list of elements, the faces, the material"

        ! pour chaque element, les noeuds, les faces, les sommets, le materiau
        ! attention decalage de la numerotation, elle demarre a zero
        do j=1,sem%n_elem
            do i=1,sem%typeElem
                sem%control_nodes(i,j)=sem%control_nodes(i,j)-1
            enddo
        enddo

        do j=1,sem%n_elem
            do i=1,4
                sem%near_vertex(i,j)=sem%near_vertex(i,j)-1
            enddo
        enddo


        do j=1,sem%n_elem
            do i=0,3
                sem%Near_Face(i,j)=sem%Near_Face(i,j)-1
            enddo
        enddo



        do j=1,sem%n_elem

            if (sem%typeElem==4) then
                write(sem%fileId,'(13(1X,I7))') (sem%control_nodes(i,j),i=1,4),&
                    & (sem%Near_Face(i,j),i=0,3)  ,&
                    & (sem%near_vertex(i,j),i=1,4)  ,&
                    & 0 ! + numero du materiau
            elseif (sem%typeElem==8) then
                print*,"ecriture des faces",sem%Near_Face(0,j),sem%Near_Face(1,j),sem%Near_Face(2,j),sem%Near_Face(3,j)
                write(sem%fileId,'(17(1X,I7))') (sem%control_nodes(i,j),i=1,8),&
                    &  (sem%Near_Face(i,j),i=0,3)  ,&
                    &  (sem%near_vertex(i,j),i=1,4)  ,&
                    &  0 ! + numero du materiau
            else
                print*,"ecrireFichierSem : type d element non traite"
                stop
            endif

        enddo


        ! une ligne de commentaire
        write(sem%fileId,*)"!!!!!!!!!!!!!!!"

        ! nombre de faces
        write(sem%fileId,*) sem%n_face

        ! une ligne de commentaire
        write(sem%fileId,'(A)')"! For any face, listed to the elements, what face of the element and the vertices"

        ! pour chaque face, les 2 elements, les 2 numeros locaux de face, les 2 noeuds
        ! attention decalage de la numerotation, elle demarre a zero
        do j=1,sem%n_face
            do i=0,1
                if (sem%face_Element(i,j).ne.-1) sem%face_Element(i,j)=sem%face_Element(i,j)-1
                if (sem%face_vertex(i,j).ne.-1) sem%face_vertex(i,j)=sem%face_vertex(i,j)-1
            enddo


            ! il faut ranger les elements de la face en ordre croissant sauf pour el=-1 qui reste en 2nd position
            if (sem%face_Element(0,j).gt.sem%face_Element(1,j).and.sem%face_Element(1,j).ne.-1) then
                temp = sem%face_Element(0,j)
                sem%face_Element(0,j) = sem%face_Element(1,j)
                sem%face_Element(1,j) = temp

                ! les numero locaux des faces sont donnes dans le meme ordre que les elements donc on les rearrange
                ! en meme temps
                temp = sem%which_face(0,j)
                sem%which_face(0,j) = sem%which_face(1,j)
                sem%which_face(1,j) = temp

            endif

        enddo


        do j=1,sem%n_face
            write(sem%fileId,'(6(1X,I7))') (sem%face_Element(i,j),i=0,1),&
                (sem%which_face(i,j),i=0,1),&
                (sem%face_vertex(i,j),i=0,1)
        enddo

        ! une ligne de commentaire
        write(sem%fileId,*)"!!!!!!!!!!!!!!!"


        ! nombre de sommet
        write(sem%fileId,*) sem%n_vertex

        ! une ligne de commentaire
        write(sem%fileId,'(A)')"! For any vertex, address to the global numbering"



        !!

!!!!!  open(UNIT=77,IOSTAT=CodeErreur,FILE=trim(sem%fichier)//"corres",FORM='formatted',STATUS='new',ACTION='write')

        ! numerotation globale des noeuds
        do i=1,sem%n_vertex
            write(sem%fileId,*) sem%glob_numbering(i)-1
        enddo

!!!!  close(77)



        ! une ligne de commentaire
        write(sem%fileId,*)"!!!!!!!!!!!!!!!"

        write(sem%fileId,*)"Define super_objects for communications"

        ! nombre de processeurs (forcement 1)
        write(sem%fileId,*)"       1          Number of processors"


        ! nombre de communications (forcement 0)
        write(sem%fileId,*)"       0          Number of effective processors in communication"
        ! la suite ne concerne que les cas // et est donc sans objet ici

        close(sem%fileId)


    end subroutine ecrireFichierSem
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !


    !>
    !! \brief
    !!
    !<


    subroutine lister_faces_sem ()

        implicit none

        integer :: iel,n1,n2,n3,n4

        ! tableau de correspondance num local num globale des faces pour chaque element
        allocate(sem%Near_face(0:3,sem%n_elem))


        !numero de faces initialement -1, la 1ere face est en 0
        sem%n_face = 0

        !liste initialement vide
        nullify(listeFace)


        do iel=1,sem%n_elem

            if (sem%typeElem==4) then

                n1 = sem%Near_Vertex(1,iel)
                n2 = sem%Near_Vertex(2,iel)
                n3 = sem%Near_Vertex(3,iel)
                n4 = sem%Near_Vertex(4,iel)

                call updateListeFace(n1,n2,iel,0)
                call updateListeFace(n2,n3,iel,1)
                call updateListeFace(n3,n4,iel,2)
                call updateListeFace(n4,n1,iel,3)

            elseif (sem%typeElem==8) then

                n1 = sem%Near_Vertex(1,iel)
                n2 = sem%Near_Vertex(2,iel)
                n3 = sem%Near_Vertex(3,iel)
                n4 = sem%Near_Vertex(4,iel)

                call updateListeFace(n1,n2,iel,0)
                call updateListeFace(n2,n3,iel,1)
                call updateListeFace(n3,n4,iel,2)
                call updateListeFace(n4,n1,iel,3)

            endif

        enddo


    end subroutine lister_faces_sem

    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !


    subroutine updateListeFace(noeud1,noeud2,iel,numlocalFace)

        integer :: noeud1,noeud2 ! numero globaux de 2 noeuds composant une face
        integer :: iel ! numero global d element
        integer :: numlocalFace ! numero local de la face de element iel

        logical :: trouveFace,trouveElem
        integer :: n1,n2 ! numero global de noeud tq n1 < n2
        type(Tface),pointer :: nouvelleFace,chercheFace
        type(TelementDeFace),pointer :: nouvelElement, chercheElement


        ! les numeros de noeuds sont ranges en ordre croissant dans la liste des faces, pour faciliter l'algo de recherche
        if (noeud1.lt.noeud2) then
            n1 = noeud1
            n2 = noeud2   !! n1 = MIN(n1, n2) ; n2 = MAX(n1,n2)
        else
            n1 = noeud2
            n2 = noeud1
        endif



        trouveFace = .FALSE.
        trouveElem = .FALSE.

        chercheFace => listeFace !listeFace est un pointeur sur un Tface (dont les champs sont globnum, n1, n2, listelement (pointeur sur Telementdeface), suivant (pointeur sur Tface) donc chercheface est comme une instance de Tface

        ! parcourir les elements de listeFace pour savoir si les noeuds n1 et n2 ont deja ete stockes
        do
            if (.not.associated(chercheFace)) exit ! on est en fin de liste
            if (chercheFace%n1.eq.n1.and.chercheFace%n2.eq.n2) then ! face deja prise en compte
                trouveFace=.TRUE.
                ! si elle a ete prise en compte pour iel, rien de plus a faire
                chercheElement => chercheFace%listeElement
                do
                    if (.not.associated(chercheElement)) exit
                    if (chercheElement%numElem.eq.iel) then
                        trouveElem=.TRUE.
                        exit
                    endif
                    chercheElement=>chercheElement%suivant
                enddo
                if (.not.trouveElem) then
                    ! elle a ete prise en compte mais pas pour iel, il faut completer listeElement
                    allocate(nouvelElement)
                    nouvelElement%numElem=iel
                    nouvelElement%numlocalFace=numlocalFace
                    sem%Near_face(numlocalFace,iel)=chercheFace%globnum
                    nouvelElement%suivant=>chercheFace%listeElement
                    chercheFace%listeElement=>nouvelElement

                endif
                exit
            endif
            chercheFace => chercheFace%suivant
        enddo

        if (trouveFace) return

        ! nouvelle face
        sem%n_face=sem%n_face+1
        allocate(nouvelleFace)
        nouvelleFace%n1=n1
        nouvelleFace%n2=n2
        nouvelleFace%globnum = sem%n_face
        nullify(nouvelleFace%listeElement)
        nouvelleFace%suivant=>listeFace
        listeFace=>nouvelleFace

        allocate(nouvelElement)
        nouvelElement%numElem=iel
        nouvelElement%numlocalFace=numlocalFace
        nouvelElement%suivant=>listeFace%listeElement
        listeFace%listeElement=>nouvelElement

        sem%Near_face(numlocalFace,iel)=sem%n_face

    end subroutine updateListeFace


end module ipg2
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
