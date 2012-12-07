!>
!!\file convertirUnv.f90
!!\brief Convertit un maillage universel en maillage ipg2.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \fn subroutine convertirUnv(tDomain)
!! \brief
!!
!! \param type (domain), intent (INOUT) Tdomain
!<


subroutine convertirUnv(tDomain, rg)


    ! convertit un maillage universel en maillage ipg2
    !

    use sdomain
    use unv
    use ipg2

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: rg

    ! Initialisation des fid
    sunv%fileId = 10
    sem%fileId  = 11

    ! Récupération du nom du fichier unv
    write (sunv%fichier, "(a,a5,I4.4)") trim(adjustl(Tdomain%mesh_file)), ".unv.", rg
    PRINT*,'sunv%fichier=',sunv%fichier

    ! Récupération du nom du fichier sem
    write (sem%fichier, "(a,a1,I2.2)") trim(adjustl(Tdomain%mesh_file)), ".", rg

    call lireFichierUnv
    write(6,*) 'End of reading the .unv mesh file'

    call convertir
    write(6,*) 'End of the conversion from unv to sem -> no writing'

    call ecrireFichierSem
    write(6,*) 'End of the conversion from unv to sem -> writing done'

    Tdomain%mesh_file = sem%fichier
    !  if (Tdomain%n_proc.gt.1) call ecrireFichierSemFrontiere

end subroutine convertirUnv
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
!>
!! \fn subroutine convertir
!! \brief
!!
!<


subroutine convertir
    use unv
    use ipg2

    implicit none

    type(Tface),pointer :: chercheFace
    type(TelementDeFace),pointer :: chercheElement
    type(Tedge),pointer :: chercheEdge
    type(TelementDArete), pointer :: chercheElementdArete

    integer i,j,indice,iNumLocal
    integer, dimension(:), allocatable :: global2local

    real P(4,3), Q(4,3)
    integer index0(2)
    integer iref, jref, iel, ii, ij
    integer TAB(0:5,1:4), TAB2(0:11,1:2)
    real P1(3), Q1(3), Q2(3), normeP
    real, parameter :: epsilon = 1e-7
    ! dimension du pb : 3D en dur, a voir
    sem%n_dime = 3

    ! type element : 8, 20, 27
    sem%typeElem = sunv%typeElem


    ! recuperation du nombre de noeuds
    sem%n_glob_nodes = sunv%n_noeud

    ! recuperation des coordonnees 2D des noeuds
    allocate(sem%Coord_nodes(3,sem%n_glob_nodes))
    sem%Coord_nodes = sunv%noeudCoord

    ! recuperation du nombre de materiaux ! A modifier GS
    sem%n_mat = 1 ! pour l'instant


    ! recuperation du nombre d'elements
    sem%n_elem = sunv%n_elem


    ! recuperation du nombre de sommets par element
    sem%n_nodes = sunv%typeElem


    ! recuperation des connectivites noeuds elements
    ! a priori, pas de double numerotation sur les sommets
    ! dans le fichier UNV
    allocate(sem%control_nodes(sem%n_nodes,sem%n_elem))
    !allocate(sem%near_vertex(sem%n_nodes,sem%n_elem))
    allocate(sem%near_vertex(8,sem%n_elem))

    allocate (global2local(sem%n_glob_nodes))
    global2local(:)=-1 ! par defaut un faux numero local
    iNumLocal = 1 ! numero du premier vertex

    do j=1,sem%n_elem

        if (sem%typeElem==8) then ! les sommets st dans le meme ordre : 1 2 3 4
            do i=1,8
                sem%control_nodes(i,j) = sunv%elemConnec(i,j)
                sem%near_vertex(i,j)   = sunv%elemConnec(i,j)
            enddo


            ! MISE EN COMMENTAIRE SI sem%typeElem==20 (8 en 2D)
!!!        elseif (sem%typeElem==8) then  ! les sommets ne sont pas dans le meme ordre :
!!!           ! 1 2 3 4 5 6 7 8 pour unv
!!!           ! 1 3 5 7 2 4 6 8 pour sem2d
            !!
            !!
!!!        sem%control_nodes(1,j)=sunv%elemConnec(1,j)
!!!        sem%control_nodes(2,j)=sunv%elemConnec(3,j)
!!!        sem%control_nodes(3,j)=sunv%elemConnec(5,j)
!!!        sem%control_nodes(4,j)=sunv%elemConnec(7,j)
!!!        sem%control_nodes(5,j)=sunv%elemConnec(2,j)
!!!        sem%control_nodes(6,j)=sunv%elemConnec(4,j)
!!!        sem%control_nodes(7,j)=sunv%elemConnec(6,j)
!!!        sem%control_nodes(8,j)=sunv%elemConnec(8,j)
            !!
!!!        if (global2local(sem%control_nodes(1,j)).eq.-1) then
!!!           global2local(sem%control_nodes(1,j)) = iNumLocal
!!!           iNumLocal=iNumLocal+1
!!!        endif
!!!        sem%near_vertex(1,j) = global2local(sem%control_nodes(1,j))
            !!
            !!
!!!        if (global2local(sem%control_nodes(2,j)).eq.-1) then
!!!           global2local(sem%control_nodes(2,j)) = iNumLocal
!!!           iNumLocal=iNumLocal+1
!!!        endif
!!!        sem%near_vertex(2,j) = global2local(sem%control_nodes(2,j))
            !!
!!!        if (global2local(sem%control_nodes(3,j)).eq.-1) then
!!!           global2local(sem%control_nodes(3,j)) = iNumLocal
!!!           iNumLocal=iNumLocal+1
!!!        endif
!!!        sem%near_vertex(3,j) = global2local(sem%control_nodes(3,j))
            !!
            !!
            !!
!!!        if (global2local(sem%control_nodes(4,j)).eq.-1) then
!!!           global2local(sem%control_nodes(4,j)) = iNumLocal
!!!           iNumLocal=iNumLocal+1
!!!        endif
!!!        sem%near_vertex(4,j) = global2local(sem%control_nodes(4,j))
            !!
            ! FIN MISE EN COMMENTAIRE SI sem%typeElem==20 (8 en 2D)
        endif


    enddo

    ! nombre de vertex
    if (sem%typeElem==8) sem%n_vertex=sem%n_glob_nodes

!!! MISE EN COMMENTAIRE SI sem%typeElem==20 (8 en 2D)
!!!  if (sem%typeElem==8) sem%n_vertex=iNumLocal-1
    ! FIN MISE EN COMMENTAIRE SI sem%typeElem==20 (8 en 2D)


    ! creer les faces et leurs numerotations (locale et globale)
    call lister_faces_sem

    ! creer les aretes et leurs numerotations (locale et globale)
    call lister_edges_sem

    ! recuperation du numero de materiau associe a chaque element


    !! MISE EN COMMENTAIRE DE CE PASSAGE CAR PAS FACILE
    ! recuperation des numeros globaux des 2 elements constituant chaque face
    ! recuperation des numeros locaux des 2 elements constituant chaque face
    ! recuperation des numeros globaux des 2 sommets constituant chaque face

    !!  allocate(sem%which_face(0:1,sem%n_face))
    allocate(sem%Orient_face(0:5,sem%n_elem))
    allocate(sem%Orient_edges(0:11,sem%n_elem))
    allocate(sem%face_Element(0:1,sem%n_face))
    allocate(sem%face_vertex(0:3,sem%n_face))
    allocate(sem%edge_vertex(0:1,sem%n_edge))

    ! recuperation du nombre de faces
    !sem%n_face = sunv%n_face

    sem%face_Element(:,:) = -1
    !!  sem%which_face(:,:) = -1

    allocate (chercheFace)
    chercheFace => ListeFace ! on liste les faces
    do
        if (.not.associated(chercheFace)) exit !  fin de parcours de la liste des faces
        ! pour chaque face, on parcourt sa liste des elements
        indice=0
        allocate(chercheElement)
        chercheElement => chercheFace%listeElement
        do
            if (.not.associated(chercheElement)) exit ! fin de parcours de la liste des elements
            sem%face_Element(indice,chercheFace%globnum)=chercheElement%numElem
            !        sem%which_face(indice,chercheFace%globnum)=chercheElement%numlocalFace

            !! GS: Mettre ici la definition de Orient_faces en utilisant numlocalface ?
!!!        if(indice==1) then
!!!           select case (chercheElement%numlocalface)
!!!              case(0)
!!!           if (chercheElement%suivant%numlocalface ==2) then
!!!           et si sens different
!!!                 Orient_Faces(0,chercheElement%numElem)=1
!!!           end select
!!!        endif
            if (.not.associated(chercheElement%suivant)) exit !  fin de parcours de la liste des elements
            chercheElement => chercheElement%suivant
            indice = indice+1

        enddo

        sem%face_vertex(0,chercheFace%globnum) = chercheFace%n1
        sem%face_vertex(1,chercheFace%globnum) = chercheFace%n2
        sem%face_vertex(2,chercheFace%globnum) = chercheFace%n3
        sem%face_vertex(3,chercheFace%globnum) = chercheFace%n4

        if (.not.associated(chercheFace%suivant)) exit !  fin de parcours de la liste des faces
        chercheFace => chercheFace%suivant
    enddo

    sem%Orient_face = 0

    TAB(0,1) =1; TAB(0,2)=2; TAB(0,3)=3;  TAB(0,4)=4 !face 0 - k=0
    TAB(1,1) =1; TAB(1,2)=2; TAB(1,3)=6;  TAB(1,4)=5 !face 1 - j=0
    TAB(2,1) =2; TAB(2,2)=3; TAB(2,3)=7;  TAB(2,4)=6 !face 2 - i=n-1
    TAB(3,1) =4; TAB(3,2)=3; TAB(3,3)=7;  TAB(3,4)=8 !face 3 - j=n-1
    TAB(4,1) =1; TAB(4,2)=4; TAB(4,3)=8;  TAB(4,4)=5 !face 4 - i=0
    TAB(5,1) =5; TAB(5,2)=6; TAB(5,3)=7;  TAB(5,4)=8 !face 5 - k=n-1

    do i=1,sem%n_face
        iref = sem%face_element(0,i)
        iel = sem%face_element(1,i)
        ! iel donne le numero de l'elt de la face i pr lequel i est la face d'etude.
        ! si face frontiere (iel=-1) on ne fait rien
        if(iel/=-1) then
            do j=0,5
                if (sem%Near_faces(j,iref)==i) then
                    jref = j
                    exit
                endif
            enddo

            !        print*,'iref,iel',iref,iel,jref

            do j=0,5

                if (sem%Near_faces(j,iel)==i) then
                    !           sem%Orient_faces(j,iel) pas forcement zero



                    do ij=1,4
                        do ii=1,3
                            P(ij,ii) = sem%Coord_nodes(ii,sem%control_nodes(TAB(jref,ij),iref))
                            Q(ij,ii) = sem%Coord_nodes(ii,sem%control_nodes(TAB(j,ij),iel))
                        enddo
                    enddo
                    !              P(1:4,1:3) = transpose(sem%Coord_nodes(1:3,sem%control_nodes(TAB(jref,1:4),iref)))
                    !              Q(1:4,1:3) = transpose(sem%Coord_nodes(1:3,sem%control_nodes(TAB(j,1:4),iel)))
                    call calcorient(P, Q, index0)
                    !
                    if ((index0(1)==1) .AND. (index0(2)==2)) then
                        sem%Orient_Face(j,iel) = 0
                    elseif ((index0(1)==2) .AND. (index0(2)==3)) then
                        sem%Orient_Face(j,iel) = 5
                    elseif ((index0(1)==3) .AND. (index0(2)==4)) then
                        sem%Orient_Face(j,iel) = 3
                    elseif ((index0(1)==4) .AND. (index0(2)==1)) then
                        sem%Orient_Face(j,iel) = 6
                    elseif ((index0(1)==2) .AND. (index0(2)==1)) then
                        sem%Orient_Face(j,iel) = 1
                    elseif ((index0(1)==3) .AND. (index0(2)==2)) then
                        sem%Orient_Face(j,iel) = 7
                    elseif ((index0(1)==4) .AND. (index0(2)==3)) then
                        sem%Orient_Face(j,iel) = 2
                    elseif ((index0(1)==1) .AND. (index0(2)==4)) then
                        sem%Orient_Face(j,iel) = 4
                    else
                        write(6,*) 'Pb in the definition of Orient_faces - Stop'
                        stop
                    endif
                    exit
                endif
            enddo
        endif
    enddo

    PRint*,'Fin Orient_Faces'
    do i=1,sem%n_face
        WRITE(15,*) 'Face', i, 'Elts ',sem%face_element(0:1,i)
        WRITE(15,'(A11,12(1X,E12.5))') 'Coord', (sem%Coord_nodes(1:3,sem%face_vertex(j,i)),j=0,3)
        WRITE(15,*) 'Noeuds Fac1', sem%face_vertex(0:3,i)
        WRITE(15,'(A11,8(1X,I4))') 'Noeuds Elt1', sem%Control_nodes(1:8,sem%face_element(0,i))
        WRITE(15,'(A11,8(1X,I4))') 'Noeuds Elt2', sem%Control_nodes(1:8,sem%face_element(1,i))
        WRITE(15,*) '-------------------------'
        WRITE(15,*)
    enddo

    chercheEdge => ListeEdge ! on liste les faces
    do
        if (.not.associated(chercheEdge)) exit !  fin de parcours de la liste des faces
        ! pour chaque arete, on parcourt sa liste des elements
        indice = 0
        allocate(chercheElementdArete)
        chercheElementdArete => chercheEdge%listeElement
        do
            if (.not.associated(chercheElementdArete)) exit ! fin de parcours de la liste des elements

            if (.not.associated(chercheElementdArete%suivant)) exit !  fin de parcours de la liste des elements
            chercheElementdArete => chercheElementdArete%suivant
            indice = indice+1
        enddo

        sem%edge_vertex(0,chercheEdge%globnum) = chercheEdge%n1
        sem%edge_vertex(1,chercheEdge%globnum) = chercheEdge%n2

        if (.not.associated(chercheEdge%suivant)) exit !  fin de parcours de la liste des faces
        chercheEdge => chercheEdge%suivant
    enddo

    PRint*,'Avt debut Orient_Edges'
    sem%Orient_Edges = 0

    ! Mise en commentaire GS - A finir
    TAB2(0,1) =1; TAB2(0,2)=2
    TAB2(1,1) =2; TAB2(1,2)=3
    TAB2(2,1) =4; TAB2(2,2)=3  !2 sens possibles
    TAB2(3,1) =1; TAB2(3,2)=4  !2 sens possibles
    TAB2(4,1) =2; TAB2(4,2)=6  !2 sens possibles
    TAB2(5,1) =5; TAB2(5,2)=6  !2 sens possibles
    TAB2(6,1) =1; TAB2(6,2)=5  !pas logique avec le sens de parcours des faces
    TAB2(7,1) =3; TAB2(7,2)=7
    TAB2(8,1) =6; TAB2(8,2)=7  !2 sens possibles
    TAB2(9,1) =8; TAB2(9,2)=7  !pas logique avec le sens de parcours des faces
    TAB2(10,1) =4; TAB2(10,2)=8 !2 sens possibles
    TAB2(11,1) =5; TAB2(11,2)=8 !pas logique avec le sens de parcours des faces

    do i=1,sem%n_elem !on parcourt tous les elements
        do j=0,11   !on parcourt toutes les aretes de l'element i
            iref = sem%Near_Edges(j,i)
            !        print *,'iref=',iref
            do ii=1,sem%n_elem !on parcourt tous les elements
                do ij=0,11 !on parcourt toutes les aretes de l'element ii
                    if(sem%Near_Edges(ij,ii) == iref) then !on cherche la premiere occurrence de l'arete iref
                        !                 print *,'trouve',ij,ii
                        normeP = sqrt(P1(1)**2 + P1(2)**2 + P1(3)**2)
                        P1(1) = sem%Coord_nodes(1,sem%control_nodes(TAB2(ij,1),ii))
                        P1(2) = sem%Coord_nodes(2,sem%control_nodes(TAB2(ij,1),ii))
                        P1(3) = sem%Coord_nodes(3,sem%control_nodes(TAB2(ij,1),ii))
                        Q1(1) = sem%Coord_nodes(1,sem%control_nodes(TAB2(j,1),i))
                        Q1(2) = sem%Coord_nodes(2,sem%control_nodes(TAB2(j,1),i))
                        Q1(3) = sem%Coord_nodes(3,sem%control_nodes(TAB2(j,1),i))
                        Q2(1) = sem%Coord_nodes(1,sem%control_nodes(TAB2(j,2),i))
                        Q2(2) = sem%Coord_nodes(2,sem%control_nodes(TAB2(j,2),i))
                        Q2(3) = sem%Coord_nodes(3,sem%control_nodes(TAB2(j,2),i))
                        if((abs(P1(1)-Q1(1))<=epsilon*normeP) .AND. (abs(P1(2)-Q1(2))<=epsilon*normeP) .AND. (abs(P1(3)-Q1(3))<=epsilon*normeP)) then
                            sem%Orient_edges(j,i) = 0
                        elseif((abs(P1(1)-Q2(1))<=epsilon*normeP) .AND. (abs(P1(2)-Q2(2))<=epsilon*normeP) .AND. (abs(P1(3)-Q2(3))<=epsilon*normeP)) then
                            sem%Orient_edges(j,i) = 1
                        endif
                        exit
                    endif
                enddo

            enddo
        enddo

    enddo

    PRint*,'Fin Orient_Edge'

    do i=1,sem%n_edge
        WRITE(75,*) 'Arete', i
        WRITE(75,'(A17,2(1X,I5))') 'Numero des sommets', sem%edge_vertex(0:1,i)
        WRITE(75,'(A11,6(1X,E12.5))') 'Coord', (sem%Coord_nodes(1:3,sem%edge_vertex(j,i)),j=0,1)
        WRITE(75,*) '-------------------------'
        WRITE(75,*)
    enddo
    ! recuperation de la numerotation globale des noeuds
    allocate (sem%glob_numbering(sem%n_vertex))


    if (sem%typeElem==8) then
        ! les numeros globaux st forcement en ordre croissant
        do i=1,sem%n_glob_nodes
            sem%glob_numbering(i)=i
        enddo
    endif
    !! FIN MISE EN COMMENTAIRE PASSAGE PAS FACILE

!!! MISE EN COMMENTAIRE SI sem%typeElem==20 (8 en 2D)
!!!  if (sem%typeElem==8) then
!!!     do i=1,sem%n_glob_nodes
!!!        if  (global2local(i).ne.-1) then
!!!           sem%glob_numbering(global2local(i))=i
!!!        endif
!!!     enddo
!!!  endif
!!! FIN MISE EN COMMENTAIRE SI sem%typeElem==20 (8 en 2D)

end subroutine convertir
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!

!>
!! \fn function vertexDejaPresent(tab,dimTab,vertex)
!! \brief
!!
!! \param integer dimTab
!! \param integer tab(dimTab)
!! \param logical vertexDejaPresent
!<


function vertexDejaPresent(tab,dimTab,vertex)

    integer :: i,dimTab,vertex
    integer :: tab(dimTab)
    logical :: vertexDejaPresent

    vertexDejaPresent = .false.

    do i=1,dimTab

        if (tab(i).eq.vertex) then
            vertexDejaPresent = .true.
            print*,"vertex deja present", vertex
            exit
        endif

    enddo


end function vertexDejaPresent


!---------------------------------
!---------------------------------
subroutine calcorient(P, Q, index0)
    implicit none
    real, intent(in) :: P(4,3), Q(4,3)
    integer, intent(out) :: index0(2)
    integer i, j
    real, parameter :: epsilon = 1e-7
    real normeP
    ! INPUT: P, Q
    ! P contient trois coordonnees de 4 points formant une face pour la configuration de reference
    ! Q contient trois coordonnees de 4 points formant une face pour la configuration d'etude
    ! OUTPUT: index0 est un tableau qui renvoie le numero local des 2 1ers points de la configuration de reference dans la configuration d'etude
    ! Objectif: permet de reperer l'orientation
    !Les faces sont:
    ! 1,2,3,4 - face 0
    ! 1,2,6,5 - face 1
    ! 2,3,7,6 - face 2
    ! 4,3,7,8 - face 3
    ! 1,4,8,5 - face 4
    ! 5,6,7,8 - face 5
    do i=1,2
        normeP = sqrt(P(i,1)**2 + P(i,2)**2 + P(i,3)**2)
        do j=1,4
            if((abs(P(i,1)-Q(j,1))<=epsilon*normeP) .AND. (abs(P(i,2)-Q(j,2))<=epsilon*normeP) .AND. (abs(P(i,3)-Q(j,3))<=epsilon*normeP)) then
                !Le point i de la face de reference est le meme que le point j de la face d'etude
                !2 sommets de la face de reference suffisent pour determiner la disposition des sommets de la face d'etude
                index0(i) = j
                exit
            endif
        enddo
    enddo

end subroutine calcorient
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
