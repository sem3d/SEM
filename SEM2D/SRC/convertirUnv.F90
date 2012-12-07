!>
!!\file convertirUnv.F90
!!\brief Convertit un maillage universel en maillage ipg2.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type (domain), intent (INOUT) Tdomain
!<


subroutine convertirUnv(tDomain)


    ! convertit un maillage universel en maillage ipg2
    !

    use sdomain
    use unv
    use ipg2
    use mailleFictive



    implicit none
    type (domain), intent (INOUT) :: Tdomain


    ! Initialisation des fid
    sunv%fileId = 10
    sem%fileId  = 11
    ajout%fileId = 12

    ! Récupération du nom du fichier unv
    write (sunv%fichier, "(a,a5,I4.4)") trim(adjustl(Tdomain%mesh_file)), ".unv.", Tdomain%Mpi_var%my_rank


    ! Récupération du nom du fichier sem
    write (sem%fichier, "(a,a1,I2.2)") trim(adjustl(Tdomain%mesh_file)), ".", Tdomain%Mpi_var%my_rank




    call lireFichierUnv

    call convertir


    call ecrireFichierSem

    if (Tdomain%Mpi_var%n_proc.gt.1) call ecrireFichierSemFrontiere

end subroutine convertirUnv
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!
!>
!! \brief
!!
!<


subroutine convertir
    use unv
    use ipg2
    use mailleFictive

    implicit none

    type(Tface),pointer :: chercheFace
    type(TelementDeFace),pointer :: chercheElement

    integer i,j,indice,iNumLocal
    integer, dimension(:), allocatable :: global2local


    ! dimension du pb : 2D en dur, a voir
    sem%n_dime = 2

    ! type element : 4 pour quad4 ou 8 pour quad8
    sem%typeElem = sunv%typeElem


    ! recuperation du nombre de noeuds
    sem%n_glob_nodes = sunv%n_noeud
    !sem%n_glob_nodes = sunv%n_noeud+ajout%n_noeud


    ! recuperation des coordonnees 2D des noeuds
    allocate(sem%Coord_nodes(2,sem%n_glob_nodes))
    ! ici, passage du 3D de UNV a du 2D, on suppose que Z est constant
    ! a verifier !
    do i=1,sunv%n_noeud
        sem%Coord_nodes(1,i)=sunv%noeudCoord(1,i)
        sem%Coord_nodes(2,i)=sunv%noeudCoord(2,i)
    enddo

    !  ! ajout des noeuds du fichier ajout
    !  do i=1,ajout%n_noeud
    !    numGlobal=ajout%noeudIndice(i)
    !    sem%Coord_nodes(1,numGlobal)=ajout%noeudCoord(1,i)
    !    sem%Coord_nodes(2,numGlobal)=ajout%noeudCoord(2,i)
    !  enddo


    ! recuperation du nombre de materiaux
    sem%n_mat = 1 ! pour l'instant


    ! recuperation du nombre d'elements
    !sem%n_elem = sunv%n_elem+ajout%n_elem
    sem%n_elem = sunv%n_elem


    ! recuperation du nombre de sommets par element
    sem%n_nodes = sunv%typeElem


    ! recuperation des connectivites noeuds elements
    ! a priori, pas de double numerotation sur les sommets
    ! dans le fichier UNV
    allocate(sem%control_nodes(sem%n_nodes,sem%n_elem))
    !allocate(sem%near_vertex(sem%n_nodes,sem%n_elem))
    allocate(sem%near_vertex(4,sem%n_elem))

    allocate (global2local(sem%n_glob_nodes))
    global2local(:)=-1 ! par defaut un faux numero local
    iNumLocal = 1 ! numero du premier vertex

    do j=1,sem%n_elem

        if (sem%typeElem==4) then ! les sommets st dans le meme ordre : 1 2 3 4
            do i=1,4
                sem%control_nodes(i,j) = sunv%elemConnec(i,j)
                sem%near_vertex(i,j)   = sunv%elemConnec(i,j)
            enddo

        elseif (sem%typeElem==8) then  ! les sommets ne sont pas dans le meme ordre :
            ! 1 2 3 4 5 6 7 8 pour unv
            ! 1 3 5 7 2 4 6 8 pour sem2d


            sem%control_nodes(1,j)=sunv%elemConnec(1,j)
            sem%control_nodes(2,j)=sunv%elemConnec(3,j)
            sem%control_nodes(3,j)=sunv%elemConnec(5,j)
            sem%control_nodes(4,j)=sunv%elemConnec(7,j)
            sem%control_nodes(5,j)=sunv%elemConnec(2,j)
            sem%control_nodes(6,j)=sunv%elemConnec(4,j)
            sem%control_nodes(7,j)=sunv%elemConnec(6,j)
            sem%control_nodes(8,j)=sunv%elemConnec(8,j)

            if (global2local(sem%control_nodes(1,j)).eq.-1) then
                global2local(sem%control_nodes(1,j)) = iNumLocal
                iNumLocal=iNumLocal+1
            endif
            sem%near_vertex(1,j) = global2local(sem%control_nodes(1,j))


            if (global2local(sem%control_nodes(2,j)).eq.-1) then
                global2local(sem%control_nodes(2,j)) = iNumLocal
                iNumLocal=iNumLocal+1
            endif
            sem%near_vertex(2,j) = global2local(sem%control_nodes(2,j))

            if (global2local(sem%control_nodes(3,j)).eq.-1) then
                global2local(sem%control_nodes(3,j)) = iNumLocal
                iNumLocal=iNumLocal+1
            endif
            sem%near_vertex(3,j) = global2local(sem%control_nodes(3,j))



            if (global2local(sem%control_nodes(4,j)).eq.-1) then
                global2local(sem%control_nodes(4,j)) = iNumLocal
                iNumLocal=iNumLocal+1
            endif
            sem%near_vertex(4,j) = global2local(sem%control_nodes(4,j))

        endif

    enddo

    ! nombre de vertex
    if (sem%typeElem==4) sem%n_vertex=sem%n_glob_nodes
    if (sem%typeElem==8) sem%n_vertex=iNumLocal-1

    ! ajout des connectivites du maillage ajout
    !  do i=1,sem%n_nodes
    !  do j=1,ajout%n_elem
    !    numGlobal=ajout%elemIndice(j)
    !    sem%control_nodes(i,numGlobal)=ajout%elemConnec(i,j)
    !    sem%near_vertex(i,numGlobal)=ajout%elemConnec(i,j)
    !  enddo
    !  enddo



    ! recuperation des numeros globaux des faces
    !allocate(sem%Near_Face(0:sem%n_nodes-1,sem%n_elem))
    !sem%Near_Face(:,:)=sunv%face(:,:)
    ! creer les faces et leurs numerotations (locale et globale)
    call lister_faces_sem


    ! recuperation du numero de materiau associe a chaque element


    ! recuperation du nombre de faces
    !sem%n_face = sunv%n_face

    ! recuperation des numeros globaux des 2 elements constituant chaque face
    ! recuperation des numeros locaux des 2 elements constituant chaque face
    ! recuperation des numeros globaux des 2 sommets constituant chaque face
    allocate(sem%face_Element(0:1,sem%n_face))
    allocate(sem%which_face(0:1,sem%n_face))
    allocate(sem%face_vertex(0:1,sem%n_face))

    sem%face_Element(:,:) = -1
    sem%which_face(:,:) = -1

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
            sem%which_face(indice,chercheFace%globnum)=chercheElement%numlocalFace

            if (.not.associated(chercheElement%suivant)) exit !  fin de parcours de la liste des elements
            chercheElement => chercheElement%suivant
            indice = indice+1
        enddo

        sem%face_vertex(0,chercheFace%globnum) = chercheFace%n1
        sem%face_vertex(1,chercheFace%globnum) = chercheFace%n2

        if (.not.associated(chercheFace%suivant)) exit !  fin de parcours de la liste des faces
        chercheFace => chercheFace%suivant
    enddo


    ! recuperation de la numerotation globale des noeuds
    allocate (sem%glob_numbering(sem%n_vertex))


    if (sem%typeElem==4) then
        ! les numeros globaux st forcement en ordre croissant
        do i=1,sem%n_glob_nodes
            sem%glob_numbering(i)=i
        enddo
    endif

    if (sem%typeElem==8) then
        do i=1,sem%n_glob_nodes
            if  (global2local(i).ne.-1) then
                sem%glob_numbering(global2local(i))=i
            endif
        enddo
    endif


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
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
