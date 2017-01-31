!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file unv.F90
!!\brief Contient le module unv.
!!\author
!!\version 1.0
!!\date 10/03/2009
!! Assure la gestion des donnees rattachees au maillage universel.
!<

module unv

    use semdatafiles

    type Tunv

       integer :: fileId
       integer :: n_elem, n_face, n_noeud,n_sommet,typeElem
       integer ::  n_mat

       !! integer, dimension (:), allocatable :: noeudIndice
       !! integer, dimension (:), allocatable :: elemIndice
       !! integer, dimension (:), allocatable :: RenumeroteNoeudEnSem
       !! integer, dimension (:), allocatable :: RenumeroteElemEnSem
       !! integer, dimension (:), allocatable :: tabNoeudUtile
       !!
       !! integer, dimension (:,:), allocatable :: elemConnec
       !! integer, dimension (:,:), allocatable :: face
       !! real   , dimension (:,:), allocatable :: noeudCoord

       integer, dimension (:), pointer :: noeudIndice
       integer, dimension (:), pointer :: elemIndice
       integer, dimension (:), pointer :: RenumeroteNoeudEnSem
       integer, dimension (:), pointer :: RenumeroteElemEnSem
       integer, dimension (:), pointer :: tabNoeudUtile

       integer, dimension (:,:), pointer :: elemConnec
       integer, dimension (:,:), pointer :: face
       real(fpp),dimension(:,:), pointer :: noeudCoord


       character (len=110) :: fichier
       character(Len=MAX_FILE_SIZE) :: fnamef

    end type Tunv



    type(Tunv) :: sunv



    character(LEN=6),parameter ::repere="-1"
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !

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


    subroutine lireFichierUnv

        integer :: CodeErreur
        logical :: status
        character(LEN=20) :: ligne
        character(Len=MAX_FILE_SIZE) :: fnamef

        !controle d'existence du fichier
        INQUIRE(File=trim(sunv%fichier),Exist=status)

        if ( .not.status ) then
            write (*,*)"file not found :",trim(sunv%fichier)
            stop
        endif

        call semname_unv_fichier(sunv%fichier,fnamef)
        open(UNIT=sunv%fileId,IOSTAT=CodeErreur,FILE=fnamef,FORM='formatted',STATUS='old',ACTION='read')
        if (CodeErreur .ne.0 ) print*,'Ouverture du fichier UNV  :CodeErreur=',CodeErreur


        ligne="vide"

        do while (trim(ligne).ne.repere)

            ligne=repere
            read(sunv%fileId,*,end=100) ligne

            if (trim(ligne).ne.repere) cycle

            do while (trim(ligne).eq.repere)
                read(sunv%fileId,*,END=100) ligne
            enddo


            if (trim(ligne).ne.repere) then
                select case (trim(ligne))
                case ("15","781","2411"); call lire_noeuds(sunv%fileId)
                case ("780","2412");    call lire_mailles(sunv%fileId)
                end select
            endif

        enddo

100     continue ! fin du fichier

        call ValidationNoeud


    end subroutine lireFichierUnv
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !
    !>
    !! \brief
    !!
    !! \param integer,intent(IN) fid
    !<


    subroutine lire_noeuds (fid)

        integer,intent(IN) :: fid
        integer :: i,valeurLue
        character(LEN=10) ligne
        integer :: plusGrandeValeur,plusPetiteValeur


        sunv%n_noeud = 0
        read (fid,*)valeurLue
        do while (valeurLue.ne.-1)
            sunv%n_noeud = sunv%n_noeud+1
            read(fid,*) ligne
            read (fid,*)valeurLue
        enddo

        allocate (sunv%noeudIndice(sunv%n_noeud))
        allocate (sunv%noeudCoord(3,sunv%n_noeud))

        do i=1,2*sunv%n_noeud+1
            backspace(fid)
        enddo



        plusPetiteValeur=2**30
        plusGrandeValeur=0

        do i=1,sunv%n_noeud
            read (fid,*) sunv%noeudIndice(i)
            plusGrandeValeur=max(plusGrandeValeur,sunv%noeudIndice(i))
            plusPetiteValeur=min(plusPetiteValeur,sunv%noeudIndice(i))
            !if (sunv%noeudIndice(i).gt.sunv%n_noeud) then
            !write(*,*)"lireFichierUnv : cas non prevu : numero de noeud > nombre total de noeud"
            !stop
            !endif
            read (fid,*) sunv%noeudCoord(1,i),&
                sunv%noeudCoord(2,i),&
                sunv%noeudCoord(3,i)
        enddo



        ! ce tableau contiendra la correspondance de numerotation entre un noeud unv et ce meme noeud sem
        allocate(sunv%RenumeroteNoeudEnSem(plusPetiteValeur:plusGrandeValeur))



    end subroutine lire_noeuds
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !

    !>
    !! \brief
    !!
    !! \param integer,intent(IN) fid
    !<


    subroutine lire_mailles (fid)

        integer,intent(IN) :: fid
        integer :: i,valeurLue
        integer :: plusGrandeValeur,plusPetiteValeur


        sunv%n_elem = 0

        read (fid,*)valeurLue


        sunv%n_sommet = 8 ! quad8



        do while (valeurLue.ne.-1)
            sunv%n_elem = sunv%n_elem+1
            read(fid,*)
            read (fid,*)valeurLue
        enddo

        sunv%n_sommet = 8

        allocate (sunv%elemIndice(sunv%n_elem))
        allocate (sunv%elemConnec(sunv%n_sommet,sunv%n_elem))

        do i=1,2*sunv%n_elem+1
            backspace(fid)
        enddo


        plusPetiteValeur=2**30
        plusGrandeValeur=0

        do i=1,sunv%n_elem
            read (fid,*) sunv%elemIndice(i),valeurLue


            plusGrandeValeur=max(plusGrandeValeur,sunv%elemIndice(i))
            plusPetiteValeur=min(plusPetiteValeur,sunv%elemIndice(i))

            if (sunv%elemIndice(i).gt.sunv%n_elem) then
                !write(*,*)"lireFichierUnv : cas non prevu : numero element > nombre total element"
                !stop
            endif

            select case (valeurLue)
            case (44,54,64,71,84,94) ! elements quad4
                sunv%typeElem=4
                read (fid,*) sunv%elemConnec(1,i),&
                    sunv%elemConnec(2,i),&
                    sunv%elemConnec(3,i),&
                    sunv%elemConnec(4,i)

            case (95) ! elements quad8
                sunv%typeElem=8
                read (fid,*) sunv%elemConnec(1,i),&
                    sunv%elemConnec(2,i),&
                    sunv%elemConnec(3,i),&
                    sunv%elemConnec(4,i),&
                    sunv%elemConnec(5,i),&
                    sunv%elemConnec(6,i),&
                    sunv%elemConnec(7,i),&
                    sunv%elemConnec(8,i)

            case default
                write(*,*)"lireFichierUnv : cas non prevu : descripteur element=",valeurLue
                stop
            end select

            call orienteMaille(i)

        enddo

        ! la numerotation unv peut etre "fantaisiste", celle de sem continue, demarramt  a 0
        allocate(sunv%RenumeroteElemEnSem(plusPetiteValeur:plusGrandeValeur))
        do i=1,sunv%n_elem
            sunv%RenumeroteElemEnSem(sunv%elemIndice(i))=i
        enddo



    end subroutine lire_mailles
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !

    !>
    !! \brief
    !!
    !<


    subroutine ValidationNoeud

        integer :: i,j,k,noeud

        ! il peut y avoir plus de noeuds que necessaire lus dans unv, donc a purger pour sem

        allocate(sunv%tabNoeudUtile(sunv%n_noeud))
        k=0
        do j=1,sunv%n_elem
            do i=1,sunv%typeElem
                noeud=sunv%elemConnec(i,j)
                if (.not.NoeudPresent(noeud,k)) then
                    k=k+1
                    sunv%tabNoeudUtile(k)=noeud
                endif
            enddo
        enddo


        if (k>sunv%n_noeud) then
            write(*,*)"nombre de noeuds insuffisant !"  ! nbre de noeuds definis dans les connectivite > nbre de noeuds definis par leurs coordonnees
        elseif (k<sunv%n_noeud) then
            write(*,*)"Le maillage unv contient plus de noeuds que necessaire => suppression des noeuds inutiles pour sem"
            call purgeNoeud(k)
        else
            write(*,*)"Le nombre de noeuds du maillage unv est coherent avec les connectivites de mailles"
        endif

        deallocate(sunv%tabNoeudUtile)


        ! correspondance entre la numerotation unv et la numerotation sem
        do i=1,sunv%n_noeud
            sunv%RenumeroteNoeudEnSem(sunv%noeudIndice(i))=i
        enddo




        ! les connectivites sont mises en numerotation sem
        do j=1,sunv%n_elem
            do i=1,sunv%typeElem
                noeud=sunv%elemConnec(i,j)
                sunv%elemConnec(i,j)=sunv%RenumeroteNoeudEnSem(noeud)
            enddo
        enddo



    end subroutine ValidationNoeud


    !>
    !! \brief
    !!
    !! \param integer dim
    !<


    subroutine purgeNoeud(dim)

        integer :: dim
        integer, dimension (:,:),allocatable :: dummy

        integer i,k,noeud



        allocate(dummy(dim,0:3))

        k=0
        do i=1,sunv%n_noeud
            noeud=sunv%noeudIndice(i)
            if (NoeudPresent(noeud,dim)) then
                k=k+1
                dummy(k,0)=noeud
                dummy(k,1)=sunv%noeudCoord(1,i)
                dummy(k,2)=sunv%noeudCoord(2,i)
                dummy(k,3)=sunv%noeudCoord(3,i)
            endif
        enddo

        deallocate(sunv%noeudIndice)
        deallocate(sunv%noeudCoord)

        sunv%n_noeud=dim
        allocate(sunv%noeudIndice(dim))
        allocate(sunv%noeudCoord(3,dim))

        do i=1,sunv%n_noeud
            sunv%noeudIndice(i)=dummy(i,0)
            sunv%noeudCoord(1,i)=dummy(i,1)
            sunv%noeudCoord(2,i)=dummy(i,2)
            sunv%noeudCoord(3,i)=dummy(i,3)
        enddo



        deallocate(dummy)


    end subroutine purgeNoeud
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !

    !>
    !! \fn function NoeudPresent(noeud,dim)
    !! \brief
    !!
    !! \param integer noeud
    !! \param integer dim
    !<


    function NoeudPresent(noeud,dim)

        integer :: noeud,dim,k
        logical :: NoeudPresent

        NoeudPresent=.false.
        do k=1,dim
            if (sunv%tabNoeudUtile(k)==noeud) then
                NoeudPresent=.true.
                exit
            endif
        enddo

    end function NoeudPresent
    !
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !

    subroutine orienteMaille(i)


        ! orientation de la maille dans le sens inverse des aiguilles d'une montre

        implicit none

        integer :: i ! numero de la maille a orienter
        integer :: s1, s2, s3,s4,s5,s6,s7,s8 ! numero des sommets de la maille
        real(fpp) :: z

        ! calcul de la composante z du produit vectoriel s1s2^s1s4
        s1=sunv%elemConnec(1,sunv%elemIndice(i))
        s2=sunv%elemConnec(2,sunv%elemIndice(i))
        s3=sunv%elemConnec(3,sunv%elemIndice(i))
        s4=sunv%elemConnec(4,sunv%elemIndice(i))
        s5=sunv%elemConnec(5,sunv%elemIndice(i))
        s6=sunv%elemConnec(6,sunv%elemIndice(i))
        s7=sunv%elemConnec(7,sunv%elemIndice(i))
        s8=sunv%elemConnec(8,sunv%elemIndice(i))


        z=(sunv%noeudCoord(1,s2)-sunv%noeudCoord(1,s1))* &
            (sunv%noeudCoord(2,s4)-sunv%noeudCoord(2,s1))- &
            (sunv%noeudCoord(2,s2)-sunv%noeudCoord(2,s1))* &
            (sunv%noeudCoord(1,s4)-sunv%noeudCoord(1,s1))

        if (z.lt.0) then
            write(*,*)"reorientation de la maille i"

            if (sunv%typeElem.eq.4) then
                sunv%elemConnec(1,sunv%elemIndice(i))=s1
                sunv%elemConnec(2,sunv%elemIndice(i))=s4
                sunv%elemConnec(3,sunv%elemIndice(i))=s3
                sunv%elemConnec(4,sunv%elemIndice(i))=s2

            elseif (sunv%typeElem.eq.8) then
                sunv%elemConnec(1,sunv%elemIndice(i))=s1
                sunv%elemConnec(2,sunv%elemIndice(i))=s8
                sunv%elemConnec(3,sunv%elemIndice(i))=s7
                sunv%elemConnec(4,sunv%elemIndice(i))=s6
                sunv%elemConnec(5,sunv%elemIndice(i))=s5
                sunv%elemConnec(6,sunv%elemIndice(i))=s4
                sunv%elemConnec(7,sunv%elemIndice(i))=s3
                sunv%elemConnec(8,sunv%elemIndice(i))=s2
            endif

        endif


    end subroutine orienteMaille


    !>
    !! \brief
    !!
    !>
    !! \brief
    !!
    !<


    subroutine lister_elements_v1 ()

        implicit none

        type contenu
           integer :: numElem
           type(contenu), pointer :: suivant
        end type contenu

        type pointeur
           type(contenu),pointer :: pt
        end type pointeur

        type(pointeur),dimension(:),allocatable :: listeElement
        type(contenu), pointer :: courant

        integer :: i,iel,inoeud


        allocate (listeElement(sunv%n_noeud))

        ! le pointeur liste%noeud est nullifie
        do i=1,sunv%n_noeud
            nullify (listeElement(i)%pt) ! nullify du pointeur associe a chaque liste d element
        enddo


        do iel=1,sunv%n_elem
            do i=1,4
                inoeud=sunv%elemConnec(i,iel) ! recuperation du numero de noeud dans le tableau des connectivites
                allocate(courant)    ! creation d un nouveau element de type contenu
                courant%numElem=iel  ! courant recele le numero d element
                courant%suivant=>listeElement(inoeud)%pt ! et pointe sur le pointeur de listeElement
                listeElement(inoeud)%pt=>courant ! le pointeur de listeElement pointe sur le nouvel element cree
            enddo
        enddo


        do i=1,sunv%n_noeud
            write(12,*)"Le noeud",i
            do
                if (.not.associated(listeElement(i)%pt)) exit
                write(12,*)"     contient element=",listeElement(i)%pt%numElem
                listeElement(i)%pt=>listeElement(i)%pt%suivant
            enddo
        enddo





    end subroutine lister_elements_v1



    !>
    !! \brief
    !!
    !<


    subroutine lister_faces ()

        implicit none

        integer :: i,iel,n1,n2,n3,n4

        ! tableau de correspondance num local num globale des faces pour chaque element
        allocate(sunv%face(0:3,sunv%n_elem))


        !nombre de faces initialement nul
        sunv%n_face = 0

        !liste initialement vide
        !nullify(listeFace)


        do i=1,sunv%n_elem
            iel=sunv%elemIndice(i) ! recuperation du numero de l element
            n1 = sunv%elemConnec(1,iel)
            n2 = sunv%elemConnec(2,iel)
            n3 = sunv%elemConnec(3,iel)
            n4 = sunv%elemConnec(4,iel)

            !call updateListeFace(n1,n2,iel,0)
            !call updateListeFace(n2,n3,iel,1)
            !call updateListeFace(n3,n4,iel,2)
            !call updateListeFace(n4,n1,iel,3)

        enddo


    end subroutine lister_faces




end module unv

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
