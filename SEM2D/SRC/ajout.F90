!>
!!\file ajout.F90
!!\brief Gestion du fichier AJOUT.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module mailleFictive

    type :: Tajout

       integer :: fileId
       integer :: n_elem, n_noeud, n_sommet

       !!modif Gsa allocatable en pointer a cause de gfortran
       !! integer, dimension (:), allocatable :: noeudIndice
       !! integer, dimension (:), allocatable :: elemIndice
       !!
       !! integer, dimension (:,:), allocatable :: elemConnec
       !! real   , dimension (:,:), allocatable :: noeudCoord

       integer, dimension (:), pointer :: noeudIndice
       integer, dimension (:), pointer :: elemIndice

       integer, dimension (:,:), pointer :: elemConnec
       real   , dimension (:,:), pointer :: noeudCoord

       character (len=110) :: fichier

    end type Tajout

    type(Tajout) :: ajout

contains



    !>
    !! \brief
    !!
    !<


    subroutine lireFichierAjout

        integer :: CodeErreur
        logical :: status



        !controle d'existence du fichier
        INQUIRE(File=trim(ajout%fichier),Exist=status)


        if ( .not.status ) then
            write (*,*)"file not found :",trim(ajout%fichier)
            stop
        endif

        open(UNIT=ajout%fileId,IOSTAT=CodeErreur,FILE=trim(ajout%fichier),FORM='formatted',STATUS='old',ACTION='read')
        if (CodeErreur .ne.0 ) print*,'Ouverture du fichier AJOUT  :CodeErreur=',CodeErreur

        ! lecture du nombre de noeuds
        read(ajout%fileId,*)ajout%n_noeud

        allocate(ajout%noeudIndice(ajout%n_noeud))
        allocate(ajout%noeudCoord(3,ajout%n_noeud))


        ! lecture du numero de chaque noeud et leurs coordonnees
        do i=1,ajout%n_noeud
            read(ajout%fileId,*) ajout%noeudIndice(i),ajout%noeudCoord(1,i),ajout%noeudCoord(2,i),ajout%noeudCoord(3,i)
        enddo

        ! lecture du nombre d'elements
        read(ajout%fileId,*)ajout%n_elem


        ajout%n_sommet=4 ! on ne traite que les quadrangles

        allocate (ajout%elemIndice(ajout%n_elem))
        allocate (ajout%elemConnec(ajout%n_sommet,ajout%n_elem))

        ! materiau unique
        numMat=1

        ! lecture du numero de chaque element et de ses connectivites
        do i=1,ajout%n_elem
            read(ajout%fileId,*) ajout%elemIndice(i),ajout%elemConnec(1,i),ajout%elemConnec(2,i),ajout%elemConnec(3,i),ajout%elemConnec(4,i)
        enddo



    end subroutine lireFichierAjout



end module mailleFictive
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
