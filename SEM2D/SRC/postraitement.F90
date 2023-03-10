!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file postraitement.F90
!!\brief Assure le post traitement des donnees.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module Postraitement

    use sdomain
    use semdatafiles
    use mpi
    implicit none

    public  :: POST_Ensight
    private :: POST_EnsightGeo,POST_Ensight1Geo,POST_EnsightCase
    private :: POST_FichierGEOBinaire ,POST_FichierGEOAscii, POST_FichierGEOEntete


    ! Longueur de chaine de caractere
    integer,parameter :: GLO_MAXSTRING = 110

    ! Nom du fichier de geometrie Ensight
    character(Len=MAX_FILE_SIZE) :: POST_NomFichierGeo
    character(Len=MAX_FILE_SIZE) :: POST_NomFichierCase

    ! numero logique du fichier geo
    integer :: POST_IdFichierGeo,POST_IdFichierCase


    ! format d'ecriture des fichiers ensight (ascii:0 ou binaire:1)
    integer, parameter :: ENS_Binaire=1, ENS_Ascii=0
    integer :: ENS_format

    ! code erreur pour allocation des tableaux
    integer :: CodeErreur

    type tPART
       integer :: NbNodeId, PremNodeId
       integer :: NbElemId, PremElemId
       integer :: NbConnec
       integer :: Number

       character(GLO_MAXSTRING) :: Titre
       character(GLO_MAXSTRING) :: TypeElem

       real(fpp), dimension(:,:) , pointer :: Coord

       integer, dimension(:,:), pointer :: Connec

    end type tPART


contains

    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    !>
    !! \brief
    !!
    !! \param type (Domain), intent (IN) Tdomain
    !<


    subroutine POST_Ensight(Tdomain)

        implicit none
        type (Domain), intent (IN) :: Tdomain
        integer :: ierr

        ENS_Format = 0 ! type du fichier .geo : ascii

        POST_IdFichierGeo=200+Tdomain%MPI_var%my_rank
        POST_IdFichierCase=100+Tdomain%MPI_var%my_rank

        call semname_posttraitement_case(Tdomain%MPI_var%my_rank,POST_NomFichierCase)
        call semname_posttraitement_geo(Tdomain%MPI_var%my_rank,POST_NomFichierGeo)

        call POST_EnsightCase()
        call POST_EnsightGeo(Tdomain)

        ! important, chaque proc doit avoir fini sa tache avant d'entamer POST_Ensight1Geo
        call MPI_Barrier (Tdomain%communicateur, ierr)
        if (Tdomain%MPI_var%my_rank == 0) call POST_Ensight1Geo(Tdomain)

    end subroutine POST_Ensight


    subroutine POST_EnsightCase ()


        character(GLO_MAXSTRING) :: Cts, Cfs, Cns, Cfi,Ctv
        integer :: idF, ts, ns, fs, fi,tv




        ! temps initial
        ts=1
        write(Cts,*) ts

        ! nombre de pas
        ns=1
        write(Cns,*) ns

        ! index initial de fichier
        fs=0
        write(Cfs,*) fs

        ! increment de fichier
        fi=1
        write(Cfi,*) fi

        ! Id du fichier case
        idF=POST_IdFichierCase

        ! Ouverture du fichier case a constituer
        open(UNIT=idF,FILE=POST_NomFichierCase,FORM='formatted',STATUS='replace',ACTION='write')



        write(idF,'(A)') 'FORMAT'
        write(idF,'(A)') 'type:     ensight gold'
        write(idF,'(A)') ''
        write(idF,'(A)') 'GEOMETRY'
        write(idF,'(A)') 'model: 1 '//trim(POST_NomFichierGeo)
        write(idF,'(A)') ''
        !write(idF,'(A)') 'VARIABLE'
        !write(idF,'(A)') 'vector per node:  1    Deplacement        deplacement-****.depla'
        !write(idF,'(A)') 'vector per element:  1    Vitesse         vitesse-****.vitesse'
        !write(idF,'(A)') 'scalar per element:  1    AltMaxEau       hauteurMaxEau-****.milieu'
        !write(idF,'(A)') 'scalar per element:  1    Nature          nature-****.milieu'
        !write(idF,'(A)') 'scalar per element:  1    Cote            cote.milieu'
        !write(idF,'(A)') 'scalar per element:  1    Submersion      submersion-****.milieu'
        write(idF,'(A)') ''
        write(idF,'(A)') 'TIME'
        write(idF,'(A)') 'time set:             '//trim(Cts)
        write(idF,'(A)') 'number of steps:      '//trim(Cns)
        write(idF,'(A)') 'filename start number:'//trim(Cfs)
        write(idF,'(A)') 'filename increment:   '//trim(Cfi)

        write(idF,'(A)') 'time values:'
        tv=0
        write(Ctv,*) tv
        write(idF,*) Ctv


        close (idF)


    end subroutine POST_EnsightCase


    !>
    !! \brief
    !!
    !! \param type (Domain), intent (IN) Tdomain
    !<


    subroutine POST_Ensight1Geo(Tdomain)

        implicit none
        type (Domain), intent (IN) :: Tdomain
        integer :: i,j,k,NbreNoeuds,NbreElem,fid,typeElem,DernierNoeud,DernierElem
        character*50 :: NomDeLaPart
        character*5 :: ctypeElem

        integer,dimension(:),allocatable   :: indexNoeud,indexElem
        integer,dimension(:,:),allocatable :: connec
        real(fpp),dimension(:),allocatable :: coordX,coordY,coordZ


        ! ouverture du fichier geo du proc0
        i=0
        fid=200

        call semname_posttraitement_geo(i,POST_NomFichierGeo)

        if (ENS_Format .eq. ENS_Ascii) then
            open(UNIT=fid,FILE=POST_NomFichierGeo,FORM='formatted',ACTION='write',position='append')
        else
            open(UNIT=fid,FILE=POST_NomFichierGeo,FORM='unformatted',ACTION='write',position='append')
        endif

        DernierNoeud = Tdomain%n_glob_nodes
        DernierElem  = Tdomain%n_elem

        do i=1,Tdomain%Mpi_var%n_proc-1

            ! lecture du fichier ensight.geo du proc i
            POST_IdFichierGeo=200+i
            call semname_posttraitement_geo(i,POST_NomFichierGeo)

            if (ENS_Format .eq. ENS_Ascii) then
                open(UNIT=POST_IdFichierGeo,FILE=POST_NomFichierGeo,FORM='formatted',ACTION='read')
            else
                open(UNIT=POST_IdFichierGeo,FILE=POST_NomFichierGeo,FORM='unformatted',ACTION='read')
            endif
            do j=1,6
                read(POST_IdFichierGeo,*)
            enddo
            read(POST_IdFichierGeo,*) NomDeLaPart
            read(POST_IdFichierGeo,*)
            read(POST_IdFichierGeo,*) NbreNoeuds
            allocate(indexNoeud(NbreNoeuds),coordX(NbreNoeuds),coordY(NbreNoeuds),coordZ(NbreNoeuds))
            do j=1,NbreNoeuds
                read(POST_IdFichierGeo,*)indexNoeud(j)
            enddo
            do j=1,NbreNoeuds
                read(POST_IdFichierGeo,*)coordX(j)
            enddo
            do j=1,NbreNoeuds
                read(POST_IdFichierGeo,*)coordY(j)
            enddo
            do j=1,NbreNoeuds
                read(POST_IdFichierGeo,*)coordZ(j)
            enddo

            !ecriture sur le fichier ensight.geo0
            write(fid,'(A)') 'part'
            write(fid,'(I10)') 1+i ! ecriture du numero de la part
            write(fid,'(A,I4.4)') trim(NomDeLaPart),i
            write(fid,'(A)') 'coordinates'
            write(fid,'(I10)') NbreNoeuds
            do j=1,NbreNoeuds
                write(fid,'(I10)') indexNoeud(j)+DernierNoeud
            enddo
            do j=1,NbreNoeuds
                write(fid,'(E12.5)') coordX(j)
            enddo
            do j=1,NbreNoeuds
                write(fid,'(E12.5)') coordY(j)
            enddo
            do j=1,NbreNoeuds
                write(fid,'(E12.5)') coordZ(j)
            enddo

            ! lecture des elem par type
            do
                read(POST_IdFichierGeo,*,end=555)ctypeElem
                if (ctypeElem.eq."quad4") typeElem=4
                if (ctypeElem.eq."quad8") typeElem=8
                read(POST_IdFichierGeo,*) NbreElem
                allocate(indexElem(NbreElem),connec(NbreElem,typeElem))
                do j=1,NbreElem
                    read(POST_IdFichierGeo,*)indexElem(j)
                enddo
                do j=1,NbreElem
                    read(POST_IdFichierGeo,*)(connec(j,k),k=1,typeElem)
                enddo

                !ecriture
                write(fid,'(A)') ctypeElem
                write(fid,'(I10)') NbreElem
                do j=1,NbreElem
                    write(fid,'(I10)') indexElem(j)+DernierElem
                enddo
                do j=1,NbreElem
                    write(fid,'(4I10)') (connec(j,k),k=1,typeElem)


                enddo

                deallocate(indexElem)
                deallocate(connec)
                DernierElem=DernierElem+NbreElem

            enddo  ! fin de boucle sur les differents types d elements
555         continue ! fin de lecture des differents types
            deallocate(indexNoeud)
            deallocate(coordX,coordY,coordZ)

            close(POST_IdFichierGeo)

        enddo ! fin de traitement du ensight geo pour le proc i


        close(fid)

    end subroutine POST_Ensight1Geo




    !>
    !! \brief
    !!
    !! \param type (Domain), intent (IN) Tdomain
    !<


    subroutine POST_EnsightGeo(Tdomain)


        type (Domain), intent (IN) :: Tdomain

        type(tPART) :: Part
        integer :: i,j

        if (Tdomain%n_nodes==4) Part%TypeElem = "quad4"
        if (Tdomain%n_nodes==8) Part%TypeElem = "quad8" ! type d'element (on suppose que le maillage est constitue d'un seul type d'element)

        Part%NbConnec = Tdomain%n_nodes ! nombre de connectivite par element

        Part%NbNodeId = Tdomain%n_glob_nodes ! nombre de noeud constituant le maillage
        Part%NbElemId = Tdomain%n_elem       ! nombre d'element constituant le maillage

        Part%PremNodeId = 1 ! + petit numero de noeud

        Part%PremElemId = 1 ! ! + petit numero d'element

        Part%Number = 1 ! numero de la part

        Part%Titre  = 'part' ! titre de la part



        ! Constitution du tableau des coordonnees pour ecriture dans le fichier geo
        ! allocation du tableau des coordonnees
        allocate (Part%Coord(3,Part%NbNodeId), STAT=CodeErreur)

        if (CodeErreur.ne.0) then
            write(*,*) 'Probleme allocation tableau Part%Coord'
            stop
        endif

        do i=1,2
            do j=1,Part%NbNodeId
                Part%Coord(i,j) = Tdomain%coord_nodes(i-1,j-1)
            enddo
        enddo
        Part%Coord(3,1:Part%NbNodeId)=0.   ! on est en 2D, les z sont mis a zero

        ! Constitution du tableau des connectivites pour ecriture dans le fichier geo
        ! allocation du tableau des connectivites
        allocate (Part%Connec(Part%NbElemId,Part%NbConnec), STAT=CodeErreur)

        if (CodeErreur.ne.0) then
            write(*,*) 'Probleme allocation tableau Part%Connec'
            stop
        endif


        do i=1,Part%NbElemId
            do j=1,Part%NbConnec
                Part%Connec(i,j)= Tdomain%specel(i-1)%control_nodes(j-1)+1
            enddo
        enddo


        ! Ouverture du fichier de geometrie ensight a constituer
        if (ENS_Format .eq. ENS_Ascii) then
            open(UNIT=POST_IdFichierGeo,FILE=POST_NomFichierGeo,FORM='formatted',STATUS='replace',ACTION='write')
        else
            open(UNIT=POST_IdFichierGeo,FILE=POST_NomFichierGeo,FORM='unformatted',STATUS='replace',ACTION='write')
        endif

        call POST_FichierGEOEntete (POST_IdFichierGeo)


        ! Ecriture dans le fichier geo
        if (ENS_Format .eq. ENS_Ascii) call POST_FichierGEOAscii (POST_IdFichierGeo, Part)
        if (ENS_Format .eq. ENS_Binaire) call POST_FichierGEOBinaire (POST_IdFichierGeo, Part)


        close (POST_IdFichierGeo)

    end subroutine POST_EnsightGeo




    !>
    !! \brief
    !!
    !! \param integer, intent(in) id_Fichier
    !<


    subroutine POST_FichierGEOEntete (id_Fichier)
        integer, intent(in) :: id_Fichier
        character*80 :: ligne

        if (ENS_format .eq. ENS_Ascii) then
            write(id_Fichier,'(A)') '#commentaire'
            write(id_Fichier,'(A)') '#commentaire'

            write(id_Fichier,'(A)') 'node id given'
            write(id_Fichier,'(A)') 'element id given'
        else

            ligne = 'Fortran Binary'
            write(id_Fichier) ligne
            ligne = '#commentaire'
            write(id_Fichier) ligne
            write(id_Fichier) ligne
            ligne = 'node id given'
            write(id_Fichier) ligne
            ligne = 'element id given'
            write(id_Fichier) ligne
        endif


    end subroutine POST_FichierGEOEntete


    !>
    !! \brief
    !!
    !! \param integer, intent(in) id_Fichier
    !! \param type(tPART), intent(in) Part
    !<


    subroutine POST_FichierGEOAscii (id_Fichier, Part)


        integer, intent(in)    :: id_Fichier

        type(tPART), intent(in) :: Part


        integer :: inum, jnum


        write(id_Fichier,'(A)') 'part'
        write(id_Fichier,'(I10)') Part%Number

        ! Ecriture du libelle de la PART
        write(id_Fichier,'(A)') Part%Titre
        write(id_Fichier,'(A)') 'coordinates'

        ! Ecriture du nombre d'id de la PART
        write(id_Fichier,'(I10)') Part%NbNodeId

        ! Ecriture des id (=numeros) des noeuds
        do inum=Part%PremNodeId,Part%PremNodeId+Part%NbNodeId-1
            write(id_Fichier,'(I10)') inum
        enddo

        ! Ecriture des coordonnees
        do inum=Part%PremNodeId,Part%PremNodeId+Part%NbNodeId-1
            write(id_Fichier,'(E12.5)') Part%Coord(1,inum)
        enddo
        do inum=Part%PremNodeId,Part%PremNodeId+Part%NbNodeId-1
            write(id_Fichier,'(E12.5)') Part%Coord(2,inum)
        enddo
        do inum=Part%PremNodeId,Part%PremNodeId+Part%NbNodeId-1
            write(id_Fichier,'(E12.5)') Part%Coord(3,inum)
        enddo


        ! Ecriture du type d'element
        write(id_Fichier,'(A)') Part%TypeElem

        ! Ecriture du nombre d'elements
        write(id_Fichier,'(I10)') Part%NbElemId

        ! Ecriture des id (=numeros) des elements
        do inum=Part%PremElemId,Part%PremElemId+Part%NbElemId-1
            write(id_Fichier,'(I10)') inum
        enddo

        ! Ecriture des connectivites pour chaque element
        do inum=1,Part%NbElemId
            do jnum=1,Part%NbConnec -1
                write(UNIT=id_Fichier,ADVANCE='no',FMT='(I10)') Part%Connec(inum,jnum)
            enddo
            jnum = Part%NbConnec
            write(UNIT=id_Fichier,ADVANCE='yes',FMT='(I10)') Part%Connec(inum,jnum)
        enddo


    end subroutine POST_FichierGeoAscii

    !>
    !! \brief
    !!
    !! \param integer, intent(in) id_Fichier
    !! \param type(tPART), intent(in) Part
    !<


    subroutine POST_FichierGEOBinaire (id_Fichier, Part)


        integer, intent(in)    :: id_Fichier

        type(tPART), intent(in) :: Part

        character*80 :: ligne
        real(kind=4), dimension(:,:), allocatable :: TabR4
        integer :: inum, jnum


        ligne = 'part'
        write(id_Fichier) ligne
        write(id_Fichier) Part%Number

        ! Ecriture du libelle de la PART
        ligne = trim(Part%Titre)
        write(id_Fichier) ligne
        ligne = 'coordinates'
        write(id_Fichier) ligne

        ! Ecriture du nombre d'id de la PART
        write(id_Fichier) Part%NbNodeId

        ! Ecriture des id (=numeros) des noeuds
        write(id_Fichier) (inum,inum=Part%PremNodeId,Part%PremNodeId+Part%NbNodeId-1)


        ! passage en real*4 du tableau coord
        allocate (TabR4(2,Part%NbNodeId))
        TabR4 = Part%Coord

        ! Ecriture des coordonnees
        write(UNIT=id_Fichier) (TabR4(1,inum),inum=Part%PremNodeId,Part%PremNodeId+Part%NbNodeId-1)
        write(UNIT=id_Fichier) (TabR4(2,inum),inum=Part%PremNodeId,Part%PremNodeId+Part%NbNodeId-1)
        write(UNIT=id_Fichier) (TabR4(3,inum),inum=Part%PremNodeId,Part%PremNodeId+Part%NbNodeId-1)
        deallocate (TabR4)

        ! Ecriture du type d'element
        ligne = trim(Part%TypeElem)
        write(id_Fichier) ligne

        ! Ecriture du nombre d'elements
        write(id_Fichier) Part%NbElemId

        ! Ecriture des id (=numeros) des elements
        write(id_Fichier) (inum,inum=Part%PremElemId,Part%PremElemId+Part%NbElemId-1)

        ! Ecriture des connectivites pour tous les elements
        write(id_Fichier) ((Part%Connec(inum,jnum),jnum=1,Part%NbConnec),inum=1,Part%NbElemId)


    end subroutine POST_FichierGeoBinaire



end module Postraitement

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
