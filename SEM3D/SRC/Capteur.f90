!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Capteur.f90
!!\brief Permet de manipuler les quantités associées aux capteurs.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module mCapteur

    use sdomain
    use semdatafiles
    use sem_hdf5
    use sem_c_config
    use constants
    use mshape8
    use mshape27
    use msnapshots, only : integrate_on_element
    use lagrange_prop
#include "index.h"
    implicit none

    public :: save_capteur, evalueSortieCapteur, flushAllCapteurs, create_capteurs
    private ::  flushCapteur

    type :: tCapteur
        type(tCapteur),pointer :: suivant ! pour passer au capteur suivant
        integer :: periode          ! frequence de captation des grandeur
        real(fpp), dimension (3) :: Coord  ! localisation du capteur
        character(LEN=20) :: nom      ! nom du capteur
        integer :: n_el ! numero de la maille dans laquelle se trouve le capteur
        ! si le capteur est partage entre plusieurs mailles, une seule suffit (type_calcul=1)
        real(fpp) :: xi, eta, zeta ! abscisses curvilignes pour le capteur en cas d'interpolation (type_calcul=1)
        integer :: numproc               ! numero du proc localisant le capteur
        integer :: icache
        real(fpp), dimension(:,:), allocatable :: valuecache
        integer :: type
        integer :: domain
        integer :: lnum
        !
        integer :: ngll
        real(fpp), dimension(:), allocatable :: outx,outy,outz,doutx,douty,doutz
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: dUdX
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev_pl
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
    end type tCapteur

    integer           :: dimCapteur        ! nombre total de capteurs

    type(Tcapteur), pointer :: listeCapteur
    type(Tcapteur), pointer :: capt_En_PS

    type(Tcapteur), dimension(:), allocatable, target :: localCapteurs
    integer :: nCapteursOnRank
    integer,parameter :: fileIdCapteur=200  ! id fichier capteur

    logical :: traces_h5_created
contains

    subroutine create_capteurs(Tdomain)
        use sem_mpi
        implicit none
        type(domain), intent (inout) :: Tdomain
        !
        type(Tcapteur),pointer :: capteur, oldcapteur
        type(C_PTR) :: station_next;
        type(sem_station), pointer :: station_ptr
        character(Len=MAX_FILE_SIZE) :: nom
        real(fpp) :: xc, yc, zc, xi, eta, zeta
        real(fpp) :: xc0, yc0, zc0
        character(len=MAX_FILE_SIZE) :: fnamef
        integer :: numproc, numproc_max, ierr, n_el, n_eln, i, n_out, ngll
        real(fpp) :: dmin, glob_dmin
        real(fpp),dimension(:),allocatable :: gllc

        real(fpp), dimension(0:2, 0:Tdomain%n_nodes-1) :: coordl
        integer :: periodeRef
        logical :: flag
        integer :: c

        flag = .false.
        station_next = Tdomain%config%stations
        nullify(listeCapteur)
        periodeRef = -1

        ! En reprise on ne recree pas le fichier capteur
        if (Tdomain%TimeD%NtimeMin==0) then
            traces_h5_created = .false.
        else
            traces_h5_created = .true.
        endif

        Tdomain%has_station = .false. !Stations other than Total Energy


        Tdomain%out_var_offset = 0
        do i = 0,size(Tdomain%out_var_capt)-2
            if (Tdomain%out_var_capt(i) == 1) then
                Tdomain%out_var_offset(i+1) = Tdomain%out_var_offset(i) + OUT_VAR_DIMS_3D(i)
            else
                Tdomain%out_var_offset(i+1) = Tdomain%out_var_offset(i)
            end if
        end do
        nCapteursOnRank = 0

        do while (C_ASSOCIATED(station_next))
            call c_f_pointer(station_next, station_ptr)
            xc = station_ptr%coords(1)
            yc = station_ptr%coords(2)
            zc = station_ptr%coords(3)
            numproc = -1
            nom = fromcstr(station_ptr%name)
            !if (trim(nom(1:20))=="04x29") flag = .true.
            call trouve_capteur(Tdomain, xc, yc, zc, n_el, dmin, xi, eta, zeta, flag)
            ! Cas ou le capteur est dans le maillage
            call MPI_AllReduce(dmin, glob_dmin, 1, MPI_REAL_FPP, &
                MPI_MIN, Tdomain%communicateur, ierr)
            if (dmin==glob_dmin) then
                numproc = Tdomain%rank
            end if
            call MPI_AllReduce(numproc, numproc_max, 1, MPI_INTEGER, &
                MPI_MAX, Tdomain%communicateur, ierr)
            ! Cas ou le capteur est completement en dehors du maillage
            if (numproc_max==-1) then
                if (Tdomain%rank==0) then
                    nom = fromcstr(station_ptr%name)
                    write(*,*) "Something wrong happened..."
                    write(*,*) "The station ", trim(nom), " doesn't appear to be on any processor"
                    write(*,*) "Please verify that the station location is within the computation domain"
                end if
                stop 1
            end if


            ! attention si le capteur est partage par plusieurs procs. On choisit le proc de num max
            if(Tdomain%rank==numproc_max) then
                allocate(capteur)
                Tdomain%has_station = .true.
                n_out = Tdomain%nReqOut

                if (glob_dmin>0) then
                    do i = 0, Tdomain%n_nodes-1
                        coordl(0:2, i) = Tdomain%Coord_Nodes(0:2, Tdomain%specel(n_el)%Control_Nodes(i))
                    enddo
                    if (Tdomain%n_nodes==8) then
                        call shape8_local2global(coordl, xi, eta, zeta, xc0, yc0, zc0)
                    else
                        call shape27_local2global(coordl, xi, eta, zeta, xc0, yc0, zc0)
                    end if
                    write(*,*) "The station",trim(nom)," is outside. Moved from ", xc, yc, zc, " to ", xc0, yc0, zc0
                else
                    xc0 = xc
                    yc0 = yc
                    zc0 = zc
                end if
                nom = fromcstr(station_ptr%name)
                ! Initialisation.
                ngll = domain_ngll(Tdomain, Tdomain%specel(n_el)%domain)

                capteur%nom = nom(1:20)     ! ses caracteristiques par defaut
                capteur%type = CPT_INTERP
                capteur%periode = station_ptr%period
                periodeRef = capteur%periode
                capteur%coord(1) = xc0
                capteur%coord(2) = yc0
                capteur%coord(3) = zc0
                capteur%xi = xi
                capteur%eta = eta
                capteur%zeta = zeta
                capteur%n_el = n_el
                capteur%numproc = numproc_max
                capteur%domain = Tdomain%specel(n_el)%domain
                capteur%lnum = Tdomain%specel(n_el)%lnum
                capteur%icache = 1
                capteur%suivant => listeCapteur
                listeCapteur => capteur
                nCapteursOnRank = nCapteursOnRank + 1
                write(*,"(A,A,A,I5,A,I6,A,F8.4,A,F8.4,A,F8.4,A,I1)") "Capteur:", trim(capteur%nom), &
                    " on proc ", Tdomain%rank, " in elem ", n_el, " at ", xi, ",", eta, ",", zeta, &
                    " in domain ", Tdomain%specel(n_el)%domain

                ! si c'est un nouveau run, suppression de l'eventuel fichier de sortie des capteurs
                if (Tdomain%traces_format == 1) then
                    if ( .not.Tdomain%logicD%run_restart) then
                        call semname_capteur_type(capteur%nom,".txt",fnamef)
                        open(123,file=trim(fnamef),status="replace",form="formatted")
                        close(123)
                    end if
                end if
            end if

            station_next = station_ptr%next
        end do

        ! Create array for GPU computing
        if (nCapteursOnRank.gt.0) then
            capteur => listeCapteur
            allocate(localCapteurs(0:nCapteursOnRank-1))
            do c = 0,nCapteursOnRank-1
                localCapteurs(c) = capteur
                if (c.lt.nCapteursOnRank-2) then
                    localCapteurs(c)%suivant => localCapteurs(c+1)
                else
                    nullify(localCapteurs(c)%suivant)
                endif
                n_el = localCapteurs(c)%n_el
                localCapteurs(c)%ngll = domain_ngll(Tdomain, Tdomain%specel(n_el)%domain)

                allocate(localCapteurs(c)%fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
                allocate(localCapteurs(c)%fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                allocate(localCapteurs(c)%fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                allocate(localCapteurs(c)%fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                allocate(localCapteurs(c)%eps_vol(0:ngll-1,0:ngll-1,0:ngll-1))
                allocate(localCapteurs(c)%P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                allocate(localCapteurs(c)%S_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                allocate(localCapteurs(c)%eps_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
                allocate(localCapteurs(c)%dUdX(0:ngll-1,0:ngll-1,0:ngll-1,0:8))
                allocate(localCapteurs(c)%eps_dev_pl(0:ngll-1,0:ngll-1,0:ngll-1,0:6))
                allocate(localCapteurs(c)%sig_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
                allocate(localCapteurs(c)%outx(0:ngll-1))
                allocate(localCapteurs(c)%outy(0:ngll-1))
                allocate(localCapteurs(c)%outz(0:ngll-1))
                allocate(localCapteurs(c)%doutx(0:ngll-1))
                allocate(localCapteurs(c)%douty(0:ngll-1))
                allocate(localCapteurs(c)%doutz(0:ngll-1))
                allocate(localCapteurs(c)%valuecache(1:n_out+1,NCAPT_CACHE))
                call domain_gllc(Tdomain, Tdomain%specel(n_el)%domain, GLLc)

                do i = 0,ngll - 1
                    call  pol_lagrange(ngll,GLLc,i,capteur%xi  ,localCapteurs(c)%outx(i))
                    call  pol_lagrange(ngll,GLLc,i,capteur%eta ,localCapteurs(c)%outy(i))
                    call  pol_lagrange(ngll,GLLc,i,capteur%zeta,localCapteurs(c)%outz(i))
                    call  der_lagrange(ngll,GLLc,i,capteur%xi  ,localCapteurs(c)%doutx(i))
                    call  der_lagrange(ngll,GLLc,i,capteur%eta ,localCapteurs(c)%douty(i))
                    call  der_lagrange(ngll,GLLc,i,capteur%zeta,localCapteurs(c)%doutz(i))
                end do
                oldcapteur => capteur
                capteur => capteur%suivant
                deallocate(oldcapteur)
            enddo
            listeCapteur => localCapteurs(0)
        endif
        !$acc enter data    copyin(localCapteurs)
        do c = 0,nCapteursOnRank-1
            !$acc enter data &
            !$acc&      copyin(localCapteurs(c)) &
            !$acc&      copyin(localCapteurs(c)%fieldP) &
            !$acc&      copyin(localCapteurs(c)%fieldU) &
            !$acc&      copyin(localCapteurs(c)%fieldV) &
            !$acc&      copyin(localCapteurs(c)%fieldA) &
            !$acc&      copyin(localCapteurs(c)%eps_vol) &
            !$acc&      copyin(localCapteurs(c)%P_energy) &
            !$acc&      copyin(localCapteurs(c)%S_energy) &
            !$acc&      copyin(localCapteurs(c)%eps_dev) &
            !$acc&      copyin(localCapteurs(c)%dUdX) &
            !$acc&      copyin(localCapteurs(c)%eps_dev_pl) &
            !$acc&      copyin(localCapteurs(c)%sig_dev) &
            !$acc&      copyin(localCapteurs(c)%valuecache) &
            !$acc&      copyin(localCapteurs(c)%outx) &
            !$acc&      copyin(localCapteurs(c)%outy) &
            !$acc&      copyin(localCapteurs(c)%outz) &
            !$acc&      copyin(localCapteurs(c)%doutx) &
            !$acc&      copyin(localCapteurs(c)%douty) &
            !$acc&      copyin(localCapteurs(c)%doutz) &
            !$acc&
        end do
        if(periodeRef < 1) periodeRef = 1

        ! Energy outputs
        if(Tdomain%out_var_capt(OUT_TOTAL_ENERGY) == 1) then

            if(Tdomain%rank == 0) write(*,*) "CREATING ENERGY SENSORS"

            allocate(capt_En_PS)

            n_out = 5
            if (.not.allocated(capt_En_PS%valuecache)) allocate(capt_En_PS%valuecache(1:n_out+1,NCAPT_CACHE))

            capt_En_PS%nom = "En_PS"
            capt_En_PS%type = CPT_ENERGY
            !capt_En_PS%periode = periodeRef !TODO
            capt_En_PS%periode = 1 !TODO
            capt_En_PS%coord(1) = -1111
            capt_En_PS%coord(2) = -1111
            capt_En_PS%coord(3) = -1111
            capt_En_PS%xi = -1111
            capt_En_PS%eta = -1111
            capt_En_PS%zeta = -1111
            capt_En_PS%n_el = -1
            capt_En_PS%numproc = Tdomain%rank
            capt_En_PS%icache = 1
            capt_En_PS%suivant => listeCapteur
            listeCapteur => capt_En_PS
            write(*,"(A,A,A,I5,A,I6,A,F8.4,A,F8.4,A,F8.4)") "Capteur:", trim(capt_En_PS%nom), &
                " on proc ", Tdomain%rank, " in elem ", n_el, " at ", xi, ",", eta, ",", zeta

        end if

    end subroutine create_capteurs


    subroutine evalueSortieCapteur(it, sortie_capteur)
        implicit none
        integer, intent(in) :: it
        logical, intent(out) :: sortie_capteur
        type(tCapteur),pointer :: capteur

        sortie_capteur = .FALSE.
        ! boucle sur les capteurs

        !write(*,*)  "Before pointer"
        capteur=>listeCapteur
        !write(*,*)  "Before boucle"
        do while (associated(capteur))
            !write(*,*)  "capteur%nom =", capteur%nom
            !write(*,*)  "capteur%periode=", capteur%periode
            if(capteur%periode < 1) stop "ERROR, station with period smaller than 1"
            if (mod(it,capteur%periode)==0) then ! on fait la sortie
                sortie_capteur = .TRUE.
            endif

            !write(*,*)  "Before capteur suivant"
            capteur=>capteur%suivant
            !write(*,*)  "After capteur suivant"
        enddo
    end subroutine evalueSortieCapteur

    !---------------------------------------------------------------------
    !---------------------------------------------------------------------

    !>
    !! \fn subroutine save_capteur(Tdomain)
    !! \brief
    !!
    !! \param type (domain) TDomain
    !<
    subroutine save_capteur(Tdomain, ntime)


        implicit none

        integer :: ntime, c, ngll, icache
        type (domain) :: TDomain

        ! type(tCapteur),pointer :: capteur
        logical :: do_flush
        real(fpp) :: rtime
        rtime = TDomain%timeD%rtime
        do_flush = .false.
        ! boucle sur les capteurs kernel si openacc
        !$acc parallel loop gang async(2) wait(1)
        do c = 0,nCapteursOnRank-1
            if (mod(ntime, localCapteurs(c)%periode)==0) then
                if (localCapteurs(c)%type == CPT_INTERP) then
                    call sortieGrandeurCapteur_interp(Tdomain, rtime, localCapteurs(c)%ngll, localCapteurs(c))
                else if (localCapteurs(c)%type == CPT_ENERGY) then
!                    call sortieGrandeurCapteur_energy(Tdomain, ngll, localCapteurs(c))
                end if
            endif
        enddo
        ! check counters (always on cpu)
        do c = 0,nCapteursOnRank-1
            if (localcapteurs(c)%icache==NCAPT_CACHE) then
                do_flush = .true.
            end if
        enddo
        if (do_flush) then
            call flushAllCapteurs(Tdomain)
        endif
        !! XXX Bug : flushAllCapteurs resets icache to 1
        !! but we increment it there, so it starts at 2 after a flush
        ! update counters (always on cpu)
        do c = 0,nCapteursOnRank-1
            if (mod(ntime, localCapteurs(c)%periode)==0) then
                localcapteurs(c)%icache = localcapteurs(c)%icache + 1
                if (localcapteurs(c)%icache>NCAPT_CACHE) then
                    localcapteurs(c)%icache = 1
                endif
                !$acc update device(localCapteurs(c)%icache) async(2)
            endif
        enddo
    end subroutine save_capteur

    function dset_capteur_name(capteur)
        implicit none
        type(tCapteur),pointer :: capteur
        character(len=40) :: dset_capteur_name
        dset_capteur_name = trim(adjustl(capteur%nom)) !//"_"//trim(adjustl(capteur%grandeur))
    end function dset_capteur_name

    function dset_capteur_posname(capteur)
        implicit none
        type(tCapteur),pointer :: capteur
        character(len=40) :: dset_capteur_posname
        dset_capteur_posname = trim(adjustl(capteur%nom)) //"_pos"
    end function dset_capteur_posname

    subroutine create_traces_h5_skel(Tdomain)
        use HDF5
        implicit none
        type (domain), intent(inout) :: TDomain
        type(tCapteur),pointer :: capteur
        character (len=MAX_FILE_SIZE) :: fnamef
        character (len=40) :: dname
        integer(HID_T) :: fid, dset_id
        integer :: hdferr, n_out

        call init_hdf5()

        call semname_tracefile_h5(Tdomain%rank, fnamef)
        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
        call create_capteur_descriptions(Tdomain, fid)

        n_out = Tdomain%nReqOut+1

        capteur=>listeCapteur
        do while (associated(capteur))
            dname = dset_capteur_name(capteur)
            call create_dset_2d(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, &
                int(n_out,HSIZE_T), int(H5S_UNLIMITED_F,HSIZE_T), dset_id)
            call h5dclose_f(dset_id, hdferr)
            dname = dset_capteur_posname(capteur)
            call write_dataset(fid, trim(adjustl(dname)), capteur%Coord)
            capteur=>capteur%suivant
        enddo

        call h5fclose_f(fid, hdferr)
    end subroutine create_traces_h5_skel

    ! Creates a string dataset describing each columns of the trace file
    subroutine create_capteur_descriptions(Tdomain, fid)
        use HDF5
        use constants, only : OUT_VAR_NAMES, OUT_VAR_DIMS_3D
        type (domain), intent(inout) :: TDomain
        integer(HID_T), intent(in) :: fid
        !
        integer(HID_T) :: tid, dsetid, spaceid
        integer :: hdferr
        character(len=12), dimension(:), allocatable :: varnames
        character(len=12), dimension(6) :: energy_varnames = ["Time       1", &
                                                              "Eng_P      1", &
                                                              "Eng_S      1", &
                                                              "Eng_Resid  1", &
                                                              "Eng_Cine   1", &
                                                              "Eng_Total  1"]
        character(len=12) :: temp
        integer :: d,k,dim,dimtot
        integer(HSIZE_T), dimension(1) :: dims

        dimtot = Tdomain%nReqOut
        allocate(varnames(0:dimtot))
        varnames(0) = "Time       1"
        d = 1
        do k=0,OUT_LAST
            if (Tdomain%out_var_capt(k)==1) then
                if(k == OUT_TOTAL_ENERGY) cycle
                do dim=1,OUT_VAR_DIMS_3D(k)
                    write(temp,"(A,I2)") OUT_VAR_NAMES(k),dim
                    varnames(d) = temp
                    d = d+1
                end do
            end if
        end do
        !
        if(Tdomain%has_station) then
            dims(1) = d
            call H5Tcopy_f(H5T_FORTRAN_S1, tid, hdferr)
            call H5Tset_size_f(tid, 12_HSIZE_T, hdferr)
            call H5Screate_simple_f(1, dims, spaceid, hdferr)
            call H5Dcreate_f(fid, "Variables", tid, spaceid, dsetid, hdferr)
            call H5Dwrite_f(dsetid, tid, varnames, dims, hdferr, spaceid, spaceid)
            call H5Dclose_f(dsetid, hdferr)
            call H5Sclose_f(spaceid, hdferr)
            call H5Tclose_f(tid, hdferr)
        end if
        !
        if(Tdomain%out_var_capt(OUT_TOTAL_ENERGY) == 1) then
            dims(1) = size(energy_varnames)
            call H5Tcopy_f(H5T_FORTRAN_S1, tid, hdferr)
            call H5Tset_size_f(tid, 12_HSIZE_T, hdferr)
            call H5Screate_simple_f(1, dims, spaceid, hdferr)
            call H5Dcreate_f(fid, "En_PS_Variables", tid, spaceid, dsetid, hdferr)
            call H5Dwrite_f(dsetid, tid, energy_varnames, dims, hdferr, spaceid, spaceid)
            call H5Dclose_f(dsetid, hdferr)
            call H5Sclose_f(spaceid, hdferr)
            call H5Tclose_f(tid, hdferr)
        end if
    end subroutine create_capteur_descriptions
    
    subroutine append_traces_h5(Tdomain)
        implicit none
        type (domain), intent(inout) :: TDomain
        type(tCapteur),pointer :: capteur
        character (len=40) :: dname
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: dset_id, fid
        integer :: hdferr

        call semname_tracefile_h5(Tdomain%rank, fnamef)

        call h5fopen_f(fnamef, H5F_ACC_RDWR_F, fid, hdferr)

        capteur=>listeCapteur
        do while (associated(capteur))
            !write(*,*) "Capteur:", capteur%nom
            dname = dset_capteur_name(capteur)
            if (capteur%icache==1) then
                capteur=>capteur%suivant
                cycle
            endif
            call h5dopen_f(fid, trim(dname), dset_id, hdferr)
            call append_dataset_2d(dset_id, capteur%valuecache(:,1:capteur%icache-1), hdferr)
            call h5dclose_f(dset_id, hdferr)
            capteur%icache=1
            capteur=>capteur%suivant
        enddo

        call h5fclose_f(fid, hdferr)
    end subroutine append_traces_h5

    subroutine flushAllCapteurs(Tdomain)
        implicit none
        type (domain), intent(inout) :: TDomain
        type(tCapteur),pointer :: capteur
        integer :: c
        do c = 0,nCapteursOnRank-1
            !$acc update host(localCapteurs(c)%valuecache) async(2)
        enddo
        !$acc wait(2)

        ! Default unspecified value is 'text'
        if (Tdomain%traces_format == 0) Tdomain%traces_format = 1

        if (Tdomain%traces_format == 1) then
            ! boucle sur les capteurs
            capteur=>listeCapteur
            do while (associated(capteur))
                call flushCapteur(capteur)
                capteur=>capteur%suivant
            enddo
        else
            ! Sauvegarde au format hdf5
            if (associated(listeCapteur)) then
                ! On ne fait rien sur ce proc si on n'a pas de capteur
                if (.not. traces_h5_created) then

                    call create_traces_h5_skel(Tdomain)
                    traces_h5_created = .true.
                end if
                call append_traces_h5(Tdomain)
            end if
        end if
    end subroutine flushAllCapteurs

    subroutine flushCapteur(capteur)
        implicit none
        type(tCapteur),pointer :: capteur
        !
        integer, parameter :: fileId=123
        integer :: j
        character(len=MAX_FILE_SIZE) :: fnamef
        character(len=20) :: sizeChar

        if (capteur%icache==1) return

        call semname_capteur_type(capteur%nom,".txt",fnamef)

        open(fileId,file=trim(fnamef),status="unknown",form="formatted",position="append")
        do j=1,capteur%icache-1
            ! start modifs
            write(sizeChar, *) size(capteur%valuecache)
            write(fileId,'('//trim(sizeChar)//'(1X,E16.8E3))') capteur%valuecache(:,j)
            ! end modifs
        end do
        close(fileId)
        capteur%icache = 1
    end subroutine flushCapteur


    !---------------------------------------------------------------------
    !---------------------------------------------------------------------


    !! effectue l'interpolation des grandeurs dans la maille dans laquelle se trouve le capteur
    !! la maille se trouve dans un seul proc
    !! seul le proc gere l'ecriture
    !!
    subroutine sortieGrandeurCapteur_interp(Tdomain, rtime, ngll, capteur)
        !$acc routine worker
        use constants
        use dom_solid
        use dom_solid_dg
        use dom_fluid
        use dom_solidpml
        use dom_fluidpml
        implicit none
        !
        type(domain),intent(inout)    :: TDomain
        type(tCapteur),intent(inout)  :: capteur
        integer, intent(in)           :: ngll
        real(fpp), intent(in)         :: rtime
        !
        integer                       :: i, j, k, ioff, ni, c
        real(fpp)                     :: weight
        ! Evaluation of derivative of lagrange polynomial d/dx at xi, d/dy at eta d/dz at zeta
        logical :: nl_flag
        integer :: nComp

        ! Verification : le capteur est il gere par le proc. ?

        nl_flag = Tdomain%nl_flag
        ! On recupere les variables de l'element associe au capteur.
        select case(capteur%domain)
            case (DM_SOLID_DG)
              !call get_solid_dg_dom_var(Tdomain%sdomdg, capteur%lnum, out_variables, &
              !  fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev, &
              !  dUdX)
            case (DM_SOLID_CG)
              call get_solid_dom_var(Tdomain%sdom, capteur%lnum, Tdomain%out_var_capt, &
                capteur%fieldU, capteur%fieldV, capteur%fieldA, capteur%fieldP, &
                capteur%P_energy, capteur%S_energy, capteur%eps_vol, capteur%eps_dev, capteur%sig_dev, &
                capteur%dUdX, nl_flag, capteur%eps_dev_pl)
            case (DM_FLUID_CG)
              call get_fluid_dom_var(Tdomain%fdom, capteur%lnum, ngll, Tdomain%out_var_capt, &
                capteur%fieldU, capteur%fieldV, capteur%fieldA, capteur%fieldP, &
                capteur%P_energy, capteur%S_energy, capteur%eps_vol, &
                capteur%eps_dev, capteur%sig_dev, capteur%dUdX)
            case (DM_SOLID_CG_PML)
              call get_solidpml_dom_var(Tdomain%spmldom, capteur%lnum, Tdomain%out_var_capt, &
                capteur%fieldU, capteur%fieldV, capteur%fieldA, capteur%fieldP, capteur%P_energy, &
                capteur%S_energy, capteur%eps_vol, capteur%eps_dev, capteur%sig_dev)
            case (DM_FLUID_CG_PML)
              call get_fluidpml_dom_var(Tdomain%fpmldom, capteur%lnum, Tdomain%out_var_capt, &
                capteur%fieldU, capteur%fieldV, capteur%fieldA, capteur%fieldP, &
                capteur%P_energy, capteur%S_energy, capteur%eps_vol, capteur%eps_dev, capteur%sig_dev)
            case default
              stop "unknown domain"
        end select

        ni = capteur%icache
        ! Sauvegarde des valeurs dans le capteur.
        capteur%valuecache(1,ni) = rtime
        capteur%valuecache(2:Tdomain%nReqOut+1,ni) = 0

       ! On interpole le DOF a la position du capteur.
       !$acc loop worker vector collapse(3)
       do i = 0,ngll - 1
           do j = 0,ngll - 1
               do k = 0,ngll - 1
                   weight = capteur%outx(i)*capteur%outy(j)*capteur%outz(k)
                   if (Tdomain%out_var_capt(OUT_DEPLA) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_DEPLA)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_DEPLA)-1
                       do c=ioff,ioff+nComp
                           !$acc atomic update
                           capteur%valuecache(c,ni) &
                               = capteur%valuecache(c,ni) + weight*capteur%fieldU(i,j,k,c-ioff)
                       enddo
                   end if

                   if (Tdomain%out_var_capt(OUT_VITESSE) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_VITESSE)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_VITESSE)-1
                       do c=ioff,ioff+nComp
                           !$acc atomic update
                           capteur%valuecache(c,ni) &
                               = capteur%valuecache(c,ni) + weight*capteur%fieldV(i,j,k,c-ioff)
                           enddo
                   end if

                   if (Tdomain%out_var_capt(OUT_ACCEL) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_ACCEL)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_ACCEL)-1
                       capteur%valuecache(ioff:ioff+nComp,ni) &
                           = capteur%valuecache(ioff:ioff+nComp,ni) + weight*capteur%fieldA(i,j,k,:)
                   end if

                   if (Tdomain%out_var_capt(OUT_PRESSION) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_PRESSION)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_PRESSION)-1
                       capteur%valuecache(ioff+nComp,ni) &
                           = capteur%valuecache(ioff+nComp,ni) + weight*capteur%fieldP(i,j,k)
                   end if

                   if (Tdomain%out_var_capt(OUT_ENERGYP) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_ENERGYP)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_ENERGYP)-1
                       capteur%valuecache(ioff+nComp,ni) = capteur%valuecache(ioff+nComp,ni) + weight*capteur%P_energy(i,j,k)
                   end if

                   if (Tdomain%out_var_capt(OUT_ENERGYS) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_ENERGYS)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_ENERGYS)-1
                       capteur%valuecache(ioff+nComp,ni) = capteur%valuecache(ioff+nComp,ni) + weight*capteur%S_energy(i,j,k)
                   end if

                   if (Tdomain%out_var_capt(OUT_DUDX) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_DUDX)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_DUDX)-1
                       capteur%valuecache(ioff:ioff+nComp,ni) = capteur%valuecache(ioff:ioff+nComp,ni)+weight*capteur%dUdX(i,j,k,:)
                   end if

                   if (Tdomain%out_var_capt(OUT_EPS_VOL) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_EPS_VOL)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_EPS_VOL)-1
                       capteur%valuecache(ioff+nComp,ni) = capteur%valuecache(ioff+nComp,ni) + weight*capteur%eps_vol(i,j,k)
                   end if

                   if (Tdomain%out_var_capt(OUT_EPS_DEV) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_EPS_DEV)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_EPS_DEV)-1
                       capteur%valuecache(ioff:ioff+nComp,ni) = capteur%valuecache(ioff:ioff+nComp,ni)+weight*capteur%eps_dev(i,j,k,:)
                   end if

                   if (Tdomain%out_var_capt(OUT_EPS_DEV_PL) == 1) then
                       ioff=Tdomain%out_var_offset(OUT_EPS_DEV_PL)
                       nComp = OUT_VAR_DIMS_3D(OUT_EPS_DEV_PL)-1
                       capteur%valuecache(ioff:ioff+nComp,ni) = capteur%valuecache(ioff:ioff+nComp,ni)+weight*capteur%eps_dev_pl(i,j,k,:)
                   end if

                   if (Tdomain%out_var_capt(OUT_STRESS_DEV) == 1) then
                       ioff = Tdomain%out_var_offset(OUT_STRESS_DEV)+1
                       nComp = OUT_VAR_DIMS_3D(OUT_STRESS_DEV)-1
                       capteur%valuecache(ioff:ioff+nComp,ni) = capteur%valuecache(ioff:ioff+nComp,ni) &
                       + (/weight*capteur%sig_dev(i,j,k,0), weight*capteur%sig_dev(i,j,k,1), &
                           weight*capteur%sig_dev(i,j,k,2), &
                           weight*capteur%sig_dev(i,j,k,3), weight*capteur%sig_dev(i,j,k,4), & 
                           weight*capteur%sig_dev(i,j,k,5)/)
                   end if
               enddo
           enddo
       enddo
!       write(*,*) "NI=", ni, Tdomain%out_var_capt(OUT_VITESSE)
!       do i=1,Tdomain%nReqOut+1
!           write(*,*) "CC", capteur%valuecache(i+1,ni)
!       end do
    end subroutine sortieGrandeurCapteur_interp

    subroutine sortieGrandeurCapteur_energy(Tdomain, ngll, capteur)
        !$acc routine worker
        use sdomain
        use dom_solid
        use dom_fluid
        use dom_solidpml
        use dom_fluidpml
        implicit none
        !
        type(domain)   :: TDomain
        type(tCapteur) :: capteur
        integer, intent(in) :: ngll
        !
        integer :: domain_type
        integer                                    :: i, j, k, n
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, R_energy, C_energy
        real(fpp) :: local_sum_P_energy, local_sum_S_energy, local_sum_R_energy, local_sum_C_energy
        real(fpp) :: global_sum_P_energy, global_sum_S_energy, global_sum_R_energy, global_sum_C_energy
        real(fpp), dimension(:), allocatable :: GLLw
        integer :: bnum, ee
        real(fpp), dimension(:,:,:), allocatable :: jac
        real(fpp) :: elem_P_En, elem_S_En, elem_R_En, elem_C_En
        type(Element), pointer :: el
        type(subdomain), pointer :: sub_dom_mat

        ! Verification : is an Energy Captor?
        if(capteur%type /= CPT_ENERGY) return
        !print *, "ENERGY CAPTOR"

        local_sum_P_energy = 0d0
        local_sum_S_energy = 0d0
        local_sum_R_energy = 0d0
        local_sum_C_energy = 0d0
        !count_P = 0
        !count_S = 0


        do n = 0,Tdomain%n_elem-1
            el => Tdomain%specel(n)
            domain_type = el%domain

            select case(domain_type)
                case (DM_SOLID_CG_PML)
                  cycle !We don't want the energy on PMLs
                case (DM_FLUID_CG_PML)
                  cycle !We don't want the energy on PMLs
                case (DM_SOLID_CG)
                    !We continue calculations
                case (DM_FLUID_CG)
                    !We continue calculations
                case default
                    stop "unknown domain"
            end select

            sub_dom_mat => Tdomain%sSubdomain(el%mat_index)
            bnum = el%lnum/VCHUNK
            ee = mod(el%lnum,VCHUNK)

            !print *, "BEFORE Jac"

            if(allocated(jac)) then
                if(size(jac) /= ngll*ngll*ngll) deallocate(jac)
            end if

            !print *, "INIT allocated(jac) = ", allocated(jac)
            if(.not. allocated(jac)) allocate(jac(0:ngll-1,0:ngll-1,0:ngll-1))
            jac (:,:,:) = 0.0d0

            !print *, "END allocated(jac) = ", allocated(jac)

            select case (domain_type)
                case (DM_SOLID_CG)
                    do k = 0, ngll-1
                        do j = 0, ngll-1
                            do i = 0, ngll-1
                                jac(i,j,k) = Tdomain%sdom%Jacob_   (i,j,k,bnum,ee)
                            enddo
                        enddo
                    enddo

                    call get_solid_dom_elem_energy(Tdomain%sdom, el%lnum, &
                                                   P_energy, S_energy, &
                                                   R_energy, C_energy)

                case (DM_FLUID_CG)
                    do k = 0, ngll-1
                        do j = 0, ngll-1
                            do i = 0, ngll-1
                                jac(i,j,k) = Tdomain%fdom%Jacob_   (i,j,k,bnum,ee)
                            enddo
                        enddo
                    enddo

                    call get_fluid_dom_elem_energy(Tdomain%fdom, el%lnum, P_energy, S_energy) !TODO

            end select

            if(allocated(GLLw)) deallocate(GLLw)
            call domain_gllw(Tdomain, domain_type, GLLw)

            call integrate_on_element(ngll, jac, GLLw, P_energy, elem_P_En)
            call integrate_on_element(ngll, jac, GLLw, S_energy, elem_S_En)
            call integrate_on_element(ngll, jac, GLLw, R_energy, elem_R_En)
            call integrate_on_element(ngll, jac, GLLw, C_energy, elem_C_En)

            local_sum_P_energy = local_sum_P_energy + elem_P_En
            local_sum_S_energy = local_sum_S_energy + elem_S_En
            local_sum_R_energy = local_sum_R_energy + elem_R_En
            local_sum_C_energy = local_sum_C_energy + elem_C_En

        enddo

        ! !TOTO, take out this part and put only local values (total values on post-processing)
        ! call MPI_ALLREDUCE(local_sum_S_energy, global_sum_S_energy, 1, MPI_DOUBLE_PRECISION, &
        !                    MPI_SUM, Tdomain%communicateur_global, ierr)
        ! call MPI_ALLREDUCE(local_sum_P_energy, global_sum_P_energy, 1, MPI_DOUBLE_PRECISION, &
        !                    MPI_SUM, Tdomain%communicateur_global, ierr)
        ! call MPI_ALLREDUCE(local_sum_R_energy, global_sum_R_energy, 1, MPI_DOUBLE_PRECISION, &
        !                    MPI_SUM, Tdomain%communicateur_global, ierr)
        ! call MPI_ALLREDUCE(local_sum_C_energy, global_sum_C_energy, 1, MPI_DOUBLE_PRECISION, &
        !                    MPI_SUM, Tdomain%communicateur_global, ierr)

        global_sum_P_energy = local_sum_P_energy
        global_sum_S_energy = local_sum_S_energy
        global_sum_R_energy = local_sum_R_energy
        global_sum_C_energy = local_sum_C_energy

        ! Sauvegarde des valeurs dans le capteur.
        i = capteur%icache
        capteur%valuecache(1,i) = Tdomain%TimeD%rtime
        capteur%valuecache(2,i) = global_sum_P_energy
        capteur%valuecache(3,i) = global_sum_S_energy
        capteur%valuecache(4,i) = global_sum_R_energy
        capteur%valuecache(5,i) = global_sum_C_energy
        capteur%valuecache(6,i) = global_sum_P_energy + global_sum_S_energy &
                                + global_sum_R_energy + global_sum_C_energy

        ! Deallocation.
        if(allocated(jac)) deallocate(jac)
        if(allocated(GLLw)) deallocate(GLLw)
        if(allocated(fieldU))   deallocate(fieldU)
        if(allocated(P_energy)) deallocate(P_energy)
        if(allocated(S_energy)) deallocate(S_energy)
        if(allocated(R_energy)) deallocate(R_energy)
        if(allocated(C_energy)) deallocate(C_energy)

    end subroutine sortieGrandeurCapteur_energy

    !!on identifie la maille dans laquelle se trouve le capteur. Il peut y en avoir plusieurs,
    !! alors le capteur est sur une face, arete ou sur un sommet
    !!
    subroutine trouve_capteur(Tdomain, xc, yc, zc, n_el, dmin, xi, eta, zeta, flag)
        use mshape8
        use mshape27
        use mlocations3d
        implicit none
        type (domain), INTENT(INOUT)  :: Tdomain
        real(fpp), intent(in) :: xc, yc, zc
        integer, intent(out) :: n_el
        real(fpp), intent(out) :: xi, eta, zeta, dmin
        logical, intent(in) :: flag
        !
        integer :: i, iel, emin, j
        logical :: inside
        integer :: nmax
        integer, parameter :: NMAXEL=100
        integer, dimension(NMAXEL) :: elems
        real(fpp), dimension(0:2,NMAXEL) :: coordloc
        real(fpp), parameter :: EPSD = 1D-10, EPSN = 0.1D0
        real(fpp) :: EPS
        real(fpp), dimension(0:2, 0:Tdomain%n_nodes-1) :: coordl
        real(fpp) :: dist, xc0, yc0, zc0

        EPS=EPSD*Tdomain%dxmax
        nmax = NMAXEL
        call find_location(Tdomain, xc, yc, zc, nmax, elems, coordloc)
        n_el = -1
        ! On calcule la distance entre la station et la maille (0 si interieur)
        dmin = 1e20
        emin = -1
        do i=1,nmax
            xi   = coordloc(0,i)
            eta  = coordloc(1,i)
            zeta = coordloc(2,i)
            iel = elems(i)
            do j = 0, Tdomain%n_nodes-1
                coordl(0:2, j) = Tdomain%Coord_Nodes(0:2, Tdomain%specel(iel)%Control_Nodes(j))
            enddo
            inside = .true.
            if ((xi<-1_fpp) .or. eta<(-1_fpp) .or. zeta<(-1_fpp)) inside = .false.
            if ((xi>+1_fpp) .or. eta>(+1_fpp) .or. zeta>(+1_fpp)) inside = .false.
            ! On projette sur le bord de l'element si besoin
            if (xi<-1D0) xi = -1D0
            if (xi>+1D0) xi = +1D0
            if (eta<-1D0) eta = -1D0
            if (eta>+1D0) eta = +1D0
            if (zeta<-1D0) zeta = -1D0
            if (zeta>+1D0) zeta = +1D0
            ! calcule less coord du nouveau pts
            if (.not. inside) then
                if (Tdomain%n_nodes==8) then
                    call shape8_local2global(coordl, xi, eta, zeta, xc0, yc0, zc0)
                else
                    call shape27_local2global(coordl, xi, eta, zeta, xc0, yc0, zc0)
                end if
                dist = (xc-xc0)**2 + (yc-yc0)**2 + (zc-zc0)**2
            else
                dist = 0.
            endif
            if (dist<dmin) then
                dmin = dist
                emin = iel
            endif
            if (flag .and. Tdomain%rank<40) then
                write(*,*) "CPT:", i, Tdomain%rank, xi, eta, zeta, inside, dist
            end if
        end do
        n_el = emin
        return
    end subroutine trouve_capteur

end module Mcapteur

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
