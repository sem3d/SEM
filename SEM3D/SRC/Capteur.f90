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
    use mpi
    use sem_hdf5
    use sem_c_config
    use constants
    use mshape8
    use mshape27
    use msnapshots, only : integrate_on_element
#include "index.h"
    implicit none

    public :: save_capteur, evalueSortieCapteur, flushAllCapteurs, create_capteurs
    private ::  flushCapteur
    ! start modifs
    integer, parameter :: CAPT_DIM=32
    ! end modifs
    type :: tCapteur
        type(tCapteur),pointer :: suivant ! pour passer au capteur suivant
        integer :: periode          ! frequence de captation des grandeur
        real, dimension (3) :: Coord  ! localisation du capteur
        character(LEN=20) :: nom      ! nom du capteur
        integer :: n_el ! numero de la maille dans laquelle se trouve le capteur
        ! si le capteur est partage entre plusieurs mailles, une seule suffit (type_calcul=1)
        real :: xi, eta, zeta ! abscisses curvilignes pour le capteur en cas d'interpolation (type_calcul=1)
        integer :: numproc               ! numero du proc localisant le capteur
        integer :: icache
        real, dimension(:,:), allocatable :: valuecache
        integer :: type
    end type tCapteur

    integer           :: dimCapteur        ! nombre total de capteurs

    type(Tcapteur), pointer :: listeCapteur
    type(Tcapteur), pointer :: capt_En_PS

    integer,parameter :: fileIdCapteur=200  ! id fichier capteur

    logical :: traces_h5_created
contains

    subroutine create_capteurs(Tdomain)
        implicit none
        type(domain), intent (inout) :: Tdomain
        !
        type(Tcapteur),pointer :: capteur
        type(C_PTR) :: station_next;
        type(sem_station), pointer :: station_ptr
        character(Len=MAX_FILE_SIZE) :: nom
        double precision :: xc, yc, zc, xi, eta, zeta
        character(len=MAX_FILE_SIZE) :: fnamef
        integer :: numproc, numproc_max, ierr, n_el, n_eln, i, n_out
        double precision, allocatable, dimension(:,:) :: coordl
        integer :: periodeRef


        station_next = Tdomain%config%stations
        nullify(listeCapteur)
        periodeRef = -1

        ! En reprise on ne recree pas le fichier capteur
        if (Tdomain%TimeD%NtimeMin==0) then
            traces_h5_created = .false.
        else
            traces_h5_created = .true.
        endif



        do while (C_ASSOCIATED(station_next))
            call c_f_pointer(station_next, station_ptr)
            xc = station_ptr%coords(1)
            yc = station_ptr%coords(2)
            zc = station_ptr%coords(3)
            numproc = -1
            call trouve_capteur(Tdomain, xc, yc, zc, n_el, n_eln, xi, eta, zeta)
            ! Cas ou le capteur est dans le maillage
            if (n_el/=-1) then
                numproc = Tdomain%rank
            end if
            call MPI_AllReduce(numproc, numproc_max, 1, MPI_INTEGER, &
                MPI_MAX, Tdomain%communicateur, ierr)
            ! Cas ou le capteur est en dehors du maillage mais < abs(1.1)
            if (numproc_max==-1) then
                n_el = n_eln
                if (n_el/=-1) then
                    numproc = Tdomain%rank
                end if
                call MPI_AllReduce(numproc, numproc_max, 1, MPI_INTEGER, &
                    MPI_MAX, Tdomain%communicateur, ierr)
                ! Re-calcul des coordonnees globales de la station apres
                ! modification des coordonnes locales
                if (n_el/=-1) then
                    allocate(coordl(0:2, 0:Tdomain%n_nodes-1))
                    do i = 0, Tdomain%n_nodes-1
                        coordl(0:2, i) = Tdomain%Coord_Nodes(0:2, Tdomain%specel(n_el)%Control_Nodes(i))
                    enddo
                    if (Tdomain%n_nodes==8) then
                        call shape8_local2global(coordl, xi, eta, zeta, xc, yc, zc)
                    else
                        call shape27_local2global(coordl, xi, eta, zeta, xc, yc, zc)
                    end if
                    deallocate(coordl)
                end if
            end if
            ! Cas ou le capteur est completement en dehors du maillage
            if (numproc_max==-1) then
                if (Tdomain%rank==0) then
                    write(*,*) "One of the station doesn't appear to be on any processor"
                    write(*,*) "Please verify that the station location is within the computation domain"
                end if
                stop 1
            end if


            ! attention si le capteur est partage par plusieurs procs. On choisit le proc de num max
            if(Tdomain%rank==numproc_max) then
                allocate(capteur)

                n_out = Tdomain%nReqOut
                if (.not.allocated(capteur%valuecache)) allocate(capteur%valuecache(1:n_out+1,NCAPT_CACHE))

                nom = fromcstr(station_ptr%name)
                capteur%nom = nom(1:20)     ! ses caracteristiques par defaut
                capteur%type = CPT_INTERP
                capteur%periode = station_ptr%period
                periodeRef = capteur%periode
                capteur%coord(1) = xc
                capteur%coord(2) = yc
                capteur%coord(3) = zc
                capteur%xi = xi
                capteur%eta = eta
                capteur%zeta = zeta
                capteur%n_el = n_el
                capteur%numproc = numproc_max
                capteur%icache = 0
                capteur%suivant => listeCapteur
                listeCapteur => capteur
                write(*,"(A,A,A,I5,A,I6,A,F8.4,A,F8.4,A,F8.4)") "Capteur:", trim(capteur%nom), &
                    " on proc ", Tdomain%rank, " in elem ", n_el, " at ", xi, ",", eta, ",", zeta

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

        if(periodeRef < 1) periodeRef = 1

        ! Energy outputs
        if(Tdomain%out_variables(OUT_TOTAL_ENERGY) == 1) then

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
            capt_En_PS%icache = 0
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

        integer :: ntime
        type (domain) :: TDomain

        type(tCapteur),pointer :: capteur
        logical :: do_flush

        do_flush = .false.
        ! boucle sur les capteurs


        capteur=>listeCapteur
        do while (associated(capteur))
            if (mod(ntime, capteur%periode)==0) then ! on fait la sortie
                if (capteur%type == CPT_INTERP) then
                    call sortieGrandeurCapteur_interp(Tdomain, capteur)
                else if (capteur%type == CPT_ENERGY) then
                    !print*, "BEFORE sortieGrandeurCapteur_energy"
                    call sortieGrandeurCapteur_energy(Tdomain, capteur)
                    !print*, "AFTER sortieGrandeurCapteur_energy"
                end if
                if (capteur%icache==NCAPT_CACHE) do_flush = .true.
            end if
            capteur=>capteur%suivant
        enddo

        if (do_flush) call flushAllCapteurs(Tdomain)

    end subroutine save_capteur

    function dset_capteur_name(capteur)
        implicit none
        type(tCapteur),pointer :: capteur
        character(len=40) :: dset_capteur_name
        dset_capteur_name = trim(adjustl(capteur%nom)) !//"_"//trim(adjustl(capteur%grandeur))
    end function dset_capteur_name

    subroutine create_traces_h5_skel(Tdomain)
        use HDF5
        implicit none
        type (domain), intent(inout) :: TDomain
        type(tCapteur),pointer :: capteur
        character (len=MAX_FILE_SIZE) :: fnamef
        character (len=40) :: dname
        integer(HID_T) :: fid, dset_id
        integer :: hdferr

        
        call init_hdf5()

        call semname_tracefile_h5(Tdomain%rank, fnamef)
        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
        call create_capteur_descriptions(Tdomain, fid)
        

        capteur=>listeCapteur
        do while (associated(capteur))
            dname = dset_capteur_name(capteur)
            call create_dset_2d(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, &
                int(CAPT_DIM,HSIZE_T), int(H5S_UNLIMITED_F,HSIZE_T), dset_id)
            call h5dclose_f(dset_id, hdferr)
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
        character(len=12) :: temp
        integer :: d,k,dim,dimtot
        integer(HSIZE_T), dimension(1) :: dims
        
        dimtot = sum(OUT_VAR_DIMS_3D)-OUT_VAR_DIMS_3D(OUT_TOTAL_ENERGY)
        allocate(varnames(0:dimtot))
        varnames(0) = "Time"
        d = 1
        do k=0,dimtot-1
            if (Tdomain%out_variables(k)==1) then
                if(k == OUT_TOTAL_ENERGY) cycle 
                do dim=1,OUT_VAR_DIMS_3D(k)
                    write(temp,"(A,I2)") OUT_VAR_NAMES(k),dim
                    varnames(d) = temp
                    d = d+1
                end do
            end if
        end do
        dims(1) = d
        call H5Tcopy_f(H5T_FORTRAN_S1, tid, hdferr)
        call H5Tset_size_f(tid, 12_HSIZE_T, hdferr)
        call H5Screate_simple_f(1, dims, spaceid, hdferr)
        call H5Dcreate_f(fid, "Variables", tid, spaceid, dsetid, hdferr)
        call H5Dwrite_f(dsetid, tid, varnames, dims, hdferr, spaceid, spaceid)
        call H5Dclose_f(dsetid, hdferr)
        call H5Sclose_f(spaceid, hdferr)
        call H5Tclose_f(tid, hdferr)
        !
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
            if (capteur%icache==0) then
                capteur=>capteur%suivant
                cycle
            endif
            call h5dopen_f(fid, trim(dname), dset_id, hdferr)
            call append_dataset_2d(dset_id, capteur%valuecache(:,1:capteur%icache), hdferr)
            call h5dclose_f(dset_id, hdferr)
            capteur%icache=0
            capteur=>capteur%suivant
        enddo

        call h5fclose_f(fid, hdferr)
    end subroutine append_traces_h5

    subroutine flushAllCapteurs(Tdomain)
        implicit none
        type (domain), intent(inout) :: TDomain
        type(tCapteur),pointer :: capteur

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

        if (capteur%icache==0) return

        call semname_capteur_type(capteur%nom,".txt",fnamef)

        open(fileId,file=trim(fnamef),status="unknown",form="formatted",position="append")
        do j=1,capteur%icache
            ! start modifs
            write(sizeChar, *) size(capteur%valuecache)
            write(fileId,'('//trim(sizeChar)//'(1X,E16.8E3))') capteur%valuecache(:,j)
            ! end modifs
        end do
        close(fileId)
        capteur%icache = 0
    end subroutine flushCapteur


    !---------------------------------------------------------------------
    !---------------------------------------------------------------------


    !! effectue l'interpolation des grandeurs dans la maille dans laquelle se trouve le capteur
    !! la maille se trouve dans un seul proc
    !! seul le proc gere l'ecriture
    !!
    subroutine sortieGrandeurCapteur_interp(Tdomain, capteur)
        use constants
        use dom_solid
        use dom_fluid
        use dom_solidpml
        use dom_fluidpml
        implicit none
        !
        type(domain)   :: TDomain
        type(tCapteur) :: capteur
        !
        integer                                    :: i, j, k, ioff
        integer                                    :: n_el, ngll
        real(fpp)                                  :: weight
        real(fpp), dimension(:), allocatable       :: outx, outy, outz
        real(fpp), dimension(:), allocatable       :: grandeur
        integer, dimension(0:size(Tdomain%out_variables)-1):: out_variables, offset
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev_pl
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        real, dimension(:), allocatable :: GLLc
        logical :: nl_flag
        integer :: nComp

        ! Verification : le capteur est il gere par le proc. ?

        n_el = capteur%n_el
        if((n_el==-1) .OR. (capteur%numproc/=Tdomain%rank)) return

        ! Initialisation.
        ngll = domain_ngll(Tdomain, Tdomain%specel(n_el)%domain)
        call domain_gllc(Tdomain, Tdomain%specel(n_el)%domain, GLLc)

        allocate(outx(0:ngll-1))
        allocate(outy(0:ngll-1))
        allocate(outz(0:ngll-1))
        do i = 0,ngll - 1
            call  pol_lagrange(ngll,GLLc,i,capteur%xi,outx(i))
        end do
        do j = 0,ngll - 1
            call  pol_lagrange(ngll,GLLc,j,capteur%eta,outy(j))
        end do
        do k = 0,ngll - 1
            call  pol_lagrange(ngll,GLLc,k,capteur%zeta,outz(k))
        end do
        deallocate(GLLc)
        
        allocate(grandeur(0:Tdomain%nReqOut-1))
        grandeur(:) = 0. ! si maillage vide donc pas de pdg, on fait comme si il y en avait 1
        
        out_variables(:) = Tdomain%out_variables(:)
        nl_flag = Tdomain%nl_flag
        offset = 0
        do i = 0,size(out_variables)-2
            if (out_variables(i) == 1) then
                offset(i+1) = offset(i) + OUT_VAR_DIMS_3D(i)
            else
                offset(i+1) = offset(i)
            end if
        end do

        ! On recupere les variables de l'element associe au capteur.
        
        select case(Tdomain%specel(n_el)%domain)
            case (DM_SOLID)
              call get_solid_dom_var(Tdomain%sdom, Tdomain%specel(n_el)%lnum, out_variables, &
                fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev, &
                nl_flag, eps_dev_pl)
            case (DM_FLUID)
              call get_fluid_dom_var(Tdomain%fdom, Tdomain%specel(n_el)%lnum, out_variables, &
                fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
            case (DM_SOLID_PML)
              call get_solidpml_dom_var(Tdomain%spmldom, Tdomain%specel(n_el)%lnum, out_variables, &
                fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
            case (DM_FLUID_PML)
              call get_fluidpml_dom_var(Tdomain%fpmldom, Tdomain%specel(n_el)%lnum, out_variables, &
                fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
            case default
              stop "unknown domain"
        end select

        ! On interpole le DOF a la position du capteur.

        do i = 0,ngll - 1
            do j = 0,ngll - 1
                do k = 0,ngll - 1
                    weight = outx(i)*outy(j)*outz(k)

                    if (out_variables(OUT_DEPLA) == 1 .AND. allocated(fieldU)) then
                        ioff = offset(OUT_DEPLA)
                        nComp = OUT_VAR_DIMS_3D(OUT_DEPLA)-1
                        grandeur(ioff:ioff+nComp) &
                            = grandeur(ioff:ioff+nComp) + weight*fieldU(i,j,k,:)
                    end if

                    if (out_variables(OUT_VITESSE) == 1 .AND. allocated(fieldV)) then
                        ioff = offset(OUT_VITESSE)
                        nComp = OUT_VAR_DIMS_3D(OUT_VITESSE)-1
                        grandeur(ioff:ioff+nComp) &
                            = grandeur(ioff:ioff+nComp) + weight*fieldV(i,j,k,:)
                    end if

                    if (out_variables(OUT_ACCEL) == 1 .AND. allocated(fieldA)) then
                        ioff = offset(OUT_ACCEL)
                        nComp = OUT_VAR_DIMS_3D(OUT_ACCEL)-1
                        grandeur(ioff:ioff+nComp) &
                            = grandeur(ioff:ioff+nComp) + weight*fieldA(i,j,k,:)
                    end if

                    if (out_variables(OUT_PRESSION) == 1 .AND. allocated(fieldP)) then
                        ioff = offset(OUT_PRESSION)
                        nComp = OUT_VAR_DIMS_3D(OUT_PRESSION)-1
                        grandeur(ioff+nComp) &
                            = grandeur(ioff+nComp) + weight*fieldP(i,j,k)
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        ioff = offset(OUT_ENERGYP)
                        nComp = OUT_VAR_DIMS_3D(OUT_ENERGYP)-1
                        grandeur (ioff+nComp) = grandeur (ioff+nComp) + weight*P_energy(i,j,k)
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        ioff = offset(OUT_ENERGYS)
                        nComp = OUT_VAR_DIMS_3D(OUT_ENERGYS)-1
                        grandeur (ioff+nComp) = grandeur (ioff+nComp) + weight*S_energy(i,j,k)
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        ioff = offset(OUT_EPS_VOL)
                        nComp = OUT_VAR_DIMS_3D(OUT_EPS_VOL)-1
                        grandeur (ioff+nComp) = grandeur (ioff+nComp) + weight*eps_vol(i,j,k)
                    end if
                     
                    if (out_variables(OUT_EPS_DEV) == 1) then
                        ioff = offset(OUT_EPS_DEV)
                        nComp = OUT_VAR_DIMS_3D(OUT_EPS_DEV)-1
                        grandeur (ioff:ioff+nComp) = grandeur(ioff:ioff+nComp)+weight*eps_dev(i,j,k,:)
                    end if

                    if (out_variables(OUT_EPS_DEV_PL) == 1) then
                        ioff=offset(OUT_EPS_DEV_PL)
                        nComp = OUT_VAR_DIMS_3D(OUT_EPS_DEV_PL)-1
                        grandeur (ioff:ioff+nComp) = grandeur(ioff:ioff+nComp)+weight*eps_dev_pl(i,j,k,:)
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        ioff = offset(OUT_STRESS_DEV)
                        nComp = OUT_VAR_DIMS_3D(OUT_STRESS_DEV)-1
                        grandeur (ioff:ioff+nComp) = grandeur (ioff:ioff+nComp) &
                        + (/weight*sig_dev(i,j,k,0), weight*sig_dev(i,j,k,1), weight*sig_dev(i,j,k,2), &
                            weight*sig_dev(i,j,k,3), weight*sig_dev(i,j,k,4), weight*sig_dev(i,j,k,5)/)
                    end if
                enddo
            enddo
        enddo
        ! Sauvegarde des valeurs dans le capteur.

        i = capteur%icache+1
        capteur%valuecache(1,i) = Tdomain%TimeD%rtime
        capteur%valuecache(2:Tdomain%nReqOut+1,i) = grandeur(:)
        capteur%icache = i

        ! Deallocation.

        if(allocated(fieldU))   deallocate(fieldU)
        if(allocated(fieldV))   deallocate(fieldV)
        if(allocated(fieldA))   deallocate(fieldA)
        if(allocated(fieldP))   deallocate(fieldP)
        if(allocated(P_energy)) deallocate(P_energy)
        if(allocated(S_energy)) deallocate(S_energy)
        if(allocated(eps_vol))  deallocate(eps_vol)
        if(allocated(eps_dev))  deallocate(eps_dev)
        if(allocated(eps_dev_pl))  deallocate(eps_dev_pl)
        if(allocated(sig_dev))  deallocate(sig_dev)
        if(allocated(grandeur)) deallocate(grandeur)
        deallocate(outx)
        deallocate(outy)
        deallocate(outz)
    end subroutine sortieGrandeurCapteur_interp

    subroutine sortieGrandeurCapteur_energy(Tdomain, capteur)
        use sdomain
        use dom_solid
        use dom_fluid
        use dom_solidpml
        use dom_fluidpml
        implicit none
        !
        type(domain)   :: TDomain
        type(tCapteur) :: capteur
        !
        integer :: domain_type
        integer                                    :: i, j, k, n, ioff
        integer                                    :: n_el, ngll
        real(fpp)                                  :: weight
        real(fpp), dimension(:), allocatable       :: grandeur
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, R_energy, C_energy
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        real, dimension(:), allocatable :: GLLc
        real(fpp) :: local_sum_P_energy, local_sum_S_energy, local_sum_R_energy, local_sum_C_energy
        real(fpp) :: global_sum_P_energy, global_sum_S_energy, global_sum_R_energy, global_sum_C_energy
        real(fpp) :: Whei, mult
        real, dimension(:), allocatable :: GLLw
        integer :: bnum, ee
        real(fpp), dimension(:,:,:), allocatable :: jac
        real(fpp) :: elem_P_En, elem_S_En, elem_R_En, elem_C_En
        type(Element), pointer :: el
        type(subdomain), pointer :: sub_dom_mat
        integer :: ierr


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
                case (DM_SOLID_PML)
                  cycle !We don't want the energy on PMLs
                case (DM_FLUID_PML)
                  cycle !We don't want the energy on PMLs
                case (DM_SOLID)
                    !We continue calculations
                case (DM_FLUID)
                    !We continue calculations
                case default
                    stop "unknown domain"
            end select

            sub_dom_mat => Tdomain%sSubdomain(el%mat_index)
            ngll = domain_ngll(Tdomain, el%domain)
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
                case (DM_SOLID)
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

                case (DM_FLUID)
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

        !TOTO, take out this part and put only local values (total values on post-processing)
        call MPI_ALLREDUCE(local_sum_S_energy, global_sum_S_energy, 1, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, Tdomain%communicateur_global, ierr)
        call MPI_ALLREDUCE(local_sum_P_energy, global_sum_P_energy, 1, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, Tdomain%communicateur_global, ierr)
        call MPI_ALLREDUCE(local_sum_R_energy, global_sum_R_energy, 1, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, Tdomain%communicateur_global, ierr)
        call MPI_ALLREDUCE(local_sum_C_energy, global_sum_C_energy, 1, MPI_DOUBLE_PRECISION, &
                           MPI_SUM, Tdomain%communicateur_global, ierr)

        ! Sauvegarde des valeurs dans le capteur.
        i = capteur%icache+1
        capteur%valuecache(1,i) = Tdomain%TimeD%rtime
        capteur%valuecache(2,i) = global_sum_P_energy
        capteur%valuecache(3,i) = global_sum_S_energy
        capteur%valuecache(4,i) = global_sum_R_energy
        capteur%valuecache(5,i) = global_sum_C_energy
        capteur%valuecache(6,i) = global_sum_P_energy + global_sum_S_energy &
                                + global_sum_R_energy + global_sum_C_energy
        capteur%icache = i

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
    subroutine trouve_capteur(Tdomain, xc, yc, zc, n_el, n_eln, xi, eta, zeta)
        use mshape8
        use mshape27
        use mlocations3d
        implicit none
        type (domain), INTENT(INOUT)  :: Tdomain
        double precision, intent(in) :: xc, yc, zc
        integer, intent(out) :: n_el, n_eln
        double precision, intent(out) :: xi, eta, zeta
        !
        integer :: i
        logical :: inside
        integer :: nmax
        integer, parameter :: NMAXEL=20
        integer, dimension(NMAXEL) :: elems
        double precision, dimension(0:2,NMAXEL) :: coordloc
        double precision, parameter :: EPS = 1D-13, EPSN = 0.1D0

        nmax = NMAXEL
        call find_location(Tdomain, xc, yc, zc, nmax, elems, coordloc)
        n_el = -1
        n_eln = -1
        ! Cas ou la station est dans le maillage
        do i=1,nmax
            inside = .true.
            xi   = coordloc(0,i)
            eta  = coordloc(1,i)
            zeta = coordloc(2,i)
            if (xi<(-1-EPS) .or. eta<(-1-EPS) .or. zeta<(-1-EPS)) inside = .false.
            if (xi>(1+EPS) .or. eta>(1+EPS) .or. zeta>(1+EPS)) inside = .false.
            if (inside) then
                n_el = elems(i)
                return
            end if
        end do
        ! Cas ou la station est en dehors du maillage mais < abs(1.1)
        do i=1,nmax
            xi   = coordloc(0,i)
            eta  = coordloc(1,i)
            zeta = coordloc(2,i)
            if (abs(xi)<(1+EPS) .and. abs(eta)<(1+EPS) .and. abs(zeta)<(1+EPSN)) then
                if (zeta<0) zeta = -1.D0
                if (zeta>0) zeta = 1.D0
                n_eln = elems(i)
                return
            else if (abs(xi)<(1+EPS) .and. abs(eta)<(1+EPSN) .and. abs(zeta)<(1+EPS)) then
                if (eta<0) eta = -1.D0
                if (eta>0) eta = 1.D0
                n_eln = elems(i)
                return
            else if (abs(xi)<(1+EPSN) .and. abs(eta)<(1+EPS) .and. abs(zeta)<(1+EPS)) then
                if (xi<0) xi = -1.D0
                if (xi>0) xi = 1.D0
                n_eln = elems(i)
                return
            end if
        end do
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
