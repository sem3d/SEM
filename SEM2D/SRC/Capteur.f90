!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Capteur.f90
!!\brief Permet de manipuler les quantités associées aux capteurs en 2D.
!!\version 1.0
!!\date 2026-07-11
!!
!<

module mCapteur

    use, intrinsic :: ISO_C_BINDING, only : C_PTR
    use sdomain, only : domain
    use semdatafiles
    use mpi
    use sem_hdf5
    use sem_c_config, only : sem_station, fromcstr
    use constants
    use shape_lin
    use shape_quad
    use orientation
    implicit none

    public :: save_capteur, evalueSortieCapteur, flushAllCapteurs, create_capteurs
    private ::  flushCapteur

    type :: tCapteur
        type(tCapteur),pointer :: suivant ! pour passer au capteur suivant
        integer :: periode          ! frequence de captation des grandeur
        real(fpp), dimension (2) :: Coord  ! localisation du capteur (X, Z)
        character(LEN=20) :: nom      ! nom du capteur
        integer :: n_el ! numero de la maille dans laquelle se trouve le capteur
        real(fpp) :: xi, eta ! abscisses curvilignes pour le capteur en cas d'interpolation
        integer :: numproc               ! numero du proc localisant le capteur
        integer :: icache
        real(fpp), dimension(:,:), allocatable :: valuecache
        integer :: type
    end type tCapteur

    integer           :: dimCapteur        ! nombre total de capteurs

    type(Tcapteur), pointer :: listeCapteur
    type(Tcapteur), pointer :: capt_Energy

    integer,parameter :: fileIdCapteur=200  ! id fichier capteur

    logical :: traces_h5_created
contains

    subroutine create_capteurs(Tdomain)
        implicit none
        type(domain), intent (inout) :: Tdomain
        !
        type(Tcapteur),pointer :: capteur
        type(C_PTR) :: station_next
        type(sem_station), pointer :: station_ptr
        character(Len=MAX_FILE_SIZE) :: nom
        real(fpp) :: xc, zc, xi, eta
        real(fpp) :: xc0, zc0
        character(len=MAX_FILE_SIZE) :: fnamef
        integer :: numproc, numproc_max, ierr, n_el, i, n_out
        real(fpp) :: dmin, glob_dmin
        real(fpp), dimension(0:1, 0:Tdomain%n_nodes-1) :: coordl
        integer :: periodeRef
        logical :: flag

        flag = .false.
        station_next = Tdomain%stations
        nullify(listeCapteur)
        periodeRef = -1

        ! En reprise on ne recree pas le fichier capteur
        if (Tdomain%TimeD%NtimeMin==0) then
            traces_h5_created = .false.
        else
            traces_h5_created = .true.
        endif

        Tdomain%has_station = .false. !Stations other than Total Energy

        do while (C_ASSOCIATED(station_next))
            call c_f_pointer(station_next, station_ptr)
            xc = station_ptr%coords(1)
            zc = station_ptr%coords(2)
            numproc = -1
            nom = fromcstr(station_ptr%name)
            call trouve_capteur(Tdomain, xc, zc, n_el, dmin, xi, eta, flag)
            ! Cas ou le capteur est dans le maillage
            call MPI_AllReduce(dmin, glob_dmin, 1, MPI_DOUBLE, &
                MPI_MIN, Tdomain%communicateur, ierr)
            if ((n_el >= 0) .and. (dmin==glob_dmin)) then
                numproc = Tdomain%Mpi_var%my_rank
            end if
            call MPI_AllReduce(numproc, numproc_max, 1, MPI_INTEGER, &
                MPI_MAX, Tdomain%communicateur, ierr)
            ! Cas ou le capteur est completement en dehors du maillage
            if (numproc_max==-1) then
                if (Tdomain%Mpi_var%my_rank==0) then
                    nom = fromcstr(station_ptr%name)
                    write(*,*) "Something wrong happened..."
                    write(*,*) "The station ", trim(nom), " doesn't appear to be on any processor"
                    write(*,*) "Please verify that the station location is within the computation domain"
                end if
                stop 1
            end if

            ! attention si le capteur est partage par plusieurs procs. On choisit le proc de num max
            if(Tdomain%Mpi_var%my_rank==numproc_max) then
                if (n_el < 0) then
                    write(*,*) "Internal error while creating station ", trim(nom)
                    write(*,*) "Selected rank has no containing element (n_el = ", n_el, ")"
                    write(*,*) "Please check station coordinates and mesh decomposition"
                    stop 1
                end if
                allocate(capteur)
                Tdomain%has_station = .true.
                n_out = Tdomain%nReqOut
                if (.not.allocated(capteur%valuecache)) allocate(capteur%valuecache(1:n_out+1,NCAPT_CACHE))

                if (glob_dmin>0) then
                    do i = 0, Tdomain%n_nodes-1
                        coordl(0:1, i) = Tdomain%Coord_Nodes(0:1, Tdomain%specel(n_el)%Control_Nodes(i))
                    enddo
                    if (Tdomain%n_nodes==4) then
                        call shape4_local2global(coordl, xi, eta, xc0, zc0)
                    else
                        call shape8_local2global(coordl, xi, eta, xc0, zc0)
                    end if
                    write(*,*) "The station",trim(nom)," is outside. Moved from ", xc, zc, " to ", xc0, zc0
                else
                    xc0 = xc
                    zc0 = zc
                end if
                nom = fromcstr(station_ptr%name)
                capteur%nom = nom(1:20)     ! ses caracteristiques par defaut
                capteur%type = CPT_INTERP
                capteur%periode = station_ptr%period
                periodeRef = capteur%periode
                capteur%coord(1) = xc0
                capteur%coord(2) = zc0
                capteur%xi = xi
                capteur%eta = eta
                capteur%n_el = n_el
                capteur%numproc = numproc_max
                capteur%icache = 0
                capteur%suivant => listeCapteur
                listeCapteur => capteur
                write(*,"(A,A,A,I5,A,I6,A,F8.4,A,F8.4,A,I1)") "Capteur:", trim(capteur%nom), &
                    " on proc ", Tdomain%Mpi_var%my_rank, " in elem ", n_el, " at ", xi, ",", eta, &
                    " in domain ", Tdomain%specel(n_el)%mat_index

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
        if(Tdomain%out_var_capt(OUT_TOTAL_ENERGY) == 1) then

            if(Tdomain%Mpi_var%my_rank == 0) write(*,*) "CREATING ENERGY SENSORS"

            allocate(capt_Energy)

            n_out = 5
            if (.not.allocated(capt_Energy%valuecache)) allocate(capt_Energy%valuecache(1:n_out+1,NCAPT_CACHE))

            capt_Energy%nom = "Energy"
            capt_Energy%type = CPT_ENERGY
            capt_Energy%periode = 1 
            capt_Energy%coord(1) = -1111
            capt_Energy%coord(2) = -1111
            capt_Energy%xi = -1111
            capt_Energy%eta = -1111
            capt_Energy%n_el = -1
            capt_Energy%numproc = Tdomain%Mpi_var%my_rank
            capt_Energy%icache = 0
            capt_Energy%suivant => listeCapteur
            listeCapteur => capt_Energy
            write(*,"(A,A,A,I5,A,I6,A,F8.4,A,F8.4)") "Capteur:", trim(capt_Energy%nom), &
                " on proc ", Tdomain%Mpi_var%my_rank, " in elem ", n_el, " at ", xi, ",", eta

        end if

        Tdomain%has_station = associated(listeCapteur)

    end subroutine create_capteurs


    subroutine evalueSortieCapteur(it, sortie_capteur)
        implicit none
        integer, intent(in) :: it
        logical, intent(out) :: sortie_capteur
        type(tCapteur),pointer :: capteur

        sortie_capteur = .FALSE.
        capteur=>listeCapteur
        do while (associated(capteur))
            if(capteur%periode < 1) stop "ERROR, station with period smaller than 1"
            if (mod(it,capteur%periode)==0) then ! on fait la sortie
                sortie_capteur = .TRUE.
            endif
            capteur=>capteur%suivant
        enddo
    end subroutine evalueSortieCapteur


    subroutine save_capteur(Tdomain, ntime)
        implicit none

        integer :: ntime
        type (domain) :: TDomain

        type(tCapteur),pointer :: capteur
        logical :: do_flush

        do_flush = .false.
        capteur=>listeCapteur
        do while (associated(capteur))
            if (mod(ntime, capteur%periode)==0) then ! on fait la sortie
                if (capteur%type == CPT_INTERP) then
                    call sortieGrandeurCapteur_interp(Tdomain, capteur)
                else if (capteur%type == CPT_ENERGY) then
                    call sortieGrandeurCapteur_energy(Tdomain, capteur)
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
        dset_capteur_name = trim(adjustl(capteur%nom))
    end function dset_capteur_name

    function dset_capteur_posname(capteur)
        implicit none
        type(tCapteur),pointer :: capteur
        character(len=40) :: dset_capteur_posname
        dset_capteur_posname = trim(adjustl(capteur%nom)) //"_pos"
    end function dset_capteur_posname

    subroutine create_traces_h5_skel(Tdomain)
        use HDF5
        use sem_git_version, only : SEM_GIT_HASH, SEM_GIT_DIRTY
        implicit none
        type (domain), intent(inout) :: TDomain
        type(tCapteur),pointer :: capteur
        character (len=MAX_FILE_SIZE) :: fnamef
        character (len=40) :: dname
        integer(HID_T) :: fid, dset_id
        integer :: hdferr, n_out

        call init_hdf5()

        call semname_tracefile_h5(Tdomain%Mpi_var%my_rank, fnamef)
        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
        call write_attr_string(fid, "Code", "SEM2D")
        call write_attr_string(fid, "GitHash", SEM_GIT_HASH)
        call write_attr_string(fid, "GitStatus", SEM_GIT_DIRTY)
        call create_capteur_descriptions(Tdomain, fid)

        n_out = Tdomain%nReqOut+1

        capteur=>listeCapteur
        do while (associated(capteur))
            dname = dset_capteur_name(capteur)
            if (capteur%type == CPT_ENERGY) then
                n_out = 6
            else
                n_out = Tdomain%nReqOut+1
            endif
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
        use constants, only : OUT_VAR_NAMES, OUT_VAR_DIMS_2D
        type (domain), intent(inout) :: TDomain
        integer(HID_T), intent(in) :: fid
        !
        integer(HID_T) :: tid, dsetid, spaceid
        integer :: hdferr
        character(len=12), dimension(:), allocatable :: varnames
        character(len=12), dimension(6) :: energy_varnames = ["Time       1", &
                                                              "EnergyP    1", &
                                                              "EnergyK    1", &
                                                              "EnergyL    1", &
                                                              "EnergyS    1", &
                                                              "EnergyR    1"]
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
                do dim=1,OUT_VAR_DIMS_2D(k)
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
            call H5Dcreate_f(fid, "Energy_Variables", tid, spaceid, dsetid, hdferr)
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

        call semname_tracefile_h5(Tdomain%Mpi_var%my_rank, fnamef)

        call h5fopen_f(fnamef, H5F_ACC_RDWR_F, fid, hdferr)

        capteur=>listeCapteur
        do while (associated(capteur))
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
            write(sizeChar, *) size(capteur%valuecache, 1)
            write(fileId,'('//trim(sizeChar)//'(1X,E16.8E3))') capteur%valuecache(:,j)
        end do
        close(fileId)
        capteur%icache = 0
    end subroutine flushCapteur


    subroutine sortieGrandeurCapteur_interp(Tdomain, capteur)
        implicit none
        !
        type(domain)   :: TDomain
        type(tCapteur) :: capteur
        !
        integer                                    :: i, j, ioff
        integer                                    :: n_el, ngllx, ngllz, mat
        real(fpp)                                  :: weight
        ! Evaluation of lagrange polynomial at xi/eta capteur
        real(fpp), dimension(:), allocatable       :: outx, outz
        real(fpp), dimension(:), allocatable       :: grandeur
        integer, dimension(0:OUT_LAST)             :: out_variables, offset
        real(fpp), dimension(:,:,:), allocatable   :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:), allocatable     :: fieldP
        real(fpp), dimension(:,:), allocatable     :: P_energy, K_energy, eps_vol
        real(fpp), dimension(:,:,:), allocatable   :: D_energy
        real(fpp), dimension(:,:,:), allocatable   :: dUdX
        real(fpp), dimension(:,:,:), allocatable   :: eps_dev
        real(fpp), dimension(:,:,:), allocatable   :: eps_dev_pl
        real(fpp), dimension(:,:,:), allocatable   :: sig_dev
        real(fpp), dimension(:), allocatable       :: GLLcx, GLLcz
        integer :: nComp, r
        real(fpp), dimension(0:1,0:1)              :: invgrad_ij
        real(fpp)                                  :: dUx_dxi, dUx_deta, dUz_dxi, dUz_deta
        real(fpp)                                  :: DXX, DXZ, DZX, DZZ
        real(fpp)                                  :: eps_xx, eps_zz, eps_xz, eps_v
        real(fpp)                                  :: sig_xx, sig_zz, sig_xz, sig_mean
        real(fpp), dimension(3)                    :: strain_v, stress_v

        n_el = capteur%n_el
        if((n_el==-1) .OR. (capteur%numproc/=Tdomain%Mpi_var%my_rank)) return

        ! Initialisation
        ngllx = Tdomain%specel(n_el)%ngllx
        ngllz = Tdomain%specel(n_el)%ngllz
        mat   = Tdomain%specel(n_el)%mat_index

        GLLcx = Tdomain%sSubdomain(mat)%GLLcx
        GLLcz = Tdomain%sSubdomain(mat)%GLLcz

        allocate(fieldP(0:ngllx-1,0:ngllz-1))
        allocate(fieldU(0:ngllx-1,0:ngllz-1,0:1))
        allocate(fieldV(0:ngllx-1,0:ngllz-1,0:1))
        allocate(fieldA(0:ngllx-1,0:ngllz-1,0:1))
        allocate(eps_vol(0:ngllx-1,0:ngllz-1))
        allocate(P_energy(0:ngllx-1,0:ngllz-1))
        allocate(K_energy(0:ngllx-1,0:ngllz-1))
        allocate(D_energy(0:ngllx-1,0:ngllz-1,0:2))
        allocate(eps_dev(0:ngllx-1,0:ngllz-1,0:2))
        allocate(dUdX(0:ngllx-1,0:ngllz-1,0:3))
        allocate(eps_dev_pl(0:ngllx-1,0:ngllz-1,0:2))
        allocate(sig_dev(0:ngllx-1,0:ngllz-1,0:2))
        allocate(outx(0:ngllx-1))
        allocate(outz(0:ngllz-1))

        do i = 0,ngllx - 1
            call  pol_lagrange(ngllx,GLLcx,i,capteur%xi,outx(i))
        end do
        do j = 0,ngllz - 1
            call  pol_lagrange(ngllz,GLLcz,j,capteur%eta,outz(j))
        end do

        allocate(grandeur(0:Tdomain%nReqOut-1))
        grandeur(:) = 0.

        out_variables(:) = Tdomain%out_var_capt(:)
        offset = 0
        do i = 0,size(out_variables)-2
            if (out_variables(i) == 1) then
                offset(i+1) = offset(i) + OUT_VAR_DIMS_2D(i)
            else
                offset(i+1) = offset(i)
            end if
        end do

        ! Gather variables
        if (Tdomain%type_timeInteg == TIME_INTEG_RK4) then
            fieldU = Tdomain%specel(n_el)%Displ
            fieldV = Tdomain%specel(n_el)%Veloc
            fieldA = Tdomain%specel(n_el)%Accel
        else
            call gather_elem_displ(Tdomain, n_el, fieldU)
            call gather_elem_veloc(Tdomain, n_el, fieldV, .false.)
            call gather_elem_accel(Tdomain, n_el, fieldA)
        endif

        ! On compute les contraintes/déformations et énergies locales sur l'élément GLL
        fieldP = 0.0_fpp
        eps_vol = 0.0_fpp
        P_energy = 0.0_fpp
        K_energy = 0.0_fpp
        D_energy = 0.0_fpp
        eps_dev = 0.0_fpp
        eps_dev_pl = 0.0_fpp
        sig_dev = 0.0_fpp
        dUdX = 0.0_fpp

        do j = 0,ngllz-1
            do i = 0,ngllx-1
                ! Calcul de l'énergie cinétique
                K_energy(i,j) = 0.5_fpp * Tdomain%specel(n_el)%Density(i,j) * (fieldV(i,j,0)**2 + fieldV(i,j,1)**2)

                if (Tdomain%specel(n_el)%acoustic) then
                    ! Cas Fluide (Potentiel ou déplacement)
                    if (allocated(Tdomain%specel(n_el)%IDensTensor2d)) then
                        ! Fluide-anisotrope (Velocity potential field): p = -VelPhi = -veloc(0)
                        fieldP(i,j) = -fieldV(i,j,0)
                    else
                        ! Fluide-isotrope (displacement-based): p = -lambda * div(U)
                        dUx_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * fieldU(:,j,0))
                        dUx_deta = sum(fieldU(i,:,0) * Tdomain%sSubdomain(mat)%hprimez(:,j))
                        dUz_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * fieldU(:,j,1))
                        dUz_deta = sum(fieldU(i,:,1) * Tdomain%sSubdomain(mat)%hprimez(:,j))

                        invgrad_ij = Tdomain%specel(n_el)%InvGrad(i,j,:,:)
                        DXX = invgrad_ij(0,0)*dUx_dxi + invgrad_ij(0,1)*dUx_deta
                        DZZ = invgrad_ij(1,0)*dUz_dxi + invgrad_ij(1,1)*dUz_deta
                        eps_v = DXX + DZZ
                        eps_vol(i,j) = eps_v
                        fieldP(i,j) = -Tdomain%specel(n_el)%Lambda(i,j) * eps_v
                    endif
                else
                    ! Cas Solide (Élastique)
                    dUx_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * fieldU(:,j,0))
                    dUx_deta = sum(fieldU(i,:,0) * Tdomain%sSubdomain(mat)%hprimez(:,j))
                    dUz_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * fieldU(:,j,1))
                    dUz_deta = sum(fieldU(i,:,1) * Tdomain%sSubdomain(mat)%hprimez(:,j))

                    invgrad_ij = Tdomain%specel(n_el)%InvGrad(i,j,:,:)
                    DXX = invgrad_ij(0,0)*dUx_dxi + invgrad_ij(0,1)*dUx_deta
                    DXZ = invgrad_ij(1,0)*dUx_dxi + invgrad_ij(1,1)*dUx_deta
                    DZX = invgrad_ij(0,0)*dUz_dxi + invgrad_ij(0,1)*dUz_deta
                    DZZ = invgrad_ij(1,0)*dUz_dxi + invgrad_ij(1,1)*dUz_deta

                    dUdX(i,j,0) = DXX
                    dUdX(i,j,1) = DXZ
                    dUdX(i,j,2) = DZX
                    dUdX(i,j,3) = DZZ

                    eps_xx = DXX
                    eps_zz = DZZ
                    eps_xz = 0.5_fpp * (DXZ + DZX)
                    eps_v = eps_xx + eps_zz
                    eps_vol(i,j) = eps_v

                    eps_dev(i,j,0) = eps_xx - 0.5_fpp * eps_v
                    eps_dev(i,j,1) = eps_zz - 0.5_fpp * eps_v
                    eps_dev(i,j,2) = eps_xz

                    if (allocated(Tdomain%specel(n_el)%Cij2d)) then
                        strain_v = (/ eps_xx, eps_zz, 2.0_fpp * eps_xz /)
                        stress_v = matmul(Tdomain%specel(n_el)%Cij2d(:,:,i,j), strain_v)
                        sig_xx = stress_v(1)
                        sig_zz = stress_v(2)
                        sig_xz = stress_v(3)
                    else
                        sig_xx = Tdomain%specel(n_el)%Lambda(i,j) * eps_v + 2.0_fpp * Tdomain%specel(n_el)%Mu(i,j) * eps_xx
                        sig_zz = Tdomain%specel(n_el)%Lambda(i,j) * eps_v + 2.0_fpp * Tdomain%specel(n_el)%Mu(i,j) * eps_zz
                        sig_xz = 2.0_fpp * Tdomain%specel(n_el)%Mu(i,j) * eps_xz
                    endif

                    sig_mean = 0.5_fpp * (sig_xx + sig_zz)
                    sig_dev(i,j,0) = sig_xx - sig_mean
                    sig_dev(i,j,1) = sig_zz - sig_mean
                    sig_dev(i,j,2) = sig_xz

                    fieldP(i,j) = -sig_mean
                    P_energy(i,j) = 0.5_fpp * (sig_xx * eps_xx + sig_zz * eps_zz + 2.0_fpp * sig_xz * eps_xz)
                endif
            end do
        end do

        ! Interpolation à la position curviligne du capteur
        do i = 0,ngllx - 1
            do j = 0,ngllz - 1
                weight = outx(i)*outz(j)

                if (out_variables(OUT_DEPLA) == 1) then
                    ioff = offset(OUT_DEPLA)
                    nComp = OUT_VAR_DIMS_2D(OUT_DEPLA)-1
                    grandeur(ioff:ioff+nComp) = grandeur(ioff:ioff+nComp) + weight*fieldU(i,j,:)
                end if

                if (out_variables(OUT_VITESSE) == 1) then
                    ioff = offset(OUT_VITESSE)
                    nComp = OUT_VAR_DIMS_2D(OUT_VITESSE)-1
                    grandeur(ioff:ioff+nComp) = grandeur(ioff:ioff+nComp) + weight*fieldV(i,j,:)
                end if

                if (out_variables(OUT_ACCEL) == 1) then
                    ioff = offset(OUT_ACCEL)
                    nComp = OUT_VAR_DIMS_2D(OUT_ACCEL)-1
                    grandeur(ioff:ioff+nComp) = grandeur(ioff:ioff+nComp) + weight*fieldA(i,j,:)
                end if

                if (out_variables(OUT_PRESSION) == 1) then
                    ioff = offset(OUT_PRESSION)
                    nComp = OUT_VAR_DIMS_2D(OUT_PRESSION)-1
                    grandeur(ioff+nComp) = grandeur(ioff+nComp) + weight*fieldP(i,j)
                end if

                if (out_variables(OUT_ENERGYP) == 1) then
                    ioff = offset(OUT_ENERGYP)
                    nComp = OUT_VAR_DIMS_2D(OUT_ENERGYP)-1
                    grandeur (ioff+nComp) = grandeur (ioff+nComp) + weight*P_energy(i,j)
                end if

                if (out_variables(OUT_ENERGYK) == 1) then
                    ioff = offset(OUT_ENERGYK)
                    nComp = OUT_VAR_DIMS_2D(OUT_ENERGYK)-1
                    grandeur (ioff+nComp) = grandeur (ioff+nComp) + weight*K_energy(i,j)
                end if

                if (out_variables(OUT_DUDX) == 1) then
                    ioff = offset(OUT_DUDX)
                    nComp = OUT_VAR_DIMS_2D(OUT_DUDX)-1
                    grandeur (ioff:ioff+nComp) = grandeur (ioff:ioff+nComp)+weight*dUdX(i,j,:)
                end if

                if (out_variables(OUT_EPS_VOL) == 1) then
                    ioff = offset(OUT_EPS_VOL)
                    nComp = OUT_VAR_DIMS_2D(OUT_EPS_VOL)-1
                    grandeur (ioff+nComp) = grandeur (ioff+nComp) + weight*eps_vol(i,j)
                end if

                if (out_variables(OUT_EPS_DEV) == 1) then
                    ioff = offset(OUT_EPS_DEV)
                    nComp = OUT_VAR_DIMS_2D(OUT_EPS_DEV)-1
                    grandeur (ioff:ioff+nComp) = grandeur(ioff:ioff+nComp)+weight*eps_dev(i,j,:)
                end if

                if (out_variables(OUT_EPS_DEV_PL) == 1) then
                    ioff=offset(OUT_EPS_DEV_PL)
                    nComp = OUT_VAR_DIMS_2D(OUT_EPS_DEV_PL)-1
                    grandeur (ioff:ioff+nComp) = grandeur(ioff:ioff+nComp)+weight*eps_dev_pl(i,j,:)
                end if

                if (out_variables(OUT_STRESS_DEV) == 1) then
                    ioff = offset(OUT_STRESS_DEV)
                    nComp = OUT_VAR_DIMS_2D(OUT_STRESS_DEV)-1
                    grandeur (ioff:ioff+nComp) = grandeur(ioff:ioff+nComp)+weight*sig_dev(i,j,:)
                end if

                if (out_variables(OUT_ENERGYD) == 1) then
                    ioff = offset(OUT_ENERGYD)
                    nComp = OUT_VAR_DIMS_2D(OUT_ENERGYD)-1
                    grandeur(ioff:ioff+nComp) = grandeur (ioff:ioff+nComp) + weight*D_energy(i,j,:)
                end if
            end do
        end do

        ! Sauvegarde dans le cache
        r = capteur%icache+1
        capteur%valuecache(1,r) = Tdomain%TimeD%rtime
        capteur%valuecache(2:Tdomain%nReqOut+1,r) = grandeur(:)
        capteur%icache = r

        deallocate(fieldU)
        deallocate(fieldV)
        deallocate(fieldA)
        deallocate(fieldP)
        deallocate(P_energy)
        deallocate(K_energy)
        deallocate(D_energy)
        deallocate(eps_vol)
        deallocate(eps_dev)
        deallocate(eps_dev_pl)
        deallocate(sig_dev)
        deallocate(grandeur)
        deallocate(dUdX)
        deallocate(outx)
        deallocate(outz)
    end subroutine sortieGrandeurCapteur_interp

    subroutine integrate_on_element(ngllx, ngllz, jac, GLLwx, GLLwz, input_field, output_integral)
        implicit none
        integer, intent(in) :: ngllx, ngllz
        real(fpp), dimension(0:ngllx-1), intent(in) :: GLLwx
        real(fpp), dimension(0:ngllz-1), intent(in) :: GLLwz
        real(fpp), dimension(0:ngllx-1,0:ngllz-1), intent(in) :: jac
        real(fpp), dimension(0:ngllx-1,0:ngllz-1), intent(in) :: input_field
        real(fpp), intent(out) :: output_integral
        integer :: i, j
        output_integral = 0d0
        do j = 0, ngllz-1
            do i = 0, ngllx-1
                output_integral = output_integral + input_field(i,j) * jac(i,j) * GLLwx(i) * GLLwz(j)
            end do
        end do
    end subroutine integrate_on_element

    subroutine sortieGrandeurCapteur_energy(Tdomain, capteur)
        use sdomain
        implicit none
        type(domain)   :: TDomain
        type(tCapteur) :: capteur
        !
        integer :: n, mat, ngllx, ngllz
        real(fpp), dimension(:,:,:), allocatable   :: fieldU, fieldV
        real(fpp), dimension(:,:), allocatable     :: P_energy, K_energy, eps_vol, L_energy, S_energy, R_energy
        real(fpp), dimension(:,:,:), allocatable   :: D_energy
        real(fpp), dimension(:,:,:), allocatable   :: eps_dev
        real(fpp), dimension(:,:,:), allocatable   :: sig_dev
        real(fpp) :: local_sum_P_energy, local_sum_K_energy
        real(fpp) :: local_sum_L_energy, local_sum_S_energy, local_sum_R_energy
        real(fpp) :: global_sum_P_energy, global_sum_K_energy
        real(fpp) :: global_sum_L_energy, global_sum_S_energy, global_sum_R_energy
        real(fpp) :: elem_P_En, elem_K_En, elem_L_En, elem_S_En, elem_R_En
        real(fpp), dimension(0:1) :: elem_D_En
        type(element), pointer :: el
        integer :: i, j, ierr
        real(fpp) :: dUx_dxi, dUx_deta, dUz_dxi, dUz_deta
        real(fpp) :: DXX, DXZ, DZX, DZZ
        real(fpp) :: eps_xx, eps_zz, eps_xz, eps_v
        real(fpp) :: sig_xx, sig_zz, sig_xz, sig_mean
        real(fpp), dimension(3) :: strain_v, stress_v
        real(fpp), dimension(0:1,0:1) :: invgrad_ij

        if(capteur%type /= CPT_ENERGY) return
        local_sum_P_energy = 0d0
        local_sum_K_energy = 0d0
        local_sum_L_energy = 0d0
        local_sum_S_energy = 0d0
        local_sum_R_energy = 0d0

        do n = 0,Tdomain%n_elem-1
            el => Tdomain%specel(n)
            if (el%PML) cycle ! Pas d'énergie dans les PMLs

            ngllx = el%ngllx
            ngllz = el%ngllz
            mat   = el%mat_index

            allocate(fieldU(0:ngllx-1,0:ngllz-1,0:1))
            allocate(fieldV(0:ngllx-1,0:ngllz-1,0:1))
            allocate(P_energy(0:ngllx-1,0:ngllz-1))
            allocate(K_energy(0:ngllx-1,0:ngllz-1))
            allocate(L_energy(0:ngllx-1,0:ngllz-1))
            allocate(S_energy(0:ngllx-1,0:ngllz-1))
            allocate(R_energy(0:ngllx-1,0:ngllz-1))

            if (Tdomain%type_timeInteg == TIME_INTEG_RK4) then
                fieldU = el%Displ
                fieldV = el%Veloc
            else
                call gather_elem_displ(Tdomain, n, fieldU)
                call gather_elem_veloc(Tdomain, n, fieldV, .false.)
            endif

            L_energy = 0.0_fpp
            S_energy = 0.0_fpp
            R_energy = 0.0_fpp

            do j = 0,ngllz-1
                do i = 0,ngllx-1
                    K_energy(i,j) = 0.5_fpp * el%Density(i,j) * (fieldV(i,j,0)**2 + fieldV(i,j,1)**2)

                    if (el%acoustic) then
                        P_energy(i,j) = 0d0 ! Pas d'énergie potentielle standard
                        L_energy(i,j) = 0d0
                        S_energy(i,j) = 0d0
                        R_energy(i,j) = 0d0
                    else
                        dUx_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * fieldU(:,j,0))
                        dUx_deta = sum(fieldU(i,:,0) * Tdomain%sSubdomain(mat)%hprimez(:,j))
                        dUz_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * fieldU(:,j,1))
                        dUz_deta = sum(fieldU(i,:,1) * Tdomain%sSubdomain(mat)%hprimez(:,j))

                        invgrad_ij = el%InvGrad(i,j,:,:)
                        DXX = invgrad_ij(0,0)*dUx_dxi + invgrad_ij(0,1)*dUx_deta
                        DXZ = invgrad_ij(1,0)*dUx_dxi + invgrad_ij(1,1)*dUx_deta
                        DZX = invgrad_ij(0,0)*dUz_dxi + invgrad_ij(0,1)*dUz_deta
                        DZZ = invgrad_ij(1,0)*dUz_dxi + invgrad_ij(1,1)*dUz_deta

                        eps_xx = DXX
                        eps_zz = DZZ
                        eps_xz = 0.5_fpp * (DXZ + DZX)
                        eps_v = eps_xx + eps_zz

                        if (allocated(el%Cij2d)) then
                            strain_v = (/ eps_xx, eps_zz, 2.0_fpp * eps_xz /)
                            stress_v = matmul(el%Cij2d(:,:,i,j), strain_v)
                            sig_xx = stress_v(1)
                            sig_zz = stress_v(2)
                            sig_xz = stress_v(3)
                        else
                            sig_xx = el%Lambda(i,j) * eps_v + 2.0_fpp * el%Mu(i,j) * eps_xx
                            sig_zz = el%Lambda(i,j) * eps_v + 2.0_fpp * el%Mu(i,j) * eps_zz
                            sig_xz = 2.0_fpp * el%Mu(i,j) * eps_xz
                        endif

                        L_energy(i,j) = el%Mu(i,j)/2.0_fpp * (DXZ - DZX)**2
                        S_energy(i,j) = (0.5_fpp * el%Lambda(i,j) + el%Mu(i,j)) * eps_v**2
                        R_energy(i,j) = 2.0_fpp * el%Mu(i,j) * DXZ * DZX - 2.0_fpp * el%Mu(i,j) * eps_xx * eps_zz

                        P_energy(i,j) = 0.5_fpp * (sig_xx * eps_xx + sig_zz * eps_zz + 2.0_fpp * sig_xz * eps_xz)
                    endif
                end do
            end do

            call integrate_on_element(ngllx, ngllz, el%Jacob, Tdomain%sSubdomain(mat)%GLLwx, Tdomain%sSubdomain(mat)%GLLwz, P_energy, elem_P_En)
            call integrate_on_element(ngllx, ngllz, el%Jacob, Tdomain%sSubdomain(mat)%GLLwx, Tdomain%sSubdomain(mat)%GLLwz, K_energy, elem_K_En)
            call integrate_on_element(ngllx, ngllz, el%Jacob, Tdomain%sSubdomain(mat)%GLLwx, Tdomain%sSubdomain(mat)%GLLwz, L_energy, elem_L_En)
            call integrate_on_element(ngllx, ngllz, el%Jacob, Tdomain%sSubdomain(mat)%GLLwx, Tdomain%sSubdomain(mat)%GLLwz, S_energy, elem_S_En)
            call integrate_on_element(ngllx, ngllz, el%Jacob, Tdomain%sSubdomain(mat)%GLLwx, Tdomain%sSubdomain(mat)%GLLwz, R_energy, elem_R_En)

            local_sum_P_energy = local_sum_P_energy + elem_P_En
            local_sum_K_energy = local_sum_K_energy + elem_K_En
            local_sum_L_energy = local_sum_L_energy + elem_L_En
            local_sum_S_energy = local_sum_S_energy + elem_S_En
            local_sum_R_energy = local_sum_R_energy + elem_R_En

            deallocate(fieldU, fieldV, P_energy, K_energy, L_energy, S_energy, R_energy)
        enddo

        call MPI_AllReduce(local_sum_P_energy, global_sum_P_energy, 1, MPI_DOUBLE, MPI_SUM, Tdomain%communicateur, ierr)
        call MPI_AllReduce(local_sum_K_energy, global_sum_K_energy, 1, MPI_DOUBLE, MPI_SUM, Tdomain%communicateur, ierr)
        call MPI_AllReduce(local_sum_L_energy, global_sum_L_energy, 1, MPI_DOUBLE, MPI_SUM, Tdomain%communicateur, ierr)
        call MPI_AllReduce(local_sum_S_energy, global_sum_S_energy, 1, MPI_DOUBLE, MPI_SUM, Tdomain%communicateur, ierr)
        call MPI_AllReduce(local_sum_R_energy, global_sum_R_energy, 1, MPI_DOUBLE, MPI_SUM, Tdomain%communicateur, ierr)

        if(Tdomain%Mpi_var%my_rank == 0) then
            i = capteur%icache+1
            capteur%valuecache(1,i) = Tdomain%TimeD%rtime
            capteur%valuecache(2,i) = global_sum_P_energy
            capteur%valuecache(3,i) = global_sum_K_energy
            capteur%valuecache(4,i) = global_sum_L_energy
            capteur%valuecache(5,i) = global_sum_S_energy
            capteur%valuecache(6,i) = global_sum_R_energy
            capteur%icache = i
        endif

    end subroutine sortieGrandeurCapteur_energy


    subroutine trouve_capteur(Tdomain, xc, zc, n_el, dmin, xi, eta, flag)
        use shape_lin
        use shape_quad
        use mlocations2d
        implicit none
        type (domain), INTENT(INOUT)  :: Tdomain
        real(fpp), intent(in) :: xc, zc
        integer, intent(out) :: n_el
        real(fpp), intent(out) :: xi, eta, dmin
        logical, intent(in) :: flag
        !
        integer :: i, iel, emin, j
        logical :: inside
        integer :: nmax
        integer, parameter :: NMAXEL=100
        integer, dimension(NMAXEL) :: elems
        real(fpp), dimension(0:1,NMAXEL) :: coordloc
        real(fpp), parameter :: EPSD = 1D-10
        real(fpp) :: EPS
        real(fpp), dimension(0:1, 0:Tdomain%n_nodes-1) :: coordl
        real(fpp) :: dist, xc0, zc0

        EPS = EPSD
        nmax = NMAXEL
        call find_location(Tdomain, xc, zc, nmax, elems, coordloc)
        n_el = -1
        dmin = 1e20
        emin = -1
        do i=1,nmax
            xi   = coordloc(0,i)
            eta  = coordloc(1,i)
            iel = elems(i)
            do j = 0, Tdomain%n_nodes-1
                coordl(0:1, j) = Tdomain%Coord_Nodes(0:1, Tdomain%specel(iel)%Control_Nodes(j))
            enddo
            inside = .true.
            if ((xi<-1_fpp) .or. eta<(-1_fpp)) inside = .false.
            if ((xi>+1_fpp) .or. eta>(+1_fpp)) inside = .false.
            ! On projette sur le bord de l'element si besoin
            if (xi<-1D0) xi = -1D0
            if (xi>+1D0) xi = +1D0
            if (eta<-1D0) eta = -1D0
            if (eta>+1D0) eta = +1D0
            ! calcule less coord du nouveau pts
            if (.not. inside) then
                if (Tdomain%n_nodes==4) then
                    call shape4_local2global(coordl, xi, eta, xc0, zc0)
                else
                    call shape8_local2global(coordl, xi, eta, xc0, zc0)
                end if
                dist = (xc-xc0)**2 + (zc-zc0)**2
            else
                dist = 0.
            endif
            if (dist<dmin) then
                dmin = dist
                emin = iel
            endif
        end do
        n_el = emin
        return
    end subroutine trouve_capteur

end module mCapteur
