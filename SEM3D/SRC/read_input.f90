! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module semconfig
    use semdatafiles
    use mesh3d
    use iso_c_binding
    use sem_c_config

    implicit none

    public :: read_input

contains

    subroutine read_material_file(Tdomain)
        use sdomain
        use orientation
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: i, mat, dom

        call read_material_file_v2(Tdomain)
        ! Complete the material definition with optional material.spec
        call read_material_spec(Tdomain)

        Tdomain%any_sdom = .false.
        Tdomain%any_fdom = .false.
        Tdomain%any_spml = .false.
        Tdomain%any_fpml = .false.
        Tdomain%any_sdomdg = .false.
        Tdomain%sdom%ngll    = 0
        Tdomain%sdomdg%ngll  = 0
        Tdomain%fdom%ngll    = 0
        Tdomain%spmldom%ngll = 0
        Tdomain%fpmldom%ngll = 0

        do mat = 0, Tdomain%n_mat-1
            dom = Tdomain%sSubDomain(mat)%dom
            select case (dom)
            case (DM_SOLID_CG)
                Tdomain%sdom%ngll = Tdomain%sSubDomain(mat)%NGLL
                Tdomain%any_sdom = .true.
            case (DM_SOLID_DG)
                Tdomain%sdomdg%ngll = Tdomain%sSubDomain(mat)%NGLL
                Tdomain%any_sdomdg = .true.
            case (DM_FLUID_CG)
                Tdomain%fdom%ngll = Tdomain%sSubDomain(mat)%NGLL
                Tdomain%any_fdom = .true.
            case (DM_SOLID_CG_PML)
                Tdomain%spmldom%ngll = Tdomain%sSubDomain(mat)%NGLL
                Tdomain%any_spml = .true.
            case (DM_FLUID_CG_PML)
                Tdomain%fpmldom%ngll = Tdomain%sSubDomain(mat)%NGLL
                Tdomain%any_fpml = .true.
            case default
                stop " Fatal Error : unknown domain"
            end select
        end do
        !- GLL properties in elements, on faces, edges.
        do i = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(i)%mat_index
            Tdomain%specel(i)%domain = Tdomain%sSubDomain(mat)%dom
        end do

        call apply_mat_to_faces(Tdomain)
        call apply_mat_to_edges(Tdomain)
        call apply_mat_to_vertices(Tdomain)
        call apply_interface(Tdomain, Tdomain%intSolPml, DM_SOLID_CG, DM_SOLID_CG_PML, .false.)
        call apply_interface(Tdomain, Tdomain%intFluPml, DM_FLUID_CG, DM_FLUID_CG_PML, .false.)
        call apply_interface(Tdomain, Tdomain%SF%intSolFlu, DM_SOLID_CG, DM_FLUID_CG, .false.)
        call apply_interface(Tdomain, Tdomain%SF%intSolFluPml, DM_SOLID_CG_PML, DM_FLUID_CG_PML, .false.)
        call apply_interface(Tdomain, Tdomain%intSolPml, DM_SOLID_CG, DM_SOLID_CG_PML, .true.)
        call apply_interface(Tdomain, Tdomain%intFluPml, DM_FLUID_CG, DM_FLUID_CG_PML, .true.)
        call apply_interface(Tdomain, Tdomain%SF%intSolFlu, DM_SOLID_CG, DM_FLUID_CG, .true.)
        call apply_interface(Tdomain, Tdomain%SF%intSolFluPml, DM_SOLID_CG_PML, DM_FLUID_CG_PML, .true.)

    end subroutine read_material_file

    subroutine read_material_spec(Tdomain)
        use sdomain
        use semdatafiles
        use iso_c_binding, only: c_loc, c_associated
        implicit none
        type(domain), intent(inout)                  :: Tdomain
        !
        type(sem_material), pointer :: matdesc
        integer :: i, code, num, nprop

        call read_sem_materials(Tdomain%material_list, Tdomain%rank, "material.spec"//C_NULL_CHAR, code)

        call c_f_pointer(Tdomain%material_list%head, matdesc)

        do while(associated(matdesc))
            num = matdesc%num
            ! init global
            Tdomain%sSubdomain(num)%is_sph = .false.
            Tdomain%sSubdomain(num)%n_prop = 0
            select case(matdesc%defspatial)
            case (0) ! CONSTANT
                Tdomain%sSubdomain(num)%material_definition = MATERIAL_CONSTANT
            case (1) ! FILE
                Tdomain%sSubdomain(num)%material_definition = MATERIAL_FILE
                nprop = 3
                select case(matdesc%deftype)
                case(MATDEF_VP_VS_RHO_D,MATDEF_E_NU_RHO_D,MATDEF_LAMBDA_MU_RHO_D,&
                    MATDEF_KAPPA_MU_RHO_D,MATDEF_NLKP_VS_RHO_D,&
                    MATDEF_NU_VS_RHO_D)
                    nprop = 5
                case (MATDEF_VTI_ANISO)
                    Tdomain%aniso=.true.
                    nprop = 8
                case (MATDEF_HOOKE_ANISO)
                    nprop = 22
                    Tdomain%aniso=.true.
                end select

                ! check for spherical material
                if (matdesc%is_sph==1) then
                    Tdomain%sSubdomain(num)%is_sph = .true.
                    Tdomain%sSubdomain(num)%sph_args%theta_chk = 90.0_fpp-matdesc%lat_center
                    Tdomain%sSubdomain(num)%sph_args%phi_chk = matdesc%lon_center
                end if
                Tdomain%sSubdomain(num)%n_prop = nprop
                allocate(Tdomain%sSubdomain(num)%prop_field(nprop))
                if (nprop==3) then
                    Tdomain%sSubdomain(num)%prop_field(1)%propFilePath = trim(fromcstr(matdesc%filename0))
                    Tdomain%sSubdomain(num)%prop_field(2)%propFilePath = trim(fromcstr(matdesc%filename1))
                    Tdomain%sSubdomain(num)%prop_field(3)%propFilePath = trim(fromcstr(matdesc%filename2))
                elseif (nprop==5) then
                    Tdomain%sSubdomain(num)%prop_field(1)%propFilePath = trim(fromcstr(matdesc%filename0))
                    Tdomain%sSubdomain(num)%prop_field(2)%propFilePath = trim(fromcstr(matdesc%filename1))
                    Tdomain%sSubdomain(num)%prop_field(3)%propFilePath = trim(fromcstr(matdesc%filename2))
                    Tdomain%sSubdomain(num)%prop_field(4)%propFilePath = trim(fromcstr(matdesc%filename3))
                    Tdomain%sSubdomain(num)%prop_field(5)%propFilePath = trim(fromcstr(matdesc%filename4))
                elseif (nprop>5) then
                    do i = 1,nprop
                        Tdomain%sSubdomain(num)%prop_field(i)%propFilePath = trim(fromcstr(matdesc%filename0))
                    end do
                endif
               !MC MC
               ! do i = 6,nprop
               !     Tdomain%sSubdomain(num)%prop_field(i)%propFilePath = trim(fromcstr(matdesc%filename0))
               ! end do
            end select
            Tdomain%sSubdomain(num)%deftype  = matdesc%deftype
            Tdomain%sSubdomain(num)%Ddensity = matdesc%rho
            Tdomain%sSubdomain(num)%Pspeed   = matdesc%Vp
            Tdomain%sSubdomain(num)%Sspeed   = matdesc%Vs
            Tdomain%sSubdomain(num)%DE       = matdesc%E
            Tdomain%sSubdomain(num)%DNu      = matdesc%nu
            Tdomain%sSubdomain(num)%DLambda  = matdesc%lambda
            Tdomain%sSubdomain(num)%DMu      = matdesc%mu
            Tdomain%sSubdomain(num)%DKappa   = matdesc%kappa
            Tdomain%sSubdomain(num)%DRinf    = matdesc%rinf
            Tdomain%sSubdomain(num)%DBiso    = matdesc%biso
            Tdomain%sSubdomain(num)%DNlkp    = matdesc%nlkp
            Tdomain%sSubdomain(num)%Qp       = matdesc%Qp
            Tdomain%sSubdomain(num)%Qs       = matdesc%Qs

            !!! TODO Add Qp/Qmu, PML
            call c_f_pointer(matdesc%next, matdesc)
        end do
        if (code<=0) then
            write(*,*) "Erreur reading material.spec"
            stop 1
        endif
    end subroutine read_material_spec

    subroutine read_material_file_v2(Tdomain)
        use sdomain
        use semdatafiles
        use sem_mpi
        implicit none

        type(domain), intent(inout)   :: Tdomain
        character(Len=MAX_FILE_SIZE)  :: fnamef, buffer
        integer                       :: i, n_aus, npml
        integer                       :: rg, fid, assocMat
        character                     :: material_type

        rg = Tdomain%rank
        npml = 0
        fid = 13

        call semname_read_inputmesh_parametrage(Tdomain%material_file,fnamef)
        if(rg==0) write(*,*) "read material file : ", trim(adjustL(fnamef))
        open (fid, file=fnamef, status="old", form="formatted")

        buffer = getLine (fid, "#")
        read(buffer,*) n_aus

        if(n_aus /= Tdomain%n_mat) then
            write(*,*) trim(fnamef), n_aus, Tdomain%n_mat
            stop "Incompatibility between the mesh file and the material file for n_mat"
        endif

        !!! GB GB
        !!! if (Tdomain%aniso) then
        !!!     print *,"The code can't put anisotropy in a homogeneous media"
        !!!     stop
        !!! endif

        do i = 0,Tdomain%n_mat-1

            buffer = getLine (fid, "#")
            read(buffer,*) material_type, &
                Tdomain%sSubDomain(i)%Pspeed,               &
                Tdomain%sSubDomain(i)%Sspeed,               &
                Tdomain%sSubDomain(i)%dDensity,             &
                Tdomain%sSubDomain(i)%Qpression,            &
                Tdomain%sSubDomain(i)%Qmu

            Tdomain%sSubDomain(i)%NGLL = Tdomain%ngll
            Tdomain%sSubDomain(i)%dom = domain_from_type_char(material_type)
            Tdomain%sSubdomain(i)%material_definition = MATERIAL_CONSTANT
            Tdomain%sSubDomain(i)%deftype = MATDEF_VP_VS_RHO
            Tdomain%sSubDomain(i)%pml_width = 0.
            Tdomain%sSubDomain(i)%is_sph = .false.

            ! call Lame_coefficients (Tdomain%sSubDomain(i)) ! XXX Make sure its done in definearrays

            if (is_pml(Tdomain%sSubDomain(i)))  then
                npml = npml + 1
            endif
        enddo

        if(npml > 0) then
            do i = 0,Tdomain%n_mat-1
                if (is_pml(Tdomain%sSubDomain(i)))  then
                    buffer = getLine (fid, "#")
                    read(buffer,*) Tdomain%sSubdomain(i)%npow,  &
                        Tdomain%sSubdomain(i)%Apow,         &
                        Tdomain%sSubdomain(i)%pml_pos(0), &
                        Tdomain%sSubdomain(i)%pml_width(0), &
                        Tdomain%sSubdomain(i)%pml_pos(1), &
                        Tdomain%sSubdomain(i)%pml_width(1), &
                        Tdomain%sSubdomain(i)%pml_pos(2), &
                        Tdomain%sSubdomain(i)%pml_width(2), &
                        assocMat
                endif
            enddo
        endif

        close(13)

    end subroutine read_material_file_v2


    subroutine create_sem_sources(Tdomain, config)
        use sdomain
        implicit none
        type(domain), intent(inout)  :: Tdomain
        type(sem_config), intent(in) :: config
        type(sem_source), pointer :: src
        integer :: nsrc
        real(fpp) :: ndir

        call c_f_pointer(config%source, src)
        nsrc = 0
        do while(associated(src))
            Tdomain%Ssource(nsrc)%Xsource = src%coords(1)
            Tdomain%Ssource(nsrc)%Ysource = src%coords(2)
            Tdomain%Ssource(nsrc)%Zsource = src%coords(3)
            Tdomain%Ssource(nsrc)%i_type_source = src%type
            Tdomain%Ssource(nsrc)%amplitude_factor = src%amplitude
            if (src%func .eq. 5) then
                Tdomain%Ssource(nsrc)%time_file = trim(fromcstr(src%time_file))
            end if
            if (src%func .eq. 15) then
                Tdomain%Ssource(nsrc)%time_file = trim(fromcstr(src%time_file))
            endif
            ! Comportement temporel
            Tdomain%Ssource(nsrc)%i_time_function = src%func
            Tdomain%Ssource(nsrc)%cutoff_freq = src%freq ! func=2,4
            Tdomain%Ssource(nsrc)%tau_b = src%tau ! func=1,2,3,4,5
            Tdomain%Ssource(nsrc)%fh = src%band  ! func=3
            Tdomain%Ssource(nsrc)%gamma = src%gamma ! func=4
            Tdomain%Ssource(nsrc)%ts = src%ts   ! func=4
            Tdomain%Ssource(nsrc)%Q = src%Q
            Tdomain%Ssource(nsrc)%Y = src%Y
            Tdomain%Ssource(nsrc)%X = src%X
            Tdomain%Ssource(nsrc)%L = src%L
            Tdomain%Ssource(nsrc)%v = src%v
            Tdomain%Ssource(nsrc)%a = src%a
            Tdomain%Ssource(nsrc)%d = src%d

            ! Comportement Spacial
            ! i_type_source==1
            ndir = sqrt(src%dir(1)**2 + src%dir(2)**2 + src%dir(3)**2)
            if (ndir>0) Tdomain%Ssource(nsrc)%dir(0:2) = src%dir(1:3)/ndir

            ! i_type_source==2
            Tdomain%Ssource(nsrc)%Moment(0,0) = src%moments(1)
            Tdomain%Ssource(nsrc)%Moment(1,1) = src%moments(2)
            Tdomain%Ssource(nsrc)%Moment(2,2) = src%moments(3)
            Tdomain%Ssource(nsrc)%Moment(0,1) = src%moments(4)
            Tdomain%Ssource(nsrc)%Moment(1,0) = src%moments(4)
            Tdomain%Ssource(nsrc)%Moment(0,2) = src%moments(5)
            Tdomain%Ssource(nsrc)%Moment(2,0) = src%moments(5)
            Tdomain%Ssource(nsrc)%Moment(1,2) = src%moments(6)
            Tdomain%Ssource(nsrc)%Moment(2,1) = src%moments(6)

            nsrc = nsrc + 1
            Tdomain%logicD%any_source = .true.
            call c_f_pointer(src%next, src)

        end do


    end subroutine create_sem_sources


    subroutine create_sem_extended_sources(Tdomain, config)
        use sdomain
        use sem_hdf5
        implicit none
        type(domain), intent(inout)  :: Tdomain
        type(sem_config), intent(in) :: config
        type(sem_extended_source), pointer :: ext_src
        integer :: nsrc, ntotal, nextsrc, hdferr
        integer(HID_T) :: fid

        nsrc = config%nsources ! Number of isolated point sources
        Tdomain%n_extsource = config%nextended_sources
        allocate (Tdomain%sExtendSource(0:Tdomain%n_extsource-1))
        call init_hdf5()

        ! Calcul nombre total de sources pour allouer Tdomain%Ssource
        call c_f_pointer(config%extended_source, ext_src)
        ntotal = nsrc ; nextsrc = 0
        do while(associated(ext_src))
            Tdomain%sExtendSource(nextsrc)%kine_file = trim(fromcstr(ext_src%kine_file))
            Tdomain%sExtendSource(nextsrc)%slip_file = trim(fromcstr(ext_src%slip_file))
            Tdomain%sExtendSource(nextsrc)%is_force = .false.
            if(ext_src%is_force == 1) Tdomain%sExtendSource(nextsrc)%is_force = .true.
            call h5fopen_f(Tdomain%sExtendSource(nextsrc)%kine_file, H5F_ACC_RDONLY_F, fid, hdferr)
            call read_attr_int (fid, "Ns",    Tdomain%sExtendSource(nextsrc)%Ns)
            call read_attr_int (fid, "Nd",    Tdomain%sExtendSource(nextsrc)%Nd)
            call read_attr_int (fid, "Nt",    Tdomain%sExtendSource(nextsrc)%Nt)
            call read_attr_real(fid, "dt",    Tdomain%sExtendSource(nextsrc)%Dt)

            Tdomain%sExtendSource(nextsrc)%Npt = Tdomain%sExtendSource(nextsrc)%Ns * Tdomain%sExtendSource(nextsrc)%Nd
            ntotal = ntotal + Tdomain%sExtendSource(nextsrc)%Npt
            call h5fclose_f(fid, hdferr)
            call c_f_pointer(ext_src%next, ext_src)
            nextsrc = nextsrc+1
        end do

        ! Allocation du tableau de sources
        Tdomain%n_source = ntotal
        allocate (Tdomain%Ssource(0:Tdomain%n_source-1))
        if (ntotal>0)  Tdomain%logicD%any_source = .true.

        ! Creation en premier des sources individuelles
        call create_sem_sources(Tdomain, config)

        ! Nouvelle boucle pour creation des sources etendues
        nextsrc = 0
        call c_f_pointer(config%extended_source, ext_src)
        do while(associated(ext_src))
            call h5fopen_f(trim(fromcstr(ext_src%kine_file)), H5F_ACC_RDONLY_F, fid, hdferr)
            call create_point_sources_from_fault(Tdomain,Tdomain%sExtendSource(nextsrc),fid,nsrc)
            call h5fclose_f(fid, hdferr)
            call c_f_pointer(ext_src%next, ext_src)
            nextsrc = nextsrc+1
        end do

    end subroutine create_sem_extended_sources
    !


    subroutine create_point_sources_from_fault(Tdomain, extsrc, fid, nsrc)
        use sdomain
        use constants
        use sem_hdf5
        use shape_geom_3d
        implicit none
        type(domain), intent(inout)       :: Tdomain
        type(Extended_source), intent(inout) :: extsrc
        integer(HID_T), intent(in)        :: fid
        integer, intent(inout)            :: nsrc

        if (.NOT.  extsrc%is_force) &
            call create_point_sources_from_fault_moment(Tdomain, extsrc, fid, nsrc)

        if (extsrc%is_force) &
            call create_point_sources_from_fault_force(Tdomain, extsrc, fid, nsrc)



    end subroutine create_point_sources_from_fault
    !


    subroutine create_point_sources_from_fault_moment(Tdomain, extsrc, fid, nsrc)
        use sdomain
        use constants
        use sem_hdf5
        use shape_geom_3d
        implicit none
        type(domain), intent(inout)       :: Tdomain
        ! type(Extended_source), intent(in) :: extsrc
        type(Extended_source), intent(inout) :: extsrc
        integer(HID_T), intent(in)        :: fid
        integer, intent(inout)            :: nsrc
        real(fpp), dimension(0:2)              :: normal, U, v1, v2, dSn
        real(fpp), dimension(0:2,0:2)          :: Moment
        real(fpp), allocatable, dimension(:,:) :: tmp, Xtemp, Ytemp, Ztemp
        integer :: i, j



        ! Lire les vecteurs dans le fichier kine
        call read_attr_real_vec(fid, "Vnormal",extsrc%Normal(:))
        call read_attr_real_vec(fid, "Vslip",  extsrc%Uslip(:))



        ! Recuperation du vecteur normal et du vecteur de glissement
        normal = extsrc%Normal
        U      = extsrc%Uslip

        ! Construction du moment correspondant au slip
        do i=0,2
            do j=0,2
                Moment(i,j) = normal(i)*U(j) + normal(j)*U(i)
            enddo
        enddo


        ! Lecture des datasets, avec (x,y,z) repere direct.
        ! Allocation des tableaux de coordonnes
        allocate (Xtemp(0:extsrc%Ns-1,0:extsrc%Nd-1))
        allocate (Ytemp(0:extsrc%Ns-1,0:extsrc%Nd-1))
        allocate (Ztemp(0:extsrc%Ns-1,0:extsrc%Nd-1))

        call read_dataset(fid, "x", tmp)
        do i=0,extsrc%Ns-1
            Xtemp(i,:) = tmp(:,i+1)
        enddo
        deallocate(tmp)
        call read_dataset(fid, "y", tmp)
        do i=0,extsrc%Ns-1
            Ytemp(i,:) = tmp(:,i+1)
        enddo
        deallocate(tmp)
        call read_dataset(fid, "z", tmp)
        do i=0,extsrc%Ns-1
            Ztemp(i,:) = tmp(:,i+1)
        enddo
        deallocate(tmp)

        ! Surface elementaire
        v1(0) = Xtemp(1,0)-Xtemp(0,0) ; v2(0) = Xtemp(0,1)-Xtemp(0,0)
        v1(1) = Ytemp(1,0)-Ytemp(0,0) ; v2(1) = Ytemp(0,1)-Ytemp(0,0)
        v1(2) = Ztemp(1,0)-Ztemp(0,0) ; v2(2) = Ztemp(0,1)-Ztemp(0,0)
        call cross_prod(v1,v2,dSn)

        do i=0,extsrc%Ns-1
            do j=0,extsrc%Nd-1
                ! Source position
                Tdomain%Ssource(nsrc)%Xsource = Xtemp(i,j)
                Tdomain%Ssource(nsrc)%Ysource = Ytemp(i,j)
                Tdomain%Ssource(nsrc)%Zsource = Ztemp(i,j)

                ! Position in the grid
                Tdomain%Ssource(nsrc)%ind_i = i
                Tdomain%Ssource(nsrc)%ind_j = j

                ! Moment-type source for all the points
                Tdomain%Ssource(nsrc)%i_type_source = 2
                Tdomain%Ssource(nsrc)%amplitude_factor = 1.

                ! Comportement temporel
                Tdomain%Ssource(nsrc)%time_file = extsrc%slip_file ! Slip-rate history file
                Tdomain%Ssource(nsrc)%i_time_function = 15 ! flag pour source file...
                Tdomain%Ssource(nsrc)%ts = extsrc%Dt
                Tdomain%Ssource(nsrc)%Nt = extsrc%Nt

                Tdomain%Ssource(nsrc)%moment(:,:) = Moment(:,:)

                nsrc = nsrc+1
            enddo
        enddo

        deallocate(Xtemp, Ytemp, Ztemp)

    end subroutine create_point_sources_from_fault_moment
    !


    subroutine create_point_sources_from_fault_force(Tdomain, extsrc, fid, nsrc)

        use sdomain
        use constants
        use sem_hdf5
        use shape_geom_3d
        implicit none
        type(domain), intent(inout)            :: Tdomain
        type(Extended_source), intent(in)      :: extsrc
        integer(HID_T), intent(in)             :: fid
        integer, intent(inout)                 :: nsrc
        real(fpp), dimension(0:2)              :: v1, v2, dSn, dirvec
        real(fpp), allocatable, dimension(:,:) :: tmp, Xtemp, Ytemp, Ztemp
        integer :: i, j

        ! Lecture des datasets, avec (x,y,z) repere direct.
        ! Allocation des tableaux de coordonnes
        allocate (Xtemp(0:extsrc%Ns-1,0:extsrc%Nd-1))
        allocate (Ytemp(0:extsrc%Ns-1,0:extsrc%Nd-1))
        allocate (Ztemp(0:extsrc%Ns-1,0:extsrc%Nd-1))

        call read_dataset(fid, "x", tmp)
        do i=0,extsrc%Ns-1
            Xtemp(i,:) = tmp(:,i+1)
        enddo
        deallocate(tmp)
        call read_dataset(fid, "y", tmp)
        do i=0,extsrc%Ns-1
            Ytemp(i,:) = tmp(:,i+1)
        enddo
        deallocate(tmp)
        call read_dataset(fid, "z", tmp)
        do i=0,extsrc%Ns-1
            Ztemp(i,:) = tmp(:,i+1)
        enddo
        deallocate(tmp)

        ! Vecteur de direction
        call read_attr_real_vec(fid, "Dir", dirvec)


        ! Surface elementaire
        v1(0) = Xtemp(1,0)-Xtemp(0,0) ; v2(0) = Xtemp(0,1)-Xtemp(0,0)
        v1(1) = Ytemp(1,0)-Ytemp(0,0) ; v2(1) = Ytemp(0,1)-Ytemp(0,0)
        v1(2) = Ztemp(1,0)-Ztemp(0,0) ; v2(2) = Ztemp(0,1)-Ztemp(0,0)
        call cross_prod(v1,v2,dSn)

        do i=0,extsrc%Ns-1
            do j=0,extsrc%Nd-1
                ! Source position
                Tdomain%Ssource(nsrc)%Xsource = Xtemp(i,j)
                Tdomain%Ssource(nsrc)%Ysource = Ytemp(i,j)
                Tdomain%Ssource(nsrc)%Zsource = Ztemp(i,j)

                ! Position in the grid
                Tdomain%Ssource(nsrc)%ind_i = i
                Tdomain%Ssource(nsrc)%ind_j = j

                ! Moment-type source for all the points
                Tdomain%Ssource(nsrc)%i_type_source = 1            ! CHANGED
                Tdomain%Ssource(nsrc)%amplitude_factor = 1.0

                ! Comportement temporel
                Tdomain%Ssource(nsrc)%time_file = extsrc%slip_file ! Slip-rate history file
                Tdomain%Ssource(nsrc)%i_time_function = 15         ! flag pour source file...
                Tdomain%Ssource(nsrc)%ts = extsrc%Dt
                Tdomain%Ssource(nsrc)%Nt = extsrc%Nt

                Tdomain%Ssource(nsrc)%dir = dirvec


                nsrc = nsrc+1
            enddo
        enddo
        deallocate(Xtemp, Ytemp, Ztemp)

    end subroutine create_point_sources_from_fault_force
    !


    function is_in_box(pos, box)
        real(fpp), dimension(3), intent(in) :: pos
        double precision, dimension(6), intent(in) :: box
        logical :: is_in_box
        !
        integer :: i

        is_in_box = .true.
        do i=1,3
            if (pos(i)<box(i)) is_in_box = .false.
            if (pos(i)>box(3+i)) is_in_box = .false.
        end do
    end function is_in_box
    !
    ! Renvoie true si le point est a une distance <1 du plane
    function is_in_plane(pos, plane)
        real(fpp), dimension(3), intent(in) :: pos
        double precision, dimension(4), intent(in) :: plane
        logical :: is_in_plane
        !
        real(fpp) :: w

        is_in_plane = .false.
        w = plane(1)*pos(1)+plane(2)*pos(2)+plane(3)*pos(3)+plane(4)
        if (abs(w)<1) is_in_plane = .true.
    end function is_in_plane
    !>
    ! Selectionne les elements pour les inclure ou non dans les snapshots
    !<
    subroutine select_output_elements(Tdomain, config)
        use sdomain
        implicit none
        type(domain), intent(inout)  :: Tdomain
        type(sem_config), intent(in) :: config

        type(sem_snapshot_cond), pointer :: selection
        integer :: n, i, ipoint
        real(fpp), dimension(3) :: pos
        logical :: sel

        do n = 0, Tdomain%n_elem-1
            call c_f_pointer(config%snapshot_selection, selection)
            do while(associated(selection))
                if (selection%include==1) then
                    sel = .true.
                else
                    sel = .false.
                end if
                pos = 0.
                do i=0,7
                    ipoint = Tdomain%specel(n)%Control_Nodes(i)
                    pos = pos + Tdomain%Coord_Nodes(:,ipoint)
                end do
                pos = pos/8
                select case (selection%type)
                case (1)
                    ! All
                    Tdomain%specel(n)%output = sel
                case (2)
                    ! Material
                    if (Tdomain%specel(n)%mat_index == selection%material) Tdomain%specel(n)%output = sel
                case (3)
                    ! Box
                    if (is_in_box(pos, selection%box)) Tdomain%specel(n)%output = sel
                case (4)
                    ! Plane
                    if (is_in_plane(pos, selection%plane)) Tdomain%specel(n)%output = sel
                end select

                call c_f_pointer(selection%next, selection)
            end do
        end do

    end subroutine select_output_elements

    subroutine read_input (Tdomain, code)
        use sdomain
        use semdatafiles
        use sem_mpi
        use constants
        use mcapteur
        use surface_input !, only : read_surface_input,

        implicit none

        type(domain), intent(inout)  :: Tdomain
        integer, intent(out)         :: code
        character(Len=MAX_FILE_SIZE) :: fnamef
        logical                      :: logic_scheme
        integer                      :: rg
        integer                      :: i
        integer                      :: outflag

        rg = Tdomain%rank

        call semname_file_input_spec(fnamef)

        call read_sem_config(Tdomain%config, Tdomain%rank, 3, trim(fnamef)//C_NULL_CHAR, code)

        if (code/=1) then
            stop 1
        endif
        if (rg==0) call dump_config(Tdomain%config) !Print of configuration on the screen
        ! On copie les parametres renvoyes dans la structure C
        Tdomain%Title_simulation          = fromcstr(Tdomain%config%run_name)
        Tdomain%TimeD%acceleration_scheme = Tdomain%config%accel_scheme .ne. 0
        Tdomain%TimeD%velocity_scheme     = Tdomain%config%veloc_scheme .ne. 0
        Tdomain%TimeD%duration            = Tdomain%config%sim_time
        Tdomain%TimeD%alpha               = Tdomain%config%alpha
        Tdomain%TimeD%beta                = Tdomain%config%beta
        Tdomain%TimeD%gamma               = Tdomain%config%gamma
        Tdomain%TimeD%fmax                = Tdomain%config%fmax
        Tdomain%TimeD%type_timeinteg      = Tdomain%config%type_timeinteg
        select case(Tdomain%TimeD%type_timeinteg)
        case (0)  ! Newmark
            Tdomain%TimeD%nsubsteps = 1
        case (1) ! RK4
            Tdomain%TimeD%nsubsteps = 5
        case (2) ! LDDRK64
            Tdomain%TimeD%nsubsteps = 2
        case (3) ! Unused
            Tdomain%TimeD%nsubsteps = 1
        case (4) ! Unused
            Tdomain%TimeD%nsubsteps = 1
        end select
        Tdomain%nl_flag = .false.
        Tdomain%use_avg = .false.
        Tdomain%prot_at_time = Tdomain%config%prot_at_time
        if(Tdomain%config%nl_flag == 1) Tdomain%nl_flag = .true.
        if(Tdomain%config%use_avg == 1) Tdomain%use_avg = .true.
        if (rg==0) then
            if (Tdomain%TimeD%alpha /= 0.5 .or. Tdomain%TimeD%beta /= 0.5 .or. Tdomain%TimeD%gamma /= 1.) then
                write(*,*) "***WARNING*** : Les parametres alpha,beta,gamma sont ignores dans cette version"
                write(*,*) "***WARNING*** : on prend: alpha=0.5, beta=0.5, gamma=1"
            endif
        end if

        Tdomain%TimeD%alpha = 0.5
        Tdomain%TimeD%beta = 0.5
        Tdomain%TimeD%gamma = 1.
        ! OUTPUT FIELDS
        Tdomain%nReqOut = 0
        do i = 0, OUT_LAST
            outflag = Tdomain%config%out_variables(i+1)
            Tdomain%out_var_snap(i) = 0
            Tdomain%out_var_capt(i) = 0
            if (outflag==1 .or. outflag==3) then
                ! For snapshots
                Tdomain%out_var_snap(i) = 1
            endif
            if (outflag==1 .or. outflag==2) then
                ! For traces
                Tdomain%out_var_capt(i) = 1
                Tdomain%nReqOut = Tdomain%nReqOut + OUT_VAR_DIMS_3D(i)
            endif
        end do


        Tdomain%TimeD%courant             = Tdomain%config%courant
        Tdomain%mesh_file                 = fromcstr(Tdomain%config%mesh_file)
        call semname_read_input_meshfile(rg,Tdomain%mesh_file,fnamef) !indicates the path to the mesh file for this proc"
        Tdomain%mesh_file             = fnamef
        Tdomain%aniso                 = Tdomain%config%anisotropy .ne. 0
        Tdomain%material_file         = fromcstr(Tdomain%config%mat_file)
        Tdomain%logicD%save_trace     = Tdomain%config%save_traces .ne. 0
        Tdomain%logicD%save_snapshots = Tdomain%config%save_snap .ne. 0
        Tdomain%ngroup                = Tdomain%config%n_group_outputs
        Tdomain%logicD%save_restart   = Tdomain%config%prorep_iter .ne. 0
        Tdomain%logicD%MPML           = .false.
        Tdomain%MPML_coeff            = Tdomain%config%mpml
        if (Tdomain%config%mpml/=0) then
            Tdomain%logicD%MPML = .true.
        end if

        Tdomain%logicD%run_restart = Tdomain%config%prorep .ne. 0
        Tdomain%TimeD%iter_reprise = Tdomain%config%prorep_restart_iter
        Tdomain%TimeD%ncheck       = Tdomain%config%prorep_iter ! frequence de sauvegarde

        Tdomain%TimeD%ntrace         = Tdomain%config%traces_interval ! XXX
        Tdomain%traces_format        = Tdomain%config%traces_format
        Tdomain%TimeD%time_snapshots = Tdomain%config%snap_interval
        logic_scheme                 = Tdomain%TimeD%acceleration_scheme .neqv. Tdomain%TimeD%velocity_scheme
        if(.not. logic_scheme) then
            stop "Both acceleration and velocity schemes: no compatibility, chose only one."
        end if
        ! Amortissement
        Tdomain%n_sls     = Tdomain%config%nsolids
        Tdomain%T1_att    = Tdomain%config%atn_band(1)
        Tdomain%T2_att    = Tdomain%config%atn_band(2)
        Tdomain%T0_modele = Tdomain%config%atn_period
        if (rg==0) then
            write(*,*) "Attenuation SLS =", Tdomain%n_sls
            write(*,*) "         period =", Tdomain%T0_modele
            write(*,*) "           band =", Tdomain%T1_att, Tdomain%T2_att
        end if

        ! boundary conditions? If yes: geometrical properties read in the mesh files.
        Tdomain%logicD%surfBC = Tdomain%config%surface_find /= 0
        !! Add by Mtaro
        call init_surface_source(Tdomain)
        if (Tdomain%logicD%surfBC) then
            call read_surface_input(Tdomain, Tdomain%config)
        endif

        ! Gestion des Sources
        if (Tdomain%config%nextended_sources==0) then
            ! Only a few ponctual sources
            Tdomain%n_source = Tdomain%config%nsources

            !write(*,*) 'TOTAL SOURCE NUMBER IS   ', Tdomain%n_source
            allocate (Tdomain%Ssource(0:Tdomain%n_source-1))
            ! Create point sources from C structures
            call create_sem_sources(Tdomain, Tdomain%config)
        else
            ! Create extended sources from data files
            call create_sem_extended_sources(Tdomain, Tdomain%config)
        endif

        ! NGLL points
        Tdomain%ngll = Tdomain%config%ngll
        ! Mirror
        Tdomain%use_mirror = .false.
        if (Tdomain%config%use_mirror/=0) Tdomain%use_mirror = .true.
        Tdomain%mirror_type = Tdomain%config%mirror_type
        if (Tdomain%config%use_mirror==0) Tdomain%mirror_type = -1
        Tdomain%mirror_expl = .false.
        if (Tdomain%config%mirror_expl/=0) Tdomain%mirror_expl = .true.
        Tdomain%mirror_recalc = .false.
        if (Tdomain%config%mirror_recalc/=0) Tdomain%mirror_recalc = .true.

        ! XXX Unused flags
        Tdomain%logicD%Neumann = .false.

        !---   Reading mesh file
        call read_mesh_file_h5(Tdomain)

        !---
        if (Tdomain%logicD%surfBC) then
            call surface_in_list(Tdomain)
        endif

        !---   Properties of materials.
        call read_material_file(Tdomain)
        call compute_material_boundaries(Tdomain)

        ! Material Earthchunk

        Tdomain%earthchunk_isInit=0
        if( Tdomain%config%material_type == MATERIAL_EARTHCHUNK) then
            Tdomain%earthchunk_isInit=1

        endif


        call select_output_elements(Tdomain, Tdomain%config)
    end subroutine read_input

    subroutine init_surface_source(Tdomain)

        use sdomain

        implicit none
        type(domain), intent(inout) :: Tdomain

        Tdomain%nsurface = 0
        Tdomain%n_NEBC =0
        Tdomain%n_PWBC =0
        Tdomain%n_FTBC =0
        Tdomain%n_DIRIC=0

    end subroutine init_surface_source

    function getLine (fid, comment_Tag) result(nextLine)

        integer,          intent(in) :: fid
        character(len=1), intent(in) :: comment_Tag
        character(len=MAX_FILE_SIZE) :: nextLine
        integer :: lineCount = 200, i, stat

        do i = 1, lineCount
            read(fid, fmt="(A)",IOSTAT = stat) nextLine
            nextLine = adjustL(nextLine)
            if(stat /= 0) then
                nextLine = " "
                exit
            else if(nextLine(1:1) /= comment_Tag) then
                exit
            end if
        end do

    end function getLine

end module semconfig

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
