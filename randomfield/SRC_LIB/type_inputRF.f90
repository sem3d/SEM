module type_inputRF

    use mpi
    use readFile_RF
    use readUNV_RF
    use systemUt_RF

    implicit none

    type :: IPT_RF
        !LOCAL
        integer :: log_ID = -1
        integer :: comm = -1
        integer :: rang = -1
        integer :: nb_procs = -1
        logical :: init=.false.
        logical :: alloc=.false.
        !MESH
        integer :: nDim_mesh
        integer :: meshMod
        double precision, dimension(:), allocatable :: xMaxGlob, xMinGlob;
        double precision, dimension(:), allocatable :: xMaxGlob_in, xMinGlob_in;
        integer         , dimension(:), allocatable :: pointsPerCorrL;
        integer         , dimension(:), allocatable :: procPerDim
        integer         , dimension(:), allocatable :: coords
        !UNV
        double precision, dimension(:,:), allocatable :: coordList_local
        integer         , dimension(:,:), allocatable :: connectList_local
        double precision, dimension(:,:), pointer :: coordList
        integer         , dimension(:,:), pointer :: connectList
        logical :: monotype
        logical :: unv = .false.
        character(len=buf_RF) :: unv_path, outputFolder, outputName

        !GENERATION
        integer :: nDim_gen
        double precision   :: fieldAvg = -1, fieldVar = -1;
        double precision, dimension(:), allocatable :: corrL_in, overlap_in
        integer :: corrMod = -1 !1 for Gaussian
        integer :: margiFirst = -1 !1 for Gaussian, 2 for Lognormal
        integer :: method = -1 !1 for Isotropic, 2 for Shinozuka, 3 for Randomization, 4 for FFT
        integer :: Nmc = -1, seedStart
        integer, dimension(:), allocatable :: nFields
        integer :: localizationLevel = -1, nTotalFields = -1
        logical :: global
        !FILE MANAGER
        logical :: writeDataSet = .true.
        logical :: sameFolder = .false.
        integer :: outputStyle = 1 !1: parallel hdf5, 2: hdf5 per proc
        logical :: write_intermediate_files = .false.
        !integer :: ignoreTillLocLevel = 0 !<1 doesn't affetct the behaviour of the program (for restarts)
        logical :: sampleFields = .true.
        logical :: writeUNVinterpolation

        !PROCESS
        integer :: nDim = -1
        double precision, dimension(:), allocatable :: corrL, overlap, xStep
        double precision, dimension(:), allocatable :: stepProc, procExtent
        integer :: gen_groupMax, gen_group, gen_comm, gen_nbProcs, gen_rang
        integer :: loc_groupMax, loc_group, loc_comm, loc_nbProcs, loc_rang
        logical :: extLoc
        integer, dimension(:), allocatable :: neigh, op_neigh
        integer, dimension(:,:), allocatable :: neighShift

    end type IPT_RF

contains
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_IPT_RF(&
                                IPT, &
                                log_ID, &
                                comm, &
                                rang, &
                                nb_procs, &
                                nDim, &
                                meshMod, &
                                xMinGlob_in, &
                                xMaxGlob_in, &
                                pointsPerCorrL, &
                                procPerDim, &
                                fieldAvg, &
                                fieldVar, &
                                corrL_in, &
                                overlap_in, &
                                corrMod, &
                                margiFirst, &
                                method, &
                                Nmc, &
                                seedStart, &
                                nFields, &
                                localizationLevel, &
                                writeDataSet, &
                                sameFolder, &
                                outputStyle, &
                                write_intermediate_files, &
                                sampleFields, &
                                writeUNVinterpolation, &
                                outputFolder, &
                                outputName, &
                                unv, &
                                unv_path, &
                                monotype, &
                                coordList_local, &
                                connectList_local)

            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT

            !INPUT----------------------------------
            integer, intent(in) :: log_ID
            integer, intent(in) :: comm
            integer, intent(in) :: rang
            integer, intent(in) :: nb_procs
            integer, intent(in) :: nDim
            character(len=*), intent(in) :: outputFolder, outputName

            !MESH
            integer, intent(in) :: meshMod
            double precision, dimension(:), intent(in) :: xMaxGlob_in, xMinGlob_in;
            integer         , dimension(:), intent(in) :: pointsPerCorrL;
            integer         , dimension(:), intent(in) :: procPerDim
            !UNV
            double precision, dimension(:,:), intent(in), target, optional :: coordList_local
            integer         , dimension(:,:), intent(in), target, optional :: connectList_local
            logical, intent(in), optional :: monotype
            logical, intent(in), optional :: unv
            character(len=1024), intent(in), optional :: unv_path

            !GENERATION
            double precision, intent(in)   :: fieldAvg, fieldVar
            double precision, dimension(:), intent(in) :: corrL_in, overlap_in
            integer, intent(in) :: corrMod!1 for Gaussian
            integer, intent(in) :: margiFirst!1 for Gaussian, 2 for Lognormal
            integer, intent(in) :: method!1 for Isotropic, 2 for Shinozuka, 3 for Randomization, 4 for FFT
            integer, intent(in) :: Nmc, seedStart
            integer, dimension(:), intent(in) :: nFields
            integer, intent(in) :: localizationLevel

            !FILE MANAGER
            logical, intent(in) :: writeDataSet
            logical, intent(in) :: sameFolder
            integer, intent(in) :: outputStyle!1: parallel hdf5, 2: hdf5 per proc
            logical, intent(in) :: write_intermediate_files
            logical, intent(in) :: sampleFields
            logical, intent(in) :: writeUNVinterpolation

            call create_folder(outputFolder, ".", rang, comm)
            !create xmf and h5 folders
            call create_folder("xmf", outputFolder, rang, comm)
            call create_folder("h5", outputFolder, rang, comm)

            call allocate_IPT_RF(IPT, nDim, log_ID, rang, comm, nb_procs)

            IPT%nDim_mesh = nDim
            IPT%nDim_gen  = nDim
            IPT%meshMod   = meshMod
            IPT%xMaxGlob_in = xMaxGlob_in
            IPT%xMinGlob_in = xMinGlob_in
            IPT%pointsPerCorrL = pointsPerCorrL
            if(present(coordList_local)) IPT%coordList => coordList_local
            if(present(connectList_local)) IPT%connectList => connectList_local
            if(present(monotype)) IPT%monotype = monotype
            if(present(unv)) IPT%unv = unv
            if(present(unv_path)) IPT%unv_path = unv_path
            IPT%fieldAvg = fieldAvg
            IPT%fieldVar = fieldVar
            IPT%corrL_in    = corrL_in
            IPT%overlap_in  = overlap_in
            IPT%corrMod  = corrMod
            IPT%method   = method
            IPT%Nmc      = Nmc
            IPT%nFields  = nFields
            IPT%margiFirst    = margiFirst
            IPT%seedStart     = seedStart
            IPT%procPerDim = procPerDim
            IPT%localizationLevel = localizationLevel
            IPT%writeDataSet = writeDataSet
            IPT%sameFolder = sameFolder
            IPT%outputStyle = outputStyle
            IPT%write_intermediate_files = write_intermediate_files
            IPT%sampleFields = sampleFields
            IPT%writeUNVinterpolation = writeUNVinterpolation
            IPT%outputFolder = outputFolder
            IPT%outputName = outputName

            IPT%xMinGlob = IPT%xMinGlob_in
            IPT%xMaxGlob = IPT%xMaxGlob_in


        end subroutine init_IPT_RF


        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_IPT_RF_std(&
                                    IPT, &
                                    comm, &
                                    nDim, &
                                    xMinGlob_in, &
                                    xMaxGlob_in, &
                                    fieldAvg, &
                                    fieldVar, &
                                    corrL_in, &
                                    corrMod, &
                                    margiFirst, &
                                    seedStart, &
                                    outputFolder, &
                                    outputName)

            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT

            !INPUT----------------------------------
            integer, intent(in) :: comm
            integer, intent(in) :: nDim
            character(len=*), intent(in) :: outputFolder, outputName

            !MESH
            double precision, dimension(:), intent(in) :: xMaxGlob_in, xMinGlob_in;

            !GENERATION
            double precision, intent(in)   :: fieldAvg, fieldVar
            double precision, dimension(:), intent(in) :: corrL_in
            integer, intent(in) :: corrMod!1 for Gaussian
            integer, intent(in) :: margiFirst!1 for Gaussian, 2 for Lognormal
            integer, intent(in) :: seedStart

            !LOCAL
            integer :: rang
            integer :: nb_procs
            integer :: code
            integer, dimension(nDim) :: pointsPerCorrL, procPerDim, nFields
            double precision, dimension(nDim) :: overlap_in

            pointsPerCorrL = 5
            procPerDim = 1
            overlap_in = 0.0D0
            nFields = 1
#ifdef MAKELOG
            stop("init_IPT_RF_std cannot be used if MAKELOG is activated")
#endif
            call MPI_COMM_RANK(comm, rang, code)
            call MPI_COMM_SIZE(comm, nb_procs, code)

            call init_IPT_RF(&
                            IPT=IPT, &
                            log_ID=0, &
                            comm=comm, &
                            rang=rang, &
                            nb_procs=nb_procs, &
                            nDim=nDim, &
                            meshMod=1, &
                            xMaxGlob_in=xMaxGlob_in, &
                            xMinGlob_in=xMinGlob_in, &
                            pointsPerCorrL=pointsPerCorrL, &
                            procPerDim=procPerDim, &
                            fieldAvg=fieldAvg, &
                            fieldVar=fieldVar, &
                            corrL_in=corrL_in, &
                            overlap_in=overlap_in, &
                            corrMod=corrMod, &
                            margiFirst=margiFirst, &
                            method=4, &
                            Nmc=1, &
                            seedStart=seedStart, &
                            nFields=nFields, &
                            localizationLevel=1, &
                            writeDataSet=.true., &
                            sameFolder=.true., &
                            outputStyle=1, &
                            write_intermediate_files=.false., &
                            sampleFields=.true., &
                            writeUNVinterpolation=.false., &
                            outputFolder=outputFolder, &
                            outputName=outputName)

        end subroutine init_IPT_RF_std

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine allocate_IPT_RF(IPT, nDim, log_file_RF_ID, rang, comm, nb_procs)

            implicit none

            !INPUT
            integer, intent(in) :: nDim
            integer, intent(in) :: log_file_RF_ID, rang, comm, nb_procs
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT

            if(.not. allocated(IPT%xMaxGlob)) allocate(IPT%xMaxGlob(nDim))
            if(.not. allocated(IPT%xMinGlob)) allocate(IPT%xMinGlob(nDim))
            if(.not. allocated(IPT%xMaxGlob_in)) allocate(IPT%xMaxGlob_in(nDim))
            if(.not. allocated(IPT%xMinGlob_in)) allocate(IPT%xMinGlob_in(nDim))
            if(.not. allocated(IPT%pointsPerCorrL)) allocate(IPT%pointsPerCorrL(nDim))
            if(.not. allocated(IPT%corrL_in)) allocate(IPT%corrL_in(nDim))
            if(.not. allocated(IPT%overlap_in)) allocate(IPT%overlap_in(nDim))
            if(.not. allocated(IPT%procPerDim)) allocate(IPT%procPerDim(nDim))
            if(.not. allocated(IPT%nFields)) allocate(IPT%nFields(nDim))
            !Process
            if(.not. allocated(IPT%corrL)) allocate(IPT%corrL(nDim))
            if(.not. allocated(IPT%overlap)) allocate(IPT%overlap(nDim))
            if(.not. allocated(IPT%xStep)) allocate(IPT%xStep(nDim))
            if(.not. allocated(IPT%stepProc)) allocate(IPT%stepProc(nDim))
            if(.not. allocated(IPT%procExtent)) allocate(IPT%procExtent(nDim))
            if(.not. allocated(IPT%coords)) allocate(IPT%coords(nDim))
            if(.not. allocated(IPT%neigh)) allocate(IPT%neigh((3**nDim)-1))
            if(.not. allocated(IPT%op_neigh)) allocate(IPT%op_neigh((3**nDim)-1))
            if(.not. allocated(IPT%neighShift)) allocate(IPT%neighShift(nDim,(3**nDim)-1))

            IPT%log_ID = log_file_RF_ID
            IPT%rang   = rang
            IPT%comm   = comm
            IPT%nb_procs = nb_procs
            IPT%alloc   = .true.
            IPT%xMaxGlob = -1.0D0
            IPT%xMinGlob = -1.0D0
            IPT%xMaxGlob_in = -1.0D0
            IPT%xMinGlob_in = -1.0D0
            IPT%pointsPerCorrL = -1
            IPT%corrL   = -1.0D0
            IPT%overlap = -1.0D0

        end subroutine allocate_IPT_RF


        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine finalize_IPT_RF(IPT)

            implicit none

            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT

            if(allocated(IPT%xMaxGlob)) deallocate(IPT%xMaxGlob)
            if(allocated(IPT%xMinGlob)) deallocate(IPT%xMinGlob)
            if(allocated(IPT%pointsPerCorrL)) deallocate(IPT%pointsPerCorrL)
            if(allocated(IPT%corrL_in)) deallocate(IPT%corrL_in)
            if(allocated(IPT%overlap_in)) deallocate(IPT%overlap_in)
            if(allocated(IPT%procPerDim)) deallocate(IPT%procPerDim)
            if(allocated(IPT%nFields)) deallocate(IPT%nFields)
            if(allocated(IPT%xMaxGlob_in)) deallocate(IPT%xMaxGlob_in)
            if(allocated(IPT%xMinGlob_in)) deallocate(IPT%xMinGlob_in)

            if(allocated(IPT%coordList_local)) deallocate(IPT%coordList_local)
            if(allocated(IPT%connectList_local)) deallocate(IPT%connectList_local)
            if(associated(IPT%coordList)) nullify(IPT%coordList)
            if(associated(IPT%connectList)) nullify(IPT%connectList)

            if(allocated(IPT%corrL)) deallocate(IPT%corrL)
            if(allocated(IPT%overlap)) deallocate(IPT%overlap)
            if(allocated(IPT%xStep)) deallocate(IPT%xStep)
            if(allocated(IPT%stepProc)) deallocate(IPT%stepProc)
            if(allocated(IPT%procExtent)) deallocate(IPT%procExtent)
            if(allocated(IPT%coords)) deallocate(IPT%coords)

            if(allocated(IPT%neigh)) deallocate(IPT%neigh)
            if(allocated(IPT%op_neigh)) deallocate(IPT%op_neigh)
            if(allocated(IPT%neighShift)) deallocate(IPT%neighShift)

        end subroutine finalize_IPT_RF

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine copy_IPT_RF(IPT, IPT_orig)

            implicit none

            !INPUT
            type(IPT_RF), intent(in)  :: IPT_orig
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT

            IPT%init   = IPT_orig%init
            IPT%xMaxGlob = IPT_orig%xMaxGlob
            IPT%xMinGlob = IPT_orig%xMinGlob
            IPT%xMinGlob_in = IPT_orig%xMinGlob_in
            IPT%xMinGlob_in = IPT_orig%xMinGlob_in
            IPT%pointsPerCorrL = IPT_orig%pointsPerCorrL
            IPT%corrL_in   = IPT_orig%corrL_in
            IPT%overlap_in = IPT_orig%overlap_in
            IPT%procPerDim = IPT_orig%procPerDim
            IPT%nFields = IPT_orig%nFields
            IPT%log_ID = IPT_orig%log_ID
            IPT%rang   = IPT_orig%rang
            IPT%overlap = IPT_orig%overlap
            IPT%comm = IPT_orig%comm
            IPT%nb_procs = IPT_orig%nb_procs
            IPT%writeDataSet = IPT_orig%writeDataSet

            IPT%nDim_mesh = IPT_orig%nDim_mesh
            IPT%meshMod = IPT_orig%meshMod

            IPT%coordList   => IPT_orig%coordList
            IPT%connectList => IPT_orig%connectList

            !UNV
            !if(allocated(IPT_orig%coordList)) then
            !    if(.not. allocated(IPT%coordList)) then
            !        allocate(IPT%coordList(size(IPT_orig%connectList,1), size(IPT_orig%connectList,2)))
            !    end if
            !    IPT%coordList = IPT_orig%coordList
            !end if
            !if(allocated(IPT_orig%connectList)) then
            !   if(.not. allocated(IPT%connectList)) then
            !        allocate(IPT%connectList(size(IPT_orig%connectList,1), size(IPT_orig%connectList,2)))
            !    end if
            !    IPT%connectList = IPT_orig%connectList
            !end if

            IPT%monotype = IPT_orig%monotype
            IPT%unv = IPT_orig%unv
            IPT%unv_path = IPT_orig%unv_path

            !GENERATION
            IPT%nDim_gen = IPT_orig%nDim_gen
            IPT%fieldAvg = IPT_orig%fieldAvg
            IPT%fieldVar = IPT_orig%fieldVar
            IPT%margiFirst = IPT_orig%margiFirst
            IPT%corrL = IPT_orig%corrL
            IPT%overlap = IPT_orig%overlap
            IPT%corrMod = IPT_orig%corrMod
            IPT%method = IPT_orig%method
            IPT%Nmc = IPT_orig%Nmc
            IPT%seedStart = IPT_orig%seedStart
            IPT%nFields = IPT_orig%nFields
            IPT%localizationLevel = IPT_orig%localizationLevel

            IPT%corrL = IPT_orig%corrL
            IPT%overlap = IPT_orig%overlap
            IPT%xStep = IPT_orig%xStep
            IPT%stepProc = IPT_orig%stepProc
            IPT%procExtent = IPT_orig%procExtent

        end subroutine copy_IPT_RF

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_mesh_input(path, IPT)

            implicit none
            !INPUT
            character(len=*), intent(in) :: path
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT
            !LOCAL
            character(len=1024) , dimension(:,:), allocatable :: dataTable;

            call set_DataTable(path, dataTable)

            call wLog("    set_DataTable")

            call read_DataTable(dataTable, "nDim", IPT%nDim_mesh)
            call read_DataTable(dataTable, "meshMod", IPT%meshMod)

            if(.not. IPT%alloc) call allocate_IPT_RF(IPT, IPT%nDim_mesh, IPT%log_ID, IPT%rang, IPT%comm, IPT%nb_procs)

            select case (IPT%meshMod)
                case(1)
                    if(IPT%rang==0) write(*,*) "   Mesh automatic"
                    call wLog("    Mesh automatic")
                    call read_DataTable(dataTable, "Min", IPT%xMinGlob_in)
                    call read_DataTable(dataTable, "Max", IPT%xMaxGlob_in)
                    call read_DataTable(dataTable, "pointsPerCorrL", IPT%pointsPerCorrL)
                    IPT%xMinGlob = IPT%xMinGlob_in
                    IPT%xMaxGlob = IPT%xMaxGlob_in
                case(2)
                    !stop("Inside read_mesh_input UNV not updated")
                    if(IPT%rang==0) write(*,*) "   Mesh UNV"
                    call wLog("    Mesh UNV")
                    call wLog("        file: ")
                    !call wLog(trim(IPT%unv_path))
                    call read_DataTable(dataTable, "unv_path", IPT%unv_path)
                    call read_DataTable(dataTable, "pointsPerCorrL", IPT%pointsPerCorrL)
                    IPT%unv = .true.
                    call readUNV(IPT%unv_path, IPT%nDim_mesh, IPT%coordList_local, IPT%connectList_local, IPT%monotype, &
                                 IPT%rang, IPT%nb_procs, IPT%coordList, IPT%connectList)
                    !call readUNV_many([IPT%unv_path, IPT%unv_path], IPT%nDim_mesh, IPT%coordList, &
                    !                  IPT%connectList, IPT%monotype, &
                    !                  IPT%rang, IPT%nb_procs, IPT%comm)
                    call wLog("-> defining_UNV_extremes")
                    if(IPT%rang==0) write(*,*) "-> defining_UNV_extremes"
                    call get_Global_Extremes_Mesh(IPT%coordList, IPT%comm, IPT%xMinGlob_in, IPT%xMaxGlob_in)
                    if(IPT%rang==0) write(*,*) " IPT%xMinGlob_in = ", IPT%xMinGlob_in
                    if(IPT%rang==0) write(*,*) " IPT%xMaxGlob_in = ", IPT%xMaxGlob_in
!                    call DispCarvalhol(transpose(IPT%connectList), "transpose(IPT%connectList)", &
!                                       nColumns=8, unit_in=IPT%log_ID)
!                    call DispCarvalhol(transpose(IPT%coordList), "transpose(IPT%coordList)", &
!                                       unit_in=IPT%log_ID)
                case default
                    write(*,*) "meshMod not accepted: ", IPT%meshMod
                    stop(" ")
            end select


            if(allocated(dataTable)) deallocate(dataTable)

        end subroutine read_mesh_input


        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_generation_input(path, IPT)

            implicit none
            !INPUT
            character(len=*), intent(in) :: path
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT
            !LOCAL
            character(len=1024) , dimension(:,:), allocatable :: dataTable;


            call set_DataTable(path, dataTable)
            call read_DataTable(dataTable, "nDim", IPT%nDim_gen)

            if(.not. IPT%alloc) call allocate_IPT_RF(IPT, IPT%nDim_gen, IPT%log_ID, IPT%rang, IPT%comm, IPT%nb_procs)

            call read_DataTable(dataTable, "Nmc", IPT%Nmc)
            call read_DataTable(dataTable, "corrMod"    , IPT%corrMod)
            call read_DataTable(dataTable, "margiFirst" , IPT%margiFirst)
            call read_DataTable(dataTable, "fieldAvg"   , IPT%fieldAvg)
            call read_DataTable(dataTable, "fieldVar"   , IPT%fieldVar)
            call read_DataTable(dataTable, "method"     , IPT%method)
            call read_DataTable(dataTable, "seedStart"  , IPT%seedStart)
            call read_DataTable(dataTable, "corrL"      , IPT%corrL_in)
            call read_DataTable(dataTable, "nFields"    , IPT%nFields)
            call read_DataTable(dataTable, "localizationLevel", IPT%localizationLevel)
            call read_DataTable(dataTable, "overlap", IPT%overlap_in)



            deallocate(dataTable)

        end subroutine read_generation_input

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine validate_input(IPT)

            !OUTPUT
            type(IPT_RF), intent(in)  :: IPT
            !LOCAL
            integer :: i

            !Input Validation
            if(IPT%Nmc < 1) then
               write(*,*) ""
               write(*,*) "ERROR - Number of events should be a positive integer"
               write(*,*) "IPT%Nmc = ", IPT%Nmc
               stop(" ")
            end if

            if((IPT%nDim_gen < 1) .or. (IPT%nDim_mesh < 1) .or. (IPT%nDim_gen /= IPT%nDim_mesh)) then
               write(*,*) ""
               write(*,*) "ERROR - nDim should be a positive integer and should be the same in both mesh and generation files"
               write(*,*) "IPT%nDim_gen  = ", IPT%nDim_gen
               write(*,*) "IPT%nDim_mesh = ", IPT%nDim_mesh
               stop(" ")
            end if

            do i = 1, IPT%nDim_gen
               if(IPT%corrL_in(i) <= 0.0d0) then
                   write(*,*) ""
                   write(*,*) "ERROR - corrL_in should be a positive number greater than 0.0"
                   stop(" ")
               end if
            end do

            if(any(IPT%nFields < 1)) then
               write(*,*) ""
               write(*,*) "ERROR - IPT%nFields should be a positive non-zero value vector"
               stop(" ")
            end if

            if(IPT%localizationLevel == 0 .and. any(IPT%nFields /=1)) then
               write(*,*) ""
               write(*,*) "ERROR - If locLevel = 0 IPT%nFields should have only ones"
               stop(" ")
            end if

        end subroutine validate_input

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_IPT_RF(IPT, name, forLog_in, unit_in)
            !INPUT
            type(IPT_RF), intent(in) :: IPT
            character(len=*), intent(in), optional :: name
            integer, intent(in), optional :: unit_in
            logical, intent(in), optional :: forLog_in
            !LOCAL
            !character(len = 20) :: dblFmt, matDblFmt, intFmt, matIntFmt
            integer :: unit
            logical :: forLog
            logical :: active

            active = .true.
            unit = 6 !Screen
            if(present(unit_in)) unit = unit_in

            forLog = .false.
            if(present(forLog_in)) then
#ifdef MAKELOG
                if(forLog_in) unit = IPT%log_ID
                forLog = forLog_in
#else
                if(forLog_in) write(*,*) "WARNING!!! Inside show_IPT_RF, forLog_in = .true. but MAKELOG was not defined"
                active = .false.
#endif
            end if

            if(active) then

                write(unit,*) "INPUT------------------------------------------------------------"
                if(present(name)) write(unit,*) "  ", name
                write(unit,*) " "
                write(unit,*) " log_ID = ", IPT%log_ID
                write(unit,*) " rang = ", IPT%rang
                write(unit,*) " init = ", IPT%init
                write(unit,*) " comm = ", IPT%comm
                write(unit,*) " "
                write(unit,*) " MESH -----------------"
                write(unit,*) " nDim_mesh = ", IPT%nDim_mesh
                write(unit,*) " meshMod = ", IPT%meshMod
                write(unit,*) " xMaxGlob_in = ", IPT%xMaxGlob_in
                write(unit,*) " xMinGlob_in = ", IPT%xMinGlob_in
                write(unit,*) " xMaxGlob = ", IPT%xMaxGlob
                write(unit,*) " xMinGlob = ", IPT%xMinGlob
                write(unit,*) " pointsPerCorrL = ", IPT%pointsPerCorrL
                write(unit,*) " unv = ", IPT%unv
                if(IPT%unv) then
                    write(unit,*) " unv_path = "//trim(adjustL(IPT%unv_path))
                    write(unit,*) " shape(coordList) = "  , shape(IPT%coordList)
                    write(unit,*) " shape(connectList) = ", shape(IPT%connectList)
                    write(unit,*) " monotype = ", IPT%monotype
                end if
                write(unit,*) " "
                write(unit,*) " GENERATION -----------------"
                write(unit,*) " nDim_gen = ", IPT%nDim_gen
                write(unit,*) " fieldAvg = ", IPT%fieldAvg
                write(unit,*) " fieldVar = ", IPT%fieldVar
                write(unit,*) " overlap_in = ", IPT%overlap_in
                write(unit,*) " corrL_in = ", IPT%corrL_in
                write(unit,*) " overlap = ", IPT%overlap
                write(unit,*) " corrL = ", IPT%corrL
                write(unit,*) " corrMod = ", IPT%corrMod
                write(unit,*) " margiFirst = ", IPT%margiFirst
                write(unit,*) " method = ", IPT%method
                write(unit,*) " Nmc = ", IPT%Nmc
                write(unit,*) " seedStart = ", IPT%seedStart
                write(unit,*) " nFields = ", IPT%nFields
                write(unit,*) " "

                write(unit,*) " "
                write(unit,*) " PROCESS -----------------"
                write(unit,*) " xStep      = ", IPT%xStep
                write(unit,*) " stepProc   = ", IPT%stepProc
                write(unit,*) " procExtent = ", IPT%procExtent
                write(unit,*) " procExtent = ", IPT%procExtent
                write(unit,*) " outputFolder = ", trim(adjustL(IPT%outputFolder))
                write(unit,*) " outputName   = ", trim(adjustL(IPT%outputName))

            end if

        end subroutine show_IPT_RF

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine get_Global_Extremes_Mesh(coordList, comm, xMinGlob, xMaxGlob)
        implicit none

        !INPUT
        double precision, dimension(:,:), intent(in) :: coordList
        integer, intent(in) :: comm
        !OUTPUT
        double precision, dimension(:), intent(out) :: xMinGlob, xMaxGlob
        !LOCAL
        integer :: nDim
        double precision, dimension(:), allocatable :: xMinLoc, xMaxLoc
        integer :: code, i

        !call DispCarvalhol(transpose(coordList), "transpose(coordList)")

        nDim = size(coordList,1)
        allocate(xMinLoc(nDim))
        allocate(xMaxLoc(nDim))

        xMinLoc = minval(coordList, 2)
        xMaxLoc = maxval(coordList, 2)

        call wLog("xMinLoc = ")
        call wLog(xMinLoc)
        call wLog("xMaxLoc = ")
        call wLog(xMaxLoc)

        do i = 1, nDim
            call MPI_ALLREDUCE (xMinLoc(i), xMinGlob(i), 1, MPI_DOUBLE_PRECISION, MPI_MIN, comm,code)
            call MPI_ALLREDUCE (xMaxLoc(i), xMaxGlob(i), 1, MPI_DOUBLE_PRECISION, MPI_MAX, comm,code)
        end do

        call wLog(" ")
        call wLog("xMinGlob = ")
        call wLog(xMinGlob)
        call wLog("xMaxGlob = ")
        call wLog(xMaxGlob)

        if(allocated(xMinLoc)) deallocate(xMinLoc)
        if(allocated(xMaxLoc)) deallocate(xMaxLoc)

    end subroutine get_Global_Extremes_Mesh

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_IPTneigh(IPT, name, onlyExisting, forLog, unit_in)
            implicit none
            !INPUT
            type(IPT_RF), intent(in) :: IPT
            character(len=*), intent(in), optional :: name
            logical, intent(in) :: onlyExisting
            logical, intent(in), optional :: forLog
            integer, intent(in), optional :: unit_in
            !LOCAL
            character(len = 3) :: nDim, space
            character(len = 50) :: fmtNum, fmtChar
            character(len = 1) :: nbChar
            integer :: i
            integer :: unit
            logical :: active
            integer :: numb_space = 9

            active = .true.
            unit = 6 !Screen
            if(present(unit_in)) unit = unit_in

            if(present(forLog)) then
#ifdef MAKELOG
                if(forLog) unit = IPT%log_ID
#else
                !if(forLog) write(*,*) "WARNING!!! Inside show_MESHneigh, forLog = .true. but MAKELOG was not defined"
                if(forLog) active = .false.
#endif
            end if

            if(active) then

                if(unit <= 0) then
                    write(*,*) "ERROR!!! Inside show_MESHneigh unit = ", unit
                    stop("")
                end if

                nDim   = trim(numb2String(IPT%nDim))
                space  = trim(numb2String(IPT%nDim*numb_space + 1))
                nbChar = trim(numb2String(numb_space))

                if(present(name)) write(unit,*) "|  ", name

                fmtChar = "(A7, A"//space//")"
                fmtNum = "(I6, A1, "//nDim//"I"//nbChar//",A1)"

                write(unit,fmtChar) "Neigh|","Shift|"

                do i = 1, size(IPT%neigh)
                    if(onlyExisting .and. IPT%neigh(i)<0) cycle
                    write(unit,fmtNum) IPT%neigh(i), &
                                       "|", IPT%neighShift(:,i), "|"
                end do


                write(unit,*) ""

            end if

        end subroutine show_IPTneigh

end module type_inputRF
