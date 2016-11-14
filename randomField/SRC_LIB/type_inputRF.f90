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
        integer :: timeFolder
        integer :: nSamples
        integer :: application !1 for library, 2 for SEM
        character(len=buf_RF), dimension(:), allocatable :: out_folders, out_names
        character(len=buf_RF), dimension(:), allocatable :: mesh_inputs, gen_inputs
        character(len=buf_RF) :: appFolder
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
        !STATISTICS
        integer :: calculateCorrL
        integer :: deleteSampleAfterStatistics
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

    type :: property_RF
        integer :: mat
        character(len=buf_RF) :: name
        double precision :: avg
        integer :: corrMod
        double precision, dimension(3) :: corrL
        integer :: margiF
        double precision :: CV
        integer :: seedStart
        double precision, dimension(3) :: bbox_min, bbox_max

    end type property_RF

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
                                connectList_local, &
                                calculateCorrL, &
                                deleteSampleAfterStatistics)

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

            !STATISTICS
            integer :: calculateCorrL, deleteSampleAfterStatistics

            !FILE MANAGER
            logical, intent(in) :: writeDataSet
            logical, intent(in) :: sameFolder
            integer, intent(in) :: outputStyle!1: parallel hdf5, 2: hdf5 per proc
            logical, intent(in) :: write_intermediate_files
            logical, intent(in) :: sampleFields
            logical, intent(in) :: writeUNVinterpolation

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
            IPT%calculateCorrL = calculateCorrL
            IPT%deleteSampleAfterStatistics = deleteSampleAfterStatistics

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
                            outputName=outputName, &
                            calculateCorrL=0, &
                            deleteSampleAfterStatistics = 0)

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

            if(allocated(IPT%out_folders)) deallocate(IPT%out_folders)
            if(allocated(IPT%out_names)) deallocate(IPT%out_names)
            if(allocated(IPT%mesh_inputs)) deallocate(IPT%mesh_inputs)
            if(allocated(IPT%gen_inputs)) deallocate(IPT%gen_inputs)

        end subroutine finalize_IPT_RF

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine read_main_input(IPT, auto)

            implicit none
            !INPUT
            logical, intent(in), optional :: auto
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT
            !LOCAL
            character(len=buf_RF) , dimension(:,:), allocatable :: dataTable;
            character(len=buf_RF) :: SEM_gen_path, path
            integer :: i
            integer :: code
            logical :: fExist
            logical :: secondStep

            path = "./random.spec"
            fExist = fileExist(path)
            if(IPT%rang == 0) write(*,*) "      ", trim(adjustL(path)), " exist? ", fExist

            secondStep = .false.
            if(present(auto)) then
                if(auto) fExist = .false.
                secondStep = .true.
            end if

            if(fExist) then
                path = "./domains.txt"
                if(.not. fileExist("./domains.txt")) stop "ERROR: ./random.spec file found but domains.txt not (you need to make the mesh before the generating the properties"
                IPT%application = SEM
                SEM_gen_path = "./mat"
            else
                path = "./RF_main_input"
                fExist = fileExist(path)
                if(fExist) then
                    IPT%application = NATIVE
                else
                    stop "ERROR: nor ./RF_main_input  or ./random.spec files found"
                end if
                if(fileExist("./domains.txt") .and. (.not. secondStep)) stop "To use the randomfield library with SEM the inputs should be given by random.spec (there's no more RF_main_input file)" 
                !write(*,*) "Before Datatable"
                if(IPT%rang == 0) write(*,*) "      path = ", trim(adjustL(path))
                call set_DataTable(path, dataTable)

           end if

            if(IPT%application == NATIVE) then
                if(IPT%rang == 0) write(*,*) "      Native lecture"
                call read_DataTable(dataTable, "nSamples", IPT%nSamples)
                call read_DataTable(dataTable, "timeFolder", IPT%timeFolder)
                call read_DataTable(dataTable, "calculateCorrL", IPT%calculateCorrL)
                call read_DataTable(dataTable, "deleteSampleAfterStatistics", IPT%deleteSampleAfterStatistics)
                allocate(IPT%out_folders(IPT%nSamples))
                allocate(IPT%out_names(IPT%nSamples))
                allocate(IPT%mesh_inputs(IPT%nSamples))
                allocate(IPT%gen_inputs(IPT%nSamples))
                do i =1, IPT%nSamples
                    call read_DataTable(dataTable, stringNumb_join("mesh_input_",i), IPT%mesh_inputs(i))
                    call read_DataTable(dataTable, stringNumb_join("gen_input_",i), IPT%gen_inputs(i))
                    call read_DataTable(dataTable, stringNumb_join("out_folder_",i), IPT%out_folders(i))
                    call read_DataTable(dataTable, stringNumb_join("out_name_",i), IPT%out_names(i))
                end do

                IPT%outputName = IPT%out_names(1)

            else if (IPT%application == SEM) then
                if(IPT%rang == 0) write(*,*) "     SEM files lecture"
                IPT%appFolder = SEM_gen_path
                if(IPT%rang == 0) call generateMain_inputSEM(SEM_gen_path)
                if(IPT%rang == 0) write(*,*) "     END of SEM main_input generation"
                call MPI_BARRIER(IPT%comm, code)
            end if

            if(allocated(dataTable)) deallocate(dataTable)

        end subroutine read_main_input

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine check_constantField(IPT)

            implicit none
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT
            !LOCAL
            double precision :: tolerance

            tolerance = 1e-6

            if(IPT%fieldVar < tolerance) then
                write(*,*) "CONSTANT FIELD"
                IPT%pointsPerCorrL = 2
                IPT%corrL_in = (IPT%xMaxGlob_in - IPT%xMinGlob_in)
                IPT%corrMod = cm_GAUSSIAN
                IPT%margiFirst = fom_GAUSSIAN

            end if

        end subroutine check_constantField


        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine generateMain_inputSEM(SEM_gen_path)

            implicit none
            !INPUT
            character(len=*), intent(in) :: SEM_gen_path
            !LOCAL
            integer :: i, j
            integer :: fid, fid_2
            integer :: nMat
            double precision, dimension(:,:), allocatable :: bbox_min, bbox_max
            character(len=buf_RF) :: mesh_path, gen_path, absPath
            character(len=buf_RF) :: out_folder
            character(len=buf_RF) :: buffer
            integer, dimension(:), allocatable :: nb_Mat
            integer, dimension(:), allocatable :: nProp_Mat, label_Mat
            integer :: mat, nSamples, prop_count, mat_Nb, assocMat
            type(property_RF), dimension(:), allocatable :: prop
            integer :: VP_VS_RHO_flag
            double precision :: VP, VS, RHO, VAR
            double precision, dimension(0:2) :: bb_min, bb_max
            logical :: found

            fid_2 = 19

            !READING random.spec
            open (unit = fid_2 , file = "./random.spec", action = 'read', status="old", form="formatted")

                !!1)  Counting number of samples to be generated
                mat = -1
                buffer = getLine(fid_2, '#') !number of materials
                read(buffer,*) nMat
                allocate(nProp_Mat(0:nMat-1))
                allocate(nb_Mat(0:nMat-1))
                print*, "nMat = ", nMat
                do mat = 0, nMat - 1
                    buffer = getLine(fid_2, '#') !Material Number
                    read(buffer,*) nb_Mat(mat)
                    buffer = getLine(fid_2, '#') !VP_VS_RHO and value
                    buffer = getLine(fid_2, '#') !number of properties
                    read(buffer,*) nProp_Mat(mat)
                    do j = 1, nProp_Mat(mat)
                        buffer = getLine(fid_2, '#')
                    end do
                end do

                nSamples = sum(nProp_Mat(:))

                allocate(prop(0:nSamples-1))
                
                print*, "size(prop) = ", size(prop)

                !!2) Reading Samples Properties
                rewind(fid_2)
                prop_count = 0

                buffer = getLine(fid_2, '#') !number of materials
                do mat = 0, nMat - 1
                    buffer = getLine(fid_2, '#') !Material Number
                    buffer = getLine(fid_2, '#') !VP_VS_RHO flag
                    read(buffer,*) VP_VS_RHO_flag
                    if(VP_VS_RHO_flag == 1) then
                        read(buffer,*) VP_VS_RHO_flag, VP, VS, RHO
                        print*, "VP = ", VP, "VS = ", VS, "RHO = ", RHO
                    end if
                    buffer = getLine(fid_2, '#') !Number of Properties
                    do j = 1, nProp_Mat(mat)
                        buffer = getLine(fid_2, '#')
                        prop(prop_count)%mat = nb_Mat(mat)

                        read(buffer,*) prop(prop_count)%name, &
                                       prop(prop_count)%avg, &
                                       prop(prop_count)%corrMod,   &
                                       prop(prop_count)%corrL(1),  &
                                       prop(prop_count)%corrL(2),  &
                                       prop(prop_count)%corrL(3),  &
                                       prop(prop_count)%margiF,    &
                                       prop(prop_count)%CV,        &
                                       prop(prop_count)%seedStart

                        if(VP_VS_RHO_flag == 1) then
                            if(trim(adjustl(prop(prop_count)%name)) == "Lambda") then
                                prop(prop_count)%avg = 2d0*RHO*VS**2d0
                            else if(trim(adjustl(prop(prop_count)%name)) == "Kappa") then
                                prop(prop_count)%avg = RHO*(VP**2d0 - 4d0*(VS**2d0)/3d0)
                            else if(trim(adjustl(prop(prop_count)%name)) == "Mu") then
                                prop(prop_count)%avg = RHO*VS**2d0
                            else if(trim(adjustl(prop(prop_count)%name)) == "Density") then
                                prop(prop_count)%avg = RHO 
                            end if
                        end if

                        prop_count = prop_count + 1
                    end do
                end do

            close(fid_2)

            !READING domains.txti
            allocate(bbox_min(0:2,0:nMat-1))
            allocate(bbox_max(0:2,0:nMat-1))

            bbox_min(:,:) = MAX_DOUBLE
            bbox_max(:,:) = MIN_DOUBLE

            open (unit = fid_2 , file = "./domains.txt", action = 'read', status="old", form="formatted")
                
                buffer = getLine(fid_2, '#')
                write(*,*) "buffer = ", trim(adjustL(buffer))
                
                do while (trim(adjustL(buffer)) /= "eof_gl")
                   
                    
                    read(buffer,*) mat_Nb, bb_min(0), bb_min(1), bb_min(2), &
                                          bb_max(0), bb_max(1), bb_max(2), assocMat
                    buffer = getLine(fid_2, '#')
                    !write(*,*) "buffer = ", trim(adjustL(buffer))
                    
                    found = .false.
                    do i = 0, size(nb_Mat)-1
                        if (nb_Mat(i) .eq. assocMat) then
                            found = .true.
                            exit
                        end if
                    end do

                    if(.not. found) cycle
                    
                    where(bbox_min(:,i) > bb_min(:)) bbox_min(:,i) = bb_min
                    where(bbox_max(:,i) < bb_max(:)) bbox_max(:,i) = bb_max
                    !buffer = "eof_gl" !TEST 
                end do


            close(fid_2)

            prop_Count = 0
            do mat = 0, nMat-1
                do j = 0, nProp_Mat(mat)-1
                    mat_Nb = prop(prop_count)%mat
                    prop(prop_Count)%bbox_min = bbox_min(:, mat)
                    prop(prop_Count)%bbox_max = bbox_max(:, mat)
                    prop_Count = prop_count + 1
                end do
            end do


            !CREATING FOLDER
            call create_folder(".", SEM_gen_path, 0, 0, singleProc = .true.)
            call create_folder(".", string_join_many(SEM_gen_path,"/input"), 0, 0, singleProc = .true.)

            !WRITING new RF_main_input

            fid = 18

            open (unit = fid , file = string_join_many("RF_main_input"), action = 'write')

            write(fid,"(A)")  "$application 1"
            write(fid,"(A)")  "$timeFolder 0"
            write(fid,*) " "
            write(fid,"(A)")  "$calculateCorrL 1"
            write(fid,"(A)")  "$deleteSampleAfterStatistics 0"
            write(fid,*) " "
            
            call getcwd(absPath)

            write(*,*) "absPath = ", trim(adjustL(absPath))

            write(fid,"(A)") "$nSamples "//numb2String(size(prop))

            do prop_Count = 0, size(prop) - 1
                mat_Nb = prop(prop_count)%mat
                VAR = (prop(prop_Count)%CV*prop(prop_Count)%avg)**2d0
                mesh_path = trim(string_join_many(SEM_gen_path,"/input/", &
                             stringNumb_join("Mat_", mat_Nb),"_",prop(prop_Count)%name,"_mesh"))
                gen_path  = trim(string_join_many(SEM_gen_path,"/input/", &
                             stringNumb_join("Mat_", mat_Nb),"_",prop(prop_Count)%name,"_gen"))
                write(*,*) "mesh_path = ", mesh_path
                write(*,*) "gen_path = ", gen_path
  
                write(fid,"(A)") trim(string_join_many(stringNumb_join("$mesh_input_", prop_Count+1)))//' "'//&
                                 trim(string_join_many(absPath,"/",mesh_path,'"'))
                call write_mesh_file(3, prop(prop_Count)%bbox_min(:), prop(prop_Count)%bbox_max(:), [5, 5, 5], &
                                     mesh_path)
  
  
                write(fid,"(A)") trim(string_join_many(stringNumb_join("$gen_input_", prop_Count+1)))//' "'//&
                             trim(string_join_many(absPath,'/',gen_path,'"'))
                call write_gen_file(3, 1, prop(prop_Count)%corrMod, prop(prop_Count)%margiF, prop(prop_Count)%corrL, &
                                    prop(prop_Count)%avg, VAR, 4, &
                                    prop(prop_Count)%seedStart, [5d0, 5d0, 5d0], &
                                    gen_path,  &
                                    1, [0, 0, 0])
  
                out_folder = trim(string_join_many(stringNumb_join("$out_folder_", prop_Count+1)))//&
                                  ' "'//trim(adjustL(absPath))//"/"//trim(adjustL(SEM_gen_path))//'"'
                write(fid,"(A)") trim(adjustL(out_folder))
  
  
                write(fid,"(A)") trim(string_join_many(stringNumb_join("$out_name_", prop_Count+1)))//' "'//&
                             trim(stringNumb_join("Mat_", mat_Nb))//"_"//&
                             trim(prop(prop_Count)%name)//'"'
                write(fid,"(A)") " "
  
            end do
                write(*,*) "AFTER 5"

            close(fid)

            if(allocated(prop)) deallocate(prop)
            if(allocated(nProp_Mat)) deallocate(nProp_Mat)
            if(allocated(bbox_max)) deallocate(bbox_max)
            if(allocated(bbox_min)) deallocate(bbox_min)


        end subroutine generateMain_inputSEM

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
                    if(IPT%rang==0) write(*,*) "   Mesh UNV"
                    call wLog("    Mesh UNV")
                    call wLog("        file: ")
                    call read_DataTable(dataTable, "unv_path", IPT%unv_path)
                    call read_DataTable(dataTable, "pointsPerCorrL", IPT%pointsPerCorrL)
                    !call wLog(trim(IPT%unv_path))
                    if(IPT%rang==0) write(*,*) "   unv_path = ", IPT%unv_path
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
        subroutine estimate_nFields(IPT)

            implicit none
            !OUTPUT
            type(IPT_RF), intent(inout)  :: IPT
            !LOCAL
            integer, dimension(IPT%nDim_gen) :: nPoints, nBlocks
            integer, dimension(IPT%nDim_gen) :: nFieldsIdeal, nBlocksIdeal, nFieldsChosen
            integer, dimension(IPT%nDim_gen) :: nPointsBase, nPointsOvlp, ratio
            !double precision, dimension(IPT%nDim_gen) :: nBlocsDble
            logical :: nFieldsOK
            integer, dimension(100) :: factors
            integer :: i, pos, np, np_total, np_start, np_end
            double precision :: vol_surf_factor, vol_surf_factor_temp

            if(any(IPT%nFields <=0)) then

                if(IPT%rang == 0)  write(*,*) "IPT%nFields will be decided automatically"

                nPointsBase = IPT%pointsPerCorrL * &
                              ceiling((IPT%xMaxGlob_in - IPT%xMinGlob_in)/IPT%corrL_in)
                nPointsOvlp = IPT%pointsPerCorrL * &
                              ceiling(IPT%overlap_in/IPT%corrL_in)
                nBlocks = 1

                nFieldsOK = .false.

                !write(*,*) "nBlocks = ", nBlocks
                !write(*,*) "nPointsOvlp = ", nPointsOvlp

                !write(*,*) "BEFORE WHILE"

                do while (.not. nFieldsOK)
                    nPoints = nPointsBase + (nBlocks - 1)*2*nPointsOvlp
                    !write(*,*) "nPoints = ", nPoints

                    !nBlocks = ceiling(nPoints/200)

                    ratio = ceiling(dble(nPoints)/dble(nBlocks))
                    !write(*,*) "ratio = ", ratio

                    where(nPoints/nBlocks > 200) nBlocks = nBlocks + 1
                    !write(*,*) "nBlocks = ", nBlocks

                    ratio = ceiling(dble(nPoints)/dble(nBlocks))

                    if(all(ratio < 200)) nFieldsOK = .true.
                end do

                nBlocksIdeal = nBlocks

                if(IPT%rang == 0) then
                    if(IPT%rang == 0)  write(*,*) "nBlocks (ideal) = ", nBlocksIdeal
                end if

                np_total = IPT%nb_procs
                vol_surf_factor_temp = 0D0

                np_start = ceiling(0.9*np_total)
                np_end   = np_total

                do np = np_start, np_end


                    nBlocks = nBlocksIdeal

                    factors = 0
                    nFieldsIdeal  = 1
                    !call find_factors(IPT%nb_procs, factors)
                    call find_factors(np, factors)
                    !if(IPT%rang == 0) write(*,*) "factors = ", factors
                    !Adapting to our number of processors


                    do i = size(factors), 1, -1
                        !if(IPT%rang == 0) write(*,*) "i = ", i

                        if(factors(i) == 0) cycle

                        if(all(nBlocks == 1)) exit

                        pos = MAXLOC(nBlocks,1)
                        nFieldsIdeal(pos) = nFieldsIdeal(pos)*factors(i)
                        nBlocks(pos) = ceiling(dble(nBlocks(pos))/dble(factors(i)))

                    end do

                    vol_surf_factor_temp = dble(product(nFieldsIdeal))/dble(2*sum(nFieldsIdeal*CSHIFT(nFieldsIdeal, shift=1)))

                    !if(IPT%rang == 0) write(*,*) "nBlocks (ideal) = ", nBlocks
                    !if(IPT%rang == 0) write(*,*) "nFieldsIdeal (ideal) = ", nFieldsIdeal
                    !if(IPT%rang == 0) write(*,*) "product(nFieldsIdeal)", product(nFieldsIdeal)
                    !if(IPT%rang == 0) write(*,*) "sum(nFieldsIdeal**2)", sum(nFieldsIdeal*CSHIFT(nFieldsIdeal, shift=1))
                    !if(IPT%rang == 0) write(*,*) "Vol/Surf = ", vol_surf_factor_temp

                    if(vol_surf_factor_temp >= vol_surf_factor) then
                        vol_surf_factor = vol_surf_factor_temp
                        nFieldsChosen   = nFieldsIdeal
                    end if
                end do

                IPT%nFields = nFieldsChosen

                if(IPT%rang == 0) write(*,*) "vol_surf_factor = ", vol_surf_factor
            end if

            if(IPT%rang == 0) write(*,*) "IPT%nFields = ", IPT%nFields

        end subroutine estimate_nFields

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
                !write(unit,*) " outputFolder = ", trim(adjustL(IPT%outputFolder))
                !write(unit,*) " outputName   = ", trim(adjustL(IPT%outputName))

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

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_mesh_file(nDim, xMinGlob, xMaxGlob, pointsPerCorrL, mesh_path)

        implicit none
        !INPUT
        integer, intent(in) :: nDim
        double precision, dimension(:), intent(in) :: xMinGlob, xMaxGlob
        integer, dimension(:), intent(in) :: pointsPerCorrL
        character(len=*), intent(in) :: mesh_path
        !LOCAL
        integer :: i
        integer :: fileId

        fileID = 19

        open (unit = fileId , file = mesh_path, action = 'write')

        write(fileId,*) "$$nDim ", nDim
        write(fileId,*) "$$meshMod 1"
        write(fileId,*) "          $Min            $Max           $pointsPerCorrL"
        do i = 1, nDim
            write(fileId, "(2(F15.5, A), (I15))") xMinGlob(i), " ", xMaxGlob(i), " ", pointsPerCorrL(i)
        end do

        close(fileId)

        call system("chmod a+r "//trim(mesh_path))

    end subroutine write_mesh_file

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine write_gen_file(nDim, Nmc, corrMod, margiFirst, corrL, fieldAvg, fieldVar, method, &
                              seedStart, overlap, gen_path,                         &
                              localizationLevel, nFields)

        implicit none
        !INPUT
        integer, intent(in) :: nDim, Nmc, corrMod, margiFirst, method, seedStart
        double precision, dimension(:), intent(in) :: corrL, overlap
        double precision, intent(in) :: fieldAvg, fieldVar
        character(len=*), intent(in) :: gen_path
        integer, intent(in) :: localizationLevel
        integer, dimension(:), intent(in) :: nFields
        !LOCAL
        integer :: fileId

        fileID = 19

        open (unit = fileId , file = gen_path, action = 'write')

        write(fileId,*) "$$nDim ", nDim
        write(fileId,*) "$$Nmc ", Nmc
        write(fileId,*) "$$corrMod ", corrMod
        write(fileId,*) "$$margiFirst ", margiFirst
        write(fileId,*) "$$localizationLevel ", localizationLevel
        write(fileId,*) "$nFields "
        write(fileId,*) nFields
        write(fileId,*) "$corrL "
        write(fileId,*) corrL
        write(fileId,*) "$$fieldAvg "
        write(fileId,*) fieldAvg
        write(fileId,*) "$$fieldVar "
        write(fileId,*) fieldVar
        write(fileId,*) "$$method "
        write(fileId,*) method
        write(fileId,*) "$$seedStart"
        write(fileId,*) seedStart
        !write(fileId,*) "$$independent"
        !write(fileId,*) independent
        write(fileId,*) "$overlap"
        write(fileId,*) overlap

        close(fileId)

        call system("chmod a+r "//trim(gen_path))

    end subroutine write_gen_file

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine find_factors(n, d)
        integer, intent(in) :: n
        integer, dimension(:), intent(out) :: d

        integer :: div, next, rest
        integer :: i

        i = 1
        div = 2; next = 3; rest = n

        do while ( rest /= 1 )
           do while ( mod(rest, div) == 0 )
              d(i) = div
              i = i + 1
              rest = rest / div
           end do
           div = next
           next = next + 2
        end do

  end subroutine find_factors

end module type_inputRF
