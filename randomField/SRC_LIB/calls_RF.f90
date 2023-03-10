module calls_RF

    use displayCarvalhol
    use write_Log_File
    use math_RF
    use constants_RF
    use mpi
    use writeResultFile_RF
    use type_RF
    use localization_RF
    use type_MESH
    use common_variables_RF
    use randomFieldND
    use mesh_RF
    use type_inputRF
    use sample_RF

    implicit none

    interface createRandomField
       module procedure create_RF_Unstruct_Init
    end interface createRandomField

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine make_random_field (IPT, times)
        implicit none
        !INPUT
        type(IPT_RF), intent(inout)  :: IPT
        !OUTPUT
        double precision, dimension(:), intent(inout) :: times
        !LOCAL
        double precision :: t_initial
        integer :: code

        call CPU_TIME(t_initial)

        call wLog("-> Processing INPUTS----------------------------------------")
        if(IPT%rang == 0) write(*,*) "-> Processing INPUTS---------------------------"

        call validate_input(IPT)

        !Initialization
        IPT%corrL   = IPT%corrL_in
        IPT%overlap = IPT%overlap_in

        call process_input(IPT, &
                           IPT%stepProc, IPT%procExtent, IPT%overlap, &
                           IPT%gen_groupMax, IPT%gen_group, IPT%gen_Comm, IPT%gen_nbProcs, &
                           IPT%gen_rang, &
                           IPT%loc_groupMax, IPT%loc_group, IPT%loc_Comm, IPT%loc_nbProcs, &
                           IPT%loc_rang, IPT%extLoc, IPT%nDim, &
                           IPT%xMinGlob, IPT%xMaxGlob, IPT%xStep, IPT%localizationLevel, &
                           IPT%nTotalFields, IPT%coords, IPT%neigh, IPT%op_neigh, IPT%neighShift, &
                           IPT%global)

        call wLog("-----LOCALIZATION---------------")
        call wLog("     IPT%loc_group = ")
        call wLog(IPT%loc_group)
        call wLog("     IPT%loc_groupMax = ")
        call wLog(IPT%loc_groupMax)
        call wLog("     IPT%loc_Comm = ")
        call wLog(IPT%loc_Comm)
        call wLog("     IPT%loc_nbProcs = ")
        call wLog(IPT%loc_nbProcs)
        call wLog("     IPT%extLoc = ")
        call wLog(IPT%extLoc)
        call wLog("     IPT%nTotalFields = ")
        call wLog(IPT%nTotalFields)
        call wLog(" ")
        call wLog("-----GENERATION---------------")
        call wLog("     IPT%gen_group = ")
        call wLog(IPT%gen_group)
        call wLog("     IPT%gen_groupMax = ")
        call wLog(IPT%gen_groupMax)
        call wLog("     IPT%gen_Comm = ")
        call wLog(IPT%gen_Comm)
        call wLog("     IPT%gen_nbProcs = ")
        call wLog(IPT%gen_nbProcs)
        call wLog("     IPT%gen_rang = ")
        call wLog(IPT%gen_rang)
        call wLog(" ")
        call wLog("-----ROUNDING---------------")
        call wLog("     IPT%xMinGlob = ")
        call wLog(IPT%xMinGlob)
        call wLog("     IPT%xMaxGlob = ")
        call wLog(IPT%xMaxGlob)
        call wLog("     IPT%localizationLevel = ")
        call wLog(IPT%localizationLevel)
        call wLog("-----TOPOLOGY---------------")
        call wLog("     IPT%coords = ")
        call wLog(IPT%coords)
        call wLog("     IPT%nDim = ")
        call wLog(IPT%nDim)
!        call wLog("     IPT%neigh = ")
!        call wLog(IPT%neigh)
!        call wLog("     IPT%neighShift = ")
!        call wLog(IPT%neighShift)
        call show_IPTneigh(IPT, "IPT-Neighbours", .false., forLog=.true.)

        if(IPT%rang == 0) call show_IPT_RF(IPT, "IPT")
        call build_random_field (IPT, times, t_initial)

        call MPI_COMM_FREE (IPT%gen_Comm, code)
        call MPI_COMM_FREE (IPT%loc_Comm, code)

    end subroutine make_random_field


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine build_random_field (IPT, times, t_initial)
        !INPUT
        type(IPT_RF), intent(in)  :: IPT
        double precision, dimension(:), intent(inout) :: times
        double precision, intent(in) :: t_initial

        !LOCAL
        double precision, dimension(:,:), allocatable :: UNV_randField
        integer               :: code
        double precision, dimension(IPT%nDim) :: gen_GroupRange
        double precision, dimension(IPT%nDim, IPT%nTotalFields) :: subdivisionCoords
        integer         , dimension(IPT%nDim, IPT%nTotalFields) :: subdivisionId
        double precision      :: t_bef
        integer               :: fieldNumber
        character(len=buf_RF) :: BBoxPath, XMFPath, MONO_FileName
        character(len=buf_RF), dimension(:), allocatable :: MONO_FileNames
        character(len=buf_RF) :: HDF5_RePath_out
        double precision, dimension(IPT%nTotalFields) :: gen_times, temp_gen_times
        integer :: i, d, countFields
        integer :: nSamplesInProc
        double precision, dimension(:,:), allocatable, target :: randField_Gen
        double precision, dimension(:,:), allocatable, target :: randField_Group
        double precision, dimension(:,:), allocatable :: randField_Local
        double precision, dimension(:), allocatable ::unityPartition
        integer(kind=8) :: xNTotal_Proc, xNTotal_Group
        double precision, dimension(:,:), allocatable :: xMinFiles, xMaxFiles
        double precision, dimension(IPT%nDim) :: xMin_Group, xMax_Group
        integer, dimension(IPT%nDim) :: xNStep_Proc, xNStep_Group, origin_Group
        double precision, dimension(:, :), pointer :: RF_2D_Group
        double precision, dimension(:, :, :), pointer :: RF_3D_Group
        double precision :: gen_WALL_Time
        double precision, dimension(8) :: build_times, BT_sum, BT2_sum
        double precision, dimension(8) :: BT_avg
        double precision, dimension(8) :: BT_stdDev, BT_max, BT_min
        double precision :: t_final
        double precision, dimension(IPT%nDim) :: kMax_out
        integer, dimension(IPT%nDim) :: kNStep_out
        integer(kind=8) :: file_bytes_size
        double precision :: file_mb_size
        character(len=buf_RF) :: absPath_HDF5, absPath_XMF

        if(IPT%rang == 0) write(*,*) "  "
        if(IPT%rang == 0) write(*,*) " Inside 'build_random_field'"
        if(IPT%rang == 0) write(*,*) "  "
        write(*,*) "-> Inside 'build_random_field---------------------------rang = ",IPT%rang

        build_times(:) = 0.0D0
        BT_sum = 0.0D0
        BT2_sum = 0.0D0
        BT_avg = 0.0D0
        BT_stdDev = 0.0D0
        BT_max = 0.0D0
        BT_min = 0.0D0
        gen_times = 0.0D0
        temp_gen_times = 0.0D0
        gen_WALL_Time = 0.0D0

        build_times(1) = MPI_Wtime() !Reference

        !Discovering number of fields in each proc
        nSamplesInProc = 0
        xNTotal_Proc   = 0

        call find_nSamples_in_proc(IPT, nSamplesInProc, xNStep_Proc, xNTotal_Proc)

        if(IPT%loc_group == 0) then
            allocate(randField_Gen(xNTotal_Proc, IPT%Nmc))
            if(IPT%write_intermediate_files) allocate(MONO_FileNames(nSamplesInProc))
            allocate(xMinFiles(IPT%nDim, nSamplesInProc))
            allocate(xMaxFiles(IPT%nDim, nSamplesInProc))
        end if

        call prepare_Localization(IPT, subdivisionCoords, subdivisionId, &
                                    xMinFiles, xMaxFiles, xMin_Group, xMax_Group, &
                                    gen_GroupRange, xNStep_Group, origin_Group, &
                                    xNTotal_Group)

        if(IPT%loc_group == 0) then
            allocate(randField_Group(xNTotal_Group, IPT%Nmc))
            randField_Group = 0.0D0
            if(IPT%nDim_gen == 2) then
                RF_2D_Group(1:xNStep_Group(1),1:xNStep_Group(2)) => randField_Group
            else if(IPT%nDim_gen == 3) then
                RF_3D_Group(1:xNStep_Group(1),1:xNStep_Group(2),1:xNStep_Group(3)) => randField_Group
            end if

            if(.not. IPT%global) then
                allocate(unityPartition(xNTotal_Proc))
                call wLog("Generating Partition of Unity")
                call generateUnityPartition_Matrix(xNStep_Proc, IPT%overlap, IPT%corrL, IPT%xStep,&
                                             1, unityPartition, IPT%nDim)
                if(IPT%write_intermediate_files) then
                    MONO_FileName = string_join_many("PofUnit_L0_Group",numb2String(IPT%gen_group))
                    call write_MONO_proc_result(IPT%procExtent*0, IPT%procExtent, &
                                                IPT%xStep, IPT%nDim, &
                                                unityPartition, MONO_FileName, IPT%outputFolder)
                end if
            end if
        end if


        !call MPI_BARRIER(IPT%comm, code)
        build_times(2) = MPI_Wtime() !Organizing Localization
        !times(3) = MPI_Wtime() !Organizing Localization

        !MAKING ALL REALIZATIONS--------------------------------------------------
        gen_times(:) = 0.0D0
        countFields  = 0
        !if(.false.)then
        if(IPT%sampleFields)then

            if(IPT%rang == 0) write(*,*) " "
            if(IPT%rang == 0) write(*,*) "-> SAMPLING----------------------------------------"
            call wLog("-> SAMPLING----------------------------------------")
            do i = 1, IPT%nTotalFields
            !do i = 1, 1 !FOR TESTS
                if(all(subdivisionId(:,i)/(IPT%nFields**(IPT%localizationLevel - 1)) == IPT%coords) &
                   .or. (.not. IPT%extLoc)) then

                    if(IPT%gen_rang == 0) write(*,*)  "-> Gen_Group ", IPT%gen_group, " making Field ", i
                    call wLog("-> Making Field")
                    call wLog(i)

                    t_bef = MPI_Wtime()
                    fieldNumber = i;
                    countFields = countFields + 1
                    write(*,*) "-> SINGLE REALIZATION proc ", IPT%rang
                    call single_realization(IPT, &
                                            IPT%gen_Comm, fieldNumber, subdivisionCoords(:,i), &
                                            randField_Local, randField_Gen, kMax_out, kNStep_out)

                    if(IPT%gen_rang == 0) then
                        if(IPT%write_intermediate_files) then
                            call wLog("Writing intermediate generation file")
                            MONO_FileNames(countFields) = "GEN"

                            do d = 1, IPT%nDim
                                MONO_FileNames(countFields) = &
                                string_join_many(MONO_FileNames(countFields),"_",numb2String(subdivisionId(d,i)+1,3))
                            end do

                            call write_MONO_proc_result(xMinFiles(:, countFields), xMaxFiles(:, countFields), &
                                                        IPT%xStep, IPT%nDim, &
                                                        randField_Gen(:,1), &
                                                        MONO_FileNames(countFields), IPT%outputFolder)
                        end if

                        if(.not. IPT%global) then
                            if(IPT%nDim == 2) then
                                call add_RF_to_Group(IPT, randField_Gen, xNStep_Proc, &
                                                unityPartition, &
                                                xMinFiles(:, countFields), &
                                                xMin_Group, RF_2D_Group=RF_2D_Group)
                           else if(IPT%nDim == 3) then
                               call add_RF_to_Group(IPT, randField_Gen, xNStep_Proc, &
                                                unityPartition, &
                                                xMinFiles(:, countFields), &
                                                xMin_Group, RF_3D_Group=RF_3D_Group)
                           end if


                           if(IPT%write_intermediate_files) then
                                call wLog("Writing intermediate localized generation file")
                                MONO_FileName = MONO_FileNames(countFields)
                                MONO_FileName = string_join_many("LOC_L0_P", numb2String(IPT%rang), "-",MONO_FileName)
                                if(IPT%gen_rang == 0) call write_MONO_proc_result( &
                                    xMinFiles(:, countFields), xMaxFiles(:, countFields), &
                                    IPT%xStep, IPT%nDim, &
                                    randField_Gen(:,1), MONO_FileName, &
                                    IPT%outputFolder)
                            end if
                        else
                            randField_Group = randField_Gen
                        end if
                   end if

                    !if(allocated(randField_Local)) deallocate(randField_Local)
                    !write(*,*) "After single"
                end if
            end do
        end if

        build_times(3) = MPI_Wtime() !Sampling

        write(*,*) "-> AFTER SAMPLING proc ", IPT%rang

        if(.not. IPT%global) then

            if(IPT%write_intermediate_files) then
                MONO_FileName = string_join_many("LOC_L1_P", numb2String(IPT%rang), &
                                                 "BEF_Cor-GROUP",numb2String(IPT%gen_group))
                if(IPT%gen_rang == 0) call write_MONO_proc_result(xMin_Group, xMax_Group, &
                                                              IPT%xStep, IPT%nDim, &
                                                              randField_Group(:,1), MONO_FileName, &
                                                              IPT%outputFolder, &
                                                              HDF5_RePath_out)
            end if

            !Correcting Borders (From internal localization)
            if(allocated(unityPartition)) deallocate(unityPartition)

            if(IPT%loc_group == 0) then
                allocate(unityPartition(xNTotal_Group))

                call wLog("shape(randField_Group(:,1)) = ")
                call wLog(shape(randField_Group(:,1)))
                call wLog("shape(unityPartition) = ")
                call wLog(shape(unityPartition))

                if(any(shape(unityPartition) /= shape(randField_Group(:,1)))) then
                    write(*,*) "shape(randField_Group(:,1)) = ", shape(randField_Group(:,1))
                    write(*,*) "shape(unityPartition) = ", shape(unityPartition)
                    write(*,*) "ERROR in correcting borders, unityPartition and randField_Group don't have the same sizes"
                end if

                call generateUnityPartition_Matrix(xNStep_Group, IPT%overlap, IPT%corrL, IPT%xStep,&
                                                   1, unityPartition, IPT%nDim, &
                                                   IPT%neigh, IPT%neighShift, reverse = .true.)

                if(IPT%write_intermediate_files) then
                    MONO_FileName = string_join_many("PofUnit_L1_Group",numb2String(IPT%gen_group))
                    if(IPT%gen_rang == 0) call write_MONO_proc_result(xMin_Group, xMax_Group, &
                                                                  IPT%xStep, IPT%nDim, &
                                                                  unityPartition, MONO_FileName, &
                                                                  IPT%outputFolder)
                end if

                call wLog("BEFmaxval(randField_Group(:,:)) = ")
                call wLog(maxval(randField_Group(:,:)))
                call wLog("BEFminval(randField_Group(:,:)) = ")
                call wLog(minval(randField_Group(:,:)))

                randField_Group(:,1) = randField_Group(:,1)/unityPartition

                if(IPT%write_intermediate_files) then
                    MONO_FileName = string_join_many("LOC_L1_P", numb2String(IPT%rang), "AFT_Cor-GROUP",numb2String(IPT%gen_group))
                    if(IPT%gen_rang == 0) call write_MONO_proc_result(xMin_Group, xMax_Group, &
                                                                  IPT%xStep, IPT%nDim, &
                                                                  randField_Group(:,1), MONO_FileName, &
                                                                  IPT%outputFolder)
                end if

                call wLog("AFTmaxval(randField_Group(:,:)) = ")
                call wLog(maxval(randField_Group(:,:)))
                call wLog("AFTminval(randField_Group(:,:)) = ")
                call wLog(minval(randField_Group(:,:)))
            end if

        end if

        if(allocated(unityPartition)) deallocate(unityPartition)

        !call MPI_BARRIER(IPT%loc_Comm, code)
        !build_times(4) = MPI_Wtime() !Internal Localization Time

         ! EXTERNAL LOCALIZATION-----------------------------------
         if(IPT%rang == 0) write(*,*) " "
         if(IPT%rang == 0) then
            write(*,*) "-> EXTERNAL LOCALIZATION----------------------------------------"
            if(IPT%extLoc) then
                write(*,*) "YES"
            else
                write(*,*) "NO"
            end if
         end if

        if(IPT%loc_group == 0 .and. IPT%extLoc) then
            !Combining realizations (communicating between procs)

            if(IPT%rang == 0) write(*,*) " "
            if(IPT%rang == 0) write(*,*) "-> COMBINING----------------------------------------"
            call wLog("-> COMBINING----------------------------------------")
             call addNeighboursFieldsV3(randField_Group, xNStep_Group, IPT%overlap, IPT%corrL, IPT%xStep,&
                                        IPT%nDim, IPT%neigh, IPT%op_neigh, IPT%neighShift, xNTotal_Group, &
                                        IPT%rang, IPT%loc_comm)

            if(IPT%write_intermediate_files) then
                MONO_FileName = string_join_many("LOC_L1_P", numb2String(IPT%rang), "AFT_Sum-GROUP",numb2String(IPT%gen_group))
                if(IPT%gen_rang == 0) call write_MONO_proc_result(xMin_Group, xMax_Group, &
                                                                  IPT%xStep, IPT%nDim, &
                                                                  randField_Group(:,1), MONO_FileName, &
                                                                  IPT%outputFolder)
            end if
        end if

        if(IPT%write_intermediate_files) then
            MONO_FileName = string_join_many("EXT_LOC_L1_P", numb2String(IPT%rang), "-GROUP",numb2String(IPT%gen_group))
            if(IPT%gen_rang == 0) call write_MONO_proc_result(xMin_Group, xMax_Group, &
                                                          IPT%xStep, IPT%nDim, &
                                                          randField_Group(:,1), MONO_FileName, &
                                                          IPT%outputFolder)
        end if

        !call MPI_BARRIER(IPT%loc_Comm, code)
        build_times(4:8) = MPI_Wtime() !External Localization Time
        times(6) = build_times(4) !External Localization Time

        ! TRANSFORMATION AND OUTPUT WRITING

        if(IPT%loc_group == 0) then
            if(IPT%rang == 0) write(*,*) "-> TRANFORMING AND WRITING OUTPUT--------------------"
            call wLog("-> TRANFORMING AND WRITING OUTPUT--------------------")
            !Normalizing and Writing files
            call transform_and_write_output(randField_Group, xNStep_Group, origin_Group, &
                                            IPT, build_times, BBoxPath, XMFPath)
            if(IPT%rang == 0) write(*,*) "BBoxPath = ", trim(BBoxPath)
            if(IPT%rang == 0) write(*,*) "fileExist = ", fileExist (BBoxPath)

            if(IPT%rang == 0) write(*,*) "-> WRITING ATTRIBUTES--------------------"
            call wLog("-> WRITING ATTRIBUTES--------------------")
            if(IPT%rang == 0) call write_HDF5_attributes(BBoxPath, &
                           IPT%nb_procs, IPT%nDim, IPT%Nmc, IPT%method, IPT%seedStart, &
                           IPT%corrMod, IPT%margiFirst, &
                           IPT%localizationLevel, IPT%nFields, &
                           IPT%xMinGlob, IPT%xMaxGlob, IPT%xStep, IPT%corrL, IPT%overlap, &
                           IPT%procExtent, kMax_out, kNStep_out, .false.)
        end if

        if(IPT%write_intermediate_files) then
            MONO_FileName = string_join_many("RF_",numb2String(IPT%gen_group))
            if(IPT%gen_rang == 0) call write_MONO_proc_result(xMin_Group, xMax_Group, &
                                                          IPT%xStep, IPT%nDim, &
                                                          randField_Group(:,1), MONO_FileName, &
                                                          IPT%outputFolder)


        write(*,*) "IPT%gen_groupMax = ", IPT%gen_groupMax

        if(IPT%rang == 0) call write_XMF_Elements_per_proc(trim(string_join(MONO_FileName,".h5")), xMin_Group, xMax_Group, IPT%xStep, &
                                     IPT%nDim, trim(MONO_FileName), trim(string_join(IPT%outputFolder,"/xmf")), &
                                     "../h5", "RF",IPT%gen_groupMax)
        end if

        if(IPT%rang == 0) call write_stat_input("./stat_input", BBoxPath, IPT%calculateCorrL, IPT%deleteSampleAfterStatistics)

        if(allocated(unityPartition))   deallocate(unityPartition)
        if(allocated(randField_Gen)) deallocate(randField_Gen)
        if(allocated(randField_Local))  deallocate(randField_Local)
        if(allocated(randField_Group))  deallocate(randField_Group)
        if(associated(RF_2D_Group)) nullify(RF_2D_Group)
        if(associated(RF_3D_Group)) nullify(RF_3D_Group)

        !Build Times
        !1 - Reference
        !2 - Organizing
        !3 - Sampling
        !4 - Internal Localization
        !5 - External Localization
        !6 - Normalizing
        !7 - Transformation
        !8 - Writing Files
        !9 - UNV Interpolation


        ! WRITING INTERPOLATION FILE-----------------------------------------------
        if(IPT%unv .and. IPT%writeUNVinterpolation .and. IPT%outputStyle == 1) then
            call MPI_BCAST (BBoxPath, len(BBoxPath), MPI_CHARACTER, 0, IPT%comm, code)
            if(IPT%rang == 0) write(*,*) "-> WRITING INTERPOLATION FILE----------------------------"
            if(IPT%rang == 0) write(*,*) "-> Writing 'UNV' XMF and hdf5 files for"
            if(IPT%rang == 0) write(*,*) trim(adjustL(IPT%unv_path))
            allocate(UNV_randField(size(IPT%coordList,2),1))
            if(IPT%rang == 0) write(*,*) "  Source:"
            if(IPT%rang == 0) write(*,*) trim(adjustL(BBoxPath))
            if(IPT%rang == 0) write(*,*) "-> INTERPOLATING TO GIVEN MESH----------------------------------------"
            call wLog("-> INTERPOLATING TO GIVEN MESH----------------------------------------")
            call interpolateToMesh(BBoxPath, IPT%coordList, UNV_randField, IPT%rang)
            call write_UNV_XMF_h5(UNV_randField, IPT%coordList, IPT%connectList, &
                                  "UNV_", IPT%rang, IPT%outputFolder, &
                                  IPT%comm, 0)
            if(allocated(UNV_randField)) deallocate(UNV_randField)
        end if

        build_times(8) = MPI_Wtime() !Interpolation Time

        times(3:9) = build_times(2:8)

        BT_sum  = 0.0D0 !Just a temporary variable
        BT2_sum = 0.0D0 !Just a temporary variable
        BT_sum(2:8)  = build_times(2:8)
        BT2_sum(2:8) = build_times(1:7)
        build_times = BT_sum - BT2_sum

        !BT_stdDev = -1.0D0 !TEST
        !BT_avg = -1.0D0 !TEST

        call MPI_REDUCE (build_times, BT_sum, size(build_times), MPI_DOUBLE_PRECISION, MPI_SUM, &
                         0, IPT%comm,code)
        call MPI_REDUCE (build_times**2.0D0, BT2_sum, size(build_times), MPI_DOUBLE_PRECISION, MPI_SUM, &
                         0, IPT%comm,code)
        call MPI_REDUCE (temp_gen_times, gen_times, size(gen_times), MPI_DOUBLE_PRECISION, MPI_SUM, &
                         0, IPT%comm,code)

        !BT_min = -1.0D0 !TEST
        !BT_max = -1.0D0 !TEST
        do i = 1, size(build_times)
            call MPI_REDUCE (build_times(i), BT_min(i), 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
                             0, IPT%comm,code)
            call MPI_REDUCE (build_times(i), BT_max(i), 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
                             0, IPT%comm,code)
        end do

        call CPU_TIME(t_final)
        gen_WALL_Time = t_final - t_initial

        call wLog("build_times = ")
        call wLog(build_times)
        call wLog("t_initial = ")
        call wLog(t_initial)
        call wLog("t_final = ")
        call wLog(t_final)
        call wLog("gen_WALL_Time = ")
        call wLog(gen_WALL_Time)

        if(IPT%rang == 0) then

            BT_avg    = BT_sum/dble(IPT%nb_procs)
            BT_stdDev = sqrt(BT2_sum/dble(IPT%nb_procs) - BT_avg**2.0D0)

            call wLog("Build Times (BT)")
            call wLog("|1 - Reference")
            call wLog("|2 - Organizing")
            call wLog("|3 - Sampling")
            call wLog("|4 - Internal Localization")
            call wLog("|5 - External Localization")
            call wLog("|6 - Normalizing")
            call wLog("|7 - Transformation")
            call wLog("|8 - Writing Files")
            call wLog("|9 - UNV Interpolation")

            call wLog("BT_min = ")
            call wLog(BT_min)
            call wLog("BT_max = ")
            call wLog(BT_max)
            call wLog("BT_avg    = ")
            call wLog(BT_avg)
            call wLog("BT_stdDev = ")
            call wLog(BT_stdDev)

            call write_HDF5_time_attributes(BBoxPath, &
                                            BT_avg, BT_stdDev, BT_min, BT_max, &
                                            gen_times, gen_WALL_Time)

            call getcwd(MONO_FileName)
            !write(*,*) "MONO_FileName(len(trim(MONO_FileName)):len(trim(MONO_FileName))) "
            !write(*,*) MONO_FileName(len(trim(MONO_FileName)):len(trim(MONO_FileName)))
            write(*,*) "PATH: ", trim(MONO_FileName)

            if(BBoxPath(1:1) == "/") then
                absPath_HDF5 = BBoxPath
            else
                absPath_HDF5 = string_join_many(MONO_FileName, "/", BBoxPath)
            end if

            if(XMFPath(1:1) == "/") then
                absPath_XMF = XMFPath
            else
                absPath_XMF = string_join_many(MONO_FileName, "/", XMFPath)
            end if

            write(*,*) "OUTPUT HDF5 ON: ", trim(absPath_HDF5)
            inquire(FILE=trim(absPath_HDF5), SIZE=file_bytes_size)
            file_mb_size = dble(file_bytes_size)/dble(1024.0D0 ** 2.0D0)
            write(*,*) "    file_mb_size: ", file_mb_size
            write(*,*) "OUTPUT XMF  ON: ", trim(absPath_XMF)
        end if


        if(allocated(xMinFiles))        deallocate(xMinFiles)
        if(allocated(xMaxFiles))        deallocate(xMaxFiles)
        if(allocated(unityPartition))   deallocate(unityPartition)
        if(allocated(randField_Gen)) deallocate(randField_Gen)
        !if(allocated(randField_inProc)) deallocate(randField_inProc)
        if(allocated(randField_Local))  deallocate(randField_Local)
        if(allocated(randField_Group))  deallocate(randField_Group)
        if(associated(RF_2D_Group)) nullify(RF_2D_Group)
        if(associated(RF_3D_Group)) nullify(RF_3D_Group)
        !if(associated(RF_2D_Gen))  nullify(RF_2D_Gen)
        !if(associated(RF_3D_Gen))  nullify(RF_3D_Gen)

    end subroutine build_random_field

end module calls_RF
