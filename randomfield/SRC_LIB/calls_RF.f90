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

        call CPU_TIME(t_initial)

        if(IPT%rang == 0) write(*,*) "-> Processing INPUTS----------------------------------------"
        call wLog("-> Processing INPUTS----------------------------------------")

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
                           IPT%nTotalFields, IPT%coords, IPT%neigh, IPT%op_neigh, IPT%neighShift)

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
        double precision      :: t_bef, t_aft
        integer               :: fieldNumber
        character(len=200) :: BBoxPath, XMFPath, MONO_FileName
        character(len=200), dimension(:), allocatable :: MONO_FileNames
        double precision, dimension(IPT%nTotalFields) :: gen_times, temp_gen_times
        integer :: i, d, countFields, j
        integer :: nSamplesInProc, nSamplesInAllProc, rest, sum_SamplesInProc
        double precision, dimension(:,:,:), allocatable, target :: randField_inProc
        double precision, dimension(:,:), allocatable, target :: randField_Group
        double precision, dimension(:,:), allocatable :: randField_Local
        double precision, dimension(:), allocatable ::unityPartition
        integer(kind=8) :: xNTotal_Proc, xNTotal_Group
        double precision, dimension(:,:), allocatable :: xMinFiles, xMaxFiles
        double precision, dimension(IPT%nDim) :: ones, xMin_Group, xMax_Group
        integer, dimension(IPT%nDim) :: xNStep_Proc, xNStep_Group, origin_Group
        integer, dimension(IPT%nDim) :: locStep, minP, maxP
        double precision, dimension(:, :), pointer :: RF_2D_Proc, RF_2D_Group
        double precision, dimension(:, :, :), pointer :: RF_3D_Proc, RF_3D_Group
        double precision :: gen_WALL_Time
        double precision, dimension(9) :: build_times, BT_sum, BT2_sum
        double precision, dimension(9) :: BT_avg
        double precision, dimension(9) :: BT_stdDev, BT_max, BT_min
        double precision :: t_final
        double precision, dimension(IPT%nDim) :: kMax_out
        integer, dimension(IPT%nDim) :: kNStep_out

        if(IPT%rang == 0) write(*,*) "  "
        if(IPT%rang == 0) write(*,*) " Inside 'build_random_field'"
        if(IPT%rang == 0) write(*,*) "  "

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
        !BT_stdDev = -1.0D0 !TEST

        !Init
        ones = 1.0D0
        locStep = IPT%nFields**(IPT%localizationLevel)
        call setGrid(subdivisionCoords, IPT%xMinGlob, IPT%stepProc, locStep, inverse=.true.)
        do i = 1, size(subdivisionCoords, 2)
            subdivisionId(:,i) = nint((subdivisionCoords(:,i)-IPT%xMinGlob)/IPT%stepProc)
        end do
        if(IPT%rang == 0) call DispCarvalhol(subdivisionCoords, "subdivisionCoords")
        !if(IPT%rang == 0) call DispCarvalhol(subdivisionId, "subdivisionId")
        if(IPT%rang == 0) write(*,*) "Max Coord = ", subdivisionCoords(:, size(subdivisionCoords,2)) + IPT%procExtent

        !Discovering number of fields in each proc
        nSamplesInProc = 0
        xNTotal_Proc   = 0

        if(IPT%gen_rang == 0) then
            nSamplesInAllProc = int(IPT%nTotalFields/(IPT%loc_nbProcs))
            rest              = IPT%nTotalFields - (nSamplesInAllProc*IPT%loc_nbProcs)
            nSamplesInProc    = nSamplesInAllProc
            if(IPT%rang < rest) nSamplesInProc = nSamplesInProc + 1

            xNStep_Proc = find_xNStep(xMaxExt=IPT%procExtent, xStep=IPT%xStep)
            xNTotal_Proc = product(int(xNStep_Proc, 8))

            allocate(randField_inProc(xNTotal_Proc, IPT%Nmc, nSamplesInProc))

            if(IPT%write_intermediate_files) allocate(MONO_FileNames(nSamplesInProc))
            allocate(xMinFiles(IPT%nDim, nSamplesInProc))
            allocate(xMaxFiles(IPT%nDim, nSamplesInProc))

        end if


        !Verification
        call MPI_ALLREDUCE (nSamplesInProc, sum_SamplesInProc, 1, MPI_INTEGER, MPI_SUM, IPT%comm, code)
        if(IPT%nTotalFields /= sum_SamplesInProc) then
            write(*,*) "ERROR in 'build_random_field' IPT%nTotalFields and sum_SamplesInProc are different"
            write(*,*) "IPT%nTotalFields  = ", IPT%nTotalFields
            write(*,*) "sum_SamplesInProc = ", sum_SamplesInProc
            call wLog("ERROR in 'build_random_field' IPT%nTotalFields and sum_SamplesInProc are different")
            call wLog("nSamplesInProc = ")
            call wLog(nSamplesInProc)
            call wLog("IPT%nTotalFields = ")
            call wLog(IPT%nTotalFields)
            call wLog("sum_SamplesInProc = ")
            call wLog(sum_SamplesInProc)
            stop (" ")
        end if

        call wLog("IPT%gen_rang = ")
        call wLog(IPT%gen_rang)
        call wLog("nSamplesInProc = ")
        call wLog(nSamplesInProc)
        call wLog("xNTotal_Proc = ")
        call wLog(xNTotal_Proc)
        call wLog("IPT%coords = ")
        call wLog(IPT%coords)


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
                !if(mod(i, gen_groupMax) == gen_group) then
                if(all(subdivisionId(:,i)/(IPT%nFields**(IPT%localizationLevel - 1)) == IPT%coords) &
                   .or. (.not. IPT%extLoc)) then

                    if(IPT%gen_rang == 0) write(*,*)  "-> Gen_Group ", IPT%gen_group, " making Field ", i
                    call wLog("-> Making Field")
                    call wLog(i)

                    t_bef = MPI_Wtime()
                    fieldNumber = i;
                    countFields = countFields + 1
                    call single_realization(IPT, &
                                            IPT%gen_Comm, fieldNumber, subdivisionCoords(:,i), &
                                            IPT%stepProc, randField_Local, kMax_out, kNStep_out)
                    call wLog("Gathering Sample")
                    call gather_sample(randField_Local, randField_inProc(:,:,countFields), &
                                       IPT%gen_rang, IPT%gen_nbProcs, IPT%gen_comm)

                    if(IPT%gen_rang == 0) then
                        xMinFiles(:, countFields) = subdivisionCoords(:,i)
                        xMaxFiles(:, countFields) = xMinFiles(:, countFields) + IPT%procExtent
                        if(IPT%write_intermediate_files) then
                            call wLog("Writing intermediate generation file")
                            MONO_FileNames(countFields) = "GEN"

                            do d = 1, IPT%nDim
                                MONO_FileNames(countFields) = &
                                string_join_many(MONO_FileNames(countFields),"_",numb2String(subdivisionId(d,i)+1,3))
                            end do

                            call write_MONO_proc_result(xMinFiles(:, countFields), xMaxFiles(:, countFields), &
                                                        IPT%xStep, IPT%nDim, &
                                                        randField_inProc(:,1,countFields), &
                                                        MONO_FileNames(countFields), single_path)
                        end if
                    end if
                    if(allocated(randField_Local)) deallocate(randField_Local)
                    t_aft = MPI_Wtime()
                    temp_gen_times(i) = t_aft-t_bef
                    !write(*,*) "After single"
                end if
            end do
        end if

        build_times(3) = MPI_Wtime() !Sampling

        ! INTERNAL LOCALIZATION-----------------------------------
        if(IPT%rang == 0) write(*,*) " "
        if(IPT%rang == 0) write(*,*) "-> INTERNAL LOCALIZATION----------------------------------------"

        if(IPT%loc_group == 0) then

            !Localization Inside Group
            if(IPT%rang == 0) write(*,*) "Internal Localization"
            xMin_Group    = minval(xMinFiles(:, :),2)
            xMax_Group    = maxval(xMaxFiles(:, :),2)
            gen_GroupRange = xMax_Group - xMin_Group
            xNStep_Group  = find_xNStep(xMaxExt=gen_GroupRange, xStep=IPT%xStep)
            origin_Group  = find_xNStep(xMinExt=IPT%xMinGlob, xMaxExt=xMin_Group, xStep=IPT%xStep)
            xNTotal_Group = product(int(xNStep_Group,8))
            call wLog("xMin_Group = ")
            call wLog(xMin_Group)
            call wLog("xMax_Group = ")
            call wLog(xMax_Group)
            call wLog("gen_GroupRange = ")
            call wLog(gen_GroupRange)
            call wLog("xNStep_Group = ")
            call wLog(xNStep_Group)
            call wLog("xNTotal_Group = ")
            call wLog(xNTotal_Group)

            allocate(randField_Group(xNTotal_Group, IPT%Nmc))
            randField_Group = 0.0D0

            if(IPT%nDim_gen == 2) then
                RF_2D_Group(1:xNStep_Group(1),1:xNStep_Group(2)) => randField_Group
            else if(IPT%nDim_gen == 3) then
                RF_3D_Group(1:xNStep_Group(1),1:xNStep_Group(2),1:xNStep_Group(3)) => randField_Group
            end if

            allocate(unityPartition(xNTotal_Proc))

            call wLog("Generating Partition of Unity")
            call generateUnityPartition_Matrix(xNStep_Proc, IPT%overlap, IPT%corrL, IPT%xStep,&
                                         1, unityPartition, IPT%nDim)
            if(IPT%write_intermediate_files) then
                MONO_FileName = string_join_many("PofUnit_L0_Group",numb2String(IPT%gen_group))
                call write_MONO_proc_result(IPT%procExtent*0, IPT%procExtent, &
                                            IPT%xStep, IPT%nDim, &
                                            unityPartition, MONO_FileName, single_path)
            end if

            if(any(shape(unityPartition) /= shape(randField_inProc(:,1,1)))) then
                write(*,*) "ERROR in internal localization, unityPartition and randField_inProc don't have the same sizes"
            end if

            do i = 1, nSamplesInProc
                !Multiplication
                do j = 1, IPT%Nmc
                    randField_inProc(:,j,i) = randField_inProc(:,j,i)*unityPartition
                    if(IPT%write_intermediate_files) then
                        MONO_FileName = MONO_FileNames(i)
                        MONO_FileName = string_join_many("LOC_L0_P", numb2String(IPT%rang), "-",MONO_FileName)
                        if(IPT%gen_rang == 0) call write_MONO_proc_result(xMinFiles(:, i), xMaxFiles(:, i), &
                                                                      IPT%xStep, IPT%nDim, &
                                                                      randField_inProc(:,j,i), MONO_FileName, &
                                                                      single_path)

                    end if
                end do

                !Sum
                minP = find_xNStep(xMin_Group, xMinFiles(:, i), IPT%xStep)
                maxP = minP + xNStep_Proc - 1

                if(IPT%nDim_gen == 2) then
                    RF_2D_Proc(1:xNStep_Proc(1),1:xNStep_Proc(2)) => randField_inProc(:,1,i)

                    RF_2D_Group(minP(1):maxP(1),minP(2):maxP(2)) = RF_2D_Proc &
                                                       + RF_2D_Group(minP(1):maxP(1),minP(2):maxP(2))
                else if(IPT%nDim_gen == 3) then
                    RF_3D_Proc(1:xNStep_Proc(1),1:xNStep_Proc(2),1:xNStep_Proc(3)) => randField_inProc(:,1,i)

                    RF_3D_Group(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3)) = RF_3D_Proc &
                                              + RF_3D_Group(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3))
                end if

                if(associated(RF_2D_Proc)) nullify(RF_2D_Proc)
                if(associated(RF_3D_Proc)) nullify(RF_3D_Proc)

            end do

            if(allocated(unityPartition)) deallocate(unityPartition)

            if(IPT%write_intermediate_files) then
                MONO_FileName = string_join_many("LOC_L1_P", numb2String(IPT%rang), "BEF_Cor-GROUP",numb2String(IPT%gen_group))
                if(IPT%gen_rang == 0) call write_MONO_proc_result(xMin_Group, xMax_Group, &
                                                              IPT%xStep, IPT%nDim, &
                                                              randField_Group(:,1), MONO_FileName, &
                                                              single_path)
            end if

            !Correcting Borders
            allocate(unityPartition(xNTotal_Group))
            call generateUnityPartition_Matrix(xNStep_Group, IPT%overlap, IPT%corrL, IPT%xStep,&
                                               1, unityPartition, IPT%nDim, &
                                               IPT%neigh, IPT%neighShift, reverse = .true.)

            if(IPT%write_intermediate_files) then
                MONO_FileName = string_join_many("PofUnit_L1_Group",numb2String(IPT%gen_group))
                if(IPT%gen_rang == 0) call write_MONO_proc_result(xMin_Group, xMax_Group, &
                                                              IPT%xStep, IPT%nDim, &
                                                              unityPartition, MONO_FileName, &
                                                              single_path)
            end if

            call wLog("shape(randField_Group(:,1)) = ")
            call wLog(shape(randField_Group(:,1)))
            call wLog("shape(unityPartition) = ")
            call wLog(shape(unityPartition))

            if(any(shape(unityPartition) /= shape(randField_Group(:,1)))) then
                write(*,*) "ERROR in internal localization, unityPartition and randField_Group don't have the same sizes"
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
                                                              single_path)
            end if

            call wLog("AFTmaxval(randField_Group(:,:)) = ")
            call wLog(maxval(randField_Group(:,:)))
            call wLog("AFTminval(randField_Group(:,:)) = ")
            call wLog(minval(randField_Group(:,:)))

            if(allocated(unityPartition)) deallocate(unityPartition)

        end if !END INTERNAL LOCALIZATION

        !call MPI_BARRIER(IPT%loc_Comm, code)
        build_times(4) = MPI_Wtime() !Internal Localization Time

         ! EXTERNAL LOCALIZATION-----------------------------------
         if(IPT%rang == 0) write(*,*) " "
         if(IPT%rang == 0) write(*,*) "-> EXTERNAL LOCALIZATION----------------------------------------"

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
                                                                  single_path)
            end if
        end if

        !call MPI_BARRIER(IPT%loc_Comm, code)
        build_times(5:9) = MPI_Wtime() !External Localization Time
        times(6) = build_times(5) !External Localization Time

        ! TRANSFORMATION AND OUTPUT WRITING
        if(IPT%rang == 0) write(*,*) "-> TRANFORMING AND WRITING OUTPUT--------------------"
        call wLog("-> TRANFORMING AND WRITING OUTPUT--------------------")

        if(IPT%loc_group == 0) then
            !Normalizing and Writing files
            call transform_and_write_output(randField_Group, xNStep_Group, origin_Group, &
                                            IPT, build_times, BBoxPath, XMFPath)
        end if


        if(IPT%rang == 0) call write_stat_input("./stat_input", BBoxPath)

        !call MPI_BARRIER(IPT%loc_Comm, code)

        if(allocated(unityPartition))   deallocate(unityPartition)
        if(allocated(randField_inProc)) deallocate(randField_inProc)
        if(allocated(randField_Local))  deallocate(randField_Local)
        if(allocated(randField_Group))  deallocate(randField_Group)
        if(associated(RF_2D_Proc))  nullify(RF_2D_Proc)
        if(associated(RF_2D_Group)) nullify(RF_2D_Group)
        if(associated(RF_3D_Proc))  nullify(RF_3D_Proc)
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
            if(IPT%rang == 0) write(*,*) BBoxPath
            if(IPT%rang == 0) write(*,*) "-> INTERPOLATING TO GIVEN MESH----------------------------------------"
            call wLog("-> INTERPOLATING TO GIVEN MESH----------------------------------------")
            call interpolateToMesh(BBoxPath, IPT%coordList, UNV_randField, IPT%rang)
            call write_UNV_XMF_h5(UNV_randField, IPT%coordList, IPT%connectList, &
                                  "UNV_", IPT%rang, single_path, &
                                  IPT%comm, 0)
            if(allocated(UNV_randField)) deallocate(UNV_randField)
        end if

        build_times(9) = MPI_Wtime() !Interpolation Time

        times(3:10) = build_times(2:9)

        BT_sum  = 0.0D0 !Just a temporary variable
        BT2_sum = 0.0D0 !Just a temporary variable
        BT_sum(2:9)  = build_times(2:9)
        BT2_sum(2:9) = build_times(1:8)
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

            call getcwd(MONO_FileName)
            !write(*,*) "PATH: ", MONO_FileName
            write(*,*) "OUTPUT HDF5 ON: ", trim(string_join_many(MONO_FileName, BBoxPath(2:)))
            write(*,*) "OUTPUT XMF  ON: ", trim(string_join_many(MONO_FileName, XMFPath(2:)))

            call write_HDF5_attributes(BBoxPath, &
                IPT%nb_procs, IPT%nDim, IPT%Nmc, IPT%method, IPT%seedStart, &
                IPT%corrMod, IPT%margiFirst, &
                BT_avg, BT_stdDev, BT_min, BT_max, gen_times, gen_WALL_Time, &
                IPT%localizationLevel, IPT%nFields, &
                IPT%xMinGlob, IPT%xMaxGlob, IPT%xStep, IPT%corrL, IPT%overlap, &
                IPT%procExtent, kMax_out, kNStep_out)
        end if


        if(allocated(xMinFiles))        deallocate(xMinFiles)
        if(allocated(xMaxFiles))        deallocate(xMaxFiles)
        if(allocated(unityPartition))   deallocate(unityPartition)
        if(allocated(randField_inProc)) deallocate(randField_inProc)
        if(allocated(randField_Local))  deallocate(randField_Local)
        if(allocated(randField_Group))  deallocate(randField_Group)
        if(associated(RF_2D_Proc))  nullify(RF_2D_Proc)
        if(associated(RF_2D_Group)) nullify(RF_2D_Group)
        if(associated(RF_3D_Proc))  nullify(RF_3D_Proc)
        if(associated(RF_3D_Group)) nullify(RF_3D_Group)

    end subroutine build_random_field
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine create_RF_Unstruct_noInit (xPoints, corrL, corrMod, Nmc,   &
!                                          randField, method, seedStart,   &
!                                          margiFirst, fieldAvg, fieldVar, &
!                                          comm, rang, nb_procs, calculate, MSH)
!        !INPUT
!        double precision, dimension(1:, 1:), intent(in), target :: xPoints;
!        double precision, dimension(1:)    , intent(in) :: corrL;
!        integer                            , intent(in) :: corrMod;
!        integer                            , intent(in) :: Nmc;
!        integer                            , intent(in) :: method
!        integer                            , intent(in) :: seedStart
!        integer                            , intent(in) :: margiFirst;
!        double precision                   , intent(in) :: fieldAvg
!        double precision                   , intent(in) :: fieldVar;
!        integer                            , intent(in) :: comm, rang, nb_procs
!        logical, dimension(1:), optional   , intent(in) :: calculate
!        type(MESH), intent(inout) :: MSH
!
!        !OUTPUT
!        double precision, dimension(:, :), intent(out), target :: randField;
!
!        !LOCAL
!        type(RF) :: RDF
!
!        write(*,*) "Inside create_RF_Unstruct_noInit"
!
!        !Initializing RF
!        call init_RF(RDF, size(corrL), Nmc, comm, rang, nb_procs)
!        RDF%xPoints   => xPoints
!        RDF%randField => randField
!        RDF%xNTotal    = size(RDF%xPoints, 2)
!        RDF%corrL      = corrL
!        RDF%corrMod    = corrMod
!        RDF%Nmc        = Nmc
!        RDF%method     = method
!        RDF%seedStart  = seedStart
!        RDF%margiFirst = margiFirst
!        RDF%fieldAvg   = fieldAvg
!        RDF%fieldVar   = fieldVar
!        if(present(calculate)) RDF%calculate  = calculate
!
!        call create_RF_Unstruct_Init(RDF, MSH)
!
!    end subroutine create_RF_Unstruct_noInit

end module calls_RF
