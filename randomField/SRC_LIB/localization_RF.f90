module localization_RF

    use displayCarvalhol
    use math_RF
    use constants_RF
    use mpi
    use write_Log_File
    use type_RF
    use type_MESH
    use common_variables_RF
    use randomFieldND
    use type_inputRF

    implicit none

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine getNeighIndexRange(MSH, minIndexNeigh, maxIndexNeigh)
        !INPUT
        type(MESH), intent(in) :: MSH
        !OUTPUT
        integer, intent(out) :: minIndexNeigh, maxIndexNeigh
        integer :: i
        !logical, dimension(:), intent(out) :: considerNeighbour

        !considerNeighbour = .true.
        !where(MSH%neigh < 0) considerNeighbour = .false.
        !where(minval(MSH%neighShift,1) < 0) considerNeighbour = .false.

        !Global Min and Max positions
        !minIndexNeigh = minval(pack(MSH%indexNeigh(1,:), MSH%considerNeighbour))
        !maxIndexNeigh = maxval(pack(MSH%indexNeigh(2,:), MSH%considerNeighbour))

        minIndexNeigh = -1
        maxIndexNeigh = -1
        do i = 1, size(MSH%indexNeigh,2)
            if(MSH%considerNeighbour(i)) then
                if(minIndexNeigh == -1) then
                    minIndexNeigh = MSH%indexNeigh(1,i)
                    maxIndexNeigh = MSH%indexNeigh(2,i)
                end if
                if(minIndexNeigh > MSH%indexNeigh(1,i)) minIndexNeigh = MSH%indexNeigh(1,i)
                if(maxIndexNeigh < MSH%indexNeigh(2,i)) maxIndexNeigh = MSH%indexNeigh(2,i)
            end if
        end do

    end subroutine getNeighIndexRange

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine find_nSamples_in_proc(IPT, nSamplesInProc, xNStep_Proc, xNTotal_Proc)

        implicit none

        !INPUT
        type(IPT_RF), intent(in)  :: IPT

        !OUTPUT
        integer, intent(out) :: nSamplesInProc
        integer(kind=8), intent(out) :: xNTotal_Proc
        integer, dimension(IPT%nDim), intent(out) :: xNStep_Proc

        !LOCAL
        integer :: rest, code
        integer :: sum_SamplesInProc, nSamplesInAllProc

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

            if(IPT%rang == 0) then
                write(*,*) 'nSamplesInProc = ', nSamplesInProc
                write(*,*) 'xNStep_Proc    = ', xNStep_Proc
                write(*,*) 'xNTotal_Proc   = ', xNTotal_Proc
            end if

            call wLog("nSamplesInProc = ")
            call wLog(nSamplesInProc)
            call wLog("xNStep_Proc = ")
            call wLog(xNStep_Proc)
            call wLog("xNTotal_Proc = ")
            call wLog(xNTotal_Proc)
        end if

        !Verification
        call MPI_ALLREDUCE (nSamplesInProc, sum_SamplesInProc, 1, MPI_INTEGER, MPI_SUM, IPT%comm, code)
        if(IPT%rang == 0) then
            write(*,*) 'sum_SamplesInProc = ', sum_SamplesInProc
            write(*,*) 'IPT%nTotalFields  = ', IPT%nTotalFields
            write(*,*) 'Check = ', (IPT%nTotalFields == sum_SamplesInProc)
        end if
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

    end subroutine find_nSamples_in_proc

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine prepare_Localization(IPT, subdivisionCoords, subdivisionId, &
                                    xMinFiles, xMaxFiles, xMin_Group, xMax_Group, &
                                    gen_GroupRange, xNStep_Group, origin_Group, &
                                    xNTotal_Group)

        implicit none

        !INPUT
        type(IPT_RF), intent(in)  :: IPT

        !OUTPUT
        double precision, dimension(:,:), intent(out) :: subdivisionCoords
        integer, dimension(:,:), intent(out) :: subdivisionId
        double precision, dimension(:,:), intent(out) :: xMinFiles, xMaxFiles
        double precision, dimension(:), intent(out) :: xMin_Group, xMax_Group
        double precision, dimension(:), intent(out) :: gen_GroupRange
        integer, dimension(IPT%nDim), intent(out) :: xNStep_Group, origin_Group
        integer(kind=8), intent(out) :: xNTotal_Group

        !LOCAL
        integer :: countFields, i
        integer, dimension(IPT%nDim) :: locStep

        !Init
        locStep = IPT%nFields**(IPT%localizationLevel)
        call setGrid(subdivisionCoords, IPT%xMinGlob, IPT%stepProc, locStep, inverse=.true.)
        do i = 1, size(subdivisionCoords, 2)
            subdivisionId(:,i) = nint((subdivisionCoords(:,i)-IPT%xMinGlob)/IPT%stepProc)
        end do
        if(IPT%rang == 0) write(*,*) "Max Coord = ", subdivisionCoords(:, size(subdivisionCoords,2)) + IPT%procExtent

        ! PREPARING INTERNAL LOCALIZATION-----------------------------------
        if(IPT%rang == 0) write(*,*) " "
        if(IPT%rang == 0) write(*,*) "-> PREPARING INTERNAL LOCALIZATION----------------------------------------"
        if(IPT%loc_group == 0) then

            countFields  = 0
            do i = 1, IPT%nTotalFields
                if(all(subdivisionId(:,i)/(IPT%nFields**(IPT%localizationLevel - 1)) == IPT%coords) &
                   .or. (.not. IPT%extLoc)) then
                    countFields = countFields + 1
                    xMinFiles(:, countFields) = subdivisionCoords(:,i)
                    xMaxFiles(:, countFields) = xMinFiles(:, countFields) + IPT%procExtent
                end if
            end do

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

        end if

    end subroutine prepare_Localization


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine add_RF_to_Group(IPT, randField_Gen, xNStep_Proc, &
                               unityPartition, xMinFile, &
                               xMin_Group, RF_2D_Group, RF_3D_Group)

        implicit none

        !INPUT
        type(IPT_RF), intent(in)  :: IPT
        integer, dimension(:), intent(in) :: xNStep_Proc
        double precision, dimension(:), intent(in) :: unityPartition
        double precision, dimension(:), intent(in) :: xMinFile, xMin_Group

        !OUTPUT
        double precision, dimension(:,:), intent(inout), target :: randField_Gen
        double precision, dimension(:, :), pointer, intent(inout), optional :: RF_2D_Group
        double precision, dimension(:, :, :), pointer, intent(inout), optional :: RF_3D_Group


        !LOCAL
        double precision, dimension(:, :), pointer :: RF_2D_Gen
        double precision, dimension(:, :, :), pointer :: RF_3D_gen
        integer, dimension(IPT%nDim) :: minP, maxP
        integer :: j

        if(IPT%gen_rang == 0) then

            do j = 1, IPT%Nmc

                !Multiplication
                randField_Gen(:,j) = randField_Gen(:,j)*unityPartition

                !Sum
                minP = find_xNStep(xMin_Group, xMinFile, IPT%xStep)
                maxP = minP + xNStep_Proc - 1

                if(IPT%nDim_gen == 2) then
                    if(.not. present(RF_2D_Group)) stop("RF_2D_Group not present inside add_RF_to_Group")
                    RF_2D_Gen(1:xNStep_Proc(1),1:xNStep_Proc(2)) => randField_Gen(:,j)

                    if(size(randField_Gen(:,j)) /= product(shape(RF_2D_Group(minP(1):maxP(1),&
                                                                   minP(2):maxP(2)))))&
                 stop("RF_2D_Gen and randField_Gen don't have consistent sizes inside add_RF_to_Group")

                    RF_2D_Group(minP(1):maxP(1),minP(2):maxP(2)) = RF_2D_Gen &
                        + RF_2D_Group(minP(1):maxP(1),minP(2):maxP(2))
                else if(IPT%nDim_gen == 3) then
                    if(.not. present(RF_3D_Group)) stop("RF_3D_Group not present inside add_RF_to_Group")
                    RF_3D_Gen(1:xNStep_Proc(1),1:xNStep_Proc(2),1:xNStep_Proc(3)) => randField_Gen(:,j)

                    if(size(randField_Gen(:,j)) /= product(shape(RF_3D_Group(minP(1):maxP(1),&
                                                                   minP(2):maxP(2),&
                                                                   minP(3):maxP(3)))))&
                 stop("RF_3D_Gen and randField_Gen don't have consistent sizes inside add_RF_to_Group")

                    RF_3D_Group(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3)) = RF_3D_Gen &
                        + RF_3D_Group(minP(1):maxP(1),minP(2):maxP(2),minP(3):maxP(3))
                end if
            end do

        end if

        if(associated(RF_2D_Gen))  nullify(RF_2D_Gen)
        if(associated(RF_3D_Gen))  nullify(RF_3D_Gen)

    end subroutine add_RF_to_Group

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine addNeighboursFieldsV3(randFieldGroup, xNStep, overlap, corrL, xStep,&
                                     nDim, neigh, op_neigh, neighShift, xNTotal, &
                                     rang, loc_comm)
        implicit none

        !INPUT OUTPUT
        double precision, dimension(:,:), intent(inout) , target :: randFieldGroup

        !INPUT
        integer, dimension(:), intent(in) ::  xNStep
        double precision, dimension(:)  , intent(in) :: overlap, corrL, xStep
        integer, intent(in) :: nDim
        integer, dimension(:), intent(in) :: neigh, op_neigh
        integer, dimension(:,:), intent(in) :: neighShift
        integer, intent(in) :: rang, loc_comm
        integer(kind=8), intent(in) :: xNTotal

        !LOCAL
        integer :: direction, myRank
        integer :: neighRank
        integer, dimension(nDim) :: dirShift
        integer :: op_direction
        integer(kind=8) :: nOvlpMax
        integer :: code
        integer, dimension(nDim) :: minPos, maxPos
        integer, dimension(nDim) :: overlapNPoints
        integer :: double_size
        integer(kind=8) :: totalSize, overHead, overEst, bufferSize
        double precision, dimension(:), allocatable :: buffer
        !integer :: testrank = 1, testrank2 = 11
        double precision, dimension(1:xNTotal), target :: tempRandField, tempRandField_Now
        integer :: request
        integer, dimension(MPI_STATUS_SIZE) :: status
        integer :: tag
        logical :: snd, rcv
        integer :: BufDT_size
        double precision, dimension(:, :), pointer :: TRF_2D, TRF_2D_Now, RF_2D
        double precision, dimension(:, :, :), pointer :: TRF_3D, TRF_3D_Now, RF_3D


        if(nDim == 2) RF_2D(1:xNStep(1),1:xNStep(2)) => randFieldGroup(:,1)
        if(nDim == 3) RF_3D(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => randFieldGroup(:,1)


        call wLog ("Inside addNeighboursFieldsV3")

        !Initialization
        tempRandField = 0 !Contributions over the random field
        tag = 0
        myRank = rang
        overlapNPoints = nint((overlap*corrL - 2*xStep)/xStep) + 1
        where(overlapNPoints < 1) overlapNPoints = 0

        nOvlpMax = 0
        !Set Max Number of Points on Overlap (buffer size)
        do direction = 1, size(neigh)
            dirShift   = neighShift(:, direction)
            !if(neigh(direction) < 0) cycle
            minPos(:) = 1
            maxPos(:) = xNStep
            where(dirShift == 1 ) minPos = xNStep - overlapNPoints + 1
            where(dirShift == -1) maxPos = overlapNPoints
            totalSize = (product(int(maxPos - minPos + 1,8)))
            nOvlpMax = nOvlpMax + totalSize
        end do
        !nOvlpMax = maxval(overlapNPoints)*&
        !           nint((maxval(dble(xNStep))**(dble(nDim-1)))) !Max hypercube to send

        !write(*,*) "xStep = ", xStep
        !write(*,*) "corrL = ", corrL
        !write(*,*) "overlap = ", overlap
        !write(*,*) "overlapNPoints = ", overlapNPoints
        !write(*,*) "maxval(overlapNPoints) = ", maxval(overlapNPoints)
        !write(*,*) "maxval(dble(xNStep))**(dble(nDim-1) = ", maxval(dble(xNStep))**(dble(nDim-1))
        !write(*,*) "nOvlpMax = ", nOvlpMax

        !Buffer allocation
        call wLog ("Allocating buffer")
        overEst = 2
        !overEst = nint(2.0D0**(dble(nDim)))*10 !<- Parameter not fully understood
        call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,double_size,code)
        overHead = int(1+(MPI_BSEND_OVERHEAD)/double_size)
        bufferSize = overEst*(nOvlpMax+overHead)
        call wLog ("bufferSize = ")
        call wLog (bufferSize)
        allocate(buffer(bufferSize))
        call MPI_BUFFER_ATTACH(buffer, int(double_size*bufferSize),code)
        call wLog ("buffer allocation code= ")
        call wLog (code)

        !Allocation
        if(nDim == 2) TRF_2D(1:xNStep(1),1:xNStep(2)) => tempRandField
        if(nDim == 2) TRF_2D_Now(1:xNStep(1),1:xNStep(2)) => tempRandField_Now
        if(nDim == 3) TRF_3D(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => tempRandField
        if(nDim == 3) TRF_3D_Now(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => tempRandField_Now

        if(nDim == 2) then
            call wLog("    shape(TRF_2D) ")
            call wLog(shape(TRF_2D))
            call wLog("    shape(RF_2D) ")
            call wLog(shape(RF_2D))
        else if(nDim == 3) then
            call wLog("    shape(TRF_3D) ")
            call wLog(shape(TRF_3D))
            call wLog("    shape(RDF%RF_3D) ")
            call wLog(shape(RF_3D))
        end if

        do direction = 1, size(neigh)

            !write(*,*) "              FB RANG ", RDF%rang

            call MPI_BARRIER(loc_comm, code) !This is necessary so we can reduce the size of the buffer

            call wLog(" ")
            call wLog(" ")
            call wLog(" DIRECTION =========================================")
            call wLog(direction)
            call wLog("    Shift")
            call wLog(neighShift(:, direction))

            !SENDING---------------------------------------------------------------

            snd = .true.
            dirShift   = neighShift(:, direction)
            neighRank  = neigh(direction)
            minPos(:) = 1
            maxPos(:) = xNStep
            where(dirShift == 1 ) minPos = xNStep - overlapNPoints + 1
            where(dirShift == -1) maxPos = overlapNPoints
            totalSize = (product(maxPos - minPos + 1))



            if(neigh(direction) < 0) snd = .false. !Check if this direction exists

            if(snd) then
                call wLog(" ")
                call wLog(" SENDING ==============")
                !Defining slice that should be sent
                !call wLog("DIMENSIONING GOING SLICE!!!")

                !xMinDir = MSH%xMinNeigh(:, direction)
                !xMaxDir = MSH%xMaxNeigh(:, direction)
                !minPos = nint((xMinDir - MSH%xMinBound(:))/MSH%xStep) + 1
                !maxPos = nint((xMaxDir - MSH%xMinBound(:))/MSH%xStep) + 1

                tag = neighRank

                call wLog("totalSize")
                call wLog(totalSize)
                call wLog("minPos")
                call wLog(minPos)
                call wLog("maxPos")
                call wLog(maxPos)
                call wLog("dirShift")
                call wLog(dirShift)
                call wLog(" ")
                call wLog(" TAG ")
                call wLog(tag)
                call wLog("    TO  rang ")
                call wLog(neighRank)

                if(nDim == 2) then
                    call MPI_IBSEND (RF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
                            int(totalSize), MPI_DOUBLE_PRECISION, &
                            neighRank, tag, loc_comm, request, code)

                end if
                if(nDim == 3) then
                    call MPI_IBSEND (RF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)), &
                        int(totalSize), MPI_DOUBLE_PRECISION, &
                        neighRank, tag, loc_comm, request, code)
                end if
            else

                call wLog(" not sending ===========")

            end if

            !RECEIVING---------------------------------------------------------------
            rcv = .true.
            dirShift   = -dirShift
            neighRank  = op_neigh(direction)
            minPos(:) = 1
            maxPos(:) = xNStep
            where(dirShift == 1 ) minPos = xNStep - overlapNPoints + 1
            where(dirShift == -1) maxPos = overlapNPoints
            totalSize = (product(maxPos - minPos + 1))

            !write(*,*) "              FC RANG ", RDF%rang

            if(op_neigh(direction) < 0) rcv = .false. !Check if this direction exists

            !Finding opposite direction index
            do op_direction = 1, size(neigh)
                if(all(dirShift == neighShift(:, op_direction))) exit
            end do

            call wLog("op_direction = ")
            call wLog(op_direction)

            if(rcv) then
                call wLog(" ")
                call wLog(" RECEIVING ==============")
                !Defining slice that should be sent
                tag       = myRank

                call wLog("totalSize")
                call wLog(totalSize)
                call wLog("minPos")
                call wLog(minPos)
                call wLog("maxPos")
                call wLog(maxPos)
                call wLog("dirShift")
                call wLog(dirShift)
                call wLog(" ")
                call wLog(" TAG ")
                call wLog(tag)
                call wLog("    FROM  rang ")
                call wLog(neighRank)

                tempRandField_Now = 0.0D0

                if(nDim == 2) then
                    call MPI_RECV (TRF_2D_Now(minPos(1):maxPos(1),minPos(2):maxPos(2)), &
                        int(totalSize), MPI_DOUBLE_PRECISION, &
                        neighRank, tag, loc_comm, status, code)
                        TRF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)) = &
                        TRF_2D(minPos(1):maxPos(1),minPos(2):maxPos(2))   &
                        + TRF_2D_Now(minPos(1):maxPos(1),minPos(2):maxPos(2))
                end if
                if(nDim == 3) then
                    call MPI_RECV (TRF_3D_Now(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)), &
                        int(totalSize), MPI_DOUBLE_PRECISION, &
                        neighRank, tag, loc_comm, status, code)
                        TRF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3)) = &
                        TRF_3D(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3))   &
                        + TRF_3D_Now(minPos(1):maxPos(1),minPos(2):maxPos(2),minPos(3):maxPos(3))
                end if
            else
                call wLog(" not receiving ==============")
            end if
        end do

        randFieldGroup(:,1) = randFieldGroup(:,1) + tempRandField

        BufDT_size = int(double_size*bufferSize)
        call MPI_BUFFER_DETACH (buffer,BufDT_size,code)
        if(allocated(buffer)) deallocate(buffer)

        if(associated(TRF_2D_Now)) nullify(TRF_2D_Now)
        if(associated(TRF_3D_Now)) nullify(TRF_3D_Now)
        if(associated(TRF_2D)) nullify(TRF_2D)
        if(associated(TRF_3D)) nullify(TRF_3D)


    end subroutine addNeighboursFieldsV3


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine applyWeightingFunctions_OnMatrix(RDF, MSH, partitionType)
    !INPUT
    type(RF), intent(in) :: RDF
    type(MESH), intent(in) :: MSH
    integer, intent(in) :: partitionType

    !LOCAL
    integer :: direction, i
    double precision, dimension(1:MSH%xNTotal) :: unityPartition
    double precision, dimension(MSH%nDim) :: originCorner
    integer, dimension(MSH%nDim) :: minPos, maxPos

    !call wLog("shape(RDF%RF_2D) = ")
    !call wLog(shape(RDF%RF_2D))
    !call wLog("shape(RDF%xPoints_2D) = ")
    !call wLog(shape(RDF%xPoints_2D))
    !call wLog("RDF%xPoints_2D(:,1,1) = ")
    !call wLog(RDF%xPoints_2D(:,1,1))
    !call wLog("RDF%xPoints_2D(:,30,29) = ")
    !call wLog(RDF%xPoints_2D(:,30,29))

    call wLog("Inside 'applyWeightingFunctions_OnMatrix'")


    !Modify extremes of local Random Field-------------------------------------------------------
    unityPartition(:) = 1

    !Building Shape Functions in all directions
    do direction = 1, size(MSH%neigh)

        !if(.not. MSH%considerNeighbour(direction)) cycle !Don't consider Neighbours in this direction
        if(MSH%neigh(direction) < 0) cycle
        if(sum(abs(MSH%neighShift(:, direction))) /= 1) cycle !We need only the "pure directions"

        !if(direction/=4) cycle !FOR TESTS

        !Finding range
        originCorner = MSH%xOrNeigh(:, direction)
        minPos = nint((MSH%xMinNeigh(:, direction) - MSH%xMinBound(:))/MSH%xStep) + 1
        maxPos = nint((MSH%xMaxNeigh(:, direction) - MSH%xMinBound(:))/MSH%xStep) + 1

        call wLog("DIRECTION")
        call wLog(direction)
        call wLog("originCorner")
        call wLog(originCorner)
        call wLog("MSH%neighShift(:, direction)")
        call wLog(MSH%neighShift(:, direction))
        call wLog("minPos = ")
        call wLog(minPos)
        call wLog("maxPos")
        call wLog(maxPos)

        !Shape Function Generation
        call generateUnityPartition_OnMatrix(RDF, MSH%xNStep, originCorner, MSH%overlap, &
            MSH%corrL, MSH%neighShift(:, direction), partitionType, &
            unityPartition, minPos, maxPos)

    end do !Direction

    !if(RDF%rang == 0) call DispCarvalhol(unityPartition, "unityPartition")

    !if(RDF%rang == 0) write(*,*) "Showing Unity Partition -------------------------------"

    do i = 1, size(RDF%randField(:,1),1)
        RDF%randField(i,1) = RDF%randField(i,1) * sqrt(unityPartition(i))
        !if(RDF%rang == 0 .and. (unityPartition(i) /= 1.0D0))  write(*,*) unityPartition(i)
    end do

    !if(RDF%rang == 0) write(*,*) "That's ALL -------------------------------"
    !RDF%randField(:,1) = RDF%randField(:,1) * sqrt(unityPartition(:)) !This syntax demands too much memory for Occygen

    !RDF%randField(:,1) = unityPartition(:) !TEST
    !RDF%randField(:,1) = sqrt(unityPartition(:)) !TEST

end subroutine applyWeightingFunctions_OnMatrix

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
subroutine generateUnityPartition_OnMatrix(RDF, xNStep, originCorner, overlap, corrL, neighShift, &
    partitionType, unityPartition, minPos, maxPos)

    implicit none

    !INPUT
    type(RF), intent(in) ::  RDF
    integer, dimension(:), intent(in) ::  xNStep
    double precision, dimension(:)  , intent(in) :: originCorner, overlap, corrL
    integer, dimension(:)  , intent(in) :: neighShift
    integer, intent(in) :: partitionType
    integer, dimension(:), intent(in) :: minPos, maxPos

    !OUTPUT
    double precision, dimension(:), intent(out), target :: unityPartition

    !LOCAL
    integer :: i
    double precision, dimension(:, :), pointer :: UP_2D
    double precision, dimension(:, :, :), pointer :: UP_3D

    if(RDF%nDim == 2) UP_2D(1:xNStep(1),1:xNStep(2)) => unityPartition
    if(RDF%nDim == 3) UP_3D(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => unityPartition

    !call wLog("shape(unityPartition) = ")
    !call wLog(shape(unityPartition))
    !call wLog("shape(UP_2D) = ")
    !call wLog(shape(UP_2D))
    !call wLog("shape(RDF%xPoints_2D) = ")
    !call wLog(shape(RDF%xPoints_2D))
    !call wLog("shape(UP_3D) = ")
    !call wLog(shape(UP_3D))
    !call wLog("shape(RDF%xPoints_3D) = ")
    !call wLog(shape(RDF%xPoints_3D))

    !call wLog("Point in minimal position = ")
    !call wLog(RDF%xPoints_3D(:,22,1,1))
    !call wLog(RDF%xPoints_3D(:,minPos(1),minPos(2),minPos(3)))



    do i = 1, RDF%nDim
        !do i = 1, 1 !For Tests
        if(neighShift(i) == 0) cycle
        !            call wLog("Dimension = ")
        !call wLog(i)

        if(partitionType == 1) then
            if(RDF%nDim == 2) then
                !                    call wLog("Point in minimal position = ")
                !                    call wLog(RDF%xPoints_2D(i,minPos(1),minPos(2)))
                !                    call wLog("Point in maximal position = ")
                !                    call wLog(RDF%xPoints_2D(i,maxPos(1),maxPos(2)))
                !                    call wLog("originCorner(i) = ")
                !                    call wLog(originCorner(i))
                !UP_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)) = 1 !For Tests
                UP_2D(minPos(1):maxPos(1),minPos(2):maxPos(2)) = &
                    ((1.0D0 + cos(PI*(RDF%xPoints_2D(&
                    i,minPos(1):maxPos(1),minPos(2):maxPos(2)) &
                    - originCorner(i))/(overlap(i)*corrL(i))))&
                    / 2.0D0) &
                    * UP_2D(minPos(1):maxPos(1),minPos(2):maxPos(2))
            else if (RDF%nDim == 3) then
                !                    call wLog("Point in minimal position = ")
                !                    call wLog(RDF%xPoints_3D(i,minPos(1),minPos(2),minPos(3)))
                !                    call wLog("Point in maximal position = ")
                !                    call wLog(RDF%xPoints_3D(i,maxPos(1),maxPos(2),maxPos(3)))
                !                    call wLog("originCorner(i) = ")
                !                    call wLog(originCorner(i))
                UP_3D(minPos(1):maxPos(1),minPos(2):maxPos(2), minPos(3):maxPos(3)) = &
                    ((1.0D0 + cos(PI*(RDF%xPoints_3D(&
                    i,minPos(1):maxPos(1),minPos(2):maxPos(2), minPos(3):maxPos(3)) &
                    - originCorner(i))/(overlap(i)*corrL(i))))&
                    / 2.0D0) &
                    * UP_3D(minPos(1):maxPos(1),minPos(2):maxPos(2), minPos(3):maxPos(3))
            else
                stop("Unity Partition not implemented in this Dimension (generateUnityPartition_OnMatrix)")
            end if
        else
        stop('ERROR!! Inside "generateUnityPartition_OnMatrix" - partition Type not defined')
    end if
end do

if(associated(UP_2D)) nullify(UP_2D)
if(associated(UP_3D)) nullify(UP_3D)

end subroutine generateUnityPartition_OnMatrix

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
subroutine generateUnityPartition_Matrix(xNStep, overlap, corrL, xStep,&
                                         partitionType, unityPartition, nDim, &
                                         neigh, neighShift, reverse)

    implicit none

    !INPUT
    integer, dimension(:), intent(in) ::  xNStep
    double precision, dimension(:)  , intent(in) :: overlap, corrL, xStep
    integer, intent(in) :: partitionType
    integer, intent(in) :: nDim
    integer, dimension(:), intent(in), optional :: neigh
    integer, dimension(:,:), intent(in), optional :: neighShift
    logical, intent(in), optional :: reverse
    !OUTPUT
    double precision, dimension(:), intent(out), target :: unityPartition

    !LOCAL
    integer :: i, j, k, s
    double precision, dimension(:, :), pointer :: UP_2D
    double precision, dimension(:, :, :), pointer :: UP_3D
    integer, dimension(nDim) :: overlapNPoints, overlapSizeInt
    double precision, dimension(nDim) :: originCorner
    double precision, dimension(:), allocatable :: pattern
    integer, dimension(nDim) :: U_Lim, D_Lim
    logical :: considerNeighbour

    if(nDim == 2) UP_2D(1:xNStep(1),1:xNStep(2)) => unityPartition
    if(nDim == 3) UP_3D(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => unityPartition

    !call wLog("shape(unityPartition) = ")
    !call wLog(shape(unityPartition))
    !call wLog("shape(UP_2D) = ")
    !call wLog(shape(UP_2D))
    !call wLog("shape(RDF%xPoints_2D) = ")
    !call wLog(shape(RDF%xPoints_2D))
    !call wLog("shape(UP_3D) = ")
    !call wLog(shape(UP_3D))
    !call wLog("shape(RDF%xPoints_3D) = ")
    !call wLog(shape(RDF%xPoints_3D))

    overlapNPoints = nint((overlap*corrL - 2*xStep)/xStep) + 1
    overlapSizeInt = overlapNPoints + 2
    unityPartition = 1.0D0

    call wLog("xNStep =")
    call wLog(xNStep)
    call wLog("size(unityPartition) =")
    call wLog(size(unityPartition))

    do j = 1, 2

        call wLog("j =================")
        call wLog(j)

        if(j == 1) then
            D_Lim = 1
            U_Lim = overlapNPoints
            originCorner = dble(D_Lim - 1)
        else
            U_Lim = xNStep
            D_Lim = xNStep - overlapNPoints + 1
            originCorner = dble(U_Lim + 1)
        end if

        call wLog("D_Lim = ")
        call wLog(D_Lim)
        call wLog("U_Lim = ")
        call wLog(U_Lim)
        call wLog("originCorner = ")
        call wLog(originCorner)

        do i = 1, nDim

            if(overlapNPoints(i) < 1) cycle

            considerNeighbour = .true.

            if(present(neigh)) then
                considerNeighbour = .false.
                if(j==1 .and. any(neighShift(i,:) == -1 .and. neigh >=0)) then
                    considerNeighbour = .true.
                else if(j==2 .and. any(neighShift(i,:) == 1 .and. neigh >=0)) then
                    considerNeighbour = .true.
                end if

                if(reverse) considerNeighbour = .not.considerNeighbour
            end if

            if(.not. considerneighbour) cycle

            if(allocated(pattern)) deallocate(pattern)
            allocate(pattern(overlapNPoints(i)))
            call wLog("size(pattern) = ")
            call wLog(size(pattern))
            pattern = [(dble(k), k=1, size(pattern))]
            !call wLog("Pattern Before 1")
            !call wLog(pattern)
            pattern = abs((pattern)/dble(overlapSizeInt(i)))
            !call wLog("Pattern Before 1.5")
            !call wLog(pattern)
            if(j == 1) pattern = pattern(size(pattern):1:-1) !Pattern inversion

            !call wLog("Pattern Before 2 =")
            !call wLog(pattern)

            if(partitionType == 1) then
                pattern = ((1.0D0 + cos(PI*(pattern)))&
                            / 2.0D0)
            end if

            call wLog("Pattern Developed =")
            call wLog(pattern)

            if(nDim == 2) then

                if(i==1) then
                    do s = 1, size(pattern)
                        UP_2D(D_Lim(i)+s-1,:) = UP_2D(D_Lim(i)+s-1,:) * pattern(s)
                    end do
                else if(i==2) then
                    do s = 1, size(pattern)
                        UP_2D(:,D_Lim(i)+s-1) = UP_2D(:,D_Lim(i)+s-1) * pattern(s)
                    end do
                end if



            else if(nDim == 3) then

                if(i==1) then
                    do s = 1, size(pattern)
                        UP_3D(D_Lim(i)+s-1,:,:) = UP_3D(D_Lim(i)+s-1,:,:) * pattern(s)
                    end do
                else if(i==2) then
                    do s = 1, size(pattern)
                        UP_3D(:,D_Lim(i)+s-1,:) = UP_3D(:,D_Lim(i)+s-1,:) * pattern(s)
                    end do
                else if(i==3) then
                    do s = 1, size(pattern)
                        UP_3D(:,:,D_Lim(i)+s-1) = UP_3D(:,:,D_Lim(i)+s-1) * pattern(s)
                    end do
                end if

            end if
        end do
    end do

    if(associated(UP_2D)) nullify(UP_2D)
    if(associated(UP_3D)) nullify(UP_3D)
    if(allocated(pattern)) deallocate(pattern)

end subroutine generateUnityPartition_Matrix

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
subroutine generateUnityPartition(xPoints, originCorner, overlap, neighShift, partitionType, unityPartition)

    implicit none

    !INPUT
    double precision, dimension(:,:), intent(in) ::  xPoints
    double precision, dimension(:)  , intent(in) :: originCorner, overlap
    integer, dimension(:)  , intent(in) :: neighShift
    integer, intent(in) :: partitionType

    !OUTPUT
    double precision, dimension(:), intent(out) :: unityPartition

    !LOCAL
    integer :: i, nDim

    nDim = size(originCorner)


    unityPartition = 1.0D0
    do i = 1, nDim
        if(neighShift(i) == 0) cycle

        if(partitionType == 1) then
            unityPartition = ((1.0D0 + cos(PI*(xPoints(i,:) - originCorner(i))/overlap(i)))&
                / 2.0D0) &
                * unityPartition
        else
            stop('ERROR!! Inside "generateUnityPartition" - partition Type not defined')
        end if
    end do

end subroutine generateUnityPartition

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
subroutine reorderRandomFieldStruct(RDF, MSH)
    implicit none

    type(RF)  , intent(inout) :: RDF
    type(MESH), intent(in)    :: MSH

    !LOCAL
    integer, dimension(MSH%nDim) :: ind3D, offset
    double precision, dimension(MSH%nDim) :: orig, coordValue, tempCoordValue, xStep
    double precision, dimension(RDF%Nmc)  :: tempRFValue, RFValue
    integer :: ind1D
    integer :: i

    orig  = (dble(MSH%origin-1)*MSH%xStep + MSH%xMinGlob)
    xStep = MSH%xStep

    call wLog("xStep")
    call wLog(xStep)
    call wLog("orig")
    call wLog(orig)
    call wLog("MSH%xNStep")
    call wLog(MSH%xNStep)

    offset(1) = 1
    do i = 2, MSH%nDim
        offset(i) = product(MSH%xNStep(1:i-1))
    end do

    call wLog("offset")
    call wLog(offset)
    call wLog("shape(RDF%randField)")
    call wLog(shape(RDF%randField))

    do i = 1, size(RDF%randField,1)
        coordValue = RDF%xPoints(:,i)
        RFValue    = RDF%randField(i,:)
        ind3D = nint((coordValue-orig)/xStep)
        ind1D = sum(ind3D*offset) + 1

        !The point is not where it was supposed to be
        do while (ind1D /= i)
            !Finding index
            ind3D = nint((coordValue-orig)/xStep)
            ind1D = sum(ind3D*offset) + 1
            !Saving temp data
            tempRFValue    = RDF%randField(ind1D,:)
            tempCoordValue = RDF%xPoints(:, ind1D)
            !Replacing data
            RDF%randField(ind1D,:) = RFvalue
            RDF%xPoints(:, ind1D)  = coordValue
            !Going to the next coordinate (the one that was in the index we took)
            RFValue     = tempRFValue
            coordValue  = tempCoordValue
        end do


    end do

end subroutine reorderRandomFieldStruct

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
function findTag(MSH, neighPos, direction, send) result(tag)

    implicit none
    !INPUT
    type(MESH), intent(in) :: MSH
    integer, intent(in) :: neighPos, direction
    logical, intent(in) :: send
    integer :: sinal

    !OUTPUT
    integer :: tag

    !LOCAL
    integer :: i


    if(MSH%neigh(neighPos) < 0 .or. MSH%neigh(direction) < 0) then
        !write(*,*) "Inside findTag , Invalid Neighbour"
        tag = -1
    else
        tag = 0

        if(send) then
            do i = 1, MSH%nDim
                tag = tag + MSH%neighShift(i, direction)*3**(MSH%nDim - i)
            end do
        else
            do i = 1, MSH%nDim
                sinal = 1
                if(MSH%neighShift(i, neighPos) /= 0) sinal = -1
                tag = tag + sinal*MSH%neighShift(i, direction)*3**(MSH%nDim - i)
            end do
        end if

        tag = tag + (3**(MSH%nDim)-1)/2

    end if

end function findTag

!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    !-----------------------------------------------------------------------------------------------
!    subroutine copy_RF_xPoints(MSH, orRDF, destRDF, destXPoints, direction)
!        implicit none
!
!        !INPUT AND OUTPUT
!        type(MESH), intent(in) :: MSH
!        type(RF), intent(in)   :: orRDF
!        type(RF) ::destRDF
!        integer, intent(in) :: direction
!        double precision, dimension(:, :), allocatable, intent(out), target :: destXPoints;
!
!        !LOCAL
!        integer :: i, j, totalSize
!
!        totalSize = MSH%indexNeigh(2,direction) - MSH%indexNeigh(1,direction) + 1
!
!        if(allocated(destXPoints)) then
!            if(size(destXPoints) /= totalSize) then
!                nullify(destRDF%xPoints)
!                deallocate(destXPoints)
!            end if
!        end if
!
!        if(.not.allocated(destXPoints)) allocate(destXPoints(MSH%nDim, 1:totalSize))
!
!        destRDF%xPoints => destXPoints
!
!        destRDF%xPoints(:,:) = orRDF%xPoints(:, MSH%indexNeigh(1,direction):MSH%indexNeigh(2,direction))
!        destRDF%xNTotal = totalSize
!
!    end subroutine copy_RF_xPoints

end module localization_RF
