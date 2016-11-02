module readFile_RF
    use displayCarvalhol

    interface read_DataTable
       module procedure read_DataTable_DbleVec,  &
           read_DataTable_DbleScal, &
           read_DataTable_IntVec,   &
           read_DataTable_IntScal,  &
           read_DataTable_CharVec,  &
           read_DataTable_CharScal, &
           read_DataTable_LogicalScal
    end interface read_DataTable

contains
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_DataTable(path, dataTable, commentDef, tagIDDef, wordsMax, tagPatternMax)
        implicit none

        !Fill up data table with the data (double) of path.

        !INPUT
        character (len=*)  , intent(in) :: path;
        character, optional, intent(in) :: commentDef;
        character, optional, intent(in) :: tagIDDef;
        integer,   optional, intent(in) :: wordsMax;
        integer,   optional, intent(in) :: tagPatternMax;
        !OUTPUT
        character(len=1024), dimension(:,:), allocatable, intent(out) :: dataTable;

        !LOCAL VARIABLES
        integer            :: fileID, contentSize, unitTags;
        character (len=1024) :: empty='';
        integer            :: i, j, stat, wdMax, conStart, commentCount,   &
            tagCount, tagTotal, blankTotal,ratioTagData,         &
            dataRow, dataColumn, dataCount,dataTotal, tagPatMax;
        logical            :: posTaken, dataPassed, labelPassed;
        character          :: comment, tagID;
        character (len=1024), dimension(500)               :: foundedTags
        integer,            dimension(:),   allocatable :: tagPattern, lastLine;
        character (len=1024), dimension(:),   allocatable :: contentVector;


        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Treating optional arguments
        if(.not.present(commentDef))    comment   = '!'
        if     (present(commentDef))    comment   = commentDef
        if(.not.present(tagIDDef))      tagID     = '$'
        if     (present(tagIDDef))      tagID     = tagIDDef
        if     (present(tagPatternMax)) tagPatMax = tagPatternMax
        if(.not.present(tagPatternMax)) tagPatMax = 20
        if     (present(wordsMax))      wdMax     = wordsMax
        if(.not.present(wordsMax))      wdMax     = 1000

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Putting all the data in a vector (for accessibility)
        fileID      = 11;
        stat        = 0;

        !write(*,*) "FIRST LOOP"
        !write(*,*) "path = ", path
        open (unit = fileID , file = path, action = 'read')
        !write(*,*) "AFTER OPEN"
        contentSize = 0
        !wdMax       = 20 !For Tests
        allocate (contentVector(wdMax)) !Limitation about the number of words the file can have
        contentVector = "notUsed"
        do while(stat.eq.0)
            contentSize = contentSize + 1
            !write(*,*) "contentVector = ", contentVector
            read(fileID, fmt=*,IOSTAT = stat) contentVector(1:contentSize)
            rewind(fileID);
            !write(*,*) "stat", stat
            !if(contentSize > 15) stop("ERROR READING FILE (truncated at 200 words)")
        end do

        !write(*,*) "AFTER"
        contentSize = contentSize - 1
        deallocate (contentVector)
        allocate (contentVector(contentSize))
        rewind(fileID);
        read(fileID, fmt=*,IOSTAT = stat) contentVector(:)
        !write(*,*) "contentVector = ", contentVector
        close(fileID)

        conStart     = 1
        commentCount = 0

        !write(*,*) "contentVector = ", contentVector
        do i = 1, size(contentVector) !Treating comments between (!comment!)
            if(contentVector(i)(1:1) == comment) then
                commentCount = commentCount + 1
                if (mod(commentCount,2) == 0) then
                    contentVector(conStart:i-1) = ""
                    contentVector(i) = contentVector(i)(2:)
                else
                    conStart = i
                end if
            end if
            j = len(trim(contentVector(i)))
            if(contentVector(i)(j:j) == comment .and. j>1) then
                commentCount = commentCount + 1
                if (mod(commentCount,2) == 0) then
                    contentVector(conStart:i) = ""
                else
                    conStart = i
                end if
            end if
        end do
        if (mod(commentCount,2) == 1) contentVector(conStart:) = "" !Treating final comment
        !call DispCarvalhol(contentVector, "contentVector BEFORE")
        !write(*,*) " size(contentVector) = ", size(contentVector)

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Starting data treatment
        blankTotal  = 0;
        dataTotal   = 0;
        foundedTags = 'notUsed';
        posTaken    = .FALSE.;
        stat        = 0;
        tagTotal    = 0;
        unitTags    = 0;

        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>First loop: finding tags and counting data
        do i = 1, size(contentVector)
            !            nLines = nLines+1;
            if (contentVector(i)(1:1).eq.tagID) then !If a tag i founded
                posTaken = .FALSE.
                do j = 1, size(foundedTags) !Founding tag position
                    if(contentVector(i) == foundedTags(j)) then
                        posTaken = .TRUE.
                        exit
                    endif
                enddo
                if(posTaken .eqv. .FALSE.) then
                    tagTotal = tagTotal +1
                    if(tagTotal > size(foundedTags)) stop("Too many tags, max is 500")
                    foundedTags(tagTotal) = contentVector(i)
                    posTaken = .TRUE.
                endif
                if (contentVector(i)(2:2) == tagID) unitTags = unitTags + 1
                !write(*,*) "tag = ", tagTotal;
            else if (.not.((contentVector(i)(1:1) == comment) .or. (contentVector(i)(1:1) == empty))) then
                dataTotal = dataTotal+1
                !write(*,*) "data = ", dataTotal;
            else
                blankTotal = blankTotal+1
                !write(*,*) "Comment or blank";
            endif
        enddo

        !call DispCarvalhol(contentVector, "contentVector")
        !write(*,*) "dataTotal            = ", dataTotal;
        !write(*,*) "tagTotal             = ", tagTotal;
        !write(*,*) "blank/comment Total  = ", blankTotal;

        if(dataTotal == 0 .or. tagTotal == 0) then
            write(*,*) "ERROR! In set_DataTable, no tags and/or data bad conditioned "
            stop "ERROR! In set_DataTable, no tags and/or data bad conditioned ";
            write(*,*) "dataTotal = ", dataTotal
            write(*,*) " tagTotal = ", tagTotal
            call dispCarvalhol(foundedTags,"foundedTags")
            call dispCarvalhol(contentVector, "contentVector")
        else if((tagTotal-unitTags) == 0) then
            !            write(*,*) "Only unit data"
            ratioTagData = 1;
        else if(mod((dataTotal-unitTags),(tagTotal-unitTags)) == 0) then
            ratioTagData = (dataTotal-unitTags)/(tagTotal-unitTags);
            !write(*,*) "ratioTagData         = ", ratioTagData;
        else
            write(*,*) "ERROR! In set_DataTable, Tags and data dimensions don't match "
            stop "ERROR! In set_DataTable, Tags and data dimensions don't match ";
            write(*,*) "dataTotal = ", dataTotal
            write(*,*) " tagTotal = ", tagTotal
            call dispCarvalhol(foundedTags,"foundedTags")
            call dispCarvalhol(contentVector, "contentVector")
        endif

        !write(*,*) "SECOND LOOP"
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>Second loop: filling data table

        !Initializing
        allocate(dataTable(ratioTagData+1, tagTotal));
        allocate(tagPattern(tagTotal));
        allocate(lastLine(size(dataTable,2)));


        dataRow        = 0;
        dataColumn     = 0;
        dataCount      = 0;
        dataPassed     = .FALSE.;
        dataTable      = "    ---";
        labelPassed    = .FALSE.;
        lastLine       = 1; !To take into account the headers
        tagCount       = 0;
        tagPattern     = 0;


        do i = 1, size(contentVector)
            if (contentVector(i)(1:1).eq.tagID) then !If a tag i founded
                labelPassed = .true.;  !says we have put a label (so now we can put data)

                if(dataPassed.eqv..TRUE.) then !If it's a new pattern we reset the variables
                    tagPattern = 0;       !Vector indicating all tags columns position
                    tagCount   = 0;       !Integer indicating this tag position in this pattern
                    dataPassed = .FALSE.; !Says that it'll be a new data packet
                    dataCount  = 0;       !Integer saying in which position of the tag pattern the data should go
                end if

                tagCount = tagCount +1;

                do j = 1, tagTotal !Founding tag position
                    if(contentVector(i).eq.foundedTags(j)) then
                        tagPattern(tagCount) = j;
                        exit
                    endif
                enddo

            else if ((.not.((contentVector(i)(1:1) == comment) .or. (contentVector(i)(1:1) == empty)))&
                .and.labelPassed.eqv..TRUE.) then

                dataPassed = .TRUE.;
                dataCount  = dataCount + 1;
                dataColumn = tagPattern(modCyclic(dataCount,tagCount));
                lastLine(dataColumn) = lastLine(dataColumn) +1;
                if(lastLine(dataColumn) > size(dataTable, 1)) write(*,*) "ERROR! Data is not equally distributed over the columns"
                if(lastLine(dataColumn) > size(dataTable, 1)) stop "ERROR! Data is not equally distributed over the columns"
                dataTable(lastLine(dataColumn), dataColumn) = contentVector(i)

                !write(*,*) "modCyclic(dataCount,tagCount)= ", modCyclic(dataCount,tagCount);
                !write(*,*) "tagPattern(modCyclic(dataCount,tagCount)) = ", tagPattern(modCyclic(dataCount,tagCount));
                !write(*,*) "Data line = ", dataLine;
                !write(*,*) "dataRow = ", lastLine(dataColumn);
            end if
        enddo

        do i = 1, tagTotal !Putting the headers
            if(foundedTags(i)(2:2) == tagID) then
                dataTable(1,i) = foundedTags(i)(3:); !Jumping the "$$"
            else
                dataTable(1,i) = foundedTags(i)(2:); !Jumping the "$"
            end if
        end do

        if(allocated(tagPattern))    deallocate(tagPattern)
        if(allocated(lastLine))      deallocate(lastLine)
        if(allocated(contentVector)) deallocate(contentVector)

    end subroutine set_DataTable

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_DataTable_DbleVec(dataTable, tagName, dataDestination)
        implicit none
        !INPUT
        character(len=*), dimension(:,:), intent(in)  :: dataTable
        character(len=*),                 intent(in)  :: tagName
        !OUTPUT
        double precision, dimension(:),   intent(out), allocatable :: dataDestination
        !LOCAL VARIABLES
        integer :: i, stat
        logical :: tagFounded

        allocate(dataDestination(size(dataTable,1)-1))

        tagFounded = .FALSE.
        stat       = -1

        do i = 1, size(dataTable,2)
            if (trim(dataTable(1,i)) == tagName) then
                read(dataTable(2:,i), fmt=*, IOSTAT = stat) dataDestination
                tagFounded = .TRUE.
                exit
            end if
        end do

        if(stat /= 0) then
            write(*,*) "ERROR! in -read_DataTable_DbleVec- read failed (check types and Tag name)"
            write(*,*) "tagName = ", tagName
            stop
        end if

        if(.not.tagFounded) then
            write(*,*) "ERROR! in -read_DataTable- TAG not founded"
            write(*,*) "tagName = ", tagName
            stop
        end if

    end subroutine read_DataTable_DbleVec

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_DataTable_IntVec(dataTable, tagName, dataDestination)
        implicit none
        !INPUT
        character(len=*), dimension(:,:), intent(in)  :: dataTable
        character(len=*),                 intent(in)  :: tagName
        !OUTPUT
        integer,          dimension(:),   intent(out), allocatable :: dataDestination
        !LOCAL VARIABLES
        integer :: i, stat
        logical :: tagFounded

        allocate(dataDestination(size(dataTable,1)-1))

        tagFounded = .FALSE.
        stat       = -1

        do i = 1, size(dataTable,2)
            if (trim(dataTable(1,i)) == tagName) then
                read(dataTable(2:,i), fmt=*, IOSTAT = stat) dataDestination
                tagFounded = .TRUE.
                exit
            end if
        end do

        if(stat /= 0) then
            write(*,*) "ERROR! in -read_DataTable_DbleVec- read failed (check types and Tag name)"
            write(*,*) "tagName = ", tagName
            stop
        end if

        if(.not.tagFounded) then
            write(*,*) "ERROR! in -read_DataTable- TAG not founded"
            write(*,*) "tagName = ", tagName
            stop
        end if

    end subroutine read_DataTable_IntVec

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_DataTable_CharVec(dataTable, tagName, dataDestination)
        implicit none
        !INPUT
        character(len=*), dimension(:,:), intent(in)  :: dataTable
        character(len=*),                 intent(in)  :: tagName
        !OUTPUT
        character(len=*), dimension(:),   intent(out), allocatable :: dataDestination
        !LOCAL VARIABLES
        integer :: i, stat
        logical :: tagFounded

        allocate(dataDestination(size(dataTable,1)-1))

        tagFounded = .FALSE.
        stat       = -1

        do i = 1, size(dataTable,2)
            if (trim(dataTable(1,i)) == tagName) then
                read(dataTable(2:,i), fmt=*, IOSTAT = stat) dataDestination
                tagFounded = .TRUE.
                exit
            end if
        end do

        if(stat /= 0) then
            write(*,*) "ERROR! in -read_DataTable_DbleVec- read failed (check types and Tag name)"
            write(*,*) "tagName = ", tagName
            stop
        end if

        if(.not.tagFounded) then
            write(*,*) "ERROR! in -read_DataTable- TAG not founded"
            write(*,*) "tagName = ", tagName
            stop
        end if

    end subroutine read_DataTable_CharVec

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_DataTable_DbleScal(dataTable, tagName, dataDestination)
        implicit none
        !INPUT
        character(len=*), dimension(:,:), intent(in)  :: dataTable
        character(len=*),                 intent(in)  :: tagName
        !OUTPUT
        double precision,                intent(out)  :: dataDestination
        !LOCAL VARIABLES
        integer :: i, stat
        logical :: tagFounded

        tagFounded = .FALSE.
        stat       = -1

        do i = 1, size(dataTable,2)
            if (trim(dataTable(1,i)) == tagName) then
                read(dataTable(2,i), fmt=*, IOSTAT = stat) dataDestination
                tagFounded = .TRUE.
                exit
            end if
        end do

        if(stat /= 0) then
            write(*,*) "ERROR! in -read_DataTable_DbleVec- read failed (check types and Tag name)"
            write(*,*) "tagName = ", tagName
            stop
        end if

        if(.not.tagFounded) then
            write(*,*) "ERROR! in -read_DataTable- TAG not founded"
            write(*,*) "tagName = ", tagName
            stop
        end if

    end subroutine read_DataTable_DbleScal

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_DataTable_LogicalScal(dataTable, tagName, dataDestination)
        implicit none
        !INPUT
        character(len=*), dimension(:,:), intent(in)  :: dataTable
        character(len=*),                 intent(in)  :: tagName
        !OUTPUT
        logical, intent(out)  :: dataDestination
        !LOCAL VARIABLES
        integer :: i, stat
        logical :: tagFounded

        tagFounded = .FALSE.
        stat       = -1

        do i = 1, size(dataTable,2)
            if (trim(dataTable(1,i)) == tagName) then
                read(dataTable(2,i), fmt=*, IOSTAT = stat) dataDestination
                tagFounded = .TRUE.
                exit
            end if
        end do

        if(stat /= 0) then
            write(*,*) "ERROR! in -read_DataTable_DbleVec- read failed (check types and Tag name)"
            write(*,*) "tagName = ", tagName
            stop
        end if

        if(.not.tagFounded) then
            write(*,*) "ERROR! in -read_DataTable- TAG not founded"
            write(*,*) "tagName = ", tagName
            stop
        end if

    end subroutine read_DataTable_LogicalScal

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_DataTable_IntScal(dataTable, tagName, dataDestination)
        implicit none
        !INPUT
        character(len=*), dimension(:,:), intent(in)  :: dataTable
        character(len=*),                 intent(in)  :: tagName
        !OUTPUT
        integer,                         intent(out)  :: dataDestination
        !LOCAL VARIABLES
        integer :: i, stat
        logical :: tagFounded

        tagFounded = .FALSE.
        stat       = -1

        do i = 1, size(dataTable,2)
            if (trim(dataTable(1,i)) == tagName) then
                read(dataTable(2,i), fmt=*, IOSTAT = stat) dataDestination
                tagFounded = .TRUE.
                exit
            end if
        end do

        if(stat /= 0) then
            write(*,*) "ERROR! in -read_DataTable_IntScal- read failed (check types and Tag name)"
            write(*,*) "tagName = ", tagName
            stop
        end if

        if(.not.tagFounded) then
            write(*,*) "ERROR! in -read_DataTable_IntScal- TAG not founded"
            write(*,*) "tagName = ", tagName
            stop
        end if

    end subroutine read_DataTable_IntScal

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine read_DataTable_CharScal(dataTable, tagName, dataDestination)
        implicit none
        !INPUT
        character(len=*), dimension(:,:), intent(in)  :: dataTable
        character(len=*),                 intent(in)  :: tagName
        !OUTPUT
        character(len=*),                 intent(out)  :: dataDestination
        !LOCAL VARIABLES
        integer :: i, stat
        logical :: tagFounded

        tagFounded = .FALSE.
        stat       = -1

        do i = 1, size(dataTable,2)
            if (trim(dataTable(1,i)) == tagName) then
                !read(dataTable(2,i), fmt=*, IOSTAT = stat) dataDestination
                dataDestination = adjustL(dataTable(2,i))
                !write(*,*) "dataDestination = ", dataDestination
                tagFounded = .TRUE.
                exit
            end if
        end do

!        if(stat /= 0) then
!            write(*,*) "ERROR! in -read_DataTable_CharScal- read failed (check types and Tag name)"
!            write(*,*) "tagName = ", tagName
!            stop
!        end if

        if(.not.tagFounded) then
            write(*,*) "ERROR! in -read_DataTable_CharScal- TAG not founded"
            write(*,*) "tagName = ", tagName
            stop
        end if

    end subroutine read_DataTable_CharScal

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function modCyclic(dividend,divisor) result(cyclicRest)
        implicit none
        integer :: dividend, divisor, cyclicRest

        if (mod(dividend, divisor).eq.0) then
            cyclicRest = divisor;
        else
            cyclicRest = mod(dividend, divisor);
        endif
    end function modCyclic

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function getLine (fid, comment_Tag) result(nextLine)

        integer,          intent(in) :: fid
        character(len=1), intent(in) :: comment_Tag
        character(len=buf_RF) :: nextLine
        integer :: lineCount = 200, i, stat

        do i = 1, lineCount
            read(fid, fmt="(A)",IOSTAT = stat) nextLine
            nextLine = adjustL(nextLine)
            !print*, "nextLine = ", trim(nextLine)
            if(stat /= 0) then
                nextLine = "eof_gl"
                exit
            else if(nextLine(1:1) /= comment_Tag) then
                !write(*,*) "nextLine = ", trim(nextLine)
                exit
            end if
        end do

    end function getLine

end module readFile_RF
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
