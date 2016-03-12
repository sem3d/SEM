module type_MESH

    use mpi
    use charFunctions
    use type_inputRF

    implicit none

    type :: MESH
        !MPI VARIABLES
        integer :: comm = -1
        integer :: rang = -1
        integer :: nb_procs = -1
        integer :: topComm = -1
        integer :: log_ID = -1
        logical :: validProc

        !MESH VARIABLES
            !nDim independent
        integer :: meshMod, method
        integer :: nDim = -1
        integer(kind=8) :: xNTotal
        integer :: nOvlpPoints, nOvlpMax
        !integer :: xNTotal
        !integer :: xNTotal = -1, xNInit = -1, xNEnd = -1, xNGlob = -1
        !logical :: independent
            !nDim dependent
        integer         , dimension(:), allocatable :: xNStep, procPerDim;
        integer         , dimension(:), allocatable :: neigh, coords, op_neigh;
        double precision, dimension(:), allocatable :: xMaxGlob, xMinGlob;
        double precision, dimension(:), allocatable :: xStep, corrL
        double precision, dimension(:), allocatable :: procStart, procExtent
        integer         , dimension(:), allocatable :: pointsPerCorrL;
        double precision, dimension(:), allocatable :: xMaxInt, xMinInt; !Non-overlapping area
        double precision, dimension(:), allocatable :: xMaxExt, xMinExt; !Interface with neighbours proc
        double precision, dimension(:), allocatable :: xMaxBound, xMinBound; !Bounding box for generation in this proc
        double precision, dimension(:,:), allocatable :: xMaxNeigh, xMinNeigh; !Overlapping area in each direction
        double precision, dimension(:,:), allocatable :: xOrNeigh; !Origin for Shape Functions
        integer         , dimension(2,1) :: indexLocal
        integer         , dimension(:,:), allocatable :: indexNeigh, neighShift
        logical, dimension(:), allocatable :: considerNeighbour
        double precision, dimension(:), allocatable :: overlap !Size of the overlap (in corrL)
        integer, dimension(:)  , allocatable :: origin
        integer, dimension(:)  , allocatable :: mappingFromShift
        logical :: init = .false.

        !coords - Coordinates vector on the comunications Carte
        !procPerDim - number of processor in each direction
        !neigh - rang of the neighbour in a given direction (if there is not it should be -1)
        !op_neigh - rang of the neighbour in the oposite direction (if there is not it should be -1)
        !procStart - minimum coordinate that this processor will write on the HDF5 (defines origin)
        !xMaxBound, xMinBound - Bounding box for generation in this proc. It limitates the values of xPoints

    end type MESH

    contains
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_MESH(MESH_a, IPT, comm, rang, newNbProcs)
            implicit none
            !INPUT
            type(IPT_RF), intent(in)  :: IPT
            integer, intent(in) :: comm, rang
            integer, intent(in), optional :: newNbProcs
            !OUTPUT
            type(MESH) :: MESH_a
            !LOCAL
            integer :: nDim

            nDim = IPT%nDim_mesh

            allocate(MESH_a%xStep(nDim))
            allocate(MESH_a%xNStep(nDim))
            allocate(MESH_a%xMaxGlob(nDim))
            allocate(MESH_a%xMinGlob(nDim))
            allocate(MESH_a%xMaxInt(nDim))
            allocate(MESH_a%xMinInt(nDim))
            allocate(MESH_a%procPerDim(nDim))
            allocate(MESH_a%coords(nDim))
            allocate(MESH_a%xMaxExt(nDim))
            allocate(MESH_a%xMinExt(nDim))
            allocate(MESH_a%xMaxBound(nDim))
            allocate(MESH_a%xMinBound(nDim))
            allocate(MESH_a%neigh((3**nDim)-1))
            allocate(MESH_a%op_neigh((3**nDim)-1))
            allocate(MESH_a%xMaxNeigh(nDim,(3**nDim)-1))
            allocate(MESH_a%xMinNeigh(nDim,(3**nDim)-1))
            allocate(MESH_a%xOrNeigh(nDim,(3**nDim)-1))
            allocate(MESH_a%indexNeigh(2,(3**nDim)-1))
            allocate(MESH_a%neighShift(nDim,(3**nDim)-1))
            allocate(MESH_a%considerNeighbour((3**nDim)-1))
            allocate(MESH_a%mappingFromShift((3**nDim)-1))
            allocate(MESH_a%overLap(nDim))
            allocate(MESH_a%pointsPerCorrL(nDim))
            allocate(MESH_a%corrL(nDim))
            allocate(MESH_a%origin(nDim))
            allocate(MESH_a%procExtent(nDim))
            allocate(MESH_a%procStart(nDim))

            MESH_a%log_ID   = IPT%log_ID
            MESH_a%nDim     = IPT%nDim_mesh
            MESH_a%comm     = comm
            if(present(newNbProcs)) then
                MESH_a%nb_procs = newNbProcs
                MESH_a%procPerDim = -1
            else
                MESH_a%nb_procs = product(IPT%nFields)
                MESH_a%procPerDim = IPT%nFields
            end if
            MESH_a%rang     = rang
            MESH_a%meshMod  = IPT%meshMod
            MESH_a%xMaxGlob = IPT%xMaxGlob
            MESH_a%xMinGlob = IPT%xMinGlob
            MESH_a%corrL    = IPT%corrL
            MESH_a%method   = IPT%method
            MESH_a%pointsPerCorrL = IPT%pointsPerCorrL
            MESH_a%overlap(:)     = IPT%overlap
            !MESH_a%independent    = IPT%independent

            MESH_a%xMaxBound     = -1
            MESH_a%xMinBound     = -1
            MESH_a%xMaxExt     = -1
            MESH_a%xMinExt     = -1
            MESH_a%xStep    = -1
            MESH_a%xNStep   = -1
            MESH_a%xMaxInt  = -1
            MESH_a%xMinInt  = -1
            MESH_a%coords(:)  = -1
            MESH_a%xMaxNeigh(:,:) = 0
            MESH_a%xMinNeigh(:,:) = 0
            MESH_a%xOrNeigh(:,:) = 0
            MESH_a%indexNeigh(:,:) = -1
            MESH_a%mappingFromShift(:) = -1
            MESH_a%neigh(:) = -2 !-1 is already the default when the proc is in the topology border
            MESH_a%neighShift(:,:) = 0
            MESH_a%considerNeighbour = .true.
            MESH_a%init = .true.
            MESH_a%validProc = .true.

        end subroutine init_MESH

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_MESH(MESH_a, name, fmt, forLog_in, unit_in)
            !INPUT
            type(MESH), intent(in) :: MESH_a
            character(len=*), intent(in), optional :: name
            character(len = 20), intent(in), optional :: fmt
            integer, intent(in), optional :: unit_in
            logical, intent(in), optional :: forLog_in
            !LOCAL
            character(len = 200) :: dblFmt, matDblFmt, intFmt, matIntFmt
            integer :: unit
            logical :: forLog
            logical :: active

            active = .true.
            unit = 6 !Screen
            if(present(unit_in)) unit = unit_in

            forLog = .false.
            if(present(forLog_in)) then
#ifdef MAKELOG
                if(forLog_in) unit = MESH_a%log_ID
                forLog = forLog_in
#else
                if(forLog_in) write(*,*) "WARNING!!! Inside show_MESH, forLog_in = .true. but MAKELOG was not defined"
                active = .false.
#endif
            end if

            if(active) then

                if(unit <= 0) then
                    write(*,*) "ERROR!!! Inside show_MESH unit = ", unit
                    stop("")
                end if

                dblFmt = "T20,F15.5"
                intFmt = "T20,I15"
                matDblFmt = string_join_many("T20,",trim(numb2String(MESH_a%nDim)),"F15.5")
                matIntFmt = string_join_many("T20,",trim(numb2String(MESH_a%nDim)),"I15")
                if(present(fmt)) dblFmt = fmt

                write(unit,*) "MESH------------------------------------------------------------"
                if(present(name)) write(unit,*) "|  ", name

                if(MESH_a%init) then
                    write(unit,*) "|  init     = ", MESH_a%init
                    write(unit,*) "|"
                    write(unit,*) "|  MPI---"
                    write(unit,*) "|  |rang     = ", MESH_a%rang
                    write(unit,*) "|  |comm     = ", MESH_a%comm
                    write(unit,*) "|  |nb_procs = ", MESH_a%nb_procs
                    write(unit,*) "|  |topComm  = ", MESH_a%topComm
                    write(unit,*) "|"
                    write(unit,*) "|  Input---"
                    write(unit,*) "|  |nDim     = ", MESH_a%nDim
                    !write(unit,*) "|  |independent = ", MESH_a%independent
                    write(unit,*) "|  |meshMod  = ", MESH_a%meshMod
                    write(unit,"(A,("//dblFmt//"))") " |  |xMinGlob = ", MESH_a%xMinGlob
                    write(unit,"(A,("//dblFmt//"))") " |  |xMaxGlob = ", MESH_a%xMaxGlob
                    write(unit,"(A,("//dblFmt//"))") " |  |xStep    = ", MESH_a%xStep
                    write(unit,"(A,F15.5)")" |  |overLap  = ", MESH_a%overLap
                    write(unit,*) "|"
                    write(unit,*) "|  |Process--"
                    write(unit,*) "|  |xNStep     = ", MESH_a%xNStep
                    !write(unit,"(A,("//dblFmt//"))") " |  |xMinExt       = ", MESH_a%xMinExt
                    !write(unit,"(A,("//dblFmt//"))") " |  |xMaxExt       = ", MESH_a%xMaxExt
                    write(unit,*) "|  |xNTotal    = ", MESH_a%xNTotal
                    write(unit,*) "|  |procPerDim = ", MESH_a%procPerDim
                    write(unit,*) "|  |coords     = ", MESH_a%coords
                    write(unit,"(A,("//dblFmt//"))") " |  |xMinInt    = ", MESH_a%xMinInt
                    write(unit,"(A,("//dblFmt//"))") " |  |xMaxInt    = ", MESH_a%xMaxInt
                    write(unit,"(A,("//dblFmt//"))") " |  |xMinExt  = ", MESH_a%xMinExt
                    write(unit,"(A,("//dblFmt//"))") " |  |xMaxExt  = ", MESH_a%xMaxExt
                    call show_MESHneigh(MESH_a, onlyExisting = .false., forLog = forLog, unit_in = unit)

                else
                    write(unit,*) "|  init     = ", MESH_a%init
                    write(unit,*) "|  MESH has not been initialized----"
                end if
                write(unit,*) "|---------------------------------------------------------------"
                write(unit,*) ""

            end if

        end subroutine show_MESH

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_MESHneigh(MESH_a, name, onlyExisting, forLog, unit_in)
            implicit none
            !INPUT
            type(MESH), intent(in) :: MESH_a
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
                if(forLog) unit = MESH_a%log_ID
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

                nDim   = trim(numb2String(MESH_a%nDim))
                space  = trim(numb2String(MESH_a%nDim*numb_space + 1))
                nbChar = trim(numb2String(numb_space))

                if(present(name)) write(unit,*) "|  ", name

                fmtChar = "(A15, A"//space//", A"//space//", A"//space//")"
                fmtNum = "(I6, A3, I5, A1, "//nDim//"I"//nbChar//", A1, "//nDim//"F"//nbChar//".2, A1, "//nDim//"F"//nbChar//".2)"

                if(MESH_a%init) then

                    write(unit,fmtChar) "Neighbour/ op|","Shift|", "xMin|", "xMax|"

                    do i = 1, size(MESH_a%neigh)
                        if(onlyExisting .and. MESH_a%neigh(i)<0) cycle
                        write(unit,fmtNum) MESH_a%neigh(i), " / ", MESH_a%op_neigh(i), &
                                           "|", MESH_a%neighShift(:,i), "|", &
                                           MESH_a%xMinNeigh(:,i), "|", MESH_a%xMaxNeigh(:,i)
                    end do

                else
                    write(unit,*) "|  init     = ", MESH_a%init
                    write(unit,*) "|  MESH has not been initialized----"
                end if

                write(unit,*) ""

            end if

        end subroutine show_MESHneigh

        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        function findLabel(neighShift) result(label)

            implicit none
            !INPUT
            integer, dimension(:), intent(in) :: neighShift

            !OUTPUT
            integer :: label

            !LOCAL
            integer :: i, nDim

            nDim = size(neighShift)

            label = 0

            do i = 1, nDim
                label = label + neighShift(i)*3**(nDim - i)
            end do

            if(label > 0) label = label - 1 !Because no neighbour has shift (0,0,0)

            label = label + 1 +(3**(nDim)-1)/2

            !call wLog(" label = ")
            !call wLog(label)

            !In 3D labels go from 1 (-1,-1,-1) to 26 (1, 1, 1)


        end function findLabel

        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        !-----------------------------------------------------------------------------------------------
        function findShiftRang(MESH_a, neighShift) result(rang)

            implicit none
            !INPUT
            !INPUT
            type(MESH), intent(in) :: MESH_a
            integer, dimension(:), intent(in) :: neighShift
            !OUTPUT
            integer :: rang

            rang = MESH_a%mappingFromShift(findLabel(neighShift))
            call wLog(" rang = ")
            call wLog(rang)

        end function findShiftRang

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine finalize_MESH(MESH_a)
            type(MESH) :: MESH_a

            !if (allocated(MESH_a%xMaxExt))       deallocate(MESH_a%xMaxExt)
            !if (allocated(MESH_a%xMinExt))       deallocate(MESH_a%xMinExt)
            if (allocated(MESH_a%xStep))      deallocate(MESH_a%xStep)
            if (allocated(MESH_a%xNStep))     deallocate(MESH_a%xNStep)
            if (allocated(MESH_a%xMaxGlob))   deallocate(MESH_a%xMaxGlob)
            if (allocated(MESH_a%xMinGlob))   deallocate(MESH_a%xMinGlob)
            if (allocated(MESH_a%procPerDim)) deallocate(MESH_a%procPerDim)
            if (allocated(MESH_a%coords))     deallocate(MESH_a%coords)
            if (allocated(MESH_a%neigh))      deallocate(MESH_a%neigh)
            if (allocated(MESH_a%xOrNeigh))      deallocate(MESH_a%xOrNeigh)
            if (allocated(MESH_a%xMaxNeigh))  deallocate(MESH_a%xMaxNeigh)
            if (allocated(MESH_a%xMinNeigh))  deallocate(MESH_a%xMinNeigh)
            if (allocated(MESH_a%xMaxInt))    deallocate(MESH_a%xMaxInt)
            if (allocated(MESH_a%xMinInt))    deallocate(MESH_a%xMinInt)
            if (allocated(MESH_a%indexNeigh)) deallocate(MESH_a%indexNeigh)
            if (allocated(MESH_a%xMaxExt))  deallocate(MESH_a%xMaxExt)
            if (allocated(MESH_a%xMinExt))  deallocate(MESH_a%xMinExt)
            if (allocated(MESH_a%xMaxBound))  deallocate(MESH_a%xMaxBound)
            if (allocated(MESH_a%xMinBound))  deallocate(MESH_a%xMinBound)
            if (allocated(MESH_a%neighShift)) deallocate(MESH_a%neighShift)
            if (allocated(MESH_a%overlap))    deallocate(MESH_a%overlap)
            if (allocated(MESH_a%pointsPerCorrL)) deallocate(MESH_a%pointsPerCorrL)
            if (allocated(MESH_a%considerNeighbour)) deallocate(MESH_a%considerNeighbour)
            if (allocated(MESH_a%corrL))  deallocate(MESH_a%corrL)
            if (allocated(MESH_a%origin)) deallocate(MESH_a%origin)
            if (allocated(MESH_a%procExtent)) deallocate(MESH_a%procExtent)
            if (allocated(MESH_a%procStart)) deallocate(MESH_a%procStart)
            if (allocated(MESH_a%mappingFromShift)) deallocate(MESH_a%mappingFromShift)
            if (allocated(MESH_a%op_neigh)) deallocate(MESH_a%op_neigh)


            MESH_a%init = .false.

        end subroutine finalize_MESH

end module type_MESH
