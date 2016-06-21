module stat
    use constants, only : fpp
    ! TYPE DE STATISTIQUES
    integer, parameter :: STAT_COUNT = 12
    integer, parameter :: STAT_GIVE  = 0 !
    integer, parameter :: STAT_TAKE  = 1 !
    integer, parameter :: STAT_WAIT  = 2 !
    integer, parameter :: STAT_FSOL  = 3 ! Forces      solide
    integer, parameter :: STAT_FFLU  = 4 ! Forces      fluide
    integer, parameter :: STAT_PSOL  = 5 ! Forces  PML solide
    integer, parameter :: STAT_PFLU  = 6 ! Forces  PML fluide
    integer, parameter :: STAT_FEXT  = 7 !
    integer, parameter :: STAT_FULL  = 8 !
    integer, parameter :: STAT_TSTEP = 9 ! Timestep synchronisation
    integer, parameter :: STAT_IO    =10 ! IO
    integer, parameter :: STAT_START =11 ! Initialisation/startup time

    character(len=5), dimension(0:STAT_COUNT-1) :: stat_labels = (/ &
        "GIVE ", "TAKE ", "WAIT ", "FSOL ", "FFLU ", "PSOL ", "PFLU ", "FEXT ", "FULL ", "TSTEP","IO   ", "INIT " /)

    integer, private :: clockRate, maxPeriod

    integer, private :: startFullTick, stopFullTick
    real(fpp)   , private :: fullTime

    integer     , private :: startTick, stopTick, deltaTick
    real(fpp), dimension(0:STAT_COUNT-1), private :: statTimes
    contains

    subroutine stat_init()
        implicit none

        call system_clock(COUNT_RATE=clockRate)
        call system_clock(COUNT_MAX=maxPeriod)

        call system_clock(count=startFullTick)

        statTimes(:) = 0d0
    end subroutine stat_init

    subroutine stat_finalize()
        use mpi
        implicit none
        integer :: ierr, rank, sz, r, i
        integer :: status(MPI_STATUS_SIZE)
        real(fpp) :: calctime

        call system_clock(count=stopFullTick)
        deltaTick = stopFullTick-startFullTick
        if (deltaTick < 0) deltaTick = deltaTick+maxPeriod
        fullTime = real(deltaTick)/clockRate

        statTimes(STAT_FULL) = fullTime
        call MPI_Comm_Rank (MPI_COMM_WORLD, rank, ierr)
        if (rank .gt. 0) then
            call MPI_SEND(statTimes,STAT_COUNT,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierr)
        else
            open (123, file="stat.log", status="replace", action="write")
            call MPI_Comm_Size (MPI_COMM_WORLD,   sz, ierr)
            do r = 0, sz-1
                if (r .gt. 0) then
                    call MPI_RECV(statTimes,STAT_COUNT,MPI_DOUBLE_PRECISION,r,0,MPI_COMM_WORLD,status,ierr)
                end if

                do i=0,STAT_COUNT-1
                    write(123, '(a6,I4,f10.3)') stat_labels(i), r, statTimes(i)
                end do
            end do
            close (123)
        end if
    end subroutine stat_finalize

    subroutine stat_starttick()
        implicit none

        call system_clock(count=startTick)
    end subroutine stat_starttick

    subroutine stat_stoptick(step)
        implicit none
        integer :: step

        real(fpp) :: time

        call system_clock(count=stopTick)
        deltaTick = stopTick-startTick
        if (deltaTick < 0) deltaTick = deltaTick+maxPeriod
        time = real(deltaTick)/clockRate
        statTimes(step) = statTimes(step) + time
    end subroutine stat_stoptick

end module stat
