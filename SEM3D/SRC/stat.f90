module stat

    use constants
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
        integer :: ierr, rank, sz, r
        integer :: status(MPI_STATUS_SIZE)
        real(fpp) :: calctime

        call system_clock(count=stopFullTick)
        deltaTick = stopFullTick-startFullTick
        if (deltaTick < 0) deltaTick = deltaTick+maxPeriod
        fullTime = real(deltaTick)/clockRate

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

                write (123,'(a,i4,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a)') "TIMING - stat comm : rank ", r, &
                ", comm time ", statTimes(STAT_GIVE) + statTimes(STAT_WAIT) + statTimes(STAT_TAKE), &
                " sec [give ", statTimes(STAT_GIVE), " sec, wait ", &
                statTimes(STAT_WAIT), " sec, take ", &
                statTimes(STAT_TAKE), " sec]"
                calctime = statTimes(STAT_FSOL)+statTimes(STAT_FFLU)+statTimes(STAT_PSOL)+statTimes(STAT_PFLU)
                write (123,'(a,i4,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a)') &
                    "TIMING - stat calc : rank ", r,         &
                    ", calc time ",  calctime + statTimes(STAT_FEXT), &
                    " sec [fsol:", statTimes(STAT_FSOL), &
                    " sec fflu:", statTimes(STAT_FFLU), &
                    " sec psol:", statTimes(STAT_PSOL), &
                    " sec pflu:", statTimes(STAT_PFLU), &
                    " sec, fext ", statTimes(STAT_FEXT), " sec]"

                write (123,'(a,i4,a,f10.3,a)') "TIMING - stat full : rank ", r, &
                    ", full time ", statTimes(STAT_FULL), " sec"
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
