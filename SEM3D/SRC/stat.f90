module stat

    integer*8, private :: clockRate, maxPeriod

    integer*8, private :: startFullTick, stopFullTick
    real*8   , private :: fullTime

    integer*8, private :: startTick, stopTick, deltaTick
    real*8   , private :: waitTime
    real*8   , private :: giveTime
    real*8   , private :: takeTime
    real*8   , private :: fintTime
    real*8   , private :: fextTime

    contains

    subroutine stat_init()
        implicit none

        call system_clock(COUNT_RATE=clockRate)
        call system_clock(COUNT_MAX=maxPeriod)

        fullTime=0.
        call system_clock(count=startFullTick)

        waitTime=0.
        giveTime=0.
        takeTime=0.
        fintTime=0.
        fextTime=0.
    end subroutine stat_init

    subroutine stat_finalize()
        use mpi
        implicit none
        integer :: ierr, rank
        character(len=5) :: rankName

        call MPI_Comm_Rank (MPI_COMM_WORLD, rank, ierr)
        write (rankName, '(I5)') rank
        open (123, file="stat."//trim(adjustl(rankName))//".log", status="replace", action="write")

        write (123,'(a,i4,a,f10.3,a,f10.3,a,f10.3,a,f10.3,a)') "TIMING - stat comm : rank ", rank, &
        ", comm time ", giveTime+waitTime+takeTime,                                                  &
        " sec [give ", giveTime, " sec, wait ", waitTime, " sec, take ", takeTime, " sec]"

        write (123,'(a,i4,a,f10.3,a,f10.3,a,f10.3,a)') "TIMING - stat calc : rank ", rank,         &
        ", calc time ", fintTime+fextTime, " sec [fint ", fintTime, " sec, fext ", fextTime, " sec]"

        call system_clock(count=stopFullTick)
        deltaTick = stopFullTick-startFullTick
        if (deltaTick < 0) deltaTick = deltaTick+maxPeriod
        fullTime = real(deltaTick)/clockRate
        write (123,'(a,i4,a,f10.3,a)') "TIMING - stat full : rank ", rank, ", full time ", fullTime, " sec"

        close (123)
    end subroutine stat_finalize

    subroutine stat_starttick()
        implicit none

        call system_clock(count=startTick)
    end subroutine stat_starttick

    subroutine stat_stoptick(step)
        implicit none
        character(len=4) :: step

        real*8 :: time

        call system_clock(count=stopTick)
        deltaTick = stopTick-startTick
        if (deltaTick < 0) deltaTick = deltaTick+maxPeriod
        time = real(deltaTick)/clockRate
        if (step == 'give') giveTime = giveTime+time
        if (step == 'take') takeTime = takeTime+time
        if (step == 'wait') waitTime = waitTime+time
        if (step == 'fint') fintTime = fintTime+time
        if (step == 'fext') fextTime = fextTime+time
        if (step == 'full') fullTime = fullTime+time
    end subroutine stat_stoptick

end module stat
