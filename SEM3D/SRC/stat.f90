module stat
    use constants, only : fpp, DM_MAX, OUT_DOM_NAMES
    use sem_c_bindings, only : sem_mkdir
    implicit none
    ! TYPE DE STATISTIQUES
    integer, parameter :: STAT_COUNT = 14
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
    integer, parameter :: STAT_START  =11 ! Initialisation/startup time
    integer, parameter :: STAT_FSOL_DG=12 ! Forces      solide
    integer, parameter :: STAT_ITER=13    !

    character(len=5), dimension(0:STAT_COUNT-1) :: stat_labels = (/ &
        "GIVE ", "TAKE ", "WAIT ", "FSOL ", "FFLU ", "PSOL ", "PFLU ", "FEXT ", "FULL ", "TSTEP","IO   ", "INIT ", "FSDG ", "ITER " /)

    integer*8, private :: clockRate, maxPeriod

    integer*8, private :: startFullTick, stopFullTick, deltaTick
    real(fpp), private :: fullTime

    real(fpp), dimension(0:STAT_COUNT-1), private :: statTimes
    integer*8, dimension(0:STAT_COUNT-1), private :: statStart, statStop
    integer, dimension(1:DM_MAX) :: elemCounts
    integer, dimension(1:DM_MAX) :: ngllCounts
    logical :: log_traces
    integer :: ntraces, itrace, rank
    integer*8, allocatable, dimension(:,:) :: traces
contains


    subroutine add_start_trace(step)
        implicit none
        integer, intent(in) :: step

        traces(itrace,1) = step
        traces(itrace,2) = statStart(step)-startFullTick
        itrace = itrace + 1
    end subroutine add_start_trace

    subroutine add_stop_trace(step)
        implicit none
        integer, intent(in) :: step

        traces(itrace,1) = step+1000
        traces(itrace,2) = statStop(step)-startFullTick
        itrace = itrace + 1
        if (itrace>ntraces) then
            write(124) traces(1:itrace-1,1)
            write(125) traces(1:itrace-1,2)
            itrace = 1
        end if
    end subroutine add_stop_trace

    subroutine stat_init(comm, rg, nprocs, dotraces)
        use mpi
        implicit none
        integer, intent(in) :: comm, rg, nprocs
        logical, intent(in) :: dotraces
        character(Len=1000) :: fname
        integer :: ierr

        call system_clock(COUNT_RATE=clockRate)
        call system_clock(COUNT_MAX=maxPeriod)
        call system_clock(count=startFullTick)

        log_traces = dotraces
        statStart(STAT_FULL)=startFullTick
        if (log_traces) call add_start_trace(STAT_FULL)

        statTimes(:) = 0d0
        statStart(:) = 0d0
        statStop(:) = 0d0

        rank = rg

        if (log_traces) then
            ntraces = 2000
            itrace = 1
            if (rg==0) then
                ierr = sem_mkdir("timings")
            end if
            call MPI_Barrier(comm, ierr)
            write(fname,"(A,I6.6,A)") "timings/type.", rank,".bin"
            open(124, file=trim(adjustl(fname)),form="unformatted", access="stream")
            write(fname,"(A,I6.6,A)") "timings/time.", rank,".bin"
            open(125, file=trim(adjustl(fname)),form="unformatted", access="stream")

            allocate(traces(ntraces+4,2))
        end if

    end subroutine stat_init


    subroutine stat_finalize(Tdomain)
        use mpi
        use sdomain, only : domain_nelems, domain
        implicit none
        type(domain), intent(in) :: Tdomain
        integer :: ierr, rank, sz, r, i
        integer :: status(MPI_STATUS_SIZE)

        call system_clock(count=stopFullTick)
        statStop(STAT_FULL)=stopFullTick
        if (log_traces) then
            ntraces=1 ! force output
            call add_stop_trace(STAT_FULL)
            close(124)
            close(125)
        end if

        deltaTick = stopFullTick-startFullTick
        if (deltaTick < 0) deltaTick = deltaTick+maxPeriod
        fullTime = real(deltaTick,fpp)/real(clockRate,fpp)

        do i=1, DM_MAX
            call domain_nelems(Tdomain, i, elemCounts(i), ngllCounts(i))
        end do
        statTimes(STAT_FULL) = fullTime
        call MPI_Comm_Rank (MPI_COMM_WORLD, rank, ierr)
        if (rank .gt. 0) then
            call MPI_SEND(elemCounts,DM_MAX,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
            call MPI_SEND(ngllCounts,DM_MAX,MPI_INTEGER,0,0,MPI_COMM_WORLD,ierr)
            call MPI_SEND(statTimes,STAT_COUNT,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierr)
        else
            open (123, file="stat.log", status="replace", action="write")
            call MPI_Comm_Size (MPI_COMM_WORLD,   sz, ierr)
            do r = 0, sz-1
                if (r .gt. 0) then
                    call MPI_RECV(elemCounts,DM_MAX,MPI_INTEGER,r,0,MPI_COMM_WORLD,status,ierr)
                    call MPI_RECV(ngllCounts,DM_MAX,MPI_INTEGER,r,0,MPI_COMM_WORLD,status,ierr)
                    call MPI_RECV(statTimes,STAT_COUNT,MPI_DOUBLE_PRECISION,r,0,MPI_COMM_WORLD,status,ierr)
                end if
                do i=1,DM_MAX
                    write(123, '(a12,I6,I8,I8)') OUT_DOM_NAMES(i), r, elemCounts(i), ngllCounts(i)
                end do
                do i=0,STAT_COUNT-1
                    write(123, '(a8,I6,f12.5)') stat_labels(i), r, statTimes(i)
                end do
            end do
            close (123)
        end if
    end subroutine stat_finalize

    subroutine stat_starttick(step)
        implicit none
        integer, intent(in) :: step

        call system_clock(count=statStart(step))

        if (log_traces) call add_start_trace(step)

    end subroutine stat_starttick

    subroutine stat_stoptick(step)
        implicit none
        integer, intent(in) :: step

        real(fpp) :: time

        call system_clock(count=statStop(step))
        if (log_traces) call add_stop_trace(step)
        deltaTick = statStop(step)-statStart(step)
        if (deltaTick < 0) deltaTick = deltaTick+maxPeriod
        time = real(deltaTick,fpp)/real(clockRate,fpp)
        statTimes(step) = statTimes(step) + time
    end subroutine stat_stoptick

end module stat

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
