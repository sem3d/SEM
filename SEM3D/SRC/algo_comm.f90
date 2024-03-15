!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module scomm
    use constants, only : fpp
#if OPENACC
    use openacc
#endif
    implicit none
    interface comm_give_data
        module procedure comm_give_data_1, comm_give_data_2, comm_give_data_3
    end interface comm_give_data

    interface comm_take_data
        module procedure comm_take_data_1, comm_take_data_2, comm_take_data_3
    end interface comm_take_data
    logical :: startup_init_done

contains
    subroutine exchange_proc(Tdomain, tag, vector, i, give, take)
        use sdomain
        use sem_mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: tag, i
        type(comm_vector), intent(inout) :: vector
        real(fpp), intent(inout), dimension(0:vector%Data(i)%ndata-1) :: give
        real(fpp), intent(inout), dimension(0:vector%Data(i)%ndata-1) :: take
        !
#ifdef __MPI
        ! Won't work with stub mpi module and nvhpc/openacc because
        ! somehow the compiler changes the datatype of give/take
        integer :: dest, ierr
        dest = vector%Data(i)%dest
        !$acc host_data use_device(give) if (startup_init_done)
        call MPI_Isend(give, vector%Data(i)%ndata, &
            MPI_REAL_FPP, dest, tag, Tdomain%communicateur, &
            vector%send_reqs(i), ierr)
        !$acc end host_data
        !$acc host_data use_device(take) if (startup_init_done)
        call MPI_Irecv(take, vector%Data(i)%ndata, &
            MPI_REAL_FPP, dest, tag, Tdomain%communicateur, &
            vector%recv_reqs(i), ierr)
        !$acc end host_data
#endif
    end subroutine exchange_proc
    subroutine exchange_sem_var(Tdomain, tag, vector)
        use sdomain
        use sem_mpi
        use stat, only : stat_starttick, stat_stoptick, STAT_WAIT


        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: tag
        type(comm_vector), intent(inout) :: vector

        type(MPI_Status), dimension(vector%ncomm) :: statuses
        integer :: dest, src, ierr, i

        !- now we can exchange (communication global arrays)
        vector%send_reqs = MPI_REQUEST_NULL
        vector%recv_reqs = MPI_REQUEST_NULL

        call stat_starttick(STAT_WAIT)
        !$acc wait(1)   if (startup_init_done)
        do i = 0,vector%ncomm-1
            if (vector%Data(i)%ndata>0) then
                call exchange_proc(Tdomain, tag, vector, i, vector%Data(i)%Give, vector%Data(i)%Take)
            end if
        enddo

        call MPI_Waitall(vector%ncomm, vector%recv_reqs, statuses, ierr)
        call MPI_Waitall(vector%ncomm, vector%send_reqs, statuses, ierr)
        call stat_stoptick(STAT_WAIT)

        if (startup_init_done .and. .false.) then
            open(111,file="comm.nok."//strrank(Tdomain%rank),status="unknown",form="formatted",position="append")
            write(111,*) "----------------------------------------------------"
            write(111,*) "ITERX:", Tdomain%TimeD%rtime
            do i = 0,vector%ncomm-1
                write(111,*) "COMM:", i, ":", vector%Data(i)%src, " -> ", vector%Data(i)%dest
                if (vector%Data(i)%ndata>0) then
                    write(111,"(*(6(F10.5,3X),/))") vector%Data(0)%Give
                end if
                write(111,*)
            end do
            close(111)
        end if
    end subroutine exchange_sem_var

    subroutine comm_give_data_1(give, igive, data, pos)
        real(fpp), dimension(0:), intent(inout) :: give
        real(fpp), dimension(0:), intent(in)  :: data
        integer, dimension(0:), intent(in) :: igive
        integer, intent(inout) :: pos
        !
        integer :: i, p0
        integer :: ni
        ni = size(igive)
        p0 = pos
        !$acc parallel loop async(1) if (startup_init_done)
        do i=0,ni-1
            give(p0+i) = data(igive(i))
        end do
        pos = p0+ni
    end subroutine comm_give_data_1

    subroutine comm_give_data_2(give, igive, data, pos)
        real(fpp), dimension(0:), intent(inout) :: give
        real(fpp), dimension(0:,0:), intent(in)  :: data
        integer, dimension(0:), intent(in) :: igive
        integer, intent(inout) :: pos
        !
        integer :: i,j,p0
        integer :: ni, nj
        ni = size(igive)
        nj = size(data,2)
        p0 = pos
        !$acc parallel loop collapse(2) async(1)  if (startup_init_done)
        do i=0,ni-1
            do j=0,nj-1
                give(p0+j+nj*i) = data(igive(i),j)
            end do
        end do
        pos = p0 + ni*nj
    end subroutine comm_give_data_2

    subroutine comm_give_data_3(give, igive, data, pos)
        real(fpp), dimension(0:), intent(inout) :: give
        real(fpp), dimension(0:,0:,0:), intent(in)  :: data
        integer, dimension(0:), intent(in) :: igive
        integer, intent(inout) :: pos
        !
        integer :: i,j,k,p0
        integer :: ni, nj, nk
        ni = size(igive)
        nj = size(data,2)
        nk = size(data,3)
        p0 = pos
        !$acc parallel loop collapse(3) async(1)  if (startup_init_done)
        do i=0,ni-1
            do j=0,nj-1
                do k=0,nk-1
                    give(p0 + k + nk*(j + nj*i)) = data(igive(i),j,k)
                end do
            end do
        end do
        pos = pos+ni*nj*nk
    end subroutine comm_give_data_3

    subroutine comm_take_data_1(take, itake, data, pos)
        real(fpp), dimension(0:), intent(in) :: take
        integer, dimension(0:), intent(in) :: itake
        real(fpp), dimension(0:), intent(inout)  :: data
        integer, intent(inout) :: pos
        !
        integer :: i, idx, p0
        integer :: ni
        ni = size(itake)
        p0 = pos
        !$acc parallel loop async(1)  if (startup_init_done)
        do i=0,ni-1
            idx = itake(i)
            data(idx) = data(idx) + take(p0+i)
        end do
        pos = p0+ni
    end subroutine comm_take_data_1

    subroutine comm_take_data_2(take, itake, data, pos)
        real(fpp), dimension(0:), intent(in) :: take
        integer, dimension(0:), intent(in) :: itake
        real(fpp), dimension(0:,0:), intent(inout)  :: data
        integer, intent(inout) :: pos
        !
        integer :: i, j, idx, p0
        integer :: ni, nj
        ni = size(itake)
        nj = size(data,2)
        p0 = pos
        !$acc parallel loop collapse(2) async(1)  if (startup_init_done)
        do i=0,ni-1
            do j=0,nj-1
                idx = itake(i)
                data(idx,j) = data(idx,j) + take(p0+j+nj*i)
            end do
        end do
        pos = p0+ni*nj
    end subroutine comm_take_data_2

    subroutine comm_take_data_3(take, itake, data, pos)
        real(fpp), dimension(0:), intent(in) :: take
        integer, dimension(0:), intent(in) :: itake
        real(fpp), dimension(0:,0:,0:), intent(inout)  :: data
        integer, intent(inout) :: pos
        !
        integer :: i, j, k, idx, p0
        integer :: ni, nj, nk
        ni = size(itake)
        nj = size(data,2)
        nk = size(data,3)
        p0 = pos
        !$acc parallel loop collapse(3) async(1)  if (startup_init_done)
        do i=0,ni-1
            do j=0,nj-1
                do k=0,nk-1
                    idx = itake(i)
                    data(idx,j,k) = data(idx,j,k) + take(p0 + k + nk*(j + nj*i))
                end do
            end do
        end do
        pos = p0 +ni*nj*nk
    end subroutine comm_take_data_3

    subroutine exchange_sf_normals(Tdomain)
        use sdomain
        type(domain), intent(inout) :: Tdomain
        !
        integer :: n, k, p
        if(Tdomain%Comm_solflu%ncomm <= 0) return

        do n = 0,Tdomain%Comm_solflu%ncomm-1
            ! Domain SOLID
            ! Domain SOLID PML
            ! ignored
            k = 0

            ! Domain FLUID
            do p=0,2
                call comm_give_data(Tdomain%Comm_solflu%Data(n)%Give, &
                    Tdomain%Comm_SolFlu%Data(n)%IGiveF, Tdomain%SF%SF_BtN(p,:), k)
            end do
            ! Domain FLUID PML
            if (Tdomain%Comm_SolFlu%Data(n)%nflupml>0) then
                do p=0,2
                    call comm_give_data(Tdomain%Comm_SolFlu%Data(n)%Give, &
                        Tdomain%Comm_SolFlu%Data(n)%IGiveFPML, Tdomain%SF%SFpml_BtN(p,:), k)
                end do
            end if
            Tdomain%Comm_SolFlu%Data(n)%nsend = k
        end do

        ! Exchange
        call exchange_sem_var(Tdomain, 108, Tdomain%Comm_SolFlu)

        ! Take
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            ! Domain SOLID
            ! Domain SOLID PML
            ! .. ignored
            ! Domain FLUID
            k = 0
            do p=0,2
                call comm_take_data(Tdomain%Comm_SolFlu%Data(n)%Take, &
                    Tdomain%Comm_SolFlu%Data(n)%IGiveF, Tdomain%SF%SF_BtN(p,:), k)
            end do
            ! Domain FLUID PML
                if (Tdomain%Comm_SolFlu%Data(n)%nflupml>0) then
                    do p=0,2
                        call comm_take_data(Tdomain%Comm_SolFlu%Data(n)%Take, &
                            Tdomain%Comm_SolFlu%Data(n)%IGiveFPML, Tdomain%SF%SFpml_BtN(p,:), k)
                    end do
            end if
        end do

    end subroutine exchange_sf_normals

end module scomm

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
