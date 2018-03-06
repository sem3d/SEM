!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module scomm
    use constants, only : fpp
    implicit none
    interface comm_give_data
        module procedure comm_give_data_1, comm_give_data_2, comm_give_data_3
    end interface comm_give_data

    interface comm_take_data
        module procedure comm_take_data_1, comm_take_data_2, comm_take_data_3
    end interface comm_take_data

contains
#ifdef SINGLEPRECISION
#define MPI_REAL_FPP MPI_FLOAT
#else
#define MPI_REAL_FPP MPI_DOUBLE_PRECISION
#endif
    subroutine exchange_sem_var(Tdomain, tag, vector)
        use sdomain
        use mpi
        use stat, only : stat_starttick, stat_stoptick, STAT_WAIT


        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: tag
        type(comm_vector), intent(inout) :: vector

        integer, dimension(MPI_STATUS_SIZE,vector%ncomm) :: statuses
        integer :: dest, src, ierr, i

        !- now we can exchange (communication global arrays)
        vector%send_reqs = MPI_REQUEST_NULL
        vector%recv_reqs = MPI_REQUEST_NULL

        call stat_starttick(STAT_WAIT)
        do i = 0,vector%ncomm-1
            dest = vector%Data(i)%dest
            src = vector%Data(i)%src
            if (vector%Data(i)%ndata>0) then
                call MPI_Isend(vector%Data(i)%Give, vector%Data(i)%ndata, &
                    MPI_REAL_FPP, dest, tag, Tdomain%communicateur, &
                    vector%send_reqs(i), ierr)

                call MPI_Irecv(vector%Data(i)%Take, vector%Data(i)%ndata, &
                    MPI_REAL_FPP, dest, tag, Tdomain%communicateur, &
                    vector%recv_reqs(i), ierr)
            end if
        enddo

        call MPI_Waitall(vector%ncomm, vector%recv_reqs, statuses, ierr)
        call MPI_Waitall(vector%ncomm, vector%send_reqs, statuses, ierr)
        call stat_stoptick(STAT_WAIT)

    end subroutine exchange_sem_var

    subroutine comm_give_data_1(give, igive, data, pos)
        real(fpp), dimension(0:), intent(inout) :: give
        real(fpp), dimension(0:), intent(in)  :: data
        integer, dimension(0:), intent(in) :: igive
        integer, intent(inout) :: pos
        !
        integer :: i
        do i=0,size(igive)-1
            give(pos) = data(igive(i))
            pos = pos + 1
        end do
    end subroutine comm_give_data_1

    subroutine comm_give_data_2(give, igive, data, pos)
        real(fpp), dimension(0:), intent(inout) :: give
        real(fpp), dimension(0:,0:), intent(in)  :: data
        integer, dimension(0:), intent(in) :: igive
        integer, intent(inout) :: pos
        !
        integer :: i,j
        do i=0,size(igive)-1
            do j=0,size(data,2)-1
                give(pos) = data(igive(i),j)
                pos = pos + 1
            end do
        end do
    end subroutine comm_give_data_2

    subroutine comm_give_data_3(give, igive, data, pos)
        real(fpp), dimension(0:), intent(inout) :: give
        real(fpp), dimension(0:,0:,0:), intent(in)  :: data
        integer, dimension(0:), intent(in) :: igive
        integer, intent(inout) :: pos
        !
        integer :: i,j,k
        do i=0,size(igive)-1
            do j=0,size(data,2)-1
                do k=0,size(data,3)-1
                    give(pos) = data(igive(i),j,k)
                    pos = pos + 1
                end do
            end do
        end do
    end subroutine comm_give_data_3

    subroutine comm_take_data_1(take, itake, data, pos)
        real(fpp), dimension(0:), intent(in) :: take
        integer, dimension(0:), intent(in) :: itake
        real(fpp), dimension(0:), intent(inout)  :: data
        integer, intent(inout) :: pos
        !
        integer :: i, idx
        do i=0,size(itake)-1
            idx = itake(i)
            data(idx) = data(idx) + take(pos)
            pos = pos + 1
        end do
    end subroutine comm_take_data_1

    subroutine comm_take_data_2(take, itake, data, pos)
        real(fpp), dimension(0:), intent(in) :: take
        integer, dimension(0:), intent(in) :: itake
        real(fpp), dimension(0:,0:), intent(inout)  :: data
        integer, intent(inout) :: pos
        !
        integer :: i, j, idx
        do i=0,size(itake)-1
            idx = itake(i)
            do j=0,size(data,2)-1
                data(idx,j) = data(idx,j) + take(pos)
                pos = pos + 1
            end do
        end do
    end subroutine comm_take_data_2

    subroutine comm_take_data_3(take, itake, data, pos)
        real(fpp), dimension(0:), intent(in) :: take
        integer, dimension(0:), intent(in) :: itake
        real(fpp), dimension(0:,0:,0:), intent(inout)  :: data
        integer, intent(inout) :: pos
        !
        integer :: i, j, k, idx
        do i=0,size(itake)-1
            idx = itake(i)
            do j=0,size(data,2)-1
                do k=0,size(data,3)-1
                    data(idx,j,k) = data(idx,j,k) + take(pos)
                    pos = pos + 1
                end do
            end do
        end do
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
