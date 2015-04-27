!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module scomm

contains

    subroutine exchange_sem_var(Tdomain, tag, vector)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: tag
        type(comm_vector), intent(inout) :: vector

        integer, dimension(MPI_STATUS_SIZE,vector%ncomm) :: statuses
        integer :: dest, src, ierr, i

        !- now we can exchange (communication global arrays)
        vector%send_reqs = MPI_REQUEST_NULL
        vector%recv_reqs = MPI_REQUEST_NULL

        do i = 0,vector%ncomm-1
            dest = vector%Data(i)%dest
            src = vector%Data(i)%src
            call MPI_Isend(vector%Data(i)%Give, vector%Data(i)%ndata, &
                           MPI_DOUBLE_PRECISION, dest, tag, Tdomain%communicateur, &
                           vector%send_reqs(i), ierr)

            call MPI_Irecv(vector%Data(i)%Take, vector%Data(i)%ndata, &
                           MPI_DOUBLE_PRECISION, dest, tag, Tdomain%communicateur, &
                           vector%recv_reqs(i), ierr)
        enddo

        call MPI_Waitall(vector%ncomm, vector%recv_reqs, statuses, ierr)
        call MPI_Waitall(vector%ncomm, vector%send_reqs, statuses, ierr)

    end subroutine exchange_sem_var

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
