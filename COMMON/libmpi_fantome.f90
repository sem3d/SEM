!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!     Fichier contenant une emulation des fonctions     !!!!!!!!!!
!!!!!!!!!!     MPI pour permettre un fonctionnement sequentiel   !!!!!!!!!!
!!!!!!!!!!     liste des fonctions :                             !!!!!!!!!!
!!!!!!!!!!     MPI_INIT                                          !!!!!!!!!!
!!!!!!!!!!     MPI_FINALIZE                                      !!!!!!!!!!
!!!!!!!!!!     MPI_COMM_RANK                                     !!!!!!!!!!
!!!!!!!!!!     MPI_COMM_SIZE                                     !!!!!!!!!!
!!!!!!!!!!     MPI_GROUP_SIZE                                    !!!!!!!!!!
!!!!!!!!!!     MPI_GROUP_RANK                                    !!!!!!!!!!
!!!!!!!!!!     MPI_COMM_GROUP                                    !!!!!!!!!!
!!!!!!!!!!     MPI_SEND                                          !!!!!!!!!!
!!!!!!!!!!     MPI_RECV                                          !!!!!!!!!!
!!!!!!!!!!     MPI_BARRIER                                       !!!!!!!!!!
!!!!!!!!!!     MPI_ALLREDUCE                                     !!!!!!!!!!
!!!!!!!!!!     MPI_GATHER                                        !!!!!!!!!!
!!!!!!!!!!     MPI_GATHERV                                       !!!!!!!!!!
!!!!!!!!!!     MPI_SCATTER                                       !!!!!!!!!!
!!!!!!!!!!     MPI_SCATTERV                                      !!!!!!!!!!
!!!!!!!!!!     MPI_ALLGATHER                                     !!!!!!!!!!
!!!!!!!!!!     MPI_SENDRECV                                      !!!!!!!!!!
!!!!!!!!!!     MPI_ISEND                                         !!!!!!!!!!
!!!!!!!!!!     MPI_ALLGATHERV                                    !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mpi
    implicit none

    INTERFACE MPI_Allreduce
       MODULE PROCEDURE MPI_Allreduce_rv, MPI_Allreduce_r1, MPI_Allreduce_i1, MPI_Allreduce_iv
    END INTERFACE
    INTERFACE MPI_Send
       MODULE PROCEDURE MPI_Send_rv, MPI_Send_r1, MPI_Send_iv, MPI_Send_im, MPI_Send_i1, MPI_Send_rm, MPI_Send_rm3
    END INTERFACE
    INTERFACE MPI_Recv
       MODULE PROCEDURE MPI_Recv_rv, MPI_Recv_r1, MPI_Recv_iv, MPI_Recv_i1, MPI_Recv_rm, MPI_Recv_rm3
    END INTERFACE
    INTERFACE MPI_Reduce
       MODULE PROCEDURE MPI_Reduce_rv, MPI_Reduce_r1, MPI_Reduce_i1, MPI_Reduce_iv
    END INTERFACE
    INTERFACE MPI_Isend
       MODULE PROCEDURE MPI_Isend_rv, MPI_Isend_r1
    END INTERFACE
    INTERFACE MPI_Gatherv
       MODULE PROCEDURE MPI_Gatherv_rv, MPI_Gatherv_r1, MPI_Gatherv_iv, MPI_Gatherv_i1
    END INTERFACE
    INTERFACE MPI_Scatterv
       MODULE PROCEDURE MPI_Scatterv_rv, MPI_Scatterv_r1, MPI_Scatterv_iv, MPI_Scatterv_i1, MPI_Scatterv_im, MPI_Scatterv_im3
    END INTERFACE
    INTERFACE MPI_Bcast
       MODULE PROCEDURE MPI_Bcast_rv, MPI_Bcast_r1, MPI_Bcast_iv, MPI_Bcast_i1
    END INTERFACE
    INTERFACE MPI_Sendrecv
       MODULE PROCEDURE MPI_Sendrecv_rv, MPI_Sendrecv_r1, MPI_Sendrecv_iv, MPI_Sendrecv_i1, MPI_Sendrecv_rm
    END INTERFACE
    INTERFACE MPI_Allgather
       MODULE PROCEDURE MPI_Allgather_rv, MPI_Allgather_r1, MPI_Allgather_iv, MPI_Allgather_i1
    END INTERFACE

    include 'mpif.h'

contains


    subroutine MPI_INIT (ierr)
        implicit none
        integer, intent(inout) :: ierr

        ierr = 0

    end subroutine MPI_INIT

    subroutine MPI_FINALIZE (ierr)
        implicit none
        integer, intent(inout) :: ierr

        ierr = 0

    end subroutine MPI_FINALIZE

    subroutine MPI_COMM_RANK (comm, rank, ierr)
        implicit none
        integer, intent(in) :: comm
        integer, intent(out) :: rank
        integer, intent(inout) :: ierr

        rank = 0
        ierr = 0

    end subroutine MPI_COMM_RANK

    subroutine MPI_COMM_SIZE (comm, size, ierr)
        implicit none
        integer, intent(in) :: comm
        integer, intent(out):: size
        integer, intent(inout) :: ierr

        size = 1
        ierr = 0

    end subroutine MPI_COMM_SIZE

    subroutine MPI_GROUP_SIZE (group, size, ierr)
        implicit none
        integer, intent(in) :: group
        integer, intent(out):: size
        integer, intent(inout) :: ierr

        size = 1
        ierr = 0

    end subroutine MPI_GROUP_SIZE

    subroutine MPI_GROUP_RANK (group, rank, ierr)
        implicit none
        integer, intent(in) :: group
        integer, intent(out):: rank
        integer, intent(inout) :: ierr

        rank = 0
        ierr = 0

    end subroutine MPI_GROUP_RANK

    subroutine MPI_COMM_GROUP (comm, group, ierr)
        implicit none
        integer, intent(in) :: comm
        integer, intent(out):: group
        integer, intent(inout) :: ierr

        group = 0
        ierr = 0

    end subroutine MPI_COMM_GROUP


    !! ----------------------- MPI_Bcast ---------------------------
    subroutine MPI_BCAST_rv (buffer, count, type, root, comm, ierr)
        implicit none
        integer, intent(in) :: count, type, root, comm
        integer, intent(inout) :: ierr
        real, intent(inout) :: buffer(count)
        ierr = 0
    end subroutine MPI_BCAST_rv

    subroutine MPI_BCAST_r1 (buffer, count, type, root, comm, ierr)
        implicit none
        integer, intent(in) :: count, type, root, comm
        integer, intent(inout) :: ierr
        real, intent(inout) :: buffer
        ierr = 0
    end subroutine MPI_BCAST_r1

    subroutine MPI_BCAST_iv (buffer, count, type, root, comm, ierr)
        implicit none
        integer, intent(in) :: count, type, root, comm
        integer, intent(inout) :: ierr
        integer, intent(inout) :: buffer(count)
        ierr = 0
    end subroutine MPI_BCAST_iv

    subroutine MPI_BCAST_i1 (buffer, count, type, root, comm, ierr)
        implicit none
        integer, intent(in) :: count, type, root, comm
        integer, intent(inout) :: ierr
        integer, intent(inout) :: buffer
        ierr = 0
    end subroutine MPI_BCAST_i1

    !! ----------------------- MPI_Send ---------------------------
    subroutine MPI_SEND_r1 (val, long, type, dest, tag, comm, ierr)
        implicit none
        real, intent(in) :: val
        integer, intent(in) :: long, type, dest, tag, comm
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine MPI_SEND_r1

    subroutine MPI_SEND_rv (val, long, type, dest, tag, comm, ierr)
        implicit none
        real, intent(in) :: val(long)
        integer, intent(in) :: long, type, dest, tag, comm
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine MPI_SEND_rv

    subroutine MPI_SEND_rm (val, long, type, dest, tag, comm, ierr)
        implicit none
        real, intent(in), dimension(:,:) :: val
        integer, intent(in) :: long, type, dest, tag, comm
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine MPI_SEND_rm

    subroutine MPI_SEND_rm3 (val, long, type, dest, tag, comm, ierr)
        implicit none
        real, intent(in), dimension(:,:,:) :: val
        integer, intent(in) :: long, type, dest, tag, comm
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine MPI_SEND_rm3

    subroutine MPI_SEND_i1 (val, long, type, dest, tag, comm, ierr)
        implicit none
        integer, intent(in) :: val
        integer, intent(in) :: long, type, dest, tag, comm
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine MPI_SEND_i1

    subroutine MPI_SEND_iv (val, long, type, dest, tag, comm, ierr)
        implicit none
        integer, intent(in) :: val(long)
        integer, intent(in) :: long, type, dest, tag, comm
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine MPI_SEND_iv

    subroutine MPI_SEND_im (val, long, type, dest, tag, comm, ierr)
        implicit none
        integer, intent(in), dimension(:,:) :: val
        integer, intent(in) :: long, type, dest, tag, comm
        integer, intent(inout) :: ierr
        ierr = 0
    end subroutine MPI_SEND_im

    !! ----------------------- MPI_Recv ---------------------------
    subroutine MPI_RECV_r1 (val, long, type, src, tag, comm, status, ierr)
        implicit none
        real, intent(in) :: val
        integer, intent(in) :: long, type, src, tag, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_RECV_r1

    subroutine MPI_RECV_rv (val, long, type, src, tag, comm, status, ierr)
        implicit none
        real, intent(in) :: val(long)
        integer, intent(in) :: long, type, src, tag, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_RECV_rv

    subroutine MPI_RECV_rm (val, long, type, src, tag, comm, status, ierr)
        implicit none
        real, intent(in), dimension(:,:) :: val
        integer, intent(in) :: long, type, src, tag, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_RECV_rm

    subroutine MPI_RECV_rm3 (val, long, type, src, tag, comm, status, ierr)
        implicit none
        real, intent(in), dimension(:,:,:) :: val
        integer, intent(in) :: long, type, src, tag, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_RECV_rm3

    subroutine MPI_RECV_i1 (val, long, type, src, tag, comm, status, ierr)
        implicit none
        integer, intent(in) :: val
        integer, intent(in) :: long, type, src, tag, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_RECV_i1

    subroutine MPI_RECV_iv (val, long, type, src, tag, comm, status, ierr)
        implicit none
        integer, intent(in) :: val(long)
        integer, intent(in) :: long, type, src, tag, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_RECV_iv

    !! ----------------------- MPI_Barrier ---------------------------
    subroutine MPI_BARRIER (comm, ierr)
        implicit none
        integer, intent(in) :: comm
        integer, intent(inout) ::ierr

        ierr = 0

    end subroutine MPI_BARRIER


    !! ----------------------- MPI_Allreduce ---------------------------
    subroutine MPI_Allreduce_rv (sendbuf, recvbuf, long, type, op, comm, ierr)
        implicit none

        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: op
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr
        real, intent(inout) :: recvbuf(long)
        real, intent(in) :: sendbuf(long)

        if ( (op == MPI_MAX) .OR. (op == MPI_MIN) .OR. (op == MPI_SUM) ) then
            recvbuf = sendbuf
        else
            print*,'MPI_ALLREDUCE : type d operation non traite par la librairie mpi fantome',op
            stop
        endif
        ierr = 0
    end subroutine MPI_Allreduce_rv

    subroutine MPI_Allreduce_r1 (sendbuf, recvbuf, long, type, op, comm, ierr)
        implicit none

        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: op
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr
        real, intent(inout) :: recvbuf
        real, intent(in) :: sendbuf

        if ( (op == MPI_MAX) .OR. (op == MPI_MIN) .OR. (op == MPI_SUM) ) then
            recvbuf = sendbuf
        else
            print*,'MPI_ALLREDUCE : type d operation non traite par la librairie mpi fantome',op
            stop
        endif
        ierr = 0
    end subroutine MPI_Allreduce_r1


    subroutine MPI_Allreduce_i1 (sendbuf, recvbuf, long, type, op, comm, ierr)
        implicit none

        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: op
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr
        integer, intent(inout) :: recvbuf
        integer, intent(in) :: sendbuf

        if ( (op == MPI_MAX) .OR. (op == MPI_MIN) .OR. (op == MPI_SUM) ) then
            recvbuf = sendbuf
        else
            print*,'MPI_ALLREDUCE : type d operation non traite par la librairie mpi fantome',op
            stop
        endif
        ierr = 0
    end subroutine MPI_Allreduce_i1

    subroutine MPI_Allreduce_iv (sendbuf, recvbuf, long, type, op, comm, ierr)
        implicit none

        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: op
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr
        integer, intent(inout) :: recvbuf(long)
        integer, intent(in) :: sendbuf(long)

        if ( (op == MPI_MAX) .OR. (op == MPI_MIN) .OR. (op == MPI_SUM) ) then
            recvbuf = sendbuf
        else
            print*,'MPI_ALLREDUCE : type d operation non traite par la librairie mpi fantome',op
            stop
        endif
        ierr = 0
    end subroutine MPI_Allreduce_iv

    !! ----------------------- MPI_Reduce ---------------------------
    subroutine MPI_REDUCE_rv (sendbuf, recvbuf, long, type, op, dest, comm, ierr)
        implicit none

        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: op
        integer, intent(in) :: dest
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr

        real, intent(inout) :: recvbuf(long)
        real, intent(in) :: sendbuf(long)

        recvbuf = sendbuf
        ierr = 0

    end subroutine MPI_REDUCE_rv

    subroutine MPI_REDUCE_r1 (sendbuf, recvbuf, long, type, op, dest, comm, ierr)
        implicit none

        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: op
        integer, intent(in) :: dest
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr

        real, intent(inout) :: recvbuf
        real, intent(in) :: sendbuf

        recvbuf = sendbuf
        ierr = 0

    end subroutine MPI_REDUCE_r1

    subroutine MPI_REDUCE_iv (sendbuf, recvbuf, long, type, op, dest, comm, ierr)
        implicit none

        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: op
        integer, intent(in) :: dest
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr

        integer, intent(inout) :: recvbuf(long)
        integer, intent(in) :: sendbuf(long)

        recvbuf = sendbuf
        ierr = 0

    end subroutine MPI_REDUCE_iv

    subroutine MPI_REDUCE_i1 (sendbuf, recvbuf, long, type, op, dest, comm, ierr)
        implicit none

        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: op
        integer, intent(in) :: dest
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr

        integer, intent(inout) :: recvbuf
        integer, intent(in) :: sendbuf

        recvbuf = sendbuf
        ierr = 0

    end subroutine MPI_REDUCE_i1

    !! ----------------------- MPI_Gather ---------------------------
    subroutine MPI_GATHER(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr )
        implicit none

        integer, intent(in) :: sendtype, recvtype, sendcount, recvcount, root, comm
        integer, intent(inout) :: ierr
        real, intent(inout) :: recvbuf(recvcount)
        real, intent(in) :: sendbuf(sendcount)
        recvbuf = sendbuf
    end subroutine MPI_GATHER

    !! ----------------------- MPI_Gatherv ---------------------------
    subroutine MPI_GATHERV_rv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,displs, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, sendcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: recvcounts(*)
        integer, intent(in) :: displs(*)
        real, intent(inout) :: recvbuf(sendcount)
        real, intent(in) :: sendbuf(sendcount)
        recvbuf = sendbuf
    end subroutine MPI_GATHERV_rv

    subroutine MPI_GATHERV_r1(sendbuf, sendcount, sendtype, recvbuf, recvcounts,displs, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, sendcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: recvcounts(*)
        integer, intent(in) :: displs(*)
        real, intent(inout) :: recvbuf(sendcount)
        real, intent(in) :: sendbuf
        recvbuf = sendbuf
    end subroutine MPI_GATHERV_r1

    subroutine MPI_GATHERV_iv(sendbuf, sendcount, sendtype, recvbuf, recvcounts,displs, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, sendcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: recvcounts(*)
        integer, intent(in) :: displs(*)
        integer, intent(inout) :: recvbuf(sendcount)
        integer, intent(in) :: sendbuf(sendcount)
        recvbuf = sendbuf
    end subroutine MPI_GATHERV_iv

    subroutine MPI_GATHERV_i1(sendbuf, sendcount, sendtype, recvbuf, recvcounts,displs, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, sendcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: recvcounts(*)
        integer, intent(in) :: displs(*)
        integer, intent(inout) :: recvbuf(sendcount)
        integer, intent(in) :: sendbuf
        recvbuf = sendbuf
    end subroutine MPI_GATHERV_i1

    !! ----------------------- MPI_Scatter ---------------------------
    subroutine MPI_SCATTER(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr )
        implicit none

        integer, intent(in) :: sendtype
        integer, intent(in) :: recvtype
        integer, intent(in) :: sendcount
        integer, intent(in) :: recvcount
        integer, intent(in) :: root
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr

        real, intent(inout) :: recvbuf(recvcount)
        real, intent(in) :: sendbuf(sendcount)

        recvbuf = sendbuf


    end subroutine MPI_SCATTER

    !! ----------------------- MPI_Scatterv ---------------------------
    subroutine MPI_SCATTERV_rv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, recvcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: sendcounts(*)
        integer, intent(in) :: displs(*)
        real, intent(inout) :: recvbuf(recvcount)
        real, intent(in), dimension(:) :: sendbuf
        ierr = 0
        recvbuf=sendbuf
    end subroutine MPI_SCATTERV_rv

    subroutine MPI_SCATTERV_r1(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, recvcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: sendcounts(*)
        integer, intent(in) :: displs(*)
        real, intent(inout) :: recvbuf
        real, intent(in) :: sendbuf
        ierr = 0
        recvbuf=sendbuf
    end subroutine MPI_SCATTERV_r1

    subroutine MPI_SCATTERV_iv(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, recvcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: sendcounts(*)
        integer, intent(in) :: displs(*)
        integer, intent(inout) :: recvbuf(recvcount)
        integer, intent(in) :: sendbuf(recvcount)
        ierr = 0
        recvbuf = sendbuf
    end subroutine MPI_SCATTERV_iv

    subroutine MPI_SCATTERV_im(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, recvcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: sendcounts(*)
        integer, intent(in) :: displs(*)
        integer, intent(inout), dimension(:,:) :: recvbuf
        integer, intent(in), dimension(:,:) :: sendbuf
        ierr = 0
        recvbuf = sendbuf
    end subroutine MPI_SCATTERV_im

    subroutine MPI_SCATTERV_im3(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, recvcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: sendcounts(*)
        integer, intent(in) :: displs(*)
        integer, intent(inout), dimension(:,:,:) :: recvbuf
        integer, intent(in), dimension(:,:,:) :: sendbuf
        ierr = 0
        recvbuf = sendbuf
    end subroutine MPI_SCATTERV_im3


    subroutine MPI_SCATTERV_i1(sendbuf, sendcounts, displs, sendtype, recvbuf, recvcount, recvtype, root, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, recvcount, root, comm
        integer, intent(inout) :: ierr
        integer, intent(in) :: sendcounts(*)
        integer, intent(in) :: displs(*)
        integer, intent(inout) :: recvbuf
        integer, intent(in) :: sendbuf
        ierr = 0
        recvbuf = sendbuf
    end subroutine MPI_SCATTERV_i1

    !! ----------------------- MPI_Allgather ---------------------------
    subroutine MPI_ALLGATHER_rv(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, sendcount, recvcount, comm
        integer, intent(inout) :: ierr
        real, intent(inout) :: recvbuf(recvcount)
        real, intent(in) :: sendbuf(sendcount)
        recvbuf = sendbuf
    end subroutine MPI_ALLGATHER_rv

    subroutine MPI_ALLGATHER_r1(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, sendcount, recvcount, comm
        integer, intent(inout) :: ierr
        real, intent(inout) :: recvbuf(recvcount)
        real, intent(in) :: sendbuf
        ierr = 0
        recvbuf = sendbuf
    end subroutine MPI_ALLGATHER_r1

    subroutine MPI_ALLGATHER_iv(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, sendcount, recvcount, comm
        integer, intent(inout) :: ierr
        integer, intent(inout) :: recvbuf(recvcount)
        integer, intent(in) :: sendbuf(sendcount)
        ierr = 0
        recvbuf = sendbuf
    end subroutine MPI_ALLGATHER_iv

    subroutine MPI_ALLGATHER_i1(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm, ierr )
        implicit none
        integer, intent(in) :: sendtype, recvtype, sendcount, recvcount, comm
        integer, intent(inout) :: ierr
        integer, intent(inout) :: recvbuf
        integer, intent(in) :: sendbuf
        ierr = 0
        recvbuf = sendbuf
    end subroutine MPI_ALLGATHER_i1

    !! ----------------------- MPI_Sendrecv ---------------------------
    subroutine MPI_SENDRECV_r1(val1, long1, type1, dest1, tag1, val2, long2, type2, dest2, tag2, comm, status, ierr)
        implicit none
        real, intent(in) :: val1, val2
        integer, intent(in) :: long1, type1, dest1, tag1, long2, type2, dest2, tag2, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_SENDRECV_r1

    subroutine MPI_SENDRECV_rv(val1, long1, type1, dest1, tag1, val2, long2, type2, dest2, tag2, comm, status, ierr)
        implicit none
        real, intent(in), dimension(long1) :: val1
        real, intent(in), dimension(long2) :: val2
        integer, intent(in) :: long1, type1, dest1, tag1, long2, type2, dest2, tag2, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_SENDRECV_rv

    subroutine MPI_SENDRECV_rm(val1, long1, type1, dest1, tag1, val2, long2, type2, dest2, tag2, comm, status, ierr)
        implicit none
        real, intent(in), dimension(:,:) :: val1, val2
        integer, intent(in) :: long1, type1, dest1, tag1, long2, type2, dest2, tag2, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_SENDRECV_rm

    subroutine MPI_SENDRECV_i1(val1, long1, type1, dest1, tag1, val2, long2, type2, dest2, tag2, comm, status, ierr)
        implicit none
        integer, intent(in) :: val1, val2
        integer, intent(in) :: long1, type1, dest1, tag1, long2, type2, dest2, tag2, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_SENDRECV_i1

    subroutine MPI_SENDRECV_iv(val1, long1, type1, dest1, tag1, val2, long2, type2, dest2, tag2, comm, status, ierr)
        implicit none
        integer, intent(in), dimension(long1) :: val1
        integer, intent(in), dimension(long2) :: val2
        integer, intent(in) :: long1, type1, dest1, tag1, long2, type2, dest2, tag2, comm
        integer, intent(out), dimension(MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr
        ierr = 0
        status = 0
    end subroutine MPI_SENDRECV_iv

    !! ----------------------- MPI_Isend ---------------------------
    subroutine MPI_ISEND_r1 (val, long, type, dest, tag, comm, requete, ierr)
        implicit none
        real, intent(in) :: val
        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: dest
        integer, intent(in) :: tag
        integer, intent(in) :: comm
        integer, intent(out) :: requete
        integer, intent(inout) :: ierr

        ierr = 0
        requete = 0
    end subroutine MPI_ISEND_r1

    subroutine MPI_ISEND_rv (val, long, type, dest, tag, comm, requete, ierr)
        implicit none
        real, intent(in), dimension(long) :: val
        integer, intent(in) :: long
        integer, intent(in) :: type
        integer, intent(in) :: dest
        integer, intent(in) :: tag
        integer, intent(in) :: comm
        integer, intent(out) :: requete
        integer, intent(inout) :: ierr

        ierr = 0
        requete = 0
    end subroutine MPI_ISEND_rv

    subroutine MPI_ALLGATHERV(sendbuf, sendcount, sendtype, recvbuf, recvcounts,displs, recvtype, comm, ierr )
        implicit none

        integer, intent(in) :: sendtype
        integer, intent(in) :: recvtype
        integer, intent(in) :: sendcount
        integer, intent(in) :: comm
        integer, intent(inout) :: ierr

        integer, intent(in) :: recvcounts(*)
        integer, intent(in) :: displs(*)

        real, intent(inout), dimension(sendcount) :: recvbuf
        real, intent(in), dimension(sendcount) :: sendbuf

        recvbuf = sendbuf


    end subroutine MPI_ALLGATHERV

    subroutine MPI_WAIT(request, status, ierr)
        implicit none
        integer, intent(in) :: request
        integer, intent(out), dimension  (MPI_STATUS_SIZE) :: status
        integer, intent(inout) :: ierr

        status = 0
        ierr = 0

    end subroutine MPI_WAIT
end module mpi
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
