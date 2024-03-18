!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! stub mpi module that either use the real mpi module or a fake one
!!
!<

module sem_mpi
#ifdef __MPI
    use mpi_f08
#ifdef SINGLEPRECISION
    type(MPI_datatype), parameter :: MPI_REAL_FPP=MPI_REAL
#else
    type(MPI_datatype), parameter :: MPI_REAL_FPP=MPI_DOUBLE_PRECISION
#endif
#else
    use fake_mpi
    integer, parameter :: MPI_REAL_FPP=3
#endif

contains

    subroutine acc_init()
#if defined( _OPENACC) && defined(__MPI)
        USE openacc
        character(len=16) :: local_rank_env
        integer :: status, local_rank
        ! OpenACC initialisation
        !$acc init
        ! Use the environment variable set by Slurm, as MPI_Comm_rank cannot be used here
        ! since this routine is used BEFORE the initialisation of MPI
        call get_environment_variable(name="SLURM_LOCALID", value=local_rank_env, status=status)
        if (status /= 0) then
            call get_environment_variable(name="OMPI_COMM_WORLD_LOCAL_RANK", value=local_rank_env, status=status)
        end if
        if (status == 0) then
            read(local_rank_env, *) local_rank
            ! Define the GPU to use via OpenACC
            call acc_set_device_num(local_rank, acc_get_device_type())
        else
            print *, "Error: impossible to determine the local rank of the process"
            stop 1
        end if
#endif
    end subroutine acc_init
end module sem_mpi
