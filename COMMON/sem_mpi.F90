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
    use mpi
#else
    use fake_mpi
#endif
#ifdef SINGLEPRECISION
    parameter (MPI_REAL_FPP=MPI_REAL)
#else
    parameter (MPI_REAL_FPP=MPI_DOUBLE_PRECISION)
#endif
end module sem_mpi
