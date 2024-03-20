!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
program SEM3D
    use sem_mpi
    use semdatafiles
    use drive_sem

    character(Len=MAX_FILE_SIZE),parameter :: p_param = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_traces = "./traces"
    character(Len=MAX_FILE_SIZE),parameter :: p_results = "./res"
    character(Len=MAX_FILE_SIZE),parameter :: p_data = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_prot = "./prot"
    character(Len=MAX_FILE_SIZE),parameter :: p_mat = "./mat"
    character(Len=MAX_FILE_SIZE),parameter :: p_mirror = "./mirror"

    integer :: ierr, rg

    call init_sem_path(p_param, p_traces, p_results, p_data, p_prot, p_mat, p_mirror)
    call ACC_Init ()
    call MPI_Init (ierr)

    call sem(-1, MPI_COMM_WORLD, MPI_COMM_WORLD)

    ! synchro pour s'assurer que le message s'affiche lorsque tout le monde a bien termine
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Comm_Rank (MPI_COMM_WORLD, rg, ierr)
    if (rg==0) then
        write (*,*) "fin du calcul sur processeurs "
    end if
    call mpi_finalize(ierr)

end program SEM3D

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
