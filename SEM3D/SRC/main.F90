program drive_sem
    use mpi
    use semdatafiles

    character(Len=MAX_FILE_SIZE),parameter :: p_param = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_traces = "./traces"
    character(Len=MAX_FILE_SIZE),parameter :: p_results = "./res"
    character(Len=MAX_FILE_SIZE),parameter :: p_data = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_prot = "./prot"

    integer :: ierr, rg

    call init_sem_path(p_param, p_traces, p_results, p_data, p_prot)

    call sem(-1, MPI_COMM_WORLD, MPI_COMM_WORLD)

    if (rg==0) then
        ! synchro pour s'assurer que le message s'affiche lorsque tout le monde a bien termine
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
        write (*,*) "fin du calcul sur processeurs "
    end if
    call mpi_finalize(ierr)

end program drive_sem
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
