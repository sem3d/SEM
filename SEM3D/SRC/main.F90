program drive_sem

    use semdatafiles

    character(Len=MAX_FILE_SIZE),parameter :: p_param = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_traces = "./traces"
    character(Len=MAX_FILE_SIZE),parameter :: p_results = "./res"
    character(Len=MAX_FILE_SIZE),parameter :: p_data = "."
    character(Len=MAX_FILE_SIZE),parameter :: p_prot = "./prot"
    character(Len=MAX_FILE_SIZE) :: fnamef

    integer :: code, rg

    call init_sem_path(p_param, p_traces, p_results, p_data, p_prot)

    call MPI_Init (code)
    call MPI_Comm_Rank (MPI_COMM_WORLD, rg, code)

    call sem(-1, MPI_COMM_WORLD, MPI_COMM_WORLD)

    if (rg == 0) write (*,*) "fin du calcul sur processeurs ",rg
    call mpi_finalize(code)

end program drive_sem
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
