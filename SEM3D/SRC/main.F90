program drive_sem

    integer :: code, rg


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
