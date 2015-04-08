program mesher_mat

    !- Building of mesh files compatible with SEM 3D

    use mesh2spec
    implicit none

    integer                   :: n_parts
    logical, parameter        :: VRAI = .true., FAUX = .false.
    integer, parameter        :: NMAX_PROCS = 8192


    write(*,*) "-----------------------------------------------------------------------"
    write(*,*) "-----------------------------------------------------------------------"
    write(*,*) "----------                                                 ------------"
    write(*,*) "----------  Construction of input files for SEM runs (3D)  ------------"
    write(*,*) "----------                                                 ------------"
    write(*,*) "-----------------------------------------------------------------------"
    write(*,*) "-----------------------------------------------------------------------"
    write(*,*)
    write(*,*) "  --> How many procs for the run?"
    read(*,*) n_parts
    if(n_parts < 1 .or. n_parts > NMAX_PROCS)   &
        stop "In mesher_mat: sure, for the proc number?"
    write(*,*)
    call gen_mesh(n_parts)


end program mesher_mat
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
