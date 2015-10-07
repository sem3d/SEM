!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
subroutine read_restart (Tdomain,isort)
    use protrep
    use sdomain
    use semdatafiles
    use mpi

    implicit none

    type (domain), intent (inout):: Tdomain
    integer,intent (inout)::isort

    ! local variables
    character (len=MAX_FILE_SIZE) :: file_prot

    Tdomain%TimeD%NtimeMin = Tdomain%TimeD%iter_reprise
    call init_restart(Tdomain%communicateur,Tdomain%Mpi_var%my_rank,Tdomain%TimeD%NtimeMin,file_prot)
    call read_write_prot(Tdomain,.true.,Tdomain%TimeD%rtime,Tdomain%TimeD%dtmin,Tdomain%TimeD%NtimeMin,isort)

    ! init des iterations de protection
    Tdomain%TimeD%prot_m0 = Tdomain%TimeD%NtimeMin !!GSa 11/2009
    Tdomain%TimeD%prot_m1 = Tdomain%TimeD%NtimeMin !!GSa 11/2009

    ! preparation pour le pas de temps suivant
    Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin
    Tdomain%TimeD%NtimeMin=Tdomain%TimeD%NtimeMin+1

    if (Tdomain%Mpi_var%my_rank == 0) then
        write (*,'(A,I8,A,f10.6,A)') "SEM : REPRISE (iteration ", Tdomain%TimeD%NtimeMin,", temps ",Tdomain%TimeD%rtime, ")"
    endif

    call clean_prot(Tdomain%TimeD%prot_m0, Tdomain%MPI_var%my_rank)

    return
end subroutine read_restart

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
