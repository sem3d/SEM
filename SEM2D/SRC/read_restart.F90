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
    integer :: n,ngllx,ngllz,i,j,ngll
    character (len=100) :: commande
    character (len=6) :: sit !Gsa

    call init_restart(Tdomain%communicateur,Tdomain%Mpi_var%my_rank,Tdomain%TimeD%iter_reprise,file_prot)
    open (61,file=file_prot,status="unknown",form="formatted")
    read(61,*) Tdomain%TimeD%rtime,Tdomain%TimeD%dtmin
    read(61,*) Tdomain%TimeD%NtimeMin,isort

    Tdomain%TimeD%prot_m0 = Tdomain%TimeD%NtimeMin !!GSa 11/2009 !!init des iterations de protection
    Tdomain%TimeD%prot_m1 = Tdomain%TimeD%NtimeMin !!GSa 11/2009

    !preparation pour le pas de temps suivant
    Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin
    Tdomain%TimeD%NtimeMin=Tdomain%TimeD%NtimeMin+1

    ! Save Fields for Elements
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        if ( .not. Tdomain%specel(n)%PML ) then
            do j = 1,ngllz-2
                do i = 1,ngllx-2
                    read(61,*) Tdomain%specel(n)%Veloc(i,j,0)
                    read(61,*) Tdomain%specel(n)%Veloc(i,j,1)
                    read(61,*) Tdomain%specel(n)%Displ(i,j,0)
                    read(61,*) Tdomain%specel(n)%Displ(i,j,1)
                enddo
            enddo
        else
            do j = 1,ngllz-2
                do i = 1,ngllx-2
                    read(61,*) Tdomain%specel(n)%Veloc(i,j,0)
                    read(61,*) Tdomain%specel(n)%Veloc(i,j,1)
                    read(61,*) Tdomain%specel(n)%Veloc1(i,j,0)
                    read(61,*) Tdomain%specel(n)%Veloc1(i,j,1)
                    read(61,*) Tdomain%specel(n)%Veloc2(i,j,0)
                    read(61,*) Tdomain%specel(n)%Veloc2(i,j,1)
                enddo
            enddo

            do j = 0,ngllz-1
                do i = 0,ngllx-1
                    read(61,*) Tdomain%specel(n)%Stress1(i,j,0)
                    read(61,*) Tdomain%specel(n)%Stress1(i,j,1)
                    read(61,*) Tdomain%specel(n)%Stress1(i,j,2)
                    read(61,*) Tdomain%specel(n)%Stress2(i,j,0)
                    read(61,*) Tdomain%specel(n)%Stress2(i,j,1)
                    read(61,*) Tdomain%specel(n)%Stress2(i,j,2)
                enddo
            enddo

        endif
    enddo

    ! Save Fields for Faces
    do n = 0,Tdomain%n_face-1
        ngll = Tdomain%sFace(n)%ngll
        if (.not. Tdomain%sFace(n)%PML ) then
            do i = 1,ngll-2
                read(61,*) Tdomain%sFace(n)%Veloc(i,0)
                read(61,*) Tdomain%sFace(n)%Veloc(i,1)
                read(61,*) Tdomain%sFace(n)%Displ(i,0)
                read(61,*) Tdomain%sFace(n)%Displ(i,1)
#ifdef MKA3D
                read(61,*) Tdomain%sFace(n)%Forces(i,0)
                read(61,*) Tdomain%sFace(n)%Forces(i,1)
#endif
            enddo
        else
            do i = 1,ngll-2
                read(61,*) Tdomain%sFace(n)%Veloc(i,0)
                read(61,*) Tdomain%sFace(n)%Veloc(i,1)
                read(61,*) Tdomain%sFace(n)%Veloc1(i,0)
                read(61,*) Tdomain%sFace(n)%Veloc1(i,1)
                read(61,*) Tdomain%sFace(n)%Veloc2(i,0)
                read(61,*) Tdomain%sFace(n)%Veloc2(i,1)
            enddo
        endif
    enddo



    ! Save Fields for Vertices
    do n = 0,Tdomain%n_vertex-1
        if (.not. Tdomain%sVertex(n)%PML ) then
            do i = 0,1
                read(61,*) Tdomain%sVertex(n)%Veloc(i)
                read(61,*) Tdomain%sVertex(n)%Displ(i)
            enddo
        else
            do i = 0,1
                read(61,*) Tdomain%sVertex(n)%Veloc(i)
                read(61,*) Tdomain%sVertex(n)%Veloc1(i)
                read(61,*) Tdomain%sVertex(n)%Veloc2(i)
            enddo
        endif
    enddo
    close(61)

    if (Tdomain%MPI_var%my_rank == 0) then  !!GSa 11/2009
        write(sit,'(I6)') Tdomain%TimeD%prot_m0
        commande='find ./ProRep/sem -name "Prot*" ! -name "Prot*'//trim(adjustl(sit))//'*" -exec rm -fr  {} \; '
        ! on supprime les fichiers et repertoire de protection autres que celui contenant prot_m0
        !!   write(6,*) 'commande',commande
        call system(commande)
    endif

    return
end subroutine read_restart
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
