subroutine read_restart (Tdomain,isort)


    use sdomain
    use semdatafiles
    use mpi

    implicit none

    type (domain), intent (inout):: Tdomain
    integer,intent (inout)::isort

    ! local variables
    integer :: n,ngllx,ngllz,i,j,ngll,ierr
    character (len=MAX_FILE_SIZE) :: fnamef,fnamer,fnamec
    character (len=100) :: commande
    character (len=6) :: sit !Gsa

#ifdef MKA3D

    call semname_couplage_iter(Tdomain%TimeD%iter_reprise,Tdomain%Mpi_var%my_rank,fnamef)
    call semname_couplage_iterr(Tdomain%TimeD%iter_reprise,fnamer)

    if (Tdomain%MPI_var%my_rank == 0) then

        ! copie du fichier temps.dat dans le rep de Resultat
        call semname_couplage_commandecpt(fnamer,fnamec)
        commande="cp "//trim(adjustl(fnamec)) !!modif 09/11
        call system(commande)

        ! copie du repertoire des sorties capteurs sem dans le rep de resultats
        call semname_couplage_commanderm(fnamec)
        commande="rm -Rf "//trim(adjustl(fnamec))
        call system(commande)
        call semname_couplage_commandecp(fnamer,fnamec)
        commande="cp -r "//trim(adjustl(fnamec))
        call system(commande)

    endif

#else
    call semname_read_restart_save_checkpoint_rank(Tdomain%Mpi_var%my_rank,fnamef)
#endif

    ! pour s'assurer que le proc 0 a bien eu le temps de remettre en place tous les fichiers proteges
    call MPI_Barrier(Tdomain%communicateur, ierr)

    open (61,file=fnamef,status="unknown",form="formatted")

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
