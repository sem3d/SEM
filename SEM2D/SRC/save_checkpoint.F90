!>
!!\file save_checkpoint.F90
!!\brief Contient la subroutine save_checkpoint.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief La routine save_checkpoint() assure la sauvegarde des champs (force,vitesse,...). Cette sauvegarde permet une reprise éventuelle.
!!
!! \param type (domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) it
!! \param integer, intent (IN) isort
!! \param real, intent (IN) rtime
!! \param real, intent (IN) dtmin
!<


subroutine save_checkpoint (Tdomain,rtime,dtmin,it,isort)

    ! rtime : checkpoint time
    ! it    : iteration



    use sdomain
    use semdatafiles
    use mpi

    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: it,isort
    real, intent (IN) :: rtime,dtmin

    ! local variables
    integer :: n,ngllx,ngllz,i,j,ngll,ierr
    character (len=MAX_FILE_SIZE) :: fnamef,fnamer,fnamec

#ifdef MKA3D
    character (len=100) :: commande
#endif

    if (Tdomain%MPI_var%my_rank == 0) then
        write (*,'(A40,I8,A1,f10.6)') "SEM : protection a iteration et tps:",it," ",rtime
    endif

    ! on doit s'assurer que les processeurs ont fini d ecrire leurs resultats=> synchronisation
    call MPI_Barrier (Tdomain%communicateur, ierr)


#ifdef MKA3D

    call semname_couplage_iter(it,Tdomain%Mpi_var%my_rank,fnamef)
    call semname_couplage_iterr(it,fnamer)
    ! creation du repertoire data/sem/Protection_<it> (par le proc le plus rapide)
    commande="mkdir -p "//trim(fnamer)
    call system(commande)


    ! recherche et destruction au fur et a mesure des anciennes prots
    if (Tdomain%MPI_var%my_rank == 0) then

        Tdomain%TimeD%prot_m2 = Tdomain%TimeD%prot_m1
        Tdomain%TimeD%prot_m1 = Tdomain%TimeD%prot_m0
        Tdomain%TimeD%prot_m0 = it

        if (Tdomain%TimeD%prot_m2>0) then
            call semname_couplage_iterr(Tdomain%TimeD%prot_m2,fnamec)
            commande="rm -Rf "//trim(adjustl(fnamec))
            call system(commande)        ! suppression de l avant avant derniere protection
        endif


        ! copie du fichier temps.dat dans le rep de protection
        call semname_save_checkpoint_cpt(fnamer,fnamec)
        commande="cp "//trim(adjustl(fnamec)) !!modif 09/11
        call system(commande)

        ! copie du repertoire des sorties capteurs sem dans le rep de protection
        call semname_save_checkpoint_cp(fnamer,fnamec)
        commande="cp -r "//trim(adjustl(fnamec))
        call system(commande)
    endif

#else
    call semname_read_restart_save_checkpoint_rank(Tdomain%Mpi_var%my_rank,fnamef)
#endif



    open (61,file=fnamef,status="unknown",form="formatted")
    write(61,*) rtime,dtmin
    write(61,*) it,isort

    ! Save Fields for Elements
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        if ( .not. Tdomain%specel(n)%PML ) then
            do j = 1,ngllz-2
                do i = 1,ngllx-2
                    write(61,*) Tdomain%specel(n)%Veloc(i,j,0)
                    write(61,*) Tdomain%specel(n)%Veloc(i,j,1)
                    write(61,*) Tdomain%specel(n)%Displ(i,j,0)
                    write(61,*) Tdomain%specel(n)%Displ(i,j,1)
                enddo
            enddo
        else
            do j = 1,ngllz-2
                do i = 1,ngllx-2
                    write(61,*) Tdomain%specel(n)%Veloc(i,j,0)
                    write(61,*) Tdomain%specel(n)%Veloc(i,j,1)
                    write(61,*) Tdomain%specel(n)%Veloc1(i,j,0)
                    write(61,*) Tdomain%specel(n)%Veloc1(i,j,1)
                    write(61,*) Tdomain%specel(n)%Veloc2(i,j,0)
                    write(61,*) Tdomain%specel(n)%Veloc2(i,j,1)
                enddo
            enddo

            do j = 0,ngllz-1
                do i = 0,ngllx-1
                    write(61,*) Tdomain%specel(n)%Stress1(i,j,0)
                    write(61,*) Tdomain%specel(n)%Stress1(i,j,1)
                    write(61,*) Tdomain%specel(n)%Stress1(i,j,2)
                    write(61,*) Tdomain%specel(n)%Stress2(i,j,0)
                    write(61,*) Tdomain%specel(n)%Stress2(i,j,1)
                    write(61,*) Tdomain%specel(n)%Stress2(i,j,2)
                enddo
            enddo

        endif
    enddo

    ! Save Fields for Faces
    do n = 0,Tdomain%n_face-1
        ngll = Tdomain%sFace(n)%ngll
        if (.not. Tdomain%sFace(n)%PML ) then
            do i = 1,ngll-2
                write(61,*) Tdomain%sFace(n)%Veloc(i,0)
                write(61,*) Tdomain%sFace(n)%Veloc(i,1)
                write(61,*) Tdomain%sFace(n)%Displ(i,0)
                write(61,*) Tdomain%sFace(n)%Displ(i,1)
#ifdef MKA3D
                write(61,*) Tdomain%sFace(n)%Forces(i,0)
                write(61,*) Tdomain%sFace(n)%Forces(i,1)
#endif


            enddo
        else
            do i = 1,ngll-2
                write(61,*) Tdomain%sFace(n)%Veloc(i,0)
                write(61,*) Tdomain%sFace(n)%Veloc(i,1)
                write(61,*) Tdomain%sFace(n)%Veloc1(i,0)
                write(61,*) Tdomain%sFace(n)%Veloc1(i,1)
                write(61,*) Tdomain%sFace(n)%Veloc2(i,0)
                write(61,*) Tdomain%sFace(n)%Veloc2(i,1)
            enddo
        endif
    enddo





    ! Save Fields for Vertices
    do n = 0,Tdomain%n_vertex-1
        if (.not. Tdomain%sVertex(n)%PML ) then
            do i = 0,1
                write(61,*) Tdomain%sVertex(n)%Veloc(i)
                write(61,*) Tdomain%sVertex(n)%Displ(i)
            enddo
        else
            do i = 0,1
                write(61,*) Tdomain%sVertex(n)%Veloc(i)
                write(61,*) Tdomain%sVertex(n)%Veloc1(i)
                write(61,*) Tdomain%sVertex(n)%Veloc2(i)
            enddo
        endif
    enddo
    close(61)



    return
end subroutine save_checkpoint
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
