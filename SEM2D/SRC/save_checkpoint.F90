!> -*- coding: utf-8 -*-
!!\file save_checkpoint.F90
!!\brief Contient la subroutine save_checkpoint.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief La routine save_checkpoint() assure la sauvegarde des champs
!! (force,vitesse,...). Cette sauvegarde permet une reprise
!! éventuelle.
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
    character (len=MAX_FILE_SIZE) :: prot_file
    character (len=MAX_FILE_SIZE) :: dir_prot, dir_prot_prev, times_file, dir_traces
    character (len=100) :: commande

    integer :: rg

    rg = Tdomain%MPI_var%my_rank

    if (rg == 0) then
        write (*,'(A40,I8,A1,f10.6)') "SEM : protection a iteration et tps:",it," ",rtime
    endif

    ! on doit s'assurer que les processeurs ont fini d ecrire leurs resultats=> synchronisation
    call MPI_Barrier (Tdomain%communicateur, ierr)




    ! recherche et destruction au fur et a mesure des anciennes prots
    if (rg == 0) then

        call semname_protection_iter_dir(it,dir_prot)
        ierr = sem_mkdir(trim(adjustl(dir_prot)))
        if (ierr/=0) then
            write(*,*) "Erreur creating dir:", trim(adjustl(dir_prot))
        endif
        Tdomain%TimeD%prot_m2 = Tdomain%TimeD%prot_m1
        Tdomain%TimeD%prot_m1 = Tdomain%TimeD%prot_m0
        Tdomain%TimeD%prot_m0 = it

        if (Tdomain%TimeD%prot_m2>0) then
            call semname_protection_iter_dir(Tdomain%TimeD%prot_m2,dir_prot_prev)
            commande="rm -Rf "//trim(adjustl(dir_prot_prev))
            call system(commande)        ! suppression de l avant avant derniere protection
        endif

        ! copie du fichier temps.dat dans le rep de protection
        call semname_results_temps_sem(times_file)
        commande="cp "//trim(adjustl(times_file))//" "//trim(adjustl(dir_prot))
        call system(trim(commande))

        ! copie du repertoire des sorties capteurs sem dans le rep de protection
        call semname_dir_capteurs(dir_traces)
        commande="cp -R "//trim(adjustl(dir_traces))//" "//dir_prot
        call system(trim(commande))
    endif
    call MPI_Barrier (Tdomain%communicateur, ierr)

    call semname_protection_iter_rank_file(it,rg,prot_file)
    open (61,file=prot_file,status="unknown",form="formatted")
    write(61,*) rtime,dtmin
    write(61,*) it,isort

    ! Save Fields for Elements
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        if ( (.not.Tdomain%specel(n)%PML) .or. Tdomain%specel(n)%CPML) then
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
        if ((.not.Tdomain%sFace(n)%PML) .or. Tdomain%sFace(n)%CPML) then
            do i = 1,ngll-2
                write(61,*) Tdomain%sFace(n)%Veloc(i,0)
                write(61,*) Tdomain%sFace(n)%Veloc(i,1)
                write(61,*) Tdomain%sFace(n)%Displ(i,0)
                write(61,*) Tdomain%sFace(n)%Displ(i,1)
                write(61,*) Tdomain%sFace(n)%Forces(i,0)
                write(61,*) Tdomain%sFace(n)%Forces(i,1)
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
        if ((.not.Tdomain%sVertex(n)%PML) .or. Tdomain%sVertex(n)%CPML) then
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
