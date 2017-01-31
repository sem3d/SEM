!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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
!! \param real(fpp), intent (IN) rtime
!! \param real(fpp), intent (IN) dtmin
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
    real(fpp), intent (IN) :: rtime,dtmin

    ! local variables
    integer :: ierr
    character (len=MAX_FILE_SIZE) :: dir_prot, dir_prot_prev, times_file
    character (len=MAX_FILE_SIZE) :: dir_traces, dir_prot_traces
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
        call semname_protection_iter_dir_capteurs(it,dir_prot_traces)
        commande="cp -R "//trim(adjustl(dir_traces))//" "//dir_prot_traces
        call system(trim(commande))
    endif
    call MPI_Barrier (Tdomain%communicateur, ierr)

    call read_write_prot(Tdomain,.false.,rtime,dtmin,it,isort)
    return
end subroutine save_checkpoint

! read_write_prot :
! 1. make sure read / write are done the exact same way (save_chekpoint // read_restart)
! 2. use "read(61,*) tab(:)" instead of "do i=1,n read(61,*) tab(i)" : this reduces file size, and, improves IO time
subroutine read_write_prot(Tdomain,readprot,rtime,dtmin,it,isort)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    logical,      intent(in)     :: readprot
    real(fpp)                    :: rtime,dtmin
    integer                      :: it,isort

    character(len=MAX_FILE_SIZE) :: prot_file
    integer                      :: n, ngll, ngllx, ngllz

    call semname_protection_iter_rank_file(it,Tdomain%MPI_var%my_rank,prot_file)
    open (61,file=prot_file,status="unknown",form="formatted")
    if (readprot .eqv. .true.) then
        read (61,*) rtime,dtmin
        read (61,*) it,isort
    else
        write(61,*) rtime,dtmin
        write(61,*) it,isort
    endif

    ! Save Fields for Elements
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        if ( (.not.Tdomain%specel(n)%PML) .or. Tdomain%specel(n)%CPML .or. Tdomain%specel(n)%ADEPML) then
            if (readprot .eqv. .true.) then
                read (61,*) Tdomain%specel(n)%Veloc(1:ngllx-2,1:ngllz-2,0:1)
                read (61,*) Tdomain%specel(n)%Displ(1:ngllx-2,1:ngllz-2,0:1)
            else
                write(61,*) Tdomain%specel(n)%Veloc(1:ngllx-2,1:ngllz-2,0:1)
                write(61,*) Tdomain%specel(n)%Displ(1:ngllx-2,1:ngllz-2,0:1)
            endif
        else
            if (readprot .eqv. .true.) then
                read (61,*) Tdomain%specel(n)%Veloc  (1:ngllx-2,1:ngllz-2,0:1)
                read (61,*) Tdomain%specel(n)%Veloc1 (1:ngllx-2,1:ngllz-2,0:1)
                read (61,*) Tdomain%specel(n)%Veloc2 (1:ngllx-2,1:ngllz-2,0:1)
                read (61,*) Tdomain%specel(n)%Stress1(1:ngllx-2,1:ngllz-2,0:2)
                read (61,*) Tdomain%specel(n)%Stress2(1:ngllx-2,1:ngllz-2,0:2)
            else
                write(61,*) Tdomain%specel(n)%Veloc  (1:ngllx-2,1:ngllz-2,0:1)
                write(61,*) Tdomain%specel(n)%Veloc1 (1:ngllx-2,1:ngllz-2,0:1)
                write(61,*) Tdomain%specel(n)%Veloc2 (1:ngllx-2,1:ngllz-2,0:1)
                write(61,*) Tdomain%specel(n)%Stress1(1:ngllx-2,1:ngllz-2,0:2)
                write(61,*) Tdomain%specel(n)%Stress2(1:ngllx-2,1:ngllz-2,0:2)
            endif
        endif
    enddo

    ! Save Fields for Faces
    do n = 0,Tdomain%n_face-1
        ngll = Tdomain%sFace(n)%ngll
        if ((.not.Tdomain%sFace(n)%PML) .or. Tdomain%sFace(n)%CPML .or. Tdomain%sFace(n)%ADEPML) then
            if (readprot .eqv. .true.) then
                read (61,*) Tdomain%sFace(n)%Veloc (1:ngll-2,0:1)
                read (61,*) Tdomain%sFace(n)%Displ (1:ngll-2,0:1)
                read (61,*) Tdomain%sFace(n)%Forces(1:ngll-2,0:1)
            else
                write(61,*) Tdomain%sFace(n)%Veloc (1:ngll-2,0:1)
                write(61,*) Tdomain%sFace(n)%Displ (1:ngll-2,0:1)
                write(61,*) Tdomain%sFace(n)%Forces(1:ngll-2,0:1)
            endif
        else
            if (readprot .eqv. .true.) then
                read (61,*) Tdomain%sFace(n)%Veloc (1:ngll-2,0:1)
                read (61,*) Tdomain%sFace(n)%Veloc1(1:ngll-2,0:1)
                read (61,*) Tdomain%sFace(n)%Veloc2(1:ngll-2,0:1)
            else
                write(61,*) Tdomain%sFace(n)%Veloc (1:ngll-2,0:1)
                write(61,*) Tdomain%sFace(n)%Veloc1(1:ngll-2,0:1)
                write(61,*) Tdomain%sFace(n)%Veloc2(1:ngll-2,0:1)
            endif
        endif
    enddo

    ! Save Fields for Vertices
    do n = 0,Tdomain%n_vertex-1
        if ((.not.Tdomain%sVertex(n)%PML) .or. Tdomain%sVertex(n)%CPML .or. Tdomain%sVertex(n)%ADEPML) then
            if (readprot .eqv. .true.) then
                read (61,*) Tdomain%sVertex(n)%Veloc(0:1)
                read (61,*) Tdomain%sVertex(n)%Displ(0:1)
            else
                write(61,*) Tdomain%sVertex(n)%Veloc(0:1)
                write(61,*) Tdomain%sVertex(n)%Displ(0:1)
            endif
        else
            if (readprot .eqv. .true.) then
                read (61,*) Tdomain%sVertex(n)%Veloc (0:1)
                read (61,*) Tdomain%sVertex(n)%Veloc1(0:1)
                read (61,*) Tdomain%sVertex(n)%Veloc2(0:1)
            else
                write(61,*) Tdomain%sVertex(n)%Veloc (0:1)
                write(61,*) Tdomain%sVertex(n)%Veloc1(0:1)
                write(61,*) Tdomain%sVertex(n)%Veloc2(0:1)
            endif
        endif
    enddo

    ! Save n_quad for xdmf results visualisation
    if (readprot .eqv. .true.) then
        read (61,*) Tdomain%n_quad
    else
        write(61,*) Tdomain%n_quad
    endif

    close(61)
end subroutine read_write_prot

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
