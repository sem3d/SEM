!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module protrep
    use semdatafiles
    use sem_mpi
contains
    !>
    !! \brief Assure la reprise par la lecture des fichiers de protection.
    !!
    !<
    subroutine init_restart(comm, rg, iter, file_prot)
        implicit none
        integer, intent (in):: comm, rg, iter
        character (len=MAX_FILE_SIZE), intent(out) :: file_prot
        character (len=MAX_FILE_SIZE) :: dir_prot,  file_temps_sem
        character (len=MAX_FILE_SIZE) :: dir_prot_capteurs, dir_capteurs, file_prot_temps_sem
        character (len=MAX_FILE_SIZE) :: commande
        integer :: ierr

        call semname_protection_iter_rank_file(iter,rg,file_prot)

        if (rg == 0) then

            call semname_protection_iter_dir(iter,dir_prot)

            ! copie du fichier temps.dat dans le rep de Resultat
            call semname_results_temps_sem(file_temps_sem)
            call semname_protection_temps_sem(iter,file_prot_temps_sem)
            commande="cp "//trim(adjustl(file_prot_temps_sem))//" "//trim(adjustl(file_temps_sem))
            call system(commande)

            ! copie du repertoire des sorties capteurs sem dans le rep de resultats
            !! Suppression du repertoire existant
            call semname_dir_capteurs(dir_capteurs)
            commande="rm -Rf "//trim(adjustl(dir_capteurs))
            call system(commande)
            !! Copie du repertoire des protections vers le repertoire dest
            call semname_protection_iter_dir_capteurs(iter,dir_prot_capteurs)
            commande="cp -r "//trim(adjustl(dir_prot_capteurs))//" "//trim(adjustl(dir_capteurs))
            call system(commande)

        endif

        ! pour s'assurer que le proc 0 a bien eu le temps de remettre en place tous les fichiers proteges
        call MPI_Barrier(comm, ierr)
    end subroutine init_restart

    subroutine clean_prot(protit, rg)
        implicit none
        integer, intent (IN)       :: protit
        integer, intent (IN)             :: rg
        character (len=MAX_FILE_SIZE)    :: commande
        character (len=6) :: sit

        if (rg == 0) then
            write(sit,'(I6)') protit
            commande='find ./ProRep/sem -name "Prot*" ! -name "Prot*'//trim(adjustl(sit))//'*" -exec rm -fr  {} \; '
            call system(commande)
        endif
    end subroutine clean_prot



end module protrep

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
