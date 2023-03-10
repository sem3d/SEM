!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

subroutine test_path()
    use semdatafiles
    
    character(Len=MAX_FILE_SIZE) :: fnamef
        
    call semname_dir_capteurs(fnamef)
    write(*,*) "dir capteurs='", trim(adjustl(fnamef)),"'"
    call semname_capteur_type("cap00", ".vel", fnamef)
    write(*,*) "capteur_pos:", trim(adjustl(fnamef))
        
end subroutine test_path


program test_sem_path

    use semdatafiles

    character(Len=MAX_FILE_SIZE),parameter :: p_param = "./param"
    character(Len=MAX_FILE_SIZE),parameter :: p_traces = "./traces"
    character(Len=MAX_FILE_SIZE),parameter :: p_results = "./results"
    character(Len=MAX_FILE_SIZE),parameter :: p_data = "./data"
    character(Len=MAX_FILE_SIZE),parameter :: p_prot = "./prot"
    character(Len=MAX_FILE_SIZE),parameter :: p_mat = "./mat"
    character(Len=MAX_FILE_SIZE),parameter :: p_mirror = "./mirror"



    write(*,*) "pjoin: ", trim(adjustl(pjoin("abc","def")))

    call init_sem_path(p_param, p_traces, p_results, p_data, p_prot, p_mat, p_mirror)
    write(*,*) "CHEMINS SEM"
    write(*,*) "==========="
    call test_path()

end program test_sem_path

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
