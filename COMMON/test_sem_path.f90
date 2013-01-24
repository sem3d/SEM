
subroutine test_path()
    use semdatafiles
    
    character(Len=MAX_FILE_SIZE) :: fnamef
        
    call semname_dir_capteurs(fnamef)
    write(*,*) "dir capteurs='", trim(adjustl(fnamef)),"'"
    call semname_capteur_pos("cap00", fnamef)
    write(*,*) "capteur_pos:", trim(adjustl(fnamef))
    call semname_capteur_fichiercapteur(fnamef)
    write(*,*) "capteur.dat:", trim(adjustl(fnamef))
        
end subroutine test_path


program test_sem_path

    use semdatafiles

    character(Len=MAX_FILE_SIZE),parameter :: p_param = "./param"
    character(Len=MAX_FILE_SIZE),parameter :: p_traces = "./traces"
    character(Len=MAX_FILE_SIZE),parameter :: p_results = "./results"
    character(Len=MAX_FILE_SIZE),parameter :: p_data = "./data"
    character(Len=MAX_FILE_SIZE),parameter :: p_prot = "./prot"



    write(*,*) "pjoin: ", trim(adjustl(pjoin("abc","def")))

    call init_sem_path(p_param, p_traces, p_results, p_data, p_prot)
    write(*,*) "CHEMINS SEM"
    write(*,*) "==========="
    call test_path()

    call init_mka3d_path()
    write(*,*) "CHEMINS MKA3D"
    write(*,*) "============="
    call test_path()

end program test_sem_path
