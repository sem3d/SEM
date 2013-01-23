!>
!!\file semdatafiles.F90
!!\brief Permet de gerer les noms de fichier des repertoire SEM
!!\author Coupey Quentin
!!\version 1.0
!!\date 23/08/2011
!!
!<

#if 0
#define DEBUG(f) write(*,"(a)") f
#else
#define DEBUG(f)
#endif

module semdatafiles
    integer, parameter :: MAX_FILE_SIZE=1024
    !character(Len=20) :: datadir
    character(Len=MAX_FILE_SIZE) :: path_param
    character(Len=MAX_FILE_SIZE) :: path_traces
    character(Len=MAX_FILE_SIZE) :: path_results
    character(Len=MAX_FILE_SIZE) :: path_data
    character(Len=MAX_FILE_SIZE) :: path_prot
contains

    function pjoin(s1, s2)
        implicit none
        character(Len=*) :: s1, s2
        character(Len=MAX_FILE_SIZE) :: pjoin
        integer :: l1
        l1 = len_trim(s1)
        if (s1(l1:l1)=="/") then
            pjoin = trim(adjustl(s1))//trim(adjustl(s2))
        else
            pjoin = trim(adjustl(s1))//"/"//trim(adjustl(s2))
        end if
    end function pjoin

    subroutine init_mka3d_path()
        path_param = "./Parametrage/sem"
        path_traces = "./Capteurs/sem"
        path_results = "./Resultats"
        path_data = "./data"
        path_prot = "./ProRep/sem"
    end subroutine init_mka3d_path

    subroutine init_sem_path(param, traces, results, data, prorep)
        character(Len=MAX_FILE_SIZE), intent(in) :: param
        character(Len=MAX_FILE_SIZE), intent(in) :: traces
        character(Len=MAX_FILE_SIZE), intent(in) :: results
        character(Len=MAX_FILE_SIZE), intent(in) :: data
        character(Len=MAX_FILE_SIZE), intent(in) :: prorep

        path_param = param
        path_traces = traces
        path_results = results
        path_data = data
        path_prot = prorep
    end subroutine init_sem_path


    subroutine semname_dir_capteurs(dirname)
        implicit none
        character(Len=MAX_FILE_SIZE), intent(out) :: dirname
        dirname = path_traces
    end subroutine semname_dir_capteurs

    !!fichier capteur 2d 3d
    subroutine semname_capteur_pos (name,fnamef)
        !SEMFILE 98 RW ./Capteurs/sem/capteurs_XXX.position (MKA) & "capteur_XXX.position (NOMKA)
        implicit none
        character(Len=*),intent(in) :: name
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        fnamef = pjoin(path_traces, "capteurs_"//trim(adjustl(name))//".position")
    end subroutine semname_capteur_pos

    subroutine semname_capteur_type (name,type,fnamef)
        !SEMFILE 99 W ./Capteurs/sem/XXX_YYY (MKA) & XXX_YYY (NOMKA)
        implicit none
        character(Len=*), intent(in) :: name
        character(Len=*),intent(in) :: type
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: dirname

        call semname_dir_capteurs(dirname)

        write(fnamef,"(a,a,a,a)") trim(adjustl(dirname)),"/",trim(adjustl(name)),trim(adjustl(type))

        DEBUG(fnamef)
    end subroutine semname_capteur_type

    subroutine semname_capteur_fichiercapteur (fnamef)
        !SEMFILE 200 R ./Parametrage/sem/capteurs.dat
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        fnamef = pjoin(path_param,"capteurs.dat")
    end subroutine semname_capteur_fichiercapteur
    !!end fichier capteur 2d 3d

    !!fichier couplage 2d 3d
    subroutine semname_couplage_iter (iter,rank,fnamef) !! egalement read restart et save checkpoint
        !SEMFILE 61 R ./ProRep/sem/Protection_III/Protection_III.JJJ
        implicit none
        integer,intent(in) :: iter
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a24,I8.8,a12,I8.8,a1,I4.4)")"./ProRep/sem/Protection_",iter,"/Protection_",iter,".",rank

        DEBUG(fnamef)
    end subroutine semname_couplage_iter
    subroutine semname_couplage_iterr (iter,fnamef)
        !SEMFILE 61 R ./ProRep/sem/Protection_III
        implicit none
        integer,intent(in) :: iter
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a24,I8.8)")"./ProRep/sem/Protection_",iter

        DEBUG(fnamef)
    end subroutine semname_couplage_iterr
    subroutine semname_couplage_commandecpt(fnamer,fnamef)
        !SEMFILE * C XXX/temps_sem.dat ./Resultats/temps_sem.dat
        implicit none
        character(Len=*),intent(in) :: fnamer
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,a)") trim(adjustl(fnamer)),"/temps_sem.dat ./Resultats/temps_sem.dat"

        DEBUG(fnamef)
    end subroutine semname_couplage_commandecpt
    subroutine semname_couplage_commanderm(fnamef)
        !SEMFILE * C ./Capteurs/sem
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a)")"./Capteurs/sem"

        DEBUG(fnamef)
    end subroutine semname_couplage_commanderm
    subroutine semname_couplage_commandecp(fnamer,fnamef)
        !SEMFILE * C XXX/Capteurs ./Capteurs/sem
        implicit none
        character(Len=*),intent(in) :: fnamer
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,a)") trim(adjustl(fnamer)),"/Capteurs ./Capteurs/sem"

        DEBUG(fnamef)
    end subroutine semname_couplage_commandecp
    subroutine semname_couplage_listepts (rank,fnamef)
        !SEMFILE 75 W liste_pts_interp.dat.III
        implicit none
        integer,intent(in) :: rank
        !character(Len=20) :: crank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        !write(crank,"(I20)") rank
        write(fnamef,"(a,I1)")"liste_pts_interp.dat.",rank

        DEBUG(fnamef)
    end subroutine semname_couplage_listepts

    subroutine semname_couplage_listeproj (rank,fnamef)
        !SEMFILE 77 W liste_proj_mka.dat.III
        implicit none
        integer,intent(in) :: rank
        !character(Len=20) :: crank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        !write(crank,"(I20)") rank
        write(fnamef,"(a,I1)")"liste_proj_mka.dat.",rank

        DEBUG(fnamef)
    end subroutine semname_couplage_listeproj

    subroutine semname_couplage_verif1 (rank,fnamef)
        !SEMFILE 85 W verifW.dat.III
        implicit none
        character(Len=*),intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        write(fnamef,"(a,a1,a)")"verifW.dat",'.',trim(adjustl(rank))

        DEBUG(fnamef)
    end subroutine semname_couplage_verif1

    subroutine semname_couplage_verifglob (fnamef)
        !SEMFILE 93 W verif_mka_glob.dat
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        write(fnamef,"(a18)")"verif_mka_glob.dat"

        DEBUG(fnamef)
    end subroutine semname_couplage_verifglob

    subroutine semname_couplage_NCdat (fnamef) !3d
        !SEMFILE 51 W tableau_NC.dat
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        write(fnamef,"(a)")"tableau_NC.dat"

        DEBUG(fnamef)
    end subroutine semname_couplage_NCdat

    subroutine semname_couplage_imatdat (fnamef) !3d
        !SEMFILE 58 W tableau_imat.dat
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        write(fnamef,"(a)")"tableau_imat.dat"

        DEBUG(fnamef)
    end subroutine semname_couplage_imatdat
    !!end fichier couplage 2d 3d

    !!fichier define_fault_properties 2d
    subroutine semname_define_fault_data (domain,fnamef)
        !SEMFILE 24 R ./data/sem/XXX (MKA) & XXX (NOMKA)
        implicit none
        character(Len=*),intent(in) :: domain
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        fnamef = pjoin(path_data, trim(adjustl(domain)))
    end subroutine semname_define_fault_data

    subroutine semname_define_fault_rankl (rank,fnamef)
        !SEMFILE 23 W ./Resultats/initial.III (MKA) initial.III (NOMKA)
        implicit none
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=20) :: temp

        write(temp,"(a,I4.4)") "initial.",rank
        fnamef = pjoin(path_results, temp)
    end subroutine semname_define_fault_rankl

    subroutine semname_define_fault_rankn (rank,fnamef)
        !SEMFILE 24 W "./Resultats/initian.III (MKA) initian.III (NOMKA)
        implicit none
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
#ifdef MKA3D
        write(fnamef,"(a,I4.4)")"./Resultats/initian.",rank
#else
        write(fnamef,"(a,I4.4)") "initian.",rank
#endif

        DEBUG(fnamef)
    end subroutine semname_define_fault_rankn
    !!end fichier define_fault_properties 2d

    !!fichier define_neu_properties
    subroutine semname_define_neu_properties_file (file,fnamef)
        !SEMFILE 24 R XXX
        implicit none
        character(Len=*),intent(in) :: file
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        fnamef = pjoin(path_param, file)
    end subroutine semname_define_neu_properties_file
    !!end fichier define_neu_properties

    !!fichier define_planew_properties
    subroutine semname_define_planew_properties_file (file,fnamef)
        !SEMFILE 24 R XXX
        implicit none
        character(Len=*),intent(in) :: file
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        fnamef = pjoin(path_param, file)
    end subroutine semname_define_planew_properties_file
    !!end fichier define_planew_properties

    !!fichier dumptrace 2d
    subroutine semname_dumptrace_tracex (i,fnamef)
        !SEMFILE 31 W ./Capteurs/sem/tracexIII (MKA) tracexIII (NOMKA)
        implicit none
        integer,intent(in) :: i
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=20) :: temp

        write(temp,"(a,I3.3)") "tracex",i
        DEBUG(fnamef)
        fnamef = pjoin(path_traces, temp)
    end subroutine semname_dumptrace_tracex

    subroutine semname_dumptrace_tracez (i,fnamef)
        !SEMFILE 32 W ./Capteurs/sem/tracezIII (MKA) tracezIII (NOMKA)
        implicit none
        integer,intent(in) :: i
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=20) :: temp

        write(temp,"(a,I3.3)") "tracez",i
        DEBUG(fnamef)
        fnamef = pjoin(path_traces, temp)
    end subroutine semname_dumptrace_tracez
    !!end fichier dumptrace 2d

    !!fichier main 2d
    subroutine semname_couplage_dat (fnamef)
        !SEMFILE 50 R data/sem/couplage.dat (COUPLAGE)
        ! main 2d, driver_sem3d
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        write(fnamef,"(a)")"Parametrage/voronoi/couplage.dat"

        DEBUG(fnamef)
    end subroutine semname_couplage_dat

    subroutine semname_main_temp (fnamef)
        !SEMFILE 78 W ./Resultats/temps_sem.dat (MKA)
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        write(fnamef,"(a)")"./Resultats/temps_sem.dat"

        DEBUG(fnamef)
    end subroutine semname_main_temp

    subroutine semname_main_result(cit,fnamef)
        !SEMFILE * C ./Resultats/RsemIII
        implicit none
        character(Len=*),intent(in) :: cit
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp
        write(temp,"(a,a)") "Rsem",trim(adjustl(cit))
        fnamef = pjoin(path_results, temp)
    end subroutine semname_main_result
    !!end fichier main 2d

    !!fichier plot_grid 2d
    subroutine semname_plot_grid_grid (crank,fnamef)
        !SEMFILE 22 W ./Resultats/grid.gnu.III (MKA) & "grid.gnu.III (NOMKA)
        implicit none
        integer,intent(in) :: crank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp

        write(temp,"(a,I4.4)") "grid.gnu.",crank
        fnamef = pjoin(path_results, temp)
    end subroutine semname_plot_grid_grid

    subroutine semname_plot_grid_surface (crank,fnamef)
        !SEMFILE 23 W ./data/sem/surface.gnu.III (MKA) & surface.gnu.III (NOMKA)
        implicit none
        integer,intent(in) :: crank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
#ifdef MKA3D
        write(fnamef,"(a,I4.4)")"./data/sem/surface.gnu.",crank
#else
        write(fnamef,"(a,I4.4)")"surface.gnu.",crank
#endif

        DEBUG(fnamef)
    end subroutine semname_plot_grid_surface
    !!end fichier plot_grid 2d

    subroutine semname_snap_geom_file (srank,fnamef)
        implicit none
        integer,intent(in) :: srank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp

        write(temp,"(a,I4.4,a)") "geometry",srank,".h5"
        fnamef = pjoin(path_results, temp)
    end subroutine semname_snap_geom_file

    subroutine semname_snap_result_file (srank,isort,fnamef)
        implicit none
        integer,intent(in) :: srank, isort
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp

        write(temp,"(a,I4.4,a,I4.4,a)") "Rsem",isort,"/sem_field.",srank,".h5"
        fnamef = pjoin(path_results, temp)
    end subroutine semname_snap_result_file


    !! Nom du repertoire de sortie d'un pas de temps
    subroutine semname_snap_result_dir(isort,fnamef)
        implicit none
        integer,intent(in) :: isort
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp

        write(temp,"(a,I4.4)") "Rsem",isort
        fnamef = pjoin(path_results, temp)
    end subroutine semname_snap_result_dir

    !!fichier posttraitement 2d 3d
    subroutine semname_posttraitement_geo (srank,fnamef)
        !SEMFILE 200+i W ./Resultats/ensightIII.geo (MKA) & ensightIII.geo (NOMKA)
        implicit none
        integer,intent(in) :: srank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp

        write(temp,"(a,I4.4,a)")"ensight",srank,".geo"
        fnamef = pjoin(path_results, temp)
    end subroutine semname_posttraitement_geo



    subroutine semname_posttraitement_case (srank,fnamef)
        !SEMFILE 100 W ./Resultats/ensightIII.case (MKA) & ensightIII.case (NOMKA)
        implicit none
        integer,intent(in) :: srank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp

        write(fnamef,"(a,I4.4,a)")"ensight",srank,".case"
        fnamef = pjoin(path_results, temp)
    end subroutine semname_posttraitement_case
    !!end fichier postraitement 2d 3d

    !!fichier read_input 2d 3d
    subroutine semname_read_input_input (fnamef)
        !SEMFILE 11 R ./Parametrage/sem/input.spec (MKA) & input.spec (NOMKA)
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        fnamef = pjoin(path_param, "input.spec")
    end subroutine semname_read_input_input

    subroutine semname_read_input_source (fnamef)
        !SEMFILE 21 R ./Parametrage/sem/source.dat
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        fnamef = pjoin(path_param, "source.dat")
    end subroutine semname_read_input_source

    subroutine semname_read_input_amortissement (fnamef) !3d
        !SEMFILE 22 R ./Parametrage/sem/amortissement.dat (MKA)
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        fnamef = pjoin(path_param, "amortissement.dat")
    end subroutine semname_read_input_amortissement

    subroutine semname_read_input_spec (fnamef)
        !SEMFILE 91 W ./data/sem/input_spec_echo (MKA) & input_spec_echo (NOMKA)
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
#ifdef MKA3D
        write(fnamef,"(a)")"./data/sem/input_spec_echo"
#else
        write(fnamef,"(a)")"input_spec_echo"
#endif

        DEBUG(fnamef)
    end subroutine semname_read_input_spec

    subroutine semname_read_input_meshfile (rg,meshfile,fnamef)
        !SEMFILE 12 R ./data/sem/XXX.III (MKA)
        implicit none
        integer,intent(in) :: rg
        character(Len=*),intent(in) :: meshfile
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp
        write(temp,"(a,a1,I4.4)") trim(adjustl(meshfile)),".",rg
        fnamef = pjoin(path_data, temp)
    end subroutine semname_read_input_meshfile

    subroutine semname_read_inputmesh_parametrage (file,fnamef) !!valable egalement pour le fichier read_mesh
        !SEMFILE 13 R ./Parametrage/sem/XXX (MKA) & XXX (NOMKA)
        implicit none
        character(Len=*),intent(in) :: file
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        fnamef = pjoin( path_param, file)
    end subroutine semname_read_inputmesh_parametrage

    !!end fichier read_input 2d 3d

    !!fichier read_mesh 2d
    subroutine semname_read_mesh_rank (mesh,rank,fnamef)
        !SEMFILE 12 R ./data/sem/XXX.III (MKA) & XXX.III (NOMKA)
        implicit none
        integer,intent(in) :: rank
        character(Len=*),intent(in) :: mesh
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
#ifdef MKA3D
        write(fnamef,"(a11,a,a1,I4.4)")"./data/sem/",trim(adjustl(mesh)),".",rank
#else
        write(fnamef,"(a,a1,I4.4)")trim(adjustl(mesh)),".",rank
#endif

        DEBUG(fnamef)
    end subroutine semname_read_mesh_rank

    subroutine semname_read_mesh_echo (rank,fnamef)
        !SEMFILE 92 W ./data/sem/mesh_echoIII (MKA) & mesh_echoIII (NOMKA)
        implicit none
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
#ifdef MKA3D
        write(fnamef,"(a,I4.4)")"./data/sem/mesh_echo",rank
#else
        write(fnamef,"(a,I4.4)")"mesh_echo",rank
#endif

        DEBUG(fnamef)
    end subroutine semname_read_mesh_echo

    subroutine semname_read_mesh_material_echo (fnamef)
        !SEMFILE 93 W ./data/sem/material_echo (MKA) & material_echo (NOMKA)
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
#ifdef MKA3D
        write(fnamef,"(a,I4.4)")"./data/sem/material_echo"
#else
        write(fnamef,"(a,I4.4)")"material_echo"
#endif

        DEBUG(fnamef)
    end subroutine semname_read_mesh_material_echo

    subroutine semname_read_mesh_station_echo (fnamef)
        !SEMFILE 94 W ./data/sem/station_file_echo (MKA) & station_file_echo (NOMKA)
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
#ifdef MKA3D
        write(fnamef,"(a,I4.4)")"./data/sem/station_file_echo"
#else
        write(fnamef,"(a,I4.4)")"station_file_echo"
#endif

        DEBUG(fnamef)
    end subroutine semname_read_mesh_station_echo
    !!end fichier read_mesh 2d

    !!fichier read_restart 2d 3d & save_checkpoint 2d 3d
    subroutine semname_read_restart_save_checkpoint_protectionf (iter,rank,fnamef)
        !SEMFILE: 61  R W  ./ProRep/sem/Protection_III/Protection_III.JJJ (MKA) : donnee protections
        implicit none
        integer,intent(in) :: rank
        integer,intent(in) :: iter
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a24,I8.8,a12,I8.8,a1,I4.4)") "./ProRep/sem/Protection_",iter,"/Protection_",iter,".",rank

        DEBUG(fnamef)
    end subroutine semname_read_restart_save_checkpoint_protectionf
    subroutine semname_read_restart_save_checkpoint_protectionr (iter,fnamef)
        !SEMFILE: 61  R W  ./ProRep/sem/Protection_III (MKA) : donnee protections
        implicit none
        integer,intent(in) :: iter
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a24,I8.8)") "./ProRep/sem/Protection_",iter

        DEBUG(fnamef)
    end subroutine semname_read_restart_save_checkpoint_protectionr
    subroutine semname_read_restart_save_checkpoint_rank (rank,fnamef)
        !SEMFILE 61 R W backupIII (NOMKA)
        implicit none
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a6,I4.4)") "backup",rank

        DEBUG(fnamef)
    end subroutine semname_read_restart_save_checkpoint_rank

    subroutine semname_save_checkpoint_rank (rank,fnamef) !seulement pour save_checkpoint 3d
        !SEMFILE 61 W ./data/SBackup/backupIII (NOMKA)
        implicit none
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write (fnamef,"(a,I4.4)") "./data/SBackup/backup",rank

        DEBUG(fnamef)
    end subroutine semname_save_checkpoint_rank
    subroutine semname_save_checkpoint_cpt(fnamer,fnamef)
        character,intent(in) :: fnamer
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,a)") "./Resultats/temps_sem.dat ",trim(adjustl(fnamer))

        DEBUG(fnamef)
    end subroutine semname_save_checkpoint_cpt
    subroutine semname_save_checkpoint_cp(fnamer,fnamef)
        character(Len=*),intent(in) :: fnamer
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,a,a)") "./Capteurs/sem ",trim(adjustl(fnamer)),"/Capteurs"

        DEBUG(fnamef)
    end subroutine semname_save_checkpoint_cp
    subroutine semname_save_checkpoint_cp3(fnamer,fnamef)
        character(Len=*),intent(in) :: fnamer
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,a,a)")trim(adjustl(fnamer)),";cp ./Resultats/temps_sem.dat ",trim(adjustl(fnamer))

        DEBUG(fnamef)
    end subroutine semname_save_checkpoint_cp3
    subroutine semname_save_checkpoint_cp2(fnamer,fnamef)
        character(Len=*),intent(in) :: fnamer
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,a,a,a)")trim(adjustl(fnamer)),";cp -r ./Capteurs/sem ",trim(adjustl(fnamer)),"/Capteurs"

        DEBUG(fnamef)
    end subroutine semname_save_checkpoint_cp2
    subroutine semname_read_restart_commandedl(sit,fnamef)
        character(Len=*),intent(in) :: sit
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,a,a)")'./ProRep/sem -name "Prot*" ! -name "Prot*',trim(adjustl(sit)),'*" -exec rm -fr  {} \; '

        DEBUG(fnamef)
    end subroutine semname_read_restart_commandedl
    !!end fichier read_restart 2d 3d

    !!fichier save_deformation 2d
    subroutine semname_save_deformation_datdef (sorties,it,rank,fnamef)
        !SEMFILE 71 W XXX/devol_III.dat.JJJ (MKA)
        implicit none
        character(Len=*),intent(in) :: sorties
        integer,intent(in) :: it
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=20) :: Ait
        write(Ait,"(I20)") it
        write(fnamef,"(a,a,a,a,I4.4)")trim(adjustl(sorties)),"/defvol_",trim(adjustl(Ait)),".dat.",rank
        DEBUG(fnamef)
    end subroutine semname_save_deformation_datdef

    subroutine semname_save_deformation_datepsilon (sorties,it,rank,fnamef)
        !SEMFILE 72 W XXX/epsilon_minus_III.dat.JJJ (MKA)
        implicit none
        character(Len=*),intent(in) :: sorties
        integer,intent(in) :: rank
        integer,intent(in) :: it
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=20) :: Ait
        write(Ait,"(I20)") it
        write(fnamef,"(a,a,a,a,I4.4)")trim(adjustl(sorties)),"/epsilon_minus_",trim(adjustl(Ait)),".dat.",rank

        DEBUG(fnamef)
    end subroutine semname_save_deformation_datepsilon

    subroutine semname_save_deformation_epsilonp (rank,it,fnamef)
        !SEMFILE 71 W III.JJJ (NOMKA)
        implicit none
        integer ,intent(in) :: rank
        integer ,intent(in) :: it
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        !character(Len=20) :: Ait
        !write(Ait,"(I20)") it
        write(fnamef,"(I2.2,a1,I4)")rank,".",it

        DEBUG(fnamef)
    end subroutine semname_save_deformation_epsilonp

    subroutine semname_save_deformation_epsilonm (rank,it,fnamef) !!maybe a degager
        !SEMFILE 72 W III.JJJ (NOMKA)
        implicit none
        integer ,intent(in) :: rank
        integer ,intent(in) :: it
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        !character(Len=20) :: Ait
        !write(Ait,"(I20)") it
        write(fnamef,"(I2.2,a1,I4)")rank,".",it

        DEBUG(fnamef)
    end subroutine semname_save_deformation_epsilonm
    !!end fichier save_deformation 2d

    !!fichier save_fault_trace 2d
    subroutine semname_save_fault_trace_rankit (path,rank,it,fnamef)
        !SEMFILE 62->66 W XXXIII.JJJ (MKA & NOMKA)
        implicit none
        integer ,intent(in) :: rank
        integer ,intent(in) :: it
        character(Len=*),intent(in) :: path
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,I4.4,a1,I5.5)")trim(adjustl(path)),rank,".", it

        DEBUG(fnamef)
    end subroutine semname_save_fault_trace_rankit
    !!end fichier save_fault_trace 2d

    !!fichier save_vorticity 2d
    subroutine semname_save_vorticity_epsidat (sorties,it,rank,fnamef)
        !SEMFILE 73 W XXX/epsi_point_III.dat.JJJ (MKA)
        implicit none
        character(Len=*),intent(in) :: sorties
        integer,intent(in) :: it
        integer,intent(in) :: rank
        character(Len=20) :: Ait
        !character(Len=20) :: crank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(Ait,"(I20)") it
        !write(crank,"(I20)") rank
        write(fnamef,"(a,a,a,a,I4.4)")trim(adjustl(sorties)),"/epsi_point_",trim(adjustl(Ait)),".dat.",rank

        DEBUG(fnamef)
    end subroutine semname_save_vorticity_epsidat

    subroutine semname_save_vorticity_curldat (sorties,it,rank,fnamef)
        !SEMFILE 74 W XXX/epsi_point_III.dat.JJJ (MKA)
        implicit none
        character(Len=*),intent(in) :: sorties
        integer,intent(in) :: it
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=20) :: Ait
        !character(Len=20) :: crank
        write(Ait,"(I20)") it
        !write(crank,"(I20)") rank
        write(fnamef,"(a,a,a,a,I4.4)")trim(adjustl(sorties)),"/curl_v_",trim(adjustl(Ait)),".dat.",rank

        DEBUG(fnamef)
    end subroutine semname_save_vorticity_curldat

    subroutine semname_save_vorticity_epsilon (rank,it,fnamef)
        !SEMFILE 73 W epsi_pointIII.JJJ (NOMKA)
        implicit none
        integer ,intent(in) :: rank
        integer ,intent(in) :: it
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,I4.4,a1,I5.5)")"epsi_point",rank, ".", it

        DEBUG(fnamef)
    end subroutine semname_save_vorticity_epsilon

    subroutine semname_save_vorticity_curl (rank,it,fnamef)
        !SEMFILE 74 W curl_vIII.JJJ (NOMKA)
        implicit none
        integer ,intent(in) :: rank
        integer ,intent(in) :: it
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,I4.4,a1,I5.5)")"curl_v",rank, ".", it

        DEBUG(fnamef)
    end subroutine semname_save_vorticity_curl
    !!end fichier save_vorticity 2d

    !!fichier savefield 2d 3d
    subroutine semname_savefield_dat (sorties,grandeur,rank,it,fnamef) !2d
        !SEMFILE 61 W XXX/YYY_III.dat.JJJ (MKA)
        implicit none
        integer,intent(in) :: rank
        integer,intent(in) :: it
        character(Len=20) :: cit
        character(Len=*),intent(in) :: sorties
        character(Len=*),intent(in) :: grandeur
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(cit,"(I4)") it
        !character(Len=20) :: crank
        !write(crank,"(I20)")rank
        write(fnamef,"(a,a,a,a,a,a,I4.4)")trim(adjustl(sorties)),"/",trim(adjustl(grandeur)),"_",trim(adjustl(cit)),".dat.",rank
        DEBUG(fnamef)
    end subroutine semname_savefield_dat

    subroutine semname_savefield_fieldt (rank,it,fnamef)!2d
        !SEMFILE 61 W III.JJJ (NOMKA)
        implicit none
        integer ,intent(in) :: rank
        integer ,intent(in) :: it
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(I4.4,a,I5.5)")rank,".",it

        DEBUG(fnamef)
    end subroutine semname_savefield_fieldt

    subroutine semname_savefield_results (it,rank,fnamef)
        implicit none
        integer ,intent(in) :: rank
        integer ,intent(in) :: it
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp1, temp2

        call semname_snap_result_dir(it, temp1)
        write(temp2,"(a,I4.4,a,I4.4)") "vel_",it,".dat.",rank

        fnamef = pjoin(temp1,temp2)
    end subroutine semname_savefield_results

    subroutine semname_savefield_datfields (rank,field,count,fnamef) !3d savefield_disp aussi
        !SEMFILE 61->63 W data/SField/ProcIIIfield(x~y~z)JJJ (NOMKA)
        implicit none
        integer ,intent(in) :: rank
        integer ,intent(in) :: count
        character(Len=*),intent(in) :: field
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,I2.2,a,I3.3)")"data/SField/Proc",rank,trim(adjustl(field)),count

        DEBUG(fnamef)
    end subroutine semname_savefield_datfields
    !!end fichier savefield 2d 3d

    !!fichier savefield_disp 3d
    subroutine semname_savefield_disp_datsorties (rank,it,fnamef) !3d
        implicit none
        integer ,intent(in) :: rank
        integer ,intent(in) :: it
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp1, temp2

        call semname_snap_result_dir(it, temp1)
        write(temp2,"(a,I4.4,a,I4.4)") "displ_",it,".dat.",rank

        fnamef = pjoin(temp1,temp2)
    end subroutine semname_savefield_disp_datsorties

    !!end fichier savefield_disp 3d

    !!fichier savetrace 3d
    subroutine semname_savetrace_trace (n,coord,k,fnamef) !3d
        !SEMFILE 61->63 W ./data/Trace(X~Y~Z)/traceIII(X~Y~Z)JJJ
        implicit none
        integer ,intent(in) :: n
        integer ,intent(in) :: k
        character(Len=*),intent(in) :: coord
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,a,a,I4.4,a,I3.3)")"./data/Trace",trim(adjustl(coord)),"/trace",n,trim(adjustl(coord)),k

        DEBUG(fnamef)
    end subroutine semname_savetrace_trace

    !!end fichier savetrace 3d

    !!fichier sem 2d 3d
    subroutine semname_sem_frontiere(fichier,fnamef) !2d
        !SEMFILE ?(sem%fileId) W XXXfrontiere
        implicit none
        character(Len=*),intent(in) :: fichier
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,a)")trim(adjustl(fichier)),"frontiere"

        DEBUG(fnamef)
    end subroutine semname_sem_frontiere

    subroutine semname_sem_fichier (fichier,fnamef)
        !SEMFILE ? (sem%fileId) W XXX
        implicit none
        character(Len=*),intent(in) :: fichier
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a)")trim(adjustl(fichier))

        DEBUG(fnamef)
    end subroutine semname_sem_fichier
    !!end fichier sem 2d 3d

    !!fichier shape4 2d shape8 2d 3d
    subroutine semname_posptg (rank,fnamef)
        !SEMFILE 88 W ./Resultats/posPTG.III (MKA)
        implicit none
        integer,intent(in) :: rank
        !character(Len=20) :: crank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        !write(crank,"(I20)") rank
        write(fnamef,"(a,I4.4)")"./Resultats/posPtG.",rank

        DEBUG(fnamef)
    end subroutine semname_posptg

    subroutine semname_post_data (rank,fnamef)
        !SEMFILE 88 W ./Resultats/posPTG.III (MKA)
        implicit none
        integer,intent(in) :: rank
        !character(Len=20) :: crank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        !write(crank,"(I20)") rank
        write(fnamef,"(a,I4.4,a)")"./Resultats/post_data.",rank,".h5"

        DEBUG(fnamef)
    end subroutine semname_post_data

    subroutine semname_connecptg (rank,fnamef)
        !SEMFILE 89 W ./Resultats/connecPTG.III (MKA)
        implicit none
        integer,intent(in) :: rank
        !character(Len=20) :: crank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        !write(crank,"(I20)") rank
        write(fnamef,"(a,I4.4)")"./Resultats/connecPtG.",rank

        DEBUG(fnamef)
    end subroutine semname_connecptg

    subroutine semname_shape8_norm (rank,fnamef) !3d
        !SEMFILE 60 W PWNormalIII
        implicit none
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,I2.2)")"PWNormale",rank

        DEBUG(fnamef)
    end subroutine semname_shape8_norm

    subroutine semname_shape8_coord (rank,fnamef) !3d
        !SEMFILE 62 W CoordPWIII
        implicit none
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,I2.2)")"CoordPW",rank

        DEBUG(fnamef)
    end subroutine semname_shape8_coord

    subroutine semname_shape8_neu (rank,fnamef) !3d
        !SEMFILE 60 W NeuNormaleIII
        implicit none
        integer,intent(in) :: rank
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a,I2.2)")"NeuNormale",rank

        DEBUG(fnamef)
    end subroutine semname_shape8_neu
    !!end fichier shape4 2d shape8 2d 3d

    !!fichier drive_sem 3d
    ! semname_couplage_dat

    subroutine semname_drive_sem_listing (rg, fnamef)
        !SEMFILE 50 ? sem.listingIII (COUPLAGE)
        implicit none
        integer, intent(in) :: rg
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=2) :: sit
        !! fichier d'ecriture de la sortie
        if(rg.LE.9) then
            write(sit,'(I1)') rg
            sit='0'//sit
        elseif(rg.LE.99) then
            write(sit,'(I2)') rg
        endif

        write(fnamef,"(a,a)")"./listings/sem.listing.",sit

        DEBUG(fnamef)
    end subroutine semname_drive_sem_listing

    subroutine semname_drive_sem_resulttemp (fnamef)
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef

        fnamef = pjoin(path_results, "temps_sem.dat")
    end subroutine semname_drive_sem_resulttemp

    subroutine semname_xdmf(isort, fnamef)
        implicit none
        integer, intent(in) :: isort
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        character(Len=MAX_FILE_SIZE) :: temp
        write(temp,"(a,I4.4,a)") "mesh.",isort,".xmf"
        fnamef = pjoin(path_results, temp)
    end subroutine semname_xdmf

    !! Nom du fichier contenant le nombre de processeurs ayant genere une sortie
    subroutine semname_nb_proc(isort,fnamef)
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        integer, intent(in) :: isort
        character(Len=MAX_FILE_SIZE) :: temp
        call semname_snap_result_dir(isort, temp)
        fnamef = pjoin(temp, "Nb_proc")
    end subroutine semname_nb_proc


    !!fichier ondelette 3d
    subroutine semname_ondelette_sgn (fnamef)
        !SEMFILE 60 W sgn_src
        implicit none
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a)")"sgn_src"

        DEBUG(fnamef)
    end subroutine semname_ondelette_sgn
    !!end fichier ondelette 3d

    !!fichier unv 2d 3d
    subroutine semname_unv_fichier (fichier,fnamef)
        !SEMFILE ? sunv%fileId R XXX
        implicit none
        character(Len=*),intent(in) :: fichier
        character(Len=MAX_FILE_SIZE),intent(out) :: fnamef
        write(fnamef,"(a)")trim(adjustl(fichier))

        DEBUG(fnamef)
    end subroutine semname_unv_fichier
    !!end fichier unv 2d 3d
end module semdatafiles
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
