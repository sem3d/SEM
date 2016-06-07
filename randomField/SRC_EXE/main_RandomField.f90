program main_RandomField

    use mpi
    use constants_RF
    use randomFieldND
    use mesh_RF
    use writeResultFile_RF
    use displayCarvalhol
    use charFunctions
    use write_Log_File
    use systemUt_RF
    use common_variables_RF
    use type_RF
    use type_MESH
    use type_inputRF
    use calls_RF
    use sample_RF

    implicit none

    !INPUTS

    !LOCAL VARIABLES
    character(len=buf_RF) :: path
    double precision, dimension(9) :: times, all_times
    integer :: code, i

    !INPUT VARIABLES
    type(IPT_RF)  :: IPT_Temp !Only for sake of practicity
    type(IPT_RF)  :: IPT !The one that shoud be initialized when calling from an external program

    !call MPI_INIT(code)
    !write(*,*) "HERE"
    !call MPI_FINALIZE(code)

    !Initializing MPI
    call init_communication(MPI_COMM_WORLD, IPT_Temp%comm, IPT_Temp%rang, IPT_Temp%nb_procs)

    !Options
    IPT_Temp%writeDataSet = .true.
    IPT_Temp%writeUNVinterpolation = .true.
    IPT_Temp%sameFolder = .true.
    IPT_Temp%outputStyle = 1 !1: parallel hdf5, 2: hdf5 per proc
    IPT_Temp%write_intermediate_files = .false.
    IPT_Temp%sampleFields = .true.

    times(1) = MPI_Wtime() !Initial Time

    if(IPT_Temp%rang == 0)then
        write(*,*)
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "************************                             ***********************"
        write(*,*) "************************  RANDOM FIELD LIBRARY TEST  ***********************"
        write(*,*) "************************                             ***********************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*) "****************************************************************************"
        write(*,*)
    end if

    if(IPT_Temp%rang == 0) write(*,*) "-> MPI_communications started"
    if(IPT_Temp%rang == 0) write(*,*) "         running on: "
    if(IPT_Temp%rang == 0) call system("pwd")
    if(IPT_Temp%rang == 0) write(*,*) "         nb_procs    = ", IPT_Temp%nb_procs
    if(IPT_Temp%rang == 0) write(*,*) "         outputStyle = ", IPT_Temp%outputStyle

    !Initializing folders
    if(IPT_Temp%rang == 0) write(*,*)  "-> Initialize Folders"
    call init_basic_folders(IPT_Temp%comm, IPT_Temp)

    !READING INPUTS--------------------------------------------------
    !----------------------------------------------------------------
    if(IPT_Temp%rang == 0) write(*,*)  "-> Reading inputs"
    call wLog("-> Reading inputs")
    !Reading Main---------------------------------------------------
    if(IPT_Temp%rang == 0) write(*,*)  "     -> Reading Main Input"
    call wLog("     -> Reading Main Input")
    call read_main_input("./RF_main_input", IPT_Temp)

    call MPI_BARRIER(IPT_Temp%comm, code)

    if(IPT_Temp%application /= 1) then
        if(IPT_Temp%rang == 0) write(*,*)  "     SEM generation"
        call read_main_input("./TEMP_RF_main_input", IPT_Temp)
        call MPI_BARRIER(IPT_Temp%comm, code)
        if(IPT_Temp%rang == 0) call system("mv TEMP_RF_main_input "//IPT_Temp%appFolder)
    end if

    !Initial allocation---------------------------------------------
    call allocate_init()

    do i = 1, IPT_Temp%nSamples

        !Reading Mesh---------------------------------------------------
        if(IPT_Temp%rang == 0) write(*,*)  "     -> Reading Mesh Input"
        call wLog("     -> Reading Mesh Input")
        path = IPT_Temp%mesh_inputs(i)
        path = adjustL(path)
        !call wLog("        file: "//trim(path))
        call read_mesh_input(path, IPT_Temp)
        !Reading Generation Input---------------------------------------
        if(IPT_Temp%rang == 0) write(*,*)  "     -> Reading Generation Input"
        call wLog("     -> Reading Generation Input")
        path = IPT_Temp%gen_inputs(i)
        path = adjustL(path)
        !call wLog("        file: "//trim(path))
        call read_generation_input(path, IPT_Temp)

        !Estimating ideal number of fields
        call estimate_nFields(IPT_Temp)

        !Validating Inputs----------------------------------------------
        if(IPT_Temp%rang == 0) write(*,*)  "     -> Validating Input (IPT_Temp)"
        call wLog("    Validating Inputs (IPT_Temp)")
        call validate_input(IPT_Temp)
        if(IPT_Temp%rang == 0) call show_IPT_RF(IPT_Temp, "IPT_Temp")
#ifdef MAKELOG
        call show_IPT_RF(IPT_Temp, forLog_in=.true.)
#endif

        times(2) = MPI_Wtime() !Reading Inputs

        !Initialize Inputs
        if(IPT_Temp%rang == 0) write(*,*)  " "
        if(IPT_Temp%rang == 0) write(*,*)  " "
        if(IPT_Temp%rang == 0) write(*,*)  "-> Initializing Input (IPT)"

        if(.true.)then
            call init_IPT_RF(&
                IPT, &
                log_ID = IPT_Temp%log_ID, &
                comm = IPT_Temp%comm, &
                rang = IPT_Temp%rang, &
                nb_procs = IPT_Temp%nb_procs, &
                nDim = IPT_Temp%nDim_gen, &
                meshMod = IPT_Temp%meshMod, &
                xMinGlob_in = IPT_Temp%xMinGlob_in, &
                xMaxGlob_in = IPT_Temp%xMaxGlob_in, &
                pointsPerCorrL = IPT_Temp%pointsPerCorrL, &
                procPerDim = IPT_Temp%procPerDim, &
                fieldAvg = IPT_Temp%fieldAvg, &
                fieldVar = IPT_Temp%fieldVar, &
                corrL_in = IPT_Temp%corrL_in, &
                overlap_in = IPT_Temp%overlap_in, &
                corrMod = IPT_Temp%corrMod, &
                margiFirst = IPT_Temp%margiFirst, &
                method = IPT_Temp%method, &
                Nmc = IPT_Temp%Nmc, &
                seedStart = IPT_Temp%seedStart, &
                nFields = IPT_Temp%nFields, &
                localizationLevel = IPT_Temp%localizationLevel, &
                writeDataSet = IPT_Temp%writeDataSet, &
                sameFolder = IPT_Temp%sameFolder, &
                outputStyle = IPT_Temp%outputStyle, &
                write_intermediate_files = IPT_Temp%write_intermediate_files, &
                sampleFields = IPT_Temp%sampleFields, &
                writeUNVinterpolation = IPT_Temp%writeUNVinterpolation, &
                outputFolder = IPT_Temp%out_folders(i), &
                outputName = IPT_Temp%out_names(i), &
                unv = IPT_Temp%unv, &
                unv_path = IPT_Temp%unv_path, &
                monotype = IPT_Temp%monotype, &
                coordList_local = IPT_Temp%coordList_local, &
                connectList_local = IPT_Temp%connectList_local)

        else
            call init_IPT_RF_std(&
                             IPT, &
                             comm = IPT_Temp%comm, &
                             nDim = IPT_Temp%nDim_gen, &
                             xMinGlob_in = IPT_Temp%xMinGlob_in, &
                             xMaxGlob_in = IPT_Temp%xMaxGlob_in, &
                             fieldAvg = IPT_Temp%fieldAvg, &
                             fieldVar = IPT_Temp%fieldVar, &
                             corrL_in = IPT_Temp%corrL_in, &
                             corrMod = IPT_Temp%corrMod, &
                             margiFirst = IPT_Temp%margiFirst, &
                             seedStart = IPT_Temp%seedStart, &
                             outputFolder = IPT_Temp%out_folders(i), &
                             outputName = IPT_Temp%out_names(i))
        end if

        !Generating random fields
        call make_random_field(IPT, times)

        call MPI_ALLREDUCE (times, all_times, size(times), MPI_DOUBLE_PRECISION, MPI_SUM, IPT%comm,code)

        if(IPT_Temp%rang == 0) write(*,*) ""
        if(IPT_Temp%rang == 0) write(*,*) ""
        if(IPT_Temp%rang == 0) write(*,*) "AVERAGE TIMES (WALL)------------------------ "
        if(IPT_Temp%rang == 0) write(*,*) "Reading Inputs        = ", (all_times(2) - all_times(1))/dble(IPT_Temp%nb_procs)
        if(IPT_Temp%rang == 0) write(*,*) "Pre Organization      = ", (all_times(3) - all_times(2))/dble(IPT_Temp%nb_procs)
        if(IPT_Temp%rang == 0) write(*,*) "Generation            = ", (all_times(4) - all_times(3))/dble(IPT_Temp%nb_procs)
        !if(IPT_Temp%rang == 0) write(*,*) "Localization Int      = ", (all_times(5) - all_times(4))/dble(IPT_Temp%nb_procs)
        if(IPT_Temp%rang == 0) write(*,*) "Localization Ext      = ", (all_times(5) - all_times(4))/dble(IPT_Temp%nb_procs)
        if(IPT_Temp%rang == 0) write(*,*) "Writing Normalization = ", (all_times(6) - all_times(5))/dble(IPT_Temp%nb_procs)
        if(IPT_Temp%rang == 0) write(*,*) "Transforming          = ", (all_times(7) - all_times(6))/dble(IPT_Temp%nb_procs)
        if(IPT_Temp%rang == 0) write(*,*) "Writing Files         = ", (all_times(8) - all_times(7))/dble(IPT_Temp%nb_procs)
        if(IPT_Temp%rang == 0) write(*,*) "UNV Interpolation     = ", (all_times(9) - all_times(8))/dble(IPT_Temp%nb_procs)
        if(IPT_Temp%rang == 0) write(*,*) ""

        !4 - Sampling
        !5 - Internal Localization
        !6 - External Localization
        !7 - Normalizing
        !8 - Transformation
        !9 - Writing Files
        !10 - UNV Interpolation

        !Deallocating
        call deallocate_all(IPT)

        if(IPT_Temp%rang == 0) then
            write(*,*) ""
            write(*,*) "---------------------------------------------------------------------";
            write(*,*) "-----------------END RANDOM FIELD LIBRARY TEST-----------------------";
            write(*,*) "---------------------------------------------------------------------";
            write(*,*) ""
        end if

        call MPI_BARRIER(IPT_Temp%comm, code)
    end do

    call finalize_IPT_RF(IPT_Temp)
    !Finalizing MPI
    call end_communication()
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
        !----------------------------------------------------------------------------------------------------
    contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_communication(comm_local, comm, rang, nb_procs)
            implicit none
            !INPUT
            integer, intent(in) :: comm_local
            !OUTPUT
            integer, intent(out) :: comm, rang, nb_procs
            !LOCAL
            integer :: code

            !call MPI_Init_thread(MPI_THREAD_MULTIPLE, &provided
            call MPI_INIT(code)
            call MPI_COMM_RANK(comm_local, rang, code)
            call MPI_COMM_SIZE(comm_local, nb_procs, code)

            comm = comm_local

        end subroutine init_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_basic_folders(comm, IPT_Temp)
            implicit none
            !INPUT
            integer, intent(in) :: comm
            type(IPT_RF), intent(inout)  :: IPT_Temp !Only for sake ok practicity
            !LOCAL
            integer, dimension(8) :: date_time
            integer :: code
            character(len=10), dimension(3) :: strings
            character(len=buf_RF) :: results_folder_name
            !    !LOCAL VARIABLES
            character(len=buf_RF) :: logFilePath, log_folder_name
            logical :: folderExist1


            if(IPT_Temp%sameFolder) then
                results_folder_name = "res"
            else
                !date_time_label
                if(IPT_Temp%rang == 0) then
                    folderExist1 = .true.
                    do while(folderExist1)
                        call date_and_time(strings(1), strings(2), strings(3), date_time)
                        results_folder_name = strings(1)(3:8)//"_"//strings(2)(1:6)//"_res"
                        folderExist1 = folderExist(results_folder_name, results_path)
                        if(folderExist1) write(*,*) "Repeated result folder, changing name"
                    end do
                end if
                call MPI_BCAST (results_folder_name, len(results_folder_name), &
                                MPI_CHARACTER, 0, comm, code)
            end if

            if(IPT_Temp%rang == 0) write(*,*) "-> Setting folder path"
            IPT_Temp%outputFolder = string_join_many(results_path,"/",results_folder_name)
            if(IPT_Temp%rang == 0) write(*,*) "     single_path = "//trim(IPT_Temp%outputFolder)

#ifdef MAKELOG
            if(IPT_Temp%rang == 0) write(*,*) "IFDEF MAKELOG DEFINED"

            log_folder_name     = trim(adjustL(results_folder_name))//"/log"
            if(IPT_Temp%sameFolder) then
                log_folder_name     = "."
                logFilePath = trim(string_join_many("./",log_filename))
            else
                call create_folder(log_folder_name, results_path, IPT_Temp%rang, comm)
                logFilePath = trim(adjustL(&
                                  string_join_many(results_path,"/",log_folder_name,"/",log_filename)))
            end if

            !Initializing logFiles
            if(IPT_Temp%rang == 0) write(*,*)  "-> Initialize logFiles"
            if(IPT_Temp%rang == 0) write(*,*)  " logFilePath = ", trim(adjustL(logFilePath)), "<RANK>"
            call init_log_file(trim(adjustL(logFilePath)), IPT_Temp%rang, IPT_Temp%log_ID, IPT_Temp%nb_procs)
#else
            if(IPT_Temp%rang == 0) write(*,*) "IFDEF MAKELOG NOT DEFINED"
#endif

            !create xmf and h5 folders
            !path = string_join_many(results_path,"/",results_folder_name)
            !call create_folder("xmf", path, IPT_Temp%rang, comm)
            !call create_folder("h5", path, IPT_Temp%rang, comm)

        end subroutine init_basic_folders


        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine allocate_init()

        end subroutine allocate_init


        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine end_communication()
            implicit none
            !LOCAL
            integer :: code

            !call finalize_log_file()
            call MPI_FINALIZE(code)

        end subroutine end_communication

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine deallocate_all(IPT)
            implicit none
            !LOCAL
            type(IPT_RF), intent(inout)  :: IPT

            call finalize_IPT_RF(IPT)

        end subroutine deallocate_all

end program main_RandomField
