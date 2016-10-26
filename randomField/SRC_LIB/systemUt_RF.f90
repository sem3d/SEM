module systemUt_RF
    !use mpi
    use write_Log_File

contains
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine create_folder(folder, path, rang, comm, singleProc)

        implicit none
        !INPUT
        character(len = *), intent(in) :: folder, path
        integer, intent(in) :: rang, comm
        logical, intent(in), optional :: singleProc
        !LOCAL
        character(len = 200) command, fullName
        integer :: code
        logical :: effectSP

        effectSP = .false.
        if(present(singleProc)) effectSP = singleProc

        !if(rang==0) write(*,*) "creating Folder"

        !if(.not. folderExist (folder, path, compiler)) then
            fullName = trim(adjustL(path)) // "/" // trim(adjustL(folder))
            !write(*,*) "fullName = ", fullName
            !write(*,*) "Directory is being created: ", fullName
            if(rang==0) then
                !write(*,*) "fullName = ", fullName
                !call system("ls "//trim(adjustL(path)))
                command = 'mkdir -pv '// trim(adjustL(fullName))
                !write(*,*) "command = ", command
                call system(command)
            end if
        !end if

        if(.not. effectSP) call MPI_BARRIER (comm ,code)

    end subroutine create_folder

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine delete_folder(folder, path, comm)
        implicit none
        !INPUT
        character(len = *), intent(in) :: folder, path
        integer, intent(in) :: comm
        !LOCAL
        character(len = 200) command, fullName

        integer :: code

        fullName = trim(adjustL(path)) // "/" // trim(adjustL(folder))

        !if(folderExist (folder, path, compiler)) then
        write (*,*) "WARNING!!! Directory: '", trim(fullName) ,"' will be deleted"
        command = 'rm -r '// trim(adjustL(fullName))
        call system(command)
        !end if

        call MPI_BARRIER (comm ,code)

    end subroutine delete_folder

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine rename_folder(folder, path, newName, comm)
        implicit none
        !INPUT
        character(len = *), intent(in) :: folder, path, newName
        integer, intent(in) :: comm
        !LOCAL
        character(len = 200) command, fullName_new, fullName_old
        integer :: code

        fullName_old = trim(adjustL(path)) // "/" // trim(adjustL(folder))
        fullName_new = trim(adjustL(path)) // "/" // trim(adjustL(newName))

        !if(folderExist (folder, path, compiler)) then
        write (*,*) "WARNING!!! Directory: '", trim(fullName_old) ,"' will be renamed to: ", trim(fullName_new)
        command = 'mv '// trim(adjustL(fullName_old))//" "//trim(adjustL(fullName_new))
        call system(command)
        !end if

        call MPI_BARRIER (comm ,code)

    end subroutine rename_folder

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function folderExist (folder, path) result (dirExists)
        implicit none
        !INPUT
        character(len = *), intent(in) :: folder, path
        !OUTPUT
        logical :: dirExists
        !LOCAL
        character(len=1024) :: filePath
        integer :: fileId=56, error

        filePath = trim(adjustL(path))//"/"//trim(adjustL(folder))//"/INQUIRE_file"
        !write(*,*) "filePath = ", trim(filePath)

        dirExists = .false.
        open (unit = fileId , file = filePath, action = 'write', iostat=error)
        !write(*,*) "error = ", error
        if(error==0) then
            close(fileId)
            dirExists = .true.
            call system("rm "//filePath)
        end if


    end function folderExist

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    function fileExist (filePath) result (fileExists)
        implicit none
        !INPUT
        character(len = *), intent(in) :: filePath
        !OUTPUT
        logical :: fileExists
        !LOCAL
        integer :: fileId=56, error

        fileExists = .false.
        open (unit = fileId , file = filePath, action = 'read', iostat=error)
        !write(*,*) "error = ", error
        if(error==0) then
            close(fileId)
            fileExists = .true.
        end if


    end function fileExist

end module systemUt_RF


!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
