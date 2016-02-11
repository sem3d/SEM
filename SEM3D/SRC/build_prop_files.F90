module build_prop_files

    use sdomain
    use mpi
    use constants
    use scomm
    use scommutils
    use randomFieldND
    use displayCarvalhol
    use writeResultFile_RF
    use define_random

    implicit none
#include "index.h"

    character(len=15) :: procFileName = "prop"
    character(len=50) :: h5folder  = "./prop/h5", &
        XMFfolder = "./prop", &
        h5_to_xmf = "./h5"

contains
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine create_prop_files(Tdomain, rg)

        !INPUTS
        type (domain), intent (INOUT), target :: Tdomain
        integer      , intent(IN) :: rg

        !LOCAL
        integer :: mat, assocMat
        integer :: randMethod = 1 !1 for Isotropic method, 2 for Shinozuka's
        real               , dimension(:)   , allocatable :: avgProp;
        integer            , dimension(:)   , allocatable :: nSubDPoints;
        double precision   , dimension(:, :), allocatable :: prop !Properties
        character(len=110) , dimension(:)   , allocatable :: HDF5NameList

        allocate(avgProp(0:nProp-1)) !Obs: should have nProp declared before
        allocate(prop(0:size(Tdomain%GlobCoord,2)-1, 0:nProp-1)) !Subdomain properties Matrix ((:,0) = Dens, (:,1) = Lambda, (:,2) = Mu) per proc
        allocate(HDF5NameList(0:Tdomain%n_mat-1))
        allocate(nSubDPoints(0:Tdomain%n_mat-1))
        HDF5NameList(:) = "not_Used"

        !Writing hdf5 files
        if(rg == 0) write(*,*)
        do mat = 0, Tdomain%n_mat - 1
            if(propOnFile(Tdomain, mat)) then
                if(rg == 0) write(*,*) "  Material ", mat, " will have properties on file"
                assocMat = Tdomain%sSubdomain(mat)%assocMat
                avgProp  = [Tdomain%sSubDomain(mat)%Ddensity, &
                    Tdomain%sSubDomain(mat)%DLambda,  &
                    Tdomain%sSubDomain(mat)%DMu]

                if(Tdomain%sSubDomain(assocMat)%material_type == "S".or. &
                    Tdomain%sSubDomain(assocMat)%material_type == "P") then
                    prop(:,0) = avgProp(0)
                    prop(:,1) = avgProp(1)
                    prop(:,2) = avgProp(2)
                else if(Tdomain%sSubDomain(assocMat)%material_type == "R") then
                    if(rg == 0) write(*,*) "  Generating Random Properties"
                    call build_random_properties(Tdomain, rg, mat, prop, randMethod)
                end if
                if(rg == 0) write(*,*) "  Writing hdf5 files"
                call write_ResultHDF5Unstruct_MPI(Tdomain%GlobCoord, prop, trim(procFileName)//"_read", &
                    rg, trim(h5folder), Tdomain%communicateur,   &
                    ["_proc", "_subD"], [rg, mat], HDF5NameList(mat))
            end if
        end do

        !Writing XMF File
        if(rg == 0) write(*,*) "  Writing XMF file"
        nSubDPoints(:) = size(Tdomain%GlobCoord,2)
        call writeXMF_RF_MPI(nProp, HDF5NameList, nSubDPoints, Tdomain%subD_exist, Tdomain%n_dime, &
            trim(string_join(procFileName,"-TO_READ")), rg, trim(XMFfolder),     &
            Tdomain%communicateur, trim(h5_to_xmf),                               &
            ["Density","Lambda ","Mu     "])

        !Deallocating
        if(allocated(prop))         deallocate(prop)
        if(allocated(HDF5NameList)) deallocate(HDF5NameList)
        if(allocated(nSubDPoints))  deallocate(nSubDPoints)
        if(allocated(avgProp))      deallocate(avgProp)

        !Deallocating arrays associated to random properties
        do mat = 0, Tdomain%n_mat - 1
            if (allocated(Tdomain%sSubdomain(mat)%varProp))       deallocate(Tdomain%sSubdomain(mat)%varProp)
            if (allocated(Tdomain%sSubDomain(mat)%corrL))         deallocate(Tdomain%sSubDomain(mat)%corrL)
            if (allocated(Tdomain%sSubDomain(mat)%margiFirst))    deallocate(Tdomain%sSubDomain(mat)%margiFirst)
            if (allocated(Tdomain%sSubDomain(mat)%chosenSeed))    deallocate(Tdomain%sSubDomain(mat)%chosenSeed)
        end do

    end subroutine create_prop_files

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine apply_prop_files(Tdomain, rg)
        !INPUTS
        type (domain), intent (INOUT), target :: Tdomain
        integer      , intent(IN) :: rg

        !LOCAL
        integer :: mat, n, ipoint, i, j, k, lnum
        integer :: ngllx,nglly,ngllz
        double precision, dimension(:, :, :), allocatable :: propMatrix !Properties
        allocate(propMatrix(0:size(Tdomain%GlobCoord,2)-1, 0:nProp-1, 0:Tdomain%n_mat-1))


        !Reading files and putting it on a matrix
        do mat = 0, Tdomain%n_mat-1
            if(propOnFile(Tdomain, mat)) then

                call read_properties_from_file(rg, propMatrix(:,:, mat), &
                    trim(procFileName)//"_read", trim(h5folder),     &
                    ["_proc", "_subD"], [rg, mat])
            end if
        end do

        !Applying properties to elements
        do n = 0, Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            lnum = Tdomain%specel(n)%lnum
            if(propOnFile(Tdomain, mat)) then

                ngllx    = Tdomain%sSubDomain(mat)%NGLLx
                nglly    = Tdomain%sSubDomain(mat)%NGLLy
                ngllz    = Tdomain%sSubDomain(mat)%NGLLz

                do i = 0, ngllx-1
                    do j = 0, nglly-1
                        do k = 0, ngllz-1
                            ipoint  = Tdomain%specel(n)%Iglobnum(i,j,k)
                            select case (Tdomain%specel(n)%domain)
                                case (DM_SOLID)
                                    Tdomain%sdom%Density_(i,j,k,lnum) = propMatrix(ipoint, 0, mat)
                                    Tdomain%sdom%Lambda_ (i,j,k,lnum) = propMatrix(ipoint, 1, mat)
                                    Tdomain%sdom%Mu_     (i,j,k,lnum) = propMatrix(ipoint, 2, mat)
                                case (DM_FLUID)
                                    Tdomain%fdom%Density_(i,j,k,lnum) = propMatrix(ipoint, 0, mat)
                                    Tdomain%fdom%Lambda_ (i,j,k,lnum) = propMatrix(ipoint, 1, mat)
                                case (DM_SOLID_PML)
                                    Tdomain%spmldom%Density_(i,j,k,lnum) = propMatrix(ipoint, 0, mat)
                                    Tdomain%spmldom%Lambda_ (i,j,k,lnum) = propMatrix(ipoint, 1, mat)
                                    Tdomain%spmldom%Mu_     (i,j,k,lnum) = propMatrix(ipoint, 2, mat)
                                case (DM_FLUID_PML)
                                    Tdomain%fpmldom%Density_(i,j,k,lnum) = propMatrix(ipoint, 0, mat)
                                    Tdomain%fpmldom%Lambda_ (i,j,k,lnum) = propMatrix(ipoint, 1, mat)
                                    Tdomain%fpmldom%Mu_     (i,j,k,lnum) = propMatrix(ipoint, 2, mat)
                            end select
                        end do
                    end do
                end do

            end if
        end do

        deallocate(propMatrix)

    end subroutine apply_prop_files

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine read_properties_from_file(rg, prop,  &
        fileName, folderPath, labels, indexes)

        use sem_hdf5
        use hdf5

        implicit none
        !INPUT
        integer         , intent(in)            :: rg
        character(len=*)                 , intent(in) :: filename;
        character(len=*)                 , intent(in) :: folderPath
        character(len=*) , dimension(1:), optional  , intent(in) :: labels
        integer          , dimension(1:), optional  , intent(in) :: indexes
        !OUTPUT
        real        , intent(out), dimension(0:, 0:) :: prop !Properties
        !LOCAL
        character (len=12) :: rangStr;
        character(len=110) :: fileHDF5Name, fullPath !File name
        integer            :: i, error
        integer(HID_T)     :: file_id
        real,    allocatable, dimension(:,:) :: rtemp2

        if(.not. present(labels)) then
            write(rangStr,'(I100)'  ) rg
            rangStr      = adjustL(rangStr)
            fileHDF5Name = trim(fileName)//"-proc"//trim(rangStr)//".h5"
        else
            fileHDF5Name = fileName
            do i = 1, size(labels)
                fileHDF5Name =  string_join(fileHDF5Name,stringNumb_join(labels(i), indexes(i)))
            end do
        end if

        fileHDF5Name = string_join(fileHDF5Name,".h5")
        fullPath     = string_join(folderPath,"/"//fileHDF5Name)

        !write(*,*) "fileHDF5Name =",fileHDF5Name
        !write(*,*) "fullPath =",fullPath

        file_id = 12
        call h5open_f(error) ! Initialize FORTRAN interface.
        call h5fopen_f(trim(fullPath), H5F_ACC_RDONLY_F, file_id, error)
        do i = 0, size(prop, 2)-1
            call read_dataset(file_id, trim(stringNumb_join("RF_",i+1)), rtemp2)
            prop(:,i:i) = rtemp2(:,:)

            deallocate(rtemp2)
        end do
        call h5fclose_f(file_id, error) ! Close the file.
        call h5close_f(error) ! Close FORTRAN interface.

    end subroutine read_properties_from_file

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine create_prop_visu_files (Tdomain, rg)

        !INPUTS
        type (domain), intent (INOUT), target :: Tdomain
        integer      , intent(IN) :: rg

        !LOCAL
        integer :: n, ipoint, i, j, k, mat, lnum
        integer :: ngllx,nglly,ngllz
        logical         , dimension(:,:) , allocatable :: globCoordMask, propMask
        double precision, dimension(:, :), allocatable :: prop !Properties
        integer         , dimension(:)   , allocatable :: nSubDPoints;
        character(len=110) , dimension(:), allocatable :: HDF5NameList

        allocate(prop(0:size(Tdomain%GlobCoord,2)-1, 0:nProp-1)) !Subdomain properties Matrix ((:,0) = Dens, (:,1) = Lambda, (:,2) = Mu) per proc

        !Creating  visualization files (by subdomain)
        allocate(HDF5NameList(0:Tdomain%n_mat - 1))
        allocate(nSubDPoints(0:Tdomain%n_mat - 1))
        allocate(globCoordMask(0:size(Tdomain%GlobCoord,1)-1, 0:size(Tdomain%GlobCoord,2)-1))
        allocate(propMask(0:size(Tdomain%GlobCoord,2)-1, 0:nProp-1))

        HDF5NameList(:) = "not_Used"
        nSubDPoints(:)  = 0

        if(rg == 0) write(*,*) ">>>>Writing visualization hdf5 files"
        do mat = 0, Tdomain%n_mat - 1
            globCoordMask(:,:) = .false.
            propMask(:,:)      = .false.
            ngllx = Tdomain%sSubDomain(mat)%NGLLx
            nglly = Tdomain%sSubDomain(mat)%NGLLy
            ngllz = Tdomain%sSubDomain(mat)%NGLLz

            do n = 0, Tdomain%n_elem-1
                lnum = Tdomain%specel(n)%lnum
                if(Tdomain%specel(n)%mat_index == mat) then
                    !Building masks for this subdomain in this proc
                    do i = 0, ngllx-1
                        do j = 0, nglly-1
                            do k = 0, ngllz-1
                                !write(*,*) "i = ", i, "j = ", j, "k = ", k
                                ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)
                                globCoordMask(:,ipoint) = .true.
                                propMask(ipoint,:)      = .true.
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID)
                                        prop(ipoint, 0) = Tdomain%sdom%Density_(i,j,k,lnum)
                                        prop(ipoint, 1) = Tdomain%sdom%Lambda_ (i,j,k,lnum)
                                        prop(ipoint, 2) = Tdomain%sdom%Mu_     (i,j,k,lnum)
                                    case (DM_FLUID)
                                        prop(ipoint, 0) = Tdomain%fdom%Density_(i,j,k,lnum)
                                        prop(ipoint, 1) = Tdomain%fdom%Lambda_ (i,j,k,lnum)
                                        prop(ipoint, 2) = 0
                                    case (DM_SOLID_PML)
                                        prop(ipoint, 0) = Tdomain%spmldom%Density_(i,j,k,lnum)
                                        prop(ipoint, 1) = Tdomain%spmldom%Lambda_ (i,j,k,lnum)
                                        prop(ipoint, 2) = Tdomain%spmldom%Mu_     (i,j,k,lnum)
                                    case (DM_FLUID_PML)
                                        prop(ipoint, 0) = Tdomain%fpmldom%Density_(i,j,k,lnum)
                                        prop(ipoint, 1) = Tdomain%fpmldom%Lambda_ (i,j,k,lnum)
                                        prop(ipoint, 2) = Tdomain%fpmldom%Mu_     (i,j,k,lnum)
                                end select
                            end do
                        end do
                    end do !END Loop over GLLs
                end if
            end do

            nSubDPoints(mat) = count(globCoordMask(0,:))

            call write_ResultHDF5Unstruct_MPI(                                      &
                reshape(pack(Tdomain%GlobCoord(:,:),       &
                mask = globCoordMask(:,:)),        &
                shape = [3, nSubDPoints(mat)]),    &
                reshape(pack(prop(:,:),                    &
                mask = propMask(:,:)),             &
                shape = [nSubDPoints(mat), 3]),    &
                trim(procFileName)//"_view",                &
                rg, trim(h5folder), Tdomain%communicateur, &
                ["_proc", "_subD"], [rg, mat], HDF5NameList(mat))
        end do

        !Writing XMF File
        if(rg == 0) write(*,*) "  Writing visualization XMF file"
        call writeXMF_RF_MPI(nProp, HDF5NameList, nSubDPoints, Tdomain%subD_exist, Tdomain%n_dime, &
            trim(string_join(procFileName,"-TO_VIEW")), rg, trim(XMFfolder),     &
            Tdomain%communicateur, trim(h5_to_xmf),                               &
            ["Density","Lambda ","Mu     "])

        deallocate(prop)
        deallocate(HDF5NameList)
        deallocate(globCoordMask)
        deallocate(propMask)
        deallocate(nSubDPoints)

    end subroutine create_prop_visu_files

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    function propOnFile(Tdomain, mat) result(authorization)

        !INPUTS
        type (domain), intent (in), target :: Tdomain
        integer, intent(in) :: mat

        !OUTPUT
        logical :: authorization

        !LOCAL
        integer :: assocMat

        assocMat = Tdomain%sSubdomain(mat)%assocMat
        authorization = .false.

        !if(Tdomain%subD_exist(mat)) then
        if(Tdomain%sSubDomain(assocMat)%material_type == "R") then
            authorization = .true.
        end if
        !end if

    end function

end module build_prop_files
