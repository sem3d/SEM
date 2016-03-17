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
        integer :: mat, assocMat, code

        !Writing hdf5 files
        if(rg == 0) write(*,*)
        do mat = 0, Tdomain%n_mat - 1
            assocMat = Tdomain%sSubdomain(mat)%assocMat
            if(propOnFile(Tdomain, assocMat)) then
                !if(rg == 0) write(*,*) "  Material ", mat, " will have properties on file"
                if(is_rand(Tdomain%sSubdomain(assocMat))) then
                !write(*,*) " rang ", rg," Flag 1---------------------", "mat = ", mat
                call MPI_BARRIER(Tdomain%communicateur, code)
                    !write(*,*) "FOR SURE"
                    call build_random_properties(Tdomain, rg, mat)
                end if

            end if
        end do

    end subroutine create_prop_files

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine apply_prop_files(Tdomain, rg)
        use sample_RF

        implicit none
        !INPUTS
        type (domain), intent (INOUT), target :: Tdomain
        integer      , intent(IN) :: rg

        !LOCAL
        integer :: mat, n, ipoint, i, j, k, lnum
        integer :: ngll, elem_mat
        double precision, dimension(:, :), allocatable :: interpolatedRF !Properties
        double precision, dimension(:), allocatable :: kappa
        integer :: assocMat, propId

        !Putting properties on elements
        if(Tdomain%any_Random) then
            allocate(interpolatedRF(0:size(Tdomain%GlobCoord, 2) - 1, 0:nProp - 1))
            allocate(kappa(0:size(Tdomain%GlobCoord, 2) - 1))
        end if

        do mat = 0, Tdomain%n_mat - 1

            if(.not. Tdomain%subD_exist(mat)) cycle

            if(propOnFile(Tdomain, mat)) then
                if(rg == 0) write(*,*) "  Material ", mat, " have properties on file"
                assocMat = Tdomain%sSubdomain(mat)%assocMat
                if(is_rand(Tdomain%sSubdomain(assocMat))) then

                    if(rg == 0) write(*,*) "  Reading and Interpolating Random Properties"

                    do propId = 0, nProp - 1
                        call interpolateToMesh(BBoxFileName=Tdomain%sSubDomain(mat)%propFilePath(propId),  &
                                               coordList=Tdomain%GlobCoord, &
                                               UNV_randField=interpolatedRF(:,propId:propId), &
                                               rang=Tdomain%rank, &
                                               xMinLoc_In = Tdomain%sSubDomain(mat)%MinBound_Loc,  &
                                               xMaxLoc_In = Tdomain%sSubDomain(mat)%MaxBound_Loc)
                    end do
                end if

                kappa(:) = interpolatedRF(:, 1) + 2.*interpolatedRF(:, 2) /3.
                !S%DMu = S%Sspeed**2 * S%Ddensity
                !S%DLambda = (S%Pspeed**2 - 2 * S%Sspeed **2 ) * S%Ddensity
                !S%DKappa = S%DLambda + 2.*S%DMu /3.



                !Applying properties to elements
                do n = 0, Tdomain%n_elem-1
                    elem_mat = Tdomain%specel(n)%mat_index

                    if(elem_mat /= mat) cycle

                    lnum = Tdomain%specel(n)%lnum
                    ngll = Tdomain%sSubDomain(mat)%NGLL

                    !Properties by element
                    select case (Tdomain%specel(n)%domain)
                        case (DM_SOLID)
                            if (Tdomain%sdom%n_sls>0)  then
                                if (Tdomain%sdom%aniso) then
                                    Tdomain%sdom%Q_(:,:,:,lnum) = Tdomain%sSubDomain(mat)%Qmu
                                else
                                    Tdomain%sdom%Qs_(:,:,:,lnum) = Tdomain%sSubDomain(mat)%Qmu
                                    Tdomain%sdom%Qp_(:,:,:,lnum) = Tdomain%sSubDomain(mat)%Qpression
                                endif
                            endif
                        case (DM_FLUID)
                            !Nothing to do, all definitions are by point
                        case (DM_SOLID_PML)
                            !Nothing to do, all definitions are by point
                        case (DM_FLUID_PML)
                            !Nothing to do, all definitions are by point
                    end select

                    !Properties by GLL
                    do i = 0, ngll-1
                        do j = 0, ngll-1
                            do k = 0, ngll-1
                                ipoint  = Tdomain%specel(n)%Iglobnum(i,j,k)
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID)
                                        Tdomain%sdom%Density_(i,j,k,lnum) = interpolatedRF(ipoint, 0)
                                        Tdomain%sdom%Lambda_ (i,j,k,lnum) = interpolatedRF(ipoint, 1)
                                        Tdomain%sdom%Mu_     (i,j,k,lnum) = interpolatedRF(ipoint, 2)
                                        Tdomain%sdom%Kappa_  (i,j,k,lnum) = kappa(ipoint)
                                    case (DM_FLUID)
                                        Tdomain%fdom%IDensity_(i,j,k,lnum) = 1d0/interpolatedRF(ipoint, 0)
                                        Tdomain%fdom%Lambda_  (i,j,k,lnum) = interpolatedRF(ipoint, 1)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,lnum) = interpolatedRF(ipoint, 0)
                                        Tdomain%spmldom%Lambda_ (i,j,k,lnum) = interpolatedRF(ipoint, 1)
                                        Tdomain%spmldom%Mu_     (i,j,k,lnum) = interpolatedRF(ipoint, 2)
                                    case (DM_FLUID_PML)
                                        Tdomain%fpmldom%Density_(i,j,k,lnum) = interpolatedRF(ipoint, 0)
                                        Tdomain%fpmldom%Lambda_ (i,j,k,lnum) = interpolatedRF(ipoint, 1)
                                end select
                            end do
                        end do
                    end do
                end do

            end if
        end do


        if(allocated(interpolatedRF)) deallocate(interpolatedRF)
        if(allocated(kappa))          deallocate(kappa)

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
        integer :: ngll
        logical         , dimension(:,:) , allocatable :: globCoordMask, propMask
        double precision, dimension(:, :), allocatable :: prop !Properties
        integer         , dimension(:)   , allocatable :: nSubDPoints;
        character(len=110) , dimension(:), allocatable :: HDF5NameList

!        allocate(prop(0:size(Tdomain%GlobCoord,2)-1, 0:nProp-1)) !Subdomain properties Matrix ((:,0) = Dens, (:,1) = Lambda, (:,2) = Mu) per proc
!
!        !Creating  visualization files (by subdomain)
!        allocate(HDF5NameList(0:Tdomain%n_mat - 1))
!        allocate(nSubDPoints(0:Tdomain%n_mat - 1))
!        allocate(globCoordMask(0:size(Tdomain%GlobCoord,1)-1, 0:size(Tdomain%GlobCoord,2)-1))
!        allocate(propMask(0:size(Tdomain%GlobCoord,2)-1, 0:nProp-1))
!
!        HDF5NameList(:) = "not_Used"
!        nSubDPoints(:)  = 0
!
!        if(rg == 0) write(*,*) ">>>>Writing visualization hdf5 files"
!        do mat = 0, Tdomain%n_mat - 1
!            globCoordMask(:,:) = .false.
!            propMask(:,:)      = .false.
!            ngll = Tdomain%sSubDomain(mat)%NGLL
!            do n = 0, Tdomain%n_elem-1
!                lnum = Tdomain%specel(n)%lnum
!                if(Tdomain%specel(n)%mat_index == mat) then
!                    !Building masks for this subdomain in this proc
!                    do i = 0, ngll-1
!                        do j = 0, ngll-1
!                            do k = 0, ngll-1
!                                !write(*,*) "i = ", i, "j = ", j, "k = ", k
!                                ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)
!                                globCoordMask(:,ipoint) = .true.
!                                propMask(ipoint,:)      = .true.
!                                select case (Tdomain%specel(n)%domain)
!                                    case (DM_SOLID)
!                                        prop(ipoint, 0) = Tdomain%sdom%Density_(i,j,k,lnum)
!                                        prop(ipoint, 1) = Tdomain%sdom%Lambda_ (i,j,k,lnum)
!                                        prop(ipoint, 2) = Tdomain%sdom%Mu_     (i,j,k,lnum)
!                                    case (DM_FLUID)
!                                        prop(ipoint, 0) = Tdomain%fdom%IDensity_(i,j,k,lnum)
!                                        prop(ipoint, 1) = Tdomain%fdom%Lambda_ (i,j,k,lnum)
!                                        prop(ipoint, 2) = 0
!                                    case (DM_SOLID_PML)
!                                        prop(ipoint, 0) = Tdomain%spmldom%Density_(i,j,k,lnum)
!                                        prop(ipoint, 1) = Tdomain%spmldom%Lambda_ (i,j,k,lnum)
!                                        prop(ipoint, 2) = Tdomain%spmldom%Mu_     (i,j,k,lnum)
!                                    case (DM_FLUID_PML)
!                                        prop(ipoint, 0) = Tdomain%fpmldom%Density_(i,j,k,lnum)
!                                        prop(ipoint, 1) = Tdomain%fpmldom%Lambda_ (i,j,k,lnum)
!                                        prop(ipoint, 2) = 0
!                                end select
!                            end do
!                        end do
!                    end do !END Loop over GLLs
!                end if
!            end do
!
!            nSubDPoints(mat) = count(globCoordMask(0,:))
!
!            call write_ResultHDF5Unstruct_MPI(                                      &
!                reshape(pack(Tdomain%GlobCoord(:,:),       &
!                mask = globCoordMask(:,:)),        &
!                shape = [3, nSubDPoints(mat)]),    &
!                reshape(pack(prop(:,:),                    &
!                mask = propMask(:,:)),             &
!                shape = [nSubDPoints(mat), 3]),    &
!                trim(procFileName)//"_view",                &
!                rg, trim(h5folder), Tdomain%communicateur, &
!                ["_proc", "_subD"], [rg, mat], HDF5NameList(mat))
!        end do
!
!        !Writing XMF File
!        if(rg == 0) write(*,*) "  Writing visualization XMF file"
!        call writeXMF_RF_MPI(nProp, HDF5NameList, nSubDPoints, Tdomain%subD_exist, Tdomain%n_dime, &
!            trim(string_join(procFileName,"-TO_VIEW")), rg, trim(XMFfolder),     &
!            Tdomain%communicateur, trim(h5_to_xmf),                               &
!            ["Density","Lambda ","Mu     "])
!
!        deallocate(prop)
!        deallocate(HDF5NameList)
!        deallocate(globCoordMask)
!        deallocate(propMask)
!        deallocate(nSubDPoints)

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

        if(Tdomain%sSubdomain(assocMat)%initial_material_type == 'R') then
            authorization = .true.
        end if
        !end if

    end function

end module build_prop_files
