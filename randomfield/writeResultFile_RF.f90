!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module writeResultFile_RF

    use displayCarvalhol
    use statistics_RF
    use math_RF
    use hdf5
    use mpi

    interface write_ResultHDF5
        module procedure write_ResultHDF5Unstruct_MPI
    end interface write_ResultHDF5

contains

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine write_ResultHDF5Unstruct_MPI(xPoints, randField, fileName, rang, folderPath, &
        communicator, labels, indexes, HDF5Name)
        implicit none

        !INPUTS
        double precision,  dimension(1:,1:), intent(in) :: xPoints, randField;
        character (len=*)                , intent(in) :: filename;
        integer                          , intent(in) :: rang;
        character(len=*)                 , intent(in) :: folderPath
        integer                          , intent(in) :: communicator
        character(len=*) , dimension(1:), optional  , intent(in) :: labels
        integer          , dimension(1:), optional  , intent(in) :: indexes

        !OUTPUTS
        character(len=110) , optional  , intent(out) ::HDF5Name

        !HDF5 VARIABLES
        character(len=110)             :: fileHDF5Name, fullPath !File name
        character(len=30)              :: eventName, coordName;        !Dataset names
        integer(HID_T)                 :: file_id       !File identifier
        integer(HID_T)                 :: dset_id       !Dataset identifier
        integer(HID_T)                 :: dspace_id     !Dataspace identifier
        integer                        :: rank = 2      !Dataset rank (number of dimensions)
        integer(HSIZE_T), dimension(2) :: dims !Dataset dimensions
        integer                        :: error !Error flag

        !LOCAL VARIABLES
        integer :: nDim, Nmc, nPoints, i
        integer :: effectComm
        character (len=12) :: numberStr, rangStr;

        effectComm = communicator
        nDim       = size(xPoints , 1)
        nPoints    = size(randField, 1)
        Nmc        = size(randField, 2)

        if(.not. present(labels)) then
            write(rangStr,'(I8)'  ) rang
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

        if (nDim > 3) then
            write(*,*) "Dimension exceeds 3, HDF file won't be created"
        else

            call h5open_f(error) ! Initialize FORTRAN interface.
            call h5fcreate_f(fullPath, H5F_ACC_TRUNC_F, file_id, error) ! Create a new file using default properties.

            dims = shape(xPoints)
            write(coordName,'(A)') "XYZ"
            call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
            call h5dcreate_f(file_id, coordName, H5T_NATIVE_DOUBLE,  &
                dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, &
                xPoints, dims, error) ! Write the dataset.
            call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
            call h5sclose_f(dspace_id, error) ! Terminate access to the data space.

            dims(1) = size(randField,1)
            dims(2) = 1 !One random field in each dataset

            do i = 1, Nmc
                write(numberStr,'(I8)'  ) i
                numberStr = adjustL(numberStr)
                write(eventName,'(2A)') "RF_", trim(numberStr)

                call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
                call h5dcreate_f(file_id, eventName, H5T_NATIVE_DOUBLE,  &
                    dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
                call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, randField(:,i), dims, error) ! Write the dataset.
                call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
                call h5sclose_f(dspace_id, error) ! Terminate access to the data space.
            end do

            call h5fclose_f(file_id, error) ! Close the file.
            call h5close_f(error) ! Close FORTRAN interface.

        end if

        if(present(HDF5Name)) then
            HDF5Name = trim(adjustL(fileHDF5Name))
            HDF5Name = adjustL(HDF5Name)
        end if

    end subroutine write_ResultHDF5Unstruct_MPI


    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine writeXMF_RF_MPI(nSamples, HDF5nameList, nPointList, mask, nDim, fileName, rang, folderPath, &
        communicator, HDF5relativePath, attName, byProc)
        implicit none

        !INPUTS
        integer                           , intent(in) :: nSamples, nDim;
        logical         , dimension(:)    , intent(in) :: mask
        character(len=*), dimension(1:)   , intent(in) :: HDF5nameList;
        integer         , dimension(1:)   , intent(in) :: nPointList;
        character(len=*)                  , intent(in) :: filename;
        integer                           , intent(in) :: rang;
        character(len=*)                  , intent(in) :: folderPath
        integer                           , intent(in) :: communicator
        character(len=*), optional        , intent(in) :: HDF5relativePath
        character(len=*), dimension(1:)   , optional  , intent(in) :: attName;
        logical                           , optional  , intent(in) :: byProc !To divise the mesh by proc and not by subdomain

        !LOCAL VARIABLES
        integer             :: Nmc, i, j, file, nb_procs, code;
        integer             :: effectComm
        character (len=110) :: fileXMFName, fullPathXMF, effecHDF5path;
        character (len=35)  :: meshName;
        character (len=50) , dimension(:), allocatable :: effectAttName;
        integer            , dimension(:), allocatable :: all_nPointList
        character (len=110), dimension(:), allocatable :: all_HDF5nameList
        logical            , dimension(:), allocatable :: all_mask

        effectComm = communicator

        call MPI_COMM_SIZE(effectComm, nb_procs, code)

        !Common parameters
        Nmc         = nSamples

        if(rang == 0) then
            allocate(all_nPointList(nb_procs*size(HDF5nameList)))
            allocate(all_HDF5nameList(nb_procs*size(HDF5nameList)))
            allocate(all_mask(nb_procs*size(mask)))
            allocate(effectAttName(Nmc))
        end if

        call MPI_GATHER(nPointList    , size(nPointList), MPI_INTEGER,     &
            all_nPointList, size(nPointList), MPI_INTEGER,     &
            0         , effectComm    , code)

        call MPI_GATHER(HDF5nameList    , len(HDF5nameList)*size(HDF5nameList), MPI_CHARACTER,     &
            all_HDF5nameList, len(HDF5nameList)*size(HDF5nameList), MPI_CHARACTER,     &
            0              , effectComm       , code)

        call MPI_GATHER(mask    , size(mask), MPI_LOGICAL,     &
            all_mask, size(mask), MPI_LOGICAL,     &
            0      , effectComm       , code)

        if(rang == 0 .and. (nDim == 2 .or. nDim == 3)) then

            fileXMFName = string_join(fileName,".xmf")
            fullPathXMF = string_join(folderPath, "/"//fileXMFName)
            !if(rang == 0) write(*,*) "fileXMFName = ", fileXMFName
            !if(rang == 0) write(*,*) "fullPathXMF = ", fullPathXMF
            !Optional inputs
            if(present(attName)) then
                if (size(attName)==Nmc) then
                    effectAttName = attName
                end if
            else
                do i = 1, Nmc
                    effectAttName(i) = stringNumb_join("RF_", i)
                end do
            end if

            if(present(HDF5relativePath)) then
                effecHDF5path = string_join(HDF5relativePath, "/")
            else
                effecHDF5path = "./"
            end if

            !Building file
            file=21;
            open (unit = file , file = fullPathXMF, action = 'write')

            write (file,'(A)'      )'<?xml version="1.0" ?>'
            write (file,'(A)'      )'<Xdmf Version="2.0">'
            write (file,'(A)'      )' <Domain>'
            write (file,'(A)'      )'  <Grid CollectionType="Spatial" GridType="Collection">' !Opens the Collection

            do j = 1, size(all_HDF5nameList)
                if(all_mask(j)) then
                    !write(*,*) trim(all_HDF5nameList(j))
                    meshName = trim(stringNumb_join("SubD_", mod(j-1,size(HDF5nameList))))//&
                        trim(stringNumb_join("Proc_", int((j-1)/size(HDF5nameList))))
                    if(present(byProc)) then
                        if(byProc) meshName = stringNumb_join("Proc_", int((j-1)/size(HDF5nameList)))
                    end if
                    write (file,'(3A)'     )'   <Grid Name="',trim(meshName),'" GridType="Uniform">' !START Writing the data of one subdomain
                    write (file,'(3A)'     )'    <Topology Type="Polyvertex" NodesPerElements="1" NumberOfElements="',trim(numb2String(all_nPointList(j))),'">'
                    write (file,'(A)'      )'    </Topology>'
                    if(nDim == 2) write (file,'(A)'      )'     <Geometry GeometryType="XY">'
                    if(nDim == 3) write (file,'(A)'      )'     <Geometry GeometryType="XYZ">'
                    write (file,'(5A)'     )'      <DataItem Name="Coordinates" Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(numb2String(all_nPointList(j))), ' ',trim(numb2String(nDim))  ,'">'
                    write (file,'(4A)'     )'           ',trim(effecHDF5path),trim(all_HDF5nameList(j)),':/XYZ'
                    write (file,'(A)'      )'      </DataItem>'
                    write (file,'(A)'      )'    </Geometry>'

                    do i = 1, Nmc
                        write (file,'(3A)'     )'     <Attribute Name="',trim(effectAttName(i)),'" Center="Node" AttributeType="Scalar">'
                        write (file,'(3A)'     )'      <DataItem Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(numb2String(all_nPointList(j))),'">'
                        write (file,'(5A)'     )'          ',trim(effecHDF5path),trim(all_HDF5nameList(j)),":/", trim(stringNumb_join("RF_", i))
                        write (file,'(A)'      )'       </DataItem>'
                        write (file,'(A)'      )'     </Attribute>'
                    end do

                    write (file,'(A)'      )'   </Grid>' !END Writing the data of one subdomain
                end if
            end do !END Loop Over HDF5 names

            write (file,'(A)'      )'  </Grid>' !Closes the Collection
            write (file,'(A)'      )' </Domain>'
            write (file,'(A)'      )'</Xdmf>'

            close(file)
        else if(rang == 0) then
            write(*,*) "No XDMF Model prepared for this dimension"
            write(*,*) "The file won't be created"
        end if

        if(rang == 0) then
            deallocate(all_nPointList)
            deallocate(all_HDF5nameList)
            deallocate(effectAttName)
        end if

    end subroutine writeXMF_RF_MPI

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine write_MatlabTable(randField, filename)

        implicit none
        !INPUT
        double precision,  dimension(:,:), intent(in) :: randField;
        character (len=*)              ,   intent(in) :: filename;

        !OUTPUT - Result File

        !LOCAL VARIABLES
        integer            :: i, j, file, nColumns;
        character (len=40) :: doubleFmt;
        character (len=50) :: path

        !write(*,*) "";
        !write(*,*) "------------START Writing result file-----------------------";
        !write(*,*) "";

        path = trim(adjustL(filename));

        nColumns = size(randField,2);

        write(doubleFmt, fmt="(A,i10,A)") "(", nColumns, "F25.16)";

        !>>>>>>>>>>>>>>>>>>>>
        file=11;
        open (unit = file , file = trim(path)//"MATLAB", action = 'write')
        write (file,"("//doubleFmt//")") ((randField(i,j),j=1, size(randField,2)), i=1, size(randField,1));
        close(file)

    end subroutine write_MatlabTable

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function string_join(string1, string2) result(stringTot)

        implicit none

        !INPUT
        character (len=*)  , intent(in) :: string1, string2;

        !OUTPUT
        character (len=100) :: stringTot;

        stringTot = trim(adjustL(string1))//trim(adjustL(string2))

    end function string_join

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine jump_To_Line(fileUnit, line)

        implicit none

        !INPUT
        integer, intent(in) :: fileUnit, line;
        !LOCAL
        integer :: i

        do i = 1, line-1
            read(fileUnit, *)
        end do

    end subroutine jump_To_Line

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function stringNumb_join(string, number) result(stringTot)

        implicit none

        !INPUT
        character (len=*)  , intent(in) :: string
        integer            , intent(in) :: number

        !OUTPUT
        character (len=100) :: stringTot;

        !LOCAL
        character (len=30)  :: nmbString

        write(nmbString, fmt='(I8)') number
        stringTot = string_join(string, nmbString)

    end function

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    function numb2String(number) result(stringTot)

        implicit none

        !INPUT
        integer, intent(in) :: number

        !OUTPUT
        character (len=30) :: stringTot;

        !LOCAL
        character (len=30)  :: nmbString

        write(nmbString, fmt='(I8)') number
        stringTot = trim(adjustL(nmbString))

    end function

end module writeResultFile_RF

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
