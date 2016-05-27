module interpolation

    use mpi
    use hdf5
    use sdomain
    use mpi
    use constants
    use sem_hdf5

   implicit none

contains

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine interpolateToMesh(BBoxFileName, coordList, UNV_randField, rang, &
                                     xMinLoc_In, xMaxLoc_In)
            implicit none

            !INPUT
            character(len=*), intent(in)      :: BBoxFileName
            double precision, dimension(:,:), intent(in) :: coordList
            integer, intent(in) :: rang
            double precision, dimension(:), intent(in), optional :: xMinLoc_In, xMaxLoc_In
            !OUTPUT
            double precision, dimension(:,:), intent(out) :: UNV_randField
            !LOCAL
            integer :: nDim, Nmc
            !logical :: independent
            character(len=50) :: attr_Name, dset="samples"
            integer :: hdferr
            integer(HID_T) :: file_id, space_id, dset_id, mem_id
            integer(HSIZE_T), dimension(size(coordList,1)) :: offset, locDims
            integer, dimension(size(coordList,1)) :: xNStep, coordPosInt
            integer, dimension(size(coordList,1), 2**size(coordList,1)) :: neighCoord
            double precision, dimension(size(coordList,1)) :: coordPos
            double precision, dimension(size(coordList,1)) :: distance
            double precision, dimension(size(coordList,1)) :: xMinGlob, xMaxGlob, xStep
            double precision, dimension(:,:), allocatable, target :: BB_randField
            double precision, dimension(size(coordList,1)) :: xMin_Loc_UNV, xMax_Loc_UNV
            integer, dimension(size(coordList,1)) :: minPos, maxPos, extent
            integer(HSIZE_T), dimension(2) :: locShape, zero2D
            integer :: i, j
            double precision, dimension(:,:)    , pointer :: BB_2D
            double precision, dimension(:,:,:)  , pointer :: BB_3D
            double precision :: weight


            call h5open_f(hdferr) ! Initialize FORTRAN interface.
            call h5fopen_f(trim(BBoxFileName), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
            if(hdferr /= 0) stop("ERROR OPENING FILE")
            !write(*,*) "hdferr = ", hdferr

            !READING SCALARS----------------------------
            !BOOL
            !attr_name = "independent"
            !call read_h5attr_bool(file_id, trim(adjustL(attr_name)), independent)

            !INTEGERS
            attr_name = "nDim"
            call read_attr_int(file_id, trim(adjustL(attr_name)), nDim)
            attr_name = "Nmc"
            call read_attr_int(file_id, trim(adjustL(attr_name)), Nmc)

            !DOUBLE VEC
            attr_name = "xMinGlob"
            call read_attr_real_vec(file_id, attr_name, xMinGlob)
            attr_name = "xMaxGlob"
            call read_attr_real_vec(file_id, attr_name, xMaxGlob)
            attr_name = "xStep"
            call read_attr_real_vec(file_id, attr_name, xStep)

            xNStep = nint((xMaxGlob-xMinGlob)/xStep) +1

            !DEFINE LOCAL BOUNDING BOX
            if(present(xMinLoc_In) .and. present(xMaxLoc_In)) then
                xMin_Loc_UNV = xMinLoc_In
                xMax_Loc_UNV = xMaxLoc_In
            else
                do i = 1, nDim
                    xMin_Loc_UNV(i) = minval(coordList(i,:))
                    xMax_Loc_UNV(i) = maxval(coordList(i,:))
                end do
            end if

            minPos = floor((xMin_Loc_UNV-xMinGlob)/xStep) + 1
            maxPos = ceiling((xMax_Loc_UNV-xMinGlob)/xStep) + 1
            where(minPos < 1) minPos = 1
            where(maxPos > xNStep) maxPos = xNStep

            extent = maxPos - minPos + 1

            allocate(BB_randField(product(extent),1))

            !READING MATRIX BLOCK----------------------------------------------------------------
            locDims = extent
            call h5dopen_f(file_id, trim(dset), dset_id, hdferr)! Open Dataset
            call h5dget_space_f(dset_id, space_id, hdferr) !Open Dataset Space
            !call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr) !Measure Dataset Space

            !allocate(STA%randField(product(dims),1))
            offset = minPos-1
            locShape = shape(BB_randField)
            zero2D = 0
           !For hyperslab lecture

            !IN
            call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, offset, locDims, hdferr) !Select Hyperslab IN

            !OUT
            call h5screate_simple_f(2, locShape, mem_id, hdferr) !Create memory dataspace
            call h5sselect_hyperslab_f(mem_id, H5S_SELECT_SET_F, zero2D, locShape, hdferr) !Select Hyperslab OUT
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, BB_randField, locShape, hdferr, mem_id, space_id) !Read Dataset Hyperslab


            call h5dclose_f(dset_id, hdferr) !Close Dataset
            call h5sclose_f(space_id, hdferr) !Close Dataset Space

            call h5fclose_f(file_id, hdferr) ! Close the file.
            call h5close_f(hdferr) ! Close FORTRAN interface.


            !INTERPOLATION
            if(nDim == 2) then
                BB_2D(minPos(1):maxPos(1), minPos(2):maxPos(2)) => BB_randField
                neighCoord(:,1) = [0, 0]
                neighCoord(:,2) = [1, 0]
                neighCoord(:,3) = [0, 1]
                neighCoord(:,4) = [1, 1]

            else if(nDim == 3) then
                BB_3D(minPos(1):maxPos(1), minPos(2):maxPos(2), minPos(3):maxPos(3)) => BB_randField
                neighCoord(:,1) = [0, 0, 0]
                neighCoord(:,2) = [1, 0, 0]
                neighCoord(:,3) = [0, 1, 0]
                neighCoord(:,4) = [0, 0, 1]
                neighCoord(:,5) = [1, 1, 0]
                neighCoord(:,6) = [0, 1, 1]
                neighCoord(:,7) = [1, 0, 1]
                neighCoord(:,8) = [1, 1, 1]

            end if

            UNV_randField(:,:) = 0

            do i = 1, size(coordList, 2)

                !MAPPING COORD In BB
                coordPos = ((coordList(:,i)-xMinGlob)/xStep) + 1.0D0
                coordPosInt = floor(coordPos)
                where(coordPosInt == maxPos) coordPosInt = coordPosInt - 1 !Dealing with points on the positive border
                if(any(coordPosInt<0)) stop("coordPosInt smaller than 1")

                !Applying Values
                do j = 1, size(neighCoord, 2)

                    distance(:) = abs(coordPos - dble(coordPosInt+neighCoord(:,j)))
                    weight      = product(1.0D0 - distance)
                    !weight      = 1.0D0 - sqrt(sum(distance**2))/sqrt(dble(nDim))

                    !write(*,*) "weight = ", weight

                    if(any(coordPosInt(:)+neighCoord(:,j) > maxPos)) then
                        stop(" ERROR! UNV TRIED POSITION OUT OF RANGE")

                    end if

                    if(any(coordPosInt(:)+neighCoord(:,j) < minPos)) then
                        stop(" ERROR! UNV TRIED POSITION OUT OF RANGE")
                    end if

                    if(nDim == 2) then
                        UNV_randField(i,1) = UNV_randField(i,1) +                  &
                            (                                     &
                            BB_2D(coordPosInt(1)+neighCoord(1,j), &
                            coordPosInt(2)+neighCoord(2,j)) &
                            * weight                              &
                            )
                    else if (nDim == 3) then
                        UNV_randField(i,1) = UNV_randField(i,1) +                  &
                            (                                     &
                            BB_3D(coordPosInt(1)+neighCoord(1,j), &
                            coordPosInt(2)+neighCoord(2,j), &
                            coordPosInt(3)+neighCoord(3,j)) &
                            * weight                              &
                            )
                    end if

                end do


!                if(nDim == 2) then
!                    UNV_randField(i,1) = BB_2D(coordPosInt(1), coordPosInt(2))
!                else if (nDim == 3) then
!                    UNV_randField(i,1) = BB_3D(coordPosInt(1), coordPosInt(2), coordPosInt(3))
!                end if

!                coordPosInt = nint(coordPos)
!
!                if(nDim == 2) then
!                    UNV_randField(i,1) = BB_2D(coordPosInt(1), coordPosInt(2))
!                else if (nDim == 3) then
!                    UNV_randField(i,1) = BB_3D(coordPosInt(1), coordPosInt(2), coordPosInt(3))
!                end if

            end do




            if (allocated(BB_randField))   deallocate(BB_randField)

        end subroutine interpolateToMesh


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine interpolateToMesh_V2(BBoxFileName, coordList, UNV_randField, rang, &
                                 xMinLoc_In, xMaxLoc_In, mat, Tdomain)
        implicit none

        !INPUT
        character(len=*), intent(in)      :: BBoxFileName
        double precision, dimension(1:,1:), intent(in) :: coordList
        integer, intent(in) :: rang
        double precision, dimension(1:), intent(in), optional :: xMinLoc_In, xMaxLoc_In
        type (domain), intent (INOUT), target :: Tdomain
        integer, intent(in) :: mat
        !OUTPUT
        double precision, dimension(1:,1:), intent(out) :: UNV_randField
        !LOCAL
        integer :: nDim, Nmc
        !logical :: independent
        character(len=50) :: attr_Name, dset="samples"
        integer :: hdferr
        integer(HID_T) :: file_id, space_id, dset_id, mem_id
        integer(HSIZE_T), dimension(size(coordList,1)) :: offset, locDims
        integer, dimension(size(coordList,1)) :: xNStep, coordPosInt
        integer, dimension(size(coordList,1), 2**size(coordList,1)) :: neighCoord
        integer, dimension(size(coordList,2)) :: nContrib
        double precision, dimension(size(coordList,1)) :: coordPos
        double precision, dimension(size(coordList,1)) :: distance
        double precision, dimension(size(coordList,1)) :: xMinGlob, xMaxGlob, xStep
        double precision, dimension(:,:), allocatable, target :: BB_randField
        double precision, dimension(1:size(coordList,1)) :: xMin_Loc_UNV, xMax_Loc_UNV
        integer, dimension(1:size(coordList,1)) :: minPos, maxPos, extent
        integer(HSIZE_T), dimension(2) :: locShape, zero2D
        integer(kind=8) :: extTotal
        integer :: ii, jj, kk
        integer :: i, j
        double precision, dimension(:,:)    , pointer :: BB_2D
        double precision, dimension(:,:,:)  , pointer :: BB_3D
        double precision :: weight
        logical :: findPropertyonTable
        real :: t_0, t_f
        integer :: n, lnum, elem_mat, ngll, ipoint


        call h5open_f(hdferr) ! Initialize FORTRAN interface.
        call h5fopen_f(trim(BBoxFileName), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
        if(hdferr /= 0) then
            write(*,*) " "
            write(*,*) "-----------------------------------------------------------------------"
            write(*,*) "Your simulation requires random properties. Sample them and re-run SEM"
            write(*,*) "   Tried to open ", trim(BBoxFileName)
            write(*,*) "-----------------------------------------------------------------------"
            write(*,*) " "
            stop(" ")
        end if

        !INTEGERS
        attr_name = "nDim"
        call read_attr_int(file_id, trim(adjustL(attr_name)), nDim)
        attr_name = "Nmc"
        call read_attr_int(file_id, trim(adjustL(attr_name)), Nmc)

        !DOUBLE VEC
        attr_name = "xMinGlob"
        call read_attr_real_vec(file_id, attr_name, xMinGlob)
        attr_name = "xMaxGlob"
        call read_attr_real_vec(file_id, attr_name, xMaxGlob)
        attr_name = "xStep"
        call read_attr_real_vec(file_id, attr_name, xStep)

        xNStep = nint((xMaxGlob-xMinGlob)/xStep) +1

        !DEFINE LOCAL BOUNDING BOX
        if(present(xMinLoc_In) .and. present(xMaxLoc_In)) then
            xMin_Loc_UNV = xMinLoc_In
            xMax_Loc_UNV = xMaxLoc_In
        else
            do i = 1, nDim
                xMin_Loc_UNV(i) = minval(coordList(i,:))
                xMax_Loc_UNV(i) = maxval(coordList(i,:))
            end do
        end if

        !write(*,*) "xMinGlob = ", xMinGlob
        !write(*,*) "xMaxGlob = ", xMaxGlob
        !write(*,*) "xMin_Loc_UNV = ", xMin_Loc_UNV
        !write(*,*) "xMax_Loc_UNV = ", xMax_Loc_UNV
        !write(*,*) "shape(coordList) = ", shape(coordList)

        minPos = floor((xMin_Loc_UNV-xMinGlob)/xStep) + 1
        maxPos = ceiling((xMax_Loc_UNV-xMinGlob)/xStep) + 1
        where(minPos < 1) minPos = 1
        where(maxPos > xNStep) maxPos = xNStep

        extent = maxPos - minPos + 1
        extTotal = product(int(extent,8))
        write(*,*) "minPos = ", minPos
        write(*,*) "minPos Coord = ", xMin_Loc_UNV
        write(*,*) "maxPos = ", maxPos
        write(*,*) "maxPos Coord = ", xMax_Loc_UNV
        write(*,*) "extent = ", extent
        write(*,*) "extTotal = ", extTotal

        allocate(BB_randField(extTotal,1))

        !READING MATRIX BLOCK----------------------------------------------------------------
        locDims = extent
        !write(*,*) "Open hdf5 file"
        !call cpu_time(t_0)
        call h5dopen_f(file_id, trim(dset), dset_id, hdferr)! Open Dataset
        call h5dget_space_f(dset_id, space_id, hdferr) !Open Dataset Space
        !call h5sget_simple_extent_dims_f(space_id, dims, maxdims, hdferr) !Measure Dataset Space

        offset = minPos-1
        locShape = shape(BB_randField)
        zero2D = 0

        !For hyperslab lecture

        !IN
        call h5sselect_hyperslab_f(space_id, H5S_SELECT_SET_F, offset, locDims, hdferr) !Select Hyperslab IN

        !OUT
        call h5screate_simple_f(2, locShape, mem_id, hdferr) !Create memory dataspace
        call h5sselect_hyperslab_f(mem_id, H5S_SELECT_SET_F, zero2D, locShape, hdferr) !Select Hyperslab OUT
        !call cpu_time(t_f)
        !write(*,*) "time = ", t_f - t_0
        !call cpu_time(t_0)

        !READING
        !write(*,*) "Reading hdf5 file"
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, BB_randField, locShape, hdferr, mem_id, space_id) !Read Dataset Hyperslab

        !call cpu_time(t_f)
        !write(*,*) "time = ", t_f - t_0
        !call cpu_time(t_0)

        call h5dclose_f(dset_id, hdferr) !Close Dataset
        call h5sclose_f(space_id, hdferr) !Close Dataset Space

        call h5fclose_f(file_id, hdferr) ! Close the file.
        call h5close_f(hdferr) ! Close FORTRAN interface.


        !INTERPOLATION
        if(nDim == 2) then
            BB_2D(minPos(1):maxPos(1), minPos(2):maxPos(2)) => BB_randField
            neighCoord(:,1) = [0, 0]
            neighCoord(:,2) = [1, 0]
            neighCoord(:,3) = [0, 1]
            neighCoord(:,4) = [1, 1]

        else if(nDim == 3) then
            BB_3D(minPos(1):maxPos(1), minPos(2):maxPos(2), minPos(3):maxPos(3)) => BB_randField
            neighCoord(:,1) = [0, 0, 0]
            neighCoord(:,2) = [1, 0, 0]
            neighCoord(:,3) = [0, 1, 0]
            neighCoord(:,4) = [0, 0, 1]
            neighCoord(:,5) = [1, 1, 0]
            neighCoord(:,6) = [0, 1, 1]
            neighCoord(:,7) = [1, 0, 1]
            neighCoord(:,8) = [1, 1, 1]

        end if

        UNV_randField(:,:) = 0
        nContrib(:) = 0

        !write(*,*) "Interpolating"

        !Applying properties to elements
        do n = 0, Tdomain%n_elem-1
            elem_mat = Tdomain%specel(n)%mat_index

            if(elem_mat /= mat) cycle

            lnum = Tdomain%specel(n)%lnum
            ngll = Tdomain%sSubDomain(mat)%NGLL

            !Properties by GLL
            do ii = 0, ngll-1
                do jj = 0, ngll-1
                    do kk = 0, ngll-1
                        ipoint  = Tdomain%specel(n)%Iglobnum(ii,jj,kk) + 1

                        !MAPPING COORD In BB
                        coordPos = ((coordList(:,ipoint)-xMinGlob)/xStep) + 1.0D0
                        coordPosInt = floor(coordPos)
                        where(coordPosInt == maxPos) coordPosInt = coordPosInt - 1 !Dealing with points on the positive border
                        where(coordPosInt < minPos) coordPosInt = coordPosInt + 1
                        if(any(coordPosInt > maxPos)) stop("coordPosInt bigger than maxPos")
                        if(any(coordPosInt < minPos)) stop("coordPosInt bigger than minPos")
                        if(any(coordPosInt<0)) stop("coordPosInt smaller than 1")

                        nContrib(ipoint) = nContrib(ipoint) + 1

                        !Applying Values
                        do j = 1, size(neighCoord, 2)
                            findPropertyOnTable = .true.
                            distance(:) = abs(coordPos - dble(coordPosInt+neighCoord(:,j)))
                            weight      = product(1.0D0 - distance)

                            !write(*,*) "weight = ", weight

                            if(any(coordPosInt(:)+neighCoord(:,j) > maxPos)) then
                                findPropertyOnTable = .false.
                                write(*,*) "Error in rang ", rang
                                write(*,*) "  coordList(:,ipoint) = ", coordList(:,ipoint)
                                write(*,*) "   coordPos = ", coordPos
                                write(*,*) "          j = ", j
                                write(*,*) "coordPosInt(:)+neighCoord(:,j) = ", coordPosInt(:)+neighCoord(:,j)
                                write(*,*) "maxPos = ", maxPos
                                stop(" ERROR! UNV TRIED POSITION OUT OF RANGE (>Max)")

                            end if

                            if(any(coordPosInt(:)+neighCoord(:,j) < minPos)) then
                                findPropertyOnTable = .false.
                                write(*,*) "Error in rang ", rang
                                write(*,*) "  coordList(:,i) = ", coordList(:,ipoint)
                                write(*,*) "   coordPos = ", coordPos
                                write(*,*) "          j = ", j
                                write(*,*) "coordPosInt(:)+neighCoord(:,j) = ", coordPosInt(:)+neighCoord(:,j)
                                write(*,*) "minPos = ", minPos
                                stop(" ERROR! UNV TRIED POSITION OUT OF RANGE (<Min)")
                            end if

                            if(findPropertyOnTable) then
                                if(nDim == 2) then
                                    UNV_randField(ipoint,1) = UNV_randField(ipoint,1) +      &
                                        (                                     &
                                        BB_2D(coordPosInt(1)+neighCoord(1,j), &
                                        coordPosInt(2)+neighCoord(2,j)) &
                                        * weight                              &
                                        )

                                else if (nDim == 3) then
                                    UNV_randField(ipoint,1) = UNV_randField(ipoint,1) +     &
                                        (                                     &
                                        BB_3D(coordPosInt(1)+neighCoord(1,j), &
                                        coordPosInt(2)+neighCoord(2,j), &
                                        coordPosInt(3)+neighCoord(3,j)) &
                                        * weight                              &
                                        )
                                end if
                            end if
                        end do
                    end do
                end do
            end do

        end do

        where(nContrib > 0) UNV_randField(:,1) = UNV_randField(:,1)/dble(nContrib)

        !call cpu_time(t_f)
        !write(*,*) "time = ", t_f - t_0
        !call cpu_time(t_0)

        if (allocated(BB_randField))   deallocate(BB_randField)

    end subroutine interpolateToMesh_V2

end module interpolation
