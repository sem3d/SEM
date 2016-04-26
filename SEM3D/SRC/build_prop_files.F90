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


        call build_random_properties(Tdomain, rg)


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
        double precision, dimension(:), allocatable :: lambda
        integer :: assocMat, propId

        !Putting properties on elements
        if(Tdomain%any_Random) then
            allocate(interpolatedRF(0:size(Tdomain%GlobCoord, 2) - 1, 0:nProp - 1))
            allocate(lambda(0:size(Tdomain%GlobCoord, 2) - 1))
        end if

        do mat = 0, Tdomain%n_mat - 1
        !do mat = 0, 0 !For Tests

            write(*,*) "  "
            write(*,*) "  Analyzing Material ", mat, "------------------- in rg ", rg
            !write(*,*) "  Tdomain%subD_exist ", Tdomain%subD_exist, "------------------- in rg ", rg

            if(.not. Tdomain%subD_exist(mat)) cycle

            !write(*,*) "  Exists in proc ", rg

            assocMat = Tdomain%sSubdomain(mat)%assocMat

            if(propOnFile(Tdomain, assocMat)) then
                !write(*,*) "  Material ", mat, " have properties on file"
                if(is_rand(Tdomain%sSubdomain(assocMat))) then

                    do propId = 0, nProp - 1
                    !write(*,*) "  Reading and Interpolating Random Properties on:"
                    !do propId = 0, 0 !For Tests
                        !write(*,*) " prop = ", propId, "rg = ", rg
                        write(*,*) trim(Tdomain%sSubDomain(mat)%propFilePath(propId))
                        call interpolateToMesh_V2(BBoxFileName=Tdomain%sSubDomain(mat)%propFilePath(propId),  &
                                               coordList=Tdomain%GlobCoord, &
                                               UNV_randField=interpolatedRF(:,propId:propId), &
                                               rang=rg, &
                                               xMinLoc_In = Tdomain%sSubDomain(mat)%MinBound_Loc,  &
                                               xMaxLoc_In = Tdomain%sSubDomain(mat)%MaxBound_Loc,  &
                                               mat = mat, &
                                               Tdomain  = Tdomain)
                        !write(*,*) " AFTER prop = ", propId, "rg = ", rg
                    end do

                    !write(*,*) " AFTER 2 rg = ", rg

                    !INFO
                    !interpolatedRF(:, 0) --Density
                    !interpolatedRF(:, 1) --Kappa
                    !interpolatedRF(:, 2) --Mu

                    lambda(:) = interpolatedRF(:, 1) - 2.0D0*interpolatedRF(:, 2) /3.0D0
                    !S%DMu = S%Sspeed**2 * S%Ddensity
                    !S%DLambda = (S%Pspeed**2 - 2 * S%Sspeed **2 ) * S%Ddensity
                    !S%DKappa = S%DLambda + 2.*S%DMu /3.
                    !S%DLambda = S%DKappa - 2.*S%DMu /3.

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
                                !Nothing to do, all definitions are on GLLs
                            case (DM_SOLID_PML)
                                !Nothing to do, all definitions are on GLLs
                            case (DM_FLUID_PML)
                                !Nothing to do, all definitions are on GLLs
                        end select

                        !Properties by GLL
                        do i = 0, ngll-1
                            do j = 0, ngll-1
                                do k = 0, ngll-1
                                    ipoint  = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    select case (Tdomain%specel(n)%domain)
                                        case (DM_SOLID)
                                            Tdomain%sdom%Density_(i,j,k,lnum) = interpolatedRF(ipoint, 0)
                                            Tdomain%sdom%Lambda_ (i,j,k,lnum) = lambda(ipoint)
                                            Tdomain%sdom%Mu_     (i,j,k,lnum) = interpolatedRF(ipoint, 2)
                                            Tdomain%sdom%Kappa_  (i,j,k,lnum) = interpolatedRF(ipoint, 1)
                                        case (DM_SOLID_PML)
                                            Tdomain%spmldom%Density_(i,j,k,lnum) = interpolatedRF(ipoint, 0)
                                            Tdomain%spmldom%Lambda_ (i,j,k,lnum) = lambda(ipoint)
                                            Tdomain%spmldom%Mu_     (i,j,k,lnum) = interpolatedRF(ipoint, 2)
                                        case (DM_FLUID)
                                            Tdomain%fdom%IDensity_(i,j,k,lnum) = 1d0/interpolatedRF(ipoint, 0)
                                            Tdomain%fdom%Lambda_  (i,j,k,lnum) = lambda(ipoint)
                                        case (DM_FLUID_PML)
                                            Tdomain%fpmldom%Density_(i,j,k,lnum) = interpolatedRF(ipoint, 0)
                                            Tdomain%fpmldom%Lambda_ (i,j,k,lnum) = lambda(ipoint)
                                    end select
                                end do
                            end do
                        end do
                    end do
                    !write(*,*) " AFTER 3 rg = ", rg
                end if
                !write(*,*) " AFTER 4 rg = ", rg
            end if
            !write(*,*) " AFTER 5 rg = ", rg
        end do

        !write(*,*) " AFTER 6 rg = ", rg

        !Adjusting PML values


        if(allocated(interpolatedRF)) deallocate(interpolatedRF)
        if(allocated(lambda))         deallocate(lambda)


        write(*,*) "  Propagating properties to PML in rg ", rg
        call propagate_PML_properties(Tdomain)

        !write(*,*) " AFTER 7 rg = ", rg

    end subroutine apply_prop_files

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
            write(*,*) "-----------------------------------------------------------------------"
            write(*,*) " "
            stop(" ")
        end if

        !INTEGERS
        attr_name = "nDim"
        call read_h5attr_int(file_id, trim(adjustL(attr_name)), nDim)
        attr_name = "Nmc"
        call read_h5attr_int(file_id, trim(adjustL(attr_name)), Nmc)

        !DOUBLE VEC
        attr_name = "xMinGlob"
        call read_h5attr_real_vec(file_id, attr_name, xMinGlob)
        attr_name = "xMaxGlob"
        call read_h5attr_real_vec(file_id, attr_name, xMaxGlob)
        attr_name = "xStep"
        call read_h5attr_real_vec(file_id, attr_name, xStep)

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

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine propagate_PML_properties(Tdomain)
        !This subroutine modifies the "prop" table.
        !In the PML elements it extrudes the properties according to the PML direction

        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain

        !LOCAL
        integer :: ngll, lnum
        integer :: iPoint
        integer :: i, j, k, n !counters
        integer :: LimPML1, LimPML2, LimPML3
        integer :: dir, mat_index
        double precision :: pointDens, pointLambda, pointMu;
        logical :: verbose = .false.

        !Propagating Properties over the PML

        ! With the new PML definition, we can do better:
        ! We can project x,y,z, on the pml plane, and interpolate the properties from the
        ! projected coordinates

        write(*,*) "Tdomain%not_PML_List = ", Tdomain%not_PML_List

        do n = 0, Tdomain%n_elem-1
            mat_index = Tdomain%specel(n)%mat_index

            if(Tdomain%not_PML_List(mat_index)) cycle

            ngll = Tdomain%sSubDomain(mat_index)%NGLL
            lnum = Tdomain%specel(n)%lnum
            !On the lower bound initialization
            LimPML1 = 0
            LimPML2 = 0
            LimPML3 = 0
            !To the upper bound transformation
            if (Tdomain%sSubDomain(mat_index)%pml_width(0)<0) LimPML1 = ngll-1
            if (Tdomain%sSubDomain(mat_index)%pml_width(1)<0) LimPML2 = ngll-1
            if (Tdomain%sSubDomain(mat_index)%pml_width(2)<0) LimPML3 = ngll-1

            dir = read_PML_Direction(Tdomain, mat_index)
            if(verbose) write(*,*) "dir = ", dir
            if(verbose) write(*,*) "Tdomain%specel(n)%domain = ", Tdomain%specel(n)%domain
            pointDens = 0

            !FOR TESTS
            !select case (Tdomain%specel(n)%domain)
            !    case (DM_SOLID_PML)
            !        Tdomain%spmldom%Density_(:,:,:,lnum) = 0
            !    case (DM_FLUID_PML)
            !
            !end select



            select case(dir)
                !Face X oriented
                case(0)
                    if(verbose) write(*,*) "Face X oriented"
                    do j = 0, ngll-1
                        do k = 0, ngll-1
                            select case (Tdomain%specel(n)%domain)
                                case (DM_SOLID_PML)
                                    pointDens   = Tdomain%spmldom%Density_(LimPML1,j,k,lnum)
                                    pointLambda = Tdomain%spmldom%Lambda_(LimPML1,j,k,lnum)
                                    pointMu     = Tdomain%spmldom%Mu_(LimPML1,j,k,lnum)
                                case (DM_FLUID_PML)

                            end select

                            do i = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,lnum) = pointDens
                                        Tdomain%spmldom%Lambda_(i,j,k,lnum)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,lnum)      = pointMu
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        enddo
                    enddo
                !
                !Face Y oriented
                case(1)
                    if(verbose) write(*,*) "Face Y oriented"
                    do i = 0, ngll-1
                        do k = 0, ngll-1
                            select case (Tdomain%specel(n)%domain)
                                case (DM_SOLID_PML)
                                    pointDens   = Tdomain%spmldom%Density_(i, LimPML2,k,lnum)
                                    pointLambda = Tdomain%spmldom%Lambda_(i, LimPML2,k,lnum)
                                    pointMu     = Tdomain%spmldom%Mu_(i, LimPML2,k,lnum)
                                case (DM_FLUID_PML)

                            end select

                            do j = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,lnum) = pointDens
                                        Tdomain%spmldom%Lambda_(i,j,k,lnum)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,lnum)      = pointMu
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        enddo
                    enddo

                !Face Z oriented
                case(2)
                    if(verbose) write(*,*) "Face Z oriented"
                    do i = 0, ngll-1
                        do j = 0, ngll-1
                            select case (Tdomain%specel(n)%domain)
                                case (DM_SOLID_PML)
                                    pointDens   = Tdomain%spmldom%Density_(i,j,LimPML3,lnum)
                                    pointLambda = Tdomain%spmldom%Lambda_(i,j,LimPML3,lnum)
                                    pointMu     = Tdomain%spmldom%Mu_(i,j,LimPML3,lnum)
                                case (DM_FLUID_PML)

                            end select

                            do k = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,lnum) = pointDens
                                        Tdomain%spmldom%Lambda_(i,j,k,lnum)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,lnum)      = pointMu
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        enddo
                    enddo

                !Edge in XY
                case(3)
                    if(verbose) write(*,*) "Edge in XY"
                    do k = 0, ngll-1
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_PML)
                                pointDens   = Tdomain%spmldom%Density_(LimPML1,LimPML2,k,lnum)
                                pointLambda = Tdomain%spmldom%Lambda_(LimPML1,LimPML2,k,lnum)
                                pointMu     = Tdomain%spmldom%Mu_(LimPML1,LimPML2,k,lnum)
                            case (DM_FLUID_PML)

                        end select

                        do i = 0, ngll-1
                            do j = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,lnum) = pointDens
                                        Tdomain%spmldom%Lambda_(i,j,k,lnum)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,lnum)      = pointMu
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        end do
                    enddo

                !Edge in YZ
                case(4)
                    if(verbose) write(*,*) "Edge in YZ"
                    do i = 0, ngll-1
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_PML)
                                pointDens   = Tdomain%spmldom%Density_(i, LimPML2,LimPML3,lnum)
                                pointLambda = Tdomain%spmldom%Lambda_(i, LimPML2,LimPML3,lnum)
                                pointMu     = Tdomain%spmldom%Mu_(i, LimPML2,LimPML3,lnum)
                            case (DM_FLUID_PML)

                        end select

                        do j = 0, ngll-1
                            do k = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,lnum) = pointDens
                                        Tdomain%spmldom%Lambda_(i,j,k,lnum)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,lnum)      = pointMu
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        end do
                    enddo

                !Edge in ZX
                case(5)
                    if(verbose) write(*,*) "Edge in ZX"
                    do j = 0, ngll-1
                        select case (Tdomain%specel(n)%domain)
                            case (DM_SOLID_PML)
                                pointDens   = Tdomain%spmldom%Density_(LimPML1, j,LimPML3,lnum)
                                pointLambda = Tdomain%spmldom%Lambda_(LimPML1, j,LimPML3,lnum)
                                pointMu     = Tdomain%spmldom%Mu_(LimPML1, j,LimPML3,lnum)
                            case (DM_FLUID_PML)

                        end select

                        do i = 0, ngll-1
                            do k = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,lnum) = pointDens
                                        Tdomain%spmldom%Lambda_(i,j,k,lnum)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,lnum)      = pointMu
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        end do
                    enddo

                !Vertex
                case(6)
                    if(verbose) write(*,*) "Vertex"
                    select case (Tdomain%specel(n)%domain)
                        case (DM_SOLID_PML)
                            pointDens   = Tdomain%spmldom%Density_(LimPML1, LimPML2,LimPML3,lnum)
                            pointLambda = Tdomain%spmldom%Lambda_(LimPML1, LimPML2,LimPML3,lnum)
                            pointMu     = Tdomain%spmldom%Mu_(LimPML1, LimPML2,LimPML3,lnum)
                        case (DM_FLUID_PML)

                    end select

                    do i = 0, ngll-1
                        do j = 0, ngll-1
                            do k = 0, ngll-1
                                select case (Tdomain%specel(n)%domain)
                                    case (DM_SOLID_PML)
                                        Tdomain%spmldom%Density_(i,j,k,lnum) = pointDens
                                        Tdomain%spmldom%Lambda_(i,j,k,lnum)  = pointLambda
                                        Tdomain%spmldom%Mu_(i,j,k,lnum)      = pointMu
                                    case (DM_FLUID_PML)

                                end select
                            end do
                        end do
                    end do
            end select
        end do

    end subroutine propagate_PML_properties

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    function read_PML_Direction(Tdomain, mat) result(dir)

        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: mat

        !OUTPUT
        integer :: dir

        !LOCAL
        integer :: error, code

        real :: vx, vy, vz

        !!! XXX: this wont work always, read comments from the call site
        vx = Tdomain%sSubDomain(mat)%pml_width(0)
        vy = Tdomain%sSubDomain(mat)%pml_width(1)
        vz = Tdomain%sSubDomain(mat)%pml_width(2)

        !/////////////Defining PML orientation
        dir = 0
        !Face X oriented
        if  (vx /= 0 .and. vy == 0d0 .and. vz == 0d0) then
            dir = 0
        !Face Y oriented
        elseif (vx == 0 .and. vy /= 0d0 .and. vz == 0d0) then
            dir = 1
        !Face Z oriented
        elseif (vx == 0 .and. vy == 0d0 .and. vz /= 0d0) then
            dir = 2
        !Edge in XY
        elseif (vx /= 0 .and. vy /= 0d0 .and. vz == 0d0) then
            dir = 3

        !Edge in YZ
        elseif (vx == 0 .and. vy /= 0d0 .and. vz /= 0d0) then
            dir = 4

        !Edge in ZX
        elseif (vx /= 0 .and. vy == 0d0 .and. vz /= 0d0) then
            dir = 5

        !Vertex in XYZ
        elseif (vx /= 0 .and. vy /= 0d0 .and. vz /= 0d0) then
            dir = 6

        !Undefined PML
        else
            write(*,*) "ERROR in mat ", mat, " (PML) definition (directions), check 'material.input'"
            call MPI_ABORT(Tdomain%communicateur, error, code)
        end if

    end function read_PML_Direction


end module build_prop_files
                     
