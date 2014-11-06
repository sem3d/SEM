module writeResultFile_RF

    use displayCarvalhol
    use statistics_RF
    use math_RF
    use hdf5
    use mpi

    interface write_ResultHDF5
		module procedure write_ResultHDF5Structured,   &
		                 write_ResultHDF5Unstruct_MPI
		                 !write_ResultHDF5Unstruct
	end interface write_ResultHDF5

contains


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine write_ResultHDF5Structured(xMin, xMax, xNStep, randField, fileName, rang)

        implicit none

        !INPUTS
        double precision,  dimension(:),   intent(in) :: xMax, xMin;
        integer,           dimension(:),   intent(in) :: xNStep;
        double precision,  dimension(:,:), intent(in) :: randField;
        character (len=*)                , intent(in) :: filename;
        integer                          , intent(in) :: rang;

		!LOCAL VARIABLES
        double precision, dimension(:,:), allocatable :: xPoints
        integer :: nDim, nPoints, i

        nDim         = size(xNStep , 1)
        nPoints      = size(randField, 1)

        allocate(xPoints(nPoints, nDim))

		do i = 1, nPoints
				call get_Permutation(i, xMax, xNStep, xPoints(i,:), xMin);
		end do

		!call write_ResultHDF5Unstruct_MPI(xPoints, randField, fileName, rang)


    end subroutine write_ResultHDF5Structured

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
        integer              , optional  , intent(in) :: communicator
        character(len=*) , dimension(1:), optional  , intent(in) :: labels
        integer          , dimension(1:), optional  , intent(in) :: indexes

        !OUTPUTS
		character(len=110) , optional  , intent(out) ::HDF5Name

        !HDF5 VARIABLES
        character(len=110)             :: fileHDF5Name, fullPath !File name
        character(len=30)              :: eventName, coordName;  !Dataset names
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
        double precision, dimension(:,:), allocatable :: grid_data

!		if(rang == 0) then
!      	  	write(*,*) "";
!       		write(*,*) "------------START Writing result HDF5 file (MPI)-----------------------";
!      	 	write(*,*) "";
!      	end if

      	!if(rang == 0) then
      		!write (*,*) "lbound(xPoints) = ", lbound(xPoints)
      		!write (*,*) "lbound(randField) = ", lbound(randField)
      		!write (*,*) "size(xPoints,1) = ", size(xPoints,1)
      		!write (*,*) "size(xPoints,2) = ", size(xPoints,2)
      		!write (*,*) "xPoints(1, 1) = ", xPoints(1, 1)
      		!write (*,*) "xPoints(1:10, :) = ", xPoints(1:10, :)
      		!call dispCarvalhol(xPoints(:,:)  , "xPoints(:,:)"  , "F30.5")
      		!call dispCarvalhol(randField(:,:), "randField(:,:)", "F30.5")
      	!end if


		effectComm = communicator
	    nDim       = size(xPoints , 2)
	    nPoints    = size(randField, 1)
	    Nmc        = size(randField, 2)

		!Creating file name
		if(.not. present(labels)) then
			write(rangStr,'(I)'  ) rang
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

			!if(rang == 0) write(*,*) ">>>>>>>>> Opening file";
        	allocate (grid_data(3, nPoints)) !3 lines to put X, Y and Z
        	grid_data = 0;
        	grid_data(1:nDim, :) = transpose(xPoints)
	        call h5open_f(error) ! Initialize FORTRAN interface.
	        call h5fcreate_f(fullPath, H5F_ACC_TRUNC_F, file_id, error) ! Create a new file using default properties.

			!if(rang == 0) write(*,*) ">>>>>>>>> Creating Coordinates dataset 'XYZ table'";

	        dims = shape(grid_data)
			write(coordName,'(A)') "XYZ"
			!if(Nmc < 11) write(*,*) "coordName = ", coordName

	        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
	        call h5dcreate_f(file_id, coordName, H5T_NATIVE_DOUBLE,  &
	                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).

	 		call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, grid_data, dims, error) ! Write the dataset.
	        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
	        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.

	        deallocate (grid_data)

			!if(rang == 0) write(*,*) ">>>>>>>>> Creating Quantities dataset 'random field'";
	        dims(1) = size(randField,1)
	        dims(2) = 1 !One random field in each dataset

			do i = 1, Nmc
				write(numberStr,'(I)'  ) i
				numberStr = adjustL(numberStr)
				write(eventName,'(2A)') "RF_", trim(numberStr)
				!if(Nmc < 11) write(*,*) "eventName = ", eventName

		        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
		        call h5dcreate_f(file_id, eventName, H5T_NATIVE_DOUBLE,  &
		                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
		        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, randField(:,i), dims, error) ! Write the dataset.
		        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
		        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.
		    end do

			!if(rang == 0) write(*,*) ">>>>>>>>> Closing file";
	        call h5fclose_f(file_id, error) ! Close the file.
	        call h5close_f(error) ! Close FORTRAN interface.

        end if

		if(present(HDF5Name)) HDF5Name = fileHDF5Name

!		if(rang == 0) then
!        	write(*,*) "";
!        	write(*,*) "------------END Writing result HDF5 file (MPI)-----------------------";
!       		write(*,*) "";
!       	end if

    end subroutine write_ResultHDF5Unstruct_MPI


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine writeXMF_RF_MPI(nSamples, HDF5nameList, nPointList, fileName, rang, folderPath, &
    						   communicator, HDF5relativePath, attName, byProc)
        implicit none

        !INPUTS
        integer                           , intent(in) :: nSamples;
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
        character (len=35)  :: eventName, meshName;
        character (len=50) , dimension(:), allocatable :: effectAttName;
        integer            , dimension(:), allocatable :: all_nPointList
        character (len=110), dimension(:), allocatable :: all_HDF5nameList

        effectComm = communicator

!        if(rang == 0) then
!        	write(*,*) "";
!        	write(*,*) "------------START Writing result XMF file-----------------------";
!        	write(*,*) "";
!        end if

		call MPI_COMM_SIZE(effectComm, nb_procs, code)

		if(rang == 0) then
			allocate(all_nPointList(nb_procs*size(HDF5nameList)))
			allocate(all_HDF5nameList(nb_procs*size(HDF5nameList)))
		end if

		call MPI_GATHER(nPointList    , size(nPointList), MPI_INTEGER,     &
		                all_nPointList, size(nPointList), MPI_INTEGER,     &
		                 0         , effectComm    , code)

		call MPI_GATHER(HDF5nameList    , len(HDF5nameList)*size(HDF5nameList), MPI_CHARACTER,     &
		                all_HDF5nameList, len(HDF5nameList)*size(HDF5nameList), MPI_CHARACTER,     &
		                 0              , effectComm       , code)

		!write(*,*) "all_HDF5nameList = ", all_HDF5nameList

		if(rang == 0) then
			!Common parameters
			Nmc         = nSamples
			fileXMFName = string_join(fileName,".xmf")
			fullPathXMF = string_join(folderPath, "/"//fileXMFName)
			!if(rang == 0) write(*,*) "all_nPointList in rang 0 = ", all_nPointList
			!if(rang == 0) write(*,*) "all_HDF5nameList in rang 0 = ", all_HDF5nameList
			if(rang == 0) write(*,*) "fileXMFName = ", fileXMFName
			if(rang == 0) write(*,*) "fullPathXMF = ", fullPathXMF

			!Optional inputs
			allocate(effectAttName(Nmc))
        	if(present(attName) .and. (size(attName)==Nmc)) then
        		effectAttName = attName
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
					!write(numberStr,'(I)'  ) j
					!numberStr = adjustL(numberStr)
					!write(meshName,'(2A)' ) "meshRF_", trim(adjustL(numberStr))
					meshName = stringNumb_join("Subdomain_", j-1)
					if(present(byProc)) then
						if(byProc) meshName = stringNumb_join("Proc_", int((j-1)/size(HDF5nameList)))
						!write(*,*) "                       j-1 = ",j-1
						!write(*,*) "        size(HDF5nameList) = ",size(HDF5nameList)
						!write(*,*) "int((j-1)/size(HDF5nameList) = ",int((j-1)/size(HDF5nameList))
					end if
				write (file,'(3A)'     )'   <Grid Name="',trim(meshName),'" GridType="Uniform">' !START Writing the data of one subdomain
				write (file,'(3A)'     )'    <Topology Type="Polyvertex" NodesPerElements="1" NumberOfElements="',trim(numb2String(all_nPointList(j))),'">'
				write (file,'(A)'      )'    </Topology>'
				write (file,'(A)'      )'     <Geometry GeometryType="XYZ">'
				write (file,'(3A)'     )'      <DataItem Name="Coordinates" Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(numb2String(all_nPointList(j))), ' 3">'
				write (file,'(4A)'     )'     	  ',trim(effecHDF5path),trim(all_HDF5nameList(j)),':/XYZ'
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
				end do !END Loop Over Subdomains

				write (file,'(A)'      )'  </Grid>' !Closes the Collection
				write (file,'(A)'      )' </Domain>'
				write (file,'(A)'      )'</Xdmf>'

	        close(file)
        end if

		if(rang == 0) then
			deallocate(all_nPointList)
			deallocate(all_HDF5nameList)
			deallocate(effectAttName)

!        	write(*,*) "";
!        	write(*,*) "------------END Writing result XMF file-----------------------";
!        	write(*,*) "";
        end if

    end subroutine writeXMF_RF_MPI

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine writeResultFileND(corrMod, corrL, xMax, xPeriod,   &
							     kMax, xNStep, kNStep, randField, &
							     average, stdDeviation, averageCorrL, fileName, time)

        implicit none
        !INPUT
        character (len=*)                             :: corrMod;
        double precision,  dimension(:),   intent(in) :: corrL, xMax, xPeriod, kMax;
        integer,           dimension(:),   intent(in) :: xNStep, kNStep;
        double precision,  dimension(:,:), intent(in) :: randField;
        double precision,  dimension(:),   intent(in) :: average, stdDeviation, averageCorrL;
        character (len=*)              ,   intent(in) :: filename;
        double precision,  optional,       intent(in) :: time;

        !OUTPUT - Result File

        !LOCAL VARIABLES
        integer            :: i, j, file, nLines, nColumns, nDim;
        character (len=40) :: titleFmt, stringFmt, doubleFmt, intFmt, dimFmt;
        character (len=50) :: path
        character (len=20), dimension(:), allocatable :: eventLabel;

        write(*,*) "";
        write(*,*) "------------START Writing result file-----------------------";
        write(*,*) "";

		path = trim(adjustL(filename));

        nLines   = size(randField,1);
        nColumns = size(randField,2);
        nDim     = size(xMax,     1);
        titleFmt = "A18"

        allocate(eventLabel(nColumns))

        write(stringFmt, fmt="(A,i10,A)") "(t25,", nColumns, "A25)";
        write(intFmt,    fmt="(A,i10,A)") "(t25,", nColumns, "I25)";
        write(doubleFmt, fmt="(A,i10,A)") "(t25,", nColumns, "F25.16)";
        write(dimFmt,    fmt="(A,i10,A)") "(t25,", nDim, "F25.16)";

        do i = 1, nColumns
            write(eventLabel(i), fmt="(A10,i10)") "Event", i;
        end do

        !>>>>>>>>>>>>>>>>>>>>
        file=10;
        open (unit = file , file = path, action = 'write')

        write (file, *) "DIMENSION = ", size(xMax,1);
        write (file, *) "";

        write (file,"("//titleFmt//","//stringFmt//")") "", eventLabel;

        write (file, *) "INPUTS";
        write (file,"("//titleFmt//","//stringFmt//")") "corrMod", adjustr(corrMod);
        write (file,"("//titleFmt//","//doubleFmt//")") "corrL",   corrL;
        write (file,"("//titleFmt//","//doubleFmt//")") "xMax",    xMax;
        write (file,"("//titleFmt//","//doubleFmt//")") "xPeriod", xPeriod;

        write (file, *) "";
        write (file, *) "COMPUTED";
        write (file,"("//titleFmt//","//doubleFmt//")") "kMax",   kMax;
        write (file,"("//titleFmt//","//intFmt   //")") "xNStep", xNStep;
        write (file,"("//titleFmt//","//intFmt   //")") "kNStep", kNStep;

		write(*,*) ">>>>>>>>> Writing Random Field";
        write (file, *) "";
        write (file, *) "RANDOM FIELD";
        write (file,"("//doubleFmt//")") ((randField(i,j),j=1, size(randField,2)),i=1, size(randField,1));;

		write(*,*) ">>>>>>>>> Writing Statistics";
        write (file, *) "";
        write (file, *) "STATISTICS - POINT BY POINT";
        write (file,"("//titleFmt//","//stringFmt//","//stringFmt//")") "", "average", "stdDeviation";
        do i = 1, size(average)
            write (file,"("//titleFmt//","//"(t25,2F25.16)"//")") "", average(i), stdDeviation(i);
        end do

!        write (file, *) "STATISTICS - GLOBAL";
!        write (file, *) "";
!        write (file,"("//titleFmt//","//"(t25, F25.16)"//")") "average = ",   calculateAverage(average)
!        write (file,"("//titleFmt//","//"(t25, F25.16)"//")") "stdDeviation = ", calculateAverage(stdDeviation);
!        write (file,"("//titleFmt//","//dimFmt//")") "avgCorrL = ", averageCorrL;

        if(present(time)) write (file,"("//titleFmt//","//"(t25, F25.16)"//")") "Total time (s) = ", time;
        close(file)

    end subroutine writeResultFileND

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

        write(*,*) "";
        write(*,*) "------------START Writing result file-----------------------";
        write(*,*) "";

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

    	!write(*,*) "WRITE Flag string_join"

    	stringTot = trim(adjustL(string1))//trim(adjustL(string2))

    end function string_join

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

    	!write(*,*) "WRITE Flag stringNumb_join"

    	write(nmbString, fmt='(I)') number
    	stringTot = string_join(string, nmbString)

    	!write(*,*) "WRITE Flag 2 stringNumb_join"

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

    	write(nmbString, fmt='(I)') number
    	stringTot = trim(adjustL(nmbString))

    end function
!TRASH

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    subroutine write_ResultHDF5Unstruct(xPoints, randField, fileName)
!        implicit none
!
!        !INPUTS
!        double precision,  dimension(:,:), intent(in) :: randField, xPoints;
!        character (len=*),                 intent(in) :: filename;
!
!        !HDF5 VARIABLES
!        character(len=50)              :: fileHDF5Name !File name
!        character(len=4),  parameter   :: dsetXYZ = "XYZ", dsetRF  = "RF"  !Dataset name
!        integer(HID_T)                 :: file_id       !File identifier
!        integer(HID_T)                 :: dset_id       !Dataset identifier
!        integer(HID_T)                 :: dspace_id     !Dataspace identifier
!        integer                        :: rank = 2      !Dataset rank (number of dimensions)
!        integer(HSIZE_T), dimension(2) :: dims !Dataset dimensions
!        integer                        :: error !Error flag
!
!		!LOCAL VARIABLES
!        integer :: nDim, Nmc, nPoints, i
!        character (len=50) :: eventName;
!        character (len=12) :: numberStr;
!        double precision, dimension(:,:), allocatable :: grid_data
!
!        write(*,*) "";
!        write(*,*) "------------START Writing result HDF5 file-----------------------";
!        write(*,*) "";
!
!        fileHDF5Name = trim(fileName)//".h5"
!        nDim         = size(xPoints , 2)
!        nPoints      = size(randField, 1)
!        Nmc          = size(randField, 2)
!
!        write(*,*) "fileHDF5Name", fileHDF5Name
!
!
!        if (nDim > 3) then
!        	write(*,*) "Dimension exceeds 3, HDF file won't be created"
!
!        else
!
!        	allocate (grid_data(3, nPoints)) !3 lines to put X, Y and Z
!        	grid_data = 0;
!        	grid_data(1:nDim, :) = transpose(xPoints)
!
!			write(*,*) ">>>>>>>>> Opening file";
!	        call h5open_f(error) ! Initialize FORTRAN interface.
!	        call h5fcreate_f(fileHDF5Name, H5F_ACC_TRUNC_F, file_id, error) ! Create a new file using default properties.
!
!			write(*,*) ">>>>>>>>> Creating Coordinates dataset 'XYZ table'";
!
!	        dims = shape(grid_data)
!	        !write(*,*) "dims   = ", dims
!	        !write(*,*) dsetXYZ
!
!	        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
!	        call h5dcreate_f(file_id, dsetXYZ, H5T_NATIVE_DOUBLE,  &
!	                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
!
!	 		call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, grid_data, dims, error) ! Write the dataset.
!	        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
!	        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.
!	        deallocate (grid_data)
!
!			write(*,*) ">>>>>>>>> Creating Quantities dataset 'random field'";
!	        dims(1) = size(randField,1)
!	        dims(2) = 1 !One random field in each dataset
!			!write(*,*) "dims   = ", dims;
!
!			do i = 1, Nmc
!				write(numberStr,'(I)'  ) i
!				numberStr = adjustL(numberStr)
!				write(eventName,'(A,A)') "RF_", numberStr
!
!		        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
!		        call h5dcreate_f(file_id, eventName, H5T_NATIVE_DOUBLE,  &
!		                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
!		        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, randField(:,i), dims, error) ! Write the dataset.
!		        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
!		        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.
!		    end do
!
!			write(*,*) ">>>>>>>>> Closing file";
!	        call h5fclose_f(file_id, error) ! Close the file.
!	        call h5close_f(error) ! Close FORTRAN interface.
!	        call writeXMF_RF(randField, fileHDF5Name, fileName)
!
!        end if
!
!        write(*,*) "";
!        write(*,*) "------------END Writing result HDF5 file-----------------------";
!        write(*,*) "";
!
!    end subroutine write_ResultHDF5Unstruct


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    subroutine writeXMF_RF(randField, fileHDF5Name, fileName)
!        implicit none
!
!        !INPUTS
!        double precision,  dimension(:,:), intent(in) :: randField;
!        character (len=*),                 intent(in) :: filename;
!        character (len=*),                 intent(in) :: fileHDF5Name;
!
!		!LOCAL VARIABLES
!        integer            :: Nmc, i, file, nElem
!        character (len=35) :: fileXMFName;
!        character (len=35) :: eventName;
!        character (len=12) :: numberStr, nElemStr, NmcStr;
!
!        write(*,*) "";
!        write(*,*) "------------START Writing result XMF file-----------------------";
!        write(*,*) "";
!
!		fileXMFName = trim(fileName)//".xmf"
!        Nmc     = size(randField, 2)
!        nElem   = size(randField, 1)
!        write(NmcStr,'(I)'  ) Nmc
!		write(nElemStr,'(I)') nElem
!		NmcStr   = adjustL(NmcStr)
!		nElemStr = adjustL(nElemStr)
!
!        write(*,*) "fileXMFName", fileXMFName
!        !write(*,*) "nElem = ", nElemStr
!		!write(*,*) "Nmc   = ", NmcStr
!
!        file=21;
!
!        open (unit = file , file = trim(fileXMFName), action = 'write')
!
!			write (file,'(A)'    )'<?xml version="1.0" ?>'
!			write (file,'(A)'    )'<Xdmf Version="2.0">'
!			write (file,'(A)'    )' <Domain>'
!			write (file,'(A)'    )'  <Grid Name="meshRF" GridType="Uniform">'
!			write (file,'(A,A,A)')'   <Topology Type="Polyvertex" NodesPerElements="1" NumberOfElements="',trim(nElemStr),'">'
!			write (file,'(A)'    )'   </Topology>'
!			write (file,'(A)'    )'    <Geometry GeometryType="XYZ">'
!			write (file,'(A,A,A)')'     <DataItem Name="Coordinates" Format="HDF" DataType="Float" Precision="8" Dimensions="',trim(nElemStr), ' 3">'
!			write (file,'(A,A,A)')'     	 ',trim(fileHDF5Name),':/XYZ'
!			write (file,'(A)'    )'     </DataItem>'
!			write (file,'(A)'    )'    </Geometry>'
!			do i = 1, Nmc
!
!
!				write(numberStr,'(I)'  ) i
!				numberStr = adjustL(numberStr)
!				write(eventName,'(A,A)') "RF_", numberStr
!
!				write (file,'(A,A,A)'  )'    <Attribute Name="',trim(eventName),'" Center="Node" AttributeType="Scalar">'
!				write (file,'(A,A,A)'  )'      <DataItem Dimensions="',trim(nElemStr),'" Format="HDF" DataType="Float" Precision="8">'
!				write (file,'(A,A,A,A)')'         ',trim(fileHDF5Name),":/", eventName
!				write (file,'(A)'      )'      </DataItem>'
!				write (file,'(A)'      )'    </Attribute>'
!			end do
!
!			write (file,'(A)'    )'  </Grid>'
!			write (file,'(A)'    )' </Domain>'
!			write (file,'(A)'    )'</Xdmf>'
!
!        close(file)
!
!
!
!        write(*,*) "";
!        write(*,*) "------------END Writing result XMF file-----------------------";
!        write(*,*) "";
!
!    end subroutine writeXMF_RF


!iN "writeXMF_RF"
! START - VERSION TESTED FOR Nmc = 1
!        open (unit = file , file = trim(fileXMFName), action = 'write')
!
!			write (file,'(A)'    )'<?xml version="1.0" ?>'
!			write (file,'(A)'    )'<Xdmf Version="2.0">'
!			write (file,'(A)'    )' <Domain>'
!			write (file,'(A)'    )'  <Grid Name="meshRF" GridType="Uniform">'
!			write (file,'(A,I,A)')'   <Topology Type="Polyvertex" NodesPerElements="1" NumberOfElements="',nElem,'">'
!			write (file,'(A)'    )'   </Topology>'
!			write (file,'(A)'    )'    <Geometry GeometryType="XYZ">'
!			write (file,'(A,I,A)')'     <DataItem Name="Coordinates" Format="HDF" DataType="Float" Precision="8" Dimensions="',nElem, ' 3">'
!			write (file,'(A,A,A)')'     	 ',trim(fileHDF5Name),':/XYZ'
!			write (file,'(A)'    )'     </DataItem>'
!			write (file,'(A)'    )'    </Geometry>'
!			write (file,'(A)'    )'    <Attribute Name="randomField" Center="Node" AttributeType="Scalar">'
!			write (file,'(A,I,A)')'     <DataItem Dimensions="',nElem,'" Format="HDF" DataType="Float" Precision="8">'
!			write (file,'(A,A,A)')'        ',trim(fileHDF5Name),':/RF'
!			write (file,'(A)'    )'     </DataItem>'
!			write (file,'(A)'    )'    </Attribute>'
!			write (file,'(A)'    )'  </Grid>'
!			write (file,'(A)'    )' </Domain>'
!			write (file,'(A)'    )'</Xdmf>'
!
!        close(file)
! END - VERSION TESTED FOR Nmc = 1

!iN "writeRsultHDF5"
! START - VERSION TESTED FOR Nmc = 1
!		if (nDim > 3) then
!        	write(*,*) "Dimension exceeds 3, HDF file won't be created"
!
!        else
!
!			write(*,*) ">>>>>>>>> Opening file";
!	        call h5open_f(error) ! Initialize FORTRAN interface.
!	        call h5fcreate_f(fileHDF5Name, H5F_ACC_TRUNC_F, file_id, error) ! Create a new file using default properties.
!
!			write(*,*) ">>>>>>>>> Creating Coordinates dataset 'XYZ table'";
!
!			allocate (grid_data(3, nPoints)) !3 lines to put X, Y and Z
!	        grid_data = 0D0
!
!	        dims = shape(grid_data)
!	        write(*,*) "dims   = ", dims
!
!	        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
!	        call h5dcreate_f(file_id, dsetXYZ, H5T_NATIVE_DOUBLE,  &
!	                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
!
!
!
!			do i = 1, nPoints
!				call get_Permutation(i, xMax, xNStep, grid_data(1:nDim,i));
!			end do
!
!	 		call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, grid_data, dims, error) ! Write the dataset.
!	        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
!	        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.
!
!			deallocate (grid_data)
!
!			write(*,*) ">>>>>>>>> Creating Quantities dataset 'random field'";
!	        dims = shape(randField)
!
!			write(*,*) "dims   = ", dims;
!
!	        call h5screate_simple_f(rank, dims, dspace_id, error) ! Create the dataspace (dspace_id).
!	        call h5dcreate_f(file_id, dsetRF, H5T_NATIVE_DOUBLE,  &
!	                         dspace_id, dset_id, error) ! Create the dataset with default properties (dset_id).
!	        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, randField, dims, error) ! Write the dataset.
!	        call h5dclose_f(dset_id, error) ! End access to the dataset and release resources used by it.
!	        call h5sclose_f(dspace_id, error) ! Terminate access to the data space.
!
!			write(*,*) ">>>>>>>>> Closing file";
!	        call h5fclose_f(file_id, error) ! Close the file.
!	        call h5close_f(error) ! Close FORTRAN interface.
!
!	        call writeXMF_RF(randField, fileHDF5Name, fileName)
!
!        end if
! END - VERSION TESTED FOR Nmc = 1


end module writeResultFile_RF
