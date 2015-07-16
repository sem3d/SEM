program main_RandomField

	use randomFieldND
	use statistics_RF
	use obsolete_RF
	use writeResultFile_RF
	use readFile_RF
	use displayCarvalhol
	use mpi

    implicit none

    !INPUTS
    integer                                        :: Nmc, nDim;
    character (len=30), parameter                  :: mainInput        = "input_main"
    character (len=50), parameter                  :: sampleSpecFolder = "../sampleSpec"
    character (len=50), parameter                  :: meshFolder       = "../mesh"
    character (len=50), parameter                  :: outputFolder     = "../output"
    character (len=50)                             :: sampleSpecName, meshName, outputName;
    character (len=15)                             :: corrMod, margiFirst, meshType, meshMod;
    double precision,   dimension(:),  allocatable :: corrL;
    double precision                               :: fieldAvg, fieldVar;

	!OUTPUTS (in this case variables to be filled in in the proces)
	integer         , dimension(:)   , allocatable :: xNStep, kNStep;
    integer         , dimension(:, :), allocatable :: all_xNStep;
	double precision, dimension(:)   , allocatable :: kMax, xMax, xMin;
    double precision, dimension(:, :), allocatable :: randField, xPoints;
    double precision, dimension(:, :), allocatable :: all_RandField, all_xPoints;
    double precision, dimension(:, :), allocatable :: all_xMin, all_xMax;
    double precision, dimension(:)   , allocatable :: evntAvg, evntStdDev, procCorrL, globalCorrL;
    double precision, dimension(:)   , allocatable :: compEvntAvg, compEvntStdDev, compGlobCorrL;
    double precision                               :: globalAvg, globalStdDev, compGlobAvg, compGlobStdDev;

	!LOCAL VARIABLES
    logical           :: structured = .false., calcComp = .false.
    integer           :: i, baseStep, pointsPerCorrL;
    integer           :: code, rang, error, nb_procs;
    double precision  :: pi = 3.1415926535898;
    character(len=30) :: rangChar;
    character(len=30) , dimension(:,:), allocatable :: dataTable;
    character(len=100) :: path
    double precision  , dimension(:)  , allocatable :: sumRF, sumRFsquare, &
                                                      totalSumRF, totalSumRFsquare;

	call MPI_INIT(code)
	call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
	call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

	if(rang == 0) then
	    write(*,*) "------------START main-----------------------------------------";
		write(*,*) "";
	    write(*,*) "------------START MPI processing-----------------------------------------";
		write(*,*) "Number of procs = ", nb_procs
	    write(*,*) "";
		write(*,*) ">>>>>>>>> Reading file input";
	end if

	!Reading Filenames Input
	call set_DataTable("../"//mainInput, dataTable)
	!if(rang == 0) call DispCarvalhol (dataTable, "input_main");
	call read_DataTable(dataTable, "sampleSpecName" , sampleSpecName)
	call read_DataTable(dataTable, "meshName"       , meshName)
	call read_DataTable(dataTable, "outputName"     , outputName)
	deallocate(dataTable)

	!Reading Input
	path = trim(adjustL(sampleSpecFolder))//"/"//sampleSpecName
	call set_DataTable(path, dataTable)
	!if(rang == 0) call DispCarvalhol (dataTable, sampleSpecName);
	call read_DataTable(dataTable, "Nmc"        ,  Nmc)
	call read_DataTable(dataTable, "corrMod"    , corrMod)
	call read_DataTable(dataTable, "margiFirst" , margiFirst)
	call read_DataTable(dataTable, "corrL"      , corrL)
	call read_DataTable(dataTable, "fieldAvg"   , fieldAvg)
	call read_DataTable(dataTable, "fieldVar"   , fieldVar)
	deallocate(dataTable)

	nDim = size(corrL)

	!Reading Mesh
	path = trim(adjustL(meshFolder))//"/"//meshName
	call set_DataTable(path, dataTable)
	!if(rang == 0) call DispCarvalhol (dataTable, meshName);
	call read_DataTable(dataTable, "meshType", meshType)
	call read_DataTable(dataTable, "meshMod", meshMod)
	if (meshType == "structured" .or. meshType == "unstructured") then
		select case (meshMod)
			case("manual")
				write(rangChar,'(I7)') rang
				call read_DataTable(dataTable, "Max"//trim(adjustl(rangChar)), xMax)
				call read_DataTable(dataTable, "Min"//trim(adjustl(rangChar)), xMin)
			case("automatic")
				baseStep = nint(dble(nb_procs)**(1.0d0/nDim))
				if (nb_procs == baseStep**nDim) then
					call read_DataTable(dataTable, "Max", xMax)
					call read_DataTable(dataTable, "Min", xMin)
					call set_ProcMaxMin(xMin, xMax, rang, nb_procs)
				else
					write(*,*) "In meshing automatic mode the number of processors should be an integer like (N)**nDim, where N is an integer and nDim is the number of dimensions"
					call MPI_ABORT(MPI_COMM_WORLD, error, code)
				end if
		end select
		nDim = size(xMax)
	end if
	deallocate(dataTable)

	!Showing read inputs
	if(rang == 0) then
		write(*,*) ""
		write(*,*) ">>>>>>>>INPUT Read Data"
		write(*,*) "##IN/OUT"
		write(*,*) "sampleSpecName = ", sampleSpecName
		write(*,*) "outputName     = ", outputName
		write(*,*) ""
		write(*,*) "##MESH"
		write(*,*) "meshName       = ", meshName
		write(*,*) "meshType       = ", meshType
		write(*,*) "meshMod        = ", meshMod
		write(*,*) "xMin           = ", xMin
		write(*,*) "xMax           = ", xMax
		write(*,*) ""
		write(*,*) "##EVENTS"
		write(*,*) "nDim           = ", nDim
		write(*,*) "Nmc            = ", Nmc !number of Monte-Carlo experiments
		write(*,*) "corrMod        = ", corrMod
		write(*,*) "margiFirst     = ", margiFirst
		write(*,*) "corrL          = ", corrL
		write(*,*) "fieldAvg       = ", fieldAvg
		write(*,*) "fieldVar       = ", fieldVar
	end if

	!Input Validation
	if(Nmc < 1) then
		write(*,*) ""
		write(*,*) "ERROR - Number of events should be a positive integer"
		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if
	if(nDim < 1) then
		write(*,*) ""
		write(*,*) "ERROR - nDim should be a positive integer"
		write(*,*) "nDim           = ", nDim
		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if
	do i = 1, nDim
		if(corrL(i) < 0.001) then
			write(*,*) ""
			write(*,*) "ERROR - corrL should be a positive number greater than 0.001"
			write(*,*) "corrL ", i, "of proc ", rang, " is ", corrL(i)
			call MPI_ABORT(MPI_COMM_WORLD, error, code)
		end if
	end do

	!--------------------------------------------------------------------------------------
	!------------------------------UNSTRUCTURED--------------------------------------------
	!--------------------------------------------------------------------------------------

	if(meshType == "unstructured") then

		!Creating x points for each processor (in the future this would be one of the inputs)
		write(*,*) ""
		if(rang == 0) write(*,*) ">>>>>>>>> Creating xPoints for unstructured mesh";
		pointsPerCorrL = 10
		call set_XPoints(corrL, xMin, xMax, xPoints, pointsPerCorrL)
		!write(*,*) " Rang = ", rang
		!call dispCarvalhol(xPoints, "xPoints")

		!Calculating a part of the random field in each processor
		if(rang == 0) write(*,*) ">>>>>>>>> Calculating Random field (unstructured)";
		call createRandomField(xPoints, corrMod, margiFirst, corrL, fieldAvg, fieldVar, Nmc, randField);
		!call dispCarvalhol(xPoints, "xPoints", mpi = .true.)
		!call dispCarvalhol(randField, "randField", mpi = .true.)

		if(rang == 0) write(*,*) ">>>>>>>>> Writing Result File";
		path = trim(adjustL(outputFolder))//"/"//outputName
		call write_ResultHDF5(xPoints, randField, path, rang) !HDF5 of this processor

		!Calculating Statistics
		if(rang == 0) write(*,*) ">>>>>>>>> Calculating Statistics MPI (unstructured)";
!		call set_Statistics_MPI(randField, xPoints, rang,             &
!							 	evntAvg, evntStdDev, procCorrL,       &
!							 	globalAvg, globalStdDev, globalCorrL)

		!Writing results and comparison statistics (To be deleted)
		if(rang == 0) write(*,*) ">>>>>>>>> Calculating Comparison Statistics (unstructured)";
		call set_allRandField(randField, xPoints, rang,    &
						      all_RandField, all_xPoints)
		!if(rang == 0) call DispCarvalhol(all_xPoints, "all_xPoints")
		!if(rang == 0) call DispCarvalhol(all_RandField, "all_RandField")
		if(rang == 0) then
			call set_CompStatistics(all_RandField, all_xPoints,                   &
			                        compEvntAvg, compEvntStdDev,                  &
			                        compGlobAvg, compGlobStdDev, compGlobCorrL)
		end if

		if (rang == 0 .and. nDim == 1) then
			path = trim(adjustL(outputFolder))//"/"//outputName
			call write_MatlabTable (all_RandField, path)
		end if

	!--------------------------------------------------------------------------------------
	!--------------------------------STRUCTURED (need to be updated)--------------------------------------------
	!--------------------------------------------------------------------------------------

	else if(meshType == "structured") then
		!Calculating a part of the random field in each processor
		if(rang == 0) write(*,*) ">>>>>>>>> Calculating Random field (structured)";
		call createRandomField(xMin, xMax, corrMod, corrL, Nmc, xNStep, randField);
		!Calculating Statistics
		if(rang == 0) write(*,*) ">>>>>>>>> Calculating Statistics MPI (structured)";
		call set_Statistics_MPI(randField, xMin, xMax, xNStep, rang, &
							 	evntAvg, evntStdDev, procCorrL,      &
							 	globalAvg, globalStdDev, globalCorrL)

		if(rang == 0) write(*,*) ">>>>>>>>> Writing Result File";
		path = trim(adjustL(outputFolder))//"/"//outputName
		call write_ResultHDF5(xMin, xMax, xNStep, randField, path, rang) !HDF5 of this processor

		!Writing results and comparison statistics (To be deleted)
		if (calcComp) then
			if(rang == 0) write(*,*) ">>>>>>>>> Calculating Comparison Statistics (structured)";
			call set_allRandField(randField, xMin, xMax, xNStep, rang,    &
							      all_xMin, all_xMax, all_xNStep,         &
							      all_RandField, all_xPoints)
			if(rang == 0) then
				call set_CompStatistics(all_RandField, all_xPoints,                   &
										all_xMin, all_xMax, all_xNStep,               &
				                        compEvntAvg, compEvntStdDev,                  &
				                        compGlobAvg, compGlobStdDev, compGlobCorrL)
			end if
		end if

	!--------------------------------------------------------------------------------------
	!--------------------------------------------------------------------------------------
	!--------------------------------------------------------------------------------------

	else
		write(*,*) "ERROR - meshtype not accepted"
		write(*,*) "meshtype = ", meshtype
		call MPI_ABORT(MPI_COMM_WORLD, error, code)
	end if

	if(rang == 0) then
		write(*,*) "Showing Statistics------------------------------- "
		write(*,*) ""
		if(Nmc < 11) then
			write(*,*) ">>BY EVENT"
			call DispCarvalhol(evntAvg, "    MPI Avg = ")
			if (calcComp) call DispCarvalhol(compEvntAvg, "   Comp Avg = ")
			call DispCarvalhol(evntStdDev**2, "    MPI Var = ")
			if (calcComp) call DispCarvalhol(compEvntStdDev**2, "   Comp Var = ")
		end if

		write(*,*) ">>GLOBAL"
		write(*,*) "  INPUT Avg = ", fieldAvg
		if (calcComp) write(*,*) "   Comp Avg = ", compGlobAvg;
		write(*,*) "    MPI Avg = ", globalAvg;
		write(*,*) " "
		write(*,*) "  INPUT Var = ", fieldVar
		if (calcComp) write(*,*) "   Comp Var = ", compGlobStdDev**2;
		write(*,*) "    MPI Var = ", globalStdDev**2;
		write(*,*) " "
		write(*,*) "INPUT CorrL = ", corrL
		if (calcComp) write(*,*) " Comp CorrL = ", compGlobCorrL;
		write(*,*) "  MPI CorrL = ", globalCorrL;
!		call DispCarvalhol(globalCorrL, "  MPI CorrL ")
!		call DispCarvalhol(compGlobCorrL, " Comp CorrL = ")
	end if


	!Deallocating
		if(allocated(corrL))   deallocate(corrL);
		if(allocated(xMax))    deallocate(xMax);
		if(allocated(xMin))    deallocate(xMin);
		if(allocated(kMax))      deallocate(kMax);
		if(allocated(xNStep))    deallocate(xNStep);
		if(allocated(kNStep))    deallocate(kNStep);
		if(allocated(xPoints))   deallocate(xPoints);
		if(allocated(randField)) deallocate(randField);

		if(allocated(evntAvg))      deallocate(evntAvg);
		if(allocated(evntStdDev))   deallocate(evntStdDev);
		if(allocated(procCorrL))    deallocate(procCorrL);
		if(allocated(globalCorrL))  deallocate(globalCorrL);
		if(allocated(sumRF))        deallocate(sumRF)
		if(allocated(sumRFsquare))  deallocate(sumRFsquare)
		if(allocated(totalSumRF))        deallocate(totalSumRF)
		if(allocated(totalSumRFsquare))  deallocate(totalSumRFsquare)
		if(allocated(all_RandField))     deallocate(all_RandField);
		if(allocated(all_xPoints))       deallocate(all_xPoints);
		if(allocated(all_xNStep)) deallocate(all_xNStep)
		if(allocated(all_xMin))   deallocate(all_xMin)
		if(allocated(all_xMax))   deallocate(all_xMax)

		if(allocated(compEvntAvg))    deallocate(compEvntAvg);
		if(allocated(compEvntStdDev)) deallocate(compEvntStdDev);
		if(allocated(compGlobCorrL))  deallocate(compGlobCorrL);

	if(rang == 0) then
		write(*,*) "";
	    write(*,*) "------------END main-----------------------";
		write(*,*) "";
	end if

	call MPI_FINALIZE(code)

end program main_RandomField
