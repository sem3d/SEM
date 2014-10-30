module obsolete_RF

	use displayCarvalhol
	use math_RF
	use mpi

	interface set_allRandField
		module procedure set_allRandFieldStructured,   &
		                 set_allRandFieldUnstruct
	end interface set_allRandField

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine get_sizes_MPI(xNStep, sizeLoc, sizeUnif, start, end)

		implicit none
    	!INPUT
    	integer, dimension(:), intent(in) :: xNStep;
    	!OUTPUT
        integer, intent(out), optional :: sizeLoc, sizeUnif;
        integer, intent(out), optional :: start, end;

        !LOCAL VARIABLES
        integer :: rang, nb_procs, code;
        integer :: xStart, xEnd, xNStepTotal;

    	call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

		xNStepTotal = product(xNStep)

    	if(rang == nb_procs-1) then
			xStart = (nb_procs-1)*ceiling(dble(xNStepTotal)/dble(nb_procs)) + 1
			xEnd   = xNStepTotal
		else
			xStart = rang*ceiling(dble(xNStepTotal)/dble(nb_procs)) + 1
			xEnd   = (rang+1)*ceiling(dble(xNStepTotal)/dble(nb_procs))
		end if

		if(present(sizeLoc))  sizeLoc  = xEnd - xStart + 1 !Used to escape the not-used places of sizeUnif
		if(present(sizeUnif)) sizeUnif = ceiling(dble(xNStepTotal)/dble(nb_procs)) !Used to allow the MPI_GATHER afterwards
		if(present(start))    start    = xStart
		if(present(end))      end      = xEnd

	end subroutine get_sizes_MPI
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine set_allRandFieldUnstruct(randField, xPoints, rang,    &
						      			all_RandField, all_xPoints)
    	implicit none
    	!INPUT
    	double precision, dimension(:, :), intent(in) :: randField
		double precision, dimension(:, :), intent(in) :: xPoints;
    	integer                          , intent(in) :: rang

    	!OUTPUT
    	double precision, dimension(:, :), allocatable, intent(out) :: all_RandField;
    	double precision, dimension(:, :), allocatable, intent(out) :: all_xPoints;

    	!LOCAL VARIABLES
    	integer :: nb_procs, code, Nmc, nDim, xNTotalProc;
    	integer :: xNTotalGlob, i, j, globCount, start;
		double precision, dimension(:, :), allocatable :: all_RandFieldTransp, all_xPointsTransp;
    	integer         , dimension(:)   , allocatable :: all_xNTotal, procDpl;
    	integer         , dimension(:)   , allocatable :: mapping;


    	call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

		Nmc         = size(randField, 2);
		nDim        = size(xPoints, 2);
		xNTotalProc = size(xPoints, 1)


		allocate (all_xNTotal(nb_procs))
		allocate (procDpl(nb_procs))

		!Gathering range and steps information
		call MPI_ALLGATHER(xNTotalProc , 1, MPI_INTEGER, &
		                   all_xNTotal , 1, MPI_INTEGER, &
		                   MPI_COMM_WORLD,  code)

		xNTotalGlob = sum(all_xNTotal)

		if(rang == 0) then
			allocate(all_xPoints(xNTotalGlob, nDim))
			allocate(all_xPointsTransp(nDim, xNTotalGlob))
			allocate(all_RandField(xNTotalGlob, Nmc))
			allocate(all_RandFieldTransp(Nmc, xNTotalGlob))
			all_xPoints         = 0;
			all_xPointsTransp   = 0;
			all_RandField       = 0;
			all_RandFieldTransp = 0;
		end if

		!Gathering points coordinates
		procDpl         = 0
		do i = 2, nb_procs
			procDpl(i) = procDpL(i-1) + nDim*all_xNTotal(i-1) !Calculating deplacements to gather the random fields
		end do
		call MPI_GATHERV(transpose(xPoints), size(xPoints), MPI_DOUBLE_PRECISION,              &
		                 all_xPointsTransp, all_xNTotal*nDim, procDpl, MPI_DOUBLE_PRECISION, &
		                 0            , MPI_COMM_WORLD , code)
		if(rang == 0) all_xPoints = transpose(all_xPointsTransp)

		!Gathering Random Fields
		procDpl         = 0
		do i = 2, nb_procs
			procDpl(i) = procDpL(i-1) + Nmc*all_xNTotal(i-1) !Calculating deplacements to gather the random fields
		end do
		call MPI_GATHERV(transpose(randField), size(randField), MPI_DOUBLE_PRECISION,              &
		                 all_RandFieldTransp, all_xNTotal*Nmc, procDpl, MPI_DOUBLE_PRECISION, &
		                 0            , MPI_COMM_WORLD , code)
		if(rang == 0) all_RandField = transpose(all_RandFieldTransp)

		!if(rang == 0) call DispCarvalhol(all_RandField, "all_RandField")

		!Reordering
		if(rang == 0) call reorderToGlobal(all_RandField, all_xPoints, mapping)

		!!Printing
!		start = 0
!		if(rang /= 0) start = sum(xNTotalProc(:rang))
!		write(*,*) "xPoints and randField. Proc = ", rang
!		do i = 1, product(xNStep)
!			write(*,'(A, I10, A, 3F10.5, A, F10.5)') "Point ", start + i, ">>", xPoints(i,:), "  RF = ", RandField(i,1)
!		end do
!		if(rang == 0) then
!			write(*,*) "all_xPoints and all_RandField"
!			do i = 1, xNTotalGlob
!				write(*,'(I5, A, I5, A, 3F10.5, A, F10.5)') i, " from ", mapping(i), ">>", all_xPoints(i,:), "  RF = ", all_RandField(i,1)
!			end do
!		end if
		!!/Printing


		if(allocated(procDpl))             deallocate(procDpl)
		if(allocated(all_xNTotal))         deallocate(all_xNTotal)
		if(allocated(all_xPointsTransp))   deallocate(all_xPointsTransp)
		if(allocated(all_RandFieldTransp)) deallocate(all_RandFieldTransp)
		if(allocated(mapping))             deallocate(mapping)

    end subroutine set_allRandFieldUnstruct

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine set_allRandFieldStructured(randField, xMin, xMax, xNStep, rang,    &
					            	      all_xMin, all_xMax, all_xNStep,         &
					           			  all_RandField, all_xPoints)
    	implicit none
    	!INPUT
    	double precision, dimension(:, :), intent(in) :: randField
    	double precision, dimension(:)   , intent(in) :: xMax, xMin;
    	integer         , dimension(:)   , intent(in) :: xNStep;
    	integer                          , intent(in) :: rang

    	!OUTPUT
    	double precision, dimension(:, :), allocatable, intent(out) :: all_RandField;
    	double precision, dimension(:, :), allocatable, intent(out) :: all_xPoints;
    	integer         , dimension(:, :), allocatable, intent(out) :: all_xNStep;
    	double precision, dimension(:, :), allocatable, intent(out) :: all_xMin, all_xMax;

    	!LOCAL VARIABLES
    	integer :: nb_procs, code, Nmc, nDim, xNTotalProc;
    	integer :: xNTotalGlob, i, j, start, globCount;
		double precision, dimension(:, :), allocatable :: all_RandFieldTransp;
    	integer         , dimension(:)   , allocatable :: all_xNTotal, procDpl;
    	integer         , dimension(:)   , allocatable :: mapping;
!    	double precision, dimension(:, :), allocatable :: xPoints; !Only for tests

		Nmc  = size(randField, 2);
		nDim = size(xNStep);
    	call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

		allocate (all_xMin   (nDim, nb_procs))
		allocate (all_xMax   (nDim, nb_procs))
		allocate (all_xNStep (nDim, nb_procs))
		allocate (all_xNTotal(nb_procs))
		allocate (procDpl    (nb_procs))

		!Gathering range and steps information
		call MPI_ALLGATHER(xMin     , nDim, MPI_DOUBLE_PRECISION, &
		                   all_xMin , nDim, MPI_DOUBLE_PRECISION, &
		                   MPI_COMM_WORLD,  code)
		call MPI_ALLGATHER(xMax     , nDim, MPI_DOUBLE_PRECISION, &
		                   all_xMax , nDim, MPI_DOUBLE_PRECISION, &
		                   MPI_COMM_WORLD,  code)
		call MPI_ALLGATHER(xNStep     , nDim, MPI_INTEGER, &
		                   all_xNStep , nDim, MPI_INTEGER, &
		                   MPI_COMM_WORLD,  code)

		if(rang == 0) then
			call DispCarvalhol(all_xMin, "all_xMin")
			call DispCarvalhol(all_xMax, "all_xMax")
			call DispCarvalhol(all_xNStep, "all_xNStep")
		end if

		xNTotalProc = product(xNStep)
		all_xNTotal = product(all_xNStep, 1)
		xNTotalGlob = sum(all_xNTotal)

		if(rang == 0) then
			allocate(all_xPoints(xNTotalGlob, nDim))
			allocate(all_RandField(xNTotalGlob, Nmc))
			allocate(all_RandFieldTransp(Nmc, xNTotalGlob))
			all_xPoints   = 0;
			all_RandField = 0;
			all_RandFieldTransp = 0;
		end if

		!Gathering Random Fields
		procDpl = 0
		do i = 2, nb_procs
			procDpl(i) = procDpL(i-1) + Nmc*all_xNTotal(i-1)
		end do
		call MPI_GATHERV(transpose(randField), size(randField), MPI_DOUBLE_PRECISION,              &
		                 all_RandFieldTransp, all_xNTotal*Nmc, procDpl, MPI_DOUBLE_PRECISION, &
		                 0            , MPI_COMM_WORLD , code)
		if(rang == 0) all_RandField = transpose(all_RandFieldTransp)

		!Calculating global vector of coordinates (all_xPoints)
		if(rang == 0) then
			globCount = 1
			do i = 1, nb_procs
				!write(*,*) "Permutation X (Nmc = 1), rang ", i
				do j = 1, all_xNTotal(i)
					call get_Permutation(j, all_xMax(:, i), all_xNStep(:, i), all_xPoints(globCount,:), all_xMin(:, i));
					globCount = globCount + 1
				end do
			end do
		end if

		!Reordering
		if(rang == 0) call reorderToGlobal(all_RandField, all_xPoints, mapping)

		!!Printing
!		start = 0
!		if(rang /= 0) start = sum(xNTotalProc(:rang))
!		write(*,*) "xPoints and randField. Proc = ", rang
!		do i = 1, product(xNStep)
!			write(*,'(A, I10, A, 3F10.5, A, F10.5)') "Point ", start + i, ">>", xPoints(i,:), "  RF = ", RandField(i,1)
!		end do
!		if(rang == 0) then
!			write(*,*) "all_xPoints and all_RandField"
!			do i = 1, xNTotalGlob
!				write(*,'(I5, A, I5, A, 3F10.5, A, F10.5)') i, " from ", mapping(i), ">>", all_xPoints(i,:), "  RF = ", all_RandField(i,1)
!			end do
!		end if
		!!/Printing


!		if(allocated(xPoints))             deallocate(xPoints)
		if(allocated(procDpl))             deallocate(procDpl)
		if(allocated(all_xNTotal))         deallocate(all_xNTotal)
		if(allocated(all_RandFieldTransp)) deallocate(all_RandFieldTransp)
		if(allocated(mapping))             deallocate(mapping)

    end subroutine set_allRandFieldStructured

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!TRASH

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!	function calculateAverage(vector) result(average)
!
!		implicit none
!		!INPUT
!		double precision, dimension(:), intent(in) :: vector;
!		!OUTPUT
!		double precision :: average;
!		!LOCAL VARIABLES
!		integer :: i;
!
!
!		average = 0d0;
!
!		do i = 1, size(vector)
!			average = average + vector(i);
!		end do
!
!		average = average/dble(size(vector))
!
!	end function calculateAverage

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!	function calculateStdDeviation(vector, average) result(stdDeviation)
!		implicit none
!		!INPUT
!		double precision, dimension(:), intent(in) :: vector;
!		double precision, intent(in), optional :: average;
!		!OUTPUT
!		double precision :: stdDeviation;
!		!LOCAL VARIABLES
!		integer :: i;
!		double precision :: avg;
!
!
!		if(present(average)) then
!			avg = average;
!		else
!			avg = calculateAverage(vector);
!		endif
!
!		stdDeviation = (calculateAverage(vector**2d0) - avg**2d0)**(0.5d0);
!
!	end function calculateStdDeviation
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!	subroutine set_StatisticsND_MPIOLD(randField, xNStep, xMax, rang,                  &
!									sumRF, sumRFsquare, ptAvg, ptStdDev, procCorrL, &
!									totalSumRF, totalSumRFsquare,                   &
!									evntAvg, evntStdDev,                            &
!									globalAvg, globalStdDev)
!
!		implicit none
!		!INPUT
!		integer                          , intent(in) :: rang;
!		integer         , dimension(:)   , intent(in) :: xNStep;
!		double precision, dimension(:, :), intent(in) :: randField;
!		double precision, dimension(:)   , intent(in) :: xMax;
!		!OUTPUT
!		double precision, dimension(:), allocatable, intent(out) :: sumRF, sumRFsquare, totalSumRF, totalSumRFsquare;
!		double precision, dimension(:), allocatable, intent(out) :: evntAvg, evntStdDev, procCorrL;
!    	double precision, dimension(:), allocatable, intent(out) :: ptAvg, ptStdDev;
!    	double precision                           , intent(out) :: globalAvg, globalStdDev;
!
!		!LOCAL VARIABLES
!		integer :: i, Nmc, nPermut,  nDim;
!		integer :: sizeLoc;
!		integer :: code, xNTotal;

!		write(*,*) "";
!		write(*,*) "------------START set_StatisticsNDmpi (proc = ", rang, ")-----------------------------------------";
!		write(*,*) "";

!		!write (*,*) ">>>>get_sizes_MPI Proc = ", rang
!		call get_sizes_MPI (xNStep, sizeLoc)
!		nDim        = size(xNStep)
!		xNTotal = product(xNStep)
!		Nmc         = size(randField,2)
!
!		!write (*,*) ">>>>Allocating Proc = ", rang
!
!		!Allocating
!		allocate(sumRF(Nmc))
!		allocate(sumRFsquare(Nmc))
!		allocate(ptAvg(sizeLoc))
!		allocate(ptStdDev(sizeLoc))
!		allocate(procCorrL(nDim));
!		if(rang == 0) then
!			allocate (totalSumRF(Nmc))
!			allocate (totalSumRFsquare(Nmc))
!			totalSumRF = -1.0d0
!			totalSumRFsquare = -1.0d0
!			allocate (evntAvg(Nmc))
!			allocate (evntStdDev(Nmc))
!		end if
!
!
!		!write (*,*) ">>>>sizeLoc = ", sizeLoc, " Proc = ", rang
!
!		sumRF(:)       = sum( randField(1:sizeLoc, :)    , dim = 1)
!		sumRFsquare(:) = sum((randField(1:sizeLoc, :))**2, dim = 1)
!
!		!call DispCarvalhol(sumRF, "sumRF");
!		!call DispCarvalhol(sumRFsquare, "sumRFsquare", "F15.8");
!
!		!Calculating events statistics
!		!write (*,*) ">>>>MPI_REDUCE Proc = ", rang
!
!		call MPI_REDUCE(sumRF, totalSumRF, Nmc, MPI_DOUBLE_PRECISION, &
!					    MPI_SUM, 0, MPI_COMM_WORLD, code)
!		call MPI_REDUCE(sumRFsquare, totalSumRFsquare, Nmc, MPI_DOUBLE_PRECISION, &
!					MPI_SUM, 0, MPI_COMM_WORLD, code)
!
!		!write (*,*) ">>>>after MPI_REDUCE Proc = ", rang
!
!		if(rang == 0) then
!			evntAvg      = totalSumRF/dble(xNTotal);
!			evntStdDev   = sqrt(totalSumRFsquare/dble(xNTotal) &
!		             	   - (totalSumRF/dble(xNTotal))**2)
!			globalAvg    = sum(totalSumRF)/dble(xNTotal*Nmc);
!			globalStdDev = sqrt(sum(totalSumRFsquare)/dble(xNTotal*Nmc) &
!		              	   - (sum(totalSumRF)/dble(xNTotal*Nmc))**2)

!		    call DispCarvalhol(evntAvg, "evntAvg")
!		    call DispCarvalhol(evntStdDev, "evntStdDev")
!		    write(*,*) "globalAvg    = ", globalAvg
!		    write(*,*) "globalStdDev = ", globalStdDev

!		write(*,*) ">>>>>>>>> Showing Statistics of each event";
!		write(*,*) "totalSumRF = ", totalSumRF;
!		write(*,*) "totalSumRFsquare = ", totalSumRFsquare;
!		write(*,*) " Average = ", totalSumRF/dble(xNTotal);
!		write(*,*) "Variance = ", totalSumRFsquare/dble(xNTotal) - (totalSumRF/dble(xNTotal))**2;
!		write(*,*) " Std Dev = ", sqrt(totalSumRFsquare/dble(xNTotal) - (totalSumRF/dble(xNTotal))**2)
!		write(*,*) ">>>>>>>>> Showing Global Statistics (all events)";
!		write(*,*) " Average = ", sum(totalSumRF)/dble(xNTotal*Nmc);
!		write(*,*) "Variance = ", sum(totalSumRFsquare)/dble(xNTotal*Nmc) - (sum(totalSumRF)/dble(xNTotal*Nmc))**2;
!		write(*,*) " Std Dev = ", sqrt(sum(totalSumRFsquare)/dble(xNTotal*Nmc) - (sum(totalSumRF)/dble(xNTotal*Nmc))**2)
!		end if
!
!		do i = 1, sizeLoc
!			ptAvg(i)    = calculateAverage(randField(i,:));
!			ptStdDev(i) = calculateStdDeviation(randField(i,:), ptAvg(i));
!		end do

		!write (*,*) "Rang = ", rang
		!call DispCarvalhol(ptAvg, "ptAvg")
		!if(rang == 0) call DispCarvalhol(ptStdDev, "ptStdDev")

		!call calculateAverageCorrL_MPI(randField, xMax, xNStep, ptAvg, ptStdDev, procCorrL)

		!write(*,*) "      sumRF = ", sumRF
		!write(*,*) "sumRFsquare = ", sumRFsquare

		!write(*,*) "";
		!write(*,*) "------------END set_StatisticsNDmpi (proc = ", rang, ")-----------------------------------------";
		!write(*,*) "";

!	end subroutine set_StatisticsND_MPIOLD

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!	subroutine set_StatisticsND(totalRandField, xNStep, xMax, &
!					 			totalAverage, totalStdDeviation, totalAverageCorrL)
!
!		implicit none
!		!INPUT
!		double precision, dimension(:, :),           intent(in) :: totalRandField;
!		double precision, dimension(:),              intent(in) :: xMax;
!		integer,          dimension(:),              intent(in) :: xNStep;
!		!OUTPUT
!		double precision, dimension(:), allocatable, intent(out) :: totalAverage, totalStdDeviation, &
!																	totalAverageCorrL;
!		!LOCAL VARIABLES
!		integer :: i, Nmc, nPermut,  nDim;
!
!		write(*,*) "";
!		write(*,*) "------------START set_StatisticsND-----------------------------------------";
!		write(*,*) "";
!
!		nPermut = size(totalRandField, 1);
!		Nmc     = size(totalRandField, 2);
!		nDim    = size(xNStep   , 1)
!
!		allocate(totalAverage(nPermut));
!		allocate(totalStdDeviation(nPermut));
!		allocate(totalAverageCorrL(nDim   ));
!		totalAverage = -1
!		totalStdDeviation = -1
!		totalAverageCorrL = -1
!		write(*,*) ">>>>>>>>> Calculating Average and Stantard Deviation";
!		do i = 1, nPermut
!			totalAverage(i)      = calculateAverage(totalRandField(i,:));
!			totalStdDeviation(i) = calculateStdDeviation(totalRandField(i,:), totalAverage(i));
!		end do
!
!		if(Nmc > 1) then
!			write(*,*) ">>>>>>>>> Calculating Correlation Length";
!			call calculateAverageCorrL(totalRandField, xMax, xNStep, totalAverageCorrL)
!		else
!			write(*,*) ">>>>>>>>> WARNING! Correlation Length coudn't be computed (only one event)";
!		end if
!
!		write(*,*) "";
!		write(*,*) "------------END set_StatisticsND-----------------------------------------";
!		write(*,*) "";
!
!	end subroutine set_StatisticsND

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    subroutine set_normalisation(xMin, xMax, corrL, xNorm)
!
!    	implicit none
!    	!INPUT
!    	double precision, dimension(:), intent(in) :: xMin, xMax, corrL;
!    	!LOCAL VARIABLES
!    	double precision, dimension(:, :), allocatable, intent(out) :: all_xMin, all_xMax;
!		integer :: nb_procs, code, nDim;
!
!		nDim = size(corrL);-
!    	call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
!
!		allocate (all_xMin(nDim, nb_procs))
!		allocate (all_xMax(nDim, nb_procs))
!
!		!Gathering range and steps information
!		call MPI_ALLGATHER(xMin     , nDim, MPI_DOUBLE_PRECISION, &
!		                   all_xMin , nDim, MPI_DOUBLE_PRECISION, &
!		                   MPI_COMM_WORLD,  code)
!		call MPI_ALLGATHER(xMax     , nDim, MPI_DOUBLE_PRECISION, &
!		                   all_xMax , nDim, MPI_DOUBLE_PRECISION, &
!		                   MPI_COMM_WORLD,  code)
!
!
!
!		deallocate (all_xMin)
!		deallocate (all_xMax)
!
!   	end subroutine set_normalisation

!TRASH
		!Assembling all the fields in one matrix";

!		!Type all_lineRF (global)
!		call MPI_TYPE_VECTOR(Nmc, 1, xNTotal,                   &
!		                     MPI_DOUBLE_PRECISION, type_Temp, code) !type_Temp and code are outputs
!		call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dblBitSize, code)
!		call MPI_TYPE_CREATE_RESIZED(type_Temp, 0, dblBitSize, all_lineRF, code) !Changing the starting point to the next "slice"
!		call MPI_TYPE_COMMIT(all_lineRF, code)
!		write(*,*) "Type line global rang ", rang, " id = ", all_lineRF


!		write(*,*) "Type line local rang ", rang, " id = ", lineRF(rang+1)

!		call MPI_TYPE_FREE(all_lineRF, code)
!		call MPI_TYPE_FREE(lineRF, code)
!
!		call MPI_TYPE_VECTOR(size(randField,2), size(randField,1), size(randField,1)*nb_procs, &
!		                     MPI_DOUBLE_PRECISION, type_Temp, code) !type_Temp and code are outputs
!		call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION, dblBitSize, code)
!		newExtent = size(randField,1)*dblBitSize
!		call MPI_TYPE_CREATE_RESIZED(type_Temp, 0, newExtent, type_RF, code) !Changing the starting point to the next "slice"
!		call MPI_TYPE_COMMIT(type_RF, code)
!		call MPI_GATHER(randField    , size(randField), MPI_DOUBLE_PRECISION, &
!		                all_RandField, 1              , type_RF,               &
!		                0            , MPI_COMM_WORLD , code)
!		call MPI_TYPE_FREE(type_RF, code)
!		write (*,*) "Proc = ", rang
!		call DispCarvalhol(randField, "randField")

		!!!!!TEST
!		do i = 1, size(randField, 2)
!			do j = 1, size(randField, 1)
!				randField(j,i) = 100*rang + (i-1)*size(randField, 1) + j; !Testing where the fields go in the global matrix
!			end do
!		end do
		!!!!!TEST

!		if(rang == 0) call DispCarvalhol(procDpl, "procDpl")
!		if(rang == 0) call DispCarvalhol(transpose(randField), "transpose(randField)")
!		if(rang == 0) call DispCarvalhol(all_RandFieldTransp, "all_RandFieldTransp")
!		if(rang == 0) write(*,*) "xNTotalProc*nDim = ", xNTotalProc*nDim

end module obsolete_RF
