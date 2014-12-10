module statistics_RF

	use displayCarvalhol
	use math_RF
	use mpi

    interface set_Statistics_MPI
		module procedure set_StatisticsStructured_MPI,   &
		                 set_StatisticsUnstruct_MPI
	end interface set_Statistics_MPI

	interface set_CompStatistics
		module procedure set_CompStatisticsStructured,   &
		                 set_CompStatisticsUnstruct
	end interface set_CompStatistics

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_StatisticsUnstruct_MPI(randField, xPoints, rang, &
							 	           evntAvg, evntStdDev, procCorrL,             &
							 	           globalAvg, globalStdDev, globalCorrL)
		implicit none

		!INPUT
		double precision, dimension(:, :), intent(in) :: randField, xPoints;
		integer                          , intent(in) :: rang;

		!OUTPUT
		double precision, dimension(:), allocatable, intent(out) :: evntAvg, evntStdDev, procCorrL, globalCorrL;
    	double precision                           , intent(out) :: globalAvg, globalStdDev;

		!LOCAL VARIABLES
		double precision, dimension(:)   , allocatable :: avg;
		double precision, dimension(:),    allocatable :: sumRF, sumRFsquare, &
		                                                  totalSumRF, totalSumRFsquare;
		integer :: Nmc, nPoints, nDim, xNTotal;
		integer :: i, j, code, nb_procs;

		Nmc          = size(randField, 2)
		nPoints      = size(randField, 1)
		nDim         = size(xPoints, 1)
		xNTotal  	 = size(xPoints, 2)
		globalAvg    = -1.0d0
		globalStdDev = -1.0d0

		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

		!Allocating
		allocate(sumRF(Nmc))
		allocate(sumRFsquare(Nmc))
		allocate(procCorrL(nDim));
		if(rang == 0) then
			allocate (totalSumRF(Nmc))
			allocate (totalSumRFsquare(Nmc))
			totalSumRF = -1.0d0
			totalSumRFsquare = -1.0d0
			allocate (evntAvg(Nmc))
			allocate (evntStdDev(Nmc))
			allocate(globalCorrL(nDim));
		end if

		!Calculating Correlation Length (should be reformulated to take advantage of matrix symmetry)
!		call set_CorrelationLengthUnstruct(randField, xPoints, procCorrL)

		!Setting variables to calculate Average and Standard Deviation (by event and global)
		sumRF(:)       = sum( randField    , dim = 1)
		sumRFsquare(:) = sum((randField)**2, dim = 1)

		call MPI_REDUCE(sumRF, totalSumRF, Nmc, MPI_DOUBLE_PRECISION, &
					    MPI_SUM, 0, MPI_COMM_WORLD, code)
		call MPI_REDUCE(sumRFsquare, totalSumRFsquare, Nmc, MPI_DOUBLE_PRECISION, &
					    MPI_SUM, 0, MPI_COMM_WORLD, code)
		call MPI_REDUCE(procCorrL, globalCorrL, nDim, MPI_DOUBLE_PRECISION, &
					    MPI_SUM, 0, MPI_COMM_WORLD, code)

		if(rang == 0) then
			!by Event
			evntAvg      = totalSumRF/dble(xNTotal*nb_procs);
			evntStdDev   = sqrt(totalSumRFsquare/dble(xNTotal*nb_procs) &
		             	   - (evntAvg)**2)
		    !Global
			globalAvg    = sum(totalSumRF)/dble(xNTotal*nb_procs*Nmc);
			globalStdDev = sqrt(sum(totalSumRFsquare)/dble(xNTotal*Nmc*nb_procs) &
		              	   - (globalAvg)**2)
		    !globalCorrL  = globalCorrL / nb_procs

		    !call DispCarvalhol(evntAvg   , "evntAvg")
		    !call DispCarvalhol(evntStdDev, "evntStdDev")
		    !write(*,*) "globalAvg    = ", globalAvg
		    !write(*,*) "globalStdDev = ", globalStdDev

		end if

		if(allocated(sumRF))       deallocate(sumRF)
		if(allocated(sumRFsquare)) deallocate(sumRFsquare)
		if(allocated(totalSumRF))       deallocate (totalSumRF)
		if(allocated(totalSumRFsquare)) deallocate (totalSumRFsquare)

	end subroutine set_StatisticsUnstruct_MPI

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_CorrelationLengthUnstruct(randField, xPoints, procCorrL)
		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in)  :: randField, xPoints;
		!OUTPUT
		double precision, dimension(:), intent(out) :: procCorrL;
		!LOCAL VARIABLES
		double precision, dimension(:, :)   , allocatable :: covMatrix, sumOthers;
		double precision, dimension(:, :)   , allocatable :: deltaMatrix;
		double precision, dimension(:, :, :), allocatable :: distMatrix;
		integer          :: Nmc, nPoints, nDim, i, j, nFactors;
		double precision :: tolerance = 1.0d-6;

		Nmc     = size(randField, 2)
		nPoints = size(xPoints, 2)
		nDim    = size(xPoints, 1)

		allocate(covMatrix(nPoints, nPoints))
		allocate(sumOthers(nPoints, nPoints))
		allocate(deltaMatrix(nDim, nPoints))
		allocate(distMatrix(nPoints, nPoints, nDim))

		call set_CovMatrix(randField, covMatrix)
		call set_DeltaMatrix(xPoints, deltaMatrix)
		call set_DistMatrix(xPoints, distMatrix)

!		call DispCarvalhol(xPoints,"xPoints")
!		call DispCarvalhol(covMatrix,"covMatrix", nColumns = 15)
!		call DispCarvalhol(deltaMatrix,"deltaMatrix", nColumns = 15)
!		call DispCarvalhol(distMatrix,"distMatrix", nColumns = 15)

		do i = 1, nDim
			!write(*,*) "nDim = ", i, "-------------------------"
			sumOthers = sum(distMatrix, 3) - distMatrix(:,:,i)
			!call DispCarvalhol(sumOthers,"sumOthers BEFORE")
			nFactors = count(sumOthers < tolerance);
			!write(*,*) "nFactors = ", nFactors
			where(sumOthers < tolerance)
				sumOthers(:,:) = covMatrix (:,:)
			elsewhere
				sumOthers(:,:) = 0
			end where

			!call DispCarvalhol(sumOthers,"sumOthers BEFORE")
			do j = 1, nPoints
				sumOthers(:,j) = deltaMatrix(i, j) * sumOthers(:,j)
			end do
			!call DispCarvalhol(sumOthers,"sumOthers AFTER")

			procCorrL(i) = sum(sumOthers)/sqrt(dble(nFactors))
		end do
!
!		call DispCarvalhol(procCorrL,"procCorrL")

		if(allocated(distMatrix))  deallocate(distMatrix)
		if(allocated(covMatrix))   deallocate(covMatrix)
		if(allocated(deltaMatrix)) deallocate(deltaMatrix)
		if(allocated(sumOthers))   deallocate(sumOthers)

	end subroutine set_CorrelationLengthUnstruct

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_StatisticsStructured_MPI(randField, xMin, xMax, xNStep, rang, &
							 	    		evntAvg, evntStdDev, procCorrL,      &
							 	    		globalAvg, globalStdDev, globalCorrL)
		implicit none

		!INPUT
		double precision, dimension(:, :), intent(in) :: randField;
		double precision, dimension(:)   , intent(in) :: xMin,xMax;
		integer         , dimension(:)   , intent(in) :: xNStep;
		integer                          , intent(in) :: rang;

		!OUTPUT
		double precision, dimension(:), allocatable, intent(out) :: evntAvg, evntStdDev, procCorrL, globalCorrL;
    	double precision                           , intent(out) :: globalAvg, globalStdDev;

		!LOCAL VARIABLES
		double precision, dimension(:)   , allocatable :: avg;
		double precision, dimension(:),    allocatable :: sumRF, sumRFsquare, &
		                                                  totalSumRF, totalSumRFsquare;
		integer :: Nmc, nPoints, nDim, xNTotal;
		integer :: i, j, code, nb_procs;

		Nmc          = size(randField, 2)
		nPoints      = size(randField, 1)
		nDim         = size(xNStep)
		xNTotal  = product(xNStep)
		globalAvg    = -1.0d0
		globalStdDev = -1.0d0

		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)

		!Allocating
		allocate(sumRF(Nmc))
		allocate(sumRFsquare(Nmc))
		allocate(procCorrL(nDim));
		if(rang == 0) then
			allocate (totalSumRF(Nmc))
			allocate (totalSumRFsquare(Nmc))
			totalSumRF = -1.0d0
			totalSumRFsquare = -1.0d0
			allocate (evntAvg(Nmc))
			allocate (evntStdDev(Nmc))
			allocate(globalCorrL(nDim));
		end if

		!Calculating Correlation Length (should be reformulated to take advantage of matrix symmetry)
		call set_CorrelationLengthStructured(randField, xMin, xMax, xNStep, procCorrL)

		!Setting variables to calculate Average and Standard Deviation (by event and global)
		sumRF(:)       = sum( randField    , dim = 1)
		sumRFsquare(:) = sum((randField)**2, dim = 1)

		call MPI_REDUCE(sumRF, totalSumRF, Nmc, MPI_DOUBLE_PRECISION, &
					    MPI_SUM, 0, MPI_COMM_WORLD, code)
		call MPI_REDUCE(sumRFsquare, totalSumRFsquare, Nmc, MPI_DOUBLE_PRECISION, &
					    MPI_SUM, 0, MPI_COMM_WORLD, code)
		call MPI_REDUCE(procCorrL, globalCorrL, nDim, MPI_DOUBLE_PRECISION, &
					    MPI_SUM, 0, MPI_COMM_WORLD, code)

		if(rang == 0) then
			!by Event
			evntAvg      = totalSumRF/dble(xNTotal*nb_procs);
			evntStdDev   = sqrt(totalSumRFsquare/dble(xNTotal*nb_procs) &
		             	   - (evntAvg)**2)
		    !Global
			globalAvg    = sum(totalSumRF)/dble(xNTotal*nb_procs*Nmc);
			globalStdDev = sqrt(sum(totalSumRFsquare)/dble(xNTotal*Nmc*nb_procs) &
		              	   - (globalAvg)**2)
		    globalCorrL  = globalCorrL / nb_procs

		    !call DispCarvalhol(evntAvg   , "evntAvg")
		    !call DispCarvalhol(evntStdDev, "evntStdDev")
		    !write(*,*) "globalAvg    = ", globalAvg
		    !write(*,*) "globalStdDev = ", globalStdDev

		end if

		if(allocated(sumRF))       deallocate(sumRF)
		if(allocated(sumRFsquare)) deallocate(sumRFsquare)
		if(allocated(totalSumRF))       deallocate (totalSumRF)
		if(allocated(totalSumRFsquare)) deallocate (totalSumRFsquare)

	end subroutine set_StatisticsStructured_MPI

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_CorrelationLengthStructured(randField, xMin, xMax, xNStep, procCorrL)
		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in) :: randField;
		double precision, dimension(:)   , intent(in) :: xMin,xMax;
		integer         , dimension(:)   , intent(in) :: xNStep;
		!OUTPUT
		double precision, dimension(:), intent(out) :: procCorrL;
		!LOCAL VARIABLES
		double precision, dimension(:, :), allocatable :: covMatrix;
		integer :: Nmc, nPoints, nDim, xNTotal;
		integer :: beg, step, end, plane, radius, rStart, rEnd, elemStep, patternStep;
		integer :: i, j, code, nb_procs;

		Nmc     = size(randField, 2)
		nPoints = size(randField, 1)
		nDim    = size(xNStep)

		allocate(covMatrix(nPoints, nPoints))

		call set_CovMatrix(randField, covMatrix)

		do i = 1, nDim
			!if (rang == 0) write(*,*) "Dim = ",   i,"----------------------"
			!nPlanes = xNTotal/xNStep(i)
			call get_SequenceParam(i, 1, xNStep, beg, elemStep, end) !establishing some reference parameters
			patternStep = end - beg + elemStep
			do j = 1, xNTotal
				plane = cyclicMod(j, elemStep) +  int_roundDown(j, patternStep)*elemStep !Obs: j and patternRange are integers so we are rounding down
				call get_SequenceParam(i, plane, xNStep, beg, step, end)

				procCorrL(i) = procCorrL(i) +                   &
				               sum(covMatrix(j, beg:end:step))
				!Averaging the double terms
				radius = min((j - beg)/step, (end - j)/step)*step
				rStart = j - radius
				rEnd   = j + radius
				procCorrL(i) =   procCorrL(i)                                      &
							   - sum(covMatrix(j, rStart:rEnd:step))/2 &
							   + covMatrix(j,j)/2

				!if (rang == 0) write(*,*) ""
				!if (rang == 0) write(*,*) "Plane = ", plane, "patternStep = ", patternStep
				!if (rang == 0) write(*,*) j, ">> beg = ", beg, "step  = ", step, "end = ", end
				!if (rang == 0) write(*,*) "radius ", radius, "Doubles from = ", rStart, "    to ", rEnd

			end do
			procCorrL(i) = ((xMax(i)-xMin(i))/xNStep(i))*procCorrL(i)/xNTotal
		end do

	end subroutine set_CorrelationLengthStructured

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_CovMatrix(randField, covMatrix)

		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in)  :: randField;
		!OUTPUT
		double precision, dimension(:, :), intent(out) :: covMatrix;
		!LOCAL VARIABLES
		double precision, dimension(:)   , allocatable :: avg;
		double precision, dimension(:, :), allocatable :: centeredRF;
		double precision :: diagStdDev
		integer          :: Nmc, nPoints, i;

		Nmc     = size(randField, 2)
		nPoints = size(randField, 1)

		allocate(avg        (nPoints))
		allocate(centeredRF (nPoints, Nmc))

		avg = sum(randField, 2) / Nmc

		do i = 1, nPoints
			centeredRF(i,:) = randField(i,:) - avg(i)
		end do

		covMatrix = matmul(centeredRF, transpose(centeredRF)) / Nmc

		do i = 1, nPoints
			!Normalising the Covariance Matrix
			diagStdDev = sqrt(covMatrix(i, i))
			covMatrix(i, :) = covMatrix(i, :)/diagStdDev
			covMatrix(:, i) = covMatrix(:, i)/diagStdDev
		end do

		deallocate(avg)
		deallocate(centeredRF)

	end subroutine set_CovMatrix

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_CompStatisticsUnstruct(all_RandField, all_xPoints,                   &
			                        	  compEvntAvg, compEvntStdDev,                  &
			                        	  compGlobAvg, compGlobStdDev, compGlobCorrL)
		implicit none
		!INPUT
		double precision, dimension(:, :), intent(in) :: all_RandField, all_xPoints;
		!OUTPUT
		double precision, dimension(:), allocatable, intent(out) :: compEvntAvg, compEvntStdDev, compGlobCorrL;
		double precision, intent(out) :: compGlobAvg, compGlobStdDev;
		!LOCAL VARIABLES
		integer :: i, Nmc, nPoints,  nDim;

		nPoints = size(all_RandField, 1);
		Nmc     = size(all_RandField, 2);
		nDim    = size(all_xPoints  , 2);

		allocate(compEvntAvg(Nmc));
		allocate(compEvntStdDev(Nmc));
		allocate(compGlobCorrL(nDim));

		compEvntAvg    = sum(all_RandField, 1)/nPoints;
		compEvntStdDev = sqrt(sum(all_RandField**2, 1)/nPoints - compEvntAvg**2);
		compGlobAvg    = sum(all_RandField)/(nPoints*Nmc);
		compGlobStdDev = sqrt(sum(all_RandField**2)/(nPoints*Nmc) - compGlobAvg**2);

		call set_CorrelationLengthUnstruct(all_RandField, all_xPoints, compGlobCorrL)

	end subroutine set_CompStatisticsUnstruct

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine set_CompStatisticsStructured(all_RandField, all_xPoints,                 &
									   		all_xMin, all_xMax, all_xNStep,               &
		                    			    compEvntAvg, compEvntStdDev,                  &
		                   				    compGlobAvg, compGlobStdDev, compGlobCorrL)
		implicit none
		!INPUT
		integer         , dimension(:, :), intent(in) :: all_xNStep;
		double precision, dimension(:, :), intent(in) :: all_RandField, all_xPoints;
		double precision, dimension(:, :), intent(in) :: all_xMin, all_xMax;
		!OUTPUT
		double precision, dimension(:), allocatable, intent(out) :: compEvntAvg, compEvntStdDev, compGlobCorrL;
		double precision, intent(out) :: compGlobAvg, compGlobStdDev;
		!LOCAL VARIABLES
		integer :: i, Nmc, nPoints,  nDim, nbProcs;
		integer, dimension(:), allocatable :: totalXNStep;
		double precision, dimension(:), allocatable :: totalXRange, totalXMax, totalXMin;

		nPoints = size(all_RandField, 1);
		Nmc     = size(all_RandField, 2);
		nDim    = size(all_xPoints  , 2);
		nbProcs = size(all_xMax, 2)

		allocate(compEvntAvg(Nmc));
		allocate(compEvntStdDev(Nmc));
		allocate(compGlobCorrL(nDim));
		allocate(totalXRange(nDim));
		allocate(totalXNStep(nDim));
		allocate(totalXMax(nDim))
		allocate(totalXMin(nDim))

		totalXNStep    = sum(all_xNStep,2);
		compEvntAvg    = sum(all_RandField, 1)/nPoints;
		compEvntStdDev = sqrt(sum(all_RandField**2, 1)/nPoints - compEvntAvg**2);
		compGlobAvg    = sum(all_RandField)/(nPoints*Nmc);
		compGlobStdDev = sqrt(sum(all_RandField**2)/(nPoints*Nmc) - compGlobAvg**2);
		totalXMax      = maxval(all_xMax, 2)
		totalXMin      = minval(all_xMin, 2)
		totalXRange    = totalXMax - totalXMin

		call calculateAverageCorrL(all_RandField, totalXRange, totalXNStep, compGlobCorrL) !To change by set_CorrelationLengthStructured
		!call set_CorrelationLengthStructured(all_RandField, totalXMax, totalXMin, totalXNStep, compGlobCorrL)

		if(allocated(totalXRange)) deallocate(totalXRange);
		if(allocated(totalXNStep)) deallocate(totalXNStep);
		if(allocated(totalXMax))   deallocate(totalXMax);
		if(allocated(totalXMin))   deallocate(totalXMin);

	end subroutine set_CompStatisticsStructured

end module statistics_RF
