module randomFieldND

    use displayCarvalhol
    use spectra_RF
    use math_RF
    use obsolete_RF
	use mpi

    interface createRandomField
		module procedure createRandomFieldUnstruct,   &
		                 createRandomFieldStructured
	end interface createRandomField

	double precision :: periodMult = 1.1 !"range" multiplier
	double precision :: kAdjust = 5 !"kNStep minimum" multiplier

contains

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine createRandomFieldUnstruct(xPoints, corrMod, margiFirst, corrL, fieldAvg, fieldVar, Nmc, randField);

        implicit none

        !INPUT
        double precision, dimension(:, :), intent(in)    :: xPoints;
        double precision, dimension(:)   , intent(in)    :: corrL;
        character (len=*)                , intent(in)    :: corrMod, margiFirst;
        integer                          , intent(in)    :: Nmc;
        double precision                 , intent(in)    :: fieldAvg, fieldVar;

        !OUTPUT
        double precision, dimension(:, :), allocatable, intent(out) :: randField;

		call createStandardGaussianFieldUnstruct (xPoints, corrL, corrMod, Nmc, randField)
		call multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

    end subroutine createRandomFieldUnstruct

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

        implicit none

        !INPUT
        character (len=*), intent(in) :: margiFirst;
        double precision , intent(in) :: fieldAvg, fieldVar;

        !OUTPUT (IN)
        double precision, dimension(1:, 1:), intent(inout) :: randField;

        !LOCAL VARIABLES
        double precision :: normalVar, normalAvg
        integer          :: error, code, i

		select case (margiFirst)
			case("gaussian")
				normalVar = fieldVar
				normalAvg = fieldAvg
			case("lognormal")
				if(fieldAvg <= 0) then
					write(*,*) ""
					write(*,*) "ERROR - when using lognormal fieldAvg should be a positive number greater than 0.001"
					call MPI_ABORT(MPI_COMM_WORLD, error, code)
				end if
				normalVar = log(1 + fieldVar/(fieldAvg**2))
				normalAvg = log(fieldAvg) - normalVar/2
		end select

		randField(:,:) = randField(:,:) * sqrt(normalVar) &
						 + normalAvg;

		if (margiFirst == "lognormal") then
			randField(:,:) = exp(randField(:,:))
		end if

	end subroutine

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine createStandardGaussianFieldUnstruct (xPoints, corrL, corrMod, Nmc, randField, &
	                                                chosenSeed, MinBound, MaxBound, communicator, calculate)

		!INPUT
        double precision, dimension(1:, 1:), intent(in) :: xPoints;
        double precision, dimension(1:)    , intent(in) :: corrL;
        integer                            , intent(in) :: Nmc;
        character (len=*)                  , intent(in) :: corrMod;
        integer, dimension(1:), optional   , intent(in) :: chosenSeed
        real   , dimension(1:), optional   , intent(in) :: MinBound, MaxBound
        integer               , optional   , intent(in) :: communicator
        logical, dimension(1:), optional   , intent(in) :: calculate

        !OUTPUT
        double precision, dimension(:, :), intent(out) :: randField;

        !LOCAL VARIABLES
        integer         , dimension(:)  , allocatable :: kNStep;
        double precision, dimension(:)  , allocatable :: kMax;
        double precision, dimension(:,:), allocatable :: kSign, phiN, xPointsNorm;
        double precision, dimension(:)  , allocatable :: kVec, kVecUnsigned; !Allocated in function
        double precision, dimension(:)  , allocatable :: deltaK, kDelta;
        double precision, dimension(:)  , allocatable :: xMaxGlob, xMinGlob;
        logical         , dimension(:)  , allocatable :: effectCalc;
        integer          :: effectComm
        integer          :: i, j, k, m, nDim;
        integer          :: xNTotal, kNTotal;
        integer          :: nb_procs, rang, code, error;
        double precision :: Sk, deltaKprod;
        double precision :: pi = 3.1415926535898, zero = 0d0;

		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
		call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)

		nDim    = size(xPoints, 2);
		xNTotal = size(xPoints, 1);

		allocate(effectCalc(Nmc))
		effectCalc(:) = .true.
		if(present(calculate)) effectCalc = calculate

		if(present(communicator)) then
			effectComm = communicator
		else
			effectComm = MPI_COMM_WORLD
		end if

		!Normalization
		allocate (xPointsNorm (size(xPoints,1), size(xPoints,2)))
		xPointsNorm(:,:) = xPoints(:,:)
		do i = 1, nDim
			if(corrL(i) /= 1) xPointsNorm(:,i) = xPointsNorm (:,i)/corrL(i)
		end do

		!if(rang == 0) write(*,*) ">>>>>>>>> Variables initialization: kMax, kDelta, kNStep, xMinGlob, xMaxGlob";
		!Allocating
		allocate(kMax   (nDim));
		allocate(kNStep (nDim));
		allocate(kDelta (nDim));


		if((present(MinBound) .and. (.not.present(MaxBound))) .or. &
		   (present(MaxBound) .and. (.not.present(MinBound)))) then
			write(*,*) "ERROR - In 'createStandardGaussianFieldUnstruct': both or none of the bounds should be defined"
			if(.not.present(MinBound)) write(*,*) "MinBound not present"
			if(.not.present(MaxBound)) write(*,*) "MaxBound not present"
			call MPI_ABORT(effectComm  , error, code)
		endif

		!Communicating normalized extremes
        allocate(xMinGlob(nDim))
        allocate(xMaxGlob(nDim))
		if(.not.present(MinBound)) then
			write(*,*) "Min bound not present, the extremes will be calculated automatically"
			call set_Extremes(xPointsNorm, xMinGlob, xMaxGlob, effectComm)
		else
			!write(*,*) "lbound(MinBound) = ", lbound(MinBound)
			!write(*,*) "ubound(MinBound) = ", ubound(MinBound)
			!write(*,*) "lbound(MaxBound) = ", lbound(MaxBound)
			!write(*,*) "ubound(MaxBound) = ", ubound(MaxBound)
			!write(*,*) "lbound(corrL)    = ", lbound(corrL)
			!write(*,*) "ubound(corrL)    = ", ubound(corrL)

			do i = 1, nDim
				!write(*,*) "i = ", i
				xMinGlob(i) = MinBound(i)/corrL(i)
				xMaxGlob(i) = MaxBound(i)/corrL(i)
			end do
		end if

		call set_kMaxND(corrMod, kMax) !Defining kMax according to corrMod
		kDelta  = 2*pi/(periodMult*(xMaxGlob - xMinGlob)) !Delta min in between two wave numbers to avoid periodicity
		kNStep  = kAdjust*(ceiling(kMax/kDelta) + 1);
		kNTotal = product(kNStep);

		!if(rang == 0) write(*,*) "Nmc     = ", Nmc
		!if(rang == 0) write(*,*) "kNTotal = ", kNTotal
		!if(rang == 0) write(*,*) "kDelta  = ", kDelta
		!if(rang == 0) write(*,*) "kNStep  = ", kNStep

		if(kNTotal < 1) then
			write(*,*) "ERROR - In 'createStandardGaussianFieldUnstruct': kNTotal should be a positive integer (possibly a truncation problem)"
			call MPI_ABORT(effectComm, error, code)
		endif

		!Random Field
		allocate(deltaK      (nDim));
		allocate(kVec        (nDim));
		allocate(kVecUnsigned(nDim));
		allocate(kSign       (2**(nDim-1), nDim));
		allocate(phiN        (size(kSign,1), kNTotal));

		!write(*,*) "Flag5 inside"


		if (size(randField, 1) /= xNTotal .or. size(randField, 2) /= Nmc) then
			write(*,*) "ERROR - In 'createStandardGaussianFieldUnstruct': randfield dimensions are imcompatible with the coordinates (xPoints)"
			write(*,*) "shape(randfield(:,:)) = ", size(randField, 1), size(randField, 2)
			call MPI_ABORT(effectComm, error, code)
		end if

		randField(:,:) = 0;
		deltaK(:)      = 0;
	    deltaK(:)      = (kMax)/(kNStep-1); !Defines deltaK
		call set_kSign(kSign) !Set the sign permutations for kVec

		!write(*,*) "Flag6 inside"
		!Initializing the seed
		if(present(chosenSeed)) then
			call init_random_seed(chosenSeed)
		else
			call init_random_seed()
		end if

		!write(*,*) "Flag7 inside"
		!Generating random field samples
		do k = 1, Nmc
			if(effectCalc(k)) then
				call random_number(phiN(:,:))
				do j = 1, kNTotal
					call get_Permutation(j, kMax, kNStep, kVecUnsigned);
				    do m = 1, size(kSign,1)
				    	kVec           = kVecUnsigned * kSign(m, :)
						Sk             = get_SpectrumND(kVec, corrMod);
						randField(:,k) = sqrt(Sk) * cos(matmul(xPoints(:,:), kVec(:)) &
										 + 2*pi*phiN(m, j)) &
										 + randField(:,k)
	!					if(rang == 0 .and. (k*j*m<50)) write(*,*) "Sk  = ", Sk
	!					if(rang == 0 .and. (k*j*m<50)) write(*,*) "cos(matmul(xPoints(:,:), kVec(:))  = ", matmul(xPoints(:,:), kVec(:))
	!					if(rang == 0 .and. (k*j*m<50)) write(*,*) "2*pi*phiN(m, j)  = ", 2*pi*phiN(m, j)
	!					if(rang == 0 .and. (k*j*m<50)) write(*,*) "Contrib  = ", sqrt(Sk) * cos(matmul(xPoints(:,:), kVec(:)) &
	!									 										+ 2*pi*phiN(m, j))
					end do
				end do
			else
				randField(:,k) = 0.0
			end if
		end do

		if(rang == 0) write(*,*) "Spectra (Sk) cut in: ", Sk

		randField(:,:) = 2*sqrt(product(deltaK)/((2*pi)**(nDim)))&
						 * randField(:,:) !Obs: sqrt(product(corrL)) is not needed because of normalization

		if(allocated(deltaK))       deallocate(deltaK);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
		if(allocated(kVec))         deallocate(kVec);
		if(allocated(kMax))         deallocate(kMax);
		if(allocated(kNStep))       deallocate(kNStep);
		if(allocated(kDelta))       deallocate(kDelta);
		if(allocated(kSign))        deallocate(kSign);
		if(allocated(phiN))         deallocate(phiN);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
		if(allocated(xMinGlob))     deallocate(xMinGlob);
        if(allocated(xMaxGlob))     deallocate(xMaxGlob);
        if(allocated(xPointsNorm))  deallocate (xPointsNorm);

	end subroutine createStandardGaussianFieldUnstruct

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine createStandardGaussianFieldUnstructVictor (xPoints, corrL, corrMod, Nmc, randField, &
	                                                      chosenSeed, MinBound, MaxBound, communicator, calculate)

		!INPUT
        double precision, dimension(:, :), intent(in) :: xPoints; !inout because of normalization
        double precision, dimension(:)   , intent(in) :: corrL;
        integer                          , intent(in) :: Nmc;
        character (len=*)                , intent(in) :: corrMod;
        integer, dimension(1:), optional , intent(in) :: chosenSeed
        real   , dimension(1:), optional , intent(in) :: MinBound, MaxBound
        integer               , optional , intent(in) :: communicator
        logical, dimension(1:), optional , intent(in) :: calculate
        !OUTPUT
        double precision, dimension(:, :), intent(out) :: randField;

        !LOCAL VARIABLES
        integer         , dimension(:)  , allocatable :: kNStep;
        double precision, dimension(:)  , allocatable :: kMax;
        double precision, dimension(:,:), allocatable :: kSign,  xPointsNorm;
        double precision, dimension(:)  , allocatable :: gammaN, phiN, thetaN, psiN;
        double precision, dimension(:)  , allocatable :: kVec, kVecUnsigned, rVec; !Allocated in function
        double precision, dimension(:)  , allocatable :: deltaK, kDelta;
        double precision, dimension(:)  , allocatable :: xMaxGlob, xMinGlob;
        double precision, dimension(1)                 :: rMax
        logical         , dimension(:)  , allocatable :: effectCalc;
        integer          :: effectComm
        integer          :: i, j, k, m, nDim, N;
        integer          :: xNTotal, kNTotal, rNTotal;
        integer          :: nb_procs, rang, code, error;
        integer          :: pointsPerCorrl
        double precision :: Sk, deltaKprod, step;
        double precision :: pi = 3.1415926535898, zero = 0d0;

		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
		call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)
		nDim    = size(xPoints, 2);
		xNTotal = size(xPoints, 1);

		allocate (xPointsNorm (xNTotal, nDim))
		allocate(kMax   (nDim));
		allocate(kNStep (nDim));
		allocate(kDelta (nDim));
		allocate(kVec        (nDim));
		allocate(rVec        (nDim));

		allocate(effectCalc(Nmc))
		effectCalc(:) = .true.
		if(present(calculate)) effectCalc = calculate

		if(present(communicator)) then
			effectComm = communicator
		else
			effectComm = MPI_COMM_WORLD
		end if

		!Normalization
		xPointsNorm = xPoints
		do i = 1, nDim
			if(corrL(i) /= 1) xPointsNorm(:,i) = xPointsNorm (:,i)/corrL(i)
		end do


		!!
		if((present(MinBound) .and. (.not.present(MaxBound))) .or. &
		   (present(MaxBound) .and. (.not.present(MinBound)))) then
			write(*,*) "ERROR - In 'createStandardGaussianFieldUnstruct': both or none of the bounds should be defined"
			if(.not.present(MinBound)) write(*,*) "MinBound not present"
			if(.not.present(MaxBound)) write(*,*) "MaxBound not present"
			call MPI_ABORT(effectComm  , error, code)
		endif

		!Communicating normalized extremes
        allocate(xMinGlob(nDim))
        allocate(xMaxGlob(nDim))
		if(.not.present(MinBound)) then
			write(*,*) "Min bound not present, the extremes will be calculated automatically"
			call set_Extremes(xPointsNorm, xMinGlob, xMaxGlob, effectComm)
		else
			!write(*,*) "lbound(MinBound) = ", lbound(MinBound)
			!write(*,*) "ubound(MinBound) = ", ubound(MinBound)
			!write(*,*) "lbound(MaxBound) = ", lbound(MaxBound)
			!write(*,*) "ubound(MaxBound) = ", ubound(MaxBound)
			!write(*,*) "lbound(corrL)    = ", lbound(corrL)
			!write(*,*) "ubound(corrL)    = ", ubound(corrL)

			do i = 1, nDim
				!write(*,*) "i = ", i
				xMinGlob(i) = MinBound(i)/corrL(i)
				xMaxGlob(i) = MaxBound(i)/corrL(i)
			end do
		end if
		!!

		!if(rang == 0) write(*,*) ">>>>>>>>> Variables initialization: kMax, kDelta, kNStep, xMinGlob, xMaxGlob";
		!Allocating


		call set_kMaxND(corrMod, kMax) !Defining kMax according to corrMod
		kDelta  = 2*pi/(periodMult*(xMaxGlob - xMinGlob)) !Delta min in between two wave numbers to avoid periodicity
		kNStep  = kAdjust*(ceiling(kMax/kDelta) + 1);
		kNTotal = product(kNStep);

		!N         = size(xPoints,1)
		pointsPerCorrl = 10
		rMax(1)        = sqrt(sum((xMaxGlob - xMinGlob)**2))
		N              = ceiling(rMax(1) * dble(pointsPerCorrl))
		rNTotal        = ceiling(sqrt(dble(N)))

!		if(rang == 0) write(*,*) "N       = ", N;
!		if(rang == 0) write(*,*) "rMax(1) = ",rMax(1);
!		if(rang == 0) write(*,*) "rNTotal = ",rNTotal;
!		if(rang == 0) write(*,*) "size(xPoints,1) = ",size(xPoints,1);
!		if(rang == 0) write(*,*) "size(xPoints,2) = ",size(xPoints,2);
!		call dispCarvalhol(xPoints, "xPoints")
		!Random Field

		randField  = 0;
		step       = rMax(1)/sqrt(dble(N))

		!Initializing the seed
		if(present(chosenSeed)) then
			call init_random_seed(chosenSeed)
		else
			call init_random_seed()
		end if

		if (nDim == 2) then
			allocate(psiN        (rNTotal)); !Out of phase
			allocate(thetaN      (rNTotal));
			allocate(gammaN      (rNTotal));
			do k = 1, Nmc
				if(effectCalc(k)) then
					!if(rang == 0) write(*,*) "k = ",k;
					!write(*,*) "rNTotal = ",rNTotal;
					call random_number(psiN(:))
					call random_number(thetaN(:))
					call random_number(gammaN(:))
					psiN   = 2*pi*psiN
					thetaN = 2*pi*psiN
					gammaN = 2*pi*gammaN

					do j = 1, rNTotal

						rVec           = [cos(thetaN(j)) * j*step, &
						                  sin(thetaN(j)) * j*step]
						Sk             = get_SpectrumND([j*step], corrMod);
						randField(:,k) = sqrt(Sk*j*(step**2)) * gammaN(j) &
						                 * cos(matmul(xPoints(:,:), rVec(:)) + psiN(j)) &
										 + randField(:,k)

	!					if(rang == 0) write(*,*) "rN = ",j, "rVec = ", rVec;
	!					if(rang == 0) write(*,*) "randField(1,k) = ",randField(1,k), "Sk = ", Sk;
					end do
				else
					randField(:,k) = 0.0
				end if
			end do

		else if (nDim == 3) then
			!write(*,*) "nDim = 3 !!!"
			!write(*,*) "k = ",k;
			allocate(psiN   (rNTotal));
			allocate(thetaN (rNTotal));
			allocate(phiN   (rNTotal));
			allocate(gammaN (rNTotal));
			do k = 1, Nmc
				if(effectCalc(k)) then
					!write(*,*) "k = ",k;
					!write(*,*) "rNTotal = ",rNTotal;
					call random_number(phiN(:))
					call random_number(thetaN(:))
					call random_number(gammaN(:))
					call random_number(psiN(:))

					psiN   = 2*pi*psiN
					thetaN = 2*pi*psiN
					phiN   = pi*phiN
					gammaN = 2*pi*gammaN

					do j = 1, rNTotal
						!write(*,*) "j = ", j
						rVec           = [cos(thetaN(j))*sin(phiN(j)) * j*step, &
						                  sin(thetaN(j))*sin(phiN(j)) * j*step, &
						                  cos(phiN(j))                * j*step]
						Sk             = get_SpectrumND([j*step], corrMod);
						randField(:,k) = sqrt(Sk*sin(phiN(j))*step*(j*step)**2) * gammaN(j) &
						                 * cos(matmul(xPoints(:,:), rVec(:)) + psiN(j)) &
										 + randField(:,k)

	!					if(rang == 0) write(*,*) "rN = ",j, "rVec = ", rVec;
	!					if(rang == 0) write(*,*) "randField(1,k) = ",randField(1,k), "Sk = ", Sk;
					end do
				else
					randField(:,k) = 0.0
				end if
			end do
		else
			write(*,*) "The number of dimensions is not accepted in this method (Victor). nDim = ", nDim;
			call MPI_ABORT(MPI_COMM_WORLD, error, code)
		end if

		if(rang == 0) write(*,*) "Spectra (Sk) cut in: ", Sk

		if(allocated(deltaK))       deallocate(deltaK);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
		if(allocated(kVec))         deallocate(kVec);
		if(allocated(kMax))         deallocate(kMax);
		if(allocated(kNStep))       deallocate(kNStep);
		if(allocated(kDelta))       deallocate(kDelta);
		if(allocated(kSign))        deallocate(kSign);
		if(allocated(phiN))         deallocate(phiN);
		if(allocated(psiN))         deallocate(psiN);
		if(allocated(thetaN))       deallocate(thetaN);
		if(allocated(gammaN))       deallocate(gammaN);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
		if(allocated(xMinGlob))     deallocate(xMinGlob);
        if(allocated(xMaxGlob))     deallocate(xMaxGlob);
        if(allocated(xPointsNorm))  deallocate (xPointsNorm);
        if(allocated(rVec))         deallocate(rVec);

	end subroutine createStandardGaussianFieldUnstructVictor

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	subroutine createRandomFieldStructured(xMin, xMax, corrMod, corrL, Nmc, &
										   xNStep, randField);

        implicit none

        !INPUT
        integer                       , intent(in) :: Nmc;
        character (len=*)             , intent(in) :: corrMod;
        double precision, dimension(:), intent(in) :: corrL, xMax, xMin;

        !OUTPUT
        integer         , dimension(:)   , allocatable, intent(out) :: xNStep;
        double precision, dimension(:, :), allocatable, intent(out) :: randField;

        !LOCAL VARIABLES
        integer         , dimension(:)  , allocatable :: kNStep;
        double precision, dimension(:)  , allocatable :: kMax;
        double precision, dimension(:,:), allocatable :: kSign, phiN;
        double precision, dimension(:),   allocatable :: xVec, kVec, kVecUnsigned; !Allocated in function
        double precision, dimension(:),   allocatable :: deltaK, angleVec, kDelta;
        double precision, dimension(:)  , allocatable :: xMaxGlob, xMinGlob;
        integer          :: i, j, k, m, nDim;
        integer          :: xNTotal, kNTotal;
        integer          :: testDim; !Only for tests
        integer          :: nb_procs, rang, code !, xStart, xEnd, sizeUnif, sizeLoc;
        double precision :: xPerCorrL = 5; !number of points per Correlation Length
        double precision :: Sk;
        double precision :: pi = 3.1415926535898, zero = 0d0;

		call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, code)
		call MPI_COMM_RANK(MPI_COMM_WORLD, rang, code)

		!if(rang == 0) write(*,*) "";
        !if(rang == 0) write(*,*) "------------START randomFieldND (proc = ", rang, ")---------";
		!if(rang == 0) write(*,*) "";

		nDim = size(xMax,1);

		!if(rang == 0) write(*,*) ">>>>>>>>> Variables initialization: kMax, xNStep and kNStep";
		!Allocating
		allocate(kMax   (nDim));
		allocate(xNStep (nDim));
		allocate(kNStep (nDim));
		allocate(kDelta (nDim));

		!Communicating extremes
		call set_Extremes(xMin, xMax, xMinGlob, xMaxGlob)

		!Initializing (obs: parameters FOR EACH PROC)
		call set_kMaxND(corrMod, kMax)
		xNStep  = xPerCorrL*ceiling((xMax-xMin)/corrL);
		kDelta  = 2*pi/(periodMult*(xMaxGlob - xMinGlob)) !Delta max in between two wave numbers to avoid periodicity
		kNStep  = kAdjust*(ceiling(kMax/kDelta) + 1);
		xNTotal = product(xNStep);
		kNTotal = product(kNStep);

!		if (rang == 0) then
!			call DispCarvalhol (xMax, "xMax")
!			call DispCarvalhol (xMin, "xMin")
!			call DispCarvalhol (corrL, "corrL")
!			call DispCarvalhol ((xMax-xMin), "(xMax-xMin)")
!			call DispCarvalhol (ceiling((xMax-xMin)/corrL), "ceiling((xMax-xMin)/corrL)")
!			call DispCarvalhol (kMax, "kMax")
!			call DispCarvalhol (xPeriod, "xPeriod")
!			call DispCarvalhol (kNStep, "kNStep")
!		end if



		!Random Field
		allocate(randField((xNTotal),Nmc));
		allocate(angleVec   (xNTotal));
		allocate(deltaK    (nDim));
		allocate(kSign     (2**(nDim-1), nDim));
		allocate(phiN      (size(kSign,1), kNTotal));
		allocate(kVecUnsigned(nDim));
		allocate(kVec(nDim));
		allocate(xVec(nDim));

		if(rang == 0) then
			write(*,*) "";
			write(*,*) ">>>>>>>>> Random Field Creation (only showing proc 0)";
			write(*,*) "";
			write(*,*) "Number of dimensions = ", nDim;
			write(*,*) "Number of events     = ", Nmc;
			write(*,*) "Number x points ";
			write(*,*) "      by Dimension   = ", xNStep;
			write(*,*) "      Total          = ", xNTotal;
			write(*,*) "Number k points ";
			write(*,*) "      by Dimension   = ", kNStep;
			write(*,*) "      Total          = ", kNTotal;
			write(*,*) "";
		end if

		randField = 0;
		angleVec   = 0;
		deltaK    = 0;
	    deltaK    = (kMax)/(kNStep-1); !Defines deltaK
		call set_kSign(kSign)

		do k = 1, Nmc
			call random_number(phiN(:,1:kNTotal))
			do j = 1, kNTotal
				call get_Permutation(j, kMax, kNStep, kVecUnsigned);
			    do m = 1, size(kSign,1)
			    	kVec    = kVecUnsigned * kSign(m, :)
					Sk      = get_SpectrumND(kVec, corrMod, corrL);
					angleVec(:) = 0;
					angleVec(:) = 2 * pi * phiN(m, j) !TO VERIFY - Random part of the angle matrix for this k permutation
					do i = 1, xNTotal
						call get_Permutation(i, xMax, xNStep, xVec, xMin);
						angleVec(i) = angleVec(i) + dot_product(kVec, xVec); !Not-random part of the angle matrix
					end do
					angleVec(:)     = sqrt(Sk*2*product(deltaK)) &
										     *cos(angleVec(:))
					randField(:,k) = randField(:,k) + angleVec(:)
				end do
			end do
			randField(:,k) = product(corrL)*randField(:,k)
			!if(rang == 0) write(*,*) "Event ",k, "of", Nmc, "completed (counting only in proc 0)";
		end do

		randField(:,:) = sqrt(2d0)*randField(:,:)

		!Only printing------------------------------------------

		if(rang == 0) then
			write (*,*) "randField = ", rang
			do i = 1, 10! size(randField, 1)
				write(*,*) randField(i,:)
			end do
		end if

!		call DispCarvalhol(xNStep, "xNStep");
!		call DispCarvalhol(kNStep, "kNStep");
!
!
!		if (rang == 0) then
!			write(*,*) "Permutation X (Nmc = 1), rang ", rang
!			do i = 1, xNTotal
!				call get_Permutation(i, xMax, xNStep, xVec, xMin);
!				write(*,'(I,A,3F10.5)') i, ">", xVec;
!			end do
!		end if
!
!		if (rang == 0) then
!			call DispCarvalhol(kMax  , "kMax")
!			call DispCarvalhol(kNStep, "kNStep")
!			write(*,*) "Permutation K (Nmc = 1)"
!			do i = 1, kNTotal
!			call get_Permutation(i, kMax, kNStep, kVecUnsigned);
!				write(*,*) i, ">", kVecUnsigned;
!			end do
!		end if
!
!		!Spectrum
!		write(*,*) "Spectrum (Nmc = 1)"
!		do i = 1, kNTotal
!			call get_Permutation(j, kMax, kNStep, kVecUnsigned);
!		    do m = 1, size(kSign,1)
!		    	kVec           = kVecUnsigned * kSign(m, :)
!				Sk             = get_SpectrumND(kVec, corrMod, corrL);
!			end do
!		end do
!
!		call DispCarvalhol(randField)

		!---------------------------------------------------

		if(allocated(deltaK))       deallocate(deltaK);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
		if(allocated(kVec))         deallocate(kVec);
		if(allocated(xVec))         deallocate(xVec)
		if(allocated(kMax))         deallocate(kMax);
		if(allocated(kNStep))       deallocate(kNStep);
		if(allocated(angleVec))     deallocate(angleVec);
		if(allocated(kSign))        deallocate(kSign);
		if(allocated(phiN))         deallocate(phiN);
		if(allocated(kVecUnsigned)) deallocate(kVecUnsigned);
        if(allocated(xMaxGlob))     deallocate(xMaxGlob);
        if(allocated(xMinGlob))     deallocate(xMinGlob);
        if(allocated(kDelta))       deallocate(kDelta);


		!if(rang == 0) write(*,*) "";
        !if(rang == 0) write(*,*) "------------END randomFieldND (proc = ", rang, "---------";
		!if(rang == 0) write(*,*) "";
    end subroutine createRandomFieldStructured

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine set_XPoints(corrL, xMin, xMax, xPoints, nPoinsPerCorrL)
    	implicit none

    	!INPUT
    	double precision, dimension(:), intent(in) :: corrL, xMin, xMax;
    	integer                       , intent(in) :: nPoinsPerCorrL !Number of points by correlation length
    	!OUTPUT
    	double precision, dimension(:,:), allocatable, intent(OUT) :: xPoints;
    	!LOCAL VARIABLES
    	integer :: nDim, i, xNTotal;
    	integer , dimension(:) , allocatable :: xNStep;

    	nDim    = size(corrL)
    	allocate(xNStep(nDim))
		xNStep  = nPoinsPerCorrL*ceiling((xMax-xMin)/corrL);
		xNTotal = product(xNStep)
		allocate(xPoints(xNTotal, nDim))

		!call DispCarvalhol(xNStep,"xNStep")

		do i = 1, xNTotal
			call get_Permutation(i, xMax, xNStep, xPoints(i,:), xMin);
		end do

		deallocate(xNStep)

    end subroutine set_XPoints

end module randomFieldND
