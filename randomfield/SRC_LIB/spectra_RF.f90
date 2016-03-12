module spectra_RF
    use displayCarvalhol
    use math_RF
    use write_Log_File
    use constants_RF
    use type_RF

    implicit none
contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kMaxND(corrMod, kMax);
        implicit none
        !INPUT
        integer, intent(in) :: corrMod;
        !double precision,   dimension(:),  intent(in), optional :: corrL;

        !OUTPUT
        double precision, dimension(:),   intent(out) :: kMax;

        !LOCAL VARIABLES
        integer          :: nDim

        nDim = size(kMax)

        select case(corrMod)
            case(cm_GAUSSIAN)
                select case(nDim)
                    case(1)
                        kMax(:) = 6.457D0; !Value to cover 99% Spectra area
                    case(2)
                        kMax(:) = 7.035D0; !Value to cover 99% Spectra area
                    case(3)
                        kMax(:) = 7.355D0; !Value to cover 99% Spectra area
                end select
        end select

        !write(get_fileId(),*) "kMax = ", kMax


    end subroutine set_kMaxND

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kPoints(RDF, xStep);
        implicit none

        !INPUT OUTPUT
        type(RF) :: RDF
        !INPUT
        double precision, dimension(:), intent(in) :: xStep

        !LOCAL
        integer :: i
        integer(kind=8) :: i_long, kNLocal
        double precision :: kAdjust    = 1.0D0 !"kNStep minimum" multiplier
        double precision :: periodMult = 1.1D0 !"range" multiplier
        double precision :: rAdjust    = 5.0D0 !"kNStep minimum" multiplier (isotropic)

        if(allocated(RDF%kPoints)) deallocate(RDF%kPoints)

        call set_kMaxND(RDF%corrMod, RDF%kMax) !Defining kMax according to corrMod

        select case (RDF%method)
            case(ISOTROPIC)
                RDF%kDelta(1) = 2.0D0*PI/(periodMult*sqrt(sum(RDF%xRange**2))) !Diagonal
                RDF%kNStep(1) = 1 + int(rAdjust*(ceiling(maxval(RDF%kMax)/RDF%kDelta(1)))); !Number of points in k
                RDF%kDelta(1) = maxval(RDF%kMax)/(RDF%kNStep(1)-1); !Redefining kDelta after ceiling and adjust
                RDF%kNTotal   = RDF%kNStep(1);

                allocate(RDF%kPoints(1, RDF%kNTotal))

                do i_long = 1, RDF%kNTotal
                    call get_Permutation(i, [RDF%kMax(1)], [RDF%kNStep(1)], RDF%kPoints(:, i), snapExtremes = .true.);
                end do

            case(SHINOZUKA)
                RDF%kNStep(:)   = 1 + int(kAdjust*(ceiling(RDF%kMax/RDF%kDelta(:)))); !Number of points in k
                RDF%kDelta(:) = (RDF%kMax)/(RDF%kNStep-1); !Redefining kDelta after ceiling and adjust
                RDF%kNTotal = product(int(RDF%kNStep,8));

                !!write(get_fileId(),*) "RDF%kNStep = ", RDF%kNStep
                !!write(get_fileId(),*) "RDF%kNTotal = ", RDF%kNTotal

                allocate(RDF%kPoints(RDF%nDim, RDF%kNTotal))

                do i_long = 1, RDF%kNTotal
                    call get_Permutation(i, RDF%kMax, RDF%kNStep, RDF%kPoints(:, i), snapExtremes = .true.);
                end do

            case(RANDOMIZATION)
                RDF%kNStep(:)   = 1 + int(kAdjust*(ceiling(RDF%kMax/RDF%kDelta(:)))); !Number of points in k
                RDF%kDelta(:) = (RDF%kMax)/(RDF%kNStep-1); !Redefining kDelta after ceiling and adjust
                RDF%kNTotal = product(int(RDF%kNStep,8));

                allocate(RDF%kPoints(RDF%nDim, RDF%kNTotal))
                call random_number(RDF%kPoints(:,:))
                do i = 1, RDF%nDim
                    RDF%kPoints(i,:) = RDF%kPoints(i,:) * RDF%kMax(i)
                end do

            case(FFT)
                RDF%kNStep(:) = find_xNStep(xMaxExt = RDF%xRange, xStep=xStep)
                !RDF%kMax(:)   = (dble(RDF%kNStep(:) - 1)/((RDF%xRange)**1.0D0) * RDF%kDelta(:) !Redefinition of kMax (divided by 2 because of symmetric plane)
                !RDF%kMax(:)   = 2.0D0*PI/(periodMult*RDF%xRange)
                !RDF%kDelta(:) = (RDF%kMax)/(RDF%kNStep-1);
                RDF%kDelta(:) = 2.0D0*PI/(periodMult*RDF%xRange)
                RDF%kMax(:)   = ((RDF%kNStep(:) - 1) * RDF%kDelta(:))/1.0D0
                !RDF%kDelta(:) = RDF%kDelta(:)
                !RDF%kDelta(:) = RDF%kMax(:)/dble(RDF%kNStep(:))
                !RDF%kNStep(:) = RDF%xNStep(:);
                RDF%kNTotal   = product(int(RDF%kNStep,8));
                !RDF%kMax(:)   = (dble(RDF%kNStep(:) - 1)/(2.0D0)) * RDF%kDelta(:)/RDF%corrL(:)!Redefinition of kMax (divided by 2 because of symmetric plane)
                !RDF%kMax(:)   = (dble(RDF%kNStep(:) - 1)/(1.0D0)) * RDF%kDelta(:)!Redefinition of kMax (divided by 2 because of symmetric plane)

                !if(RDF%independent) then
                !    call wLog("RDF%kNStep = ")
                !    call wLog(RDF%kNStep)
                !    call wLog("RDF%kMax = ")
                !    call wLog(RDF%kMax)
                !    call wLog("RDF%kNTotal = ")
                !    allocate(RDF%kPoints(RDF%nDim, RDF%kNTotal))
                !    call wLog("shape(RDF%kPoints) = ")
                !    call wLog(shape(RDF%kPoints))
                !    call wLog(RDF%kNTotal)
                !    !call wLog("RDF%kPoints = ")
                !    do i = 1, RDF%kNTotal
                !        call get_Permutation(i, RDF%kMax, RDF%kNStep, RDF%kPoints(:, i), snapExtremes = .true.);
                !        !call wLog(RDF%kPoints(:, i))
                !    end do
                !else
                    call wLog("RDF%kNInit = ")
                    call wLog(RDF%kNInit)
                    call wLog("RDF%kNEnd = ")
                    call wLog(RDF%kNEnd)
                    kNLocal = RDF%kNEnd - RDF%kNInit + 1
                    call wLog("kNLocal = ")
                    call wLog(kNLocal)
                    call wLog("HERE !!!!!!!!!!!!")
                    allocate(RDF%kPoints(RDF%nDim, kNLocal))
                    call wLog("shape(RDF%kPoints) = ")
                    call wLog(shape(RDF%kPoints))
                    call wLog("Making kPoints= ")
                    !call setGrid(RDF%kPoints, dble(RDF%kNInit-1)*RDF%kDelta, RDF%kDelta, RDF%kNStep)
                    do i_long = RDF%kNInit, RDF%kNEnd
                        call get_Permutation(int(i_long), RDF%kMax, RDF%kNStep, &
                                             RDF%kPoints(:, i_long-RDF%kNInit+1), &
                                             snapExtremes = .true.)
                    end do
                !end if

        end select

    end subroutine set_kPoints

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_SkVec(RDF, corrL_in);
        implicit none

        !OBS: corrL is supposed = 1. The complete formula is RDF%SkVec = product(corrL) * exp(-dot_product(kVector**2, corrL_effec**2)/(4.0d0*pi))

        !INPUT OUTPUT
        type(RF) :: RDF
        double precision, dimension(:), optional, intent(in) :: corrL_in
        !LOCAL
        integer :: i, freqK = 6
        double precision, dimension(RDF%nDim) :: corrL

        corrL(:) = 1.0D0
        if(present(corrL_in)) corrL = corrL_in

        if(allocated(RDF%SkVec)) deallocate(RDF%SkVec)
        allocate(RDF%SkVec(size(RDF%kPoints,2)))


        call wLog("RDF%corrMod")
        call wLog(RDF%corrMod)

        !write(*,*) "Inside set_SkVec, corrL = ", corrL

        select case(RDF%corrMod)

            case(cm_GAUSSIAN)
                RDF%SkVec(:) = 1.0D0
                !write(*,*) "Gaussian Correlation Model"
                call wLog("cm_GAUSSIAN")
                do i = 1, RDF%nDim
                    RDF%SkVec(:) = RDF%SkVec(:) * corrL(i) * exp(-((RDF%kPoints(i,:)**2.0D0) * (corrL(i)**2.0D0))/(4.0d0*pi))
                end do
                !call DispCarvalhol(RDF%SkVec, "RDF%SkVec")
            case(cm_COS)
                call wLog("cm_COS")
                RDF%SkVec = 0
                if(size(RDF%SkVec)>freqK) then
                    RDF%SkVec(freqK) = 1
                    call wLog("kPoint used in cosinus CM: ")
                    call wLog(RDF%kPoints(:,freqK))
                    if(RDF%rang == 0) write(*,*) "kPoint used in cosinus CM: ", RDF%kPoints(:,freqK)
                else
                    stop("When using the cosinus correlation Model you should have more than 'freqK' kPoints")
                end if
        end select

    end subroutine set_SkVec

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_rMax(corrMod, rMax, corrL);
        implicit none
        !INPUT
        integer,                intent(in) :: corrMod;
        double precision,   dimension(:),  intent(in), optional :: corrL;

        !OUTPUT
        double precision, dimension(:),   intent(out) :: rMax;

        !LOCAL VARIABLES
        double precision :: pi = 3.1415926535898
        double precision, dimension(:), allocatable:: corrL_effec

        allocate(corrL_effec (size(rMax)))

        if (present(corrL)) corrL_effec = corrL
        if (.not. present(corrL)) corrL_effec = 1
        !        do i = 1, 100
        !            kMax = i/10.0 * corrL(:)
        !            write(*,*) "kMax = ", kMax
        !            write(*,*) "Spectrum = ", get_SpectrumND(kMax, corrMod, corrL)
        !            call DispCarvalhol (kMax, "kMax")
        !            call DispCarvalhol (get_SpectrumND(kMax, corrMod, corrL), "Spectrum")
        !        end do

        select case(corrMod)
        case(cm_GAUSSIAN)
            rMax(:) = 2*pi*corrL_effec(:); !CRITERIA STILL TO BE TESTED
        end select

        deallocate(corrL_effec)

    end subroutine set_rMax

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kSign(kSign);
        implicit none
        !INPUT and OUTPUT
        double precision,   dimension(:, :), intent(inout) :: kSign;

        !LOCAL VARIABLES
        integer:: i, j, pos, seedStep, nDim;

        nDim = size(kSign, 2);

        kSign(:,:) = 1d0

        do i = 2, nDim
            do j = 1, size(kSign,1)
                seedStep = 2**(nDim-i);
                pos      = cyclicMod(int((j-0.9)/seedStep), 2)
                if (mod(pos,2) == 1) kSign(j,i) = -1d0
            end do
        end do
    end subroutine set_kSign

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kArray_rand(corrMod, kArray_rand, nDim);
        implicit none
        !INPUT
        integer, intent(in) :: corrMod;
        integer :: nDim

        !OUTPUT
        double precision, dimension(:),  intent(out) :: kArray_rand;

        !LOCAL VARIABLES
        integer :: i
        double precision ::  bound,mean,cdf_x,q,sd,uniRand
        integer :: st, which
        !real :: av = 0.0
        !real :: gennor
        !real :: std = 1.0
        !real :: normRand

!        mean = 0.0D0
!        sd = 1.0D0
!        test2 = gennor (real(mean), real(sd))
!        write(*,*) "gennor (0.0, 1.0) = ", test2

!        call random_number(kArray_rand)
!        kArray_rand = kMaxScal*kArray_rand

!TODO RANDOMIZATION
        which = 2 !Find x, for a given cdf_x

        !write(*,*) "corrMod = ", corrMod

        select case(corrMod)
            case(cm_GAUSSIAN)
                !write(*,*) "Gaussian kArray"
                mean = 0.0D0
                sd = SQRT_2PI
                !sd = 1.0D0/SQRT_2PI
                do i = 1, size(kArray_rand)
                    call random_number(uniRand)
                    cdf_x = (0.5D0 + uniRand/2.0D0)
                    q = 1-cdf_x
                    call cdfnor(which,cdf_x,q,kArray_rand(i),mean,sd,st,bound) !Normal CDF
                    kArray_rand = kArray_rand**(dble(nDim))
                end do

!                sd = SQRT_2PI
!                !sd = 1.0D0/SQRT_2PI
!                do i = 1, size(kArray_rand)
!                    cdf_x = 1.0D0
!                    do j = 1, nDim
!                        call random_number(uniRand)
!                        cdf_x = cdf_x*(0.5D0 + uniRand/2.0D0)
!                    end do
!                    q = 1-cdf_x
!                    call cdfnor(which,cdf_x,q,kArray_rand(i),mean,sd,st,bound) !Normal CDF
!                end do


                !write(*,*) "gennor (mean, sd) = ",gennor (0.5, 1.0)
                ! PROBLEM, how to put the seed
!                mean = 0.0D0
!                sd = 1.0D0/SQRT_2PI
!
!                kArray_rand(:) = 1
!                do i = 1, size(kArray_rand)
!                    do j = 1, nDim
!                        normRand = gennor (real(mean), real(sd))
!                        kArray_rand(i) = kArray_rand(i)*abs(dble(normRand))
!                    end do
!                end do
        end select



        !call reorder_vector(kArray_rand)

        !kArray_rand = kArray_rand**(dble(nDim))

        !call dispCarvalhol(kArray, "kArray")

    end subroutine set_kArray_rand

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_kDelta_rand(kArray_rand, kDelta_rand);
        implicit none
        !INPUT
        double precision, dimension(:),  intent(in) :: kArray_rand;

        !OUTPUT
        double precision, dimension(:),  intent(out) :: kDelta_rand;

        !LOCAL VARIABLES
        integer :: i
        double precision :: inc

        kDelta_rand = 0.0D0

        !Correction for the first term
        inc = (kArray_rand(1) - 0.0D0)
        kDelta_rand(1) = kDelta_rand(1) + inc

        !Calculating distance between numbers (we supose they're ordered)
        do i = 1, size(kArray_rand)-1
            inc = (kArray_rand(i+1) - kArray_rand(i))/2.0D0
            kDelta_rand(i)   = kDelta_rand(i)   + inc
            kDelta_rand(i+1) = kDelta_rand(i+1) + inc
        end do

    end subroutine set_kDelta_rand

end module spectra_RF
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
