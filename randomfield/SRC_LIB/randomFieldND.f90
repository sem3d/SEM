module randomFieldND

    use displayCarvalhol
    use spectra_RF
    use math_RF
    use constants_RF
    use mesh_RF
    use mpi
    use write_Log_File
    use type_RF
    use type_MESH
    use common_variables_RF
    use writeResultFile_RF
    use fftw3
    !use blas
    implicit none




contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine create_RF_Unstruct_Init (RDF, MSH)
        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        type(MESH), intent(inout) :: MSH
        !LOCAL

        !if(RDF%rang == 0) write(*,*) "Inside create_RF_Unstruct_Init"

        !Generating standard Gaussian Field
        call gen_Std_Gauss(RDF, MSH)

    end subroutine create_RF_Unstruct_Init

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss (RDF, MSH)
        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        type(MESH), intent(inout) :: MSH

        !LOCAL VARIABLES
        integer :: i;

        !logical, dimension(size(MSH%neigh)) :: considerNeighbour
        !integer, dimension(16) :: testVec
        !integer :: partitionType = 1


        !testVec = [(i, i = 1, 16)]

        !Normalization
        call wLog(" ")
        call wLog("->Normalizing Coordinates")
        call wLog(" ")
        do i = 1, RDF%nDim
            if(associated(RDF%xPoints)) RDF%xPoints(i,:)   = RDF%xPoints(i,:)/RDF%corrL(i)
            MSH%xStep(i)       = MSH%xStep(i) /RDF%corrL(i)
            MSH%xMinInt(i)     = MSH%xMinInt(i)/RDF%corrL(i)
            MSH%xMaxInt(i)     = MSH%xMaxInt(i)/RDF%corrL(i)
            MSH%xMinExt(i)     = MSH%xMinExt(i)/RDF%corrL(i)
            MSH%xMaxExt(i)     = MSH%xMaxExt(i)/RDF%corrL(i)
            MSH%xMinGlob(i)    = MSH%xMinGlob(i)/RDF%corrL(i)
            MSH%xMaxGlob(i)    = MSH%xMaxGlob(i)/RDF%corrL(i)
            MSH%xMaxNeigh(i,:) = MSH%xMaxNeigh(i,:)/RDF%corrL(i)
            MSH%xMinNeigh(i,:) = MSH%xMinNeigh(i,:)/RDF%corrL(i)
            MSH%xMaxBound(i)   = MSH%xMaxBound(i)/RDF%corrL(i)
            MSH%xMinBound(i)   = MSH%xMinBound(i)/RDF%corrL(i)
            MSH%xOrNeigh(i,:)  = MSH%xOrNeigh(i,:)/RDF%corrL(i)
            RDF%xRange(i)      = RDF%xRange(i)/RDF%corrL(i)
        end do

        !Generating Standard Gaussian Field
        call wLog("")
        call wLog("GENERATING RANDOM FIELDS")
        call wLog("-------------------------------")
        call wLog("")

        select case (RDF%method)
            case(ISOTROPIC)
                call wLog(" ISOTROPIC")
                !if(RDF%rang == 0) write(*,*)"ISOTROPIC"
                call gen_Std_Gauss_Isotropic(RDF, MSH)
            case(SHINOZUKA)
                call wLog(" SHINOZUKA")
                !if(RDF%rang == 0) write(*,*)"SHINOZUKA"
                call gen_Std_Gauss_Shinozuka(RDF, MSH)
            case(RANDOMIZATION)
                call wLog(" RANDOMIZATION")
                if(RDF%rang == 0) write(*,*)"RANDOMIZATION"
                call gen_Std_Gauss_Randomization(RDF, MSH)
            case(FFT)
                call wLog(" FFT")
                write(*,*) "-> FFT ", RDF%rang
                !if(RDF%rang == 0) write(*,*)"FFT"
                call gen_Std_Gauss_FFT(RDF, MSH)
                write(*,*) "-> AFTER FFT ", RDF%rang
        end select

        !RDF%randField = 1.0 ! For Tests

        call wLog("minval(RDF%randField,1) =")
        call wLog(minval(RDF%randField,1))
        call wLog("maxval(RDF%randField,1) =")
        call wLog(maxval(RDF%randField,1))

        !Reverting Normalization
        call wLog(" ")
        call wLog("->Reverting Normalization")
        do i = 1, RDF%nDim
            if(associated(RDF%xPoints)) RDF%xPoints(i,:)   = RDF%xPoints(i,:)*RDF%corrL(i)
            RDF%xRange(i)      = RDF%xRange(i)*RDF%corrL(i)
            MSH%xStep(i)       = MSH%xStep(i)*RDF%corrL(i)
            MSH%xMinInt(i)     = MSH%xMinInt(i)*RDF%corrL(i)
            MSH%xMaxInt(i)     = MSH%xMaxInt(i)*RDF%corrL(i)
            MSH%xMinExt(i)     = MSH%xMinExt(i)*RDF%corrL(i)
            MSH%xMaxExt(i)     = MSH%xMaxExt(i)*RDF%corrL(i)
            MSH%xMinGlob(i)    = MSH%xMinGlob(i)*RDF%corrL(i)
            MSH%xMaxGlob(i)    = MSH%xMaxGlob(i)*RDF%corrL(i)
            MSH%xMaxNeigh(i,:) = MSH%xMaxNeigh(i,:)*RDF%corrL(i)
            MSH%xMinNeigh(i,:) = MSH%xMinNeigh(i,:)*RDF%corrL(i)
            MSH%xMaxBound(i)   = MSH%xMaxBound(i)*RDF%corrL(i)
            MSH%xMinBound(i)   = MSH%xMinBound(i)*RDF%corrL(i)
            MSH%xOrNeigh(i,:)  = MSH%xOrNeigh(i,:)*RDF%corrL(i)
            RDF%xRange(i)      = RDF%xRange(i)*RDF%corrL(i)
        end do

        !RDF%randField = RDF%rang ! For Tests

    end subroutine gen_Std_Gauss


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine multiVariateTransformation (margiFirst, fieldAvg, fieldVar, randField)

        implicit none

        !INPUT
        integer          , intent(in) :: margiFirst;
        double precision , intent(in) :: fieldAvg, fieldVar;

        !OUTPUT (IN)
        double precision, dimension(1:, 1:), intent(inout) :: randField;

        !LOCAL VARIABLES
        double precision :: normalVar, normalAvg
        integer          :: error, code

        select case (margiFirst)
        case(fom_GAUSSIAN)
            normalVar = fieldVar
            normalAvg = fieldAvg
        case(fom_LOGNORMAL)
            if(fieldAvg <= 0.0D0) then
                write(*,*) ""
                write(*,*) "ERROR - when using lognormal fieldAvg should be a positive number"
                call MPI_ABORT(MPI_COMM_WORLD, error, code)
            end if
            normalVar = log(1 + fieldVar/(fieldAvg**2))
            normalAvg = log(fieldAvg) - normalVar/2
        end select

        randField(:,:) = randField(:,:) * sqrt(normalVar) &
            + normalAvg;

        if (margiFirst == fom_LOGNORMAL) then
            randField(:,:) = exp(randField(:,:))
        end if

    end subroutine multiVariateTransformation

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss_Shinozuka(RDF, MSH, randomK_in)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        logical, intent(in), optional ::randomK_in
        !INPUT
        type(MESH), intent(in) :: MSH

        !LOCAL
        double precision, dimension(:, :), allocatable :: phiK, kVec;
        double precision, dimension(:)   , allocatable :: dgemm_mult;
        double precision, dimension(:,:) , allocatable :: k_x_phi, kSign;
        double precision :: ampMult

        integer :: n, i, m
        integer(kind=8) :: j
        logical :: randomK
        integer(kind=8) :: xNTotal

        !!write(get_fileId(),*) "Inside Shinozuka"

        xNTotal = product(MSH%xNStep)
        randomK = .false.
        if(present(randomK_in)) randomK = randomK_in
        call init_random_seed(RDF%seed)

        !!write(get_fileId(),*) "Defining kPoints and SkVec"
        call set_kPoints(RDF, MSH%xStep)
        call set_SkVec(RDF)

        !!write(get_fileId(),*) "Calculating Fields"
        allocate(k_x_phi (xNTotal, 1))
        allocate(kSign (2**(RDF%nDim-1), RDF%nDim));
        allocate(kVec(RDF%nDim, 1))

        call set_kSign(kSign) !Set the sign permutations for kVec

        allocate(phiK (RDF%kNTotal, size(kSign,1)));

        if(randomK) then
            !write(get_fileId(),*) "-----Shinozuka, k random-----"
            ampMult = 2.0d0*sqrt(1/(RDF%kNTotal*(2.0d0*PI)**(dble(RDF%nDim))))
        else
            !write(get_fileId(),*) "-----Shinozuka, k discrete-----"
            ampMult = 2.0d0*sqrt(product(RDF%kDelta)/((2.0d0*PI)**(dble(RDF%nDim))))
        end if

        !write(get_fileId(),*) "     kNStep  = ", RDF%kNStep
        !write(get_fileId(),*) "     kNTotal = ", RDF%kNTotal
        !write(get_fileId(),*) "     xNTotal = ", size(RDF%xPoints, 2)

        RDF%randField(:,:) = 0.0d0;

        do n = 1, RDF%Nmc
            !write(get_fileId(),*) "  --Generating Field Number ", n

            if(.not. RDF%calculate(n)) cycle

            call random_number(phiK(:,:))
            phiK(:,:) = 2.0D0*pi*phiK(:,:)

            !write(get_fileId(),*) "     First PhiK = ", phiK(1,1)
            !write(get_fileId(),*) "     Last PhiK  = ", phiK(size(phiK,1), size(phiK,2))

            !Loop on k sign
            do m = 1, size(kSign,1)

                !Changing kPoints Sign
                do i = 1, RDF%nDim
                    if(kSign(m, i) == -1) RDF%kPoints(i,:) = -RDF%kPoints(i,:)
                end do

                !Loop on k
                do j = 1, RDF%kNTotal

                    call DGEMM_simple(RDF%xPoints, RDF%kPoints(:,j:j), k_x_phi(:,:), "T", "N") !x*k

                    RDF%randField(:,n) = sqrt(RDF%SkVec(j)) * &
                                         cos(k_x_phi(:,1) + phiK(j, m)) &
                                         + RDF%randField(:,n)

                end do !END Loop on k

                !Reverting kPoints Sign
                do i = 1, RDF%nDim
                    if(kSign(m, i) == -1) RDF%kPoints(i,:) = - RDF%kPoints(i,:)
                end do

            end do !END Loop on k sign

        end do !END Loop on Nmc

        RDF%randField(:,:) = ampMult * RDF%randField(:,:);

        if(allocated(dgemm_mult))       deallocate(dgemm_mult)
        if(allocated(phiK))             deallocate(phiK);
        if(allocated(k_x_phi))          deallocate(k_x_phi)
        if(allocated(kSign))            deallocate(kSign)

    end subroutine gen_Std_Gauss_Shinozuka

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss_Randomization(RDF, MSH)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        !INPUT
        type(MESH), intent(in) :: MSH

        !LOCAL
        logical :: randomK;

        !write(get_fileId(),*) "-----Randomization-----"

        randomK = .true.

        call gen_Std_Gauss_Shinozuka(RDF, MSH, randomK)

    end subroutine gen_Std_Gauss_Randomization


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss_Isotropic(RDF, MSH)

        implicit none

        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        !INPUT
        type(MESH), intent(in) :: MSH

        !LOCAL
        double precision, dimension(:)  , allocatable :: gammaN, phiN, thetaN, psiN;
        double precision, dimension(RDF%nDim) :: rVec;

        !double precision :: rMax, Sk
        !double precision, dimension(1) :: rMaxVec
        integer          :: k
        integer(kind=8)  :: j
        !integer          :: rNTotal;
        double precision :: step;
        double precision, dimension(MSH%xNTotal) :: dgemm_mult;

        call wLog("-----Inside Isotropic-----")

        write(*,*) "Setting kPoints"
        call set_kPoints(RDF, MSH%xStep)
        write(*,*) "Setting SkVec"
        call set_SkVec(RDF)

        call wLog("RDF%kPoints------------")
        call wLog(RDF%kPoints)
        call wLog("RDF%SkVec--------------")
        call wLog(RDF%SkVec)

!        !r Definition
!        call set_kMaxND(RDF%corrMod, rMaxVec)
!        call set_kPoints(RDF)
!        rMax = rMaxVec(1)
!        rDelta  = maxval(RDF%kDelta(:))/5.0D0 !Delta min in between two wave numbers to avoid periodicity
!        rNTotal = ceiling(rMax/rDelta) + 1;
!
!        !Generating random field samples
        step = RDF%kPoints(1,2)-RDF%kPoints(1,1)
        call wLog("step = ")
        call wLog(step)
!        RDF%randField(:,:) = 0.0D0;
!
        call init_random_seed(RDF%seed)

        if (RDF%nDim == 2) then
            allocate(psiN   (RDF%kNTotal));
            allocate(thetaN (RDF%kNTotal));
            allocate(gammaN (RDF%kNTotal));
            do k = 1, RDF%Nmc
                if(RDF%calculate(k)) then
                    write(*,*) "k = ",k;
                    !write(*,*) "dgemm_mult"
                    call random_number(psiN(:))
                    call random_number(thetaN(:))
                    call random_number(gammaN(:))
                    psiN   = 2d0*pi*psiN
                    thetaN = 2d0*pi*psiN
                    gammaN = 2d0*pi*gammaN

                    do j = 1, RDF%kNTotal
                        rVec           = [cos(thetaN(j)) * RDF%kPoints(1,j), &
                                          sin(thetaN(j)) * RDF%kPoints(1,j)]
                        call DGEMM ( "T", "N", size(RDF%randField,1), 1, RDF%nDim, &
                            1.0d0, RDF%xPoints, RDF%nDim, rVec, RDF%nDim, 0.0d0, dgemm_mult, size(RDF%randField,1))
                        call wLog(dgemm_mult)
                        RDF%randField(:,k) = sqrt(RDF%SkVec(j)*RDF%kPoints(1,j)) * gammaN(j) &
                                            * cos(                           &
                                            dgemm_mult                &
                                            + psiN(j)                 &
                                            )                          &
                                            + RDF%randField(:,k)
                    end do
                else
                    RDF%randField(:,k) = 0.0
                end if
            end do

        else if (RDF%nDim == 3) then
            !write(*,*) "nDim = 3 !!!"
            !write(*,*) "k = ",k;
            allocate(psiN   (RDF%kNTotal));
            allocate(thetaN (RDF%kNTotal));
            allocate(phiN   (RDF%kNTotal));
            allocate(gammaN (RDF%kNTotal));
            do k = 1, RDF%Nmc
                if(RDF%calculate(k)) then
                    !write(*,*) "k = ",k, "-------------------------";
                    !write(*,*) "rNTotal = ",rNTotal;
                    call random_number(phiN(:))
                    call random_number(thetaN(:))
                    call random_number(gammaN(:))
                    call random_number(psiN(:))

                    psiN   = 2*pi*psiN
                    thetaN = 2*pi*psiN
                    phiN   = pi*phiN
                    gammaN = sqrt(12.0)*(gammaN -0.5d0)



                    do j = 1, RDF%kNTotal
                        !write(*,*) "j = ", j
                        rVec           = [cos(thetaN(j))*sin(phiN(j)) * RDF%kPoints(1,j), &
                                          sin(thetaN(j))*sin(phiN(j)) * RDF%kPoints(1,j), &
                                          cos(phiN(j))                * RDF%kPoints(1,j)]
                        !call wLog(rVec)
                        call DGEMM ( "T", "N", size(RDF%randField,1), 1, RDF%nDim, &
                                    1.0d0, RDF%xPoints, RDF%nDim, rVec, RDF%nDim, 0.0d0, dgemm_mult, size(RDF%randField,1))
                        RDF%randField(:,k) = sqrt(RDF%SkVec(j)*sin(phiN(j))*(RDF%kPoints(1,j))**2) * gammaN(j) &
                                              * cos(                                             &
                                              dgemm_mult                                   &
                                              + psiN(j)                                    &
                                              )                                            &
                                              + RDF%randField(:,k)
                    end do
                else
                    RDF%randField(:,k) = 0.0
                end if
            end do

        else
!            write(*,*) "ERROR The number of dimensions is not accepted in this method (Isotropic)";
!            write(*,*) "RDF%nDim = ", RDF%nDim;
!            stop
        end if
!
!        !if(rang == 0) write(*,*) "Spectra (Sk) cut in: ", Sk
!
        RDF%randField(:,:) = sqrt((step)/((2.0d0*pi)**(RDF%nDim)))&
                             * RDF%randField(:,:)

        call wLog("RDF%randField--------------")
        call wLog(RDF%randField)
!
!        RDF%randField = 1.0 ! For Tests
        !RDF%randField = RDF%rang ! For Tests

        if(allocated(phiN))         deallocate(phiN);
        if(allocated(psiN))         deallocate(psiN);
        if(allocated(thetaN))       deallocate(thetaN);
        if(allocated(gammaN))       deallocate(gammaN);

    end subroutine gen_Std_Gauss_Isotropic

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine gen_Std_Gauss_FFT(RDF, MSH)

        implicit none
        !INPUT OUTPUT
        type(RF), intent(inout) :: RDF
        !INPUT
        type(MESH), intent(in) :: MSH
        !LOCAL
        integer(C_INTPTR_T) :: L, M, N
        integer(C_INTPTR_T) :: local_LastDim
        integer(C_INTPTR_T) :: local_LD_offset
        real(C_DOUBLE), pointer :: data_real_2D(:,:), data_real_3D(:,:,:)
        type(C_PTR) :: cdata, plan
        integer(C_INTPTR_T) :: alloc_local
        integer :: sliceSize
        integer :: kNLocal

        double precision, dimension(:), allocatable :: gammaK, phiK
        integer(kind=8) :: i_long, kNCumulated_Init, kNCumulated_End
        double precision :: trashNumber
        integer, dimension(RDF%nDim) :: xNStepGlob
        double precision :: ampMult

        call wLog("gen_Std_Gauss_FFT_init")
        call wLog(" ")


        L = MSH%xNStep(1)
        if(RDF%nDim >= 2) M = MSH%xNStep(2)
        if(RDF%nDim >= 3) N = MSH%xNStep(3)

        call wLog("    GLOBAL")
        call fftw_mpi_init()

        xNStepGlob = find_xNStep(MSH%xMinGlob, MSH%xMaxGlob, MSH%xStep)
        call wLog("MSH%xMinGlob = ")
        call wLog(MSH%xMinGlob)
        call wLog("MSH%xMaxGlob = ")
        call wLog(MSH%xMaxGlob)
        call wLog("MSH%xStep = ")
        call wLog(MSH%xStep)
        call wLog("xNStepGlob = ")
        call wLog(xNStepGlob)

        if(RDF%nDim == 2) then
            L = xNStepGlob(1)
            M = xNStepGlob(2)
            alloc_local = fftw_mpi_local_size_2d(M, L, RDF%comm, &
                                                 local_LastDim, local_LD_offset) !FOR MPI
            cdata = fftw_alloc_real(alloc_local)
            call c_f_pointer(cdata, data_real_2D, [L, local_LastDim])
            call wLog("L = ")
            call wLog(L)
            call wLog("M = ")
            call wLog(L)

        else if(RDF%nDim == 3) then
            L = xNStepGlob(1)
            M = xNStepGlob(2)
            N = xNStepGlob(3)
            alloc_local = fftw_mpi_local_size_3d(N, M, L, RDF%comm, &
                                                 local_LastDim, local_LD_offset) !FOR MPI
            cdata = fftw_alloc_real(alloc_local)
            call c_f_pointer(cdata, data_real_3D, [L, M, local_LastDim])
            call wLog("L = ")
            call wLog(L)
            call wLog("M = ")
            call wLog(L)
            call wLog("N = ")
            call wLog(L)
            write(*,*) "In rang ", RDF%rang, " L=",L, "M=",M, "N=",N &
                       , "local_LastDim=", local_LastDim

        else
            stop("Inside gen_Std_Gauss_FFT dimension not yet implemented for this generation method")

        end if

        !Defining kInit and kEnd
        RDF%kNInit = int(local_LD_offset) + 1
        RDF%kNEnd  = RDF%kNInit + int(local_LastDim) - 1
        call wLog("local_LD_offset")
        call wLog(local_LD_offset)
        call wLog("local_LastDim")
        call wLog(local_LastDim)

        if(local_LastDim < 20) then
            write(*,*) "WARNING, local_LastDim = ", local_LastDim&
                       , "problems have been observed when a FFTW slab is this thin"
            call wLog("WARNING, local_LastDim = ")
            call wLog(local_LastDim)
            call wLog("problems have been observed when a FFTW slab is this thin")
        end if

        call wLog("RDF%kNInit")
        call wLog(RDF%kNInit)
        call wLog("RDF%kNEnd")
        call wLog(RDF%kNEnd)

        sliceSize = 1
        if(RDF%nDim > 1) sliceSize = product(MSH%xNStep(1:RDF%nDim -1))
        kNCumulated_Init = (RDF%kNInit - 1) * sliceSize + 1
        kNCumulated_End = RDF%kNEnd *sliceSize

        !call wLog("kNCumulated_Init")
        !call wLog(kNCumulated_Init)
        !call wLog("kNCumulated_End")
        !call wLog(kNCumulated_End)

        !RDF%origin  = [1, int(local_LD_offset) + 1]
        !RDF%kExtent = [L , local_M]
        !RDF%origin = 1
        !RDF%origin(RDF%nDim) = int(local_LD_offset) + 1
        RDF%kExtent = MSH%xNStep
        RDF%kExtent(RDF%nDim) = int(local_LastDim)

        call wLog("MSH%origin")
        call wLog(MSH%origin)
        call wLog("MSH%origin (IDEAL)")
        if(MSH%nDim==2) call wLog([1, int(local_LD_offset) + 1])
        if(MSH%nDim==3) call wLog([1, 1, int(local_LD_offset) + 1])
        call wLog("RDF%kExtent")
        call wLog(RDF%kExtent)
        call wLog("RDF%kExtent (IDEAL)")
        if(MSH%nDim==2) call wLog([int(L), int(local_LastDim)])
        if(MSH%nDim==3) call wLog([int(L), int(M), int(local_LastDim)])

        call wLog("Before set_kPoints")
        call set_kPoints(RDF, MSH%xStep)
        call set_SkVec(RDF)

        !call wLog("RDF%kPoints")
        !call DispCarvalhol(RDF%kPoints, unit_in = RDF%log_ID)
        !call DispCarvalhol(RDF%kPoints)
        !call wLog("RDF%SkVec")
        !call DispCarvalhol(RDF%SkVec, unit_in = RDF%log_ID)

        kNLocal = size(RDF%kPoints,2)
        !kNStepLocal = RDF%kNStep
        !kNStepLocal(RDF%nDim) = kNLocal/sliceSize

        call wLog("MSH%origin")
        call wLog(MSH%origin)
        call wLog("RDF%kExtent")
        call wLog(RDF%kExtent)
        call wLog("kNLocal")
        call wLog(kNLocal)
        allocate(gammaK(kNLocal))
        allocate(phik(kNLocal))

        !Putting away the random numbers from others k (that are in others procs)
        do i_long = 1, kNCumulated_Init-1
            call random_number(trashNumber)
        end do
        call random_number(gammaK(:))
        do i_long = kNCumulated_End+1, product(RDF%kNStep)
            call random_number(trashNumber)
        end do
        do i_long = 1, kNCumulated_End-1
            call random_number(trashNumber)
        end do
        call random_number(phiK(:))

        call wLog("shape(gammaK)")
        call wLog(shape(gammaK))
        call wLog("shape(phiK)")
        call wLog(shape(phiK))
        call wLog("shape(RDF%SkVec)")
        call wLog(shape(RDF%SkVec))

        gammaK       = gammaK -0.5
        RDF%SkVec(:) =  gammak*sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);
        !RDF%SkVec(:) =  sqrt(RDF%SkVec)*cos(2.0D0*PI*phik);

        if(allocated(gammaK)) deallocate(gammaK)
        if(allocated(phik))   deallocate(phik)

        !cdata = fftw_alloc_real(alloc_local)
        !call c_f_pointer(cdata, data_real_2D, [L,local_M])

        ampMult = 2.0d0*sqrt(product(MSH%xStep)/((2.0d0)**(dble(RDF%nDim))))

        if(RDF%nDim == 2) then
            plan = fftw_mpi_plan_r2r(RDF%nDim, [M, L], data_real_2D, data_real_2D, &
                                     RDF%comm, [FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
            data_real_2D(:,:) = reshape(RDF%SkVec, [L, local_LastDim])
            call fftw_mpi_execute_r2r(plan, data_real_2D, data_real_2D)
            data_real_2D = data_real_2D*sqrt(product(MSH%xStep))
            !RDF%randField(:,1) = pack(data_real_2D(1:L,1:local_LastDim), .true.)
            RDF%randField(:,1) = reshape(data_real_2D, [L*local_LastDim])

        else if(RDF%nDim == 3) then
            call wLog("L = ")
            call wLog(L)
            call wLog("M = ")
            call wLog(L)
            call wLog("N = ")
            call wLog(L)
            call wLog("local_LastDim")
            call wLog(local_LastDim)
            call wLog("shape(data_real_3D) = ")
            call wLog(shape(data_real_3D))
            write(*,*) "RDF%comm = ", RDF%comm

            plan = fftw_mpi_plan_r2r(RDF%nDim, [N, M, L], data_real_3D, data_real_3D, &
                                     RDF%comm, [FFTW_REDFT01, FFTW_REDFT01, FFTW_REDFT01], FFTW_ESTIMATE)
            data_real_3D(:,:,:) = reshape(RDF%SkVec, [L, M, local_LastDim])
            write(*,*) "Calculating FFT In rang ", RDF%rang
            !write(*,*) "plan = ", plan
            write(*,*) "shape(data_real_3D) = ", shape(data_real_3D)
            call fftw_mpi_execute_r2r(plan, data_real_3D, data_real_3D)
            write(*,*) "AFTER Calculating FFT In rang ", RDF%rang
            data_real_3D = data_real_3D*(2.0D0)*sqrt(product(MSH%xStep))
            !data_real_3D = data_real_3D*(2.0D0**((RDF%nDim-1))/2.0D0))*sqrt(product(MSH%xStep))
            !RDF%randField(:,1) = pack(data_real_3D(1:L, 1:M, 1:local_LastDim), .true.)
            RDF%randField(:,1) = reshape(data_real_3D, [L*M*local_LastDim])
        end if

        call wLog("shape(RDF%SkVec)")
        call wLog(shape(RDF%SkVec))
        call wLog("shape(RDF%randField)")
        call wLog(shape(RDF%randField))


        call fftw_destroy_plan(plan)
        call fftw_free(cdata)

        call wLog("L")
        call wLog(L)
        call wLog("M")
        call wLog(M)
        call wLog("local_LastDim")
        call wLog(local_LastDim)
        call wLog("local_LD_offset")
        call wLog(local_LD_offset)
        call wLog("RDF%kNInit")
        call wLog(RDF%kNInit)
        call wLog("RDF%kNEnd")
        call wLog(RDF%kNEnd)

        if(allocated(gammaK)) deallocate(gammaK)
        if(allocated(phik)) deallocate(phik)

    end subroutine gen_Std_Gauss_FFT

    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    subroutine allocate_randField(RDF, xNstep, randField)
        implicit none
        !INPUT AND OUTPUT
        type(RF)   :: RDF
        double precision, dimension(:,:), allocatable, target :: randField
        !INPUT
        integer, dimension(:), intent(in) :: xNstep
        !LOCAL
        integer(kind=8) :: xNTotal

        xNTotal = product(int(xNStep,8))

        if(allocated(randField)) then
            if(.not.(size(randField,1) == xNTotal .and. size(randField,1) == RDF%Nmc)) then
                nullify(RDF%randField)
                deallocate(randField)
            end if
        end if

        if(.not.allocated(randField)) then
            allocate(randField(xNTotal, RDF%Nmc))
        end if

        call associate_randField(RDF, xNstep, randField)

    end subroutine allocate_randField

    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    subroutine associate_randField(RDF, xNstep, randField)
        implicit none
        !INPUT!
        double precision, dimension(:,:), intent(in), target :: randField
        integer, dimension(:), intent(in) :: xNstep
        !OUTPUT
        type(RF)   :: RDF

        RDF%randField => randField
        call wLog(" Inside allocate_randField")
        call wLog("     IN xNStep = ")
        call wLog(xNStep)
        if(RDF%nDim == 2) RDF%RF_2D(1:xNStep(1),1:xNStep(2)) => randField(:,1)
        if(RDF%nDim == 3) RDF%RF_3D(1:xNStep(1),1:xNStep(2),1:xNStep(3)) => randField(:,1)
    end subroutine associate_randField

    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------
    subroutine normalize_randField(randField, minPos, maxPos, &
                                   nDim, Nmc, loc_Comm, RF_2D, RF_3D)
        implicit none
        !INPUT OUTPUT
        double precision, dimension(:,:), intent(inout) :: randField
        !INPUT
        integer, dimension(:), intent(in) :: minPos, maxPos
        integer, intent(in) :: nDim, loc_Comm, Nmc !OBS Nmc, was setted manually in the caller
        double precision, dimension(:,:), intent(in), optional :: RF_2D
        double precision, dimension(:,:,:), intent(in), optional :: RF_3D
        !LOCAL
        double precision, dimension(Nmc) :: sumRF, sumRFsquare
        double precision, dimension(Nmc) :: totalSumRF, totalSumRFsquare;
        double precision, dimension(Nmc) :: evntAvg, evntStdDev;
        !integer, dimension(:), allocatable :: xNTotal_Vec, deplacement
        integer :: code
        integer :: i
        integer(kind=8) :: xNTotal, xNLocal

        if(nDim == 2 .and. (.not. present(RF_2D))) then
            write(*,*) "ERROR on normalize_randField nDim = 2 and RF_2D not present"
        else if(nDim == 3 .and. (.not. present(RF_3D))) then
            write(*,*) "ERROR on normalize_randField nDim = 3 and RF_3D not present"
        end if

        !Total Number of Points
        call wLog("Discovering Total Number of Points")
        xNLocal = product(int(maxPos - minPos + 1, 8))
        call MPI_ALLREDUCE (xNLocal, xNTotal, 1,MPI_INTEGER8, &
                            MPI_SUM,loc_Comm,code)
        call wLog("xNTotal = ")
        call wLog(xNTotal)


        call wLog("Calculating Average and stdVar")


        !AVERAGE---------------------------------------------------------
        if(nDim == 2) sumRF(:) = sum(RF_2D(minPos(1):maxPos(1), &
                                           minPos(2):maxPos(2)))
        if(nDim == 3) sumRF(:) = sum(RF_3D(minPos(1):maxPos(1), &
                                           minPos(2):maxPos(2), &
                                           minPos(3):maxPos(3)))

        !sumRF(:)       = sum( RDF%randField    , dim = 1)
        call wLog("sumRF(1) = ")
        call wLog(sumRF(1))
        call MPI_ALLREDUCE (sumRF,totalSumRF,Nmc,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,loc_Comm,code)
        evntAvg      = totalSumRF/dble(xNTotal);
        call wLog("Initial Average = ")
        call wLog(evntAvg)

        do i = 1, Nmc
            randField(:,i) = randField(:,i) - evntAvg(i)
        end do

        !Verifying Average
        if(nDim == 2) sumRF(:) = sum(RF_2D(minPos(1):maxPos(1), &
                                           minPos(2):maxPos(2)))
        if(nDim == 3) sumRF(:) = sum(RF_3D(minPos(1):maxPos(1), &
                                           minPos(2):maxPos(2), &
                                           minPos(3):maxPos(3)))
        !sumRF(:)       = sum( RDF%randField    , dim = 1)
        call MPI_ALLREDUCE (sumRF,totalSumRF,Nmc,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,loc_Comm,code)
        evntAvg      = totalSumRF/dble(xNTotal);
        call wLog("Final Average = ")
        call wLog(evntAvg)


        !STANDARD DEVIATION-----------------------------------------------------
        if(nDim == 2) sumRFsquare(:) = sum((RF_2D(minPos(1):maxPos(1), &
                                                  minPos(2):maxPos(2)) &
                                                  )**2.0D0)
        if(nDim == 3) sumRFsquare(:) = sum((RF_3D(minPos(1):maxPos(1), &
                                                  minPos(2):maxPos(2), &
                                                  minPos(3):maxPos(3)) &
                                                  )**2.0D0)
        !sumRFsquare(:) = sum((RDF%randField)**2, dim = 1)
        call wLog("sumRFsquare(1) = ")
        call wLog(sumRFsquare(1))
        call MPI_ALLREDUCE (sumRFsquare,totalSumRFsquare,Nmc,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,loc_Comm,code)
        evntStdDev   = sqrt(totalSumRFsquare/dble(xNTotal)) !Mean of random field is supposed to be 0
        call wLog("Initial StdDev = ")
        call wLog(evntStdDev)
        do i = 1, Nmc
            randField(:,i) = randField(:,i)/evntStdDev(i)
        end do

        !Verifying Standard Deviation
        if(nDim == 2) sumRFsquare(:) = sum((RF_2D(minPos(1):maxPos(1), &
                                                  minPos(2):maxPos(2)) &
                                                  )**2.0D0)
        if(nDim == 3) sumRFsquare(:) = sum((RF_3D(minPos(1):maxPos(1), &
                                                  minPos(2):maxPos(2), &
                                                  minPos(3):maxPos(3)) &
                                                  )**2.0D0)
        !sumRFsquare(:) = sum((RDF%randField)**2.0D0, dim = 1)
        call MPI_ALLREDUCE (sumRFsquare,totalSumRFsquare,Nmc,MPI_DOUBLE_PRECISION, &
                            MPI_SUM,loc_Comm,code)
        evntStdDev   = sqrt(totalSumRFsquare/dble(xNTotal)) !Mean of random field is supposed to be 0
        call wLog("Final StdDev = ")
        call wLog(evntStdDev)

    end subroutine normalize_randField


end module randomFieldND
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
