module type_STAT

    use mpi

    implicit none

    type :: STAT
        !MPI VARIABLES
        integer :: comm = -1, valid_comm = -1
        integer :: rang = -1
        integer :: nb_procs = -1, valid_nb_procs = -1
        integer :: log_ID = -1

        !GENERATION VARIABLES
            !nDim independent
        integer :: nDim = -1
        integer :: corrMod = -1 !1 for Gaussian
        integer :: margiFirst = -1 !1 for Gaussian, 2 for Lognormal
        integer :: method = -1 !1 for Isotropic, 2 for Shinozuka, 3 for Randomization, 4 for FFT
        integer :: Nmc = -1
        logical :: init = .false.
        logical :: valid_proc
        logical :: independent
            !nDim dependent
        double precision, dimension(:)   , allocatable :: corrL, corrL_out
        double precision, dimension(:, :), allocatable :: randField
        double precision, dimension(:), allocatable :: xMaxGlob, xMinGlob, xStep;
        double precision, dimension(:), allocatable :: overlap;
        integer(kind=8) , dimension(:), allocatable :: xNStep, xNStep_Loc!, kNStep
        integer,          dimension(:), allocatable :: procPerDim, coords
        integer(kind=8),  dimension(:,:), allocatable :: localRange
        integer(kind=8),  dimension(:,:), allocatable :: Sk_Ind, SkTot_Ind
        double precision, dimension(:), allocatable :: Sk_Dir, SkTot_Dir
        !integer(kind=8) :: xNTotal, kNTotal
        !integer(kind=8) :: sum_xNTotal, sum_kNTotal

        double precision, dimension(:), allocatable :: evntAvg, evntStdDev;
        double precision, dimension(:), allocatable :: pointAvg, pointStdDev;
        double precision :: globalAvg, globalStdDev;

    end type STAT

    contains
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine init_STAT(STAT_a, nDim, Nmc, method, corrMod, margiFirst, independent)

            type(STAT) :: STAT_a
            integer, intent(in) :: nDim, Nmc, method, corrMod, margiFirst
            logical, intent(in) :: independent

            !Obs: comm, rang and nb_procs sould be defined before
            STAT_a%nDim        = nDim
            STAT_a%Nmc         = Nmc
            STAT_a%method      = method
            STAT_a%corrMod     = corrMod
            STAT_a%margiFirst  = margiFirst
            STAT_a%independent = independent
            allocate(STAT_a%corrL(nDim))
            allocate(STAT_a%xMinGlob(nDim))
            allocate(STAT_a%xMaxGlob(nDim))
            allocate(STAT_a%xStep(nDim))
            allocate(STAT_a%overlap(nDim))
            allocate(STAT_a%xNStep(nDim))
            allocate(STAT_a%xNStep_Loc(nDim))
            allocate(STAT_a%procPerDim(nDim))
            allocate(STAT_a%coords(nDim))
            allocate(STAT_a%localRange(nDim,2))
            allocate(STAT_a%Sk_Ind(nDim,2))
            allocate(STAT_a%SkTot_Ind(nDim,2))
            allocate(STAT_a%corrL_out(nDim))
            allocate(STAT_a%evntAvg(Nmc))
            allocate(STAT_a%evntStdDev(Nmc))
            STAT_a%corrL  = -1
            STAT_a%xMinGlob = -1
            STAT_a%xMaxGlob = -1
            STAT_a%xStep = -1
            STAT_a%overlap = -1
            STAT_a%xNStep = -1
            STAT_a%globalAvg = -1
            STAT_a%globalStdDev = -1
            STAT_a%procPerDim = -1
            STAT_a%localRange = -1
            STAT_a%coords = -1
            STAT_a%corrL_out = -1
            STAT_a%init  = .true.

        end subroutine init_STAT

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine show_STAT(STAT_a, name, unit)

            type(STAT), intent(in) :: STAT_a
            character(len=*) :: name
            integer :: unit

            !Obs: comm, rang and nb_procs sould be defined before
            write(unit,*) " ---------------------------- "
            write(unit,*) "  NAME = ", name
            write(unit,*) "  "
            write(unit,*) "  comm           = ", STAT_a%comm
            write(unit,*) "  rang           = ", STAT_a%rang
            write(unit,*) "  nb_procs       =  ", STAT_a%nb_procs
            write(unit,*) "  nDim           = ", STAT_a%nDim
            write(unit,*) "  Nmc            = ", STAT_a%Nmc
            write(unit,*) "  method         = ", STAT_a%method
            write(unit,*) "  corrMod        = ", STAT_a%corrMod
            write(unit,*) "  margiFirst     = ", STAT_a%margiFirst
            write(unit,*) "  independent    = ", STAT_a%independent
            write(unit,*) "  corrL (IN)     = ", STAT_a%corrL
            write(unit,*) "  xMinGlob       = ", STAT_a%xMinGlob
            write(unit,*) "  xMaxGlob       = ", STAT_a%xMaxGlob
            write(unit,*) "  xStep          = ", STAT_a%xStep
            write(unit,*) "  xNStep         = ", STAT_a%xNStep
            write(unit,*) "  xNStep_Loc     = ", STAT_a%xNStep_Loc
            write(unit,*) "  overlap        = ", STAT_a%overlap
            write(unit,*) "  init           = ", STAT_a%init
            write(unit,*) "  allocated RF   = ", allocated(STAT_a%randField)
            write(unit,*) "  shape(RF)      = ", shape(STAT_a%randField)
            write(unit,*) "  procPerDim     = ", STAT_a%procPerDim
            write(unit,*) "  localRange INF = ", STAT_a%localRange(:,1)
            write(unit,*) "  localRange SUP = ", STAT_a%localRange(:,2)
            write(unit,*) "  allocated(Sk_Dir) = ", allocated(STAT_a%Sk_Dir)
            if(allocated(STAT_a%Sk_Dir)) write(unit,*) "  shape(Sk_Dir)      = ", shape(STAT_a%Sk_dir)
            write(unit,*) "  Sk_Ind INF     = ", STAT_a%Sk_Ind(:,1)
            write(unit,*) "  Sk_Ind SUP     = ", STAT_a%Sk_Ind(:,2)
            write(unit,*) "  coords         = ", STAT_a%coords
            write(unit,*) "  STATISTICS----------------- "
            if(allocated(STAT_a%evntAvg)) then
                write(unit,*) "  evntAvg      = ", STAT_a%evntAvg
                write(unit,*) "  evntStdDev   = ", STAT_a%evntStdDev
            end if
            if(allocated(STAT_a%pointAvg)) then
                write(unit,*) "  pointAvg     = ", STAT_a%pointAvg(1:3)
                write(unit,*) "  pointStdDev  = ", STAT_a%pointStdDev(1:3)
            end if
            write(unit,*) "  globalAvg    = ", STAT_a%globalAvg
            write(unit,*) "  globalStdDev = ", STAT_a%globalStdDev
            write(unit,*) "  corrL (OUT)  = ", STAT_a%corrL_out
            write(unit,*) "  "

        end subroutine show_STAT

        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        !---------------------------------------------------------------------------------
        subroutine finalize_STAT(STAT_a)
            type(STAT) :: STAT_a

            STAT_a%nDim     = -1
            STAT_a%Nmc      = -1
            STAT_a%comm     = -1
            STAT_a%rang     = -1
            STAT_a%nb_procs = -1
            if(allocated(STAT_a%corrL))       deallocate(STAT_a%corrL)
            if(allocated(STAT_a%xMinGlob))    deallocate(STAT_a%xMinGlob)
            if(allocated(STAT_a%xMaxGlob))    deallocate(STAT_a%xMaxGlob)
            if(allocated(STAT_a%xStep))       deallocate(STAT_a%xStep)
            if(allocated(STAT_a%overlap))     deallocate(STAT_a%overlap)
            if(allocated(STAT_a%randField))   deallocate(STAT_a%randField)
            if(allocated(STAT_a%xNStep))      deallocate(STAT_a%xNStep)
            if(allocated(STAT_a%evntAvg))     deallocate(STAT_a%evntAvg)
            if(allocated(STAT_a%evntStdDev))  deallocate(STAT_a%evntStdDev)
            if(allocated(STAT_a%pointAvg))    deallocate(STAT_a%pointAvg)
            if(allocated(STAT_a%pointStdDev)) deallocate(STAT_a%pointStdDev)
            if(allocated(STAT_a%procPerDim))  deallocate(STAT_a%procPerDim)
            if(allocated(STAT_a%localRange))  deallocate(STAT_a%localRange)
            if(allocated(STAT_a%coords))      deallocate(STAT_a%coords)
            if(allocated(STAT_a%Sk_Ind))      deallocate(STAT_a%Sk_Ind)
            if(allocated(STAT_a%Sk_Dir))      deallocate(STAT_a%Sk_Dir)
            if(allocated(STAT_a%xNStep_Loc))  deallocate(STAT_a%xNStep_Loc)
            if(allocated(STAT_a%SkTot_Dir))   deallocate(STAT_a%SkTot_Dir)
            if(allocated(STAT_a%SkTot_Ind))   deallocate(STAT_a%SkTot_Ind)
            if(allocated(STAT_a%corrL_out))   deallocate(STAT_a%corrL_out)
            if(allocated(STAT_a%evntAvg))     deallocate(STAT_a%evntAvg)
            if(allocated(STAT_a%evntStdDev))  deallocate(STAT_a%evntStdDev)

            STAT_a%init = .false.

        end subroutine finalize_STAT

end module type_STAT
