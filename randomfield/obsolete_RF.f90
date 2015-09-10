!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module obsolete_RF

    use displayCarvalhol
    use math_RF
    use mpi

    interface set_allRandField
        module procedure set_allRandFieldStructured,   &
            set_allRandFieldUnstruct
    end interface set_allRandField

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine get_sizes_MPI(xNStep, sizeLoc, sizeUnif, start, ends)

        implicit none
        !INPUT
        integer, dimension(:), intent(in) :: xNStep;
        !OUTPUT
        integer, intent(out), optional :: sizeLoc, sizeUnif;
        integer, intent(out), optional :: start, ends;

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
        if(present(ends))     ends     = xEnd

    end subroutine get_sizes_MPI
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
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
        integer :: xNTotalGlob, i;
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

        !Reordering
        if(rang == 0) call reorderToGlobal(all_RandField, all_xPoints, mapping)

        if(allocated(procDpl))             deallocate(procDpl)
        if(allocated(all_xNTotal))         deallocate(all_xNTotal)
        if(allocated(all_xPointsTransp))   deallocate(all_xPointsTransp)
        if(allocated(all_RandFieldTransp)) deallocate(all_RandFieldTransp)
        if(allocated(mapping))             deallocate(mapping)

    end subroutine set_allRandFieldUnstruct

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
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
        integer :: xNTotalGlob, i, j, globCount;
        double precision, dimension(:, :), allocatable :: all_RandFieldTransp;
        integer         , dimension(:)   , allocatable :: all_xNTotal, procDpl;
        integer         , dimension(:)   , allocatable :: mapping;
        !        double precision, dimension(:, :), allocatable :: xPoints; !Only for tests

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

        if(allocated(procDpl))             deallocate(procDpl)
        if(allocated(all_xNTotal))         deallocate(all_xNTotal)
        if(allocated(all_RandFieldTransp)) deallocate(all_RandFieldTransp)
        if(allocated(mapping))             deallocate(mapping)

    end subroutine set_allRandFieldStructured

end module obsolete_RF

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
