module mesh_RF

    !use mpi
    use math_RF
    use write_Log_File
    use type_RF
    use type_MESH
    use fftw3

    implicit none

contains
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine set_XPoints(MSH, RDF, xPoints)
        implicit none

        !INPUT AND OUTPUT
        type(MESH) :: MSH
        type(RF)   :: RDF
        double precision, dimension(:, :), allocatable, intent(out), target :: xPoints;

        !LOCAL

        !write(*,*) "Allocating xPoints"
        allocate(xPoints(MSH%nDim, MSH%xNTotal))

        call wLog("MSH%xMinBound = ")
        call wLog(MSH%xMinBound)
        call wLog("MSH%xStep = ")
        call wLog(MSH%xStep)
        call wLog("MSH%xNStep = ")
        call wLog(MSH%xNStep)
        call wLog("MSH%xNTotal = ")
        call wLog(MSH%xNTotal)

        !write(*,*) "Setting Grid"
        !write(*,*) "shape(xPoints) = ", shape(xPoints)
        !write(*,*) "MSH%xStep      = ", MSH%xStep
        !write(*,*) "MSH%xMinBound  = ", MSH%xMinBound
        !write(*,*) "MSH%xNStep = ", MSH%xNStep
        call setGrid(xPoints, MSH%xMinBound, MSH%xStep, MSH%xNStep)
        !write(*,*) "AFTER Setting Grid"

        !write(*,*) "Pointing xPoints"
        RDF%xPoints => xPoints

        if(RDF%nDim == 2) then
            RDF%xPoints_2D(1:MSH%nDim, 1:MSH%xNStep(1), 1:MSH%xNStep(2)) => xPoints
            !call DispCarvalhol(RDF%xPoints_2D(1,:,:), "RDF%xPoints_2D 1", unit_in=RDF%log_ID)
            !call DispCarvalhol(RDF%xPoints_2D(2,:,:), "RDF%xPoints_2D 2", unit_in=RDF%log_ID)
        else if(RDF%nDim == 3) then
            RDF%xPoints_3D(1:MSH%nDim, 1:MSH%xNStep(1), 1:MSH%xNStep(2), 1:MSH%xNStep(3)) => xPoints
            !call wLog("Point in minimal position = ")
            !call wLog(RDF%xPoints_3D(:,22,1,1))
        end if

        !write(*,*) "AFTER Pointing xPoints"

        !call DispCarvalhol(transpose(RDF%xPoints), "transpose(RDF%xPoints)", unit_in=RDF%log_ID)

    end subroutine set_XPoints


end module mesh_RF
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
