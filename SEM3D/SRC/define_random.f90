module define_random

    use sdomain
    use mpi
    use constants
    use randomFieldND
    use displayCarvalhol
    use writeResultFile_RF

contains
    !---------------------------------------------------------------------------
    subroutine define_random_subdomains(Tdomain, rg)
        implicit none
        !INPUT
        type(domain), intent(inout), target :: Tdomain
        integer     , intent(in)  :: rg

        !LOCAL
        integer :: mat, assocMat
        integer :: fileId, code


        !TEST
        integer :: rgSubD, error
        integer(kind= MPI_OFFSET_KIND ) :: filePos
        integer, dimension( MPI_STATUS_SIZE ) :: status
        integer :: nbBitesInt, sizeRand
        integer :: i
        integer, parameter :: nb_valeurs=10
        integer, dimension(nb_valeurs) :: valeurs

        !Non PMLs
        !        call MPI_FILE_OPEN ( Tdomain%communicateur ,"randSeed.txt", MPI_MODE_RDWR + MPI_MODE_CREATE , &
        !                                     MPI_INFO_NULL ,fileId, code)
        call MPI_FILE_OPEN ( Tdomain%communicateur ,"./prop/randSeed.txt", MPI_MODE_RDWR + MPI_MODE_CREATE , &
            MPI_INFO_NULL ,fileId, code)
        do mat = 0, Tdomain%n_mat - 1
            if(Tdomain%sSubDomain(mat)%material_type == "R" .and. &
                Tdomain%subD_exist(mat)) then
                call build_mat_extremes(Tdomain, rg, mat)
                call define_random_seed(Tdomain, rg, mat, fileId)
            end if
        end do
        call MPI_FILE_CLOSE (fileId,code)

        !        !START TEST
        !        call MPI_FILE_OPEN ( Tdomain%communicateur ,"./prop/randSeed.txt", MPI_MODE_RDWR + MPI_MODE_CREATE , &
        !                                     MPI_INFO_NULL ,fileId, code)
        !        Tdomain%logicD%run_restart = .true.
        !        do mat = 0, Tdomain%n_mat - 1
        !            if(Tdomain%sSubDomain(mat)%material_type == "R" .and. &
        !                  Tdomain%subD_exist(mat)) then
        !                   write(*,*) "1) proc ", rg, " mat ", mat, "%chosenSeed(:) = ", Tdomain%sSubdomain(mat)%chosenSeed(:)
        !                   !Tdomain%sSubdomain(mat)%chosenSeed(:) = -1
        !                   !write(*,*) "2) proc ", rg, " mat ", mat, "%chosenSeed(:) = ", Tdomain%sSubdomain(mat)%chosenSeed(:)
        !                   call define_random_seed(Tdomain, rg, mat, fileId)
        !                   write(*,*) "3) proc ", rg, " mat ", mat, "%chosenSeed(:) = ", Tdomain%sSubdomain(mat)%chosenSeed(:)
        !            end if
        !        end do
        !        call MPI_FILE_CLOSE (fileId,code)
        !        Tdomain%logicD%run_restart = .false.
        !        !END TEST

        !PMLS
        do mat = 0, Tdomain%n_mat - 1
            assocMat = Tdomain%sSubDomain(mat)%assocMat
            if(.not. (Tdomain%not_PML_List(mat))                        .and. &
                Tdomain%sSubdomain(assocMat)%material_type == "R" .and. &
                Tdomain%subD_exist(mat)) then !PMLs associated to random subdomains

                Tdomain%sSubdomain(mat)%chosenSeed(:) = Tdomain%sSubdomain(assocMat)%chosenSeed(:)
                Tdomain%sSubdomain(mat)%MinBound(:)   = Tdomain%sSubdomain(assocMat)%MinBound(:)
                Tdomain%sSubdomain(mat)%MaxBound(:)   = Tdomain%sSubdomain(assocMat)%MaxBound(:)
            end if
        end do

    end subroutine define_random_subdomains

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine build_mat_extremes(Tdomain, rg, mat)
        !WARNING: *The subdomains communicators should've been created before this subroutine
        implicit none
        !INPUT
        type(domain), intent(inout), target :: Tdomain
        integer, intent(in)  :: rg, mat

        !LOCAL
        integer              :: ngllx, nglly, ngllz, ipoint
        integer              :: m, n, i, j, k, coord
        integer              :: code, error
        real, dimension(0:2) :: tempMin, tempMax

        !Verification
        if(.not.(Tdomain%sSubDomain(mat)%material_type == "R" .and. Tdomain%subD_exist(mat))) then
            write(*,*) "!!!WARNING:'build_mat_extremes' was for a non-random material_type"
            write(*,*) "Rang =", rg, "mat = ", mat
            write(*,*) "material_type =", Tdomain%sSubDomain(mat)%material_type, ", existence = ", Tdomain%subD_exist(mat)
            !call MPI_ABORT(Tdomain%communicateur, error, code)
        end if

        !Initialize extremes
        if(.not.allocated(Tdomain%sSubDomain(mat)%MinBound)) allocate(Tdomain%sSubDomain(mat)%MinBound(0:2))
        if(.not.allocated(Tdomain%sSubDomain(mat)%MaxBound)) allocate(Tdomain%sSubDomain(mat)%MaxBound(0:2))
        n = Tdomain%sSubDomain(mat)%elemList(0)
        ipoint = Tdomain%specel(n)%Iglobnum(0,0,0)
        Tdomain%sSubDomain(mat)%MinBound = [Tdomain%GlobCoord(0,ipoint), &
            Tdomain%GlobCoord(1,ipoint), &
            Tdomain%GlobCoord(2,ipoint)]
        Tdomain%sSubDomain(mat)%MaxBound = Tdomain%sSubDomain(mat)%MinBound

        !Calculate local extremes
        ngllx = Tdomain%sSubDomain(mat)%NGLLx
        nglly = Tdomain%sSubDomain(mat)%NGLLy
        ngllz = Tdomain%sSubDomain(mat)%NGLLz

        do m = 0, Tdomain%sSubDomain(mat)%nElem - 1
            n = Tdomain%sSubDomain(mat)%elemList(m)
            do i = 0, ngllx-1
                do j = 0, nglly-1
                    do k = 0, ngllz-1
                        ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)
                        !write(*,*) "Is R material and ipoint = ", ipoint
                        do coord = 0, 2
                            if(Tdomain%GlobCoord(coord,ipoint) < Tdomain%sSubDomain(mat)%MinBound(coord)) &
                                Tdomain%sSubDomain(mat)%MinBound(coord) = Tdomain%GlobCoord(coord,ipoint)
                            if(Tdomain%GlobCoord(coord,ipoint) > Tdomain%sSubDomain(mat)%MaxBound(coord)) &
                                Tdomain%sSubDomain(mat)%MaxBound(coord) = Tdomain%GlobCoord(coord,ipoint)
                        end do
                    end do
                end do
            end do !END Loop over GLLs
        end do !END Loop over subdomain elements

        !Establishing the global extremes
        tempMin = Tdomain%sSubDomain(mat)%MinBound
        tempMax = Tdomain%sSubDomain(mat)%MaxBound

        call MPI_ALLREDUCE(tempMin,                                &
            Tdomain%sSubDomain(mat)%MinBound,       &
            size(Tdomain%sSubDomain(mat)%MinBound), &
            MPI_DOUBLE_PRECISION, MPI_MIN,          &
            Tdomain%subDComm(mat) ,code)
        call MPI_ALLREDUCE(tempMax,                                &
            Tdomain%sSubDomain(mat)%MaxBound,       &
            size(Tdomain%sSubDomain(mat)%MaxBound), &
            MPI_DOUBLE_PRECISION, MPI_MAX,          &
            Tdomain%subDComm(mat) ,code)
    end subroutine build_mat_extremes

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine define_random_seed(Tdomain, rg, mat, fileId)
        !WARNING: *The subdomains communicators should've been created before this subroutine
        !         *"allocate_random_subdomain" should be called before this subroutine
        implicit none
        !INPUT
        type(domain), intent(inout), target :: Tdomain
        integer     , intent(in)  :: rg, mat, fileId
        !LOCAL
        integer :: rgSubD, code, error
        integer(kind= MPI_OFFSET_KIND ) :: filePos
        integer, dimension( MPI_STATUS_SIZE ) :: status
        integer :: nbBitesInt, sizeRand

        if(.not.(Tdomain%sSubDomain(mat)%material_type == "R" .and. Tdomain%subD_exist(mat))) then
            write(*,*) "!!!ERROR:'define_random_seed' was called wrongly"
            write(*,*) "Rang =", rg, "mat = ", mat
            write(*,*) "material_type =", Tdomain%sSubDomain(mat)%material_type, ", existence = ", Tdomain%subD_exist(mat)
            call MPI_ABORT(Tdomain%communicateur, error, code)
        end if

        !Choosing the seed for properties field creation
        rgSubD  = -1
        sizeRand = size(Tdomain%sSubdomain(mat)%chosenSeed)
        call MPI_COMM_RANK  (Tdomain%subDComm(mat) ,rgSubD, code)

        if(Tdomain%logicD%run_restart) then
            write(*,*) "This is a reprise"
            if(rgSubD == 0) then
                write(*,*) "And rgSubD is in proc ", rg
                call MPI_TYPE_SIZE (MPI_INTEGER, nbBitesInt, code)
                filePos = mat*nbBitesInt*sizeRand
                call MPI_FILE_READ_AT (fileId, filePos,                                           &
                    Tdomain%sSubdomain(mat)%chosenSeed,sizeRand, MPI_INTEGER , &
                    status,code)
            end if
        else !This is not reprise
            if (rgSubD == 0) then
                if(Tdomain%sSubdomain(mat)%seedStart >= 0) then
                    call calculate_random_seed(Tdomain%sSubdomain(mat)%chosenSeed, Tdomain%sSubdomain(mat)%seedStart)
                else
                    call calculate_random_seed(Tdomain%sSubdomain(mat)%chosenSeed)
                end if
                call MPI_TYPE_SIZE (MPI_INTEGER, nbBitesInt, code)
                filePos = mat*nbBitesInt*sizeRand
                call MPI_FILE_WRITE_AT (fileId, filePos,                                           &
                    Tdomain%sSubdomain(mat)%chosenSeed,sizeRand, MPI_INTEGER , &
                    status,code)
                !write(*,*) "mat        = ", mat
                !write(*,*) "nbBitesInt = ", nbBitesInt
                !write(*,*) "sizeRand   = ", sizeRand
                !write(*,*) "mat*nbBitesInt*sizeRand = ", mat*nbBitesInt*sizeRand
                !write(*,*) "...%chosenSeed          = ", Tdomain%sSubdomain(mat)%chosenSeed
            end if
        end if

        call MPI_BCAST (Tdomain%sSubdomain(mat)%chosenSeed,             &
            size(Tdomain%sSubdomain(mat)%chosenSeed),       &
            MPI_INTEGER, 0, Tdomain%subDComm(mat), code)
        !write(*,*)"..%chosenSeed   = ", Tdomain%sSubdomain(mat)%chosenSeed
    end subroutine define_random_seed

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine build_random_properties(Tdomain, rg, mat, xPoints, prop, method)
        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: rg, mat
        real            , intent(in),   dimension(0:, 0:) :: xPoints;
        integer         , intent(in),   optional          :: method
        !OUTPUT
        real        , intent(out), dimension(0:, 0:) :: prop !Properties
        !LOCAL
        integer :: ngllx,nglly,ngllz, assocMat
        integer :: iPoint, RFpoint
        integer :: i, j, k, m, n !counters
        integer :: error, code
        integer :: LimPML1, LimPML2, LimPML3
        integer :: nProp = 3
        integer :: effecMethod
        logical, dimension(:), allocatable :: calculate
        logical :: face_X, face_Y, face_Z, edge_XY, edge_YZ, edge_ZX, vertex_XYZ
        real, dimension(0:2) :: pointProp, avgProp;
        double precision, dimension(:, :), allocatable :: xPointsPML;
        real            , dimension(:, :), allocatable :: propPML !Properties
        logical         , dimension (:,:), allocatable :: PML_mask

        allocate(calculate(0:nProp-1))
        calculate(:) = .true.
        do i = 0, nProp - 1
            if(Tdomain%sSubDomain(mat)%varProp(i) <= 0) calculate(i) = .false.
        end do

        effecMethod = 1;
        if (present(method)) then
            if((method > 0) .and. (method<3)) then
                effecMethod = method
            else
                effecMethod = 1
                write(*,*) "WARNING! The chosen method is not an avaiable choice - method = ", method
                write(*,*) "         The method was automatically setted - method = ", effecMethod
            end if
        end if

        avgProp = [Tdomain%sSubDomain(mat)%Ddensity, &
            Tdomain%sSubDomain(mat)%DLambda,  &
            Tdomain%sSubDomain(mat)%DMu]
        assocMat = Tdomain%sSubdomain(mat)%assocMat
        ngllx    = Tdomain%sSubDomain(mat)%NGLLx
        nglly    = Tdomain%sSubDomain(mat)%NGLLy
        ngllz    = Tdomain%sSubDomain(mat)%NGLLz

        !            if(rg == 0) write(*,*) ">>>>Creating Standard Gaussian Field"
        !            if(rg == 0) write(*,*) "corrL                = ", Tdomain%sSubDomain(mat)%corrL
        !            if(rg == 0) write(*,*) "corrMod              = ", Tdomain%sSubDomain(mat)%corrMod
        !            if(rg == 0) write(*,*) "Nmc (nProp)          = ", nProp
        !            if(rg == 0) write(*,*) "chosenSeed           = ", Tdomain%sSubDomain(mat)%chosenSeed
        !            if(rg == 0) write(*,*) "MinBound             = ", Tdomain%sSubDomain(mat)%MinBound
        !            if(rg == 0) write(*,*) "MaxBound             = ", Tdomain%sSubDomain(mat)%MaxBound
        !            !if(rg == 0) write(*,*) "allocated(xPoints)   = ", allocated(xPoints)
        !            !if(rg == 0) write(*,*) "allocated(prop)      = ", allocated(prop)

        if(Tdomain%not_PML_List(mat)) then
            !write(*,*) "Material ", mat, " is not a PML"
            !write(*,*) "prop BEFORE RANDOM FIELD CREATION = ", prop(:,0)

            select case(effecMethod)
            case( 1 ) !Victor
                write(*,*) "Victor's method"
                call createStandardGaussianFieldUnstructVictor(&
                    xPoints(:, :),                            &
                    Tdomain%sSubDomain(mat)%corrL,            &
                    Tdomain%sSubDomain(mat)%corrMod,          &
                    nProp,                                    &
                    prop(:, :),                               &
                    Tdomain%sSubDomain(mat)%chosenSeed,       &
                    Tdomain%sSubDomain(mat)%MinBound,         &
                    Tdomain%sSubDomain(mat)%MaxBound,         &
                    Tdomain%communicateur,                    &
                    calculate)

            case( 2 ) !Shinozuka
                write(*,*) "Shinozuka's method"
                call createStandardGaussianFieldUnstructShinozuka (&
                    xPoints(:, :),                                 &
                    Tdomain%sSubDomain(mat)%corrL,                   &
                    Tdomain%sSubDomain(mat)%corrMod,                 &
                    nProp,                                           &
                    prop(:, :),                                      &
                    Tdomain%sSubDomain(mat)%chosenSeed,              &
                    Tdomain%sSubDomain(mat)%MinBound,                &
                    Tdomain%sSubDomain(mat)%MaxBound,                &
                    Tdomain%communicateur,                           &
                    calculate)

            case default
                write(*,*) "ERROR! The chosen method is not an avaiable choice"
                call MPI_ABORT(Tdomain%communicateur, error, code)
            end select

            !call dispCarvalhol(prop, "prop GAUSS", "F30.10")
            !////////Transfoming Stantard Gaussian Field
            if(rg == 0) write(*,*) ">>>>Transforming Standard Gaussian Field"

            i = 0
            if(rg == 0) write(*,*) "Dens------------ "
            if(rg == 0) write(*,*) "margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
            if(rg == 0) write(*,*) "average      = ", avgProp(i)
            if(rg == 0) write(*,*) "variance     = ", Tdomain%sSubDomain(mat)%varProp(i)
            i = 1
            if(rg == 0) write(*,*) "Lambda----------- "
            if(rg == 0) write(*,*) "margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
            if(rg == 0) write(*,*) "average      = ", avgProp(i)
            if(rg == 0) write(*,*) "variance     = ", Tdomain%sSubDomain(mat)%varProp(i)
            i = 2
            if(rg == 0) write(*,*) "Mu--------------- "
            if(rg == 0) write(*,*) "margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
            if(rg == 0) write(*,*) "average      = ", avgProp(i)
            if(rg == 0) write(*,*) "variance     = ", Tdomain%sSubDomain(mat)%varProp(i)

            do i = 0, nProp - 1
                call multiVariateTransformation (          &
                    Tdomain%sSubDomain(mat)%margiFirst(i), &
                    avgProp(i),                            &
                    Tdomain%sSubDomain(mat)%varProp(i),    &
                    prop(:, i:i))
            end do

            !call dispCarvalhol(prop, "prop AFTER", "F30.10")


        else !Random PML
            allocate(PML_Mask(0:size(Tdomain%GlobCoord,1)-1, 0:size(Tdomain%GlobCoord,2)-1))
            face_X     = .false.
            face_Y     = .false.
            face_Z     = .false.
            edge_XY    = .false.
            edge_YZ    = .false.
            edge_ZX    = .false.
            vertex_XYZ = .false.
            PML_Mask(:,:) = .false.

            !write(*,*) ">>>>Defining PML orientation"
            !/////////////Defining PML orientation
            !Face X oriented
            if  (        Tdomain%sSubDomain(mat)%Px   .and. &
                (.not.Tdomain%sSubDomain(mat)%Py)  .and. &
                (.not.Tdomain%sSubDomain(mat)%Pz)) then
                face_X = .true.
                !write(*,*) "face_X"
                !Face Y oriented
            elseif ((.not.Tdomain%sSubDomain(mat)%Px)  .and. &
                Tdomain%sSubDomain(mat)%Py   .and. &
                (.not.Tdomain%sSubDomain(mat)%Pz)) then
                face_Y = .true.
                !write(*,*) "face_Y"
                !Face Z oriented
            elseif ((.not.Tdomain%sSubDomain(mat)%Px) .and. &
                (.not.Tdomain%sSubDomain(mat)%Py) .and. &
                Tdomain%sSubDomain(mat)%Pz) then
                face_Z = .true.
                !write(*,*) "face_Z"

                !Edge in XY
            elseif (     Tdomain%sSubDomain(mat)%Px  .and. &
                (     Tdomain%sSubDomain(mat)%Py) .and. &
                (.not.Tdomain%sSubDomain(mat)%Pz)) then
                edge_XY = .true.
                !write(*,*) "edge_XY"

                !Edge in YZ
            elseif ((.not.Tdomain%sSubDomain(mat)%Px) .and. &
                (      Tdomain%sSubDomain(mat)%Py) .and. &
                (      Tdomain%sSubDomain(mat)%Pz)) then
                edge_YZ = .true.
                !write(*,*) "edge_YZ"

                !Edge in ZX
            elseif (      Tdomain%sSubDomain(mat)%Px   .and. &
                ((.not.Tdomain%sSubDomain(mat)%Py)) .and. &
                (      Tdomain%sSubDomain(mat)%Pz)) then
                edge_ZX = .true.
                !write(*,*) "edge_ZX"

                !Vertex in XYZ
            elseif (  Tdomain%sSubDomain(mat)%Px   .and. &
                (  Tdomain%sSubDomain(mat)%Py)  .and. &
                (  Tdomain%sSubDomain(mat)%Pz)) then
                vertex_XYZ = .true.
                !write(*,*) "vertex_XYZ"
                !Undefined PML
            else
                write(*,*) "ERROR in mat ", mat, " (PML) definition (directions), check 'material.input'"
                write(*,*) "Tdomain%sSubDomain(", mat, ")%Px = ", Tdomain%sSubDomain(mat)%Px
                write(*,*) "Tdomain%sSubDomain(", mat, ")%Py = ", Tdomain%sSubDomain(mat)%Py
                write(*,*) "Tdomain%sSubDomain(", mat, ")%Pz = ", Tdomain%sSubDomain(mat)%Pz
                call MPI_ABORT(Tdomain%communicateur, error, code)
            end if

            !/////////////Creating PML Mask
            !On the lower bound initialization
            LimPML1         = 0
            LimPML2         = 0
            LimPML3         = 0
            !To the upper bound transformation
            if (Tdomain%sSubDomain(mat)%Left)    LimPML1 = ngllx-1
            if (Tdomain%sSubDomain(mat)%Forward) LimPML2 = nglly-1
            if (Tdomain%sSubDomain(mat)%Down)    LimPML3 = ngllz-1

            !                    if (Tdomain%sSubDomain(mat)%Left) then
            !                        write(*,*) "Left"
            !                    else
            !                        write(*,*) "Right"
            !                    end if
            !                    if (Tdomain%sSubDomain(mat)%Forward) then
            !                        write(*,*) "Forward"
            !                    else
            !                        write(*,*) "Backward"
            !                    end if
            !                    if (Tdomain%sSubDomain(mat)%Down) then
            !                        write(*,*) "Down"
            !                    else
            !                        write(*,*) "Up"
            !                    end if

            do m = 0, Tdomain%sSubDomain(mat)%nElem - 1
                n = Tdomain%sSubDomain(mat)%elemList(m)
                !Face X oriented
                if  (face_X) then
                    do j = 0, nglly-1
                        do k = 0, ngllz-1
                            !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                            ipoint             = Tdomain%specel(n)%Iglobnum(LimPML1,j,k)
                            RFpoint            = count(Tdomain%sSubDomain(mat)%globCoordMask(0,0:ipoint)) - 1
                            PML_Mask(:,ipoint) = .true.
                        enddo
                    enddo
                    !Face Y oriented
                elseif (face_Y) then
                    do i = 0, ngllx-1
                        do k = 0, ngllz-1
                            !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                            ipoint             = Tdomain%specel(n)%Iglobnum(i, LimPML2, k)
                            RFpoint            = count(Tdomain%sSubDomain(assocMat)%globCoordMask(0,0:ipoint)) - 1
                            PML_Mask(:,ipoint) = .true.
                        enddo
                    enddo
                    !Face Z oriented
                elseif (face_Z) then
                    do i = 0, ngllx-1
                        do j = 0, nglly-1
                            !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                            ipoint             = Tdomain%specel(n)%Iglobnum(i,j,LimPML3)
                            RFpoint            = count(Tdomain%sSubDomain(assocMat)%globCoordMask(0,0:ipoint)) - 1
                            PML_Mask(:,ipoint) = .true.
                        enddo
                    enddo
                    !Edge in XY
                elseif (edge_XY) then
                    do k = 0, ngllz-1
                        !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                        ipoint             = Tdomain%specel(n)%Iglobnum(LimPML1,LimPML2,k)
                        RFpoint            = count(Tdomain%sSubDomain(assocMat)%globCoordMask(0,0:ipoint)) - 1
                        PML_Mask(:,ipoint) = .true.
                    enddo
                    !Edge in YZ
                elseif (edge_YZ) then
                    do i = 0, ngllx-1
                        !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                        ipoint             = Tdomain%specel(n)%Iglobnum(i, LimPML2,LimPML3)
                        RFpoint            = count(Tdomain%sSubDomain(assocMat)%globCoordMask(0,0:ipoint)) - 1
                        PML_Mask(:,ipoint) = .true.
                    enddo
                    !Edge in ZX
                elseif (edge_ZX) then
                    do j = 0, nglly-1
                        ipoint          = Tdomain%specel(n)%Iglobnum(LimPML1, j,LimPML3)
                        RFpoint         = count(Tdomain%sSubDomain(assocMat)%globCoordMask(0,0:ipoint)) - 1
                        PML_Mask(:,ipoint) = .true.
                    enddo
                    !Vertex
                elseif (vertex_XYZ) then
                    ipoint          = Tdomain%specel(n)%Iglobnum(LimPML1, LimPML2,LimPML3)
                    RFpoint         = count(Tdomain%sSubDomain(assocMat)%globCoordMask(0,0:ipoint)) - 1
                    PML_Mask(:,ipoint) = .true.
                end if
                !Loop over subdomain elements
            end do

            if(rg == 0) write(*,*) ">>>>Creating Stantard Random Field (PML)"
            !Creating PML properties
            allocate(xPointsPML (0:size(Tdomain%GlobCoord, 1)-1, 0:count(PML_Mask(0,:))-1))
            allocate(propPML    (0:count(PML_Mask(0,:))-1, 0:nProp-1))
            propPML    = -1
            xPointsPML = (reshape(pack(Tdomain%GlobCoord(:,:), mask = PML_Mask(:,:)), &
                shape = [3, count(PML_Mask(0,:))]))

            !write(*,*) "propPML BEFORE RANDOM FIELD CREATION = ", propPML(:,0)
            select case(effecMethod)
            case( 1 ) !Victor
                write(*,*) "Victor's method"
                call createStandardGaussianFieldUnstructVictor(&
                    xPointsPML(:, :),                         &
                    Tdomain%sSubDomain(mat)%corrL,            &
                    Tdomain%sSubDomain(mat)%corrMod,          &
                    nProp,                                    &
                    propPML(:, :),                            &
                    Tdomain%sSubDomain(mat)%chosenSeed,       &
                    Tdomain%sSubDomain(mat)%MinBound,         &
                    Tdomain%sSubDomain(mat)%MaxBound,         &
                    Tdomain%communicateur,                    &
                    calculate)

            case( 2 ) !Shinozuka
                write(*,*) "Shinozuka's method"
                call createStandardGaussianFieldUnstructShinozuka(&
                    xPointsPML(:, :),                             &
                    Tdomain%sSubDomain(mat)%corrL,                &
                    Tdomain%sSubDomain(mat)%corrMod,              &
                    nProp,                                        &
                    propPML(:, :),                                &
                    Tdomain%sSubDomain(mat)%chosenSeed,           &
                    Tdomain%sSubDomain(mat)%MinBound,             &
                    Tdomain%sSubDomain(mat)%MaxBound,             &
                    Tdomain%communicateur,                           &
                    calculate)

            case default
                write(*,*) "ERROR! The chosen method is not an avaiable choice"
                call MPI_ABORT(Tdomain%communicateur, error, code)
            end select
            deallocate(xPointsPML)
            !write(*,*) "propPML BEFORE TRANSFORMATION = ", propPML(:,0)

            !////////Transfoming Stantard Gaussian Field
            if(rg == 0) write(*,*) ">>>>Transforming Stantard Random Field (PML)"
            !i = 0
            !if(rg == 0) write(*,*) "Dens------------ "
            !if(rg == 0) write(*,*) "margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
            !if(rg == 0) write(*,*) "average      = ", avgProp(i)
            !if(rg == 0) write(*,*) "variance     = ", Tdomain%sSubDomain(mat)%varProp(i)
            !i = 1
            !if(rg == 0) write(*,*) "Lambda----------- "
            !if(rg == 0) write(*,*) "margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
            !if(rg == 0) write(*,*) "average      = ", avgProp(i)
            !if(rg == 0) write(*,*) "variance     = ", Tdomain%sSubDomain(mat)%varProp(i)
            !i = 2
            !if(rg == 0) write(*,*) "Mu--------------- "
            !if(rg == 0) write(*,*) "margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
            !if(rg == 0) write(*,*) "average      = ", avgProp(i)
            !if(rg == 0) write(*,*) "variance     = ", Tdomain%sSubDomain(mat)%varProp(i)

            !            do i = 0, nProp - 1
            !                if(Tdomain%sSubDomain(mat)%varProp(i) > 0) then
            !                    call multiVariateTransformation (          &
            !                        Tdomain%sSubDomain(mat)%margiFirst(i), &
            !                        avgProp(i),                            &
            !                        Tdomain%sSubDomain(mat)%varProp(i),    &
            !                        propPML(:, i:i))
            !                else
            !                    propPML(:,i) = avgProp(i)
            !                end if
            !            end do

            do i = 0, nProp - 1
                call multiVariateTransformation (          &
                    Tdomain%sSubDomain(mat)%margiFirst(i), &
                    avgProp(i),                            &
                    Tdomain%sSubDomain(mat)%varProp(i),    &
                    propPML(:, i:i))
            end do
            !write(*,*) "propPML AFTER TRANSFORMATION = ", propPML(:,0)

            !// Propagating Properties over the PML

            if(rg == 0) write(*,*) ">>>>Propagating Properties over the PML"

            do m = 0, Tdomain%sSubDomain(mat)%nElem - 1
                n = Tdomain%sSubDomain(mat)%elemList(m)
                !Face X oriented
                if(face_X) then

                    do j = 0, nglly-1
                        do k = 0, ngllz-1
                            ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1,j,k)
                            RFpoint      = count(PML_Mask(0,0:ipoint)) - 1
                            pointProp(:) = propPML(RFpoint, :)
                            do i = 0, ngllx-1
                                ipoint           = Tdomain%specel(n)%Iglobnum(i,j,k)
                                RFpoint          = count(Tdomain%sSubDomain(mat)%globCoordMask(0,0:ipoint)) - 1
                                prop(RFpoint, :) = pointProp(:)
                            end do
                        enddo
                    enddo
                    !Face Y oriented
                elseif (face_Y) then
                    do i = 0, ngllx-1
                        do k = 0, ngllz-1
                            !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                            ipoint       = Tdomain%specel(n)%Iglobnum(i, LimPML2, k)
                            RFpoint      = count(PML_Mask(0,0:ipoint)) - 1
                            pointProp(:) = propPML(RFpoint, :)
                            do j = 0, nglly-1
                                ipoint           = Tdomain%specel(n)%Iglobnum(i,j,k)
                                RFpoint          = count(Tdomain%sSubDomain(mat)%globCoordMask(0,0:ipoint)) - 1
                                prop(RFpoint, :) = pointProp(:)
                            end do
                        enddo
                    enddo

                    !Face Z oriented
                elseif (face_Z) then
                    do i = 0, ngllx-1
                        do j = 0, nglly-1
                            !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                            ipoint       = Tdomain%specel(n)%Iglobnum(i,j,LimPML3)
                            RFpoint      = count(PML_Mask(0,0:ipoint)) - 1
                            pointProp(:) = propPML(RFpoint, :)
                            do k = 0, ngllz-1
                                ipoint           = Tdomain%specel(n)%Iglobnum(i,j,k)
                                RFpoint          = count(Tdomain%sSubDomain(mat)%globCoordMask(0,0:ipoint)) - 1
                                prop(RFpoint, :) = pointProp(:)
                            end do
                        enddo
                    enddo

                    !Edge in XY
                elseif (edge_XY) then
                    do k = 0, ngllz-1
                        !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                        ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1,LimPML2,k)
                        RFpoint      = count(PML_Mask(0,0:ipoint)) - 1
                        pointProp(:) = propPML(RFpoint, :)
                        do i = 0, ngllx-1
                            do j = 0, nglly-1
                                ipoint           = Tdomain%specel(n)%Iglobnum(i,j,k)
                                RFpoint          = count(Tdomain%sSubDomain(mat)%globCoordMask(0,0:ipoint)) - 1
                                prop(RFpoint, :) = pointProp(:)
                            end do
                        end do
                    enddo

                    !Edge in YZ
                elseif (edge_YZ) then
                    do i = 0, ngllx-1
                        !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                        ipoint       = Tdomain%specel(n)%Iglobnum(i, LimPML2,LimPML3)
                        RFpoint      = count(PML_Mask(0,0:ipoint)) - 1
                        pointProp(:) = propPML(RFpoint, :)
                        do j = 0, nglly-1
                            do k = 0, ngllz-1
                                ipoint           = Tdomain%specel(n)%Iglobnum(i,j,k)
                                RFpoint          = count(Tdomain%sSubDomain(mat)%globCoordMask(0,0:ipoint)) - 1
                                prop(RFpoint, :) = pointProp(:)
                            end do
                        end do
                    enddo

                    !Edge in ZX
                elseif (edge_ZX) then
                    do j = 0, nglly-1
                        ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1, j,LimPML3)
                        RFpoint      = count(PML_Mask(0,0:ipoint)) - 1
                        pointProp(:) = propPML(RFpoint, :)
                        do i = 0, ngllx-1
                            do k = 0, ngllz-1
                                ipoint           = Tdomain%specel(n)%Iglobnum(i,j,k)
                                RFpoint          = count(Tdomain%sSubDomain(mat)%globCoordMask(0,0:ipoint)) - 1
                                prop(RFpoint, :) = pointProp(:)
                            end do
                        end do
                    enddo

                    !Vertex
                elseif (vertex_XYZ) then
                    ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1, LimPML2,LimPML3)
                    RFpoint      = count(PML_Mask(0,0:ipoint)) - 1
                    pointProp(:) = propPML(RFpoint, :)
                    do i = 0, ngllx-1
                        do j = 0, nglly-1
                            do k = 0, ngllz-1
                                ipoint           = Tdomain%specel(n)%Iglobnum(i,j,k)
                                RFpoint          = count(Tdomain%sSubDomain(mat)%globCoordMask(0,0:ipoint)) - 1
                                prop(RFpoint, :) = pointProp(:)
                            end do
                        end do
                    end do
                end if
                !Loop over subdomain elements
            end do

            deallocate(PML_Mask)
            deallocate(propPML)

            !write(*,*) "prop AFTER PROPAGATION (PML) = ", prop(:,0)
        end if
        !END Non-PML/PML if

        deallocate(calculate)

    end subroutine build_random_properties

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine read_properties_from_file(Tdomain, rg, mat, prop,  &
        fileName, folderPath, labels, indexes)

        use sem_hdf5
        use hdf5

        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: rg, mat
        character (len=*)                , intent(in) :: filename;
        character(len=*)                 , intent(in) :: folderPath
        character(len=*) , dimension(1:), optional  , intent(in) :: labels
        integer          , dimension(1:), optional  , intent(in) :: indexes
        !OUTPUT
        real        , intent(out), dimension(0:, 0:) :: prop !Properties
        !LOCAL
        character (len=12) :: numberStr, rangStr;
        character(len=110) :: fileHDF5Name, fullPath !File name
        integer            :: i, error
        integer(HID_T)     :: file_id
        real,    allocatable, dimension(:,:) :: rtemp2

        write(*,*) "Reading properties from file"

        if(.not. present(labels)) then
            write(rangStr,'(I8)'  ) rg
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

        write(*,*) "fileHDF5Name =",fileHDF5Name
        write(*,*) "fullPath =",fullPath

        file_id = 12
        call h5open_f(error) ! Initialize FORTRAN interface.
        call h5fopen_f(trim(fullPath), H5F_ACC_RDONLY_F, file_id, error)
        do i = 0, size(prop, 2)-1
            write(*,*) "i = ", i
            write(*,*) "fullPath = ", fullPath
            write(*,*) "trim(stringNumb_join('RF_',i+1)) = ", trim(stringNumb_join("RF_",i+1))
            call read_dataset(file_id, trim(stringNumb_join("RF_",i+1)), rtemp2)
            !call dispCarvalhol(rtemp2, "rtemp2", "F30.10")
            prop(:,i:i) = rtemp2(:,:)

            deallocate(rtemp2)
        end do
        call h5fclose_f(file_id, error) ! Close the file.
        call h5close_f(error) ! Close FORTRAN interface.

    end subroutine read_properties_from_file


    !TRASH-------------------------

    !    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !    subroutine allocate_random_subdomain(Tdomain, rg, mat)
    !        implicit none
    !
    !        !INPUT & OUTPUT
    !        type (domain), intent(inout), target :: Tdomain
    !           integer      , intent(in)            :: rg, mat
    !
    !           !LOCAL
    !           integer :: n !Counters
    !           integer :: ipoint
    !
    !        call random_seed(size = n)
    !        allocate(Tdomain%sSubdomain(mat)%chosenSeed(n))
    !        allocate(Tdomain%sSubDomain(mat)%MinBound(0:2))
    !        allocate(Tdomain%sSubDomain(mat)%MaxBound(0:2))
    !
    !    end subroutine allocate_random_subdomain
    !    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !    subroutine define_random_PML(Tdomain, rg, mat)
    !        implicit none
    !        !INPUT
    !        type(domain), intent(inout), target :: Tdomain
    !        integer     , intent(in)  :: rg, mat
    !
    !        !LOCAL
    !        integer :: assocMat, n
    !        integer :: error, code
    !
    !        assocMat = Tdomain%sSubDomain(mat)%assocMat
    !        if(.not.(.not. (Tdomain%not_PML_List(mat))            .and. &
    !           Tdomain%sSubdomain(assocMat)%material_type == "R" .and. &
    !           Tdomain%subD_exist(mat))) then !PMLs associated to random subdomains
    !            write(*,*) "!!!ERROR: 'define_random_PML' was called wrongly"
    !            write(*,*) "Rang =", rg, "mat = ", mat
    !            write(*,*) "material_type =", Tdomain%sSubDomain(mat)%material_type, ", existence = ", Tdomain%subD_exist(mat)
    !            call MPI_ABORT(Tdomain%communicateur, error, code)
    !        end if
    !
    !        !write(*,*) "allocated(Tdomain%sSubdomain(mat)%chosenSeed(n)) = ", allocated(Tdomain%sSubdomain(mat)%chosenSeed(n))
    !        !write(*,*) "allocated(Tdomain%sSubDomain(mat)%MinBound) = ", allocated(Tdomain%sSubDomain(mat)%MinBound)
    !        !write(*,*) "allocated(Tdomain%sSubDomain(mat)%MaxBound) = ", allocated(Tdomain%sSubDomain(mat)%MaxBound)
    !        !call random_seed(size = n)
    !        !allocate(Tdomain%sSubdomain(mat)%chosenSeed(n))
    !        !allocate(Tdomain%sSubDomain(mat)%MinBound(0:2))
    !        !allocate(Tdomain%sSubDomain(mat)%MaxBound(0:2))
    !
    !        Tdomain%sSubdomain(mat)%chosenSeed(:) = Tdomain%sSubdomain(assocMat)%chosenSeed(:)
    !        Tdomain%sSubdomain(mat)%MinBound(:)   = Tdomain%sSubdomain(assocMat)%MinBound(:)
    !        Tdomain%sSubdomain(mat)%MaxBound(:)   = Tdomain%sSubdomain(assocMat)%MaxBound(:)
    !
    !    end subroutine define_random_PML


end module define_random
