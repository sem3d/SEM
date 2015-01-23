module define_random

    use sdomain
    use mpi
    use constants
    use randomFieldND
    use displayCarvalhol
    use writeResultFile_RF

contains
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    subroutine define_random_subdomains(Tdomain, rg)
        !This routine define the seed for random subdomains
        !and spreads the seed and Max/Min limits to PMLs interfacin random subdomains
        implicit none
        !INPUT
        type(domain), intent(inout), target :: Tdomain
        integer     , intent(in)  :: rg

        !LOCAL
        integer :: mat, assocMat
        integer :: seedSize

        call random_seed(size = seedSize)

        !Non PMLs
        do mat = 0, Tdomain%n_mat - 1
            if(Tdomain%sSubDomain(mat)%material_type == "R") then
                allocate(Tdomain%sSubdomain(mat)%chosenSeed(seedSize))
                call define_random_seed(Tdomain, rg, mat)
            end if
        end do

        !PMLS
        do mat = 0, Tdomain%n_mat - 1
            assocMat = Tdomain%sSubDomain(mat)%assocMat

            if(.not. (Tdomain%not_PML_List(mat))                    .and. &
                Tdomain%sSubdomain(assocMat)%material_type == "R" ) then !PMLs associated to random subdomains
                allocate(Tdomain%sSubdomain(mat)%chosenSeed(seedSize))
                Tdomain%sSubdomain(mat)%chosenSeed(:) = Tdomain%sSubdomain(assocMat)%chosenSeed(:)
                Tdomain%sSubdomain(mat)%MinBound(:)   = Tdomain%sSubdomain(assocMat)%MinBound(:)
                Tdomain%sSubdomain(mat)%MaxBound(:)   = Tdomain%sSubdomain(assocMat)%MaxBound(:)
            end if
            !write(*,*) "Chosen seed mat ", mat, "= ", Tdomain%sSubdomain(mat)%chosenSeed(:)
        end do


    end subroutine define_random_subdomains

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    subroutine build_random_properties(Tdomain, rg, mat, prop, method)
        !This routine fills up the properties table (prop) using chosen method to create the samples
        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: rg, mat
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
        real, dimension(0:2) :: pointProp, avgProp;

        !Defining which properties will be calculated
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

!        if(rg == 0) write(*,*) ">>>>Creating Standard Gaussian Field mat = ", mat
!        if(rg == 0) write(*,*) "corrL                = ", Tdomain%sSubDomain(mat)%corrL
!        if(rg == 0) write(*,*) "corrMod              = ", Tdomain%sSubDomain(mat)%corrMod
!        if(rg == 0) write(*,*) "Nmc (nProp)          = ", nProp
!        if(rg == 0) write(*,*) "chosenSeed           = ", Tdomain%sSubDomain(mat)%chosenSeed
!        if(rg == 0) write(*,*) "MinBound             = ", Tdomain%sSubDomain(mat)%MinBound
!        if(rg == 0) write(*,*) "MaxBound             = ", Tdomain%sSubDomain(mat)%MaxBound

        !call dispCarvalhol(prop(1:20,:), "prop(1:20,:) RAW", "F30.10")

        select case(effecMethod)
        case( 1 ) !Victor
            write(*,*) "Victor's method"
            call createStandardGaussianFieldUnstructVictor(&
                Tdomain%GlobCoord(:, :),                   &
                Tdomain%sSubDomain(mat)%corrL,             &
                Tdomain%sSubDomain(mat)%corrMod,           &
                nProp,                                     &
                prop(:, :),                                &
                Tdomain%sSubDomain(mat)%chosenSeed,        &
                Tdomain%sSubDomain(mat)%MinBound,          &
                Tdomain%sSubDomain(mat)%MaxBound,          &
                Tdomain%communicateur,                     &
                calculate)

        case( 2 ) !Shinozuka
            write(*,*) "Shinozuka's method"
            call createStandardGaussianFieldUnstructShinozuka (&
                Tdomain%GlobCoord,                             &
                Tdomain%sSubDomain(mat)%corrL,                 &
                Tdomain%sSubDomain(mat)%corrMod,               &
                nProp,                                         &
                prop(:, :),                                    &
                Tdomain%sSubDomain(mat)%chosenSeed,            &
                Tdomain%sSubDomain(mat)%MinBound,              &
                Tdomain%sSubDomain(mat)%MaxBound,              &
                Tdomain%communicateur,                         &
                calculate)

        case default
            write(*,*) "ERROR! The chosen method is not an avaiable choice"
            call MPI_ABORT(Tdomain%communicateur, error, code)
        end select

        !call dispCarvalhol(prop(1:20,:), "prop(1:20,:) GAUSS", "F30.10")

!        !////////Transfoming Stantard Gaussian Field
!        if(rg == 0) write(*,*) ">>>>Transforming Standard Gaussian Field mat = ", mat
!
!        i = 0
!        if(rg == 0) write(*,*) "Dens------------ "
!        if(rg == 0) write(*,*) "margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
!        if(rg == 0) write(*,*) "average      = ", avgProp(i)
!        if(rg == 0) write(*,*) "variance     = ", Tdomain%sSubDomain(mat)%varProp(i)
!        i = 1
!        if(rg == 0) write(*,*) "Lambda----------- "
!        if(rg == 0) write(*,*) "margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
!        if(rg == 0) write(*,*) "average      = ", avgProp(i)
!        if(rg == 0) write(*,*) "variance     = ", Tdomain%sSubDomain(mat)%varProp(i)
!        i = 2
!        if(rg == 0) write(*,*) "Mu--------------- "
!        if(rg == 0) write(*,*) "margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
!        if(rg == 0) write(*,*) "average      = ", avgProp(i)
!        if(rg == 0) write(*,*) "variance     = ", Tdomain%sSubDomain(mat)%varProp(i)

        do i = 0, nProp - 1
            call multiVariateTransformation (          &
                Tdomain%sSubDomain(mat)%margiFirst(i), &
                avgProp(i),                            &
                Tdomain%sSubDomain(mat)%varProp(i),    &
                prop(:, i:i))
        end do

        !call dispCarvalhol(prop(1:20,:), "prop(1:20,:) TRANSF", "F30.10")

        if(.not. Tdomain%not_PML_List(mat)) then !Random PML
            call propagate_PML_properties(Tdomain, rg, prop)
        end if

        deallocate(calculate)

    end subroutine build_random_properties


    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine propagate_PML_properties(Tdomain, rg, prop)
        !This subroutine modifies the "prop" table.
        !In the PML elements it extrudes the properties according to the PML direction

        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: rg

        !OUTPUT
        real        , intent(inout), dimension(0:, 0:) :: prop !Properties that should be modified

        !LOCAL
        integer :: ngllx,nglly,ngllz, assocMat
        integer :: iPoint
        integer :: i, j, k, m, n !counters
        integer :: LimPML1, LimPML2, LimPML3
        integer :: dir, mat_index
        real, dimension(:), allocatable :: pointProp;
        logical :: verbose = .false.

        allocate(pointProp(0:size(prop,2)-1))

        !Propagating Properties over the PML
        if(rg == 0) write(*,*) ">>>>Propagating Properties over the PML"

        do n = 0, Tdomain%n_elem-1
            mat_index = Tdomain%specel(n)%mat_index
            if(.not. Tdomain%not_PML_List(mat_index)) then

                ngllx = Tdomain%sSubDomain(mat_index)%NGLLx
                nglly = Tdomain%sSubDomain(mat_index)%NGLLy
                ngllz = Tdomain%sSubDomain(mat_index)%NGLLz
                !On the lower bound initialization
                LimPML1 = 0
                LimPML2 = 0
                LimPML3 = 0
                !To the upper bound transformation
                if (Tdomain%sSubDomain(mat_index)%Left)    LimPML1 = ngllx-1
                if (Tdomain%sSubDomain(mat_index)%Forward) LimPML2 = nglly-1
                if (Tdomain%sSubDomain(mat_index)%Down)    LimPML3 = ngllz-1

                dir = read_PML_Direction(Tdomain, mat_index)
                if(verbose) write(*,*) "dir = ", dir

                select case(dir)
                    !Face X oriented
                    case(0)
                        if(verbose) write(*,*) "Face X oriented"
                        do j = 0, nglly-1
                            do k = 0, ngllz-1
                                ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1,j,k)
                                pointProp(:) = prop(ipoint, :)
                                do i = 0, ngllx-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            enddo
                        enddo
!
                    !Face Y oriented
                    case(1)
                        if(verbose) write(*,*) "Face Y oriented"
                        do i = 0, ngllx-1
                            do k = 0, ngllz-1
                                ipoint       = Tdomain%specel(n)%Iglobnum(i, LimPML2, k)
                                pointProp(:) = prop(ipoint, :)
                                do j = 0, nglly-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            enddo
                        enddo

                    !Face Z oriented
                    case(2)
                        if(verbose) write(*,*) "Face Z oriented"
                        do i = 0, ngllx-1
                            do j = 0, nglly-1
                                !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                                ipoint       = Tdomain%specel(n)%Iglobnum(i,j,LimPML3)
                                pointProp(:) = prop(ipoint, :)
                                do k = 0, ngllz-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            enddo
                        enddo

                    !Edge in XY
                    case(3)
                        if(verbose) write(*,*) "Edge in XY"
                        do k = 0, ngllz-1
                            !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                            ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1,LimPML2,k)
                            pointProp(:) = prop(ipoint, :)
                            do i = 0, ngllx-1
                                do j = 0, nglly-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            end do
                        enddo

                    !Edge in YZ
                    case(4)
                        if(verbose) write(*,*) "Edge in YZ"
                        do i = 0, ngllx-1
                            !if(ipoint > size(randMu, 1)) write (*,*) "ERROR ipoint = ", ipoint, "and size(randMu. 1) = ", size(randMu, 1)
                            ipoint       = Tdomain%specel(n)%Iglobnum(i, LimPML2,LimPML3)
                            pointProp(:) = prop(ipoint, :)
                            do j = 0, nglly-1
                                do k = 0, ngllz-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            end do
                        enddo

                    !Edge in ZX
                    case(5)
                        if(verbose) write(*,*) "Edge in ZX"
                        do j = 0, nglly-1
                            ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1, j,LimPML3)
                            pointProp(:) = prop(ipoint, :)
                            do i = 0, ngllx-1
                                do k = 0, ngllz-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            end do
                        enddo

                    !Vertex
                    case(6)
                        if(verbose) write(*,*) "Vertex"
                        ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1, LimPML2,LimPML3)
                        pointProp(:) = prop(ipoint, :)
                        do i = 0, ngllx-1
                            do j = 0, nglly-1
                                do k = 0, ngllz-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            end do
                        end do
                end select
            end if
        end do

        deallocate(pointProp)

    end subroutine propagate_PML_properties

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    subroutine define_random_seed(Tdomain, rg, mat)
        !WARNING: *The subdomains communicators should've been created before this subroutine
        !         *"allocate_random_subdomain" should be called before this subroutine
        implicit none
        !INPUT
        type(domain), intent(inout), target :: Tdomain
        integer     , intent(in)  :: rg, mat
        !LOCAL
        integer :: code, error
        integer(kind= MPI_OFFSET_KIND ) :: filePos
        integer, dimension( MPI_STATUS_SIZE ) :: status
        integer :: nbBitesInt, seedSize

        if(.not.(Tdomain%sSubDomain(mat)%material_type == "R" .and. Tdomain%subD_exist(mat))) then
            write(*,*) "!!!ERROR:'define_random_seed' was called wrongly"
            write(*,*) "Rang =", rg, "mat = ", mat
            write(*,*) "material_type =", Tdomain%sSubDomain(mat)%material_type, ", existence = ", Tdomain%subD_exist(mat)
            call MPI_ABORT(Tdomain%communicateur, error, code)
        end if


        if(rg == 0) then
            if(Tdomain%sSubdomain(mat)%seedStart >= 0) then
                call calculate_random_seed(Tdomain%sSubdomain(mat)%chosenSeed, Tdomain%sSubdomain(mat)%seedStart)
            else
                call calculate_random_seed(Tdomain%sSubdomain(mat)%chosenSeed)
            end if
            !write(*,*) "mat        = ", mat
            !write(*,*) "nbBitesInt = ", nbBitesInt
            !write(*,*) "sizeRand   = ", sizeRand
            !write(*,*) "mat*nbBitesInt*sizeRand = ", mat*nbBitesInt*sizeRand
            !write(*,*) "...%chosenSeed          = ", Tdomain%sSubdomain(mat)%chosenSeed
        end if

        call MPI_BCAST (Tdomain%sSubdomain(mat)%chosenSeed,             &
            size(Tdomain%sSubdomain(mat)%chosenSeed),       &
            MPI_INTEGER, 0, Tdomain%communicateur, code)
        write(*,*)"rg = ", rg, "..%chosenSeed   = ", Tdomain%sSubdomain(mat)%chosenSeed
    end subroutine define_random_seed

    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    function read_PML_Direction(Tdomain, mat) result(dir)

        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: mat

        !OUTPUT
        integer :: dir

        !LOCAL
        integer :: error, code

        !write(*,*) ">>>>Defining PML orientation"
        !/////////////Defining PML orientation
        !Face X oriented
        if  (        Tdomain%sSubDomain(mat)%Px   .and. &
            (.not.Tdomain%sSubDomain(mat)%Py)  .and. &
            (.not.Tdomain%sSubDomain(mat)%Pz)) then
            dir = 0
            !write(*,*) "face_X"
        !Face Y oriented
        elseif ((.not.Tdomain%sSubDomain(mat)%Px)  .and. &
            Tdomain%sSubDomain(mat)%Py   .and. &
            (.not.Tdomain%sSubDomain(mat)%Pz)) then
            dir = 1
            !write(*,*) "face_Y"
        !Face Z oriented
        elseif ((.not.Tdomain%sSubDomain(mat)%Px) .and. &
            (.not.Tdomain%sSubDomain(mat)%Py) .and. &
            Tdomain%sSubDomain(mat)%Pz) then
            dir = 2
            !write(*,*) "face_Z"

        !Edge in XY
        elseif (     Tdomain%sSubDomain(mat)%Px  .and. &
            (     Tdomain%sSubDomain(mat)%Py) .and. &
            (.not.Tdomain%sSubDomain(mat)%Pz)) then
            dir = 3
            !write(*,*) "edge_XY"

        !Edge in YZ
        elseif ((.not.Tdomain%sSubDomain(mat)%Px) .and. &
            (      Tdomain%sSubDomain(mat)%Py) .and. &
            (      Tdomain%sSubDomain(mat)%Pz)) then
            dir = 4
            !write(*,*) "edge_YZ"

       !Edge in ZX
        elseif (      Tdomain%sSubDomain(mat)%Px   .and. &
            ((.not.Tdomain%sSubDomain(mat)%Py)) .and. &
            (      Tdomain%sSubDomain(mat)%Pz)) then
            dir = 5
            !write(*,*) "edge_ZX"

       !Vertex in XYZ
        elseif (  Tdomain%sSubDomain(mat)%Px   .and. &
            (  Tdomain%sSubDomain(mat)%Py)  .and. &
            (  Tdomain%sSubDomain(mat)%Pz)) then
            dir = 6
            !write(*,*) "vertex_XYZ"
        !Undefined PML
        else
            write(*,*) "ERROR in mat ", mat, " (PML) definition (directions), check 'material.input'"
            write(*,*) "Tdomain%sSubDomain(", mat, ")%Px = ", Tdomain%sSubDomain(mat)%Px
            write(*,*) "Tdomain%sSubDomain(", mat, ")%Py = ", Tdomain%sSubDomain(mat)%Py
            write(*,*) "Tdomain%sSubDomain(", mat, ")%Pz = ", Tdomain%sSubDomain(mat)%Pz
            call MPI_ABORT(Tdomain%communicateur, error, code)
        end if

    end function read_PML_Direction


end module define_random
