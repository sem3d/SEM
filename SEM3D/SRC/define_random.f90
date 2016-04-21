!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module define_random

    use sdomain
    use mpi
    use constants
    use randomFieldND
    use displayCarvalhol
    use writeResultFile_RF
    use RF_fileCreator

    implicit none


contains
!    !---------------------------------------------------------------------------
!    !---------------------------------------------------------------------------
!    !---------------------------------------------------------------------------
!    !---------------------------------------------------------------------------
!    subroutine define_random_subdomains(Tdomain, rg)
!        !This routine define the seed for random subdomains
!        !and spreads the seed and Max/Min limits to PMLs interfacin random subdomains
!        implicit none
!        !INPUT
!        type(domain), intent(inout), target :: Tdomain
!        integer     , intent(in)  :: rg
!
!        !LOCAL
!        integer :: mat, assocMat, i
!        integer :: seedSize
!
!        call random_seed(size = seedSize)
!
!        !Non PMLs
!        do mat = 0, Tdomain%n_mat - 1
!            if(Tdomain%sSubdomain(mat)%initial_material_type == 'R') then
!                allocate(Tdomain%sSubdomain(mat)%chosenSeed(seedSize))
!                call define_random_seed(Tdomain, rg, mat)
!                if(rg == 0) then
!                    write(*,*) ""
!                    write(*,*) "Material ", mat, " (random)"
!
!                    write(*,*) " INPUTS:"
!                    write(*,*) "  corrL        = ", Tdomain%sSubDomain(mat)%corrL
!                    write(*,*) "  corrMod      = ", Tdomain%sSubDomain(mat)%corrMod
!                    i = 0
!                    write(*,*) "  Dens------------ "
!                    write(*,*) "   average      = ", Tdomain%sSubDomain(mat)%Ddensity
!                    write(*,*) "   variance     = ", Tdomain%sSubDomain(mat)%varProp(i)
!                    write(*,*) "   margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
!                    i = 1
!                    write(*,*) "  Lambda----------- "
!                    write(*,*) "   average      = ", Tdomain%sSubDomain(mat)%DLambda
!                    write(*,*) "   variance     = ", Tdomain%sSubDomain(mat)%varProp(i)
!                    write(*,*) "   margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
!                    i = 2
!                    write(*,*) "  Mu--------------- "
!                    write(*,*) "   average      = ", Tdomain%sSubDomain(mat)%DMu
!                    write(*,*) "   variance     = ", Tdomain%sSubDomain(mat)%varProp(i)
!                    write(*,*) "   margiFirst   = ", Tdomain%sSubDomain(mat)%margiFirst(i)
!
!                    write(*,*) " COMPUTED:"
!                    write(*,*) "  chosenSeed   = ", Tdomain%sSubDomain(mat)%chosenSeed
!                    write(*,*) "  MinBound     = ", Tdomain%sSubDomain(mat)%MinBound
!                    write(*,*) "  MaxBound     = ", Tdomain%sSubDomain(mat)%MaxBound
!                    write(*,*) ""
!                end if
!            end if
!        end do
!
!        !Propagating properties to PMLS
!        do mat = 0, Tdomain%n_mat - 1
!            assocMat = Tdomain%sSubDomain(mat)%assocMat
!
!            if(.not. (Tdomain%not_PML_List(mat)) .and. &
!               Tdomain%sSubdomain(assocMat)%initial_material_type == 'R') then !PMLs associated to random subdomains
!                allocate(Tdomain%sSubdomain(mat)%varProp(0:nProp-1))
!                allocate(Tdomain%sSubdomain(mat)%corrL(0:nProp-1))
!                allocate(Tdomain%sSubdomain(mat)%margiFirst(0:nProp-1))
!                allocate(Tdomain%sSubdomain(mat)%chosenSeed(seedSize))
!                !Min/Max Bounds already allocated in mesh3d.f90
!                Tdomain%sSubDomain(mat)%Ddensity      = Tdomain%sSubDomain(assocMat)%Ddensity
!                Tdomain%sSubDomain(mat)%DLambda       = Tdomain%sSubDomain(assocMat)%DLambda
!                Tdomain%sSubDomain(mat)%DMu           = Tdomain%sSubDomain(assocMat)%DMu
!                Tdomain%sSubDomain(mat)%varProp       = Tdomain%sSubDomain(assocMat)%varProp
!                Tdomain%sSubdomain(mat)%corrMod       = Tdomain%sSubdomain(assocMat)%corrMod
!                Tdomain%sSubdomain(mat)%corrL(:)      = Tdomain%sSubdomain(assocMat)%corrL(:)
!                Tdomain%sSubdomain(mat)%margiFirst(:) = Tdomain%sSubdomain(assocMat)%margiFirst(:)
!                Tdomain%sSubdomain(mat)%chosenSeed(:) = Tdomain%sSubdomain(assocMat)%chosenSeed(:)
!                Tdomain%sSubdomain(mat)%MinBound(:)   = Tdomain%sSubdomain(assocMat)%MinBound(:)
!                Tdomain%sSubdomain(mat)%MaxBound(:)   = Tdomain%sSubdomain(assocMat)%MaxBound(:)
!
!                if(rg == 0) write(*,*) "Material ", mat, " is a random PML linked to material ", assocMat
!            end if
!        end do
!
!        if(rg == 0) write(*,*) ""
!
!    end subroutine define_random_subdomains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine build_random_properties(Tdomain, rg)
        use calls_RF
        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: rg
        !LOCAL
        type(IPT_RF) :: IPT
        integer :: mat,assocMat
        integer :: i, propCounter
        integer :: error, code
        integer :: nProp = 3
        integer :: effecMethod
        logical, dimension(:), allocatable :: calculate
        real, dimension(0:2) :: avgProp;
        integer :: propId
        double precision, dimension(10) :: times
        character(len=50), dimension(0:2) :: propName
        integer :: seedStart
        character(len=1024) :: runPath, fileNameBase
        integer :: main_fId
        integer :: nProcsPerChunk, nChunks
        double precision :: Density, Lambda, Kappa, Mu
        double precision :: P_Speed, S_Speed
        

        propName(0) = "Density"
        propName(1) = "Kappa"
        propName(2) = "Mu"
        propCounter = 0

        main_fId = 25
        if(rg == 0) open (unit = main_fId , file = "./mat/main_input", action = 'write')

        do mat = 0, Tdomain%n_mat - 1

            Density = Tdomain%sSubDomain(mat)%Ddensity
            P_Speed = Tdomain%sSubDomain(mat)%Pspeed
            S_Speed = Tdomain%sSubDomain(mat)%Sspeed

            Mu = S_Speed**2 * Density
            Lambda = (P_Speed**2 - 2 * S_Speed**2 ) * Density
            Kappa = Lambda + 2.0D0*Mu /3.0D0

            assocMat = Tdomain%sSubdomain(mat)%assocMat

            avgProp  = [Density, &
                        Kappa,  &
                        Mu]

            !stop("Random Properties not yet functional on SEM")


            nProcsPerChunk = Tdomain%nb_procs
            if(nProcsPerChunk > 12) nProcsPerChunk = 12
            nChunks = ceiling(dble(Tdomain%nb_procs)/dble(nProcsPerChunk))

            do propId = 0, nProp - 1
            !do propId = 0, 0 !FOR TESTS
                Tdomain%sSubDomain(mat)%propFilePath(propId) = string_join_many("./mat/h5/",propName(propId), "_Mat_", numb2String(assocMat),".h5")
                if(.not. is_rand(Tdomain%sSubdomain(mat))) cycle

                if(rg == 0) write(*,*) "  Material ", mat, "is Random"
                if(rg == 0) write(*,*) "  avgProp = ", avgProp

                if(is_rand(Tdomain%sSubdomain(mat)) .and. rg == 0) then

                    write(*,*) " rang ", rg," in build_random_properties"
                    if(rg == 0) write(*,*) "  Generating Random Properties"
                    seedStart = Tdomain%sSubDomain(mat)%seedStart
                    if(seedStart >= 0) seedStart = seedStart + 7*(propId+1)

                    write(*,*) "  Generating Random Properties Files"
                    propCounter = propCounter + 1
                    write(*,*) "  propCounter = ", propCounter, "rank = ", rg

                    fileNameBase = string_join_many(propName(propId), "_Mat_", numb2String(mat))

                    call makeCase(&
                            nDim=3, &
                            Nmc=1, &
                            corrMod=Tdomain%sSubdomain(mat)%corrMod, &
                            margiFirst=Tdomain%sSubdomain(mat)%margiFirst(propId), &
                            corrL=Tdomain%sSubdomain(mat)%corrL, &
                            fieldAvg=avgProp(propId), &
                            fieldVar=Tdomain%sSubdomain(mat)%varProp(propId), &
                            method=4, &
                            seedStart=seedStart, &
                            overlap=[0.0D0, 0.0D0, 0.0D0], &
                            xMinGlob=Tdomain%sSubdomain(mat)%MinBound, &
                            xMaxGlob=Tdomain%sSubdomain(mat)%MaxBound, &
                            pointsPerCorrL=[5, 5, 5], &
                            nProcsTotal=Tdomain%nb_procs, &
                            nProcsPerChunk=nProcsPerChunk, &
                            localizationLevel=1, &
                            nFields=[1, 1, 1], &
                            nChunks=nChunks, &
                            memPerChunk=12000, &
                            queue="iceq", &
                            wallTime="2:00:00", &
                            cluster=1, &
                            folderPath="./mat/input", &
                            runPath=runPath, &
                            rfPath=string_join_many(Tdomain%random_library_path,"/randomField.exe"), &
                            statPath=string_join_many(Tdomain%random_library_path,"/statistics.exe"), &
                            gen_input_name="gen_"//fileNameBase, &
                            mesh_input_name="mesh_"//fileNameBase)

                    write(main_fId, *) trim(string_join_many("$mesh_input_",numb2String(propCounter)))//' '//&
                                       trim(string_join_many('"./input/mesh_'//fileNameBase,'"'))
                    write(main_fId, *) trim(string_join_many("$gen_input_",numb2String(propCounter)))//' '//&
                                       trim(string_join_many('"./input/gen_'//fileNameBase,'"'))
                    write(main_fId, *) trim(string_join_many("$out_folder_",numb2String(propCounter)))//' '//'"."'
                    write(main_fId, *) trim(string_join_many("$out_name_",numb2String(propCounter)))//' '//&
                                       trim(fileNameBase)
                end if


            end do

        end do

        if(rg == 0) write(main_fId,*) "$nSamples ", propCounter
        if(rg == 0) close(main_fId)

        call MPI_BARRIER(Tdomain%communicateur, code)

    end subroutine build_random_properties


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine propagate_PML_properties(Tdomain, prop)
        !This subroutine modifies the "prop" table.
        !In the PML elements it extrudes the properties according to the PML direction

        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain

        !OUTPUT
        real        , intent(inout), dimension(0:, 0:) :: prop !Properties that should be modified

        !LOCAL
        integer :: ngll
        integer :: iPoint
        integer :: i, j, k, n !counters
        integer :: LimPML1, LimPML2, LimPML3
        integer :: dir, mat_index
        real, dimension(:), allocatable :: pointProp;
        logical :: verbose = .false.

        allocate(pointProp(0:size(prop,2)-1))

        !Propagating Properties over the PML

        ! With the new PML definition, we can do better:
        ! We can project x,y,z, on the pml plane, and interpolate the properties from the
        ! projected coordinates

        do n = 0, Tdomain%n_elem-1
            mat_index = Tdomain%specel(n)%mat_index
            if(.not. Tdomain%not_PML_List(mat_index)) then

                ngll = Tdomain%sSubDomain(mat_index)%NGLL
                !On the lower bound initialization
                LimPML1 = 0
                LimPML2 = 0
                LimPML3 = 0
                !To the upper bound transformation
                if (Tdomain%sSubDomain(mat_index)%pml_width(0)<0) LimPML1 = ngll-1
                if (Tdomain%sSubDomain(mat_index)%pml_width(1)<0) LimPML2 = ngll-1
                if (Tdomain%sSubDomain(mat_index)%pml_width(2)<0) LimPML3 = ngll-1

                dir = read_PML_Direction(Tdomain, mat_index)
                if(verbose) write(*,*) "dir = ", dir

                select case(dir)
                    !Face X oriented
                    case(0)
                        if(verbose) write(*,*) "Face X oriented"
                        do j = 0, ngll-1
                            do k = 0, ngll-1
                                ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1,j,k)
                                pointProp(:) = prop(ipoint, :)
                                do i = 0, ngll-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            enddo
                        enddo
                    !
                    !Face Y oriented
                    case(1)
                        if(verbose) write(*,*) "Face Y oriented"
                        do i = 0, ngll-1
                            do k = 0, ngll-1
                                ipoint       = Tdomain%specel(n)%Iglobnum(i, LimPML2, k)
                                pointProp(:) = prop(ipoint, :)
                                do j = 0, ngll-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            enddo
                        enddo

                    !Face Z oriented
                    case(2)
                        if(verbose) write(*,*) "Face Z oriented"
                        do i = 0, ngll-1
                            do j = 0, ngll-1
                                ipoint       = Tdomain%specel(n)%Iglobnum(i,j,LimPML3)
                                pointProp(:) = prop(ipoint, :)
                                do k = 0, ngll-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            enddo
                        enddo

                    !Edge in XY
                    case(3)
                        if(verbose) write(*,*) "Edge in XY"
                        do k = 0, ngll-1
                            ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1,LimPML2,k)
                            pointProp(:) = prop(ipoint, :)
                            do i = 0, ngll-1
                                do j = 0, ngll-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            end do
                        enddo

                    !Edge in YZ
                    case(4)
                        if(verbose) write(*,*) "Edge in YZ"
                        do i = 0, ngll-1
                            ipoint       = Tdomain%specel(n)%Iglobnum(i, LimPML2,LimPML3)
                            pointProp(:) = prop(ipoint, :)
                            do j = 0, ngll-1
                                do k = 0, ngll-1
                                    ipoint          = Tdomain%specel(n)%Iglobnum(i,j,k)
                                    prop(ipoint, :) = pointProp(:)
                                end do
                            end do
                        enddo

                    !Edge in ZX
                    case(5)
                        if(verbose) write(*,*) "Edge in ZX"
                        do j = 0, ngll-1
                            ipoint       = Tdomain%specel(n)%Iglobnum(LimPML1, j,LimPML3)
                            pointProp(:) = prop(ipoint, :)
                            do i = 0, ngll-1
                                do k = 0, ngll-1
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
                        do i = 0, ngll-1
                            do j = 0, ngll-1
                                do k = 0, ngll-1
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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine define_random_seed(Tdomain, rg, mat)
        !WARNING: *The subdomains communicators should've been created before this subroutine
        !         *"allocate_random_subdomain" should be called before this subroutine
        implicit none
        !INPUT
        type(domain), intent(inout), target :: Tdomain
        integer     , intent(in)  :: rg, mat
        !LOCAL
        integer :: code, error

        if(Tdomain%sSubdomain(mat)%initial_material_type /= 'R') then
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
        end if

        call MPI_BCAST (Tdomain%sSubdomain(mat)%chosenSeed,             &
            size(Tdomain%sSubdomain(mat)%chosenSeed),       &
            MPI_INTEGER, 0, Tdomain%communicateur, code)
    end subroutine define_random_seed

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    function read_PML_Direction(Tdomain, mat) result(dir)

        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: mat

        !OUTPUT
        integer :: dir

        !LOCAL
        integer :: error, code

        real :: vx, vy, vz

        !!! XXX: this wont work always, read comments from the call site
        vx = Tdomain%sSubDomain(mat)%pml_width(0)
        vy = Tdomain%sSubDomain(mat)%pml_width(1)
        vz = Tdomain%sSubDomain(mat)%pml_width(2)

        !/////////////Defining PML orientation
        dir = 0
        !Face X oriented
        if  (vx /= 0 .and. vy == 0d0 .and. vz == 0d0) then
            dir = 0
        !Face Y oriented
        elseif (vx == 0 .and. vy /= 0d0 .and. vz == 0d0) then
            dir = 1
        !Face Z oriented
        elseif (vx == 0 .and. vy == 0d0 .and. vz /= 0d0) then
            dir = 2
        !Edge in XY
        elseif (vx /= 0 .and. vy /= 0d0 .and. vz == 0d0) then
            dir = 3

        !Edge in YZ
        elseif (vx == 0 .and. vy /= 0d0 .and. vz /= 0d0) then
            dir = 4

        !Edge in ZX
        elseif (vx /= 0 .and. vy == 0d0 .and. vz /= 0d0) then
            dir = 5

        !Vertex in XYZ
        elseif (vx /= 0 .and. vy /= 0d0 .and. vz /= 0d0) then
            dir = 6

        !Undefined PML
        else
            write(*,*) "ERROR in mat ", mat, " (PML) definition (directions), check 'material.input'"
            call MPI_ABORT(Tdomain%communicateur, error, code)
        end if

    end function read_PML_Direction


end module define_random

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
