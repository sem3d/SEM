!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module define_random

    use sdomain
    use mpi
    use constants
    use RF_fileCreator

    implicit none


contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine build_random_properties(Tdomain, rg)
!        use calls_RF
        implicit none
        !INPUT
        type(domain)    , intent(inout), target :: Tdomain
        integer         , intent(in)            :: rg
        !LOCAL
!        type(IPT_RF) :: IPT
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
        double precision :: varProp
        

        propName(0) = "Density"
        propName(1) = "Kappa"
        propName(2) = "Mu"
        propCounter = 0

        main_fId = 25
!        if(rg == 0) open (unit = main_fId , file = "./mat/main_input", action = 'write')

        do mat = 0, Tdomain%n_mat - 1

!            Density = Tdomain%sSubDomain(mat)%Ddensity
!            P_Speed = Tdomain%sSubDomain(mat)%Pspeed
!            S_Speed = Tdomain%sSubDomain(mat)%Sspeed
!
!            Mu = S_Speed**2 * Density
!            Lambda = (P_Speed**2 - 2 * S_Speed**2 ) * Density
!            Kappa = Lambda + 2.0D0*Mu /3.0D0

            assocMat = Tdomain%sSubdomain(mat)%assocMat

!            avgProp  = [Density, &
!                        Kappa,  &
!                        Mu]

            !stop("Random Properties not yet functional on SEM")


!            nProcsPerChunk = Tdomain%nb_procs
!            if(nProcsPerChunk > 12) nProcsPerChunk = 12
!            nChunks = ceiling(dble(Tdomain%nb_procs)/dble(nProcsPerChunk))

            do propId = 0, nProp - 1
            !do propId = 0, 0 !FOR TESTS
                Tdomain%sSubDomain(mat)%propFilePath(propId) = string_join_many("./mat/h5/", "Mat_", numb2String(assocMat),"_",propName(propId),".h5")
!                if(.not. is_rand(Tdomain%sSubdomain(mat))) cycle
!
!                if(rg == 0) write(*,*) "  Material ", mat, "is Random"
!                if(rg == 0) write(*,*) "  avgProp = ", avgProp
!
!                if(is_rand(Tdomain%sSubdomain(mat)) .and. rg == 0) then
!
!                    write(*,*) " rang ", rg," in build_random_properties"
!                    if(rg == 0) write(*,*) "  Generating Random Properties"
!                    seedStart = Tdomain%sSubDomain(mat)%seedStart
!                    if(seedStart >= 0) seedStart = seedStart + 7*(propId+1)
!
!                    write(*,*) "  Generating Random Properties Files"
!                    propCounter = propCounter + 1
!                    write(*,*) "  propCounter = ", propCounter, "rank = ", rg
!                    varProp = (avgProp(propId)*Tdomain%sSubdomain(mat)%varCoef(propId))**2D0
!
!                    fileNameBase = string_join_many(propName(propId), "_Mat_", numb2String(mat))
!
!                    call makeCase(&
!                            nDim=3, &
!                            Nmc=1, &
!                            corrMod=Tdomain%sSubdomain(mat)%corrMod, &
!                            margiFirst=Tdomain%sSubdomain(mat)%margiFirst(propId), &
!                            corrL=Tdomain%sSubdomain(mat)%corrL, &
!                            fieldAvg=avgProp(propId), &
!                            fieldVar=varProp, &
!                            method=4, &
!                            seedStart=seedStart, &
!                            overlap=[0.0D0, 0.0D0, 0.0D0], &
!                            xMinGlob=Tdomain%sSubdomain(mat)%MinBound, &
!                            xMaxGlob=Tdomain%sSubdomain(mat)%MaxBound, &
!                            pointsPerCorrL=[5, 5, 5], &
!                            nProcsTotal=Tdomain%nb_procs, &
!                            nProcsPerChunk=nProcsPerChunk, &
!                            localizationLevel=1, &
!                            nFields=[1, 1, 1], &
!                            nChunks=nChunks, &
!                            memPerChunk=12000, &
!                            queue="iceq", &
!                            wallTime="2:00:00", &
!                            cluster=1, &
!                            folderPath="./mat/input", &
!                            runPath=runPath, &
!                            rfPath=string_join_many(Tdomain%random_library_path,"/randomField.exe"), &
!                            statPath=string_join_many(Tdomain%random_library_path,"/statistics.exe"), &
!                            gen_input_name="gen_"//fileNameBase, &
!                            mesh_input_name="mesh_"//fileNameBase)
!
!                    write(main_fId, *) trim(string_join_many("$mesh_input_",numb2String(propCounter)))//' '//&
!                                       trim(string_join_many('"./input/mesh_'//fileNameBase,'"'))
!                    write(main_fId, *) trim(string_join_many("$gen_input_",numb2String(propCounter)))//' '//&
!                                       trim(string_join_many('"./input/gen_'//fileNameBase,'"'))
!                    write(main_fId, *) trim(string_join_many("$out_folder_",numb2String(propCounter)))//' '//'"."'
!                    write(main_fId, *) trim(string_join_many("$out_name_",numb2String(propCounter)))//' '//&
!                                       trim(fileNameBase)
!                end if


            end do

        end do

!        if(rg == 0) write(main_fId,*) "$nSamples ", propCounter
!        if(rg == 0) close(main_fId)

!        call MPI_BARRIER(Tdomain%communicateur, code)

    end subroutine build_random_properties




!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    subroutine define_random_seed(Tdomain, rg, mat)
!        !WARNING: *The subdomains communicators should've been created before this subroutine
!        !         *"allocate_random_subdomain" should be called before this subroutine
!        implicit none
!        !INPUT
!        type(domain), intent(inout), target :: Tdomain
!        integer     , intent(in)  :: rg, mat
!        !LOCAL
!        integer :: code, error
!
!        if(Tdomain%sSubdomain(mat)%initial_material_type /= 'R') then
!            write(*,*) "!!!ERROR:'define_random_seed' was called wrongly"
!            write(*,*) "Rang =", rg, "mat = ", mat
!            write(*,*) "material_type =", Tdomain%sSubDomain(mat)%material_type, ", existence = ", Tdomain%subD_exist(mat)
!            call MPI_ABORT(Tdomain%communicateur, error, code)
!        end if
!
!
!        if(rg == 0) then
!            if(Tdomain%sSubdomain(mat)%seedStart >= 0) then
!                call calculate_random_seed(Tdomain%sSubdomain(mat)%chosenSeed, Tdomain%sSubdomain(mat)%seedStart)
!            else
!                call calculate_random_seed(Tdomain%sSubdomain(mat)%chosenSeed)
!            end if
!        end if
!
!        call MPI_BCAST (Tdomain%sSubdomain(mat)%chosenSeed,             &
!            size(Tdomain%sSubdomain(mat)%chosenSeed),       &
!            MPI_INTEGER, 0, Tdomain%communicateur, code)
!    end subroutine define_random_seed

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
