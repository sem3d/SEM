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

            assocMat = Tdomain%sSubdomain(mat)%assocMat
            propName(1) = "Kappa"
            if(Tdomain%sSubdomain(assocMat)%lambdaSwitch == 1) propName(1) = "Lambda"
            !write(*,*) "AFTER READ Tdomain%sSubdomain(assocMat)%lambdaSwitch = ", Tdomain%sSubdomain(assocMat)%lambdaSwitch


            do propId = 0, nProp - 1
            !do propId = 0, 0 !FOR TESTS
                Tdomain%sSubDomain(mat)%propFilePath(propId) = string_join_many("./mat/h5/", "Mat_", numb2String(assocMat),"_",propName(propId),".h5")

            end do

        end do

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
