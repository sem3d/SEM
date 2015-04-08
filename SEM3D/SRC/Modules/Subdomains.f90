!>
!!\file Subdomains.f90
!!\brief Calcul des coefficients de Lame.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssubdomains

    type Subdomain

        logical :: Filtering, Px, Py, Pz, Left, Forward, Down

        integer :: NGLLx, NGLLy, NGLLz, npow

        !  modif mariotti fevrier 2007 cea
        ! Qmu en plus
        real :: Pspeed, Sspeed, Ddensity, Dt, Apow, freq, DLambda, DMu, Qmu, Q
        real :: DKappa, Qpression
        real, dimension (:), pointer :: GLLcx, GLLpolx, GLLwx
        real, dimension (:,:), pointer :: hprimex, hTprimex
        real, dimension (:), pointer :: GLLcy, GLLpoly, GLLwy
        real, dimension (:,:), pointer :: hprimey, hTprimey
        real, dimension (:), pointer :: GLLcz, GLLpolz, GLLwz
        real, dimension (:,:), pointer :: hprimez, hTprimez

        character(len=1) :: material_type
        integer          :: material_definition
        !Modification to accept random media
        character(len = 15) :: corrMod
        integer             :: assocMat = -1
        integer             :: seedStart
        integer             :: nElem = 0 !number of elements in each subdomain (by proc) - mesh3d.f90(362)
        integer            , dimension(:)   , allocatable :: elemList !List of elements in "Tdomain%specel(:)" that belong to this subdomain (by proc)
        character(len = 30), dimension(:)   , allocatable :: margiFirst
        real               , dimension(:)   , allocatable :: varProp
        integer            , dimension(:)   , allocatable :: chosenSeed
        logical            , dimension(:,:) , allocatable :: globCoordMask
        real               , dimension(:)   , allocatable :: corrL
        real               , dimension(:)   , allocatable :: MinBound, MaxBound

    end type Subdomain

contains

    !>
    !! \fn subroutine Lame_coefficients (S)
    !! \brief
    !!
    !! \param type (Subdomain) S
    !<
    subroutine Lame_coefficients (S)

        type (Subdomain) :: S

        S%DMu = S%Sspeed**2 * S%Ddensity
        S%DLambda = (S%Pspeed**2 - 2 * S%Sspeed **2 ) * S%Ddensity
        S%DKappa = S%DLambda + 2.*S%DMu /3.

    end subroutine Lame_coefficients

end module ssubdomains
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
