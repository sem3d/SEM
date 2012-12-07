!>
!!\file Subdomains.f90
!!\brief Calcul des coefficients de Lame.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssubdomains
    ! Gaetano Festa 01/02/2005

    type Subdomain

       logical :: Filtering, Px, Py, Pz, Left, Forward, Down

       integer :: NGLLx, NGLLy, NGLLz, wpml, npow

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

        !!print*,' lame ',S%DMu,S%DLambda ,S%DKappa
    end subroutine Lame_coefficients

end module ssubdomains
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
