!>
!!\file Subdomains.F90
!!\brief Calcul des coefficients de Lame.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssubdomains
    ! Gaetano Festa 18/10/2004

    type Subdomain

       integer :: NGLLx, NGLLz, n_loc_dim, wpml, npow, pml_type

       real :: Pspeed, Sspeed, Ddensity, DT, DLambda, DMu, Apow, freq, k
       real, dimension (:), pointer :: GLLcx, GLLpolx, GLLwx
       real, dimension (:,:), pointer :: hprimex, hTprimex
       real, dimension (:), pointer :: GLLcz, GLLpolz, GLLwz
       real, dimension (:,:), pointer :: hprimez, hTprimez

       logical :: Filtering, Px, Pz, Left, Down

       character(len=1) :: material_type

    end type Subdomain

contains

    !>
    !! \brief
    !!
    !! \param type (Subdomain) S
    !<


    subroutine Lame_coefficients (S)

        type (Subdomain) :: S

        S%DMu = S%Sspeed**2 * S%Ddensity
        S%DLambda = (S%Pspeed**2 - 2 * S%Sspeed **2 ) * S%Ddensity

    end subroutine Lame_coefficients

end module ssubdomains
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
