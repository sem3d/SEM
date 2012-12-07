module ssubdomains

    type Subdomain

       logical :: Filtering, Px, Py, Pz, Left, Forward, Down

       integer :: NGLLx, NGLLy, NGLLz, wpml, npow

       real :: Pspeed, Sspeed, Ddensity, Dt, Apow, freq, DLambda, DMu
       real, dimension (:), pointer :: GLLcx, GLLpolx, GLLwx
       real, dimension (:,:), pointer :: hprimex, hTprimex
       real, dimension (:), pointer :: GLLcy, GLLpoly, GLLwy
       real, dimension (:,:), pointer :: hprimey, hTprimey
       real, dimension (:), pointer :: GLLcz, GLLpolz, GLLwz
       real, dimension (:,:), pointer :: hprimez, hTprimez

       character(len=1) :: material_type

    end type Subdomain

contains

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
