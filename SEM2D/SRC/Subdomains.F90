!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Subdomains.F90
!!\brief Calcul des coefficients de Lame.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssubdomains

    use constants

    type Subdomain

       integer :: NGLLx, NGLLz,  n_loc_dim, wpml, npow, pml_type
       integer :: type_DG,type_Flux

       real(fpp) :: Pspeed, Sspeed, Ddensity, DT, DLambda, DMu, Apow, freq, k
       real(fpp), dimension (:), pointer :: GLLcx, GLLpolx, GLLwx
       real(fpp), dimension (:,:), pointer :: hprimex, hTprimex
       real(fpp), dimension (:), pointer :: GLLcz, GLLpolz, GLLwz
       real(fpp), dimension (:,:), pointer :: hprimez, hTprimez

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

        ! Check Case Fluid (Acoustic)
        if ((S%material_type == "F") .AND. (S%Sspeed .NE. 0.)) then
           STOP "Error in material input file : Fluids must have S-speed equal to 0. "
        endif

    end subroutine Lame_coefficients

end module ssubdomains

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
