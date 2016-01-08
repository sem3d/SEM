!>
!!\file champs_fluid.f90
!!\brief Contient la définition du type champs pour un domaine fluide
!!
!<

module champs_fluid

    implicit none

    type champsfluid

        !! Fluide
        real, dimension(:), allocatable :: ForcesFl
        real, dimension(:), allocatable :: Phi
        real, dimension(:), allocatable :: VelPhi

    end type champsfluid

    contains

end module champs_fluid

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
