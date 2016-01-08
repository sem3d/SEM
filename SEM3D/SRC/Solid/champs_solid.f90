!>
!!\file champs_solid.f90
!!\brief Contient la définition du type champs pour un domaine solide
!!
!<

module champs_solid

    implicit none

    type :: champssolid

        !! Solide
        real, dimension(:,:), allocatable :: Forces
        real, dimension(:,:), allocatable :: Depla
        real, dimension(:,:), allocatable :: Veloc

    end type champssolid

    contains

end module champs_solid

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
