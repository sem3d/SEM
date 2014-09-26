!>
!!\file Champs.f90
!!\brief Contient le définition du type champs
!!
!<

module schamps
    implicit none

    type :: champs

        !! Solide
        real, dimension(:,:), allocatable :: Forces
        real, dimension(:,:), allocatable :: Depla !! V0
        real, dimension(:,:), allocatable :: Veloc

        !! Solide PML
        real, dimension(:,:), allocatable :: ForcesPML
        real, dimension(:,:), allocatable :: VelocPML
        real, dimension(:,:), allocatable :: DumpV

        !! Fluide
        real, dimension(:), allocatable :: ForcesFl
        real, dimension(:), allocatable :: Phi
        real, dimension(:), allocatable :: VelPhi

        ! Champs rempli de 0. et de 1. pour annuler VelPhi sur faces libres
        real, dimension(:), allocatable :: Fluid_dirich

        !! Fluide PML

        !! Couplage solide / fluide
        real, dimension(:,:), allocatable :: Save_forces
        real, dimension(:,:), allocatable :: Save_depla

    end type champs

end module schamps
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
