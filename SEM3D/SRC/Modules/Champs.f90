!>
!!\file Champs.f90
!!\brief Contient le d�finition du type champs
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
        real, dimension(:,:,:), allocatable :: ForcesPML
        real, dimension(:,:,:), allocatable :: VelocPML
        real, dimension(:,:,:), allocatable :: DumpV

        !! Fluide
        real, dimension(:), allocatable :: ForcesFl
        real, dimension(:), allocatable :: Phi
        real, dimension(:), allocatable :: VelPhi

        !! Fluide PML
        real, dimension(:,:), allocatable   :: fpml_Forces
        real, dimension(:,:), allocatable   :: fpml_Phi
        real, dimension(:,:), allocatable   :: fpml_VelPhi
        real, dimension(:,:,:), allocatable :: fpml_DumpV
!        !! Non Linear
!        real, dimension(:,:), allocatable :: Stress
!        real, dimension(:,:), allocatable :: Epsilon_pl
!        real, dimension(:,:), allocatable :: Xkin
!        real, dimension(:),   allocatable :: Riso
!        integer, dimension(:), allocatable :: element_connectivity
    end type champs

end module schamps
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
