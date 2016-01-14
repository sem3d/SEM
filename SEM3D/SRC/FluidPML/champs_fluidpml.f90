!>
!!\file champs_fluidpml.f90
!!\brief Contient la définition du type champs pour un domaine fluide PML
!!
!<

module champs_fluidpml

    use constants
    implicit none

    type :: champsfluidpml

        !! Fluide PML
        real(fpp), dimension(:,:), allocatable   :: fpml_Forces
        real(fpp), dimension(:,:), allocatable   :: fpml_Phi
        real(fpp), dimension(:,:), allocatable   :: fpml_VelPhi
        real(fpp), dimension(:,:,:), allocatable :: fpml_DumpV

    end type champsfluidpml

    contains

end module champs_fluidpml

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
