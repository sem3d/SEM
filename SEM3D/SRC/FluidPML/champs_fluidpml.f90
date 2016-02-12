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

    type domain_fluidpml
        ! D'abord, les données membres qui ne sont pas modifiées

        ! Nombre de gll dans chaque element du domaine
        integer :: ngllx
        integer :: nglly
        integer :: ngllz

        ! Nombre total de gll du domaine (assembles)
        integer :: ngll

        ! Nombre d'elements dans le domaine
        integer :: nbelem

        ! MassMat pour elements solide, fluide, solide pml et fluide pml
        real(fpp), dimension(:), allocatable :: MassMat
        real(fpp), dimension(:,:), allocatable :: DumpMass

        real(fpp), dimension (:,:,:,:), allocatable :: m_Lambda, m_Density

        real(fpp), dimension(:,:,:,:),     allocatable :: m_Jacob
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_InvGrad

        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! Index of a gll node within the physical domain
        integer, dimension (:,:,:,:), allocatable :: m_Idom ! Idom copied from element

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champsfluidpml) :: champs0
        type(champsfluidpml) :: champs1
        real(fpp), dimension(:,:,:,:,:), allocatable :: Veloc
        real(fpp), dimension(:,:,:,:,:), allocatable :: PMLDumpSx,PMLDumpSy,PMLDumpSz
    end type domain_fluidpml

    contains

end module champs_fluidpml

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
