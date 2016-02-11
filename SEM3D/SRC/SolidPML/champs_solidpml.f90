!>
!!\file champs_solidpml.f90
!!\brief Contient la définition du type champs pour un domaine solide PML
!!
!<

module champs_solidpml

    use constants
    implicit none

    type :: champssolidpml

        !! Solide PML
        real(fpp), dimension(:,:,:), allocatable :: ForcesPML
        real(fpp), dimension(:,:,:), allocatable :: VelocPML
        real(fpp), dimension(:,:,:), allocatable :: DumpV

    end type champssolidpml

    type domain_solidpml
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

        real(fpp), dimension (:,:,:,:), allocatable :: m_Lambda, m_Mu, m_Kappa, m_Density

        real(fpp), dimension(:,:,:,:),     allocatable :: m_Jacob
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_InvGrad

        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champssolidpml) :: champs0
        type(champssolidpml) :: champs1
        real(fpp), dimension(:,:,:,:,:), allocatable :: Diagonal_Stress1, Diagonal_Stress2, Diagonal_Stress3
        real(fpp), dimension(:,:,:,:,:), allocatable :: Residual_Stress1, Residual_Stress2, Residual_Stress3
        real(fpp), dimension(:,:,:,:,:), allocatable :: Diagonal_Stress, Residual_Stress
        real(fpp), dimension(:,:,:,:,:), allocatable :: PMLDumpSx,PMLDumpSy,PMLDumpSz
    end type domain_solidpml

    contains

end module champs_solidpml

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
