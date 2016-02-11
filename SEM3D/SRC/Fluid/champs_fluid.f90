!>
!!\file champs_fluid.f90
!!\brief Contient la définition du type champs pour un domaine fluide
!!
!<

module champs_fluid

    use constants
    implicit none

    type champsfluid

        !! Fluide
        real(fpp), dimension(:), allocatable :: ForcesFl
        real(fpp), dimension(:), allocatable :: Phi
        real(fpp), dimension(:), allocatable :: VelPhi

    end type champsfluid

    type domain_fluid
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

        real(fpp), dimension (:,:,:,:), allocatable :: m_Lambda, m_Mu, m_Kappa, m_Density

        real(fpp), dimension(:,:,:,:),     allocatable :: m_Jacob
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_InvGrad

        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champsfluid) :: champs0
        type(champsfluid) :: champs1
    end type domain_fluid

    contains

end module champs_fluid

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
