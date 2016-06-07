!>
!!\file champs_solidcpml.f90
!!\brief Contient la définition du type champs pour un domaine solide CPML
!!
!<

module champs_solidpml

    use constants
    implicit none

    type :: champssolidpml
        ! Inconnues du problème : déplacement, vitesse
        real(fpp), dimension(:,:), allocatable :: Depla
        real(fpp), dimension(:,:), allocatable :: Veloc
    end type champssolidpml

    !! ATTENTION: voir index.h en ce qui concerne les champs dont les noms commencent par m_
    type domain_solidpml
        ! D'abord, les données membres qui ne sont pas modifiées

        ! Nombre de gll dans chaque element du domaine
        integer :: ngll

        ! Nombre total de gll du domaine (assembles)
        integer :: nglltot

        ! Nombre d'elements dans le domaine
        integer :: nbelem
        integer :: nblocks

        ! Nombre d'elements alloues dans le domaine (>=nbelem)
        integer :: nbelem_alloc
        integer :: nb_chunks ! nbelem_alloc == nb_chunks*VCHUNK

        ! Points, poids de gauss et derivees
        real(fpp), dimension (:), allocatable :: GLLc
        real(fpp), dimension (:), allocatable :: GLLw
        real(fpp), dimension (:,:), allocatable :: hprime
        real(fpp), dimension (:,:), allocatable :: hTprime

        ! MassMat pour elements solide, fluide, solide pml et fluide pml
        real(fpp), dimension(:), allocatable :: MassMat

        ! PML is "like" an anisotropic material : do not use Lambda/Mu but use Cij
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_Density
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_Cij

        real(fpp), dimension(:,:,:,:,:),     allocatable :: m_Jacob
        real(fpp), dimension(:,:,:,:,:,:,:), allocatable :: m_InvGrad

        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! Index of a gll node within the physical domain
        integer, dimension (:,:,:,:,:), allocatable :: m_Idom ! Idom copied from element

        ! Copy of node global coords : mandatory to compute distances in the PML
        ! TODO : replace copy by handle (pointer) on Tdomain%GlobCoord ?
        !        pointer on allocatable seems to be a (huge) mess in Fortran : giving up !...
        !        a copy of node coord should work for a very first version.
        real(fpp), allocatable :: GlobCoord(:,:)

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champssolidpml) :: champs0
        type(champssolidpml) :: champs1
        real(fpp), dimension(:,:), allocatable :: Forces
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_R1 ! Convolutional term R1 (19a) from Ref1
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_R2 ! Convolutional term R2 (19b) from Ref1
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_R3 ! Convolutional term R3 (19c) from Ref1

        ! CPML parameters
        real(fpp) :: c_x, c_y, c_z
        integer   :: n_x, n_y, n_z
        real(fpp) :: r_c
        integer   :: kappa_0, kappa_1
        real(fpp) :: L_x, L_y, L_z
        real(fpp) :: bpp_x, bpp_y, bpp_z ! bpp : begin PML position
        real(fpp) :: alphamax

    end type domain_solidpml

    contains

end module champs_solidpml

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
