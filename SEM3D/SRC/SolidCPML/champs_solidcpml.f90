!>
!!\file champs_solidcpml.f90
!!\brief Contient la définition du type champs pour un domaine solide CPML
!!
!<

module champs_solidpml

    use constants
    use ssubdomains
    use selement
    use mdombase
    implicit none

    ! Inconnues du problème : déplacement, vitesse
    type :: champssolidpml
        ! Le déplacement et la vitesse sont "staggered" (leap frog) : ceci est imposé par Solid/Fluid
        real(fpp), dimension(:,:), allocatable :: Depla ! U_n+3/2
        real(fpp), dimension(:,:), allocatable :: Veloc ! V_n+1
        real(fpp), dimension(:,:), allocatable :: Forces
    end type champssolidpml

    !! ATTENTION: voir index.h en ce qui concerne les champs dont les noms commencent par m_
    type, extends(dombase) :: domain_solidpml
        ! D'abord, les données membres qui ne sont pas modifiées

        ! Nombre d'elements alloues dans le domaine (>=nbelem)
        integer :: nbelem_alloc
        integer :: nb_chunks ! nbelem_alloc == nb_chunks*VCHUNK

        ! Materiaux
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda, m_Mu, m_Density

        ! Masse pour elements solide cpml
        real(fpp), dimension(:), allocatable :: DumpMat ! Delta 1st derivative term in (12a) from Ref1
        real(fpp), dimension(:), allocatable :: MasUMat ! M^U <=> Delta term in (12a) from Ref1

        ! Element materials
        type(subdomain), dimension (:), pointer :: sSubDomain ! Point to Tdomain%sSubDomain

        ! Material index
        integer, dimension(:,:), allocatable :: mat_index

        ! Copy of node global coords : mandatory to compute distances in the PML
        real(fpp), allocatable, dimension(:,:) :: GlobCoord
        real(fpp), allocatable, dimension(:,:,:,:,:,:) :: Kappa

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champssolidpml) :: champs0 ! Etat courant
        type(champssolidpml) :: champs1 ! Prediction (à partir de l'état courant)

        real(fpp), dimension(:,:), allocatable :: R1 ! Convolutional term R1 (19a) from Ref1
        real(fpp), dimension(:,:), allocatable :: R2 ! Convolutional term R2 (19b) from Ref1
        real(fpp), dimension(:,:), allocatable :: R3 ! Convolutional term R3 (19c) from Ref1

        ! CPML parameters
        real(fpp) :: c(0:2)
        integer   :: n(0:2)
        real(fpp) :: rc
        real(fpp) :: kappa_0, kappa_1
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
