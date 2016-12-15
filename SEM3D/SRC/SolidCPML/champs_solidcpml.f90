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

        ! Copy of node global coords : mandatory to compute distances in the PML
        real(fpp), allocatable, dimension(:,:) :: GlobCoord
        ! We keep separate variables for the cases where we have 1, 2 or 3 attenuation directions
        ! since we always have at least one, xxx_0 are indexed by ee,bnum
        real(fpp), allocatable, dimension(:,:,:,:,:) :: Alpha_0, dxi_k_0, Kappa_0 ! dxi_k = dxi/kappa
        ! for the other two cases we have much less elements (12*N and 8 in case of cube NxNxN)
        ! so we maintain two separate indirection indices I1/I2
        ! I1 = -1 means we have only one dir, I1>=0 and I2==-1 we have two directions, 3 otherwise
        ! D0 : index of first direction with pml!=0
        ! D1 : index of second direction with pml!=0
        ! There is no D2 because with 3 dirs!=0 we have D0=0 D1=1 D2=2
        integer,   allocatable, dimension(:,:) :: I1, I2, D0, D1 ! ee, bnum
        real(fpp), allocatable, dimension(:,:,:,:) :: Alpha_1, Kappa_1, dxi_k_1
        real(fpp), allocatable, dimension(:,:,:,:) :: Alpha_2, Kappa_2, dxi_k_2
        ! the number of elements with 2 (resp. 3) attenuation direction
        integer :: dir1_count, dir2_count

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champssolidpml), dimension(0:1) :: champs ! Etat courant

        ! Convolutional terms R
        ! First dimension : 0, 1, 2 <=> u, t*u, t^2*u
        ! Last  dimension : 0, 1, 2 <=> x,   y,     z
        ! Convolutional terms R from Ref1 for L and Lijk with only one direction of attenuation
        ! R1 = L*u,  R2=exp(-a0t)*du/dx
        ! Dimensions : R1_0(VS,3,N,N,N,NB) R2_0(VS,9,N,N,N,NB)
        ! Dimensions : R1_1(3,N,N,N,N1) R2_0(9,N,N,N,N1)
        ! Dimensions : R1_2(3,N,N,N,N2) R2_0(9,N,N,N,N2)
        ! with N=ngll, VS:vector size, NB : Nelements/VS
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: R1_0, R2_0
        real(fpp), dimension(:,:,:,:,:)  , allocatable :: R1_1, R2_1
        real(fpp), dimension(:,:,:,:,:)  , allocatable :: R1_2, R2_2
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: DUDVold, Uold

        ! Save forces contributions for snapshots
        real(fpp), dimension(:,:), allocatable :: FDump, FMasU, Fint

        ! Solid - Fluid coupling

        real(fpp), dimension(:,:), allocatable :: Kappa_SF, Alpha_SF, dxi_k_SF
        real(fpp), dimension(:,:), allocatable :: R2_SF

        ! CPML parameters
        real(fpp) :: cpml_c
        real(fpp) :: cpml_n
        real(fpp) :: cpml_rc
        real(fpp) :: cpml_kappa_0, cpml_kappa_1
        real(fpp) :: alphamax
        ! Integration Rxx
        real(fpp) :: dt
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
