!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Domain.f90
!!\brief Contient le d�finition du type domain
!!
!<

module mdombase
    use constants
    implicit none
#include "index.h"

    ! Ce type regroupe les donnees communes a tous les domaines physiques (sol/flu/pml)
    type :: dombase
        ! Nombre de gll dans chaque element du domaine
        integer :: ngll

        ! Nombre total de gll du domaine (assembles)
        integer :: nglltot

        ! Nombre d'elements dans le domaine
        integer :: nbelem

        ! Nombre d'elements alloues dans le domaine (>=nbelem)
        integer :: nblocks ! nbelem_alloc == nblocks*VCHUNK

        ! Pas de temps d'integration (normalement le meme pour tout domaine)
        real(fpp) :: dt

        ! Points, poids de gauss et derivees
        real(fpp), dimension (:), allocatable :: GLLc
        real(fpp), dimension (:), allocatable :: GLLw
        real(fpp), dimension (:,:), allocatable :: hprime
        real(fpp), dimension (:,:), allocatable :: hTprime

        ! MassMat pour elements solide, fluide, solide pml et fluide pml
        real(fpp), dimension(:), allocatable :: MassMat

        real(fpp), dimension(:,:,:,:,:),     allocatable :: m_Jacob
        real(fpp), dimension(:,:,:,:,:,:,:), allocatable :: m_InvGrad

        ! Condition de dirichlet : liste des noeuds à mettre à 0 pour chaque domaine
        integer :: n_dirich
        integer, dimension(:), allocatable :: dirich

        ! Index of a gll node within the physical domain
        integer, dimension (:,:,:,:,:), allocatable :: m_Idom ! Idom copied from element
    end type dombase

    type, extends(dombase) :: dombase_cpml
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
        ! CPML parameters
        real(fpp) :: cpml_c
        real(fpp) :: cpml_n
        real(fpp) :: cpml_rc
        real(fpp) :: cpml_kappa_0, cpml_kappa_1
        real(fpp) :: alphamax

        ! Solid - Fluid coupling
        real(fpp), dimension(:,:), allocatable :: Kappa_SF, Alpha_SF, dxi_k_SF
        integer, dimension(:), allocatable :: D0_SF, D1_SF
    end type dombase_cpml
contains

    subroutine init_dombase(bz)
        use gll3d
        class(dombase), intent(inout) :: bz
        !
        integer :: nbelem, ngll, nblocks
        !
        ngll   = bz%ngll
        nbelem = bz%nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call compute_gll_data(ngll, bz%gllc, bz%gllw, bz%hprime, bz%htprime)

        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)
            nblocks = ((nbelem+VCHUNK-1)/VCHUNK)
            bz%nblocks = nblocks

            allocate (bz%Jacob_  (        0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate (bz%InvGrad_(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate (bz%Idom_(0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            ! Initialize to point after the last gll
            ! Since we allocate nblock*VCHUNK elements we could have up to VCHUNK-1 extra (fake)
            ! elements those will point to this extra gll so we can ignore them in most
            ! vectorized loops
            bz%m_Idom = bz%nglltot
        end if
        if (bz%nglltot /= 0) then
            ! Allocation de MassMat
            allocate(bz%MassMat(0:bz%nglltot))
            bz%MassMat = 0d0
            bz%MassMat(bz%nglltot) = 1d0
        end if
    end subroutine init_dombase

    subroutine deallocate_dombase(bz)
        class(dombase), intent(inout) :: bz
        !
        if(allocated(bz%m_Jacob  )) deallocate(bz%m_Jacob  )
        if(allocated(bz%m_InvGrad)) deallocate(bz%m_InvGrad)
        if(allocated(bz%m_Idom))    deallocate(bz%m_Idom)
        if(allocated(bz%gllc))    deallocate(bz%gllc)
        if(allocated(bz%gllw))    deallocate(bz%gllw)
        if(allocated(bz%hprime))  deallocate(bz%hprime)
        if(allocated(bz%htprime)) deallocate(bz%htprime)
        if(allocated(bz%MassMat)) deallocate(bz%MassMat)
    end subroutine deallocate_dombase

    subroutine allocate_dombase_cpml(dom)
        class(dombase_cpml) :: dom
    end subroutine allocate_dombase_cpml

    subroutine deallocate_dombase_cpml(bz)
        class(dombase_cpml), intent(inout) :: bz
        !
        if(allocated(bz%Alpha_0)) deallocate(bz%Alpha_0)
        if(allocated(bz%dxi_k_0)) deallocate(bz%dxi_k_0)
        if(allocated(bz%Kappa_0)) deallocate(bz%Kappa_0)

        if(allocated(bz%Alpha_1)) deallocate(bz%Alpha_1)
        if(allocated(bz%dxi_k_1)) deallocate(bz%dxi_k_1)
        if(allocated(bz%Kappa_1)) deallocate(bz%Kappa_1)

        if(allocated(bz%Alpha_2)) deallocate(bz%Alpha_2)
        if(allocated(bz%dxi_k_2)) deallocate(bz%dxi_k_2)
        if(allocated(bz%Kappa_2)) deallocate(bz%Kappa_2)

        if(allocated(bz%I1)) deallocate(bz%I1)
        if(allocated(bz%I2)) deallocate(bz%I2)
        if(allocated(bz%D0)) deallocate(bz%D0)
        if(allocated(bz%D1)) deallocate(bz%D1)

        if(allocated(bz%GlobCoord)) deallocate(bz%GlobCoord)
    end subroutine deallocate_dombase_cpml
end module mdombase

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
