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
            bz%m_Idom = 0
        end if
        if (bz%nglltot /= 0) then
            ! Allocation de MassMat
            allocate(bz%MassMat(0:bz%nglltot-1))
            bz%MassMat = 0d0
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
