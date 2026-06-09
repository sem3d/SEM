!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

!>
!!\file champs_fluid_aniso.f90
!!\brief Type definitions for the anisotropic acoustic domain.
!!
!! Physical model: scalar density rho, anisotropic bulk modulus tensor Kij (3x3 symmetric).
!! Governing equation: rho d2p/dt2 = d/dxi( Kij dp/dxj )
!! Primary field variable: pressure p.
!<

module champs_fluid_aniso

    use constants
    use mdombase
    implicit none

    type champsfluid_aniso
        real(fpp), dimension(:), allocatable   :: ForcesP  ! force residual (rhs of weak form)
        real(fpp), dimension(:), allocatable   :: P        ! pressure (primary unknown)
        real(fpp), dimension(:), allocatable   :: VelP     ! dp/dt
        real(fpp), dimension(:,:), allocatable :: Vel      ! particle velocity (0:2, 0:nglltot)
    end type champsfluid_aniso

    type, extends(dombase) :: domain_fluid_aniso
        ! 6 independent components of the symmetric acoustic stiffness tensor Kij:
        !   index 0 = K11, 1 = K22, 2 = K33, 3 = K12, 4 = K13, 5 = K23
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_Kij
        ! Scalar density
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_Rho

        ! Time fields (indexed 0..nsubsteps)
        type(champsfluid_aniso), dimension(:), allocatable :: champs
    end type domain_fluid_aniso

contains

end module champs_fluid_aniso

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
!! vim: set sw=4 ts=8 et tw=80 smartindent :
