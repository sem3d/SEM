!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

!>
!!\file champs_fluid_aniso.f90
!!\brief Type definitions for the anisotropic-density acoustic domain.
!!
!! Physical model (Capdeville & Cance 2015): the anisotropy is carried by the
!! DENSITY, not by the bulk modulus. Velocity-potential formulation, mirroring
!! the regular fluid domain (champs_fluid.f90):
!!
!!     (1/kappa) phi_tt = d_i( rho^{-1}_ij d_j phi ) ,   v_i = rho^{-1}_ij d_j phi
!!
!! - scalar bulk modulus kappa (m_Lambda) in the inertial/mass term,
!! - symmetric inverse-density tensor rho^{-1}_ij (m_IDensTensor, 6 indep comps)
!!   generalising the regular fluid scalar IDensity = 1/rho.
!! Primary unknown: velocity potential phi (champs%Phi).
!<

module champs_fluid_aniso

    use constants
    use mdombase
    implicit none

    type champsfluid_aniso
        real(fpp), dimension(:), allocatable :: ForcesFl  ! rhs of the weak form
        real(fpp), dimension(:), allocatable :: Phi       ! velocity potential (primary unknown)
        real(fpp), dimension(:), allocatable :: VelPhi    ! d(phi)/dt
    end type champsfluid_aniso

    type, extends(dombase) :: domain_fluid_aniso
        ! Inverse-density tensor rho^{-1}_ij (symmetric, 6 independent comps):
        !   index 0 = 11, 1 = 22, 2 = 33, 3 = 12, 4 = 13, 5 = 23
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: m_IDensTensor
        ! Scalar bulk modulus kappa
        real(fpp), dimension(:,:,:,:,:),   allocatable :: m_Lambda

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
