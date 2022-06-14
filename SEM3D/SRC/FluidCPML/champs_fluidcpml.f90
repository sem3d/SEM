!>
!!\file champs_fluidpml.f90
!!\brief Contient la définition du type champs pour un domaine fluide
!!
!<

module champs_fluidpml

    use constants
    use ssubdomains
    use mdombase
    implicit none

    type champsfluidpml

        !! Fluide
        real(fpp), dimension(:), allocatable :: ForcesFl
        real(fpp), dimension(:), allocatable :: Phi
        real(fpp), dimension(:), allocatable :: VelPhi

    end type champsfluidpml

    type, extends(dombase_cpml) :: domain_fluidpml
        ! D'abord, les données membres qui ne sont pas modifiées

        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_IDensity ! Inverse of density

        ! Element materials
        type(subdomain), dimension (:), pointer :: sSubDomain ! Point to Tdomain%sSubDomain

        ! A partir de là, les données membres sont modifiées en cours de calcul

        ! Champs
        type(champsfluidpml), dimension(:), allocatable :: champs

        ! Masse pour elements solide cpml
        real(fpp), dimension(:), allocatable :: DumpMat ! Delta 1st derivative term in (12a) from Ref1
        real(fpp), dimension(:), allocatable :: MasUMat ! M^U <=> Delta term in (12a) from Ref1

        ! Convolutional terms R
        ! First dimension : 0, 1, 2 <=> u, t*u, t^2*u
        ! Last  dimension : 0, 1, 2 <=> x,   y,     z
        ! Convolutional terms R from Ref1 for L and Lijk with only one direction of attenuation
        ! R1 = L*u,  R2=exp(-a0t)*du/dx
        ! Dimensions : R1_0(VS,N,N,N,NB) R2_0(VS,3,N,N,N,NB)
        ! Dimensions : R1_1(N,N,N,N1) R2_0(3,N,N,N,N1)
        ! Dimensions : R1_2(N,N,N,N2) R2_0(3,N,N,N,N2)
        ! with N=ngll, VS:vector size, NB : Nelements/VS
        real(fpp), dimension(:,:,:,:,:)  , allocatable :: R1_0
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: R2_0
        real(fpp), dimension(:,:,:,:)    , allocatable :: R1_1
        real(fpp), dimension(:,:,:,:,:)  , allocatable :: R2_1
        real(fpp), dimension(:,:,:,:)    , allocatable :: R1_2
        real(fpp), dimension(:,:,:,:,:)  , allocatable :: R2_2

        ! Solid - Fluid coupling
        real(fpp), dimension(:), allocatable :: R_0_SF  ! exp(-a0.t)*X
        real(fpp), dimension(:), allocatable :: R_1_SF  ! exp(-a1.t)*X or (if a0=a1) t.exp(-a0.t)*X
    end type domain_fluidpml

    contains

end module champs_fluidpml

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
