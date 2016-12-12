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
    type, extends(dombase_cpml) :: domain_solidpml
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

        ! Integration Rxx
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
