!>
!!\file champs_solid_dg.f90
!!\brief Contient la définition du type champs pour un domaine solide
!!
!<

module champs_solid_dg

    use constants
    use mdombase
    implicit none

    type :: champssolid_dg
        real(fpp), dimension(:,:,:,:,:,:), allocatable :: Q       ! VCH*N*N*N*9*(NE/VCH)       
        real(fpp), dimension(:,:,:,:,:),   allocatable :: trace_Q ! VCH*N*N*9*(NF/VCH)
    end type champssolid_dg

    !! ATTENTION: voir index.h en ce qui concerne les champs dont les noms commencent par m_
    type, extends(dombase) :: domain_solid_dg
        ! D'abord, les données membres qui ne sont pas modifiées
        logical :: aniso
        real(fpp), dimension (:,:,:,:,:), allocatable :: m_Lambda, m_Mu, m_Density

        ! Champs
        type(champssolid_dg), dimension(:), allocatable :: champs
        real(fpp), dimension(:,:), allocatable :: Qasm ! NN*9        
        
        ! VCHUNK,F,N,EB) : N: 0:
        ! N:
        ! 0: domain-local face number,
        ! 1-3 : I0
        ! 4-6: DI
        ! 7-9: DJ  face(i,j) = elem(I0+i*DI+j*DJ)
        ! 10: side: 0 or 1
        integer, dimension (:,:,:,:), allocatable :: m_Itrace ! Itrace copied from element

    end type domain_solid_dg

    contains

end module champs_solid_dg

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
