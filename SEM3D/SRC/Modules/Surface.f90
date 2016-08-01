!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file Surface.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

module ssurf

    use sinterface
    
    type elastic_
       real(kind=8)                       :: Mu, Lambda
       real(kind=8)                       :: Sspeed, Pspeed
       real(kind=8)                       :: PWspeed, density
       integer                            :: mat_index
    end type elastic_

    type SurfaceParam
        real(kind=8)                       :: f0, amplitude, Rickertau, size, Speed
        real(kind=8)                       :: dir(0:2)
        real(kind=8)                       :: Kdir(0:2)
        real(kind=8)                       :: scoord(0:2)
        character                          :: wtype
        character(len=2)                   :: what_bc
        character(len=1500)                :: funcx, funcy, funcz
        character(len=1500)                :: funcxy, funcxz, funcyz
        character(len=12)                  :: varia
        character                          :: source
        integer                            :: dim, mat_index, shape, wave_type
        integer, allocatable               :: index(:), indexes(:)
        real(kind=8), allocatable          :: paravalue(:)
        character(len=2), dimension(1:100) :: paramname
        integer                            :: nparamvar, paramvar
    end type SurfaceParam

    type SurfaceT
        type(surf_num) :: surf_sl
        type(surf_num) :: surf_fl
        type(surf_num) :: surf_spml
        type(surf_num) :: surf_fpml
        type(elastic_) :: Elastic
        character(len=100) :: name
        integer            :: domain
        integer :: cond_type ! from constants.F90 COND_*
        real(kind=8), dimension(:,:), allocatable :: Surf_BtN
        real(kind=8), dimension(:,:), allocatable :: coord
        real(kind=8), dimension(:)  , allocatable :: source
    end type SurfaceT

end module ssurf

module sbassin
    implicit none

    type :: Bassin
       !       real     :: ymin,ymax
       !     n_layer nombre de couches
       !     n_colonne nombre de colonnes en x ici uniquement
       !     x_type == 0 on remet des materiaux  homogenes dans chaque bloc
       !     x_type == 1 on met des gradients pour chaque colonne en interpolant suivant z
       !     x_type == 2 on met des gradients pour chaque colonne en interpolant suivant z
       !      uniquement dans les couches entre zmin et zmax
       integer  :: n_colonne, n_layer, x_type
       real     :: zmin,zmax
       !    x_coord correspond aux abscisses des colonnes
       real, pointer, dimension(:) :: x_coord
       !      z_layer profondeur de  linterface pour chaque x de colonne
       !      on definit egalement le materiaux par rho, Cp , Cs
       real, pointer, dimension(:,:) :: z_layer, z_rho, z_Cp, z_Cs
    end type Bassin


end module sbassin

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
