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
    type SurfaceT
        type(surf_num) :: surf_sl
        type(surf_num) :: surf_fl
        type(surf_num) :: surf_spml
        type(surf_num) :: surf_fpml
        character(len=100) :: name
        integer            :: domain
        integer :: cond_type ! from constants.F90 COND_*
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
