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
    use constants, only : fpp
    use sinterface
    
    type elastic_
       real(fpp)                       :: Mu, Lambda
       real(fpp)                       :: Sspeed, Pspeed
       real(fpp)                       :: PWspeed, density
       integer                         :: mat_index
    end type elastic_

    type SurfaceParam
        real(fpp)                          :: f0, amplitude, Rickertau, size, Speed
        real(fpp)                          :: dir(0:2)
        real(fpp)                          :: Kdir(0:2)
        real(fpp)                          :: scoord(0:2)
        character                          :: wtype
        character(len=2)                   :: what_bc
        character(len=1500)                :: funcx, funcy, funcz
        character(len=1500)                :: funcxy, funcxz, funcyz
        character(len=12)                  :: varia
        character                          :: source
        integer                            :: dim, mat_index, shape, wave_type
        integer, allocatable               :: index(:), indexes(:)
        real(fpp), allocatable             :: paravalue(:)
        character(len=2), dimension(1:100) :: paramname
        integer                            :: nparamvar, paramvar
    end type SurfaceParam

    type SurfaceT
        type(surf_num) :: surf_sl
        type(surf_num) :: surf_fl
        type(surf_num) :: surf_spml
        type(surf_num) :: surf_fpml
        type(surf_num) :: surf_sldg
        type(elastic_) :: Elastic
        character(len=100) :: name
        integer            :: domain
        integer :: cond_type ! from constants.F90 COND_*
        real(fpp), dimension(:,:), allocatable :: Surf_BtN
        real(fpp), dimension(:,:), allocatable :: coord
        real(fpp), dimension(:)  , allocatable :: source
    end type SurfaceT

end module ssurf

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
