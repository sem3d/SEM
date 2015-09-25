!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file Edge.f90
!! \brief
!!
!<

module sedges
    implicit none
    type :: edge

        integer :: ngll
        integer :: domain
        integer, dimension (:), allocatable :: Iglobnum_Edge
        integer, dimension (:), allocatable :: Idom
        integer, dimension(0:1) :: inodes

        ! TODO REMOVE
        logical :: PML, Abs, FPML ! redondant avec mat_index
        ! Lien entre ngll et numérotation des champs globaux
        integer, dimension (:), allocatable :: Renum ! renommage en Igll
        ! solid-fluid
        logical  :: solid, fluid_dirich ! redondant avec mat_index


       !! Couplage Externe
!       real, dimension (:,:), allocatable :: ForcesExt
!       real, dimension (:), allocatable :: tsurfsem

    end type edge

contains

    ! ###########################################################
    subroutine init_edge(ed)
        type(Edge), intent(inout) :: ed
        !
        ed%ngll = 0
        ed%domain = -1
        !
        ed%PML = .false.
        ed%Abs = .false.
        ed%FPML = .false.
        ed%solid = .true.
    end subroutine init_edge

end module sedges

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
