!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Element.f90
!!\brief contient les méthodes qui assure la gestion du type Element.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module selement
    use deriv3d
    use blas
    implicit none

    type :: element
        integer :: mat_index
        integer :: lnum ! local number of element within its domain
        integer :: domain ! Type de domaine, voir constants : DOM_SOLID, DOM_FLUID, ...
        real, dimension (:,:,:), allocatable :: MassMat

        ! Whether this element will be part of snapshot outputs
        logical :: OUTPUT

        ! These should not be used during the simulation, only at init time
        integer, dimension (:), allocatable :: Control_nodes
        integer, dimension (0:5) :: Near_Faces
        integer, dimension (0:11) :: Near_Edges
        integer, dimension (0:7) :: Near_Vertices

        ! Index of a gll node within the global nodes array
        integer, dimension (:,:,:), allocatable :: Iglobnum
        ! Index of a gll node within the physical domain
        integer, dimension (:,:,:), allocatable :: Idom ! deallocated in define_array + copied into domain_XXX

        ! Integrated Energy inside the Element
        double precision :: En_S_int
        double precision :: En_P_int
        ! Averaged Energy inside the subElements (defined by GLLs division)
        double precision, dimension (:,:,:), allocatable :: En_S_avg
        double precision, dimension (:,:,:), allocatable :: En_P_avg

    end type element

contains

    subroutine init_element(el)
        type(element), intent(inout) :: el

        el%mat_index=-1
        el%domain = -1
        el%En_S_int = -1
        el%En_P_int = -1
    end subroutine init_element

end module selement


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
