!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Face.f90
!!\brief Gère les faces des éléments.
!!
!<

module sfaces
    implicit none
    ! Faces are unique between elements of the same domain.
    ! If two elements of a different domain share a geometric face
    ! then two faces exists.
    type :: face
        integer :: ngll
        integer :: domain
        integer :: lnum ! domain-local face number
        !
        integer :: elem, elem_0, elem_1
        ! True if this face doesn't have an associated element on this cpu
        logical :: orphan
        ! Index dans Tdomain%GlobalCoord des coordonnees du pt de gauss
        ! Long term : on peut avoir un champ GlobalCoord par domaine et virer les Iglobnum_*
        integer, dimension (:,:), allocatable :: Iglobnum_Face
        ! Index dans un domaine des valeurs associées au gll
        integer, dimension (:,:), allocatable :: Idom
        ! Index dans Tdomain%Coord_nodes des 4 sommets de la face (y compris pour les Hex27)
        integer, dimension(0:3) :: inodes
        !
        logical                              :: PML
    end type face

contains
    ! ############################################################

    subroutine init_face(fc)
        type(Face), intent(inout) :: fc
        fc%domain = -1
        fc%ngll = 0
    end subroutine init_face

end module sfaces

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
