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
    type :: face
        integer :: ngll1, ngll2
        integer :: domain
        !
        integer :: elem
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
        real, dimension (:,:,:), allocatable :: Forces, Forces3
        logical                              :: PML
        !! Couplage Externe
!       real, dimension (:,:,:), allocatable :: ForcesExt
!       real, dimension (:,:), allocatable :: tsurfsem
    end type face

contains
    ! ############################################################

    subroutine init_face(fc)
        type(Face), intent(inout) :: fc
        fc%domain = -1
        fc%ngll1 = 0
        fc%ngll2 = 0
    end subroutine init_face

    subroutine allocate_face_force(fc)
     
     implicit none
     type(Face), intent(inout) :: fc

        allocate(fc%Forces(1:fc%ngll1-2,1:fc%ngll2-2,0:2))
        allocate(fc%Forces3(1:fc%ngll1-2,1:fc%ngll2-2,0:2))
        fc%Forces = 0.0
        fc%Forces3= 0.0

    end subroutine allocate_face_force

   subroutine free_face_force(fc)
    
     implicit none
     type(Face), intent(inout) :: fc
     
       deallocate(fc%Forces)
       deallocate(fc%Forces3)
   end subroutine free_face_force

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
