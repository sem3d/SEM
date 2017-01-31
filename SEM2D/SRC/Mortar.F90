!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Face.F90
!!\brief Gere les faces des elements.
!!\author S. TERRANA
!!\version 1.0
!!\date 01/05/2016
!!
!<

module smortars
    use constants

    implicit none

    type :: mortar

       integer :: ngllmin, ngllmax
       integer, dimension (0:1) :: Near_Face
       integer, dimension (0:1) :: Near_Element

       real(fpp), dimension (:), allocatable   :: Coeff_Integr
       real(fpp), dimension (:,:), allocatable :: MatReinterp, MatProj



    end type mortar

contains

end module smortars

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
