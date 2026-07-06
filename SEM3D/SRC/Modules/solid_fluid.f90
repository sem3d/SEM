!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module solid_fluid
    use constants, only : fpp
    use sinterface
    implicit none

    ! in general, when refering to both sides of faces, edges or vertices:
    !       index 0 for fluid part, 1 for the solid one


    ! general SF object
    type :: SF_object
        type(inter_num) :: intSolFlu  ! 0: Solid 1: fluid
        type(inter_num) :: intSolFluPml ! 0: spml 1: fpml
        real(fpp), allocatable, dimension(:,:) :: SF_BtN
        real(fpp), allocatable, dimension(:,:) :: SFPml_BtN
    end type SF_object

end module solid_fluid

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
