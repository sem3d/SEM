!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Fault.f90
!!\brief
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module sfault
    implicit none
    type Elem_fault

       integer :: nUp,nDown,nFleft
       logical :: CommLeft
       real :: X0,X1
       real, dimension (:), pointer :: Xco
       real, dimension (:), pointer :: Cgammam1, Zeta, ZetaHF,Bt,ds,ReactionV
       real, dimension (:,:) , pointer :: Reaction, ReactionBF,Dislop, DislopBF
       real, dimension (:), pointer :: Dislo
       real, dimension (:), pointer :: Second_or,IdisloHF,A_sec,DisloHF

    end type Elem_fault


    type Faul

       integer :: nelem
       real :: tau0,tauf,vr,dx,Xnucl,deltaTau,tauu,Dc
       type (Elem_fault), dimension (:), pointer :: specel1

    end type Faul

end module sfault

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
