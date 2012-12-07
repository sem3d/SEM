module sfault

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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
