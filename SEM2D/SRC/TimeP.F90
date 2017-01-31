!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file TimeP.F90
!!\brief
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module stimeparam

    use constants

    type :: time
       integer :: ntime
       integer :: ntimeMax,nsnap,ncheck,NtimeMin,iter_reprise
       integer :: prot_m0,prot_m1,prot_m2
       logical :: acceleration_scheme, velocity_scheme
       real(fpp) :: alpha, beta, gamma
       real(fpp) :: duration
       real(fpp) :: dtmin,rtime,Time_snapshots
       real(fpp) :: courant
    end type time

end module stimeparam

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
