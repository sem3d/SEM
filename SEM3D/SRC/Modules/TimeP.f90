!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file TimeP.f90
!!\brief
!!
!<
module stimeparam
    use constants
    implicit none
    type :: time

       logical :: acceleration_scheme, velocity_scheme
       integer :: ntimeMax, NtimeMin, nSnap, ntrace, ncheck
       real(fpp) :: alpha, beta, gamma, duration, Time_snapshots, dtmin, rtime
       integer :: iter_reprise
       integer :: type_timeinteg
       integer :: prot_m0, prot_m1, prot_m2
       real(fpp) :: courant, fmax
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
