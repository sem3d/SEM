!>
!!\file TimeP.f90
!!\brief
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module stimeparam

    ! Modified by Gaetano 10/02/2005

    type :: time

       logical :: acceleration_scheme, velocity_scheme
       integer :: ntimeMax, NtimeMin, nSnap, ntrace, ncheck
       real :: alpha, beta, gamma, duration, Time_snapshots, dtmin, rtime
       integer :: iter_reprise !Gsa Ipsis
       integer :: prot_m0, prot_m1, prot_m2 !Gsa Ipsis
    end type time

end module stimeparam
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
