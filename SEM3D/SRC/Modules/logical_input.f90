!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file logical_input.f90
!!\brief Contient la definition du type Logical_array.
!!
!<

module logical_input
    implicit none
    type :: Logical_array

       logical :: save_trace,save_snapshots,save_restart,run_restart
       logical :: any_source
       logical :: Neumann, Neumann_local_present
       logical :: surfBC
       ! solid-fluid
       logical :: SF_local_present
       ! MPML
       logical :: MPML
    end type Logical_array

end module logical_input

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
