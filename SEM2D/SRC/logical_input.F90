!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file logical_input.F90
!!\brief Contient la definition du type Logical_array.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module logical_input

    ! Modified by Gaetano 07/04/05

    type :: Logical_array
       logical :: save_trace,save_snapshots,save_deformation, save_fault_trace
       logical :: run_restart, save_restart
       logical :: any_source, super_object, super_object_local_present
       logical :: compEnerg, Lamb_test, post_proc
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
