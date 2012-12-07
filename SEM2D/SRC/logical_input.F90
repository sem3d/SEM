!>
!!\file logical_input.F90
!!\brief Contient la définition du type Logical_array.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module logical_input

    ! Modified by Gaetano 07/04/05

    type :: Logical_array
       logical :: save_trace,save_snapshots,save_energy,plot_grid, save_deformation, save_fault_trace
       logical :: run_exec, run_debug, run_echo
       logical :: run_restart, save_restart
       logical :: any_source, super_object, super_object_local_present
    end type Logical_array


end module logical_input
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
