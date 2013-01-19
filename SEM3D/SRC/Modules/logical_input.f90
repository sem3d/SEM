!>
!!\file logical_input.f90
!!\brief Contient la définition du type Logical_array.
!!
!<

module logical_input

    type :: Logical_array

       logical :: save_trace,save_snapshots,save_energy,plot_grid,save_restart,run_restart
       logical :: run_exec, run_debug, run_echo
       logical :: any_source
       logical :: super_object, Neumann, super_object_local_present, Neumann_local_present,Save_Surface
       ! solid-fluid
       logical :: solid_fluid, all_fluid, SF_local_present
       ! MPML
       logical :: MPML
       logical :: grad_bassin
    end type Logical_array

end module logical_input
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
