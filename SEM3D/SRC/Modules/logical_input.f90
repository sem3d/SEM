module logical_input


type :: Logical_array 

logical :: save_trace,save_snapshots,save_energy,plot_grid,save_restart,run_restart
logical :: run_exec, run_debug, run_echo
logical :: any_source
logical :: super_object, Neumann, super_object_local_present, Neumann_local_present,Save_Surface
! solid-fluid
logical :: solid_fluid, all_fluid, SF_local_present
! MPML
logical  :: MPML
end type 

end module logical_input
