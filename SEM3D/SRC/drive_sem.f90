program  drive_sem


use sdomain

implicit none

include 'mpif.h'

type(domain), target :: Tdomain

integer :: code, rg, nb_procs, ntime, i_snap, n,i,j, icount, icountc
character (len=100) :: fnamef


call mpi_init(code)
call mpi_comm_rank(mpi_comm_world, rg, code)
call mpi_comm_size(mpi_comm_world, nb_procs, code)

! ##############  Begin of the program  #######################
write (*,*) "Define mesh properties ",rg
call read_input (Tdomain,rg)

if (Tdomain%logicD%super_object_local_present) then
 if (Tdomain%super_object_type == "P") then
   write(*,*) "Define Plane Wave properties",rg
   call define_planew_properties (Tdomain)
 endif
endif
if (Tdomain%logicD%neumann_local_present) then
  write(*,*) "Define Neumann properties",rg
  call define_neu_properties (Tdomain)
endif

write (*,*) "Compute Gauss-Lobatto-Legendre weights and zeroes ",rg
call compute_GLL (Tdomain)

write (*,*) "Define a global numbering for the collocation points ",rg
call global_numbering (Tdomain)

write (*,*) "Computing shape functions within their derivatives ",rg
if (Tdomain%n_nodes == 8) then
      call shape8(TDomain)   ! Linear interpolation
else if (Tdomain%n_nodes == 27) then
      write (*,*) "Sorry, 27 control points not yet implemented in the code... Wait for an upgrade ",rg
      stop
else
      write (*,*) "Bad number of nodes for hexaedral shape ",rg
      stop
endif

write (*,*) "Compute Courant parameter ",rg
call compute_Courant (Tdomain)

if (Tdomain%any_PML)   then
    write (*,*) "Attribute PML properties ",rg
    call PML_definition (Tdomain) 
endif

if (Tdomain%logicD%any_source) then
    write (*,*) "Computing point-source parameters ",rg
    call SourcePosition(Tdomain,nb_procs,rg)
endif
if (Tdomain%logicD%save_trace) then
    write (*,*) "Computing receiver locations ",rg
    call ReceiverExactPosition(Tdomain,rg)
endif

write (*,*) "Allocate fields ",rg
call allocate_domain (Tdomain)

write (*,*) "Compute mass matrix and internal forces coefficients ",rg
call define_arrays (Tdomain,rg)

if ( Tdomain%logicD%run_exec ) then
write (*,*) "Entering the time evolution ",rg
Tdomain%TimeD%rtime = 0
Tdomain%TimeD%NtimeMin = 0

if (Tdomain%logicD%run_restart) then
  call read_restart(Tdomain,rg)
  write (*,*) "Restart done",rg
endif
 
if (Tdomain%logicD%save_snapshots) then
    Tdomain%timeD%nsnap = Tdomain%TimeD%time_snapshots / Tdomain%TimeD%dtmin
    icount = 0
endif

icountc = 0

do ntime = Tdomain%TimeD%NtimeMin, Tdomain%TimeD%NtimeMax-1
    call Newmark (Tdomain, rg, ntime)

    Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin

    if (Tdomain%logicD%save_snapshots) then
        i_snap = mod (ntime+1, Tdomain%TimeD%nsnap) 
        if (i_snap == 0) then
            icount = icount+1
            call savefield (Tdomain,ntime,rg,icount)
            write (*,*) "A new snapshot is recorded ",rg
        endif
    endif

    if (Tdomain%logicD%save_trace) then
      if (mod (ntime+1,Tdomain%TimeD%ntrace) == 0) then
        call savetrace (Tdomain,rg,int(ntime/Tdomain%TimeD%ntrace))
        do n = 0, Tdomain%n_receivers-1
          if (rg == Tdomain%sReceiver(n)%proc) deallocate (Tdomain%sReceiver(n)%StoreTrace)
        enddo
      endif
    endif

    ! Checkpoint restart
    if (Tdomain%logicD%save_restart)  then
      if ( mod (ntime+1,Tdomain%TimeD%ncheck) == 0 ) then
         icountc = icountc + 1
         call save_checkpoint(Tdomain,Tdomain%TimeD%rtime,ntime+1,rg,icountc)
      endif
    endif
enddo

write (*,*) "Simulation is finished ",rg

endif

write (*,*) "Deallocate fields ",rg
call deallocate_domain(Tdomain)

write (*,*) "END ",rg

call mpi_finalize(code)

end program drive_sem
