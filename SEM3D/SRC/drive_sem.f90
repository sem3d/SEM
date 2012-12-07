program drive_sem

    use sdomain
    use mpi
    implicit none

    type(domain) :: Tdomain
    integer :: code, rg, n_proc, ntime, i_snap, n,i,j, icount, icountc
    character(len=100) :: fnamef


!- Communicators (MPI) initialization 
call MPI_INIT(code)
call MPI_COMM_RANK(MPI_COMM_WORLD,rg,code)
call MPI_COMM_SIZE(MPI_COMM_WORLD,n_proc,code)

!------------  Starting the code --------------------
if(rg == 0)then
    write(*,*) "-----------------------------------------"
    write(*,*) "----         SEM - 3d version      ------"
    write(*,*) "-----------------------------------------"
    write(*,*)
end if
call MPI_BARRIER(MPI_COMM_WORLD,code)

write(*,*) "  --> Reading run properties, proc. # ",rg
call read_input(Tdomain,rg)
call MPI_BARRIER(MPI_COMM_WORLD, code)

if(Tdomain%logicD%neumann_local_present)then
  write(*,*) "  --> Defining Neumann properties, proc. # ",rg
  call define_Neumann_properties(Tdomain,rg)
endif
call MPI_BARRIER(MPI_COMM_WORLD, code)

write(*,*) "  --> Computing Gauss-Lobatto-Legendre weights and zeroes, proc. # ",rg
call compute_GLL(Tdomain)
call MPI_BARRIER(MPI_COMM_WORLD, code)

write(*,*) "  --> Defining global numbering for collocation points, proc. # ",rg
call global_numbering(Tdomain,rg)
call MPI_BARRIER(MPI_COMM_WORLD, code)

write(*,*) "  --> Computing shape functions within their derivatives, proc. # ",rg
if(Tdomain%n_nodes == 8) then
    call shape8(Tdomain,rg)   ! Linear interpolation
else if(Tdomain%n_nodes == 27) then
    write(*,*) "    Sorry, 27 control points not yet implemented in the code... Wait for an upgrade ",rg
    stop
else
    write(*,*) "    Bad number of nodes for hexaedral shape ",rg
    stop
endif
call MPI_BARRIER(MPI_COMM_WORLD,code)


write(*,*) "  --> Computing Courant parameter, proc. #",rg
call compute_Courant(Tdomain,rg)
call MPI_BARRIER(MPI_COMM_WORLD,code)

if(Tdomain%any_PML)then
    write(*,*) "  --> Attribute PML properties, proc. # ",rg
    call PML_definition(Tdomain) 
endif
call MPI_BARRIER(MPI_COMM_WORLD,code)

if(Tdomain%logicD%any_source)then
    write(*,*) "  --> Computing point-source parameters",rg
    call SourcePosition(Tdomain,n_proc,rg)
    call source_excit(Tdomain,rg)
endif
if(Tdomain%logicD%save_trace)then
    write(*,*) "  --> Computing receivers'locations, proc. # ",rg
    call ReceiverExactPosition(Tdomain,rg)
endif
call MPI_BARRIER(MPI_COMM_WORLD,code)

write(*,*) "  --> Allocating fields, proc. # ",rg
call allocate_domain(Tdomain,rg)
call MPI_BARRIER(MPI_COMM_WORLD,code)

write(*,*) "  --> Computing mass matrix and internal forces coefficients, proc. # ",rg
call define_arrays(Tdomain,rg)
call MPI_BARRIER(MPI_COMM_WORLD,code)

if(Tdomain%logicD%run_exec)then
    write(*,*) "  --> Entering the time evolution, proc. # ",rg
    Tdomain%TimeD%rtime = 0
    Tdomain%TimeD%NtimeMin = 0
    call MPI_BARRIER(MPI_COMM_WORLD,code)

    if(Tdomain%logicD%run_restart)then
        call read_restart(Tdomain,rg)
        write(*,*) "  --> Restarting : OK.",rg
    endif
 
    if(Tdomain%logicD%save_snapshots)then
        Tdomain%timeD%nsnap = Tdomain%TimeD%time_snapshots/Tdomain%TimeD%dtmin
        icount = 0
    endif

    icountc = 0

    do ntime = Tdomain%TimeD%NtimeMin, Tdomain%TimeD%NtimeMax-1
        call Newmark(Tdomain,rg,ntime)

        Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin

        if(Tdomain%logicD%save_snapshots)then
            i_snap = mod(ntime+1,Tdomain%TimeD%nsnap) 
            if(i_snap == 0)then
                icount = icount+1
                call savefield(Tdomain,ntime,rg,icount)
                write (*,*) "    --> New snapshot recorded ",rg
            endif
        endif

        if(Tdomain%logicD%save_trace) then
            if(mod(ntime+1,Tdomain%TimeD%ntrace) == 0)then
                call savetrace(Tdomain,rg,int(ntime/Tdomain%TimeD%ntrace))
                do n = 0, Tdomain%n_receivers-1
                    if(rg == Tdomain%sReceiver(n)%proc)  &
                          deallocate(Tdomain%sReceiver(n)%StoreTrace)
                enddo
            endif
        endif

    ! Checkpoint restart
        if(Tdomain%logicD%save_restart)then
            if(mod(ntime+1,Tdomain%TimeD%ncheck) == 0)then
                icountc = icountc + 1
                call save_checkpoint(Tdomain,Tdomain%TimeD%rtime,ntime+1,rg,icountc)
            endif
        endif
    enddo

    write(*,*) "  --> Simulation is over.",rg
endif

write(*,*) "  --> Deallocating fields, ",rg
call deallocate_domain(Tdomain)

write(*,*) "  --> Ok for proc. #",rg

call MPI_FINALIZE(code)

write(*,*) "------------------------------------------"
write(*,*) "---------------    END  ------------------"
write(*,*) "------------------------------------------"

end program drive_sem
