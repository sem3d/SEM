module partition_mesh

    implicit none

contains

    !-----------------------------------------------------
    subroutine part_mesh_3D(n_elem, n_points, Ipoint, nproc, dxadj, dxadjncy, procs)
        implicit none
        integer, intent(in)                            :: n_elem, nproc, n_points
        integer, intent(in), dimension(0:7,0:n_elem-1) :: Ipoint
        integer, dimension(0:8*n_elem-1)               :: elmnts
        integer   :: i,n,idummy
        ! Metis'variables
        integer  :: edgecut
        integer, parameter  ::  etype = 3, numflag = 0, wgtflag = 0, options = 0
        integer, allocatable :: vwgt(:), adjwgt(:)
        integer, allocatable, dimension(:), intent(out)  :: procs, dxadj, dxadjncy


        write(*,*)
        write(*,*) "****************************************"
        write(*,*) "  --> Partitioning"
        write(*,*) "****************************************"
        !- recasting of nodes (vertices), for use in Metis routines.
        do n = 0,n_elem-1
            do i = 0, 7   ! just the nodes at vertices of elements
                idummy = 8*n+i
                elmnts(idummy) = Ipoint(i,n)
            end do
        end do
        !- Metis'routines called
        write(*,*) "  --> Number of procs: ",nproc
        allocate(vwgt(0:n_elem-1)) ; vwgt = 0
        allocate(dxadj(0:n_elem),dxadjncy(0:8*n_elem-1))
        call METIS_MeshToDual(n_elem,n_points,elmnts,etype,numflag,dxadj,dxadjncy)
        allocate(adjwgt(0:dxadj(n_elem)-1)) ; adjwgt = 0
        allocate(procs(0:n_elem-1))
        if(nproc == 1)then
            do n = 0, n_elem-1
                procs(n) = 0
            end do
        else
            call METIS_PartGraphKway(n_elem,dxadj,dxadjncy(0:dxadj(n_elem)-1),    &
                vwgt, adjwgt, wgtflag, numflag,nproc,options,edgecut,procs)
        end if

    end subroutine part_mesh_3D
    !------------------------------------------------------------------
end module partition_mesh
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
