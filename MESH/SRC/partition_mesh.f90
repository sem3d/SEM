module partition_mesh
    use iso_c_binding
    implicit none

    integer, parameter :: METIS_NOPTIONS = 40

    interface
        subroutine METIS_MeshToDual(ne,nn,eptr,eind,ncommon,numflag,xadj,dxadjncy) bind(c)
            use iso_c_binding
            integer(C_INT), intent(in) :: ne
            integer(C_INT), intent(in) :: nn
            type(C_PTR), intent(in), value :: eptr
            type(C_PTR), intent(in), value :: eind
            integer(C_INT), intent(in) :: ncommon
            integer(C_INT), intent(in) :: numflag
            type(c_ptr), intent(out) :: xadj
            type(c_ptr), intent(out) :: dxadjncy
        end subroutine METIS_MeshToDual

    end interface
contains

    !-----------------------------------------------------
    subroutine part_mesh_3D(n_elem, n_points, Ipoint, nproc, dxadj, dxadjncy, procs)
        implicit none
        integer, intent(in)                            :: n_elem, nproc, n_points
        integer, intent(in), dimension(0:7,0:n_elem-1) :: Ipoint
        integer, dimension(0:8*n_elem-1),target        :: eind
        integer, dimension(0:n_elem),target            :: eptr
        integer, dimension(0:n_elem-1)                 :: vwgt, vsize
        integer   :: i,n,idummy
        ! Metis'variables
        integer  :: edgecut
        integer, parameter  ::  etype = 3, numflag = 0, wgtflag = 0, ncon = 1
        integer, dimension(METIS_NOPTIONS) :: options
        integer, allocatable, dimension(:) :: adjwgt
        integer, allocatable, dimension(:), intent(out)  :: procs
        type(c_ptr) :: dxadj_p, dxadjncy_p
        integer, pointer, dimension(:) :: dxadj, dxadjncy
        real(kind=4), dimension(nproc) :: tpwgts
        real(kind=4), dimension(ncon) :: ubvec

        write(*,*)
        write(*,*) "****************************************"
        write(*,*) "  --> Partitioning"
        write(*,*) "****************************************"
        !- recasting of nodes (vertices), for use in Metis routines.
        do n = 0,n_elem-1
            do i = 0, 7   ! just the nodes at vertices of elements
                idummy = 8*n+i
                eind(idummy) = Ipoint(i,n)
                eptr(n)=n*8
            end do
        end do
        eptr(n_elem)=8*n_elem
        !- Metis'routines called
        write(*,*) "  --> Number of procs: ",nproc
        !allocate(dxadj(0:n_elem),dxadjncy(0:8*n_elem-1))
        call METIS_SetDefaultOptions(options)
        call METIS_MeshToDual(n_elem,n_points,c_loc(eptr),c_loc(eind),1,numflag,dxadj_p,dxadjncy_p)
        call c_f_pointer(dxadj_p, dxadj, [n_elem+1])
        call c_f_pointer(dxadjncy_p, dxadjncy, [dxadj(n_elem+1)])
        allocate(procs(0:n_elem-1))
        if(nproc == 1)then
            do n = 0, n_elem-1
                procs(n) = 0
            end do
        else
            allocate(adjwgt(0:dxadj(n_elem)-1)) ; adjwgt = 1
            vwgt = 1
            vsize = 1
            tpwgts = 1./nproc
            ubvec = 1.001
            call METIS_PartGraphKway(n_elem,ncon,dxadj,dxadjncy,    &
                vwgt, vsize, adjwgt, nproc, tpwgts, ubvec, options, edgecut, procs)
        end if

    end subroutine part_mesh_3D
    !------------------------------------------------------------------
end module partition_mesh
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
