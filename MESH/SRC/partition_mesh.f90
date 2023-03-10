!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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
    subroutine part_mesh_3D(n_nodes, n_elem, n_points, Ipoint, nproc, dxadj, dxadjncy, procs, &
        material, tabmat)
        implicit none
        integer, intent(in)                            :: n_nodes,n_elem, nproc, n_points
        integer, intent(in), dimension(0:n_nodes-1,0:n_elem-1) :: Ipoint
        integer, allocatable, dimension(:), intent(out)  :: procs
        integer, intent(in), dimension(0:) :: material
        character, intent(in), dimension(0:) :: tabmat
        !
        integer, dimension(:), allocatable, target :: eind
        integer, dimension(:), allocatable, target :: eptr
        integer, dimension(:), allocatable :: vwgt, vsize
        integer   :: i,n,idummy
        character :: etype
        ! Metis'variables
        integer  :: edgecut
        integer, parameter  ::  numflag = 0, wgtflag = 0, ncon = 1
        integer, dimension(METIS_NOPTIONS) :: options
        integer, allocatable, dimension(:) :: adjwgt
        type(c_ptr) :: dxadj_p, dxadjncy_p
        integer, pointer, dimension(:) :: dxadj, dxadjncy
        real(kind=4), dimension(nproc) :: tpwgts
        real(kind=4), dimension(ncon) :: ubvec
        integer :: nsol, nflu, nspml, nfpml, nrand

        nsol = 0
        nflu = 0
        nspml = 0
        nfpml = 0
        nrand = 0
        allocate(eind(0:8*n_elem-1))
        allocate(eptr(0:n_elem))
        allocate(vwgt(0:n_elem))
        allocate(vsize(0:n_elem))

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
        call METIS_MeshToDual(n_elem,n_points,c_loc(eptr),c_loc(eind),4,numflag,dxadj_p,dxadjncy_p)
        call c_f_pointer(dxadj_p, dxadj, [n_elem+1])
        call c_f_pointer(dxadjncy_p, dxadjncy, [dxadj(n_elem+1)])
        allocate(procs(0:n_elem-1))
        if(nproc == 1)then
            do n = 0, n_elem-1
                procs(n) = 0
            end do
        else
            allocate(adjwgt(0:dxadj(n_elem+1)-1))
            adjwgt = 1
            do n = 0, n_elem-1
                etype = tabmat(Material(n))
                vwgt(n)=0
                select case (etype)
                case('S')
                    vwgt(n) = 3
                    nsol = nsol+1
                case('F')
                    vwgt(n) = 1
                    nflu = nflu+1
                case ('P')
                    vwgt(n) = 9
                    nspml = nspml+1
                case ('L')
                    vwgt(n) = 3
                    nfpml = nfpml+1
                case ('R')
                    vwgt(n) = 9
                    nrand = nrand+1
                end select
                if (vwgt(n)==0) then
                    write(*,*) "Erreur: elem",n, Material(n), tabmat(Material(n))
                end if
            end do
            vsize = 1
            tpwgts = real(1./nproc,kind=4)
            ubvec = real(1.001,kind=4)
            call METIS_PartGraphKway(n_elem,ncon,dxadj,dxadjncy,    &
                vwgt, vsize, adjwgt, nproc, tpwgts, ubvec, options, edgecut, procs)
        end if
        write(*,*) "Found:"
        write(*,*) nsol, "Sol"
        write(*,*) nflu, "Flu"
        write(*,*) nrand, "Rand"
        write(*,*) nspml, "Sol PML"
        write(*,*) nfpml, "Flu PML"
        deallocate(eind, eptr, vwgt, vsize)
    end subroutine part_mesh_3D
    !------------------------------------------------------------------
end module partition_mesh

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
