!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP

module msnapdata

    use constants
    implicit none

    ! Type assurant la gestion des sorties snapshots
    type :: output_var_t
        ! Communicator used for output grouping. Only rank 0 of this comm produces outputs
        integer :: comm
        ! Nombre de processeur dans le groupe de communication associe a comm_output
        integer :: nprocs
        integer :: rank
        ! Numero du groupe de sauvegarde
        integer :: group

!        integer, dimension(:), allocatable :: output_nodes, output_nodes_offset, output_elems
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer, allocatable, dimension(:) :: domains
        ! Nombre de noeuds elements et cellules par proc
        integer :: nnodes
        integer :: ncells
        integer :: nelems
        ! Nodes fields
        real(fpp), dimension(:,:), allocatable :: displ, veloc, accel
        real(fpp), dimension(:)  , allocatable :: press_n
        ! Cell fields
        real(fpp), dimension(:)  , allocatable   :: press_c, eps_vol
        real(fpp), dimension(:,:), allocatable   :: eps_dev, sig_dev
        real(fpp), dimension(:), allocatable     :: P_energy, S_energy
        real(fpp), dimension(:,:), allocatable   :: eps_dev_pl

        ! Storage for communications
        integer, dimension(:),   allocatable :: displs_n, displs_c
        integer, dimension(:),   allocatable :: counts_n, counts_c
        real(fpp), dimension(:), allocatable :: all_data_1d_n, all_data_1d_c
        integer, dimension(:),   allocatable :: displs2d_n ! so far no vectors on cells
        integer, dimension(:),   allocatable :: counts2d_n
        real(fpp), dimension(:,:), allocatable :: all_data_2d_n
        ! Sizes (group)
        integer :: ntot_nodes, ntot_cells
    end type output_var_t

end module msnapdata
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
