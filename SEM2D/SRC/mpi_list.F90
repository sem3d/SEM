!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file mpi_list.F90
!!\brief Contient la définition du type mpi_objects.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module mpi_list

    ! Written by Gaetano Festa 12/10/2005


    type :: mpi_objects

       integer :: my_rank, n_proc, number_of_communications
       integer, dimension (:), pointer :: communication_list

    end type mpi_objects


end module mpi_list

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
