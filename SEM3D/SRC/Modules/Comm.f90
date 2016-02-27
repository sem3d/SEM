!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module scomms
    implicit none

    type :: exchange_vector
       integer :: ndata ! Taille allouee de Give/Take
       integer :: nsend ! Taille a echanger pour une communication (<=ndata)
       integer :: nsol, nsolpml, nflu, nflupml ! Nombre de point de gauss a echanger pour chaque domaine
       ! si on echange 4 ddl pour solpml et 2 pour sol alors ndata=2*nsol+4*nsolpml...
       integer :: src, dest
       integer :: ncomm ! bookeeping numero de la structure comm associee
       real, dimension(:), allocatable :: Give, Take
       integer, dimension(:), allocatable :: IGiveS, IGiveSPML, IGiveF, IGiveFPML
    end type exchange_vector

    type :: comm_vector
       integer :: ncomm
       type(exchange_vector), dimension(:), allocatable :: Data ! dimmension : 0,ncomm-1
       integer, dimension(:), allocatable :: send_reqs
       integer, dimension(:), allocatable :: recv_reqs
    end type comm_vector

    type :: comm
       ! Numero du proc avec qui on communique
       integer :: dest

       integer :: nb_faces, nb_edges, nb_vertices
       integer, dimension(:), allocatable :: faces, edges, vertices
    end type comm

contains

end module scomms

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
