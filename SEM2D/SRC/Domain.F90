!>
!!\file Domain.F90
!!\brief Contient le définition du type domain
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module sdomain

    ! Modified by Gaetano 12/10/05

    use selement
    use sfaces
    use svertices
    use ssources
    use stimeparam
    use logical_input
    use ssubdomains
    use sreceivers
    use sfaults
    use mpi_list
    use communication_object



    type :: domain

       integer :: communicateur
#ifdef MKA3D
       integer :: communicateur_global
       integer :: master_superviseur
#endif

       integer :: n_elem, n_face, n_vertex, n_source,n_dime,n_glob_nodes, n_line ,n_receivers
       integer :: n_nodes, n_mat,n_glob_points, n_super_object, n_fault, n_communications
       integer, dimension (:), pointer :: Line_index, Communication_list

       real, dimension (:,:), pointer :: Coord_nodes, GlobCoord
       real, dimension(:,:,:), pointer :: Store_Trace


       real, dimension(:,:), pointer :: GrandeurVitesse     ! tableau conservatoire des grandeurs
       real, dimension(:,:), pointer :: GrandeurDepla     ! tableau conservatoire des grandeurs
       real, dimension(:), pointer :: GrandeurDeformation ! pour tous les points de gauss a une iteration donnee (utilise pour les sorties capteurs)


       character (len=30) :: Title_simulation, mesh_file, station_file, material_file
       character (len=1), dimension (:), pointer ::  Name_Line, Super_object_type
       character (len=30), dimension (:), pointer :: Super_object_file

       logical :: any_PML,bMailUnv,bCapteur

       type (time) :: TimeD
       type (logical_array) :: logicD
       type (mpi_objects) :: Mpi_var
       type(element), dimension(:), pointer :: specel
       type(face), dimension (:), pointer :: sface
       type (vertex), dimension (:), pointer :: svertex
       type (source), dimension (:), pointer :: sSource
       type (subdomain), dimension (:), pointer :: sSubDomain
       type (receiver), dimension (:), pointer :: sReceiver
       type (fault), dimension (:), pointer :: sFault
       type (communicating_wall), dimension (:), pointer :: sWall



    end type domain

contains

    subroutine dist_max_elem(Tdomain)

        type (Domain), intent (INOUT) :: Tdomain
        integer ipoint, jpoint
        integer n
        real coor_i(0:1), coor_j(0:1)
        !!real dist_max


        do n = 0,Tdomain%n_elem-1
            dist_max = 0.

            do i=0,Tdomain%n_nodes-1
                ipoint = Tdomain%specel(n)%Control_Nodes(i)
                coor_i = Tdomain%Coord_nodes(0:1,ipoint)
                do j=i+1,Tdomain%n_nodes-1
                    jpoint = Tdomain%specel(n)%Control_Nodes(j)
                    coor_j = Tdomain%Coord_nodes(0:1,jpoint)
                    dist_max = max(dist_max, sqrt((coor_i(0)-coor_j(0))**2+(coor_i(1)-coor_j(1))**2))
                enddo
            enddo
            Tdomain%specel(n)%dist_max = dist_max
        enddo

    end subroutine dist_max_elem



end module sdomain
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
