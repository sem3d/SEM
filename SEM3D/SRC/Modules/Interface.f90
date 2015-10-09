!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Domain.f90
!!\brief Contient les definition des surfaces et interfaces
!!
!<

module sinterface
    implicit none

    ! Definition d'une surface (liste des éléments faces, edges, vertex)
    type :: surf_num
        integer :: n_faces, n_edges, n_vertices
        integer :: nbtot ! nombre total de points de gauss de l'interface
        ! A list of gll points from the interface
        integer, dimension(:), allocatable :: map ! dimension(0:nbtot-1)
        ! List of faces, edges, vertices
        integer, dimension(:), allocatable :: if_faces
        integer, dimension(:), allocatable :: if_edges
        integer, dimension(:), allocatable :: if_vertices
        ! mapping between face index and a sign to apply whether the local face has
        ! the same orientation as the surface outside normal
        integer, dimension(:), allocatable :: if_norm
    end type surf_num

    type :: inter_num
        type(surf_num) :: surf0
        type(surf_num) :: surf1
    end type inter_num
contains
    ! simple initialisation
    subroutine init_surface(surf)
        type(surf_num), intent(inout) :: surf

        surf%nbtot = 0
        surf%n_faces = 0
        surf%n_edges = 0
        surf%n_vertices = 0
    end subroutine init_surface

    ! simple initialisation
    subroutine init_interface(inter)
        type(inter_num), intent(inout) :: inter

        call init_surface(inter%surf0)
        call init_surface(inter%surf1)
    end subroutine init_interface

    ! Free temporary memory (ie faces, edges, vertices) that we won't need during computation
    subroutine free_temp_surface(surf)
        type(surf_num), intent(inout) :: surf
        if (surf%n_faces>0) then
            deallocate(surf%if_faces)
            deallocate(surf%if_norm)
        end if
        if (surf%n_edges>0) then
            deallocate(surf%if_edges)
        end if
        if (surf%n_vertices>0) then
            deallocate(surf%if_vertices)
        end if
    end subroutine free_temp_surface

    subroutine allocate_surface(surf)
        type(surf_num), intent(inout) :: surf

        if (surf%n_faces>0) then
            allocate(surf%if_faces(0:surf%n_faces-1))
            allocate(surf%if_norm(0:surf%n_faces-1))
            ! We Don't allocate map  since we can't know yet the number of glls
        end if
        if (surf%n_edges>0) then
            allocate(surf%if_edges(0:surf%n_edges-1))
        end if
        if (surf%n_vertices>0) then
            allocate(surf%if_vertices(0:surf%n_vertices-1))
        end if
    end subroutine allocate_surface

    subroutine allocate_interface(inter)
        type(inter_num), intent(inout) :: inter

        call allocate_surface(inter%surf0)
        call allocate_surface(inter%surf1)
    end subroutine allocate_interface
end module sinterface

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
