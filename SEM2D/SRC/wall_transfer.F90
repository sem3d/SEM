!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file wall_transfer.F90
!!\brief Contient la subroutine wall_transfer().
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

subroutine wall_transfer (Tdomain)

    ! Subroutine to evaluate the amount of data
    ! transferred in the communications
    ! Gaetano Festa 14/12/05

    use sdomain
    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer :: n, n_points, nf, nnf, ngll


    do n = 0, Tdomain%n_communications - 1
        n_points = Tdomain%sWall(n)%n_vertices
        do nf = 0, Tdomain%sWall(n)%n_faces-1
            nnf = Tdomain%sWall(n)%Face_List(nf)
            ngll = Tdomain%sFace(nnf)%ngll
            n_points = n_points + ngll - 2
        enddo
        Tdomain%sWall(n)%n_points = n_points
    enddo


    if (Tdomain%any_PML) then
        do n = 0, Tdomain%n_communications - 1
            n_points = 0
            do nf = 0, Tdomain%sWall(n)%n_pml_faces-1
                nnf = Tdomain%sWall(n)%FacePML_List(nf)
                ngll = Tdomain%sFace(nnf)%ngll
                n_points = n_points + ngll - 2
            enddo
            Tdomain%sWall(n)%n_points_pml = n_points
        enddo
    endif

    do n = 0, Tdomain%n_communications - 1
        do n_points = 0, Tdomain%sWall(n)%n_vertices-1
            nf = Tdomain%sWall(n)%Vertex_List(n_points)
            if (.not. allocated(Tdomain%sVertex(nf)%Double_Value)) then
                allocate (Tdomain%sVertex(nf)%Double_Value(0:1))
            end if
        enddo
    enddo

    return
end subroutine wall_transfer

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
