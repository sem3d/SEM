!>
!!\file mesh_deformations.F90
!!\brief Contains some routine that can deform the initial mesh, for testing purposes.
!!\author Sebastien Terrana
!!\version 1.0
!!\date 01/09/2014
!!
!<

subroutine rotate_mesh(Tdomain)

    use sdomain

    implicit none
    type (Domain), intent (INOUT) :: Tdomain

    ! Local declarations
    integer :: i
    real    :: xc, yc, theta, x, y, xnew, ynew, PI

    ! Parameters for the rotation (center coordinates and angle)
    xc = 1250.
    yc = 1250.
    PI = 4.D0*DATAN(1.D0)
    theta = PI/4.

    do i=0,Tdomain%n_glob_nodes-1
        x = Tdomain%Coord_nodes(0,i)
        y = Tdomain%Coord_nodes(1,i)
        xnew = cos(theta)*(x - xc) - sin(theta)*(y - yc) + xc
        ynew = sin(theta)*(x - xc) + cos(theta)*(y - yc) + yc
        Tdomain%Coord_nodes(0,i) = xnew
        Tdomain%Coord_nodes(1,i) = ynew
    enddo

end subroutine rotate_mesh

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
