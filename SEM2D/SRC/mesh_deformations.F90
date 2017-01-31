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
    integer  :: i
    real(fpp):: xc, yc, theta, x, y, xnew, ynew, PI

    ! Parameters for the rotation (center coordinates and angle)
    xc = 1250.
    yc = 1250.
    PI = 4.D0*DATAN(1.D0)
    theta = PI/6.

    do i=0,Tdomain%n_glob_nodes-1
        x = Tdomain%Coord_nodes(0,i)
        y = Tdomain%Coord_nodes(1,i)
        xnew = cos(theta)*(x - xc) - sin(theta)*(y - yc) + xc
        ynew = sin(theta)*(x - xc) + cos(theta)*(y - yc) + yc
        Tdomain%Coord_nodes(0,i) = xnew
        Tdomain%Coord_nodes(1,i) = ynew
    enddo

end subroutine rotate_mesh

!####################################################

subroutine random_mesh_deformation(Tdomain)

    use sdomain

    implicit none
    type (Domain), intent (INOUT) :: Tdomain

    ! Local declarations
    integer  :: i, n1, n2
    real(fpp):: Lc, d, x, y, xnew, ynew, harv1, harv2

    ! Computes caracteristical length from the first element of the domain
    n1 = Tdomain%specel(0)%Control_Nodes(0)
    n2 = Tdomain%specel(0)%Control_Nodes(1)
    Lc = sqrt( (Tdomain%Coord_nodes(0,n1) - Tdomain%Coord_nodes(0,n2))**2 &
              +(Tdomain%Coord_nodes(1,n1) - Tdomain%Coord_nodes(1,n2))**2)
    ! Parameter for the random deformation
    d  = 0.50

    do i=0,Tdomain%n_glob_nodes-1
        x = Tdomain%Coord_nodes(0,i)
        y = Tdomain%Coord_nodes(1,i)
        call random_number(harv1)
        call random_number(harv2)
        xnew = x + Lc * d * 0.5 * (2.*harv1 - 1.)
        ynew = y + Lc * d * 0.5 * (2.*harv2 - 1.)
        Tdomain%Coord_nodes(0,i) = xnew
        Tdomain%Coord_nodes(1,i) = ynew
    enddo

end subroutine random_mesh_deformation

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
