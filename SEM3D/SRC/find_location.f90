!>
!!\file find_location.f90
!!\brief Utilitaire permettant de localiser un point dans le maillage
!!
!!
!<
module mlocations
    use sdomain
    use mshape8
    use mshape27
    implicit none
contains

    ! Finds the element(s) containing the point (x0,y0,z0)
    ! returns xi,eta,zeta coordinates within each elements
    ! nmax is the maximum number of elements returned, (can be 1) and
    ! is changed by the subroutine to reflect the number of elements found
    ! multiple elements are returned only if the point lies on elements boudaries
    ! The local coords are always within the element except when the point is outside
    ! the mesh in which case nmax is -1 on return and only one candidate is selected
    subroutine find_location(Tdomain, x0,y0,z0, nmax, elems, localcoord)
        type(domain), intent(in) :: Tdomain
        double precision, intent(in) :: x0, y0, z0
        integer, intent(inout) :: nmax
        integer, dimension(nmax) :: elems
        integer, dimension(0:2,nmax) :: localcoord
        !
        integer :: n, nnodes, i, k
        double precision :: mindist, dist
        double precision :: best_node
        double precision, dimension(0:2) :: p
        !! Find the closest node
        nnodes = Tdomain%n_nodes
        best_node = 0
        mindist = 1d100
        do n = 0, Tdomain%n_glob_nodes-1
            p = Tdomain%GlobCoord(:,n)
            dist = (p(0)-x0)**2 + (p(1)-y0)**2 + (p(2)-z0)**2
            if (dist<mindist) then
                best_node = n
                dist = mindist
            end if
        end do
        !! Now find all the elements in contact with the best_node
        k = 0
        do n = 0, Tdomain%n_elem-1
            do i = 0,nnodes-1
                if (Tdomain%specel(n)%Control_Nodes(i)==best_node) then
                    k = k+1
                    elems(k) = n
                    if (k>=nmax) exit
                end if
            end do
        end do
        
    end subroutine find_location

end module mlocations
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
