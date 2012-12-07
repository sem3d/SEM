!>
!! \file shape27.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

!>
!! shape27: alloue et calcule la jacobienne et l'inverse du gradient de ??
!<
subroutine shape27(Tdomain)

    ! Modified by Paul Cupillard 26/06/2006


    use sdomain

    implicit none

    type(domain), target, intent (INOUT) :: Tdomain

    integer :: n, ngllx,nglly,ngllz, mat, i,j,k, i_count,j_count, ipoint
    real :: xi,eta,zeta, xp,yp,zp, f, Jac
    real, dimension(:,:), allocatable :: coord
    real, dimension(0:2,0:2) :: LocInvGrad


    allocate (Tdomain%GlobCoord(0:2,0:Tdomain%n_glob_points-1))

    do n = 0,Tdomain%n_elem - 1

        allocate (coord(0:Tdomain%n_nodes-1,0:2))
        do i = 0,Tdomain%n_nodes-1
            j = Tdomain%specel(n)%Control_Nodes(i)
            coord(i,0:2) = Tdomain%Coord_Nodes(0:2,j)
        enddo

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        mat = Tdomain%specel(n)%mat_index

        allocate (Tdomain%specel(n)%Jacob(0:ngllx-1,0:nglly-1,0:ngllz-1) )
        allocate (Tdomain%specel(n)%InvGrad(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2) )

        do k = 0,ngllz - 1
            zeta =  Tdomain%sSubdomain(mat)%GLLcz (k)
            do j = 0,nglly - 1
                eta =  Tdomain%sSubdomain(mat)%GLLcy (j)
                do i = 0,ngllx - 1
                    xi = Tdomain%sSubdomain(mat)%GLLcx (i)

                    xp = 0;   yp = 0;   zp = 0;
                    do i_count = 0,Tdomain%n_nodes-1
                        f = Comp_shapefunc(i_count,xi,eta,zeta)
                        xp = xp + coord(i_count,0)*f
                        yp = yp + coord(i_count,1)*f
                        zp = zp + coord(i_count,2)*f
                    enddo
                    ipoint = Tdomain%specel(n)%Iglobnum(i,j,k)
                    Tdomain%GlobCoord (0,ipoint) = xp
                    Tdomain%GlobCoord (1,ipoint) = yp
                    Tdomain%GlobCoord (2,ipoint) = zp

!!! Computation of the derivative matrix, dx_(jj)/dxi_(ii) !!!
                    LocInvGrad = 0.
                    do i_count = 0,Tdomain%n_nodes-1
                        do j_count = 0,2
                            f = Comp_derivshapefunc(i_count,xi,eta,zeta,j_count)
                            LocInvGrad(j_count,0) = LocInvGrad(j_count,0) + coord(i_count,0)*f
                            LocInvGrad(j_count,1) = LocInvGrad(j_count,1) + coord(i_count,1)*f
                            LocInvGrad(j_count,2) = LocInvGrad(j_count,2) + coord(i_count,2)*f
                        enddo
                    enddo
                    call invert_3d (LocInvGrad, Jac)
                    Tdomain%specel(n)%Jacob(i,j,k) = Jac
                    Tdomain%specel(n)%InvGrad (i,j,k,0:2,0:2) = LocInvGrad(0:2,0:2)
                enddo
            enddo
        enddo

        deallocate (coord)

    enddo

end subroutine shape27
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
