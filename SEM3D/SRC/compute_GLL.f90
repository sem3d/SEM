!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file compute_GLL.f90
!!\brief Contient la routine compute_GLL().
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<
module gll3d
    implicit none
contains
    subroutine compute_gll_data(ngll, gllc, gllw, hprime, htprime)
        use constants
        use splib, only : zelegl, welegl, dmlegl
        implicit none
        integer, intent(in) :: ngll
        real(fpp), dimension(:), allocatable, intent(out) :: gllc, gllw
        real(fpp), dimension(:,:), allocatable, intent(out) :: hprime, htprime
        !
        real(fpp), dimension(:), allocatable :: gllpol

        allocate(GLLc(0:ngll-1))
        allocate(GLLpol(0:ngll-1))
        allocate(GLLw(0:ngll-1))
        allocate(hprime(0:ngll-1,0:ngll-1))
        allocate(hTprime(0:ngll-1,0:ngll-1))

        ! USING FUNARO SUBROUTINES
        ! ZELEGL computes the coordinates of GLL points
        ! WELEGL computes the respective weights
        ! DMLEGL compute the matrix of the first derivatives in GLL points

        call zelegl (ngll-1, GLLc, GLLpol)
        call welegl (ngll-1, GLLc, GLLpol, GLLw)
        call dmlegl (ngll-1, ngll-1, GLLc, GLLpol, hTprime)
        hprime = TRANSPOSE(hTprime)
        deallocate (GLLpol)
    end subroutine compute_gll_data

    subroutine compute_GLL(Tdomain)
        use sdomain
        implicit none
        type(domain), intent(inout) :: Tdomain

        integer ::  ngll, i, ndomains

        ndomains = Tdomain%n_mat
        do i = 0, ndomains-1
            !- x-part
            ngll = Tdomain%sSubdomain(i)%NGLL
            call compute_gll_data(ngll, Tdomain%sSubdomain(i)%GLLc, Tdomain%sSubdomain(i)%GLLw, &
                Tdomain%sSubdomain(i)%hprime, Tdomain%sSubdomain(i)%htprime)
        enddo
    end subroutine compute_GLL
end module gll3d
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
