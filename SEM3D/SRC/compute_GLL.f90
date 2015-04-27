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
subroutine compute_GLL(Tdomain)

    use sdomain
    use splib, only : zelegl, welegl, dmlegl

    implicit none

    type(domain), intent(inout) :: Tdomain

    integer ::  ngll, i, ndomains

    ndomains = Tdomain%n_mat

    do i = 0, ndomains-1
        !- x-part
        ngll = Tdomain%sSubdomain(i)%NGLLx
        allocate(Tdomain%sSubdomain(i)%GLLcx(0:ngll-1))
        allocate(Tdomain%sSubdomain(i)%GLLpolx(0:ngll-1))
        allocate(Tdomain%sSubdomain(i)%GLLwx(0:ngll-1))
        allocate(Tdomain%sSubdomain(i)%hprimex(0:ngll-1,0:ngll-1))
        allocate(Tdomain%sSubdomain(i)%hTprimex(0:ngll-1,0:ngll-1))

        ! USING FUNARO SUBROUTINES
        ! ZELEGL computes the coordinates of GLL points
        ! WELEGL computes the respective weights
        ! DMLEGL compute the matrix of the first derivatives in GLL points

        call zelegl (ngll-1,Tdomain%sSubdomain(i)%GLLcx,Tdomain%sSubdomain(i)%GLLpolx)
        call welegl (ngll-1, Tdomain%sSubdomain(i)%GLLcx, Tdomain%sSubdomain(i)%GLLpolx, Tdomain%sSubdomain(i)%GLLwx)
        call dmlegl (ngll-1, ngll-1, Tdomain%sSubdomain(i)%GLLcx, Tdomain%sSubdomain(i)%GLLpolx, Tdomain%sSubdomain(i)%hTprimex)

        Tdomain%sSubdomain(i)%hprimex =  TRANSPOSE ( Tdomain%sSubdomain(i)%hTprimex )

        deallocate (Tdomain%sSubdomain(i)%GLLpolx)

        if (Tdomain%sSubDomain(i)%NGLLy == Tdomain%sSubDomain(i)%NGLLx ) then
            Tdomain%sSubdomain(i)%GLLcy => Tdomain%sSubdomain(i)%GLLcx
            Tdomain%sSubdomain(i)%GLLwy => Tdomain%sSubdomain(i)%GLLwx
            Tdomain%sSubdomain(i)%hprimey => Tdomain%sSubdomain(i)%hprimex
            Tdomain%sSubdomain(i)%hTprimey => Tdomain%sSubdomain(i)%hTprimex
        else
            ngll = Tdomain%sSubdomain(i)%NGLLy
            allocate (Tdomain%sSubdomain(i)%GLLcy (0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%GLLpoly (0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%GLLwy (0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%hprimey (0:ngll-1,0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%hTprimey (0:ngll-1,0:ngll-1))
            call zelegl (ngll-1,Tdomain%sSubdomain(i)%GLLcy,Tdomain%sSubdomain(i)%GLLpoly)
            call welegl (ngll-1, Tdomain%sSubdomain(i)%GLLcy, Tdomain%sSubdomain(i)%GLLpoly, Tdomain%sSubdomain(i)%GLLwy)
            call dmlegl (ngll-1, ngll-1, Tdomain%sSubdomain(i)%GLLcy, Tdomain%sSubdomain(i)%GLLpoly, Tdomain%sSubdomain(i)%hTprimey)
            Tdomain%sSubdomain(i)%hprimey =  TRANSPOSE ( Tdomain%sSubdomain(i)%hTprimey )
            deallocate (Tdomain%sSubdomain(i)%GLLpoly)
        endif

        if (Tdomain%sSubDomain(i)%NGLLz == Tdomain%sSubDomain(i)%NGLLx ) then
            Tdomain%sSubdomain(i)%GLLcz => Tdomain%sSubdomain(i)%GLLcx
            Tdomain%sSubdomain(i)%GLLwz => Tdomain%sSubdomain(i)%GLLwx
            Tdomain%sSubdomain(i)%hprimez => Tdomain%sSubdomain(i)%hprimex
            Tdomain%sSubdomain(i)%hTprimez => Tdomain%sSubdomain(i)%hTprimex
        else if (Tdomain%sSubDomain(i)%NGLLz == Tdomain%sSubDomain(i)%NGLLy ) then
            Tdomain%sSubdomain(i)%GLLcz => Tdomain%sSubdomain(i)%GLLcy
            Tdomain%sSubdomain(i)%GLLwz => Tdomain%sSubdomain(i)%GLLwy
            Tdomain%sSubdomain(i)%hprimez => Tdomain%sSubdomain(i)%hprimey
            Tdomain%sSubdomain(i)%hTprimez => Tdomain%sSubdomain(i)%hTprimey
        else
            ngll = Tdomain%sSubdomain(i)%NGLLz
            allocate (Tdomain%sSubdomain(i)%GLLcz (0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%GLLpolz (0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%GLLwz (0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%hprimez (0:ngll-1,0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%hTprimez (0:ngll-1,0:ngll-1))
            call zelegl (ngll-1,Tdomain%sSubdomain(i)%GLLcz,Tdomain%sSubdomain(i)%GLLpolz)
            call welegl (ngll-1, Tdomain%sSubdomain(i)%GLLcz, Tdomain%sSubdomain(i)%GLLpolz, Tdomain%sSubdomain(i)%GLLwz)
            call dmlegl (ngll-1, ngll-1, Tdomain%sSubdomain(i)%GLLcz, Tdomain%sSubdomain(i)%GLLpolz, Tdomain%sSubdomain(i)%hTprimez)
            Tdomain%sSubdomain(i)%hprimez =  TRANSPOSE ( Tdomain%sSubdomain(i)%hTprimez )
            deallocate (Tdomain%sSubdomain(i)%GLLpolz)
        endif

    enddo

    return
end subroutine compute_GLL

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
