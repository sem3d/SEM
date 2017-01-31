!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file compute_GLL.F90
!!\brief Contient la routine compute_GLL.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

subroutine compute_GLL(Tdomain)

    use sdomain
    use splib, only : zelegl, welegl, dmlegl

    implicit none
    type (Domain), intent (INOUT) :: Tdomain

    ! Local declarations
    integer :: i,ndomains
    integer ::  ngll

    ! Gaetano Festa 26/05/05

    ndomains = Tdomain%n_mat

    do i = 0, ndomains-1

        ngll = Tdomain%sSubdomain(i)%NGLLx

        allocate (Tdomain%sSubdomain(i)%GLLcx (0:ngll-1))
        allocate (Tdomain%sSubdomain(i)%GLLpolx (0:ngll-1))
        allocate (Tdomain%sSubdomain(i)%GLLwx (0:ngll-1))
        allocate (Tdomain%sSubdomain(i)%hprimex (0:ngll-1,0:ngll-1))
        allocate (Tdomain%sSubdomain(i)%hTprimex (0:ngll-1,0:ngll-1))

        ! USING FUNARO SUBROUTINES
        ! ZELEGL computes the coordinates of GLL points
        ! WELEGL computes the respective weights
        ! DMLEGL compute the matrix of the first derivatives in GLL points

        call zelegl (ngll-1,Tdomain%sSubdomain(i)%GLLcx,Tdomain%sSubdomain(i)%GLLpolx)
        call welegl (ngll-1, Tdomain%sSubdomain(i)%GLLcx, Tdomain%sSubdomain(i)%GLLpolx, Tdomain%sSubdomain(i)%GLLwx)
        call dmlegl (ngll-1, ngll-1, Tdomain%sSubdomain(i)%GLLcx, Tdomain%sSubdomain(i)%GLLpolx, Tdomain%sSubdomain(i)%hTprimex)

        Tdomain%sSubdomain(i)%hprimex =  TRANSPOSE ( Tdomain%sSubdomain(i)%hTprimex )

        deallocate (Tdomain%sSubdomain(i)%GLLpolx)
        if ( Tdomain%sSubdomain(i)%n_loc_dim == 2) then

            ngll = Tdomain%sSubdomain(i)%NGLLz

            allocate (Tdomain%sSubdomain(i)%GLLcz (0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%GLLpolz (0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%GLLwz (0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%hprimez(0:ngll-1,0:ngll-1))
            allocate (Tdomain%sSubdomain(i)%hTprimez (0:ngll-1,0:ngll-1))

            call zelegl (ngll-1,Tdomain%sSubdomain(i)%GLLcz,Tdomain%sSubdomain(i)%GLLpolz)
            call welegl (ngll-1, Tdomain%sSubdomain(i)%GLLcz, Tdomain%sSubdomain(i)%GLLpolz, Tdomain%sSubdomain(i)%GLLwz)
            call dmlegl (ngll-1, ngll-1, Tdomain%sSubdomain(i)%GLLcz, Tdomain%sSubdomain(i)%GLLpolz, Tdomain%sSubdomain(i)%hTprimez)

            Tdomain%sSubdomain(i)%hprimez =  TRANSPOSE ( Tdomain%sSubdomain(i)%hTprimez )
            deallocate (Tdomain%sSubdomain(i)%GLLpolz)
        else
            Tdomain%sSubdomain(i)%GLLcz => Tdomain%sSubdomain(i)%GLLcx
            Tdomain%sSubdomain(i)%GLLwz => Tdomain%sSubdomain(i)%GLLwx
            Tdomain%sSubdomain(i)%hprimez => Tdomain%sSubdomain(i)%hprimex
            Tdomain%sSubdomain(i)%hTprimez => Tdomain%sSubdomain(i)%hTprimex
        endif
    enddo

    ! Addition S. Terrana May 2016 for mortars :
    do i=0,Tdomain%n_mortar-1
        call build_mortar_interpol_matrices(Tdomain%sMortar(i))
    enddo
    return
end subroutine compute_GLL



subroutine build_mortar_interpol_matrices(Mor)

    use smortars
    use splib, only : zelegl

    implicit none
    type (Mortar), intent (INOUT) :: Mor

    ! Local declarations
    integer  :: i, j, ngmin, ngmax
    real(fpp):: xpos, res
    real(fpp), dimension (0:Mor%ngllmin-1) :: Xi_ngmin, PolMin
    real(fpp), dimension (0:Mor%ngllmax-1) :: Xi_ngmax, PolMax

    ngmin = Mor%ngllmin
    ngmax = Mor%ngllmax
    allocate(Mor%MatReinterp(0:ngmax-1,0:ngmin-1))
    allocate(Mor%MatProj(0:ngmin-1,0:ngmax-1))

    call zelegl(ngmin-1,Xi_ngmin,PolMin)
    call zelegl(ngmax-1,Xi_ngmax,PolMax)

    do i=0,ngmax-1
        xpos = Xi_ngmax(i)
        do j=0,ngmin-1
            call pol_lagrange(ngmin,Xi_ngmin,j,xpos,res)
            Mor%MatReinterp(i,j) = res
        enddo
    enddo

    do i=0,ngmin-1
        xpos = Xi_ngmin(i)
        do j=0,ngmax-1
            call pol_lagrange(ngmax,Xi_ngmax,j,xpos,res)
            Mor%MatProj(i,j) = res
        enddo
    enddo
    return

end subroutine build_mortar_interpol_matrices


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
