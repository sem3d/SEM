!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file pol_lagrange.F90
!!\brief Evalue une fonction et sa dérivée via les polynômes de Lagrange.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

subroutine pol_lagrange(n,GLLc,k,x,y)
    ! value of Lagrange polynomial at x, on the basis defined by GLLc
    implicit none
    integer, intent(in) :: n,k
    double precision, dimension(0:n-1), intent(in) :: GLLc
    double precision,intent(in) :: x
    double precision, intent(out) :: y
    integer :: i

    y = 1
    if(n == 0)then
        stop "Bad number n = 0. It should be a constant!"
    endif
    if(n == 1) return
    do i =0,n-1
        if(i /= k)  y = y * (x-GLLc(i))/(GLLc(k)-GLLc(i))
    enddo

    return
end subroutine pol_lagrange
!------------------------------------------------------------------------------
subroutine derivlag(GLLc,n,ip,x0,dy0)

    ! computes derivatives of Lagrange polynomials given by GLLc, at x0
    implicit none
    integer, intent(in) :: n,ip
    double precision, dimension(0:n-1), intent(in) :: GLLc
    double precision, intent(in) :: x0
    double precision, intent(out):: dy0
    integer :: i,k
    double precision :: y0


    dy0 = 0
    do i = 0 ,n-1
        y0 = 1
        if(i /= ip)then
            do k = 0,n-1
                if(k /= i .and. k /= ip)then
                    y0 = y0*(x0-GLLc(k))/(GLLc(ip)-GLLc(k))
                endif
            enddo
            y0 = y0/(GLLc(ip)-GLLc(i))
            dy0 = dy0 + y0
        endif
    enddo

    return

end subroutine derivlag

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
