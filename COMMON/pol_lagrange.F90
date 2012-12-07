!>
!!\file pol_lagrange.F90
!!\brief Evalue une fonction et sa dérivée via les polynômes de Lagrange.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param integer, intent (IN) n
!! \param integer, intent (IN) k
!! \param real, dimension (0:n-1), intent (IN) GLLc
!! \param real,intent (IN) x
!! \param real, intent (OUT) y
!<


subroutine pol_lagrange (n,GLLc,k,x,y)
    implicit none
    integer, intent (IN) :: n,k
    real, dimension  (0:n-1), intent (IN) :: GLLc
    real,intent (IN) :: x
    real, intent (OUT) :: y

    integer :: i

    y = 1
    if (n ==0 ) then
        write (*,*) "Bad number n = 0. It should be a constant!"
        stop
    endif
    if (n ==1 ) return
    do i =0,n-1
        if (i /= k)  y = y * (x-GLLc(i))/(GLLc(k)-GLLc(i))
    enddo
    return
end subroutine pol_lagrange


! ################################################################
!>
!! \brief
!!
!! \param integer, intent (IN) n
!! \param integer, intent (IN) ip
!! \param real, dimension (0:n-1), intent (IN) GLLc
!! \param real, intent (IN) x0
!! \param real, intent (OUT) dy0
!<


subroutine DERIVAL (GLLc,n,ip,x0,dy0)

    implicit none
    integer, intent (IN) :: n,ip
    real, dimension  (0:n-1), intent (IN) :: GLLc
    real, intent (IN) :: x0
    real, intent (OUT):: dy0

    integer :: i,k
    real :: y0

    ! Compute derivatives of Lagrange polynomials

    dy0 = 0
    do i = 0 ,n-1
        y0 = 1
        if (i /= ip) then
            do k = 0,n-1
                if (k /= i .and. k /= ip)  then
                    y0  = y0 * (x0 - GLLc(k) ) /(GLLc(ip)-GLLc(k) )
                endif
            enddo
            y0 = y0 / (GLLc(ip)-GLLc(i) )
            dy0 = dy0 + y0
        endif
    enddo
    return
end subroutine DERIVAL
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
