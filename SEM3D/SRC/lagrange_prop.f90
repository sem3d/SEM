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
!------------------------------------------------------------------------
real function interp_lag(Tdomain,nel,ngllx,nglly,ngllz,mat,xi,eta,zeta,func)
 ! gives the Lagrange interpolation at (xi,eta,zeta) for a function func
 !    whose values are known at GLL points
    use sdomain
    implicit none
    
    type(domain), intent(in)   :: Tdomain
    integer, intent(in)   :: nel,mat,ngllx,nglly,ngllz
    real, intent(in)  :: xi,eta,zeta
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: func 
    integer :: i,j,k
    real  :: f,wxi,weta,wzeta

    f = 0d0

    do k = 0,ngllz-1
        call pol_lagrange(ngllz,Tdomain%sSubdomain(mat)%GLLcz,k,zeta,wzeta)
        do j = 0,nglly-1
            call pol_lagrange(nglly,Tdomain%sSubdomain(mat)%GLLcy,j,eta,weta)
            do i = 0,ngllx-1
                call pol_lagrange(ngllx,Tdomain%sSubdomain(mat)%GLLcx,i,xi,wxi)
                f = f+func(i,j,k)*wxi*weta*wzeta
            end do
        end do
    end do
    
    interp_lag = f

    return
end function interp_lag
