!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module sspline
    use constants, only : fpp

    public :: bmn,bspln,bmdiag,bchfac
    private

contains

    subroutine bmn(b,m,n)
        implicit none
        integer :: m,n,i
        real(fpp) :: b((n+1)*m+1)

        do i = 1,(n+1)*m+1
            b(i) = bspln(0,n,dble(i-1)/dble((n+1)*m)*dble(n+1))
        enddo
        return

    end subroutine bmn

    recursive function bspln(i,k,x) result(b)
        implicit none

        real(fpp) :: x,b
        integer :: k,i

        b=0.d0
        if(k+1>0)then
            if(k+1==1)then
                if(x>=i.and.x<i+1)b=1.d0
            else
                b=(x-i)*bspln(i,k-1,x)/(k+1-1)+(i+k+1-x)*bspln(i+1,k-1,x)/(k+1-1)
                if(k==0)b=0
            endif
        endif
    end function bspln

    subroutine bmdiag(m,n,nsp,nt,spmat)
        implicit none
        integer, intent(in) :: m,n,nsp,nt
        real(fpp) :: spmat(n+1,nsp),b((n+1)*m+1)
        integer :: i,j,i1,i2,i3,i4

        call bmn(b,m,n)
        do j = 1,nsp
            do i = 1,n+1
                i1 = max(1+(i-1)*m,(n+1-j)*m+1)
                i2 = min(1+(n+1)*m,(n-j+1)*m+nt)
                i3 = max(1,(n+2-j-i)*m+1)
                i4 = min(1+(n+2-i)*m,(n+2-j-i)*m+nt)
                if (j+i-1<=nsp)then
                    spmat(i,j) = dot_product(b(i1:i2),b(i3:i4))
                else
                    spmat(i,j) = 0.
                endif
            enddo
        enddo
        return

    end subroutine bmdiag

    subroutine bchfac(w,nbands,nrow,diag)
        implicit none
        real(fpp), intent(inout) :: w(nbands,nrow)
        integer, intent(in) :: nbands,nrow
        real(fpp), intent(out) :: diag(nrow)
        integer :: i,imax,j,jmax,n
        real(fpp) :: ratio

        if (nrow<=1) then
            if (w(1,1)>0.) then
                w(1,1) = 1./w(1,1)
            endif
            return
        endif
        diag(1:nrow) = w(1,1:nrow)
        do n = 1,nrow
            if (w(1,n)+diag(n)<=diag(n)) then
                w(1:nbands,n) = 0.
            else
                w(1,n) = 1./w(1,n)
                imax = min(nbands-1,nrow-n)
                jmax = imax
                do i = 1,imax
                    ratio = w(i+1,n)*w(1,n)
                    do j = 1,jmax
                        w(j,n+i) = w(j,n+i)-w(j+i,n)*ratio
                    enddo
                    jmax = jmax-1
                    w(i+1,n) = ratio
                enddo
            endif
        enddo
        return

    end subroutine bchfac

end module sspline
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
