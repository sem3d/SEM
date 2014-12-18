
module mleastsq
    implicit none
    public :: minimize_cg
    ! if we need these put them in a separate module
    private :: vnorm, vdot, mnorm, line_search
contains

    double precision function vnorm(dim, x)
        integer, intent(in) :: dim
        double precision, dimension(0:dim-1), intent(in) :: x
        !
        double precision :: sum
        integer :: i
        sum = 0.0
        do i=0, dim-1
            sum = sum + x(i)*x(i)
        end do
        vnorm = sqrt(sum)
    end function vnorm

    double precision function vdot(dim, x, y)
        integer, intent(in) :: dim
        double precision, dimension(0:dim-1), intent(in) :: x, y
        !
        integer :: i
        vdot = 0.0
        do i=0, dim-1
            vdot = vdot + x(i)*y(i)
        end do
    end function vdot

    double precision function mnorm(dim, m)
        integer, intent(in) :: dim
        double precision, dimension(0:dim-1,0:dim-1), intent(in) :: m
        !
        double precision :: sum
        integer :: i,j
        sum = 0.0
        do j=0, dim-1
            do i=0, dim-1
                sum = sum + m(i,j)*m(i,j)
            end do
        end do
        mnorm = sqrt(sum)
    end function mnorm

    subroutine line_search(dim, nn, args0, args1, fun, x0, dir, alpha)
        integer, intent(in) :: dim, nn
        double precision, dimension(0:dim-1), intent(in) :: args1, x0, dir
        double precision, dimension(0:dim-1,0:nn-1), intent(in) :: args0
        double precision, intent(inout) :: alpha
        !
        interface
           double precision function fun(dim, nn, x, args0, args1)
               integer, intent(in) :: dim, nn
               double precision, dimension(0:dim-1), intent(in) :: x, args1
               double precision, dimension(0:dim-1,0:nn-1), intent(in) :: args0
           end function fun
        end interface
        double precision :: r, val
        double precision :: fmin, alphamin
        ! Very simple line search that chooses between 0.1 x alpha and 2 x alpha
        alphamin = alpha
        fmin = 1d100
        r = 0.1
        do while (r<2.0)
            val = fun(dim, nn, x0+r*alpha*dir, args0, args1)
            if (val<fmin) then
                fmin = val
                alphamin = r*alpha
            end if
            r = r + 0.1
        end do
        alpha = alphamin
    end subroutine line_search
    
    subroutine minimize_cg(dim, nn, x0, args0, args1, fun, gradfun, alpha0, xout, niter)
        ! Minimize fun(x) with gradfun using conjugate gradient method
        ! starting with initial guess x0
        ! 
        integer, intent(in) :: dim, nn
        double precision, dimension(0:dim-1), intent(in) :: x0, args1
        double precision, dimension(0:dim-1), intent(out) :: xout
        double precision, dimension(0:dim-1,0:nn-1), intent(in) :: args0
        double precision, intent(in) :: alpha0
        integer, intent(out) :: niter
        interface
           double precision function fun(dim, nn, x, args0, args1)
               integer, intent(in) :: dim, nn
               double precision, dimension(0:dim-1), intent(in) :: x, args1
               double precision, dimension(0:dim-1,0:nn-1), intent(in) :: args0
           end function fun
           subroutine gradfun(dim, nn, x, args0, args1, grad)
               integer, intent(in) :: dim, nn
               double precision, dimension(0:dim-1), intent(in) :: x, args1
               double precision, dimension(0:dim-1,0:nn-1), intent(in) :: args0
               double precision, dimension(0:dim-1), intent(out) :: grad
           end subroutine gradfun
        end interface
        !
        double precision, dimension(0:dim-1) :: grad, pgrad, sdir, gdiff
        double precision :: gnorm, old_gnorm, alpha, beta
        integer :: i
        !
        integer, parameter :: maxiter=1000
        double precision, parameter :: gtol=1e-20
        
        alpha = alpha0
        xout = x0
        call gradfun(dim, nn, xout, args0, args1, grad)
        pgrad = -grad
        sdir = pgrad
        gnorm = vdot(dim, grad, grad)
        old_gnorm = gnorm
        do i=1, maxiter
            write(*,*) gnorm, xout, grad
            if (gnorm<gtol) exit
            call line_search(dim, nn, args0, args1, fun, xout, sdir, alpha)
            xout = xout + alpha * pgrad
            call gradfun(dim, nn, xout, args0, args1, grad)
            !sdir = -grad
            ! Conjugate
            gdiff = grad+pgrad
            ! Polak-Ribiere
            beta = vdot(dim, grad, gdiff)/gnorm
            if (beta<0) beta = 0d0
            ! Pure gradient
            !beta = 0
            sdir = -grad + beta*sdir
            pgrad = -grad
            gnorm = vdot(dim, grad, grad)
        end do
        niter = i
    end subroutine minimize_cg
end module mleastsq
