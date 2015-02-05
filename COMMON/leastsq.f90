
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
        r = 0.05
        do while (r<1.1)
            val = fun(dim, nn, x0+r*alpha*dir, args0, args1)
            if (val<fmin) then
                fmin = val
                alphamin = r*alpha
            end if
            r = r + 0.05
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

    subroutine cg_inv(m,n,Amat,dataval,model_val)
        !-> inversion of the stiffness matrix, and calculation
        !   of the corrector's coefficient
        implicit none
        integer, intent(in)  :: m,n
        real, dimension(0:m-1,0:n-1),intent(in)  :: Amat
        real, dimension(0:m-1), intent(in)  :: dataval
        real, dimension(0:n-1), intent(out)  :: model_val
        real, allocatable   :: Atrans(:,:), Acp(:,:), Dmat(:,:)
        real, allocatable   :: U(:,:),VT(:,:),SI(:),sigma(:,:),WORK(:)
        integer             :: info,lwork,i

        !---------------------------------------------------
        !---------------------------------------------------
        allocate(Acp(m,n),Atrans(n,m))
        allocate(Dmat(n,m))

        Acp(:,:) = real(Amat(0:,0:),8)

        lwork = 2*(3*min(M,N)*min(M,N) + max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N)))
        allocate(WORK(lwork))
        allocate(U(m,m),VT(n,n),SI(min(m,n)),sigma(m,n))
        ! svd
        call dgesvd('A','A',m,n,Acp,m,SI,U,m,VT,n,WORK,lwork,info)

        sigma=0d0
        do i=1,size(si)
            sigma(i,i)=si(i)
        end do

        !! verif
        !    allocate(mat_verif(m,n))
        !    mat_verif=matmul(U,matmul(sigma,VT))
        !   do i=1,m
        !     do j=1,n
        !      print*,i,j,Amat(i-1,j-1),mat_verif(i,j)
        !     end do
        !   end do

        !! list of singular values and ratio to the highest
        open(10,file="singval.out",status="replace",action="write")
        write(10,*) "-----------------------------------------"
        write(10,*) "Singular values and ratio to the highest:"
        write(10,*) "----------------------------------------"
        do i=1,size(si)
            write(10,*) si(i),si(i)/maxval(si)
        end do
        close(10)

        where(abs(si/maxval(si))<1d-10)
            si = 0d0
        elsewhere
            si = 1d0/si
        end where

        sigma = 0d0
        do i=1,size(si)
            sigma(i,i) = si(i)
        end do

        !- calculation of "inverse"
        Dmat = MATMUL(TRANSPOSE(sigma),TRANSPOSE(U))

        Atrans = MATMUL(TRANSPOSE(VT),Dmat)

        !- coefficients for the correctors
        model_val(0:) = MATMUL(Atrans(:,:),dataval(0:))

        deallocate(WORK,VT,U,SI,sigma,Atrans,Acp,Dmat)
        !---
    end subroutine cg_inv

end module mleastsq
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
