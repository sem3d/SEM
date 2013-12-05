

SUBROUTINE DGEMM2 (M, N, K, A, B, C)
    !     .. Scalar Arguments ..
    INTEGER, INTENT(IN)             :: M, N, K
    !     .. Array Arguments ..
    DOUBLE PRECISION, INTENT(IN)    :: A(M,K), B(K,N)
    DOUBLE PRECISION, INTENT(INOUT) :: C(M,N)
    !     ..
    !
    !  Purpose
    !  =======
    !
    !  DGEMM2  performs one of the matrix-matrix operations
    !
    !     C := A * B
    !
    !  Parameters
    !  ==========
    !
    !
    !  M      - INTEGER.
    !           On entry,  M  specifies  the number  of rows  of the  matrix
    !           op( A )  and of the  matrix  C.  M  must  be at least  zero.
    !           Unchanged on exit.
    !
    !  N      - INTEGER.
    !           On entry,  N  specifies the number  of columns of the matrix
    !           op( B ) and the number of columns of the matrix C. N must be
    !           at least zero.
    !           Unchanged on exit.
    !
    !  K      - INTEGER.
    !           On entry,  K  specifies  the number of columns of the matrix
    !           op( A ) and the number of rows of the matrix op( B ). K must
    !           be at least  zero.
    !           Unchanged on exit.
    !
    !  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
    !           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
    !           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
    !           part of the array  A  must contain the matrix  A,  otherwise
    !           the leading  k by m  part of the array  A  must contain  the
    !           matrix A.
    !           Unchanged on exit.
    !
    !  LDA    - INTEGER.
    !           On entry, LDA specifies the first dimension of A as declared
    !           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
    !           LDA must be at least  max( 1, m ), otherwise  LDA must be at
    !           least  max( 1, k ).
    !           Unchanged on exit.
    !
    !  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
    !           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
    !           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
    !           part of the array  B  must contain the matrix  B,  otherwise
    !           the leading  n by k  part of the array  B  must contain  the
    !           matrix B.
    !           Unchanged on exit.
    !
    !  LDB    - INTEGER.
    !           On entry, LDB specifies the first dimension of B as declared
    !           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
    !           LDB must be at least  max( 1, k ), otherwise  LDB must be at
    !           least  max( 1, n ).
    !           Unchanged on exit.
    !
    !  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
    !           Before entry, the leading  m by n  part of the array  C must
    !           contain the matrix  C,  except when  beta  is zero, in which
    !           case C need not be set on entry.
    !           On exit, the array  C  is overwritten by the  m by n  matrix
    !           ( alpha*op( A )*op( B ) + beta*C ).
    !
    !  LDC    - INTEGER.
    !           On entry, LDC specifies the first dimension of C as declared
    !           in  the  calling  (sub)  program.   LDC  must  be  at  least
    !           max( 1, m ).
    !           Unchanged on exit.
    !
    !
    !  Level 3 Blas routine.
    !
    !  -- Written on 8-February-1989.
    !     Jack Dongarra, Argonne National Laboratory.
    !     Iain Duff, AERE Harwell.
    !     Jeremy Du Croz, Numerical Algorithms Group Ltd.
    !     Sven Hammarling, Numerical Algorithms Group Ltd.
    !
    !
    !     .. Intrinsic Functions ..
    INTRINSIC          MAX
    !     .. Local Scalars ..
    INTEGER            :: I, J, L
    DOUBLE PRECISION   :: TEMP
    !     .. Parameters ..
    DOUBLE PRECISION, PARAMETER   :: ZERO = 0.0D+0
    !
    IF( ( M.EQ.0 ).OR.( N.EQ.0 ) ) RETURN
    !
    !     Start the operations.
    !
    !
    !           Form  C := alpha*A*B + beta*C.
    !
    DO J = 1, N
        DO I = 1, M
            C( I, J ) = ZERO
        END DO
        DO L = 1, K
            TEMP = B( L, J )
            DO I = 1, M
                C( I, J ) = C( I, J ) + TEMP*A( I, L )
            END DO
        END DO
    END DO
    !
    RETURN
    !
    !     End of DGEMM .
    !
END SUBROUTINE DGEMM2

#if 1
    subroutine elem_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,Scalp,    &
                      dScalp_dxi,dScalp_deta,dScalp_dzeta)
    !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dScalp_dxi,dScalp_deta,dScalp_dzeta
        integer :: i,j,k,l
        DOUBLE PRECISION, PARAMETER   :: ZERO = 0.0D+0
        DO K = 0, ngllz-1
            ! d(Scalp)_dxi
            DO J = 0, nglly-1
                DO I = 0, ngllx-1
                    dScalp_dxi(I,J,K) = ZERO
                END DO
                DO L = 0, ngllx-1
                    DO I = 0, ngllx-1
                        dScalp_dxi(I,J,K) = dScalp_dxi(I,J,K) + Scalp(L,J,K)*hTprimex(I,L)
                    END DO
                END DO
            END DO

            ! d(Scalp)_deta
            DO J = 0, nglly-1
                DO I = 0, ngllx-1
                    dScalp_deta( I, J, K ) = ZERO
                END DO
                DO L = 0, nglly-1
                    DO I = 0, ngllx-1
                        Dscalp_deta( I, J, K ) = dScalp_deta( I, J, K ) + Scalp( I, L, K )*hprimey(L,J)
                    END DO
                END DO
            END DO

            ! d(Scalp)_dzeta
            DO J = 0,nglly-1
                DO I = 0,ngllx-1
                    dScalp_dzeta( I, J, K ) = ZERO
                END DO
            END DO
            DO L = 0,ngllz-1
                DO J = 0,nglly-1
                    DO I = 0,ngllx-1
                        dScalp_dzeta(I,J,K) = dScalp_dzeta(I,J,K) + Scalp(I,J,L)*hprimez(L,K)
                    END DO
                END DO
            END DO
        END DO
    end subroutine elem_part_deriv
#else
    subroutine elem_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,Scalp,    &
                      dScalp_dxi,dScalp_deta,dScalp_dzeta)
    !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dScalp_dxi,dScalp_deta,dScalp_dzeta
        integer   :: n_z

        ! d(Scalp)_dxi
        call DGEMM('N','N',ngllx,nglly*ngllz,ngllx,1.0,htprimex,ngllx,Scalp(:,:,:),ngllx,0.0,dScalp_dxi,ngllx)
        ! d(Scalp)_deta
        do n_z = 0,ngllz-1
            call DGEMM('N','N',ngllx,nglly,nglly,1.0,Scalp(:,:,n_z),ngllx,hprimey,nglly,0.0,dScalp_deta(:,:,n_z),ngllx)
        enddo
        ! d(Scalp)_dzeta
        call DGEMM('N','N',ngllx*nglly,ngllz,ngllz,1.0,Scalp(:,:,:),ngllx*nglly,hprimez,ngllz,0.0,dScalp_dzeta,ngllx*nglly)


    end subroutine elem_part_deriv
#endif
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

#if 1
    subroutine physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad, &
        Scalp, dS_dx,dS_dy,dS_dz)
        !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dS_dx,dS_dy,dS_dz

        if (ngllx==5 .and. nglly==5 .and. ngllz==5) then
            call physical_part_deriv_555(hTprimex,hprimey,hprimez,InvGrad, &
                Scalp, dS_dx,dS_dy,dS_dz)
            return
        else if (ngllx==7 .and. nglly==7 .and. ngllz==7) then
            call physical_part_deriv_777(hTprimex,hprimey,hprimez,InvGrad, &
                Scalp, dS_dx,dS_dy,dS_dz)
            return
        endif
        call physical_part_deriv_nnn(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad, &
            Scalp, dS_dx,dS_dy,dS_dz)
    end subroutine physical_part_deriv

    subroutine physical_part_deriv_nnn(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad, &
        Scalp, dS_dx,dS_dy,dS_dz)
    !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dS_dx,dS_dy,dS_dz
        real :: dS_dxi, dS_deta, dS_dzeta
        integer :: i,j,k,l
        DOUBLE PRECISION, PARAMETER   :: ZERO = 0.0D+0

        DO K = 0, ngllz-1
            ! d(Scalp)_dxi
            DO J = 0, nglly-1
                DO I = 0, ngllx-1
                    dS_dxi = ZERO
                    dS_deta = ZERO
                    dS_dzeta = ZERO
                    DO L = 0, ngllx-1
                        dS_dxi   = dS_dxi   + Scalp(L,J,K)*hTprimex(I,L)
                        DS_deta  = dS_deta  + Scalp(I,L,K)*hprimey(L,J)
                        dS_dzeta = dS_dzeta + Scalp(I,J,L)*hprimez(L,K)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dxi*InvGrad(I,J,K,0,0) + dS_deta*InvGrad(I,J,K,0,1) +  &
                        dS_dzeta*InvGrad(I,J,K,0,2)
                    dS_dy(I,J,K) = dS_dxi*InvGrad(I,J,K,1,0) + dS_deta*InvGrad(I,J,K,1,1) +  &
                        dS_dzeta*InvGrad(I,J,K,1,2)
                    dS_dz(I,J,K) = dS_dxi*InvGrad(I,J,K,2,0) + dS_deta*InvGrad(I,J,K,2,1) +  &
                        dS_dzeta*InvGrad(I,J,K,2,2)
                END DO
            END DO
        END DO
    end subroutine physical_part_deriv_nnn

    subroutine physical_part_deriv_555(hTprimex,hprimey,hprimez,InvGrad, &
        Scalp, dS_dx,dS_dy,dS_dz)
    !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, parameter  :: ngllx=5,nglly=5,ngllz=5
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dS_dx,dS_dy,dS_dz
        real :: dS_dxi, dS_deta, dS_dzeta
        integer :: i,j,k,l
        DOUBLE PRECISION, PARAMETER   :: ZERO = 0.0D+0

        DO K = 0, ngllz-1
            ! d(Scalp)_dxi
            DO J = 0, nglly-1
                DO I = 0, ngllx-1
                    dS_dxi = ZERO
                    dS_deta = ZERO
                    dS_dzeta = ZERO
                    DO L = 0, ngllx-1
                        dS_dxi   = dS_dxi   + Scalp(L,J,K)*hTprimex(I,L)
                        DS_deta  = dS_deta  + Scalp(I,L,K)*hprimey(L,J)
                        dS_dzeta = dS_dzeta + Scalp(I,J,L)*hprimez(L,K)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dxi*InvGrad(I,J,K,0,0) + dS_deta*InvGrad(I,J,K,0,1) +  &
                        dS_dzeta*InvGrad(I,J,K,0,2)
                    dS_dy(I,J,K) = dS_dxi*InvGrad(I,J,K,1,0) + dS_deta*InvGrad(I,J,K,1,1) +  &
                        dS_dzeta*InvGrad(I,J,K,1,2)
                    dS_dz(I,J,K) = dS_dxi*InvGrad(I,J,K,2,0) + dS_deta*InvGrad(I,J,K,2,1) +  &
                        dS_dzeta*InvGrad(I,J,K,2,2)
                END DO
            END DO
        END DO
    end subroutine physical_part_deriv_555

    subroutine physical_part_deriv_777(hTprimex,hprimey,hprimez,InvGrad, &
        Scalp, dS_dx,dS_dy,dS_dz)
    !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, parameter  :: ngllx=7,nglly=7,ngllz=7
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dS_dx,dS_dy,dS_dz
        real :: dS_dxi, dS_deta, dS_dzeta
        integer :: i,j,k,l
        DOUBLE PRECISION, PARAMETER   :: ZERO = 0.0D+0

        DO K = 0, ngllz-1
            ! d(Scalp)_dxi
            DO J = 0, nglly-1
                DO I = 0, ngllx-1
                    dS_dxi = ZERO
                    dS_deta = ZERO
                    dS_dzeta = ZERO
                    DO L = 0, ngllx-1
                        dS_dxi   = dS_dxi   + Scalp(L,J,K)*hTprimex(I,L)
                        DS_deta  = dS_deta  + Scalp(I,L,K)*hprimey(L,J)
                        dS_dzeta = dS_dzeta + Scalp(I,J,L)*hprimez(L,K)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dxi*InvGrad(I,J,K,0,0) + dS_deta*InvGrad(I,J,K,0,1) +  &
                        dS_dzeta*InvGrad(I,J,K,0,2)
                    dS_dy(I,J,K) = dS_dxi*InvGrad(I,J,K,1,0) + dS_deta*InvGrad(I,J,K,1,1) +  &
                        dS_dzeta*InvGrad(I,J,K,1,2)
                    dS_dz(I,J,K) = dS_dxi*InvGrad(I,J,K,2,0) + dS_deta*InvGrad(I,J,K,2,1) +  &
                        dS_dzeta*InvGrad(I,J,K,2,2)
                END DO
            END DO
        END DO
    end subroutine physical_part_deriv_777

#else
    subroutine physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad, &
        Scalp, dS_dx,dS_dy,dS_dz)
    !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dS_dx,dS_dy,dS_dz
        real :: dS_dxi, dS_deta, dS_dzeta
        integer :: i,j,k,l
        DOUBLE PRECISION, PARAMETER   :: ZERO = 0.0D+0

        DO K = 0, ngllz-1
            ! d(Scalp)_dxi
            DO J = 0, nglly-1
                DO I = 0, ngllx-1
                    dS_dxi = ZERO
                    DO L = 0, ngllx-1
                        dS_dxi   = dS_dxi   + Scalp(L,J,K)*hTprimex(I,L)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dxi*InvGrad(I,J,K,0,0)
                    dS_dy(I,J,K) = dS_dxi*InvGrad(I,J,K,1,0)
                    dS_dz(I,J,K) = dS_dxi*InvGrad(I,J,K,2,0)
                END DO
            END DO
        END DO

        DO K = 0, ngllz-1
            ! d(Scalp)_dxi
            DO J = 0, nglly-1
                DO I = 0, ngllx-1
                    dS_deta = ZERO
                    DO L = 0, ngllx-1
                        dS_deta  = dS_deta  + Scalp(I,L,K)*hprimey(L,J)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dx(I,J,K) + dS_deta*InvGrad(I,J,K,0,1)
                    dS_dy(I,J,K) = dS_dy(I,J,K) + dS_deta*InvGrad(I,J,K,1,1)
                    dS_dz(I,J,K) = dS_dz(I,J,K) + dS_deta*InvGrad(I,J,K,2,1)
                END DO
            END DO
        END DO

        DO K = 0, ngllz-1
            ! d(Scalp)_dxi
            DO J = 0, nglly-1
                DO I = 0, ngllx-1
                    dS_dzeta = ZERO
                    DO L = 0, ngllx-1
                        dS_dzeta = dS_dzeta + Scalp(I,J,L)*hprimez(L,K)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dx(I,J,K) + dS_dzeta*InvGrad(I,J,K,0,2)
                    dS_dy(I,J,K) = dS_dy(I,J,K) + dS_dzeta*InvGrad(I,J,K,1,2)
                    dS_dz(I,J,K) = dS_dz(I,J,K) + dS_dzeta*InvGrad(I,J,K,2,2)
                END DO
            END DO
        END DO



    end subroutine physical_part_deriv
#endif
