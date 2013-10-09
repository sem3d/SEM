

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
        call DGEMM2(ngllx,nglly*ngllz,ngllx,htprimex,Scalp(:,:,:),dScalp_dxi)
        ! d(Scalp)_deta
        do n_z = 0,ngllz-1
            call DGEMM2(ngllx,nglly,nglly,Scalp(:,:,n_z),hprimey,dScalp_deta(:,:,n_z))
        enddo
        ! d(Scalp)_dzeta
        call DGEMM2(ngllx*nglly,ngllz,ngllz,Scalp(:,:,:),hprimez,dScalp_dzeta)


    end subroutine elem_part_deriv
#endif
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
    subroutine physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
                        Scalp,dScalp_dx,dScalp_dy,dScalp_dz)
    !- partial derivatives in the (x,y,z) space, of the scalar Scalp
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dScalp_dx,dScalp_dy,dScalp_dz
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: dScalp_dxi,dScalp_deta,dScalp_dzeta
  
    !- in the reference element
        call elem_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,Scalp,    &
                             dScalp_dxi,dScalp_deta,dScalp_dzeta)

    !- in the physical domain
    ! dScalp_dX
        dScalp_dx(:,:,:) = dScalp_dxi(:,:,:)*InvGrad(:,:,:,0,0) + dScalp_deta(:,:,:)*InvGrad(:,:,:,0,1) +  &
                  dScalp_dzeta(:,:,:)*InvGrad(:,:,:,0,2)
    ! dScalp_dY
        dScalp_dy(:,:,:) = dScalp_dxi(:,:,:)*InvGrad(:,:,:,1,0) + dScalp_deta(:,:,:)*InvGrad(:,:,:,1,1) +  &
                  dScalp_dzeta(:,:,:)*InvGrad(:,:,:,1,2)
    ! dScalp_dZ
        dScalp_dz(:,:,:) = dScalp_dxi(:,:,:)*InvGrad(:,:,:,2,0) + dScalp_deta(:,:,:)*InvGrad(:,:,:,2,1) +  &
                  dScalp_dzeta(:,:,:)*InvGrad(:,:,:,2,2)



    end subroutine physical_part_deriv
