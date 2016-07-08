      SUBROUTINE DGER ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
!Source: http://www.netlib.no/netlib/maspar/f90/dblas/dger.f
!***************************************************************************
!                                                                          *
!   DATA PARALLEL BLAS based on MPL                                        *
!                                                                          *
!   Version 1.0   1/9-92 ,                                                 *
!   For MasPar MP-1 computers                                              *
!                                                                          *
!   para//ab, University of Bergen, NORWAY                                 *
!                                                                          *
!   These programs must be called using F90 style array syntax.            *
!   Note that the F77 style calling sequence has been retained             *
!   in this version for compatibility reasons, be aware that               *
!   parameters related to the array dimensions and shape therefore may     *
!   be redundant and without any influence.                                *
!   The calling sequence may be changed in a future version.               *
!   Please report any BUGs, ideas for improvement or other                 *
!   comments to                                                            *
!                    adm@parallab.uib.no                                   *
!                                                                          *
!   Future versions may then reflect your suggestions.                     *
!   The most current version of this software is available                 *
!   from netlib@nac.no , send the message `send index from maspar'         *
!                                                                          *
!   REVISIONS:                                                             *
!                                                                          *
!***************************************************************************

    implicit none

!     .. Scalar Arguments ..
      DOUBLE PRECISION   :: ALPHA
      INTEGER            :: INCX, INCY, LDA, M, N
!     .. Array Arguments ..
      double precision, dimension(1:,1:) :: A
      double precision, dimension(1:) :: X
      double precision, dimension(1:) :: Y
!     ..
!
!  Purpose
!  =======
!
!  DGER   performs the rank 1 operation
!
!     A := alpha*x*y' + A,
!
!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.
!
!  Parameters
!  ==========
!
!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.
!
!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.
!
!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.
!
!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.
!
!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.
!
!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.
!
!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.
!
!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.
!
!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.
!
!
!  Level 2 Blas routine.
!
!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.
!
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
!     .. Local Arrays ..
      double precision, dimension(m) :: xloc
      double precision, dimension(n) :: yloc
!     .. Local Scalars ..
      INTEGER            INFO, KY, KX
!     .. External Subroutines ..
      EXTERNAL           XERBLA
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
      INTRINSIC          spread
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF

      print *, "Here 1"
      print *, "ALPHA = ", ALPHA
      print *, "INCX = ", INCX
      print *, "INCY = ", INCY
      print *, "LDA = ", LDA
      print *, "M = ", M
      print *, "N = ", N
      write(*,*) "A    = ", A
      write(*,*) "Y    = ", Y
      write(*,*) "X    = ", X
!
!     Quick return if possible.
!
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) ) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( M - 1 )*INCX
      END IF

      print *, "Here 2"
!
      IF( INCX.EQ.1 )THEN
       print *, "Here 1.125"
       print *, "XLOC(1:M) = ", XLOC(1:M)
       !print *, "X(:)    = ", X(:)
       !print *, "X(1:M)    = ", X(1:M)
       XLOC(1:M) = X(1:M)
      ELSE
       XLOC(1:M) = X(KX : KX+INCX*(M-1) : INCX)
      END IF
!
      print *, "Here 2.25"
      xloc = xloc * alpha
      print *, "Here 2.75"
!
      IF( INCY.EQ.1 )THEN
       YLOC(1:N) = Y(1:N)
      ELSE
       YLOC(1:N) = Y(KY : KY+INCY*(N-1) : INCY)
      END IF

      print *, "Here 3"
!
      a(1:m,1:n) = a(1:m,1:n) + spread(xloc,2,n) * spread(yloc,1,m)

      print *, "Here 4"
!
      RETURN
!
!     End of DGER  .
!
      END
