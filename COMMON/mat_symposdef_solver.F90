!>
!!\file mat_symposdef_solver.F90
!!\brief Algorithmes de resolutions de systemes matriciels
!!\version 1.0
!!\date 05/10/2014
!! Ces routines ont ete copiees de la librairie LAPACK, version 3.5
!! Les deux routines du module ci dessous qui seront appelees dans notre code seront :
!! SPPTRF (pour la factorisation de Cholesky), et SPPTRS (pour la resolution)
!<

module smat_solver

implicit none

contains


      SUBROUTINE SPPTRF( UPLO, N, AP, INFO )
!
!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, N
!     ..
!     .. Array Arguments ..
      REAL               AP( * )
!     ..
!
!  Purpose
!  =======
!
!  SPPTRF computes the Cholesky factorization of a real symmetric
!  positive definite matrix A stored in packed format.
!
!  The factorization has the form
!     A = U**T * U,  if UPLO = 'U', or
!     A = L  * L**T,  if UPLO = 'L',
!  where U is an upper triangular matrix and L is lower triangular.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER!1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  AP      (input/output) REAL array, dimension (N*(N+1)/2)
!          On entry, the upper or lower triangle of the symmetric matrix
!          A, packed columnwise in a linear array.  The j-th column of A
!          is stored in the array AP as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!          See below for further details.
!
!          On exit, if INFO = 0, the triangular factor U or L from the
!          Cholesky factorization A = U**T*U or A = L*L**T, in the same
!          storage format as A.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the leading minor of order i is not
!                positive definite, and the factorization could not be
!                completed.
!
!  Further Details
!  ======= =======
!
!  The packed storage scheme is illustrated by the following example
!  when N = 4, UPLO = 'U':
!
!  Two-dimensional storage of the symmetric matrix A:
!
!     a11 a12 a13 a14
!         a22 a23 a24
!             a33 a34     (aij = aji)
!                 a44
!
!  Packed storage of the upper triangle of A:
!
!  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
!
!  =====================================================================
!
!     .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            J, JC, JJ
      REAL               AJJ
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
      REAL               sdot_lapack
!      EXTERNAL           LSAME, SDOT
!     ..
!     .. External Subroutines ..
      !EXTERNAL           SSCAL, SSPR, STPSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0

      IF( UPLO .EQ. 'U') THEN
         UPPER = .TRUE.
     else
         UPPER = .FALSE.
         INFO = -1
     endif
!
!     Quick return if possible
!
      IF( N.EQ.0 )  RETURN
!
      IF( UPPER ) THEN
!
!        Compute the Cholesky factorization A = U!!T!U.
!
         JJ = 0
         do J = 1, N
            JC = JJ + 1
            JJ = JJ + J
!
!           Compute elements 1:J-1 of column J.
!
!            IF( J.GT.1 )    CALL STPSV( 'Upper', 'Transpose', 'Non-unit', J-1, AP, AP( JC ), 1 )
            IF( J.GT.1 )    CALL STPSV( 'Upper', 'T', 'Non-unit', J-1, AP, AP( JC ), 1 )
!
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            AJJ = AP( JJ ) - sdot_lapack( J-1, AP( JC ), 1, AP( JC ), 1 )
            IF( AJJ.LE.ZERO ) THEN
               AP( JJ ) = AJJ
               GO TO 30
            END IF
            AP( JJ ) = SQRT( AJJ )
        enddo
      ELSE
!
!        Compute the Cholesky factorization A = L*L**T.
!
!         JJ = 1
!         DO  J = 1, N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
!            AJJ = AP( JJ )
!            IF( AJJ.LE.ZERO ) THEN
!               AP( JJ ) = AJJ
!               GO TO 30
!            END IF
!            AJJ = SQRT( AJJ )
!            AP( JJ ) = AJJ
!
!           Compute elements J+1:N of column J and update the trailing
!           submatrix.
!
!            IF( J.LT.N ) THEN
!               CALL SSCAL( N-J, ONE / AJJ, AP( JJ+1 ), 1 )
!               CALL SSPR( 'Lower', N-J, -ONE, AP( JJ+1 ), 1, AP( JJ+N-J+1 ) )
!               JJ = JJ + N - J + 1
!            END IF
!        enddo
      END IF
      GO TO 40
!
   30 CONTINUE
      INFO = J
!
   40 CONTINUE
      RETURN
!
!     End of SPPTRF
!
  END SUBROUTINE SPPTRF



      SUBROUTINE SPPTRS( UPLO, N, NRHS, AP, B, LDB, INFO )

!  -- LAPACK routine (version 3.3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!  -- April 2011                                                      --
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDB, N, NRHS
!     ..
!     .. Array Arguments ..
      REAL               AP( * ), B( LDB, * )
!     ..
!
!  Purpose
!  =======
!
!  SPPTRS solves a system of linear equations A*X = B with a symmetric
!  positive definite matrix A in packed storage using the Cholesky
!  factorization A = U**T*U or A = L*L**T computed by SPPTRF.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER!1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.
!
!  AP      (input) REAL array, dimension (N*(N+1)/2)
!          The triangular factor U or L from the Cholesky factorization
!          A = U**T*U or A = L*L**T, packed columnwise in a linear
!          array.  The j-th column of U or L is stored in the array AP
!          as follows:
!          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
!
!  B       (input/output) REAL array, dimension (LDB,NRHS)
!          On entry, the right hand side matrix B.
!          On exit, the solution matrix X.
!
!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
      !EXTERNAL           STPSV, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0

      IF( UPLO .EQ. 'U') THEN
         UPPER = .TRUE.
     else
         UPPER = .FALSE.
         INFO = -1
     endif
!
!     Quick return if possible
!
      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN
!
      IF( UPPER ) THEN
!
!        Solve A*X = B where A = U**T * U.
!
         DO 10 I = 1, NRHS
!
!           Solve U**T *X = B, overwriting B with X.
!
!            CALL STPSV( 'Upper', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
            CALL STPSV( 'Upper', 'T', 'Non-unit', N, AP, B( 1, I ), 1 )
!
!           Solve U*X = B, overwriting B with X.
!
!            CALL STPSV( 'Upper', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
            CALL STPSV( 'Upper', 'N', 'Non-unit', N, AP, B( 1, I ), 1 )
   10    CONTINUE
      ELSE
!
!        Solve A*X = B where A = L * L**T.
!
         DO 20 I = 1, NRHS
!
!           Solve L*Y = B, overwriting B with X.
!
            CALL STPSV( 'Lower', 'No transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
!
!           Solve L**T *X = Y, overwriting B with X.
!
            CALL STPSV( 'Lower', 'Transpose', 'Non-unit', N, AP, B( 1, I ), 1 )
   20    CONTINUE
      END IF
!
      RETURN
!
!     End of SPPTRS
!
      END SUBROUTINE SPPTRS

!> \brief \b STPSV
!
! =========== DOCUMENTATION ===========
!
! Online html documentation available at
! http://www.netlib.org/lapack/explore-html/
!
! Definition:
! ===========

! SUBROUTINE STPSV(UPLO,TRANS,DIAG,N,AP,X,INCX)
!
! .. Scalar Arguments ..
! INTEGER INCX,N
! CHARACTER DIAG,TRANS,UPLO
! ..
! .. Array Arguments ..
! REAL AP(*),X(*)
! ..


!> \par Purpose:
! =============
!>
!> \verbatim
!>
!> STPSV solves one of the systems of equations
!>
!> A*x = b, or A**T*x = b,
!>
!> where b and x are n element vectors and A is an n by n unit, or
!> non-unit, upper or lower triangular matrix, supplied in packed form.
!>
!> No test for singularity or near-singularity is included in this
!> routine. Such tests must be performed before calling this routine.
!> \endverbatim

! Arguments:
! ==========

!> \param[in] UPLO
!> \verbatim
!> UPLO is CHARACTER*1
!> On entry, UPLO specifies whether the matrix is an upper or
!> lower triangular matrix as follows:
!>
!> UPLO = 'U' or 'u' A is an upper triangular matrix.
!>
!> UPLO = 'L' or 'l' A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!> TRANS is CHARACTER*1
!> On entry, TRANS specifies the equations to be solved as
!> follows:
!>
!> TRANS = 'N' or 'n' A*x = b.
!>
!> TRANS = 'T' or 't' A**T*x = b.
!>
!> TRANS = 'C' or 'c' A**T*x = b.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!> DIAG is CHARACTER*1
!> On entry, DIAG specifies whether or not A is unit
!> triangular as follows:
!>
!> DIAG = 'U' or 'u' A is assumed to be unit triangular.
!>
!> DIAG = 'N' or 'n' A is not assumed to be unit
!> triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!> N is INTEGER
!> On entry, N specifies the order of the matrix A.
!> N must be at least zero.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!> AP is REAL array of DIMENSION at least
!> ( ( n*( n + 1 ) )/2 ).
!> Before entry with UPLO = 'U' or 'u', the array AP must
!> contain the upper triangular matrix packed sequentially,
!> column by column, so that AP( 1 ) contains a( 1, 1 ),
!> AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
!> respectively, and so on.
!> Before entry with UPLO = 'L' or 'l', the array AP must
!> contain the lower triangular matrix packed sequentially,
!> column by column, so that AP( 1 ) contains a( 1, 1 ),
!> AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
!> respectively, and so on.
!> Note that when DIAG = 'U' or 'u', the diagonal elements of
!> A are not referenced, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!> X is REAL array of dimension at least
!> ( 1 + ( n - 1 )*abs( INCX ) ).
!> Before entry, the incremented array X must contain the n
!> element right-hand side vector b. On exit, X is overwritten
!> with the solution vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!> INCX is INTEGER
!> On entry, INCX specifies the increment for the elements of
!> X. INCX must not be zero.
!> \endverbatim
!*
!* Authors:
!* ========
!*
!*> \author Univ. of Tennessee
!*> \author Univ. of California Berkeley
!*> \author Univ. of Colorado Denver
!*> \author NAG Ltd.
!*
!*> \date November 2011
!*
!*> \ingroup single_blas_level2
!*
!*> \par Further Details:
!* =====================
!*>
!*> \verbatim
!*>
!*> Level 2 Blas routine.
!*>
!*> -- Written on 22-October-1986.
!*> Jack Dongarra, Argonne National Lab.
!*> Jeremy Du Croz, Nag Central Office.
!*> Sven Hammarling, Nag Central Office.
!*> Richard Hanson, Sandia National Labs.
!*> \endverbatim
!*>
!* =====================================================================
 SUBROUTINE stpsv(UPLO,TRANS,DIAG,N,AP,X,INCX)
!*
!* -- Reference BLAS level2 routine (version 3.4.0) --
!* -- Reference BLAS is a software package provided by Univ. of Tennessee, --
!* -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!* November 2011
!*
!* .. Scalar Arguments ..
 INTEGER incx,n
 CHARACTER diag,trans,uplo
!* ..
!* .. Array Arguments ..
 REAL ap(*),x(*)
!* ..
!*
!* =====================================================================
!*
!* .. Parameters ..
  REAL zero
  parameter(zero=0.0e+0)
 ! ..
 ! .. Local Scalars ..
  REAL temp
  INTEGER i,info,ix,j,jx,k,kk,kx
  LOGICAL nounit


 ! Quick return if possible.
 
  IF (n.EQ.0) RETURN
 
  !nounit = lsame(diag,'N')
  nounit = .true.

 ! Set up the start point in X if the increment is not unity. This
 ! will be ( N - 1 )*INCX too small for descending loops.
 
  IF (incx.LE.0) THEN
      kx = 1 - (n-1)*incx
  ELSE IF (incx.NE.1) THEN
      kx = 1
  END IF
 
 ! Start the operations. In this version the elements of AP are
 !accessed sequentially with one pass through AP.
 !
  IF (trans .EQ. 'N') THEN
 !
 ! Form x := inv( A )*x.
 !
  IF (uplo .EQ. 'Upper') THEN
  kk = (n* (n+1))/2
  IF (incx.EQ.1) THEN
  DO 20 j = n,1,-1
  IF (x(j).NE.zero) THEN
  IF (nounit) x(j) = x(j)/ap(kk)
  temp = x(j)
  k = kk - 1
  DO 10 i = j - 1,1,-1
  x(i) = x(i) - temp*ap(k)
  k = k - 1
  10 CONTINUE
  END IF
  kk = kk - j
  20 CONTINUE
  ELSE
  jx = kx + (n-1)*incx
  DO 40 j = n,1,-1
  IF (x(jx).NE.zero) THEN
  IF (nounit) x(jx) = x(jx)/ap(kk)
  temp = x(jx)
  ix = jx
  DO 30 k = kk - 1,kk - j + 1,-1
  ix = ix - incx
  x(ix) = x(ix) - temp*ap(k)
  30 CONTINUE
  END IF
  jx = jx - incx
  kk = kk - j
  40 CONTINUE
  END IF
  ELSE
  kk = 1
  IF (incx.EQ.1) THEN
  DO 60 j = 1,n
  IF (x(j).NE.zero) THEN
  IF (nounit) x(j) = x(j)/ap(kk)
  temp = x(j)
  k = kk + 1
  DO 50 i = j + 1,n
  x(i) = x(i) - temp*ap(k)
  k = k + 1
  50 CONTINUE
  END IF
  kk = kk + (n-j+1)
  60 CONTINUE
  ELSE
  jx = kx
  DO 80 j = 1,n
  IF (x(jx).NE.zero) THEN
  IF (nounit) x(jx) = x(jx)/ap(kk)
  temp = x(jx)
  ix = jx
  DO 70 k = kk + 1,kk + n - j
  ix = ix + incx
  x(ix) = x(ix) - temp*ap(k)
  70 CONTINUE
  END IF
  jx = jx + incx
  kk = kk + (n-j+1)
  80 CONTINUE
  END IF
  END IF
  ELSE
 !
 ! Form x := inv( A**T )*x.
 !
  IF (uplo .EQ. 'Upper') THEN
  kk = 1
  IF (incx.EQ.1) THEN
  DO 100 j = 1,n
  temp = x(j)
  k = kk
  DO 90 i = 1,j - 1
  temp = temp - ap(k)*x(i)
  k = k + 1
  90 CONTINUE
  IF (nounit) temp = temp/ap(kk+j-1)
  x(j) = temp
  kk = kk + j
  100 CONTINUE
  ELSE
  jx = kx
  DO 120 j = 1,n
  temp = x(jx)
  ix = kx
  DO 110 k = kk,kk + j - 2
  temp = temp - ap(k)*x(ix)
  ix = ix + incx
  110 CONTINUE
  IF (nounit) temp = temp/ap(kk+j-1)
  x(jx) = temp
  jx = jx + incx
  kk = kk + j
  120 CONTINUE
  END IF
  ELSE
  kk = (n* (n+1))/2
  IF (incx.EQ.1) THEN
  DO 140 j = n,1,-1
  temp = x(j)
  k = kk
  DO 130 i = n,j + 1,-1
  temp = temp - ap(k)*x(i)
  k = k - 1
  130 CONTINUE
  IF (nounit) temp = temp/ap(kk-n+j)
  x(j) = temp
  kk = kk - (n-j+1)
  140 CONTINUE
  ELSE
  kx = kx + (n-1)*incx
  jx = kx
  DO 160 j = n,1,-1
  temp = x(jx)
  ix = kx
  DO 150 k = kk,kk - (n- (j+1)),-1
  temp = temp - ap(k)*x(ix)
  ix = ix - incx
  150 CONTINUE
  IF (nounit) temp = temp/ap(kk-n+j)
  x(jx) = temp
  jx = jx - incx
  kk = kk - (n-j+1)
  160 CONTINUE
  END IF
  END IF
  END IF
 
  RETURN
 
 ! End of STPSV .
 
END SUBROUTINE stpsv

end module smat_solver

 REAL FUNCTION sdot_lapack(N,SX,INCX,SY,INCY)

! -- Reference BLAS level1 routine (version 3.4.0) --
! -- Reference BLAS is a software package provided by Univ. of Tennessee, --
! -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
! November 2011
!
! .. Scalar Arguments ..
 INTEGER incx,incy,n
! ..
! .. Array Arguments ..
 REAL sx(*),sy(*)
! ..
! =====================================================================
!
! .. Local Scalars ..
 REAL stemp
 INTEGER i,ix,iy,m,mp1
! ..
! .. Intrinsic Functions ..
INTRINSIC mod
! ..
 stemp = 0.0e0
 sdot_lapack = 0.0e0
 IF (n.LE.0) RETURN
 IF (incx.EQ.1 .AND. incy.EQ.1) THEN

! code for both increments equal to 1


! clean-up loop
!
 m = mod(n,5)
 IF (m.NE.0) THEN
     DO i = 1,m
         stemp = stemp + sx(i)*sy(i)
     END DO
     IF (n.LT.5) THEN
         sdot_lapack=stemp
         RETURN
     END IF
 END IF
 mp1 = m + 1
 DO i = mp1,n,5
     stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) + &
         sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
 END DO
  ELSE
 
 ! code for unequal increments or equal increments
 ! not equal to 1
 
      ix = 1
      iy = 1
      IF (incx.LT.0) ix = (-n+1)*incx + 1
      IF (incy.LT.0) iy = (-n+1)*incy + 1
      DO i = 1,n
          stemp = stemp + sx(ix)*sy(iy)
          ix = ix + incx
          iy = iy + incy
      END DO
  END IF
  sdot_lapack = stemp
  RETURN
END FUNCTION sdot_lapack



