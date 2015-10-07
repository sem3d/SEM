!! The original code in this module was written by Daniele Funaro (see below)
!! It has been reworked as F90 module
MODULE splib
    USE constants, only : M_PI, M_PI_2
    IMPLICIT none

    PUBLIC ZELEGL, WELEGL, DMLEGL
CONTAINS

    !! NAG ROUTINES

    SUBROUTINE C06EBF(X, N, IFAIL)
        !C06EBF calculates the discrete Fourier transform of a Hermitian sequence of n complex data values. (No
        !extra workspace required.)

        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(OUT) :: IFAIL
        REAL(kind=8), DIMENSION(N), INTENT(INOUT) :: X
        X=N+1
        IFAIL=1
    END SUBROUTINE C06EBF

    SUBROUTINE C06EAF(X, N, IFAIL)
        !C06EAF calculates the discrete Fourier transform of a sequence of n real data values.
        !workspace required.)

        INTEGER, INTENT(IN) :: N
        INTEGER, INTENT(OUT) :: IFAIL
        REAL(kind=8), DIMENSION(N), INTENT(INOUT) :: X
        X=N+1
        IFAIL=1
    END SUBROUTINE C06EAF


    !*****************************************************************
    !*     FORTRAN ROUTINES FOR SPECTRAL METHODS
    !*     BY  DANIELE FUNARO
    !*     DEPARTMENT OF MATHEMATICS
    !*     UNIVERSITY OF MODENA
    !*     VIA CAMPI 213/B, 41100 MODENA, ITALY
    !*     E-MAIL: funaro@unimo.it
    !*****************************************************************
    SUBROUTINE GAMMAF(X,GX)
        !*****************************************************************
        !*     COMPUTES THE GAMMA FUNCTION AT A GIVEN POINT
        !*     X = ARGUMENT GREATER THAN 1.E-75 AND SMALLER THAN 57.
        !*     GX= VALUE OF GAMMA IN X
        !*****************************************************************
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), INTENT(OUT) :: GX
        INTEGER :: IND, K
        REAL(kind=8), parameter :: EPS = 1.D-14
        REAL(kind=8) :: XX, PR, S, G
        REAL(kind=8), parameter, dimension(11) :: &
            C = (/-0.524741987629368444D0, 0.116154405493589130D0, &
            -0.765978624506602380D-2, 0.899719449391378898D-4, &
            -0.194536980009534621D-7, 0.199382839513630987D-10, &
            -0.204209590209541319D-11, 0.863896817907000175D-13, &
            0.152237501608472336D-13, -0.82572517527771995D-14, &
            0.29973478220522461D-14 /)

        IF (X .LE. 1.D-75 .OR. X .GE. 57.D0) THEN
            WRITE(*,*) 'ARGUMENT OUT OF RANGE IN SUBROUTINE GAMMAF'
            RETURN
        ENDIF

        XX  = X
        GX  = 1.0D0

        IF (ABS(XX-1.D0) .LT. EPS) RETURN

        DO WHILE (XX .GE. 1.D0)
            XX = XX-1.D0
            GX = GX*XX
            IF (ABS(XX-1.D0) .LT. EPS) RETURN
        ENDDO

        IND = 0
        IF (XX .LT. .5D0) THEN
            IND = 1
            GX  = GX*M_PI/SIN(M_PI*XX)
            XX  = 1.D0-XX
        ENDIF

        PR = 1.D0
        S  = 0.426401432711220868D0

        DO K=1,11
            PR = PR*(XX-real(K,kind=8))/(XX+real(K-1,kind=8))
            S  = S+C(K)*PR
        END DO

        G  = S*DEXP(1.D0-XX)*(XX+4.5D0)**(XX-.5D0)
        IF (IND .EQ. 1) THEN
            GX = GX/G
        ELSE
            GX = GX*G
        ENDIF

        RETURN
    END SUBROUTINE GAMMAF


    SUBROUTINE VAJAPO(N,A,B,X,Y,DY,D2Y)
        !************************************************************
        !*   COMPUTES THE VALUE OF THE JACOBI POLYNOMIAL OF DEGREE N
        !*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
        !*   N  = DEGREE OF THE POLYNOMIAL
        !*   A  = PARAMETER >-1
        !*   B  = PARAMETER >-1
        !*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B, X
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER :: I
        REAL(kind=8) :: AB, YP, DYP, D2YP, DI, C0, C1, C2, C3, C4, YM, DYM, D2YM

        Y   = 1.D0
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .EQ. 0) RETURN

        AB  = A+B
        Y   = .5D0*(AB+2.D0)*X+.5D0*(A-B)
        DY  = .5D0*(AB+2.D0)
        D2Y = 0.D0
        IF(N .EQ. 1) RETURN

        YP   = 1.D0
        DYP  = 0.D0
        D2YP = 0.D0
        DO I=2,N
            DI = real(I,kind=8)
            C0 = 2.D0*DI+AB
            C1 = 2.D0*DI*(DI+AB)*(C0-2.D0)
            C2 = (C0-1.D0)*(C0-2.D0)*C0
            C3 = (C0-1.D0)*(A-B)*AB
            C4 = 2.D0*(DI+A-1.D0)*C0*(DI+B-1.D0)
            YM = Y
            Y  = ((C2*X+C3)*Y-C4*YP)/C1
            YP = YM
            DYM  = DY
            DY   = ((C2*X+C3)*DY-C4*DYP+C2*YP)/C1
            DYP  = DYM
            D2YM = D2Y
            D2Y  = ((C2*X+C3)*D2Y-C4*D2YP+2.D0*C2*DYP)/C1
            D2YP = D2YM
        END DO
        RETURN
    END SUBROUTINE VAJAPO


    SUBROUTINE VALEPO(N,X,Y,DY,D2Y)
        !**************************************************************
        !*   COMPUTES THE VALUE OF THE LEGENDRE POLYNOMIAL OF DEGREE N
        !*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
        !*   N  = DEGREE OF THE POLYNOMIAL
        !*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !**************************************************************
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER, INTENT(IN) :: N
        INTEGER :: I
        REAL(kind=8) :: C1, C2, C4, YP, DYP, D2YP, YM, DYM, D2YM

        Y   = 1.D0
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .EQ. 0) RETURN

        Y   = X
        DY  = 1.D0
        D2Y = 0.D0
        IF(N .EQ. 1) RETURN

        YP   = 1.D0
        DYP  = 0.D0
        D2YP = 0.D0
        DO I=2,N
            C1 = real(I,kind=8)
            C2 = 2.D0*C1-1.D0
            C4 = C1-1.D0
            YM = Y
            Y  = (C2*X*Y-C4*YP)/C1
            YP = YM
            DYM  = DY
            DY   = (C2*X*DY-C4*DYP+C2*YP)/C1
            DYP  = DYM
            D2YM = D2Y
            D2Y  = (C2*X*D2Y-C4*D2YP+2.D0*C2*DYP)/C1
            D2YP = D2YM
        END DO

        RETURN
    END SUBROUTINE VALEPO


    SUBROUTINE VACHPO(N,X,Y,DY,D2Y)
        !***************************************************************
        !*   COMPUTES THE VALUE OF THE CHEBYSHEV POLYNOMIAL OF DEGREE N
        !*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
        !*   N  = DEGREE OF THE POLYNOMIAL
        !*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !***************************************************************
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER, INTENT(IN) :: N
        REAL(kind=8) :: YP, DYP, DYM, D2YM, D2YP, YM
        INTEGER :: K

        Y   = 1.D0
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .EQ. 0) RETURN

        Y   = X
        DY  = 1.D0
        D2Y = 0.D0
        IF (N .EQ. 1) RETURN

        YP   = 1.D0
        DYP  = 0.D0
        D2YP = 0.D0
        DO K=2,N
            YM  = Y
            Y   = 2.D0*X*Y-YP
            YP  = YM
            DYM = DY
            DY  = 2.D0*X*DY+2.D0*YP-DYP
            DYP = DYM
            D2YM= D2Y
            D2Y = 2.D0*X*D2Y+4.D0*DYP-D2YP
            D2YP= D2YM
        END DO

        RETURN
    END SUBROUTINE VACHPO


    SUBROUTINE VALAPO(N,A,X,Y,DY,D2Y)
        !**************************************************************
        !*   COMPUTES THE VALUE OF THE LAGUERRE POLYNOMIAL OF DEGREE N
        !*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
        !*   N  = DEGREE OF THE POLYNOMIAL
        !*   A  = PARAMETER >-1
        !*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !**************************************************************
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: X, A
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER, INTENT(IN) :: N
        REAL(kind=8) :: B1, B2, YP, DYP, DYM, D2YM, D2YP, DK, YM
        INTEGER :: K

        Y   = 1.D0
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .EQ. 0) RETURN

        Y   = 1.D0+A-X
        DY  = -1.D0
        D2Y = 0.D0
        IF (N .EQ. 1) RETURN

        YP  = 1.D0
        DYP = 0.D0
        D2YP= 0.D0
        DO K=2,N
            DK = real(K,kind=8)
            B1 = (2.D0*DK+A-1.D0-X)/DK
            B2 = (DK+A-1.D0)/DK
            YM = Y
            Y  = B1*Y-B2*YP
            YP = YM
            DYM = DY
            DY  = B1*DY-YP/DK-B2*DYP
            DYP = DYM
            D2YM= D2Y
            D2Y = B1*D2Y-2.D0*DYP/DK-B2*D2YP
            D2YP= D2YM
        END DO

        RETURN
    END SUBROUTINE VALAPO


    SUBROUTINE VAHEPO(N,X,Y,DY,D2Y)
        !*************************************************************
        !*   COMPUTES THE VALUE OF THE HERMITE POLYNOMIAL OF DEGREE N
        !*   AND ITS FIRST AND SECOND DERIVATIVES AT A GIVEN POINT
        !*   N  = DEGREE OF THE POLYNOMIAL
        !*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !*************************************************************
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER, INTENT(IN) :: N
        REAL(kind=8) :: DK, DN, DNN, YP, YM, YPM
        INTEGER :: K

        Y   = 1.D0
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .LE. 0) RETURN

        Y   = 2.D0*X
        DY  = 2.D0
        D2Y = 0.D0
        IF (N .EQ. 1) RETURN

        YP=1.D0
        DO K=2,N
            DK = real(K-1,kind=8)
            YM  = Y
            Y   = 2.D0*X*Y-2.D0*DK*YP
            YPM = YP
            YP  = YM
        END DO

        DN  = 2.D0*real(N,kind=8)
        DNN = 2.D0*real(N-1,kind=8)
        DY  = DN*YP
        D2Y = DN*DNN*YPM

        RETURN
    END SUBROUTINE VAHEPO


    SUBROUTINE VALASF(N,A,X,Y,DY)
        !*********************************************************************
        !*   COMPUTES THE VALUES OF THE SCALED LAGUERRE FUNCTION OF DEGREE N
        !*   AND ITS FIRST DERIVATIVE AT A GIVEN POINT
        !*   N  = DEGREE
        !*   A  = PARAMETER >-1
        !*   X  = POINT (NON NEGATIVE) IN WHICH THE COMPUTATION IS PERFORMED
        !*   Y  = VALUE OF THE FUNCTION IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*********************************************************************
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: X, A
        REAL(kind=8), INTENT(OUT) :: Y, DY
        INTEGER, INTENT(IN) :: N
        REAL(kind=8) :: DK, DK4, C0, C1, C2, C3, C4, C5, C6, C7
        REAL(kind=8) :: YP, DYP, YM, DYM
        INTEGER :: K

        Y  = 1.D0
        DY = 0.D0
        IF (N .EQ. 0) RETURN

        C0 = 1.D0/(4.D0+X)
        C1 = 4.D0*C0/(A+1.D0)
        Y  = C1*(A+1.D0-X)
        DY = -C0*C1*(A+5.D0)
        IF (N .EQ. 1) RETURN

        YP  = 1.D0
        DYP = 0.D0
        DO K=2,N
            DK  = real(K,kind=8)
            DK4 = 4.D0*DK
            C0  = 1.D0/(DK4+X)
            C1  = DK4+X-2.D0
            C2  = 1.D0/(C1-2.D0)
            C3  = DK4*C0/(DK+A)
            C4  = 2.D0*DK+A-1.D0
            C5  = C4-X
            C6  = C4+DK4
            C7  = C2*real(4*(K-1)**2,kind=8)
            DYM = DY
            DY  = C3*(C5*DY-C0*C6*Y+C7*(2.D0*C0*C1*C2*YP-DYP))
            DYP = DYM
            YM  = Y
            Y   = C3*(C5*Y-C7*YP)
            YP  = YM
        END DO

        RETURN
    END SUBROUTINE VALASF


    SUBROUTINE VAHESF(N,X,Y,DY)
        !***********************************************************
        !*   COMPUTES THE VALUES OF THE SCALED HERMITE FUNCTION OF
        !*   DEGREE N AND ITS FIRST DERIVATIVE AT A GIVEN POINT
        !*   N  = DEGREE
        !*   X  = POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   Y  = VALUE OF THE FUNCTION IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !***********************************************************
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), INTENT(OUT) :: Y, DY
        INTEGER, INTENT(IN) :: N
        INTEGER :: M
        REAL(kind=8) :: A, XX

        Y  = 1.D0
        DY = 0.D0
        IF (N .EQ. 0) RETURN

        XX = X*X
        M  = N/2
        IF(N .EQ. 2*M) THEN
            A  = -.5D0
            CALL VALASF(M,A,XX,Y,DY)
            DY = 2.D0*X*DY
        ELSE
            A  = .5D0
            CALL VALASF(M,A,XX,Y,DY)
            DY = Y+2.D0*XX*DY
            Y  = X*Y
        ENDIF

        RETURN
    END SUBROUTINE VAHESF


    SUBROUTINE ZEJAGA(N,A,B,CS,DZ)
        !***************************************************************
        !*   COMPUTES THE ZEROES OF THE JACOBI POLYNOMIAL OF DEGREE N
        !*   N  = THE NUMBER OF ZEROES
        !*   A  = PARAMETER BETWEEN -1/2 AND 1/2
        !*   B  = PARAMETER BETWEEN -1/2 AND 1/2
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !***************************************************************
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: A, B
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: CS, DZ
        INTEGER :: I, IT
        REAL(kind=8) :: AB, C, DI, CSX, Y, DY, D2Y
        REAL(kind=8), PARAMETER :: EPS=1.D-17

        IF (N .EQ. 0) RETURN

        AB = A+B
        CS(1) = (B-A)/(AB+2.D0)
        DZ(1) = .5D0*AB+1.D0
        IF(N .EQ. 1) RETURN

        C  = M_PI_2/(2.D0*real(N,kind=8)+AB+1.D0)
        DO I=1,N
            DI  = real(I,kind=8)
            CSX = -COS(C*(4.D0*DI+AB-1.D0))
            DO IT=1,8
                CALL VAJAPO(N,A,B,CSX,Y,DY,D2Y)
                IF(ABS(Y) .LT. EPS) EXIT
                CSX = CSX-Y/DY
            END DO
            IF(ABS(CSX) .LT. EPS) CSX=0.D0
            CS(I) = CSX
            DZ(I) = DY
        END DO

        RETURN
    END SUBROUTINE ZEJAGA


    SUBROUTINE ZELEGA(N,CS,DZ)
        !***************************************************************
        !*   COMPUTES THE ZEROES OF THE LEGENDRE POLYNOMIAL OF DEGREE N
        !*   N  = THE NUMBER OF ZEROES
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !***************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: CS, DZ
        INTEGER :: I, IT, N2, IN
        REAL(kind=8) :: Y, DY, D2Y, DI, CSX, C

        IF (N .EQ. 0) RETURN

        CS(1) = 0.D0
        DZ(1) = 1.D0
        IF(N .EQ. 1) RETURN

        N2 = N/2
        IN = 2*N-4*N2-1

        C  = M_PI_2/(2.D0*real(N,kind=8)+1.D0)
        DO I=1,N2
            DI  = real(I,kind=8)
            CSX = COS(C*(4.D0*DI-1.D0))
            DO IT=1,8
                CALL VALEPO(N,CSX,Y,DY,D2Y)
                CSX = CSX-Y/DY
            END DO
            CS(I) = -CSX
            CS(N-I+1) = CSX
            DZ(I) = DY*real(IN,kind=8)
            DZ(N-I+1) = DY
        END DO

        IF(IN .EQ. -1) RETURN
        CSX = 0.D0
        CS(N2+1) = CSX
        CALL VALEPO(N,CSX,Y,DY,D2Y)
        DZ(N2+1) = DY
        RETURN
    END SUBROUTINE ZELEGA


    SUBROUTINE ZECHGA(N,CS,DZ)
        !****************************************************************
        !*   COMPUTES THE ZEROES OF THE CHEBYSHEV POLYNOMIAL OF DEGREE N
        !*   N  = THE NUMBER OF ZEROES
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !****************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: CS, DZ
        INTEGER :: I, N2, IN, N4, IN2
        REAL(kind=8) :: DI, CSX, QX, SI, DN, C

        IF (N .EQ. 0) RETURN

        CS(1) = 0.D0
        DZ(1) = 1.D0
        IF(N .EQ. 1) RETURN

        N2 = N/2
        IN = 1+4*N2-2*N
        DN = real(N,kind=8)
        C  = M_PI_2/DN
        SI = -1.D0
        DO I=1,N2
            DI = real(I,kind=8)
            CSX = COS(C*(2.D0*DI-1.D0))
            CS(I) = -CSX
            CS(N-I+1) = CSX
            QX = DN/SQRT(1.D0-CSX*CSX)
            DZ(I) = QX*SI*real(IN,kind=8)
            DZ(N-I+1) = -QX*SI
            SI = -SI
        END DO

        IF(IN .EQ. 1) RETURN
        CS(N2+1) = 0.D0
        N4  = N2/2
        IN2 = 1+4*N4-2*N2
        DZ(N2+1) = DN*real(IN2,kind=8)

        RETURN
    END SUBROUTINE ZECHGA


    SUBROUTINE ZELAGA(N,A,CS,DZ)
        !************************************************************************
        !*   COMPUTES THE ZEROES OF THE LAGUERRE POLYNOMIAL OF DEGREE N
        !*   N  = THE NUMBER OF ZEROES
        !*   A  = PARAMETER >-1
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   DZ = DERIVATIVES OF THE SCALED FUNCTIONS AT THE ZEROES, DZ(I), I=1,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: CS, DZ
        INTEGER :: M, IT, IN
        REAL(kind=8) :: A1, C1, DM, C2, XN, XP, Z, ZD, CSX, CSM, DZM, DN, Y, DY

        IF (N .EQ. 0) RETURN

        A1 = A+1.D0
        CS(1) = A1
        DZ(1) = -4.D0/(A1*(A+5.D0))
        IF (N .EQ. 1) RETURN

        DN = real(N,kind=8)
        C1 = 2.D0*DN+A1
        DO M=1,N
            DM = real(M,kind=8)
            C2 = 2.D0*(DN+.75D0-DM)*M_PI/C1
            XN = (C2+M_PI)/2.D0
            DO IT=1,8
                XP = XN
                XN = (SIN(XP)-XP*COS(XP)+C2)/(1.D0-COS(XP))
            END DO
            Z  = (COS(XN/2.D0))**2
            ZD = 1.D0/(Z-1.D0)
            CSX= 2.D0*C1*Z-((1.25D0*ZD+1.D0)*ZD-1.D0+3.D0*A*A)/(6.D0*C1)
            DO IT=1,6
                CALL VALASF(N,A,CSX,Y,DY)
                CSX = CSX-Y/DY
            END DO
            CS(M) = CSX
            DZ(M) = DY
        END DO

        DO
            IN = 0
            DO M=1,N-1
                IF(CS(M) .LE. CS(M+1)) CYCLE
                CSM = CS(M)
                CS(M) = CS(M+1)
                CS(M+1) = CSM
                DZM = DZ(M)
                DZ(M) = DZ(M+1)
                DZ(M+1) = DZM
                IN = 1
            END DO
            IF (IN .NE. 1) EXIT
        END DO
        RETURN
    END SUBROUTINE ZELAGA


    SUBROUTINE ZEHEGA(N,CS,DZ)
        !************************************************************************
        !*   COMPUTES THE ZEROES OF THE HERMITE POLYNOMIAL OF DEGREE N
        !*   N  = THE NUMBER OF ZEROES
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   DZ = DERIVATIVES OF THE SCALED FUNCTIONS AT THE ZEROES, DZ(I), I=1,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: CS, DZ
        INTEGER :: I, M, IN
        REAL(kind=8) :: CSX, A

        IF (N .EQ. 0) RETURN

        M = N/2
        CS(M+1) = 0.D0
        DZ(M+1) = 1.D0
        IF(N .EQ. 1) RETURN

        IN = 2*N-4*M-1
        IF(IN .EQ. -1) THEN
            A = -.5D0
            CALL ZELAGA(M,A,CS,DZ)
        ELSE
            A = .5D0
            CALL ZELAGA(M,A,CS,DZ)
        ENDIF

        DO I=1,M
            CSX = SQRT(CS(M-I+1))
            CS(N-I+1) = CSX
            IF(IN .EQ. -1) THEN
                DZ(N-I+1) = 2.D0*CSX*DZ(M-I+1)
            ELSE
                DZ(N-I+1) = 2.D0*CSX*CSX*DZ(M-I+1)
            ENDIF
        END DO

        DO I=1,M
            CS(I) = -CS(N-I+1)
            DZ(I) = DZ(N-I+1)*real(IN,kind=8)
        END DO

        RETURN
    END SUBROUTINE ZEHEGA


    SUBROUTINE ZEJAGL(N,A,B,ET,VN)
        !********************************************************************
        !*   COMPUTES THE NODES RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   A  = PARAMETER BETWEEN -1/2 AND 1/2
        !*   B  = PARAMETER BETWEEN -1/2 AND 1/2
        !*   ET = VECTOR OF THE NODES, ET(I), I=0,N
        !*   VN = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: ET, VN
        REAL(kind=8), INTENT(IN) :: A, B
        INTEGER :: I, IT, N1
        REAL(kind=8) :: X, Y, DY, D2Y, AB, DN, A1, B1, ETX, DI, C
        REAL(kind=8), PARAMETER :: EPS=1.D-17

        IF (N .EQ. 0) RETURN

        ET(0) = -1.D0
        ET(N) = 1.D0
        X = -1.D0
        CALL VAJAPO(N,A,B,X,Y,DY,D2Y)
        VN(0) = Y
        X = 1.D0
        CALL VAJAPO(N,A,B,X,Y,DY,D2Y)
        VN(N) = Y

        IF (N .EQ. 1) RETURN

        AB = A+B
        DN = real(N,kind=8)
        C  = M_PI_2/(2.D0*DN+AB+1.D0)
        N1 = N-1
        A1 = A+1.D0
        B1 = B+1.D0
        DO I=1,N1
            DI  = real(I,kind=8)
            ETX = -COS(C*(4.D0*DI+AB+1.D0))
            DO IT=1,8
                CALL VAJAPO(N1,A1,B1,ETX,Y,DY,D2Y)
                IF(ABS(Y) .LE. EPS) EXIT
                ETX = ETX-Y/DY
            END DO
            IF(ABS(ETX) .LE. EPS) ETX=0.D0
            ET(I) = ETX
            VN(I) = -.5D0*DY*(1.D0-ETX*ETX)/DN
        END DO

        RETURN
    END SUBROUTINE ZEJAGL


    SUBROUTINE ZELEGL(N,ET,VN)
        !*********************************************************************
        !*   COMPUTES THE NODES RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   ET = VECTOR OF THE NODES, ET(I), I=0,N
        !*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: ET, VN
        INTEGER :: I, IT, N2
        REAL(kind=8) :: C, SN, X, Y, DY, D2Y, ETX

        IF (N .EQ. 0) RETURN

        N2 = (N-1)/2
        SN = real(2*N-4*N2-3,kind=8)
        ET(0) = -1.D0
        ET(N) = 1.D0
        VN(0) = SN
        VN(N) = 1.D0
        IF (N .EQ. 1) RETURN

        ET(N2+1) = 0.D0
        X = 0.D0
        CALL VALEPO(N,X,Y,DY,D2Y)
        VN(N2+1) = Y
        IF(N .EQ. 2) RETURN

        C  = M_PI/real(N,kind=8)
        DO I=1,N2
            ETX = COS(C*real(I,kind=8))
            DO IT=1,8
                CALL VALEPO(N,ETX,Y,DY,D2Y)
                ETX = ETX-DY/D2Y
            END DO
            ET(I) = -ETX
            ET(N-I) = ETX
            VN(I) = Y*SN
            VN(N-I) = Y
        END DO

        RETURN
    END SUBROUTINE ZELEGL


    SUBROUTINE ZECHGL(N,ET)
        !**********************************************************************
        !*   COMPUTES THE NODES RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   ET = VECTOR OF THE NODES, ET(I), I=0,N
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: ET
        INTEGER :: I, N2
        REAL(kind=8) :: C, ETX

        IF (N .EQ. 0) RETURN

        ET(0) = -1.D0
        ET(N) = 1.D0
        IF (N .EQ. 1) RETURN

        N2 = (N-1)/2
        ET(N2+1) = 0.D0
        IF(N .EQ. 2) RETURN

        C  = M_PI/real(N,kind=8)
        DO I=1,N2
            ETX = COS(C*real(I,kind=8))
            ET(I) = -ETX
            ET(N-I) = ETX
        END DO
        RETURN
    END SUBROUTINE ZECHGL


    SUBROUTINE ZELAGR(N,A,ET,VN)
        !************************************************************************
        !*   COMPUTES THE NODES RELATIVE TO THE LAGUERRE GAUSS-RADAU FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   A  = PARAMETER >-1
        !*   ET = VECTOR OF THE NODES, ET(I), 0=1,N-1
        !*   VN = SCALED LAGUERRE FUNCTION AT THE NODES, VN(I), I=0,N-1
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: ET, VN
        INTEGER :: M, IT, N1, IN
        REAL(kind=8) :: A1, DN, DM, C1, C2, XN, XP, Z, ZD, ETX, Y, DY, ETM, VNM

        IF (N .EQ. 0) RETURN

        ET(0) = 0.D0
        VN(0) = 1.D0
        IF (N .EQ. 1) RETURN

        A1 = A+1.D0
        N1 = N-1
        DN = real(N1,kind=8)
        C1 = 2.D0*DN+A1+1.D0
        DO M=1,N1
            DM = real(M,kind=8)
            C2 = 2.D0*(DN+.75D0-DM)*M_PI/C1
            XN = (C2+M_PI)/2.D0
            DO IT=1,8
                XP = XN
                XN = (SIN(XP)-XP*COS(XP)+C2)/(1.D0-COS(XP))
            END DO
            Z  = (COS(XN/2.D0))**2
            ZD = 1.D0/(Z-1.D0)
            ETX= 2.D0*C1*Z-((1.25D0*ZD+1.D0)*ZD-1.D0+3.D0*A1*A1)/(6.D0*C1)
            DO IT=1,6
                CALL VALASF(N1,A1,ETX,Y,DY)
                ETX = ETX-Y/DY
            END DO
            ET(M) = ETX
            CALL VALASF(N,A,ETX,Y,DY)
            VN(M) = Y
        END DO

        DO
            IN = 0
            DO M=1,N-2
                IF(ET(M) .LE. ET(M+1)) CYCLE
                ETM = ET(M)
                ET(M) = ET(M+1)
                ET(M+1) = ETM
                VNM = VN(M)
                VN(M) = VN(M+1)
                VN(M+1) = VNM
                IN = 1
            END DO
            IF(IN .NE. 1) EXIT
        END DO
        RETURN
    END SUBROUTINE ZELAGR


    SUBROUTINE WEJAGA(N,A,B,CS,DZ,WE)
        !****************************************************************
        !*   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   A  = PARAMETER > -1
        !*   B  = PARAMETER > -1
        !*   CS = ZEROES OF THE JACOBI POLYNOMIAL, CS(I), I=1,N
        !*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
        !****************************************************************
        IMPLICIT NONE
        REAL(kind=8), INTENT(IN) :: A, B
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, DZ
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: WE
        INTEGER :: I, M
        REAL(kind=8) :: AB, A2, B2, GA2, GB2, GAB, C, X, DY, DM

        IF (N .EQ. 0) RETURN

        AB = A+B+2.D0
        A2 = A+2.D0
        B2 = B+2.D0
        CALL GAMMAF(A2,GA2)
        CALL GAMMAF(B2,GB2)
        CALL GAMMAF(AB,GAB)
        C  = .5D0*(2.D0**AB)*GA2*GB2/GAB
        DO M=2,N
            DM = real(M,kind=8)
            C  = C*(DM+A)*(DM+B)/(DM*(DM+A+B))
        END DO

        DO  I=1,N
            X  = CS(I)
            DY = DZ(I)
            WE(I) = C/((1.D0-X*X)*DY*DY)
        END DO

        RETURN
    END SUBROUTINE WEJAGA


    SUBROUTINE WELEGA(N,CS,DZ,WE)
        !*****************************************************************
        !*   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   CS = ZEROES OF THE LEGENDRE POLYNOMIAL, CS(I), I=1,N
        !*   DZ = VECTOR OF THE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
        !*****************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, DZ
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: WE
        INTEGER :: I, N2
        REAL(kind=8) :: WEX, DY, X

        IF (N .EQ. 0) RETURN

        N2 = N/2
        DO  I=1,N2
            X  = CS(I)
            DY = DZ(I)
            WEX = 2.D0/((1.D0-X*X)*DY*DY)
            WE(I) = WEX
            WE(N-I+1) = WEX
        END DO

        IF(N .EQ. 2*N2) RETURN
        DY = DZ(N2+1)
        WE(N2+1) = 2.D0/(DY*DY)

        RETURN
    END SUBROUTINE WELEGA


    SUBROUTINE WECHGA(N,WE)
        !*****************************************************************
        !*   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
        !*****************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: WE
        INTEGER :: I
        REAL(kind=8) :: C

        IF (N .EQ. 0) RETURN

        C  = M_PI/real(N,kind=8)
        DO  I=1,N
            WE(I) = C
        END DO

        RETURN
    END SUBROUTINE WECHGA


    SUBROUTINE WELAGA(N,A,CS,WE)
        !****************************************************************
        !*   COMPUTES THE WEIGHTS RELATIVE TO THE LAGUERRE GAUSS FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   A  = PARAMETER > -1
        !*   CS = ZEROES OF THE LAGUERRE POLYNOMIAL, CS(I), I=1,N
        !*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
        !****************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: WE
        INTEGER :: M, I, N1
        REAL(kind=8) :: A1, DN, C, GA1, X, Y, DY, D2Y, DM

        IF (N .EQ. 0) RETURN

        A1 = A+1.D0
        CALL GAMMAF(A1,GA1)
        N1 = N+1
        DN = real(N1,kind=8)
        C = GA1/DN
        DO M=1,N
            DM = real(M,kind=8)
            C = C*(DM+A)/(DM+1.D0)
        END DO

        DO I=1,N
            X =CS(I)
            CALL VALAPO(N1,A,X,Y,DY,D2Y)
            WE(I) = C*X/(Y*Y)
        END DO

        RETURN
    END SUBROUTINE WELAGA

    SUBROUTINE WEHEGA(N,CS,WE)
        !****************************************************************
        !*   COMPUTES THE WEIGHTS RELATIVE TO THE HERMITE GAUSS FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   CS = ZEROES OF THE HERMITE POLYNOMIAL, CS(I), I=1,N
        !*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
        !****************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: WE
        INTEGER :: I, K, N2
        REAL(kind=8) :: PR, R2, C, X, YP, Y, DK, RK, QK, YM, WEX

        IF (N .EQ. 0) RETURN

        PR = 1.77245385090551588D0
        R2 = 1.41421356237309515D0
        C  = PR/real(N,kind=8)
        N2 = N/2
        DO I=1,N2
            X  = CS(I)
            YP = 1.D0
            Y  = R2*X
            DO  K=2,N-1
                DK = real(K,kind=8)
                RK = SQRT(DK)
                QK = SQRT(DK-1.D0)
                YM = Y
                Y  = (R2*X*Y-QK*YP)/RK
                YP = YM
            END DO
            WEX = C/(Y*Y)
            WE(I) = WEX
            WE(N-I+1) = WEX
        END DO

        IF(N .EQ. 2*N2) RETURN
        Y = 1.D0
        DO K=2,N-1,2
            DK = real(K,kind=8)
            Y  = Y*SQRT((DK-1.D0)/DK)
        END DO
        WE(N2+1) = C/(Y*Y)

        RETURN
    END SUBROUTINE WEHEGA

    SUBROUTINE WEJAGL(N,A,B,ET,WT)
        !*********************************************************************
        !*   COMPUTES THE WEIGHTS RELATIVE TO THE JACOBI GAUSS-LOBATTO FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   A  = PARAMETER > -1
        !*   B  = PARAMETER > -1
        !*   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
        !*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: WT
        INTEGER :: I, K, M, N1
        REAL(kind=8) :: C, C1, C2, C3, C4, A1, B1, AB, AB2, DK, SU
        REAL(kind=8) :: GA1, GB1, GAB2, Y, DY, D2Y, DN, X
        REAL(kind=8), PARAMETER :: EPS=1.D-17

        IF (N .EQ. 0) RETURN

        A1  = A+1.D0
        B1  = B+1.D0
        AB  = A+B
        AB2 = A1+B1
        CALL GAMMAF(A1,GA1)
        CALL GAMMAF(B1,GB1)
        CALL GAMMAF(AB2,GAB2)
        C = (2.D0**AB)*GA1*GB1/GAB2
        WT(0) = 2.D0*C*A1/AB2
        WT(N) = 2.D0*C*B1/AB2
        IF (N .EQ. 1) RETURN

        N1 = N-1
        DN = real(N,kind=8)
        C  = C*(2.D0*DN+AB)/(DN+AB+1.D0)
        C1 = C*A1/((B+2.D0)*AB2)
        C2 = C*B1/((A+2.D0)*AB2)
        C3 = .5D0*C*A1*B1
        DO K=1,N-2
            DK = real(K,kind=8)
            C1 = C1*(DK+A1)*DK/((DK+AB2)*(DK+B+2.D0))
            C2 = C2*(DK+B1)*DK/((DK+AB2)*(DK+A+2.D0))
            C3 = C3*(DK+A1)*(DK+B1)/((DK+2.D0)*(DK+AB+1.D0))
        END DO

        SU = 0.D0
        DO M=1,N1
            SU = SU+ET(M)
        END DO
        WT(0) = C1*(DN-1.D0-SU)
        WT(N) = C2*(DN-1.D0+SU)
        DO I=1,N1
            X = ET(I)
            CALL VAJAPO(N,A,B,X,Y,DY,D2Y)
            C4 = -C3/Y
            CALL VAJAPO(N1,A,B,X,Y,DY,D2Y)
            WT(I) = C4/DY
        END DO

        RETURN
    END SUBROUTINE WEJAGL

    SUBROUTINE WELEGL(N,ET,VN,WT)
        !***********************************************************************
        !*   COMPUTES THE WEIGHTS RELATIVE TO THE LEGENDRE GAUSS-LOBATTO FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   ET = JACOBI GAUSS-LOBATTO NODES, ET(I), I=0,N
        !*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
        !***********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, VN
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: WT
        INTEGER :: I, N2
        REAL(kind=8) :: DN, C, X, Y, WTX

        IF (N .EQ. 0) RETURN

        N2 = (N-1)/2
        DN = real(N,kind=8)
        C  = 2.D0/(DN*(DN+1.D0))
        DO I=0,N2
            X = ET(I)
            Y = VN(I)
            WTX = C/(Y*Y)
            WT(I) = WTX
            WT(N-I) = WTX
        END DO

        IF(N-1 .EQ. 2*N2) RETURN
        X = 0.D0
        Y = VN(N2+1)
        WT(N2+1) = C/(Y*Y)

        RETURN
    END SUBROUTINE WELEGL

    SUBROUTINE WECHGL(N,WT)
        !************************************************************************
        !*   COMPUTES THE WEIGHTS RELATIVE TO THE CHEBYSHEV GAUSS-LOBATTO FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: WT
        INTEGER :: I
        REAL(kind=8) :: C, C2

        IF (N .EQ. 0) RETURN

        C  = M_PI/real(N,kind=8)
        C2 = .5D0*C
        WT(0) =C2
        WT(N) =C2
        IF (N .EQ. 1) RETURN

        DO  I=1,N-1
            WT(I) = C
        END DO

        RETURN
    END SUBROUTINE WECHGL

    SUBROUTINE WELAGR(N,A,ET,WT)
        !*********************************************************************
        !*   COMPUTES THE WEIGHTS RELATIVE TO THE LAGUERRE GAUSS-RADAU FORMULA
        !*   N  = ORDER OF THE FORMULA
        !*   A  = PARAMETER > -1
        !*   ET = LAGUERRE GAUSS-RADAU NODES, ET(I), I=0,N-1
        !*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N-1
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(0:N-1), INTENT(IN) :: ET
        REAL(kind=8), DIMENSION(0:N-1), INTENT(OUT) :: WT
        INTEGER :: I, K, N1
        REAL(kind=8) :: A1, C1, GA1, C2, DK, X, C3, Y, DY, D2Y

        IF (N .EQ. 0) RETURN

        A1 = A+1.D0
        CALL GAMMAF(A1,GA1)
        C1 = GA1
        WT(0) = C1
        IF (N .EQ. 1) RETURN

        N1 = N-1
        C2 = GA1
        DO K=1,N1
            DK = real(K,kind=8)
            C1 = C1*DK/(DK+A1)
            C2 = C2*(DK+A)/(DK+1.D0)
        END DO
        WT(0) = C1
        DO I=1,N1
            X = ET(I)
            CALL VALAPO(N,A,X,Y,DY,D2Y)
            C3 = C2/Y
            CALL VALAPO(N1,A,X,Y,DY,D2Y)
            WT(I) = C3/DY
        END DO

        RETURN
    END SUBROUTINE WELAGR

    SUBROUTINE WECHCC(N,WK)
        !*********************************************************************
        !*   COMPUTES THE WEIGHTS OF THE CLENSHAW-CURTIS FORMULA OF ORDER 2*N
        !*   N  = INTEGER PARAMETER
        !*   WK = VECTOR OF THE WEIGHTS, WK(I), I=0,2*N
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:2*N), INTENT(OUT) :: WK
        INTEGER :: M, J, K
        REAL(kind=8) :: DN, DM, C, DJ, SU, DK

        IF (N .EQ. 0) RETURN

        M  = 2*N
        DN = real(N,kind=8)
        DM = real(M,kind=8)
        C  = 1.D0/(DM*DM-1.D0)
        WK(0) = C
        WK(M) = C
        WK(N) = 1.33333333333333333333D0
        IF (N .EQ. 1) RETURN

        DO J=1,N-1
            DJ = real(J,kind=8)
            SU = 1.D0-((-1.D0)**J)*C
            DO K=1,N-1
                DK = 2.D0*real(K,kind=8)
                SU = SU+2.D0*COS(DJ*DK*M_PI/DM)/(1.D0-DK*DK)
            END DO
            WK(J) = SU/DN
            WK(M-J) = SU/DN
        END DO

        SU = 1.D0-((-1.D0)**N)*C
        DO  K=1,N-1
            DK = 2.D0*real(K,kind=8)
            SU = SU+2.D0*((-1.D0)**K)/(1.D0-DK*DK)
        END DO
        WK(N) = SU/DN

        RETURN
    END SUBROUTINE WECHCC

    SUBROUTINE INJAGA(N,A,B,CS,DZ,QZ,X,QX)
        !************************************************************************
        !*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
        !*   BY THE VALUES ATTAINED AT THE ZEROES OF THE JACOBI POLYNOMIAL
        !*   N  = THE NUMBER OF ZEROES
        !*   A  = PARAMETER > -1
        !*   B  = PARAMETER > -1
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   DZ = JACOBI POLYNOMIAL DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   QX = VALUE OF THE POLYNOMIAL IN X
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B, X
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, DZ, QZ
        REAL(kind=8), INTENT(OUT) :: QX
        INTEGER :: J
        REAL(kind=8) :: ED, Y, DY, D2Y
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        CALL VAJAPO(N,A,B,X,Y,DY,D2Y)
        QX = 0.D0
        DO  J=1,N
            ED = X-CS(J)
            IF(ABS(ED) .LT. EPS) THEN
                QX = QZ(J)
                RETURN
            ELSE
                QX = QX+QZ(J)*Y/(DZ(J)*ED)
            ENDIF
        END DO

        RETURN
    END SUBROUTINE INJAGA

    SUBROUTINE INLEGA(N,CS,DZ,QZ,X,QX)
        !************************************************************************
        !*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
        !*   BY THE VALUES ATTAINED AT THE ZEROES OF THE LEGENDRE POLYNOMIAL
        !*   N  = THE NUMBER OF ZEROES
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   DZ = LEGENDRE POLYNOMIAL DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   QX = VALUE OF THE POLYNOMIAL IN X
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, DZ, QZ
        REAL(kind=8), INTENT(OUT) :: QX
        INTEGER :: J
        REAL(kind=8) :: ED, Y, DY, D2Y
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        CALL VALEPO(N,X,Y,DY,D2Y)
        QX = 0.D0
        DO J=1,N
            ED = X-CS(J)
            IF(ABS(ED) .LT. EPS) THEN
                QX = QZ(J)
                RETURN
            ELSE
                QX = QX+QZ(J)*Y/(DZ(J)*ED)
            ENDIF
        END DO

        RETURN
    END SUBROUTINE INLEGA

    SUBROUTINE INCHGA(N,DZ,QZ,X,QX)
        !************************************************************************
        !*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
        !*   BY THE VALUES ATTAINED AT THE ZEROES OF THE CHEBYSHEV POLYNOMIAL
        !*   N  = THE NUMBER OF ZEROES
        !*   DZ = CHEBYSHEV POLYNOMIAL DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   QX = VALUE OF THE POLYNOMIAL IN X
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: DZ, QZ
        REAL(kind=8), INTENT(OUT) :: QX
        INTEGER :: J
        REAL(kind=8) :: DN, C, ED, Y, DY, D2Y
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        C  = M_PI_2/DN
        CALL VACHPO(N,X,Y,DY,D2Y)
        QX = 0.D0
        DO J=1,N
            ED = X+COS(C*(2.D0*real(J,kind=8)-1.D0))
            IF(ABS(ED) .LT. EPS) THEN
                QX = QZ(J)
                RETURN
            ELSE
                QX = QX+QZ(J)*Y/(DZ(J)*ED)
            ENDIF
        END DO

        RETURN
    END SUBROUTINE INCHGA

    SUBROUTINE INLAGA(N,A,CS,DZ,QZ,X,QX)
        !*********************************************************************
        !*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
        !*   BY THE VALUES ATTAINED AT THE ZEROES OF THE LAGUERRE POLYNOMIAL
        !*   N  = THE NUMBER OF ZEROES
        !*   A  = PARAMETER > -1
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   DZ = SCALED LAGUERRE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   QX = VALUE OF THE POLYNOMIAL IN X
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, X
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, DZ, QZ
        REAL(kind=8), INTENT(OUT) :: QX
        INTEGER :: J, K
        REAL(kind=8) :: ED, Y, DY, PR, DK
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        CALL VALASF(N,A,X,Y,DY)
        QX = 0.D0
        DO J=1,N
            ED = X-CS(J)
            IF(ABS(ED) .LT. EPS) THEN
                QX = QZ(J)
                RETURN
            ELSE
                PR = 1.D0
                DO  K=1,N
                    DK = 4.D0*real(K,kind=8)
                    PR = PR*(DK+X)/(DK+CS(J))
                END DO
                QX = QX+QZ(J)*Y*PR/(DZ(J)*ED)
            ENDIF
        END DO

        RETURN
    END SUBROUTINE INLAGA


    SUBROUTINE INHEGA(N,CS,DZ,QZ,X,QX)
        !*********************************************************************
        !*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
        !*   BY THE VALUES ATTAINED AT THE ZEROES OF THE HERMITE POLYNOMIAL
        !*   N  = THE NUMBER OF ZEROES
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   DZ = SCALED HERMITE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   QX = VALUE OF THE POLYNOMIAL IN X
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, DZ, QZ
        REAL(kind=8), INTENT(OUT) :: QX
        INTEGER :: J, K
        REAL(kind=8) :: ED, Y, DY, PR, DK
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        QX = QZ(1)
        IF (N .EQ. 1) RETURN

        CALL VAHESF(N,X,Y,DY)
        QX = 0.D0
        DO J=1,N
            ED = X-CS(J)
            IF(ABS(ED) .LT. EPS) THEN
                QX = QZ(J)
                RETURN
            ELSE
                PR = 1.D0
                DO K=1,N/2
                    DK = 4.D0*real(K,kind=8)
                    PR = PR*(DK+X*X)/(DK+CS(J)**2)
                END DO
                QX = QX+QZ(J)*Y*PR/(DZ(J)*ED)
            END IF
        END DO

        RETURN
    END SUBROUTINE INHEGA

    SUBROUTINE INJAGL(N,A,B,ET,VN,QN,X,QX)
        !*********************************************************************
        !*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
        !*   BY THE VALUES ATTAINED AT THE JACOBI GAUSS-LOBATTO NODES
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   A  = PARAMETER > -1
        !*   B  = PARAMETER > -1
        !*   ET = VECTOR OF THE NODES, ET(I), I=0,N
        !*   VN = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   QX = VALUE OF THE POLYNOMIAL IN X
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B, X
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, VN, QN
        REAL(kind=8), INTENT(OUT) :: QX
        INTEGER :: J
        REAL(kind=8) :: ED, Y, DY, D2Y, DN, C
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        CALL VAJAPO(N,A,B,X,Y,DY,D2Y)
        DN = real(N,kind=8)
        C  = 1.D0/(DN*(DN+A+B+1.D0))
        QX = QN(0)*C*DY*(B+1.D0)*(X-1.D0)/VN(0)
        QX = QX+QN(N)*C*DY*(A+1.D0)*(X+1.D0)/VN(N)
        IF (N .EQ. 1) RETURN

        DO J=1,N-1
            ED = X-ET(J)
            IF(ABS(ED) .LT. EPS) THEN
                QX = QN(J)
                RETURN
            ELSE
                QX = QX+QN(J)*C*DY*(X*X-1.D0)/(VN(J)*ED)
            ENDIF
        END DO

        RETURN
    END SUBROUTINE INJAGL

    SUBROUTINE INLEGL(N,ET,VN,QN,X,QX)
        !**********************************************************************
        !*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
        !*   BY THE VALUES ATTAINED AT THE LEGENDRE GAUSS-LOBATTO NODES
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   ET = VECTOR OF THE NODES, ET(I), I=0,N
        !*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   QX = VALUE OF THE POLYNOMIAL IN X
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, VN, QN
        REAL(kind=8), INTENT(OUT) :: QX
        INTEGER :: J
        REAL(kind=8) :: ED, Y, DY, D2Y, DN, C
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        CALL VALEPO(N,X,Y,DY,D2Y)
        DN = real(N,kind=8)
        C  = 1.D0/(DN*(DN+1.D0))
        QX = QN(0)*C*DY*(X-1.D0)/VN(0)
        QX = QX+QN(N)*C*DY*(X+1.D0)/VN(N)
        IF (N .EQ. 1) RETURN

        DO J=1,N-1
            ED = X-ET(J)
            IF(ABS(ED) .LT. EPS) THEN
                QX = QN(J)
                RETURN
            ELSE
                QX = QX+QN(J)*C*DY*(X*X-1.D0)/(VN(J)*ED)
            ENDIF
        END DO

        RETURN
    END SUBROUTINE INLEGL

    SUBROUTINE INCHGL(N,QN,X,QX)
        !*********************************************************************
        !*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
        !*   BY THE VALUES ATTAINED AT THE CHEBYSHEV GAUSS-LOBATTO NODES
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   QX = VALUE OF THE POLYNOMIAL IN X
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: QN
        REAL(kind=8), INTENT(OUT) :: QX
        INTEGER :: J
        REAL(kind=8) :: ED, Y, DY, D2Y, DN, C, SN, PN, SJ
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        CALL VACHPO(N,X,Y,DY,D2Y)
        DN = real(N,kind=8)
        SN = real(1+4*(N/2)-2*N,kind=8)
        PN = M_PI/DN
        C  = 1.D0/(DN*DN)
        QX = .5D0*SN*QN(0)*C*DY*(X-1.D0)
        QX = QX+.5D0*QN(N)*C*DY*(X+1.D0)
        IF (N .EQ. 1) RETURN

        SJ = -1.D0
        DO J=1,N-1
            ED = X+COS(PN*real(J,kind=8))
            IF(ABS(ED) .LT. EPS) THEN
                QX = QN(J)
                RETURN
            ELSE
                QX = QX+QN(J)*SN*SJ*C*DY*(X*X-1.D0)/ED
            ENDIF
            SJ = -SJ
        END DO

        RETURN
    END SUBROUTINE INCHGL

    SUBROUTINE INLAGR(N,A,ET,VN,QN,X,QX)
        !*********************************************************************
        !*   COMPUTES THE VALUE AT A GIVEN POINT OF A POLYNOMIAL INDIVIDUATED
        !*   BY THE VALUES ATTAINED AT THE LAGUERRE GAUSS-RADAU NODES
        !*   N  = THE NUMBER OF NODES
        !*   A  = PARAMETER > -1
        !*   ET = VECTOR OF THE NODES, ET(I), I=0,N-1
        !*   VN = SCALED LAGUERRE FUNCTION AT THE NODES, VN(I), I=0,N-1
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N-1
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   QX = VALUE OF THE POLYNOMIAL IN X
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, X
        REAL(kind=8), DIMENSION(0:N-1), INTENT(IN) :: ET, VN, QN
        REAL(kind=8), INTENT(OUT) :: QX
        INTEGER :: J, K, M
        REAL(kind=8) :: ED, Y, DY, DN, C, PR, SU, DM, DK
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        QX = QN(0)
        IF (N .EQ. 1) RETURN

        DN = real(N,kind=8)
        C  = -1.D0/DN
        PR = 1.D0
        SU = 0.D0
        DO M=1,N
            DM = 4.D0*real(M,kind=8)
            PR = PR*(DM+X)/DM
            SU = SU+1.D0/(DM+X)
        END DO
        CALL VALASF(N,A,X,Y,DY)
        QX = QN(0)*C*(A+1.D0)*PR*(DY+SU*Y)/VN(0)
        DO J=1,N-1
            ED = X-ET(J)
            IF(ABS(ED) .LT. EPS) THEN
                QX = QN(J)
                RETURN
            ELSE
                PR = 1.D0
                DO K=1,N
                    DK = 4.D0*real(K,kind=8)
                    PR = PR*(DK+X)/(DK+ET(J))
                END DO
                QX = QX+QN(J)*C*PR*(DY+SU*Y)*X/(VN(J)*ED)
            ENDIF
        END DO

        RETURN
    END SUBROUTINE INLAGR

    SUBROUTINE NOLEGA(N,QZ,WE,QI,QM)
        !*********************************************************************
        !*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE LEGENDRE ZEROES
        !*   N  = THE NUMBER OF ZEROES
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   WE = VECTOR OF THE LEGENDRE GAUSS WEIGHTS, WE(I), I=1,N
        !*   QI = INTEGRAL NORM OF THE POLYNOMIAL
        !*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE ZEROES
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: QZ, WE
        REAL(kind=8), INTENT(OUT) :: QI, QM
        INTEGER :: J
        REAL(kind=8) :: SU, Y
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        SU  = 0.D0
        QM  = 0.D0
        DO J=1,N
            Y  = ABS(QZ(J))
            IF(Y .GT. QM) QM=Y
            IF(Y .LT. EPS) CYCLE
            SU = SU+Y*Y*WE(J)
        END DO
        QI = SQRT(SU)

        RETURN
    END SUBROUTINE NOLEGA

    SUBROUTINE NOCHGA(N,DZ,QZ,WK,QW,QI,QM)
        !**********************************************************************
        !*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE CHEBYSHEV ZEROES
        !*   N  = THE NUMBER OF ZEROES
        !*   DZ = CHEBYSHEV POLYNOMIAL DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   WK = VECTOR OF THE CLENSHAW-CURTIS WEIGHTS, WE(I), I=0,2*N
        !*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL
        !*   QI = INTEGRAL NORM OF THE POLYNOMIAL
        !*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE ZEROES
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: DZ, QZ
        REAL(kind=8), DIMENSION(0:2*N), INTENT(IN) :: WK
        REAL(kind=8), INTENT(OUT) :: QW, QI, QM
        INTEGER :: J, J2
        REAL(kind=8) :: DN, X, S1, S2, DJ, QX, Y
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        X  = -1.D0
        CALL INCHGA(N,DZ,QZ,X,QX)
        S1 = 0.D0
        S2 = QX*QX*WK(0)
        QM = 0.D0
        DO J=1,N
            DJ = real(J,kind=8)
            J2 = 2*J
            X  = -COS(M_PI*DJ/DN)
            Y  = ABS(QZ(J))
            IF(Y .GT. QM) QM=Y
            CALL INCHGA(N,DZ,QZ,X,QX)
            S1 = S1+Y*Y
            S2 = S2+Y*Y*WK(J2-1)+QX*QX*WK(J2)
        END DO
        QW = SQRT(S1*M_PI/DN)
        QI = SQRT(S2)

        RETURN
    END SUBROUTINE NOCHGA

    SUBROUTINE NOJAGL(N,A,B,VN,QN,WT,QW,QS,QM)
        !**********************************************************************
        !*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE JACOBI GAUSS-
        !*   LOBATTO NODES
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   A  = PARAMETER > -1
        !*   B  = PARAMETER > -1
        !*   VN = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   WT = VECTOR OF THE JACOBI GAUSS-LOBATTO WEIGHTS, WT(I), I=0,N
        !*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL
        !*   QS = QUADRATURE NORM OF THE POLYNOMIAL
        !*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE NODES
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: VN, QN, WT
        REAL(kind=8), INTENT(OUT) :: QW, QS, QM
        INTEGER :: K, J
        REAL(kind=8) :: A1, B1, AB, AB2, DN, C, DK, S1, S2, S3
        REAL(kind=8) :: Y1, YM, Y2, GA1, GB1, GAB2
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        A1  = A+1.D0
        B1  = B+1.D0
        AB  = A+B
        AB2 = AB+2.D0
        DN = real(N,kind=8)
        C  = ((2.D0)**(AB+1.D0))*(DN+AB+1.D0)/(2.D0*DN+AB+1.D0)
        CALL GAMMAF(A1,GA1)
        CALL GAMMAF(B1,GB1)
        CALL GAMMAF(AB2,GAB2)
        C  = C*GA1*GB1/GAB2
        DO K=1,N
            DK = real(K,kind=8)
            C  = C*(DK+A)*(DK+B)/(DK*(DK+AB+1.D0))
        END DO

        S1 = 0.D0
        S2 = 0.D0
        S3 = 0.D0
        QM = 0.D0
        DO J=0,N
            Y1 = QN(J)
            YM = ABS(Y1)
            IF(YM .GT. QM) QM=YM
            Y2 = VN(J)
            S2 = S2+Y1*Y2*WT(J)
            S3 = S3+Y2*Y2*WT(J)
            IF(YM .LT. EPS) CYCLE
            S1 = S1+Y1*Y1*WT(J)
        END DO
        QS = SQRT(S1)
        QW = SQRT(S1-(S3-C)*S2*S2/(S3*S3))

        RETURN
    END SUBROUTINE NOJAGL

    SUBROUTINE NOLEGL(N,VN,QN,WT,QI,QS,QM)
        !**********************************************************************
        !*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE LEGENDRE GAUSS-
        !*   LOBATTO NODES
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   VN = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   WT = VECTOR OF THE LEGENDRE GAUSS-LOBATTO WEIGHTS, WT(I), I=0,N
        !*   QW = INTEGRAL NORM OF THE POLYNOMIAL
        !*   QS = QUADRATURE NORM OF THE POLYNOMIAL
        !*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE NODES
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: VN, QN, WT
        REAL(kind=8), INTENT(OUT) :: QI, QS, QM
        INTEGER :: J
        REAL(kind=8) :: DN, C, S1, S2, Y1, YM, Y2
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        C  = .5D0*DN*(DN+1.D0)/(2.D0*DN+1.D0)

        S1 = 0.D0
        S2 = 0.D0
        QM = 0.D0
        DO J=0,N
            Y1 = QN(J)
            YM = ABS(Y1)
            IF(YM .GT. QM) QM=YM
            Y2 = VN(J)
            S2 = S2+Y1*Y2*WT(J)
            IF(YM .LT. EPS) CYCLE
            S1 = S1+Y1*Y1*WT(J)
        END DO
        QS = SQRT(S1)
        QI = SQRT(ABS(S1-C*S2*S2))

        RETURN
    END SUBROUTINE NOLEGL

    SUBROUTINE NOCHGL(N,QN,WK,QW,QI,QS,QM)
        !**********************************************************************
        !*   COMPUTES THE NORMS OF A POLYNOMIAL DEFINED AT THE CHEBYSHEV GAUSS-
        !*   LOBATTO NODES
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   WK = VECTOR OF THE CLENSHAW-CURTIS WEIGHTS, WE(I), I=0,2*N
        !*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL
        !*   QI = INTEGRAL NORM OF THE POLYNOMIAL
        !*   QS = QUADRATURE NORM OF THE POLYNOMIAL
        !*   QM = MAXIMUM VALUE OF THE POLYNOMIAL AT THE NODES
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: QN
        REAL(kind=8), DIMENSION(0:2*N), INTENT(IN) :: WK
        REAL(kind=8), INTENT(OUT) :: QW, QI, QS, QM
        INTEGER :: J, N2, J2
        REAL(kind=8) :: DN, SN, Y, S1, S2, S3, SJ, DJ, X, YM, QX, DD

        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        SN = real(1+4*(N/2)-2*N,kind=8)
        Y  = QN(0)
        S1 = .5D0*Y*Y
        S2 = .5D0*Y*SN
        S3 = Y*Y*WK(0)
        QM = ABS(Y)

        SJ = -1.D0
        DO J=1,N-1
            J2 = 2*J
            DJ = real(J2-1,kind=8)
            X  = -COS(M_PI_2*DJ/DN)
            Y  = QN(J)
            YM = ABS(Y)
            IF(YM .GT. QM) QM=YM
            CALL INCHGL(N,QN,X,QX)
            S1 = S1+Y*Y
            S2 = S2+Y*SN*SJ
            S3 = S3+QX*QX*WK(J2-1)+Y*Y*WK(J2)
            SJ = -SJ
        END DO
        N2 = 2*N
        DD = real(N2-1,kind=8)
        X  = -COS(M_PI_2*DD/DN)
        Y  = QN(N)
        YM = ABS(Y)
        IF(YM .GT. QM) QM=YM
        CALL INCHGL(N,QN,X,QX)
        S1 = S1+.5D0*Y*Y
        S2 = S2+.5D0*Y
        S3 = S3+QX*QX*WK(N2-1)+Y*Y*WK(N2)

        QW = SQRT(ABS(M_PI*S1/DN-M_PI_2*S2*S2/(DN*DN)))
        QI = SQRT(S3)
        QS = SQRT(M_PI*S1/DN)

        RETURN
    END SUBROUTINE NOCHGL

    SUBROUTINE COJAGA(N,A,B,CS,QZ,WE,CO)
        !**********************************************************************
        !*   COMPUTES THE JACOBI FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*   INDIVIDUATED BY THE VALUES ATTAINED AT THE JACOBI ZEROES
        !*   N  = THE NUMBER OF ZEROES
        !*   A  = PARAMETER >-1
        !*   B  = PARAMETER >-1
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, QZ, WE
        REAL(kind=8), DIMENSION(0:N-1), INTENT(OUT) :: CO
        INTEGER :: J, K
        REAL(kind=8) :: A1, B1, AB, AB2, GA1, GB1, GAB2, C, SU
        REAL(kind=8) :: X, YP, Y, DK, CC, C1, C2, C3, C4, YM

        IF (N .EQ. 0) RETURN

        A1  = A+1.D0
        B1  = B+1.D0
        AB  = A+B
        AB2 = AB+2.D0
        CALL GAMMAF(A1,GA1)
        CALL GAMMAF(B1,GB1)
        CALL GAMMAF(AB2,GAB2)
        C  = ((2.D0)**(AB+1.D0))*GA1*GB1/GAB2

        SU = 0.D0
        DO J=1,N
            SU = SU+QZ(J)*WE(J)
            CO(J-1) = 0.D0
        END DO
        CO(0) = SU/C
        IF (N .EQ. 1) RETURN

        DO J=1,N
            X  = CS(J)
            YP = QZ(J)*WE(J)
            Y  = .5D0*YP*(AB2*X+A-B)
            DO K=1,N-1
                CO(K) = CO(K)+Y
                DK = real(K+1,kind=8)
                CC = 2.D0*DK+AB
                C1 = 2.D0*DK*(DK+AB)*(CC-2.D0)
                C2 = (CC-1.D0)*(CC-2.D0)*CC
                C3 = (CC-1.D0)*(A-B)*AB
                C4 = 2.D0*(DK+A-1.D0)*CC*(DK+B-1.D0)
                YM = Y
                Y  = ((C2*X+C3)*Y-C4*YP)/C1
                YP = YM
            END DO
        END DO

        DO K=1,N-1
            DK = real(K,kind=8)
            C  = C*(DK+A)*(DK+B)/DK
            CO(K) = CO(K)*(2.D0*DK+AB+1.D0)/C
            C  = C/(DK+AB+1.D0)
        END DO

        RETURN
    END SUBROUTINE COJAGA

    SUBROUTINE COLEGA(N,CS,QZ,WE,CO)
        !**********************************************************************
        !*   COMPUTES THE LEGENDRE FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*   INDIVIDUATED BY THE VALUES ATTAINED AT THE LEGENDRE ZEROES
        !*   N  = THE NUMBER OF ZEROES
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, QZ, WE
        REAL(kind=8), DIMENSION(0:N-1), INTENT(OUT) :: CO
        INTEGER :: J, K
        REAL(kind=8) :: SU, X, YP, Y, DK, C1, C2, YM

        IF (N .EQ. 0) RETURN

        SU = 0.D0
        DO J=1,N
            SU = SU+QZ(J)*WE(J)
            CO(J-1) = 0.D0
        END DO
        CO(0) = .5D0*SU
        IF (N .EQ. 1) RETURN

        DO J=1,N
            X  = CS(J)
            YP = QZ(J)*WE(J)
            Y  = X*YP
            DO K=1,N-1
                CO(K) = CO(K)+Y
                DK = real(K+1,kind=8)
                C1 = 2.D0*DK-1.D0
                C2 = DK-1.D0
                YM = Y
                Y  = (C1*X*Y-C2*YP)/DK
                YP = YM
            END DO
        END DO

        DO K=1,N-1
            CO(K) = .5D0*CO(K)*(2.D0*real(K,kind=8)+1.D0)
        END DO

        RETURN
    END SUBROUTINE COLEGA

    SUBROUTINE COCHGA(N,QZ,CO)
        !**********************************************************************
        !*   COMPUTES THE CHEBYSHEV FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*   INDIVIDUATED BY THE VALUES ATTAINED AT THE CHEBYSHEV ZEROES
        !*   N  = THE NUMBER OF ZEROES
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: QZ
        REAL(kind=8), DIMENSION(0:N-1), INTENT(OUT) :: CO
        INTEGER :: J, K
        REAL(kind=8) :: DN, SU, SK, DK, DJ

        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        SU = 0.D0
        DO J=1,N
            SU = SU+QZ(J)
        END DO
        CO(0) = SU/DN
        IF (N .EQ. 1) RETURN

        SK = -2.D0
        DO K=1,N-1
            DK = real(K,kind=8)
            SU = 0.D0
            DO J=1,N
                DJ = 2.D0*real(J,kind=8)-1.D0
                SU = SU+QZ(J)*COS(DK*DJ*M_PI_2/DN)
            END DO
            CO(K) = SK*SU/DN
            SK = -SK
        END DO

        RETURN
    END SUBROUTINE COCHGA

    SUBROUTINE COLAGA(N,A,CS,QZ,WE,CO)
        !**********************************************************************
        !*   COMPUTES THE LAGUERRE FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*   INDIVIDUATED BY THE VALUES ATTAINED AT THE LAGUERRE ZEROES
        !*   N  = THE NUMBER OF ZEROES
        !*   A  = PARAMETER >-1
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, QZ, WE
        REAL(kind=8), DIMENSION(0:N-1), INTENT(OUT) :: CO
        INTEGER :: J, K
        REAL(kind=8) :: A1, SU, X, YP, Y, DK, C, B1, B2, YM

        IF (N .EQ. 0) RETURN

        A1  = A+1.D0
        CALL GAMMAF(A1,C)
        SU = 0.D0
        DO J=1,N
            SU = SU+QZ(J)*WE(J)
            CO(J-1) = 0.D0
        END DO
        CO(0) = SU/C
        IF (N .EQ. 1) RETURN

        DO J=1,N
            X  = CS(J)
            YP = QZ(J)*WE(J)
            Y  = (A1-X)*YP
            DO K=1,N-1
                CO(K) = CO(K)+Y
                DK = real(K+1,kind=8)
                B1 = (2.D0*DK+A-1.D0-X)/DK
                B2 = (DK+A-1.D0)/DK
                YM = Y
                Y  = B1*Y-B2*YP
                YP = YM
            END DO
        END DO

        DO K=1,N-1
            DK = real(K,kind=8)
            C  = C*(DK+A)/DK
            CO(K) = CO(K)/C
        END DO

        RETURN
    END SUBROUTINE COLAGA

    SUBROUTINE COHEGA(N,CS,QZ,WE,CO)
        !**********************************************************************
        !*   COMPUTES THE HERMITE FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*   INDIVIDUATED BY THE VALUES ATTAINED AT THE HERMITE ZEROES
        !*   N  = THE NUMBER OF ZEROES
        !*   CS = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*   QZ = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   WE = VECTOR OF THE WEIGHTS, WE(I), I=1,N
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, QZ, WE
        REAL(kind=8), DIMENSION(0:N-1), INTENT(OUT) :: CO
        INTEGER :: J, K
        REAL(kind=8) :: PR, SU, X, YP, Y, DK, YM

        IF (N .EQ. 0) RETURN

        PR = 1.77245385090551588D0
        SU = 0.D0
        DO J=1,N
            SU = SU+QZ(J)*WE(J)
            CO(J-1) = 0.D0
        END DO
        CO(0) = SU/PR
        IF (N .EQ. 1) RETURN

        DO J=1,N
            X  = CS(J)
            YP = QZ(J)*WE(J)/PR
            Y  = X*YP
            DO K=1,N-1
                CO(K) = CO(K)+Y
                DK = real(K+1,kind=8)
                YM = Y
                Y  = (X*Y-.5D0*YP)/DK
                YP = YM
            END DO
        END DO
        RETURN
    END SUBROUTINE COHEGA

    SUBROUTINE COJAGL(N,A,B,ET,VN,QN,WT,CO)
        !**********************************************************************
        !*   COMPUTES THE JACOBI FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*   INDIVIDUATED BY ITS VALUES AT THE JACOBI GAUSS-LOBATTO NODES
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   A  = PARAMETER >-1
        !*   B  = PARAMETER >-1
        !*   ET = VECTOR OF THE NODES, ET(I), I=0,N
        !*   VN = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, VN, QN, WT
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CO
        INTEGER :: J, K
        REAL(kind=8) :: A1, B1, AB, AB2, GA1, GB1, GAB2, C, SU, X, YP, Y
        REAL(kind=8) :: CN, DK, CC, C1, C2, C3, C4, YM

        CO(0) = QN(0)
        IF (N .EQ. 0) RETURN

        A1  = A+1.D0
        B1  = B+1.D0
        AB  = A+B
        AB2 = AB+2.D0
        CO(1) = (QN(1)-QN(0))/AB2
        CO(0) = .5D0*(QN(0)+QN(1)-(A-B)*CO(1))
        IF (N .EQ. 1) RETURN

        CALL GAMMAF(A1,GA1)
        CALL GAMMAF(B1,GB1)
        CALL GAMMAF(AB2,GAB2)
        C  = ((2.D0)**(AB+1.D0))*GA1*GB1/GAB2
        SU = 0.D0
        DO J=0,N
            SU = SU+QN(J)*WT(J)
            CO(J) = 0.D0
        END DO
        CO(0) = SU/C

        CN = 0.D0
        DO J=0,N
            X  = ET(J)
            YP = QN(J)*WT(J)
            Y  = .5D0*YP*(AB2*X+A-B)
            CN = CN+VN(J)*VN(J)*WT(J)
            DO K=1,N
                CO(K) = CO(K)+Y
                DK = real(K+1,kind=8)
                CC = 2.D0*DK+AB
                C1 = 2.D0*DK*(DK+AB)*(CC-2.D0)
                C2 = (CC-1.D0)*(CC-2.D0)*CC
                C3 = (CC-1.D0)*(A-B)*AB
                C4 = 2.D0*(DK+A-1.D0)*CC*(DK+B-1.D0)
                YM = Y
                Y  = ((C2*X+C3)*Y-C4*YP)/C1
                YP = YM
            END DO
        END DO
        DO K=1,N-1
            DK = real(K,kind=8)
            C  = C*(DK+A)*(DK+B)/DK
            CO(K) = CO(K)*(2.D0*DK+AB+1.D0)/C
            C  = C/(DK+AB+1.D0)
        END DO
        CO(N) = CO(N)/CN

        RETURN
    END SUBROUTINE COJAGL

    SUBROUTINE COLEGL(N,ET,QN,WT,CO)
        !**********************************************************************
        !*   COMPUTES THE LEGENDRE FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*   INDIVIDUATED BY ITS VALUES AT THE LEGENDRE GAUSS-LOBATTO NODES
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   ET = VECTOR OF THE NODES, ET(I), I=0,N
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, QN, WT
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CO
        INTEGER :: J, K
        REAL(kind=8) :: SU, X, YP, Y, DK, C1, C2, YM, DN

        CO(0) = QN(0)
        IF (N .EQ. 0) RETURN

        CO(0) = .5D0*(QN(0)+QN(1))
        CO(1) = .5D0*(QN(1)-QN(0))
        IF (N .EQ. 1) RETURN

        SU = 0.D0
        DO J=0,N
            SU = SU+QN(J)*WT(J)
            CO(J) = 0.D0
        END DO
        CO(0) = .5D0*SU

        DO J=0,N
            X  = ET(J)
            YP = QN(J)*WT(J)
            Y  = X*YP
            DO K=1,N
                CO(K) = CO(K)+Y
                DK = real(K+1,kind=8)
                C1 = 2.D0*DK-1.D0
                C2 = DK-1.D0
                YM = Y
                Y  = (C1*X*Y-C2*YP)/DK
                YP = YM
            END DO
        END DO

        DN = real(N,kind=8)
        CO(N) = .5D0*DN*CO(N)
        IF (N .EQ. 1) RETURN

        DO K=1,N-1
            CO(K) = .5D0*CO(K)*(2.D0*real(K,kind=8)+1.D0)
        END DO

        RETURN
    END SUBROUTINE COLEGL

    SUBROUTINE COCHGL(N,QN,CO)
        !**********************************************************************
        !*   COMPUTES THE CHEBYSHEV FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*   INDIVIDUATED BY ITS VALUES AT THE CHEBYSHEV GAUSS-LOBATTO NODES
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: QN
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CO
        INTEGER :: J, K
        REAL(kind=8) :: DN, DD, S0, SN, SJ, SK, DK, SU, DJ

        CO(0) = QN(0)
        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        DD = real(1+4*(N/2)-2*N,kind=8)
        CO(0) = .5D0*(QN(0)+QN(N))
        CO(N) = .5D0*(QN(0)+DD*QN(N))
        IF (N .EQ. 1) RETURN

        S0 = CO(0)
        SN = CO(N)
        SJ = -1.D0
        DO J=1,N-1
            S0 = S0+QN(J)
            SN = SN+QN(J)*SJ
            SJ = -SJ
        END DO
        CO(0) = S0/DN
        CO(N) = DD*SN/DN

        SK = -1.D0
        DO  K=1,N-1
            DK = real(K,kind=8)
            SU = .5D0*(QN(0)+QN(N)*SK)
            DO  J=1,N-1
                DJ = real(J,kind=8)
                SU = SU+QN(J)*COS(DK*DJ*M_PI/DN)
            END DO
            CO(K) = 2.D0*SK*SU/DN
            SK = -SK
        END DO

        RETURN
    END SUBROUTINE COCHGL

    SUBROUTINE COLAGR(N,A,ET,QN,WT,CO)
        !**********************************************************************
        !*   COMPUTES THE LAGUERRE FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*   INDIVIDUATED BY ITS VALUES AT THE LAGUERRE GAUSS-RADAU NODES
        !*   N  = THE NUMBER OF NODES
        !*   A  = PARAMETER >-1
        !*   ET = VECTOR OF THE NODES, ET(I), I=0,N-1
        !*   QN = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N-1
        !*   WT = VECTOR OF THE WEIGHTS, WT(I), I=0,N-1
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(0:N-1), INTENT(IN) :: ET, QN, WT
        REAL(kind=8), DIMENSION(0:N-1), INTENT(OUT) :: CO
        INTEGER :: J, K
        REAL(kind=8) :: A1, SU, X, YP, Y, DK, B1, B2, YM, C

        IF (N .EQ. 0) RETURN

        A1  = A+1.D0
        CALL GAMMAF(A1,C)
        SU = 0.D0
        DO J=0,N-1
            SU = SU+QN(J)*WT(J)
            CO(J) = 0.D0
        END DO
        CO(0) = SU/C
        IF (N .EQ. 1) RETURN

        DO J=0,N-1
            X  = ET(J)
            YP = QN(J)*WT(J)
            Y  = (A1-X)*YP
            DO K=1,N-1
                CO(K) = CO(K)+Y
                DK = real(K+1,kind=8)
                B1 = (2.D0*DK+A-1.D0-X)/DK
                B2 = (DK+A-1.D0)/DK
                YM = Y
                Y  = B1*Y-B2*YP
                YP = YM
            END DO
        END DO
        DO K=1,N-1
            DK = real(K,kind=8)
            C  = C*(DK+A)/DK
            CO(K) = CO(K)/C
        END DO

        RETURN
    END SUBROUTINE COLAGR

    SUBROUTINE PVJAEX(N,A,B,X,CO,Y,DY,D2Y)
        !**********************************************************************
        !*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND
        !*   SECOND DERIVATIVES BY KNOWING THE JACOBI FOURIER COEFFICIENTS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   A  = PARAMETER >-1
        !*   B  = PARAMETER >-1
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B, X
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER :: K
        REAL(kind=8) :: AB, P, DP, D2P, PP, DPP, D2PP, DK
        REAL(kind=8) :: CC, C1, C2, C3, C4, PM, DPM, D2PM

        Y   = CO(0)
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .EQ. 0) RETURN

        AB  = A+B
        P   = .5D0*((AB+2.D0)*X+A-B)
        DP  = .5D0*(AB+2.D0)
        D2P = 0.D0
        Y   = CO(0)+P*CO(1)
        DY  = DP*CO(1)
        D2Y = 0.D0
        IF (N .EQ. 1) RETURN

        PP  = 1.D0
        DPP = 0.D0
        D2PP = 0.D0
        DO K=2,N
            DK = real(K,kind=8)
            CC = 2.D0*DK+AB
            C1 = 2.D0*DK*(DK+AB)*(CC-2.D0)
            C2 = (CC-1.D0)*(CC-2.D0)*CC
            C3 = (CC-1.D0)*(A-B)*AB
            C4 = 2.D0*(DK+A-1.D0)*CC*(DK+B-1.D0)
            PM = P
            P  = ((C2*X+C3)*P-C4*PP)/C1
            Y  = Y+P*CO(K)
            PP = PM
            DPM = DP
            DP  = ((C2*X+C3)*DP-C4*DPP+C2*PP)/C1
            DY  = DY+DP*CO(K)
            DPP  = DPM
            D2PM = D2P
            D2P  = ((C2*X+C3)*D2P-C4*D2PP+2.D0*C2*DPP)/C1
            D2Y  = D2Y+D2P*CO(K)
            D2PP = D2PM
        END DO

        RETURN
    END SUBROUTINE PVJAEX

    SUBROUTINE PVLEEX(N,X,CO,Y,DY,D2Y)
        !**********************************************************************
        !*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND
        !*   SECOND DERIVATIVES BY KNOWING THE LEGENDRE FOURIER COEFFICIENTS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER :: K
        REAL(kind=8) :: P, DP, D2P, PP, DPP, D2PP, DK, C2, C4, PM, DPM, D2PM

        Y   = CO(0)
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .EQ. 0) RETURN

        P   = X
        DP  = 1.D0
        D2P = 0.D0
        Y   = CO(0)+P*CO(1)
        DY  = DP*CO(1)
        D2Y = 0.D0
        IF (N .EQ. 1) RETURN

        PP  = 1.D0
        DPP = 0.D0
        D2PP = 0.D0
        DO K=2,N
            DK = real(K,kind=8)
            C2 = 2.D0*DK-1.D0
            C4 = DK-1.D0
            PM = P
            P  = (C2*X*P-C4*PP)/DK
            Y  = Y+P*CO(K)
            PP = PM
            DPM = DP
            DP  = (C2*X*DP-C4*DPP+C2*PP)/DK
            DY  = DY+DP*CO(K)
            DPP  = DPM
            D2PM = D2P
            D2P  = (C2*X*D2P-C4*D2PP+2.D0*C2*DPP)/DK
            D2Y  = D2Y+D2P*CO(K)
            D2PP = D2PM
        END DO

        RETURN
    END SUBROUTINE PVLEEX

    SUBROUTINE PVCHEX(N,X,CO,Y,DY,D2Y)
        !**********************************************************************
        !*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND
        !*   SECOND DERIVATIVES BY KNOWING THE CHEBYSHEV FOURIER COEFFICIENTS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER :: K
        REAL(kind=8) :: P, DP, D2P, PP, DPP, D2PP, PM, DPM, D2PM

        Y   = CO(0)
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .EQ. 0) RETURN

        P   = X
        DP  = 1.D0
        D2P = 0.D0
        Y   = CO(0)+P*CO(1)
        DY  = DP*CO(1)
        D2Y = 0.D0
        IF (N .EQ. 1) RETURN

        PP  = 1.D0
        DPP = 0.D0
        D2PP = 0.D0
        DO K=2,N
            PM = P
            P  = 2.D0*X*P-PP
            Y  = Y+P*CO(K)
            PP = PM
            DPM = DP
            DP  = 2.D0*X*DP+2.D0*PP-DPP
            DY  = DY+DP*CO(K)
            DPP  = DPM
            D2PM = D2P
            D2P  = 2.D0*X*D2P+4.D0*DPP-D2PP
            D2Y  = D2Y+D2P*CO(K)
            D2PP = D2PM
        END DO

        RETURN
    END SUBROUTINE PVCHEX

    SUBROUTINE PVLAEX(N,A,X,CO,Y,DY,D2Y)
        !**********************************************************************
        !*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND
        !*   SECOND DERIVATIVES BY KNOWING THE LAGUERRE FOURIER COEFFICIENTS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   A  = PARAMETER >-1
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: X, A
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER :: K
        REAL(kind=8) :: P, DP, D2P, PP, DPP, D2PP, DK, B1, B2, PM, DPM, D2PM

        Y   = CO(0)
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .EQ. 0) RETURN

        P   = 1.D0+A-X
        DP  = -1.D0
        D2P = 0.D0
        Y   = CO(0)+P*CO(1)
        DY  = DP*CO(1)
        D2Y = 0.D0
        IF (N .EQ. 1) RETURN

        PP  = 1.D0
        DPP = 0.D0
        D2PP = 0.D0
        DO K=2,N
            DK = real(K,kind=8)
            B1 = (2.D0*DK+A-1.D0-X)/DK
            B2 = (DK+A-1.D0)/DK
            PM = P
            P  = B1*P-B2*PP
            Y  = Y+P*CO(K)
            PP = PM
            DPM = DP
            DP  = B1*DP-PP/DK-B2*DPP
            DY  = DY+DP*CO(K)
            DPP  = DPM
            D2PM = D2P
            D2P  = B1*D2P-2.D0*DPP/DK-B2*D2PP
            D2Y  = D2Y+D2P*CO(K)
            D2PP = D2PM
        END DO

        RETURN
    END SUBROUTINE PVLAEX

    SUBROUTINE PVHEEX(N,X,CO,Y,DY,D2Y)
        !**********************************************************************
        !*   COMPUTES THE VALUE OF A POLYNOMIAL OF DEGREE N AND ITS FIRST AND
        !*   SECOND DERIVATIVES BY KNOWING THE HERMITE FOURIER COEFFICIENTS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   X  = THE POINT IN WHICH THE COMPUTATION IS PERFORMED
        !*   CO = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   Y  = VALUE OF THE POLYNOMIAL IN X
        !*   DY = VALUE OF THE FIRST DERIVATIVE IN X
        !*   D2Y= VALUE OF THE SECOND DERIVATIVE IN X
        !**********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: X
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: Y, DY, D2Y
        INTEGER :: K
        REAL(kind=8) :: P, DP, D2P, PP, DPP, D2PP, DK, PM

        Y   = CO(0)
        DY  = 0.D0
        D2Y = 0.D0
        IF (N .EQ. 0) RETURN

        P   = 2.D0*X
        DP  = 2.D0
        D2P = 0.D0
        Y   = CO(0)+P*CO(1)
        DY  = DP*CO(1)
        D2Y = 0.D0
        IF (N .EQ. 1) RETURN

        PP  = 1.D0
        DPP = 0.D0
        D2PP = 0.D0
        DO K=2,N
            DK = real(K,kind=8)
            PM = P
            P  = 2.D0*X*P-2.D0*PP*(DK-1.D0)
            Y  = Y+P*CO(K)
            DY = DY+2.D0*DK*PM*CO(K)
            D2Y = D2Y+4.D0*DK*(DK-1.D0)*PP*CO(K)
            PP  = PM
        END DO

        RETURN
    END SUBROUTINE PVHEEX

    SUBROUTINE NOJAEX(N,A,B,CO,QW)
        !***************************************************************
        !*   COMPUTES THE INTEGRAL NORM OF A POLYNOMIAL BY KNOWING THE
        !*   FOURIER COEFFICIENTS WITH RESPECT TO THE JACOBI BASIS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   A  = PARAMETER >-1
        !*   B  = PARAMETER >-1
        !*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL
        !***************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: QW
        INTEGER :: K
        REAL(kind=8) :: A1, B1, AB, AB2, DN, C, V, SU, GA1, GB1, GAB2, DK
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        A1  = A+1.D0
        B1  = B+1.D0
        AB  = A+B
        AB2 = AB+2.D0
        DN = real(N,kind=8)
        CALL GAMMAF(A1,GA1)
        CALL GAMMAF(B1,GB1)
        CALL GAMMAF(AB2,GAB2)
        C  = ((2.D0)**(AB+1.D0))*GA1*GB1/GAB2
        V  = ABS(CO(0))
        QW = V*SQRT(C)
        IF (N .EQ. 0) RETURN

        SU = 0.D0
        IF (V .LT. EPS) GOTO 1
        SU = C*V*V
1       DO K=1,N
            DK = real(K,kind=8)
            C  = C*(DK+A)*(DK+B)/DK
            V  = ABS(CO(K))
            IF (V .LT. EPS) GOTO 3
            SU = SU+C*V*V/(2.D0*DK+AB+1.D0)
3           C  = C/(DK+AB+1.D0)
        END DO
        QW = SQRT(SU)

        RETURN
    END SUBROUTINE NOJAEX

    SUBROUTINE NOLEEX(N,CO,QI)
        !****************************************************************
        !*   COMPUTES THE INTEGRAL NORM OF A POLYNOMIAL BY KNOWING THE
        !*   FOURIER COEFFICIENTS WITH RESPECT TO THE LEGENDRE BASIS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   QI = INTEGRAL NORM OF THE POLYNOMIAL
        !****************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: QI
        INTEGER :: K
        REAL(kind=8) :: SU, DK, V
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        SU = 0.D0
        DO K=0,N
            DK = real(K,kind=8)
            V  = ABS(CO(K))
            IF (V .LT. EPS) CYCLE
            SU = SU+V*V/(2.D0*DK+1.D0)
        END DO
        QI = SQRT(2.D0*SU)

        RETURN
    END SUBROUTINE NOLEEX

    SUBROUTINE NOCHEX(N,CO,QW,QI)
        !****************************************************************
        !*   COMPUTES THE INTEGRAL NORMS OF A POLYNOMIAL BY KNOWING THE
        !*   FOURIER COEFFICIENTS WITH RESPECT TO THE CHEBYSHEV BASIS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL
        !*   QI = INTEGRAL NORM OF THE POLYNOMIAL
        !****************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: QW, QI
        INTEGER :: K, M
        REAL(kind=8) :: PR, R2, V, SU, D1, D2, C
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        PR = 1.77245385090551588D0
        R2 = 1.41421356237309515D0
        V  = ABS(CO(0))
        QW = PR*V
        QI = R2*V
        IF (N .EQ. 0) RETURN

        SU = 0.D0
        IF (V .GE. EPS) THEN
            SU = 2.D0*V*V
        END IF
        DO K=1,N
            V  = ABS(CO(K))
            IF (V .LT. EPS) CYCLE
            SU = SU+V*V
        END DO
        QW = PR*SQRT(.5D0*SU)

        SU = 0.D0
        DO K=0,N,2
            V  = CO(K)
            DO M=0,N,2
                D1 = 1.D0-real((K+M)*(K+M),kind=8)
                D2 = 1.D0-real((K-M)*(K-M),kind=8)
                C  =1.D0/D1+1.D0/D2
                SU = SU+C*V*CO(M)
            END DO
        END DO
        DO K=1,N,2
            V  = CO(K)
            DO M=1,N,2
                D1 = 1.D0-real((K+M)*(K+M),kind=8)
                D2 = 1.D0-real((K-M)*(K-M),kind=8)
                C  =1.D0/D1+1.D0/D2
                SU = SU+C*V*CO(M)
            END DO
        END DO
        QI = SQRT(SU)

        RETURN
    END SUBROUTINE NOCHEX

    SUBROUTINE NOLAEX(N,A,CO,QW)
        !*****************************************************************
        !*   COMPUTES THE INTEGRAL NORM OF A POLYNOMIAL BY KNOWING THE
        !*   FOURIER COEFFICIENTS WITH RESPECT TO THE LAGUERRE BASIS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   A  = PARAMETER >-1
        !*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL
        !*****************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: QW
        INTEGER :: K
        REAL(kind=8) :: A1, C, V, SU, DK
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        A1  = A+1.D0
        CALL GAMMAF(A1,C)
        V  = ABS(CO(0))
        QW = V*SQRT(C)
        IF (N .EQ. 0) RETURN

        SU = 0.D0
        IF (V .GE. EPS) THEN
            SU = C*V*V
        END IF
        DO K=1,N
            DK = real(K,kind=8)
            C  = C*(DK+A)/DK
            V  = ABS(CO(K))
            IF (V .LT. EPS) CYCLE
            SU = SU+C*V*V
        END DO
        QW = SQRT(SU)

        RETURN
    END SUBROUTINE NOLAEX

    SUBROUTINE NOHEEX(N,CO,QW)
        !****************************************************************
        !*   COMPUTES THE INTEGRAL NORM OF A POLYNOMIAL BY KNOWING THE
        !*   FOURIER COEFFICIENTS WITH RESPECT TO THE HERMITE BASIS
        !*   N  = THE DEGREE OF THE POLYNOMIAL
        !*   CO = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   QW = WEIGHTED INTEGRAL NORM OF THE POLYNOMIAL
        !****************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), INTENT(OUT) :: QW
        INTEGER :: K
        REAL(kind=8) :: PR, R2, V, SU, C, DK
        REAL(kind=8), PARAMETER :: EPS = 1.D-14

        PR = 1.33133536380038953D0
        R2 = 1.41421356237309515D0
        V  = ABS(CO(0))
        QW = V*PR
        IF (N .EQ. 0) RETURN

        SU = 0.D0
        IF (V .GE. EPS) THEN
            SU = V*V
        END IF
        C  = 1.D0
        DO K=1,N
            DK = real(K,kind=8)
            C  = C*R2*SQRT(DK)
            V  = ABS(CO(K))
            IF (V .LT. EPS) CYCLE
            SU = SU+(C*V)*(C*V)
        END DO
        QW = PR*SQRT(SU)

        RETURN
    END SUBROUTINE NOHEEX

    SUBROUTINE COJADE(N,G,CO,CD,CD2)
        !************************************************************************
        !*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
        !*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE JACOBI BASIS
        !*   N   = THE DEGREE OF THE POLYNOMIAL
        !*   G   = PARAMETER >-1
        !*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N
        !*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: G
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CD, CD2
        INTEGER :: K, KR
        REAL(kind=8) :: G2, DK, C1, C2, DN

        CD(N)  = 0.D0
        CD2(N) = 0.D0
        IF (N .EQ. 0) RETURN

        CD(0) = (G+1.D0)*CO(1)
        CD2(N-1) = 0.D0
        IF (N .EQ. 1) RETURN

        DN = real(N,kind=8)
        G2 = 2.D0*G
        CD(N-1) = (2.D0*DN+G2-1.D0)*(DN+G)*CO(N)/(DN+G2)
        DO K=0,N-2
            KR = N-K-2
            IF (KR .NE. 0) THEN
                DK = real(KR,kind=8)
                C1 = (2.D0*DK+G2+1.D0)*(DK+G+1.D0)/(DK+G2+1.D0)
                C2 = (DK+G+2.D0)/((2.D0*DK+G2+5.D0)*(DK+G2+2.D0))
                CD(KR)  = C1*(C2*CD(KR+2)+CO(KR+1))
                CD2(KR) = C1*(C2*CD2(KR+2)+CD(KR+1))
            ELSE
                CD(0) = .25D0*(G+2.D0)*CD(2)/(G+2.5D0)+(G+1.D0)*CO(1)
                CD2(0) = .25D0*(G+2.D0)*CD2(2)/(G+2.5D0)+(G+1.D0)*CD(1)
            ENDIF
        END DO

        RETURN
    END SUBROUTINE COJADE

    SUBROUTINE COLEDE(N,CO,CD,CD2)
        !************************************************************************
        !*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
        !*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE LEGENDRE BASIS
        !*   N   = THE DEGREE OF THE POLYNOMIAL
        !*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N
        !*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CD, CD2
        INTEGER :: K, KR
        REAL(kind=8) :: DN, DK

        CD(N)  = 0.D0
        CD2(N) = 0.D0
        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        CD(N-1)  = (2.D0*DN-1.D0)*CO(N)
        CD2(N-1) = 0.D0
        IF (N .EQ. 1) RETURN

        DO K=0,N-2
            KR = N-K-2
            DK = 2.D0*real(KR,kind=8)+1.D0
            CD(KR)  = DK*(CD(KR+2)/(DK+4.D0)+CO(KR+1))
            CD2(KR) = DK*(CD2(KR+2)/(DK+4.D0)+CD(KR+1))
        END DO

        RETURN
    END SUBROUTINE COLEDE

    SUBROUTINE COCHDE(N,CO,CD,CD2)
        !************************************************************************
        !*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
        !*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE CHEBYSHEV BASIS
        !*   N   = THE DEGREE OF THE POLYNOMIAL
        !*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N
        !*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CD, CD2
        INTEGER :: K, KR
        REAL(kind=8) :: DN, DK

        CD(N)  = 0.D0
        CD2(N) = 0.D0
        IF (N .EQ. 0) RETURN

        CD(0) = CO(1)
        CD2(N-1) = 0.D0
        IF (N .EQ. 1) RETURN

        DN = real(N,kind=8)
        CD(N-1) = 2.D0*DN*CO(N)
        DO K=0,N-2
            KR = N-K-2
            IF (KR .NE. 0) THEN
                DK = 2.D0*(real(KR,kind=8)+1.D0)
                CD(KR)  = CD(KR+2)+DK*CO(KR+1)
                CD2(KR) = CD2(KR+2)+DK*CD(KR+1)
            ELSE
                CD(0)  = .5D0*CD(2)+CO(1)
                CD2(0) = .5D0*CD2(2)+CD(1)
            ENDIF
        END DO

        RETURN
    END SUBROUTINE COCHDE

    SUBROUTINE COLADE(N,CO,CD,CD2)
        !************************************************************************
        !*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
        !*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE LAGUERRE BASIS
        !*   N   = THE DEGREE OF THE POLYNOMIAL
        !*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N
        !*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CD, CD2
        INTEGER :: K, KR

        CD(N)  = 0.D0
        CD2(N) = 0.D0
        IF (N .EQ. 0) RETURN

        CD(N-1)  = -CO(N)
        CD2(N-1) = 0.D0
        IF (N .EQ. 1) RETURN

        DO K=0,N-2
            KR = N-K-2
            CD(KR)  = CD(KR+1)-CO(KR+1)
            CD2(KR) = CD2(KR+2)-CD(KR+1)
        END DO

        RETURN
    END SUBROUTINE COLADE

    SUBROUTINE COHEDE(N,CO,CD,CD2)
        !************************************************************************
        !*   COMPUTES THE FOURIER COEFFICIENTS OF THE DERIVATIVES OF A POLYNOMIAL
        !*   FROM ITS FOURIER COEFFICIENTS WITH RESPECT TO THE HERMITE BASIS
        !*   N   = THE DEGREE OF THE POLYNOMIAL
        !*   CO  = COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*   CD  = COEFFICIENTS OF THE FIRST DERIVATIVE, CD(I), I=0,N
        !*   CD2 = COEFFICIENTS OF THE SECOND DERIVATIVE, CD2(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CD, CD2
        INTEGER :: K, KR
        REAL(kind=8) :: DN, DK

        CD(N)  = 0.D0
        CD2(N) = 0.D0
        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        CD(N-1)  = 2.D0*DN*CO(N)
        CD2(N-1) = 0.D0
        IF (N .EQ. 1) RETURN

        DO K=0,N-2
            KR = N-K-2
            DK = 2.D0*real(KR,kind=8)+2.D0
            CD(KR)  = DK*CO(KR+1)
            CD2(KR) = DK*CD(KR+1)
        END DO

        RETURN
    END SUBROUTINE COHEDE

    SUBROUTINE DEJAGA(N,A,B,CS,DZ,QZ,DQZ)
        !************************************************************************
        !*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE JACOBI ZEROES FROM
        !*   THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS
        !*   N   = THE NUMBER OF ZEROES
        !*   A   = PARAMETER >-1
        !*   B   = PARAMETER >-1
        !*   CS  = ZEROES OF THE JACOBI POLYNOMIAL, CS(I), I=1,N
        !*   DZ  = JACOBI DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*   QZ  = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   DQZ = DERIVATIVES OF THE POLYNOMIAL AT THE ZEROES, DQZ(I), I=1,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, DZ, QZ
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: DQZ
        INTEGER :: I, J
        REAL(kind=8) :: SU, CI, DI, CJ, DJ

        IF (N .EQ. 0) RETURN

        DO I=1,N
            SU = 0.D0
            CI = CS(I)
            DI = DZ(I)
            DO J=1,N
                IF (I .NE. J) THEN
                    CJ = CS(J)
                    DJ = DZ(J)
                    SU = SU+QZ(J)/(DJ*(CI-CJ))
                ELSE
                    SU = SU+.5D0*QZ(I)*((A+B+2.D0)*CI+A-B)/(DI*(1.D0-CI*CI))
                ENDIF
            END DO
            DQZ(I) = DI*SU
        END DO

        RETURN
    END SUBROUTINE DEJAGA

    SUBROUTINE DELAGA(N,A,CS,QZ,DQZ)
        !************************************************************************
        !*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE LAGUERRE ZEROES FROM
        !*   THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS
        !*   N   = THE NUMBER OF ZEROES
        !*   A   = PARAMETER >-1
        !*   CS  = ZEROES OF THE LAGUERRE POLYNOMIAL, CS(I), I=1,N
        !*   QZ  = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   DQZ = DERIVATIVES OF THE POLYNOMIAL AT THE ZEROES, DQZ(I), I=1,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS
        REAL(kind=8), DIMENSION(N), INTENT(INOUT) :: QZ
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: DQZ
        INTEGER :: I, J
        REAL(kind=8) :: Y, DY, D2Y, SU, CI, DI, CJ

        IF (N .EQ. 0) RETURN

        DO J=1,N
            CJ = CS(J)
            CALL VALAPO(N,A,CJ,Y,DY,D2Y)
            QZ(J) = QZ(J)/DY
        END DO

        DO I=1,N
            SU = 0.D0
            CI = CS(I)
            CALL VALAPO(N,A,CI,Y,DI,D2Y)
            DO J=1,N
                IF (I .NE. J) THEN
                    CJ = CS(J)
                    SU = SU+QZ(J)/(CI-CJ)
                ELSE
                    SU = SU+.5D0*QZ(I)*(CI-A-1.D0)/CI
                ENDIF
            END DO
            DQZ(I) = DI*SU
        END DO

        DO I=1,N
            CI = CS(I)
            CALL VALAPO(N,A,CI,Y,DY,D2Y)
            QZ(I) = DY*QZ(I)
        END DO

        RETURN
    END SUBROUTINE DELAGA

    SUBROUTINE DEHEGA(N,CS,QZ,DQZ)
        !************************************************************************
        !*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE HERMITE ZEROES FROM
        !*   THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS
        !*   N   = THE NUMBER OF ZEROES
        !*   CS  = ZEROES OF THE HERMITE POLYNOMIAL, CS(I), I=1,N
        !*   QZ  = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*   DQZ = DERIVATIVES OF THE POLYNOMIAL AT THE ZEROES, DQZ(I), I=1,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS
        REAL(kind=8), DIMENSION(N), INTENT(INOUT) :: QZ
        REAL(kind=8), DIMENSION(N), INTENT(OUT) :: DQZ
        INTEGER :: I, J
        REAL(kind=8) :: Y, DY, D2Y, SU, CI, DI, CJ

        IF (N .EQ. 0) RETURN

        DO J=1,N
            CJ = CS(J)
            CALL VAHEPO(N,CJ,Y,DY,D2Y)
            QZ(J) = QZ(J)/DY
        END DO

        DO I=1,N
            SU = 0.D0
            CI = CS(I)
            CALL VAHEPO(N,CI,Y,DI,D2Y)
            DO J=1,N
                IF (I .NE. J) THEN
                    CJ = CS(J)
                    SU = SU+QZ(J)/(CI-CJ)
                ELSE
                    SU = SU+CI*QZ(I)
                ENDIF
            END DO
            DQZ(I) = DI*SU
        END DO

        DO I=1,N
            CI = CS(I)
            CALL VAHEPO(N,CI,Y,DY,D2Y)
            QZ(I) = DY*QZ(I)
        END DO

        RETURN
    END SUBROUTINE DEHEGA

    SUBROUTINE DEJAGL(N,A,B,ET,VN,QN,DQN)
        !************************************************************************
        !*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE JACOBI GAUSS-LOBATTO
        !*   NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS
        !*   N   = THE DEGREE OF THE POLYNOMIAL
        !*   A   = PARAMETER >-1
        !*   B   = PARAMETER >-1
        !*   ET  = VECTOR OF THE NODES, ET(I), I=0,N
        !*   VN  = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A, B
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, VN, QN
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: DQN
        INTEGER :: I, J
        REAL(kind=8) :: A1, B1, AB, DN, C1, C2, S1, S2, C3, C4
        REAL(kind=8) :: VI, VJ, EI, EJ, SU

        DQN(0) = 0.D0
        IF (N .EQ. 0) RETURN

        A1 = A+1.D0
        B1 = B+1.D0
        AB = A+B
        DN = real(N,kind=8)
        C1 = DN*(DN+AB+1.D0)
        C2 = A1*VN(0)/(B1*VN(N))
        DQN(0) = .5D0*((A-C1)*QN(0)/(B+2.D0)-C2*QN(N))
        DQN(N) = .5D0*(QN(0)/C2+(C1-B)*QN(N)/(A+2.D0))
        IF (N .EQ. 1) RETURN

        S1 = DQN(0)
        S2 = DQN(N)
        C3 = -VN(0)/B1
        C4 = VN(N)/A1
        DO J=1,N-1
            VJ = QN(J)/VN(J)
            EJ = ET(J)
            S1 = S1+C3*VJ/(1.D0+EJ)
            S2 = S2+C4*VJ/(1.D0-EJ)
        END DO
        DQN(0) = S1
        DQN(N) = S2

        DO I=1,N-1
            VI = VN(I)
            EI = ET(I)
            SU = B1*QN(0)/((1.D0+EI)*VN(0))-A1*QN(N)/((1.D0-EI)*VN(N))
            DO J=1,N-1
                IF (I .NE. J) THEN
                    VJ = VN(J)
                    EJ = ET(J)
                    SU = SU+QN(J)/(VJ*(EI-EJ))
                ELSE
                    SU = SU+.5D0*QN(I)*(AB*EI+A-B)/(VI*(1.D0-EI*EI))
                ENDIF
            END DO
            DQN(I) = VI*SU
        END DO

        RETURN
    END SUBROUTINE DEJAGL

    SUBROUTINE DELEGL(N,ET,VN,QN,DQN)
        !************************************************************************
        !*  COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE LEGENDRE GAUSS-LOBATTO
        !*  NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS
        !*   N   = THE DEGREE OF THE POLYNOMIAL
        !*   ET  = VECTOR OF THE NODES, ET(I), I=0,N
        !*   VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, VN, QN
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: DQN
        INTEGER :: I, J
        REAL(kind=8) :: DN, C, SU, VI, VJ, EI, EJ

        DQN(0) = 0.D0
        IF (N .EQ. 0) RETURN

        DO I=0,N
            SU = 0.D0
            VI = VN(I)
            EI = ET(I)
            DO J=0,N
                IF (I .EQ. J) CYCLE
                VJ = VN(J)
                EJ = ET(J)
                SU = SU+QN(J)/(VJ*(EI-EJ))
            END DO
            DQN(I) = VI*SU
        END DO

        DN = real(N,kind=8)
        C  = .25D0*DN*(DN+1.D0)
        DQN(0) = DQN(0)-C*QN(0)
        DQN(N) = DQN(N)+C*QN(N)

        RETURN
    END SUBROUTINE DELEGL

    SUBROUTINE DECHGL(N,ET,QN,DQN)
        !************************************************************************
        !* COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE CHEBYSHEV GAUSS-LOBATTO
        !* NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS
        !*   N   = THE DEGREE OF THE POLYNOMIAL
        !*   ET  = VECTOR OF THE NODES, ET(I), I=0,N
        !*   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, QN
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: DQN
        INTEGER :: I, J
        REAL(kind=8) :: DN, CN, SN, S1, S2, SGN, SGNI, EI, SU, SGNJ, EJ, QJ

        DQN(0) = 0.D0
        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        CN = (2.D0*DN*DN+1.D0)/6.D0
        SN = real(1+4*(N/2)-2*N,kind=8)
        DQN(0) = -CN*QN(0)-.5D0*SN*QN(N)
        DQN(N) = .5D0*SN*QN(0)+CN*QN(N)
        IF (N .EQ. 1) RETURN

        S1 = DQN(0)
        S2 = DQN(N)
        SGN = -1.D0
        DO J=1,N-1
            EJ = ET(J)
            QJ = 2.D0*SGN*QN(J)
            S1 = S1-QJ/(1.D0+EJ)
            S2 = S2+QJ*SN/(1.D0-EJ)
            SGN = -SGN
        END DO
        DQN(0) = S1
        DQN(N) = S2

        SGNI = -1.D0
        DO I=1,N-1
            EI = ET(I)
            SU = .5D0*SGNI*(QN(0)/(1.D0+EI)-SN*QN(N)/(1.D0-EI))
            SGNJ = -1.D0
            DO J=1,N-1
                IF (I .NE. J) THEN
                    EJ = ET(J)
                    SU = SU+SGNI*SGNJ*QN(J)/(EI-EJ)
                ELSE
                    SU = SU-.5D0*EI*QN(I)/(1.D0-EI*EI)
                ENDIF
                SGNJ = -SGNJ
            END DO
            DQN(I) = SU
            SGNI = -SGNI
        END DO

        RETURN
    END SUBROUTINE DECHGL

    SUBROUTINE DELAGR(N,A,ET,QN,DQN)
        !************************************************************************
        !*   COMPUTES THE DERIVATIVE OF A POLYNOMIAL AT THE LAGUERRE GAUSS-RADAU
        !*   NODES FROM THE VALUES OF THE POLYNOMIAL ATTAINED AT THE SAME POINTS
        !*   N   = THE NUMBER OF NODES
        !*   A   = PARAMETER >-1
        !*   ET  = VECTOR OF THE NODES, ET(I), I=0,N-1
        !*   QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N-1
        !*   DQN = DERIVATIVES OF THE POLYNOMIAL AT THE NODES, DQZ(I), I=0,N-1
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(0:N-1), INTENT(IN) :: ET
        REAL(kind=8), DIMENSION(0:N-1), INTENT(INOUT) :: QN
        REAL(kind=8), DIMENSION(0:N-1), INTENT(OUT) :: DQN
        INTEGER :: I, J
        REAL(kind=8) :: DN, A1, SU, X, C, Y, DY, D2Y, EI, EJ

        DN = real(N,kind=8)
        DQN(0) = (1.D0-DN)*QN(0)/(A+2.D0)
        IF (N .EQ. 1) RETURN

        A1 = A+1.D0
        SU = DQN(0)
        X  = 0.D0
        CALL VALAPO(N,A,X,Y,DY,D2Y)
        C = Y
        DO J=1,N-1
            EJ = ET(J)
            CALL VALAPO(N,A,EJ,Y,DY,D2Y)
            QN(J) = QN(J)/Y
            SU = SU-C*QN(J)/(A1*EJ)
        END DO
        DQN(0) = SU

        DO I=1,N-1
            EI = ET(I)
            CALL VALAPO(N,A,EI,Y,DY,D2Y)
            SU = A1*QN(0)/(C*EI)
            DO J=1,N-1
                IF (I .NE. J) THEN
                    EJ = ET(J)
                    SU = SU+QN(J)/(EI-EJ)
                ELSE
                    SU = SU+.5D0*QN(I)*(EI-A)/EI
                ENDIF
            END DO
            DQN(I) = Y*SU
        END DO

        DO J=1,N-1
            EJ = ET(J)
            CALL VALAPO(N,A,EJ,Y,DY,D2Y)
            QN(J) = Y*QN(J)
        END DO

        RETURN
    END SUBROUTINE DELAGR

    SUBROUTINE DMJAGL(N,NM,A,B,ET,VN,DMA)
        !************************************************************************
        !*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
        !*  JACOBI GAUSS-LOBATTO NODES
        !*  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX
        !*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
        !*  A   = PARAMETER >-1
        !*  B   = PARAMETER >-1
        !*  ET  = VECTOR OF THE NODES, ET(I), I=0,N
        !*  VN  = VALUES OF THE JACOBI POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, NM
        REAL(kind=8), INTENT(IN) :: A, B
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, VN
        REAL(kind=8), DIMENSION(0:NM,0:N), INTENT(OUT) :: DMA
        INTEGER :: I, J
        REAL(kind=8) :: A1, B1, AB, DN, C1, C2, C3, C4, VJ, EJ, VI, EI


        DMA(0,0) = 0.D0
        IF (N .EQ. 0) RETURN

        A1 = A+1.D0
        B1 = B+1.D0
        AB = A+B
        DN = real(N,kind=8)
        C1 = DN*(DN+AB+1.D0)
        C2 = A1*VN(0)/(B1*VN(N))
        DMA(0,0) = .5D0*(A-C1)/(B+2.D0)
        DMA(N,N) = .5D0*(C1-B)/(A+2.D0)
        DMA(0,N) = -.5D0*C2
        DMA(N,0) = .5D0/C2
        IF (N .EQ. 1) RETURN

        C3 = VN(0)/B1
        C4 = VN(N)/A1
        DO J=1,N-1
            VJ = VN(J)
            EJ = ET(J)
            DMA(0,J) = -C3/(VJ*(1.D0+EJ))
            DMA(N,J) = C4/(VJ*(1.D0-EJ))
            DMA(J,0) = VJ/(C3*(1.D0+EJ))
            DMA(J,N) = -VJ/(C4*(1.D0-EJ))
        END DO

        DO I=1,N-1
            VI = VN(I)
            EI = ET(I)
            DO J=1,N-1
                IF (I .NE. J) THEN
                    VJ = VN(J)
                    EJ = ET(J)
                    DMA(I,J) = VI/(VJ*(EI-EJ))
                ELSE
                    DMA(I,I) = .5D0*(AB*EI+A-B)/(1.D0-EI*EI)
                ENDIF
            END DO
        END DO

        RETURN
    END SUBROUTINE DMJAGL

    SUBROUTINE DMLEGL(N,NM,ET,VN,DMA)
        !************************************************************************
        !*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
        !*  LEGENDRE GAUSS-LOBATTO NODES
        !*  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX
        !*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
        !*  ET  = VECTOR OF THE NODES, ET(I), I=0,N
        !*  VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, NM
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, VN
        REAL(kind=8), DIMENSION(0:NM,0:N), INTENT(OUT) :: DMA
        INTEGER :: I, J
        REAL(kind=8) :: DN, C, VJ, EJ, VI, EI

        DMA(0,0) = 0.D0
        IF (N .EQ. 0) RETURN

        DO I=0,N
            VI = VN(I)
            EI = ET(I)
            DO J=0,N
                IF (I .NE. J) THEN
                    VJ = VN(J)
                    EJ = ET(J)
                    DMA(I,J) = VI/(VJ*(EI-EJ))
                ELSE
                    DMA(I,I) = 0.D0
                ENDIF
            END DO
        END DO

        DN = real(N,kind=8)
        C  = .25D0*DN*(DN+1.D0)
        DMA(0,0) = -C
        DMA(N,N) = C

        RETURN
    END SUBROUTINE DMLEGL

    SUBROUTINE DMCHGL(N,NM,ET,DMA)
        !************************************************************************
        !*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
        !*  CHEBYSHEV GAUSS-LOBATTO NODES
        !*  N   = PARAMETER RELATIVE TO THE DIMENSION OF THE MATRIX
        !*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
        !*  ET  = VECTOR OF THE NODES, ET(I), I=0,N
        !*  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N  J=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, NM
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET
        REAL(kind=8), DIMENSION(0:NM,0:N), INTENT(OUT) :: DMA
        INTEGER :: I, J
        REAL(kind=8) :: DN, CN, SN, SGN, SGNI, SGNJ, EJ, EI

        DMA(0,0) = 0.D0
        IF (N .EQ. 0) RETURN

        DN = real(N,kind=8)
        CN = (2.D0*DN*DN+1.D0)/6.D0
        SN = real(1+4*(N/2)-2*N,kind=8)
        DMA(0,0) = -CN
        DMA(N,N) = CN
        DMA(0,N) = -.5D0*SN
        DMA(N,0) = .5D0*SN
        IF (N .EQ. 1) RETURN

        SGN = -1.D0
        DO J=1,N-1
            EJ = ET(J)
            DMA(0,J) = -2.D0*SGN/(1.D0+EJ)
            DMA(N,J) = 2.D0*SGN*SN/(1.D0-EJ)
            DMA(J,0) = .5D0*SGN/(1.D0+EJ)
            DMA(J,N) = -.5D0*SGN*SN/(1.D0-EJ)
            SGN = -SGN
        END DO

        SGNI = -1.D0
        DO I=1,N-1
            EI = ET(I)
            SGNJ = -1.D0
            DO J=1,N-1
                IF (I .NE. J) THEN
                    EJ = ET(J)
                    DMA(I,J) = SGNI*SGNJ/(EI-EJ)
                ELSE
                    DMA(I,I) = -.5D0*EI/(1.D0-EI*EI)
                ENDIF
                SGNJ = -SGNJ
            END DO
            SGNI = -SGNI
        END DO

        RETURN
    END SUBROUTINE DMCHGL

    SUBROUTINE DMLAGR(N,NM,A,ET,DMA)
        !************************************************************************
        !*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
        !*  LAGUERRE GAUSS-RADAU NODES
        !*  N   = DIMENSION OF THE MATRIX
        !*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
        !*  A   = PARAMETER >-1
        !*  ET  = VECTOR OF THE NODES, ET(I), I=0,N-1
        !*  DMA = DERIVATIVE MATRIX, DMA(I,J), I=0,N-1  J=0,N-1
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, NM
        REAL(kind=8), INTENT(IN) :: A
        REAL(kind=8), DIMENSION(0:N-1), INTENT(IN) :: ET
        REAL(kind=8), DIMENSION(0:NM,0:N-1), INTENT(OUT) :: DMA
        INTEGER :: I, J
        REAL(kind=8) :: DN, X, A1, Y, DY, D2Y, C, EI, EJ

        DN = real(N,kind=8)
        DMA(0,0) = (1.D0-DN)/(A+2.D0)
        IF (N .EQ. 1) RETURN

        A1 = A+1.D0
        X  = 0.D0
        CALL VALAPO(N,A,X,Y,DY,D2Y)
        C = Y
        DO J=1,N-1
            EJ = ET(J)
            CALL VALAPO(N,A,EJ,Y,DY,D2Y)
            DMA(0,J) = -C/(A1*EJ*Y)
            DMA(J,0) = A1*Y/(C*EJ)
        END DO

        DO I=1,N-1
            EI = ET(I)
            CALL VALAPO(N,A,EI,Y,DY,D2Y)
            DO J=1,N-1
                IF (I .NE. J) THEN
                    EJ = ET(J)
                    DMA(I,J) = Y/(EI-EJ)
                ELSE
                    DMA(I,I) = .5D0*(EI-A)/EI
                ENDIF
            END DO
        END DO

        DO J=1,N-1
            EJ = ET(J)
            CALL VALAPO(N,A,EJ,Y,DY,D2Y)
            DO I=1,N-1
                IF (I .EQ. J) CYCLE
                DMA(I,J) = DMA(I,J)/Y
            END DO
        END DO

        RETURN
    END SUBROUTINE DMLAGR

    SUBROUTINE FCCHGA(N,QZ,CO)
        !************************************************************************
        !*  COMPUTES USING FFT THE CHEBYSHEV FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*  INDIVIDUATED BY THE VALUES ATTAINED AT THE CHEBYSHEV ZEROES
        !*  N   = THE NUMBER OF ZEROES
        !*  QZ  = VALUES OF THE POLYNOMIAL AT THE ZEROES, QZ(I), I=1,N
        !*  CO  = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N-1
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: QZ
        REAL(kind=8), DIMENSION(0:N-1), INTENT(OUT) :: CO
        INTEGER :: I, M, N2, IFAIL
        REAL(kind=8) :: R2, DN, C, SN, SM, AR, CS, SI, V1, V2

        IF(N .EQ. 0) RETURN

        R2 = 1.41421356237309515D0
        DN = real(N,kind=8)
        CO(0) = QZ(1)/DN
        IF(N .EQ. 1) RETURN

        N2 = N/2
        C  = 2.D0/SQRT(DN)
        SN = real(1+4*N2-2*N,kind=8)
        CO(N-N2-1) = QZ(N)
        DO I=1,N2
            CO(I-1) = QZ(2*I-1)
            CO(N-I) = QZ(2*I)
        END DO

        CALL C06EAF(CO,N,IFAIL)
        IF (IFAIL .NE. 0) THEN
            WRITE(*,*) 'IFAIL IS NOT ZERO IN SUBROUTINE C06EAF'
        ENDIF

        CO(0) = .5D0*C*CO(0)
        IF (2*N2 .EQ. N) THEN
            CO(N2) = C*((-1.D0)**N2)*CO(N2)/R2
        ENDIF
        IF (N .EQ. 2) RETURN
        SM = -1.D0
        DO M=1,N-N2-1
            AR = M_PI_2*real(M,kind=8)/DN
            CS = COS(AR)
            SI = SIN(AR)
            V1 = C*SM*(CO(M)*CS+CO(N-M)*SI)
            V2 = C*SM*SN*(CO(M)*SI-CO(N-M)*CS)
            CO(M) = V1
            CO(N-M) = V2
            SM = -SM
        END DO

        RETURN
    END SUBROUTINE FCCHGA

    SUBROUTINE FCCHGL(N,QN,CO)
        !************************************************************************
        !*  COMPUTES USING FFT THE CHEBYSHEV FOURIER COEFFICIENTS OF A POLYNOMIAL
        !*  INDIVIDUATED BY ITS VALUES AT THE CHEBYSHEV GAUSS-LOBATTO NODES
        !*  N   = THE DEGREE OF THE POLYNOMIAL
        !*  QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*  CO  = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: QN
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CO
        INTEGER :: I, J, M, N2, I2, IFAIL
        REAL(kind=8) :: DN, SN, S1, S2, SJ, C, SM, V1, SI, V2, AR

        IF(N .EQ. 0) RETURN

        DN = real(N,kind=8)
        N2 = N/2
        SN = real(1+4*N2-2*N,kind=8)
        S1 = .5D0*(QN(0)+QN(N))
        S2 = .5D0*(SN*QN(0)+QN(N))
        CO(0) = S1
        CO(N) = S2
        IF(N .EQ. 1) RETURN

        CO(0) = QN(0)
        IF (2*N2 .EQ. N) THEN
            CO(N2) = QN(N)
        ENDIF
        IF(N .EQ. 2) GOTO 2

        DO I=1,N-N2-1
            I2 = 2*I
            CO(I) = QN(I2)
            CO(N-I) = QN(I2-1)-QN(I2+1)
        END DO

2       CALL C06EBF(CO,N,IFAIL)
        IF (IFAIL .NE. 0) THEN
            WRITE(*,*) 'IFAIL IS NOT ZERO IN SUBROUTINE C06EBF'
        ENDIF

        SJ = -1.D0
        DO J=1,N-1
            S1 = S1+QN(J)
            S2 = S2+SJ*SN*QN(J)
            SJ = -SJ
        END DO
        CO(0) = S1/DN
        CO(N) = S2/DN
        C  = .5D0/SQRT(DN)
        IF (2*N2 .EQ. N) THEN
            CO(N2) = 2.D0*C*((-1.D0)**N2)*CO(N2)
        ENDIF
        IF(N .EQ. 2) RETURN

        SM = -1.D0
        DO M=1,N-N2-1
            AR = M_PI*real(M,kind=8)/DN
            SI = .5D0/SIN(AR)
            V1 = CO(M)*(1.D0+SI)+CO(N-M)*(1.D0-SI)
            V2 = CO(M)*(1.D0-SI)+CO(N-M)*(1.D0+SI)
            CO(M) = C*SM*V1
            CO(N-M) = C*SM*SN*V2
            SM = -SM
        END DO

        RETURN
    END SUBROUTINE FCCHGL

    SUBROUTINE FVCHGL(N,CO,QN)
        !*******************************************************************
        !*  COMPUTES USING FFT THE VALUES OF A POLYNOMIAL AT THE CHEBYSHEV
        !*  GAUSS-LOBATTO NODES FROM ITS CHEBYSHEV FOURIER COEFFICIENTS
        !*  N   = THE DEGREE OF THE POLYNOMIAL
        !*  CO  = FOURIER COEFFICIENTS OF THE POLYNOMIAL, CO(I), I=0,N
        !*  QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*******************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: CO
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: QN
        INTEGER :: I, J, N2, I2, IFAIL, M
        REAL(kind=8) :: DN, SN, S1, S2, SJ, C, V1, SI, V2, AR

        IF(N .EQ. 0) RETURN

        DN = real(N,kind=8)
        N2 = N/2
        SN = real(1+4*N2-2*N,kind=8)
        S1 = CO(0)+SN*CO(N)
        S2 = CO(0)+CO(N)
        QN(0) = S1
        QN(N) = S2
        IF(N .EQ. 1) RETURN

        QN(0) = 2.D0*CO(0)
        IF (2*N2 .EQ. N) THEN
            QN(N2) = 2.D0*CO(N)
        ENDIF
        IF(N .EQ. 2) GOTO 2

        DO I=1,N-N2-1
            I2 = 2*I
            QN(N-I) = CO(I2+1)-CO(I2-1)
            QN(I) = CO(I2)
        END DO
        IF (2*N2 .NE. N) THEN
            QN(N2+1) = 2.D0*CO(N)-CO(N-2)
        ENDIF

2       CALL C06EBF(QN,N,IFAIL)
        IF (IFAIL .NE. 0) THEN
            WRITE(*,*) 'IFAIL IS NOT ZERO IN SUBROUTINE C06EBF'
        ENDIF

        SJ = -1.D0
        DO J=1,N-1
            S1 = S1+SJ*CO(J)
            S2 = S2+CO(J)
            SJ = -SJ
        END DO
        QN(0) = S1
        QN(N) = S2
        C  = .25D0*SQRT(DN)
        IF (2*N2 .EQ. N) THEN
            QN(N2) = 2.D0*C*QN(N2)
        ENDIF
        IF(N .EQ. 2) RETURN

        DO M=1,N-N2-1
            AR = M_PI*real(M,kind=8)/DN
            SI = .5D0/SIN(AR)
            V1 = QN(M)*(1.D0+SI)+QN(N-M)*(1.D0-SI)
            V2 = QN(M)*(1.D0-SI)+QN(N-M)*(1.D0+SI)
            QN(M) = C*V1
            QN(N-M) = C*V2
        END DO

        RETURN
    END SUBROUTINE FVCHGL

    SUBROUTINE FDCHGL(N,QN,CD,DQN)
        !*********************************************************************
        !*  COMPUTES USING FFT THE FOURIER COEFFICIENTS AND THE VALUES AT THE
        !*  CHEBYSHEV GAUSS-LOBATTO NODES OF THE DERIVATIVE OF A POLYNOMIAL
        !*  FROM THE VALUES ATTAINED BY THE POLYNOMIAL AT THE NODES
        !*  N   = THE DEGREE OF THE POLYNOMIAL
        !*  QN  = VALUES OF THE POLYNOMIAL AT THE NODES, QN(I), I=0,N
        !*  CD  = FOURIER COEFFICIENTS OF THE DERIVATIVE, CD(I), I=0,N
        !*  DQN = VALUES OF THE DERIVATIVE AT THE NODES, DQN(I), I=0,N
        !*********************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: QN
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: CD, DQN
        INTEGER :: K, KR
        REAL(kind=8) :: DN, DK

        IF(N .EQ. 0) RETURN

        CALL FCCHGL(N,QN,DQN)
        CD(N) = 0.D0
        CD(0) = DQN(1)
        IF(N .NE. 1) THEN
            DN = real(N,kind=8)
            CD(N-1) = 2.D0*DN*DQN(N)
            DO K=0,N-2
                KR = N-K-2
                IF(KR .NE. 0) THEN
                    DK = 2.D0*(real(KR,kind=8)+1.D0)
                    CD(KR) = CD(KR+2)+DK*DQN(KR+1)
                ELSE
                    CD(0) = .5D0*CD(2)+DQN(1)
                ENDIF
            END DO
        END IF

        CALL FVCHGL(N,CD,DQN)

        RETURN
    END SUBROUTINE FDCHGL

    SUBROUTINE EDLEGL(N,ET,VN,FU,S1,S2,PA,SO,DSO,D2SO)
        !************************************************************************
        !*  COMPUTES THE APPROXIMATE SOLUTION OF A DIRICHLET PROBLEM BY
        !*  COLLOCATION AT THE LEGENDRE GAUSS-LOBATTO NODES
        !*  N   = DEGREE OF THE APPROXIMATING POLYNOMIAL
        !*  ET  = VECTOR OF THE NODES, ET(I), I=0,N
        !*  VN  = VALUES OF THE LEGENDRE POLYNOMIAL AT THE NODES, VN(I), I=0,N
        !*  FU  = VALUES OF THE RIGHT-HAND SIDE AT THE NODES, FU(I), I=0,N
        !*  S1  = DIRICHLET DATUM AT X=-1
        !*  S2  = DIRICHLET DATUM AT X=1
        !*  PA  = PARAMETER >0
        !*  SO  = VALUES OF THE APPROXIMATED SOLUTION AT THE NODES, SO(I), I=0,N
        !*  DSO = DERIVATIVE OF THE SOLUTION AT THE NODES, DSO(I), I=0,N
        !*  D2SO= SECOND DERIVATIVE OF THE SOLUTION AT THE NODES, D2SO(I), I=0,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL(kind=8), INTENT(IN) :: S1, S2, PA
        REAL(kind=8), DIMENSION(0:N), INTENT(IN) :: ET, VN, FU
        REAL(kind=8), DIMENSION(0:N), INTENT(OUT) :: SO, DSO, D2SO
        INTEGER :: I, J, IR, IT, IM
        REAL(kind=8) :: C1, C2, C3, C4, TH, BE, D1, D2, D3

        IF (N .EQ. 0) RETURN

        SO(0) = S1
        SO(N) = S2
        C1 = .5D0*(S2-S1)
        DSO(0) = C1
        DSO(N) = C1
        D2SO(0) = 0.D0
        D2SO(N) = 0.D0
        IF (N .EQ. 1) RETURN

        DO J=1,N-1
            SO(J) = .5D0*(S1+S2)+C1*ET(J)
        END DO

        TH = .58D0
        IM = int(42.D0+DLOG(1.D0+PA))
        DO IT=1,IM
            CALL DELEGL(N,ET,VN,SO,DSO)
            CALL DELEGL(N,ET,VN,DSO,D2SO)
            DO I=1,N-1
                D2SO(I) = D2SO(I)+FU(I)-PA*SO(I)
            END DO

            D2SO(N) = 0.D0
            C2 = 2.D0/(ET(1)-ET(0))
            BE = C2/(ET(2)-ET(0))
            DSO(1) = D2SO(1)
            C3 = 0.D0
            IF(N .EQ. 2) GOTO 5

            DO I=1,N-2
                D1 = ET(I)-ET(I-1)
                D2 = ET(I+1)-ET(I)
                D3 = ET(I+2)-ET(I+1)
                C4 = -BE*C3+2.D0/(D1*D2)+PA
                D2SO(I) = C4
                BE = 2.D0/(C4*D2*(D2+D3))
                DSO(I+1) = D2SO(I+1)+BE*DSO(I)
                C3 = 2.D0/(D2*(D1+D2))
            END DO

5           D2SO(N-1) = -BE*C3+C2/(ET(2)-ET(1))+PA
            DSO(N) = 0.D0
            DO I=1,N-1
                IR = N-I
                D1 = ET(IR)-ET(IR-1)
                D2 = ET(IR+1)-ET(IR)
                DSO(IR) = (DSO(IR)+2.D0*DSO(IR+1)/(D2*(D1+D2)))/D2SO(IR)
            END DO

            DO I=1,N-1
                SO(I) = SO(I)+TH*DSO(I)
            END DO

        END DO

        CALL DELEGL(N,ET,VN,SO,DSO)
        CALL DELEGL(N,ET,VN,DSO,D2SO)

        RETURN
    END SUBROUTINE EDLEGL

    SUBROUTINE DMLEGA(N,NM,CS,DZ,DGA)
        !************************************************************************
        !*  COMPUTES THE ENTRIES OF THE DERIVATIVE MATRIX RELATIVE TO THE
        !*  LEGENDRE GAUSS NODES
        !*  N   = DIMENSION OF THE MATRIX
        !*  NM  = ORDER OF THE MATRIX AS DECLARED IN THE MAIN DIMENSION STATEMENT
        !*  CS  = VECTOR OF THE ZEROES, CS(I), I=1,N
        !*  DZ  = LEGENDRE DERIVATIVES AT THE ZEROES, DZ(I), I=1,N
        !*  DGA = DERIVATIVE MATRIX, DGA(I,J), I=1,N  J=1,N
        !************************************************************************
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N, NM
        REAL(kind=8), DIMENSION(N), INTENT(IN) :: CS, DZ
        REAL(kind=8), DIMENSION(NM,N), INTENT(OUT) :: DGA
        INTEGER :: I, J
        REAL(kind=8) :: VI, ZI, VJ, ZJ

        DGA(1,1) = 0.D0
        IF (N .LE. 1) RETURN

        DO  I=1,N
            VI = DZ(I)
            ZI = CS(I)
            DO J=1,N
                IF (I .NE. J) THEN
                    VJ = DZ(J)
                    ZJ = CS(J)
                    DGA(I,J) = VI/(VJ*(ZI-ZJ))
                ELSE
                    DGA(I,I) = ZI/(1.D0-ZI*ZI)
                ENDIF
            END DO
        END DO
        RETURN
    END SUBROUTINE DMLEGA
END MODULE splib
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
