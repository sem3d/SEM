!@
!@
!@
!@
FUNCTION CheckDimension(Variables)

    USE fparser,    ONLY: LowCase
    USE Alertes

    IMPLICIT NONE
    !--------------------------------------------------------------------------------
    !
    INTEGER                        :: CheckDimension
    CHARACTER(len=12), INTENT(in)  :: Variables
    CHARACTER(len=4)               :: RealChar
    INTEGER                        :: i, Length, ExistX, ExistY, ExistZ, ExistT
    CHARACTER(len=256)             :: FunctionName ='CheckDimension'
    CHARACTER(len=256)             :: SourceFile='Dimension'
    CHARACTER(len=700)             :: ErrorSMS
    !--------------------------------------------------------------------------------
    !
    Length = LEN_TRIM(Variables)
    RealChar = TRIM(ADJUSTL(Variables))
    CALL LowCase(RealChar, RealChar)
    ExistX=0
    ExistY=0
    ExistZ=0
    ExistT=0

    DO i=1,Length
        SELECT CASE(RealChar(i:i))
        CASE('x')
            ExistX = 1
        CASE('y')
            ExistY = 1
        CASE('z')
            ExistZ = 1
        CASE('t')
            ExistT = 1
        END SELECT
    END DO

    IF (((ExistX+ExistY+ExistZ).eq.1).or.(((ExistX+ExistY+ExistZ).eq.0).and.(ExistT.eq.1))) THEN
        CheckDimension = 1
    ELSEIF ((ExistX+ExistY+ExistZ).eq.2) THEN
        CheckDimension = 2
    ELSEIF ((ExistX+ExistY+ExistZ).eq.3) THEN
        CheckDimension = 3
    ELSE
        ErrorSMS = 'The variables of Analytical function are not defined.'
        CALL ErrorMessage(ErrorSMS,FunctionName,SourceFile)
    ENDIF

END FUNCTION CheckDimension

