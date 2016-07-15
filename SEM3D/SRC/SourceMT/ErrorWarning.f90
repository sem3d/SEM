!@
!@ Realized by TARO Mandikizinoyou
!@ This subroutine is used to display error message and then stop the computation
!@
MODULE Alertes

CONTAINS
!--------------------------------------------------------------------------------
 
    SUBROUTINE ErrorMessage(ErrorSMS,FunctionName,SourceFile)

    IMPLICIT NONE
    !--------------------------------------------------------------------------------
    !
    CHARACTER(len=700), INTENT(in) :: ErrorSMS
    CHARACTER(len=256), INTENT(in) :: FunctionName
    CHARACTER(len=256), INTENT(in) :: SourceFile
    !--------------------------------------------------------------------------------
    !
    WRITE(*,*) ' '
    WRITE(*,*) '__________________________________________________________'
    WRITE(*,*) ' -> Error  : '//TRIM(ErrorSMS)
    WRITE(*,*) '    In     : '//TRIM(FunctionName)//' function or subroutine'
    WRITE(*,*) '    Source : '//TRIM(SourceFile)//'.f90 file'
    WRITE(*,*) '    '
    WRITE(*,*) '    Computation will be stopped '
    STOP

    END SUBROUTINE ErrorMessage
    !!
    !!
    !!
    !!
    SUBROUTINE WarningMessage(WarningSMS,FunctionName,SourceFile)

    IMPLICIT NONE
    !--------------------------------------------------------------------------------
    !
    CHARACTER(len=700), INTENT(in) :: WarningSMS
    CHARACTER(len=256), INTENT(in) :: FunctionName
    CHARACTER(len=256), INTENT(in) :: SourceFile
    !--------------------------------------------------------------------------------
    !
    WRITE(*,*) ' '
    WRITE(*,*) '__________________________________________________________'
    WRITE(*,*) ' -> Warning  : '//TRIM(WarningSMS)
    WRITE(*,*) '    In       : '//TRIM(FunctionName)//' function or subroutine'
    WRITE(*,*) '    Source   : '//TRIM(SourceFile)//'.f90 file'
    WRITE(*,*) '    '
    WRITE(*,*) '    Computation will be continued '
    WRITE(*,*) '    '
    
    END SUBROUTINE WarningMessage
    
END MODULE Alertes
