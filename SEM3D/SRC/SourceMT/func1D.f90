!@
!@ Realized by TARO Mandikizinoyou
!@ This subroutine is used to compute the value of the function f(x,t)
!@ for given x and t values. The varaible "t" is related to time variable
!@
SUBROUTINE f1D(f, val1, valt, res)

  USE parameters, ONLY: rn, Addparametricvar
  USE fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg
  USE Alertes
  
  IMPLICIT NONE
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  !
  CHARACTER(LEN=1500), INTENT(IN)                  :: f
  REAL(rn), DIMENSION(1),   INTENT(IN)             :: val1, valt 
  REAL(rn), INTENT(OUT)                            :: res  
  
  INTEGER,                             PARAMETER   :: nfunc = 1                                                   
  INTEGER,                             PARAMETER   :: nvar = 3
  REAL(rn),          DIMENSION(nvar)               :: val 
  CHARACTER (LEN=*), DIMENSION(nvar),  PARAMETER   :: var  = (/ 'x ', &
                                                                't ', &
                                                                'pi'/)
  CHARACTER (LEN=2), ALLOCATABLE                   :: varfin(:)
  CHARACTER(len=1500)                              :: ff
  REAL(rn),          ALLOCATABLE                   :: valfin(:)
  CHARACTER(len=256)                               :: FunctionName ='f1D'
  CHARACTER(len=256)                               :: SourceFile = 'func1D'
  CHARACTER(len=700)                               :: ErrorSMS
  INTEGER                                          :: i, ierr
  REAL(rn)                                         :: pi
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  !
  
  CALL initf (nfunc)
                                   
  pi =  4.0*ATAN(1.0)
  val = (/Val1, valt, pi/)
  ALLOCATE (varfin(1:nvar+Addparametricvar%nparam))
  ALLOCATE (valfin(1:nvar+Addparametricvar%nparam))
  varfin(1:nvar)=var
  valfin(1:nvar)=val
  IF (Addparametricvar%nparam.ne.0) THEN
      DO i=1, Addparametricvar%nparam
          Varfin(nvar+i) = trim(Addparametricvar%paramname(i))
          valfin(nvar+i) =Addparametricvar%paramvalue(i)
      ENDDO
  ENDIF

  DO i=1,nfunc
     CALL stringstring(f,ff)
     CALL parsef(i,ff(1:LEN_TRIM(ff)), varfin)         
  END DO
  
  DO i=1,nfunc
     res = evalf(i, valfin)                                                     
     IF (EvalErrType > 0) THEN 
        !WRITE(*,*)'*** Error: ',EvalErrMsg()
        ErrorSMS = ' '//EvalErrMsg ()
        CALL ErrorMessage(ErrorSMS,FunctionName,SourceFile)
     ENDIF
  END DO
  DEALLOCATE(varfin,valfin)

END SUBROUTINE f1D
