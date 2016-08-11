!@
!@ Realized by TARO Mandikizinoyou
!@ This subroutine is used to compute the value of the function f(x,y,z,t)
!@ for given x, y, z and t values. The varaible "t" is related to time variable
!@
SUBROUTINE f3D(f, valx, valy, valz, valt, res)

  USE parameters, ONLY: rn, Addparametricvar
  USE fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg
  USE Alertes
  
  IMPLICIT NONE
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  !
  CHARACTER(LEN=1500), INTENT(IN)                  :: f
  REAL(rn), DIMENSION(1),   INTENT(IN)             :: valx, valy, valz, valt  
  REAL(rn), INTENT(OUT)                            :: res 
  
  INTEGER,                             PARAMETER   :: nfunc = 1
  INTEGER,                             PARAMETER   :: nvar = 5
  CHARACTER (LEN=*), DIMENSION(nvar),  PARAMETER   :: var  = (/ 'x ', &
                                                                'y ', &
                                                                'z ', &
                                                                't ', &
                                                                'pi' /)
  REAL(rn),          DIMENSION(nvar)               :: val 
  CHARACTER (LEN=2), ALLOCATABLE                   :: varfin(:)
  REAL(rn),          ALLOCATABLE                   :: valfin(:)
  CHARACTER(len=1500)                              :: ff
  INTEGER                                          :: i
  REAL(rn)                                         :: pi
  CHARACTER(len=256)                               :: FunctionName ='f3D'
  CHARACTER(len=256)                               :: SourceFile = 'func3D'
  CHARACTER(len=700)                               :: ErrorSMS
 !--------- -------- --------- --------- --------- --------- --------- --------- -----
 !
  CALL initf (nfunc)
  
  pi  = 4.0*ATAN(1.0)
  val = (/valx, &
          valy, &
          valz, &
          valt, &
          pi/)
  ALLOCATE (varfin(1:nvar+Addparametricvar%nparam))
  ALLOCATE (valfin(1:nvar+Addparametricvar%nparam))
  varfin(1:nvar)=var
  valfin(1:nvar)=val
  
  IF (Addparametricvar%nparam.ne.0) THEN
      DO i=1, Addparametricvar%nparam
         varfin(nvar+i) =trim(Addparametricvar%paramname(i))
         valfin(nvar+i) =Addparametricvar%paramvalue(i)
      ENDDO
  ENDIF
  
  DO i=1,nfunc
     CALL stringstring(f,ff)
     CALL parsef (i, ff(1:LEN_TRIM(ff)), varfin)
  END DO 
  
  DO i=1,nfunc
     res = evalf (i, valfin)  
     IF (EvalErrType > 0) THEN
        ErrorSMS = ' '//EvalErrMsg ()
        CALL ErrorMessage(ErrorSMS,FunctionName,SourceFile)
     ENDIF
  END DO
  DEALLOCATE(varfin,valfin)
 !@
 !@
END SUBROUTINE f3D


