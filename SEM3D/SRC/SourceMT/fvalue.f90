!@
!@ Realized by TARO Mandikizinoyou
!@ This subroutine is used to compute the value for a given function f(x,t)
!@ for given x and t values. The variable "t" is related to time variable
!@

SUBROUTINE ffvalue(func,XYZCoord,Valt)

USE Mathfval
USE Alertes
USE ssources
USE parameters, ONLY: rn

IMPLICIT NONE
!--------------------------------------------------------------------------------
!
TYPE(FoncValue), INTENT(inout)    :: func
REAL(rn), INTENT(in)              :: XYZCoord(3),Valt
INTEGER                           :: CheckDimension
CHARACTER(len=256)                :: FunctionName ='fvalue'
CHARACTER(len=256)                :: SourceFile = 'fvalue'
CHARACTER(len=700)                :: ErrorSMS
!--------------------------------------------------------------------------------
!

     IF (func%dim.eq.1) THEN
           SELECT CASE (CheckDimension(func%var))
               CASE (1)
                   CALL f1D(func%valuefx, XYZCoord(1), valt, func%fvalue(1))
               CASE (2)
                   CALL f2D(func%valuefx, XYZCoord(1), XYZCoord(2), valt, func%fvalue(1))
               CASE (3)
                   CALL f3D(func%valuefx, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(1))
               CASE DEFAULT
                   ErrorSMS = 'Unknown problem dimension'
                   CALL ErrorMessage(ErrorSMS,FunctionName,SourceFile)
           END SELECT
     ELSEIF (func%dim.eq.2) THEN
           IF (func%source.eq.'F') THEN
               SELECT CASE (CheckDimension(func%var))
                   CASE (1)
                         CALL f1D(func%valuefx, XYZCoord(1), valt, func%fvalue(1))
                         CALL f1D(func%valuefy, XYZCoord(1), valt, func%fvalue(2))
                   CASE (2)
                         CALL f2D(func%valuefx, XYZCoord(1), XYZCoord(2), valt, func%fvalue(1))
                         CALL f2D(func%valuefy, XYZCoord(1), XYZCoord(2), valt, func%fvalue(2))
                   CASE (3)
                         CALL f3D(func%valuefx, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(1))
                         CALL f3D(func%valuefy, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(2))
                   CASE DEFAULT
                         ErrorSMS = 'Unknown problem dimension'
                         CALL ErrorMessage(ErrorSMS,FunctionName,SourceFile)
              END SELECT
           ELSEIF (func%source.eq.'M') THEN
              SELECT CASE (CheckDimension(func%var))
                   CASE (1)
                          CALL f1D(func%valuefx, XYZCoord(1), valt, func%fvalue(1))
                          CALL f1D(func%valuefy, XYZCoord(1), valt, func%fvalue(2))
                          CALL f1D(func%valuefxy, XYZCoord(1), valt, func%fvalue(3))
                   CASE (2)
                         CALL f2D(func%valuefx, XYZCoord(1), XYZCoord(2), valt, func%fvalue(1))
                         CALL f2D(func%valuefy, XYZCoord(1), XYZCoord(2), valt, func%fvalue(2))
                         CALL f2D(func%valuefxy, XYZCoord(1), XYZCoord(2), valt, func%fvalue(3))
                   CASE (3)
                         CALL f3D(func%valuefx, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(1))
                         CALL f3D(func%valuefy, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(2))
                         CALL f3D(func%valuefxy, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(3))
                   CASE DEFAULT
                         ErrorSMS = 'Unknown problem dimension'
                         CALL ErrorMessage(ErrorSMS,FunctionName,SourceFile)
              END SELECT
           END IF
      ELSEIF (func%dim.eq.3) THEN
           IF (func%source.eq.'F') THEN
               SELECT CASE (CheckDimension(func%var))
                   CASE (1)
                         CALL f1D(func%valuefx, XYZCoord(1), valt, func%fvalue(1))
                         CALL f1D(func%valuefy, XYZCoord(1), valt, func%fvalue(2))
                         CALL f1D(func%valuefz, XYZCoord(1), valt, func%fvalue(3))
                   CASE (2)
                         CALL f2D(func%valuefx, XYZCoord(1), XYZCoord(2), valt, func%fvalue(1))
                         CALL f2D(func%valuefy, XYZCoord(1), XYZCoord(2), valt, func%fvalue(2))
                         CALL f2D(func%valuefz, XYZCoord(1), XYZCoord(2), valt, func%fvalue(3))
                   CASE (3)
                         CALL f3D(func%valuefx, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(1))
                         CALL f3D(func%valuefy, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(2))
                         CALL f3D(func%valuefz, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(3))
                   CASE DEFAULT
                         ErrorSMS = 'Unknown problem dimension'
                         CALL ErrorMessage(ErrorSMS,FunctionName,SourceFile)
                END SELECT
           ELSEIF (func%source.eq.'M') THEN
          
           	SELECT CASE (CheckDimension(func%var))
                   CASE (1)
                         CALL f1D(func%valuefx, XYZCoord(1), valt, func%fvalue(1))
                         CALL f1D(func%valuefy, XYZCoord(1), valt, func%fvalue(2))
                         CALL f1D(func%valuefz, XYZCoord(1), valt, func%fvalue(3))
                         CALL f1D(func%valuefxy, XYZCoord(1), valt, func%fvalue(4))
                         CALL f1D(func%valuefyz, XYZCoord(1), valt, func%fvalue(5))
                         CALL f1D(func%valuefxz, XYZCoord(1), valt, func%fvalue(6))
	          CASE (2)
                         CALL f2D(func%valuefx, XYZCoord(1), XYZCoord(2), valt, func%fvalue(1))
                         CALL f2D(func%valuefy, XYZCoord(1), XYZCoord(2), valt, func%fvalue(2))
                         CALL f2D(func%valuefz, XYZCoord(1), XYZCoord(2), valt, func%fvalue(3))
                         CALL f2D(func%valuefxy, XYZCoord(1),XYZCoord(2), valt, func%fvalue(4))
                         CALL f2D(func%valuefyz, XYZCoord(1), XYZCoord(2), valt, func%fvalue(5))
                         CALL f2D(func%valuefxz, XYZCoord(1), XYZCoord(2), valt, func%fvalue(6))
                  CASE (3)
                         CALL f3D(func%valuefx, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(1))
                         CALL f3D(func%valuefy, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(2))
                         CALL f3D(func%valuefz, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(3))
                         CALL f3D(func%valuefxy, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(4))
                         CALL f3D(func%valuefyz, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(5))
                         CALL f3D(func%valuefxz, XYZCoord(1), XYZCoord(2), XYZCoord(3), valt, func%fvalue(6))
                  CASE DEFAULT
                         ErrorSMS = 'Unknown problem dimension'
                         CALL ErrorMessage(ErrorSMS,FunctionName,SourceFile)
                  END SELECT
           END IF
     END IF
     
ENDSUBROUTINE ffvalue
!!
!!
SUBROUTINE stringstring(func, string)

 implicit none
 character(len=1500), intent(in)    :: func
 character(len=1500), intent(out)   :: String
 integer                            :: l,n,k,j,ns
 logical                            :: EndOfLine
 

 do n=1,len_trim(func)
    string(n:n) = ' '
 end do
 n=0; k=1; l = len_trim(func)
 EndOfLine = l-k < 0
 do while (.not.EndOfLine)
   j = index(func(k:l),"/")
   if (j == 0) then
       j=l+1
   else
       j=j+k-1
       string(j:j) = "/"
   end if
 
   if ((j.ne.k).and.(len_trim(func(k:j-1)).ne.0)) then
       read(func(k:j-1),*) string(k:j-1)
   endif
   k=j+1
   EndOfLine=l-k < 0
 enddo

END SUBROUTINE  stringstring



