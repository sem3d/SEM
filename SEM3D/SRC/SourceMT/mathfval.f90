!@
!@ Realized by TARO Mandikizinoyou
!@ This subroutine is used to compute the value for a given function f(x,y,z,t)
!@ for given x, y, z and t values. The variable "t" is related to time variable
!@
MODULE Mathfval

USE parameters, ONLY: rn

IMPLICIT NONE

    TYPE :: FoncValue
      CHARACTER(len=1500) :: valuefx
      CHARACTER(len=1500) :: valuefy
      CHARACTER(len=1500) :: valuefz
      CHARACTER(len=1500) :: valuefxy
      CHARACTER(len=1500) :: valuefyz
      CHARACTER(len=1500) :: valuefxz
      CHARACTER(len=12)   :: var
      CHARACTER           :: source
      INTEGER             :: dim
      REAL(kind=rn),ALLOCATABLE :: fvalue(:)
    END TYPE

    INTERFACE
        SUBROUTINE f2D(f, valx, valy, valt, res)
              USE parameters, ONLY: rn
              USE fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg
              USE Alertes
          IMPLICIT NONE
              CHARACTER(LEN=1500), INTENT(IN)                      :: f
              REAL(rn),            INTENT(IN)                    :: valx, valy, valt
              REAL(rn), INTENT(OUT)                              :: res
              INTEGER,                             PARAMETER     :: nfunc = 1
              INTEGER,                             PARAMETER     :: nvar = 4
              CHARACTER (LEN=*), DIMENSION(nvar),  PARAMETER     :: var  = (/ 'x ', &
                                                                              'y ', &
                                                                              't ', &
                                                                              'pi' /)
              REAL(rn),          DIMENSION(nvar)                 :: val
              CHARACTER(len=1500)                                :: ff
              INTEGER                                            :: i
              REAL(rn)                                           :: pi
              CHARACTER (LEN=2), ALLOCATABLE                     :: varfin(:)
              REAL(rn),          ALLOCATABLE                     :: valfin(:)
              CHARACTER(len=256)                                 :: FunctionName ='f2D'
              CHARACTER(len=256)                                 :: SourceFile = 'mathfval'
              CHARACTER(len=700)                                 :: ErrorSMS
        END SUBROUTINE f2D
   END INTERFACE
        !!
        !!
   INTERFACE
        SUBROUTINE f3D(f, valx, valy, valz, valt, res)
              USE parameters, ONLY: rn
              USE fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg
              USE Alertes
          IMPLICIT NONE
              CHARACTER(LEN=1500), INTENT(IN)                  :: f
              REAL(rn),            INTENT(IN)                  :: valx, valy, valz, valt
              REAL(rn), INTENT(OUT)                            :: res
              INTEGER,                             PARAMETER   :: nfunc = 1
              INTEGER,                             PARAMETER   :: nvar = 5
              CHARACTER (LEN=*), DIMENSION(nvar),  PARAMETER   :: var  = (/ 'x ', &
                                                                            'y ', &
                                                                            'z ', &
                                                                            't ', &
                                                                            'pi' /)
              REAL(rn),          DIMENSION(nvar)               :: val
              CHARACTER(len=1500)                              :: ff
              INTEGER                                          :: i
              REAL(rn)                                         :: pi
              CHARACTER (LEN=2), ALLOCATABLE                   :: varfin(:)
              REAL(rn),          ALLOCATABLE                   :: valfin(:)
              CHARACTER(len=256)                               :: FunctionName ='f3D'
              CHARACTER(len=256)                               :: SourceFile = 'mathfval'
              CHARACTER(len=700)                               :: ErrorSMS
        END SUBROUTINE f3D
   END INTERFACE
        !!
        !!
   INTERFACE
        SUBROUTINE f1D(f, valx, valt, res)
              USE parameters, ONLY: rn
              USE fparser,    ONLY: initf, parsef, evalf, EvalErrType, EvalErrMsg
              USE Alertes
          IMPLICIT NONE
              CHARACTER(LEN=1500), INTENT(IN)                  :: f
              REAL(rn),            INTENT(IN)                  :: valx, valt
              REAL(rn), INTENT(OUT)                            :: res
              INTEGER,                             PARAMETER   :: nfunc = 1
              INTEGER,                             PARAMETER   :: nvar = 3
              REAL(rn),          DIMENSION(nvar)               :: val
              CHARACTER (LEN=*), DIMENSION(nvar),  PARAMETER   :: var  = (/ 'x ', &
                                                                            't ', &
                                                                            'pi'/)
              CHARACTER(len=256)                               :: FunctionName ='f1D'
              CHARACTER(len=256)                               :: SourceFile = 'mathfval'
              CHARACTER(len=1500)                              :: ff
              CHARACTER (LEN=2), ALLOCATABLE                   :: varfin(:)
              REAL(rn),          ALLOCATABLE                   :: valfin(:)
              CHARACTER(len=700)                               :: ErrorSMS
              INTEGER                                          :: i, ierr
              REAL(rn)                                         :: pi
        END SUBROUTINE f1D
    END INTERFACE
!!
!!
END MODULE Mathfval

