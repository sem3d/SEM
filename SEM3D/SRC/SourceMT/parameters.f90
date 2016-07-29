! 
! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
MODULE parameters
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  ! Specify data types
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  IMPLICIT NONE
  INTEGER, PARAMETER :: rn = KIND(0.0)            ! Precision of real numbers (modified by M Taro)
  INTEGER, PARAMETER :: is = SELECTED_INT_KIND(1) ! Data type of bytecode
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  TYPE parametric         
    INTEGER             :: nparam
    CHARACTER(len=2)    :: paramname(1:100)
    REAL(rn)            :: paramvalue(1:100)
  END TYPE parametric
  TYPE (parametric),PUBLIC:: Addparametricvar
END MODULE parameters


