! 
! Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!
! This software is distributable under the BSD license. See the terms of the
! BSD license in the documentation provided with this software.
!
MODULE parameters
    use constants, only : fpp
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  ! Specify data types
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  IMPLICIT NONE
  INTEGER, PARAMETER :: is = SELECTED_INT_KIND(1) ! Data type of bytecode
  INTEGER, PARAMETER :: rn = fpp
  !------- -------- --------- --------- --------- --------- --------- --------- -------
  TYPE parametric         
    INTEGER            :: nparam
    CHARACTER(len=2)   :: paramname(1:100)
    REAL(fpp)          :: paramvalue(1:100)
  END TYPE parametric
  TYPE (parametric),PUBLIC:: Addparametricvar
END MODULE parameters


