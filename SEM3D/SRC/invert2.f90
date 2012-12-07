subroutine invert2 (A, det)

!-------------------------------------------------------------------------------
! Subroutine to compute a 2x2 inverse matrix
! using the definition of inverse
! 
! Modified Gaetano Festa 01/06/2005
! -------------------------------------------------------------------------------
implicit none

real, dimension (0:1,0:1), intent (INOUT) :: A
real, intent (OUT) :: det

real :: aus

! compute determinant of A using Sarrus law

det = A(0,0)*A(1,1)
det = det - A(0,1)*A(1,0)

aus = A (0,0)
A (0,0) = A(1,1)
A (1,0) = -A(1,0)
A (0,1) = -A(0,1)
A (1,1) = aus

! Compute the inverse 

A = A /det

return
end subroutine invert2
