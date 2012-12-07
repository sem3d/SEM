subroutine invert_3d (a,Det)

! Modified by Gaetano Festa 02/02/2005


implicit none

real, dimension (0:2,0:2), intent (INOUT) :: a
real, intent (OUT) :: Det 

real, dimension(0:2,0:2) :: Inverse_A


Det = a(0,0)* (a(1,1)*a(2,2) - a(2,1)*a(1,2) ) + a(1,0) * (a(0,2) * a(2,1) - a(0,1) * a(2,2) ) + a(2,0)  * (a(0,1)*a(1,2) - a(0,2)*a(1,1) ) 

Inverse_A (0,0) = a(1,1)*a(2,2) - a(2,1)*a(1,2)
Inverse_A (1,0) = a(0,2) * a(2,1) - a(0,1) * a(2,2) 
Inverse_A (2,0) = a(0,1)*a(1,2) - a(0,2)*a(1,1)

Inverse_A (0,1) = a(1,2)*a(2,0) - a(1,0)*a(2,2)
Inverse_A (1,1) = a(0,0) * a(2,2) - a(0,2) * a(2,0) 
Inverse_A (2,1) = a(0,2)*a(1,0) - a(0,0)*a(1,2)

Inverse_A (0,2) = a(1,0)*a(2,1) - a(1,1)*a(2,0)
Inverse_A (1,2) = a(0,1) * a(2,0) - a(2,1) * a(0,0) 
Inverse_A (2,2) = a(1,1)*a(0,0) - a(0,1)*a(1,0)

A = TRANSPOSE (Inverse_A)
A = A/ Det 

return
end subroutine invert_3d
