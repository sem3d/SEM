module ssources


type :: Source

integer :: i_dir, i_type_source, i_time_function, elem, proc
integer , dimension(0:2) :: gll
real :: Xsource,YSource,Zsource, tau_b,cutoff_freq

end type 

contains

! ###############################################
real function CompSource (Sour,time,np)

type (source) :: Sour
integer ::np
real :: time

select case (Sour%i_type_source)
case (1)   ! Impulse
        if (np /= Sour%i_dir) then
           CompSource = 0. 
        else   
        select case (Sour%i_time_function)
        case (1) 
                CompSource = Gaussian (time,Sour%tau_b)
        case (2)
                CompSource = Ricker (time,Sour%tau_b,Sour%cutoff_freq)
         end select
         end if
case (2) ! Explosion
        select case (Sour%i_time_function)
	case (1)
		CompSource = Gaussian (time,Sour%tau_b)
	case (2)
		CompSource = Ricker (time,Sour%tau_b,Sour%cutoff_freq)

        end select
end select
return
end function

! ################################################
real function Gaussian (time, tau)

real :: tau,time

Gaussian = -(time-tau) * exp (-(time-tau)**2/tau**2)

return
end function
 
! ################################################
real function Ricker (time,tau,f0)

use pig

real :: time, tau, f0
real :: sigma

sigma = pi * f0 * (time-tau)
sigma = sigma**2
Ricker = (1 - 2*sigma) * exp(-sigma)

return
end function

! #################################################
end module ssources
