module ssources

    type :: Source
       integer :: i_dir, i_type_source, i_time_function, elem, proc
       integer , dimension(0:2) :: gll
       real :: Xsource,YSource,Zsource,tau_b,cutoff_freq
       real, dimension(0:2)  :: RefCoord
       real, dimension(0:2,0:2)  :: Moment,InvGrad
       real, allocatable, dimension(:,:,:) :: ExtForce
       real, allocatable, dimension(:,:,:,:) :: MomForce
    end type Source

contains

    !-------------------------------------------------
    real function CompSource(Sour,time,np)

        type(source) :: Sour
        integer ::np
        real :: time

        select case(Sour%i_type_source)
        case(1)   ! Impulse
            if(np /= Sour%i_dir)then
                CompSource = 0.
            else
                select case(Sour%i_time_function)
                case(1)
                    CompSource = Gaussian(time,Sour%tau_b)
                case(2)
                    CompSource = Ricker(time,Sour%tau_b,Sour%cutoff_freq)
                end select
            end if
        case(2) ! Moment tensor source
            select case(Sour%i_time_function)
            case(1)
                CompSource = Gaussian(time,Sour%tau_b)
            case(2)
                CompSource = Ricker(time,Sour%tau_b,Sour%cutoff_freq)
            end select
        end select

        return

    end function CompSource
    !-----------------------------------------------------
    !-----------------------------------------------------
    real function Gaussian(time,tau)

        real :: tau,time

        Gaussian = -(time-tau) * exp (-(time-tau)**2/tau**2)

        return

    end function Gaussian
    !------------------------------------------------------
    !------------------------------------------------------
    real function Ricker(time,tau,f0)

        use pig
        real :: time, tau, f0
        real :: sigma

        sigma = pi * f0 * (time-tau)
        sigma = sigma**2
        Ricker = (1 - 2*sigma) * exp(-sigma)

        return

    end function Ricker
    !----------------------------------------------------
    !----------------------------------------------------
    real function Ricker_Fl(time,tau,f0)
        ! Ricker function for pressure wave: the same as the previous one, but with
        !    one time derivative to be coherent
        real  :: time, tau, f0
        real  :: sigma,pi,sigma2,sigma3

        pi = acos(-1d0)
        sigma = pi*f0*(time-tau)
        sigma2 = sigma**2
        sigma3 = pi*f0

        Ricker_Fl = -4d0*sigma*sigma3*exp(-sigma2)-2d0*sigma*sigma3*Ricker(time,tau,f0)

        return
    end function Ricker_Fl
    !----------------------------------------------------
    !----------------------------------------------------
    real function CompSource_Fl(Sour,time)
        ! only a Ricker for the time being.
        type(source) :: Sour
        real :: time

        CompSource_Fl = Ricker_fl(time,Sour%tau_b,Sour%cutoff_freq)

        return
    end function CompSource_Fl
    !----------------------------------------------------
    !----------------------------------------------------
end module ssources
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
