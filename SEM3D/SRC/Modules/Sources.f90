!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Sources.f90
!!\brief Assure la gestion des sources.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssources

    use constants
    implicit none

    type :: Source
       ! GENERAL PARAMETERS (see SourcePosition.f90-Source.f90)
       integer                     :: proc                      ! source belonging processor
       integer                     :: elem                      ! source belonging elem
       integer                     :: i_type_source             ! source type (solid pulse-moment-fluid pulse)
       ! SPATIAL PARAMETERS
       real(fpp), dimension(0:2)   :: refcoord                  ! local coordinates (master element)
       real(fpp), dimension(0:2)   :: dir                       ! source direction
       real(fpp)                   :: Xsource, YSource, Zsource ! source coordinates
       real(fpp)                   :: amplitude_factor          ! amplitude factor
       ! TIME PARAMETERS
       real(fpp)                   :: ts                        ! time shift
       integer                     :: i_time_function           ! source type (ricker, gabor, etc.)
       ! MOMENT SOURCE
       real(fpp), dimension(0:2,0:2)             :: InvGrad     ! Inverse Jacobian
       real(fpp), dimension(0:2,0:2)             :: Moment      ! Moment tensor
       real(fpp), dimension (:,:,:,:), pointer   :: coeff       ! weight coefficient
       ! SOURCE FROM EXTERNAL FILE
       character(len = 30)                  :: time_file        ! file name of external source

       real(fpp) :: tau_b,cutoff_freq,Q,X,Y,L,v,d,a
       real(fpp) :: radius,realcolat,reallong,refcolat,reflong
       !   ajout de parametres pour definir Gabor signal source de type 4
       !   ajout de gamma et ts
       real(fpp) ::  gamma

       real(fpp), dimension(0:3) :: fh

       real(fpp), dimension(:), pointer :: timefunc
       real(fpp), dimension(:), pointer :: ampli, time
       real(fpp), allocatable, dimension(:,:,:,:) :: ExtForce
    end type Source

contains

    !>
    !! \brief
    !!
    !! \param type (source) Sour
    !! \param real time
    !<
    real(fpp) function CompSource (Sour,time,ntime)
        implicit none
        type (source), intent(in) :: Sour
        integer, intent(in) :: ntime
        real(fpp), intent(in) :: time

        CompSource = 0d0

        select case (Sour%i_time_function)
        case (1)
            CompSource = Gaussian (time, Sour%ts, Sour%tau_b)
        case (2)
            CompSource = Ricker (time,Sour%tau_b,Sour%cutoff_freq)
        case (3)
            CompSource = Sour%timefunc(ntime)
        case (4)
            !   developpement de la source du benchmark can1
            CompSource = Gabor (time,Sour%tau_b,Sour%cutoff_freq,Sour%gamma,Sour%ts)
        case (5)
            !   modif pour benchmark can2
            CompSource = Source_File (time,Sour)
            !   modif pour benchmark can2
        case (6)
            ! Source benchmark spice M0*(1-(1+(t/T)**gamma)exp(-(t/T)**gamma)) avec T=1/freq
            CompSource = Source_Spice_Bench(time, Sour)
        case (7)
            ! Sinus, pour test. param : tau, cutoff_freq
            CompSource = Source_sinewave(time, Sour)
        case (8)
            ! Square. Param : tau, ts, gamma
            CompSource = Source_square(time, Sour)
        case (9)
            ! Square. Param : ts, gamma
            CompSource = Source_tanh(time, Sour)
        case (10)
            ! Square. Param : ts, gamma
            CompSource = Ricker_fl(time, Sour%tau_b, Sour%cutoff_freq)
        case (11)
            ! fonction de triangle
            CompSource = Triangle(time, Sour%tau_b)
        case (12)
            ! Heaviside Step Function (Kausel-2006)
            CompSource = HSF(time, Sour%tau_b)
        case(13)
            !Double-M wavelet (Al Shaer et al 2008)
            CompSource = DM(time, Sour%tau_b,Sour%Q,Sour%X,Sour%Y,Sour%L,Sour%v,Sour%d,Sour%a)

        end select
        CompSource = CompSource*Sour%amplitude_factor
        return
    end function CompSource


    real(fpp) function Source_Spice_Bench(time, Sour)
        implicit none
        type(source), intent(in) :: Sour
        real(fpp), intent(in) :: time
        !
        real(fpp) :: T, k, s

        if (time<Sour%ts) then
            Source_Spice_Bench = 0d0
            return
        end if

        T = 1./Sour%cutoff_freq
        k = Sour%gamma
        if (k<1d0) k=1d0

        s = ((time-Sour%ts)/T)**k

        Source_Spice_Bench = (1-(1+s)*exp(-s))
        return
    end function Source_Spice_Bench

    real(fpp) function Source_tanh(time, Sour)
        implicit none
        type(source), intent(in) :: Sour
        real(fpp), intent(in) :: time
        !
        real(fpp) :: k,t0

        t0 = Sour%ts
        k = Sour%gamma

        Source_tanh = 0.5d0*(tanh(k*(time-t0))+1d0)
        return
    end function Source_tanh

    real(fpp) function Source_square(t, Sour)
        implicit none
        ! A smoothed square
        type(source), intent(in) :: Sour
        real(fpp), intent(in) :: t
        !
        real(fpp) :: dt,t0,k,w0,winf

        dt = Sour%tau_b
        t0 = Sour%ts
        k = Sour%gamma

        ! Primitive : (log(cosh(k*(t-t0)))-log(cosh(k*(t0+dt-t))))/k
        w0 = (log(cosh(k*(-t0)))-log(cosh(k*(t0+dt))))/k
        winf = dt

        Source_square = (tanh(k*(t-t0))+tanh(k*(t0+dt-t)))/(winf-w0)
        return
    end function Source_square

    real(fpp) function Source_sinewave(time, Sour)
        implicit none
        type(source), intent(in) :: Sour
        real(fpp), intent(in) :: time
        !
        real(fpp) :: f0, t0

        f0 = Sour%cutoff_freq
        t0 = Sour%ts

        Source_sinewave = sin(2*M_PI*f0*(time-t0))
        return
    end function Source_sinewave


    real(fpp) function Source_File(tt,Sour)
        implicit none
        type(source), intent(in)  :: Sour
        real(fpp), intent(in)          :: tt
        integer :: iflag, i

        if(tt < Sour%time(0) .or. tt > Sour%time(size(Sour%time)-1))then
            Source_File = 0d0
            return
        end if
        i = 0
        iflag = 0

        do  while( iflag == 0 )
            if(tt >= Sour%time(i) .and. tt < Sour%time(i+1) )then
                Source_File = Sour%ampli(i) +            &
                    (tt-Sour%time(i))*(Sour%ampli(i+1)-Sour%ampli(i))/(Sour%time(i+1)-Sour%time(i))
                iflag = 1
                return
            else
                i = i + 1
            end if
        end do

    end function Source_File

    subroutine read_source_file(Sour)
        implicit none
        !- lecture directe d'un fichier temps-amplitude pour la source
        type(Source), intent(inout)   :: Sour
        integer                       :: nb_time_step
        integer                       :: i
        real(fpp)                          :: tr, trr

        i = 0 ; nb_time_step = 0
        ! count

        open(10,file=Sour%time_file,action="read",status="old")
        do
            read(10,*,end=100) tr, trr
            i = i + 1
        end do
100     close(10)
        ! nombre de donnees en entree pour la source
        nb_time_step = i

        allocate(Sour%time(0:nb_time_step-1),Sour%ampli(0:nb_time_step-1))

        open(10,file=Sour%time_file,action="read",status="old")
        do i=0, nb_time_step-1
            read(10,*,end=101) Sour%time(i),Sour%ampli(i)
        end do
101     close(10)

    end subroutine read_source_file
    !   modif pour benchmark can2
    !-------------------------------------------------

    real(fpp) function Gaussian (time,ts,tau)
        implicit none
        real(fpp), intent(in) :: tau, time, ts

        if ( (time-ts) < 8*tau ) then
            Gaussian = -2*(time-ts) * exp (-(time-ts)**2/tau**2)
        else
            Gaussian = 0.
        endif

        return
    end function Gaussian

    ! ################################################
    !>
    !! \fn function Ricker (time,tau,f0)
    !! \brief
    !!
    !! \param real time
    !! \param real tau
    !! \param real f0
    !<




    real(fpp) function Triangle (time,tau)
        implicit none
        !

         real(fpp), intent(in) :: time, tau

         !! tau = coefficient of pression

         if ( time < 0.005 ) then
              Triangle = - time * tau * ( 1e10 )
         elseif ( time < 0.01 ) then
              Triangle = - ( -time + 0.01 ) * tau * ( 1e10 )
         else
              Triangle = 0.
         endif

         return
    end function Triangle

   ! ###############################################################################

    real(fpp) function HSF (time, tau)
        implicit none
        ! HEAVISIDE STEP FUNCTION (KAUSEL-2006)
        real(fpp), intent(in) :: time, tau

        if ( time < 0.00000001 ) then
            HSF =0.
        endif
        if ( time == 0. ) then
            HSF = -0.5 * tau
        endif
        if ( time > 0.00000001 ) then
            HSF = -1 * tau
        endif
        return
    end function HSF



    real(fpp) function DM (time, tau,Q,X,Y,L,v,d,a)
        implicit none
        ! DoubleM (Al Shaer et al 2008)
        real(fpp), intent(in) :: time, tau,Q,X,Y,L,v,d,a

        DM = Q*Y/2*((X)**(((v*(time-tau)-a)**2/d**2))+(X)**(((v*(time-tau)-a-L)**2/d**2)))

        return
    end function DM

    ! ############################################################################

    real(fpp) function Ricker (time,tau,f0)
        implicit none
        !
        real(fpp), intent(in) :: time, tau, f0
        real(fpp) :: sigma, alpha

        alpha = -1d0*M_PI**2*f0**2
        if ( time < 2.5*tau ) then
            sigma = alpha * (time-tau)**2
            Ricker = 2d0*alpha*(1 + 2*sigma) * dexp(sigma)
        else
            Ricker = 0.
        endif

        return
    end function Ricker

    ! #################################################

    real(fpp) function Gabor (time,tau,fp,gamma,ts)
        implicit none
        !
        real(fpp), intent(in) :: time, tau, fp, gamma, ts
        !
        real(fpp) :: sigma
        real(fpp) ::  xomega,  xval1, xval2
        xomega  = M_PI*0.5

        if ( time < 32. ) then
            sigma = 2. * M_PI * fp * (time-ts)
            xval1 = cos(sigma + xomega)
            sigma = sigma/gamma
            sigma = sigma**2
            if ( sigma < 100. ) then
                xval2 = exp(-sigma)
            else
                xval2 = 0.
            endif
            Gabor = xval2*xval1*tau
        else
            Gabor = 0.
        endif

        !       print*,' Gabor ',time, Gabor
        return
    end function Gabor

    !----------------------------------------------------
    !----------------------------------------------------
    real(fpp) function Ricker_Fl(time,tau,f0)
        implicit none
        ! Ricker function for pressure wave: the same as the previous one, but with
        !    one time derivative to be coherent
        real(fpp),intent(in)  :: time, tau, f0
        !
        real(fpp)  :: sigma,sigma2,sigma3

        sigma = M_PI*f0*(time-tau)
        sigma2 = sigma**2
        sigma3 = M_PI*f0

        Ricker_Fl = -4d0*sigma*sigma3*exp(-sigma2)-2d0*sigma*sigma3*Ricker(time,tau,f0)

        return
    end function Ricker_Fl

end module ssources

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
