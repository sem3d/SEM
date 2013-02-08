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
       ! xxx
       character(len = 30)  :: time_file
       integer :: i_dir, i_type_source, i_time_function, elem, proc
       integer, dimension(0:2) :: gll
       real :: tau_b,cutoff_freq
       real :: radius,realcolat,reallong,refcolat,reflong

       !   ajout de parametres pour definir Gabor signal source de type 4
       !   ajout de gamma et ts
       real ::  gamma, ts

       real :: Xsource,YSource,Zsource

       real, dimension(0:3) :: fh
       real, dimension(0:2) :: refcoord
       real, dimension(:), pointer :: timefunc
       real, dimension(0:2,0:2) :: Moment, InvGrad
       real, dimension (:,:,:,:), pointer :: coeff
       real, dimension(:), pointer :: ampli, time
       real, allocatable, dimension(:,:,:) :: ExtForce
       real, allocatable, dimension(:,:,:,:) :: MomForce
    end type Source

contains

    !>
    !! \brief
    !!
    !! \param type (source) Sour
    !! \param real time
    !<
    real function CompSource (Sour,time,ntime)
        implicit none
        type (source), intent(in) :: Sour
        integer, intent(in) :: ntime
        real, intent(in) :: time

        select case (Sour%i_time_function)
        case (1)
            CompSource = Gaussian (time,Sour%tau_b)
        case (2)
            CompSource = Ricker (time,Sour%tau_b,Sour%cutoff_freq)
        case (3)
            CompSource = Sour%timefunc(ntime)
        case (4)
            !   developpement de la source du benchmark can1
            CompSource = Gabor (time,Sour%tau_b,Sour%cutoff_freq,Sour%gamma,Sour%ts)
        case (5)
            !   modif pour benchmark can2
            CompSource = Source_File (time,Sour%tau_b,Sour)
            !   modif pour benchmark can2
        case (6)
            ! Source benchmark spice M0*(1-(1+t/T)exp(-t/T)) avec T=1/freq
            CompSource = Source_Spice_Bench(time, Sour)
        end select

        return
    end function CompSource


    real function Source_Spice_Bench(time, Sour)
        implicit none
        ! only a Ricker for the time being.
        type(source), intent(in) :: Sour
        real, intent(in) :: time
        !
        real :: T

        T = 1./Sour%cutoff_freq

        Source_Spice_Bench = (1-(1+time/T)*exp(-time/T))
        return
    end function Source_Spice_Bench

    !>
    !! \fn function Gaussian (time, tau)
    !! \brief
    !!
    !! \param real time
    !! \param real tau
    !<
    !-------------------------------------------------
    !   modif pour benchmark can2
    real function Source_File(tt,tau,Sour)
        implicit none
        type(source), intent(in)  :: Sour
        real, intent(in)          :: tt
        real :: tau
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
                Source_File = Source_File * tau
                iflag = 1
                !print*,'source5 fichier ',Source_File,tau,' temps ',tt
                !      print*,' interval ',i,i+1
                !      print*,' '
                return
            else
                !            print*,' iflag ',iflag ,i
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
        real                          :: tr, trr

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

    real function Gaussian (time,tau)
        implicit none
        real, intent(in) :: tau,time

        if ( time < 2.5*tau ) then
            Gaussian = -(time-tau) * exp (-(time-tau)**2/tau**2)
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
    real function Ricker (time,tau,f0)
        implicit none
        !
        real, intent(in) :: time, tau, f0
        real :: sigma

        if ( time < 2.5*tau ) then
            sigma = M_PI * f0 * (time-tau)
            sigma = sigma**2
            Ricker = (1 - 2*sigma) * exp(-sigma) !version Gsa Ipsis (amplitude)
        else
            Ricker = 0.
        endif

        return
    end function Ricker

    ! #################################################

    real function Gabor (time,tau,fp,gamma,ts)
        implicit none
        !
        real, intent(in) :: time, tau, fp, gamma, ts
        !
        real :: sigma
        real ::  xomega,  xval1, xval2
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
    real function Ricker_Fl(time,tau,f0)
        implicit none
        ! Ricker function for pressure wave: the same as the previous one, but with
        !    one time derivative to be coherent
        real,intent(in)  :: time, tau, f0
        !
        real  :: sigma,sigma2,sigma3

        sigma = M_PI*f0*(time-tau)
        sigma2 = sigma**2
        sigma3 = M_PI*f0

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

end module ssources
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
