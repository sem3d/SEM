!>
!!\file Sources.f90
!!\brief Assure la gestion des sources.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssources

    ! Modified by Gaetano Festa 31/01/2005
    ! Modified by Paul Cupillard 05/06/2005

    type :: Source

       !   modif pour benchmark can2
       character(len = 30)  :: time_file
       !   modif pour benchmark can2
       integer :: i_dir, i_type_source, i_time_function, elem, proc
       integer, dimension(0:2) :: gll
       real :: tau_b,cutoff_freq, radius,realcolat,reallong,refcolat,reflong
       !   ajout de parametres pour definir Gabor signal source de type 4
       !   ajout de gama et ts
       real ::  gama, ts
       !  modif mariotti fevrier 2007 cea
       real :: Xsource,YSource,Zsource

       real, dimension(0:3) :: fh
       real, dimension(0:2) :: refcoord
       !!real, dimension(:), allocatable :: timefunc !!gfortran rale
       real, dimension(:), pointer :: timefunc
       real, dimension(0:2,0:2) :: Moment, InvGrad
       real, dimension (:,:,:,:), pointer :: coeff
       !   modif pour benchmark can2
       real, dimension(:), pointer :: ampli, time
       !   modif pour benchmark can2

    end type Source

contains

    ! ###############################################
    !>
    !! \fn function CompSource (Sour,time,np)
    !! \brief
    !!
    !! \param type (source) Sour
    !! \param integer np
    !! \param real time
    !<

    !>
    !! \fn function CompSource (Sour,time,np)
    !! \brief
    !!
    !! \param type (source) Sour
    !! \param integer np
    !! \param real time
    !<


    real function CompSource (Sour,time)

        type (source) :: Sour
        real :: time

        select case (Sour%i_time_function)
        case (1)
            CompSource = Gaussian (time,Sour%tau_b)
        case (2)
            CompSource = Ricker (time,Sour%tau_b,Sour%cutoff_freq)
        case (3)
            !    deja pris par ancienne fonction
            CompSource = 0.
        case (4)
            !   developpement de la source du benchmark can1
            CompSource = Gabor (time,Sour%tau_b,Sour%cutoff_freq,Sour%gama,Sour%ts)
        case (5)
            !   modif pour benchmark can2
            CompSource = Source_File (time,Sour%tau_b,Sour)
            !   modif pour benchmark can2
        end select

        return
    end function CompSource

    ! ################################################
    !>
    !! \fn function CompSource (Sour,time,np)
    !! \brief
    !!
    !! \param type (source) Sour
    !! \param integer np
    !! \param real time
    !<

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
        type(source), intent(in)  :: Sour
        real, intent(in)          :: tt
        real :: tau
        integer :: iflag

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

    subroutine read_source(Sour)
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

    end subroutine read_source
    !   modif pour benchmark can2
    !-------------------------------------------------

    real function Gaussian (time,tau)

        real :: tau,time
        real, parameter :: pi = 3.141592653

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


        real :: time, tau, f0
        real :: sigma
        real, parameter :: pi = 3.141592653

        if ( time < 2.5*tau ) then
            sigma = pi * f0 * (time-tau)
            sigma = sigma**2
!!!! Ricker = (1 - 2*sigma) * exp(-sigma) !version origine
            Ricker = 1e5* (1 - 2*sigma) * exp(-sigma) !version Gsa Ipsis (amplitude)
            !print*,' ricker ',pi, time, tau, f0,sigma,exp(-sigma), Ricker
        else
            Ricker = 0.
        endif

        return
    end function Ricker

    ! #################################################

    real function Gabor (time,tau,fp,gama,ts)


        real :: time, tau2, fp,gama,ts
        real :: sigma
        real ::  xomega,  xval1, xval2
        real, parameter :: pi = 3.141592653
        xomega  = pi*0.5

        if ( time < 32. ) then
            sigma = 2. * pi * fp * (time-ts)
            xval1 = cos(sigma + xomega)
            sigma = sigma/gama
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

end module ssources
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
