!>
!!\file Sources.F90
!!\brief Assure la gestion des sources.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssources

    ! Modified by Gaetano 13/9/04

    type :: elem_source
       integer :: nr
       real :: eta,xi
       real :: invE,nu
       real, dimension (0:1,0:1) :: Scoeff
       real, dimension (:,:,:), pointer :: ExtForce
    end type elem_source

    type :: Source
       integer :: i_type_source, i_time_function,ine
       real, dimension(2) :: dir
       real, dimension (0:1,0:1) :: moment
       real :: Xsource,Zsource, tau_b,cutoff_freq,amplitude,sigma
       type(elem_source), dimension(:), pointer :: Elem
    end type Source

contains

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

        CompSource = 0.
        select case (Sour%i_time_function)
        case (1)
            CompSource = Gaussian (time,Sour%tau_b,Sour%cutoff_freq)
        case (2)
            CompSource = Ricker (time,Sour%tau_b,Sour%cutoff_freq)
        case (3)
            CompSource = 1
        end select
        CompSource = Sour%amplitude*CompSource

        return
    end function CompSource

    !>
    !! \fn function Gaussian (time, tau, f0)
    !! \brief
    !!
    !! \param real time
    !! \param real tau
    !! \param real f0
    !<
    real function Gaussian (time, tau, f0)

        real :: tau,time,f0
        real :: sigma,pi

        pi = Acos(-1.)
        sigma = pi * f0 * (time - tau )
        sigma = sigma **2

        Gaussian = (time-tau) * exp (-sigma)

        return
    end function Gaussian

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
        real :: sigma,pi

        pi = Acos(-1.)
        sigma = pi * f0 * (time - tau )
        sigma = sigma **2

        Ricker = (1-2 * sigma) *exp (-sigma)
        if(sigma>50.) then   !Ajout Gsa 0508
            Ricker = 0.
        endif
        return
    end function Ricker



end module ssources
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
