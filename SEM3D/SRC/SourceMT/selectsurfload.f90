!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module Surface_prbl_type

    use constants, only : M_PI, fpp

    implicit none

contains
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine init_math_src(Param, sourcef)
        use Mathfval
        use ssurf
        use parameters, only: Addparametricvar

        implicit none
        type(SurfaceParam), intent(in) :: Param
        type(FoncValue),    intent(out):: Sourcef

        if (Param%wtype == 'A') then
            Sourcef%dim    =Param%dim
            Sourcef%source =Param%source
            Sourcef%var    =Param%varia
            Sourcef%valuefx(1:len_trim(Param%funcx))=Param%funcx(1:len_trim(Param%funcx))
            Sourcef%valuefy(1:len_trim(Param%funcy))=Param%funcy(1:len_trim(Param%funcy))
            Sourcef%valuefz(1:len_trim(Param%funcz))=Param%funcz(1:len_trim(Param%funcz))
            Sourcef%valuefxy(1:len_trim(Param%funcxy))=Param%funcxy(1:len_trim(Param%funcxy))
            Sourcef%valuefyz(1:len_trim(Param%funcyz))=Param%funcyz(1:len_trim(Param%funcyz))
            Sourcef%valuefxz(1:len_trim(Param%funcxz))=Param%funcxz(1:len_trim(Param%funcxz))
            Addparametricvar%nparam=0
            if (Param%paramvar==1) then
                Addparametricvar%nparam =Param%nparamvar
                Addparametricvar%paramname =Param%paramname
                Addparametricvar%paramvalue =Param%paravalue
            endif
            select case (Sourcef%dim)
            case (1)
                allocate(Sourcef%fvalue(1:1))
            case (2)
                if (Sourcef%source=='F') then
                    allocate(Sourcef%fvalue(1:2))
                elseif (Sourcef%source == 'M') then
                    allocate(Sourcef%fvalue(1:3))
                endif
            case(3)
                if (Sourcef%source == 'F') then
                    allocate(Sourcef%fvalue(1:3))
                elseif (Sourcef%source == 'M') then
                    allocate(Sourcef%fvalue(1:6))
                endif
            end select
            Sourcef%stat='UNIF'
            if ((Sourcef%dim==1).and.(Param%shape/=0)) Sourcef%stat='MIXT'
        endif

    end subroutine init_math_src
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine Neumanforce(Param,srcshape,coord,Btn,time,forces)

        use ssurf
        use Mathfval
        use parameters, only: Addparametricvar
        use ssources

        implicit none
        type(SurfaceParam),            intent(in) :: Param
        real(fpp),                  intent(in) :: srcshape, time
        real(fpp), dimension(0:2),  intent(in) :: coord, Btn
        real(fpp), dimension(0:2), intent(out) :: forces
        real(fpp), dimension(1:6)              :: sigma
        type(FoncValue) :: Sourcef
        type (source)   :: Sour

        forces = 0.d0
        sigma  = 0.d0
        select case(Param%wtype)
        case('R')
            !- Ricker in time, source uniformly distributed..
            Sour%i_time_function = 2
            Sour%cutoff_freq     = Param%f0
            Sour%tau_b           = Param%Rickertau
            Sour%amplitude_factor= Param%amplitude
            sigma(1:3) = Param%dir*CompSource(Sour,time,0)*srcshape
        case('G')
            !- Gaussian in time, source uniformly distributed...
            Sour%i_time_function = 1
            Sour%tau_b           = Param%Rickertau
            Sour%amplitude_factor= Param%amplitude
            sigma(1:3) = Param%dir*CompSource(Sour,time,0)*srcshape
        case ('A')
            !- analytical form
            call init_math_src(Param, sourcef)
            call ffvalue(Sourcef,coord,time)
            if ((Sourcef%dim==3).and.(Sourcef%source == 'M')) then
                sigma = Sourcef%fvalue

            elseif ((Sourcef%dim==3).and.(Sourcef%source == 'F')) then
                sigma(1:3) = Sourcef%fvalue(1:3)

            elseif ((Sourcef%dim==1).and.(Sourcef%source=='F')) then
                sigma(1:3) = Param%dir*Sourcef%fvalue(1)

                if (Sourcef%stat=='MIXT') sigma = sigma*srcshape
            endif
        end select

        forces(0) = (sigma(1)*Btn(0)+ sigma(4)*Btn(1)+sigma(6)*Btn(2))
        forces(1) = (sigma(4)*Btn(0)+ sigma(2)*Btn(1)+sigma(5)*Btn(2))
        forces(2) = (sigma(6)*Btn(0)+ sigma(5)*Btn(1)+sigma(3)*Btn(2))

    end subroutine Neumanforce
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine PlaneWaveDerive(Param, coord, time, PWspeed, deriveU, veloc0)

        use ssurf
        use Mathfval
        use parameters, only: Addparametricvar
        use Alertes

        implicit none
        type(SurfaceParam),                   intent(in) :: Param
        real(fpp),                         intent(in) :: time,PWspeed
        real(fpp), dimension(0:2)       ,  intent(in) :: coord
        real(fpp),dimension(0:2),         intent(out) :: deriveU, veloc0
        type(FoncValue)     :: Sourcef
        real(fpp)        :: dot
        character(len=256)  :: FunctionName ='PlaneWanedispl'
        character(len=256)  :: SourceFile = 'selectsurfload'
        character(len=700)  :: ErrorSMS

        deriveU = 0.d0
        dot = coord(0)*Param%dir(0)+coord(1)*Param%dir(1)+coord(2)*Param%dir(2) - PWspeed*time

        select case(Param%wtype)
        case('R')
            !- Ricker in time, source uniformly distributed..
            deriveU = Param%Kdir*RickerPW(dot, Param%f0, PWspeed, 'w')*Param%amplitude
            veloc0  = Param%Kdir*RickerPW(dot, Param%f0, PWspeed, 'v')*Param%amplitude
        case('G')
            !- Gaussian in time, source uniformly distributed...

        case('A')
            !- analytical form
            call init_math_src(Param, sourcef)
            call ffvalue(Sourcef ,coord,time)
            ErrorSMS= " Only the displacement field are needed for plane wave problem "
            if (Sourcef%source == 'M') call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
            if ((Sourcef%dim==3).and.(Sourcef%source == 'F')) then
                deriveU(0) = Sourcef%fvalue(1)*Param%Kdir(0)
                deriveU(1) = Sourcef%fvalue(2)*Param%Kdir(1)
                deriveU(2) = Sourcef%fvalue(3)*Param%Kdir(2)
            elseif ((Sourcef%dim==1).and.(Sourcef%source == 'F')) then
                deriveU = Param%Kdir*Sourcef%fvalue(1)
            endif
        end select
    end subroutine PlaneWaveDerive
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine PlaneWavedispl(Param,coord, time, PWspeed, displ, veloc, accel)

        use ssurf
        use Mathfval
        use parameters, only: Addparametricvar
        use ssources
        use Alertes

        implicit none
        type(SurfaceParam),                   intent(in) :: Param
        real(fpp),                         intent(in) :: time,PWspeed
        real(fpp), dimension(0:2)       ,  intent(in) :: coord
        real(fpp), dimension(0:2)       ,  intent(out):: displ, accel, veloc
        type(FoncValue)    :: Sourcef
        type (source)      :: Sour
        real(fpp)       :: dot, result
        character(len=256) :: FunctionName ='PlaneWanedispl'
        character(len=256) :: SourceFile = 'selectsurfload'
        character(len=700) :: ErrorSMS

        displ = 0.d0
        veloc = 0.d0
        accel = 0.d0
        dot = coord(0)*Param%dir(0)+coord(1)*Param%dir(1)+coord(2)*Param%dir(2) - PWspeed*time
        select case(Param%wtype)
        case('R')
            !- Ricker in time, source uniformly distributed..
            displ = Param%Kdir*RickerPW(dot, Param%f0, PWspeed, 'd')*Param%amplitude
            veloc = Param%Kdir*RickerPW(dot, Param%f0, PWspeed, 'v')*Param%amplitude
            accel = Param%Kdir*RickerPW(dot, Param%f0, PWspeed, 'a')*Param%amplitude
        case('G')


        case ('A')
            !- analytical form
            call init_math_src(Param, sourcef)
            call ffvalue(Sourcef ,coord,time)
            ErrorSMS= " Only the displacement field are needed for plane wave problem "
            if (Sourcef%source == 'M') call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
            if ((Sourcef%dim==3).and.(Sourcef%source == 'F')) then
                displ(0) = Sourcef%fvalue(1)*Param%Kdir(0)
                displ(1) = Sourcef%fvalue(2)*Param%Kdir(1)
                displ(2) = Sourcef%fvalue(3)*Param%Kdir(2)
            elseif ((Sourcef%dim==1).and.(Sourcef%source == 'F')) then
                displ = Param%Kdir*Sourcef%fvalue(1)
            endif
        end select
    end subroutine PlaneWavedispl
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    real(fpp) function RickerPW(xx, f0, C0, char)  result(Ricker)

        implicit none
        real(fpp), intent(in) :: f0, xx,C0
        character   , intent(in) :: char
        real(fpp)             :: sigma

        Ricker = 0.d0
        select case(char)
        case('d')
            !! expression exacte du déplacement
            sigma  = M_PI**2*f0**2*xx**2
            Ricker = - (1.d0 - 2.d0*sigma)*exp(-sigma)
        case('w')
            !! dérive spatiale
            sigma  = M_PI**2*f0**2
            Ricker = 2.d0*sigma*(3.d0*xx - 2.d0*sigma*xx**3)*exp(-sigma*xx**2)
        case('v')
            !! vitesse
            sigma  = M_PI**2*f0**2
            Ricker = -2.d0*C0*sigma*(3.d0*xx - 2.d0*sigma*xx**3)*exp(-sigma*xx**2)
        case('a')
            !! accélération
            sigma  = M_PI**2*f0**2
            Ricker = 2.d0*C0**2*sigma*(3.d0 - 12.d0*sigma*xx**2+4.d0*sigma**2*xx**4)*exp(-sigma*xx**2)
        end select

    end function RickerPW
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
end module Surface_prbl_type

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
