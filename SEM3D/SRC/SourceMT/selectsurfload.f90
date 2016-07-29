!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file Neumann.f90
!! \brief
!!
!<
!!
module Surface_prbl_type

use constants, only : M_PI

implicit none

contains
     !----------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------
     subroutine Neumanforce(Param,srcshape,Sourcef,coord,Btn,time,forces)
     
     use ssurf
     use Mathfval
     use parameters, only: Addparametricvar
     use ssources
     
     implicit none
     type(SurfaceParam),            intent(in) :: Param
     type(FoncValue),               intent(in) :: Sourcef
     real(kind=8),                  intent(in) :: srcshape, time
     real(kind=8), dimension(0:2),  intent(in) :: coord, Btn
     real(kind=8), dimension(0:2), intent(out) :: forces
     type (source)                             :: Sour
     
     forces = 0.d0
     select case(Param%wtype)
            case('R')
                !- Ricker in time, source uniformly distributed..
                Sour%i_time_function = 2
                Sour%cutoff_freq     = Param%f0
                Sour%tau_b           = Param%Rickertau
                Sour%amplitude_factor= Param%amplitude
                forces = Param%dir*CompSource(Sour,time,0)*srcshape
            case('G')
                !- Gaussian in time, source uniformly distributed...
                Sour%i_time_function = 1
                Sour%tau_b           = Param%Rickertau
                Sour%amplitude_factor= Param%amplitude
                forces = Param%dir*CompSource(Sour,time,0)*srcshape
            case ('A')
                !- analytical form
                CALL ffvalue(Sourcef,coord,time)
                if ((Sourcef%dim==3).and.(Sourcef%source == 'M')) then
                     forces(0) = (Sourcef%fvalue(1)*Btn(0)+ Sourcef%fvalue(4)*Btn(1)+Sourcef%fvalue(6)*Btn(2))
                     forces(1) = (Sourcef%fvalue(4)*Btn(0)+ Sourcef%fvalue(2)*Btn(1)+Sourcef%fvalue(5)*Btn(2))
                     forces(2) = (Sourcef%fvalue(6)*Btn(0)+ Sourcef%fvalue(5)*Btn(1)+Sourcef%fvalue(3)*Btn(2))
                elseif ((Sourcef%dim==3).and.(Sourcef%source == 'F')) then
                     forces(0) = Sourcef%fvalue(1)
                     forces(1) = Sourcef%fvalue(2)
                     forces(2) = Sourcef%fvalue(3)
                elseif ((Sourcef%dim==1).and.(Sourcef%source=='F')) then
                     forces = Param%dir*Sourcef%fvalue(1)
                     if (Sourcef%stat=='MIXT') forces = forces*srcshape
                endif
     end select
     
     end subroutine Neumanforce
     !----------------------------------------------------------------------------------
     !----------------------------------------------------------------------------------
     subroutine PlaneWavedispl(Param,Sourcef,coord,Btn,time,PWspeed,displ,char)
     
     use ssurf
     use Mathfval
     use parameters, only: Addparametricvar
     use ssources
     use Alertes
      
     implicit none
     type(SurfaceParam),            intent(in) :: Param
     type(FoncValue),               intent(in) :: Sourcef
     real(kind=8),                  intent(in) :: time,PWspeed
     real(kind=8), dimension(0:2),  intent(in) :: coord, Btn
     character(len=4)            ,  intent(in) :: char
     real(kind=8), dimension(0:2), intent(out) :: displ
     type (source)                             :: Sour
     real(kind=8)                              :: dot, result
     character(len=256)                        :: FunctionName ='PlaneWanedispl'
     character(len=256)                        :: SourceFile = 'selectsurfload'
     character(len=700)                        :: ErrorSMS

       
      displ = 0.d0
      dot = coord(0)*Param%dir(0)+coord(1)*Param%dir(1)+coord(2)*Param%dir(2) - PWspeed*time
      select case(Param%wtype)
              case('R')
                   !- Ricker in time, source uniformly distributed..
                    Sour%i_time_function = 2
                    Sour%cutoff_freq     = Param%f0
                    Sour%tau_b           = Param%Rickertau
                    Sour%amplitude_factor= Param%amplitude
                    if (char=='wave') then
                        call CompSourcePW(Sour,dot, result)
                    else
                        call CompdisplPW(Sour,dot, result)
                    endif
                    displ = Param%Kdir*result
              case('G')
                    !- Gaussian in time, source uniformly distributed...
                    Sour%i_time_function = 1
                    Sour%tau_b           = Param%Rickertau
                    Sour%amplitude_factor= Param%amplitude
                    if (char=='wave') then
                        call CompSourcePW(Sour,dot, result)
                    else
                        call CompdisplPW(Sour,dot, result)
                    endif
                    displ = Param%Kdir*result
             case ('A')
                    !- analytical form
                    CALL ffvalue(Sourcef ,coord,time)
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
     subroutine CompSourcePW(Sour,dot, result)
     
     use ssources

     implicit none
     real(kind=8), intent(in) :: dot
     type (source),intent (in):: Sour
     real(kind=8) ,intent(out):: result
     real(kind=8)             :: f0
     
     result = 0.d0 
     select case(Sour%i_time_function)
           case(1)
                
           case(2)
               f0 = Sour%cutoff_freq
               call dRickerPW(dot,f0,result)
           case(3)

     end select

     endsubroutine CompSourcePW
     !----------------------------------------------------------------------------------
     subroutine CompdisplPW(Sour,dot, result)
     
     use ssources
     
     implicit none
     real(kind=8), intent(in) :: dot
     type (source),intent (in):: Sour
     real(kind=8) ,intent(out):: result
     real(kind=8)             :: f0
     
     result = 0.d0
     select case(Sour%i_time_function)
           case(1)
     
           case(2)
               f0 = Sour%cutoff_freq
               call RickerPW(dot,f0,result)
           case(3)
     
     end select
     
     endsubroutine CompdisplPW
     !----------------------------------------------------------------------------------
     subroutine RickerPW(xx, f0, result)
     
     implicit none
                  
     real(kind=8), intent(in) :: f0, xx
     real(kind=8), intent(out):: result
     real(kind=8)             :: sigma
                                     
     result = 0.d0
     sigma = M_PI**2*f0**2*xx**2
     result = - (1 - 2.d0*sigma)*dexp(-sigma)

     endsubroutine RickerPW
     !----------------------------------------------------------------------------------
     subroutine dRickerPW (xx, f0, result)

     implicit none
     real(kind=8), intent(in) :: f0, xx
     real(kind=8), intent(out):: result
     real(kind=8)             :: sigma
        
        sigma  = M_PI**2*f0**2
        result = 2.d0*sigma*(3.d0*xx - 2.d0*sigma*xx**3)*dexp(-sigma*xx**2)

     endsubroutine dRickerPW
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
                                  
