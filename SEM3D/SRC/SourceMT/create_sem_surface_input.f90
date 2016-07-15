!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine read_surface_input(Tdomain, config)
  use sdomain
  use ssurf
  use iso_c_binding
  use sem_c_config 
  use Alertes
   
  implicit none
  type(domain), intent(inout)  :: Tdomain
  type(sem_config), intent(in) :: config
  type(sem_surfaces), pointer  :: surf
  character(Len=MAX_FILE_SIZE) :: fnamef
  character(len=800)           :: parametric_var
  character(len=12)            :: char
  character(len=70)            :: sourcename
  real(kind=8)                 :: ndir
  integer                      :: i, rg, nsurf
  logical                      :: Boolean(1:20)
  character(len=256)           :: FunctionName ='read_surface_input'
  character(len=256)           :: SourceFile = 'create_sem_surface_input'
  character(len=700)           :: ErrorSMS


Tdomain%nsurface = config%nsurface
allocate (Tdomain%nsurfsource(0:Tdomain%nsurface -1))

call c_f_pointer(config%surface, surf)
nsurf= 0
Tdomain%n_NEBC=0
Tdomain%n_PWBC=0
Tdomain%n_FTBC=0
!! Ecriture par un seul proc
rg = Tdomain%rank
 
if (rg.eq.0) then
   write(*,*) 
   write(*,*) "--> READING SURFACE INPUT"
   write(*,*)  
endif

do while(associated(surf))
    
  if (surf%surface_present/=0) then
      
     write(char,*) nsurf
   
     if  (surf%surface_whatbc == 1) then
       Tdomain%nsurfsource(nsurf)%what_bc = 'NE'
       Tdomain%n_NEBC = Tdomain%n_NEBC + 1
       sourcename = "Neumann /-> surface"//adjustl(char(1:len_trim(char)))
     elseif (surf%surface_whatbc == 2) then
       Tdomain%nsurfsource(nsurf)%what_bc = 'PW'
       Tdomain%n_PWBC = Tdomain%n_PWBC + 1
       sourcename =  "Plane Wave /-> surface"//adjustl(char(1:len_trim(char)))
     elseif (surf%surface_whatbc == 3) then
       Tdomain%nsurfsource(nsurf)%what_bc = 'FT'
       Tdomain%n_FTBC = Tdomain%n_FTBC + 1
       sourcename = "Fault /-> surface"//adjustl(char(1:len_trim(char)))
     endif
     
     if (rg.eq.0) write(*,1004) sourcename
      
     Boolean = .false.
     do i=lbound(surf%surface_list,1),ubound(surf%surface_list,1)
        Boolean(i) = surf%surface_list(i)/=0 
     enddo

     ErrorSMS = "unknown associeted surface to "//adjustl(sourcename(1:len_trim(sourcename)))
     if (Count(Boolean) == 0) call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
          
     allocate(Tdomain%nsurfsource(nsurf)%index(1:Count(Boolean)))
     Tdomain%nsurfsource(nsurf)%index = surf%surface_list(1:Count(Boolean))
     if (surf%surface_type==1) then
        Tdomain%nsurfsource(nsurf)%wtype = 'R'
        Tdomain%nsurfsource(nsurf)%f0 = surf%surface_f0
        Tdomain%nsurfsource(nsurf)%Rickertau=surf%Rtau
        sourcename = "Ricker in time"

     elseif (surf%surface_type==2) then
        Tdomain%nsurfsource(nsurf)%wtype = 'G'

     elseif (surf%surface_type==3) then
        Tdomain%nsurfsource(nsurf)%wtype = 'A'

     endif
  
     ndir = sqrt(surf%surface_L(1)**2+surf%surface_L(2)**2+surf%surface_L(3)**2)
     if (ndir>0.0) Tdomain%nsurfsource(nsurf)%dir(0:2) = surf%surface_L(1:3)/ndir
     Tdomain%nsurfsource(nsurf)%scoord(0:2) = surf%surface_C(1:3)
     Tdomain%nsurfsource(nsurf)%amplitude = surf%amplitude

     if (Tdomain%nsurfsource(nsurf)%wtype=='A') then
        Tdomain%nsurfsource(nsurf)%dim   = surf%surface_dim
        Tdomain%nsurfsource(nsurf)%varia(1:len_trim(fromcstr(surf%surface_varia))) = trim(fromcstr(surf%surface_varia))
        Tdomain%nsurfsource(nsurf)%source(1:len_trim(fromcstr(surf%surface_source)))= trim(fromcstr(surf%surface_source))
        Tdomain%nsurfsource(nsurf)%funcx(1:len_trim(fromcstr(surf%surface_funcx))) = trim(fromcstr(surf%surface_funcx))
        
        if (Tdomain%nsurfsource(nsurf)%dim.gt.1) then
            Tdomain%nsurfsource(nsurf)%funcy(1:len_trim(fromcstr(surf%surface_funcy))) = trim(fromcstr(surf%surface_funcy))
        endif
        if (Tdomain%nsurfsource(nsurf)%dim.gt.2) then
            Tdomain%nsurfsource(nsurf)%funcz(1:len_trim(fromcstr(surf%surface_funcz))) = trim(fromcstr(surf%surface_funcz))
        endif
        if (Tdomain%nsurfsource(nsurf)%source == 'M') then
            Tdomain%nsurfsource(nsurf)%funcxy(1:len_trim(fromcstr(surf%surface_funcxy)))= trim(fromcstr(surf%surface_funcxy))
            if (Tdomain%nsurfsource(nsurf)%dim.eq.3) then
                Tdomain%nsurfsource(nsurf)%funcxz(1:len_trim(fromcstr(surf%surface_funcxz)))= trim(fromcstr(surf%surface_funcxz))
                Tdomain%nsurfsource(nsurf)%funcyz(1:len_trim(fromcstr(surf%surface_funcyz)))= trim(fromcstr(surf%surface_funcyz))
            endif
        endif
        !!
        !!
        Tdomain%nsurfsource(nsurf)%paramvar = surf%surface_paramvar
        if (Tdomain%nsurfsource(nsurf)%paramvar==1) then
            Tdomain%nsurfsource(nsurf)%nparamvar = surf%surface_nparamvar
            allocate(Tdomain%nsurfsource(nsurf)%paravalue(1:surf%surface_nparamvar))
            do i=1,surf%surface_nparamvar
               Tdomain%nsurfsource(nsurf)%paravalue(i) = surf%surface_Paravalue(i)
            enddo
            parametric_var(1:len_trim(fromcstr(surf%surface_paramname))) = trim(fromcstr(surf%surface_paramname))
            call Split(parametric_var, surf%surface_nparamvar, Tdomain%nsurfsource(nsurf)%paramname)
        endif
        !!
        !!
        if (rg.eq.0) then
           write(*,*)
           write(*,*) "Surface BC Source analytical defined. "
           write(*,*) "-------------------------------------"
           write(*,1006) "Va = "// adjustr(Tdomain%nsurfsource(nsurf)%varia)
           write(*,1005) "F1 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcx)
           if (Tdomain%nsurfsource(nsurf)%dim.gt.1) then
               write(*,1005) "F2 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcy)
           endif
           if (Tdomain%nsurfsource(nsurf)%dim.gt.2) then
              write(*,1005) "F3 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcz)
           endif
           if (Tdomain%nsurfsource(nsurf)%source.eq."M") then
               write(*,1005) "F4 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcxy)
               if (surf%surface_dim.gt.2) then
                   write(*,1005) "F5 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcxz)
                   write(*,1005) "F6 = "// adjustr(Tdomain%nsurfsource(nsurf)%funcyz)
               endif
           endif
           if (Tdomain%nsurfsource(nsurf)%paramvar==1) then
               write(*,*)
               write(*,1007) "Param"," ParamVal"
               do i=1,Tdomain%nsurfsource(nsurf)%nparamvar
                  write(*,2017) Tdomain%nsurfsource(nsurf)%paramname(i), &
                                Tdomain%nsurfsource(nsurf)%paravalue(i)
               enddo
           endif
         endif
         write(*,*)
     else
         if (rg==0) then
            write(*,*) "Surface BC Source "//adjustl(sourcename(1:len_trim(sourcename)))
            write(*,2014) "Load dir : ", Tdomain%nsurfsource(nsurf)%dir(0), &
                                          Tdomain%nsurfsource(nsurf)%dir(1), &
                                          Tdomain%nsurfsource(nsurf)%dir(2)
           if (Tdomain%nsurfsource(nsurf)%f0/=0) write(*,2014) "Ricker Freqence ", &
                                                              Tdomain%nsurfsource(nsurf)%f0 
           if (Tdomain%nsurfsource(nsurf)%Rickertau/=0) write(*,2014) "Ricker Freqence ", &
                                                        Tdomain%nsurfsource(nsurf)%Rickertau
          endif
     endif

     if (rg==0) then
         write(*,2015) "Source Ref coord   : ", Tdomain%nsurfsource(nsurf)%scoord(0), &
                                                        Tdomain%nsurfsource(nsurf)%scoord(1), &
                                                        Tdomain%nsurfsource(nsurf)%scoord(2)
         write(*,2016) "Associated surface : ",(Tdomain%nsurfsource(nsurf)%index(i)-1, i=1,Count(Boolean))
         write(*,2016) "Load amplitude     : ",  Tdomain%nsurfsource(nsurf)%amplitude
     endif

     nsurf = nsurf + 1
     Tdomain%logicD%surfBC = .true.
     call c_f_pointer(surf%next, surf)
  endif
enddo

if (nsurf==0) Tdomain%logicD%surfBC = .false.
write(*,*)

include 'formats.in'

end subroutine read_surface_input
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine Split(String, nstring, SubStrings)
      implicit none
      
      character(len=*), intent(in)                         :: String
      integer, intent(in)                                  :: nstring
      character, parameter                                 :: Delimiter= " "
      character(len=*), dimension(nstring), intent(inout)  :: SubStrings
      integer                                              :: i, j, k, l, n, ns
      logical                                              :: EndOfLine
      
      do i=1,nstring
         SubStrings(i) = ' '
      end do
      n = 0
      k = 1
      l = len_trim(String)
      EndOfLine = l-k < 0
      do while (.not.EndOfLine)
         j = index(String(k:l),Delimiter)
         if (j == 0) then
            j=l+1
         else
            j=j+k-1
         end if
         n=n+1
         if ((j.ne.k).and.(len_trim(String(k:j-1)).ne.0).and.(n.le.nstring)) then
             SubStrings(n) = String(k:j-1)
             ns =n;
         else
             n  =ns;
         endif
         k=j+1
         EndOfLine=l-k < 0
      end do
         
end subroutine Split
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine define_surface_properties(Tdomain)
  
  use sdomain

  implicit none
  type(domain), intent(inout) :: Tdomain
  integer                     :: i_surf, index


  ! elastic properties
  do i_surf = 0,size(Tdomain%sSurfaces)-1
     index = Tdomain%sSurfaces(i_surf)%Elastic%mat_index
     Tdomain%sSurfaces(i_surf)%Elastic%lambda  = Tdomain%sSubDomain(index)%DLambda
     Tdomain%sSurfaces(i_surf)%Elastic%Mu      = Tdomain%sSubDomain(index)%DMu
     Tdomain%sSurfaces(i_surf)%Elastic%Sspeed  = Tdomain%sSubDomain(index)%Sspeed
     Tdomain%sSurfaces(i_surf)%Elastic%Pspeed  = Tdomain%sSubDomain(index)%Pspeed   
     !! print parameters values in the screen
     if (Tdomain%rank == 0) then
         write(*,*)
         write(*,2004) "SURFACE -> ", i_surf
         write(*,2007) "Domaim  : " , Tdomain%sSurfaces(i_surf)%Elastic%mat_index 
         write(*,2009) "Lambda  : " , Tdomain%sSurfaces(i_surf)%Elastic%lambda                                                       
         write(*,2009) "Mu      : " , Tdomain%sSurfaces(i_surf)%Elastic%Mu
         write(*,2009) "Wave Vp : ",  Tdomain%sSurfaces(i_surf)%Elastic%Pspeed
         write(*,2009) "Wave Vs : ",  Tdomain%sSurfaces(i_surf)%Elastic%Sspeed
     endif
  end do

  include 'formats.in'

end subroutine define_surface_properties

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
