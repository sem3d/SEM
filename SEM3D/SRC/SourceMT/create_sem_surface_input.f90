!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module surface_input

implicit none

contains
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
     character(len=20)            :: char
     character(len=70)            :: sourcename
     real(kind=8)                 :: ndir
     integer                      :: i, rg, nsurf, int
     integer, allocatable         :: dummylist(:)
     logical                      :: Boolean(1:20)
     character(len=800)                       :: parametric_var
     character(len=256)           :: FunctionName ='read_surface_input'
     character(len=256)           :: SourceFile = 'create_sem_surface_input'
     character(len=700)           :: ErrorSMS
   
   
     Tdomain%nsurface = config%nsurface
     allocate (Tdomain%nsurfsource(0:Tdomain%nsurface -1))
     
     call c_f_pointer(config%surface, surf)
     nsurf          =0
     Tdomain%n_NEBC =0
     Tdomain%n_PWBC =0
     Tdomain%n_FTBC =0
     Tdomain%n_DIRIC=0
     !! 
     rg = Tdomain%rank
      
     if (rg == 0) then
        write(*,*) 
        write(*,*) "--> READING SURFACE INPUT"
        write(*,*)  
     endif
     
     do while(associated(surf)) 
         
       if ( surf%surface_present /= 0 ) then
          int = nsurf
          write(char,*) int
          if (allocated(dummylist)) deallocate(dummylist)
          if  (surf%surface_whatbc == 1) then
            Tdomain%nsurfsource(nsurf)%what_bc = 'NE'
            if (Tdomain%n_NEBC-1.ge.0) allocate(dummylist(1:Tdomain%n_NEBC))
            if (allocated(Tdomain%list_NEBC)) dummylist = Tdomain%list_NEBC
            if (allocated(Tdomain%list_NEBC)) deallocate(Tdomain%list_NEBC)
            Tdomain%n_NEBC = Tdomain%n_NEBC + 1
            allocate(Tdomain%list_NEBC(1:Tdomain%n_NEBC))
            if (Tdomain%n_NEBC-1.gt.0) Tdomain%list_NEBC(1:Tdomain%n_NEBC-1) = dummylist
            Tdomain%list_NEBC(Tdomain%n_NEBC) = nsurf
            if (rg.eq.0) write(*,1004) "Neumann /-> surface"//adjustl(char(1:len_trim(char)))
     
          elseif (surf%surface_whatbc == 2) then
            Tdomain%nsurfsource(nsurf)%what_bc = 'PW'
            if (Tdomain%n_PWBC-1.ge.0) allocate(dummylist(1:Tdomain%n_PWBC))
            if (allocated(Tdomain%list_PWBC)) dummylist = Tdomain%list_PWBC
            if (allocated(Tdomain%list_PWBC)) deallocate(Tdomain%list_PWBC)
            Tdomain%n_PWBC = Tdomain%n_PWBC + 1
            allocate(Tdomain%list_PWBC(1:Tdomain%n_PWBC))
            if (Tdomain%n_PWBC-1.gt.0) Tdomain%list_PWBC(1:Tdomain%n_PWBC-1) = dummylist
            Tdomain%list_PWBC(Tdomain%n_PWBC) = nsurf
            Tdomain%nsurfsource(nsurf)%wave_type = surf%surface_wave
            Tdomain%nsurfsource(nsurf)%Speed= surf%surface_Speed
            Tdomain%nsurfsource(nsurf)%Kdir = surf%surface_dirU
            if (rg.eq.0) write(*,1004) "Plane Wave /-> surface"//adjustl(char(1:len_trim(char)))
     
          elseif (surf%surface_whatbc == 3) then
            Tdomain%nsurfsource(nsurf)%what_bc = 'FT'
            if (Tdomain%n_FTBC-1.ge.0) allocate(dummylist(1:Tdomain%n_FTBC))
            if (allocated(Tdomain%list_FTBC)) dummylist = Tdomain%list_FTBC
            if (allocated(Tdomain%list_FTBC)) deallocate(Tdomain%list_FTBC)
            Tdomain%n_FTBC = Tdomain%n_FTBC + 1
            if (Tdomain%n_FTBC-1.gt.0) Tdomain%list_FTBC(1:Tdomain%n_FTBC-1) = dummylist
            allocate(Tdomain%list_FTBC(1:Tdomain%n_FTBC))
            Tdomain%list_FTBC(Tdomain%n_FTBC) = nsurf
            if (rg.eq.0) write(*,1004)  "Fault /-> surface"//adjustl(char(1:len_trim(char)))
     
          elseif (surf%surface_whatbc == 4) then
            Tdomain%nsurfsource(nsurf)%what_bc = 'DR'
            if (Tdomain%n_DIRIC-1.ge.0) allocate(dummylist(1:Tdomain%n_DIRIC))
            if (allocated(Tdomain%list_DIRICBC)) dummylist = Tdomain%list_DIRICBC
            if (allocated(Tdomain%list_DIRICBC)) deallocate(Tdomain%list_DIRICBC)
            Tdomain%n_DIRIC = Tdomain%n_DIRIC + 1
            allocate(Tdomain%list_DIRICBC(1:Tdomain%n_DIRIC))
            if (Tdomain%n_DIRIC-1.gt.0) Tdomain%list_DIRICBC(1:Tdomain%n_DIRIC-1) = dummylist
            Tdomain%list_DIRICBC(Tdomain%n_DIRIC) = nsurf
            if (rg.eq.0) write(*,1004)  "Dirichlet /-> surface"//adjustl(char(1:len_trim(char)))
       
          else 
              ErrorSMS = "unknown surface problem (type : ??? ) /-> surface"//adjustl(char(1:len_trim(char)))
              call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
          endif
          
          Tdomain%nsurfsource(nsurf)%mat_index = surf%surface_mat
           
          Boolean = .false.
          do i=lbound(surf%surface_list,1),ubound(surf%surface_list,1)
             Boolean(i) = surf%surface_list(i) /= 0 
          enddo
     
          ErrorSMS = "unknown associeted surface to "//adjustl(sourcename(1:len_trim(sourcename)))
          if (Count(Boolean) == 0) call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
               
          allocate(Tdomain%nsurfsource(nsurf)%index(1:Count(Boolean)))
          Tdomain%nsurfsource(nsurf)%index = surf%surface_list(1:Count(Boolean))-1
          if (Tdomain%nsurfsource(nsurf)%what_bc /= 'DR') then
             if (surf%surface_type==2) then
                Tdomain%nsurfsource(nsurf)%wtype = 'R'
                Tdomain%nsurfsource(nsurf)%f0 = surf%surface_f0
                Tdomain%nsurfsource(nsurf)%Rickertau=surf%Rtau
                sourcename = "Ricker in time"
             elseif (surf%surface_type==1) then
                Tdomain%nsurfsource(nsurf)%wtype = 'G'
                sourcename = "Gaussian in time"
             elseif (surf%surface_type==14) then
                Tdomain%nsurfsource(nsurf)%wtype = 'A'
             endif
       
             ndir = sqrt(surf%surface_L(1)**2+surf%surface_L(2)**2+surf%surface_L(3)**2)
             if (ndir>0.0) Tdomain%nsurfsource(nsurf)%dir(0:2) = surf%surface_L(1:3)/ndir
             Tdomain%nsurfsource(nsurf)%scoord(0:2) = surf%surface_C(1:3)
             Tdomain%nsurfsource(nsurf)%amplitude = surf%amplitude
     
             if (Tdomain%nsurfsource(nsurf)%wtype == 'A') then
               
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
                    if (Tdomain%nsurfsource(nsurf)%dim == 3) then
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
     
                 Tdomain%nsurfsource(nsurf)%shape = surf%surface_space
                 Tdomain%nsurfsource(nsurf)%size = surf%surface_size
     
                 if (rg==0) then
                          write(*,*) "Surface BC Source  : "//adjustl(sourcename(1:len_trim(sourcename)))
                          write(*,2019) " Direction          : ", Tdomain%nsurfsource(nsurf)%dir
                     
                     sourcename = fromcstr(surf%surface_name)
                      write(*,2018) " Load shape         : ", sourcename(1:len_trim(sourcename))// &
                                                  adjustl(" (R="),Tdomain%nsurfsource(nsurf)%size,")" 
                     if (Tdomain%nsurfsource(nsurf)%f0/=0) write(*,2019) "Ricker Freqence    : ", &
                                                                   Tdomain%nsurfsource(nsurf)%f0 
                     if (Tdomain%nsurfsource(nsurf)%Rickertau/=0) write(*,2019) " Ricker tau         : ", &
                                                               Tdomain%nsurfsource(nsurf)%Rickertau
                  endif
             endif
     
             if (rg==0) then
                 write(*,2015) " Source Ref coord   : ", Tdomain%nsurfsource(nsurf)%scoord
                 write(*,2016) " Associated surface : ",(Tdomain%nsurfsource(nsurf)%index(i), i=1,Count(Boolean))
                 write(*,2015) " Load amplitude     : ", Tdomain%nsurfsource(nsurf)%amplitude
             endif
          
          else
              write(*,2016) "Associated surface : ",(Tdomain%nsurfsource(nsurf)%index(i), i=1,Count(Boolean))
              write(*,*)
          endif
          nsurf = nsurf + 1
          Tdomain%logicD%surfBC = .true.
          write(*,*)
       endif
       call c_f_pointer(surf%next, surf)
       
     enddo
     
     if (Tdomain%nsurface /= nsurf) Tdomain%nsurface = nsurf
     if (allocated(dummylist)) deallocate(dummylist)
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
      use Alertes

      implicit none
      type(domain), intent(inout) :: Tdomain
      integer                     :: ns, index, surf , s, n
      character(len=20)           :: char
      character(len=256)          :: FunctionName ='define_surface_properties'
      character(len=256)          :: SourceFile = 'create_sem_surface_input'
      character(len=700)          :: ErrorSMS
      
      ! elastic properties
      do ns = 0, Tdomain%nsurface-1
         index = Tdomain%nsurfsource(ns)%mat_index
         do s = lbound(Tdomain%nsurfsource(ns)%index,1),ubound(Tdomain%nsurfsource(ns)%index,1)
            surf = Tdomain%nsurfsource(ns)%index(s)
            write(char,*) surf
            do n=0,size(Tdomain%sSurfaces)-1
               if (Tdomain%sSurfaces(n)%name == "surface"//adjustl(char(:len_trim(char)))) then
                   Tdomain%sSurfaces(n)%Elastic%lambda  = Tdomain%sSubDomain(index)%DLambda
                   Tdomain%sSurfaces(n)%Elastic%Mu      = Tdomain%sSubDomain(index)%DMu
                   Tdomain%sSurfaces(n)%Elastic%Sspeed  = Tdomain%sSubDomain(index)%Sspeed
                   Tdomain%sSurfaces(n)%Elastic%Pspeed  = Tdomain%sSubDomain(index)%Pspeed   
                   Tdomain%sSurfaces(n)%Elastic%density = Tdomain%sSubDomain(index)%dDensity 
                   if (Tdomain%nsurfsource(ns)%what_bc=='PW') then
                       call PlaneWaveSpeed(Tdomain%nsurfsource(ns),Tdomain%sSurfaces(n))
                   endif
                   !! print parameters values in the screen
                   ErrorSMS = "unknown associeted domain to "//adjustl(Tdomain%sSurfaces(n)%name)
                   if (Tdomain%nsurfsource(ns)%mat_index<0) call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
                   if (Tdomain%rank == 0) then
                       write(*,*)
                       write(*,2004) "SURFACE -> ", surf
                       write(*,2007) "Domaim  : " , Tdomain%nsurfsource(ns)%mat_index 
                       write(*,2009) "Lambda  : " , Tdomain%sSurfaces(n)%Elastic%lambda 
                       write(*,2009) "Mu      : " , Tdomain%sSurfaces(n)%Elastic%Mu
                       write(*,2009) "Speed Vp: ",  Tdomain%sSurfaces(n)%Elastic%Pspeed
                       write(*,2009) "Speed Vs: ",  Tdomain%sSurfaces(n)%Elastic%Sspeed
                       if (Tdomain%nsurfsource(ns)%what_bc=='PW') then 
                           write(*,2009) "Speed PW: ", Tdomain%sSurfaces(n)%Elastic%PWspeed
                           write(*,2009) "PW dirK : ", Tdomain%nsurfsource(ns)%dir
                           write(*,2009) "PW dirD : ", Tdomain%nsurfsource(ns)%Kdir
                       endif
                       write(*,2009) "Density : ",  Tdomain%sSurfaces(n)%Elastic%density
                   endif
                endif
            enddo
         enddo
      end do
      write(*,*)  
      
      include 'formats.in'
   
   end subroutine define_surface_properties
   !----------------------------------------------------------------------------------
   !----------------------------------------------------------------------------------
   subroutine PlaneWaveSpeed(surfsrc,surf)
      
      use ssurf
      
      implicit none
      type(SurfaceParam),        intent(inout) :: surfsrc
      type(SurfaceT),            intent(inout) :: surf
      real(kind=8), dimension(0:2)             :: dirP, dirD
      real(kind=8)                             :: dot
      
      dirP = surfsrc%dir
      
      if (surfsrc%wave_type==1) then
             dirD = dirP
             surf%Elastic%PWspeed = surf%Elastic%Pspeed
             write(*,*)
             write(*,*) " LONGITUDINAL INCIDENT PLANE WAVE "
             write(*,*)
      elseif (surfsrc%wave_type==2) then
             dirD = DirP
             if ((dirP(1)**2 + dirP(2)**2) /=0.d0) dirD(0) = -dirP(0) /(dirP(1)**2 + dirP(2)**2 )
             if ((dirP(0)**2 + dirP(2)**2) /=0.d0) dirD(0) = -dirP(1) /(dirP(0)**2 + dirP(2)**2 )
             if ((dirP(0)**2 + dirP(1)**2) /=0.d0) dirD(0) = -dirP(0) /(dirP(0)**2 + dirP(1)**2 )
             surf%Elastic%PWspeed = surf%Elastic%Sspeed
             write(*,*)
             write(*,*) " TRANSVERSAL INCIDENT PLANE WAVE "
             write(*,*)
      elseif (surfsrc%wave_type==3) then
             dirD(0)=0.d0; dirD(1)=0.d0; dirD(2) = 1.d0
             surf%Elastic%PWspeed = surf%Elastic%Sspeed
             write(*,*)
             write(*,*) " SH-WAVE INCIDENT PLANE WAVE "
             write(*,*)
      elseif (surfsrc%wave_type==4) then
             dirD(0)=-dirP(1); dirD(1)=dirP(0); dirD(2)=0.d0
             surf%Elastic%PWspeed = surf%Elastic%Sspeed
             write(*,*)
             write(*,*) " SV-WAVE INCIDENT PLANE WAVE "
             write(*,*)
      else
          surf%Elastic%PWspeed = surf%Elastic%Sspeed
          dirD = surfsrc%Kdir
          write(*,*)
          write(*,*) " SV-WAVE INCIDENT PLANE WAVE "
          write(*,*)
      endif
      
      surfsrc%Kdir =  dirD
      
   end subroutine PlaneWaveSpeed
   !----------------------------------------------------------------------------------
   !----------------------------------------------------------------------------------
   subroutine  surface_in_list(Tdomain)
      
      use sdomain
      use Alertes

      implicit none
      type(domain), intent(in)    :: Tdomain
      integer                     :: n, ns, s, int 
      logical                     :: find
      character(len=100)          :: list_index
      character(len=256)          :: FunctionName ='surface_in_list'
      character(len=256)          :: SourceFile = 'create_sem_surface_input'
      character(len=700)          :: ErrorSMS
      character(len=20)           :: char, string
    
      list_index =' '
      do ns = 0, Tdomain%nsurface-1
         do s = lbound(Tdomain%nsurfsource(ns)%index,1),ubound(Tdomain%nsurfsource(ns)%index,1)
            find = .false.
            int = Tdomain%nsurfsource(ns)%index(s)
            write(char,*) int
            block : &
            do n=0,size(Tdomain%sSurfaces)-1
               string = Tdomain%sSurfaces(n)%name(8:len_trim(Tdomain%sSurfaces(n)%name))
               list_index = trim(list_index)//adjustl(string(1:len_trim(string)))//";"
               if (Tdomain%sSurfaces(n)%name == "surface"//adjustl(char(:len_trim(char)))) then
                  find = .true.
                  exit block
               endif
            enddo block
            ErrorSMS = "surface (index ="//adjustl(char(1:len_trim(char)))//adjustl(") is not found in tag list: ")// &
                       adjustl(list_index(:len_trim(list_index)))
            if (.not.find) call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
         enddo
      enddo
      
   end subroutine surface_in_list

end module surface_input

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
