!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file Neumann.f90
!! \brief
!!
!<

module surface_load

use Alertes

implicit none

contains
    
    function surfbool(surf, source)
      use sinterface
      use ssurf

      implicit none
      logical                        :: surfbool
      type(SurfaceT), intent(in)     :: surf
      type(SurfaceParam), intent(in) :: source
      integer                        :: surfi, i
      character(len=12)              :: char

      surfbool = .false.
      do i=1,size(source%index)
         surfi = source%index(i)-1
         write(char,*) surfi
         if (surf%name == "surface"//adjustl(char(1:len_trim(char))) ) then
             surfbool = .true.
             exit
         endif
      enddo

    end function
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine add_surface_force(Tdomain)
        use sdomain
        use constants

        implicit none
        type(domain), intent(inout)  :: Tdomain
        integer                      :: i_surf, n1, n2, n3, n4, ns
        character(len=256)                        :: FunctionName ='add_surface_force'
        character(len=256)                        :: SourceFile = 'compute_surface_BC'
        character(len=700)                        :: ErrorSMS


        if (Tdomain%n_NEBC/=0) then
           !! Add Neumann boundary condition
           do ns=0,size(Tdomain%list_NEBC)-1
              n2=Tdomain%list_NEBC(ns+1)
              if (Tdomain%nsurfsource(n2)%what_bc == 'NE') then
                 do n1 = 1,size(Tdomain%nsurfsource(n2)%index)
                    i_surf = Tdomain%nsurfsource(n2)%index(n1)-1
            
                    select case (Tdomain%sSurfaces(i_surf)%domain)
                           case(DM_SOLID)
                                call surface_force(Tdomain%sSurfaces(i_surf)%surf_sl, Tdomain%nsurfsource(n2), &
                                                   Tdomain%sSurfaces(i_surf), Tdomain)
                           case(DM_FLUID)
                                call surface_force(Tdomain%sSurfaces(i_surf)%surf_fl, Tdomain%nsurfsource(n2), &
                                                   Tdomain%sSurfaces(i_surf), Tdomain)
                           case(DM_SOLID_PML)
                                call surface_force(Tdomain%sSurfaces(i_surf)%surf_spml, Tdomain%nsurfsource(n2), &
                                                   Tdomain%sSurfaces(i_surf), Tdomain)
                           case(DM_FLUID_PML)
                                call surface_force(Tdomain%sSurfaces(i_surf)%surf_fpml, Tdomain%nsurfsource(n2), &
                                                   Tdomain%sSurfaces(i_surf), Tdomain)
                   end select
                 enddo 
              endif
           enddo
        endif

        if (Tdomain%n_PWBC /= 0) then
           ErrorSMS= "Sorry the plane wave problem is not yet implemented"
           call ErrorMessage(ErrorSMS,FunctionName,SourceFile)

           do ns=0,size(Tdomain%list_PWBC)-1
              n2=Tdomain%list_PWBC(ns+1)
              do n1 = 1,size(Tdomain%nsurfsource(n2)%index)
                 i_surf = Tdomain%nsurfsource(n2)%index(n1)-1
                 select case (Tdomain%sSurfaces(i_surf)%domain)
                       case(DM_SOLID)

                       case(DM_FLUID)

                       case(DM_SOLID_PML)

                       case(DM_FLUID_PML)

                 end select
              enddo
           enddo
        endif

        if (Tdomain%n_FTBC /= 0) then
           ErrorSMS= "Sorry the fault problem is not yet implemented"
           call ErrorMessage(ErrorSMS,FunctionName,SourceFile)

           do ns=0,size(Tdomain%list_PWBC)-1
              n2=Tdomain%list_PWBC(ns+1)
              do n1 = 1,size(Tdomain%nsurfsource(n2)%index)
                 i_surf = Tdomain%nsurfsource(n2)%index(n1)-1
                 select case (Tdomain%sSurfaces(i_surf)%domain)
                        case(DM_SOLID)
                           
                        case(DM_FLUID)
                             
                        case(DM_SOLID_PML)
              
                        case(DM_FLUID_PML)
                                  
                 end select
              enddo
            enddo
        endif

    end subroutine add_surface_force
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine get_surf_gll_coord(surf,Tdomain, Coord)
        use sdomain
        use sinterface
        use constants

        implicit none
        real(kind=8), dimension(:,:), allocatable :: coord
        type(domain), intent(in)                  :: Tdomain
        type(surf_num), intent(in)                :: surf
        integer                                   :: ngll_if, nv, ne, gll1, gll2, nf, nfs
        integer                                   :: idx0, nes, nvs, i, j, ngll1, ngll2
        character(len=256)                        :: FunctionName ='get_surf_gll_coord'
        character(len=256)                        :: SourceFile = 'compute_surface_BC'
        character(len=700)                        :: ErrorSMS
                      
        ngll_if = 0
        allocate(Coord(0:surf%nbtot-1,0:2))
        ! FACES
        do nf=0,surf%n_faces-1
           nfs = surf%if_faces(nf)
           ngll1 = Tdomain%sFace(nfs)%ngll1
           ngll2 = Tdomain%sFace(nfs)%ngll2
           do j=1,ngll2-2
              do i=1,ngll1-2
                 idx0= Tdomain%sFace(nfs)%Iglobnum_Face(i,j)
                 Coord(ngll_if,:) = Tdomain%GlobCoord(:,idx0)
                 ngll_if = ngll_if + 1
              end do
           end do
        end do
        ! EDGES
        do ne=0,surf%n_edges-1
           nes = surf%if_edges(ne)
           ngll1 = Tdomain%sEdge(nes)%ngll 
           do i=1,ngll1-2
              idx0 = Tdomain%sEdge(nes)%Iglobnum_Edge(i)
              Coord(ngll_if,:)  = Tdomain%GlobCoord(:,idx0)
              ngll_if = ngll_if + 1
           end do
        end do
        ! VERTICES
        do nv=0,surf%n_vertices-1
           nvs = surf%if_vertices(nv)
           idx0 = Tdomain%sVertex(nvs)%Iglobnum_Vertex
           Coord(ngll_if,:)  = Tdomain%GlobCoord(:,idx0)
           ngll_if = ngll_if + 1
        end do
        ! Check
        if (ngll_if/=surf%nbtot) then
           ErrorSMS= "Incoherent surface face+edge+vert != face"
           call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
        end if

    end subroutine get_surf_gll_coord
    !-------------------------------------------------------------------------------
    !------------------------------------------------------------------------------- 
    subroutine surface_force(surf, surf_source, surf_norm, Tdomain)
        use sdomain
        use sinterface
        use ssurf
        use constants

        implicit none
        type(surf_num), intent(in)                 :: surf
        type(SurfaceParam), intent(in)             :: surf_source
        type(SurfaceT), intent(in)                 :: surf_norm
        type(domain), intent(inout)                :: Tdomain
        real(kind=8), dimension(0:2)               :: BtN, force
        real(kind=8), dimension(:,:), allocatable  :: coord
        real(kind=8)                               :: xpt, ypt, zpt
        integer                                    :: i, idx
        character(len=256)                         :: FunctionName ='surface_force'
        character(len=256)                         :: SourceFile = 'compute_surface_BC'
        character(len=700)                         :: ErrorSMS
        
        do i=0,surf%nbtot-1
           idx = surf%map(i)
           BtN = surf_norm%Surf_BtN(:,i)
           xpt = surf_norm%coord(i,0)
           ypt = surf_norm%coord(i,1)
           zpt = surf_norm%coord(i,2)
            
           force = forces_on_face(xpt, ypt, zpt, BtN, surf_source, Tdomain%TimeD%dtmin, Tdomain%TimeD%rtime)

           select case (surf_norm%domain)
              case (DM_SOLID)
                 Tdomain%sdom%champs1%Forces(idx,:) = Tdomain%sdom%champs1%Forces(idx,:) + force
              case (DM_FLUID)
              
              case default 
                 ErrorSMS= "surface force computation, only solid domain is implemented "
                 call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
           end select
        enddo

    end subroutine surface_force
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function forces_on_face(xpt, ypt, zpt, Btn, Param, dt, ctime)
        ! gives the forces on faces

        use Mathfval
        use parameters, only: Addparametricvar
        use ssurf
        use ssources

        implicit none

        real(kind=8),dimension(0:2)              :: forces_on_face
        real(kind=8), intent(in)                 :: xpt, ypt, zpt, dt, ctime
        real(kind=8), dimension(0:2), intent(in) :: Btn
        type(SurfaceParam), intent(in)           :: Param
        type(FoncValue)                          :: Sourcef
        type (source)                            :: Sour
        real(kind=8)                             :: Sigma11, Sigma22 ,Sigma33, Sigma12, Sigma13, Sigma23, midtime
        
        forces_on_face = 0.0

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
                if (Sourcef%source == 'F') then
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
         endif
         
         midtime = ctime
         if(ctime /= 0.) midtime = dt/2. + ctime

         select case(Param%wtype)
                case('R')
                    !- Ricker in time, source uniformly distributed..
                    Sour%i_time_function = 2
                    Sour%cutoff_freq     = Param%f0
                    Sour%tau_b           = Param%Rickertau
                    Sour%amplitude_factor= Param%amplitude
                    
                    forces_on_face = Param%dir*CompSource(Sour,midtime,0)
                case('G')
                    !- Gaussian in time, source uniformly distributed...
                    Sour%i_time_function = 1
                    Sour%tau_b           = Param%Rickertau
                    Sour%amplitude_factor= Param%amplitude
                    
                    forces_on_face = Param%dir*CompSource(Sour,midtime,0)

                case ('A')
                      
                      CALL ffvalue(Sourcef , (/(xpt-Param%scoord(0)), (ypt-Param%scoord(1)), (zpt-Param%scoord(2))/), midtime)

                      if ((Sourcef%dim==3).and.(Sourcef%source == 'M')) then
                           forces_on_face(0) = (Sourcef%fvalue(1)*Btn(0)+ Sourcef%fvalue(4)*Btn(1)+Sourcef%fvalue(6)*Btn(2))
                           forces_on_face(1) = (Sourcef%fvalue(4)*Btn(0)+ Sourcef%fvalue(2)*Btn(1)+Sourcef%fvalue(5)*Btn(2))
                           forces_on_face(2) = (Sourcef%fvalue(6)*Btn(0)+ Sourcef%fvalue(5)*Btn(1)+Sourcef%fvalue(3)*Btn(2))
                      elseif ((Sourcef%dim==3).and.(Sourcef%source == 'F')) then
                           forces_on_face(0) = Sourcef%fvalue(1)
                           forces_on_face(1) = Sourcef%fvalue(2)
                           forces_on_face(2) = Sourcef%fvalue(3)
                      elseif ((Sourcef%dim==1).and.(Sourcef%source == 'F')) then
                           forces_on_face = Param%dir*Sourcef%fvalue(1)
                      endif
         end select

    end function forces_on_face
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
end module surface_load

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
