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
! Vitesse des ondes directes P et S
Real(kind=8) :: Velocity_P, Velocity_S
! Coefficients : propriétés matériaux
real(kind=8) :: Coef_lambda, Coef_mu, rho
! Vitesse de l'onde plane
real(kind=8) :: PWSpeed
real(kind=8), dimension(0:2) :: Velocity_PW
! amplitude spatiale
real(kind=8) :: srcshape

contains
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine add_surface_force(Tdomain)
        use sdomain
        use constants

        implicit none
        type(domain), intent(inout)  :: Tdomain
        integer                      :: n, n1, n2, n3, n4, ns
        character(len=30 )           :: char
        character(len=256)           :: FunctionName ='add_surface_force'
        character(len=256)           :: SourceFile = 'compute_surface_BC'
        character(len=700)           :: ErrorSMS


        if (Tdomain%n_NEBC/=0) then
           !! Add Neumann boundary condition
           do ns=0,size(Tdomain%list_NEBC)-1
              n2=Tdomain%list_NEBC(ns+1)
              do n1 = 1,size(Tdomain%nsurfsource(n2)%index)
                 write(char,*) Tdomain%nsurfsource(n2)%index(n1)
                 bloc : &
                 do n = 0,size(Tdomain%sSurfaces)-1
                     if (Tdomain%sSurfaces(n)%name=="surface"//adjustl(char(:len_trim(char)))) then
                        select case (Tdomain%sSurfaces(n)%domain)
                           case(DM_SOLID)
                                call surface_force(Tdomain%sSurfaces(n)%surf_sl, Tdomain%nsurfsource(n2), &
                                                   Tdomain%sSurfaces(n),Tdomain)
                           case(DM_FLUID)
                                call surface_force(Tdomain%sSurfaces(n)%surf_fl, Tdomain%nsurfsource(n2), &
                                                   Tdomain%sSurfaces(n),Tdomain)
                           case(DM_SOLID_PML)
                                call surface_force(Tdomain%sSurfaces(n)%surf_spml, Tdomain%nsurfsource(n2), &
                                                   Tdomain%sSurfaces(n),Tdomain)
                           case(DM_FLUID_PML)
                                call surface_force(Tdomain%sSurfaces(n)%surf_fpml, Tdomain%nsurfsource(n2), &
                                                   Tdomain%sSurfaces(n),Tdomain)
                        end select
                        exit bloc
                     endif
                 enddo bloc 
              enddo
           enddo
        endif

        if (Tdomain%n_PWBC /= 0) then

           do ns=0,size(Tdomain%list_PWBC)-1
              n2=Tdomain%list_PWBC(ns+1)
              do n1 = 1,size(Tdomain%nsurfsource(n2)%index)  
                 write(char,*) Tdomain%nsurfsource(n2)%index(n1)
                 blocPW : &
                 do n = 0,size(Tdomain%sSurfaces)-1
                    if (Tdomain%sSurfaces(n)%name=="surface"//adjustl(char(:len_trim(char)))) then
                        select case (Tdomain%sSurfaces(n)%domain)
                               case(DM_SOLID)
                                   call surface_force(Tdomain%sSurfaces(n)%surf_sl, Tdomain%nsurfsource(n2), &
                                                      Tdomain%sSurfaces(n),Tdomain,Tdomain%sdom%champs0%Veloc)
                               case(DM_FLUID)
                                    call surface_force(Tdomain%sSurfaces(n)%surf_fl, Tdomain%nsurfsource(n2), &
                                                       Tdomain%sSurfaces(n),Tdomain)
                        endselect
                        exit blocPW
                     endif
                 enddo blocPW
             enddo
           enddo
        endif

        if (Tdomain%n_FTBC /= 0) then
           ErrorSMS= "Sorry the fault problem is not yet implemented"
           call ErrorMessage(ErrorSMS,FunctionName,SourceFile)

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
    subroutine surface_force(surf, surf_source, surf_norm, Tdomain, veloc_field)
        use sdomain
        use sinterface
        use ssurf
        use constants

        implicit none
        type(surf_num), intent(in)                         :: surf
        type(SurfaceParam), intent(in)                     :: surf_source
        type(SurfaceT), intent(in)                         :: surf_norm
        type(domain), intent(inout)                        :: Tdomain
        real(kind=8), dimension(:,:),optional,intent(in)   :: veloc_field
        real(kind=8), dimension(0:2)                       :: BtN, force, veloc
        real(kind=8), dimension(:,:), allocatable          :: coord
        real(kind=8)                                       :: xpt, ypt, zpt
        integer                                            :: i, idx
        character(len=256)                                 :: FunctionName ='surface_force'
        character(len=256)                                 :: SourceFile = 'compute_surface_BC'
        character(len=700)                                 :: ErrorSMS
       
       
        Velocity_P = surf_norm%Elastic%Pspeed
        Velocity_S = surf_norm%Elastic%Sspeed
        if ((surf_source%what_bc=='PW').and.(present(veloc_field))) then
           PWSpeed= surf_norm%Elastic%PWspeed     
           Velocity_PW= surf_norm%Elastic%PWspeed*surf_source%dir
        endif
        Coef_lambda= surf_norm%Elastic%lambda
        Coef_mu    = surf_norm%Elastic%mu
        rho        = surf_norm%Elastic%density        
        do i=0,surf%nbtot-1
           idx      = surf%map(i)
           BtN      = surf_norm%Surf_BtN(:,i)
           xpt      = surf_norm%coord(i,0)
           ypt      = surf_norm%coord(i,1)
           zpt      = surf_norm%coord(i,2)
           srcshape = surf_norm%source(i)        
           veloc    = 0.d0
           if ((surf_source%what_bc=='PW').and.(present(veloc_field))) then 
               veloc=veloc_field(idx,:)
               force = forces_on_face(xpt, ypt, zpt, BtN, surf_source, Tdomain%TimeD%dtmin, Tdomain%TimeD%rtime, veloc)
           else
               force = forces_on_face(xpt, ypt, zpt, BtN, surf_source, Tdomain%TimeD%dtmin, Tdomain%TimeD%rtime)
           endif

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
    function forces_on_face(xpt, ypt, zpt, Btn, Param, dt, ctime, veloc)
        ! gives the forces on faces

        use Mathfval
        use parameters, only: Addparametricvar
        use ssurf
        use ssources
        use Surface_prbl_type

        implicit none

        real(kind=8),dimension(0:2)                      :: forces_on_face
        real(kind=8), intent(in)                         :: xpt, ypt, zpt, dt, ctime
        real(kind=8), dimension(0:2),optional,intent(in) :: Btn, veloc
        type(SurfaceParam), intent(in)                   :: Param
        type(FoncValue)                                  :: Sourcef
        type (source)                                    :: Sour
        real(kind=8)                                     :: midtime
        real(kind=8), dimension(0:2)                     :: coord
        character(len=256)                               :: FunctionName ='surface_force'
        character(len=256)                               :: SourceFile = 'compute_surface_BC'
        character(len=700)                               :: ErrorSMS

        
        forces_on_face = 0.d0

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
         if ((ctime /= 0.).and.(Param%what_bc /= 'PW')) midtime = dt/2. + ctime
         
         !! Nouvelle coordonnées par rapport au point de référence
         coord =  (/(xpt-Param%scoord(0)), (ypt-Param%scoord(1)), (zpt-Param%scoord(2))/)

         if (Param%what_bc == 'NE') then
            call Neumanforce(Param,srcshape,Sourcef,coord,Btn,midtime,forces_on_face)
         endif
         if ((Param%what_bc == 'PW').and.(present(veloc))) then
            call  PlaneWane_diffracted(Btn, coord, Veloc, ctime, dt, Param, Sourcef, forces_on_face)
         elseif ((Param%what_bc == 'PW').and.(.not.present(veloc))) then
            ErrorSMS= " Only the solid domain is needed for plane wave problem "
            call ErrorMessage(ErrorSMS,FunctionName,SourceFile)
         endif

    end function forces_on_face
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine PlaneWane_diffracted(Btn, coord, Veloc, time, dt, Param, Sourcef, force)
        
        use Surface_prbl_type
        use ssurf
        use Mathfval

        implicit none
        type(SurfaceParam),            intent(in)   :: Param
        real(kind=8), dimension(0:2), intent(inout) :: force
        type(FoncValue),               intent(in)   :: Sourcef
        real(kind=8), dimension(0:2), intent(in   ) :: Btn, Veloc, coord
        real(kind=8),                 intent(in   ) :: time, dt
        real(kind=8), dimension(0:2)                :: Nm, V, Traction_i, vel_i, Dirtang, Traction_d, vel_i_n, P, D
        real(kind=8)                                :: nn, VV
        real(kind=8)                                :: S11, S22 ,S33, S12, S13, S23
        

        !call PlaneWanedispl(Param,Sourcef,coord,Btn,time,vel_i_n)
        
        call PlaneWavedispl(Param,Sourcef,coord,Btn,time,PWSpeed,vel_i,'wave')
        
        ! Principales directions
        ! D = Param%Kdir
        ! P = Param%dir

        S11 = (Coef_lambda+2.*Coef_mu)*vel_i(0)+Coef_lambda*(vel_i(1)+vel_i(2))
        S22 = (Coef_lambda+2.*Coef_mu)*vel_i(1)+Coef_lambda*(vel_i(0)+vel_i(2))
        S33 = (Coef_lambda+2.*Coef_mu)*vel_i(2)+Coef_lambda*(vel_i(0)+vel_i(1))
        S12 = Coef_mu*(vel_i(0)+vel_i(1))
        S13 = Coef_mu*(vel_i(0)+vel_i(2))
        S23 = Coef_mu*(vel_i(1)+vel_i(2))
       
        Traction_i(0) = -( S11*Btn(0)+S12*Btn(1)+S13*Btn(2) )
        Traction_i(1) = -( S12*Btn(0)+S22*Btn(1)+S23*Btn(2) ) 
        Traction_i(2) = -( S13*Btn(0)+S23*Btn(1)+S33*Btn(2) )
        
        Nm = Btn
        nn = sqrt(Btn(0)**2 + Btn(1)**2 + Btn(2)**2)
        if (nn .gt. 0) Nm = Nm/nn
        V  = Veloc - Velocity_PW
        VV = Nm(0)*V(0) + Nm(1)*V(1) + Nm(2)*V(2)

        ! contraintes de traction due à la diffraction de l'onde plane
        Traction_d =  -Velocity_S*rho*V + (Velocity_P-Velocity_S)*rho*VV*Nm
        ! force de surface résultante
        force = Traction_i + Traction_d

    end subroutine PlaneWane_diffracted
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
