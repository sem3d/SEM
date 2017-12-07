!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module surface_load
    use constants, only : fpp
    implicit none
    ! Vitesse des ondes directes P et S
    real(fpp) :: Velocity_P, Velocity_S
    ! Coefficients : propriétés matériaux
    real(fpp) :: Coef_lambda, Coef_mu, rho
    ! Vitesse de l'onde plane
    real(fpp) :: PWSpeed
    real(fpp), dimension(0:2) :: Velocity_PW
    ! amplitude spatiale
    real(fpp) :: srcshape


contains
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine add_surface_force(Tdomain, f0, f1)
        use sdomain
        use constants

        implicit none
        type(domain), intent(inout)  :: Tdomain
        integer, intent(in) :: f0, f1
        integer                      :: n, n1, n2, ns
        character(len=30 )           :: char

        if (Tdomain%n_NEBC/=0) then
            !! Add Neumann boundary condition
            do ns=1,size(Tdomain%list_NEBC)
                n2=Tdomain%list_NEBC(ns)
                do n1 = 1,size(Tdomain%nsurfsource(n2)%index)
                    write(char,*) Tdomain%nsurfsource(n2)%index(n1)
                    do n = 0,size(Tdomain%sSurfaces)-1
                        if (Tdomain%sSurfaces(n)%name=="surface"//adjustl(char(:len_trim(char)))) then
                            select case (Tdomain%sSurfaces(n)%domain)
                            case(DM_SOLID)
                                call surface_force(f0, f1, Tdomain%sSurfaces(n)%surf_sl, &
                                    Tdomain%nsurfsource(n2), Tdomain%sSurfaces(n),Tdomain)
                            case(DM_FLUID)
                                call surface_force(f0, f1, Tdomain%sSurfaces(n)%surf_fl, &
                                    Tdomain%nsurfsource(n2), Tdomain%sSurfaces(n),Tdomain)
                            case(DM_SOLID_PML,DM_FLUID_PML)
                                stop "Neumann condition should not be used on PML surface"
                            end select
                            exit
                        endif
                    enddo
                enddo
            enddo
        endif

        Tdomain%sdom%PlaneW%Exist = .false.
        if (Tdomain%n_PWBC /= 0) then
            do ns=1,size(Tdomain%list_PWBC)
                n2=Tdomain%list_PWBC(ns)
                do n1 =1,size(Tdomain%nsurfsource(n2)%index)
                    write(char,*) Tdomain%nsurfsource(n2)%index(n1)

                    do n = 0,size(Tdomain%sSurfaces)-1
                        if (Tdomain%sSurfaces(n)%name=="surface"//adjustl(char(:len_trim(char)))) then
                            select case (Tdomain%sSurfaces(n)%domain)
                            case(DM_SOLID)
                                Tdomain%sdom%PlaneW%Exist = .true.
                                call surface_force(f0, f1, Tdomain%sSurfaces(n)%surf_sl, &
                                    Tdomain%nsurfsource(n2), Tdomain%sSurfaces(n), &
                                    Tdomain, Tdomain%sdom%champs(f0)%Veloc)
                            case(DM_SOLID_PML,DM_FLUID_PML,DM_FLUID)
                                stop "Plane wave problem is implemented only for solid domains"
                            endselect
                            exit
                        end if
                    end do
                enddo
            enddo
        endif

        if (Tdomain%n_FTBC /= 0) then
            stop "Sorry the fault problem is not yet implemented"
        endif

    end subroutine add_surface_force
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine get_surf_gll_coord(surf, Tdomain, Coord)
        use sdomain
        use sinterface
        use constants

        implicit none
        real(fpp), dimension(:,:), allocatable :: coord
        type(domain), intent(in)                  :: Tdomain
        type(surf_num), intent(in)                :: surf
        integer                                   :: ngll_if, nv, ne,  nf, nfs
        integer                                   :: idx0, nes, nvs, i, j, ngll

        ngll_if = 0
        allocate(Coord(0:surf%nbtot-1,0:2))
        ! FACES
        do nf=0,surf%n_faces-1
            nfs = surf%if_faces(nf)
            ngll = Tdomain%sFace(nfs)%ngll

            do j=1,ngll-2
                do i=1,ngll-2
                    idx0= Tdomain%sFace(nfs)%Iglobnum_Face(i,j)
                    Coord(ngll_if,:) = Tdomain%GlobCoord(:,idx0)
                    ngll_if = ngll_if + 1
                end do
            end do
        end do
        ! EDGES
        do ne=0,surf%n_edges-1
            nes = surf%if_edges(ne)
            ngll = Tdomain%sEdge(nes)%ngll
            do i=1,ngll-2
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
            stop "Incoherent surface face+edge+vert != face"
        end if

    end subroutine get_surf_gll_coord
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine surface_force(f0, f1, surf, surf_source, surf_norm, Tdomain, veloc_field)
        use sdomain
        use sinterface
        use ssurf
        use constants

        implicit none
        integer, intent(in) :: f0, f1
        type(surf_num),     intent(in)   :: surf
        type(SurfaceParam), intent(in)   :: surf_source
        type(SurfaceT),     intent(in)   :: surf_norm
        type(domain),       intent(inout):: Tdomain
        real(fpp), dimension(:,:), optional,intent(in)   :: veloc_field
        real(fpp), dimension(0:2)               :: BtN, force, veloc
        integer                                    :: i, idx

        Velocity_P = surf_norm%Elastic%Pspeed
        Velocity_S = surf_norm%Elastic%Sspeed
        if ((surf_source%what_bc=='PW').and.(present(veloc_field))) then
            PWSpeed = surf_norm%Elastic%PWspeed
            Velocity_PW = surf_norm%Elastic%PWspeed*surf_source%dir
        endif
        Coef_lambda= surf_norm%Elastic%lambda
        Coef_mu    = surf_norm%Elastic%mu
        rho        = surf_norm%Elastic%density
        do i=0,surf%nbtot-1
            idx      = surf%map(i)
            BtN      = surf_norm%Surf_BtN(:,i)
            srcshape = surf_norm%source(i)
            veloc    = 0.d0
            if ((surf_source%what_bc=='PW').and.(present(veloc_field))) then
                veloc=veloc_field(idx,:)
                call forces_on_face(surf_norm%coord(i,:),BtN,surf_source,Tdomain%TimeD%dtmin,Tdomain%TimeD%rtime,force,veloc)
            else
                call forces_on_face(surf_norm%coord(i,:),BtN,surf_source,Tdomain%TimeD%dtmin,Tdomain%TimeD%rtime,force)
            endif

            select case (surf_norm%domain)
            case (DM_SOLID)
                Tdomain%sdom%champs(f1)%Forces(idx,:) = Tdomain%sdom%champs(f1)%Forces(idx,:) + force
            case (DM_FLUID)

            case default
                stop "surface force computation, only solid domain is implemented "
            end select
        enddo

    end subroutine surface_force
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine forces_on_face(gllcoord, Btn, Param, dt, ctime, force, veloc)
        ! gives the forces on faces

        use Mathfval
        use ssurf
        use ssources
        use Surface_prbl_type

        implicit none
        real(fpp)                         ,intent(in ):: ctime, dt
        real(fpp),dimension(0:2)          ,intent(out):: force
        real(fpp), dimension(0:2),         intent(in) :: Btn, gllcoord
        real(fpp), dimension(0:2),optional,intent(in) :: veloc
        type(SurfaceParam),                   intent(in) :: Param
        real(fpp)                     :: midtime
        real(fpp), dimension(0:2)     :: coord

        force = 0.d0
        midtime = ctime
        if ((ctime/=0.).and.(Param%what_bc/='PW')) midtime = dt/2. + ctime

        !! Nouvelle coordonnées par rapport au point de référence
        coord = gllcoord - Param%scoord

        if (Param%what_bc == 'NE') then
            call Neumanforce(Param,srcshape,coord,Btn,midtime,force)
        endif

        if ((Param%what_bc == 'PW').and.(present(veloc))) then
            call PlaneWane_Reflected(Btn, coord, Veloc, ctime, dt, Param, force)

        elseif ((Param%what_bc == 'PW').and.(.not.present(veloc))) then
            stop " Only the solid domain is needed for plane wave problem "
        endif

    end subroutine forces_on_face
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine PlaneWane_IncidenteTraction(Btn, coord, Veloc, time, dt, Param, Traction_i)

        use Surface_prbl_type
        use ssurf
        use Mathfval

        implicit none
        type(SurfaceParam),            intent(in)   :: Param
        real(fpp), dimension(0:2), intent(inout) :: Traction_i
        real(fpp), dimension(0:2), intent(in)    :: Btn, Veloc, coord
        real(fpp),                 intent(in)    :: time, dt
        real(fpp), dimension(0:2)                :: vel_i, voloc_0
        real(fpp)                                :: S11, S22 ,S33, S12, S13, S23

        call PlaneWaveDerive(Param,coord,time,PWSpeed,vel_i, voloc_0)

        S11 = (Coef_lambda+2.d0*Coef_mu)*vel_i(0)*Param%dir(0)+Coef_lambda*(vel_i(1)*Param%dir(1)+vel_i(2)*Param%dir(2))
        S22 = (Coef_lambda+2.d0*Coef_mu)*vel_i(1)*Param%dir(1)+Coef_lambda*(vel_i(0)*Param%dir(0)+vel_i(2)*Param%dir(2))
        S33 = (Coef_lambda+2.d0*Coef_mu)*vel_i(2)*Param%dir(2)+Coef_lambda*(vel_i(0)*Param%dir(0)+vel_i(1)*Param%dir(1))
        S12 = Coef_mu*(vel_i(0)*Param%dir(1)+vel_i(1)*Param%dir(2))
        S13 = Coef_mu*(vel_i(0)*Param%dir(2)+vel_i(2)*Param%dir(0))
        S23 = Coef_mu*(vel_i(1)*Param%dir(2)+vel_i(2)*Param%dir(1))

        Traction_i(0) = ( S11*Btn(0)+S12*Btn(1)+S13*Btn(2) )
        Traction_i(1) = ( S12*Btn(0)+S22*Btn(1)+S23*Btn(2) )
        Traction_i(2) = ( S13*Btn(0)+S23*Btn(1)+S33*Btn(2) )

    end subroutine PlaneWane_IncidenteTraction
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine PlaneWane_diffracted(Btn, coord, Veloc, time, dt, Param, force)

        use Surface_prbl_type
        use ssurf
        use Mathfval

        implicit none
        type(SurfaceParam),            intent(in)   :: Param
        real(fpp), dimension(0:2), intent(inout) :: force
        real(fpp), dimension(0:2), intent(in)    :: Btn, Veloc, coord
        real(fpp),                 intent(in)    :: time, dt
        real(fpp), dimension(0:2)                :: Nm, V, voloc_0
        real(fpp)                                :: VV

        Nm = Btn
        V  = Veloc - voloc_0
        VV = Nm(0)*V(0) + Nm(1)*V(1) + Nm(2)*V(2)

        force =  Velocity_S*rho*V + (Velocity_P-Velocity_S)*rho*VV*Nm

    end subroutine PlaneWane_diffracted
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine PlaneWane_Reflected(Btn, coord, Veloc, time, dt, Param, force)
        use Surface_prbl_type
        use ssurf
        use Mathfval
        implicit none
        type(SurfaceParam),            intent(in)   :: Param
        real(fpp), dimension(0:2), intent(inout) :: force
        real(fpp), dimension(0:2), intent(in)    :: Btn, Veloc, coord
        real(fpp),                 intent(in)    :: time, dt
        real(fpp), dimension(0:2)                :: Traction_i

        call PlaneWane_IncidenteTraction(Btn, coord, Veloc, time, dt, Param, Traction_i)

        force = - Traction_i

        !write(*,*) Traction_i
        !if (Param%wave_type==1) then
        !      force = - Traction_i
        !elseif (Param%wave_type==2) then

        !elseif (Param%wave_type==3) then

        !elseif (Param%wave_type==4) then

        !else

        !endif

    end subroutine PlaneWane_Reflected
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
