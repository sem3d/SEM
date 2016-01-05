!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file PlaneW.f90
!! \brief
!!
!<

module splanew
    implicit none
    ! #####################################################################################
    ! #####################################################################################

    type Face_PlaneW

       integer :: ngll1, ngll2, mat_index, dir, Face_UP, Face_DOWN, Orient
       integer, dimension (0:3) :: Near_Edges, Near_Vertices, Orient_Edges
       real, dimension (:,:), pointer :: ds, MassMat_Up, MassMat_Down
       real, dimension (:,:,:), pointer :: normal, Btn, Coord_nodes, Forces_Up, Forces_Down

    end type Face_PlaneW


    type Edge_PlaneW

       integer :: ngll, mat_index, dir, Edge_UP, Edge_DOWN, Orient
       real, dimension (:), pointer :: MassMat_Up, MassMat_Down
       real, dimension (:,:), pointer :: Btn, Coord_nodes, Forces_Up, Forces_Down

    end type Edge_PlaneW


    type Vertex_PlaneW

       integer :: Vertex_UP, Vertex_DOWN, mat_index
       real :: MassMat_Up, MassMat_Down
       real, dimension (0:2) :: Btn, Coord_nodes, Forces_Up, Forces_Down

    end type Vertex_PlaneW


    type Param_PlaneW

       real :: Mu, Lambda, Kappa, speed, lx, ly, lz, xs, ys, zs, f0
       character (len=1) :: wtype

    end type Param_PlaneW


    type PlaneW

       integer :: n_faces, n_edges, n_vertices
       type(Face_PlaneW), dimension (:), pointer :: pFace
       type(Edge_PlaneW), dimension (:), pointer :: pEdge
       type(Vertex_PlaneW), dimension (:), pointer :: pVertex
       type(Param_PlaneW) :: pParam

    end type PlaneW

contains

    !  ########################################################################################
    subroutine compute_pforces_on_face (Face, Param, Vfree, dt, ctime)
        implicit none

        type (Face_PlaneW), intent (INOUT) :: Face
        type (Param_PlaneW), intent (IN) :: Param
        real, dimension (1:Face%ngll1-2,1:Face%ngll2-2,0:2), intent (IN) :: Vfree
        real, intent (IN) :: dt, ctime

        integer :: i,j, ngll1, ngll2
        real, dimension (:,:,:), allocatable :: Traction_i, vel_i
        real :: xpt, ypt, zpt, velixn, veliyn, velizn, Sigma11, Sigma22 ,Sigma33, Sigma12, Sigma13, Sigma23


        ngll1 = Face%ngll1
        ngll2 = Face%ngll2

        allocate ( Traction_i(1:ngll1-2,1:ngll2-2,0:2) )
        allocate ( vel_i(1:ngll1-2,1:ngll2-2,0:2) )
        Traction_i = 0
        vel_i = 0

        do j = 1,ngll2-2
            do i = 1,ngll1-2

                xpt = Face%Coord_nodes(i,j,0)
                ypt = Face%Coord_nodes(i,j,1)
                zpt = Face%Coord_nodes(i,j,2)

                call Inc_VelocPW (Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
                call Inc_VelocPW (Param,xpt,ypt,zpt,ctime+dt,vel_i(i,j,0),vel_i(i,j,1),vel_i(i,j,2))

                vel_i(i,j,0) = 0.5*(velixn+vel_i(i,j,0))
                vel_i(i,j,1) = 0.5*(veliyn+vel_i(i,j,1))
                vel_i(i,j,2) = 0.5*(velizn+vel_i(i,j,2))

                Sigma11 = (Param%Lambda+2*Param%Mu)*Param%lx*vel_i(i,j,0)+Param%Lambda*(Param%ly*vel_i(i,j,1)+Param%lz*vel_i(i,j,2))
                Sigma22 = (Param%Lambda+2*Param%Mu)*Param%ly*vel_i(i,j,1)+Param%Lambda*(Param%lx*vel_i(i,j,0)+Param%lz*vel_i(i,j,2))
                Sigma33 = (Param%Lambda+2*Param%Mu)*Param%lz*vel_i(i,j,2)+Param%Lambda*(Param%lx*vel_i(i,j,0)+Param%ly*vel_i(i,j,1))
                Sigma12 = Param%Mu*(Param%ly*vel_i(i,j,0)+Param%lx*vel_i(i,j,1))
                Sigma13 = Param%Mu*(Param%lz*vel_i(i,j,0)+Param%lx*vel_i(i,j,2))
                Sigma23 = Param%Mu*(Param%lz*vel_i(i,j,1)+Param%ly*vel_i(i,j,2))

                Traction_i(i,j,0) = -( Sigma11*Face%Btn(i,j,0) + Sigma12*Face%Btn(i,j,1) + Sigma13*Face%Btn(i,j,2) ) / Param%speed
                Traction_i(i,j,1) = -( Sigma12*Face%Btn(i,j,0) + Sigma22*Face%Btn(i,j,1) + Sigma23*Face%Btn(i,j,2) ) / Param%speed
                Traction_i(i,j,2) = -( Sigma13*Face%Btn(i,j,0) + Sigma23*Face%Btn(i,j,1) + Sigma33*Face%Btn(i,j,2) ) / Param%speed

                if ( Face%Orient == 0 ) then
                    Face%Forces_Up(i,j,0:2) = ( Vfree(i,j,0:2) -  vel_i(i,j,0:2) + dt*Face%MassMat_Down(i,j)*Traction_i(i,j,0:2) ) / &
                        ( dt * (Face%MassMat_Up(i,j)+Face%MassMat_Down(i,j)) )
                    Face%Forces_Down(i,j,0:2) = Face%Forces_Up(i,j,0:2) - Traction_i(i,j,0:2)
                else if ( Face%Orient == 1 ) then
                    Face%Forces_Up(i,j,0:2) = ( Vfree(i,j,0:2) -  vel_i(i,j,0:2) + dt*Face%MassMat_Down(ngll1-1-i,j)*Traction_i(i,j,0:2) ) / &
                        ( dt * (Face%MassMat_Up(i,j)+Face%MassMat_Down(ngll1-1-i,j)) )
                    Face%Forces_Down(ngll1-1-i,j,0:2) = Face%Forces_Up(i,j,0:2) - Traction_i(i,j,0:2)
                else if ( Face%Orient == 2 ) then
                    Face%Forces_Up(i,j,0:2) = ( Vfree(i,j,0:2) -  vel_i(i,j,0:2) + dt*Face%MassMat_Down(i,ngll2-1-j)*Traction_i(i,j,0:2) ) / &
                        ( dt * (Face%MassMat_Up(i,j)+Face%MassMat_Down(i,ngll2-1-j)) )
                    Face%Forces_Down(i,ngll2-1-j,0:2) = Face%Forces_Up(i,j,0:2) - Traction_i(i,j,0:2)
                else if ( Face%Orient == 3 ) then
                    Face%Forces_Up(i,j,0:2) = ( Vfree(i,j,0:2) -  vel_i(i,j,0:2) + dt*Face%MassMat_Down(ngll1-1-i,ngll2-1-j)&
                        *Traction_i(i,j,0:2) ) /( dt * (Face%MassMat_Up(i,j)+Face%MassMat_Down(ngll1-1-i,ngll2-1-j)) )
                    Face%Forces_Down(ngll1-1-i,ngll2-1-j,0:2) = Face%Forces_Up(i,j,0:2) - Traction_i(i,j,0:2)
                else if ( Face%Orient == 4 ) then
                    Face%Forces_Up(i,j,0:2) = ( Vfree(i,j,0:2) -  vel_i(i,j,0:2) + dt*Face%MassMat_Down(j,i)*Traction_i(i,j,0:2) ) / &
                        ( dt * (Face%MassMat_Up(i,j)+Face%MassMat_Down(j,i)) )
                    Face%Forces_Down(j,i,0:2) = Face%Forces_Up(i,j,0:2) - Traction_i(i,j,0:2)
                else if ( Face%Orient == 5 ) then
                    Face%Forces_Up(i,j,0:2) = ( Vfree(i,j,0:2) -  vel_i(i,j,0:2) + dt*Face%MassMat_Down(ngll1-1-j,i)*Traction_i(i,j,0:2) ) / &
                        ( dt * (Face%MassMat_Up(i,j)+Face%MassMat_Down(ngll1-1-j,i)) )
                    Face%Forces_Down(ngll1-1-j,i,0:2) = Face%Forces_Up(i,j,0:2) - Traction_i(i,j,0:2)
                else if ( Face%Orient == 6 ) then
                    Face%Forces_Up(i,j,0:2) = ( Vfree(i,j,0:2) -  vel_i(i,j,0:2) + dt*Face%MassMat_Down(j,ngll2-1-i)*Traction_i(i,j,0:2) ) / &
                        ( dt * (Face%MassMat_Up(i,j)+Face%MassMat_Down(j,ngll2-1-i)) )
                    Face%Forces_Down(j,ngll2-1-i,0:2) = Face%Forces_Up(i,j,0:2) - Traction_i(i,j,0:2)
                else if ( Face%Orient == 7 ) then
                    Face%Forces_Up(i,j,0:2) = ( Vfree(i,j,0:2) -  vel_i(i,j,0:2) + dt*Face%MassMat_Down(ngll1-1-j,ngll2-1-i) &
                        *Traction_i(i,j,0:2) ) /( dt * (Face%MassMat_Up(i,j)+Face%MassMat_Down(ngll1-1-j,ngll2-1-i)) )
                    Face%Forces_Down(ngll1-1-j,ngll2-1-i,0:2) = Face%Forces_Up(i,j,0:2) - Traction_i(i,j,0:2)
                endif
            enddo
        enddo
        deallocate (vel_i, Traction_i)

    end subroutine compute_pforces_on_face

    !  ########################################################################################
    subroutine compute_pforces_on_edge (Edge, Param, Vfree, dt, ctime)
        implicit none

        type (Edge_PlaneW), intent (INOUT) :: Edge
        type (Param_PlaneW), intent (IN) :: Param
        real, dimension (1:Edge%ngll-2,0:2), intent (IN) :: Vfree
        real, intent (IN) :: dt, ctime

        integer :: i, ngll
        real :: xpt, ypt, zpt, velixn, veliyn, velizn, Sigma11, Sigma22 ,Sigma33, Sigma12, Sigma13, Sigma23
        real, dimension(:,:), allocatable :: Traction_i, vel_i


        ngll = Edge%ngll
        allocate ( Traction_i(1:ngll-2,0:2) )
        allocate ( vel_i(1:ngll-2,0:2) )
        vel_i = 0.

        do i = 1,ngll-2

            xpt = Edge%Coord_nodes(i,0)
            ypt = Edge%Coord_nodes(i,1)
            zpt = Edge%Coord_nodes(i,2)

            call Inc_VelocPW (Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
            call Inc_VelocPW (Param,xpt,ypt,zpt,ctime+dt,vel_i(i,0),vel_i(i,1),vel_i(i,2))

            vel_i(i,0) = 0.5*(velixn+vel_i(i,0))
            vel_i(i,1) = 0.5*(veliyn+vel_i(i,1))
            vel_i(i,2) = 0.5*(velizn+vel_i(i,2))

            Sigma11 = (Param%Lambda+2*Param%Mu)*Param%lx*vel_i(i,0)+Param%Lambda*(Param%ly*vel_i(i,1)+Param%lz*vel_i(i,2))
            Sigma22 = (Param%Lambda+2*Param%Mu)*Param%ly*vel_i(i,1)+Param%Lambda*(Param%lx*vel_i(i,0)+Param%lz*vel_i(i,2))
            Sigma33 = (Param%Lambda+2*Param%Mu)*Param%lz*vel_i(i,2)+Param%Lambda*(Param%lx*vel_i(i,0)+Param%ly*vel_i(i,1))
            Sigma12 = Param%Mu*(Param%ly*vel_i(i,0)+Param%lx*vel_i(i,1))
            Sigma13 = Param%Mu*(Param%lz*vel_i(i,0)+Param%lx*vel_i(i,2))
            Sigma23 = Param%Mu*(Param%lz*vel_i(i,1)+Param%ly*vel_i(i,2))

            Traction_i(i,0) = -( Sigma11*Edge%Btn(i,0) + Sigma12*Edge%Btn(i,1) + Sigma13*Edge%Btn(i,2) ) / Param%speed
            Traction_i(i,1) = -( Sigma12*Edge%Btn(i,0) + Sigma22*Edge%Btn(i,1) + Sigma23*Edge%Btn(i,2) ) / Param%speed
            Traction_i(i,2) = -( Sigma13*Edge%Btn(i,0) + Sigma23*Edge%Btn(i,1) + Sigma33*Edge%Btn(i,2) ) / Param%speed

            if ( Edge%Orient == 0 ) then
                Edge%Forces_Up(i,0:2) = ( Vfree(i,0:2) -  vel_i(i,0:2) + dt*Edge%MassMat_Down(i)*Traction_i(i,0:2) ) / &
                    ( dt * (Edge%MassMat_Up(i)+Edge%MassMat_Down(i)) )
                Edge%Forces_Down(i,0:2) = Edge%Forces_Up(i,0:2) - Traction_i(i,0:2)
            else
                Edge%Forces_Up(i,0:2) = ( Vfree(i,0:2) -  vel_i(i,0:2) + dt*Edge%MassMat_Down(ngll-1-i)*Traction_i(i,0:2) ) / &
                    ( dt * (Edge%MassMat_Up(i)+Edge%MassMat_Down(ngll-1-i)) )
                Edge%Forces_Down(ngll-1-i,0:2) = Edge%Forces_Up(i,0:2) - Traction_i(i,0:2)
            endif

        enddo

        deallocate (vel_i, Traction_i)


    end subroutine compute_pforces_on_edge

    !  ########################################################################################
    subroutine compute_pforces_on_vertex (Vertex, Param, Vfree, dt, ctime)
        implicit none

        type (Vertex_PlaneW), intent (INOUT) :: Vertex
        type (Param_PlaneW), intent (IN) :: Param
        real, dimension (0:2), intent (IN) :: Vfree
        real, intent (IN) :: dt, ctime

        real :: xpt, ypt, zpt, velixn, veliyn, velizn, Sigma11, Sigma22 ,Sigma33, Sigma12, Sigma13, Sigma23
        real, dimension(0:2) :: Traction_i, vel_i

        vel_i = 0.

        xpt = Vertex%Coord_nodes(0)
        ypt = Vertex%Coord_nodes(1)
        zpt = Vertex%Coord_nodes(2)

        call Inc_VelocPW (Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
        call Inc_VelocPW (Param,xpt,ypt,zpt,ctime+dt,vel_i(0),vel_i(1),vel_i(2))

        vel_i(0) = 0.5*(velixn+vel_i(0))
        vel_i(1) = 0.5*(veliyn+vel_i(1))
        vel_i(2) = 0.5*(velizn+vel_i(2))

        Sigma11 = (Param%Lambda+2*Param%Mu)*Param%lx*vel_i(0)+Param%Lambda*(Param%ly*vel_i(1)+Param%lz*vel_i(2))
        Sigma22 = (Param%Lambda+2*Param%Mu)*Param%ly*vel_i(1)+Param%Lambda*(Param%lx*vel_i(0)+Param%lz*vel_i(2))
        Sigma33 = (Param%Lambda+2*Param%Mu)*Param%lz*vel_i(2)+Param%Lambda*(Param%lx*vel_i(0)+Param%ly*vel_i(1))
        Sigma12 = Param%Mu*(Param%ly*vel_i(0)+Param%lx*vel_i(1))
        Sigma13 = Param%Mu*(Param%lz*vel_i(0)+Param%lx*vel_i(2))
        Sigma23 = Param%Mu*(Param%lz*vel_i(1)+Param%ly*vel_i(2))

        Traction_i(0) = -( Sigma11*Vertex%Btn(0) + Sigma12*Vertex%Btn(1) + Sigma13*Vertex%Btn(2) ) / Param%speed
        Traction_i(1) = -( Sigma12*Vertex%Btn(0) + Sigma22*Vertex%Btn(1) + Sigma23*Vertex%Btn(2) ) / Param%speed
        Traction_i(2) = -( Sigma13*Vertex%Btn(0) + Sigma23*Vertex%Btn(1) + Sigma33*Vertex%Btn(2) ) / Param%speed

        !if (.not. (Traction_i(0)==0 .and. Traction_i(1)==0 .and. Traction_i(2)==0 .and. &
        !           vel_i(0)==0 .and. vel_i(1)==0 .and. vel_i(2)==0 ) ) then
        !     print*,'HELLOV',Traction_i(0),Traction_i(1),Traction_i(2),vel_i(0),vel_i(1),vel_i(2)
        !endif

        Vertex%Forces_Up(0:2) = ( Vfree(0:2) -  vel_i(0:2) + dt*Vertex%MassMat_Down*Traction_i(0:2) ) / &
            ( dt * (Vertex%MassMat_Up+Vertex%MassMat_Down) )

        Vertex%Forces_Down(0:2) = Vertex%Forces_Up(0:2) - Traction_i(0:2)

    end subroutine compute_pforces_on_vertex

    ! #########################################################################################
    subroutine Inc_VelocPW (Param,x,y,z,ctime,velx,vely,velz)

        implicit none

        type (Param_PlaneW), intent (IN) :: Param
        real, intent(IN) :: x,z,y,ctime
        real, intent(INOUT) :: velx,vely,velz

        real :: lx,ly,lz,phasex,phasey,phasez,f0,force
        real :: vel,xsour,ysour,zsour,xpt,ypt,zpt,arg,pi

        pi = Acos(-1.)
        xsour = Param%xs
        ysour = Param%ys
        zsour = Param%zs
        f0 = Param%f0
        lx = Param%lx
        ly = Param%ly
        lz = Param%lz
        vel = Param%speed
        xpt = x-vel*ctime*lx
        ypt = y-vel*ctime*ly
        zpt = z-vel*ctime*lz

        phasex = -(xpt-xsour)*lx/vel
        phasey = -(ypt-ysour)*ly/vel
        phasez = -(zpt-zsour)*lz/vel

        arg = pi*f0*(phasex + phasey + phasez)

        force = (1.-2.*arg**2)*exp(-arg**2)                   ! Plane Wave in Velocity
        !force = exp(-arg*arg)*(-6*arg+4*arg*arg*arg)*pi*f0   ! Plane Wave in Displacement
        !force = 0.

        if ( Param%wtype == 'S' ) then
            ! SH
            velx = -force
            vely = 0 !lx*force
            velz = 0 !ly*force
            !P/SV
            !  velx = lx*force
            !  vely = -lz*force
            !  velz = ly*force
        else
            velx = lx*force
            vely = ly*force
            velz = lz*force
        endif

        return
    end subroutine Inc_VelocPW

    ! ###############################################################################################


end module splanew

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
