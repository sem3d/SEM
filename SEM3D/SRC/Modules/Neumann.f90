!>
!! \file Neumann.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

module sneu

    ! Modified by Elise Delavaud 23/02/06
    ! #####################################################################################
    ! #####################################################################################

    type Face_Neu

       integer :: ngll1, ngll2, mat_index, dir, Face
       integer, dimension (0:3) :: Near_Edges, Near_Vertices, Orient_Edges
       real, dimension (:,:), pointer :: ds, MassMat
       real, dimension (:,:,:), pointer :: normal, Btn, Coord_nodes, Forces

    end type Face_Neu


    type Edge_Neu

       integer :: ngll, mat_index, dir, Edge
       real, dimension (:), pointer :: MassMat
       real, dimension (:,:), pointer :: Btn, Coord_nodes, Forces

    end type Edge_Neu


    type Vertex_Neu

       integer :: Vertex , mat_index
       real :: MassMat
       real, dimension (0:2) :: Btn, Coord_nodes, Forces

    end type Vertex_Neu

    type Param_Neu

       real :: Mu, Lambda, Kappa,speed, lx, ly, lz, xs, ys, zs, f0
       character (len=1) :: wtype

    end type Param_Neu


    type Neu

       integer :: n_faces, n_edges, n_vertices
       type(Face_Neu), dimension (:), pointer :: nFace
       type(Edge_Neu), dimension (:), pointer :: nEdge
       type(Vertex_Neu), dimension (:), pointer :: nVertex
       type(Param_Neu) :: nParam

    end type Neu

contains


    !  ########################################################################################
    subroutine compute_nforces_on_face (Face, Param, dt, ctime)

        use splanew
        implicit none

        type (Face_Neu), intent (INOUT) :: Face
        type (Param_Neu), intent (IN) :: Param
        real, intent (IN) :: dt, ctime

        integer :: i,j, ngll1, ngll2
        real, dimension (:,:,:), allocatable :: vel_i
        real :: xpt, ypt, zpt, velixn, veliyn, velizn, Sigma11, Sigma22 ,Sigma33, Sigma12, Sigma13, Sigma23


        ngll1 = Face%ngll1
        ngll2 = Face%ngll2

        allocate ( vel_i(1:ngll1-2,1:ngll2-2,0:2) )
        vel_i = 0

        do j = 1,ngll2-2
            do i = 1,ngll1-2

                xpt = Face%Coord_nodes(i,j,0)
                ypt = Face%Coord_nodes(i,j,1)
                zpt = Face%Coord_nodes(i,j,2)
                call Inc_VelocN (Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
                call Inc_VelocN (Param,xpt,ypt,zpt,ctime+dt,vel_i(i,j,0),vel_i(i,j,1),vel_i(i,j,2))

                vel_i(i,j,0) = 0.5*(velixn+vel_i(i,j,0))
                vel_i(i,j,1) = 0.5*(veliyn+vel_i(i,j,1))
                vel_i(i,j,2) = 0.5*(velizn+vel_i(i,j,2))

                Sigma11 = (Param%Lambda+2*Param%Mu)*Param%lx*vel_i(i,j,0)+Param%Lambda*(Param%ly*vel_i(i,j,1)+Param%lz*vel_i(i,j,2))
                Sigma22 = (Param%Lambda+2*Param%Mu)*Param%ly*vel_i(i,j,1)+Param%Lambda*(Param%lx*vel_i(i,j,0)+Param%lz*vel_i(i,j,2))
                Sigma33 = (Param%Lambda+2*Param%Mu)*Param%lz*vel_i(i,j,2)+Param%Lambda*(Param%lx*vel_i(i,j,0)+Param%ly*vel_i(i,j,1))
                Sigma12 = Param%Mu*(Param%ly*vel_i(i,j,0)+Param%lx*vel_i(i,j,1))
                Sigma13 = Param%Mu*(Param%lz*vel_i(i,j,0)+Param%lx*vel_i(i,j,2))
                Sigma23 = Param%Mu*(Param%lz*vel_i(i,j,1)+Param%ly*vel_i(i,j,2))

                Face%Forces(i,j,0) = -( Sigma11*Face%Btn(i,j,0) + Sigma12*Face%Btn(i,j,1) + Sigma13*Face%Btn(i,j,2) ) / Param%speed
                Face%Forces(i,j,1) = -( Sigma12*Face%Btn(i,j,0) + Sigma22*Face%Btn(i,j,1) + Sigma23*Face%Btn(i,j,2) ) / Param%speed
                Face%Forces(i,j,2) = -( Sigma13*Face%Btn(i,j,0) + Sigma23*Face%Btn(i,j,1) + Sigma33*Face%Btn(i,j,2) ) / Param%speed
                !if (Traction_i(i,j,0)>0 .or. Traction_i(i,j,1) > 0) print*,'Problem'
            enddo
        enddo


    end subroutine compute_nforces_on_face

    !  ########################################################################################
    subroutine compute_nforces_on_edge (Edge, Param, dt, ctime)

        use splanew
        implicit none

        type (Edge_Neu), intent (INOUT) :: Edge
        type (Param_Neu), intent (IN) :: Param
        real, intent (IN) :: dt, ctime

        integer :: i, ngll
        real :: xpt, ypt, zpt, velixn, veliyn, velizn, Sigma11, Sigma22 ,Sigma33, Sigma12, Sigma13, Sigma23
        real, dimension(:,:), allocatable :: vel_i

        ngll = Edge%ngll
        allocate ( vel_i(1:ngll-2,0:2) )
        vel_i = 0.

        do i = 1,ngll-2

            xpt = Edge%Coord_nodes(i,0)
            ypt = Edge%Coord_nodes(i,1)
            zpt = Edge%Coord_nodes(i,2)

            call Inc_VelocN (Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
            call Inc_VelocN (Param,xpt,ypt,zpt,ctime+dt,vel_i(i,0),vel_i(i,1),vel_i(i,2))

            vel_i(i,0) = 0.5*(velixn+vel_i(i,0))
            vel_i(i,1) = 0.5*(veliyn+vel_i(i,1))
            vel_i(i,2) = 0.5*(velizn+vel_i(i,2))

            Sigma11 = (Param%Lambda+2*Param%Mu)*Param%lx*vel_i(i,0)+Param%Lambda*(Param%ly*vel_i(i,1)+Param%lz*vel_i(i,2))
            Sigma22 = (Param%Lambda+2*Param%Mu)*Param%ly*vel_i(i,1)+Param%Lambda*(Param%lx*vel_i(i,0)+Param%lz*vel_i(i,2))
            Sigma33 = (Param%Lambda+2*Param%Mu)*Param%lz*vel_i(i,2)+Param%Lambda*(Param%lx*vel_i(i,0)+Param%ly*vel_i(i,1))
            Sigma12 = Param%Mu*(Param%ly*vel_i(i,0)+Param%lx*vel_i(i,1))
            Sigma13 = Param%Mu*(Param%lz*vel_i(i,0)+Param%lx*vel_i(i,2))
            Sigma23 = Param%Mu*(Param%lz*vel_i(i,1)+Param%ly*vel_i(i,2))

            Edge%Forces(i,0) = -( Sigma11*Edge%Btn(i,0) + Sigma12*Edge%Btn(i,1) + Sigma13*Edge%Btn(i,2) ) / Param%speed
            Edge%Forces(i,1) = -( Sigma12*Edge%Btn(i,0) + Sigma22*Edge%Btn(i,1) + Sigma23*Edge%Btn(i,2) ) / Param%speed
            Edge%Forces(i,2) = -( Sigma13*Edge%Btn(i,0) + Sigma23*Edge%Btn(i,1) + Sigma33*Edge%Btn(i,2) ) / Param%speed
        enddo

        deallocate (vel_i)


    end subroutine compute_nforces_on_edge

    !  ########################################################################################
    subroutine compute_nforces_on_vertex (Vertex, Param, dt, ctime)

        use splanew
        implicit none

        type (Vertex_Neu), intent (INOUT) :: Vertex
        type (Param_Neu), intent (IN) :: Param
        real, intent (IN) :: dt, ctime

        real :: xpt, ypt, zpt, velixn, veliyn, velizn, Sigma11, Sigma22 ,Sigma33, Sigma12, Sigma13, Sigma23
        real, dimension(0:2) :: vel_i

        vel_i = 0.

        xpt = Vertex%Coord_nodes(0)
        ypt = Vertex%Coord_nodes(1)
        zpt = Vertex%Coord_nodes(2)
        !print*,'V',rg,nf,xpt,ypt,zpt

        call Inc_VelocN (Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
        call Inc_VelocN (Param,xpt,ypt,zpt,ctime+dt,vel_i(0),vel_i(1),vel_i(2))

        vel_i(0) = 0.5*(velixn+vel_i(0))
        vel_i(1) = 0.5*(veliyn+vel_i(1))
        vel_i(2) = 0.5*(velizn+vel_i(2))

        Sigma11 = (Param%Lambda+2*Param%Mu)*Param%lx*vel_i(0)+Param%Lambda*(Param%ly*vel_i(1)+Param%lz*vel_i(2))
        Sigma22 = (Param%Lambda+2*Param%Mu)*Param%ly*vel_i(1)+Param%Lambda*(Param%lx*vel_i(0)+Param%lz*vel_i(2))
        Sigma33 = (Param%Lambda+2*Param%Mu)*Param%lz*vel_i(2)+Param%Lambda*(Param%lx*vel_i(0)+Param%ly*vel_i(1))
        Sigma12 = Param%Mu*(Param%ly*vel_i(0)+Param%lx*vel_i(1))
        Sigma13 = Param%Mu*(Param%lz*vel_i(0)+Param%lx*vel_i(2))
        Sigma23 = Param%Mu*(Param%lz*vel_i(1)+Param%ly*vel_i(2))

        Vertex%Forces(0) = -( Sigma11*Vertex%Btn(0) + Sigma12*Vertex%Btn(1) + Sigma13*Vertex%Btn(2) ) / Param%speed
        Vertex%Forces(1) = -( Sigma12*Vertex%Btn(0) + Sigma22*Vertex%Btn(1) + Sigma23*Vertex%Btn(2) ) / Param%speed
        Vertex%Forces(2) = -( Sigma13*Vertex%Btn(0) + Sigma23*Vertex%Btn(1) + Sigma33*Vertex%Btn(2) ) / Param%speed

    end subroutine compute_nforces_on_vertex

    ! #########################################################################################
    ! #########################################################################################
    subroutine Inc_VelocN (Param,x,y,z,ctime,velx,vely,velz)

        implicit none

        type (Param_Neu), intent (IN) :: Param
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

        !if ( Param%wtype == 'S' ) then
        !  velx = -lz*force
        !  vely = ly*force
        !  velz = lx*force
        !!  velx = lx*force
        !!  vely = -lz*force
        !!  velz = ly*force
        !else
        !  velx = lx*force
        !  vely = ly*force
        !  velz = lz*force
        !endif

        return
    end subroutine Inc_VelocN

    ! ###############################################################################################


end module sneu
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
