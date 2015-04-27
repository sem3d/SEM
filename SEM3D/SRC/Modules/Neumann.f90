!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file Neumann.f90
!! \brief
!!
!<

module sneu

    implicit none

    ! Neumann face properties
    type :: face_Neu
       integer :: Face,dir,ngll1,ngll2
       integer, dimension(0:3) :: Near_Edges, Near_Vertices
       ! all 3 edges and vertices for a Neumann face
       integer, dimension(0:3) :: Near_Edges_Orient
       real, dimension(:,:,:), pointer  :: normal,BtN,Forces,Coord
    end type face_Neu

    type :: edge_Neu
       integer :: Edge,ngll
       integer :: Orient_Edge
       real, dimension(:,:), pointer  :: BtN,Forces,Coord
    end type edge_Neu

    type :: vertex_Neu
       integer :: vertex
       real, dimension(0:2)  :: BtN,Forces,Coord
    end type vertex_Neu

    ! general parameters: description of the type of Neumann B.C.
    !   "R"=plane wave with a Ricker in time; "P"=plane wave with a sine
    !    in time; "F": time evolution given in a file.
    type Param_Neu
       integer  :: mat_index
       real :: Mu,Lambda,speed,lx,ly,lz,xs,ys,zs,f0
       character :: wtype,what_bc
    end type Param_Neu

    ! general Neumann object
    type :: Neu_object
       integer  :: Neu_n_faces, Neu_n_edges, Neu_n_vertices
       type(face_Neu), dimension(:), pointer  :: Neu_face
       type(edge_Neu), dimension(:), pointer  :: Neu_edge
       type(vertex_Neu), dimension(:), pointer  :: Neu_vertex
       type(Param_Neu) :: Neu_Param
    end type Neu_object

contains

    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine compute_Neu_forces_on_face(Face,Param,dt,ctime)
        ! gives the forces on faces submitted to Neumann B.C.
        implicit none

        type (Face_Neu), intent(inout) :: Face
        type(Param_Neu), intent(in) :: Param
        real, intent(in) :: dt,ctime
        integer :: i,j,ngll1,ngll2
        real, dimension(:,:,:), allocatable :: vel_i
        real :: xpt, ypt, zpt, velixn, veliyn, velizn,    &
            Sigma11, Sigma22 ,Sigma33, Sigma12, Sigma13, Sigma23


        ngll1 = Face%ngll1
        ngll2 = Face%ngll2
        allocate(vel_i(1:ngll1-2,1:ngll2-2,0:2))
        vel_i(:,:,:) = 0

        do j = 1,ngll2-2
            do i = 1,ngll1-2

                xpt = Face%Coord(i,j,0)
                ypt = Face%Coord(i,j,1)
                zpt = Face%Coord(i,j,2)

                select case(Param%what_bc)
                case('P','R')
                    !- velocities at 2 time steps..
                    call Inc_VelocN(Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
                    call Inc_VelocN(Param,xpt,ypt,zpt,ctime+dt,vel_i(i,j,0),vel_i(i,j,1),vel_i(i,j,2))
                    !-  .. to obtain the value at n+1/2
                    vel_i(i,j,0) = 0.5*(velixn+vel_i(i,j,0))
                    vel_i(i,j,1) = 0.5*(veliyn+vel_i(i,j,1))
                    vel_i(i,j,2) = 0.5*(velizn+vel_i(i,j,2))

                    Sigma11 = (Param%Lambda+2*Param%Mu)*Param%lx*vel_i(i,j,0)+    &
                        Param%Lambda*(Param%ly*vel_i(i,j,1)+Param%lz*vel_i(i,j,2))
                    Sigma22 = (Param%Lambda+2*Param%Mu)*Param%ly*vel_i(i,j,1)+    &
                        Param%Lambda*(Param%lx*vel_i(i,j,0)+Param%lz*vel_i(i,j,2))
                    Sigma33 = (Param%Lambda+2*Param%Mu)*Param%lz*vel_i(i,j,2)+    &
                        Param%Lambda*(Param%lx*vel_i(i,j,0)+Param%ly*vel_i(i,j,1))
                    Sigma12 = Param%Mu*(Param%ly*vel_i(i,j,0)+Param%lx*vel_i(i,j,1))
                    Sigma13 = Param%Mu*(Param%lz*vel_i(i,j,0)+Param%lx*vel_i(i,j,2))
                    Sigma23 = Param%Mu*(Param%lz*vel_i(i,j,1)+Param%ly*vel_i(i,j,2))

                    Face%Forces(i,j,0) = -(Sigma11*Face%Btn(i,j,0)+Sigma12*Face%Btn(i,j,1)+ &
                        Sigma13*Face%Btn(i,j,2))/Param%speed
                    Face%Forces(i,j,1) = -(Sigma12*Face%Btn(i,j,0)+Sigma22*Face%Btn(i,j,1)+ &
                        Sigma23*Face%Btn(i,j,2))/Param%speed
                    Face%Forces(i,j,2) = -(Sigma13*Face%Btn(i,j,0)+Sigma23*Face%Btn(i,j,1)+ &
                        Sigma33*Face%Btn(i,j,2))/Param%speed
                case('F')
                    !- velocities at 2 time steps..
                    call Inc_VelocN(Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
                    call Inc_VelocN(Param,xpt,ypt,zpt,ctime+dt,vel_i(i,j,0),vel_i(i,j,1),vel_i(i,j,2))
                    !-  .. to obtain the value at n+1/2
                    vel_i(i,j,0) = 0.5*(velixn+vel_i(i,j,0))
                    vel_i(i,j,1) = 0.5*(veliyn+vel_i(i,j,1))
                    vel_i(i,j,2) = 0.5*(velizn+vel_i(i,j,2))


                    Face%Forces(i,j,0) = -vel_i(i,j,0)*Face%Btn(i,j,0)
                    Face%Forces(i,j,1) = -vel_i(i,j,1)*Face%Btn(i,j,1)
                    Face%Forces(i,j,2) = -vel_i(i,j,2)*Face%Btn(i,j,2)
                
            !    case('S')
            !        ! Pour le case de source surfacique
            !        Face%Forces(i,j,0) = triangle(dt)*Face%Btn(i,j,0)
            !        Face%Forces(i,j,1) = triangle(dt)*Face%Btn(i,j,1)
            !        Face%Forces(i,j,2) = triangle(dt)*Face%Btn(i,j,2) 

                end select
            enddo
        enddo

        deallocate(vel_i)

    end subroutine compute_Neu_forces_on_face
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine compute_Neu_forces_on_edge(Edge,Param,dt,ctime)
        ! gives the forces on faces submitted to Neumann B.C.
        implicit none

        type(Edge_Neu), intent(inout) :: Edge
        type(Param_Neu), intent(in) :: Param
        real, intent(in) :: dt,ctime
        integer :: i, ngll
        real :: xpt,ypt,zpt,velixn,veliyn,velizn,    &
            Sigma11,Sigma22,Sigma33,Sigma12,Sigma13,Sigma23
        real, dimension(:,:), allocatable :: vel_i

        ngll = Edge%ngll
        allocate(vel_i(1:ngll-2,0:2))
        vel_i(:,:) = 0d0

        do i = 1,ngll-2
            xpt = Edge%Coord(i,0)
            ypt = Edge%Coord(i,1)
            zpt = Edge%Coord(i,2)

            select case(Param%what_bc)
            case('P','R')
                call Inc_VelocN(Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
                call Inc_VelocN(Param,xpt,ypt,zpt,ctime+dt,vel_i(i,0),vel_i(i,1),vel_i(i,2))

                vel_i(i,0) = 0.5*(velixn+vel_i(i,0))
                vel_i(i,1) = 0.5*(veliyn+vel_i(i,1))
                vel_i(i,2) = 0.5*(velizn+vel_i(i,2))

                Sigma11 = (Param%Lambda+2*Param%Mu)*Param%lx*vel_i(i,0)+   &
                    Param%Lambda*(Param%ly*vel_i(i,1)+Param%lz*vel_i(i,2))
                Sigma22 = (Param%Lambda+2*Param%Mu)*Param%ly*vel_i(i,1)+   &
                    Param%Lambda*(Param%lx*vel_i(i,0)+Param%lz*vel_i(i,2))
                Sigma33 = (Param%Lambda+2*Param%Mu)*Param%lz*vel_i(i,2)+   &
                    Param%Lambda*(Param%lx*vel_i(i,0)+Param%ly*vel_i(i,1))
                Sigma12 = Param%Mu*(Param%ly*vel_i(i,0)+Param%lx*vel_i(i,1))
                Sigma13 = Param%Mu*(Param%lz*vel_i(i,0)+Param%lx*vel_i(i,2))
                Sigma23 = Param%Mu*(Param%lz*vel_i(i,1)+Param%ly*vel_i(i,2))

                Edge%Forces(i,0) = -(Sigma11*Edge%Btn(i,0)+Sigma12*Edge%Btn(i,1)+  &
                    Sigma13*Edge%Btn(i,2))/Param%speed
                Edge%Forces(i,1) = -(Sigma12*Edge%Btn(i,0)+Sigma22*Edge%Btn(i,1)+  &
                    Sigma23*Edge%Btn(i,2))/Param%speed
                Edge%Forces(i,2) = -(Sigma13*Edge%Btn(i,0)+Sigma23*Edge%Btn(i,1)+  &
                    Sigma33*Edge%Btn(i,2))/Param%speed
            case('F')
                !- velocities at 2 time steps..
                call Inc_VelocN(Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
                call Inc_VelocN(Param,xpt,ypt,zpt,ctime+dt,vel_i(i,0),vel_i(i,1),vel_i(i,2))
                !-  .. to obtain the value at n+1/2
                vel_i(i,0) = 0.5*(velixn+vel_i(i,0))
                vel_i(i,1) = 0.5*(veliyn+vel_i(i,1))
                vel_i(i,2) = 0.5*(velizn+vel_i(i,2))


                Edge%Forces(i,0) = -vel_i(i,0)*Edge%Btn(i,0)
                Edge%Forces(i,1) = -vel_i(i,1)*Edge%Btn(i,1)
                Edge%Forces(i,2) = -vel_i(i,2)*Edge%Btn(i,2)
            
         !   case('S')
         !       ! Pour la source surfacique
         !       Edge%Forces(i,0) = triangle(dt)*Edge%Btn(i,0)
         !       Edge%Forces(i,1) = triangle(dt)*Edge%Btn(i,1)
         !       Edge%Forces(i,2) = triangle(dt)*Edge%Btn(i,2)

            end select
        enddo

        deallocate(vel_i)

    end subroutine compute_Neu_forces_on_edge
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    subroutine compute_Neu_forces_on_vertex(Vertex,ind,Param,dt,ctime)
        ! gives the forces on faces submitted to Neumann B.C.
        implicit none

        type(Vertex_Neu), intent(inout) :: Vertex
        type(Param_Neu), intent(in) :: Param
        real, intent(in) :: dt, ctime
        integer, intent(in)  :: ind
        real :: xpt,ypt,zpt,velixn,veliyn,velizn,   &
            Sigma11,Sigma22,Sigma33,Sigma12,Sigma13,Sigma23
        real, dimension(0:2) :: vel_i

        vel_i = 0.

        xpt = Vertex%Coord(0)
        ypt = Vertex%Coord(1)
        zpt = Vertex%Coord(2)

        select case(Param%what_bc)
        case('P','R')
            call Inc_VelocN(Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)

            call Inc_VelocN(Param,xpt,ypt,zpt,ctime+dt,vel_i(0),vel_i(1),vel_i(2))

            vel_i(0) = 0.5*(velixn+vel_i(0))
            vel_i(1) = 0.5*(veliyn+vel_i(1))
            vel_i(2) = 0.5*(velizn+vel_i(2))


            Sigma11 = (Param%Lambda+2*Param%Mu)*Param%lx*vel_i(0)+   &
                Param%Lambda*(Param%ly*vel_i(1)+Param%lz*vel_i(2))
            Sigma22 = (Param%Lambda+2*Param%Mu)*Param%ly*vel_i(1)+   &
                Param%Lambda*(Param%lx*vel_i(0)+Param%lz*vel_i(2))
            Sigma33 = (Param%Lambda+2*Param%Mu)*Param%lz*vel_i(2)+   &
                Param%Lambda*(Param%lx*vel_i(0)+Param%ly*vel_i(1))
            Sigma12 = Param%Mu*(Param%ly*vel_i(0)+Param%lx*vel_i(1))
            Sigma13 = Param%Mu*(Param%lz*vel_i(0)+Param%lx*vel_i(2))
            Sigma23 = Param%Mu*(Param%lz*vel_i(1)+Param%ly*vel_i(2))

            Vertex%Forces(0) = -(Sigma11*Vertex%Btn(0)+Sigma12*Vertex%Btn(1)+   &
                Sigma13*Vertex%Btn(2))/Param%speed
            Vertex%Forces(1) = -(Sigma12*Vertex%Btn(0)+Sigma22*Vertex%Btn(1)+   &
                Sigma23*Vertex%Btn(2))/Param%speed
            Vertex%Forces(2) = -(Sigma13*Vertex%Btn(0)+Sigma23*Vertex%Btn(1)+   &
                Sigma33*Vertex%Btn(2))/Param%speed
        case('F')
            !- velocities at 2 time steps..
            call Inc_VelocN(Param,xpt,ypt,zpt,ctime,velixn,veliyn,velizn)
            call Inc_VelocN(Param,xpt,ypt,zpt,ctime+dt,vel_i(0),vel_i(1),vel_i(2))
            !-  .. to obtain the value at n+1/2
            vel_i(0) = 0.5*(velixn+vel_i(0))
            vel_i(1) = 0.5*(veliyn+vel_i(1))
            vel_i(2) = 0.5*(velizn+vel_i(2))


            Vertex%Forces(0) = -vel_i(0)*Vertex%Btn(0)
            Vertex%Forces(1) = -vel_i(1)*Vertex%Btn(1)
            Vertex%Forces(2) = -vel_i(2)*Vertex%Btn(2)


        ! case('S')
                ! Pour la source surfacique
        !    Vertex%Forces(0) = triangle(dt)*Vertex%Btn(0)
        !    Vertex%Forces(1) = triangle(dt)*Vertex%Btn(1)
        !    Vertex%Forces(2) = triangle(dt)*Vertex%Btn(2)


        end select

    end subroutine compute_Neu_forces_on_vertex
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    subroutine Inc_VelocN(Param,x,y,z,ctime,velx,vely,velz)

        implicit none

        type(Param_Neu), intent(in) :: Param
        real, intent(in) :: x,z,y,ctime
        real, intent(out) :: velx,vely,velz
        real :: lx,ly,lz,phasex,phasey,phasez,f0,force
        real :: vel,xsour,ysour,zsour,xpt,ypt,zpt,arg,pi

        pi = Acos(-1d0)
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

        if(Param%wtype == 'R')then   ! Ricker in time
            !arg = pi*f0*(phasex + phasey + phasez)
            arg = pi*f0*(ctime-0.4)
            force = (1.-2.*arg**2)*exp(-arg**2)     ! Plane Wave in Velocity
        else if(Param%wtype == 'P')then  ! sine in time
            arg = phasex+phasey+phasez
            force = sin(2d0*pi*arg)
        end if

        if(Param%wtype == 'R')then ! SH
            !    velx = -force
            !    vely = 0d0
            !    velz = 0d0
            !else
            velx = lx*force
            vely = ly*force
            velz = lz*force
        endif


        return
    end subroutine Inc_VelocN
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------

end module sneu

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
