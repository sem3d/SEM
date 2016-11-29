!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file DimensionalShape.f90
!! \brief
!!
!<
module DimensionalShape
    use constants, only : fpp
    implicit none
    real(fpp)  :: Radius

contains
    !-------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    subroutine surfaceSource(nbtot,Tdomain,surface)

       use sdomain
       use ssurf

       implicit none
       integer,       intent(in   ):: nbtot
       type(domain),  intent(in   ):: Tdomain
       type(SurfaceT),intent(inout):: surface
       integer                     :: i, shape, ssrc
       real(fpp)                :: x, y, z, RR

       allocate(surface%source(0:nbtot-1))
       ssrc = nsurfindex(Tdomain, surface%name)
       surface%source = 1.0d0
       shape = Tdomain%nsurfsource(ssrc-1)%shape
       RR = Tdomain%nsurfsource(ssrc-1)%size
       if (ssrc.gt.0) then
          do i=0,nbtot-1
             x = surface%coord(i,0) - Tdomain%nsurfsource(ssrc)%scoord(0)
             y = surface%coord(i,1) - Tdomain%nsurfsource(ssrc)%scoord(1)
             z = surface%coord(i,2) - Tdomain%nsurfsource(ssrc)%scoord(2)
             surface%source(i) = SrcShape(x, y, z, RR, shape)
          enddo
       endif

    end subroutine surfaceSource
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function SrcShape(x, y, z, R, index)

       implicit none
       integer,      intent(in) :: index
       real(fpp), intent(in)    :: x, y, R, z
       real(fpp)                :: SrcShape

       SrcShape = 0.d0
       select case (index)
              case (1)
                   SrcShape = gaussienne(x,y,z,R)
              case (2)
                   SrcShape = Paraboloid(x,y,z,R)
              case (3)
                   SrcShape = square(x,y,z,R)
              case (4)
                   SrcShape = cylinder(x,y,z,R)
              case default
                   SrcShape = 1.0d0
       end select

    end function SrcShape
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function gaussienne(x,y,z,Radius)

       implicit none
       real(fpp)             :: gaussienne
       real(fpp), intent(in) :: x, y, z, Radius

       gaussienne = exp(-(x**2 + y**2 + z**2)/Radius**2)

    end function gaussienne
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function Paraboloid(x,y,z,Radius)

       implicit none
       real(fpp)             :: Paraboloid
       real(fpp), intent(in) :: x, y, z, Radius
       real(fpp)             :: ray

       Paraboloid = 0.d0
       ray = dsqrt(x**2 + y**2 + z**2)
       if (ray.le.Radius) Paraboloid=dsqrt(1.-ray**2/Radius**2)

    end function Paraboloid
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function square(x,y,z,Radius)

        implicit none
        real(fpp)             :: square
        real(fpp), intent(in) :: x, y, z, Radius

        square = 0
        if ((abs(x).le.radius).and.(abs(y).le.Radius)) square = 1.0d0

    end function square
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function cylinder(x,y,z,Radius)

       implicit none
       real(fpp)             :: cylinder
       real(fpp), intent(in) :: x, y, z, Radius
       real(fpp)             :: ray

       cylinder=0.d0
       ray = dsqrt(x**2 + y**2 + z**2)
       if (ray.le.Radius) cylinder = 1.0d0

    end function cylinder
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function nsurfindex(Tdomain, name)
       use sdomain

       implicit none
       integer                        :: nsurfindex
       type(domain),       intent(in) :: Tdomain
       character(len=*),   intent(in) :: name
       integer                        :: ns, s, id, n

       nsurfindex = -1
       do ns = 0, Tdomain%nsurface-1
          do s = lbound(Tdomain%nsurfsource(ns)%index,1),ubound(Tdomain%nsurfsource(ns)%index,1)
             id = Tdomain%nsurfsource(ns)%index(s)
             do n=0,size(Tdomain%sSurfaces)-1
                if ((trim(Tdomain%sSurfaces(n)%name) == name(1:len_trim(name))).and.&
                    (Tdomain%nsurfsource(ns)%what_bc=='NE')) then
                     nsurfindex = ns+1
                endif
             enddo
          enddo
       enddo
    end function nsurfindex

end module DimensionalShape

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
