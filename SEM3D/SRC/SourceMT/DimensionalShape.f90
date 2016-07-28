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

implicit none
real(kind=8)  :: Radius
  
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
       real(kind=8)                :: x, y, z, RR
    
       allocate(surface%source(0:nbtot-1))
       ssrc = nsurfindex(Tdomain, surface%name)
       surface%source = 1.0d0
       shape = Tdomain%nsurfsource(ssrc)%shape
       RR = Tdomain%nsurfsource(ssrc)%size
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
       real(kind=8), intent(in) :: x, y, R, z
       real(kind=8)             :: SrcShape
    
       Radius = R
       SrcShape = 0.d0
       select case (index)
              case (1)
                   SrcShape = gaussienne(x,y,z) 
              case (2)
                   SrcShape = Paraboloid(x,y,z)
              case (3)
                   SrcShape = square(x,y,z)
              case (4)
                   SrcShape = cylinder(x,y,z)
              case default
                   SrcShape = 1.0d0
       end select
    
    end function SrcShape 
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function gaussienne(x,y,z)
    
       implicit none
       real(kind=8)             :: gaussienne
       real(kind=8), intent(in) :: x, y, z
       
       gaussienne = exp(-(x**2 + y**2 + z**2)/Radius**2)
    
    end function gaussienne
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function Paraboloid(x,y,z)
       
       implicit none
       real(kind=8)             :: Paraboloid
       real(kind=8), intent(in) :: x, y, z
       real(kind=8)             :: ray
    
       Paraboloid = 0.d0
       ray = sqrt(x**2 + y**2 + z**2)
       if (ray.le.Radius) Paraboloid=sqrt(1.-ray**2/Radius**2) 
    
    end function Paraboloid
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function square(x,y,z)
        
        implicit none
        real(kind=8)             :: square
        real(kind=8), intent(in) :: x, y, z
    
        square = 0
        if ((abs(x).le.radius).and.(abs(y).le.Radius)) square = 1.0d0

    end function square
    !----------------------------------------------------------------------------------
    !----------------------------------------------------------------------------------
    function cylinder(x,y,z)
    
       implicit none
       real(kind=8)             :: cylinder
       real(kind=8), intent(in) :: x, y, z
       real(kind=8)             :: ray
    
       cylinder=0.d0
       ray = sqrt(x**2 + y**2 + z**2)
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
       integer                        :: ns, s, id
         
       nsurfindex = -1
       do ns = 0, Tdomain%nsurface-1
          do s = lbound(Tdomain%nsurfsource(ns)%index,1),ubound(Tdomain%nsurfsource(ns)%index,1)
             id = Tdomain%nsurfsource(ns)%index(s)
             if ((trim(Tdomain%sSurfaces(id)%name) == trim(name)).and.&
                 (Tdomain%nsurfsource(ns)%what_bc=='NE')) then
                  nsurfindex = ns
             endif
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


