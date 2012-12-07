real function f_p(t,xi,eta,zeta)
  !- gives the t-coordinate of a GLL node.
    implicit none
    real, dimension(0:7), intent(in)  :: t
    real, intent(in) :: xi,eta,zeta

    f_p = 0.125 * (t(0)*(1-xi)*(1-eta)*(1-zeta) + t(1)*(1+xi)*(1-eta)*(1-zeta) + t(2)*(1+xi)*(1+eta)*(1-zeta) + &
           t(3)*(1-xi)*(1+eta)*(1-zeta) + t(4)*(1-xi)*(1-eta)*(1+zeta) + t(5)*(1+xi)*(1-eta)*(1+zeta) + &
           t(6)*(1+xi)*(1+eta)*(1+zeta) + t(7)*(1-xi)*(1+eta)*(1+zeta))
end function f_p
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function f_xp(x,xi,eta,zeta)
  !- gives the x-coordinate of a GLL node.
    implicit none
    real, dimension(0:7), intent(in)  :: x
    real, intent(in) :: xi,eta,zeta

    f_xp = 0.125 * (x(0)*(1-xi)*(1-eta)*(1-zeta) + x(1)*(1+xi)*(1-eta)*(1-zeta) + x(2)*(1+xi)*(1+eta)*(1-zeta) + &
           x(3)*(1-xi)*(1+eta)*(1-zeta) + x(4)*(1-xi)*(1-eta)*(1+zeta) + x(5)*(1+xi)*(1-eta)*(1+zeta) + &
           x(6)*(1+xi)*(1+eta)*(1+zeta) + x(7)*(1-xi)*(1+eta)*(1+zeta))
end function f_xp
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function f_yp(y,xi,eta,zeta)
  !- gives the y-coordinate of a GLL node.
    implicit none
    real, dimension(0:7), intent(in)  :: y
    real, intent(in) :: xi,eta,zeta

    f_yp = 0.125 * (y(0)*(1-xi)*(1-eta)*(1-zeta) + y(1)*(1+xi)*(1-eta)*(1-zeta) + y(2)*(1+xi)*(1+eta)*(1-zeta) + &
           y(3)*(1-xi)*(1+eta)*(1-zeta) + y(4)*(1-xi)*(1-eta)*(1+zeta) + y(5)*(1+xi)*(1-eta)*(1+zeta) + &
           y(6)*(1+xi)*(1+eta)*(1+zeta) + y(7)*(1-xi)*(1+eta)*(1+zeta))
end function f_yp
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function f_zp(z,xi,eta,zeta)
  !- gives the z-coordinate of a GLL node.
    implicit none
    real, dimension(0:7), intent(in)  :: z
    real, intent(in) :: xi,eta,zeta

    f_zp = 0.125 * (z(0)*(1-xi)*(1-eta)*(1-zeta) + z(1)*(1+xi)*(1-eta)*(1-zeta) + z(2)*(1+xi)*(1+eta)*(1-zeta) + &
           z(3)*(1-xi)*(1+eta)*(1-zeta) + z(4)*(1-xi)*(1-eta)*(1+zeta) + z(5)*(1+xi)*(1-eta)*(1+zeta) + &
           z(6)*(1+xi)*(1+eta)*(1+zeta) + z(7)*(1-xi)*(1+eta)*(1+zeta))
end function f_zp
!-------------------------------------------------------------------------
!--   PARTIAL DERIVATIVES at GLL NODES   --
!-------------------------------------------------------------------------
real function der_dx_dxi(x,eta,zeta)
    implicit none
    real, dimension(0:7), intent(in)  :: x
    real, intent(in)  :: eta,zeta

    der_dx_dxi = 0.125 * ((x(1)-x(0))*(1-eta)*(1-zeta) + (x(2)-x(3))*(1+eta)*(1-zeta) +   &
                          (x(5)-x(4))*(1-eta)*(1+zeta) + (x(6)-x(7))*(1+eta)*(1+zeta))

end function der_dx_dxi
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function der_dx_deta(x,xi,zeta)
    implicit none
    real, dimension(0:7), intent(in)  :: x
    real, intent(in)  :: xi,zeta

    der_dx_deta = 0.125 * ((x(3)-x(0))*(1-xi)*(1-zeta) + (x(2)-x(1))*(1+xi)*(1-zeta) +   &
                          (x(7)-x(4))*(1-xi)*(1+zeta) + (x(6)-x(5))*(1+xi)*(1+zeta))

end function der_dx_deta
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function der_dx_dzeta(x,xi,eta)
    implicit none
    real, dimension(0:7), intent(in)  :: x
    real, intent(in)  :: xi,eta

    der_dx_dzeta = 0.125 * ((x(4)-x(0))*(1-xi)*(1-eta) + (x(5)-x(1))*(1+xi)*(1-eta) +   &
                          (x(7)-x(3))*(1-xi)*(1+eta) + (x(6)-x(2))*(1+xi)*(1+eta))

end function der_dx_dzeta
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function der_dy_dxi(y,eta,zeta)
    implicit none
    real, dimension(0:7), intent(in)  :: y
    real, intent(in)  :: eta,zeta

    der_dy_dxi = 0.125 * ((y(1)-y(0))*(1-eta)*(1-zeta) + (y(2)-y(3))*(1+eta)*(1-zeta) +  &
                          (y(5)-y(4))*(1-eta)*(1+zeta) + (y(6)-y(7))*(1+eta)*(1+zeta))

end function der_dy_dxi
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function der_dy_deta(y,xi,zeta)
    implicit none
    real, dimension(0:7), intent(in)  :: y
    real, intent(in)  :: xi,zeta
     
    der_dy_deta = 0.125 * ((y(3)-y(0))*(1-xi)*(1-zeta) + (y(2)-y(1))*(1+xi)*(1-zeta) +  &
                          (y(7)-y(4))*(1-xi)*(1+zeta) + (y(6)-y(5))*(1+xi)*(1+zeta))

end function der_dy_deta
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function der_dy_dzeta(y,xi,eta)
    implicit none
    real, dimension(0:7), intent(in)  :: y
    real, intent(in)  :: xi,eta

    der_dy_dzeta = 0.125 * ((y(4)-y(0))*(1-xi)*(1-eta) + (y(5)-y(1))*(1+xi)*(1-eta) +   &
                          (y(7)-y(3))*(1-xi)*(1+eta) + (y(6)-y(2))*(1+xi)*(1+eta))

end function der_dy_dzeta
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function der_dz_dxi(z,eta,zeta)
    implicit none
    real, dimension(0:7), intent(in)  :: z
    real, intent(in)  :: eta,zeta

    der_dz_dxi = 0.125 * ((z(1)-z(0))*(1-eta)*(1-zeta) + (z(2)-z(3))*(1+eta)*(1-zeta) +   &
                          (z(5)-z(4))*(1-eta)*(1+zeta) + (z(6)-z(7))*(1+eta)*(1+zeta))

end function der_dz_dxi
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function der_dz_deta(z,xi,zeta)
    implicit none
    real, dimension(0:7), intent(in)  :: z
    real, intent(in)  :: xi,zeta

    der_dz_deta = 0.125 * ((z(3)-z(0))*(1-xi)*(1-zeta) + (z(2)-z(1))*(1+xi)*(1-zeta) +   &
                          (z(7)-z(4))*(1-xi)*(1+zeta) + (z(6)-z(5))*(1+xi)*(1+zeta))

end function der_dz_deta
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
real function der_dz_dzeta(z,xi,eta)
    implicit none
    real, dimension(0:7), intent(in)  :: z
    real, intent(in)  :: xi,eta

    der_dz_dzeta = 0.125 * ((z(4)-z(0))*(1-xi)*(1-eta) + (z(5)-z(1))*(1+xi)*(1-eta) +   &
                          (z(7)-z(3))*(1-xi)*(1+eta) + (z(6)-z(2))*(1+xi)*(1+eta))

end function der_dz_dzeta
!------------------------------------------------------------------------
!------------------------------------------------------------------------
subroutine normal_face(dir,ngllx,nglly,ngllz,ngll1,ngll2,x,y,z,GLLcx,GLLcy,GLLcz,normal)
  ! determines the normal to a given face
    implicit none
    integer, intent(in)   :: dir,ngll1,ngll2,ngllx,nglly,ngllz
    real, dimension(0:7), intent(in)  :: x,y,z
    real, dimension(0:ngllx-1), intent(in)  :: GLLcx  
    real, dimension(0:nglly-1), intent(in)  :: GLLcy  
    real, dimension(0:ngllz-1), intent(in)  :: GLLcz  
    real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(out)  :: normal
    integer   :: i,j
    real      :: xi,eta,zeta
    real, dimension(0:2)  :: n1,n2

  ! loop on GLL points of the face
    do j = 0,ngll2-1
      do i = 0,ngll1-1
        select case(dir)
          case(0)
            xi = GLLcx(i) ; eta = GLLcy(j) ; zeta = -1d0 
          case(1)
            xi = GLLcx(i) ; zeta = GLLcz(j) ; eta = -1d0 
          case(2)
            eta = GLLcy(i) ; zeta = GLLcz(j) ; xi = 1d0 
          case(3)
            xi = GLLcx(i) ; zeta = GLLcz(j) ; eta = 1d0 
          case(4)
            eta = GLLcy(i) ; zeta = GLLcz(j) ; xi = -1d0 
          case(5)
            xi = GLLcx(i) ; eta = GLLcy(j) ; zeta = 1d0 
        end select
        call normal_vector_def(dir,x,y,z,xi,eta,zeta,n1,n2)
      ! normal = cross product of n1 and n2
        call cross_prod(n1,n2,normal(i,j,0:2))
      end do
    end do     

end subroutine normal_face
!--------------------------------------------------------------------
!--------------------------------------------------------------------
subroutine normal_vector_def(dir,x,y,z,xi,eta,zeta,n1,n2)
  ! determination of the 2 vectors defining a given face
    implicit none
    integer, intent(in)  :: dir
    real, dimension(0:7), intent(in)  :: x,y,z
    real, intent(in)  :: xi,eta,zeta
    real, dimension(0:2), intent(out) :: n1,n2
    real, external  :: der_dx_dxi,der_dx_deta,der_dx_dzeta,der_dy_dxi,der_dy_deta,   &
                       der_dy_dzeta,der_dz_dxi,der_dz_deta,der_dz_dzeta

   ! X = position vector 
    select case(dir)
        case(0,5)  ! n1 = dX/dxi , n2 = dX/deta
          n1(0) = der_dx_dxi(x,eta,zeta) ; n1(1) = der_dy_dxi(y,eta,zeta) ; n1(2) = der_dz_dxi(z,eta,zeta)
          n2(0) = der_dx_deta(x,xi,zeta) ; n2(1) = der_dy_deta(y,xi,zeta) ; n2(2) = der_dz_deta(z,xi,zeta)
        case(1,3)  ! n1 = dX/dxi , n2= dX/dzeta
          n1(0) = der_dx_dxi(x,eta,zeta) ; n1(1) = der_dy_dxi(y,eta,zeta) ; n1(2) = der_dz_dxi(z,eta,zeta)
          n2(0) = der_dx_dzeta(x,xi,eta) ; n2(1) = der_dy_dzeta(y,xi,eta) ; n2(2) = der_dz_dzeta(z,xi,eta)
        case(2,4)  ! n1 = dX/deta , n2= dX/dzeta
          n1(0) = der_dx_deta(x,xi,zeta) ; n1(1) = der_dy_deta(y,xi,zeta) ; n1(2) = der_dz_deta(z,xi,zeta)
          n2(0) = der_dx_dzeta(x,xi,eta) ; n2(1) = der_dy_dzeta(y,xi,eta) ; n2(2) = der_dz_dzeta(z,xi,eta)
    end select

end subroutine normal_vector_def
!------------------------------------------------------------------
!------------------------------------------------------------------
subroutine cross_prod(x,y,z)
  ! z = cross product of vectors x and y
    implicit none
    real, dimension(0:2), intent(in)  :: x,y
    real, dimension(0:2), intent(out) :: z

    z(0) = x(1)*y(2)-x(2)*y(1)
    z(1) = x(2)*y(0)-x(0)*y(2)
    z(2) = x(0)*y(1)-x(1)*y(0)

end subroutine cross_prod
