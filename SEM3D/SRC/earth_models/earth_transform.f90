!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module earth_transform
    use constants, only : fpp, M_PI
    implicit none

    real(fpp), parameter :: Earth_radius=6371000.d0, Pi180 = M_PI/180.0

contains

    subroutine getRotMat_loc2glob(lon, lat, M)
        real(fpp), intent(in) :: lon, lat
        real(fpp), intent(out), dimension(3,3) :: M

        real(fpp) :: theta, phi

        theta = (90.0-lat)*Pi180
        phi = lon*Pi180
        M(1,1) = -cos(theta)*cos(phi)
        M(1,2) = -sin(phi)
        M(1,3) = sin(theta)*cos(phi)
        M(2,1) = -cos(theta)*sin(phi)
        M(2,2) = cos(phi)
        M(2,3) = sin(theta)*sin(phi)
        M(3,1) = sin(theta)
        M(3,2) = 0
        M(3,3) = cos(theta)


       ! M(1,1) = cos(theta)*cos(phi)
       ! M(1,2) = -sin(phi)
       ! M(1,3) = sin(theta)*cos(phi)
       ! M(2,1) = cos(theta)*sin(phi)
       ! M(2,2) = cos(phi)
       ! M(2,3) = sin(theta)*sin(phi)
       ! M(3,1) = -sin(theta)
       ! M(3,2) = 0
       ! M(3,3) = cos(theta)

    end subroutine getRotMat_loc2glob

    subroutine cart2sph(x, y, z, r, theta, phi)
        real(fpp), intent(in) :: x, y, z
        real(fpp), intent(out) :: r, theta, phi

        r = sqrt(x*x+y*y+z*z)
        theta = acos(z/r)
        phi = atan2(y,x)
    end subroutine cart2sph




end module earth_transform


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
