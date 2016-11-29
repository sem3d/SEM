!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module deriv3d
    use constants
contains
    ! TODO : kill elem_part_deriv not used ?!...
    subroutine elem_part_deriv(ngllx,nglly,ngllz,hprimex,hprimey,hprimez,Scalp,    &
        dScalp_dxi,dScalp_deta,dScalp_dzeta)
        !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real(fpp), dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hprimex
        real(fpp), dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real(fpp), dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real(fpp), dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Scalp
        real(fpp), dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dScalp_dxi,dScalp_deta,dScalp_dzeta
        integer :: i,j,k,l

        dScalp_dxi = ZERO
        dScalp_deta = ZERO
        dScalp_dzeta = ZERO
        DO K = 0, ngllz-1
            ! d(Scalp)_dxi
            DO J = 0, nglly-1
                DO L = 0, ngllx-1
                    DO I = 0, ngllx-1
                        dScalp_dxi(I,J,K) = dScalp_dxi(I,J,K) + Scalp(L,J,K)*hprimex(L,I)
                    END DO
                END DO
            END DO

            ! d(Scalp)_deta
            DO J = 0, nglly-1
                DO L = 0, nglly-1
                    DO I = 0, ngllx-1
                        Dscalp_deta( I, J, K ) = dScalp_deta( I, J, K ) + Scalp( I, L, K )*hprimey(L,J)
                    END DO
                END DO
            END DO

            ! d(Scalp)_dzeta
            DO L = 0,ngllz-1
                DO J = 0,nglly-1
                    DO I = 0,ngllx-1
                        dScalp_dzeta(I,J,K) = dScalp_dzeta(I,J,K) + Scalp(I,J,L)*hprimez(L,K)
                    END DO
                END DO
            END DO
        END DO
    end subroutine elem_part_deriv
    !------------------------------------------------------------------------------------------
    !------------------------------------------------------------------------------------------
    subroutine physical_part_deriv(ngll,hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
        !- partial derivatives of the scalar property Scalp, with respect to x,y,z
        implicit none
        integer, intent(in) :: ngll
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: InvGrad
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Scalp
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: dS_dx,dS_dy,dS_dz

        ! Simple optimisation that allows the compiler to unroll loops for the most common cases
        if (ngll==5) then
            call physical_part_deriv_555(hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
            return
        else if (ngll==6) then
            call physical_part_deriv_666(hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
            return
        else if (ngll==7) then
            call physical_part_deriv_777(hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
            return
        else if (ngll==8) then
            call physical_part_deriv_888(hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
            return
        endif
        call physical_part_deriv_nnn(ngll,hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
    end subroutine physical_part_deriv

    subroutine physical_part_deriv_nnn(ngll,hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
        !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, intent(in) :: ngll
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: InvGrad
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Scalp
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: dS_dx,dS_dy,dS_dz

        integer :: i,j,k

        DO K = 0, ngll-1
            ! d(Scalp)_dxi
            DO J = 0, ngll-1
                DO I = 0, ngll-1
                    call physical_part_deriv_ijk(i,j,k,ngll,hprime,InvGrad(:,:,i,j,k),Scalp,&
                         dS_dx(i,j,k),dS_dy(i,j,k),dS_dz(i,j,k))
                END DO
            END DO
        END DO
    end subroutine physical_part_deriv_nnn

    subroutine physical_part_deriv_ijk(i,j,k,ngll,hprime,InvGrad_ijk,Scalp,dS_dx,dS_dy,dS_dz)
        implicit none
        integer :: i,j,k
        integer, intent(in) :: ngll
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:2,0:2), intent(in) :: InvGrad_ijk
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Scalp
        real(fpp), intent(out) :: dS_dx,dS_dy,dS_dz

        integer :: l
        real(fpp) :: dS_dxi, dS_deta, dS_dzeta

        dS_dxi   = 0.0D+0
        dS_deta  = 0.0D+0
        dS_dzeta = 0.0D+0
        DO L = 0, ngll-1
            dS_dxi   = dS_dxi  +Scalp(L,J,K)*hprime(L,I)
            dS_deta  = dS_deta +Scalp(I,L,K)*hprime(L,J)
            dS_dzeta = dS_dzeta+Scalp(I,J,L)*hprime(L,K)
        END DO
        !- in the physical domain
        dS_dx = dS_dxi*InvGrad_ijk(0,0)+dS_deta*InvGrad_ijk(0,1)+dS_dzeta*InvGrad_ijk(0,2)
        dS_dy = dS_dxi*InvGrad_ijk(1,0)+dS_deta*InvGrad_ijk(1,1)+dS_dzeta*InvGrad_ijk(1,2)
        dS_dz = dS_dxi*InvGrad_ijk(2,0)+dS_deta*InvGrad_ijk(2,1)+dS_dzeta*InvGrad_ijk(2,2)
    end subroutine physical_part_deriv_ijk

    subroutine physical_part_deriv_555(hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
        !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, parameter :: ngll=5
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: InvGrad
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Scalp
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: dS_dx,dS_dy,dS_dz
        real(fpp) :: dS_dxi, dS_deta, dS_dzeta
        integer :: i,j,k,l

        DO K = 0, ngll-1
            ! d(Scalp)_dxi
            DO J = 0, ngll-1
                DO I = 0, ngll-1
                    dS_dxi = ZERO
                    dS_deta = ZERO
                    dS_dzeta = ZERO
                    DO L = 0, ngll-1
                        dS_dxi   = dS_dxi   + Scalp(L,J,K)*hprime(L,I)
                        DS_deta  = dS_deta  + Scalp(I,L,K)*hprime(L,J)
                        dS_dzeta = dS_dzeta + Scalp(I,J,L)*hprime(L,K)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dxi*InvGrad(0,0,I,J,K) + dS_deta*InvGrad(0,1,I,J,K) +  &
                        dS_dzeta*InvGrad(0,2,I,J,K)
                    dS_dy(I,J,K) = dS_dxi*InvGrad(1,0,I,J,K) + dS_deta*InvGrad(1,1,I,J,K) +  &
                        dS_dzeta*InvGrad(1,2,I,J,K)
                    dS_dz(I,J,K) = dS_dxi*InvGrad(2,0,I,J,K) + dS_deta*InvGrad(2,1,I,J,K) +  &
                        dS_dzeta*InvGrad(2,2,I,J,K)
                END DO
            END DO
        END DO
    end subroutine physical_part_deriv_555

    subroutine physical_part_deriv_666(hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
        !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, parameter :: ngll=6
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: InvGrad
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Scalp
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: dS_dx,dS_dy,dS_dz
        real(fpp) :: dS_dxi, dS_deta, dS_dzeta
        integer :: i,j,k,l


        DO K = 0, ngll-1
            ! d(Scalp)_dxi
            DO J = 0, ngll-1
                DO I = 0, ngll-1
                    dS_dxi = ZERO
                    dS_deta = ZERO
                    dS_dzeta = ZERO
                    DO L = 0, ngll-1
                        dS_dxi   = dS_dxi   + Scalp(L,J,K)*hprime(L,I)
                        DS_deta  = dS_deta  + Scalp(I,L,K)*hprime(L,J)
                        dS_dzeta = dS_dzeta + Scalp(I,J,L)*hprime(L,K)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dxi*InvGrad(0,0,I,J,K) + dS_deta*InvGrad(0,1,I,J,K) +  &
                        dS_dzeta*InvGrad(0,2,I,J,K)
                    dS_dy(I,J,K) = dS_dxi*InvGrad(1,0,I,J,K) + dS_deta*InvGrad(1,1,I,J,K) +  &
                        dS_dzeta*InvGrad(1,2,I,J,K)
                    dS_dz(I,J,K) = dS_dxi*InvGrad(2,0,I,J,K) + dS_deta*InvGrad(2,1,I,J,K) +  &
                        dS_dzeta*InvGrad(2,2,I,J,K)
                END DO
            END DO
        END DO
    end subroutine physical_part_deriv_666

    subroutine physical_part_deriv_777(hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
        !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, parameter :: ngll=7
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: InvGrad
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Scalp
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: dS_dx,dS_dy,dS_dz
        real(fpp) :: dS_dxi, dS_deta, dS_dzeta
        integer :: i,j,k,l

        DO K = 0, ngll-1
            ! d(Scalp)_dxi
            DO J = 0, ngll-1
                DO I = 0, ngll-1
                    dS_dxi = ZERO
                    dS_deta = ZERO
                    dS_dzeta = ZERO
                    DO L = 0, ngll-1
                        dS_dxi   = dS_dxi   + Scalp(L,J,K)*hprime(L,I)
                        DS_deta  = dS_deta  + Scalp(I,L,K)*hprime(L,J)
                        dS_dzeta = dS_dzeta + Scalp(I,J,L)*hprime(L,K)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dxi*InvGrad(0,0,I,J,K) + dS_deta*InvGrad(0,1,I,J,K) +  &
                        dS_dzeta*InvGrad(0,2,I,J,K)
                    dS_dy(I,J,K) = dS_dxi*InvGrad(1,0,I,J,K) + dS_deta*InvGrad(1,1,I,J,K) +  &
                        dS_dzeta*InvGrad(1,2,I,J,K)
                    dS_dz(I,J,K) = dS_dxi*InvGrad(2,0,I,J,K) + dS_deta*InvGrad(2,1,I,J,K) +  &
                        dS_dzeta*InvGrad(2,2,I,J,K)
                END DO
            END DO
        END DO
    end subroutine physical_part_deriv_777

    subroutine physical_part_deriv_888(hprime,InvGrad,Scalp,dS_dx,dS_dy,dS_dz)
        !- partial derivatives of the scalar property Scalp, with respect to xi,eta,zeta
        implicit none
        integer, parameter :: ngll=8
        real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
        real(fpp), dimension(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: InvGrad
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: Scalp
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: dS_dx,dS_dy,dS_dz
        real(fpp) :: dS_dxi, dS_deta, dS_dzeta
        integer :: i,j,k,l

        DO K = 0, ngll-1
            ! d(Scalp)_dxi
            DO J = 0, ngll-1
                DO I = 0, ngll-1
                    dS_dxi = ZERO
                    dS_deta = ZERO
                    dS_dzeta = ZERO
                    DO L = 0, ngll-1
                        dS_dxi   = dS_dxi   + Scalp(L,J,K)*hprime(L,I)
                        DS_deta  = dS_deta  + Scalp(I,L,K)*hprime(L,J)
                        dS_dzeta = dS_dzeta + Scalp(I,J,L)*hprime(L,K)
                    END DO
                    !- in the physical domain
                    dS_dx(I,J,K) = dS_dxi*InvGrad(0,0,I,J,K) + dS_deta*InvGrad(0,1,I,J,K) +  &
                        dS_dzeta*InvGrad(0,2,I,J,K)
                    dS_dy(I,J,K) = dS_dxi*InvGrad(1,0,I,J,K) + dS_deta*InvGrad(1,1,I,J,K) +  &
                        dS_dzeta*InvGrad(1,2,I,J,K)
                    dS_dz(I,J,K) = dS_dxi*InvGrad(2,0,I,J,K) + dS_deta*InvGrad(2,1,I,J,K) +  &
                        dS_dzeta*InvGrad(2,2,I,J,K)
                END DO
            END DO
        END DO
    end subroutine physical_part_deriv_888
end module deriv3d

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
