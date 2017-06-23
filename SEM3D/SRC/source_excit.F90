!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!------------------------------------------------------------------------
module msource_excit
    use constants, only : fpp
    implicit none
contains
    real(fpp) function interp_lag(ngll,GLLc,xi,eta,zeta,func)
        ! gives the Lagrange interpolation at (xi,eta,zeta) for a function func
        !    whose values are known at GLL points
        implicit none

        integer, intent(in) :: ngll
        real(fpp), intent(in), dimension(0:ngll-1) :: GLLc
        real(fpp), intent(in) :: xi,eta,zeta
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1) :: func
        !
        integer :: i,j,k
        real(fpp)  :: f,wxi,weta,wzeta

        f = 0d0
        do k = 0,ngll-1
            call pol_lagrange(ngll,GLLc,k,zeta,wzeta)
            do j = 0,ngll-1
                call pol_lagrange(ngll,GLLc,j,eta,weta)
                do i = 0,ngll-1
                    call pol_lagrange(ngll,GLLc,i,xi,wxi)
                    f = f+func(i,j,k)*wxi*weta*wzeta
                end do
            end do
        end do
        interp_lag = f
    end function interp_lag

    subroutine source_excit_pulse(src,ngll,GLLc)
        use ssources
        implicit none
        type(source), intent(inout) :: src
        integer, intent(in) :: ngll
        real(fpp), dimension(0:ngll-1), intent(in) :: GLLc
        !
        integer :: np, i, j, k
        real(fpp) :: xi,eta,zeta,wxi,weta,wzeta

        xi = src%RefCoord(0)
        eta = src%RefCoord(1)
        zeta = src%RefCoord(2)

        allocate(src%ExtForce(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        do k = 0,ngll-1
            call pol_lagrange(ngll, GLLc, k, zeta, wzeta)
            do j = 0,ngll-1
                call pol_lagrange(ngll, GLLc, j, eta, weta)
                do i = 0,ngll-1
                    call pol_lagrange(ngll, GLLc, i, xi, wxi)
                    do np = 0,2
                        src%ExtForce(i,j,k,np) = wxi*weta*wzeta*src%dir(np)
                    end do
                end do
            end do
        end do
    end subroutine source_excit_pulse

    subroutine source_excit_pulse_fluid(Tdomain, nels, src, ngll, GLLc)
        use constants, only : fpp
        use ssources
        use sdomain
        implicit none
#include "index.h"
        type(domain), intent(in)   :: Tdomain
        integer, intent(in) :: nels
        type(source), intent(inout) :: src
        integer, intent(in) :: ngll
        real(fpp), dimension(0:ngll-1), intent(in) :: GLLc
        !
        integer :: i, j, k, bnum, ee
        real(fpp) :: xi,eta,zeta,wxi,weta,wzeta,lambda


        bnum = Tdomain%specel(nels)%lnum/VCHUNK
        ee = mod(Tdomain%specel(nels)%lnum,VCHUNK)

        xi = src%RefCoord(0)
        eta = src%RefCoord(1)
        zeta = src%RefCoord(2)

        allocate(src%ExtForce(0:ngll-1,0:ngll-1,0:ngll-1,0:0))
        do k = 0,ngll-1
            call pol_lagrange(ngll, GLLc, k, zeta, wzeta)
            do j = 0,ngll-1
                call pol_lagrange(ngll, GLLc, j, eta, weta)
                do i = 0,ngll-1
                    call pol_lagrange(ngll, GLLc, i, xi, wxi)
                    src%ExtForce(i,j,k,0) = wxi*weta*wzeta
                end do
            end do
        end do

        ! fluid case
        lambda = 0.
        if(Tdomain%specel(nels)%domain==DM_FLUID) then
            lambda = interp_lag(Tdomain%fdom%ngll,Tdomain%fdom%GLLc,xi,eta,zeta,&
                Tdomain%fdom%Lambda_(:,:,:,bnum,ee))
        end if
        if(Tdomain%specel(nels)%domain==DM_FLUID_PML) then
            lambda = interp_lag(Tdomain%fpmldom%ngll,Tdomain%fpmldom%GLLc,xi,eta,zeta,&
                Tdomain%fpmldom%Lambda_(:,:,:,bnum,ee))
        end if
        src%ExtForce(:,:,:,0) = -src%ExtForce(:,:,:,0)/lambda
        ! point source = moment tensor M (explosion is a special case: M(i,j) = delta _(ij))

    end subroutine source_excit_pulse_fluid
    ! TAG_MOMENT
    subroutine source_excit_moment(src, ngll, GLLc)
        use ssources
        implicit none
        type(source), intent(inout) :: src
        integer, intent(in) :: ngll
        real(fpp), dimension(0:ngll-1), intent(in) :: GLLc
        !
        integer :: i, j, k
        real(fpp) :: xi,eta,zeta,wxi,weta,wzeta,dwdxi,dwdeta,dwdzeta
        real(fpp), dimension(0:2,0:2) :: M,InvGrad

        xi = src%RefCoord(0)
        eta = src%RefCoord(1)
        zeta = src%RefCoord(2)

        allocate(src%ExtForce(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        M(:,:) = src%Moment(:,:)
        

        InvGrad(:,:) = src%InvGrad(:,:)
        do k = 0,ngll-1
            call pol_lagrange(ngll,GLLc,k,zeta,wzeta)
            call der_lagrange(ngll,GLLc,k,zeta,dwdzeta)
            do j = 0,ngll-1
                call pol_lagrange(ngll,GLLc,j,eta,weta)
                call der_lagrange(ngll,GLLc,j,eta,dwdeta)
                do i = 0,ngll-1
                    call pol_lagrange(ngll,GLLc,i,xi,wxi)
                    call der_lagrange(ngll,GLLc,i,xi,dwdxi)
                    ! final step: values at the GLL points
                    src%ExtForce(i,j,k,0) =     &
                        (InvGrad(0,0)*dwdxi*weta*wzeta+InvGrad(0,1)*wxi*dwdeta*wzeta+InvGrad(0,2)*wxi*weta*dwdzeta)*M(0,0) + &
                        (InvGrad(1,0)*dwdxi*weta*wzeta+InvGrad(1,1)*wxi*dwdeta*wzeta+InvGrad(1,2)*wxi*weta*dwdzeta)*M(0,1) + &
                        (InvGrad(2,0)*dwdxi*weta*wzeta+InvGrad(2,1)*wxi*dwdeta*wzeta+InvGrad(2,2)*wxi*weta*dwdzeta)*M(0,2)
                    src%ExtForce(i,j,k,1) =     &
                        (InvGrad(0,0)*dwdxi*weta*wzeta+InvGrad(0,1)*wxi*dwdeta*wzeta+InvGrad(0,2)*wxi*weta*dwdzeta)*M(1,0) + &
                        (InvGrad(1,0)*dwdxi*weta*wzeta+InvGrad(1,1)*wxi*dwdeta*wzeta+InvGrad(1,2)*wxi*weta*dwdzeta)*M(1,1) + &
                        (InvGrad(2,0)*dwdxi*weta*wzeta+InvGrad(2,1)*wxi*dwdeta*wzeta+InvGrad(2,2)*wxi*weta*dwdzeta)*M(1,2)
                    src%ExtForce(i,j,k,2) =     &
                        (InvGrad(0,0)*dwdxi*weta*wzeta+InvGrad(0,1)*wxi*dwdeta*wzeta+InvGrad(0,2)*wxi*weta*dwdzeta)*M(2,0) + &
                        (InvGrad(1,0)*dwdxi*weta*wzeta+InvGrad(1,1)*wxi*dwdeta*wzeta+InvGrad(1,2)*wxi*weta*dwdzeta)*M(2,1) + &
                        (InvGrad(2,0)*dwdxi*weta*wzeta+InvGrad(2,1)*wxi*dwdeta*wzeta+InvGrad(2,2)*wxi*weta*dwdzeta)*M(2,2)
                end do
            enddo
        enddo
    end subroutine source_excit_moment

    subroutine source_excit(Tdomain,rank)
        ! gives the excitation coeff. at GLL in an element, for a given source
        use mpi
        use sdomain
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer, intent(in)  :: rank
        integer :: nsour,nels,ngll
        real(fpp), dimension(:), allocatable :: GLLc


        do nsour = 0, Tdomain%n_source -1
            if(Tdomain%sSource(nsour)%proc == rank)then
                nels = Tdomain%sSource(nsour)%elem
                select case (Tdomain%specel(nels)%domain)
                case (DM_SOLID)
                    ngll = Tdomain%sdom%ngll
                    allocate(GLLc(0:ngll-1))
                    GLLc = Tdomain%sdom%GLLc
                case (DM_FLUID)
                    ngll = Tdomain%fdom%ngll
                    allocate(GLLc(0:ngll-1))
                    GLLc = Tdomain%fdom%GLLc
                case (DM_SOLID_PML)
                    ngll = Tdomain%spmldom%ngll
                    allocate(GLLc(0:ngll-1))
                    GLLc = Tdomain%spmldom%GLLc
                case (DM_FLUID_PML)
                    ngll = Tdomain%fpmldom%ngll
                    allocate(GLLc(0:ngll-1))
                    GLLc = Tdomain%fpmldom%GLLc
                case default
                    stop "unknown domain"
                end select

                ! pulse in a solid (in a given direction) or fluid (isotropic pressure source term)
                if(Tdomain%sSource(nsour)%i_type_source == 1) then
                    call source_excit_pulse(Tdomain%sSource(nsour), ngll, GLLc)
                else if (Tdomain%sSource(nsour)%i_type_source == 3) then
                    call source_excit_pulse_fluid(Tdomain, nels, Tdomain%sSource(nsour), ngll, GLLc)
                else if(Tdomain%sSource(nsour)%i_type_source == 2)then
                    call source_excit_moment(Tdomain%sSource(nsour), ngll, GLLc)

                end if  ! end i_type_source
                deallocate(GLLc)
            end if
        end do  ! end of loop over sources

    end subroutine source_excit

end module msource_excit
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
