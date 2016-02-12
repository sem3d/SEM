!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!------------------------------------------------------------------------
real function interp_lag(mat,xi,eta,zeta,func)
    ! gives the Lagrange interpolation at (xi,eta,zeta) for a function func
    !    whose values are known at GLL points
    use ssubdomains
    implicit none

    type(subdomain), intent(in)   :: mat
    real, intent(in)  :: xi,eta,zeta
    real, dimension(0:mat%ngll-1,0:mat%ngll-1,0:mat%ngll-1) :: func
    integer :: i,j,k
    real  :: f,wxi,weta,wzeta
    integer  :: ngll
    ngll = mat%ngll

    f = 0d0
    do k = 0,ngll-1
        call pol_lagrange(ngll,mat%GLLcz,k,zeta,wzeta)
        do j = 0,ngll-1
            call pol_lagrange(ngll,mat%GLLcy,j,eta,weta)
            do i = 0,ngll-1
                call pol_lagrange(ngll,mat%GLLcx,i,xi,wxi)
                f = f+func(i,j,k)*wxi*weta*wzeta
            end do
        end do
    end do
    interp_lag = f
end function interp_lag

subroutine source_excit_pulse(src,mat)
    use ssources
    use ssubdomains
    implicit none
    type(source), intent(inout) :: src
    type(subdomain), intent(in) :: mat
    integer :: ngll, np, i, j, k
    real :: xi,eta,zeta,wxi,weta,wzeta
    ngll = mat%ngll
    xi = src%RefCoord(0)
    eta = src%RefCoord(1)
    zeta = src%RefCoord(2)

    allocate(src%ExtForce(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
    do k = 0,ngll-1
        call pol_lagrange(ngll, mat%GLLcz, k, zeta, wzeta)
        do j = 0,ngll-1
            call pol_lagrange(ngll,mat%GLLcy,j,eta,weta)
            do i = 0,ngll-1
                call pol_lagrange(ngll,mat%GLLcx,i,xi,wxi)
                do np = 0,2
                    src%ExtForce(i,j,k,np) = wxi*weta*wzeta*src%dir(np)
                end do
            end do
        end do
    end do
end subroutine source_excit_pulse

subroutine source_excit_pulse_fluid(Tdomain, nels, src,mat)
    use ssources
    use sdomain
    use ssubdomains
    implicit none
#include "index.h"
    type(source), intent(inout) :: src
    type(subdomain), intent(in) :: mat
    type(domain), intent(in)   :: Tdomain
    integer, intent(in) :: nels
    integer :: ngll, i, j, k
    real :: xi,eta,zeta,wxi,weta,wzeta,lambda

    interface
       real function interp_lag(mat,xi,eta,zeta,func)
           ! gives the Lagrange interpolation at (xi,eta,zeta) for a function func
           !    whose values are known at GLL points
           use ssubdomains
           implicit none

           type(subdomain), intent(in)   :: mat
           real, intent(in)  :: xi,eta,zeta
           real, dimension(0:mat%ngll-1,0:mat%ngll-1,0:mat%ngll-1) :: func
       end function interp_lag
    end interface

    ngll = mat%ngll
    xi = src%RefCoord(0)
    eta = src%RefCoord(1)
    zeta = src%RefCoord(2)

    allocate(src%ExtForce(0:ngll-1,0:ngll-1,0:ngll-1,0:0))
    do k = 0,ngll-1
        call pol_lagrange(ngll, mat%GLLcz, k, zeta, wzeta)
        do j = 0,ngll-1
            call pol_lagrange(ngll,mat%GLLcy,j,eta,weta)
            do i = 0,ngll-1
                call pol_lagrange(ngll,mat%GLLcx,i,xi,wxi)
                src%ExtForce(i,j,k,0) = wxi*weta*wzeta
            end do
        end do
    end do

    ! fluid case
    lambda = 0.
    if(Tdomain%specel(nels)%domain==DM_FLUID) then
        lambda = interp_lag(mat,xi,eta,zeta,Tdomain%fdom%Lambda_(:,:,:,Tdomain%specel(nels)%lnum))
    end if
    if(Tdomain%specel(nels)%domain==DM_FLUID_PML) then
        lambda = interp_lag(mat,xi,eta,zeta,Tdomain%fpmldom%Lambda_(:,:,:,Tdomain%specel(nels)%lnum))
    end if
    src%ExtForce(:,:,:,0) = -src%ExtForce(:,:,:,0)/lambda
    ! point source = moment tensor M (explosion is a special case: M(i,j) = delta _(ij))

end subroutine source_excit_pulse_fluid

subroutine source_excit_moment(src,mat)
    use ssources
    use ssubdomains
    implicit none
    type(source), intent(inout) :: src
    type(subdomain), intent(in) :: mat
    integer :: ngll, i, j, k
    real :: xi,eta,zeta,wxi,weta,wzeta,dwdxi,dwdeta,dwdzeta
    real, dimension(0:2,0:2) :: M,InvGrad
    ngll = mat%ngll
    xi = src%RefCoord(0)
    eta = src%RefCoord(1)
    zeta = src%RefCoord(2)

    allocate(src%ExtForce(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
    M(:,:) = src%Moment(:,:)
    InvGrad(:,:) = src%InvGrad(:,:)
    do k = 0,ngll-1
        call pol_lagrange(ngll,mat%GLLcz,k,zeta,wzeta)
        call derivlag(mat%GLLcz,ngll,k,zeta,dwdzeta)
        do j = 0,ngll-1
            call pol_lagrange(ngll,mat%GLLcy,j,eta,weta)
            call derivlag(mat%GLLcy,ngll,j,eta,dwdeta)
            do i = 0,ngll-1
                call pol_lagrange(ngll,mat%GLLcx,i,xi,wxi)
                call derivlag(mat%GLLcx,ngll,i,xi,dwdxi)
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
    integer :: nsour,nels,mat


    do nsour = 0, Tdomain%n_source -1
        if(Tdomain%sSource(nsour)%proc == rank)then
            nels = Tdomain%sSource(nsour)%elem
            mat = Tdomain%specel(nels)%mat_index

            ! pulse in a solid (in a given direction) or fluid (isotropic pressure source term)
            if(Tdomain%sSource(nsour)%i_type_source == 1) then
                call source_excit_pulse(Tdomain%sSource(nsour), Tdomain%sSubdomain(mat))

            else if (Tdomain%sSource(nsour)%i_type_source == 3) then
                call source_excit_pulse_fluid(Tdomain, nels, Tdomain%sSource(nsour), Tdomain%sSubdomain(mat))
            else if(Tdomain%sSource(nsour)%i_type_source == 2)then
                call source_excit_moment(Tdomain%sSource(nsour), Tdomain%sSubdomain(mat))
            end if  ! end i_type_source

        end if
    end do  ! end of loop over sources

end subroutine source_excit


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
