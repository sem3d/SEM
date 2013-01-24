!------------------------------------------------------------------------
real function interp_lag(Tdomain,ngllx,nglly,ngllz,mat,xi,eta,zeta,func)
    ! gives the Lagrange interpolation at (xi,eta,zeta) for a function func
    !    whose values are known at GLL points
    use sdomain
    implicit none

    type(domain), intent(in)   :: Tdomain
    integer, intent(in)   :: mat,ngllx,nglly,ngllz
    real, intent(in)  :: xi,eta,zeta
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: func
    integer :: i,j,k
    real  :: f,wxi,weta,wzeta

    f = 0d0

    do k = 0,ngllz-1
        call pol_lagrange(ngllz,Tdomain%sSubdomain(mat)%GLLcz,k,zeta,wzeta)
        do j = 0,nglly-1
            call pol_lagrange(nglly,Tdomain%sSubdomain(mat)%GLLcy,j,eta,weta)
            do i = 0,ngllx-1
                call pol_lagrange(ngllx,Tdomain%sSubdomain(mat)%GLLcx,i,xi,wxi)
                f = f+func(i,j,k)*wxi*weta*wzeta
            end do
        end do
    end do

    interp_lag = f

    return
end function interp_lag


subroutine source_excit(Tdomain,rank)
    ! gives the excitation coeff. at GLL in an element, for a given source
    use mpi
    use sdomain
    implicit none

    interface
       real function interp_lag(Tdomain,ngllx,nglly,ngllz,mat,xi,eta,zeta,func)
           ! gives the Lagrange interpolation at (xi,eta,zeta) for a function func
           !    whose values are known at GLL points
           use sdomain
           implicit none

           type(domain), intent(in)   :: Tdomain
           integer, intent(in)   :: mat,ngllx,nglly,ngllz
           real, intent(in)  :: xi,eta,zeta
           real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: func
       end function interp_lag
    end interface

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)  :: rank
    integer :: nsour,nels,mat,ngllx,nglly,ngllz,i,j,k
    real :: xi,eta,zeta,wxi,weta,wzeta,rho,dwdxi,dwdeta,dwdzeta
    real, dimension(0:2,0:2) :: M,InvGrad


    do nsour = 0, Tdomain%n_source -1
        if(Tdomain%sSource(nsour)%proc == rank)then
            nels = Tdomain%sSource(nsour)%elem
            mat = Tdomain%specel(nels)%mat_index
            ngllx = Tdomain%specel(nels)%ngllx
            nglly = Tdomain%specel(nels)%nglly
            ngllz = Tdomain%specel(nels)%ngllz
            xi = Tdomain%sSource(nsour)%RefCoord(0)
            eta = Tdomain%sSource(nsour)%RefCoord(1)
            zeta = Tdomain%sSource(nsour)%RefCoord(2)

            ! pulse in a solid (in a given direction) or fluid (isotropic pressure source term)
            if(Tdomain%sSource(nsour)%i_type_source == 1 .or.  &
                Tdomain%sSource(nsour)%i_type_source == 3) then
                allocate(Tdomain%sSource(nsour)%ExtForce(0:ngllx-1,0:nglly-1,0:ngllz-1))
                do k = 0,ngllz-1
                    call pol_lagrange(ngllz,Tdomain%sSubdomain(mat)%GLLcz,k,zeta,wzeta)
                    do j = 0,nglly-1
                        call pol_lagrange(nglly,Tdomain%sSubdomain(mat)%GLLcy,j,eta,weta)
                        do i = 0,ngllx-1
                            call pol_lagrange(ngllx,Tdomain%sSubdomain(mat)%GLLcx,i,xi,wxi)
                            Tdomain%sSource(nsour)%ExtForce(i,j,k) = wxi*weta*wzeta
                        end do
                    end do
                end do
                ! fluid case
                if(.not. Tdomain%specel(nels)%solid)then
                    rho = interp_lag(Tdomain,ngllx,nglly,ngllz,mat,xi,eta,zeta,   &
                        Tdomain%specel(nels)%density)
                    Tdomain%sSource(nsour)%ExtForce(:,:,:) = -Tdomain%sSource(nsour)%ExtForce(:,:,:)/rho
                endif
                ! point source = moment tensor M (explosion is a special case: M(i,j) = delta _(ij))
            else if(Tdomain%sSource(nsour)%i_type_source == 2)then
                allocate(Tdomain%sSource(nsour)%MomForce(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                M(:,:) = Tdomain%sSource(nsour)%Moment(:,:)
                InvGrad(:,:) = Tdomain%sSource(nsour)%InvGrad(:,:)
                do k = 0,ngllz-1
                    call pol_lagrange(ngllz,Tdomain%sSubdomain(mat)%GLLcz,k,zeta,wzeta)
                    call derivlag(Tdomain%sSubdomain(mat)%GLLcz,ngllz,k,zeta,dwdzeta)
                    do j = 0,nglly-1
                        call pol_lagrange(nglly,Tdomain%sSubdomain(mat)%GLLcy,j,eta,weta)
                        call derivlag(Tdomain%sSubdomain(mat)%GLLcy,nglly,j,eta,dwdeta)
                        do i = 0,ngllx-1
                            call pol_lagrange(ngllx,Tdomain%sSubdomain(mat)%GLLcx,i,xi,wxi)
                            call derivlag(Tdomain%sSubdomain(mat)%GLLcx,ngllx,i,xi,dwdxi)
                            ! final step: values at the GLL points
                            Tdomain%sSource(nsour)%MomForce(i,j,k,0) =     &
                                (InvGrad(0,0)*dwdxi*weta*wzeta+InvGrad(0,1)*wxi*dwdeta*wzeta+InvGrad(0,2)*wxi*weta*dwdzeta)*M(0,0) + &
                                (InvGrad(1,0)*dwdxi*weta*wzeta+InvGrad(1,1)*wxi*dwdeta*wzeta+InvGrad(1,2)*wxi*weta*dwdzeta)*M(0,1) + &
                                (InvGrad(2,0)*dwdxi*weta*wzeta+InvGrad(2,1)*wxi*dwdeta*wzeta+InvGrad(2,2)*wxi*weta*dwdzeta)*M(0,2)
                            Tdomain%sSource(nsour)%MomForce(i,j,k,1) =     &
                                (InvGrad(0,0)*dwdxi*weta*wzeta+InvGrad(0,1)*wxi*dwdeta*wzeta+InvGrad(0,2)*wxi*weta*dwdzeta)*M(1,0) + &
                                (InvGrad(1,0)*dwdxi*weta*wzeta+InvGrad(1,1)*wxi*dwdeta*wzeta+InvGrad(1,2)*wxi*weta*dwdzeta)*M(1,1) + &
                                (InvGrad(2,0)*dwdxi*weta*wzeta+InvGrad(2,1)*wxi*dwdeta*wzeta+InvGrad(2,2)*wxi*weta*dwdzeta)*M(1,2)
                            Tdomain%sSource(nsour)%MomForce(i,j,k,2) =     &
                                (InvGrad(0,0)*dwdxi*weta*wzeta+InvGrad(0,1)*wxi*dwdeta*wzeta+InvGrad(0,2)*wxi*weta*dwdzeta)*M(2,0) + &
                                (InvGrad(1,0)*dwdxi*weta*wzeta+InvGrad(1,1)*wxi*dwdeta*wzeta+InvGrad(1,2)*wxi*weta*dwdzeta)*M(2,1) + &
                                (InvGrad(2,0)*dwdxi*weta*wzeta+InvGrad(2,1)*wxi*dwdeta*wzeta+InvGrad(2,2)*wxi*weta*dwdzeta)*M(2,2)
                        end do
                    enddo
                enddo

            end if  ! end i_type_source

        end if
    end do  ! end of loop over sources

end subroutine source_excit

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
