!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
subroutine  initialize_material_prem(Tdomain, mat, elem, coorPt, npts)
    use constants, only : DM_SOLID_CG_PML, fpp
    use selement
    use ssubdomains
    use model_prem
    use tensor_util
    use earth_transform
    use sdomain
    use dom_solid
    use dom_solidpml
    implicit none
#include "index.h"

    type (domain), intent (INOUT), target :: Tdomain
    type(subdomain), intent(in) :: mat
    type(element), intent(inout) :: elem
    real(fpp), dimension(0:2,0:npts-1), intent(in) :: coorPt
    integer, intent(in) :: npts

    integer :: i,j,k,ii,jj,ngll,idef
    real(fpp) :: x, y, z, A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu, r, theta, phi
    real(fpp), dimension(1:6,1:6) :: Cij
    real(fpp), dimension(1:6,1:6,0:mat%ngll-1,0:mat%ngll-1,0:mat%ngll-1) :: gCij
    real(fpp), dimension(0:mat%ngll-1,0:mat%ngll-1,0:mat%ngll-1) :: lambda, mu, rho

    ngll = Tdomain%sdom%ngll

    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1

                idef = elem%Iglobnum(i,j,k)

                x = coorPt(0,idef)
                y = coorPt(1,idef)
                z = coorPt(2,idef)

                call cart2sph_ML(x, y, z, r, theta, phi)

                call get_value_prem (r, rho(i,j,k),A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu)

                Cij(:,:) = 0.d0
                Cij(1,1) = C
                Cij(2,2) = A+Bc+Ec
                Cij(3,3) = A-Bc+Ec
                Cij(1,2) = F+Hc
                Cij(1,3) = F-Hc
                Cij(2,3) = A-2.d0*M-Ec
                Cij(1,4) = s2*Hs
                Cij(2,4) = s2*(Bs/2.d0+Es)
                Cij(3,4) = s2*(Bs/2.d0-Es)
                Cij(4,4) = 2.d0*(M-Ec)
                Cij(5,5) = 2.d0*(L-Gc)
                Cij(6,6) = 2.d0*(L+Gc)
                Cij(5,6) = 2.d0*Gs

                do ii = 2,6
                    do jj = 1,ii-1
                        Cij(ii,jj) = Cij(jj,ii)
                    enddo
                enddo


                lambda(i,j,k) = lambda_from_Cij(Cij)
                mu(i,j,k) = mu_from_Cij(Cij)
                if(elem%domain==DM_SOLID_CG) then
                    call c_4tensor(Cij,theta,phi)
                endif
                gCij(:,:,i,j,k) = Cij
            end do
        end do
    end do
    if(elem%domain==DM_SOLID_CG_PML) then
        call init_material_properties_solidpml(Tdomain%spmldom,elem%lnum,mat,&
            rho,lambda,mu)
    else if(elem%domain==DM_SOLID_CG) then
        call init_material_tensor_solid(Tdomain%sdom,elem%lnum,mat,rho,lambda,mu,lambda,mu,gCij)
    else
        stop "initialize material prem KO"
    endif

end subroutine initialize_material_prem



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
