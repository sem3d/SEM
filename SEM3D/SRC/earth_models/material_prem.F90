!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
subroutine  initialize_material_prem(Tdomain, elem, coorPt, npts)
    use constants, only : DM_SOLID_PML
    use selement
    use ssubdomains
    use model_prem
    use tensor_util
    use earth_transform
    use sdomain
    implicit none
#include "index.h"

    type (domain), intent (INOUT), target :: Tdomain
    type(element), intent(inout) :: elem
    real, dimension(0:2,0:npts-1), intent(in) :: coorPt
    integer, intent(in) :: npts


    integer :: i,j,k,ii,jj, ngllx, nglly, ngllz, idef
    real :: x, y, z, rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu, r, theta, phi
    real, dimension(1:6,1:6) :: Cij

    ngllx = elem%ngllx
    nglly = elem%nglly
    ngllz = elem%ngllz


    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1

                idef = elem%Iglobnum(i,j,k)

                x = coorPt(0,idef)
                y = coorPt(1,idef)
                z = coorPt(2,idef)

                call cart2sph(x, y, z, r, theta, phi)

                call get_value_prem (r, rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu)

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


                if(elem%domain==DM_SOLID_PML) then
                    Tdomain%spmldom%Lambda_ (i,j,k,elem%lnum) = lambda_from_Cij(Cij)
                    Tdomain%spmldom%Mu_     (i,j,k,elem%lnum) = mu_from_Cij(Cij)
                    Tdomain%spmldom%Density_(i,j,k,elem%lnum) = rho
                else ! DM_SOLID
                    call c_4tensor(Cij,theta,phi)
                    idef = 0
                    do ii = 1,6
                        do jj = ii,6
                            Tdomain%sdom%Cij_(idef,i,j,k,elem%lnum) = Cij(ii,jj)
                            idef = idef + 1
                        enddo
                    enddo
                    Tdomain%sdom%Density_(i,j,k,elem%lnum) = rho
                endif
            end do
        end do
    end do

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
