subroutine  initialize_material_earthchunk( elem, matInfo, coorPt, npts)
    use selement
    use ssubdomains
    use read_model_earthchunk
    use tensor_util


    type(element), intent(inout) :: elem
    type(subdomain), intent(in) :: matInfo
    real, dimension(0:2,0:npts-1), intent(in) :: coorPt
    integer, intent(in) :: npts


    integer :: i,j,k,ii,jj, ngllx, nglly, ngllz, idef
    real :: x, y, z, rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu, r, theta, phi, lon, lat
    real, dimension(1:6,1:6) :: Cij

    ngllx = elem%ngllx
    nglly = elem%nglly
    ngllz = elem%ngllz

 !   write(*,*) "--material earthchunk--"

    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1

                idef = elem%Iglobnum(i,j,k)

                x = coorPt(0,idef)
                y = coorPt(1,idef)
                z = coorPt(2,idef)

                call cart2sph(x, y, z, r, lon, lat)
                call get_value_aniso (r, lon, lat, rho,A,C,F,L,M,Gc,Gs,Hc,Hs,Bc,Bs,Ec,Es,Qmu)
                theta = (90.0-lat)*Pi180
                phi = lon*Pi180

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


                if( elem%PML) then
                    elem%Lambda(i,j,k) = lambda_from_Cij(Cij)
                    elem%Mu(i,j,k) = mu_from_Cij(Cij)
                else
                    
                    call c_4tensor(Cij,theta,phi)

                    idef = 0
                    do ii = 1,6
                        do jj = ii,6
                            elem%sl%Cij(idef,i,j,k) = Cij(ii,jj)
                            idef = idef + 1
                        enddo
                    enddo
                endif

                elem%Density(i,j,k) = rho

            end do
        end do
    end do

end subroutine initialize_material_earthchunk
