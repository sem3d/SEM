!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_solid
    use sdomain
    use constants
    use champs_solid
    use selement
    use ssubdomains
    implicit none

contains

    subroutine get_solid_dom_var(Tdomain, el, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain)                               :: TDomain
        integer, dimension(0:8)                    :: out_variables
        type(element)                              :: el
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp)                                  :: P_energy, S_energy, eps_vol
        real(fpp), dimension(0:5)                  :: eps_dev
        real(fpp), dimension(0:5)                  :: sig_dev
        !
        logical                                  :: flag_gradU
        integer                                  :: nx, ny, nz, i, j, k, ind, mat
        real(fpp), dimension(:,:,:), allocatable :: DXX, DXY, DXZ
        real(fpp), dimension(:,:,:), allocatable :: DYX, DYY, DYZ
        real(fpp), dimension(:,:,:), allocatable :: DZX, DZY, DZZ
        real(fpp), dimension (:,:), allocatable  :: htprimex, hprimey, hprimez
        real(fpp)                                :: xmu, xlambda, xkappa
        real(fpp)                                :: x2mu, xlambda2mu
        real(fpp)                                :: onemSbeta, onemPbeta

        flag_gradU = (out_variables(OUT_PRESSION)    + &
                      out_variables(OUT_ENERGYP)     + &
                      out_variables(OUT_ENERGYS)     + &
                      out_variables(OUT_EPS_VOL)     + &
                      out_variables(OUT_EPS_DEV)     + &
                      out_variables(OUT_STRESS_DEV)) /= 0

        nx = el%ngllx
        ny = el%nglly
        nz = el%ngllz
        mat = el%mat_index

        if (flag_gradU) then
            allocate(DXX(0:nx-1,0:ny-1,0:nz-1))
            allocate(DXY(0:nx-1,0:ny-1,0:nz-1))
            allocate(DXZ(0:nx-1,0:ny-1,0:nz-1))
            allocate(DYX(0:nx-1,0:ny-1,0:nz-1))
            allocate(DYY(0:nx-1,0:ny-1,0:nz-1))
            allocate(DYZ(0:nx-1,0:ny-1,0:nz-1))
            allocate(DZX(0:nx-1,0:ny-1,0:nz-1))
            allocate(DZY(0:nx-1,0:ny-1,0:nz-1))
            allocate(DZZ(0:nx-1,0:ny-1,0:nz-1))
            allocate(hTprimex(0:nx-1,0:nx-1))
            allocate(hprimey(0:ny-1,0:ny-1))
            allocate(hprimez(0:nz-1,0:nz-1))
            hTprimex=Tdomain%sSubDomain(mat)%hTprimex
            hprimey=Tdomain%sSubDomain(mat)%hprimey
            hprimez=Tdomain%sSubDomain(mat)%hprimez
        end if

        if (.not. Tdomain%aniso) then
            xmu     = el%Mu(i,j,k)
            xlambda = el%Lambda(i,j,k)
            xkappa  = el%Kappa(i,j,k)

            if (Tdomain%n_sls>0) then
                onemSbeta = el%sl%onemSbeta(i,j,k)
                onemPbeta = el%sl%onemPbeta(i,j,k)
                xmu    = xmu * onemSbeta
                xkappa = xkappa * onemPbeta
            endif
            x2mu       = 2. * xmu
            xlambda2mu = xlambda + x2mu
        end if

        do k=0,nz-1
            do j=0,ny-1
                do i=0,nx-1
                    ind = el%Idom(i,j,k)

                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        if(.not. allocated(fieldU)) allocate(fieldU(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldU(i,j,k,:) = Tdomain%sdom%champs0%Depla(ind,:)
                    end if

                enddo
            enddo
        enddo

        if (flag_gradU) then
            call physical_part_deriv(nx,ny,nz,htprimex,hprimey,hprimez,el%InvGrad,fieldU(:,:,:,0),DXX,DYX,DZX)
            call physical_part_deriv(nx,ny,nz,htprimex,hprimey,hprimez,el%InvGrad,fieldU(:,:,:,1),DXY,DYY,DZY)
            call physical_part_deriv(nx,ny,nz,htprimex,hprimey,hprimez,el%InvGrad,fieldU(:,:,:,2),DXZ,DYZ,DZZ)
        end if

        do k=0,nz-1
            do j=0,ny-1
                do i=0,nx-1
                    ind = el%Idom(i,j,k)

                    if (out_variables(OUT_VITESSE) == 1) then
                        if(.not. allocated(fieldV)) allocate(fieldV(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldV(i,j,k,:) = Tdomain%sdom%champs0%Veloc(ind,:)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldA(i,j,k,:) = Tdomain%sdom%champs0%Forces(ind,:)
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:nx-1,0:ny-1,0:nz-1))
                        fieldP(i,j,k) = -(el%lambda(i,j,k)+2d0/3d0*el%mu(i,j,k))*(DXX(i,j,k)+DYY(i,j,k)+DZZ(i,j,k))
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        eps_vol = DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k)
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        P_energy = 0.
                        if (.not. Tdomain%aniso) then
                            P_energy = .5 * xlambda2mu * eps_vol**2
                        end if
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        S_energy = 0.
                        if (.not. Tdomain%aniso) then
                            S_energy =   xmu/2 * ( DXY(i,j,k)**2 + DYX(i,j,k)**2 &
                                       +   DXZ(i,j,k)**2 + DZX(i,j,k)**2 &
                                       +   DYZ(i,j,k)**2 + DZY(i,j,k)**2 &
                                       - 2 * DXY(i,j,k) * DYX(i,j,k)     &
                                       - 2 * DXZ(i,j,k) * DZX(i,j,k)     &
                                       - 2 * DYZ(i,j,k) * DZY(i,j,k))
                        end if
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        eps_dev(0) = DXX(i,j,k) - eps_vol / 3
                        eps_dev(1) = DYY(i,j,k) - eps_vol / 3
                        eps_dev(2) = DZZ(i,j,k) - eps_vol / 3
                        eps_dev(3) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
                        eps_dev(4) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
                        eps_dev(5) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        sig_dev = 0.
                        if (.not. Tdomain%aniso) then
                            sig_dev(0) = x2mu * (DXX(i,j,k) - eps_vol * M_1_3)
                            sig_dev(1) = x2mu * (DYY(i,j,k) - eps_vol * M_1_3)
                            sig_dev(2) = x2mu * (DZZ(i,j,k) - eps_vol * M_1_3)
                            sig_dev(3) = xmu * (DXY(i,j,k) + DYX(i,j,k))
                            sig_dev(4) = xmu * (DXZ(i,j,k) + DZX(i,j,k))
                            sig_dev(5) = xmu * (DYZ(i,j,k) + DZY(i,j,k))
                        end if
                    end if
                enddo
            enddo
        enddo

        if (allocated(DXX))      deallocate(DXX)
        if (allocated(DXY))      deallocate(DXY)
        if (allocated(DXZ))      deallocate(DXZ)
        if (allocated(DYX))      deallocate(DYX)
        if (allocated(DYY))      deallocate(DYY)
        if (allocated(DYZ))      deallocate(DYZ)
        if (allocated(DZX))      deallocate(DZX)
        if (allocated(DZY))      deallocate(DZY)
        if (allocated(DZZ))      deallocate(DZZ)
        if (allocated(hTprimex)) deallocate(hTprimex)
        if (allocated(hprimey))  deallocate(hprimey)
        if (allocated(hprimez))  deallocate(hprimez)
    end subroutine get_solid_dom_var

    subroutine forces_int_solid(Elem, mat, htprimex, hprimey, htprimey, hprimez, htprimez,  &
               n_solid, aniso, champs1)

        use attenuation_solid

        type (Element), intent (INOUT) :: Elem
        type (subdomain), intent(IN) :: mat
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: htprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        integer, intent(IN) :: n_solid
        logical, intent(IN) :: aniso
        type(champssolid), intent(inout) :: champs1

        integer :: m1,m2,m3, i,j,k,i_dir
        real :: epsilon_trace_over_3
        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) ::  DXX,DXY,DXZ, &
            DYX,DYY,DYZ, &
            DZX,DZY,DZZ, &
            Fox,Foy,Foz

        real, dimension(:,:,:), allocatable :: epsilondev_xx_loc, epsilondev_yy_loc, &
            epsilondev_xy_loc, epsilondev_xz_loc, epsilondev_yz_loc
        real, dimension(:,:,:), allocatable :: epsilonvol_loc
        real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1,0:2) :: Depla


        m1 = Elem%ngllx;   m2 = Elem%nglly;   m3 = Elem%ngllz

        do i_dir = 0,2
            do k = 0,m3-1
                do j = 0,m2-1
                    do i = 0,m1-1
                        Depla(i,j,k,i_dir) = champs1%Depla(Elem%Idom(i,j,k),i_dir)
                    enddo
                enddo
            enddo
        enddo

        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Depla(:,:,:,0),dxx,dyx,dzx)
        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Depla(:,:,:,1),dxy,dyy,dzy)
        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad,Depla(:,:,:,2),dxz,dyz,dzz)

        if (n_solid>0) then
            if (aniso) then
            else
                allocate (epsilonvol_loc(0:m1-1,0:m2-1,0:m3-1))
            endif
            allocate (epsilondev_xx_loc(0:m1-1,0:m2-1,0:m3-1))
            allocate (epsilondev_yy_loc(0:m1-1,0:m2-1,0:m3-1))
            allocate (epsilondev_xy_loc(0:m1-1,0:m2-1,0:m3-1))
            allocate (epsilondev_xz_loc(0:m1-1,0:m2-1,0:m3-1))
            allocate (epsilondev_yz_loc(0:m1-1,0:m2-1,0:m3-1))
            do i = 0,m1-1
                do j = 0,m2-1
                    do k = 0,m3-1
                        epsilon_trace_over_3 = 0.333333333333333333333333333333d0 * (DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k))
                        if (aniso) then
                        else
                            epsilonvol_loc(i,j,k) = DXX(i,j,k) + DYY(i,j,k) + DZZ(i,j,k)
                        endif
                        epsilondev_xx_loc(i,j,k) = DXX(i,j,k) - epsilon_trace_over_3
                        epsilondev_yy_loc(i,j,k) = DYY(i,j,k) - epsilon_trace_over_3
                        epsilondev_xy_loc(i,j,k) = 0.5 * (DXY(i,j,k) + DYX(i,j,k))
                        epsilondev_xz_loc(i,j,k) = 0.5 * (DZX(i,j,k) + DXZ(i,j,k))
                        epsilondev_yz_loc(i,j,k) = 0.5 * (DZY(i,j,k) + DYZ(i,j,k))
                    enddo
                enddo
            enddo
        endif

        if (aniso) then
            if (n_solid>0) then
                call calcul_forces_aniso_att(Fox,Foy,Foz, &
                    Elem%Invgrad, &
                    htprimex, htprimey, htprimez, &
                    Elem%Jacob, mat%GLLwx, mat%GLLwy, mat%GLLwz, &
                    DXX,DXY,DXZ, &
                    DYX,DYY,DYZ, &
                    DZX,DZY,DZZ, &
                    Elem%Mu, Elem%Lambda, Elem%sl%Cij, &
                    
                    m1,m2,m3, n_solid, &
                    Elem%sl%onemSbeta, Elem%sl%R_xx_, Elem%sl%R_yy_, &
                    Elem%sl%R_xy_, Elem%sl%R_xz_, Elem%sl%R_yz_)
                call attenuation_aniso_update(Elem%sl%epsilondev_xx_,Elem%sl%epsilondev_yy_, &
                    Elem%sl%epsilondev_xy_,Elem%sl%epsilondev_xz_,Elem%sl%epsilondev_yz_, &
                    epsilondev_xx_loc,epsilondev_yy_loc, &
                    epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                    Elem%sl%R_xx_,Elem%sl%R_yy_,Elem%sl%R_xy_,Elem%sl%R_xz_,Elem%sl%R_yz_, &
                    Elem%sl%factor_common_3,Elem%sl%alphaval_3,Elem%sl%betaval_3,Elem%sl%gammaval_3, &
                    Elem%Mu, m1,m2,m3, n_solid)
                deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
                !      deallocate(epsilonvol_loc)
            else
                call calcul_forces_aniso(Fox,Foy,Foz,  &
                    Elem%Invgrad, &
                    htprimex, htprimey, htprimez, &
                    Elem%Jacob, mat%GLLwx, mat%GLLwy, mat%GLLwz, &
                    DXX,DXY,DXZ, &
                    DYX,DYY,DYZ, &
                    DZX,DZY,DZZ, &
                    Elem%sl%Cij, &
                    m1,m2,m3)
            endif
        else
            if (n_solid>0) then
                call calcul_forces_att(Fox,Foy,Foz, &
                    Elem%Invgrad, &
                    htprimex, htprimey, htprimez, &
                    Elem%Jacob, mat%GLLwx, mat%GLLwy, mat%GLLwz, &
                    DXX,DXY,DXZ, &
                    DYX,DYY,DYZ, &
                    DZX,DZY,DZZ, &
                    Elem%Mu, Elem%Kappa, &
                    m1,m2,m3, n_solid, &
                    Elem%sl%R_xx_, Elem%sl%R_yy_, &
                    Elem%sl%R_xy_, Elem%sl%R_xz_, Elem%sl%R_yz_, &
                    Elem%sl%R_vol_, &
                    Elem%sl%onemSbeta, &
                    Elem%sl%onemPbeta &
                    )
                !                             Elem%Mu, Elem%Lambda, &
                !                             m1,m2,m3, n_solid, &
                !                             Elem%onemSbeta, Elem%R_xx_, Elem%R_yy_, &
                !                             Elem%R_xy_, Elem%R_xz_, Elem%R_yz_)
                call attenuation_update(Elem%sl%epsilondev_xx_,Elem%sl%epsilondev_yy_, &
                    Elem%sl%epsilondev_xy_,Elem%sl%epsilondev_xz_,Elem%sl%epsilondev_yz_, &
                    epsilondev_xx_loc,epsilondev_yy_loc, &
                    epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc, &
                    Elem%sl%R_xx_,Elem%sl%R_yy_,Elem%sl%R_xy_,Elem%sl%R_xz_,Elem%sl%R_yz_, &
                    Elem%sl%factor_common_3,Elem%sl%alphaval_3,Elem%sl%betaval_3,Elem%sl%gammaval_3, &
                    Elem%Mu, m1,m2,m3, n_solid, &
                    Elem%Kappa, Elem%sl%epsilonvol_,epsilonvol_loc,Elem%sl%R_vol_, &
                    Elem%sl%factor_common_P,Elem%sl%alphaval_P,Elem%sl%betaval_P,Elem%sl%gammaval_P)
                !                              kappa_,elsilonvol_,epsilon_vol_loc,R_vol_, &
                !                              factor_common_P_,alphaval_P_,betaval_P_,gammaval_P_)
                deallocate(epsilondev_xx_loc,epsilondev_yy_loc,epsilondev_xy_loc,epsilondev_xz_loc,epsilondev_yz_loc)
                deallocate(epsilonvol_loc)
            else
                call calcul_forces(Fox,Foy,Foz,  &
                    Elem%Invgrad, &
                    htprimex, htprimey, htprimez, &
                    Elem%Jacob, mat%GLLwx, mat%GLLwy, mat%GLLwz, &
                    DXX,DXY,DXZ, &
                    DYX,DYY,DYZ, &
                    DZX,DZY,DZZ, &
                    Elem%Mu, Elem%Lambda, &
                    m1,m2,m3)
            endif
        endif

        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    champs1%Forces(Elem%Idom(i,j,k),0) = champs1%Forces(Elem%Idom(i,j,k),0)-Fox(i,j,k)
                    champs1%Forces(Elem%Idom(i,j,k),1) = champs1%Forces(Elem%Idom(i,j,k),1)-Foy(i,j,k)
                    champs1%Forces(Elem%Idom(i,j,k),2) = champs1%Forces(Elem%Idom(i,j,k),2)-Foz(i,j,k)
                enddo
            enddo
        enddo

        return
    end subroutine forces_int_solid

end module dom_solid

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
