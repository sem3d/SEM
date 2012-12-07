!>
!! \file calcul_forces_att.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

!! initial
!!subroutine calcul_forces_att(Fox,Foy,Foz, xi1,xi2,xi3,et1,et2,et3,ga1,ga2,ga3, &
!!     dx,dy,dz, jac, poidsx,poidsy,poidsz, DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
!!     mu_,ka_, ngllx,nglly,ngllz, n_solid, onemSbeta_, R_xx_,R_yy_,R_xy_,R_xz_,R_yz_, &
!!     onemPbeta_, R_vol_)
subroutine calcul_forces_att(Fox,Foy,Foz, xi1,xi2,xi3,et1,et2,et3,ga1,ga2,ga3, &
    dx,dy,dz, jac, poidsx,poidsy,poidsz, DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ, &
    mu_,ka_, ngllx,nglly,ngllz, n_solid, R_xx_,R_yy_,R_xy_,R_xz_,R_yz_,R_vol_)
    !                             mu_,la_, ngllx,nglly,ngllz, n_solid, onemSbeta_, R_xx_,R_yy_,R_xy_,R_xz_,R_yz_)

    !                             Elem%Mu, Elem%Kappa, &
    !                             m1,m2,m3, n_solid, &
    !                             Elem%onemSbeta, Elem%R_xx_, Elem%R_yy_, &
    !                             Elem%R_xy_, Elem%R_xz_, Elem%R_yz_, &
    !                             Elem%onemPbeta, Elem%R_vol_)

    !                             Elem%Mu, Elem%Lambda, &
    !                             m1,m2,m3, n_solid, &
    !                             Elem%onemSbeta, Elem%R_xx_, Elem%R_yy_, &
    !                             Elem%R_xy_, Elem%R_xz_, Elem%R_yz_)
    use sdomain

    implicit none

    integer, intent(in) :: ngllx,nglly,ngllz, n_solid
    real, dimension(0:n_solid-1,0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: R_vol_,R_xx_,R_yy_,R_xy_,R_xz_,R_yz_
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: xi1,xi2,xi3, et1,et2,et3, &
        ga1,ga2,ga3, jac, mu_,ka_, &
        DXX,DXY,DXZ,DYX,DYY,DYZ,DZX,DZY,DZZ
    !! real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: onemSbeta, onemPbeta_  !!initial
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: Fox,Foz,Foy
    real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: dx
    real, dimension(0:nglly-1,0:nglly-1), intent(in) :: dy
    real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: dz
    real, dimension(0:ngllx-1), intent(in) :: poidsx
    real, dimension(0:nglly-1), intent(in) :: poidsy
    real, dimension(0:ngllz-1), intent(in) :: poidsz

    real :: sxx,sxy,sxz,syy,syz,szz,t4,F1,F2,F3
    real :: xpression , stt
    real :: t41,t42,t43,t11,t51,t52,t53,t12,t61,t62,t63,t13
    real :: xt1,xt2,xt3,xt5,xt6,xt7,xt8,xt9,xt10
    real, parameter :: deuxtiers = 0.666666666666667, &
        quatretiers = 1.333333333333337, &
        s2 = 1.414213562373095, &
        s2o2 = 0.707106781186547, &
        zero = 0., two = 2.
    integer :: i,j,k,l, i_sls
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: xmu,x2mu,xkappa
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t1,t5,t8,t2,t6,t9
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: t3,t7,t10


    xmu(:,:,:) = mu_(:,:,:)
    !  modif mariotti fevrier 2007 cea
    !xmu(:,:,:) = xmu(:,:,:) * onemSbeta_(:,:,:)
    xkappa(:,:,:) = ka_(:,:,:)
    !  modif mariotti fevrier 2007 cea
    !xkappa(:,:,:) = xkappa(:,:,:) * onemPbeta_(:,:,:)
    x2mu(:,:,:) = two * xmu(:,:,:)

    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1


                !if ( i== 0 .and. j == 0 .and. k == 0 ) then
                !print*,' xkappa ', xkappa(i,j,k)
                !print*,' xmu ', xmu(i,j,k)
                !print*,' x2mu ', x2mu(i,j,k)
                !print*,' onemSbeta ', onemSbeta_(i,j,k)
                !print*,' xlambda ', xkappa(i,j,k)-(x2mu(i,j,k)/3.)
                !print*,' R_vol ', R_vol_(0,i,j,k)
                !if ( R_xx_(0,i,j,k)*R_xx_(0,i,j,k) > 0. ) then
                !print*,' R_xx ', R_xx_(0,i,j,k)
                !endif
                !print*,' R_yy ', R_yy_(0,i,j,k)
                !print*,' R_xy ', R_xy_(0,i,j,k)
                !print*,' R_xz ', R_xz_(0,i,j,k)
                !print*,' R_yz ', R_yz_(0,i,j,k)
                !endif
                xpression = xkappa(i,j,k) * ( DXX(i,j,k) + &
                    DYY(i,j,k) + DZZ(i,j,k) )

                sxx = x2mu(i,j,k) * DXX(i,j,k)

                sxy = xmu(i,j,k)  *( DXY(i,j,k)  +       &
                    DYX(i,j,k)  )

                sxz = xmu(i,j,k) * ( DXZ(i,j,k)  +       &
                    DZX(i,j,k)  )

                syy = x2mu(i,j,k) * DYY(i,j,k)

                syz = xmu(i,j,k) * ( DYZ(i,j,k)  +       &
                    DZY(i,j,k)  )

                szz = x2mu(i,j,k) * DZZ(i,j,k)

                !        stt = (sxx + syy + szz)/3.
                do i_sls = 0,n_solid-1
                    xpression = xpression - R_vol_(i_sls,i,j,k)
                    sxx = sxx  - R_xx_(i_sls,i,j,k)
                    syy = syy  - R_yy_(i_sls,i,j,k)
                    !         sxx = sxx - stt - R_xx_(i_sls,i,j,k)
                    !         syy = syy - stt - R_yy_(i_sls,i,j,k)
                    ! ici on utilise le fait que la trace est nulle
                    szz = szz  + R_xx_(i_sls,i,j,k) + R_yy_(i_sls,i,j,k)
                    !         szz = szz - stt + R_xx_(i_sls,i,j,k) + R_yy_(i_sls,i,j,k)
                    sxy = sxy - R_xy_(i_sls,i,j,k)
                    sxz = sxz - R_xz_(i_sls,i,j,k)
                    syz = syz - R_yz_(i_sls,i,j,k)
                enddo
                stt = (sxx + syy + szz)/3.
                sxx = sxx - stt + xpression
                syy = syy - stt + xpression
                szz = szz - stt + xpression
                !if ( i== 0 .and. j == 0 .and. k == 0 ) then
                !print*,' stt xpression ',3.*stt/x2mu(i,j,k),xpression/xkappa(i,j,k)
                !endif
                !
                !=====================
                !       FX
                !=====================
                !
                xt1 = sxx * xi1(i,j,k)
                xt2 = sxx * et1(i,j,k)
                xt3 = sxx * ga1(i,j,k)

                xt1 = xt1 + sxy * xi2(i,j,k)
                xt2 = xt2 + sxy * et2(i,j,k)
                xt3 = xt3 + sxy * ga2(i,j,k)

                xt1 = xt1 + sxz * xi3(i,j,k)
                xt2 = xt2 + sxz * et3(i,j,k)
                xt3 = xt3 + sxz * ga3(i,j,k)
                !
                !=====================
                !       FY
                !=====================
                !
                xt5 = syy * xi2(i,j,k)
                xt6 = syy * et2(i,j,k)
                xt7 = syy * ga2(i,j,k)

                xt5 = xt5 + sxy * xi1(i,j,k)
                xt6 = xt6 + sxy * et1(i,j,k)
                xt7 = xt7 + sxy * ga1(i,j,k)

                xt5 = xt5 + syz * xi3(i,j,k)
                xt6 = xt6 + syz * et3(i,j,k)
                xt7 = xt7 + syz * ga3(i,j,k)
                !
                !=====================
                !       FZ
                !=====================
                !
                xt8  = szz * xi3(i,j,k)
                xt9  = szz * et3(i,j,k)
                xt10 = szz * ga3(i,j,k)

                xt8  = xt8  + sxz * xi1(i,j,k)
                xt9  = xt9  + sxz * et1(i,j,k)
                xt10 = xt10 + sxz * ga1(i,j,k)

                xt8  = xt8  + syz * xi2(i,j,k)
                xt9  = xt9  + syz * et2(i,j,k)
                xt10 = xt10 + syz * ga2(i,j,k)


                !
                !- Multiplication par le Jacobien et le poids d'integration
                !
                t4 = jac(i,j,k) * poidsx(i)
                xt1  =  xt1 * t4
                xt5  =  xt5 * t4
                xt8  =  xt8 * t4

                t4 = jac(i,j,k) * poidsy(j)
                xt2  =  xt2 * t4
                xt6  =  xt6 * t4
                xt9  =  xt9 * t4

                t4 = jac(i,j,k) * poidsz(k)
                xt3  =  xt3 * t4
                xt7  =  xt7 * t4
                xt10 = xt10 * t4


                t1(i,j,k) = xt1
                t5(i,j,k) = xt5
                t8(i,j,k) = xt8

                t2(j,i,k) = xt2
                t6(j,i,k) = xt6
                t9(j,i,k) = xt9

                t3(k,i,j) = xt3
                t7(k,i,j) = xt7
                t10(k,i,j) = xt10


            enddo
        enddo
    enddo

    !
    !- Multiplication par la matrice de derivation puis par les poids
    !

    !=-=-=-=-=-=-=-=-=-=-
    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1
                !=-=-=-=-=-=-=-=-=-=-
                !
                t11 = poidsy(j) * poidsz(k)
                t12 = poidsx(i) * poidsz(k)
                t13 = poidsx(i) * poidsy(j)
                !
                t41 = zero
                t42 = zero
                t43 = zero
                t51 = zero
                t52 = zero
                t53 = zero
                t61 = zero
                t62 = zero
                t63 = zero
                !
                do l = 0,ngllx-1
                    t41 = t41 + dx(l,i) * t1(l,j,k)
                    t42 = t42 + dx(l,i) * t5(l,j,k)
                    t43 = t43 + dx(l,i) * t8(l,j,k)
                enddo

                do l = 0,nglly-1
                    t51 = t51 + dy(l,j) * t2(l,i,k)
                    t52 = t52 + dy(l,j) * t6(l,i,k)
                    t53 = t53 + dy(l,j) * t9(l,i,k)
                enddo
                ! FX
                F1 = t41*t11 + t51*t12
                ! FY
                F2 = t42*t11 + t52*t12
                ! FZ
                F3 = t43*t11 + t53*t12
                !
                !
                do l = 0,ngllz-1
                    t61 = t61 + dz(l,k) * t3(l,i,j)
                    t62 = t62 + dz(l,k) * t7(l,i,j)
                    t63 = t63 + dz(l,k) * t10(l,i,j)
                enddo

                ! FX
                F1 = F1 + t61*t13
                ! FY
                F2 = F2 + t62*t13
                ! FZ
                F3 = F3 + t63*t13
                !
                Fox(i,j,k) = F1
                Foy(i,j,k) = F2
                Foz(i,j,k) = F3

                !=-=-=-=-=-=-=-=-=-=-
            enddo
        enddo
    enddo
    !=-=-=-=-=-=-=-=-=-=-

end subroutine calcul_forces_att
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
