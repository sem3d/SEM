!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!



subroutine PROCNAME(dom, ngllvar, champs)
    type(domain_fluidpml), intent(inout) :: dom
    type(champsfluidpml), intent(inout) :: champs
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    !
    integer :: n
    do n = 0,dom%nblocks-1
        call PROC_PRED(dom, ngll, dom%dt, champs, n)
        call PROC_FPML(ngll, dom%nblocks, dom%nglltot, n, 0, &
            dom%GLLw, dom%hprime, dom%m_Jacob, dom%m_InvGrad, &
            dom%m_PMLVeloc, champs%fpml_Forces, dom%m_Idom)
        call PROC_FPML(ngll, dom%nblocks, dom%nglltot, n, 1, &
            dom%GLLw, dom%hprime, dom%m_Jacob, dom%m_InvGrad, &
            dom%m_PMLVeloc, champs%fpml_Forces, dom%m_Idom)
        call PROC_FPML(ngll, dom%nblocks, dom%nglltot, n, 2, &
            dom%GLLw, dom%hprime, dom%m_Jacob, dom%m_InvGrad, &
            dom%m_PMLVeloc, champs%fpml_Forces, dom%m_Idom)
    end do

end subroutine PROCNAME

subroutine PROC_FPML(ngllvar, nblocks, nglltot, bnum, d, &
    gllw, hprime, Jacob, InvGrad, PmlVeloc, fpml_Forces, Idom)
    integer, intent(in) :: ngllvar, nblocks, nglltot, bnum, d
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    real(fpp), intent(in), dimension(0:ngll-1) :: gllw
    real(fpp), intent(in), dimension(0:ngll-1,0:ngll-1) :: hprime
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: jacob
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: invgrad
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: pmlveloc
    real(fpp), intent(inout), dimension(0:nglltot,0:2) :: fpml_Forces
    integer, intent(in), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: Idom
    !
    integer :: i, j, k, l, ind,e,ee
    real(fpp) :: acoeff
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: FFl

    FFl(:,:,:,:) = 0d0
    ! F(dim=d)
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i=0,ngll-1
                LOOP_VECTORIZE
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                acoeff = -GLLw(i)*GLLw(j)*GLLw(k)*Jacob(ee,i,j,k,bnum)*PmlVeloc(ee,d,i,j,k,bnum)
                do l = 0,ngll-1
                    FFl(ee,l,j,k) = FFl(ee,l,j,k) + hprime(l,i)*acoeff*InvGrad(ee,d,0,i,j,k,bnum)
                    FFl(ee,i,l,k) = FFl(ee,i,l,k) + hprime(l,j)*acoeff*InvGrad(ee,d,1,i,j,k,bnum)
                    FFl(ee,i,j,l) = FFl(ee,i,j,l) + hprime(l,k)*acoeff*InvGrad(ee,d,2,i,j,k,bnum)
                end do
                END_SUBELEM_LOOP()
            end do
        end do
    end do

    ! Assemblage
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ind = Idom(ee,i,j,k,bnum)
                ! We should have atomic adds with openmp here
                fpml_Forces(ind,d) = fpml_Forces(ind,d) + FFl(ee,i,j,k)
                END_SUBELEM_LOOP()
            enddo
        enddo
    enddo
end subroutine PROC_FPML


subroutine PROC_PRED(dom, ngllvar, dt, champs1, bnum)
    type (domain_fluidpml), intent (INOUT) :: dom
    real(fpp), intent(in) :: dt
    type(champsfluidpml), intent(inout) :: champs1
    integer, intent(in) :: bnum
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    !
    real(fpp) :: dVelPhi_dx, dVelPhi_dy, dVelPhi_dz
    integer :: i, j, k, l, ind, e, ee
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1) :: VelPhi
    real(fpp) :: dVPhi_dxi,dVPhi_deta,dVPhi_dzeta
    real(fpp) :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3


    ! prediction in the element
    ! We do the sum V1+V2+V3 and F1+F2+F3 here
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
                !$omp simd linear(ee) safelen(VCHUNK) private(ind)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ind = dom%Idom_(i,j,k,bnum,ee)
                VelPhi(ee,i,j,k) = champs1%fpml_VelPhi(ind,0) + &
                    champs1%fpml_VelPhi(ind,1) + &
                    champs1%fpml_VelPhi(ind,2)
                END_SUBELEM_LOOP()
            enddo
        enddo
    enddo

    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
                !$omp simd linear(ee) safelen(VCHUNK) private(dVPhi_dxi,dVPhi_deta,dVPhi_dzeta,dVelPhi_dx, dVelPhi_dy, dVelPhi_dz)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
                ! d(rho*Phi)_d(xi,eta,zeta)
                dVPhi_dxi   = 0D0
                dVPhi_deta  = 0D0
                dVPhi_dzeta = 0D0

                DO L = 0, ngll-1
                    dVPhi_dxi   = dVPhi_dxi  +VelPhi(ee,L,J,K)*dom%hprime(L,I)
                    dVPhi_deta  = dVPhi_deta +VelPhi(ee,I,L,K)*dom%hprime(L,J)
                    dVPhi_dzeta = dVPhi_dzeta+VelPhi(ee,I,J,L)*dom%hprime(L,K)
                END DO

                xi1 = dom%InvGrad_(0,0,i,j,k,bnum,ee)
                xi2 = dom%InvGrad_(1,0,i,j,k,bnum,ee)
                xi3 = dom%InvGrad_(2,0,i,j,k,bnum,ee)
                et1 = dom%InvGrad_(0,1,i,j,k,bnum,ee)
                et2 = dom%InvGrad_(1,1,i,j,k,bnum,ee)
                et3 = dom%InvGrad_(2,1,i,j,k,bnum,ee)
                ga1 = dom%InvGrad_(0,2,i,j,k,bnum,ee)
                ga2 = dom%InvGrad_(1,2,i,j,k,bnum,ee)
                ga3 = dom%InvGrad_(2,2,i,j,k,bnum,ee)
                !- in the physical domain
                dVelPhi_dx = dVPhi_dxi*xi1 + dVPhi_deta*et1 + dVPhi_dzeta*ga1
                dVelPhi_dy = dVPhi_dxi*xi2 + dVPhi_deta*et2 + dVPhi_dzeta*ga2
                dVelPhi_dz = dVPhi_dxi*xi3 + dVPhi_deta*et3 + dVPhi_dzeta*ga3

                ! prediction for (physical) velocity (which is the equivalent of a stress, here)
                ! V_x^x
                dom%PMLVeloc_(i,j,k,0,bnum,ee) = dom%PMLDumpSx_(i,j,k,0,bnum,ee) * &
                    dom%PMLVeloc_(i,j,k,0,bnum,ee) + &
                    dom%PMLDumpSx_(i,j,k,1,bnum,ee) * Dt * dVelPhi_dx
                ! V_x^y = 0
                ! V_x^z = 0
                ! V_y^x = 0
                ! V_y^y
                dom%PMLVeloc_(i,j,k,1,bnum,ee) = dom%PMLDumpSy_(i,j,k,0,bnum,ee) * &
                    dom%PMLVeloc_(i,j,k,1,bnum,ee) + &
                    dom%PMLDumpSy_(i,j,k,1,bnum,ee) * Dt * dVelPhi_dy
                ! V_y^z = 0
                ! V_z^x = 0
                ! V_z^y = 0
                ! V_z^z
                dom%PMLVeloc_(i,j,k,2,bnum,ee) = dom%PMLDumpSz_(i,j,k,0,bnum,ee) * &
                    dom%PMLVeloc_(i,j,k,2,bnum,ee) + &
                    dom%PMLDumpSz_(i,j,k,1,bnum,ee) * Dt * dVelPhi_dz
                END_SUBELEM_LOOP()
            enddo
        enddo
    enddo
end subroutine PROC_PRED


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
