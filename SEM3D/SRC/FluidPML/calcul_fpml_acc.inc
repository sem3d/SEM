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
    call PROC_PRED(ngll, dom%nblocks, dom%nglltot, dom%dt, &
        dom%GLLw, dom%hprime, dom%m_InvGrad, champs%fpml_VelPhi, dom%m_Idom,  &
        dom%m_PMLVeloc, dom%m_PMLDumpSx, dom%m_PMLDumpSy, dom%m_PMLDumpSz)
    call PROC_FPML(ngll, dom%nblocks, dom%nglltot, &
        dom%GLLw, dom%hprime, dom%m_Jacob, dom%m_InvGrad, &
        dom%m_PMLVeloc, champs%fpml_Forces, dom%m_Idom)

end subroutine PROCNAME

subroutine PROC_FPML(ngllvar, nblocks, nglltot, &
    gllw, hprime, Jacob, InvGrad, PmlVeloc, fpml_Forces, Idom)
    integer, intent(in) :: ngllvar, nblocks, nglltot
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
    integer :: i, j, k, l, d, ind,e,ee, bnum
    real(fpp), dimension(0:2,0:2) :: acoeff
    real(fpp) :: acoeff_ijk

    !$acc parallel loop &
    !$acc& async(1) &
    !$acc& collapse(5) &
    !$acc& copyin(invgrad,idom,jacob,hprime,pmlveloc,gllw) &
    !$acc& present(fpml_Forces) &
    !$acc& firstprivate(nblocks) &
    !$acc& private(d,l,ind,acoeff,acoeff_ijk)
    do bnum = 0,nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    do l=0,2
                        do d=0,2
                            acoeff_ijk = -GLLw(i)*GLLw(j)*GLLw(k)*Jacob(ee,i,j,k,bnum)
                            acoeff(d,l) =  acoeff_ijk * PMLVeloc(ee,d,i,j,k,bnum)*InvGrad(ee,d,l,i,j,k,bnum)
                        end do
                    end do
                    do l = 0,ngll-1
                        ind = Idom(ee,l,j,k,bnum)
                        do d=0,2
                            OMPACC_ATOMIC_OP(update)
                            fpml_Forces(ind,d) = fpml_Forces(ind,d) + hprime(l,i)*acoeff(d,0)
                        end do
                        ind = Idom(ee,i,l,k,bnum)
                        do d=0,2
                            OMPACC_ATOMIC_OP(update)
                            fpml_Forces(ind,d) = fpml_Forces(ind,d) + hprime(l,j)*acoeff(d,1)
                        end do
                        ind = Idom(ee,i,j,l,bnum)
                        do d=0,2
                            OMPACC_ATOMIC_OP(update)
                            fpml_Forces(ind,d) = fpml_Forces(ind,d) + hprime(l,k)*acoeff(d,2)
                        end do
                    end do
                    END_SUBELEM_LOOP()
                end do
            end do
        end do
    end do
    !$acc end parallel loop
end subroutine PROC_FPML


subroutine PROC_PRED(ngllvar, nblocks, nglltot, dt, &
    gllw, hprime, InvGrad, fpml_VelPhi, Idom, &
    pmlveloc, dumpsx, dumpsy, dumpsz &
    )
    real(fpp), intent(in) :: dt
    integer, intent(in) :: nblocks, nglltot
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    real(fpp), intent(in), dimension(0:ngll-1) :: gllw
    real(fpp), intent(in), dimension(0:ngll-1,0:ngll-1) :: hprime
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: invgrad
    real(fpp), intent(inout), dimension(0:VCHUNK-1,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: pmlveloc
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: dumpsx
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: dumpsy
    real(fpp), intent(in), dimension(0:VCHUNK-1,0:1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: dumpsz
    real(fpp), intent(in), dimension(0:nglltot,0:2) :: fpml_VelPhi
    integer, intent(in), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1) :: Idom
    !
    real(fpp) :: dVelPhi_dx, dVelPhi_dy, dVelPhi_dz, VelPhi
    integer :: i, j, k, l, ind, e, ee
    real(fpp) :: dVPhi_dxi,dVPhi_deta,dVPhi_dzeta
    real(fpp) :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3
    integer :: bnum


    ! prediction in the element
    ! We do the sum V1+V2+V3 and F1+F2+F3 here
    !$acc parallel loop &
    !$acc& async(1) &
    !$acc& collapse(5) &
    !$acc& copyin(Idom,fpml_VelPhi,hprime,InvGrad) &
    !$acc& present(PMLVeloc,dumpsx,dumpsy,dumpsz) &
    !$acc& private(dVPhi_dxi,dVPhi_deta,dVPhi_dzeta,VelPhi) &
    !$acc& private(dVelPhi_dx, dVelPhi_dy, dVelPhi_dz) &
    !$acc& private(ind,l,xi1,xi3,xi2,et1,et2,et3,ga1,ga2,ga3) &
    !$acc& firstprivate(nblocks,dt)
    do bnum = 0,nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    ! d(rho*Phi)_d(xi,eta,zeta)
                    dVPhi_dxi   = 0D0
                    DO L = 0, ngll-1
                        ind = Idom(ee,L,j,k,bnum)
                        VelPhi = fpml_VelPhi(ind,0) + fpml_VelPhi(ind,1) + fpml_VelPhi(ind,2)
                        dVPhi_dxi   = dVPhi_dxi  +VelPhi*hprime(L,I)
                    END DO

                    dVPhi_deta  = 0D0
                    DO L = 0, ngll-1
                        ind = Idom(ee,i,L,k,bnum)
                        VelPhi = fpml_VelPhi(ind,0) + fpml_VelPhi(ind,1) + fpml_VelPhi(ind,2)
                        dVPhi_deta  = dVPhi_deta +VelPhi*hprime(L,J)
                    END DO

                    dVPhi_dzeta = 0D0
                    DO L = 0, ngll-1
                        ind = Idom(ee,i,j,L,bnum)
                        VelPhi = fpml_VelPhi(ind,0) + fpml_VelPhi(ind,1) + fpml_VelPhi(ind,2)
                        dVPhi_dzeta = dVPhi_dzeta+VelPhi*hprime(L,K)
                    END DO

                    xi1 = InvGrad(ee,0,0,i,j,k,bnum)
                    xi2 = InvGrad(ee,1,0,i,j,k,bnum)
                    xi3 = InvGrad(ee,2,0,i,j,k,bnum)
                    et1 = InvGrad(ee,0,1,i,j,k,bnum)
                    et2 = InvGrad(ee,1,1,i,j,k,bnum)
                    et3 = InvGrad(ee,2,1,i,j,k,bnum)
                    ga1 = InvGrad(ee,0,2,i,j,k,bnum)
                    ga2 = InvGrad(ee,1,2,i,j,k,bnum)
                    ga3 = InvGrad(ee,2,2,i,j,k,bnum)
                    !- in the physical domain
                    dVelPhi_dx = dVPhi_dxi*xi1 + dVPhi_deta*et1 + dVPhi_dzeta*ga1
                    dVelPhi_dy = dVPhi_dxi*xi2 + dVPhi_deta*et2 + dVPhi_dzeta*ga2
                    dVelPhi_dz = dVPhi_dxi*xi3 + dVPhi_deta*et3 + dVPhi_dzeta*ga3

                    ! prediction for (physical) velocity (which is the equivalent of a stress, here)
                    ! V_x^x
                    PMLVeloc(ee,0,i,j,k,bnum) = DumpSx(ee,0,i,j,k,bnum) * &
                        PMLVeloc(ee,0,i,j,k,bnum) + DumpSx(ee,1,i,j,k,bnum) * Dt * dVelPhi_dx
                    ! V_x^y = 0
                    ! V_x^z = 0
                    ! V_y^x = 0
                    ! V_y^y
                    PMLVeloc(ee,1,i,j,k,bnum) = DumpSy(ee,0,i,j,k,bnum) * &
                        PMLVeloc(ee,1,i,j,k,bnum) + DumpSy(ee,1,i,j,k,bnum) * Dt * dVelPhi_dy
                    ! V_y^z = 0
                    ! V_z^x = 0
                    ! V_z^y = 0
                    ! V_z^z
                    PMLVeloc(ee,2,i,j,k,bnum) = DumpSz(ee,0,i,j,k,bnum) * &
                        PMLVeloc(ee,2,i,j,k,bnum) + DumpSz(ee,1,i,j,k,bnum) * Dt * dVelPhi_dz
                    END_SUBELEM_LOOP()
                enddo
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
