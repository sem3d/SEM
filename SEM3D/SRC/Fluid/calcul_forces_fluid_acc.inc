!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

subroutine PROCNAME(dom, ngllvar, var, dvdt)
    use champs_fluid
    type(domain_fluid), intent (INOUT) :: dom
    type(champsfluid), intent(in) :: var
    type(champsfluid), intent(inout) :: dvdt
    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    !
    integer :: i,j,k,ee,bnum,idx
    integer :: nblocks, nglltot
    real(fpp) :: uval
    !
    real(fpp) :: val

    nblocks = dom%nblocks
    nglltot = dom%nglltot

    ! initialize fields
    call PROC_INIT(ngllvar, nblocks, nglltot, dom%m_Idom, dom%Phi, dom%ForcesFl, var%Phi)
    ! internal forces
    call PROC_MAIN(ngllvar, nblocks, dom%ForcesFl, dom%Phi, dom%m_IDensity, &
        dom%hprime, dom%gllw, dom%m_Jacob, dom%m_InvGrad)
    ! assemble fields
    call PROC_ASSM(ngllvar, nblocks, nglltot, dom%m_Idom, dom%ForcesFl, dvdt%ForcesFl)
    !
end subroutine PROCNAME

subroutine PROC_INIT(ngllvar,nblocks,nglltot,idom,Phi,FFl,vPhi)
    use champs_fluid
    implicit none

    integer, intent(in) :: ngllvar, nglltot
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    integer, intent(in) :: nblocks
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(out) :: FFl
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(out) :: Phi
    real(fpp), dimension(0:nglltot), intent(in) :: vPhi
    integer, dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(in) :: idom
    !
    integer :: idx,i,j,k,ee,bnum
    !
    !$acc parallel loop &
    !$acc& collapse(5) &
    !$acc& async(1) &
    !$acc& copyin(idom,vPhi) &
    !$acc& present(phi,ffl) &
    !$acc& private(idx) &
    !$acc& firstprivate(nblocks)
    do bnum = 0, nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        idx = Idom(ee,i,j,k,bnum)
                        Phi(ee,i,j,k,bnum) = vPhi(idx)
                        FFl(ee,i,j,k,bnum) = 0d0
                    enddo
                enddo
            enddo
        enddo
    enddo
end subroutine PROC_INIT

subroutine PROC_ASSM(ngllvar,nblocks,nglltot,idom,FFl,vForces)
    use champs_fluid
    implicit none

    integer, intent(in) :: ngllvar, nglltot
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    integer, intent(in) :: nblocks
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(out) :: FFl
    real(fpp), dimension(0:nglltot), intent(inout) :: vForces
    integer, dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(in) :: idom
    !
    integer :: idx,i,j,k,ee,bnum
    !
    !$acc parallel loop &
    !$acc& collapse(5) &
    !$acc& async(1) &
    !$acc& copyin(idom,FFl) &
    !$acc& present(vForces) &
    !$acc& private(idx) &
    !$acc& firstprivate(nblocks)
    do bnum = 0, nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        idx = Idom(i,j,k,bnum,ee)
                        if ((i==0).or.(j==0).or.(k==0).or.(i==ngll-1).or.(j==ngll-1).or.(k==ngll-1)) then
                            !$acc atomic update
                            vForces(idx) = vForces(idx) - FFl(ee,i,j,k,bnum)
                        else
                            vForces(idx) = vForces(idx) - Ffl(ee,i,j,k,bnum)
                        end if
                    enddo
                enddo
            enddo
        enddo
    enddo
    !
end subroutine PROC_ASSM

subroutine PROC_MAIN(ngllvar,nblocks,FFl,Phi,irho,hprime,gllw,jacob,invgrad)
    use champs_fluid
    implicit none

    integer, intent(in) :: ngllvar
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif
    integer, intent(in) :: nblocks
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(out) :: FFl
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(in) :: Phi
    real(fpp), dimension(0:ngll-1,0:ngll-1), intent(in) :: hprime
    real(fpp), dimension(0:ngll-1), intent(in) :: gllw
    real(fpp), dimension(0:VCHUNK-1,0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(in) :: invgrad
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(in) :: jacob
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:nblocks-1), intent(in) :: irho
    !
    real(fpp) :: dPhi_dX,dPhi_dY,dPhi_dZ
    real(fpp) :: dPhi_dxi,dPhi_deta,dPhi_dzeta
    real(fpp) :: xi1,xi2,xi3, et1,et2,et3, ga1,ga2,ga3
    integer :: i,j,k,l,ee,bnum
    real(fpp) :: sx,sy,sz,t4
    real(fpp) :: xt1,xt6,xt10
    !real(fpp), parameter :: zero = 0d0
    real(fpp) :: xdens
#ifdef NGLLVAL
    if (ngll/=ngll0) stop 1
#endif
    !$acc parallel loop &
    !$acc& collapse(5) &
    !$acc& async(1) &
    !$acc& copyin(Phi,invgrad,jacob,irho,hprime,gllw) &
    !$acc& present(FFL) &
    !$acc& private(L,dPhi_dxi,dPhi_deta,dPhi_dzeta) &
    !$acc& private(dPhi_dX,dPhi_dY,dPhi_dZ) &
    !$acc& private(xi1,xi2,xi3,et1,et2,et3,ga1,ga2,ga3) &
    !$acc& firstprivate(nblocks)
    do bnum = 0, nblocks-1
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP1(ee)
                    ! Calcul de dPhi/dx dPhi/dy dPhi/dz
                    dPhi_dxi   = 0D0
                    dPhi_deta  = 0D0
                    dPhi_dzeta = 0D0
                    DO L = 0, ngll-1
                        dPhi_dxi   = dPhi_dxi  +Phi(ee,L,J,K,bnum)*hprime(L,I)
                        dPhi_deta  = dPhi_deta +Phi(ee,I,L,K,bnum)*hprime(L,J)
                        dPhi_dzeta = dPhi_dzeta+Phi(ee,I,J,L,bnum)*hprime(L,K)
                    END DO

                    xi1 = invgrad(ee,0,0,i,j,k,bnum)
                    xi2 = invgrad(ee,1,0,i,j,k,bnum)
                    xi3 = invgrad(ee,2,0,i,j,k,bnum)
                    et1 = invgrad(ee,0,1,i,j,k,bnum)
                    et2 = invgrad(ee,1,1,i,j,k,bnum)
                    et3 = invgrad(ee,2,1,i,j,k,bnum)
                    ga1 = invgrad(ee,0,2,i,j,k,bnum)
                    ga2 = invgrad(ee,1,2,i,j,k,bnum)
                    ga3 = invgrad(ee,2,2,i,j,k,bnum)
                    !- in the physical domain
                    dPhi_dx = dPhi_dxi*xi1 + dPhi_deta*et1 + dPhi_dzeta*ga1
                    dPhi_dy = dPhi_dxi*xi2 + dPhi_deta*et2 + dPhi_dzeta*ga2
                    dPhi_dz = dPhi_dxi*xi3 + dPhi_deta*et3 + dPhi_dzeta*ga3

                    ! (fluid equivalent) stress  ( = physical velocity)
                    xdens = irho(ee,i,j,k,bnum)
                    sx = xdens*dPhi_dX
                    sy = xdens*dPhi_dY
                    sz = xdens*dPhi_dZ


                    !=====================
                    !       F1 
                    xt1 = sx*xi1 + sy*xi2 + sz*xi3

                    !=====================
                    !       F2 
                    xt6 = sx*et1 + sy*et2 + sz*et3

                    !=====================
                    !       F3 
                    xt10 = sx*ga1 + sy*ga2 + sz*ga3

                    !
                    !- Multiply par Jacobian and weight
                    !
                    t4  = Jacob(ee,i,j,k,bnum) * gllw(i)
                    xt1 = xt1 * t4

                    t4  = Jacob(ee,i,j,k,bnum) * gllw(j)
                    xt6 = xt6 * t4

                    t4   = Jacob(ee,i,j,k,bnum) * gllw(k)
                    xt10 = xt10 * t4

                    DO L = 0, ngll-1
                        !$acc atomic update
                        FFl(ee,L,j,k,bnum) = FFl(ee,L,j,k,bnum) + hprime(L,I)*xt1*gllw(j) * gllw(k)
                        !$acc atomic update
                        FFl(ee,i,L,k,bnum) = FFl(ee,i,L,k,bnum) + hprime(L,J)*xt6*gllw(i) * gllw(k)
                        !$acc atomic update
                        FFl(ee,i,j,L,bnum) = FFl(ee,i,j,L,bnum) + hprime(L,K)*xt10*gllw(i) * gllw(j)
                    ENDDO
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    enddo
    !=-=-=-=-=-=-=-=-=-=-
end subroutine PROC_MAIN

!!!! ROUTINES DE TRAITEMENT DES MIRROIRS
subroutine PROC_MIRROR(dom, ngllvar, var, dvdt, m_dump, m_expl, m_recalc)
    use champs_fluid
    type(domain_fluid), intent (INOUT) :: dom
    type(champsfluid), intent(in) :: var
    type(champsfluid), intent(inout) :: dvdt
    integer, intent(in) :: ngllvar
    logical, intent(in) :: m_dump, m_expl, m_recalc
#ifdef NGLLVAL
    integer, parameter :: ngll=NGLLVAL
#endif

end subroutine PROC_MIRROR

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
