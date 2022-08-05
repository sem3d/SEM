!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

#ifndef PROCNAME
#error "Must define PROCNAME before inclusion"
#endif

#undef ngllvar
#ifndef NGLLVAL
#define ngllvar ngll
#else
#define ngllvar ngll0
#endif


subroutine PROCNAME(dom,ngllvar,bnum,Q,dQdt)

    use champs_solid_dg
    implicit none

    type(domain_solid_dg), intent (INOUT) :: dom
    integer, intent(in) :: bnum
    integer, intent(in) :: ngllvar

#ifdef NGLLVAL
    integer, parameter :: ngll = NGLLVAL
#endif

    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:8), intent(out) :: dQdt
    real(fpp), dimension(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:8), intent(in)  :: Q

    integer   :: i, j, k, l, e, ee
    real(fpp) :: xmu, xla, tr_eps
    real(fpp) :: JAC, ww, hp
    real(fpp) :: dE11dt, dE22dt, dE33dt
    real(fpp) :: dE12dt, dE13dt, dE23dt
    real(fpp) :: dVxdt, dVydt, dVzdt
    real(fpp) :: E11, E22, E33,a
    real(fpp) :: E12, E13, E23
    real(fpp) :: S11, S22, S33
    real(fpp) :: S12, S13, S23
    real(fpp) :: Vx, Vy, Vz, tmp
    if (ngll/=ngllvar) stop 1

    dQdt = 0.d0
    
    do k = 0,ngll-1
        do j = 0,ngll-1
            do i = 0,ngll-1
#if VCHUNK>1
!$omp simd linear(e,ee)
#endif
                BEGIN_SUBELEM_LOOP(e,ee,bnum)
#define JXX_(i,j,k) dom%InvGrad_(0,0,i,j,k,bnum,ee)
#define JXY_(i,j,k) dom%InvGrad_(0,1,i,j,k,bnum,ee)
#define JXZ_(i,j,k) dom%InvGrad_(0,2,i,j,k,bnum,ee)
#define JYX_(i,j,k) dom%InvGrad_(1,0,i,j,k,bnum,ee)
#define JYY_(i,j,k) dom%InvGrad_(1,1,i,j,k,bnum,ee)
#define JYZ_(i,j,k) dom%InvGrad_(1,2,i,j,k,bnum,ee)
#define JZX_(i,j,k) dom%InvGrad_(2,0,i,j,k,bnum,ee)
#define JZY_(i,j,k) dom%InvGrad_(2,1,i,j,k,bnum,ee)
#define JZZ_(i,j,k) dom%InvGrad_(2,2,i,j,k,bnum,ee)


                JAC = dom%Jacob_(i,j,k,bnum,ee)
                xmu = dom%Mu_(i,j,k,bnum,ee)
                xla = dom%Lambda_(i,j,k,bnum,ee)

                ! TESTS
                !xla = 0.d0

                dE11dt = 0.d0
                dE22dt = 0.d0
                dE33dt = 0.d0
                dE12dt = 0.d0
                dE13dt = 0.d0
                dE23dt = 0.d0
                dVxdt  = 0.d0
                dVydt  = 0.d0
                dVzdt  = 0.d0
                ww = JAC*dom%GLLw(i)*dom%GLLw(j)*dom%GLLw(k)

                tmp = 0
                do L = 0, ngll-1
#define JXX JXX_(L,j,k)
#define JXY JXY_(L,j,k)
#define JXZ JXZ_(L,j,k)
#define JYX JYX_(L,j,k)
#define JYY JYY_(L,j,k)
#define JYZ JYZ_(L,j,k)
#define JZX JZX_(L,j,k)
#define JZY JZY_(L,j,k)
#define JZZ JZZ_(L,j,k)

                    E11 = Q(ee,L,j,k,0)
                    E22 = Q(ee,L,j,k,1)
                    E33 = Q(ee,L,j,k,2)
                    E12 = Q(ee,L,j,k,3)
                    E13 = Q(ee,L,j,k,4)
                    E23 = Q(ee,L,j,k,5)
                    Vx  = Q(ee,L,j,k,6)
                    Vy  = Q(ee,L,j,k,7)
                    Vz  = Q(ee,L,j,k,8)
                   
                    tr_eps = E11 + E22 + E33
                    S11 = xla*tr_eps + 2d0*xmu*E11
                    S22 = xla*tr_eps + 2d0*xmu*E22
                    S33 = xla*tr_eps + 2d0*xmu*E33
                    S12 = 2d0*xmu*E12
                    S13 = 2d0*xmu*E13
                    S23 = 2d0*xmu*E23
                    hp = dom%hprime(L,i)

                    dE11dt = dE11dt + Vx*JXX*hp
                    dE22dt = dE22dt + Vy*JXY*hp
                    dE33dt = dE33dt + Vz*JXZ*hp
                    
                    dE12dt = dE12dt + 0.5d0*(Vy*JXX + Vx*JXY)*hp
                    dE13dt = dE13dt + 0.5d0*(Vz*JXX + Vx*JXZ)*hp
                    dE23dt = dE23dt + 0.5d0*(Vz*JXY + Vy*JXZ)*hp

                    dVxdt = dVxdt + ( S11*JXX + S12*JXY + S13*JXZ )*hp
                    dVydt = dVydt + ( S12*JXX + S22*JXY + S23*JXZ )*hp
                    dVzdt = dVzdt + ( S13*JXX + S23*JXY + S33*JXZ )*hp                 

                end do

#undef JXX
#undef JXY
#undef JXZ
#undef JYX
#undef JYY
#undef JYZ
#undef JZX
#undef JZY
#undef JZZ

                do L = 0, ngll-1
#define JXX JXX_(i,L,k)
#define JXY JXY_(i,L,k)
#define JXZ JXZ_(i,L,k)
#define JYX JYX_(i,L,k)
#define JYY JYY_(i,L,k)
#define JYZ JYZ_(i,L,k)
#define JZX JZX_(i,L,k)
#define JZY JZY_(i,L,k)
#define JZZ JZZ_(i,L,k)

                    E11 = Q(ee,i,L,k,0)
                    E22 = Q(ee,i,L,k,1)
                    E33 = Q(ee,i,L,k,2)
                    E12 = Q(ee,i,L,k,3)
                    E13 = Q(ee,i,L,k,4)
                    E23 = Q(ee,i,L,k,5)
                    Vx  = Q(ee,i,L,k,6)
                    Vy  = Q(ee,i,L,k,7)
                    Vz  = Q(ee,i,L,k,8)
                   
                    tr_eps = E11 + E22 + E33
                    S11 = xla*tr_eps + 2d0*xmu*E11
                    S22 = xla*tr_eps + 2d0*xmu*E22
                    S33 = xla*tr_eps + 2d0*xmu*E33
                    S12 = 2d0*xmu*E12
                    S13 = 2d0*xmu*E13
                    S23 = 2d0*xmu*E23
                    hp = dom%hprime(L,j)
                    
                    dE11dt = dE11dt + Vx*JYX*hp 
                    dE22dt = dE22dt + Vy*JYY*hp
                    dE33dt = dE33dt + Vz*JYZ*hp

                    dE12dt = dE12dt + 0.5d0*(Vy*JYX + Vx*JYY)*hp
                    dE13dt = dE13dt + 0.5d0*(Vz*JYX + Vx*JYZ)*hp
                    dE23dt = dE23dt + 0.5d0*(Vz*JYY + Vy*JYZ)*hp
                    

                    dVxdt = dVxdt + ( S11*JYX + S12*JYY + S13*JYZ )*hp 
                    dVydt = dVydt + ( S12*JYX + S22*JYY + S23*JYZ )*hp
                    dVzdt = dVzdt + ( S13*JYX + S23*JYY + S33*JYZ )*hp

                    
                end do
#undef JXX
#undef JXY
#undef JXZ
#undef JYX
#undef JYY
#undef JYZ
#undef JZX
#undef JZY
#undef JZZ

                tmp = 0
                do L = 0, ngll-1
#define JXX JXX_(i,j,L)
#define JXY JXY_(i,j,L)
#define JXZ JXZ_(i,j,L)
#define JYX JYX_(i,j,L)
#define JYY JYY_(i,j,L)
#define JYZ JYZ_(i,j,L)
#define JZX JZX_(i,j,L)
#define JZY JZY_(i,j,L)
#define JZZ JZZ_(i,j,L)

                    E11 = Q(ee,i,j,L,0)
                    E22 = Q(ee,i,j,L,1)
                    E33 = Q(ee,i,j,L,2)
                    E12 = Q(ee,i,j,L,3)
                    E13 = Q(ee,i,j,L,4)
                    E23 = Q(ee,i,j,L,5)
                    Vx  = Q(ee,i,j,L,6)
                    Vy  = Q(ee,i,j,L,7)
                    Vz  = Q(ee,i,j,L,8)
                    
                    tr_eps = E11 + E22 + E33
                    S11 = xla*tr_eps + 2d0*xmu*E11
                    S22 = xla*tr_eps + 2d0*xmu*E22
                    S33 = xla*tr_eps + 2d0*xmu*E33
                    S12 = 2d0*xmu*E12
                    S13 = 2d0*xmu*E13
                    S23 = 2d0*xmu*E23
                    hp = dom%hprime(L,k)

                    dE11dt = dE11dt + Vx*JZX*hp
                    dE22dt = dE22dt + Vy*JZY*hp
                    dE33dt = dE33dt + Vz*JZZ*hp

                    dE12dt = dE12dt + 0.5d0*(Vy*JZX + Vx*JZY)*hp
                    dE13dt = dE13dt + 0.5d0*(Vz*JZX + Vx*JZZ)*hp
                    dE23dt = dE23dt + 0.5d0*(Vz*JZY + Vy*JZZ)*hp
                    
                    dVxdt = dVxdt + ( S11*JZX + S12*JZY + S13*JZZ )*hp
                    dVydt = dVydt + ( S12*JZX + S22*JZY + S23*JZZ )*hp
                    dVzdt = dVzdt + ( S13*JZX + S23*JZY + S33*JZZ )*hp

                    tmp = tmp + hp
                    
                end do
                
                !write(*,*) "i = ", i, " j = ", j, " k = ", k, "dE11dt = ", dE11dt
                !write(*,*) "i = ", i, " j = ", j, " k = ", k, "dE22dt = ", dE22dt
                !write(*,*) "i = ", i, " j = ", j, " k = ", k, "dE33dt = ", dE33dt
                !write(*,*) "i = ", i, " j = ", j, " k = ", k, "dE12dt = ", dE12dt
                !write(*,*) "i = ", i, " j = ", j, " k = ", k, "dE13dt = ", dE13dt
                !write(*,*) "i = ", i, " j = ", j, " k = ", k, "dE23dt = ", dE23dt
                !write(*,*) "i = ", i, " j = ", j, " k = ", k, "dVxdt  =  ", dVxdt
                !write(*,*) "i = ", i, " j = ", j, " k = ", k, "dVydt  = ", dVydt
                !write(*,*) "i = ", i, " j = ", j, " k = ", k, "dVzdt  = ", dVzdt
                
                dQdt(ee,i,j,k,0) = dE11dt*ww
                dQdt(ee,i,j,k,1) = dE22dt*ww
                dQdt(ee,i,j,k,2) = dE33dt*ww
                dQdt(ee,i,j,k,3) = dE12dt*ww
                dQdt(ee,i,j,k,4) = dE13dt*ww
                dQdt(ee,i,j,k,5) = dE23dt*ww
                dQdt(ee,i,j,k,6) = dVxdt*ww
                dQdt(ee,i,j,k,7) = dVydt*ww
                dQdt(ee,i,j,k,8) = dVzdt*ww
                
            END DO
        END DO
    END DO
    !=-=-=-=-=-=-=-=-=-=-

end subroutine PROCNAME

#undef JXX
#undef JXY
#undef JXZ
#undef JYX
#undef JYY
#undef JYZ
#undef JZX
#undef JZY
#undef JZZ

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
!! vim: set sw=4 ts=8 et tw=80 smartindent :e
