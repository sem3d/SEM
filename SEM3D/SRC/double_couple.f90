!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file double_couple.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine double_couple(Tdomain,rg)

    ! Modified by Paul Cupillard 21/03/2006

    use sdomain
    use constants, only : M_PI

    implicit none

    type (domain), intent(inout) :: Tdomain
    integer, intent(in) :: rg

    integer :: n, elem, x,y,z, ngll, mat, a,b, k
    real :: ct,st,cp,sp, xi,eta,zeta, coord, xa,ya,za, xb,yb,zb, num,denom,prod
    real, dimension(:), allocatable :: xpol,ypol,zpol, xdh,ydh,zdh
    real, dimension(0:2,0:2) :: Pcs, Psc, TMP, M, Rot, tRot


    do n = 0,Tdomain%n_source-1
        if (Tdomain%sSource(n)%i_type_source == 2 .and. rg == Tdomain%sSource(n)%proc) then

            M = Tdomain%sSource(n)%Moment

!!            if (Tdomain%curve) then
!!                Tdomain%sSource(n)%realcolat = M_PI*Tdomain%sSource(n)%realcolat/180
!!                Tdomain%sSource(n)%reallong = M_PI*Tdomain%sSource(n)%reallong/180
!!                ct=cos(Tdomain%sSource(n)%realcolat)  !!dcos(Tdomain%sSource(n)%realcolat)  !!Avt Gsa Ipsis
!!                st=sin(Tdomain%sSource(n)%realcolat) !!dsin(Tdomain%sSource(n)%realcolat)
!!                cp=cos(Tdomain%sSource(n)%reallong) !!dcos(Tdomain%sSource(n)%reallong)
!!                sp=sin(Tdomain%sSource(n)%reallong) !!dsin(Tdomain%sSource(n)%reallong)
!!                ! matrice de passage du systeme cartesien au systeme spherique
!!                Pcs(0,0) = st*cp; Pcs(0,1) = ct*cp; Pcs(0,2) = -sp
!!                Pcs(1,0) = st*sp; Pcs(1,1) = ct*sp; Pcs(1,2) = cp
!!                Pcs(2,0) = ct   ; Pcs(2,1) = -st  ; Pcs(2,2) = 0.0d0
!!                ! matrice de passage du systeme spherique au systeme cartesien (inverse de Pcs)
!!                Psc(0,0) = st*cp; Psc(0,1) = st*sp; Psc(0,2) = ct
!!                Psc(1,0) = ct*cp; Psc(1,1) = ct*sp; Psc(1,2) = -st
!!                Psc(2,0) = -sp  ; Psc(2,1) = cp   ; Psc(2,2) = 0.0d0
!!                ! Calcul du tenseur moment sismique dans le systeme cartesien
!!                M = Tdomain%sSource(n)%Moment
!!                do a = 0,2
!!                    do b = 0,2
!!                        TMP(a,b) = 0.0d0
!!                        do k = 0,2
!!                            TMP(a,b) = TMP(a,b) + M(a,k)*Psc(k,b)
!!                        enddo
!!                    enddo
!!                enddo
!!                do a = 0,2
!!                    do b = 0,2
!!                        M(a,b) = 0.0d0
!!                        do k = 0,2
!!                            M(a,b) = M(a,b) + Pcs(a,k)*TMP(k,b)
!!                        enddo
!!                    enddo
!!                enddo
!!                ! Rotation (du chunk reel au chunk de reference) du tenseur moment sismique
!!                Rot = Tdomain%rot
!!                tRot = transpose(Rot)
!!                do a = 0,2
!!                    do b = 0,2
!!                        TMP(a,b) = 0.0d0
!!                        do k = 0,2
!!                            TMP(a,b) = TMP(a,b) + M(a,k)*Rot(k,b)
!!                        enddo
!!                    enddo
!!                enddo
!!                do a = 0,2
!!                    do b = 0,2
!!                        M(a,b) = 0.0d0
!!                        do k = 0,2
!!                            M(a,b) = M(a,b) + tRot(a,k)*TMP(k,b)
!!                        enddo
!!                    enddo
!!                enddo
!!            endif

            do a = 0,2
                do b = 0,2
                    TMP(a,b) = 0.0d0
                    do k = 0,2
                        TMP(a,b) = TMP(a,b) + Tdomain%sSource(n)%InvGrad(k,a)*M(b,k)
                    enddo
                enddo
            enddo

            elem = Tdomain%Ssource(n)%elem
            ngll = 0
            select case (Tdomain%specel(elem)%domain)
                 case (DM_SOLID)
                     ngll = Tdomain%sdom%ngll
                 case (DM_FLUID)
                     ngll = Tdomain%fdom%ngll
                 case (DM_SOLID_PML)
                     ngll = Tdomain%spmldom%ngll
                 case (DM_FLUID_PML)
                     ngll = Tdomain%fpmldom%ngll
            end select
            allocate (xpol(0:ngll-1));   xpol = 1
            allocate (ypol(0:ngll-1));   ypol = 1
            allocate (zpol(0:ngll-1));   zpol = 1
            allocate (xdh(0:ngll-1))
            allocate (ydh(0:ngll-1))
            allocate (zdh(0:ngll-1))
            xi = Tdomain%sSource(n)%refcoord(0)
            eta = Tdomain%sSource(n)%refcoord(1)
            zeta = Tdomain%sSource(n)%refcoord(2)
            mat = Tdomain%specel(elem)%mat_index
            do x = 0,ngll-1
                coord = Tdomain%sSubdomain(mat)%GLLc(x)
                num = 0;   denom = 1
                do a = 0,ngll-1
                    if (a/=x) then
                        xa = Tdomain%sSubdomain(mat)%GLLc(a)
                        xpol(x) = xpol(x) * (xi-xa)/(coord-xa)
                        denom = denom * (coord-xa)
                        prod = 1
                        do b = 0,ngll-1
                            if ((b/=x) .and. (b/=a)) then
                                xb = Tdomain%sSubdomain(mat)%GLLc(b)
                                prod = prod * (xi-xb)
                            endif
                        enddo
                        num = num + prod
                    endif
                enddo
                xdh(x) = num/denom
            enddo
            do y = 0,ngll-1
                coord = Tdomain%sSubdomain(mat)%GLLc(y)
                num = 0;   denom = 1
                do a = 0,ngll-1
                    if (a/=y) then
                        ya = Tdomain%sSubdomain(mat)%GLLc(a)
                        ypol(y) = ypol(y) * (eta-ya)/(coord-ya)
                        denom = denom * (coord-ya)
                        prod = 1
                        do b = 0,ngll-1
                            if ((b/=y) .and. (b/=a)) then
                                yb = Tdomain%sSubdomain(mat)%GLLc(b)
                                prod = prod * (eta-yb)
                            endif
                        enddo
                        num = num + prod
                    endif
                enddo
                ydh(y) = num/denom
            enddo
            do z = 0,ngll-1
                coord = Tdomain%sSubdomain(mat)%GLLc(z)
                num = 0;   denom = 1
                do a = 0,ngll-1
                    if (a/=z) then
                        za = Tdomain%sSubdomain(mat)%GLLc(a)
                        zpol(z) = zpol(z) * (zeta-za)/(coord-za)
                        denom = denom * (coord-za)
                        prod = 1
                        do b = 0,ngll-1
                            if ((b/=z) .and. (b/=a)) then
                                zb = Tdomain%sSubdomain(mat)%GLLc(b)
                                prod = prod * (zeta-zb)
                            endif
                        enddo
                        num = num + prod
                    endif
                enddo
                zdh(z) = num/denom
            enddo

            allocate (Tdomain%sSource(n)%coeff(0:ngll-1, 0:ngll-1, 0:ngll-1, 0:2))
            do x = 0,ngll-1
                do y = 0,ngll-1
                    do z = 0,ngll-1
                        Tdomain%sSource(n)%coeff(x,y,z,:) = xdh(x)*ypol(y)*zpol(z)*TMP(0,:) + &
                            xpol(x)*ydh(y)*zpol(z)*TMP(1,:) + &
                            xpol(x)*ypol(y)*zdh(z)*TMP(2,:)
                    enddo
                enddo
            enddo

            deallocate (xpol,ypol,zpol, xdh,ydh,zdh)

        endif
    enddo

    return
end subroutine double_couple

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
