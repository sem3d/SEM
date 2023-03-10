!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file ondelette.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

module mondelette
    use constants, only : fpp
    implicit none
        !     ***************************************************************
        !     *                                                             *
        !     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
        !     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
        !     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
        !     *                                                             *
        !     ***************************************************************
    integer, parameter :: lvr=128
contains
    subroutine def_timefunc (Tdomain,rg)
        use sdomain
        use constants, only : M_PI
        use semdatafiles

        implicit none

        type (domain), intent(INOUT) :: Tdomain
        integer, intent(IN) :: rg

        character(Len=MAX_FILE_SIZE) :: fnamef
        integer :: j,nstep2,i,isrc,nstep,NBE_loc
        real(fpp) :: wt,freq,t1,t2,tmax,dt,t0,f1h,f2h,f3h,f4h
        real(fpp), parameter :: amp = 1
        complex(fpp) :: dphi
        complex(fpp), dimension(:), allocatable :: spectre,tmp


        nstep = Tdomain%TimeD%ntimeMax
        dt = Tdomain%TimeD%dtmin
        nstep2 = int(2.d0**(int(log(dble(nstep))/log(2.d0))+1))
        allocate(spectre(0:nstep2-1))
        allocate(tmp(1:nstep2))
        NBE_loc = Tdomain%n_source
        do isrc = 0,NBE_loc-1
            if (rg==Tdomain%sSource(isrc)%proc .and. Tdomain%sSource(isrc)%i_time_function==3) then
                allocate (Tdomain%sSource(isrc)%timefunc(0:nstep))
                Tdomain%sSource(isrc)%timefunc(0) = 0
                t0 = Tdomain%sSource(isrc)%tau_b
                f1h = Tdomain%sSource(isrc)%fh(0)
                f2h = Tdomain%sSource(isrc)%fh(1)
                f3h = Tdomain%sSource(isrc)%fh(2)
                f4h = Tdomain%sSource(isrc)%fh(3)
                spectre(:)=cmplx(0.d0,0.d0,fpp)
                do j = 1,nstep2
                    if (j<=nstep2/2) then
                        freq = (j-1)/(dt*nstep2)
                    else if (j==nstep2/2+1) then
                        freq = 1/(2.d0*dt)
                    else
                        freq = -(nstep2-j+1)/(dt*nstep2)
                    endif
                    dphi = exp(-2.d0*M_PI*freq*t0*cmplx(0.d0,1.d0,fpp))
                    call wtcoef(abs(freq),f1h,f2h,f3h,f4h,wt)
                    if (j/=0)   spectre(j-1) = wt*dphi
                enddo
                !le premier pas de temps correspond a t=0 avec la fft
                do j = 1,nstep2
                    tmp(j) = spectre(j-1)
                enddo
                call dfour1(tmp,nstep2,1)
                do j = 1,nstep2
                    spectre(j-1) = tmp(j)
                enddo
                Tdomain%sSource(isrc)%timefunc(1:nstep) = amp * real(spectre(0:nstep-1),fpp)/nstep2/dt
                !on met les premier pas de temps a zero:
                tmax = nstep2*dt
                t1 = 0.d0
                t2 = t0/5.d0
                do i = 0,nstep-1
                    call wtcoef(i*dt,t1,t2,tmax,tmax,wt)
                    Tdomain%sSource(isrc)%timefunc(i+1) = Tdomain%sSource(isrc)%timefunc(i+1) * wt
                enddo
                fnamef = pjoin(path_logs, "ondelette.txt")
                open (60,file=fnamef,status="unknown",form="formatted")
                do i = 0, Tdomain%TimeD%NtimeMax
                    t1 = i*dt
                    write (60,*) t1, Tdomain%sSource(isrc)%timefunc(i)
                enddo
                close (60)
            endif
        enddo

        deallocate(spectre,tmp)


    contains
        !-----------------------------------------------------
        subroutine wtcoef(f,f1,f2,f3,f4,wt)
            !-----------------------------------------------------
            implicit none

            real(fpp), intent(in) ::  f,f1,f2,f3,f4
            real(fpp), intent(out)::  wt

            if (f3.gt.f4) stop 'wtcoef: f3>f4 '
            if (f1.gt.f2) stop 'wtcoef: f1>f2 '
            if (f.le.f3.and.f.ge.f2) then
                wt = 1.
            else if (f.gt.f4.or.f.lt.f1 ) then
                wt = 0.
            else if (f.gt.f3.and.f.le.f4) then
                wt = 0.5*(1.0+cos(M_PI*(f-f3)/(f4-f3)))
            else if (f.ge.f1.and.f.lt.f2) then
                wt = 0.5*(1.0+cos(M_PI*(f-f2)/(f2-f1)))
            endif
            !-----------------------------------------------------
        end subroutine wtcoef
        !-----------------------------------------------------
    end subroutine def_timefunc

    ! ########################################################
    subroutine dfour1(u,n,isign)
        use constants, only : fpp
        integer n,isign
        complex(fpp) :: u(n)

        integer :: inc,jump,lot,i
        real(fpp) :: a(n),b(n),trigs(2*n)

        inc = 1
        jump = 1
        lot = 1

        call setgpfa(trigs,n)

        do i = 1,n
            a(i) = real(u(i),fpp)
            b(i) = aimag(u(i))
        enddo

        call gpfa(a,b,trigs,inc,jump,n,lot,isign)

        do i = 1,n
            u(i) = dcmplx(a(i),b(i))
        enddo

        return
    end subroutine dfour1

    !        SUBROUTINE 'SETGPFA'
    !        SETUP ROUTINE FOR SELF-SORTING IN-PLACE
    !            GENERALIZED PRIME FACTOR (COMPLEX) FFT [GPFA]
    !
    !        CALL SETGPFA(TRIGS,N)
    !
    !        INPUT :
    !        -----
    !        N IS THE LENGTH OF THE TRANSFORMS. N MUST BE OF THE FORM:
    !          -----------------------------------
    !            N = (2**IP) * (3**IQ) * (5**IR)
    !          -----------------------------------
    !
    !        OUTPUT:
    !        ------
    !        TRIGS IS A TABLE OF TWIDDLE FACTORS,
    !          OF LENGTH 2*IPQR (REAL) WORDS, WHERE:
    !          --------------------------------------
    !            IPQR = (2**IP) + (3**IQ) + (5**IR)
    !          --------------------------------------
    !
    !        WRITTEN BY CLIVE TEMPERTON 1990
    !
    !----------------------------------------------------------------------

    SUBROUTINE SETGPFA(TRIGS,N)
        integer :: N, NN, IFAC
        real(fpp) :: TRIGS(*)
        real(fpp) :: TWOPI,ANGLE,DEL

        integer, dimension(3) :: NJ
        integer :: LL, KK, IP, IQ, IR, NI, IROT, KINK, I, K
        !     DECOMPOSE N INTO FACTORS 2,3,5
        !     ------------------------------
        NN = N
        IFAC = 2

        DO  LL = 1 , 3
            KK = 0
10          CONTINUE
            IF (MOD(NN,IFAC).NE.0) GO TO 20
            KK = KK + 1
            NN = NN / IFAC
            GO TO 10
20          CONTINUE
            NJ(LL) = KK
            IFAC = IFAC + LL
        ENDDO

        IF (NN.NE.1) THEN
            WRITE(6,40) N
40          FORMAT(' *** WARNING!!!',I10,' IS NOT A LEGAL VALUE OF N ***')
            RETURN
        ENDIF

        IP = NJ(1)
        IQ = NJ(2)
        IR = NJ(3)

        !     COMPUTE LIST OF ROTATED TWIDDLE FACTORS
        !     ---------------------------------------
        NJ(1) = 2**IP
        NJ(2) = 3**IQ
        NJ(3) = 5**IR

        TWOPI = 4.0 * ASIN(1.0)
        I = 1

        DO LL = 1 , 3
            NI = NJ(LL)
            IF (NI.EQ.1) GO TO 60

            DEL = TWOPI / FLOAT(NI)
            IROT = N / NI
            KINK = MOD(IROT,NI)
            KK = 0

            DO K = 1 , NI
                ANGLE = FLOAT(KK) * DEL
                TRIGS(I) = COS(ANGLE)
                TRIGS(I+1) = SIN(ANGLE)
                I = I + 2
                KK = KK + KINK
                IF (KK.GT.NI) KK = KK - NI
            ENDDO
60          CONTINUE
        ENDDO

        RETURN
    END SUBROUTINE SETGPFA


    !        SUBROUTINE 'GPFA'
    !        SELF-SORTING IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT
    !
    !        *** THIS IS THE ALL-FORTRAN VERSION ***
    !            -------------------------------
    !
    !        CALL GPFA(A,B,TRIGS,INC,JUMP,N,LOT,ISIGN)
    !
    !        A IS FIRST REAL INPUT/OUTPUT VECTOR
    !        B IS FIRST IMAGINARY INPUT/OUTPUT VECTOR
    !        TRIGS IS A TABLE OF TWIDDLE FACTORS, PRECALCULATED
    !              BY CALLING SUBROUTINE 'SETGPFA'
    !        INC IS THE INCREMENT WITHIN EACH DATA VECTOR
    !        JUMP IS THE INCREMENT BETWEEN DATA VECTORS
    !        N IS THE LENGTH OF THE TRANSFORMS:
    !          -----------------------------------
    !            N = (2**IP) * (3**IQ) * (5**IR)
    !          -----------------------------------
    !        LOT IS THE NUMBER OF TRANSFORMS
    !        ISIGN = +1 FOR FORWARD TRANSFORM
    !              = -1 FOR INVERSE TRANSFORM
    !
    !        WRITTEN BY CLIVE TEMPERTON
    !        RECHERCHE EN PREVISION NUMERIQUE
    !        ATMOSPHERIC ENVIRONMENT SERVICE, CANADA
    !
    !----------------------------------------------------------------------
    !
    !        DEFINITION OF TRANSFORM
    !        -----------------------
    !
    !        X(J) = SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N))
    !
    !---------------------------------------------------------------------
    !
    !        FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
    !        SEE:
    !
    !        C TEMPERTON : "A GENERALIZED PRIME FACTOR FFT ALGORITHM
    !          FOR ANY N = (2**P)(3**Q)(5**R)",
    !          SIAM J. SCI. STAT. COMP., MAY 1992.
    !
    !----------------------------------------------------------------------

    SUBROUTINE GPFA(A,B,TRIGS,INC,JUMP,N,LOT,ISIGN)
        integer :: INC,JUMP,N,LOT,ISIGN
        real(fpp) ::  A(*), B(*), TRIGS(*)
        integer, DIMENSION(3) :: NJ
        !
        integer :: NN, IFAC, LL, KK, IP, IQ, IR, I
        !     DECOMPOSE N INTO FACTORS 2,3,5
        !     ------------------------------
        NN = N
        IFAC = 2

        DO LL = 1 , 3
            KK = 0
10          CONTINUE
            IF (MOD(NN,IFAC).NE.0) GO TO 20
            KK = KK + 1
            NN = NN / IFAC
            GO TO 10
20          CONTINUE
            NJ(LL) = KK
            IFAC = IFAC + LL
        ENDDO

        IF (NN.NE.1) THEN
            WRITE(6,40) N
40          FORMAT(' *** WARNING!!!',I10,' IS NOT A LEGAL VALUE OF N ***')
            RETURN
        ENDIF

        IP = NJ(1)
        IQ = NJ(2)
        IR = NJ(3)

        !     COMPUTE THE TRANSFORM
        !     ---------------------
        I = 1
        IF (IP.GT.0) THEN
            CALL GPFA2F(A,B,TRIGS,INC,JUMP,N,IP,LOT,ISIGN)
            I = I + 2 * ( 2**IP)
        ENDIF
        IF (IQ.GT.0) THEN
            CALL GPFA3F(A,B,TRIGS(I),INC,JUMP,N,IQ,LOT,ISIGN)
            I = I + 2 * (3**IQ)
        ENDIF
        IF (IR.GT.0) THEN
            CALL GPFA5F(A,B,TRIGS(I),INC,JUMP,N,IR,LOT,ISIGN)
        ENDIF

        RETURN
    END SUBROUTINE GPFA




    !     fortran version of *gpfa2* -
    !     radix-2 section of self-sorting, in-place, generalized pfa
    !     central radix-2 and radix-8 passes included
    !      so that transform length can be any power of 2
    !
    !-------------------------------------------------------------------
    !
    subroutine gpfa2f(a,b,trigs,inc,jump,n,mm,lot,isign)
        integer :: inc,jump,n,lot,mm,isign
        real(fpp) :: a(*), b(*), trigs(*)

        real(fpp) :: s,ss
        real(fpp) :: aja,ajc,t0,t2,ajb,ajd,t1,t3
        real(fpp) :: bja,bjc,u0,u2,bjb,bjd,u1,u3
        real(fpp) :: co1,si1,co2,si2,co3,si3
        real(fpp) :: c1,c2,c3
        real(fpp) :: aje,ajg,ajf,ajh
        real(fpp) :: bje,bjg,bjf,bjh
        real(fpp) :: co4,si4,co5,si5,co6,si6,co7,si7
        real(fpp) :: aji,ajj,ajk,ajl,ajm,ajn,ajo,ajp
        real(fpp) :: bji,bjj,bjk,bjl,bjm,bjn,bjo,bjp

        integer :: ink,inq,ipass,istart,j,jjj,jstep,jstepl,jstepx
        integer :: ja,jb,jc,jd,je,jf,jg,ji,jj,jh,jk,jl,jm,jn,jo,jp
        integer :: k,kk,l,la,laincl,left,ll,m,mh,mu,nb,nblox,ninc,nu,nvex
        integer :: m2,m8,n2

        n2 = 2**mm
        inq = n/n2
        jstepx = (n2-n) * inc
        ninc = n * inc
        ink = inc * inq

        m2 = 0
        m8 = 0
        if (mod(mm,2).eq.0) then
            m = mm/2
        else if (mod(mm,4).eq.1) then
            m = (mm-1)/2
            m2 = 1
        else if (mod(mm,4).eq.3) then
            m = (mm-3)/2
            m8 = 1
        endif
        mh = (m+1)/2

        nblox = 1 + (lot-1)/lvr
        left = lot
        s = float(isign)
        istart = 1

        !  loop on blocks of lvr transforms
        !  --------------------------------
        do nb = 1 , nblox

            if (left.le.lvr) then
                nvex = left
            else if (left.lt.(2*lvr)) then
                nvex = left/2
                nvex = nvex + mod(nvex,2)
            else
                nvex = lvr
            endif
            left = left - nvex

            la = 1

            !  loop on type I radix-4 passes
            !  -----------------------------
            mu = mod(inq,4)
            if (isign.eq.-1) mu = 4 - mu
            ss = 1.0
            if (mu.eq.3) ss = -1.0

            if (mh.eq.0) go to 200

            do ipass = 1 , mh
                jstep = (n*inc) / (4*la)
                jstepl = jstep - ninc

                !  k = 0 loop (no twiddle factors)
                !  -------------------------------
                do jjj = 0 , (n-1)*inc , 4*jstep
                    ja = istart + jjj

                    !     "transverse" loop
                    !     -----------------
                    do nu = 1 , inq
                        jb = ja + jstepl
                        if (jb.lt.istart) jb = jb + ninc
                        jc = jb + jstepl
                        if (jc.lt.istart) jc = jc + ninc
                        jd = jc + jstepl
                        if (jd.lt.istart) jd = jd + ninc
                        j = 0

                        !  loop across transforms
                        !  ----------------------
                        !cdiri$ ivdep, shortloop
                        do l = 1 , nvex
                            aja = a(ja+j)
                            ajc = a(jc+j)
                            t0 = aja + ajc
                            t2 = aja - ajc
                            ajb = a(jb+j)
                            ajd = a(jd+j)
                            t1 = ajb + ajd
                            t3 = ss * ( ajb - ajd )
                            bja = b(ja+j)
                            bjc = b(jc+j)
                            u0 = bja + bjc
                            u2 = bja - bjc
                            bjb = b(jb+j)
                            bjd = b(jd+j)
                            u1 = bjb + bjd
                            u3 = ss * ( bjb - bjd )
                            a(ja+j) = t0 + t1
                            a(jc+j) = t0 - t1
                            b(ja+j) = u0 + u1
                            b(jc+j) = u0 - u1
                            a(jb+j) = t2 - u3
                            a(jd+j) = t2 + u3
                            b(jb+j) = u2 + t3
                            b(jd+j) = u2 - t3
                            j = j + jump
                        enddo
                        ja = ja + jstepx
                        if (ja.lt.istart) ja = ja + ninc
                    enddo
                enddo

                !  finished if n2 = 4
                !  ------------------
                if (n2.eq.4) go to 490
                kk = 2 * la

                !  loop on nonzero k
                !  -----------------
                do k = ink , jstep-ink , ink
                    co1 = trigs(kk+1)
                    si1 = s*trigs(kk+2)
                    co2 = trigs(2*kk+1)
                    si2 = s*trigs(2*kk+2)
                    co3 = trigs(3*kk+1)
                    si3 = s*trigs(3*kk+2)

                    !  loop along transform
                    !  --------------------
                    do jjj = k , (n-1)*inc , 4*jstep
                        ja = istart + jjj

                        !     "transverse" loop
                        !     -----------------
                        do nu = 1 , inq
                            jb = ja + jstepl
                            if (jb.lt.istart) jb = jb + ninc
                            jc = jb + jstepl
                            if (jc.lt.istart) jc = jc + ninc
                            jd = jc + jstepl
                            if (jd.lt.istart) jd = jd + ninc
                            j = 0

                            !  loop across transforms
                            !  ----------------------
                            !cdiri$ ivdep,shortloop
                            do l = 1 , nvex
                                aja = a(ja+j)
                                ajc = a(jc+j)
                                t0 = aja + ajc
                                t2 = aja - ajc
                                ajb = a(jb+j)
                                ajd = a(jd+j)
                                t1 = ajb + ajd
                                t3 = ss * ( ajb - ajd )
                                bja = b(ja+j)
                                bjc = b(jc+j)
                                u0 = bja + bjc
                                u2 = bja - bjc
                                bjb = b(jb+j)
                                bjd = b(jd+j)
                                u1 = bjb + bjd
                                u3 = ss * ( bjb - bjd )
                                a(ja+j) = t0 + t1
                                b(ja+j) = u0 + u1
                                a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
                                b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
                                a(jc+j) = co2*(t0-t1) - si2*(u0-u1)
                                b(jc+j) = si2*(t0-t1) + co2*(u0-u1)
                                a(jd+j) = co3*(t2+u3) - si3*(u2-t3)
                                b(jd+j) = si3*(t2+u3) + co3*(u2-t3)
                                j = j + jump
                            enddo
                            !-----( end of loop across transforms )
                            ja = ja + jstepx
                            if (ja.lt.istart) ja = ja + ninc
                        enddo
                    enddo
                    !-----( end of loop along transforms )
                    kk = kk + 2*la
                enddo
                !-----( end of loop on nonzero k )
                la = 4*la
            enddo
            !-----( end of loop on type I radix-4 passes)

            !  central radix-2 pass
            !  --------------------
200         continue
            if (m2.eq.0) go to 300

            jstep = (n*inc) / (2*la)
            jstepl = jstep - ninc

            !  k=0 loop (no twiddle factors)
            !  -----------------------------
            do jjj = 0 , (n-1)*inc , 2*jstep
                ja = istart + jjj

                !     "transverse" loop
                !     -----------------
                do nu = 1 , inq
                    jb = ja + jstepl
                    if (jb.lt.istart) jb = jb + ninc
                    j = 0

                    !  loop across transforms
                    !  ----------------------
                    !cdiri$ ivdep, shortloop
                    do l = 1 , nvex
                        aja = a(ja+j)
                        ajb = a(jb+j)
                        t0 = aja - ajb
                        a(ja+j) = aja + ajb
                        a(jb+j) = t0
                        bja = b(ja+j)
                        bjb = b(jb+j)
                        u0 = bja - bjb
                        b(ja+j) = bja + bjb
                        b(jb+j) = u0
                        j = j + jump
                    enddo
                    !-----(end of loop across transforms)
                    ja = ja + jstepx
                    if (ja.lt.istart) ja = ja + ninc
                enddo
            enddo

            !  finished if n2=2
            !  ----------------
            if (n2.eq.2) go to 490

            kk = 2 * la

            !  loop on nonzero k
            !  -----------------
            do k = ink , jstep - ink , ink
                co1 = trigs(kk+1)
                si1 = s*trigs(kk+2)

                !  loop along transforms
                !  ---------------------
                do jjj = k , (n-1)*inc , 2*jstep
                    ja = istart + jjj

                    !     "transverse" loop
                    !     -----------------
                    do nu = 1 , inq
                        jb = ja + jstepl
                        if (jb.lt.istart) jb = jb + ninc
                        j = 0

                        !  loop across transforms
                        !  ----------------------
                        if (kk.eq.n2/2) then
                            !cdiri$ ivdep, shortloop
                            do l = 1 , nvex
                                aja = a(ja+j)
                                ajb = a(jb+j)
                                t0 = ss * ( aja - ajb )
                                a(ja+j) = aja + ajb
                                bjb = b(jb+j)
                                bja = b(ja+j)
                                a(jb+j) = ss * ( bjb - bja )
                                b(ja+j) = bja + bjb
                                b(jb+j) = t0
                                j = j + jump
                            enddo

                        else

                            !cdiri$ ivdep, shortloop
                            do l = 1 , nvex
                                aja = a(ja+j)
                                ajb = a(jb+j)
                                t0 = aja - ajb
                                a(ja+j) = aja + ajb
                                bja = b(ja+j)
                                bjb = b(jb+j)
                                u0 = bja - bjb
                                b(ja+j) = bja + bjb
                                a(jb+j) = co1*t0 - si1*u0
                                b(jb+j) = si1*t0 + co1*u0
                                j = j + jump
                            enddo

                        endif

                        !*-----(end of loop across transforms)
                        ja = ja + jstepx
                        if (ja.lt.istart) ja = ja + ninc
                    enddo
                enddo
                !-----(end of loop along transforms)
                kk = kk + 2 * la
            enddo
            !-----(end of loop on nonzero k)
            !-----(end of radix-2 pass)
            !
            la = 2 * la
            go to 400
            !
            !  central radix-8 pass
            !  --------------------
300         continue
            if (m8.eq.0) go to 400
            jstep = (n*inc) / (8*la)
            jstepl = jstep - ninc
            mu = mod(inq,8)
            if (isign.eq.-1) mu = 8 - mu
            c1 = 1.0
            if (mu.eq.3.or.mu.eq.7) c1 = -1.0
            c2 = sqrt(0.5)
            if (mu.eq.3.or.mu.eq.5) c2 = -c2
            c3 = c1 * c2
            !
            !  stage 1
            !  -------
            do k = 0 , jstep - ink , ink
                do jjj = k , (n-1)*inc , 8*jstep
                    ja = istart + jjj
                    !
                    !     "transverse" loop
                    !     -----------------
                    do nu = 1 , inq
                        jb = ja + jstepl
                        if (jb.lt.istart) jb = jb + ninc
                        jc = jb + jstepl
                        if (jc.lt.istart) jc = jc + ninc
                        jd = jc + jstepl
                        if (jd.lt.istart) jd = jd + ninc
                        je = jd + jstepl
                        if (je.lt.istart) je = je + ninc
                        jf = je + jstepl
                        if (jf.lt.istart) jf = jf + ninc
                        jg = jf + jstepl
                        if (jg.lt.istart) jg = jg + ninc
                        jh = jg + jstepl
                        if (jh.lt.istart) jh = jh + ninc
                        j = 0
                        !cdiri$ ivdep, shortloop
                        do l = 1 , nvex
                            aja = a(ja+j)
                            aje = a(je+j)
                            t0 = aja - aje
                            a(ja+j) = aja + aje
                            ajc = a(jc+j)
                            ajg = a(jg+j)
                            t1 = c1 * ( ajc - ajg )
                            a(je+j) = ajc + ajg
                            ajb = a(jb+j)
                            ajf = a(jf+j)
                            t2 = ajb - ajf
                            a(jc+j) = ajb + ajf
                            ajd = a(jd+j)
                            ajh = a(jh+j)
                            t3 = ajd - ajh
                            a(jg+j) = ajd + ajh
                            a(jb+j) = t0
                            a(jf+j) = t1
                            a(jd+j) = c2 * ( t2 - t3 )
                            a(jh+j) = c3 * ( t2 + t3 )
                            bja = b(ja+j)
                            bje = b(je+j)
                            u0 = bja - bje
                            b(ja+j) = bja + bje
                            bjc = b(jc+j)
                            bjg = b(jg+j)
                            u1 = c1 * ( bjc - bjg )
                            b(je+j) = bjc + bjg
                            bjb = b(jb+j)
                            bjf = b(jf+j)
                            u2 = bjb - bjf
                            b(jc+j) = bjb + bjf
                            bjd = b(jd+j)
                            bjh = b(jh+j)
                            u3 = bjd - bjh
                            b(jg+j) = bjd + bjh
                            b(jb+j) = u0
                            b(jf+j) = u1
                            b(jd+j) = c2 * ( u2 - u3 )
                            b(jh+j) = c3 * ( u2 + u3 )
                            j = j + jump
                        enddo
                        ja = ja + jstepx
                        if (ja.lt.istart) ja = ja + ninc
                    enddo
                enddo
            enddo
            !
            !  stage 2
            !  -------
            !
            !  k=0 (no twiddle factors)
            !  ------------------------
            do jjj = 0 , (n-1)*inc , 8*jstep
                ja = istart + jjj
                !
                !     "transverse" loop
                !     -----------------
                do nu = 1 , inq
                    jb = ja + jstepl
                    if (jb.lt.istart) jb = jb + ninc
                    jc = jb + jstepl
                    if (jc.lt.istart) jc = jc + ninc
                    jd = jc + jstepl
                    if (jd.lt.istart) jd = jd + ninc
                    je = jd + jstepl
                    if (je.lt.istart) je = je + ninc
                    jf = je + jstepl
                    if (jf.lt.istart) jf = jf + ninc
                    jg = jf + jstepl
                    if (jg.lt.istart) jg = jg + ninc
                    jh = jg + jstepl
                    if (jh.lt.istart) jh = jh + ninc
                    j = 0
                    !cdiri$ ivdep, shortloop
                    do l = 1 , nvex
                        aja = a(ja+j)
                        aje = a(je+j)
                        t0 = aja + aje
                        t2 = aja - aje
                        ajc = a(jc+j)
                        ajg = a(jg+j)
                        t1 = ajc + ajg
                        t3 = c1 * ( ajc - ajg )
                        bja = b(ja+j)
                        bje = b(je+j)
                        u0 = bja + bje
                        u2 = bja - bje
                        bjc = b(jc+j)
                        bjg = b(jg+j)
                        u1 = bjc + bjg
                        u3 = c1 * ( bjc - bjg )
                        a(ja+j) = t0 + t1
                        a(je+j) = t0 - t1
                        b(ja+j) = u0 + u1
                        b(je+j) = u0 - u1
                        a(jc+j) = t2 - u3
                        a(jg+j) = t2 + u3
                        b(jc+j) = u2 + t3
                        b(jg+j) = u2 - t3
                        ajb = a(jb+j)
                        ajd = a(jd+j)
                        t0 = ajb + ajd
                        t2 = ajb - ajd
                        ajf = a(jf+j)
                        ajh = a(jh+j)
                        t1 = ajf - ajh
                        t3 = ajf + ajh
                        bjb = b(jb+j)
                        bjd = b(jd+j)
                        u0 = bjb + bjd
                        u2 = bjb - bjd
                        bjf = b(jf+j)
                        bjh = b(jh+j)
                        u1 = bjf - bjh
                        u3 = bjf + bjh
                        a(jb+j) = t0 - u3
                        a(jh+j) = t0 + u3
                        b(jb+j) = u0 + t3
                        b(jh+j) = u0 - t3
                        a(jd+j) = t2 + u1
                        a(jf+j) = t2 - u1
                        b(jd+j) = u2 - t1
                        b(jf+j) = u2 + t1
                        j = j + jump
                    enddo
                    ja = ja + jstepx
                    if (ja.lt.istart) ja = ja + ninc
                enddo
            enddo
            !
            if (n2.eq.8) go to 490
            !
            !  loop on nonzero k
            !  -----------------
            kk = 2 * la
            !
            do k = ink , jstep - ink , ink
                !
                co1 = trigs(kk+1)
                si1 = s * trigs(kk+2)
                co2 = trigs(2*kk+1)
                si2 = s * trigs(2*kk+2)
                co3 = trigs(3*kk+1)
                si3 = s * trigs(3*kk+2)
                co4 = trigs(4*kk+1)
                si4 = s * trigs(4*kk+2)
                co5 = trigs(5*kk+1)
                si5 = s * trigs(5*kk+2)
                co6 = trigs(6*kk+1)
                si6 = s * trigs(6*kk+2)
                co7 = trigs(7*kk+1)
                si7 = s * trigs(7*kk+2)
                !
                do jjj = k , (n-1)*inc , 8*jstep
                    ja = istart + jjj
                    !
                    !     "transverse" loop
                    !     -----------------
                    do nu = 1 , inq
                        jb = ja + jstepl
                        if (jb.lt.istart) jb = jb + ninc
                        jc = jb + jstepl
                        if (jc.lt.istart) jc = jc + ninc
                        jd = jc + jstepl
                        if (jd.lt.istart) jd = jd + ninc
                        je = jd + jstepl
                        if (je.lt.istart) je = je + ninc
                        jf = je + jstepl
                        if (jf.lt.istart) jf = jf + ninc
                        jg = jf + jstepl
                        if (jg.lt.istart) jg = jg + ninc
                        jh = jg + jstepl
                        if (jh.lt.istart) jh = jh + ninc
                        j = 0
                        !cdiri$ ivdep, shortloop
                        do l = 1 , nvex
                            aja = a(ja+j)
                            aje = a(je+j)
                            t0 = aja + aje
                            t2 = aja - aje
                            ajc = a(jc+j)
                            ajg = a(jg+j)
                            t1 = ajc + ajg
                            t3 = c1 * ( ajc - ajg )
                            bja = b(ja+j)
                            bje = b(je+j)
                            u0 = bja + bje
                            u2 = bja - bje
                            bjc = b(jc+j)
                            bjg = b(jg+j)
                            u1 = bjc + bjg
                            u3 = c1 * ( bjc - bjg )
                            a(ja+j) = t0 + t1
                            b(ja+j) = u0 + u1
                            a(je+j) = co4*(t0-t1) - si4*(u0-u1)
                            b(je+j) = si4*(t0-t1) + co4*(u0-u1)
                            a(jc+j) = co2*(t2-u3) - si2*(u2+t3)
                            b(jc+j) = si2*(t2-u3) + co2*(u2+t3)
                            a(jg+j) = co6*(t2+u3) - si6*(u2-t3)
                            b(jg+j) = si6*(t2+u3) + co6*(u2-t3)
                            ajb = a(jb+j)
                            ajd = a(jd+j)
                            t0 = ajb + ajd
                            t2 = ajb - ajd
                            ajf = a(jf+j)
                            ajh = a(jh+j)
                            t1 = ajf - ajh
                            t3 = ajf + ajh
                            bjb = b(jb+j)
                            bjd = b(jd+j)
                            u0 = bjb + bjd
                            u2 = bjb - bjd
                            bjf = b(jf+j)
                            bjh = b(jh+j)
                            u1 = bjf - bjh
                            u3 = bjf + bjh
                            a(jb+j) = co1*(t0-u3) - si1*(u0+t3)
                            b(jb+j) = si1*(t0-u3) + co1*(u0+t3)
                            a(jh+j) = co7*(t0+u3) - si7*(u0-t3)
                            b(jh+j) = si7*(t0+u3) + co7*(u0-t3)
                            a(jd+j) = co3*(t2+u1) - si3*(u2-t1)
                            b(jd+j) = si3*(t2+u1) + co3*(u2-t1)
                            a(jf+j) = co5*(t2-u1) - si5*(u2+t1)
                            b(jf+j) = si5*(t2-u1) + co5*(u2+t1)
                            j = j + jump
                        enddo
                        ja = ja + jstepx
                        if (ja.lt.istart) ja = ja + ninc
                    enddo
                enddo
                kk = kk + 2 * la
            enddo
            !
            la = 8 * la
            !
            !  loop on type II radix-4 passes
            !  ------------------------------
400         continue
            mu = mod(inq,4)
            if (isign.eq.-1) mu = 4 - mu
            ss = 1.0
            if (mu.eq.3) ss = -1.0
            !
            do ipass = mh+1 , m
                jstep = (n*inc) / (4*la)
                jstepl = jstep - ninc
                laincl = la * ink - ninc
                !
                !  k=0 loop (no twiddle factors)
                !  -----------------------------
                do ll = 0 , (la-1)*ink , 4*jstep
                    !
                    do jjj = ll , (n-1)*inc , 4*la*ink
                        ja = istart + jjj
                        !
                        !     "transverse" loop
                        !     -----------------
                        do nu = 1 , inq
                            jb = ja + jstepl
                            if (jb.lt.istart) jb = jb + ninc
                            jc = jb + jstepl
                            if (jc.lt.istart) jc = jc + ninc
                            jd = jc + jstepl
                            if (jd.lt.istart) jd = jd + ninc
                            je = ja + laincl
                            if (je.lt.istart) je = je + ninc
                            jf = je + jstepl
                            if (jf.lt.istart) jf = jf + ninc
                            jg = jf + jstepl
                            if (jg.lt.istart) jg = jg + ninc
                            jh = jg + jstepl
                            if (jh.lt.istart) jh = jh + ninc
                            ji = je + laincl
                            if (ji.lt.istart) ji = ji + ninc
                            jj = ji + jstepl
                            if (jj.lt.istart) jj = jj + ninc
                            jk = jj + jstepl
                            if (jk.lt.istart) jk = jk + ninc
                            jl = jk + jstepl
                            if (jl.lt.istart) jl = jl + ninc
                            jm = ji + laincl
                            if (jm.lt.istart) jm = jm + ninc
                            jn = jm + jstepl
                            if (jn.lt.istart) jn = jn + ninc
                            jo = jn + jstepl
                            if (jo.lt.istart) jo = jo + ninc
                            jp = jo + jstepl
                            if (jp.lt.istart) jp = jp + ninc
                            j = 0
                            !
                            !  loop across transforms
                            !  ----------------------
                            !cdiri$ ivdep, shortloop
                            do l = 1 , nvex
                                aja = a(ja+j)
                                ajc = a(jc+j)
                                t0 = aja + ajc
                                t2 = aja - ajc
                                ajb = a(jb+j)
                                ajd = a(jd+j)
                                t1 = ajb + ajd
                                t3 = ss * ( ajb - ajd )
                                aji = a(ji+j)
                                ajc =  aji
                                bja = b(ja+j)
                                bjc = b(jc+j)
                                u0 = bja + bjc
                                u2 = bja - bjc
                                bjb = b(jb+j)
                                bjd = b(jd+j)
                                u1 = bjb + bjd
                                u3 = ss * ( bjb - bjd )
                                aje = a(je+j)
                                ajb =  aje
                                a(ja+j) = t0 + t1
                                a(ji+j) = t0 - t1
                                b(ja+j) = u0 + u1
                                bjc =  u0 - u1
                                bjm = b(jm+j)
                                bjd =  bjm
                                a(je+j) = t2 - u3
                                ajd =  t2 + u3
                                bjb =  u2 + t3
                                b(jm+j) = u2 - t3
                                !----------------------
                                ajg = a(jg+j)
                                t0 = ajb + ajg
                                t2 = ajb - ajg
                                ajf = a(jf+j)
                                ajh = a(jh+j)
                                t1 = ajf + ajh
                                t3 = ss * ( ajf - ajh )
                                ajj = a(jj+j)
                                ajg =  ajj
                                bje = b(je+j)
                                bjg = b(jg+j)
                                u0 = bje + bjg
                                u2 = bje - bjg
                                bjf = b(jf+j)
                                bjh = b(jh+j)
                                u1 = bjf + bjh
                                u3 = ss * ( bjf - bjh )
                                b(je+j) = bjb
                                a(jb+j) = t0 + t1
                                a(jj+j) = t0 - t1
                                bjj = b(jj+j)
                                bjg =  bjj
                                b(jb+j) = u0 + u1
                                b(jj+j) = u0 - u1
                                a(jf+j) = t2 - u3
                                ajh =  t2 + u3
                                b(jf+j) = u2 + t3
                                bjh =  u2 - t3
                                !----------------------
                                ajk = a(jk+j)
                                t0 = ajc + ajk
                                t2 = ajc - ajk
                                ajl = a(jl+j)
                                t1 = ajg + ajl
                                t3 = ss * ( ajg - ajl )
                                bji = b(ji+j)
                                bjk = b(jk+j)
                                u0 = bji + bjk
                                u2 = bji - bjk
                                ajo = a(jo+j)
                                ajl =  ajo
                                bjl = b(jl+j)
                                u1 = bjg + bjl
                                u3 = ss * ( bjg - bjl )
                                b(ji+j) = bjc
                                a(jc+j) = t0 + t1
                                a(jk+j) = t0 - t1
                                bjo = b(jo+j)
                                bjl =  bjo
                                b(jc+j) = u0 + u1
                                b(jk+j) = u0 - u1
                                a(jg+j) = t2 - u3
                                a(jo+j) = t2 + u3
                                b(jg+j) = u2 + t3
                                b(jo+j) = u2 - t3
                                !----------------------
                                ajm = a(jm+j)
                                t0 = ajm + ajl
                                t2 = ajm - ajl
                                ajn = a(jn+j)
                                ajp = a(jp+j)
                                t1 = ajn + ajp
                                t3 = ss * ( ajn - ajp )
                                a(jm+j) = ajd
                                u0 = bjd + bjl
                                u2 = bjd - bjl
                                bjn = b(jn+j)
                                bjp = b(jp+j)
                                u1 = bjn + bjp
                                u3 = ss * ( bjn - bjp )
                                a(jn+j) = ajh
                                a(jd+j) = t0 + t1
                                a(jl+j) = t0 - t1
                                b(jd+j) = u0 + u1
                                b(jl+j) = u0 - u1
                                b(jn+j) = bjh
                                a(jh+j) = t2 - u3
                                a(jp+j) = t2 + u3
                                b(jh+j) = u2 + t3
                                b(jp+j) = u2 - t3
                                j = j + jump
                            enddo
                            !-----( end of loop across transforms )
                            ja = ja + jstepx
                            if (ja.lt.istart) ja = ja + ninc
                        enddo
                    enddo
                enddo
                !-----( end of double loop for k=0 )
                !
                !  finished if last pass
                !  ---------------------
                if (ipass.eq.m) go to 490
                !
                kk = 2*la
                !
                !     loop on nonzero k
                !     -----------------
                do k = ink , jstep-ink , ink
                    co1 = trigs(kk+1)
                    si1 = s*trigs(kk+2)
                    co2 = trigs(2*kk+1)
                    si2 = s*trigs(2*kk+2)
                    co3 = trigs(3*kk+1)
                    si3 = s*trigs(3*kk+2)
                    !
                    !  double loop along first transform in block
                    !  ------------------------------------------
                    do ll = k , (la-1)*ink , 4*jstep
                        !
                        do jjj = ll , (n-1)*inc , 4*la*ink
                            ja = istart + jjj
                            !
                            !     "transverse" loop
                            !     -----------------
                            do nu = 1 , inq
                                jb = ja + jstepl
                                if (jb.lt.istart) jb = jb + ninc
                                jc = jb + jstepl
                                if (jc.lt.istart) jc = jc + ninc
                                jd = jc + jstepl
                                if (jd.lt.istart) jd = jd + ninc
                                je = ja + laincl
                                if (je.lt.istart) je = je + ninc
                                jf = je + jstepl
                                if (jf.lt.istart) jf = jf + ninc
                                jg = jf + jstepl
                                if (jg.lt.istart) jg = jg + ninc
                                jh = jg + jstepl
                                if (jh.lt.istart) jh = jh + ninc
                                ji = je + laincl
                                if (ji.lt.istart) ji = ji + ninc
                                jj = ji + jstepl
                                if (jj.lt.istart) jj = jj + ninc
                                jk = jj + jstepl
                                if (jk.lt.istart) jk = jk + ninc
                                jl = jk + jstepl
                                if (jl.lt.istart) jl = jl + ninc
                                jm = ji + laincl
                                if (jm.lt.istart) jm = jm + ninc
                                jn = jm + jstepl
                                if (jn.lt.istart) jn = jn + ninc
                                jo = jn + jstepl
                                if (jo.lt.istart) jo = jo + ninc
                                jp = jo + jstepl
                                if (jp.lt.istart) jp = jp + ninc
                                j = 0
                                !
                                !  loop across transforms
                                !  ----------------------
                                !cdiri$ ivdep, shortloop
                                do l = 1 , nvex
                                    aja = a(ja+j)
                                    ajc = a(jc+j)
                                    t0 = aja + ajc
                                    t2 = aja - ajc
                                    ajb = a(jb+j)
                                    ajd = a(jd+j)
                                    t1 = ajb + ajd
                                    t3 = ss * ( ajb - ajd )
                                    aji = a(ji+j)
                                    ajc =  aji
                                    bja = b(ja+j)
                                    bjc = b(jc+j)
                                    u0 = bja + bjc
                                    u2 = bja - bjc
                                    bjb = b(jb+j)
                                    bjd = b(jd+j)
                                    u1 = bjb + bjd
                                    u3 = ss * ( bjb - bjd )
                                    aje = a(je+j)
                                    ajb =  aje
                                    a(ja+j) = t0 + t1
                                    b(ja+j) = u0 + u1
                                    a(je+j) = co1*(t2-u3) - si1*(u2+t3)
                                    bjb =  si1*(t2-u3) + co1*(u2+t3)
                                    bjm = b(jm+j)
                                    bjd =  bjm
                                    a(ji+j) = co2*(t0-t1) - si2*(u0-u1)
                                    bjc =  si2*(t0-t1) + co2*(u0-u1)
                                    ajd =  co3*(t2+u3) - si3*(u2-t3)
                                    b(jm+j) = si3*(t2+u3) + co3*(u2-t3)
                                    !----------------------------------------
                                    ajg = a(jg+j)
                                    t0 = ajb + ajg
                                    t2 = ajb - ajg
                                    ajf = a(jf+j)
                                    ajh = a(jh+j)
                                    t1 = ajf + ajh
                                    t3 = ss * ( ajf - ajh )
                                    ajj = a(jj+j)
                                    ajg =  ajj
                                    bje = b(je+j)
                                    bjg = b(jg+j)
                                    u0 = bje + bjg
                                    u2 = bje - bjg
                                    bjf = b(jf+j)
                                    bjh = b(jh+j)
                                    u1 = bjf + bjh
                                    u3 = ss * ( bjf - bjh )
                                    b(je+j) = bjb
                                    a(jb+j) = t0 + t1
                                    b(jb+j) = u0 + u1
                                    bjj = b(jj+j)
                                    bjg =  bjj
                                    a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
                                    b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
                                    a(jj+j) = co2*(t0-t1) - si2*(u0-u1)
                                    b(jj+j) = si2*(t0-t1) + co2*(u0-u1)
                                    ajh =  co3*(t2+u3) - si3*(u2-t3)
                                    bjh =  si3*(t2+u3) + co3*(u2-t3)
                                    !----------------------------------------
                                    ajk = a(jk+j)
                                    t0 = ajc + ajk
                                    t2 = ajc - ajk
                                    ajl = a(jl+j)
                                    t1 = ajg + ajl
                                    t3 = ss * ( ajg - ajl )
                                    bji = b(ji+j)
                                    bjk = b(jk+j)
                                    u0 = bji + bjk
                                    u2 = bji - bjk
                                    ajo = a(jo+j)
                                    ajl =  ajo
                                    bjl = b(jl+j)
                                    u1 = bjg + bjl
                                    u3 = ss * ( bjg - bjl )
                                    b(ji+j) = bjc
                                    a(jc+j) = t0 + t1
                                    b(jc+j) = u0 + u1
                                    bjo = b(jo+j)
                                    bjl =  bjo
                                    a(jg+j) = co1*(t2-u3) - si1*(u2+t3)
                                    b(jg+j) = si1*(t2-u3) + co1*(u2+t3)
                                    a(jk+j) = co2*(t0-t1) - si2*(u0-u1)
                                    b(jk+j) = si2*(t0-t1) + co2*(u0-u1)
                                    a(jo+j) = co3*(t2+u3) - si3*(u2-t3)
                                    b(jo+j) = si3*(t2+u3) + co3*(u2-t3)
                                    !----------------------------------------
                                    ajm = a(jm+j)
                                    t0 = ajm + ajl
                                    t2 = ajm - ajl
                                    ajn = a(jn+j)
                                    ajp = a(jp+j)
                                    t1 = ajn + ajp
                                    t3 = ss * ( ajn - ajp )
                                    a(jm+j) = ajd
                                    u0 = bjd + bjl
                                    u2 = bjd - bjl
                                    a(jn+j) = ajh
                                    bjn = b(jn+j)
                                    bjp = b(jp+j)
                                    u1 = bjn + bjp
                                    u3 = ss * ( bjn - bjp )
                                    b(jn+j) = bjh
                                    a(jd+j) = t0 + t1
                                    b(jd+j) = u0 + u1
                                    a(jh+j) = co1*(t2-u3) - si1*(u2+t3)
                                    b(jh+j) = si1*(t2-u3) + co1*(u2+t3)
                                    a(jl+j) = co2*(t0-t1) - si2*(u0-u1)
                                    b(jl+j) = si2*(t0-t1) + co2*(u0-u1)
                                    a(jp+j) = co3*(t2+u3) - si3*(u2-t3)
                                    b(jp+j) = si3*(t2+u3) + co3*(u2-t3)
                                    j = j + jump
                                enddo
                                !-----(end of loop across transforms)
                                ja = ja + jstepx
                                if (ja.lt.istart) ja = ja + ninc
                            enddo
                        enddo
                    enddo
                    !-----( end of double loop for this k )
                    kk = kk + 2*la
                enddo
                !-----( end of loop over values of k )
                la = 4*la
            enddo
            !-----( end of loop on type II radix-4 passes )
            !-----( nvex transforms completed)
490         continue
            istart = istart + nvex * jump
        enddo
        !-----( end of loop on blocks of transforms )
        !
        return
    end subroutine gpfa2f




    !     fortran version of *gpfa3* -
    !     radix-3 section of self-sorting, in-place
    !        generalized PFA
    !
    !-------------------------------------------------------------------
    !
    subroutine gpfa3f(a,b,trigs,inc,jump,n,mm,lot,isign)
        integer :: inc,jump,n,lot,mm,isign
        real(fpp) ::  a(*), b(*), trigs(*)
        real(fpp), parameter ::  sin60 = 0.866025403784437d0

        real(fpp) ::  s
        real(fpp) :: aja,ajc,t2,ajb,ajd,t1,t3
        real(fpp) :: bja,bjc,u2,bjb,bjd,u1,u3
        real(fpp) :: co1,si1,co2,si2
        real(fpp) :: c1
        real(fpp) :: aje,ajg,ajf,ajh
        real(fpp) :: bje,bjg,bjf,bjh
        real(fpp) :: aji
        real(fpp) :: bji

        integer :: ink,inq,ipass,istart,j,jjj,jstep,jstepl,jstepx
        integer :: ja,jb,jc,jd,je,jf,jg,ji,jh
        integer :: k,kk,l,la,laincl,left,ll,m,mh,mu,nb,nblox,ninc,nu,nvex
        integer :: n3
        n3 = 3**mm
        inq = n/n3
        jstepx = (n3-n) * inc
        ninc = n * inc
        ink = inc * inq
        mu = mod(inq,3)
        if (isign.eq.-1) mu = 3-mu
        m = mm
        mh = (m+1)/2
        s = float(isign)
        c1 = sin60
        if (mu.eq.2) c1 = -c1
        !
        nblox = 1 + (lot-1)/lvr
        left = lot
        s = float(isign)
        istart = 1
        !
        !  loop on blocks of lvr transforms
        !  --------------------------------
        do nb = 1 , nblox
            !
            if (left.le.lvr) then
                nvex = left
            else if (left.lt.(2*lvr)) then
                nvex = left/2
                nvex = nvex + mod(nvex,2)
            else
                nvex = lvr
            endif
            left = left - nvex
            !
            la = 1
            !
            !  loop on type I radix-3 passes
            !  -----------------------------
            do ipass = 1 , mh
                jstep = (n*inc) / (3*la)
                jstepl = jstep - ninc
                !
                !  k = 0 loop (no twiddle factors)
                !  -------------------------------
                do jjj = 0 , (n-1)*inc , 3*jstep
                    ja = istart + jjj
                    !
                    !  "transverse" loop
                    !  -----------------
                    do nu = 1 , inq
                        jb = ja + jstepl
                        if (jb.lt.istart) jb = jb + ninc
                        jc = jb + jstepl
                        if (jc.lt.istart) jc = jc + ninc
                        j = 0
                        !
                        !  loop across transforms
                        !  ----------------------
                        !cdiri$ ivdep, shortloop
                        do l = 1 , nvex
                            ajb = a(jb+j)
                            ajc = a(jc+j)
                            t1 = ajb + ajc
                            aja = a(ja+j)
                            t2 = aja - 0.5 * t1
                            t3 = c1 * ( ajb - ajc )
                            bjb = b(jb+j)
                            bjc = b(jc+j)
                            u1 = bjb + bjc
                            bja = b(ja+j)
                            u2 = bja - 0.5 * u1
                            u3 = c1 * ( bjb - bjc )
                            a(ja+j) = aja + t1
                            b(ja+j) = bja + u1
                            a(jb+j) = t2 - u3
                            b(jb+j) = u2 + t3
                            a(jc+j) = t2 + u3
                            b(jc+j) = u2 - t3
                            j = j + jump
                        enddo
                        ja = ja + jstepx
                        if (ja.lt.istart) ja = ja + ninc
                    enddo
                enddo
                !
                !  finished if n3 = 3
                !  ------------------
                if (n3.eq.3) go to 490
                kk = 2 * la
                !
                !  loop on nonzero k
                !  -----------------
                do k = ink , jstep-ink , ink
                    co1 = trigs(kk+1)
                    si1 = s*trigs(kk+2)
                    co2 = trigs(2*kk+1)
                    si2 = s*trigs(2*kk+2)
                    !
                    !  loop along transform
                    !  --------------------
                    do jjj = k , (n-1)*inc , 3*jstep
                        ja = istart + jjj
                        !
                        !  "transverse" loop
                        !  -----------------
                        do nu = 1 , inq
                            jb = ja + jstepl
                            if (jb.lt.istart) jb = jb + ninc
                            jc = jb + jstepl
                            if (jc.lt.istart) jc = jc + ninc
                            j = 0
                            !
                            !  loop across transforms
                            !  ----------------------
                            !cdiri$ ivdep,shortloop
                            do l = 1 , nvex
                                ajb = a(jb+j)
                                ajc = a(jc+j)
                                t1 = ajb + ajc
                                aja = a(ja+j)
                                t2 = aja - 0.5 * t1
                                t3 = c1 * ( ajb - ajc )
                                bjb = b(jb+j)
                                bjc = b(jc+j)
                                u1 = bjb + bjc
                                bja = b(ja+j)
                                u2 = bja - 0.5 * u1
                                u3 = c1 * ( bjb - bjc )
                                a(ja+j) = aja + t1
                                b(ja+j) = bja + u1
                                a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
                                b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
                                a(jc+j) = co2*(t2+u3) - si2*(u2-t3)
                                b(jc+j) = si2*(t2+u3) + co2*(u2-t3)
                                j = j + jump
                            enddo
                            !-----( end of loop across transforms )
                            ja = ja + jstepx
                            if (ja.lt.istart) ja = ja + ninc
                        enddo
                    enddo
                    !-----( end of loop along transforms )
                    kk = kk + 2*la
                enddo
                !-----( end of loop on nonzero k )
                la = 3*la
            enddo
            !-----( end of loop on type I radix-3 passes)
            !
            !  loop on type II radix-3 passes
            !  ------------------------------
400         continue
            !
            do ipass = mh+1 , m
                jstep = (n*inc) / (3*la)
                jstepl = jstep - ninc
                laincl = la*ink - ninc
                !
                !  k=0 loop (no twiddle factors)
                !  -----------------------------
                do ll = 0 , (la-1)*ink , 3*jstep
                    !
                    do jjj = ll , (n-1)*inc , 3*la*ink
                        ja = istart + jjj
                        !
                        !  "transverse" loop
                        !  -----------------
                        do nu = 1 , inq
                            jb = ja + jstepl
                            if (jb.lt.istart) jb = jb + ninc
                            jc = jb + jstepl
                            if (jc.lt.istart) jc = jc + ninc
                            jd = ja + laincl
                            if (jd.lt.istart) jd = jd + ninc
                            je = jd + jstepl
                            if (je.lt.istart) je = je + ninc
                            jf = je + jstepl
                            if (jf.lt.istart) jf = jf + ninc
                            jg = jd + laincl
                            if (jg.lt.istart) jg = jg + ninc
                            jh = jg + jstepl
                            if (jh.lt.istart) jh = jh + ninc
                            ji = jh + jstepl
                            if (ji.lt.istart) ji = ji + ninc
                            j = 0
                            !
                            !  loop across transforms
                            !  ----------------------
                            !cdiri$ ivdep, shortloop
                            do l = 1 , nvex
                                ajb = a(jb+j)
                                ajc = a(jc+j)
                                t1 = ajb + ajc
                                aja = a(ja+j)
                                t2 = aja - 0.5 * t1
                                t3 = c1 * ( ajb - ajc )
                                ajd = a(jd+j)
                                ajb =  ajd
                                bjb = b(jb+j)
                                bjc = b(jc+j)
                                u1 = bjb + bjc
                                bja = b(ja+j)
                                u2 = bja - 0.5 * u1
                                u3 = c1 * ( bjb - bjc )
                                bjd = b(jd+j)
                                bjb =  bjd
                                a(ja+j) = aja + t1
                                b(ja+j) = bja + u1
                                a(jd+j) = t2 - u3
                                b(jd+j) = u2 + t3
                                ajc =  t2 + u3
                                bjc =  u2 - t3
                                !----------------------
                                aje = a(je+j)
                                ajf = a(jf+j)
                                t1 = aje + ajf
                                t2 = ajb - 0.5 * t1
                                t3 = c1 * ( aje - ajf )
                                ajh = a(jh+j)
                                ajf =  ajh
                                bje = b(je+j)
                                bjf = b(jf+j)
                                u1 = bje + bjf
                                u2 = bjb - 0.5 * u1
                                u3 = c1 * ( bje - bjf )
                                bjh = b(jh+j)
                                bjf =  bjh
                                a(jb+j) = ajb + t1
                                b(jb+j) = bjb + u1
                                a(je+j) = t2 - u3
                                b(je+j) = u2 + t3
                                a(jh+j) = t2 + u3
                                b(jh+j) = u2 - t3
                                !----------------------
                                aji = a(ji+j)
                                t1 = ajf + aji
                                ajg = a(jg+j)
                                t2 = ajg - 0.5 * t1
                                t3 = c1 * ( ajf - aji )
                                t1 = ajg + t1
                                a(jg+j) = ajc
                                bji = b(ji+j)
                                u1 = bjf + bji
                                bjg = b(jg+j)
                                u2 = bjg - 0.5 * u1
                                u3 = c1 * ( bjf - bji )
                                u1 = bjg + u1
                                b(jg+j) = bjc
                                a(jc+j) = t1
                                b(jc+j) = u1
                                a(jf+j) = t2 - u3
                                b(jf+j) = u2 + t3
                                a(ji+j) = t2 + u3
                                b(ji+j) = u2 - t3
                                j = j + jump
                            enddo
                            !-----( end of loop across transforms )
                            ja = ja + jstepx
                            if (ja.lt.istart) ja = ja + ninc
                        enddo
                    enddo
                enddo
                !-----( end of double loop for k=0 )
                !
                !  finished if last pass
                !  ---------------------
                if (ipass.eq.m) go to 490
                !
                kk = 2*la
                !
                !     loop on nonzero k
                !     -----------------
                do k = ink , jstep-ink , ink
                    co1 = trigs(kk+1)
                    si1 = s*trigs(kk+2)
                    co2 = trigs(2*kk+1)
                    si2 = s*trigs(2*kk+2)
                    !
                    !  double loop along first transform in block
                    !  ------------------------------------------
                    do ll = k , (la-1)*ink , 3*jstep
                        !
                        do jjj = ll , (n-1)*inc , 3*la*ink
                            ja = istart + jjj
                            !
                            !  "transverse" loop
                            !  -----------------
                            do nu = 1 , inq
                                jb = ja + jstepl
                                if (jb.lt.istart) jb = jb + ninc
                                jc = jb + jstepl
                                if (jc.lt.istart) jc = jc + ninc
                                jd = ja + laincl
                                if (jd.lt.istart) jd = jd + ninc
                                je = jd + jstepl
                                if (je.lt.istart) je = je + ninc
                                jf = je + jstepl
                                if (jf.lt.istart) jf = jf + ninc
                                jg = jd + laincl
                                if (jg.lt.istart) jg = jg + ninc
                                jh = jg + jstepl
                                if (jh.lt.istart) jh = jh + ninc
                                ji = jh + jstepl
                                if (ji.lt.istart) ji = ji + ninc
                                j = 0
                                !
                                !  loop across transforms
                                !  ----------------------
                                !cdiri$ ivdep, shortloop
                                do l = 1 , nvex
                                    ajb = a(jb+j)
                                    ajc = a(jc+j)
                                    t1 = ajb + ajc
                                    aja = a(ja+j)
                                    t2 = aja - 0.5 * t1
                                    t3 = c1 * ( ajb - ajc )
                                    ajd = a(jd+j)
                                    ajb =  ajd
                                    bjb = b(jb+j)
                                    bjc = b(jc+j)
                                    u1 = bjb + bjc
                                    bja = b(ja+j)
                                    u2 = bja - 0.5 * u1
                                    u3 = c1 * ( bjb - bjc )
                                    bjd = b(jd+j)
                                    bjb =  bjd
                                    a(ja+j) = aja + t1
                                    b(ja+j) = bja + u1
                                    a(jd+j) = co1*(t2-u3) - si1*(u2+t3)
                                    b(jd+j) = si1*(t2-u3) + co1*(u2+t3)
                                    ajc =  co2*(t2+u3) - si2*(u2-t3)
                                    bjc =  si2*(t2+u3) + co2*(u2-t3)
                                    !----------------------
                                    aje = a(je+j)
                                    ajf = a(jf+j)
                                    t1 = aje + ajf
                                    t2 = ajb - 0.5 * t1
                                    t3 = c1 * ( aje - ajf )
                                    ajh = a(jh+j)
                                    ajf =  ajh
                                    bje = b(je+j)
                                    bjf = b(jf+j)
                                    u1 = bje + bjf
                                    u2 = bjb - 0.5 * u1
                                    u3 = c1 * ( bje - bjf )
                                    bjh = b(jh+j)
                                    bjf =  bjh
                                    a(jb+j) = ajb + t1
                                    b(jb+j) = bjb + u1
                                    a(je+j) = co1*(t2-u3) - si1*(u2+t3)
                                    b(je+j) = si1*(t2-u3) + co1*(u2+t3)
                                    a(jh+j) = co2*(t2+u3) - si2*(u2-t3)
                                    b(jh+j) = si2*(t2+u3) + co2*(u2-t3)
                                    !----------------------
                                    aji = a(ji+j)
                                    t1 = ajf + aji
                                    ajg = a(jg+j)
                                    t2 = ajg - 0.5 * t1
                                    t3 = c1 * ( ajf - aji )
                                    t1 = ajg + t1
                                    a(jg+j) = ajc
                                    bji = b(ji+j)
                                    u1 = bjf + bji
                                    bjg = b(jg+j)
                                    u2 = bjg - 0.5 * u1
                                    u3 = c1 * ( bjf - bji )
                                    u1 = bjg + u1
                                    b(jg+j) = bjc
                                    a(jc+j) = t1
                                    b(jc+j) = u1
                                    a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
                                    b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
                                    a(ji+j) = co2*(t2+u3) - si2*(u2-t3)
                                    b(ji+j) = si2*(t2+u3) + co2*(u2-t3)
                                    j = j + jump
                                enddo
                                !-----(end of loop across transforms)
                                ja = ja + jstepx
                                if (ja.lt.istart) ja = ja + ninc
                            enddo
                        enddo
                    enddo
                    !-----( end of double loop for this k )
                    kk = kk + 2*la
                enddo
                !-----( end of loop over values of k )
                la = 3*la
            enddo
            !-----( end of loop on type II radix-3 passes )
            !-----( nvex transforms completed)
490         continue
            istart = istart + nvex * jump
        enddo
        !-----( end of loop on blocks of transforms )
        !
        return
    end subroutine gpfa3f




    !     fortran version of *gpfa5* -
    !     radix-5 section of self-sorting, in-place,
    !        generalized pfa
    !
    !-------------------------------------------------------------------
    !
    subroutine gpfa5f(a,b,trigs,inc,jump,n,mm,lot,isign)
        integer :: inc,jump,n,lot,mm,isign
        real(fpp) ::  a(*), b(*), trigs(*)
        real(fpp) :: sin36,sin72, qrt5
        data sin36/0.587785252292473/, sin72/0.951056516295154/, qrt5/0.559016994374947/

        real(fpp) ::  s
        real(fpp) :: t1,t2,t3,t4,t5,t6,t8,t9,t10,t11
        real(fpp) :: u1,u2,u3,u4,u5,u6,u8,u9,u10,u11
        real(fpp) :: aja,ajc,ajb,ajd
        real(fpp) :: bja,bjc,bjb,bjd
        real(fpp) :: co1,si1,co2,si2,co3,si3
        real(fpp) :: c1,c2,c3
        real(fpp) :: aje,ajg,ajf,ajh
        real(fpp) :: bje,bjg,bjf,bjh
        real(fpp) :: co4,si4
        real(fpp) :: aji,ajk,ajl,ajm,ajn,ajo,ajp
        real(fpp) :: bji,bjk,bjl,bjm,bjn,bjo,bjp
        real(fpp) :: ajj,ajq,ajr,ajs,ajt,aju,ajv,ajw,ajx,ajy
        real(fpp) :: ax,bx,t7,u7
        real(fpp) :: bjj,bjq,bjr,bjs,bjt,bju,bjv,bjw,bjx,bjy
        integer :: ink,inq,ipass,istart,j,jjj,jstep,jstepl,jstepx
        integer :: ja,jb,jc,jd,je,jf,jg,ji,jj,jh,jk,jl,jm,jn,jo,jp,jq,jr,js,jt,ju,jv,jw,jx,jy
        integer :: k,kk,l,la,laincl,left,ll,m,mh,mu,n5,nb,nblox,ninc,nu,nvex

        n5 = 5 ** mm
        inq = n / n5
        jstepx = (n5-n) * inc
        ninc = n * inc
        ink = inc * inq
        mu = mod(inq,5)
        if (isign.eq.-1) mu = 5 - mu
        !
        m = mm
        mh = (m+1)/2
        s = float(isign)
        c1 = qrt5
        c2 = sin72
        c3 = sin36
        if (mu.eq.2.or.mu.eq.3) then
            c1 = -c1
            c2 = sin36
            c3 = sin72
        endif
        if (mu.eq.3.or.mu.eq.4) c2 = -c2
        if (mu.eq.2.or.mu.eq.4) c3 = -c3
        !
        nblox = 1 + (lot-1)/lvr
        left = lot
        s = float(isign)
        istart = 1
        !
        !  loop on blocks of lvr transforms
        !  --------------------------------
        do nb = 1 , nblox
            !
            if (left.le.lvr) then
                nvex = left
            else if (left.lt.(2*lvr)) then
                nvex = left/2
                nvex = nvex + mod(nvex,2)
            else
                nvex = lvr
            endif
            left = left - nvex
            !
            la = 1
            !
            !  loop on type I radix-5 passes
            !  -----------------------------
            do ipass = 1 , mh
                jstep = (n*inc) / (5*la)
                jstepl = jstep - ninc
                kk = 0
                !
                !  loop on k
                !  ---------
                do k = 0 , jstep-ink , ink
                    !
                    if (k.gt.0) then
                        co1 = trigs(kk+1)
                        si1 = s*trigs(kk+2)
                        co2 = trigs(2*kk+1)
                        si2 = s*trigs(2*kk+2)
                        co3 = trigs(3*kk+1)
                        si3 = s*trigs(3*kk+2)
                        co4 = trigs(4*kk+1)
                        si4 = s*trigs(4*kk+2)
                    endif
                    !
                    !  loop along transform
                    !  --------------------
                    do jjj = k , (n-1)*inc , 5*jstep
                        ja = istart + jjj
                        !
                        !     "transverse" loop
                        !     -----------------
                        do nu = 1 , inq
                            jb = ja + jstepl
                            if (jb.lt.istart) jb = jb + ninc
                            jc = jb + jstepl
                            if (jc.lt.istart) jc = jc + ninc
                            jd = jc + jstepl
                            if (jd.lt.istart) jd = jd + ninc
                            je = jd + jstepl
                            if (je.lt.istart) je = je + ninc
                            j = 0
                            !
                            !  loop across transforms
                            !  ----------------------
                            if (k.eq.0) then
                                !
                                !cdiri$ ivdep, shortloop
                                do l = 1 , nvex
                                    ajb = a(jb+j)
                                    aje = a(je+j)
                                    t1 = ajb + aje
                                    ajc = a(jc+j)
                                    ajd = a(jd+j)
                                    t2 = ajc + ajd
                                    t3 = ajb - aje
                                    t4 = ajc - ajd
                                    t5 = t1 + t2
                                    t6 = c1 * ( t1 - t2 )
                                    aja = a(ja+j)
                                    t7 = aja - 0.25 * t5
                                    a(ja+j) = aja + t5
                                    t8 = t7 + t6
                                    t9 = t7 - t6
                                    t10 = c3 * t3 - c2 * t4
                                    t11 = c2 * t3 + c3 * t4
                                    bjb = b(jb+j)
                                    bje = b(je+j)
                                    u1 = bjb + bje
                                    bjc = b(jc+j)
                                    bjd = b(jd+j)
                                    u2 = bjc + bjd
                                    u3 = bjb - bje
                                    u4 = bjc - bjd
                                    u5 = u1 + u2
                                    u6 = c1 * ( u1 - u2 )
                                    bja = b(ja+j)
                                    u7 = bja - 0.25 * u5
                                    b(ja+j) = bja + u5
                                    u8 = u7 + u6
                                    u9 = u7 - u6
                                    u10 = c3 * u3 - c2 * u4
                                    u11 = c2 * u3 + c3 * u4
                                    a(jb+j) = t8 - u11
                                    b(jb+j) = u8 + t11
                                    a(je+j) = t8 + u11
                                    b(je+j) = u8 - t11
                                    a(jc+j) = t9 - u10
                                    b(jc+j) = u9 + t10
                                    a(jd+j) = t9 + u10
                                    b(jd+j) = u9 - t10
                                    j = j + jump
                                enddo
                                !
                            else
                                !
                                !cdiri$ ivdep,shortloop
                                do l = 1 , nvex
                                    ajb = a(jb+j)
                                    aje = a(je+j)
                                    t1 = ajb + aje
                                    ajc = a(jc+j)
                                    ajd = a(jd+j)
                                    t2 = ajc + ajd
                                    t3 = ajb - aje
                                    t4 = ajc - ajd
                                    t5 = t1 + t2
                                    t6 = c1 * ( t1 - t2 )
                                    aja = a(ja+j)
                                    t7 = aja - 0.25 * t5
                                    a(ja+j) = aja + t5
                                    t8 = t7 + t6
                                    t9 = t7 - t6
                                    t10 = c3 * t3 - c2 * t4
                                    t11 = c2 * t3 + c3 * t4
                                    bjb = b(jb+j)
                                    bje = b(je+j)
                                    u1 = bjb + bje
                                    bjc = b(jc+j)
                                    bjd = b(jd+j)
                                    u2 = bjc + bjd
                                    u3 = bjb - bje
                                    u4 = bjc - bjd
                                    u5 = u1 + u2
                                    u6 = c1 * ( u1 - u2 )
                                    bja = b(ja+j)
                                    u7 = bja - 0.25 * u5
                                    b(ja+j) = bja + u5
                                    u8 = u7 + u6
                                    u9 = u7 - u6
                                    u10 = c3 * u3 - c2 * u4
                                    u11 = c2 * u3 + c3 * u4
                                    a(jb+j) = co1*(t8-u11) - si1*(u8+t11)
                                    b(jb+j) = si1*(t8-u11) + co1*(u8+t11)
                                    a(je+j) = co4*(t8+u11) - si4*(u8-t11)
                                    b(je+j) = si4*(t8+u11) + co4*(u8-t11)
                                    a(jc+j) = co2*(t9-u10) - si2*(u9+t10)
                                    b(jc+j) = si2*(t9-u10) + co2*(u9+t10)
                                    a(jd+j) = co3*(t9+u10) - si3*(u9-t10)
                                    b(jd+j) = si3*(t9+u10) + co3*(u9-t10)
                                    j = j + jump
                                enddo
                                !
                            endif
                            !
                            !-----( end of loop across transforms )
                            !
                            ja = ja + jstepx
                            if (ja.lt.istart) ja = ja + ninc
                        enddo
                    enddo
                    !-----( end of loop along transforms )
                    kk = kk + 2*la
                enddo
                !-----( end of loop on nonzero k )
                la = 5*la
            enddo
            !-----( end of loop on type I radix-5 passes)
            !
            if (n.eq.5) go to 490
            !
            !  loop on type II radix-5 passes
            !  ------------------------------
400         continue
            !
            do ipass = mh+1 , m
                jstep = (n*inc) / (5*la)
                jstepl = jstep - ninc
                laincl = la * ink - ninc
                kk = 0
                !
                !     loop on k
                !     ---------
                do k = 0 , jstep-ink , ink
                    !
                    if (k.gt.0) then
                        co1 = trigs(kk+1)
                        si1 = s*trigs(kk+2)
                        co2 = trigs(2*kk+1)
                        si2 = s*trigs(2*kk+2)
                        co3 = trigs(3*kk+1)
                        si3 = s*trigs(3*kk+2)
                        co4 = trigs(4*kk+1)
                        si4 = s*trigs(4*kk+2)
                    endif
                    !
                    !  double loop along first transform in block
                    !  ------------------------------------------
                    do ll = k , (la-1)*ink , 5*jstep
                        !
                        do jjj = ll , (n-1)*inc , 5*la*ink
                            ja = istart + jjj
                            !
                            !     "transverse" loop
                            !     -----------------
                            do nu = 1 , inq
                                jb = ja + jstepl
                                if (jb.lt.istart) jb = jb + ninc
                                jc = jb + jstepl
                                if (jc.lt.istart) jc = jc + ninc
                                jd = jc + jstepl
                                if (jd.lt.istart) jd = jd + ninc
                                je = jd + jstepl
                                if (je.lt.istart) je = je + ninc
                                jf = ja + laincl
                                if (jf.lt.istart) jf = jf + ninc
                                jg = jf + jstepl
                                if (jg.lt.istart) jg = jg + ninc
                                jh = jg + jstepl
                                if (jh.lt.istart) jh = jh + ninc
                                ji = jh + jstepl
                                if (ji.lt.istart) ji = ji + ninc
                                jj = ji + jstepl
                                if (jj.lt.istart) jj = jj + ninc
                                jk = jf + laincl
                                if (jk.lt.istart) jk = jk + ninc
                                jl = jk + jstepl
                                if (jl.lt.istart) jl = jl + ninc
                                jm = jl + jstepl
                                if (jm.lt.istart) jm = jm + ninc
                                jn = jm + jstepl
                                if (jn.lt.istart) jn = jn + ninc
                                jo = jn + jstepl
                                if (jo.lt.istart) jo = jo + ninc
                                jp = jk + laincl
                                if (jp.lt.istart) jp = jp + ninc
                                jq = jp + jstepl
                                if (jq.lt.istart) jq = jq + ninc
                                jr = jq + jstepl
                                if (jr.lt.istart) jr = jr + ninc
                                js = jr + jstepl
                                if (js.lt.istart) js = js + ninc
                                jt = js + jstepl
                                if (jt.lt.istart) jt = jt + ninc
                                ju = jp + laincl
                                if (ju.lt.istart) ju = ju + ninc
                                jv = ju + jstepl
                                if (jv.lt.istart) jv = jv + ninc
                                jw = jv + jstepl
                                if (jw.lt.istart) jw = jw + ninc
                                jx = jw + jstepl
                                if (jx.lt.istart) jx = jx + ninc
                                jy = jx + jstepl
                                if (jy.lt.istart) jy = jy + ninc
                                j = 0
                                !
                                !  loop across transforms
                                !  ----------------------
                                if (k.eq.0) then
                                    !
                                    !cdiri$ ivdep, shortloop
                                    do l = 1 , nvex
                                        ajb = a(jb+j)
                                        aje = a(je+j)
                                        t1 = ajb + aje
                                        ajc = a(jc+j)
                                        ajd = a(jd+j)
                                        t2 = ajc + ajd
                                        t3 = ajb - aje
                                        t4 = ajc - ajd
                                        ajf = a(jf+j)
                                        ajb =  ajf
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        aja = a(ja+j)
                                        t7 = aja - 0.25 * t5
                                        a(ja+j) = aja + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        ajk = a(jk+j)
                                        ajc =  ajk
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        bjb = b(jb+j)
                                        bje = b(je+j)
                                        u1 = bjb + bje
                                        bjc = b(jc+j)
                                        bjd = b(jd+j)
                                        u2 = bjc + bjd
                                        u3 = bjb - bje
                                        u4 = bjc - bjd
                                        bjf = b(jf+j)
                                        bjb =  bjf
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        bja = b(ja+j)
                                        u7 = bja - 0.25 * u5
                                        b(ja+j) = bja + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        bjk = b(jk+j)
                                        bjc =  bjk
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        a(jf+j) = t8 - u11
                                        b(jf+j) = u8 + t11
                                        aje =  t8 + u11
                                        bje =  u8 - t11
                                        a(jk+j) = t9 - u10
                                        b(jk+j) = u9 + t10
                                        ajd =  t9 + u10
                                        bjd =  u9 - t10
                                        !----------------------
                                        ajg = a(jg+j)
                                        ajj = a(jj+j)
                                        t1 = ajg + ajj
                                        ajh = a(jh+j)
                                        aji = a(ji+j)
                                        t2 = ajh + aji
                                        t3 = ajg - ajj
                                        t4 = ajh - aji
                                        ajl = a(jl+j)
                                        ajh =  ajl
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        t7 = ajb - 0.25 * t5
                                        a(jb+j) = ajb + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        ajq = a(jq+j)
                                        aji =  ajq
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        bjg = b(jg+j)
                                        bjj = b(jj+j)
                                        u1 = bjg + bjj
                                        bjh = b(jh+j)
                                        bji = b(ji+j)
                                        u2 = bjh + bji
                                        u3 = bjg - bjj
                                        u4 = bjh - bji
                                        bjl = b(jl+j)
                                        bjh =  bjl
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        u7 = bjb - 0.25 * u5
                                        b(jb+j) = bjb + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        bjq = b(jq+j)
                                        bji =  bjq
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        a(jg+j) = t8 - u11
                                        b(jg+j) = u8 + t11
                                        ajj =  t8 + u11
                                        bjj =  u8 - t11
                                        a(jl+j) = t9 - u10
                                        b(jl+j) = u9 + t10
                                        a(jq+j) = t9 + u10
                                        b(jq+j) = u9 - t10
                                        !----------------------
                                        ajo = a(jo+j)
                                        t1 = ajh + ajo
                                        ajm = a(jm+j)
                                        ajn = a(jn+j)
                                        t2 = ajm + ajn
                                        t3 = ajh - ajo
                                        t4 = ajm - ajn
                                        ajr = a(jr+j)
                                        ajn =  ajr
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        t7 = ajc - 0.25 * t5
                                        a(jc+j) = ajc + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        ajw = a(jw+j)
                                        ajo =  ajw
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        bjo = b(jo+j)
                                        u1 = bjh + bjo
                                        bjm = b(jm+j)
                                        bjn = b(jn+j)
                                        u2 = bjm + bjn
                                        u3 = bjh - bjo
                                        u4 = bjm - bjn
                                        bjr = b(jr+j)
                                        bjn =  bjr
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        u7 = bjc - 0.25 * u5
                                        b(jc+j) = bjc + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        bjw = b(jw+j)
                                        bjo =  bjw
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        a(jh+j) = t8 - u11
                                        b(jh+j) = u8 + t11
                                        a(jw+j) = t8 + u11
                                        b(jw+j) = u8 - t11
                                        a(jm+j) = t9 - u10
                                        b(jm+j) = u9 + t10
                                        a(jr+j) = t9 + u10
                                        b(jr+j) = u9 - t10
                                        !----------------------
                                        ajt = a(jt+j)
                                        t1 = aji + ajt
                                        ajs = a(js+j)
                                        t2 = ajn + ajs
                                        t3 = aji - ajt
                                        t4 = ajn - ajs
                                        ajx = a(jx+j)
                                        ajt =  ajx
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        ajp = a(jp+j)
                                        t7 = ajp - 0.25 * t5
                                        ax = ajp + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        a(jp+j) = ajd
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        a(jd+j) = ax
                                        bjt = b(jt+j)
                                        u1 = bji + bjt
                                        bjs = b(js+j)
                                        u2 = bjn + bjs
                                        u3 = bji - bjt
                                        u4 = bjn - bjs
                                        bjx = b(jx+j)
                                        bjt =  bjx
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        bjp = b(jp+j)
                                        u7 = bjp - 0.25 * u5
                                        bx = bjp + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        b(jp+j) = bjd
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        b(jd+j) = bx
                                        a(ji+j) = t8 - u11
                                        b(ji+j) = u8 + t11
                                        a(jx+j) = t8 + u11
                                        b(jx+j) = u8 - t11
                                        a(jn+j) = t9 - u10
                                        b(jn+j) = u9 + t10
                                        a(js+j) = t9 + u10
                                        b(js+j) = u9 - t10
                                        !----------------------
                                        ajv = a(jv+j)
                                        ajy = a(jy+j)
                                        t1 = ajv + ajy
                                        t2 = ajo + ajt
                                        t3 = ajv - ajy
                                        t4 = ajo - ajt
                                        a(jv+j) = ajj
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        aju = a(ju+j)
                                        t7 = aju - 0.25 * t5
                                        ax = aju + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        a(ju+j) = aje
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        a(je+j) = ax
                                        bjv = b(jv+j)
                                        bjy = b(jy+j)
                                        u1 = bjv + bjy
                                        u2 = bjo + bjt
                                        u3 = bjv - bjy
                                        u4 = bjo - bjt
                                        b(jv+j) = bjj
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        bju = b(ju+j)
                                        u7 = bju - 0.25 * u5
                                        bx = bju + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        b(ju+j) = bje
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        b(je+j) = bx
                                        a(jj+j) = t8 - u11
                                        b(jj+j) = u8 + t11
                                        a(jy+j) = t8 + u11
                                        b(jy+j) = u8 - t11
                                        a(jo+j) = t9 - u10
                                        b(jo+j) = u9 + t10
                                        a(jt+j) = t9 + u10
                                        b(jt+j) = u9 - t10
                                        j = j + jump
                                    enddo
                                    !
                                else
                                    !
                                    !cdiri$ ivdep, shortloop
                                    do l = 1 , nvex
                                        ajb = a(jb+j)
                                        aje = a(je+j)
                                        t1 = ajb + aje
                                        ajc = a(jc+j)
                                        ajd = a(jd+j)
                                        t2 = ajc + ajd
                                        t3 = ajb - aje
                                        t4 = ajc - ajd
                                        ajf = a(jf+j)
                                        ajb =  ajf
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        aja = a(ja+j)
                                        t7 = aja - 0.25 * t5
                                        a(ja+j) = aja + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        ajk = a(jk+j)
                                        ajc =  ajk
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        bjb = b(jb+j)
                                        bje = b(je+j)
                                        u1 = bjb + bje
                                        bjc = b(jc+j)
                                        bjd = b(jd+j)
                                        u2 = bjc + bjd
                                        u3 = bjb - bje
                                        u4 = bjc - bjd
                                        bjf = b(jf+j)
                                        bjb =  bjf
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        bja = b(ja+j)
                                        u7 = bja - 0.25 * u5
                                        b(ja+j) = bja + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        bjk = b(jk+j)
                                        bjc =  bjk
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        a(jf+j) = co1*(t8-u11) - si1*(u8+t11)
                                        b(jf+j) = si1*(t8-u11) + co1*(u8+t11)
                                        aje =  co4*(t8+u11) - si4*(u8-t11)
                                        bje =  si4*(t8+u11) + co4*(u8-t11)
                                        a(jk+j) = co2*(t9-u10) - si2*(u9+t10)
                                        b(jk+j) = si2*(t9-u10) + co2*(u9+t10)
                                        ajd =  co3*(t9+u10) - si3*(u9-t10)
                                        bjd =  si3*(t9+u10) + co3*(u9-t10)
                                        !----------------------
                                        ajg = a(jg+j)
                                        ajj = a(jj+j)
                                        t1 = ajg + ajj
                                        ajh = a(jh+j)
                                        aji = a(ji+j)
                                        t2 = ajh + aji
                                        t3 = ajg - ajj
                                        t4 = ajh - aji
                                        ajl = a(jl+j)
                                        ajh =  ajl
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        t7 = ajb - 0.25 * t5
                                        a(jb+j) = ajb + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        ajq = a(jq+j)
                                        aji =  ajq
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        bjg = b(jg+j)
                                        bjj = b(jj+j)
                                        u1 = bjg + bjj
                                        bjh = b(jh+j)
                                        bji = b(ji+j)
                                        u2 = bjh + bji
                                        u3 = bjg - bjj
                                        u4 = bjh - bji
                                        bjl = b(jl+j)
                                        bjh =  bjl
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        u7 = bjb - 0.25 * u5
                                        b(jb+j) = bjb + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        bjq = b(jq+j)
                                        bji =  bjq
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        a(jg+j) = co1*(t8-u11) - si1*(u8+t11)
                                        b(jg+j) = si1*(t8-u11) + co1*(u8+t11)
                                        ajj =  co4*(t8+u11) - si4*(u8-t11)
                                        bjj =  si4*(t8+u11) + co4*(u8-t11)
                                        a(jl+j) = co2*(t9-u10) - si2*(u9+t10)
                                        b(jl+j) = si2*(t9-u10) + co2*(u9+t10)
                                        a(jq+j) = co3*(t9+u10) - si3*(u9-t10)
                                        b(jq+j) = si3*(t9+u10) + co3*(u9-t10)
                                        !----------------------
                                        ajo = a(jo+j)
                                        t1 = ajh + ajo
                                        ajm = a(jm+j)
                                        ajn = a(jn+j)
                                        t2 = ajm + ajn
                                        t3 = ajh - ajo
                                        t4 = ajm - ajn
                                        ajr = a(jr+j)
                                        ajn =  ajr
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        t7 = ajc - 0.25 * t5
                                        a(jc+j) = ajc + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        ajw = a(jw+j)
                                        ajo =  ajw
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        bjo = b(jo+j)
                                        u1 = bjh + bjo
                                        bjm = b(jm+j)
                                        bjn = b(jn+j)
                                        u2 = bjm + bjn
                                        u3 = bjh - bjo
                                        u4 = bjm - bjn
                                        bjr = b(jr+j)
                                        bjn =  bjr
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        u7 = bjc - 0.25 * u5
                                        b(jc+j) = bjc + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        bjw = b(jw+j)
                                        bjo =  bjw
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        a(jh+j) = co1*(t8-u11) - si1*(u8+t11)
                                        b(jh+j) = si1*(t8-u11) + co1*(u8+t11)
                                        a(jw+j) = co4*(t8+u11) - si4*(u8-t11)
                                        b(jw+j) = si4*(t8+u11) + co4*(u8-t11)
                                        a(jm+j) = co2*(t9-u10) - si2*(u9+t10)
                                        b(jm+j) = si2*(t9-u10) + co2*(u9+t10)
                                        a(jr+j) = co3*(t9+u10) - si3*(u9-t10)
                                        b(jr+j) = si3*(t9+u10) + co3*(u9-t10)
                                        !----------------------
                                        ajt = a(jt+j)
                                        t1 = aji + ajt
                                        ajs = a(js+j)
                                        t2 = ajn + ajs
                                        t3 = aji - ajt
                                        t4 = ajn - ajs
                                        ajx = a(jx+j)
                                        ajt =  ajx
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        ajp = a(jp+j)
                                        t7 = ajp - 0.25 * t5
                                        ax = ajp + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        a(jp+j) = ajd
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        a(jd+j) = ax
                                        bjt = b(jt+j)
                                        u1 = bji + bjt
                                        bjs = b(js+j)
                                        u2 = bjn + bjs
                                        u3 = bji - bjt
                                        u4 = bjn - bjs
                                        bjx = b(jx+j)
                                        bjt =  bjx
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        bjp = b(jp+j)
                                        u7 = bjp - 0.25 * u5
                                        bx = bjp + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        b(jp+j) = bjd
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        b(jd+j) = bx
                                        a(ji+j) = co1*(t8-u11) - si1*(u8+t11)
                                        b(ji+j) = si1*(t8-u11) + co1*(u8+t11)
                                        a(jx+j) = co4*(t8+u11) - si4*(u8-t11)
                                        b(jx+j) = si4*(t8+u11) + co4*(u8-t11)
                                        a(jn+j) = co2*(t9-u10) - si2*(u9+t10)
                                        b(jn+j) = si2*(t9-u10) + co2*(u9+t10)
                                        a(js+j) = co3*(t9+u10) - si3*(u9-t10)
                                        b(js+j) = si3*(t9+u10) + co3*(u9-t10)
                                        !----------------------
                                        ajv = a(jv+j)
                                        ajy = a(jy+j)
                                        t1 = ajv + ajy
                                        t2 = ajo + ajt
                                        t3 = ajv - ajy
                                        t4 = ajo - ajt
                                        a(jv+j) = ajj
                                        t5 = t1 + t2
                                        t6 = c1 * ( t1 - t2 )
                                        aju = a(ju+j)
                                        t7 = aju - 0.25 * t5
                                        ax = aju + t5
                                        t8 = t7 + t6
                                        t9 = t7 - t6
                                        a(ju+j) = aje
                                        t10 = c3 * t3 - c2 * t4
                                        t11 = c2 * t3 + c3 * t4
                                        a(je+j) = ax
                                        bjv = b(jv+j)
                                        bjy = b(jy+j)
                                        u1 = bjv + bjy
                                        u2 = bjo + bjt
                                        u3 = bjv - bjy
                                        u4 = bjo - bjt
                                        b(jv+j) = bjj
                                        u5 = u1 + u2
                                        u6 = c1 * ( u1 - u2 )
                                        bju = b(ju+j)
                                        u7 = bju - 0.25 * u5
                                        bx = bju + u5
                                        u8 = u7 + u6
                                        u9 = u7 - u6
                                        b(ju+j) = bje
                                        u10 = c3 * u3 - c2 * u4
                                        u11 = c2 * u3 + c3 * u4
                                        b(je+j) = bx
                                        a(jj+j) = co1*(t8-u11) - si1*(u8+t11)
                                        b(jj+j) = si1*(t8-u11) + co1*(u8+t11)
                                        a(jy+j) = co4*(t8+u11) - si4*(u8-t11)
                                        b(jy+j) = si4*(t8+u11) + co4*(u8-t11)
                                        a(jo+j) = co2*(t9-u10) - si2*(u9+t10)
                                        b(jo+j) = si2*(t9-u10) + co2*(u9+t10)
                                        a(jt+j) = co3*(t9+u10) - si3*(u9-t10)
                                        b(jt+j) = si3*(t9+u10) + co3*(u9-t10)
                                        j = j + jump
                                    enddo
                                    !
                                endif
                                !
                                !-----(end of loop across transforms)
                                !
                                ja = ja + jstepx
                                if (ja.lt.istart) ja = ja + ninc
                            enddo
                        enddo
                    enddo
                    !-----( end of double loop for this k )
                    kk = kk + 2*la
                enddo
                !-----( end of loop over values of k )
                la = 5*la
            enddo
            !-----( end of loop on type II radix-5 passes )
            !-----( nvex transforms completed)
490         continue
            istart = istart + nvex * jump
        enddo
        !-----( end of loop on blocks of transforms )
        !
        return
    end subroutine gpfa5f
end module mondelette
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
