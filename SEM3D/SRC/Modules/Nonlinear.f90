!>
!!\file nonlinear.f90
!!\brief Contains subroutines for nonlinear calculation
!!
!<

module nonlinear

    use sdomain
    use deriv3d
    use constants

    implicit none
    real(KIND=8), parameter :: FTOL = 0.00010000000000D0
    real(KIND=8), parameter :: LTOL = 0.0000010000000000D0
    real(KIND=8), parameter :: STOL = 0.00010000000000D0
    real(KIND=8), parameter :: PSI  = one!5.0D0
    real(KIND=8), parameter :: OMEGA= zero!1.0D6
    real(fpp), dimension(0:2),     parameter   :: veci = (/ one, zero, zero /)
    real(fpp), dimension(0:2),     parameter   :: vecj = (/ zero, one, zero /)
    real(fpp), dimension(0:2),     parameter   :: veck = (/ zero, zero, one /)
    real(fpp), dimension(0:2,0:2), parameter   :: id_matrix = reshape( (/veci,vecj,veck/), (/3,3/) )
    real(fpp), dimension(0:2,0:2), parameter   :: Mmatrix(0:2,0:2) = one
    real(fpp), dimension(0:5),     parameter   :: Mvector  = (/one,one,one,zero,zero,zero/)
    real(fpp), dimension(0:5),     parameter   :: Avector  = (/one,one,one,two,two,two/)
    real(fpp), dimension(0:5),     parameter   :: A1vector = (/one,one,one,half,half,half/)
    real(fpp), parameter, dimension(0:5,0:5) :: &
        Amatrix  = reshape((/ &
        one , zero, zero, zero, zero, zero,&
        zero, one , zero, zero, zero, zero, &
        zero, zero, one , zero, zero, zero, &
        zero, zero, zero, two , zero, zero, &
        zero, zero, zero, zero, two , zero, &
        zero, zero, zero, zero, zero, two   &
        /), (/6,6/))
    real(fpp), parameter, dimension(0:5,0:5) :: &
        A1matrix = reshape((/ &
        one , zero, zero, zero, zero, zero,&
        zero, one , zero, zero, zero, zero, &
        zero, zero, one , zero, zero, zero, &
        zero, zero, zero, half, zero, zero, &
        zero, zero, zero, zero, half, zero, &
        zero, zero, zero, zero, zero, half  &
        /), (/6,6/))

contains

    !****************************************************************************
    ! UPDATE STRESS
    !****************************************************************************

    subroutine update_stress(stress0,stress1,dincrement)
        implicit none
        ! intent IN
        real(fpp), dimension(0:5), intent(in)    :: stress0,dincrement
        ! intent INOUT
        real(fpp), dimension(0:5), intent(inout) :: stress1
        !
        ! dincrement = stress increment
        stress1 = zero
        stress1 = stress0 + dincrement
        !
        return
        !
    end subroutine update_stress

    !****************************************************************************
    ! UPDATE STRESS CRITICAL
    !****************************************************************************

    subroutine update_stress_critical(stress0,stress1,dincrement,lambda,mu)
        implicit none
        ! intent IN
        real(fpp),                 intent(in)    :: lambda,mu
        real(fpp), dimension(0:5), intent(in)    :: stress0,dincrement
        ! intent INOUT
        real(fpp), dimension(0:5), intent(inout) :: stress1
        !
        real(fpp), dimension(0:5,0:5)            :: DEL
        !
        ! dincrement = strain increment
        stress1 = zero
        call stiff_matrix_critical(stress0,dincrement,lambda,mu,DEL)
        stress1 = stress0 + matmul(DEL,dincrement)
        !
        return
        !
    end subroutine update_stress_critical

    !*********************************************************************************
    ! STIFFNESS MATRIX
    !*********************************************************************************

    subroutine stiff_matrix(lambda,mu,DEL)
        ! intent IN
        real(fpp),                     intent(in)  :: lambda,mu
        ! intent OUT
        real(fpp), dimension(0:5,0:5), intent(out) :: DEL

        DEL(0:5,0:5) = zero
        DEL(0:2,0:2) = DEL(0:2,0:2) + lambda*Mmatrix
        DEL(0:5,0:5) = DEL(0:5,0:5) + two*mu*A1matrix

        return
    end subroutine stiff_matrix

    !*********************************************************************************
    ! CRITICAL STATE STIFFNESS MATRIX
    !*********************************************************************************

    subroutine stiff_matrix_critical(stress0,dstrain,lambda,mu,DEL)
        ! intent IN
        real(fpp), intent(in) :: lambda,mu
        real(fpp), dimension(0:5) :: stress0,dstrain
        ! intent OUT
        real(fpp), intent(out), dimension(0:5,0:5) :: DEL
        !
        real(fpp) :: p0,devol,B,k,nu,spec_vol,mu_crit,lambda_crit
        spec_vol = one+0.7d0
        ! pressure
        p0 = dot_product(stress0,Mvector)/three
        ! volumetric strain increment
        devol = dot_product(dstrain,Mvector)
        ! Bulk's modulus
        B = p0/devol*(exp(spec_vol*devol/k)-one)
        ! Poisson's ratio
        nu = half*lambda/(lambda+mu)
        ! NEW shear modulus
        mu_crit = three*half*(one-two*nu)*B/(one+nu)
        ! NEW lambda
        lambda_crit = two*mu_crit*nu/(one-two*nu)
        call stiff_matrix(lambda_crit,mu_crit,DEL)
        !
        return
        !
    end subroutine stiff_matrix_critical

    !****************************************************************************
    ! MISES YIELDING LOCUS AND GRADIENT
    !****************************************************************************

    subroutine mises_yld_locus(stress,center,radius,syld,FM,gradF)
        ! intent IN
        real(fpp),                 intent(in)  :: radius,syld
        real(fpp), dimension(0:5), intent(in)  :: stress,center
        ! intent OUT
        real(fpp),                 intent(out) :: FM
        real(fpp), dimension(0:5), intent(out) :: gradF
        !
        real(fpp), dimension(0:5)              :: dev
        real(fpp)                              :: tau_eq
        ! COMPUTE STRESS COMPONENTS
        call tensor_components(stress,dev)
        ! COMPUTE MISES FUNCTION
        call tau_mises(dev-center,tau_eq)
        ! COMPUTE MISES FUNCTION GRADIENT
        FM    = tau_eq - syld - radius
        gradF = three*half*Avector*(dev-center)/tau_eq
        !
        return
        !
    end subroutine mises_yld_locus

    !****************************************************************************
    ! OCTAHEDRAL SHEAR STRESS
    !****************************************************************************

    subroutine tau_mises(dev,J2M)
        ! intent IN
        real(fpp), dimension(0:5), intent(in) :: dev
        ! intent OUT
        real(fpp),                 intent(out):: J2M
        !
        real(fpp), dimension(0:5)             :: temp
        real(fpp)                             :: J2M2
        !
        temp = three*half*Avector*dev
        J2M2 = dot_product(dev,temp)
        J2M  = sqrt(J2M2)
        !
        return
        !
    end subroutine

    !****************************************************************************
    ! TENSOR COMPONENTS (SPHERICAL & DEVIATORIC)
    !****************************************************************************

    subroutine tensor_components(stress,dev)
        ! intent IN
        real(fpp), dimension(0:5), intent(in)    :: stress
        ! intent OUT
        real(fpp), dimension(0:5), intent(inout) :: dev
        real(fpp)                                :: press
        !
        press = dot_product(stress,Mvector)/three
        dev   = stress-press*Mvector
        !
        return
        !
    end subroutine tensor_components

    !****************************************************************************
    ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)
    !****************************************************************************

    subroutine check_plasticity(dtrial,stress0,center,radius,syld,st_epl,alpha_epl)
        ! intent IN
        real(fpp),                 intent(in)    :: radius,syld
        real(fpp), dimension(0:5), intent(in)    :: center,stress0
        ! intent INOUT
        real(fpp), dimension(0:5), intent(inout) :: dtrial
        ! intent OUT
        real(fpp),                 intent(out)   :: alpha_epl
        logical,              intent(out)   :: st_epl
        !
        logical                             :: flagxit
        real(fpp)                                :: FS,FT,checkload
        real(fpp), dimension(0:5)                :: gradFS,gradFT,stress1

        !
        ! PREDICTION STRESS
        call update_stress(stress0,stress1,dtrial)
        !
        ! CHECK MISES FUNCTION
        call mises_yld_locus(stress0,center,radius,syld,FS,gradFS)
        call mises_yld_locus(stress1,center,radius,syld,FT,gradFT)

        alpha_epl = -one
        if (FT.le.FTOL) then
            alpha_epl = one
            st_epl = .false.
            flagxit = .true.
        endif

        if ((FS.lt.-FTOL).and.(FT.gt.FTOL)) then
            st_epl = .true.
            call gotoFpegasus(stress0,dtrial,center,radius,syld,1,alpha_epl)
            flagxit = .true.
        endif

        if ((abs(FS).le.FTOL).and.(FT.gt.FTOL)) then
            ! CHECK LOAD DIRECTION
            checkload = dot_product(gradFS,dtrial)/&
                sqrt(dot_product(gradFS,gradFS)*dot_product(dtrial,dtrial))

            if (checkload.ge.-LTOL) then! PLASTIC LOADING
                alpha_epl = zero
            else! PLASTIC UNLOADING
                call gotoFpegasus(stress0,dtrial,center,radius,syld,10,alpha_epl)
            endif
            st_epl = .true.
            flagxit = .true.
        endif

        if (.not.flagxit)then
            write(*,*) "ERROR IN FINDING INTERSECTION!!  F = ",FS
            stop
        endif

        ! ON-LOCUS STRESS STATE
        call update_stress(stress0,stress1,alpha_epl*dtrial)
        dtrial = stress1
        call mises_yld_locus(dtrial,center,radius,syld,FS,gradFS)
        !
        return
        !
    end subroutine check_plasticity

    !****************************************************************************
    ! CORRECT STRESS STATE
    !****************************************************************************

    subroutine plastic_corrector (dEps_alpha,stress,center,syld, &
        radius,biso,Rinf,Ckin,kkin,mu,lambda,pstrain)

        implicit none
        real(fpp), dimension(0:5), intent(inout) :: dEps_alpha,stress,center,pstrain
        real(fpp),                 intent(inout) :: radius
        real(fpp),                 intent(in)    :: syld,biso,Rinf,Ckin,kkin,mu,lambda
        real(fpp), dimension(0:5,0:5)            :: DEL
        real(fpp)                                :: Ttot,deltaTk,qq,R1,dR1,dR2
        real(fpp)                                :: err0,err1,err2,err3
        real(fpp)                                :: FM,hard1,hard2,deltaTmin
        real(fpp)                                :: Resk
        logical                             :: flag_fail
        real(fpp), dimension(0:5)                :: gradFM,S1,X1,Epl1
        real(fpp), dimension(0:5)                :: dS1,dS2,dX1,dX2,dEpl1,dEpl2
        integer                             :: counter
        call stiff_matrix(lambda,mu,DEL)
        deltaTk = one
        Ttot    = zero
        deltaTmin = 0.001d0
        flag_fail = .false.
        counter = 1
        do while ((Ttot.lt.one).and.counter.le.10)
            Resk     = zero
            dS1(0:5) = zero
            dX1(0:5) = zero
            dS2(0:5) = zero
            dX2(0:5) = zero
            dR1      = zero
            dR2      = zero
            dEpl1(0:5) = zero
            dEpl2(0:5) = zero
            Epl1(0:5)  = zero
            ! FIRST ORDER COMPUTATION
            call ep_integration(dEps_alpha*deltaTk,stress,center,radius,syld,&
                mu,lambda,biso,Rinf,Ckin,kkin,dS1,dX1,dR1,dEpl1,hard1,pstrain)

            S1 = stress + dS1
            X1 = center + dX1
            R1 = radius + dR1
            Epl1 = pstrain + dEpl1

            ! SECOND ORDER COMPUTATION
            call ep_integration(dEps_alpha*deltaTk,S1,X1,R1,syld,mu,lambda,&
                biso,Rinf,Ckin,kkin,dS2,dX2,dR2,dEpl2,hard2,Epl1)

            ! TEMPORARY VARIABLES
            S1 = stress + half*(dS1+dS2)
            X1 = center + half*(dX1+dX2)
            R1 = radius + half*(dR1+dR2)
            Epl1 = pstrain + half*(dEpl1+dEpl2)

            ! ERROR
            call tau_mises(dS2-dS1,err0)
            call tau_mises(S1,err1)
            call tau_mises(dX2-dX1,err2)
            call tau_mises(X1,err3)
            Resk = max(epsilon(Resk),half*err0/err1,half*err2/err3)

            ! CHECK CONVERGENCE
            if (Resk.le.STOL) then
                stress = S1
                center = X1
                radius = R1
                pstrain = Epl1

                call mises_yld_locus (stress, center,radius,syld,FM,gradFM)
                if (FM.gt.FTOL) then
                    call drift_corr(stress,center,radius,syld,&
                            biso,Rinf,Ckin,kkin,lambda,mu,pstrain)
                endif

                qq = min(0.9d0*sqrt(STOL/Resk),1.1d0)
                if (flag_fail) then
                    qq = min(qq,one)
                endif
                flag_fail=.false.
                counter=1
                Ttot=Ttot+deltaTk
                deltaTk=qq*deltaTk
                deltaTk=max(deltaTk,deltaTmin)
                deltaTk=min(deltaTk,one-Ttot)
            else
                counter = counter+1
                qq=max(0.90d0*sqrt(STOL/Resk),0.1d0)
                deltaTk=max(qq*deltaTk,deltaTmin)
                flag_fail=.true.
            end if
        end do
        if (counter.eq.10)then
            write(*,*) "FAILED CORRECTION"
            stop
        endif
        !
        return
        !
    end subroutine plastic_corrector

    !****************************************************************************
    ! ELASTO-PLASTIC INTEGRATOR
    !****************************************************************************

    subroutine ep_integration(dStrain,Stress,center,radius,syld,mu,lambda,biso,Rinf,&
        Ckin,kkin,dStress,dcenter,dradius,dpstrain,hard,pstrain)
        ! intent IN
        real(fpp)                , intent(in) :: radius,syld,mu,lambda,biso,Rinf,Ckin,kkin
        real(fpp), dimension(0:5), intent(in) :: dstrain,stress,center,pstrain
        ! intent INOUT
        real(fpp), dimension(0:5), intent(inout):: dstress,dcenter,dpstrain
        real(fpp),                 intent(inout):: dradius
        ! intent OUT
        real(fpp),                 intent(out)  :: hard
        real(fpp), dimension(0:5)             :: gradF
        real(fpp), dimension(0:5,0:5)         :: DEL
        real(fpp)                             :: Fmises,dPlast

        ! PREDICTION
        call mises_yld_locus (Stress,center,radius,syld,Fmises,gradF)
        call stiff_matrix(lambda,mu,DEL)

        ! PLASTIC MULTIPLIER
        call compute_plastic_modulus(dStrain,Stress,center,radius,mu,lambda,syld, &
            biso,Rinf,Ckin,kkin,dPlast,hard,pstrain)

        ! INCREMENTS
        call hardening_increments(Stress,radius,center,syld, &
            biso,Rinf,Ckin,kkin,dradius,dcenter,pstrain)

        dradius      = dPlast*dradius
        dcenter(0:5) = dPlast*dcenter(0:5)
        dpstrain(0:5) = dPlast*A1vector(0:5)*gradF(0:5)
        dstress = matmul(DEL,dstrain-dpstrain)
        !
        return
        !
    end subroutine ep_integration

    !****************************************************************************
    ! PLASTIC MULTIPLIER
    !****************************************************************************

    subroutine compute_plastic_modulus(dstrain,stress,center,radius,mu,lambda, &
        syld,biso,Rinf,Ckin,kkin,dPlastM,hard,pstrain)
        ! intent IN
        real(fpp),                 intent(in) :: mu,lambda,syld
        real(fpp),                 intent(in) :: radius,biso,Rinf,Ckin,kkin
        real(fpp), dimension(0:5), intent(in) :: dstrain,stress,center,pstrain
        ! intent OUT
        real(fpp),                 intent(out):: dPlastM,hard
        !
        real(fpp)                             :: Ah,FM,PHI,PlastM
        real(fpp), dimension(0:5)             :: tempv,gradFM
        real(fpp), dimension(0:5,0:5)         :: DEL

        call stiff_matrix(lambda,mu,DEL)

        call mises_yld_locus(stress,center,radius,syld,FM,gradFM)

        PlastM = sqrt(two/three*dot_product(pstrain,pstrain))
        PHI  = one+(PSI-one)*exp(-OMEGA*PlastM)
        hard = PHI*Ckin+biso*(Rinf-radius)
        hard = hard-kkin*dot_product(center,gradFM)

        dPlastM = zero
        tempv(0:5) = zero
        !
        tempv(0:5) = A1vector*gradFM
        tempv(0:5) = matmul(DEL,tempv)
        Ah = dot_product(tempv,gradFM)
        !
        tempv(0:5) = zero
        tempv(0:5) = matmul(DEL,dstrain)
        dPlastM = dot_product(gradFM,tempv)
        dPlastM = dPlastM/(Ah+hard)
        dPlastM = max(zero,dPlastM)
        !
        return
        !
    end subroutine compute_plastic_modulus


    !****************************************************************************
    ! HARDENING INCREMENTS
    !****************************************************************************

    subroutine hardening_increments(stress,radius,center,syld, &
                biso,rinf,ckin,kkin,dradius,dcenter,pstrain)

        ! INCREMENTS OF INTRINSIC STATIC VARIABLES
        ! intent IN
        real(fpp),                 intent(in) :: syld,radius
        real(fpp),                 intent(in) :: biso,rinf,ckin,kkin
        real(fpp), dimension(0:5), intent(in) :: stress,center,pstrain
        ! intent OUT
        real(fpp),                 intent(out):: dradius
        real(fpp), dimension(0:5), intent(out):: dcenter
        !
        real(fpp)                             :: FM,PlastM,PHI
        real(fpp), dimension(0:5)             :: gradFM

        ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
        dradius = biso*(rinf-radius)

        ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
        PlastM = sqrt(two/three*dot_product(pstrain,pstrain))
        PHI  = one+(PSI-one)*exp(-OMEGA*PlastM)
        call mises_yld_locus (stress,center,radius,syld,FM,gradFM)
        dcenter(0:5) = (two*PHI*ckin/three)*A1vector(0:5)*gradFM(0:5)-kkin*center(0:5)
        !
        return
        !
    end subroutine hardening_increments

    !****************************************************************************
    ! DRIFT CORRECTION (RADIAL RETURN)
    !****************************************************************************

    subroutine drift_corr(stress,center,radius,syld, &
        biso,Rinf,Ckin,kkin,lambda,mu,pstrain)
        ! intent IN
        real(fpp),                 intent(in)    :: lambda,mu,syld,biso,Rinf,Ckin,kkin
        ! intent INOUT
        real(fpp), dimension(0:5), intent(inout) :: stress,center,pstrain
        real(fpp),                 intent(inout) :: radius
        real(fpp)                                :: F0,F1,beta,hard,radiust,PlastM,PHI
        real(fpp), dimension(0:5)                :: gradF0,gradF1,dstress,stresst,centert,tempv
        real(fpp), dimension(0:5,0:5)            :: DEL
        integer                             :: counter
        real(fpp), parameter                     :: FTOL_DRIFT =   0.000001D0
        ! INITIAL PLASTIC CONDITION
        call mises_yld_locus(stress,center,radius,syld,F0,gradF0)
        call stiff_matrix(lambda,mu,DEL)

        do counter=0,4
            ! MISES FUNCTION
            call mises_yld_locus(stress,center,radius,syld,F0,gradF0)
            ! COMPUTE HARDENING INCREMENTS
            PlastM = sqrt(two/three*dot_product(pstrain,pstrain))
            PHI  = one+(PSI-one)*exp(-OMEGA*PlastM)
            hard = biso*(Rinf-radius)
            hard = hard + PHI*Ckin
            hard = hard - kkin*dot_product(gradF0,center)
            ! COMPUTE BETA FOR DRIFT CORRECTION
            beta = zero
            tempv(0:5) = A1vector(0:5)*gradF0(0:5)
            tempv(0:5) = matmul(DEL,gradF0)
            beta = dot_product(tempv,gradF0)
            beta=F0/(hard+beta)
            ! STRESS-STRAIN-HARDENING CORRECTION
            dstress(0:5) = zero
            dstress(0:5) = -beta*tempv
            call update_stress(stress,stresst,dstress)
            centert = center+beta*((two*PHI*Ckin/three)*A1vector*gradF0-center*kkin)
            radiust = radius+beta*(Rinf-radius)*biso

            ! CHECK DRIFT
            call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
            if (abs(F1).gt.abs(F0)) then
                beta   = F0/dot_product(gradF0,gradF0)
                dstress = -beta*gradF0
                call update_stress(stress,stresst,dstress)
                centert = center
                radiust = radius
                call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
            endif
            stress = stresst
            center = centert
            radius = radiust
            pstrain = pstrain+beta*A1vector*gradF0
            if (abs(F1).le.FTOL_DRIFT) then
                exit
            endif
        enddo
        if (abs(F1).gt.FTOL) then
            write(*,*) "DRIFT NOT CORRECTED"
        endif
        !
        return
        !
    end subroutine drift_corr

    subroutine gotoFtan(start0,dtrial0,F0,Ftrial,center,radius,s0,alpha)
        real(fpp), dimension(0:5), intent(in)    :: start0,dtrial0,center
        real(fpp),                 intent(in)    :: radius,s0
        real(fpp),                 intent(inout) :: F0,Ftrial
        real(fpp),                 intent(out)   :: alpha
        real(fpp), dimension(0:5)                :: start,gradF,dev,dev_temp,temp,dev0
        real(fpp)                                :: Fstart,err0,err1,beta
        integer                             :: counter
        call tensor_components(start0,dev0)
        alpha=F0/(F0-Ftrial)
        start(0:5)=start0(0:5)
        temp=start+alpha*dtrial0
        do counter=0,9
            call tensor_components(temp, dev_temp)
            call tensor_components(start,dev)
            call tau_mises(-dev+dev_temp,err0)
            call tau_mises(dev-dev0,err1)
            err0=err0/err1
            start(0:5)=temp(0:5)
            call mises_yld_locus(start,center,radius,s0,Fstart,gradF)
            if (abs(Fstart).le.FTOL .or. err0.lt.0.000001) then
                exit
            end if
            beta=sum(gradF(0:5)*dtrial0(0:5))
            alpha=alpha-Fstart/beta
            temp(0:5)=start(0:5)-(Fstart/beta)*dtrial0(0:5)
        end do
    end subroutine gotoFtan

    subroutine gotoFsec(start0,dtrial0,F0,Ftrial,center,radius,s0,alpha)
        real(fpp), dimension(0:5), intent(in)    :: start0,dtrial0,center
        real(fpp),                 intent(in)    :: radius,s0
        real(fpp),                 intent(inout) :: F0,Ftrial
        real(fpp),                 intent(out)   :: alpha
        real(fpp), dimension(0:5)                :: start,gradF,dev,dev_temp,temp,dev0
        real(fpp)                                :: alphanew,err0,err1,beta
        integer                             :: counter
        call tensor_components(start0,dev0)
        alpha=F0/(F0-Ftrial)
        beta=0d0
        alphanew=0d0
        start(0:5)=start0(0:5)
        temp(0:5)=start(0:5)+alpha*dtrial0(0:5)
        do counter=0,9
            call tensor_components(temp, dev_temp)
            call tensor_components(start,dev)
            call tau_mises(-dev+dev_temp,err0)
            call tau_mises(dev-dev0,err1)
            err0=err0/err1
            call mises_yld_locus(start,center,radius,s0,F0,gradF)
            call mises_yld_locus(temp,center,radius,s0,Ftrial,gradF)
            start(0:5)=temp(0:5)
            if (abs(Ftrial).le.FTOL .or. err0.lt.0.0000001) exit
            alphanew=alpha-(Ftrial/(Ftrial-F0))*(alpha-beta)
            beta=alpha
            alpha=alphanew
            temp(0:5)=start(0:5)-(alpha-beta)*dtrial0(0:5)

        enddo

    end subroutine gotoFsec

    !****************************************************************************
    ! FIND INTERSECTION
    !****************************************************************************

    subroutine gotoFpegasus(start0,dtrial,center,radius,s0,nsub,alpha)
        real(fpp), dimension(0:5), intent(in)    :: start0,dtrial,center
        real(fpp),                 intent(in)    :: radius,s0
        integer,              intent(in)    :: nsub
        real(fpp),                 intent(out)   :: alpha
        real(fpp), dimension(0:5)                :: stress0,stress1,stress,gradF
        real(fpp)                                :: dalpha,alpha0,alpha1,F0,F1,FM,Fsave
        integer                             :: counter0,counter1
        logical                             :: flagxit

        alpha0  = zero
        alpha1  = one
        call update_stress(start0,stress0,alpha0*dtrial)
        call update_stress(start0,stress1,alpha1*dtrial)

        call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
        call mises_yld_locus(stress1,center,radius,s0,F1,gradF)

        ! LOAD REVERSAL
        if (nsub.gt.1) then
            Fsave=F0
            do counter0=0,3
                dalpha = (alpha1-alpha0)/nsub
                flagxit=.false.
                do counter1=0,nsub-1
                    alpha  = alpha0+dalpha
                    call update_stress(start0,stress,alpha*dtrial)
                    call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                    if (FM.gt.FTOL) then
                        alpha1=alpha
                        if (F0.lt.-FTOL) then
                            F1=FM
                            flagxit=.true.
                        else
                            alpha0=zero
                            F0=Fsave
                        endif
                        exit
                    else
                        alpha0=alpha
                        F0=FM
                    endif
                end do
                if (flagxit) then
                    exit
                endif
            end do
            call update_stress(start0,stress0,alpha0*dtrial)
            call update_stress(start0,stress1,alpha1*dtrial)
            call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
            call mises_yld_locus(stress1,center,radius,s0,F1,gradF)

            if ((F0.lt.-FTOL).and.(F1.gt.FTOL)) then
            else
                write(*,*) "LOAD REVERSAL FAILED"
            endif
        endif

        ! ORIGINAL PEGASUS ALGORITHM
        do counter0=0,9
            alpha  = alpha1-F1*(alpha1-alpha0)/(F1-F0)
            call update_stress(start0,stress,alpha*dtrial)

            call mises_yld_locus(stress,center,radius,s0,FM,gradF)
            if (abs(FM).le.FTOL) then ! abs(FS)<=FTOL ---> INTERSECTION FOUND
                exit ! INTERSECTION FOUND
            else
                if (FM*F0.lt.zero) then
                    alpha1=alpha0
                    F1=F0
                else
                    F1 = F1*F0/(FM+F0)
                endif
                F0=FM
                alpha0=alpha
            endif
        end do
        if (FM.gt.FTOL) then
            write(*,*) "WARNING: F>TOL!!!!!!!!!"
        endif
    end subroutine gotoFpegasus

end module nonlinear

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
