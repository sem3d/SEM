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
    real(KIND=8), parameter :: FTOL = 0.0010000000000D0
    real(KIND=8), parameter :: LTOL = 0.00100000000000D0
    real(KIND=8), parameter :: STOL = 0.00000100000000D0
    real(KIND=8), parameter :: PSI  = 1.0D0!5.0D0
    real(KIND=8), parameter :: OMEGA= 0.0D0!1.0D6
    real, dimension(0:5),     parameter   :: A  = (/one,one,one,two,two,two/)
    real, dimension(0:5),     parameter   :: A1 = (/one,one,one,half,half,half/)
    real, dimension(0:2),     parameter   :: veci = (/ one, zero, zero /)
    real, dimension(0:2),     parameter   :: vecj = (/ zero, one, zero /)
    real, dimension(0:2),     parameter   :: veck = (/ zero, zero, one /)
    real, dimension(0:2,0:2), parameter   :: id_matrix = reshape( (/veci,vecj,veck/), (/3,3/) )
    real, dimension(0:2,0:2), parameter   :: M(0:2,0:2)=one
    real, dimension(0:5),     parameter   :: vect_id = (/one,one,one,zero,zero,zero/)    

contains
    ! STIFFNESS MATRIX 
    subroutine stiff_matrix(lambda,mu,DEL)
        
        real,                     intent(in)  :: lambda,mu
        real, dimension(0:5,0:5), intent(out) :: DEL
        
        DEL(:,:)       = zero
        DEL(0:2,0:2)   = DEL(0:2,0:2)+lambda*M+id_matrix*two*mu
        DEL(3:5,3:5)   = DEL(3:5,3:5)+id_matrix*mu
        return
    end subroutine stiff_matrix
    
    ! MISES YIELD LOCUS
    subroutine mises_yld_locus(stress,center,radius,syld,FM,gradF)
        
        real,                 intent(in)  :: radius,syld
        real,                 intent(out) :: FM
        real, dimension(0:5), intent(in)  :: stress,center
        real, dimension(0:5), intent(out) :: gradF
        real, dimension(0:5)              :: dev
        real                              :: tau_eq
        
        call tensor_components(stress,dev)
        call tau_mises(dev-center,tau_eq)

        FM = tau_eq - syld - radius
        gradF=three*half*A*(dev-center)/tau_eq
        return
    end subroutine mises_yld_locus
    
    ! MISES TAU
    subroutine tau_mises(stress,tau_eq)
        
        real, dimension(0:5), intent(in) :: stress
        real,                 intent(out):: tau_eq
        integer                          :: k
        
        tau_eq = zero
        tau_eq = dot_product(A,stress**2)
        tau_eq = sqrt(three*half*tau_eq)
        return
    end subroutine
    
    ! TENSOR COMPONENTS
    subroutine tensor_components(stress,dev)

        real, dimension(0:5), intent(in)  :: stress
        real, dimension(0:5), intent(out) :: dev
        real                              :: press
        integer                           :: k

        press = dot_product(stress,vect_id)/three
        dev = stress-press*vect_id
        return
    end subroutine tensor_components
    
    ! CHECK PLASTICITY
    subroutine check_plasticity (dtrial,stress0,center,radius,syld, &
        st_epl, alpha_epl)

        implicit none
        real,                 intent(in)    :: radius,syld
        real, dimension(0:5), intent(in)    :: center,stress0
        real, dimension(0:5), intent(inout) :: dtrial
        real,                 intent(out)   :: alpha_epl
        logical,              intent(out)   :: st_epl
        real, dimension(0:5)                :: gradFS,gradFT,stress1
        real                                :: FS,FT,checkload
        integer                             :: k
        logical                             :: flag
        stress1=stress0+dtrial
        call mises_yld_locus(stress0,center,radius,syld,FS,gradFS)
        call mises_yld_locus(stress1,center,radius,syld,FT,gradFT)
        checkload = dot_product(gradFS,dtrial)/&
            (dot_product(gradFS,gradFS)*dot_product(dtrial,dtrial))

        if (abs(FS).le.FTOL) then ! FS = 0
            if (checkload.ge.-LTOL) then ! PLASTIC LOADING
                alpha_epl = zero
                st_epl    = .true.
            else
                if (FT.lt.-FTOL) then ! ELASTIC UNLOADING
                    alpha_epl = one
                    st_epl    = .false.
                elseif(FT.gt.FTOL) then ! PLASTIC UNLOADING
                    call gotoFpegasus(stress0,dtrial,center,radius,syld,10,alpha_epl)  
                    st_epl    = .true.
                else
                    write(*,*) "ERROR: MOVING ON THE SURFACE"
                    stop
                endif
            end if
        elseif (FS.lt.-FTOL) then ! FS<0
            if (FT.le.FTOL) then  ! ELASTIC LOADING
                alpha_epl = one
                st_epl    = .false.
            else ! ELASTO-PLASTIC LOADING
                call gotoFpegasus(stress0,dtrial,center,radius,syld,1,alpha_epl)  
                st_epl    = .true.
            end if
        elseif (FS.gt.FTOL) then
            write(*,*) "*********************************"
            write(*,*) "Fstart:",FS,"Ftrial:",FT
            write(*,*) "Fstart: ",FS,">",FTOL,"!!!!"
            write(*,*) "Load condition:",checkload
            write(*,*) "ERROR!"
            write(*,*) "*********************************"
            write(*,*) ""
            stop
        end if
        ! RETURN TRIAL STRESS 
        dtrial=stress0+dtrial*alpha_epl
        call mises_yld_locus(dtrial,center,radius,syld,FS,gradFS)
        return
    end subroutine check_plasticity
    
    ! PLASTIC CORRECTOR
    subroutine plastic_corrector (dEps_alpha,stress,center,syld, &
        radius,biso,Rinf,Ckin,kkin,mu,lambda,pstrain)

        implicit none
        real, dimension(0:5), intent(inout) :: dEps_alpha,stress,center,pstrain 
        real,                 intent(inout) :: radius 
        real,                 intent(in)    :: syld,biso,Rinf,Ckin,kkin,mu,lambda
        real, dimension(0:5,0:5)            :: DEL
        real                                :: Ttot,deltaTk,qq,R1,R2,dR1,dR2
        real                                :: err0,err1,err2,err3
        real                                :: FM,hard1,hard2,deltaTmin
        real(8)                             :: Resk
        logical                             :: flag_fail
        real, dimension(0:5)                :: gradFM,S1,S2,X1,X2,Epl1
        real, dimension(0:5)                :: dS1,dS2,dX1,dX2,dEpl1,dEpl2
        integer                             :: counter
        call stiff_matrix(lambda,mu,DEL)
        deltaTk = one
        Ttot    = zero
        deltaTmin = 0.001d0
        flag_fail = .false.
        counter = 1

        do while (Ttot.lt.one)
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

            Resk = half*max(epsilon(Resk),err0/err1,err2/err3)
            
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
                Ttot=Ttot+deltaTk
                deltaTk=qq*deltaTk
                deltaTk=max(qq*deltaTk,deltaTmin)
                deltaTk=min(deltaTk,one-Ttot)

            else
                qq=max(0.90d0*sqrt(STOL/Resk),deltaTmin)
                deltaTk=qq*deltaTk
                flag_fail=.true.
            end if

        end do
    end subroutine plastic_corrector
    
    ! ELASTO-PLASTIC INTEGRATION
    subroutine ep_integration(dStrain,Stress,center,radius,syld,mu,lambda,biso,Rinf,&
        Ckin,kapakin,dStress,dcenter,dradius,dEplast,hard,pstrain)
        
        implicit none
        real                , intent(in) :: radius,syld,mu,lambda,biso,Rinf,Ckin,kapakin
        real, dimension(0:5), intent(in) :: dstrain,stress,center,pstrain
        real, dimension(0:5), intent(inout):: dstress,dcenter,deplast
        real,                 intent(inout):: dradius
        real,                 intent(out)  :: hard
        real, dimension(0:5)             :: gradF
        real, dimension(0:5,0:5)         :: DEL
        real                             :: Fmises,dPlast
        integer                          :: j,k
        
        ! PREDICTION
        call mises_yld_locus (Stress,center,radius,syld,Fmises,gradF)
        call stiff_matrix(lambda,mu,DEL)

        ! PLASTIC MULTIPLIER
        call compute_plastic_modulus(dStrain,Stress,center,radius,mu,lambda,syld, &
            biso,Rinf,Ckin,kapakin,dPlast,hard,pstrain)
        
        ! INCREMENTS
        call hardening_increments(Stress,radius,center,syld, &
            biso,Rinf,Ckin,kapakin,dradius,dcenter,pstrain)
        
        dradius      = dPlast*dradius
        dcenter(0:5) = dPlast*dcenter(0:5)
        deplast(0:5) = dPlast*A1*gradF(0:5)
        dstress = matmul(DEL,dstrain-deplast)
        
        return

    end subroutine ep_integration


    subroutine compute_plastic_modulus(dEps,stress,center,radius,mu,lambda, &
        syld,biso,Rinf,Ckin,kkin,dPlastMult,hard,pstrain)

        implicit none
        real,                 intent(in) :: mu,lambda,syld   
        real,                 intent(in) :: radius,biso,Rinf,Ckin,kkin
        real, dimension(0:5), intent(in) :: dEps,stress,center,pstrain
        real,                 intent(out):: dPlastMult,hard
        real                             :: temp_vec,FM
        real, dimension(0:5)             :: gradF
        real, dimension(0:5,0:5)         :: DEL
        real                             :: PHI,PlastM
        integer                          :: j,k
        
        call stiff_matrix(lambda,mu,DEL)
        call mises_yld_locus(stress,center,radius,syld,FM,gradF)
       
        PlastM = sqrt(two/three*dot_product(pstrain,pstrain))
        PHI  = one+(PSI-one)*exp(-OMEGA*PlastM)
        hard = Ckin+biso*(Rinf-radius)
        hard = hard-kkin*dot_product(center,gradF)
        
        temp_vec   = zero
        dPlastMult = zero

        do k = 0,5
            do j = 0,5
                temp_vec   = temp_vec+gradF(j)*DEL(j,k)*gradF(k)*A1(k)
                dPlastMult = dPlastMult+gradF(j)*DEL(j,k)*dEps(k)
            end do
        end do
        dPlastMult = max(zero,dPlastMult/(hard+temp_vec))
        return
    end subroutine compute_plastic_modulus


    subroutine hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
        b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij, pstrain)

        ! INCREMENTS OF INTRINSIC STATIC VARIABLES

        real, dimension(0:5), intent(in) :: Sigma_ij,pstrain! actual stress state
        real, dimension(0:5), intent(in) :: X_ij            ! actual back stress state
        real,                 intent(in) :: sigma_yld       ! first yielding limit
        real,                 intent(in) :: R               ! actual mises radius
        real,                 intent(in) :: b_lmc, Rinf_lmc ! Lamaitre and Chaboche parameters (isotropic hardening)
        real,                 intent(in) :: C_lmc, kapa_lmc ! Lamaitre and Chaboche parameters (kinematic hardening)

        real,                 intent(out):: dR              ! mises radius increment
        real, dimension(0:5), intent(out):: dX_ij           ! back stress increment
        
        real                             :: F_mises,PlastM,PHI
        real, dimension(0:5)             :: gradF_mises
        integer                          :: k

        ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
        dR = b_lmc*(Rinf_lmc-R)

        ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
        call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
        PlastM = sqrt(two/three*dot_product(pstrain,pstrain))
        PHI  = one+(PSI-one)*exp(-OMEGA*PlastM)
        dX_ij(0:5)=two*A1(0:5)*gradF_mises(0:5)*PHI*C_lmc/three-X_ij(0:5)*kapa_lmc
        return
    end subroutine hardening_increments

    subroutine drift_corr(stress,center,radius,syld, &
        biso,Rinf,Ckin,kkin,lambda,mu,pstrain)

        ! DRIFT CORRECTION (RADIAL RETURN)
        real, dimension(0:5), intent(inout) :: stress,center,pstrain
        real,                 intent(inout) :: radius
        real,                 intent(in)    :: lambda,mu,syld,biso,Rinf,Ckin,kkin
        real                                :: F0,F1,beta,hard,radiust,PlastM,PHI
        real, dimension(0:5)                :: gradF0,gradF1,dstress,stresst,centert
        real, dimension(0:5,0:5)            :: DEL
        integer                             :: k,j,counter
        real, parameter                     :: FTOL_DRIFT =   0.000001D0
        ! INITIAL PLASTIC CONDITION
        call mises_yld_locus(stress,center,radius,syld,F0,gradF0)
        call stiff_matrix(lambda,mu,DEL)
        do counter=0,4 
            ! COMPUTE HARDENING INCREMENTS
            PlastM = sqrt(two/three*dot_product(pstrain,pstrain))
            PHI  = one+(PSI-one)*exp(-OMEGA*PlastM)
            hard = biso*(Rinf-radius)
            hard = hard + PHI*Ckin
            hard = hard - kkin*dot_product(gradF0,center)
            ! COMPUTE BETA FOR DRIFT CORRECTION
            beta = zero 
            do j=0,5
                do k=0,5
                    beta=beta+gradF0(k)*DEL(k,j)*A1(j)*gradF0(j)
                end do
            end do
            beta=F0/(hard+beta)
            ! STRESS-STRAIN-HARDENING CORRECTION
            dstress=zero
            do k=0,5
                do j=0,5
                    dstress(k)=dstress(k)-beta*DEL(j,k)*A1(j)*gradF0(j)
                end do
            end do
            stresst = stress+dstress
            centert = center+beta*(two*A1*gradF0*PHI*Ckin/three-center*kkin)
            radiust = radius+beta*(Rinf-radius)*biso
            
            ! CHECK DRIFT
            call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
            if (abs(F1).gt.abs(F0)) then
                beta   = F0/dot_product(gradF0,gradF0)
                stresst = stress-beta*gradF0
                centert = center
                radiust = radius
                call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
            endif
            stress = stresst
            center = centert
            radius = radiust
            pstrain = pstrain+beta*A1*gradF0
            if (abs(F1).le.FTOL) then
                exit
            else
                F0     = F1
                gradF0 = gradF1
            endif
        enddo
        if (abs(F1).gt.FTOL) then
            write(*,*) "DRIFT NOT CORRECTED"
        endif

        return
    end subroutine drift_corr

    subroutine gotoFtan(start0,dtrial0,F0,Ftrial,center,radius,s0,alpha)
        real, dimension(0:5), intent(in)    :: start0,dtrial0,center
        real,                 intent(in)    :: radius,s0
        real,                 intent(inout) :: F0,Ftrial
        real,                 intent(out)   :: alpha
        real, dimension(0:5)                :: start,gradF,dev,dev_temp,temp,dev0
        real                                :: Fstart,err0,err1,beta
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
        real, dimension(0:5), intent(in)    :: start0,dtrial0,center
        real,                 intent(in)    :: radius,s0
        real,                 intent(inout) :: F0,Ftrial
        real,                 intent(out)   :: alpha
        real, dimension(0:5)                :: start,gradF,dev,dev_temp,temp,dev0
        real                                :: alphanew,err0,err1,beta
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

    subroutine gotoFpegasus(start0,dtrial,center,radius,s0,nsub,alpha)
        implicit none
        real, dimension(0:5), intent(in)    :: start0,dtrial,center
        real,                 intent(in)    :: radius,s0
        integer,              intent(in)    :: nsub
        real,                 intent(out)   :: alpha
        real, dimension(0:5)                :: stress0,stress1,stress,gradF
        real                                :: dalpha,alpha0,alpha1,F0,F1,FM,Fsave
        integer                             :: counter0,counter1
        logical                             :: flagxit

        alpha0  = zero 
        alpha1  = one
        stress0 = start0+alpha0*dtrial
        stress1 = start0+alpha1*dtrial
        call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
        call mises_yld_locus(stress1,center,radius,s0,F1,gradF)
        if (nsub.gt.1) then
            Fsave=F0
            do counter0=0,3
                dalpha = (alpha1-alpha0)/nsub
                flagxit=.false.
                do counter1=0,nsub-1
                    alpha  = alpha0+dalpha
                    stress = start0+alpha*dtrial
                    call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                    if (FM.gt.FTOL) then
                        alpha1=alpha
                        if (F0.lt.-FTOL) then
                            F1=FM
                            flagxit=.true.
                        else
                            alpha0=0.0d0
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
        end if

        do counter0=0,9
            alpha  = alpha1-F1*(alpha1-alpha0)/(F1-F0)
            stress = start0+alpha*dtrial
            call mises_yld_locus(stress,center,radius,s0,FM,gradF)
            if (abs(FM).le.FTOL) then ! abs(FS)<=FTOL ---> INTERSECTION FOUND
                exit ! INTERSECTION FOUND
            else
                if (FM*F1.lt.zero) then
                    alpha0=alpha1
                    F0=F1
                else
                    F0 = F0*half
                endif
                F1=FM
                alpha1=alpha
            endif
        end do
        if (FM.gt.FTOL) then
            write(*,*) "WARNING: F>TOL!!!!!!!!!"
        else
            write(*,*) "INTERCEPTED"
        endif
    end subroutine gotoFpegasus

end module nonlinear

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
