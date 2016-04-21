!>
!!\file nonlinear.f90
!!\brief Contains subroutines for nonlinear calculation
!!
!<

module nonlinear

    use sdomain
    use deriv3d
    use constants

   real(KIND=8), parameter :: FTOL = 0.0010000000000D0
   real(KIND=8), parameter :: LTOL = 0.00100000000000D0
   real(KIND=8), parameter :: STOL = 0.0000100000000D0
!   real(KIND=8), parameter :: FTOL = 0.0000010000000000D0
!   real(KIND=8), parameter :: LTOL = 0.0000010000000000D0
!   real(KIND=8), parameter :: STOL = 0.0010000000000000D0
contains
    
    subroutine stiff_matrix(lambda,mu,DEL)
        
        implicit none
        real,                     intent(in)  :: lambda,mu
        real, dimension(0:5,0:5), intent(out) :: DEL
        real, dimension(0:2),     parameter   :: veci = (/ 1.0, 0.0, 0.0 /)
        real, dimension(0:2),     parameter   :: vecj = (/ 0.0, 1.0, 0.0 /)
        real, dimension(0:2),     parameter   :: veck = (/ 0.0, 0.0, 1.0 /)
        real, dimension(0:2,0:2), parameter   :: id_matrix = reshape( (/veci,vecj,veck/), (/3,3/) )
        real, dimension(0:2,0:2), parameter   :: M(0:2,0:2)=1d0
        
        DEL(:,:)       = 0d0
        DEL(0:2,0:2)   = DEL(0:2,0:2)+lambda*M+id_matrix*2*mu
        DEL(3:5,3:5)   = DEL(3:5,3:5)+id_matrix*mu

    end subroutine stiff_matrix
    
    subroutine mises_yld_locus(stress,center,radius,syld,FM,gradF)
        
        implicit none
        real,                 intent(in)  :: radius,syld
        real,                 intent(out) :: FM
        real, dimension(0:5), intent(in)  :: stress,center
        real, dimension(0:5), intent(out) :: gradF
        real, dimension(0:5), parameter   :: A = (/1.0,1.0,1.0,2.0,2.0,2.0/)
        real, dimension(0:5)              :: dev
        real                              :: tau_eq
        integer                           :: k
        
        call tensor_components(stress,dev)
        call tau_mises(dev-center,tau_eq)

        FM = tau_eq-syld-radius
        gradF=1.5d0*A*(dev-center)/tau_eq

    end subroutine mises_yld_locus

    subroutine tau_mises(stress,tau_eq)
        
        implicit none
        real, dimension(0:5), intent(in) :: stress
        real,                 intent(out):: tau_eq
        real, dimension(0:5), parameter  :: A = (/1.0,1.0,1.0,2.0,2.0,2.0/) 
        integer                          :: k
        
        tau_eq = 0.0d0
        do k=0,5
            tau_eq = tau_eq+A(k)*(stress(k)**2)
        end do
        tau_eq = sqrt(1.5*tau_eq)

    end subroutine

    subroutine tensor_components(stress,dev)

        implicit none
        real, dimension(0:5), intent(in)  :: stress
        real, dimension(0:5), intent(out) :: dev
        real                              :: press
        integer                           :: k

        dev=stress
        press=sum(stress(0:2))/3
        dev(0:2)=dev(0:2)-press

    end subroutine tensor_components

    subroutine check_plasticity (dtrial,stress0,center,radius,syld, &
        st_elp, alpha_elp)

        implicit none
        real,                 intent(in)    :: radius,syld
        real, dimension(0:5), intent(in)    :: center,stress0
        real, dimension(0:5), intent(inout) :: dtrial
        real,                 intent(out)   :: alpha_elp
        integer,              intent(out)   :: st_elp
        real, dimension(0:5)                :: gradFS,gradFT,stress1
        real, dimension(0:5), parameter     :: A=(/1.0,1.0,1.0,2.0,2.0,2.0/)
        real                                :: FS,FT,checkload
        integer                             :: k
        logical                             :: flag
        stress1=stress0+dtrial
        call mises_yld_locus(stress0,center,radius,syld,FS,gradFS)
        call mises_yld_locus(stress1,center,radius,syld,FT,gradFT)
        checkload=sum(gradFS*dtrial)/sum(gradFS**2)/sum(dtrial**2)
        
        if (abs(FS).le.FTOL) then
            if (checkload.ge.-LTOL) then 
                alpha_elp = 0d0
                st_elp    = 1
            else
                if (FT.lt.-FTOL) then
                    alpha_elp = 1d0
                    st_elp    = 2 
                elseif(FT.gt.FTOL) then
                    write(*,*) "-------------------------------------------------"
                    write(*,*) "check load",checkload,"LTOL",LTOL
                    call gotoFpegasus(stress0,dtrial,center,radius,syld,10,alpha_elp)  
                    st_elp    = 1
                endif
            end if
        elseif (FS.lt.-FTOL) then
            if (FT.le.FTOL) then
                alpha_elp = 1d0
                st_elp    = 2
            else
                call gotoFpegasus(stress0,dtrial,center,radius,syld,1,alpha_elp)  
                st_elp    = 1
            end if
        elseif (FS.gt.FTOL) then
            write(*,*) "*********************************"
            write(*,*) "Fstart:",FS,"Ftrial:",FT
            write(*,*) "Fstart: ",FS,">",FTOL,"!!!!"
            write(*,*) "Load condition:",checkload
            write(*,*) "ERROR!"
            write(*,*) "*********************************"
            write(*,*) ""
        end if
        ! RETURN TRIAL STRESS 
        dtrial=stress0+dtrial*alpha_elp

        call mises_yld_locus(dtrial,center,radius,syld,FS,gradFS)
    end subroutine check_plasticity

    subroutine plastic_corrector (dEps_alpha,stress,center,syld, &
        radius,biso,Rinf,Ckin,kkin,mu,lambda,dEpl)

        implicit none
        real, dimension(0:5), intent(inout) :: dEps_alpha,stress,center,dEpl 
        real,                 intent(inout) :: radius 
        real,                 intent(in)    :: syld,biso,Rinf,Ckin,kkin,mu,lambda
        real, dimension(0:5,0:5)            :: DEL
        real, dimension(0:5), parameter     :: A = (/1.0,1.0,1.0,0.5,0.5,0.5/)
        real                                :: Ttot,deltaTk,qq,R1,R2,dR1,dR2,err0,err1
        real                                :: FM,hard1,hard2,deltaTmin
        real(8)                             :: Resk
        logical                             :: flag_fail
        real, dimension(0:5)                :: gradFM,S1,S2,X1,X2
        real, dimension(0:5)                :: dS1,dS2,dX1,dX2,dEpl1,dEpl2

        call stiff_matrix(lambda,mu,DEL)
        deltaTk = 1.0d0
        Ttot    = 0.0d0
        deltaTmin = 0.01d0
        do while (Ttot.lt.1d0-FTOL)
            Resk     = 0d0
            dS1(0:5) = 0d0
            dX1(0:5) = 0d0
            dS2(0:5) = 0d0
            dX2(0:5) = 0d0
            dR1      = 0d0
            dR2      = 0d0
            dEpl1(0:5) = 0d0
            dEpl2(0:5) = 0d0

            ! FIRST ORDER COMPUTATION
            call ep_integration(dEps_alpha*deltaTk,stress,center,radius,syld,&
                mu,lambda,biso,Rinf,Ckin,kkin,dS1,dX1,dR1,dEpl1,hard1)

            S1 = stress + dS1
            X1 = center + dX1 
            R1 = radius + dR1
           
            ! SECOND ORDER COMPUTATION
            call ep_integration(dEps_alpha*deltaTk,S1,X1,R1,syld,mu,lambda,&
                biso,Rinf,Ckin,kkin,dS2,dX2,dR2,dEpl2,hard2)
            
            ! TEMPORARY VARIABLES
            S1 = stress + 0.5d0*(dS1+dS2)
            X1 = center + 0.5d0*(dX1+dX2)
            R1 = radius + 0.5d0*(dR1+dR2)
            dEpl1 = 0.5d0*(dEpl1+dEpl2)

            ! ERROR
            call tau_mises(dS2-dS1,err0)
            call tau_mises(S1,err1)

            Resk=0.5d0*max(epsilon(Resk),err0/err1)
            
            if (Resk.le.STOL) then
                
                stress = S1
                center = X1
                radius = R1
                dEpl   = dEpl+dEpl1
                call mises_yld_locus (stress, center,radius,syld,FM,gradFM)
                if (FM.gt.FTOL) then
                    call drift_corr(stress,center,radius,syld,&
                            biso,Rinf,Ckin,kkin,lambda,mu,dEpl)
                endif
               
                Ttot=Ttot+deltaTk
                qq = min(0.9d0*sqrt(STOL/Resk),1.1d0) 
                if (flag_fail) then
                    qq = min(qq,1.0d0)
                endif
                flag_fail=.false.
                deltaTk=qq*deltaTk
                deltaTk=max(qq*deltaTk,deltaTmin)
                deltaTk=min(deltaTk,1d0-Ttot)

            else
                qq=max(0.9d0*sqrt(STOL/Resk),deltaTmin)
                deltaTk=qq*deltaTk
                flag_fail=.true.
            end if
        end do
    end subroutine plastic_corrector
    
    subroutine ep_integration(dStrain,Stress,center,radius,syld,mu,lambda,biso,Rinf,&
        Ckin,kapakin,dStress,dcenter,dradius,dEplast,hard)
        
        implicit none
        real                , intent(in) :: radius,syld,mu,lambda,biso,Rinf,Ckin,kapakin
        real, dimension(0:5), intent(in) :: dStrain,Stress,center
        real, dimension(0:5), intent(inout):: dStress,dcenter,dEplast
        real,                 intent(inout):: dradius
        real,                 intent(out)  :: hard
        real, dimension(0:5)             :: gradF
        real, dimension(0:5), parameter  :: A = (/1.0,1.0,1.0,0.5,0.5,0.5/)
        real, dimension(0:5,0:5)         :: DEL
        real                             :: Fmises,dPlast
        integer                          :: j,k
        
        ! PREDICTION
        call mises_yld_locus (Stress,center,radius,syld,Fmises,gradF)
        call stiff_matrix(lambda,mu,DEL)

        ! PLASTIC MULTIPLIER
        call compute_plastic_modulus(dStrain,Stress,center,radius,mu,lambda,syld, &
            biso,Rinf,Ckin,kapakin,dPlast,hard)
        
        ! INCREMENTS
        call hardening_increments(Stress,radius,center,syld, &
            biso,Rinf,Ckin,kapakin,dradius,dcenter)
        
        dradius     = dPlast*dradius
        dcenter(0:5)= dPlast*dcenter(0:5)
        do k=0,5 
            dEplast(k)=dEplast(k)+dPlast*gradF(k)*A(k)
        end do
        do j = 0,5 ! stress increment
            do k = 0,5
                dstress(j)=DEL(k,j)*(dstrain(j)-A(k)*dPlast*gradF(k))  
            end do
        end do
        
        return

    end subroutine ep_integration


    subroutine compute_plastic_modulus(dEps,stress,center,radius,mu,lambda, &
        syld,biso,Rinf,Ckin,kkin,dPlastMult,hard)

        implicit none
        real,                 intent(in) :: mu,lambda,syld   
        real,                 intent(in) :: radius,biso,Rinf,Ckin,kkin
        real, dimension(0:5), intent(in) :: dEps,stress,center
        real,                 intent(out):: dPlastMult,hard
        real                             :: temp_vec,FM
        real, dimension(0:5)             :: gradF
        real, dimension(0:5), parameter  :: A =(/1.0,1.0,1.0,0.5,0.5,0.5/)
        real, dimension(0:5,0:5)         :: DEL
        integer                          :: j,k
        
        call stiff_matrix(lambda,mu,DEL)
        call mises_yld_locus(stress,center,radius,syld,FM,gradF)
        
        hard = Ckin+biso*(Rinf-radius)
        hard = hard-kkin*sum(center*gradF)
        
        temp_vec   = 0d0
        dPlastMult = 0d0

        do k = 0,5
            do j = 0,5
                temp_vec   = temp_vec+gradF(j)*DEL(j,k)*gradF(k)*A(k)
                dPlastMult = dPlastMult+gradF(j)*DEL(j,k)*dEps(k)
            end do
        end do
        dPlastMult = max(0d0,dPlastMult/(hard+temp_vec))

    end subroutine compute_plastic_modulus


    subroutine hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
        b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)

        ! INCREMENTS OF INTRINSIC STATIC VARIABLES

        real, dimension(0:5), intent(in) :: Sigma_ij        ! actual stress state
        real, dimension(0:5), intent(in) :: X_ij            ! actual back stress state
        real,                 intent(in) :: sigma_yld       ! first yielding limit
        real,                 intent(in) :: R               ! actual mises radius
        real,                 intent(in) :: b_lmc, Rinf_lmc ! Lamaitre and Chaboche parameters (isotropic hardening)
        real,                 intent(in) :: C_lmc, kapa_lmc ! Lamaitre and Chaboche parameters (kinematic hardening)

        real,                 intent(out):: dR              ! mises radius increment
        real, dimension(0:5), intent(out):: dX_ij           ! back stress increment
        
        real                             :: F_mises
        real, dimension(0:5)             :: gradF_mises
        real, dimension(0:5), parameter  :: A = (/1.0,1.0,1.0,0.5,0.5,0.5/)
        integer                          :: k

        ! INCREMENT IN ISOTROPIC HARDENING VARIABLES (R)
        dR = b_lmc*(Rinf_lmc-R)

        ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
        call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
        dX_ij(0:5)=2*A(0:5)*gradF_mises(0:5)*C_lmc/3-X_ij(0:5)*kapa_lmc

    end subroutine hardening_increments

    subroutine drift_corr(stress,center,radius,syld, &
        biso,Rinf,Ckin,kkin,lambda,mu,dEplastic)

        ! DRIFT CORRECTION (RADIAL RETURN)
        real, dimension(0:5), intent(inout) :: stress,center,dEplastic
        real,                 intent(inout) :: radius
        real,                 intent(in)    :: lambda,mu,syld,biso,Rinf,Ckin,kkin
        real                                :: F0,F1,beta,hard,radiust
        real, dimension(0:5)                :: gradF0,gradF1,dstress,stresst,centert
        real, dimension(0:5),     parameter :: A = (/1.0,1.0,1.0,0.5,0.5,0.5/)
        real, dimension(0:5,0:5)            :: DEL
        integer                             :: k,j,counter
        
        ! INITIAL PLASTIC CONDITION
        call mises_yld_locus(stress,center,radius,syld,F0,gradF0)
        call stiff_matrix(lambda,mu,DEL)
        do counter=0,4 
            ! COMPUTE HARDENING INCREMENTS
            hard = biso*(Rinf-radius)
            hard = hard + Ckin
            hard = hard - kkin*sum(gradF0*center)
            ! COMPUTE BETA FOR DRIFT CORRECTION
            beta = 0d0
            do j=0,5
                do k=0,5
                    beta=beta+gradF0(k)*DEL(k,j)*A(j)*gradF0(j)
                end do
            end do
            beta=F0/(hard+beta)
            ! STRESS-STRAIN-HARDENING CORRECTION
            dstress=0d0
            do k=0,5
                do j=0,5
                    dstress(k)=dstress(k)-beta*DEL(j,k)*A(j)*gradF0(j)
                end do
            end do
            stresst = stress+dstress
            centert = center+beta*(2*A*gradF0*Ckin/3-center*kkin)
            radiust = radius+beta*(Rinf-radius)*biso
            
            ! CHECK DRIFT
            call mises_yld_locus(stresst,centert,radiust,syld,F1,gradF1)
            if (abs(F1).gt.abs(F0)) then
                beta   = F0/sum(gradF0*gradF0)
                stress = stress-beta*gradF0
            else
                stress = stresst
                center = centert
                radius = radiust
                dEplastic = dEplastic+beta*A*gradF0
                if (abs(F1).le.FTOL) then
                    exit
                else
                    F0     = F1
                    gradF0 = gradF1
                endif
            endif
            
        enddo
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

        alpha0  = 0d0
        alpha1  = 1d0
        stress0 = start0+alpha0*dtrial
        stress1 = start0+alpha1*dtrial
        call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
        call mises_yld_locus(stress1,center,radius,s0,F1,gradF)
        if (nsub.gt.1) then
            Fsave=F0
            write(*,*) "before pegasus-> find starting alpha0-alpha1"
            write(*,*) "F0",F0
            write(*,*) "F1",F1
            do counter0=0,2
                dalpha = (alpha1-alpha0)/nsub
                flagxit=.false.
                write(*,*) "dalpha",dalpha,"counter0",counter0
                do counter1=0,nsub-1
                    write(*,*) "counter1",counter1
                    alpha=alpha0+dalpha
                    stress=start0+alpha*dtrial
                    write(*,*) "alpha",alpha
                    write(*,*) "stress",stress
                    write(*,*) ""
                    write(*,*) "start0",start0
                    write(*,*) "dtrial",dtrial
                    call mises_yld_locus(stress,center,radius,s0,FM,gradF)
                    if (FM.gt.FTOL) then
                        alpha1=alpha
                        if (F0.lt.-FTOL) then
                            F1=FM
                            flagxit=.true.
                            write(*,*) "exit1"
                        else
                            alpha0=0d0
                            F0=Fsave
                            write(*,*) "exit2"
                        endif
                        exit
                        
                    else
                        alpha0=alpha
                        F0=FM
                    endif
                    write(*,*) "alpha0",alpha0,"alpha1",alpha1
                end do
                if (flagxit) then
                    exit
                endif
            end do
            if (.not.flagxit) then
                write(*,*) "ERROR IN FINDING F=0 (REVERSAL)"
                alpha1=1d0
                alpha0=0d0
            endif
        write(*,*) "-------------------------------------------------"
        
        end if

        do counter0=0,9
            alpha  = alpha1-F1*(alpha1-alpha0)/(F1-F0)
            stress = start0+alpha*dtrial
            call mises_yld_locus(stress,center,radius,s0,FM,gradF)
            if (abs(FM).le.FTOL) then
                exit
            else
                alpha0  =   alpha1
                alpha1  =   alpha
                F0      =   F1
                F1      =   FM
            endif

        end do
        if (FM.gt.FTOL) then
            write(*,*) "WARNING: F>TOL!!!!!!!!!"
        endif
    end subroutine gotoFpegasus

    subroutine tau_deepsoil(gamma_dp,gamma_R,gamma_rev,gamma_max,&
        flag_rev,tau_rev,mu,p1_dp,p2_dp,p3_dp,beta_dp,S_dp,tau_dp)
        real,       intent(in)  :: gamma_dp,gamma_R,gamma_max,gamma_rev
        real,       intent(in)  :: tau_rev,mu,p1_dp,p2_dp,p3_dp
        logical,    intent(in)  :: flag_rev
        real,       intent(out) :: tau_dp
        real                    :: Fred,temp

        if (.not.flag_rev) then ! LOADING
            tau_dp=(mu*gamma_dp)/(1+beta_dp*(gamma_dp/gamma_R)**S_dp)
        else                    ! UNLOADING
            tau_dp=tau_rev
            tau_dp=tau_dp+(mu*(gamma_dp-gamma_rev))/&
                (1+beta_dp*(gamma_max-gamma_R)**S_dp)
            call reduction_factor_deepsoil(gamma_max,p1_dp,p2_dp,p3_dp,mu,Fred)
            temp=(mu*(gamma_dp-gamma_rev))/&
                (1+beta_dp*(0.5*(gamma_dp-gamma_rev)*gamma_R)**S_dp)
            temp=temp-mu*(gamma_dp-gamma_rev)/(1+beta_dp*(gamma_max-gamma_R)**S_dp)
            tau_dp=tau_dp+Fred*temp
        endif
        
    end subroutine tau_deepsoil
    
    subroutine gamma_ref(stress,A_dp,C_dp,gamma_R)
        real, dimension(0:5), intent(in)    :: stress
        real,                 intent(in)    :: C_dp
        real,                 intent(out)   :: gamma_R
        real                                :: press
        press = sum(stress(0:2))/3 
        gamma_R=A_dp*(press/101325.0)**C_dp
    end subroutine gamma_ref
    
    subroutine reduction_factor_deepsoil(gamma_max,p1_dp,p2_dp,p3_dp,mu,Fred)
        real, intent(in) :: gamma_max,p1_dp,p2_dp,p3_dp,mu
        real, intent(out):: Fred
        Fred=p1_dp-p2_dp*(1-Ggamma/mu)**p3_dp
    end subroutine

end module nonlinear

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
