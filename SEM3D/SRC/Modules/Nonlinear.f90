!>
!!\file nonlinear.f90
!!\brief Contains subroutines for nonlinear calculation
!!
!<

module nonlinear

    use sdomain
    use deriv3d
    use constants
contains
    
    subroutine stiff_matrix(lambda,mu,DEL_ijhk)
        implicit none
        real,                     intent(in)  :: lambda,mu
        real, dimension(0:5,0:5), intent(out) :: DEL_ijhk
        real, dimension(0:2),     parameter   :: veci = (/ 1.0, 0.0, 0.0 /)
        real, dimension(0:2),     parameter   :: vecj = (/ 0.0, 1.0, 0.0 /)
        real, dimension(0:2),     parameter   :: veck = (/ 0.0, 0.0, 1.0 /)
        real, dimension(0:2,0:2), parameter   :: id_matrix = reshape( (/veci,vecj,veck/), (/3,3/) )
        real, dimension(0:2,0:2), parameter   :: M(0:2,0:2)=1d0
        
        ! COMPUTE ELASTIC STIFFNESS MATRIX
        DEL_ijhk(:,:)       = 0d0
        DEL_ijhk(0:2,0:2)   = DEL_ijhk(0:2,0:2)+lambda*M+id_matrix*2*mu
        DEL_ijhk(3:5,3:5)   = DEL_ijhk(3:5,3:5)+id_matrix*mu

    end subroutine stiff_matrix
    
    subroutine mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)

        ! MISES YIELDING LOCUS AND GRADIENT

        implicit none
        real, dimension(0:5), intent(in)  :: Sigma_ij    ! initial stress state
        real, dimension(0:5), intent(in)  ::     X_ij    ! current back-stress state
        real                              :: R           ! current yld locus size
        real,                 intent(in)  :: sigma_yld   ! first yield limit
        real,                 intent(out) :: F_mises     ! Mises yield locus
        real, dimension(0:5), intent(out) :: gradF_mises ! Mises yield locus gradient

        real, dimension(0:5)              :: Sigma_ij_dev
        real                              :: tau_eq_mises
        integer :: k
        real, dimension(0:5), parameter   :: A = (/1.0,1.0,1.0,2.0,2.0,2.0/)
        
        call tensor_components (Sigma_ij, Sigma_ij_dev)

        call tau_mises(Sigma_ij_dev-X_ij, tau_eq_mises)

        F_mises      =  tau_eq_mises-sigma_yld-R
        gradF_mises(0:5)=0d0

        do k=0,5
            gradF_mises(k) = 1.5*A(k)*(Sigma_ij_dev(k)-X_ij(k))/tau_eq_mises
        end do
    end subroutine mises_yld_locus

    subroutine tau_mises(tensor,tau_eq_mises)
        implicit none

        real, dimension(0:5), intent(in) :: tensor
        real,                 intent(out):: tau_eq_mises
        real, dimension(0:5), parameter  :: A = (/1.0,1.0,1.0,2.0,2.0,2.0/) 
        integer                          :: k
        
        tau_eq_mises = 0.0d0
        do k=0,5
            tau_eq_mises = tau_eq_mises+A(k)*(tensor(k)**2)
        end do
        tau_eq_mises = sqrt(1.5*tau_eq_mises)

    end subroutine

    subroutine tensor_components (Sigma_ij, Sigma_ij_dev)

        ! TENSOR COMPONENTS (SPHERICAL AND DEVIATORIC)

        implicit none
        real, dimension(0:5), intent(in)  :: Sigma_ij
        real, dimension(0:5), intent(out) :: Sigma_ij_dev
        real                              :: Sigma_P
        integer :: k

        Sigma_ij_dev(0:5) = Sigma_ij(0:5)
        Sigma_P=sum(Sigma_ij(0:2))/3
        do k=0,2
            Sigma_ij_dev(k) = Sigma_ij_dev(k)-Sigma_P
        end do

    end subroutine tensor_components

    subroutine check_plasticity (dSigma_ij_trial, Sigma_ij_start, X_ij, R, sigma_yld, &
        st_elp, alpha_elp, nx,ny,nz,nelement)

        ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)

        implicit none
        real,                 intent(in)    :: R,sigma_yld
        integer,              intent(in)    :: nx,ny,nz,nelement
        real, dimension(0:5), intent(in)    :: X_ij,Sigma_ij_start
        real, dimension(0:5), intent(inout) :: dSigma_ij_trial
        real,                 intent(out)   :: alpha_elp
        integer,              intent(out)   :: st_elp
        real, dimension(0:5)                :: gradFstart, gradFtrial, Sigma_ij_trial
        real, dimension(0:5), parameter     :: A=(/1.0,1.0,1.0,2.0,2.0,2.0/)
        real                                :: Fstart, Ftrial,checkload
        integer                             :: GLLx,GLLy,GLLz,k

        ! Yield function at Sigma_ij_start
        call mises_yld_locus (Sigma_ij_start, X_ij, R, sigma_yld, &
            Fstart, gradFstart)

        ! Yield function at Sigma_ij_trial
        Sigma_ij_trial=Sigma_ij_start+dSigma_ij_trial
        call mises_yld_locus(Sigma_ij_trial, X_ij, R, sigma_yld, &
            Ftrial, gradFtrial)
        write(*,*) "*********************************"
        write(*,*) "Fstart:",Fstart,"Ftrial:",Ftrial

        ! LOADING CONDITION    
        checkload=0d0
        do k=0,5
            checkload=checkload+10*gradFstart(k)*dSigma_ij_trial(k)
        end do
        ! KKT CONDITION
        if (abs(Fstart).le.tol_nl) then           ! ON F=0
            if (checkload.ge.0) then              ! PURE ELASTIC STEP 
                alpha_elp = 0d0
                st_elp    = 1
            else
                if (Ftrial.lt.tol_nl) then        ! REVERSE ELASTO-PLASTIC STEP
                    alpha_elp = 1d0
                    st_elp    = 2 
                else                              ! REVERSE PURE ELASTIC STEP
                    !alpha_elp = 2*sigma_yld/(2*sigma_yld+Ftrial)
                    Fstart=Fstart+sigma_yld
                    Ftrial=Ftrial+sigma_yld
                    write(*,*) "Fstart corrected:",Fstart,"Ftrial corrected:",Ftrial
                    call gotoFtan(Sigma_ij_start,dSigma_ij_trial,&
                        Fstart,Ftrial,X_ij,R,0d0,alpha_elp)
                    st_elp    = 1
                endif
            end if
        elseif (Fstart.lt.-tol_nl) then
            if (Ftrial.lt.tol_nl) then            ! PURE ELASTIC STEP
                alpha_elp = 1d0
                st_elp    = 2
            else                                  ! ELASO-PLASTIC STEP
                !call gotoFsec(Sigma_ij_start,dSigma_ij_trial,Fstart,Ftrial,X_ij,R,sigma_yld,alpha_elp)
                call gotoFpegasus(Sigma_ij_start,dSigma_ij_trial,X_ij,R,sigma_yld,alpha_elp)  
                st_elp    = 1
            end if
        elseif (Fstart .gt. tol_nl) then
            write(*,*) "Fstart:",Fstart,"Ftrial:",Ftrial
            write(*,*) "ERROR!"
            write(*,*) "Fstart: ",Fstart,">",tol_nl,"!!!!"
            write(*,*) "ERROR!"
            stop
        end if
        dSigma_ij_trial(0:5)=Sigma_ij_start(0:5)+dSigma_ij_trial(0:5)*alpha_elp
        call mises_yld_locus(dSigma_ij_trial,X_ij,R,sigma_yld,Fstart,gradFstart)
        write(*,*) "check plasticity:",Fstart,"<=",tol_nl
        write(*,*) "ALPHA:",alpha_elp
        write(*,*) "*********************************"
        write(*,*) ""
    end subroutine check_plasticity

    subroutine plastic_corrector (dEpsilon_ij_alpha, Sigma_ij, X_ij, sigma_yld, &
        R, b_lmc, Rinf_lmc, C_lmc, kapa_lmc, mu, lambda, dEpsilon_ij_pl,flag_SS,&
        nelement,ngllx,nglly,ngllz)

        ! NR ALGORITHM AND DRIFT CORRECTION

        implicit none
        real, dimension(0:5), intent(inout) :: Sigma_ij             ! starting stress state
        real, dimension(0:5), intent(inout) :: X_ij                 ! starting back stress
        real,                 intent(inout) :: R                    ! starting mises radius
        real, dimension(0:5), intent(inout) :: dEpsilon_ij_alpha    ! percentage of elastic-plastic strain
        real, dimension(0:5), intent(inout) :: dEpsilon_ij_pl       ! percentage of plastic strain
        real,                 intent(in)    :: sigma_yld            ! first yield limit
        real,                 intent(in)    :: b_lmc, Rinf_lmc      ! Lamaitre and Chaboche parameters (isotropic hardening)
        real,                 intent(in)    :: C_lmc, kapa_lmc      ! Lamaitre and Chaboche parameters (kinematic hardening)
        real,                 intent(in)    :: mu, lambda           ! elastic parameters
        integer, optional,    intent(in)    :: nelement,ngllx,nglly,ngllz
        logical,              intent(in)    :: flag_SS
        real                                :: N_incr
        real, dimension(0:5)                :: dX_ij, gradF_0, gradF_mises
        real                                :: dR,dPlastMult,F_mises_0,F_mises,Ffinal,dPlastMult1
        integer                             :: i,j,k
        real, dimension(0:5)                :: temp_vec
        real, dimension(0:5,0:5)            :: DEL_ijhk
        real, dimension(0:5), parameter     :: A = (/1.0,1.0,1.0,0.5,0.5,0.5/)
        real, parameter                     :: stol = 0.001d0
        real                                :: Ttot,deltaTk,qq,R1,R2,dR1,dR2,err0,err1
        real(8)                             :: Resk
        logical                             :: flag_fail
        real, dimension(0:5)                :: S1,S2,X1,X2
        real, dimension(0:5)                :: dS1,dS2,dX1,dX2,dEpl1,dEpl2

        call stiff_matrix(lambda,mu,DEL_ijhk)
        
        if (flag_SS) then         ! CONSTANT INCREMENT
            
            N_incr=10
            
            dEpsilon_ij_alpha(0:5) = dEpsilon_ij_alpha(0:5)/N_incr
            
            do i = 0,N_incr-1

                ! PREDICTION
                call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises_0, gradF_0)
                ! COMPUTE PLASTIC MULTIPLIER
                call compute_plastic_modulus(dEpsilon_ij_alpha, Sigma_ij, X_ij, R, mu, lambda, sigma_yld, &
                    b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dPlastMult,nelement,ngllx,nglly,ngllz)
                
                ! HARDENING INCREMENTS
                call hardening_increments(Sigma_ij, R, X_ij, sigma_yld, &
                    b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)
                dR = dPlastMult*dR
                dX_ij(0:5) = dPlastMult*dX_ij(0:5)

                ! VARIABLE UPDATE
                do k=0,5 ! plastic strains
                    dEpsilon_ij_pl(k)=dEpsilon_ij_pl(k)+dPlastMult*gradF_0(k)*A(k)
                end do
                R=R+dR                          ! isotropic hardening update
                X_ij(0:5)=X_ij(0:5)+dX_ij(0:5)  ! back-stress update

                ! stress update
                do j = 0,5
                    do k = 0,5
                        Sigma_ij(j)=Sigma_ij(j)+DEL_ijhk(k,j)*(dEpsilon_ij_alpha(k)-A(k)*dPlastMult*gradF_0(k))
                    end do
                end do

                ! DRIFT CORRECTION
                call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, Ffinal, gradF_0)
!               write(*,*) "Ffinal (BD): ",Ffinal
                if (Ffinal .gt. tol_nl) then
                    call drift_corr(1,Sigma_ij, X_ij, R, sigma_yld,&
                        b_lmc, Rinf_lmc, C_lmc, kapa_lmc, lambda, mu, dEpsilon_ij_pl)
                    call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, Ffinal, gradF_0)
                end if
                if (Ffinal.gt.tol_nl) then
                    write(*,*) "Ffinal (AD): ",Ffinal
                    write(*,*) ""
                endif

            end do
        else    ! ERROR ADAPTIVE CONTROL
            
            deltaTk = 1.0d0
            Ttot    = 0.0d0

            do while (Ttot.lt.0.99999d0)
                write(*,*) "correction rk",deltaTk 
                write(*,*) "total time",Ttot
                
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
                call ep_integration(dEpsilon_ij_alpha*deltaTk,Sigma_ij,X_ij,R,&
                    sigma_yld,mu,lambda,b_lmc,Rinf_lmc,C_lmc,kapa_lmc,        &
                    dS1,dX1,dR1,dEpl1)

                R1      = R+dR1
                X1(0:5) = X_ij(0:5) + dX1(0:5) 
                S1(0:5) = Sigma_ij(0:5) + dS1(0:5)
               
                ! SECOND ORDER COMPUTATION
                call ep_integration(dEpsilon_ij_alpha*deltaTk,S1,X1,R1,&
                    sigma_yld,mu,lambda,b_lmc,Rinf_lmc,C_lmc,kapa_lmc, &
                    dS2,dX2,dR2,dEpl2)
                
                ! TEMPORARY VARIABLES
                S1 = Sigma_ij + 0.5d0*(dS1+dS2)
                X1 = X_ij     + 0.5d0*(dX1+dX2)
                R1 = R        + 0.5d0*(dR1+dR2)
                dEpl1 = 0.5d0*(dEpl1+dEpl2)

                ! ERROR
                call tau_mises(dS2-dS1,err0)
                call tau_mises(S1,err1)

                Resk=0.5d0*max(epsilon(Resk),err0/err1)
!                call tau_mises(dX2-dX1,err0)
!                call tau_mises(X1,err1)
!
!                Resk=0.5d0*max(Resk,err0/err1)
                
                write(*,*) "residuum",Resk
                if (Resk.le.stol) then ! substep is ok
                    
                    Sigma_ij = S1
                    X_ij     = X1
                    R        = R1

!                    ! DRIFT CORRECTION
                    call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, Ffinal, gradF_0)
                    if (Ffinal .gt. tol_nl) then
                        call drift_corr(0,Sigma_ij, X_ij, R, sigma_yld,&
                                b_lmc, Rinf_lmc, C_lmc, kapa_lmc, lambda, mu, dEpsilon_ij_pl)
                        call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, Ffinal, gradF_0)
                        write(*,*) "after drift",Ffinal
                    endif
                   
                    ! UPDATE TO NEXT TIME STEP
                    Ttot=Ttot+deltaTk
                    qq = min(0.9d0*sqrt(stol/Resk),1.1d0) 
                    if (flag_fail) then
                        qq = min(qq,1.0d0)
                    endif
!                    qq=min(0.8d0*sqrt(Resk/stol),2.0d0)
                    flag_fail=.false.
                    deltaTk=qq*deltaTk
                    deltaTk=max(qq*deltaTk,0.1d0)
!                    deltaTk=min(deltaTk,1-Ttot)

                else    ! substep has failed
!                    qq=max(0.9d0*sqrt(stol/Resk),0.1d0)
                    qq=max(0.9d0*sqrt(Resk/stol),0.1d0)
                    deltaTk=qq*deltaTk
                    flag_fail=.true.
                end if
                write(*,*) "flag",flag_fail
            end do

        endif
    end subroutine plastic_corrector
    
    subroutine ep_integration(dStrain,Stress,center,radius,syld,mu,lambda,biso,Rinf,&
        Ckin,kapakin,dStress,dcenter,dradius,dEplast)
        
        implicit none
        real                , intent(in) :: radius,syld,mu,lambda,biso,Rinf,Ckin,kapakin
        real, dimension(0:5), intent(in) :: dStrain,Stress,center
        real, dimension(0:5), intent(inout):: dStress,dcenter,dEplast
        real,                 intent(inout):: dradius
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
            biso,Rinf,Ckin,kapakin,dPlast)
        
        ! INCREMENTS
        call hardening_increments(Stress,radius,center,syld, &
            biso,Rinf,Ckin,kapakin,dradius,dcenter)
        
        dradius     = dPlast*dradius
        dcenter(0:5)= dPlast*dcenter(0:5)
        dEplast(0:5)= 0d0
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


    subroutine compute_plastic_modulus(dEpsilon_ij, Sigma_ij, X_ij, R, mu, lambda, &
        sigma_yld, b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dPlastMult, nelement, ngllx, nglly, ngllz)

        ! HARDENING MODULUS AND PLASTIC MULTIPLIER INCREMENT

        implicit none
        real, dimension(0:5), intent(in) :: dEpsilon_ij         ! percentage of elastic-plastic strain
        real, dimension(0:5), intent(in) :: Sigma_ij            ! starting stress state
        real, dimension(0:5), intent(in) :: X_ij                ! starting back stress
        real,                 intent(in) :: R                   ! starting mises radius
        real,                 intent(in) :: b_lmc, Rinf_lmc     ! Lamaitre and Chaboche parameters (isotropic hardening)
        real,                 intent(in) :: C_lmc, kapa_lmc     ! Lamaitre and Chaboche parameters (kinematic hardening)
        real,                 intent(in) :: mu, lambda          ! elastic parameters
        real,                 intent(in) :: sigma_yld           ! first yielding limit
        integer, optional,       intent(in) :: nelement,ngllx,nglly,ngllz
        real,                 intent(out):: dPlastMult          ! plastic multiplier increment
        real                             :: h_iso, h_kin        ! isotropic and kinematic hardening modula
        real                             :: F_mises
        real, dimension(0:5)             :: gradF_mises
        real                             :: h_lmc               ! total hardening modulus
        real                             :: temp_vec
        real, dimension(0:5,0:5)         :: DEL_ijhk
        integer                          :: j,k
        real, dimension(0:5), parameter  :: A =(/1.0,1.0,1.0,0.5,0.5,0.5/)
        
        ! COMPUTE ELASTIC STIFFNESS MATRIX
        call stiff_matrix(lambda,mu,DEL_ijhk)

        ! COMPUTE HARDENING MODULUS
        call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)

        h_iso = b_lmc*(Rinf_lmc-R)
        h_kin = 0d0
        do k=0,5
            h_kin=h_kin+kapa_lmc*X_ij(k)*gradF_mises(k)
        end do
        h_kin=C_lmc-h_kin
        h_lmc = h_iso+h_kin
        ! COMPUTE PLASTIC MULTIPLIER
        temp_vec   = 0d0
        dPlastMult = 0d0

        do k = 0,5
            do j = 0,5
                temp_vec   = temp_vec+gradF_mises(j)*DEL_ijhk(j,k)*gradF_mises(k)*A(k)
                dPlastMult = dPlastMult+gradF_mises(j)*DEL_ijhk(j,k)*dEpsilon_ij(k)
            end do
        end do
        dPlastMult = max(0d0,dPlastMult/(h_lmc+temp_vec))

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

    subroutine drift_corr(drift,Sigma_ij, X_ij, R, sigma_yld, &
        b_lmc, Rinf_lmc, C_lmc, kapa_lmc, lambda, mu, dEpsilon_ij_pl)

        ! DRIFT CORRECTION (RADIAL RETURN)
        integer, intent(in)                 :: drift
        real, dimension(0:5), intent(inout) :: Sigma_ij, X_ij, dEpsilon_ij_pl
        real,                 intent(inout) :: R
        real,                 intent(in)    :: lambda, mu, sigma_yld, b_lmc, Rinf_lmc, C_lmc, kapa_lmc
        real, dimension(0:5)                :: gradF_mises,gradF0,Sigma_temp,Sigma_dev_temp,Sigma_dev_ij
        real, dimension(0:5),     parameter :: A = (/1.0,1.0,1.0,0.5,0.5,0.5/)
        real, dimension(0:5,0:5)            :: DEL_ijhk
        real :: Fmises,Fmises0,beta,h_kin,h_iso,h_lmc,dbeta,err0,err1
        integer :: k,j
        
        ! INITIAL PLASTIC CONDITION
        call mises_yld_locus(Sigma_ij, X_ij, R, sigma_yld, Fmises0, gradF0)
        !if (drift==1) then ! B METHOD by Potts & Gens 1985
        beta=0d0
        do 
            call mises_yld_locus(Sigma_ij,X_ij,R,sigma_yld,Fmises,gradF_mises)
            ! COMPUTE BETA FOR DRIFT CORRECTION
            dbeta=0d0
            do k=0,5
                dbeta = dbeta+gradF0(k)*gradF_mises(k)
            end do
            beta=beta-Fmises/dbeta
            ! STRESS CORRECTION (RADIAL RETURN)
            Sigma_temp(0:5)=Sigma_ij(0:5)+beta*gradF0(0:5)*A(0:5)
            call mises_yld_locus(Sigma_temp,X_ij,R,sigma_yld,Fmises,gradF_mises)
            call tensor_components(Sigma_temp, Sigma_dev_temp)
            call tensor_components(Sigma_ij,Sigma_dev_ij)
            call tau_mises(Sigma_dev_ij-Sigma_dev_temp,err0)
            call tau_mises(Sigma_dev_ij-X_ij,err1)
            err0=err0/err1
            Sigma_ij(0:5)=Sigma_temp(0:5)
            if (err0 .le. 0.0000001 .or. abs(Fmises) .le. tol_nl) then
                exit
            end if
        end do
        !elseif( ! E METHOD
            ! COMPUTE ELASTIC STIFFNESS MATRIX
!            call stiff_matrix(lambda,mu,DEL_ijhk)
!            DEL_ijhk(:,:) = 0d0
!            DEL_ijhk(0:2,0:2) = DEL_ijhk(0:2,0:2) + lambda * M + id_matrix *2*mu
!            DEL_ijhk(3:5,3:5) = DEL_ijhk(3:5,3:5) + id_matrix * mu
!            
!            do 
!                ! COMPUTE HARDENING INCREMENTS
!                h_iso = b_lmc*(Rinf_lmc-R)
!                h_kin = C_lmc
!                do k=0,5
!                    h_kin = h_kin-kapa_lmc*X_ij(k)*gradF0(k)
!                end do
!                h_lmc=h_iso+h_kin
!                write(*,*) "H: ", h_lmc
!               ! COMPUTE BETA FOR DRIFT CORRECTION
!                beta = 0d0
!                do j=0,5
!                    do k=0,5
!                        beta=beta+gradF0(k)*DEL_ijhk(k,j)*A(j)*gradF0(j)
!                    end do
!                end do
!                write(*,*) "beta",beta
!                 
!                beta=Fmises/(-h_lmc+beta)
!                write(*,*) "beta debug",beta
!                ! STRESS-STRAIN-HARDENING CORRECTION
!                do k=0,5
!                    do j=0,5
!                        Sigma_temp(k)=Sigma_ij(k)-beta*DEL_ijhk(j,k)*A(j)*gradF0(j)
!                    end do
!                end do
!                
!                call mises_yld_locus(Sigma_temp,X_ij,R,sigma_yld,Fmises,gradF_mises)
!                call tensor_components(Sigma_temp, Sigma_dev_temp)
!                call tensor_components(Sigma_ij,Sigma_dev_ij)
!                call tau_mises(Sigma_dev_ij-Sigma_dev_temp,err0)
!                call tau_mises(Sigma_dev_ij-X_ij,err1)
!                err0=err0/err1
!                Sigma_ij(0:5)=Sigma_temp(0:5)
!
!                dEpsilon_ij_pl(0:5)=dEpsilon_ij_pl(0:5)+beta*A(0:5)*gradF0(0:5)
!                X_ij(0:5)=X_ij(0:5)+beta*(2*A(0:5)*gradF0(0:5)*C_lmc/3-X_ij(0:5)*kapa_lmc)
!                R=R+beta*(Rinf_lmc-R)*b_lmc
!                call mises_yld_locus(Sigma_ij, X_ij, R, sigma_yld, Fmises0, gradF0)
!
!                if (abs(Fmises0) .lt. tol_nl .or. err0/err1.lt.tol_nl/100) exit
!            end do
!        end if
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
!            write(*,*) "err0",err0
!            write(*,*) "err1",err1
            err0=err0/err1
            start(0:5)=temp(0:5)
            call mises_yld_locus(start,center,radius,s0,Fstart,gradF)
            if (abs(Fstart).le.tol_nl .or. err0.lt.0.000001) then
                exit
            end if
            beta=sum(gradF(0:5)*dtrial0(0:5))
            alpha=alpha-Fstart/beta
            temp(0:5)=start(0:5)-(Fstart/beta)*dtrial0(0:5)
        end do
!        write(*,*) "F:",Fstart
!        write(*,*) "err0",err0
!        write(*,*) "gotoF: stress:",start
!        write(*,*) ""
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
            if (abs(Ftrial).le.tol_nl .or. err0.lt.0.0000001) exit
            alphanew=alpha-(Ftrial/(Ftrial-F0))*(alpha-beta)
            beta=alpha
            alpha=alphanew
            temp(0:5)=start(0:5)-(alpha-beta)*dtrial0(0:5)
            
        enddo

    end subroutine gotoFsec

    subroutine gotoFpegasus(start0,dtrial,center,radius,s0,alpha)
        implicit none
        real, dimension(0:5), intent(in)    :: start0,dtrial,center
        real,                 intent(in)    :: radius,s0
        real,                 intent(out)   :: alpha
        real, dimension(0:5)                :: stress0,stress1,stress,gradF
        real                                :: alpha0,alpha1,F0,F1,FM
        integer                             :: counter

        alpha1  = 1d0
        alpha0  = 0d0
        stress0 = start0+alpha0*dtrial
        stress1 = start0+alpha1*dtrial
        call mises_yld_locus(stress0,center,radius,s0,F0,gradF)
        call mises_yld_locus(stress1,center,radius,s0,F1,gradF)
        
        do counter=0,9
            alpha  = alpha1-F1*(alpha1-alpha0)/(F1-F0)
            stress = start0+alpha*dtrial
            call mises_yld_locus(stress,center,radius,s0,FM,gradF)

            if (abs(FM).le.tol_nl) then
                exit
            else
                alpha0  =   alpha1
                alpha1  =   alpha
                F0      =   F1
                F1      =   FM
            endif

        end do
        if (FM.gt.tol_nl) then
            write(*,*) "WARNING: F>TOL"
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
