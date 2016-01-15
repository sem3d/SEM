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

    subroutine mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)

        ! MISES YIELDING LOCUS AND GRADIENT

        implicit none
        real, dimension(0:5), intent(in)  :: Sigma_ij    ! initial stress state
        real, dimension(0:5), intent(in)  ::     X_ij    ! current back-stress state
        real                              :: R           ! current yld locus size
        real,                 intent(in)  :: sigma_yld   ! first yield limit
        real,                 intent(out) :: F_mises     ! Mises yield locus
        real, dimension(0:5), intent(out) :: gradF_mises ! Mises yield locus gradient

        real, dimension(0:5) :: Sigma_ij_dev
        real                 :: tau_eq_mises
        integer :: k
        integer, dimension(0:5), parameter :: A = (/1.0,1.0,1.0,2.0,2.0,2.0/)
        
        call tensor_components (Sigma_ij, Sigma_ij_dev)

        call tau_mises(Sigma_ij_dev-X_ij, tau_eq_mises)

        F_mises      =  tau_eq_mises-sigma_yld-R
        gradF_mises(0:5)=0d0

        do k=0,5
            gradF_mises(k) = 1.5*A(k)*(Sigma_ij_dev(k)-X_ij(k))/tau_eq_mises
        end do
    end subroutine mises_yld_locus

    subroutine tau_mises(A_ij,tau_eq_mises)
        implicit none

        real, dimension(0:5), intent(in) :: A_ij
        real,                 intent(out):: tau_eq_mises
        integer                          :: i

        tau_eq_mises = 0.0d0
        do i = 0,5
            if (i .le. 2) then
                tau_eq_mises = tau_eq_mises + (A_ij(i))**2
            else
                tau_eq_mises = tau_eq_mises + 2*(A_ij(i))**2
            end if
        end do
        tau_eq_mises = sqrt(1.5*tau_eq_mises)

    end subroutine

    subroutine tensor_components (Sigma_ij, Sigma_ij_dev)

        ! TENSOR COMPONENTS (SPHERICAL AND DEVIATORIC)

        implicit none
        real, dimension(0:5), intent(in)  :: Sigma_ij
        real, dimension(0:5), intent(out) :: Sigma_ij_dev
        real                              :: Sigma_P
        integer :: i

        Sigma_ij_dev(0:5) = Sigma_ij(0:5)
        Sigma_P           = 0d0
        do i = 0,2
            Sigma_P=Sigma_P+Sigma_ij(i)/3
        end do

        do i=0,2
            Sigma_ij_dev(i) = Sigma_ij_dev(i)-Sigma_P
        end do

    end subroutine tensor_components

    subroutine check_plasticity (Sigma_ij_trial, Sigma_ij_start, X_ij, R, sigma_yld, &
        st_epl, alpha_elp,nx,ny,nz,nelement)

        ! CHECK PLASTIC CONSISTENCY (KKT CONDITIONS)

        implicit none

        real, dimension(0:5), intent(in)    :: Sigma_ij_start ! initial stress state
        real, dimension(0:5), intent(in)    :: X_ij           ! current back-stress state
        real,                 intent(in)    :: R              ! current yld locus size
        real,                 intent(in)    :: sigma_yld      ! first yield limit
        real, dimension(0:5), intent(inout) :: Sigma_ij_trial ! stress trial increment

        real,                 intent(out) :: alpha_elp      ! percentage of elastic strain
        integer,              intent(out) :: st_epl         ! elasto-plastic status
        integer, intent(in)               :: nx,ny,nz,nelement
        real                                :: F_start, F_final_trial,checkload
        real, dimension(0:5)                :: gradF_start, gradF_trial
        real, dimension(0:5)                :: dSigma_ij_trial
        integer                             :: GLLx,GLLy,GLLz,k
        ! Stress trial increment
        dSigma_ij_trial = Sigma_ij_trial-Sigma_ij_start

        ! Yield function at Sigma_ij_start
        call mises_yld_locus (Sigma_ij_start, X_ij, R, sigma_yld, &
            F_start, gradF_start)

        ! Yield function at Sigma_ij_trial
        call mises_yld_locus (Sigma_ij_trial, X_ij, R, sigma_yld, &
            F_final_trial, gradF_trial)
        
        if (nelement==4 .and. nx==4 .and. ny==4 .and. nz==4) then
            write(*,*) "F_start(NL)",F_start
            write(*,*) "F_final_trial(NL)",F_final_trial
        end if
        checkload=0d0
        do k=0,5
            checkload=checkload+10*gradF_start(k)*dSigma_ij_trial(k)
        end do
        ! KKT condition
        if (abs(F_start) .le. tol_nl) then
            if (checkload .ge. 0) then
                alpha_elp = 0d0
                st_epl    = 1
            else
                if (F_final_trial .gt. 0) then
                    alpha_elp = 2*sigma_yld/(2*sigma_yld+F_final_trial)
                    st_epl    = 1
                else
                    alpha_elp = 1d0
                    st_epl    = 3 
                endif
            end if
        elseif (F_start .lt. -tol_nl) then
            if (F_final_trial .lt. tol_nl) then
                alpha_elp = 1d0
                st_epl    = 2
            else
                alpha_elp = F_start/(F_start-F_final_trial)    
                st_epl    = 1
            end if
        elseif (F_start .gt. tol_nl) then
            write(*,*) "element: ", nelement
            write(*,*) "gll: ",nx,ny,nz
            stop
        end if

        Sigma_ij_trial = Sigma_ij_start+dSigma_ij_trial*alpha_elp
        if (nelement==4 .and. nx==4 .and. ny==4 .and. nz==4) then
            call mises_yld_locus (Sigma_ij_trial, X_ij, R, sigma_yld, &
                F_start, gradF_trial)
            write(*,*) "F_after(NLend)",F_start
        end if
        
    end subroutine check_plasticity

    subroutine plastic_corrector (dEpsilon_ij_alpha, Sigma_ij_1, X_ij_1, sigma_yld, &
        R_1, b_lmc, Rinf_lmc, C_lmc, kapa_lmc, mu, lambda, dEpsilon_ij_pl,Ffinal,condition)

        ! NR ALGORITHM AND DRIFT CORRECTION

        implicit none
        real, dimension(0:5), intent(inout) :: Sigma_ij_1           ! starting stress state
        real, dimension(0:5), intent(inout) :: X_ij_1               ! starting back stress
        real,                 intent(inout) :: R_1                  ! starting mises radius
        real, dimension(0:5), intent(inout) :: dEpsilon_ij_alpha    ! percentage of elastic-plastic strain
        real, dimension(0:5), intent(inout) :: dEpsilon_ij_pl       ! percentage of plastic strain
        real,                 intent(in)    :: sigma_yld            ! first yield limit
        real,                 intent(in)    :: b_lmc, Rinf_lmc      ! Lamaitre and Chaboche parameters (isotropic hardening)
        real,                 intent(in)    :: C_lmc, kapa_lmc      ! Lamaitre and Chaboche parameters (kinematic hardening)
        real,                 intent(in)    :: mu, lambda           ! elastic parameters

        real, dimension(0:5)                :: Sigma_ij_0, X_ij_0, dX_ij_0, gradF_0, gradF_mises
        real                                :: R_0, dR_0, dPlastMult_1, F_mises_0, F_mises
        integer                             :: i,j,k
        integer, parameter                  :: N_incr = 10
        real, dimension(0:5)                :: temp_vec
        real, dimension(0:2), parameter     :: vec0 = (/ 0.0, 0.0, 0.0 /)
        real, dimension(0:2), parameter     :: veci = (/ 1.0, 0.0, 0.0 /)
        real, dimension(0:2), parameter     :: vecj = (/ 0.0, 1.0, 0.0 /)
        real, dimension(0:2), parameter     :: veck = (/ 0.0, 0.0, 1.0 /)
        real, dimension(0:5,0:5)            :: DEL_ijhk
        real, dimension(0:2,0:2)            :: M
        real, dimension(0:2,0:2), parameter :: id_matrix = reshape( (/veci,vecj,veck/), (/3,3/) )
        real, dimension(0:5), parameter     :: A = (/1.0,1.0,1.0,0.5,0.5,0.5/)
        real, intent(out) :: Ffinal
        logical, intent(in) :: condition
        DEL_ijhk(:,:)   = 0d0
        M(0:2,0:2)      = 1
        DEL_ijhk(0:2,0:2) = DEL_ijhk(0:2,0:2)+lambda*M+id_matrix*2*mu
        DEL_ijhk(3:5,3:5) = DEL_ijhk(3:5,3:5)+id_matrix*mu
        dEpsilon_ij_alpha(0:5) = dEpsilon_ij_alpha(0:5)/N_incr

        do i = 0,N_incr-1
            write(*,*) "******** INCREMENT:",i," ***********"
            WRITE(*,*) ""
            ! PREDICTION
            Sigma_ij_0 = Sigma_ij_1
            X_ij_0     = X_ij_1
            R_0        = R_1
            if (condition) then
                write(*,*) "******** PREDICTION *********"
                write(*,*) ""
                write(*,*) "sigma0",sigma_ij_0
                write(*,*) "Xij0",X_ij_0
                write(*,*) "R0", R_0
                write(*,*) ""
            end if
            call mises_yld_locus (Sigma_ij_0, X_ij_0, R_0, sigma_yld, F_mises_0, gradF_0)

            ! PLASTIC MULTIPLIER
            call compute_plastic_modulus(dEpsilon_ij_alpha, Sigma_ij_0, X_ij_0, R_0, mu, lambda, sigma_yld, &
                b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dPlastMult_1,condition)

            ! HARDENING INCREMENTS
            call hardening_increments(dPlastMult_1, Sigma_ij_0, R_0, X_ij_0, sigma_yld, &
                b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR_0, dX_ij_0)
            
            do k=0,5
                dEpsilon_ij_pl(k)=dEpsilon_ij_pl(k)+dPlastMult_1*gradF_0(k)*A(k)
            end do
            R_1=R_0+dR_0          ! isotropic hardening update
            X_ij_1=X_ij_0+dX_ij_0 ! back-stress update
            
            if (condition) then
                write(*,*) "******** HARDENING INCREMENTS *********"
                write(*,*) ""
                write(*,*) "deltaR",dR_0
                write(*,*) "deltaX",dX_ij_0
                write(*,*) "dlamba",dPlastMult_1
                write(*,*) ""
            end if
            temp_vec(0:5) = 0d0
            do j = 0,5
                do k = 0,5
                   temp_vec(j)=temp_vec(j)+DEL_ijhk(j,k)*(dEpsilon_ij_alpha(k)-A(k)*dPlastMult_1*gradF_0(k))
                end do
            end do
            Sigma_ij_1=Sigma_ij_0+temp_vec
            call mises_yld_locus (Sigma_ij_1, X_ij_1, R_1, sigma_yld, Ffinal, gradF_0)

            WRITE(*,*) "********* BEFORE DRIFT **************"
            write(*,*) "Ffinal_temporary",Ffinal
            write(*,*) ""
            ! DRIFT CORRECTION
            call drift_corr(Sigma_ij_1, X_ij_1, R_1, sigma_yld,condition) ! drift correction (radial return)
            WRITE(*,*) "********* AFTER DRIFT **************"
            write(*,*) "sigma",Sigma_ij_1
            write(*,*) "Ffinal_temporary",Ffinal
            write(*,*) ""
            write(*,*) "*********** END INCREMENT",i,"*****************"
            
        end do
        if (condition) then
            write(*,*) "*********** FINAL STRESS **********"
            write(*,*) "sigma1",Sigma_ij_1
            WRITE(*,*) ""
        end if
        call mises_yld_locus (Sigma_ij_1, X_ij_1, R_1, sigma_yld, Ffinal, gradF_0)
    end subroutine plastic_corrector

    subroutine compute_plastic_modulus(dEpsilon_ij, Sigma_ij, X_ij, R, mu, lambda, &
        sigma_yld, b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dPlastMult, condition)

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
        real,                 intent(out):: dPlastMult          ! plastic multiplier increment
        real                             :: h_iso, h_kin        ! isotropic and kinematic hardening modula
        real                             :: F_mises
        real, dimension(0:5)             :: gradF_mises
        real                             :: h_lmc               ! total hardening modulus
        real                             :: temp_vec
        real, dimension(0:2), parameter  :: vec0 = (/ 0.0, 0.0, 0.0 /), &
            veci = (/ 1.0, 0.0, 0.0 /), &
            vecj = (/ 0.0, 1.0, 0.0 /), &
            veck = (/ 0.0, 0.0, 1.0 /)
        real, dimension(0:5,0:5)         :: DEL_ijhk
        real, dimension(0:2,0:2)         :: M
        real, dimension(0:2,0:2), parameter :: id_matrix = reshape( (/veci,vecj,veck/), (/3,3/) )
        integer                          :: j,k
        integer, dimension(0:5), parameter :: A =(/1.0,1.0,1.0,0.5,0.5,0.5/)
        logical, intent(in) :: condition
        DEL_ijhk(:,:) = 0d0
        M(0:2,0:2) = 1d0

        DEL_ijhk(0:2,0:2) = DEL_ijhk(0:2,0:2) + lambda * M + id_matrix *2*mu
        DEL_ijhk(3:5,3:5) = DEL_ijhk(3:5,3:5) + id_matrix * mu

        call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)

        h_iso = b_lmc*(Rinf_lmc-R)
        h_kin = C_lmc
        do k=0,5
            h_kin=h_kin-kapa_lmc*X_ij(k)*gradF_mises(k)
        end do
        h_lmc = h_iso + h_kin
        if (condition) then
            write(*,*) ""
            write(*,*) "***************** HARDENING ****************"
            write(*,*) "hardening modulus"
            write(*,*) h_kin, h_iso
            write(*,*) ""
        end if
        ! PLASTIC MULTIPLIER
        temp_vec   = 0d0
        dPlastMult = 0d0

        do k = 0,5
            do j = 0,5
                temp_vec   = temp_vec+gradF_mises(j)*DEL_ijhk(j,k)*gradF_mises(k)*A(k)
                dPlastMult = dPlastMult+gradF_mises(j)*DEL_ijhk(j,k)*dEpsilon_ij(k)
            end do
        end do

        dPlastMult = abs(dPlastMult/(h_lmc+temp_vec))
        if (dPlastMult .lt. 0) then
            write(*,*) ""
            write(*,*) "DPLAST NEGATIVE"
            write(*,*) "SIGMA", Sigma_ij
            write(*,*) "DEPSI",dEpsilon_ij
            stop
        end if
    end subroutine compute_plastic_modulus


    subroutine hardening_increments(dPlastMult, Sigma_ij, R, X_ij, sigma_yld, &
        b_lmc, Rinf_lmc, C_lmc, kapa_lmc, dR, dX_ij)

        ! INCREMENTS OF INTRINSIC STATIC VARIABLES

        real, dimension(0:5), intent(in) :: Sigma_ij        ! actual stress state
        real, dimension(0:5), intent(in) :: X_ij            ! actual back stress state
        real,                 intent(in) :: dPlastMult      ! plastic multiplier increment
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
        dR = dPlastMult*b_lmc*(Rinf_lmc-R)

        ! INCREMENT IN KINEMATIC HARDENING VARIABLES (Xij)
        call mises_yld_locus (Sigma_ij, X_ij, R, sigma_yld, F_mises, gradF_mises)
        dX_ij(0:5)=0d0
        do k=0,5
            dX_ij(k)=dPlastMult*(2*gradF_mises(k)*C_lmc/3 - X_ij(k)*kapa_lmc)
        end do

    end subroutine hardening_increments

    subroutine drift_corr(Sigma_ij_0, X_ij, R, sigma_yld,condition)

        ! DRIFT CORRECTION (RADIAL RETURN)

        real, dimension(0:5), intent(inout) :: Sigma_ij_0 ! actual stress state
        real, dimension(0:5), intent(in)    :: X_ij       ! actual back stress state
        real,                 intent(in)    :: R          ! actual mises radius
        real,                 intent(in)    :: sigma_yld
        real, dimension(0:5)                :: Sigma_ij_1
        real                                :: F_0, alpha0, F_1, alpha1, dalpha, err0, err1
        real, dimension(0:5)                :: gradF_0, Sigma_ij_dev_0
        !real, dimension(0:5)                :: gradF_1, Sigma_ij_dev_1
        real :: beta
        integer, dimension(0:5), parameter  :: A = (/1.0,1.0,1.0,0.5,0.5,0.5/)
        integer :: k,j
        logical, intent(in) :: condition
        call mises_yld_locus(Sigma_ij_0, X_ij, R, sigma_yld, F_0, gradF_0)
        if (F_0 .gt. tol_nl) then
            beta=0d0
            do k=0,5
                beta=beta+gradF_0(k)*gradF_0(k)
            end do
            beta=-F_0/beta
            do k=0,5
                Sigma_ij_0(k)=Sigma_ij_0(k)+beta*gradF_0(k)
            end do
        end if
        !alpha0 = 0d0 ! starting iteration value (step 0)
        ! NR algorithm for drift correction (radial return) (see Sloan 1987)
        !j=0;
        !do
        !    j=j+1;
        !    call tensor_components(Sigma_ij_0, Sigma_ij_dev_0)
        !    call mises_yld_locus(Sigma_ij_0, X_ij, R, sigma_yld, F_1, gradF_1)
        !    dalpha = -F_1/DOT_PRODUCT(gradF_1, gradF_0)
        !    alpha1 = alpha0 + dalpha
        !   do k=0,5
        !        Sigma_ij_1(k) = Sigma_ij_0(k)+alpha1*A(k)*gradF_0(k)
        !    end do
        !    call tensor_components (Sigma_ij_1, Sigma_ij_dev_1)
        !    call tau_mises(Sigma_ij_dev_1 - Sigma_ij_dev_0, err0)
        !    call tau_mises(Sigma_ij_dev_0 - X_ij, err1)
        !    if (condition) then
         !       write(*,*) "error0",err0
       !         write(*,*) "error1",err1
      !      end if
      !      err0 = err0/err1
     !      Sigma_ij_0 = Sigma_ij_1
    !        alpha0 = alpha1
   !         if (err0 .lt. 0.00001) then
  !              exit
 !           end if
!        end do

    end subroutine drift_corr

end module nonlinear

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
