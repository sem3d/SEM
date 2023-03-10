
!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Face.F90
!!\brief Gere les faces des elements.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module sfaces
    use constants

    implicit none
    ! Modified by Gaetano 01/06/05
    type :: face

       integer :: ngll, mat_index, type_Flux, Type_DG, mortarID
       integer, dimension (0:1) :: Near_Element, Which_face, Near_Vertex

       real(fpp), dimension (:), allocatable :: massMat
       real(fpp), dimension (:,:), allocatable :: Veloc, Displ, Accel, V0, Forces
       real(fpp), dimension (:,:), allocatable :: Veloc1, Veloc2, Forces1, Forces2
       real(fpp), dimension (:,:), allocatable :: DumpMass, DumpVx, DumpVz
       integer, dimension (:), allocatable :: Iglobnum_Face

#ifdef MKA3D
       real(fpp), dimension (:,:), allocatable :: ForcesMka
#endif

       logical :: coherency, PML, Abs, CPML, ADEPML, freesurf, reflex

       ! DG
       real(fpp), dimension (:), allocatable :: Normal
       real(fpp), dimension (:), allocatable :: k0,k1,Zp_p,Zp_m,Zs_p,Zs_m
       real(fpp), dimension (:), allocatable :: Mu_p, Mu_m, Lambda_p, Lambda_m, Rho_m, Rho_p
       real(fpp), dimension (:,:), allocatable :: Flux, Veloc_p,Veloc_m,Strain_p,Strain_m, Flux_p
       real(fpp), dimension (:,:), allocatable :: r1, r2, r3  ! EigenVectors for DG Godunov
       real(fpp), dimension (:,:), allocatable :: Vect_RK
       ! HDG
       real(fpp), dimension (:,:), allocatable :: Normal_Nodes
       real(fpp), dimension (:,:), allocatable :: Kinv, KinvExpl, SmbrTrac
       real(fpp), dimension (0:1)    :: Coeff_integr_ends
       integer, dimension (0:1) :: pos_in_VertMat
       logical :: is_computed, changing_media, acoustic, mortar

    end type face

contains


  ! ############################################################
  !>
  !! \brief Compute the flux for Element E at the face F
  !!  Several different fluxs are available :
  !!  Type_Flux = 1 : Centered Flux
  !!  Type_Flux = 2 : Godunov Flux
  !!
  !! \param type (Face), intent (INOUT) F
  !<

  subroutine Compute_Flux_DGstrong (F)
    implicit none

    type (Face), intent (INOUT)      :: F

    ! local variables
    integer                          :: i
    real(fpp), dimension (0:F%ngll-1,0:1) :: Stress_jump
    real(fpp), dimension (0:F%ngll-1,0:1) :: Veloc_jump
    real(fpp), dimension (0:F%ngll-1,0:4) :: F_minus
    real(fpp), dimension (0:F%ngll-1)     :: coeff_p

    ! --------- CENTERED FLUX -----------
    !          (partie a revoir)          !
    if (F%Type_Flux == FLUX_CENTERED) then
        F_minus = compute_trace_F(F,.true.)
        F%Flux = 0.5* (F_minus - compute_trace_F(F,.false.))
        if(F%Abs) then
            F%Flux = compute_trace_F(F,.true.)
        endif

    ! -------- GODUNOV FLUX ----------
    else if (F%Type_Flux == FLUX_GODUNOV) then
        Stress_jump = compute_stress_jump(F)
        Veloc_jump(:,:) = F%Veloc_m(:,:) - F%Veloc_p(:,:)
        call check_r1(F,.true.)
        if (F%freesurf) then
            Veloc_jump(:,:)  = 0.
            Stress_jump(:,:) = 2. * Stress_jump(:,:)
        endif
        coeff_p(:) = Stress_jump(:,0) * F%Normal(0) + Stress_jump(:,1) * F%Normal(1) &
                   + F%Zp_p(:) * (F%Normal(0)*Veloc_jump(:,0) + F%Normal(1)*Veloc_jump(:,1))
        if (.NOT. F%acoustic) then
            F%r2 = compute_r(F,Stress_jump,.true.)
            F%r3 = compute_r(F,Veloc_jump, .true.)
            do i=0,F%ngll-1
                F%Flux(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:) - F%k1(i)*F%r2(i,:) - F%Zs_p(i)*F%k1(i)*F%r3(i,:)
            enddo
        else ! Acoustic case
            do i=0,F%ngll-1
                F%Flux(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:)
            enddo
        endif

        ! Construction du flux du cote "plus" de la face
        if (.not. F%changing_media) then
            F%Flux_p(:,:) = -F%Flux(:,:) - compute_trace_F(F,.true.) - compute_trace_F(F,.false.)
        else
            ! Compute a second flux Face%Flux_p for the other side
            Veloc_jump(:,:) = - Veloc_jump(:,:)
            call check_r1(F,.false.)
            coeff_p(:) = -Stress_jump(:,0) * F%Normal(0) - Stress_jump(:,1) * F%Normal(1) &
                - F%Zp_m(:) * (F%Normal(0)*Veloc_jump(:,0) + F%Normal(1)*Veloc_jump(:,1))
            if (F%acoustic) then ! Plus side is elastic
                F%r2 = compute_r(F,Stress_jump,.false.)
                F%r3 = compute_r(F,Veloc_jump, .false.)
                do i=0,F%ngll-1
                    F%Flux_p(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:) - F%k1(i)*F%r2(i,:) - F%Zs_m(i)*F%k1(i)*F%r3(i,:)
                enddo
            else ! Plus side is acoustic
                do i=0,F%ngll-1
                    F%Flux_p(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:)
                enddo
            endif
        endif

        if (F%Reflex) then
            F%Flux(:,:) = compute_trace_F(F,.true.)
        endif

    endif

end subroutine Compute_Flux_DGstrong

  ! ############################################################
  !>
  !! \brief Compute the flux for Element E at the face F
  !!  Several different fluxs are available :
  !!  Type_Flux = 1 : Centered Flux
  !!  Type_Flux = 2 : Godunov Flux
  !!
  !! \param type (Face), intent (INOUT) F
  !<

  subroutine Compute_Flux_DGweak (F)
    implicit none

    type (Face), intent (INOUT)      :: F

    ! local variables
    integer                          :: i
    real(fpp), dimension (0:F%ngll-1,0:1) :: Stress_jump
    real(fpp), dimension (0:F%ngll-1,0:1) :: Veloc_jump
    real(fpp), dimension (0:F%ngll-1,0:4) :: F_minus
    real(fpp), dimension (0:F%ngll-1)     :: coeff_p

    ! --------- CENTERED FLUX -----------
    if (F%Type_Flux == FLUX_CENTERED) then
        F_minus = compute_trace_F(F,.true.)
        F%Flux = 0.5* (F_minus - compute_trace_F(F,.false.))
        if(F%Abs) then
            F%Flux = compute_trace_F(F,.true.)
        endif

    ! -------- GODUNOV FLUX ----------
    else if (F%Type_Flux == FLUX_GODUNOV) then
        Stress_jump = compute_stress_jump(F)
        Veloc_jump(:,:) = F%Veloc_m(:,:) - F%Veloc_p(:,:)
        call check_r1(F,.true.)
        if (F%freesurf) then
            Veloc_jump(:,:)  = 0.
            Stress_jump(:,:) = 2. * Stress_jump(:,:)
        endif
        coeff_p(:) = Stress_jump(:,0) * F%Normal(0) + Stress_jump(:,1) * F%Normal(1) &
                   + F%Zp_p(:) * (F%Normal(0)*Veloc_jump(:,0) + F%Normal(1)*Veloc_jump(:,1))
        if (.NOT. F%acoustic) then
            F%r2 = compute_r(F,Stress_jump,.true.)
            F%r3 = compute_r(F,Veloc_jump, .true.)
            do i=0,F%ngll-1
                F%Flux(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:) - F%k1(i)*F%r2(i,:) - F%Zs_p(i)*F%k1(i)*F%r3(i,:)
            enddo
        else ! Acoustic case
            do i=0,F%ngll-1
                F%Flux(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:)
            enddo
        endif
        F_minus = compute_trace_F(F,.true.)
        F%Flux(:,:) = F%Flux(:,:) + F_minus(:,:)

        ! Comptute a second flux Face%Flux_p if the interface is elastic-acoustic
        if(F%changing_media) then
            Veloc_jump(:,:) = - Veloc_jump(:,:)
            call check_r1(F,.false.)
            coeff_p(:) = -Stress_jump(:,0) * F%Normal(0) - Stress_jump(:,1) * F%Normal(1) &
                    - F%Zp_m(:) * (F%Normal(0)*Veloc_jump(:,0) + F%Normal(1)*Veloc_jump(:,1))
            if (F%acoustic) then ! Plus side is elastic
                F%r2 = compute_r(F,Stress_jump,.false.)
                F%r3 = compute_r(F,Veloc_jump, .false.)
                do i=0,F%ngll-1
                    F%Flux_p(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:) - F%k1(i)*F%r2(i,:) - F%Zs_m(i)*F%k1(i)*F%r3(i,:)
                enddo
            else ! Plus side is acoustic
                do i=0,F%ngll-1
                    F%Flux_p(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:)
                enddo
            endif
            F_minus = compute_trace_F(F,.false.)
            F%Flux_p(:,:) = F%Flux_p(:,:) + F_minus(:,:)
        endif

        if (F%Reflex) F%Flux(:,:) = 0.

    endif

    end subroutine Compute_Flux_DGweak


  ! ############################################################
  !>
  !! \brief Compute the term of flux for the interface of coupling CG-DG.
  !! The final term computed here is " sigma(n) - tau*(V-Vhat) ". The part
  !! "sigma(n) - tau*V" has been previously computed by subroutine compute_TracFace
  !! and is stored in Face%Traction. Here is just computed the part "tau*Vhat".
  !! It will be used later as a Neumann boundary condition for the Continuous
  !! Galerkin side.
  !! Be carefull : here Face%InvMatPen is the penalization matrix NOT inverted.
  !! \param type (Face), intent (INOUT) F
  !<
  subroutine Compute_Flux_Coupling (F)

      implicit none
      type (Face), intent (INOUT) :: F

      F%Forces(:,0) = F%Forces(:,0) - (F%SmbrTrac(:,0) + (F%KinvExpl(:,0)*F%Veloc(:,0) &
                                                       + F%KinvExpl(:,2)*F%Veloc(:,1)))
      F%Forces(:,1) = F%Forces(:,1) - (F%SmbrTrac(:,1) + (F%KinvExpl(:,2)*F%Veloc(:,0) &
                                                       + F%KinvExpl(:,1)*F%Veloc(:,1)))

      F%SmbrTrac(:,:) = 0.

      !################ DESACTIVATED IN MAY 2015 #####################!
      !F%Traction(:,0) =  F%Traction(:,0) + (F%InvMatPen(:,0)*F%Veloc(:,0) + F%InvMatPen(:,2)*F%Veloc(:,1))
      !F%Traction(:,1) =  F%Traction(:,1) + (F%InvMatPen(:,2)*F%Veloc(:,0) + F%InvMatPen(:,1)*F%Veloc(:,1))

      ! Adding The traction to the forces
      !F%Forces(:,0) = F%Forces(:,0) - F%Coeff_Integr(:) * F%Traction(:,0)
      !F%Forces(:,1) = F%Forces(:,1) - F%Coeff_Integr(:) * F%Traction(:,1)

      !F%Traction(:,:) = 0.

  end subroutine Compute_Flux_Coupling



    ! ############################################################
    !>
    !! \brief
    !!
    !! \param type (Face), intent (INOUT) F
    !! \param real(fpp), intent (IN) bega
    !! \param real(fpp), intent (IN) gam1
    !! \param real(fpp), intent (IN) dt
    !! \param real(fpp), intent (IN) alpha
    !<
    !  subroutine Prediction_Face_Veloc (F,alpha,bega, dt)
    subroutine Prediction_Face_Veloc (F)
        implicit none

        type (Face), intent (INOUT) :: F
        !real(fpp), intent (IN) :: bega,  alpha ,dt

        F%Forces = F%Displ
        F%V0 = F%Veloc

        return
    end subroutine Prediction_Face_Veloc


    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param integer, intent (IN) ngll
    !! \param type (Face), intent (INOUT) F
    !! \param real(fpp), intent (IN) bega
    !! \param real(fpp), intent (IN) gam1
    !! \param real(fpp), intent (IN) dt
    !<
    !  subroutine Correction_Face_Veloc (F, ngll, bega, gam1,dt)
    subroutine Correction_Face_Veloc (F, ngll, dt)
        implicit none

        integer, intent (IN) :: ngll
        type (Face), intent (INOUT) :: F
        !real(fpp), intent (IN) :: bega, gam1
        real(fpp), intent (IN) :: dt
        real(fpp), dimension (1:ngll-2) :: masse

        integer :: i

        do i=1,ngll-2
            masse(i) = F%MassMat(i)
        enddo

        do i = 0,1
            F%Forces(:,i) = masse(:)  * F%Forces(:,i)
        enddo

        if (F%Abs .or. F%reflex) F%Forces = 0

        F%Veloc  = F%v0+ dt * F%Forces
        F%Accel  = (F%Veloc-F%V0)/dt
        F%Displ  =  F%Displ + dt * F%Veloc
        return
    end subroutine Correction_Face_Veloc


    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Face), intent (INOUT) F
    !! \param real(fpp), intent (IN) dt
    !<
    subroutine Correction_Face_PML_Veloc (F, dt)
        implicit none

        type (Face), intent (INOUT) :: F
        real(fpp), intent (IN) ::  dt

        integer :: i

        F%V0 = F%Veloc

        if  (F%Abs .or. F%reflex) then
            F%Veloc1 = 0; F%Veloc2 = 0; F%Veloc = 0
        else

            do i = 0,1
                F%Veloc1(:,i) = F%DumpVx(:,0) * F%Veloc1(:,i) + dt * F%DumpVx(:,1)*F%Forces1(:,i)
                F%Veloc2(:,i) =F%DumpVz(:,0) * F%Veloc2(:,i) + dt * F%DumpVz(:,1)*F%Forces2(:,i)
            enddo
            F%Veloc = F%Veloc1 + F%Veloc2
        endif

        return
    end subroutine Correction_Face_PML_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Face), intent (INOUT) F
    !! \param real(fpp), intent (IN) dt
    !<
    subroutine Correction_Face_CPML_Veloc (F, dt)
        implicit none

        type (Face), intent (INOUT) :: F
        real(fpp), intent (IN) ::  dt
        integer :: i,ngll

        ngll = F%ngll
        F%V0 = F%Veloc

        if  (F%Abs .or. F%reflex) then
            F%Veloc = 0
        else
            do i = 0,1
                F%Forces(:,i) = F%MassMat(1:ngll-2)  * F%Forces(:,i)
            enddo
            F%Veloc  = F%v0+ dt * F%Forces
            F%Accel  = (F%Veloc-F%V0)/dt
            F%Displ  =  F%Displ + dt * F%Veloc
        endif

        return
    end subroutine Correction_Face_CPML_Veloc

    ! ###########################################################

    !>
    !! \brief
    !!
    !! \param integer, intent (IN) ngll
    !! \param type (Face), intent (IN) F
    !! \param real(fpp), dimension (1:ngll-2,0:1), intent (INOUT) Vfree
    !! \param logical, intent (IN) logic
    !! \param logical, intent (IN) logic2
    !<
    subroutine get_vfree_face (F,Vfree, ngll, logic,logic2)
        implicit none

        integer, intent (IN) :: ngll
        type (Face), intent (IN) :: F
        real(fpp), dimension (1:ngll-2,0:1), intent (INOUT) :: Vfree
        logical, intent (IN) :: logic , logic2
        real(fpp), dimension (1:ngll-2) :: masse
        integer :: i,j

        do i=1,ngll-2
            masse(i) = F%MassMat(i)
        enddo

        if (logic) then
            if (logic2) then
                do i = 0,1
                    Vfree(1:ngll-2,i) = Vfree(1:ngll-2,i) - masse(1:ngll-2)*F%Forces(1:ngll-2,i)
                enddo
            else
                do i = 0,1
                    do j = 1, ngll-2
                        Vfree(j,i) = Vfree(j,i) - masse(ngll-1-j)*F%Forces(ngll-1-j,i)
                    enddo
                enddo
            endif
        else
            do i = 0,1
                Vfree(1:ngll-2,i) =  masse(1:ngll-2)*F%Forces(1:ngll-2,i)
            enddo
        endif
        return
    end subroutine get_vfree_face


    ! ###########################################################

    !>
    !! \brief Compute the two las componants of the vector
    !! for the P-wave according to the side the current element
    !! is from the Face.
    !!
    !! \param type (Face), intent (IN) F
    !! \param logical, intent (IN) bool_side
    !<

    subroutine check_r1 (F,bool_side)
      implicit none

      type (Face), intent (INOUT) :: F
      logical, intent(IN)      :: bool_side

      if (bool_side) then
           F%r1(:,3) = F%Normal(0) * F%Zp_m(:)
           F%r1(:,4) = F%Normal(1) * F%Zp_m(:)
      else
           F%r1(:,3) = -F%Normal(0) * F%Zp_p(:)
           F%r1(:,4) = -F%Normal(1) * F%Zp_p(:)
      endif

    end subroutine check_r1
    ! ############################################################

    function compute_r(F,jump,bool_side)

      type (Face), intent (INOUT)                 :: F
      real(fpp), dimension(0:F%ngll-1,0:1), intent(IN) :: jump
      logical, intent(IN)             :: bool_side
      real(fpp), dimension(0:F%ngll-1,0:4) :: compute_r
      real(fpp), dimension(0:F%ngll-1,0:2) :: aux

      ! Computation of n x (n x jump)
      aux(:,2) = F%Normal(0)*jump(:,1) - F%Normal(1)*jump(:,0)
      aux(:,0) = F%Normal(1)*aux(:,2)
      aux(:,1) =-F%Normal(0)*aux(:,2)

      if (bool_side) then
         compute_r(:,0) = F%Normal(0) * aux(:,0)
         compute_r(:,1) = F%Normal(1) * aux(:,1)
         compute_r(:,2) = 0.5 * (F%Normal(1) * aux(:,0) + F%Normal(0) * aux(:,1))
         compute_r(:,3) = F%Zs_m(:) * aux(:,0)
         compute_r(:,4) = F%Zs_m(:) * aux(:,1)
      else
         compute_r(:,0) = -F%Normal(0) * aux(:,0)
         compute_r(:,1) = -F%Normal(1) * aux(:,1)
         compute_r(:,2) = -0.5 * (F%Normal(1) * aux(:,0) + F%Normal(0) * aux(:,1))
         compute_r(:,3) = F%Zs_p(:) * aux(:,0)
         compute_r(:,4) = F%Zs_p(:) * aux(:,1)
      endif

    end function compute_r

    ! ############################################################

    function compute_stress_jump(F)

      type (Face), intent (INOUT)     :: F
      real(fpp), dimension(0:F%ngll-1,0:1) :: compute_stress_jump
      real(fpp), dimension(0:F%ngll-1,0:2) :: sigma

      ! For the "minus" side
      sigma = compute_stress(F,.true.)
      compute_stress_jump(:,0) = sigma(:,0)*F%Normal(0) + sigma(:,2)*F%Normal(1)
      compute_stress_jump(:,1) = sigma(:,2)*F%Normal(0) + sigma(:,1)*F%Normal(1)

      ! For the "plus" side
      sigma = compute_stress(F,.false.)
      compute_stress_jump(:,0) = compute_stress_jump(:,0) - sigma(:,0)*F%Normal(0) - sigma(:,2)*F%Normal(1)
      compute_stress_jump(:,1) = compute_stress_jump(:,1) - sigma(:,2)*F%Normal(0) - sigma(:,1)*F%Normal(1)

    end function compute_stress_jump

    ! ############################################################

    function compute_stress(F,bool_side)

      type (Face), intent (INOUT)     :: F
      logical, intent(IN)             :: bool_side
      real(fpp), dimension(0:F%ngll-1,0:2) :: compute_stress
      real(fpp), dimension(0:F%ngll-1)     :: trace

      if(bool_side) then ! Current element on the "minus side"
         trace(:) = F%Lambda_m(:)*(F%Strain_m(:,0)+F%Strain_m(:,1))
         compute_stress(:,0) = 2*F%Mu_m(:) * F%Strain_m(:,0) + trace(:)
         compute_stress(:,1) = 2*F%Mu_m(:) * F%Strain_m(:,1) + trace(:)
         compute_stress(:,2) = 2*F%Mu_m(:) * F%Strain_m(:,2)
      else ! Current element on the "plus side"
         trace(:) = F%Lambda_p(:)*(F%Strain_p(:,0)+F%Strain_p(:,1))
         compute_stress(:,0) = 2*F%Mu_p(:) * F%Strain_p(:,0) + trace(:)
         compute_stress(:,1) = 2*F%Mu_p(:) * F%Strain_p(:,1) + trace(:)
         compute_stress(:,2) = 2*F%Mu_p(:) * F%Strain_p(:,2)
      endif
    end function compute_stress

    ! ############################################################

    function compute_trace_F(F,bool_side)

      type (Face), intent (INOUT)      :: F
      logical, intent(IN)              :: bool_side
      real(fpp), dimension(0:F%ngll-1,0:4)  :: compute_trace_F
      real(fpp), dimension(0:F%ngll-1,0:2)  :: sigma

      if (bool_side) then ! Case the current element is "minus"
         sigma = compute_stress(F,bool_side)
         compute_trace_F(:,0) = - F%Normal(0) * F%Veloc_m(:,0)
         compute_trace_F(:,1) = - F%Normal(1) * F%Veloc_m(:,1)
         compute_trace_F(:,2) = - 0.5* (F%Normal(1)*F%Veloc_m(:,0) + F%Normal(0)*F%Veloc_m(:,1))
         compute_trace_F(:,3) = - (sigma(:,0)*F%Normal(0) + sigma(:,2)*F%Normal(1))
         compute_trace_F(:,4) = - (sigma(:,2)*F%Normal(0) + sigma(:,1)*F%Normal(1))
      else ! Case the current element is on "plus" side
         sigma = compute_stress(F,bool_side)
         compute_trace_F(:,0) =  F%Normal(0) * F%Veloc_p(:,0)
         compute_trace_F(:,1) =  F%Normal(1) * F%Veloc_p(:,1)
         compute_trace_F(:,2) =  0.5* (F%Normal(1)*F%Veloc_p(:,0) + F%Normal(0)*F%Veloc_p(:,1))
         compute_trace_F(:,3) =  (sigma(:,0)*F%Normal(0) + sigma(:,2)*F%Normal(1))
         compute_trace_F(:,4) =  (sigma(:,2)*F%Normal(0) + sigma(:,1)*F%Normal(1))
      endif
    end function compute_trace_F

    ! ############################################################


    function Flux_Laurent(F,bool_side)

        type (Face), intent (INOUT)      :: F
        logical, intent(IN)              :: bool_side
        real(fpp), dimension(0:F%ngll-1,0:4)  :: Flux_Laurent
        real(fpp), dimension(0:F%ngll-1,0:2)  :: sigma, sigma_m, sigma_p
        real(fpp), dimension(0:F%ngll-1,0:1)  :: veloc

        veloc = 0.5 * (F%Veloc_m + F%Veloc_p)
        sigma_m = compute_stress(F,.true.)
        sigma_p = compute_stress(F,.false.)
        sigma = 0.5 * (sigma_m + sigma_p)

        Flux_Laurent(:,0) = F%normal(0) * veloc(:,0)
        Flux_Laurent(:,1) = F%normal(1) * veloc(:,1)
        Flux_Laurent(:,2) = 0.5 * (F%normal(1) * veloc(:,0) + F%normal(0) * veloc(:,1))
        Flux_Laurent(:,3) = sigma(:,0)*F%Normal(0) + sigma(:,2)*F%Normal(1)
        Flux_Laurent(:,4) = sigma(:,2)*F%Normal(0) + sigma(:,1)*F%Normal(1)

        if (bool_side) then ! Change normals signs
            Flux_Laurent = -Flux_Laurent
        endif

    end function Flux_Laurent


    ! ############################################################

    !>
    !! \brief This subroutine computes the kinetic energy of the interior nodes
    !!  of an face (all the nodes except the verices)
    !! \param type (Face), intent (INOUT) F
    !! \param real(fpp), intent (INOUT) E_kin
    !<

    subroutine  compute_Kinetic_Energy_F (F, Dt, E_kin)
        implicit none

        type (Face), intent (IN) :: F
        real(fpp), intent (IN)    :: Dt
        real(fpp), intent (INOUT) :: E_kin
        real(fpp), dimension (1:F%ngll-2)      :: Ener_Mat
        real(fpp), dimension (1:F%ngll-2, 0:1) :: Vel_half
        integer :: ngll

        ngll = F%ngll
        Vel_half(:,:) = F%Veloc(:,:) + 0.5 * dt * F%Forces(1:ngll-2,:)
        Ener_Mat (:)  = 1./F%MassMat(:) * ( Vel_half(:,0)*Vel_half(:,0) &
                                           +Vel_half(:,1)*Vel_half(:,1))
        E_kin = 0.5 * sum(Ener_Mat)

    end subroutine compute_Kinetic_Energy_F


    ! ###########################################################
    !>
    !! \brief subroutine Invert_F is used for HDG elements only
    !! This subroutine computes the inverse of the K matrix used to find
    !! the lagrange parameters (= velocities) lying on faces, solving for
    !! the matrix system K * Lambda = Smbr.
    !! This subroutine is used only for HDG in a semi-implicit framework.
    !! \param type (Face), intent (INOUT) F
    !<
    subroutine Invert_K_face (F)
        implicit none

        type (Face), intent (INOUT) :: F
        real(fpp), dimension(0:F%ngll-1) :: Det, tmp
        integer                     :: i

        Det(:) = F%Kinv(:,0) * F%Kinv(:,1) - (F%Kinv(:,2)*F%Kinv(:,2))
        ! Check positive-definiteness of matrices on Faces
        do i=0,F%ngll-1
            if ((F%Kinv(i,0) .LE. 0.) .OR. (Det(i) .LE. 0.)) then
                write(*,*) "Matrix Kinv not positive definite on current face and for node ", i
                STOP "Matrix should be sym def pos on faces. End of computation."
            endif
        enddo
        ! Compute inverse of matrices :
        F%Kinv(:,2) =-1./Det(:) * F%Kinv(:,2)
        tmp(:) = F%Kinv(:,1)
        F%Kinv(:,1) = 1./Det(:) * F%Kinv(:,0)
        F%Kinv(:,0) = 1./Det(:) * tmp(:)

    end subroutine Invert_K_face

    ! ###########################################################
    !>
    !! \brief subroutine Invert_K_face_Expl is used for HDG elements only
    !! This subroutine computes the inverse of the K matrix used to find
    !! the lagrange parameters (= velocities) lying on faces, solving for
    !! the matrix system K * Lambda = Smbr on an EXPLICIT TIME SCHEME.
    !! This subroutine is used only for HDG in a semi-implicit framework.
    !! \param type (Face), intent (INOUT) F
    !<
    subroutine Invert_K_face_Expl (F)
        implicit none

        type (Face), intent (INOUT) :: F
        real(fpp), dimension(0:F%ngll-1) :: Det, tmp
        integer                     :: i

        if (F%Type_DG==COUPLE_CG_HDG) then
            return
        endif

        if (F%abs) then
            F%KinvExpl = 2.*F%KinvExpl
        endif

        if (F%Acoustic) then
            F%KinvExpl(:,0) = 1./F%KinvExpl(:,0)
            return
        endif

        Det(:) = F%KinvExpl(:,0) * F%KinvExpl(:,1) - (F%KinvExpl(:,2)*F%KinvExpl(:,2))
        ! Check positive-definiteness of matrices on Faces
        do i=0,F%ngll-1
            if ((F%KinvExpl(i,0) .LE. 0.) .OR. (Det(i) .LE. 0.)) then
                write(*,*) "Matrix KinvExpl not positive definite on current face and for node ", i
                STOP "Matrix should be sym def pos on faces. End of computation."
            endif
        enddo
        ! Compute inverse of matrices :
        F%KinvExpl(:,2) =-1./Det(:) * F%KinvExpl(:,2)
        tmp(:) = F%KinvExpl(:,1)
        F%KinvExpl(:,1) = 1./Det(:) * F%KinvExpl(:,0)
        F%KinvExpl(:,0) = 1./Det(:) * tmp(:)

    end subroutine Invert_K_face_Expl

    ! ###########################################################
    !>
    !! \brief subroutine compute_Vhat_face is used for HDG elements only.
    !! This subroutine computes the lagrange parameters (= velocities) lying
    !! on faces, solving for the matrix system K * Lambda = Smbr.
    !! It uses the already-inverted matrix K --> Kinv.
    !! This subroutine is used only for HDG in a semi-implicit framework.
    !! \param type (Face), intent (INOUT) F
    !<
    subroutine compute_Vhat_face (F)
        implicit none

        type (Face), intent (INOUT) :: F

        ! La second membre "smbr" du systeme K * Lambda = Smbr est homgene aux tractions
        if (F%acoustic) then
            F%Veloc(:,0) = F%Kinv(:,0)*F%SmbrTrac(:,0)
        else
            F%Veloc(:,0) = F%Kinv(:,0)*F%SmbrTrac(:,0) + F%Kinv(:,2)*F%SmbrTrac(:,1)
            F%Veloc(:,1) = F%Kinv(:,2)*F%SmbrTrac(:,0) + F%Kinv(:,1)*F%SmbrTrac(:,1)
        endif

        F%SmbrTrac = 0.

        ! Treatment of boundary faces
        if (F%reflex) then
            F%Veloc(:,:) = 0.
        endif

    end subroutine Compute_Vhat_face

    ! ###########################################################
    !>
    !! \brief subroutine compute_Vhat_face is used for HDG elements only.
    !! This subroutine computes the lagrange parameters (= velocities) lying
    !! on faces, solving for the matrix system K * Lambda = Smbr.
    !! It uses the already-inverted matrix K --> Kinv.
    !! FOR EXPLICIT TRACE COMPUTATION !!!
    !! This subroutine is used only for HDG in a semi-implicit framework.
    !! \param type (Face), intent (INOUT) F
    !<
    subroutine compute_Vhat_face_Expl (F)
        implicit none

        type (Face), intent (INOUT) :: F

        ! La second membre "smbr" du systeme K * Lambda = Smbr est homgene aux tractions
        if (F%acoustic) then
            F%Veloc(:,0) = F%KinvExpl(:,0)*F%SmbrTrac(:,0)
        else
            F%Veloc(:,0) = F%KinvExpl(:,0)*F%SmbrTrac(:,0) + F%KinvExpl(:,2)*F%SmbrTrac(:,1)
            F%Veloc(:,1) = F%KinvExpl(:,2)*F%SmbrTrac(:,0) + F%KinvExpl(:,1)*F%SmbrTrac(:,1)
        endif

        F%SmbrTrac = 0.

        ! Treatment of boundary faces
        if (F%reflex) then
            F%Veloc(:,:) = 0.
        endif

    end subroutine Compute_Vhat_face_Expl


    ! ###########################################################

end module sfaces

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
