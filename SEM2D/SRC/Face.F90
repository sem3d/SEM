!>
!!\file Face.F90
!!\brief Gère les faces des éléments.
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

       integer :: ngll, mat_index, type_Flux, Type_DG
       integer, dimension (0:1) :: Near_Element, Which_face, Near_Vertex

       real, dimension (:), allocatable :: massMat
       real, dimension (:,:), allocatable :: Veloc, Displ, Accel, V0, Forces
       real, dimension (:,:), allocatable :: Veloc1, Veloc2, Forces1, Forces2
       real, dimension (:,:), allocatable :: DumpMass, DumpVx, DumpVz

#ifdef MKA3D
       real, dimension (:,:), allocatable :: ForcesMka
#endif

       logical :: coherency, PML, Abs, CPML, ADEPML, freesurf, reflex

       ! DG
       real, dimension (:), allocatable :: Normal
       real, dimension (:), allocatable :: k0,k1,Zp_p,Zp_m,Zs_p,Zs_m, Coeff_Integr
       real, dimension (:), allocatable :: Mu_p, Mu_m, Lambda_p, Lambda_m, Rho_m, Rho_p
       real, dimension (:,:), allocatable :: Flux, Veloc_p,Veloc_m,Strain_p,Strain_m, Flux_p
       real, dimension (:,:), allocatable :: r1, r2, r3  ! EigenVectors for DG Godunov
       real, dimension (:,:), allocatable :: Vect_RK
       ! HDG
       real, dimension (:,:), allocatable :: Normal_Nodes
       real, dimension (:,:), allocatable :: Kinv, Traction, Smbr, InvMatPen
       integer, dimension (0:1) :: pos_in_VertMat
       logical :: is_computed, changing_media, acoustic

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

  subroutine Compute_Flux (F,nelem,DG_type,acoustic)
    implicit none

    type (Face), intent (INOUT)      :: F
    integer, intent (IN)             :: nelem
    integer, intent (IN)             :: DG_type
    logical, intent (IN)             :: acoustic

    ! local variables
    integer                          :: i
    real, dimension (0:F%ngll-1,0:1) :: Stress_jump
    real, dimension (0:F%ngll-1,0:1) :: Veloc_jump
    real, dimension (0:F%ngll-1,0:4) :: F_minus
    real, dimension (0:F%ngll-1)     :: coeff_p
    logical                          :: bool_side

    if(nelem==F%Near_Element(0)) then
       bool_side = .TRUE.
    else
       bool_side = .FALSE.
    endif

    ! --------- CENTERED FLUX -----------
    if (F%Type_Flux == FLUX_CENTERED) then

        if (DG_type==GALERKIN_DG_WEAK) then
            F_minus = compute_trace_F(F,bool_side)
            F%Flux = 0.5* (F_minus - compute_trace_F(F,.NOT.bool_side))
            !F%Flux = 0.5 * (compute_trace_F(F,.true.) + compute_trace_F(F,.false.))
        elseif (DG_type==GALERKIN_DG_STRONG) then
            F_minus = compute_trace_F(F,bool_side)
            F%Flux = -0.5 * (F_minus + compute_trace_F(F,.NOT.bool_side))
        endif
        F%is_computed = .TRUE.
        !endif
        ! Treating Absorbing Boundary Conditions
        ! Basee surla supposition F* = Fn-
        if(F%Abs) then
            if (DG_type==GALERKIN_DG_WEAK) then
                F%Flux = compute_trace_F(F,bool_side)
            elseif (DG_type==GALERKIN_DG_STRONG) then
                F%Flux = 0.
            endif
        endif

        ! -------- CENTERED LAURENT FLUX ---------- !
    else if (F%Type_Flux == FLUX_CUSTOM_LG) then
        if (DG_type==GALERKIN_DG_WEAK) then
            F%Flux = Flux_Laurent(F,bool_side)
        elseif (DG_type==GALERKIN_DG_STRONG) then
            F%Flux = Flux_Laurent(F,bool_side) - compute_trace_F(F,bool_side)
        endif
        ! Treating Absorbing Boundary Conditions
        if(F%Abs) then
            if (DG_type==GALERKIN_DG_WEAK) then
                F%Flux = compute_trace_F(F,bool_side)
            elseif (DG_type==GALERKIN_DG_STRONG) then
                F%Flux = 0.
            endif
        endif

        ! -------- GODUNOV FLUX ----------
    else if (F%Type_Flux == FLUX_GODUNOV) then
        if(F%is_computed .AND. .NOT. F%changing_media) then
            ! Case the flux has been already computed
            if (DG_type == GALERKIN_DG_WEAK) then
                F%Flux = -F%Flux
            elseif (DG_type == GALERKIN_DG_STRONG) then
                F%Flux = -F%Flux - compute_trace_F(F,bool_side) - compute_trace_F(F,.NOT.bool_side)
            endif
            F%is_computed = .FALSE.
        else
            ! Computation of jumps of stresses and Velocities :
            Stress_jump = compute_stress_jump(F)
            Veloc_jump(:,:) = F%Veloc_m(:,:) - F%Veloc_p(:,:)
            if (.NOT. bool_side) Veloc_jump(:,:) = - Veloc_jump(:,:)
            call check_r1(F,bool_side)
            if (F%freesurf) then
                ! Jumps for Treating Free-Surface case
                Veloc_jump(:,:)  = 0.
                Stress_jump(:,:) = 2. * Stress_jump(:,:)
            endif
            if (bool_side) then
                coeff_p(:) = Stress_jump(:,0) * F%Normal(0) + Stress_jump(:,1) * F%Normal(1) &
                    + F%Zp_p(:) * (F%Normal(0)*Veloc_jump(:,0) + F%Normal(1)*Veloc_jump(:,1))
            else
                coeff_p(:) = -Stress_jump(:,0) * F%Normal(0) - Stress_jump(:,1) * F%Normal(1) &
                    - F%Zp_m(:) * (F%Normal(0)*Veloc_jump(:,0) + F%Normal(1)*Veloc_jump(:,1))
            endif
            if (.NOT. acoustic) then
                ! Computation of eigenvectors r2 and r3 :
                F%r2 = compute_r(F,Stress_jump,bool_side)
                F%r3 = compute_r(F,Veloc_jump, bool_side)
                ! Computation of Numerical Flux :
                if(bool_side) then
                    do i=0,F%ngll-1
                        F%Flux(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:) - F%k1(i)*F%r2(i,:) - F%Zs_p(i)*F%k1(i)*F%r3(i,:)
                    enddo
                else
                    do i=0,F%ngll-1
                        F%Flux(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:) - F%k1(i)*F%r2(i,:) - F%Zs_m(i)*F%k1(i)*F%r3(i,:)
                    enddo
                endif
            else ! Acoustic case
                do i=0,F%ngll-1
                    F%Flux(i,:) = coeff_p(i)*F%k0(i)*F%r1(i,:)
                enddo
            endif
            if (DG_type==GALERKIN_DG_WEAK) then ! Forme "Faible" des DG
                F_minus = compute_trace_F(F,bool_side)
                F%Flux(:,:) = F%Flux(:,:) + F_minus(:,:)
            endif
            F%is_computed = .TRUE.
            if (F%freesurf) F%is_computed = .FALSE.
            if (F%Abs) F%is_computed = .FALSE.
            if (F%Reflex) then
                F%is_computed = .FALSE.
                if (DG_type==GALERKIN_DG_WEAK)   F%Flux(:,:) = 0.
                if (DG_type==GALERKIN_DG_STRONG) F%Flux(:,:) = F_minus(:,:)
            endif
        endif
    endif
end subroutine Compute_Flux

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
    real, dimension (0:F%ngll-1,0:1) :: Stress_jump
    real, dimension (0:F%ngll-1,0:1) :: Veloc_jump
    real, dimension (0:F%ngll-1,0:4) :: F_minus
    real, dimension (0:F%ngll-1)     :: coeff_p

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
  !! \brief Compute the traces of velocity on the Face. Stored in F%Veloc.
  !! These velocities are to be understood in a HDG framework
  !! and are the unknowns at the interface.
  !! These velocities are computec multiplying the tractions F%TracFace
  !! by the inverse penalization matrix F%InvMatPen.
  !! \param type (Face), intent (INOUT) F
  !<
  subroutine Compute_Vhat (F)

      implicit none
      type (Face), intent (INOUT) :: F

      F%Veloc(:,0) = - F%InvMatPen(:,0)*F%Traction(:,0) - F%InvMatPen(:,2)*F%Traction(:,1)
      F%Veloc(:,1) = - F%InvMatPen(:,2)*F%Traction(:,0) - F%InvMatPen(:,1)*F%Traction(:,1)
      F%Traction = 0.

      ! Treatment of boundary faces
      if (F%reflex) F%Veloc(:,:) = 0.

  end subroutine Compute_Vhat


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

      F%Traction(:,0) =  F%Traction(:,0) + (F%InvMatPen(:,0)*F%Veloc(:,0) + F%InvMatPen(:,2)*F%Veloc(:,1))
      F%Traction(:,1) =  F%Traction(:,1) + (F%InvMatPen(:,2)*F%Veloc(:,0) + F%InvMatPen(:,1)*F%Veloc(:,1))

      ! Adding The traction to the forces
      F%Forces(:,0) = F%Forces(:,0) - F%Coeff_Integr(:) * F%Traction(:,0)
      F%Forces(:,1) = F%Forces(:,1) - F%Coeff_Integr(:) * F%Traction(:,1)

      F%Traction(:,:) = 0.

  end subroutine Compute_Flux_Coupling



    ! ############################################################
    !>
    !! \brief
    !!
    !! \param type (Face), intent (INOUT) F
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !<
    !  subroutine Prediction_Face_Veloc (F,alpha,bega, dt)
    subroutine Prediction_Face_Veloc (F)
        implicit none

        type (Face), intent (INOUT) :: F
        !real, intent (IN) :: bega,  alpha ,dt

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
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !<
    !  subroutine Correction_Face_Veloc (F, ngll, bega, gam1,dt)
    subroutine Correction_Face_Veloc (F, ngll, dt)
        implicit none

        integer, intent (IN) :: ngll
        type (Face), intent (INOUT) :: F
        !real, intent (IN) :: bega, gam1
        real, intent (IN) :: dt
        real, dimension (1:ngll-2) :: masse

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
    !! \param real, intent (IN) dt
    !<
    subroutine Correction_Face_PML_Veloc (F, dt)
        implicit none

        type (Face), intent (INOUT) :: F
        real, intent (IN) ::  dt

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
    !! \param real, intent (IN) dt
    !<
    subroutine Correction_Face_CPML_Veloc (F, dt)
        implicit none

        type (Face), intent (INOUT) :: F
        real, intent (IN) ::  dt
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
    !! \param real, dimension (1:ngll-2,0:1), intent (INOUT) Vfree
    !! \param logical, intent (IN) logic
    !! \param logical, intent (IN) logic2
    !<
    subroutine get_vfree_face (F,Vfree, ngll, logic,logic2)
        implicit none

        integer, intent (IN) :: ngll
        type (Face), intent (IN) :: F
        real, dimension (1:ngll-2,0:1), intent (INOUT) :: Vfree
        logical, intent (IN) :: logic , logic2
        real, dimension (1:ngll-2) :: masse
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
      real, dimension(0:F%ngll-1,0:1), intent(IN) :: jump
      logical, intent(IN)             :: bool_side
      real, dimension(0:F%ngll-1,0:4) :: compute_r
      real, dimension(0:F%ngll-1,0:2) :: aux

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
      real, dimension(0:F%ngll-1,0:1) :: compute_stress_jump
      real, dimension(0:F%ngll-1,0:2) :: sigma

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
      real, dimension(0:F%ngll-1,0:2) :: compute_stress
      real, dimension(0:F%ngll-1)     :: trace

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
      real, dimension(0:F%ngll-1,0:4)  :: compute_trace_F
      real, dimension(0:F%ngll-1,0:2)  :: sigma

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
        real, dimension(0:F%ngll-1,0:4)  :: Flux_Laurent
        real, dimension(0:F%ngll-1,0:2)  :: sigma, sigma_m, sigma_p
        real, dimension(0:F%ngll-1,0:1)  :: veloc

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
    !! \param real, intent (INOUT) E_kin
    !<

    subroutine  compute_Kinetic_Energy_F (F, Dt, E_kin)
        implicit none

        type (Face), intent (IN) :: F
        real, intent (IN)    :: Dt
        real, intent (INOUT) :: E_kin
        real, dimension (1:F%ngll-2)      :: Ener_Mat
        real, dimension (1:F%ngll-2, 0:1) :: Vel_half
        integer :: ngll

        ngll = F%ngll
        Vel_half(:,:) = F%Veloc(:,:) + 0.5 * dt * F%Forces(1:ngll-2,:)
        Ener_Mat (:)  = 1./F%MassMat(:) * ( Vel_half(:,0)*Vel_half(:,0) &
                                           +Vel_half(:,1)*Vel_half(:,1))
        E_kin = 0.5 * sum(Ener_Mat)

    end subroutine compute_Kinetic_Energy_F

    ! ###########################################################
    !>
    !! \brief subroutine compute_InvMatPen is used for HDG elements only
    !! This subroutine computes the inverse of penalty matrix used on the
    !! traces for HDG methods.
    !! \param type (Face), intent (INOUT) F
    !<
    subroutine compute_InvMatPen (F)
        implicit none

        type (Face), intent (INOUT) :: F
        real, dimension(0:F%ngll-1) :: invDet, Zp_m, Zp_p, Zs_m, Zs_p
        logical                     :: is_acoustic

        Zp_m(:) = sqrt(F%Rho_m(:) * (F%Lambda_m(:) + 2.* F%Mu_m(:)))
        Zp_p(:) = sqrt(F%Rho_p(:) * (F%Lambda_p(:) + 2.* F%Mu_p(:)))
        Zs_m(:) = sqrt(F%Rho_m(:) * F%Mu_m(:))
        Zs_p(:) = sqrt(F%Rho_p(:) * F%Mu_p(:))

        if ((F%Mu_p(0) == 0.) .and. (F%Mu_m(0) == 0.)) then
           is_acoustic = .true.
        else
           is_acoustic = .false.
        endif

        if (.not. is_acoustic) then
           invDet(:) = 1. / ((Zp_m(:)+Zp_p(:))*(Zs_m(:)+Zs_p(:)))
           F%InvMatPen(:,0) = invDet(:) * ((Zs_m(:)+Zs_p(:)) * F%Normal_nodes(:,0)**2 &
                                          +(Zp_m(:)+Zp_p(:)) * F%Normal_nodes(:,1)**2)
           F%InvMatPen(:,1) = invDet(:) * ((Zp_m(:)+Zp_p(:)) * F%Normal_nodes(:,0)**2 &
                                          +(Zs_m(:)+Zs_p(:)) * F%Normal_nodes(:,1)**2)
           F%InvMatPen(:,2) =-invDet(:) * ((Zp_m(:)+Zp_p(:)) - (Zs_m(:)+Zs_p(:))) &
                                        * F%Normal_nodes(:,0) * F%Normal_nodes(:,1)
        else ! acoustic case
           invDet(:) = 1. / (Zp_m(:)+Zp_p(:))
           F%InvMatPen(:,0) = invDet(:)
           F%InvMatPen(:,1) = invDet(:)
           F%InvMatPen(:,2) = 0.
        endif

        ! If the face is couplic CG with HDG, the matrix InvMatPen is NOT inverted :
        if (F%Type_DG == COUPLE_CG_HDG) then
            F%InvMatPen(:,0) = (Zp_m(:)+Zp_p(:)) * F%Normal_nodes(:,0)**2 &
                             + (Zs_m(:)+Zs_p(:)) * F%Normal_nodes(:,1)**2
            F%InvMatPen(:,1) = (Zs_m(:)+Zs_p(:)) * F%Normal_nodes(:,0)**2 &
                             + (Zp_m(:)+Zp_p(:)) * F%Normal_nodes(:,1)**2
            F%InvMatPen(:,2) =((Zp_m(:)+Zp_p(:)) - (Zs_m(:)+Zs_p(:))) &
                             * F%Normal_nodes(:,0) * F%Normal_nodes(:,1)
        endif

    end subroutine compute_InvMatPen

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
        real, dimension(0:F%ngll-1) :: Det, tmp
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
        F%Veloc(:,0) = F%Kinv(:,0)*F%Smbr(:,0) + F%Kinv(:,2)*F%Smbr(:,1)
        F%Veloc(:,1) = F%Kinv(:,2)*F%Smbr(:,0) + F%Kinv(:,1)*F%Smbr(:,1)

        F%Smbr = 0.

        ! Treatment of boundary faces
        if (F%reflex) then
            F%Veloc(:,:) = 0.
        endif

    end subroutine Compute_Vhat_face


    ! ###########################################################

end module sfaces
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
