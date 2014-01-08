!>
!!\file Face.F90
!!\brief Gère les faces des éléments.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module sfaces

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

       logical :: coherency, PML, Abs, FPML
       real, dimension (:), allocatable :: Ivx, Ivz
       real, dimension (:,:), allocatable :: Iveloc1, Iveloc2

       ! DG
       real, dimension (:), allocatable :: Normal
       real, dimension (:), allocatable :: k0,k1,Zp_p,Zp_m,Zs_p,Zs_m
       real, dimension (:), allocatable :: Mu_p, Mu_m, Lambda_p, Lambda_m
       real, dimension (:,:), allocatable :: Flux, Veloc_p,Veloc_m,Strain_p,Strain_m
       real, dimension (:,:), allocatable :: r1, r2, r3  ! EigenVectors for DG Godunov
       real, dimension (:,:), allocatable :: Vect_RK
       real, dimension (:,:), allocatable :: Normal_Nodes
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

  subroutine Compute_Flux (F,nelem,DG_type)
    implicit none

    type (Face), intent (INOUT)      :: F
    integer, intent (IN)             :: nelem
    integer, intent (IN)             :: DG_type

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
    if (F%Type_Flux == 1) then

       if(F%is_computed) then
          ! Case the flux has been already computed
          if(DG_type==1) then
             F%Flux = -F%Flux
          elseif (DG_type==0) then
             F%Flux = F%Flux
          endif
          F%is_computed = .FALSE.
       else
          if (DG_type==1) then
             F_minus = compute_trace_F(F,bool_side)
             F%Flux = 0.5* (F_minus - compute_trace_F(F,.NOT.bool_side))
          elseif (DG_type==0) then
             F_minus = compute_trace_F(F,bool_side)
             F%Flux = -0.5 * (F_minus + compute_trace_F(F,.NOT.bool_side))
          endif
          F%is_computed = .TRUE.
       endif
       ! Treating Absorbing Boundary Conditions
       ! Basee surla supposition F* = Fn-
       if(F%Abs) then
          if (DG_type==1) then
             F%Flux = compute_trace_F(F,bool_side)
          elseif (DG_type==0) then
             F%Flux = 0.
          endif
       endif

    ! -------- GODUNOV FLUX ----------
    else if (F%Type_Flux == 2) then
       !if(F%is_computed .AND. .NOT. F%changing_media) then
          ! Case the flux has been already computed
          !if (DG_type ==1 ) then
          !   F%Flux = -F%Flux
          !elseif (DG_type == 0) then
          !   F%Flux = -F%Flux - compute_trace_F(F,bool_side) - compute_trace_F(F,.NOT.bool_side)
          !endif
          !F%is_computed = .FALSE.
       !else
          ! Computation of jumps of stresses and Velocities :
          Stress_jump = compute_stress_jump(F)
          Veloc_jump(:,:) = F%Veloc_m(:,:) - F%Veloc_p(:,:)
          if (.NOT. bool_side) Veloc_jump(:,:) = - Veloc_jump(:,:)
          call check_r1(F,bool_side)
          if (F%Near_Element(1) .LT. 0) then
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
          if (.NOT. F%acoustic) then
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
          if (DG_type==1) then ! Forme "Faible" des DG
             F_minus = compute_trace_F(F,bool_side)
             F%Flux(:,:) = F%Flux(:,:) + F_minus(:,:)
          endif
          F%is_computed = .TRUE.
       !endif
       ! Treating Absorbing Boundary Conditions
       ! Basee surla supposition F* = Fn-
       if(F%Abs) then
          if (DG_type==1) then
             F%Flux = compute_trace_F(F,bool_side)
          elseif (DG_type==0) then
             F%Flux = 0.
          endif
       endif
    endif
  end subroutine Compute_Flux

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

        !  test mariotti
        !     F%Forces = F%Displ + dt * F%Veloc + dt**2 * (0.5 - bega) * F%Accel
        !     F%V0 = F%Veloc
        !     F%Forces = alpha * F%Forces + (1-alpha) * F%Displ
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

        !   ancienne version
        !      do i = 0,1
        !          F%Forces(:,i) = F%MassMat(:)  * F%Forces(:,i)
        !      enddo
        !  fin test modif
        if (F%Abs) F%Forces = 0
        !  test mariotti
        !      F%Veloc  = F%v0+ dt * F%Forces
        !      F%Accel  = F%Accel + gam1 /dt * (F%Veloc-F%V0)
        !      F%Displ  =  F%Displ + bega * dt * (F%Veloc+F%V0)
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

        if  (F%Abs) then
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
    !! \param type (Face toto), intent (INOUT) F
    !! \param real, intent (IN toto) dt
    !! \param real, intent (IN toto) fil
    !<


    subroutine Correction_Face_FPML_Veloc (F, dt, fil)
        implicit none

        type (Face), intent (INOUT) :: F
        real, intent (IN) ::  dt, fil

        integer :: i
        real :: fil2
        real, dimension (1:F%ngll-2) :: Ausiliar_Velocity


        fil2 = fil**2
        F%V0 = F%Veloc

        if  (F%Abs) then
            F%Veloc1 = 0; F%Veloc2 = 0; F%Veloc = 0
        else

            do i = 0,1
                Ausiliar_Velocity = F%Veloc1(:,i)
                F%Veloc1(:,i) = F%DumpVx(:,0) * F%Veloc1(:,i) + dt * &
                    F%DumpVx(:,1)*F%Forces1(:,i) + F%Ivx * F%Iveloc1(:,i)
                F%Iveloc1(:,i) = Fil2 * F%Iveloc1(:,i) + 0.5 * (1.-Fil2) *  &
                    (Ausiliar_Velocity + F%Veloc1(:,i) )

                Ausiliar_Velocity = F%Veloc2(:,i)
                F%Veloc2(:,i) =F%DumpVz(:,0) * F%Veloc2(:,i) + Dt * &
                    F%DumpVz(:,1)*F%Forces2(:,i) + F%Ivz * F%IVeloc2(:,i)
                F%Iveloc2(:,i) = Fil2 * F%Iveloc2(:,i) + 0.5 * (1.-Fil2) * &
                    (Ausiliar_Velocity + F%Veloc2(:,i) )
            enddo
            F%Veloc = F%Veloc1 + F%Veloc2
        endif

        return
    end subroutine Correction_Face_FPML_Veloc
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
    !! \param logical, intent (IN) bool_size
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
      real, dimension(0:F%ngll-1)     :: trace

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

end module sfaces
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
