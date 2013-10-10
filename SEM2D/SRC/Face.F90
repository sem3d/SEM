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

       integer :: ngll, mat_index
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
    end type face

contains

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
    ! ############################################################

end module sfaces
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
