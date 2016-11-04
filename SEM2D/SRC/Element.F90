!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Element.F90
!!\brief contient les methodes qui assure la gestion du type Element.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module selement
    use constants
    implicit none


    type :: element

       integer :: mat_index,ngllx,ngllz
       integer, dimension (:), allocatable :: Control_nodes
       integer, dimension (0:3) :: Near_Face   !
       integer, dimension (0:3) :: Near_Vertex !
       integer, dimension (:,:), allocatable :: Iglobnum

       real, dimension (:,:), allocatable :: Jacob,Density, Lambda, Mu,MassMat
       real, dimension(:,:,:), allocatable :: Forces,Stress,Veloc,Displ,Accel,V0
       real, dimension(:,:,:), allocatable :: ACoeff
       real, dimension(:,:,:,:), allocatable :: InvGrad
       logical :: OUTPUT, is_source

       ! PML allocation
       logical :: PML
       real, dimension (:,:,:), allocatable :: Stress1,Stress2
       real, dimension (:,:,:), allocatable :: Veloc1, Veloc2,Forces1,Forces2
       real, dimension (:,:,:), allocatable :: DumpMass
       real, dimension (:,:,:), allocatable :: DumpSx,DumpSz,DumpVx,DumpVz

       ! CPML/ADEPML allocation
       logical :: CPML, ADEPML
       real, dimension (:,:), allocatable :: Ax, Bx, Az, Bz
       real, dimension (:,:), allocatable :: PsiVxx, PsiVxz, PsiVzx, PsiVzz
       real, dimension (:,:), allocatable :: PsiSxxx, PsiSzzz, PsiSxzx, PsiSxzz

       real dist_max !!Ajout Gsa 03/10 - taille caracteristique de l'element

       ! DG
       integer :: type_DG
       logical :: acoustic
       real, dimension (:),    allocatable :: Coeff_Integr_Faces
       real, dimension(:,:,:), allocatable :: Strain, Strain0
       real, dimension(:,:,:), allocatable :: Vect_RK, Psi_store
       ! HDG
       real, dimension(:,:), allocatable :: Normal_nodes
       real, dimension(:,:), allocatable :: MatPen, TracFace, Vhat, Dinv
       real, dimension(:,:,:), allocatable :: CAinv, EDinv, Jmat
       integer, dimension (0:3,0:1) :: pos_corner_in_VertMat

    end type element

contains

    subroutine init_element(el)
        implicit none
        type(Element), intent(INOUT) :: el

        el%mat_index=-1
        el%PML = .false.
        el%CPML = .false.
        el%ADEPML = .false.
        el%acoustic = .false.
        el%OUTPUT = .true.
        el%dist_max = 0.0
        el%type_DG = GALERKIN_CONT

    end subroutine init_element

    ! ############################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !<
    !subroutine Prediction_Elem_Veloc (Elem,alpha,bega, dt)
    subroutine Prediction_Elem_Veloc (Elem)
        implicit none

        type (Element), intent (INOUT) :: Elem
        !real, intent (IN) :: bega, alpha, dt

        integer :: ngllx, ngllz

        ngllx = Elem%ngllx ; ngllz = Elem%ngllz

        Elem%Forces(1:ngllx-2,1:ngllz-2,0:1) = Elem%Displ
        Elem%V0 = Elem%Veloc

        return
    end subroutine Prediction_Elem_Veloc


    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !<
    !  subroutine Correction_Elem_Veloc (Elem, bega, gam1,dt)
    subroutine Correction_Elem_Veloc (Elem, dt)
        implicit none

        type (Element), intent (INOUT) :: Elem
        !real, intent (IN) :: bega, gam1
        real, intent (IN) :: dt
        integer :: ngllx, ngllz,i

        ngllx = Elem%ngllx ; ngllz = Elem%ngllz

        do i = 0,1
            Elem%Forces(1:ngllx-2,1:ngllz-2,i) = Elem%MassMat(:,:)  * Elem%Forces(1:ngllx-2,1:ngllz-2,i)
        enddo
        Elem%Veloc(:,:,:)  = Elem%v0(:,:,:)+ dt * Elem%Forces(1:ngllx-2,1:ngllz-2,:)
        Elem%Accel  = (Elem%Veloc-Elem%V0)/dt
        Elem%Displ = Elem%Displ + dt * Elem%Veloc

        return
    end subroutine Correction_Elem_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) dt
    !<


    subroutine Correction_Elem_PML_Veloc (Elem, dt)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (IN) ::  dt
        integer :: ngllx, ngllz,i

        ngllx = Elem%ngllx; ngllz=Elem%ngllz

        do i = 0,1
            Elem%Veloc1(:,:,i) = Elem%DumpVx(:,:,0) * Elem%Veloc1(:,:,i) + dt * &
                Elem%DumpVx(:,:,1)*Elem%Forces1(1:ngllx-2,1:ngllz-2,i)
            Elem%Veloc2(:,:,i) =Elem%DumpVz(:,:,0) * Elem%Veloc2(:,:,i) + Dt * &
                Elem%DumpVz(:,:,1)*Elem%Forces2(1:ngllx-2,1:ngllz-2,i)
        enddo
        Elem%V0 = Elem%Veloc
        Elem%Veloc = Elem%Veloc1 + Elem%Veloc2

        return
    end subroutine Correction_Elem_PML_Veloc


    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTmat
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hmatz
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vxloc
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vzloc
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !<
    subroutine Prediction_Elem_PML_Veloc (Elem,alpha, bega, dt,Vxloc,Vzloc,Hmatz, HTmat)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) ::  hTmat
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hmatz
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT)  ::Vxloc, Vzloc
        real, intent (IN) :: bega, dt, alpha

        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1) :: s0,s1,s2,s3

        integer :: ngllx, ngllz

        ngllx = Elem%ngllx; ngllz = Elem%ngllz

        !Elem%Stress0 = Elem%Stress

        VxLoc(1:ngllx-2,1:ngllz-2)  = (0.5+alpha) * Elem%Veloc(:,:,0) + dt *(0.5-bega)*Elem%Accel(:,:,0) &
                                    + (0.5-alpha) * Elem%V0(:,:,0)
        VzLoc(1:ngllx-2,1:ngllz-2)  = (0.5+alpha) * Elem%Veloc(:,:,1) + dt *(0.5-bega)*Elem%Accel(:,:,1) &
                                    + (0.5-alpha) * Elem%V0(:,:,1)

        s0 = MATMUL (HTmat,VxLoc)
        s2 = MATMUL (HTmat,VzLoc)
        s1 = MATMUL (VxLoc,Hmatz)
        s3 = MATMUL (VzLoc,Hmatz)

        Elem%Stress1(:,:,0) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,0) &
                            + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,0)*s0 + Elem%Acoeff(:,:,2)*s2)
        Elem%Stress2(:,:,0) = Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,0) &
                            + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,3)*s3 + Elem%Acoeff(:,:,1)*s1)

        Elem%Stress1(:,:,1) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,1) &
                            + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,4)*s0 + Elem%Acoeff(:,:,6)*s2)
        Elem%Stress2(:,:,1) = Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,1) &
                            + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,7)*s3 + Elem%Acoeff(:,:,5)*s1)

        Elem%Stress1(:,:,2) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,2) &
                            + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,10)*s2 + Elem%Acoeff(:,:,8)*s0)
        Elem%Stress2(:,:,2) = Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,2) &
                            + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,9)*s1 + Elem%Acoeff(:,:,11)*s3)

        Elem%Stress = Elem%Stress1 + Elem%Stress2
        !Elem%Stress = alpha  * Elem%Stress + (1-alpha) * Elem%Stress0
        return
    end subroutine Prediction_Elem_PML_Veloc


    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTmat
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hmatz
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vxloc
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vzloc
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !! \param real, intent (IN) fil
    !<
    subroutine Prediction_Elem_CPML_Veloc (Elem, alpha, bega, dt, Vxloc, Vzloc, Hmatz, HTmat)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) ::  hTmat
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hmatz
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) ::Vxloc, Vzloc
        real, intent (IN) :: bega, dt, alpha

        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1) :: s0,s1,s2,s3,s4

        integer :: ngllx, ngllz

        ngllx = Elem%ngllx; ngllz = Elem%ngllz

        VxLoc(1:ngllx-2,1:ngllz-2) = (0.5+alpha) * Elem%Veloc(:,:,0) &
                                   + dt *(0.5-bega)*Elem%Accel(:,:,0) + (0.5-alpha)*Elem%V0(:,:,0)
        VzLoc(1:ngllx-2,1:ngllz-2) = (0.5+alpha) * Elem%Veloc(:,:,1) &
                                   + dt *(0.5-bega)*Elem%Accel(:,:,1) + (0.5-alpha)*Elem%V0(:,:,1)
        s0 = MATMUL (HTmat,VxLoc)
        s1 = MATMUL (VxLoc,Hmatz)
        s2 = MATMUL (HTmat,VzLoc)
        s3 = MATMUL (VzLoc,Hmatz)
        s4 = Dt * Elem%Acoeff(:,:,17)

        ! Updating convolution :
        Elem%PsiVxx(:,:) = Elem%Bx(:,:) *  Elem%PsiVxx(:,:) &
                         - Elem%Ax(:,:) * (Elem%Acoeff(:,:,12)*s0 + Elem%Acoeff(:,:,13)*s1)
        Elem%PsiVxz(:,:) = Elem%Bz(:,:) *  Elem%PsiVxz(:,:) &
                         - Elem%Az(:,:) * (Elem%Acoeff(:,:,14)*s0 + Elem%Acoeff(:,:,15)*s1)
        Elem%PsiVzx(:,:) = Elem%Bx(:,:) *  Elem%PsiVzx(:,:) &
                         - Elem%Ax(:,:) * (Elem%Acoeff(:,:,12)*s2 + Elem%Acoeff(:,:,13)*s3)
        Elem%PsiVzz(:,:) = Elem%Bz(:,:) *  Elem%PsiVzz(:,:) &
                         - Elem%Az(:,:) * (Elem%Acoeff(:,:,14)*s2 + Elem%Acoeff(:,:,15)*s3)

        ! Updating PML Stresses in the PML
        Elem%Stress(:,:,0) = Elem%Stress(:,:,0) + s4 * ((Elem%Lambda + 2*Elem%Mu) &
                           * (Elem%PsiVxx - Elem%Acoeff(:,:,12)*s0 - Elem%Acoeff(:,:,13)*s1) &
             + Elem%Lambda * (Elem%PsiVzz - Elem%Acoeff(:,:,14)*s2 - Elem%Acoeff(:,:,15)*s3))
        Elem%Stress(:,:,1) = Elem%Stress(:,:,1) + s4 * ((Elem%Lambda + 2*Elem%Mu) &
                           * (Elem%PsiVzz - Elem%Acoeff(:,:,14)*s2 - Elem%Acoeff(:,:,15)*s3) &
             + Elem%Lambda * (Elem%PsiVxx - Elem%Acoeff(:,:,12)*s0 - Elem%Acoeff(:,:,13)*s1))
        Elem%Stress(:,:,2) = Elem%Stress(:,:,2) + s4 * Elem%Mu * &
                           ( Elem%PsiVxz(:,:) - Elem%Acoeff(:,:,14)*s0 - Elem%Acoeff(:,:,15)*s1 &
                           + Elem%PsiVzx(:,:) - Elem%Acoeff(:,:,12)*s2 - Elem%Acoeff(:,:,13)*s3)

        return
    end subroutine Prediction_Elem_CPML_Veloc

    ! ###########################################################

    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hprimez
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<


    subroutine  compute_InternalForces_Elem (Elem,hprime, hTprime, hprimez,hTprimez)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime, hTprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1) :: Uxloc, Uzloc, dUx_dxi, dUx_deta, dUz_dxi, dUz_deta, s0

        Uxloc =Elem%Forces (:,:,0)
        Uzloc = Elem%Forces (:,:,1)

        dUx_dxi = MATMUL ( hTprime, Uxloc)
        dUz_dxi = MATMUL ( hTprime, Uzloc )

        dUx_deta = MATMUL ( Uxloc , hprimez )
        dUz_deta = MATMUL ( Uzloc , hprimez )

        s0 = Elem%Acoeff(:,:,0)*dUx_dxi + Elem%Acoeff(:,:,1)*dUx_deta + Elem%Acoeff(:,:,2)*dUz_dxi + Elem%Acoeff(:,:,3)*dUz_deta
        Uxloc = MATMUL ( hprime, s0 )

        s0 = Elem%Acoeff(:,:,1)*dUx_dxi + Elem%Acoeff(:,:,4)*dUx_deta + Elem%Acoeff(:,:,5)*dUz_dxi + Elem%Acoeff(:,:,6)*dUz_deta
        Uxloc = Uxloc +  MATMUL( s0 , hTprimez )

        s0 = Elem%Acoeff(:,:,2)*dUx_dxi + Elem%Acoeff(:,:,5)*dUx_deta + Elem%Acoeff(:,:,7)*dUz_dxi + Elem%Acoeff(:,:,8)*dUz_deta
        Uzloc = MATMUL ( hprime, s0  )

        s0= Elem%Acoeff(:,:,3)*dUx_dxi + Elem%Acoeff(:,:,6)*dUx_deta + Elem%Acoeff(:,:,8)*dUz_dxi + Elem%Acoeff(:,:,9)*dUz_deta
        Uzloc = Uzloc + MATMUL( s0 , hTprimez )

        Elem%Forces(:,:,0) = Uxloc
        Elem%Forces(:,:,1) = Uzloc

        return
    end subroutine compute_InternalForces_Elem



    ! ###########################################################

    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<


    subroutine  compute_InternalForces_DG_Weak (Elem,hprime,hTprimez)
      implicit none

      type (Element), intent (INOUT) :: Elem
      real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime
      real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez
      real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: aux1, aux2

      aux1 = Elem%Acoeff(:,:,0)*Elem%Veloc(:,:,0)
      aux2 = Elem%Acoeff(:,:,1)*Elem%Veloc(:,:,0)
      Elem%Forces(:,:,0) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

      aux1 = Elem%Acoeff(:,:,2)*Elem%Veloc(:,:,1)
      aux2 = Elem%Acoeff(:,:,3)*Elem%Veloc(:,:,1)
      Elem%Forces(:,:,1) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

      aux1 = Elem%Acoeff(:,:,0)*Elem%Veloc(:,:,1) + Elem%Acoeff(:,:,2)*Elem%Veloc(:,:,0)
      aux2 = Elem%Acoeff(:,:,1)*Elem%Veloc(:,:,1) + Elem%Acoeff(:,:,3)*Elem%Veloc(:,:,0)
      Elem%Forces(:,:,2) = 0.5 * (MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez))

      aux1 = (Elem%Acoeff(:,:,4) + Elem%Acoeff(:,:,5))*Elem%Strain(:,:,0) + Elem%Acoeff(:,:,4)*Elem%Strain(:,:,1) &
           + Elem%Acoeff(:,:,7)*Elem%Strain(:,:,2)
      aux2 = (Elem%Acoeff(:,:,8) + Elem%Acoeff(:,:,9))*Elem%Strain(:,:,0) + Elem%Acoeff(:,:,8)*Elem%Strain(:,:,1) &
           + Elem%Acoeff(:,:,11)*Elem%Strain(:,:,2)
      Elem%Forces(:,:,3) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

      aux1 = (Elem%Acoeff(:,:,7) + Elem%Acoeff(:,:,6))*Elem%Strain(:,:,1) + Elem%Acoeff(:,:,6)*Elem%Strain(:,:,0) &
           + Elem%Acoeff(:,:,5)*Elem%Strain(:,:,2)
      aux2 = (Elem%Acoeff(:,:,11) + Elem%Acoeff(:,:,10))*Elem%Strain(:,:,1) + Elem%Acoeff(:,:,10)*Elem%Strain(:,:,0) &
           + Elem%Acoeff(:,:,9)*Elem%Strain(:,:,2)
      Elem%Forces(:,:,4) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

      Elem%Forces(:,:,:) = - Elem%Forces(:,:,:)

      return
    end subroutine compute_InternalForces_DG_Weak


    ! ###########################################################

    !>
    !! \brief Calcul des forces internes HDG pour partie acoustique en formulation
    !! vitesse - deformation scalaire.
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<

    subroutine  compute_InternalForces_HDG_Weak (Elem,hprime,hTprimez)
      implicit none

      type (Element), intent (INOUT) :: Elem
      real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime
      real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez
      real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: aux1, aux2


      if (Elem%acoustic) then ! ACOUSTIC CASE
          aux1 = Elem%Acoeff(:,:,0)*Elem%Veloc(:,:,0) + Elem%Acoeff(:,:,2)*Elem%Veloc(:,:,1)
          aux2 = Elem%Acoeff(:,:,1)*Elem%Veloc(:,:,0) + Elem%Acoeff(:,:,3)*Elem%Veloc(:,:,1)
          Elem%Forces(:,:,0) = - MATMUL(hprime,aux1) - MATMUL(aux2,hTprimez)

          aux1 = Elem%Lambda(:,:)*Elem%Strain(:,:,0)
          Elem%Forces(:,:,1) = - MATMUL(hprime,(Elem%Acoeff(:,:,0)*aux1)) - MATMUL((Elem%Acoeff(:,:,1)*aux1),hTprimez)
          Elem%Forces(:,:,2) = - MATMUL(hprime,(Elem%Acoeff(:,:,2)*aux1)) - MATMUL((Elem%Acoeff(:,:,3)*aux1),hTprimez)

      else ! ELASTIC CASE
          aux1 = Elem%Acoeff(:,:,0)*Elem%Veloc(:,:,0)
          aux2 = Elem%Acoeff(:,:,1)*Elem%Veloc(:,:,0)
          Elem%Forces(:,:,0) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

          aux1 = Elem%Acoeff(:,:,2)*Elem%Veloc(:,:,1)
          aux2 = Elem%Acoeff(:,:,3)*Elem%Veloc(:,:,1)
          Elem%Forces(:,:,1) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

          aux1 = Elem%Acoeff(:,:,0)*Elem%Veloc(:,:,1) + Elem%Acoeff(:,:,2)*Elem%Veloc(:,:,0)
          aux2 = Elem%Acoeff(:,:,1)*Elem%Veloc(:,:,1) + Elem%Acoeff(:,:,3)*Elem%Veloc(:,:,0)
          Elem%Forces(:,:,2) = 0.5 * (MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez))

          aux1 = (Elem%Acoeff(:,:,4) + Elem%Acoeff(:,:,5))*Elem%Strain(:,:,0) + Elem%Acoeff(:,:,4)*Elem%Strain(:,:,1) &
              + Elem%Acoeff(:,:,7)*Elem%Strain(:,:,2)
          aux2 = (Elem%Acoeff(:,:,8) + Elem%Acoeff(:,:,9))*Elem%Strain(:,:,0) + Elem%Acoeff(:,:,8)*Elem%Strain(:,:,1) &
              + Elem%Acoeff(:,:,11)*Elem%Strain(:,:,2)
          Elem%Forces(:,:,3) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

          aux1 = (Elem%Acoeff(:,:,7) + Elem%Acoeff(:,:,6))*Elem%Strain(:,:,1) + Elem%Acoeff(:,:,6)*Elem%Strain(:,:,0) &
              + Elem%Acoeff(:,:,5)*Elem%Strain(:,:,2)
          aux2 = (Elem%Acoeff(:,:,11) + Elem%Acoeff(:,:,10))*Elem%Strain(:,:,1) + Elem%Acoeff(:,:,10)*Elem%Strain(:,:,0) &
              + Elem%Acoeff(:,:,9)*Elem%Strain(:,:,2)
          Elem%Forces(:,:,4) = MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez)

          Elem%Forces(:,:,:) = - Elem%Forces(:,:,:)
      endif

      return
    end subroutine compute_InternalForces_HDG_Weak


    ! ###########################################################

    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hprimez
    !<

    subroutine  compute_InternalForces_DG_Strong (Elem,hTprime,hprimez)
      implicit none

      type (Element), intent (INOUT) :: Elem
      real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
      real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
      real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: aux1, aux2

      aux1 = MATMUL(hTprime,Elem%Veloc(:,:,0))
      aux2 = MATMUL(Elem%Veloc(:,:,0),hprimez)
      Elem%Forces(:,:,0) = Elem%Acoeff(:,:,0)*aux1 + Elem%Acoeff(:,:,1)*aux2

      aux1 = MATMUL(hTprime,Elem%Veloc(:,:,1))
      aux2 = MATMUL(Elem%Veloc(:,:,1),hprimez)
      Elem%Forces(:,:,1) = Elem%Acoeff(:,:,2)*aux1 + Elem%Acoeff(:,:,3)*aux2

      aux1 = Elem%Acoeff(:,:,0) * MATMUL(hTprime,Elem%Veloc(:,:,1)) +  Elem%Acoeff(:,:,1) * MATMUL(Elem%Veloc(:,:,1),hprimez)
      aux2 = Elem%Acoeff(:,:,2) * MATMUL(hTprime,Elem%Veloc(:,:,0)) +  Elem%Acoeff(:,:,3) * MATMUL(Elem%Veloc(:,:,0),hprimez)
      Elem%Forces(:,:,2) = 0.5 * (aux1 + aux2)

      aux1 = (Elem%Acoeff(:,:,4) + Elem%Acoeff(:,:,5)) * MATMUL(hTprime,Elem%Strain(:,:,0)) &
           + Elem%Acoeff(:,:,4) * MATMUL(hTprime,Elem%Strain(:,:,1)) &
           + Elem%Acoeff(:,:,7) * MATMUL(hTprime,Elem%Strain(:,:,2))
      aux2 = (Elem%Acoeff(:,:,8) + Elem%Acoeff(:,:,9)) * MATMUL(Elem%Strain(:,:,0),hprimez) &
           + Elem%Acoeff(:,:,8) * MATMUL(Elem%Strain(:,:,1),hprimez) &
           + Elem%Acoeff(:,:,11)* MATMUL(Elem%Strain(:,:,2),hprimez)
      Elem%Forces(:,:,3) = aux1 + aux2

      aux1 = (Elem%Acoeff(:,:,7) + Elem%Acoeff(:,:,6)) * MATMUL(hTprime,Elem%Strain(:,:,1)) &
           + Elem%Acoeff(:,:,6) * MATMUL(hTprime,Elem%Strain(:,:,0)) &
           + Elem%Acoeff(:,:,5) * MATMUL(hTprime,Elem%Strain(:,:,2))
      aux2 = (Elem%Acoeff(:,:,11) + Elem%Acoeff(:,:,10)) * MATMUL(Elem%Strain(:,:,1),hprimez) &
           + Elem%Acoeff(:,:,10)* MATMUL(Elem%Strain(:,:,0),hprimez) &
           + Elem%Acoeff(:,:,9) * MATMUL(Elem%Strain(:,:,2),hprimez)
      Elem%Forces(:,:,4) = aux1 + aux2

      return
    end subroutine compute_InternalForces_DG_Strong


    ! ###########################################################
    !>
    !! \brief This subroutine computes the "internal forces" terms for the equations
    !! of evolution of strains (in a "weak" form), and the equation of evolution of
    !! velocities (in a "strong" form i.e after anothe integration by part).
    !! It is suitable for Hybridizable Discontinuous Galerkin elements only,
    !! and used in a framework of semi-implicit time scheme.
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hprimez
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<
    subroutine  compute_InternalForces_HDG (Elem,hprime,hTprime,hprimez,hTprimez)
      implicit none

      type (Element), intent (INOUT) :: Elem
      real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime , hTprime
      real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
      real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: aux1, aux2

      aux1 = Elem%Acoeff(:,:,0)*Elem%Veloc(:,:,0)
      aux2 = Elem%Acoeff(:,:,1)*Elem%Veloc(:,:,0)
      Elem%Forces(:,:,0) = - MATMUL(hprime,aux1) - MATMUL(aux2,hTprimez)

      aux1 = Elem%Acoeff(:,:,2)*Elem%Veloc(:,:,1)
      aux2 = Elem%Acoeff(:,:,3)*Elem%Veloc(:,:,1)
      Elem%Forces(:,:,1) = - MATMUL(hprime,aux1) - MATMUL(aux2,hTprimez)

      aux1 = Elem%Acoeff(:,:,0)*Elem%Veloc(:,:,1) + Elem%Acoeff(:,:,2)*Elem%Veloc(:,:,0)
      aux2 = Elem%Acoeff(:,:,1)*Elem%Veloc(:,:,1) + Elem%Acoeff(:,:,3)*Elem%Veloc(:,:,0)
      Elem%Forces(:,:,2) = - 0.5 * (MATMUL(hprime,aux1) + MATMUL(aux2,hTprimez))

      aux1 = (Elem%Acoeff(:,:,4) + Elem%Acoeff(:,:,5)) * MATMUL(hTprime,Elem%Strain(:,:,0)) &
           + Elem%Acoeff(:,:,4) * MATMUL(hTprime,Elem%Strain(:,:,1)) &
           + Elem%Acoeff(:,:,7) * MATMUL(hTprime,Elem%Strain(:,:,2))
      aux2 = (Elem%Acoeff(:,:,8) + Elem%Acoeff(:,:,9)) * MATMUL(Elem%Strain(:,:,0),hprimez) &
           + Elem%Acoeff(:,:,8) * MATMUL(Elem%Strain(:,:,1),hprimez) &
           + Elem%Acoeff(:,:,11)* MATMUL(Elem%Strain(:,:,2),hprimez)
      Elem%Forces(:,:,3) = aux1 + aux2

      aux1 = (Elem%Acoeff(:,:,7) + Elem%Acoeff(:,:,6)) * MATMUL(hTprime,Elem%Strain(:,:,1)) &
           + Elem%Acoeff(:,:,6) * MATMUL(hTprime,Elem%Strain(:,:,0)) &
           + Elem%Acoeff(:,:,5) * MATMUL(hTprime,Elem%Strain(:,:,2))
      aux2 = (Elem%Acoeff(:,:,11) + Elem%Acoeff(:,:,10)) * MATMUL(Elem%Strain(:,:,1),hprimez) &
           + Elem%Acoeff(:,:,10)* MATMUL(Elem%Strain(:,:,0),hprimez) &
           + Elem%Acoeff(:,:,9) * MATMUL(Elem%Strain(:,:,2),hprimez)
      Elem%Forces(:,:,4) = aux1 + aux2

    end subroutine compute_InternalForces_HDG


    ! ###########################################################
    !>
    !! \brief Thus subroutine is used in a Runge-Kutta 4 framework :
    !! It updates the variables for each of the 5 steps of the LSERK4 time scheme.
    !! For continuous Galerkin elements, velocities and displacements are updated ;
    !! For discontinuous Galerkin elements, velocities and strains are updated.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) coeff1
    !! \param real, intent (IN) coeff2
    !! \param real, intent (IN) Dt
    !<
    subroutine  update_Elem_RK4 (Elem,coeff1,coeff2,Dt)
        implicit none

        type(Element), intent (INOUT) :: Elem
        real, intent(IN) :: coeff1
        real, intent(IN) :: coeff2
        real, intent(IN) :: Dt
        integer          :: ngx, ngz

        if (Elem%type_DG==GALERKIN_CONT) then
            ngx = Elem%ngllx ; ngz = Elem%ngllz
            Elem%Vect_RK(:,:,0:1) = coeff1 * Elem%Vect_RK(:,:,0:1) &
                                  + Dt * Elem%Forces(1:ngx-2,1:ngz-2,0:1)
            Elem%Vect_RK(:,:,2:3) = coeff1 * Elem%Vect_RK(:,:,2:3)  + Dt * Elem%Veloc(:,:,:)
            Elem%Veloc = Elem%Veloc + coeff2 * Elem%Vect_RK(:,:,0:1)
            Elem%Displ = Elem%Displ + coeff2 * Elem%Vect_RK(:,:,2:3)
        else if (Elem%Acoustic) then
            Elem%Vect_RK = coeff1 * Elem%Vect_RK + Dt * Elem%Forces
            Elem%Strain(:,:,0)  = Elem%Strain(:,:,0) + coeff2 * Elem%Vect_RK(:,:,0)
            Elem%Veloc   = Elem%Veloc  + coeff2 * Elem%Vect_RK(:,:,1:2)
        else ! Elastic case
            Elem%Vect_RK = coeff1 * Elem%Vect_RK + Dt * Elem%Forces
            Elem%Strain  = Elem%Strain + coeff2 * Elem%Vect_RK(:,:,0:2)
            Elem%Veloc   = Elem%Veloc  + coeff2 * Elem%Vect_RK(:,:,3:4)
        endif

    end subroutine update_Elem_RK4

    ! ###########################################################
    !>
    !! \brief Thus subroutine adds the contribution of the memory PML terms
    !! for the ADE-PML in case of DG or HDG.
    !! \param type (Element), intent (INOUT) Elem
    !<
    subroutine  add_Psi4PML (Elem)
        implicit none

        type (Element), intent (INOUT) :: Elem

        if (Elem%acoustic) then
            Elem%Forces(:,:,0) = Elem%Forces(:,:,0) + Elem%PsiVxx(:,:) + Elem%PsiVzz(:,:)
            Elem%Forces(:,:,1) = Elem%Forces(:,:,1) + Elem%PsiSxxx(:,:)
            Elem%Forces(:,:,2) = Elem%Forces(:,:,2) + Elem%PsiSzzz(:,:)
        else ! Elastic Case
            Elem%Forces(:,:,0) = Elem%Forces(:,:,0) + Elem%PsiVxx(:,:)
            Elem%Forces(:,:,1) = Elem%Forces(:,:,1) + Elem%PsiVzz(:,:)
            Elem%Forces(:,:,2) = Elem%Forces(:,:,2) + 0.5 * (Elem%PsiVxz(:,:) + Elem%PsiVzx(:,:))
            Elem%Forces(:,:,3) = Elem%Forces(:,:,3) + Elem%PsiSxxx(:,:) + Elem%PsiSxzz(:,:)
            Elem%Forces(:,:,4) = Elem%Forces(:,:,4) + Elem%PsiSxzx(:,:) + Elem%PsiSzzz(:,:)
        endif

        return
    end subroutine add_Psi4PML

    ! ###########################################################
    !>
    !! \brief Initialize the memory variables Psi for ADE-PML in a context
    !! of a predictor-corrector-like time scheme.
    !! USED ONLY for the ADE-PML in case of DG or HDG.
    !! \param type (Element), intent (INOUT) Elem
    !<
    subroutine  initialize_Psi (Elem)
        implicit none

        type (Element), intent (INOUT) :: Elem

        if (Elem%acoustic) then
            Elem%Psi_store(:,:,0) = Elem%PsiSxxx(:,:)
            Elem%Psi_store(:,:,1) = Elem%PsiSzzz(:,:)
            Elem%Psi_store(:,:,2) = Elem%PsiVxx(:,:)
            Elem%Psi_store(:,:,3) = Elem%PsiVzz(:,:)
        else ! Elastic Case
            Elem%Psi_store(:,:,0) = Elem%PsiSxxx(:,:)
            Elem%Psi_store(:,:,1) = Elem%PsiSzzz(:,:)
            Elem%Psi_store(:,:,2) = Elem%PsiSxzx(:,:)
            Elem%Psi_store(:,:,3) = Elem%PsiSxzz(:,:)
            Elem%Psi_store(:,:,4) = Elem%PsiVxx(:,:)
            Elem%Psi_store(:,:,5) = Elem%PsiVxz(:,:)
            Elem%Psi_store(:,:,6) = Elem%PsiVzx(:,:)
            Elem%Psi_store(:,:,7) = Elem%PsiVzz(:,:)
        endif

        return
    end subroutine initialize_Psi

    ! ###########################################################
    !>
    !! \brief Thus subroutine is used in a ADE-PML framework :
    !! the memory terms PsiS*** and PsiV*** are updated for the next iteration.
    !! In a RK4 framework, the ADE for the Psi variables is solved using RK4,
    !! like the others variables.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
    !! \param real, intent (IN) coeff1
    !! \param real, intent (IN) coeff2
    !! \param real, intent (IN) Dt
    !<
    subroutine  update_Psi_ADEPML_RK4 (Elem,hTprime,hprimez,coeff1,coeff2,Dt)
        implicit none

        type(Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
        real, intent(IN) :: coeff1
        real, intent(IN) :: coeff2
        real, intent(IN) :: Dt
        real, dimension(0:Elem%ngllx-1,0:Elem%ngllz-1,0:7) :: smbr

        if (Elem%Acoustic) then
            ! Compute second member of memory variables evolution equation
            call compute_smbr_acou_ADEPML (Elem,hTprime,hprimez,smbr)
            ! UPDATING in time using usual LSERK4
            Elem%Psi_store(:,:,:)  = coeff1 * Elem%Psi_store(:,:,:) + Dt*smbr(:,:,0:3)
            Elem%PsiSxxx(:,:) = Elem%PsiSxxx(:,:) + coeff2 * Elem%Psi_store(:,:,0)
            Elem%PsiSzzz(:,:) = Elem%PsiSzzz(:,:) + coeff2 * Elem%Psi_store(:,:,1)
            Elem%PsiVxx(:,:)  = Elem%PsiVxx (:,:) + coeff2 * Elem%Psi_store(:,:,2)
            Elem%PsiVzz(:,:)  = Elem%PsiVzz (:,:) + coeff2 * Elem%Psi_store(:,:,3)

        else ! Elastic Case
            ! Compute second member of memory variables evolution equation
            call compute_smbr_ADEPML (Elem,hTprime,hprimez,smbr)
            ! UPDATING in time using usual LSERK4
            Elem%Psi_store(:,:,:)  = coeff1 * Elem%Psi_store(:,:,:) + Dt*smbr(:,:,:)
            Elem%PsiSxxx(:,:) = Elem%PsiSxxx(:,:) + coeff2 * Elem%Psi_store(:,:,0)
            Elem%PsiSzzz(:,:) = Elem%PsiSzzz(:,:) + coeff2 * Elem%Psi_store(:,:,1)
            Elem%PsiSxzx(:,:) = Elem%PsiSxzx(:,:) + coeff2 * Elem%Psi_store(:,:,2)
            Elem%PsiSxzz(:,:) = Elem%PsiSxzz(:,:) + coeff2 * Elem%Psi_store(:,:,3)
            Elem%PsiVxx(:,:)  = Elem%PsiVxx (:,:) + coeff2 * Elem%Psi_store(:,:,4)
            Elem%PsiVxz(:,:)  = Elem%PsiVxz (:,:) + coeff2 * Elem%Psi_store(:,:,5)
            Elem%PsiVzx(:,:)  = Elem%PsiVzx (:,:) + coeff2 * Elem%Psi_store(:,:,6)
            Elem%PsiVzz(:,:)  = Elem%PsiVzz (:,:) + coeff2 * Elem%Psi_store(:,:,7)
        endif

        return
    end subroutine update_Psi_ADEPML_RK4

    ! ###########################################################
    !>
    !! \brief Thus subroutine is used in a ADE-PML framework :
    !! the memory terms PsiS*** and PsiV*** are updated for the next iteration.
    !! In a RK4 framework, the ADE for the Psi variables is solved using RK4,
    !! like the others variables.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
    !! \param real, intent (IN) Dt
    !<
    subroutine  update_Psi_ADEPML (Elem,hTprime,hprimez,Dt)
        implicit none

        type(Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
        real, intent(IN) :: Dt
        real, dimension(0:Elem%ngllx-1,0:Elem%ngllz-1,0:7) :: smbr

        if (Elem%ADEPML) then
            if (Elem%acoustic) then
                ! Compute second member of memory variables evolution equation
                call compute_smbr_acou_ADEPML (Elem,hTprime,hprimez,smbr)
                ! UPDATING in time for midpoint schemes
                Elem%PsiSxxx(:,:) = Elem%Psi_store(:,:,0) + Dt * smbr(:,:,0)
                Elem%PsiSzzz(:,:) = Elem%Psi_store(:,:,1) + Dt * smbr(:,:,1)
                Elem%PsiVxx(:,:)  = Elem%Psi_store(:,:,2) + Dt * smbr(:,:,2)
                Elem%PsiVzz(:,:)  = Elem%Psi_store(:,:,3) + Dt * smbr(:,:,3)

            else ! Elastic Case
                ! Compute second member of memory variables evolution equation
                call compute_smbr_ADEPML (Elem,hTprime,hprimez,smbr)
                ! UPDATING in time for midpoint schemes
                Elem%PsiSxxx(:,:) = Elem%Psi_store(:,:,0) + Dt * smbr(:,:,0)
                Elem%PsiSzzz(:,:) = Elem%Psi_store(:,:,1) + Dt * smbr(:,:,1)
                Elem%PsiSxzx(:,:) = Elem%Psi_store(:,:,2) + Dt * smbr(:,:,2)
                Elem%PsiSxzz(:,:) = Elem%Psi_store(:,:,3) + Dt * smbr(:,:,3)
                Elem%PsiVxx(:,:)  = Elem%Psi_store(:,:,4) + Dt * smbr(:,:,4)
                Elem%PsiVxz(:,:)  = Elem%Psi_store(:,:,5) + Dt * smbr(:,:,5)
                Elem%PsiVzx(:,:)  = Elem%Psi_store(:,:,6) + Dt * smbr(:,:,6)
                Elem%PsiVzz(:,:)  = Elem%Psi_store(:,:,7) + Dt * smbr(:,:,7)
            endif
        else if (Elem%CPML) then
            smbr(:,:,:) = -smbr(:,:,:)
            Elem%PsiSxxx(:,:) = smbr(:,:,0)
            Elem%PsiSzzz(:,:) = smbr(:,:,1)
            Elem%PsiSxzx(:,:) = smbr(:,:,2)
            Elem%PsiSxzz(:,:) = smbr(:,:,3)
            Elem%PsiVxx(:,:)  = smbr(:,:,4)
            Elem%PsiVxz(:,:)  = smbr(:,:,5)
            Elem%PsiVzx(:,:)  = smbr(:,:,6)
            Elem%PsiVzz(:,:)  = smbr(:,:,7)
        endif

        ! Adding PML terms to the Elem%Forces
        call add_Psi4PML (Elem)

    end subroutine update_Psi_ADEPML


    ! ###########################################################
    !>
    !! \brief Thus subroutine is used in a ADE-PML framework : it computes the second member
    !!  of memory variables evolution equation for an ACOUSTIC element.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1, 0:7), intent (IN) :: smbr
    !<
    subroutine  compute_smbr_acou_ADEPML (Elem,hTprime,hprimez,smbr)
        implicit none

        type(Element), intent (IN) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
        real, dimension (0:Elem%ngllx-1,0:Elem%ngllz-1,0:7), intent (INOUT) :: smbr
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1) :: aux

        ! Defining Second Member of time evolution equation for the Psi
        ! Second member for Pressure memory variables
        aux(:,:) = Elem%Lambda(:,:) * Elem%Strain(:,:,0)
        smbr(:,:,0) = - Elem%Bx(:,:) * Elem%PsiSxxx(:,:) - Elem%Ax(:,:) * &
                      ( Elem%Acoeff(:,:,0) * MATMUL(HTprime,aux) &
                      + Elem%Acoeff(:,:,1) * MATMUL(aux,Hprimez))
        smbr(:,:,1) = - Elem%Bz(:,:) * Elem%PsiSzzz(:,:) - Elem%Az(:,:) * &
                      ( Elem%Acoeff(:,:,2) * MATMUL(HTprime,aux) &
                      + Elem%Acoeff(:,:,3) * MATMUL(aux,Hprimez))
        ! Second member for Velocities memory variables
        smbr(:,:,2) = - Elem%Bx(:,:) * Elem%PsiVxx(:,:) - Elem%Ax(:,:) * &
                      ( Elem%Acoeff(:,:,0) * MATMUL(HTprime,Elem%Veloc(:,:,0)) &
                      + Elem%Acoeff(:,:,1) * MATMUL(Elem%Veloc(:,:,0),Hprimez))
        smbr(:,:,3) = - Elem%Bz(:,:) * Elem%PsiVzz(:,:) - Elem%Az(:,:) * &
                      ( Elem%Acoeff(:,:,2) * MATMUL(HTprime,Elem%Veloc(:,:,1)) &
                      + Elem%Acoeff(:,:,3) * MATMUL(Elem%Veloc(:,:,1),Hprimez))

        return
    end subroutine compute_smbr_acou_ADEPML


    ! ###########################################################
    !>
    !! \brief Thus subroutine is used in a ADE-PML framework :
    !! It computes the second member of memory variables evolution equation.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1, 0:7), intent (IN) :: smbr
    !<
    subroutine  compute_smbr_ADEPML (Elem,hTprime,hprimez,smbr)
        implicit none

        type(Element), intent (IN) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
        real, dimension (0:Elem%ngllx-1,0:Elem%ngllz-1,0:7), intent (INOUT) :: smbr

        ! Defining Second Member of time evolution equation for the Psi
        ! Second member for Stresses memory variables
        smbr(:,:,0) = - Elem%Bx(:,:) * Elem%PsiSxxx(:,:) - Elem%Ax(:,:) * &
                      ((Elem%Acoeff(:,:,4) + Elem%Acoeff(:,:,5)) * MATMUL(HTprime,Elem%Strain(:,:,0)) &
                      + Elem%Acoeff(:,:,4) * MATMUL(HTprime,Elem%Strain(:,:,1)) &
                      +(Elem%Acoeff(:,:,8) + Elem%Acoeff(:,:,9)) * MATMUL(Elem%Strain(:,:,0),Hprimez) &
                      + Elem%Acoeff(:,:,8) * MATMUL(Elem%Strain(:,:,1),Hprimez))
        smbr(:,:,1) = - Elem%Bz(:,:) * Elem%PsiSzzz(:,:) - Elem%Az(:,:) * &
                      ((Elem%Acoeff(:,:,6) + Elem%Acoeff(:,:,7)) * MATMUL(HTprime,Elem%Strain(:,:,1)) &
                      + Elem%Acoeff(:,:,6) * MATMUL(HTprime,Elem%Strain(:,:,0)) &
                      +(Elem%Acoeff(:,:,10) + Elem%Acoeff(:,:,11)) * MATMUL(Elem%Strain(:,:,1),Hprimez) &
                      + Elem%Acoeff(:,:,10) * MATMUL(Elem%Strain(:,:,0),Hprimez))
        smbr(:,:,2) = - Elem%Bx(:,:) * Elem%PsiSxzx(:,:) - Elem%Ax(:,:) * &
                      ( Elem%Acoeff(:,:,5) * MATMUL(HTprime,Elem%Strain(:,:,2)) &
                      + Elem%Acoeff(:,:,9) * MATMUL(Elem%Strain(:,:,2),Hprimez))
        smbr(:,:,3) = - Elem%Bz(:,:) * Elem%PsiSxzz(:,:) - Elem%Az(:,:) * &
                      ( Elem%Acoeff(:,:,7) * MATMUL(HTprime,Elem%Strain(:,:,2)) &
                      + Elem%Acoeff(:,:,11) * MATMUL(Elem%Strain(:,:,2),Hprimez))
        ! Second member for Velocities memory variables
        smbr(:,:,4) = - Elem%Bx(:,:) * Elem%PsiVxx(:,:) - Elem%Ax(:,:) * &
                      ( Elem%Acoeff(:,:,0) * MATMUL(HTprime,Elem%Veloc(:,:,0)) &
                      + Elem%Acoeff(:,:,1) * MATMUL(Elem%Veloc(:,:,0),Hprimez))
        smbr(:,:,5) = - Elem%Bz(:,:) * Elem%PsiVxz(:,:) - Elem%Az(:,:) * &
                      ( Elem%Acoeff(:,:,2) * MATMUL(HTprime,Elem%Veloc(:,:,0)) &
                      + Elem%Acoeff(:,:,3) * MATMUL(Elem%Veloc(:,:,0),Hprimez))
        smbr(:,:,6) = - Elem%Bx(:,:) * Elem%PsiVzx(:,:) - Elem%Ax(:,:) * &
                      ( Elem%Acoeff(:,:,0) * MATMUL(HTprime,Elem%Veloc(:,:,1)) &
                      + Elem%Acoeff(:,:,1) * MATMUL(Elem%Veloc(:,:,1),Hprimez))
        smbr(:,:,7) = - Elem%Bz(:,:) * Elem%PsiVzz(:,:) - Elem%Az(:,:) * &
                      ( Elem%Acoeff(:,:,2) * MATMUL(HTprime,Elem%Veloc(:,:,1)) &
                      + Elem%Acoeff(:,:,3) * MATMUL(Elem%Veloc(:,:,1),Hprimez))

        return
    end subroutine compute_smbr_ADEPML

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<
    subroutine  compute_InternalForces_PML_Elem (Elem,hprime, hTprimez)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hTprimez
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: s0,s1

        s0 = Elem%Acoeff(:,:,12) * Elem%Stress(:,:,0) + Elem%Acoeff(:,:,14) * Elem%Stress(:,:,2)
        s1= MATMUL (hprime,s0)
        Elem%Forces1(:,:,0) = s1
        s0 = Elem%Acoeff(:,:,15) * Elem%Stress (:,:,2)  + Elem%Acoeff(:,:,13) * Elem%Stress (:,:,0)
        s1  =   MATMUL (s0,hTprimez)

        Elem%Forces2(:,:,0) = s1
        s0 = Elem%Acoeff(:,:,12) * Elem%Stress(:,:,2) + Elem%Acoeff(:,:,14)* Elem%Stress(:,:,1)
        s1 = MATMUL (hprime,s0)
        Elem%Forces1(:,:,1) = s1
        s0 = Elem%Acoeff(:,:,15) * Elem%Stress(:,:,1) + Elem%Acoeff(:,:,13) * Elem%Stress (:,:,2)
        s1  =   MATMUL (s0,hTprimez)
        Elem%Forces2(:,:,1) = s1

        Elem%Forces = Elem%Forces1 + Elem%Forces2

        return
    end subroutine compute_InternalForces_PML_Elem



    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hprime
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hprimez
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hTprimez
    !<

    subroutine  compute_InternalForces_CPML_Elem (Elem, hprime, hTprime, hprimez, hTprimez)
        implicit none

        type (Element), intent (INOUT) :: Elem

        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: dAxdx,dAzdz,s0,s1
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hprime, hTprime
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        logical :: withIntegrationByPart  = .false.

        dAxdx = Elem%Acoeff(:,:,12)*MATMUL(hTprime,Elem%Ax) &
              + Elem%Acoeff(:,:,13)*MATMUL(Elem%Ax,hprimez)
        dAzdz = Elem%Acoeff(:,:,14)*MATMUL(hTprime,Elem%Az) &
              + Elem%Acoeff(:,:,15)*MATMUL(Elem%Az,hprimez)

        ! Updating convolution :
        if (withIntegrationByPart) then ! With Integration by part
        Elem%PsiSxxx(:,:) = Elem%Bx(:,:)*Elem%PsiSxxx(:,:) + Elem%Stress(:,:,0)*dAxdx(:,:) &
                          + MATMUL(hprime,Elem%Stress(:,:,0)*Elem%Acoeff(:,:,12)*Elem%Ax(:,:)) &
                          + MATMUL(Elem%Stress(:,:,0)*Elem%Acoeff(:,:,13)*Elem%Ax(:,:),hTprimez)
        Elem%PsiSzzz(:,:) = Elem%Bz(:,:)*Elem%PsiSzzz(:,:) + Elem%Stress(:,:,1)*dAzdz(:,:) &
                          + MATMUL(hprime,Elem%Stress(:,:,1)*Elem%Acoeff(:,:,14)*Elem%Az(:,:)) &
                          + MATMUL(Elem%Stress(:,:,1)*Elem%Acoeff(:,:,15)*Elem%Az(:,:),hTprimez)
        Elem%PsiSxzx(:,:) = Elem%Bx(:,:)*Elem%PsiSxzx(:,:) + Elem%Stress(:,:,2)*dAxdx(:,:) &
                          + MATMUL(hprime,Elem%Stress(:,:,2)*Elem%Acoeff(:,:,12)*Elem%Ax(:,:)) &
                          + MATMUL(Elem%Stress(:,:,2)*Elem%Acoeff(:,:,13)*Elem%Ax(:,:),hTprimez)
        Elem%PsiSxzz(:,:) = Elem%Bz(:,:)*Elem%PsiSxzz(:,:) + Elem%Stress(:,:,2)*dAzdz(:,:) &
                          + MATMUL(hprime,Elem%Stress(:,:,2)*Elem%Acoeff(:,:,14)*Elem%Az(:,:)) &
                          + MATMUL(Elem%Stress(:,:,2)*Elem%Acoeff(:,:,15)*Elem%Az(:,:),hTprimez)
        else ! Without Integration by part
        Elem%PsiSxxx(:,:) = Elem%Bx(:,:)*Elem%PsiSxxx(:,:) - Elem%Ax(:,:) * &
                          ( Elem%Acoeff(:,:,12) * MATMUL(hTprime,Elem%Stress(:,:,0)) &
                          + Elem%Acoeff(:,:,13) * MATMUL(Elem%Stress(:,:,0),hprimez))
        Elem%PsiSzzz(:,:) = Elem%Bz(:,:)*Elem%PsiSzzz(:,:) - Elem%Az(:,:) * &
                          ( Elem%Acoeff(:,:,14) * MATMUL(hTprime,Elem%Stress(:,:,1)) &
                          + Elem%Acoeff(:,:,15) * MATMUL(Elem%Stress(:,:,1),hprimez))
        Elem%PsiSxzx(:,:) = Elem%Bx(:,:)*Elem%PsiSxzx(:,:) - Elem%Ax(:,:) * &
                          ( Elem%Acoeff(:,:,12) * MATMUL(hTprime,Elem%Stress(:,:,2)) &
                          + Elem%Acoeff(:,:,13) * MATMUL(Elem%Stress(:,:,2),hprimez))
        Elem%PsiSxzz(:,:) = Elem%Bz(:,:)*Elem%PsiSxzz(:,:) - Elem%Az(:,:) * &
                          ( Elem%Acoeff(:,:,14) * MATMUL(hTprime,Elem%Stress(:,:,2)) &
                          + Elem%Acoeff(:,:,15) * MATMUL(Elem%Stress(:,:,2),hprimez))
        endif

        ! Updating Forces :
        s0 = Elem%Acoeff(:,:,12) * Elem%Stress(:,:,0) + Elem%Acoeff(:,:,14) * Elem%Stress(:,:,2)
        s1 = Elem%Acoeff(:,:,13) * Elem%Stress(:,:,0) + Elem%Acoeff(:,:,15) * Elem%Stress(:,:,2)
        Elem%Forces(:,:,0) = MATMUL(hprime,s0) + MATMUL(s1,hTprimez) + Elem%PsiSxxx + Elem%PsiSxzz

        s0 = Elem%Acoeff(:,:,12) * Elem%Stress(:,:,2) + Elem%Acoeff(:,:,14) * Elem%Stress(:,:,1)
        s1 = Elem%Acoeff(:,:,13) * Elem%Stress(:,:,2) + Elem%Acoeff(:,:,15) * Elem%Stress(:,:,1)
        Elem%Forces(:,:,1) = MATMUL(hprime,s0) + MATMUL(s1,hTprimez) + Elem%PsiSxzx + Elem%PsiSzzz

        return
    end subroutine compute_InternalForces_CPML_Elem


    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !<
    subroutine inversion_massmat(Elem)

      implicit none
      type (Element), intent (INOUT) :: Elem
      integer  :: ngllx,ngllz,i

      if(Elem%type_DG==GALERKIN_CONT) then
         ngllx = Elem%ngllx ; ngllz = Elem%ngllz
         do i = 0,1
            Elem%Forces(1:ngllx-2,1:ngllz-2,i) = Elem%MassMat(:,:)  * Elem%Forces(1:ngllx-2,1:ngllz-2,i)
         enddo
      else if (Elem%Acoustic) then
         Elem%Forces(:,:,0) = 1./Elem%Acoeff(:,:,4) * Elem%Forces(:,:,0)
         do i = 1,2
            Elem%Forces(:,:,i) = Elem%MassMat(:,:) * Elem%Forces(:,:,i)
         enddo
      else ! Elastic Case
         do i = 0,2
            Elem%Forces(:,:,i) = 1./Elem%Acoeff(:,:,12) * Elem%Forces(:,:,i)
         enddo
         do i = 3,4
            Elem%Forces(:,:,i) = Elem%MassMat(:,:) * Elem%Forces(:,:,i)
         enddo
      endif

    end subroutine inversion_massmat


    ! ###########################################################
    !>
    !! \brief This subroutine computes the energy of the elastic deformation of the element
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTprime
    !! \param real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hprimez
    !! \param real, intent (INOUT) E_elas
    !<

    subroutine  compute_Elastic_Energy (Elem, hTprime, hprimez, E_elas)
        implicit none

        type (Element), intent (IN) :: Elem
        real, intent (INOUT) :: E_elas
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: hTprime
        real, dimension ( 0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: dUx_dxi, dUx_deta,  dUz_dxi, dUz_deta
        real, dimension ( 0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: s0, EMat, Uxloc, Uzloc

        ! This subroutine is called outside the Newmark scheme, so the local
        ! displacements are store in Forces
        Uxloc =Elem%Forces (:,:,0)
        Uzloc = Elem%Forces (:,:,1)

        dUx_dxi = MATMUL ( hTprime, Uxloc)
        dUz_dxi = MATMUL ( hTprime, Uzloc )
        dUx_deta = MATMUL ( Uxloc , hprimez )
        dUz_deta = MATMUL ( Uzloc , hprimez )

         s0 =  Elem%Acoeff(:,:,0)*dUx_dxi + Elem%Acoeff(:,:,1)*dUx_deta &
             + Elem%Acoeff(:,:,2)*dUz_dxi + Elem%Acoeff(:,:,3)*dUz_deta
         EMat = dUx_dxi * s0

         s0 =  Elem%Acoeff(:,:,2)*dUx_dxi + Elem%Acoeff(:,:,5)*dUx_deta &
             + Elem%Acoeff(:,:,7)*dUz_dxi + Elem%Acoeff(:,:,8)*dUz_deta
         EMat = EMat + dUz_dxi * s0

         s0 =  Elem%Acoeff(:,:,1)*dUx_dxi + Elem%Acoeff(:,:,4)*dUx_deta &
             + Elem%Acoeff(:,:,5)*dUz_dxi + Elem%Acoeff(:,:,6)*dUz_deta
         EMat = EMat + dUx_deta * s0

         s0 =  Elem%Acoeff(:,:,3)*dUx_dxi + Elem%Acoeff(:,:,6)*dUx_deta &
             + Elem%Acoeff(:,:,8)*dUz_dxi + Elem%Acoeff(:,:,9)*dUz_deta
         EMat = EMat + dUz_deta * s0

         E_elas = -0.5 * sum(EMat)

    end subroutine compute_Elastic_Energy


    ! ###########################################################
    !>
    !! \brief This subroutine computes the energy of the elastic deformation of the element.
    !!  Suitable for DG elements only.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (INOUT) E_elas
    !<

    subroutine  compute_Elastic_Energy_DG (Elem, E_elas)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (INOUT) :: E_elas
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: sigma11, sigma22, sigma12, EMat

        if (Elem%Acoustic) then
            Emat = Elem%Lambda * Elem%Strain(:,:,0) * Elem%Strain(:,:,0)
            EMat = EMat * Elem%Acoeff(:,:,4)
        else ! Elastic element
            sigma11 = (Elem%Lambda + 2*Elem%Mu) * Elem%Strain(:,:,0) + Elem%Lambda * Elem%Strain(:,:,1)
            sigma22 = (Elem%Lambda + 2*Elem%Mu) * Elem%Strain(:,:,1) + Elem%Lambda * Elem%Strain(:,:,0)
            sigma12 = 2*Elem%Mu * Elem%Strain(:,:,2)

            EMat = Elem%Strain(:,:,0)*sigma11 + 2.*Elem%Strain(:,:,2)*sigma12 &
                + Elem%Strain(:,:,1)*sigma22
            EMat = EMat * Elem%Acoeff(:,:,12)  ! Pour multiplier par Jac*wx*wz
        endif

        E_elas = 0.5 * sum(EMat)

  end subroutine compute_Elastic_Energy_DG


    ! ###########################################################
    !>
    !! \brief This subroutine computes the kinetic energy of the inner nodes
    !!  of an element
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (INOUT) E_kin
    !<

    subroutine  compute_Kinetic_Energy (Elem, Dt, E_kin)
        implicit none

        type (Element), intent (IN) :: Elem
        real, intent (IN)    :: Dt
        real, intent (INOUT) :: E_kin
        real, dimension (1:Elem%ngllx-2, 1:Elem%ngllz-2)      :: Ener_Mat
        real, dimension (1:Elem%ngllx-2, 1:Elem%ngllz-2, 0:1) :: Vel_half
        integer :: ngllx, ngllz

        ngllx = Elem%ngllx ; ngllz = Elem%ngllz

        Vel_half(:,:,:) = Elem%Veloc(:,:,:) + 0.5 * dt * Elem%Forces(1:ngllx-2,1:ngllz-2,:)
        Ener_Mat (:,:)  = 1./Elem%MassMat(:,:) * ( Vel_half(:,:,0)*Vel_half(:,:,0) &
                                                  +Vel_half(:,:,1)*Vel_half(:,:,1))
        E_kin = 0.5 * sum(Ener_Mat)

    end subroutine compute_Kinetic_Energy


    ! ###########################################################
    !>
    !! \brief This subroutine computes the kinetic energy using all the nodes
    !!  of an element. It suitable for Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (INOUT) E_kin
    !<

    subroutine  compute_Kinetic_Energy_DG (Elem, E_kin)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real, intent (INOUT) :: E_kin
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1)  :: Ener_Mat

        Ener_Mat = 1./Elem%MassMat(:,:) * ( Elem%Veloc(:,:,0)*Elem%Veloc(:,:,0) &
                                           +Elem%Veloc(:,:,1)*Elem%Veloc(:,:,1))
        E_kin = 0.5 * sum(Ener_Mat)

    end subroutine compute_Kinetic_Energy_DG

! ###########################################################
    !>
    !! \brief This subroutine redirects to subroutines computing the numerical Tractions
    !! according to the kind of media (acoustic or elastic).
    !! It suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_TracFace (Elem)
        implicit none

        type (Element), intent (INOUT)   :: Elem

        if (Elem%Acoustic) then
            call compute_TracFace_Acou (Elem)
        else
            call compute_TracFace_Elas (Elem)
        endif

    end subroutine compute_TracFace


! ###########################################################
    !>
    !! \brief This subroutine computes the numerical Tractions at the
    !! face of an element on each external node of the element. For ACOUSTIC media.
    !! It is suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_TracFace_Acou (Elem)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        integer    :: imin, imax, ngx, ngz
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! For the Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        Elem%TracFace(imin:imax,0) =  Elem%Lambda(0:ngx-1,0) * Elem%Strain(0:ngx-1,0,0) &
            - Elem%MatPen(imin:imax,0) * (Elem%Veloc(0:ngx-1,0,0) * Elem%Normal_nodes(imin:imax,0) &
                                       +  Elem%Veloc(0:ngx-1,0,1) * Elem%Normal_nodes(imin:imax,1))
        ! For the Right Face :
        call get_iminimax(Elem,1,imin,imax)
        Elem%TracFace(imin:imax,0) =  Elem%Lambda(ngx-1,0:ngz-1) * Elem%Strain(ngx-1,0:ngz-1,0) &
            - Elem%MatPen(imin:imax,0) * (Elem%Veloc(ngx-1,0:ngz-1,0) * Elem%Normal_nodes(imin:imax,0) &
                                       +  Elem%Veloc(ngx-1,0:ngz-1,1) * Elem%Normal_nodes(imin:imax,1))
        ! For the Top Face :
        call get_iminimax(Elem,2,imin,imax)
        Elem%TracFace(imin:imax,0) =  Elem%Lambda(0:ngx-1,ngz-1) * Elem%Strain(0:ngx-1,ngz-1,0) &
            - Elem%MatPen(imin:imax,0) * (Elem%Veloc(0:ngx-1,ngz-1,0) * Elem%Normal_nodes(imin:imax,0) &
                                       +  Elem%Veloc(0:ngx-1,ngz-1,1) * Elem%Normal_nodes(imin:imax,1))
        ! For the Left Face :
        call get_iminimax(Elem,3,imin,imax)
        Elem%TracFace(imin:imax,0) =  Elem%Lambda(0,0:ngz-1) * Elem%Strain(0,0:ngz-1,0) &
            - Elem%MatPen(imin:imax,0) * (Elem%Veloc(0,0:ngz-1,0) * Elem%Normal_nodes(imin:imax,0) &
                                       +  Elem%Veloc(0,0:ngz-1,1) * Elem%Normal_nodes(imin:imax,1))
    end subroutine compute_TracFace_Acou

! ###########################################################
    !>
    !! \brief This subroutine computes the numerical Tractions at the
    !! face of an element on each external node of the element. For ELASTIC media
    !! It is suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_TracFace_Elas (Elem)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        integer    :: imin, imax, ngx, ngz
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! For the Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        Elem%TracFace(imin:imax,0) = (Elem%Lambda(0:ngx-1,0) + 2*Elem%Mu(0:ngx-1,0)) &
                                    * Elem%Normal_nodes(imin:imax,0) * Elem%Strain(0:ngx-1,0,0) &
                                    + Elem%Lambda(0:ngx-1,0) * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0:ngx-1,0,1) &
                                    + 2*Elem%Mu(0:ngx-1,0)  * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0:ngx-1,0,2)
        Elem%TracFace(imin:imax,1) =  Elem%Lambda(0:ngx-1,0) * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0:ngx-1,0,0) &
                                    +(Elem%Lambda(0:ngx-1,0) + 2*Elem%Mu(0:ngx-1,0)) &
                                    * Elem%Normal_nodes(imin:imax,1) * Elem%Strain(0:ngx-1,0,1) &
                                    + 2*Elem%Mu(0:ngx-1,0)  * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0:ngx-1,0,2)
        Elem%TracFace(imin:imax,0) =  Elem%TracFace(imin:imax,0) &
                                    - Elem%MatPen(imin:imax,0) * Elem%Veloc(0:ngx-1,0,0) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,0,1)
        Elem%TracFace(imin:imax,1) =  Elem%TracFace(imin:imax,1) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,0,0) &
                                    - Elem%MatPen(imin:imax,1) * Elem%Veloc(0:ngx-1,0,1)
        ! For the Right Face :
        call get_iminimax(Elem,1,imin,imax)
        Elem%TracFace(imin:imax,0) = (Elem%Lambda(ngx-1,0:ngz-1) + 2*Elem%Mu(ngx-1,0:ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,0) * Elem%Strain(ngx-1,0:ngz-1,0) &
                                    + Elem%Lambda(ngx-1,0:ngz-1) * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(ngx-1,0:ngz-1,1) &
                                    + 2*Elem%Mu(ngx-1,0:ngz-1)  * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(ngx-1,0:ngz-1,2)
        Elem%TracFace(imin:imax,1) =  Elem%Lambda(ngx-1,0:ngz-1) * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(ngx-1,0:ngz-1,0) &
                                    +(Elem%Lambda(ngx-1,0:ngz-1) + 2*Elem%Mu(ngx-1,0:ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,1) * Elem%Strain(ngx-1,0:ngz-1,1) &
                                    + 2*Elem%Mu(ngx-1,0:ngz-1)  * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(ngx-1,0:ngz-1,2)
        Elem%TracFace(imin:imax,0) =  Elem%TracFace(imin:imax,0) &
                                    - Elem%MatPen(imin:imax,0) * Elem%Veloc(ngx-1,0:ngz-1,0) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(ngx-1,0:ngz-1,1)
        Elem%TracFace(imin:imax,1) =  Elem%TracFace(imin:imax,1) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(ngx-1,0:ngz-1,0) &
                                    - Elem%MatPen(imin:imax,1) * Elem%Veloc(ngx-1,0:ngz-1,1)
        ! For the Top Face :
        call get_iminimax(Elem,2,imin,imax)
        Elem%TracFace(imin:imax,0) = (Elem%Lambda(0:ngx-1,ngz-1) + 2*Elem%Mu(0:ngx-1,ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,0) * Elem%Strain(0:ngx-1,ngz-1,0) &
                                    + Elem%Lambda(0:ngx-1,ngz-1) * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0:ngx-1,ngz-1,1) &
                                    + 2*Elem%Mu(0:ngx-1,ngz-1)  * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0:ngx-1,ngz-1,2)
        Elem%TracFace(imin:imax,1) =  Elem%Lambda(0:ngx-1,ngz-1) * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0:ngx-1,ngz-1,0) &
                                    +(Elem%Lambda(0:ngx-1,ngz-1) + 2*Elem%Mu(0:ngx-1,ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,1) * Elem%Strain(0:ngx-1,ngz-1,1) &
                                    + 2*Elem%Mu(0:ngx-1,ngz-1)  * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0:ngx-1,ngz-1,2)
        Elem%TracFace(imin:imax,0) =  Elem%TracFace(imin:imax,0) &
                                    - Elem%MatPen(imin:imax,0) * Elem%Veloc(0:ngx-1,ngz-1,0) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,ngz-1,1)
        Elem%TracFace(imin:imax,1) =  Elem%TracFace(imin:imax,1) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,ngz-1,0) &
                                    - Elem%MatPen(imin:imax,1) * Elem%Veloc(0:ngx-1,ngz-1,1)
        ! For the Left Face :
        call get_iminimax(Elem,3,imin,imax)
        Elem%TracFace(imin:imax,0) = (Elem%Lambda(0,0:ngz-1) + 2*Elem%Mu(0,0:ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,0) * Elem%Strain(0,0:ngz-1,0) &
                                    + Elem%Lambda(0,0:ngz-1) * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0,0:ngz-1,1) &
                                    + 2*Elem%Mu(0,0:ngz-1)  * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0,0:ngz-1,2)
        Elem%TracFace(imin:imax,1) =  Elem%Lambda(0,0:ngz-1) * Elem%Normal_nodes(imin:imax,1) &
                                    * Elem%Strain(0,0:ngz-1,0) &
                                    +(Elem%Lambda(0,0:ngz-1) + 2*Elem%Mu(0,0:ngz-1)) &
                                    * Elem%Normal_nodes(imin:imax,1) * Elem%Strain(0,0:ngz-1,1) &
                                    + 2*Elem%Mu(0,0:ngz-1)  * Elem%Normal_nodes(imin:imax,0) &
                                    * Elem%Strain(0,0:ngz-1,2)
        Elem%TracFace(imin:imax,0) =  Elem%TracFace(imin:imax,0) &
                                    - Elem%MatPen(imin:imax,0) * Elem%Veloc(0,0:ngz-1,0) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0,0:ngz-1,1)
        Elem%TracFace(imin:imax,1) =  Elem%TracFace(imin:imax,1) &
                                    - Elem%MatPen(imin:imax,2) * Elem%Veloc(0,0:ngz-1,0) &
                                    - Elem%MatPen(imin:imax,1) * Elem%Veloc(0,0:ngz-1,1)
    end subroutine compute_TracFace_Elas

! ###########################################################
    !>
    !! \brief This subroutine computes the strain or pressure traces at the
    !! on each external node of the element, and add these strains to the
    !! Elem%Forces after weighting for integrals computation.
    !! It is suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_Traces (Elem)
        implicit none

        type (Element), intent (INOUT)   :: Elem

        if (Elem%Acoustic) then
            call compute_Traces_Acou (Elem)
        else
            call compute_Traces_Elas (Elem)
        endif

    end subroutine compute_Traces

! ###########################################################
    !>
    !! \brief This subroutine computes the pressure traces at the
    !! on each external node of an ACOUSTIC element, and add these strains to the
    !! Elem%Forces after weighting for integrals computation.
    !! It is suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_Traces_Acou (Elem)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        real, dimension(0:2*(Elem%ngllx+Elem%ngllz)-1) :: Vhatn
        integer    :: imin, imax, ngx, ngz
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! Adding contribution of Vhat on the trace Traction :
        Elem%TracFace(:,0) = Elem%TracFace(:,0) + Elem%Coeff_Integr_Faces(:)*Elem%MatPen(:,0)*Elem%Vhat(:,0)

        ! Preparing traces for the strain equation :
        Vhatn(:) = Elem%Vhat(:,0) * Elem%Coeff_Integr_Faces(:)

        ! Adding Strain and velocities traces to Elem%Forces :
        ! For the Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        Elem%Forces(0:ngx-1,0,0) = Elem%Forces(0:ngx-1,0,0) + Vhatn (imin:imax)
        Elem%Forces(0:ngx-1,0,1) = Elem%Forces(0:ngx-1,0,1) &
                                 + Elem%TracFace(imin:imax,0)*Elem%Normal_Nodes(imin:imax,0)
        Elem%Forces(0:ngx-1,0,2) = Elem%Forces(0:ngx-1,0,2) &
                                 + Elem%TracFace(imin:imax,0)*Elem%Normal_Nodes(imin:imax,1)
        ! For the Right Face :
        call get_iminimax(Elem,1,imin,imax)
        Elem%Forces(ngx-1,0:ngz-1,0) = Elem%Forces(ngx-1,0:ngz-1,0) + Vhatn (imin:imax)
        Elem%Forces(ngx-1,0:ngz-1,1) = Elem%Forces(ngx-1,0:ngz-1,1) &
                                     + Elem%TracFace(imin:imax,0)*Elem%Normal_Nodes(imin:imax,0)
        Elem%Forces(ngx-1,0:ngz-1,2) = Elem%Forces(ngx-1,0:ngz-1,2) &
                                     + Elem%TracFace(imin:imax,0)*Elem%Normal_Nodes(imin:imax,1)
        ! For the Top Face :
        call get_iminimax(Elem,2,imin,imax)
        Elem%Forces(0:ngx-1,ngz-1,0) = Elem%Forces(0:ngx-1,ngz-1,0) + Vhatn (imin:imax)
        Elem%Forces(0:ngx-1,ngz-1,1) = Elem%Forces(0:ngx-1,ngz-1,1) &
                                     + Elem%TracFace(imin:imax,0)*Elem%Normal_Nodes(imin:imax,0)
        Elem%Forces(0:ngx-1,ngz-1,2) = Elem%Forces(0:ngx-1,ngz-1,2) &
                                     + Elem%TracFace(imin:imax,0)*Elem%Normal_Nodes(imin:imax,1)
        ! For the Left Face :
        call get_iminimax(Elem,3,imin,imax)
        Elem%Forces(0,0:ngz-1,0) = Elem%Forces(0,0:ngz-1,0) + Vhatn (imin:imax)
        Elem%Forces(0,0:ngz-1,1) = Elem%Forces(0,0:ngz-1,1) &
                                 + Elem%TracFace(imin:imax,0)*Elem%Normal_Nodes(imin:imax,0)
        Elem%Forces(0,0:ngz-1,2) = Elem%Forces(0,0:ngz-1,2) &
                                 + Elem%TracFace(imin:imax,0)*Elem%Normal_Nodes(imin:imax,1)

    end subroutine compute_Traces_Acou

! ###########################################################
    !>
    !! \brief This subroutine computes the strain traces at the
    !! on each external node of an ELASTIC element, and add these strains to the
    !! Elem%Forces after weighting for integrals computation.
    !! It is suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_Traces_Elas (Elem)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        real, dimension(0:2*(Elem%ngllx+Elem%ngllz)-1) :: Vhatx,Vhatz,Vhatxz
        integer    :: imin, imax, ngx, ngz
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! Adding contribution of Vhat on the trace Traction :
        Elem%TracFace(:,0) = Elem%TracFace(:,0) + Elem%Coeff_Integr_Faces(:) * &
                            (Elem%MatPen(:,0)*Elem%Vhat(:,0) + Elem%MatPen(:,2)*Elem%Vhat(:,1))
        Elem%TracFace(:,1) = Elem%TracFace(:,1) + Elem%Coeff_Integr_Faces(:) * &
                            (Elem%MatPen(:,2)*Elem%Vhat(:,0) + Elem%MatPen(:,1)*Elem%Vhat(:,1))

        ! Preparing traces for the strain equation :
        Vhatxz(:) = 0.5 * ( Elem%Vhat(:,0) * Elem%Normal_Nodes(:,1) &
                          + Elem%Vhat(:,1) * Elem%Normal_Nodes(:,0))* Elem%Coeff_Integr_Faces(:)
        Vhatx(:) = Elem%Vhat(:,0) * Elem%Normal_Nodes(:,0) * Elem%Coeff_Integr_Faces(:)
        Vhatz(:) = Elem%Vhat(:,1) * Elem%Normal_Nodes(:,1) * Elem%Coeff_Integr_Faces(:)

        ! Adding Strain and velocities traces to Elem%Forces :
        ! For the Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        Elem%Forces(0:ngx-1,0,0) = Elem%Forces(0:ngx-1,0,0) + Vhatx (imin:imax)
        Elem%Forces(0:ngx-1,0,1) = Elem%Forces(0:ngx-1,0,1) + Vhatz (imin:imax)
        Elem%Forces(0:ngx-1,0,2) = Elem%Forces(0:ngx-1,0,2) + Vhatxz(imin:imax)
        Elem%Forces(0:ngx-1,0,3) = Elem%Forces(0:ngx-1,0,3) + Elem%TracFace(imin:imax,0)
        Elem%Forces(0:ngx-1,0,4) = Elem%Forces(0:ngx-1,0,4) + Elem%TracFace(imin:imax,1)
        ! For the Right Face :
        call get_iminimax(Elem,1,imin,imax)
        Elem%Forces(ngx-1,0:ngz-1,0) = Elem%Forces(ngx-1,0:ngz-1,0) + Vhatx (imin:imax)
        Elem%Forces(ngx-1,0:ngz-1,1) = Elem%Forces(ngx-1,0:ngz-1,1) + Vhatz (imin:imax)
        Elem%Forces(ngx-1,0:ngz-1,2) = Elem%Forces(ngx-1,0:ngz-1,2) + Vhatxz(imin:imax)
        Elem%Forces(ngx-1,0:ngz-1,3) = Elem%Forces(ngx-1,0:ngz-1,3) + Elem%TracFace(imin:imax,0)
        Elem%Forces(ngx-1,0:ngz-1,4) = Elem%Forces(ngx-1,0:ngz-1,4) + Elem%TracFace(imin:imax,1)
        ! For the Top Face :
        call get_iminimax(Elem,2,imin,imax)
        Elem%Forces(0:ngx-1,ngz-1,0) = Elem%Forces(0:ngx-1,ngz-1,0) + Vhatx (imin:imax)
        Elem%Forces(0:ngx-1,ngz-1,1) = Elem%Forces(0:ngx-1,ngz-1,1) + Vhatz (imin:imax)
        Elem%Forces(0:ngx-1,ngz-1,2) = Elem%Forces(0:ngx-1,ngz-1,2) + Vhatxz(imin:imax)
        Elem%Forces(0:ngx-1,ngz-1,3) = Elem%Forces(0:ngx-1,ngz-1,3) + Elem%TracFace(imin:imax,0)
        Elem%Forces(0:ngx-1,ngz-1,4) = Elem%Forces(0:ngx-1,ngz-1,4) + Elem%TracFace(imin:imax,1)
        ! For the Left Face :
        call get_iminimax(Elem,3,imin,imax)
        Elem%Forces(0,0:ngz-1,0) = Elem%Forces(0,0:ngz-1,0) + Vhatx (imin:imax)
        Elem%Forces(0,0:ngz-1,1) = Elem%Forces(0,0:ngz-1,1) + Vhatz (imin:imax)
        Elem%Forces(0,0:ngz-1,2) = Elem%Forces(0,0:ngz-1,2) + Vhatxz(imin:imax)
        Elem%Forces(0,0:ngz-1,3) = Elem%Forces(0,0:ngz-1,3) + Elem%TracFace(imin:imax,0)
        Elem%Forces(0,0:ngz-1,4) = Elem%Forces(0,0:ngz-1,4) + Elem%TracFace(imin:imax,1)

    end subroutine compute_Traces_Elas

    ! ###########################################################
    !>
    !! \brief This subroutine adds to the second member (i.e the forces) of the
    !!  system we solve, the terms comming from Int(tau*v*w dS)
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  add_tau_v (Elem)
        implicit none

        type (Element), intent (INOUT) :: Elem
        integer                        :: ngx, ngz, imin, imax
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        call get_iminimax(Elem,0,imin,imax)
        Elem%Forces(0:ngx-1,0,3) = Elem%Forces(0:ngx-1,0,3) - Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,0) * Elem%Veloc(0:ngx-1,0,0) &
                                 + Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,0,1))
        Elem%Forces(0:ngx-1,0,4) = Elem%Forces(0:ngx-1,0,4) - Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,0,0) &
                                 + Elem%MatPen(imin:imax,1) * Elem%Veloc(0:ngx-1,0,1))

        call get_iminimax(Elem,1,imin,imax)
        Elem%Forces(ngx-1,0:ngz-1,3) = Elem%Forces(ngx-1,0:ngz-1,3) - Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,0) * Elem%Veloc(ngx-1,0:ngz-1,0) &
                                 + Elem%MatPen(imin:imax,2) * Elem%Veloc(ngx-1,0:ngz-1,1))
        Elem%Forces(ngx-1,0:ngz-1,4) = Elem%Forces(ngx-1,0:ngz-1,4) - Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,2) * Elem%Veloc(ngx-1,0:ngz-1,0) &
                                 + Elem%MatPen(imin:imax,1) * Elem%Veloc(ngx-1,0:ngz-1,1))

        call get_iminimax(Elem,2,imin,imax)
        Elem%Forces(0:ngx-1,ngz-1,3) = Elem%Forces(0:ngx-1,ngz-1,3) - Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,0) * Elem%Veloc(0:ngx-1,ngz-1,0) &
                                 + Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,ngz-1,1))
        Elem%Forces(0:ngx-1,ngz-1,4) = Elem%Forces(0:ngx-1,ngz-1,4) - Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,2) * Elem%Veloc(0:ngx-1,ngz-1,0) &
                                 + Elem%MatPen(imin:imax,1) * Elem%Veloc(0:ngx-1,ngz-1,1))

        call get_iminimax(Elem,3,imin,imax)
        Elem%Forces(0,0:ngz-1,3) = Elem%Forces(0,0:ngz-1,3) - Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,0) * Elem%Veloc(0,0:ngz-1,0) &
                                 + Elem%MatPen(imin:imax,2) * Elem%Veloc(0,0:ngz-1,1))
        Elem%Forces(0,0:ngz-1,4) = Elem%Forces(0,0:ngz-1,4) - Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,2) * Elem%Veloc(0,0:ngz-1,0) &
                                 + Elem%MatPen(imin:imax,1) * Elem%Veloc(0,0:ngz-1,1))

    end subroutine add_tau_v

        ! ###########################################################
    !>
    !! \brief This subroutine adds to the second member (i.e the forces) of the
    !!  system we solve, the terms comming from Int(tau*v*w dS)
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  add_tau_v0 (Elem)
        implicit none

        type (Element), intent (INOUT) :: Elem
        integer                        :: ngx, ngz, imin, imax
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        call get_iminimax(Elem,0,imin,imax)
        Elem%Forces(0:ngx-1,0,3) = Elem%Forces(0:ngx-1,0,3) - 0.5*Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,0) * Elem%V0(0:ngx-1,0,0) &
                                 + Elem%MatPen(imin:imax,2) * Elem%V0(0:ngx-1,0,1))
        Elem%Forces(0:ngx-1,0,4) = Elem%Forces(0:ngx-1,0,4) - 0.5*Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,2) * Elem%V0(0:ngx-1,0,0) &
                                 + Elem%MatPen(imin:imax,1) * Elem%V0(0:ngx-1,0,1))

        call get_iminimax(Elem,1,imin,imax)
        Elem%Forces(ngx-1,0:ngz-1,3) = Elem%Forces(ngx-1,0:ngz-1,3) - 0.5*Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,0) * Elem%V0(ngx-1,0:ngz-1,0) &
                                 + Elem%MatPen(imin:imax,2) * Elem%V0(ngx-1,0:ngz-1,1))
        Elem%Forces(ngx-1,0:ngz-1,4) = Elem%Forces(ngx-1,0:ngz-1,4) - 0.5*Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,2) * Elem%V0(ngx-1,0:ngz-1,0) &
                                 + Elem%MatPen(imin:imax,1) * Elem%V0(ngx-1,0:ngz-1,1))

        call get_iminimax(Elem,2,imin,imax)
        Elem%Forces(0:ngx-1,ngz-1,3) = Elem%Forces(0:ngx-1,ngz-1,3) - 0.5*Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,0) * Elem%V0(0:ngx-1,ngz-1,0) &
                                 + Elem%MatPen(imin:imax,2) * Elem%V0(0:ngx-1,ngz-1,1))
        Elem%Forces(0:ngx-1,ngz-1,4) = Elem%Forces(0:ngx-1,ngz-1,4) - 0.5*Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,2) * Elem%V0(0:ngx-1,ngz-1,0) &
                                 + Elem%MatPen(imin:imax,1) * Elem%V0(0:ngx-1,ngz-1,1))

        call get_iminimax(Elem,3,imin,imax)
        Elem%Forces(0,0:ngz-1,3) = Elem%Forces(0,0:ngz-1,3) - 0.5*Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,0) * Elem%V0(0,0:ngz-1,0) &
                                 + Elem%MatPen(imin:imax,2) * Elem%V0(0,0:ngz-1,1))
        Elem%Forces(0,0:ngz-1,4) = Elem%Forces(0,0:ngz-1,4) - 0.5*Elem%Coeff_Integr_Faces(imin:imax) * &
                                 ( Elem%MatPen(imin:imax,2) * Elem%V0(0,0:ngz-1,0) &
                                 + Elem%MatPen(imin:imax,1) * Elem%V0(0,0:ngz-1,1))

    end subroutine add_tau_v0


    ! ###########################################################
    !>
    !! \brief This subroutine adds to the second member (i.e the forces) of the
    !!  system we solve, the terms comming from the previous time-step.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) Dt
    !!
    !<
    subroutine  add_previous_state2forces (Elem,Dt)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        real, intent (IN)                :: Dt

        if (Elem%acoustic) then
            Elem%Forces(:,:,0) = Elem%Forces(:,:,0) + 1./Dt * Elem%Acoeff(:,:,4) * Elem%Strain0(:,:,0)
            Elem%Forces(:,:,1) = Elem%Forces(:,:,1) + 1./Dt * Elem%Acoeff(:,:,4) * Elem%Density(:,:)*Elem%V0(:,:,0)
            Elem%Forces(:,:,2) = Elem%Forces(:,:,2) + 1./Dt * Elem%Acoeff(:,:,4) * Elem%Density(:,:)*Elem%V0(:,:,1)
        else ! Elastic Case
            Elem%Forces(:,:,0) = Elem%Forces(:,:,0) + 1./Dt * Elem%Acoeff(:,:,12) * Elem%Strain0(:,:,0)
            Elem%Forces(:,:,1) = Elem%Forces(:,:,1) + 1./Dt * Elem%Acoeff(:,:,12) * Elem%Strain0(:,:,1)
            Elem%Forces(:,:,2) = Elem%Forces(:,:,2) + 1./Dt * Elem%Acoeff(:,:,12) * Elem%Strain0(:,:,2)
            Elem%Forces(:,:,3) = Elem%Forces(:,:,3) + 1./Dt * Elem%Acoeff(:,:,12) * Elem%Density(:,:)*Elem%V0(:,:,0)
            Elem%Forces(:,:,4) = Elem%Forces(:,:,4) + 1./Dt * Elem%Acoeff(:,:,12) * Elem%Density(:,:)*Elem%V0(:,:,1)
        endif

    end subroutine add_previous_state2forces

    ! ###########################################################
    !>
    !! \brief This subroutine updates Strains and velocities from forces in an
    !! explicit forward Euler time scheme.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) Dt
    !!
    !<
    subroutine  update_Elem_Forward_Euler (Elem,Dt)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        real, intent (IN)                :: Dt

        if (Elem%acoustic) then
            Elem%Strain(:,:,0) = Dt * Elem%Forces(:,:,0)
            Elem%Veloc (:,:,:) = Dt * Elem%Forces(:,:,1:2)
        else ! Elastic Case
            Elem%Strain(:,:,:) = Dt * Elem%Forces(:,:,0:2)
            Elem%Veloc (:,:,:) = Dt * Elem%Forces(:,:,3:4)
        endif

    end subroutine update_Elem_Forward_Euler

    ! ###########################################################
    !>
    !! \brief This subroutine computes the second members "R" of the system
    !! K*Lambda = R on the Lagrange multiplicators (= Vhat) living on the faces.
    !! R is the computed on each external node of the current element, and is stored
    !! on the variable Elem%TracFace because indeed, R is homogeneous to a traction.
    !! This Subroutine is suitable for Hybridizable Discontinuous Galerkin elements
    !! only, and only in a framework of semi-implicit time-scheme.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) Dt
    !!
    !<
    subroutine  compute_smbr_R (Elem,Dt)
        implicit none

        type (Element), intent (INOUT)   :: Elem
        real, intent (IN) :: Dt
        real, dimension(0:2*(Elem%ngllx+Elem%ngllz)-1,0:4) :: smbr
        integer    :: imin, imax, ngx, ngz
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! Storing the Forces in the array smbr :
        call get_iminimax(Elem,0,imin,imax)
        smbr(imin:imax,0:4) = Elem%Forces(0:ngx-1,0,0:4)
        call get_iminimax(Elem,1,imin,imax)
        smbr(imin:imax,0:4) = Elem%Forces(ngx-1,0:ngz-1,0:4)
        call get_iminimax(Elem,2,imin,imax)
        smbr(imin:imax,0:4) = Elem%Forces(0:ngx-1,ngz-1,0:4)
        call get_iminimax(Elem,3,imin,imax)
        smbr(imin:imax,0:4) = Elem%Forces(0,0:ngz-1,0:4)

        ! Building R = L - CAinv * Forces_strain + EDinv * Forces_Veloc
        ! L is for Neumann boundary conditions. So far L = 0.
        Elem%TracFace(:,0) = -Elem%CAinv(:,0,0)*smbr(:,0)-Elem%CAinv(:,0,1)*smbr(:,1)-Elem%CAinv(:,0,2)*smbr(:,2)&
                             +Elem%EDinv(:,0,0)*smbr(:,3)+Elem%EDinv(:,0,1)*smbr(:,4)
        Elem%TracFace(:,1) = -Elem%CAinv(:,1,0)*smbr(:,0)-Elem%CAinv(:,1,1)*smbr(:,1)-Elem%CAinv(:,1,2)*smbr(:,2)&
                             +Elem%EDinv(:,1,0)*smbr(:,3)+Elem%EDinv(:,1,1)*smbr(:,4)
        Elem%TracFace(:,:) = Dt * Elem%TracFace(:,:)

    end subroutine compute_smbr_R

    ! ###########################################################
    !>
    !! \brief This subroutine performs the local inversion :
    !! | A  0 | (Strain)     ( Smbr_e )
    !! |      | (      )  =  (        )
    !! | 0  D | (Veloc )     ( Smbr_v )
    !! This Subroutine is suitable for Hybridizable Discontinuous Galerkin elements
    !! only, and only in a framework of semi-implicit time-scheme.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) Dt
    !!
    !<
    subroutine  inversion_local_solver (Elem,Dt)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real,           intent (IN)    :: Dt
        real, dimension(0:1)           :: aux
        !real, dimension(0:2*(Elem%ngllx+Elem%ngllz)-1,0:1) :: smbr, res
        integer                        :: ngx, ngz, imin, imax, i
        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! Inversion for strains :
        Elem%Forces(:,:,0) = Dt / Elem%Acoeff(:,:,12) * Elem%Forces(:,:,0)
        Elem%Forces(:,:,1) = Dt / Elem%Acoeff(:,:,12) * Elem%Forces(:,:,1)
        Elem%Forces(:,:,2) = Dt / Elem%Acoeff(:,:,12) * Elem%Forces(:,:,2)
        ! Inversion for velocities :
!        Elem%Forces(:,:,3) = Dt * Elem%MassMat(:,:) * Elem%Forces(:,:,3)
!        Elem%Forces(:,:,4) = Dt * Elem%MassMat(:,:) * Elem%Forces(:,:,4)

        ! FOLLOWING PART USED IF THE TERM int(tau*v*w) IS TAKEN IMPLICIT
        ! Inversion for velocities at inner GLL nodes:
        Elem%Forces(1:ngx-2,1:ngz-2,3) = Dt * Elem%MassMat(1:ngx-2,1:ngz-2) * Elem%Forces(1:ngx-2,1:ngz-2,3)
        Elem%Forces(1:ngx-2,1:ngz-2,4) = Dt * Elem%MassMat(1:ngx-2,1:ngz-2) * Elem%Forces(1:ngx-2,1:ngz-2,4)

        ! Inversion for velocity at outer GLL nodes :
        ! Bottom Face :
        do i=0,ngx-1
            aux = MATMUL(Elem%Jmat(i,:,:),Elem%Forces(i,0,0:2))
            Elem%Forces(i,0,3) = Dt * ( Elem%Dinv(i,0)*Elem%Forces(i,0,3) &
                                      + Elem%Dinv(i,2)*Elem%Forces(i,0,4))&
                                      + aux(0)
            Elem%Forces(i,0,4) = Dt * ( Elem%Dinv(i,2)*Elem%Forces(i,0,3) &
                                      + Elem%Dinv(i,1)*Elem%Forces(i,0,4))&
                                      + aux(1)
        enddo
        ! Right Face :
        imin = ngx ; imax = ngx+ngz-2
        do i=1,ngz-2
            aux = MATMUL(Elem%Jmat(imin+i,:,:),Elem%Forces(ngx-1,i,0:2))
            Elem%Forces(ngx-1,i,3) = Dt * ( Elem%Dinv(imin+i,0)*Elem%Forces(ngx-1,i,3) &
                                          + Elem%Dinv(imin+i,2)*Elem%Forces(ngx-1,i,4))&
                                          + aux(0)
            Elem%Forces(ngx-1,i,4) = Dt * ( Elem%Dinv(imin+i,2)*Elem%Forces(ngx-1,i,3) &
                                          + Elem%Dinv(imin+i,1)*Elem%Forces(ngx-1,i,4))&
                                          + aux(1)
        enddo
        ! Top Face :
        imin = ngx+ngz ; imax = 2*ngx+ngz-1
        do i=0,ngx-1
            aux = MATMUL(Elem%Jmat(imin+i,:,:),Elem%Forces(i,ngz-1,0:2))
            Elem%Forces(i,ngz-1,3) = Dt * ( Elem%Dinv(imin+i,0)*Elem%Forces(i,ngz-1,3) &
                                          + Elem%Dinv(imin+i,2)*Elem%Forces(i,ngz-1,4))&
                                          + aux(0)
            Elem%Forces(i,ngz-1,4) = Dt * ( Elem%Dinv(imin+i,2)*Elem%Forces(i,ngz-1,3) &
                                          + Elem%Dinv(imin+i,1)*Elem%Forces(i,ngz-1,4))&
                                          + aux(1)
        enddo
        ! Left Face :
        imin = 2*ngx+ngz ; imax = 2*ngx+2*ngz-2
        do i=1,ngz-2
            aux = MATMUL(Elem%Jmat(imin+i,:,:),Elem%Forces(0,i,0:2))
            Elem%Forces(0,i,3) = Dt * ( Elem%Dinv(imin+i,0)*Elem%Forces(0,i,3) &
                                      + Elem%Dinv(imin+i,2)*Elem%Forces(0,i,4))&
                                      + aux(0)
            Elem%Forces(0,i,4) = Dt * ( Elem%Dinv(imin+i,2)*Elem%Forces(0,i,3) &
                                      + Elem%Dinv(imin+i,1)*Elem%Forces(0,i,4))&
                                      + aux(1)
        enddo

    end subroutine inversion_local_solver


    ! ###########################################################
    !>
    !! \brief This subroutine computes computes the local fields of the current element
    !! Elem%Veloc, and Elem%Strain from the previous state V0 and Strain0, and from
    !! the traces Vhat (= Lagrange multiplicators).
    !! This Subroutine is suitable for Hybridizable Discontinuous Galerkin elements
    !! only, and only in a framework of semi-implicit time-scheme.
    !! \param type (Element), intent (INOUT) Elem
    !! \param real, intent (IN) Dt
    !!
    !<
    subroutine  local_solver (Elem, Dt)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real,           intent (IN)    :: Dt

        ! First step : traces terms are computed using the subroutine compute_traces
        ! and setting the former tractions TracFace to 0
        Elem%TracFace = 0.
        call compute_Traces (Elem)

        ! second step : the Mass matrices are inverted
        call inversion_local_solver(Elem, Dt)

        ! Last step : updates of Velocities and strains
        !Elem%Forces(:,:,:) = Dt * Elem%Forces(:,:,:)
        Elem%Strain(:,:,:) = Elem%Forces(:,:,0:2)
        Elem%Veloc (:,:,:) = Elem%Forces(:,:,3:4)

    end subroutine local_solver

    ! ###########################################################
    !>
    !!\brief Subroutine that computes the "imin" and "imax" indexes  which corresponds to the
    !! begining and the end of the arrays that contain data for the local face "w_face" of element "Elem".
    !!\version 1.0
    !!\date 03/04/2014
    !! This subroutine is used with DG and HDG elements.
    !<
    subroutine get_iminimax (Elem,w_face,imin,imax)

        implicit none

        type(element), intent(IN) :: Elem
        integer, intent (IN)      :: w_face
        integer, intent(INOUT)    :: imin
        integer, intent(INOUT)    :: imax

        select case (w_face)
        case(0)
            imin = 0
            imax = Elem%ngllx-1
        case(1)
            imin = Elem%ngllx
            imax = Elem%ngllx   + Elem%ngllz-1
        case(2)
            imin = Elem%ngllx   + Elem%ngllz
            imax = 2*Elem%ngllx + Elem%ngllz-1
        case(3)
            imin = 2*Elem%ngllx + Elem%ngllz
            imax = 2*Elem%ngllx + 2*Elem%ngllz -1
        end select

    end subroutine get_iminimax


    ! ###########################################################
    !>
    !! \brief This subroutine gets the local numerotation of gll nodes
    !! around a given corner nc of an element Elem.
    !! Note that the "corner" nc corresponds to the vertex Near_Vertx(nc)
    !! \param type (element), intent (INOUT) Elem
    !! \param integer, intent (IN)  nc
    !! \param integer, intent (OUt) n1, n2
    !<
    subroutine get_gll_arround_corner(Elem,nc,n1,n2)

        implicit none
        type (Element), intent (IN) :: Elem
        integer, intent(IN)  :: nc
        integer, intent(OUT) :: n1, n2

        select case (nc)
        case(0)
            n1 = 2*Elem%ngllx + Elem%ngllz
            n2 = 0
        case(1)
            n1 = Elem%ngllx-1
            n2 = Elem%ngllx
        case(2)
            n1 = Elem%ngllx+Elem%ngllz-1
            n2 = 2*Elem%ngllx+Elem%ngllz-1
        case(3)
            n1 = Elem%ngllx+Elem%ngllz
            n2 = 2*(Elem%ngllx+Elem%ngllz)-1
        end select

    end subroutine get_gll_arround_corner

! ###########################################################

end module selement

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
