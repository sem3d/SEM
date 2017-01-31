
!>
!!\file FPML_bak.F90
!!\brief This file contains all the subroutines needed for the FPML
!! (Frequency-Dependant PML).
!!\version 1.0
!!\date 20/03/2015
!! This module is neither implemented, neither called. Here is some
!! dead code that should be used in case FPML should be re-implemented.
!! Before each subroutine, a comment states in which file the subroutine
!! should be re-introduced.
!<

module FPML
    use sdomain
    use constants

    implicit none
contains


!###########################################################!
!! Morceaux de code a rajouter dans allocate_domain.F90 :
!###########################################################!
    if (Tdomain%specel(n)%FPML) then
        allocate (Tdomain%specel(n)%Isx(0:ngllx-1,0:ngllz-1))
        allocate (Tdomain%specel(n)%Isz(0:ngllx-1,0:ngllz-1))
        allocate (Tdomain%specel(n)%Ivx(0:ngllx-1,0:ngllz-1))
        allocate (Tdomain%specel(n)%Ivz(0:ngllx-1,0:ngllz-1))
        allocate (Tdomain%specel(n)%IStress1(0:ngllx-1,0:ngllz-1,0:2))
        allocate (Tdomain%specel(n)%IStress2(0:ngllx-1,0:ngllz-1,0:2))
        allocate (Tdomain%specel(n)%IVeloc1(1:ngllx-2,1:ngllz-2,0:1))
        allocate (Tdomain%specel(n)%IVeloc2(1:ngllx-2,1:ngllz-2,0:1))
        allocate (Tdomain%specel(n)%DumpMass(0:ngllx-1,0:ngllz-1,0:3))
        Tdomain%specel(n)%IStress1 = 0.
        Tdomain%specel(n)%IStress2 = 0.
        Tdomain%specel(n)%IVeloc1 = 0.
        Tdomain%specel(n)%IVeloc2 = 0.
    endif

    if (Tdomain%sFace(n)%FPML) then
        allocate (Tdomain%sFace(n)%Ivx(1:ngll-2))
        allocate (Tdomain%sFace(n)%Ivz(1:ngll-2))
        allocate (Tdomain%sFace(n)%IVeloc1(1:ngll-2,0:1))
        allocate (Tdomain%sFace(n)%IVeloc2(1:ngll-2,0:1))
        allocate (Tdomain%sFace(n)%DumpMass(1:ngll-2,0:3))
        Tdomain%sFace(n)%IVeloc1 = 0.
        Tdomain%sFace(n)%IVeloc2 = 0.
        Tdomain%sFace(n)%Ivx = 0.
        Tdomain%sFace(n)%Ivz = 0.
    endif

    if (Tdomain%sFace(n)%FPML) then
        allocate (Tdomain%sVertex(n)%Ivx(0:0))
        allocate (Tdomain%sVertex(n)%Ivz(0:0))
        allocate (Tdomain%sVertex(n)%IVeloc1(0:1))
        allocate (Tdomain%sVertex(n)%IVeloc2(0:1))
        allocate (Tdomain%sVertex(n)%DumpMass(0:1))
        Tdomain%sVertex(n)%IVeloc1 = 0.
        Tdomain%sVertex(n)%IVeloc2 = 0.
    endif


!###########################################################!
!! Morceaux de code a rajouter dans Element.F90 :
!###########################################################!
       ! FPML allocation
       logical :: FPML
       real(fpp), dimension (:,:), allocatable ::  Isx, Isz, Ivx, Ivz
       real(fpp), dimension (:,:,:), allocatable :: Istress1, IStress2, Iveloc1, Iveloc2


    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real(fpp), intent (IN) dt
    !! \param real(fpp), intent (IN) fil
    !<
    subroutine Correction_Elem_FPML_Veloc (Elem, dt, fil)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real(fpp), intent (IN) ::  dt, fil

        integer   :: ngllx, ngllz,i
        real(fpp) :: fil2
        real(fpp), dimension (1:Elem%ngllx-2,1:Elem%ngllz-2) :: Ausiliar_Velocity

        ngllx = Elem%ngllx; ngllz=Elem%ngllz
        fil2 = fil**2

        do i = 0,1
            Ausiliar_Velocity = Elem%Veloc1(:,:,i)
            Elem%Veloc1(:,:,i) = Elem%DumpVx(:,:,0) * Elem%Veloc1(:,:,i) + dt * &
                Elem%DumpVx(:,:,1)*Elem%Forces1(1:ngllx-2,1:ngllz-2,i) + Elem%Ivx * Elem%Iveloc1(:,:,i)
            Elem%Iveloc1(:,:,i) = Fil2 * Elem%Iveloc1(:,:,i) + 0.5 * (1.-Fil2) *  &
                (Ausiliar_Velocity + Elem%Veloc1(:,:,i) )

            Ausiliar_Velocity = Elem%Veloc2(:,:,i)
            Elem%Veloc2(:,:,i) =Elem%DumpVz(:,:,0) * Elem%Veloc2(:,:,i) + Dt * &
                Elem%DumpVz(:,:,1)*Elem%Forces2(1:ngllx-2,1:ngllz-2,i) + Elem%Ivz * Elem%IVeloc2(:,:,i)
            Elem%Iveloc2(:,:,i) = Fil2 * Elem%Iveloc2(:,:,i) + 0.5 * (1.-Fil2) * &
                (Ausiliar_Velocity + Elem%Veloc2(:,:,i) )
        enddo

        Elem%V0 = Elem%Veloc
        Elem%Veloc = Elem%Veloc1 + Elem%Veloc2

        return
    end subroutine Correction_Elem_FPML_Veloc


    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Element), intent (INOUT) Elem
    !! \param real(fpp), dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) hTmat
    !! \param real(fpp), dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) hmatz
    !! \param real(fpp), dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vxloc
    !! \param real(fpp), dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT) Vzloc
    !! \param real(fpp), intent (IN) bega
    !! \param real(fpp), intent (IN) dt
    !! \param real(fpp), intent (IN) alpha
    !! \param real(fpp), intent (IN) fil
    !<
    subroutine Prediction_Elem_FPML_Veloc (Elem,alpha, bega, dt,Vxloc,Vzloc,Hmatz, HTmat,fil)
        implicit none

        type (Element), intent (INOUT) :: Elem
        real(fpp), dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) ::  hTmat
        real(fpp), dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hmatz
        real(fpp), dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1), intent (INOUT)  ::Vxloc, Vzloc
        real(fpp), intent (IN) :: bega, dt, alpha, fil

        real(fpp) :: fil2
        real(fpp), dimension (0:Elem%ngllx-1, 0:Elem%ngllz-1) :: s0,s1,s2,s3, Stress_ausiliar

        integer :: ngllx, ngllz

        ngllx = Elem%ngllx; ngllz = Elem%ngllz
        fil2 = fil**2

        !Elem%Stress0 = Elem%Stress

        VxLoc(1:ngllx-2,1:ngllz-2)  = (0.5+alpha) * Elem%Veloc(:,:,0) + dt *(0.5-bega) *Elem%Accel(:,:,0) + (0.5-alpha) * Elem%V0(:,:,0)
        VzLoc(1:ngllx-2,1:ngllz-2)  =(0.5+alpha) * Elem%Veloc(:,:,1) + dt *(0.5-bega) *Elem%Accel(:,:,1)+ (0.5-alpha) * Elem%V0(:,:,1)

        s0 = MATMUL (HTmat,VxLoc)
        s2 = MATMUL (HTmat,VzLoc)
        s1 = MATMUL (VxLoc,Hmatz)
        s3 = MATMUL (VzLoc,Hmatz)

        Stress_Ausiliar = Elem%Stress1(:,:,0)
        Elem%Stress1(:,:,0) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,0) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,0) * s0 + &
            Elem%Acoeff(:,:,2) * s2 ) + Elem%Isx * Elem%Istress1(:,:,0)
        Elem%Istress1 (:,:,0) =  Fil2 * Elem%Istress1 (:,:,0) + 0.5 * (1.-Fil2) * (Elem%Stress1(:,:,0) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress2(:,:,0)
        Elem%Stress2(:,:,0) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,0) + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,3) * s3 + &
            Elem%Acoeff(:,:,1)* s1 )  + Elem%Isz * Elem%Istress2(:,:,0)
        Elem%Istress2 (:,:,0) =  Fil2 * Elem%Istress2 (:,:,0) + 0.5 * (1.-Fil2) * (Elem%Stress2(:,:,0) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress1(:,:,1)
        Elem%Stress1(:,:,1) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,1) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,4) * s0 + &
            Elem%Acoeff(:,:,6) * s2 ) + Elem%Isx * Elem%Istress1(:,:,1)
        Elem%Istress1 (:,:,1) =  Fil2 * Elem%Istress1 (:,:,1) + 0.5 * (1.-Fil2) * (Elem%Stress1(:,:,1) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress2(:,:,1)
        Elem%Stress2(:,:,1) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,1) + Elem%DumpSz(:,:,1) * Dt *( Elem%Acoeff(:,:,7) * s3 + &
            Elem%Acoeff(:,:,5) * s1 ) + Elem%Isz * Elem%Istress2(:,:,1)
        Elem%Istress2 (:,:,1) =  Fil2 * Elem%Istress2 (:,:,1) + 0.5 * (1.-Fil2) * (Elem%Stress2(:,:,1) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress1(:,:,2)
        Elem%Stress1(:,:,2) = Elem%DumpSx(:,:,0) * Elem%Stress1(:,:,2) + Elem%DumpSx(:,:,1) * Dt * (Elem%Acoeff(:,:,10) * s2 + &
            Elem%Acoeff(:,:,8) * s0) + Elem%Isx * Elem%Istress1(:,:,2)
        Elem%Istress1 (:,:,2) =  Fil2 * Elem%Istress1 (:,:,2) + 0.5 * (1.-Fil2) * (Elem%Stress1(:,:,2) + Stress_Ausiliar)

        Stress_Ausiliar = Elem%Stress2(:,:,2)
        Elem%Stress2(:,:,2) =  Elem%DumpSz(:,:,0) * Elem%Stress2(:,:,2) + Elem%DumpSz(:,:,1) * Dt * (Elem%Acoeff(:,:,9) * s1 + &
            Elem%Acoeff(:,:,11) * s3 ) + Elem%Isz * Elem%Istress2(:,:,2)
        Elem%Istress2 (:,:,2) =  Fil2 * Elem%Istress2 (:,:,2) + 0.5 * (1.-Fil2) * (Elem%Stress2(:,:,2) + Stress_Ausiliar)

        Elem%Stress = Elem%Stress1 + Elem%Stress2
        return
    end subroutine Prediction_Elem_FPML_Veloc


!###########################################################!
!! Morceaux de code a rajouter dans Face.F90 :
!###########################################################!

       real(fpp), dimension (:), allocatable :: Ivx, Ivz
       real(fpp), dimension (:,:), allocatable :: Iveloc1, Iveloc2

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Face toto), intent (INOUT) F
    !! \param real(fpp), intent (IN toto) dt
    !! \param real(fpp), intent (IN toto) fil
    !<


    subroutine Correction_Face_FPML_Veloc (F, dt, fil)
        implicit none

        type (Face), intent (INOUT) :: F
        real(fpp), intent (IN) ::  dt, fil

        integer   :: i
        real(fpp) :: fil2
        real(fpp), dimension (1:F%ngll-2) :: Ausiliar_Velocity


        fil2 = fil**2
        F%V0 = F%Veloc

        if  (F%Abs .or. F%reflex) then
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



!###########################################################!
!! Morceaux de code a rajouter dans Vertex.F90 :
!###########################################################!

    real(fpp), dimension (:), allocatable :: Ivx, Ivz,Iveloc1, Iveloc2

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real(fpp), intent (IN) dt
    !! \param real(fpp), intent (IN) fil
    !<
    subroutine Correction_Vertex_FPML_Veloc (V, dt, fil)
        implicit none

        type (Vertex), intent (INOUT) :: V
        real(fpp), intent (IN) ::  dt, fil

        integer   :: i
        real(fpp) :: fil2
        real(fpp) :: Ausiliar_Velocity

        fil2 = fil**2
        V%V0 = V%Veloc

        if  (V%Abs .or. V%reflex) then
            V%Veloc1 = 0; V%Veloc2 = 0; V%Veloc = 0
        else

            do i = 0,1
                Ausiliar_Velocity = V%Veloc1(i)
                V%Veloc1(i) = V%DumpVx(0) * V%Veloc1(i) + dt * &
                    V%DumpVx(1)*V%Forces1(i) + V%Ivx(0) * V%Iveloc1(i)
                V%Iveloc1(i) = Fil2 * V%Iveloc1(i) + 0.5 * (1.-Fil2) *  &
                    (Ausiliar_Velocity + V%Veloc1(i) )

                Ausiliar_Velocity = V%Veloc2(i)
                V%Veloc2(i) =V%DumpVz(0) * V%Veloc2(i) + Dt * &
                    V%DumpVz(1)*V%Forces2(i) + V%Ivz(0) * V%IVeloc2(i)
                V%Iveloc2(i) = Fil2 * V%Iveloc2(i) + 0.5 * (1.-Fil2) * &
                    (Ausiliar_Velocity + V%Veloc2(i) )
            enddo
            V%Veloc = V%Veloc1 + V%Veloc2
        endif

        return
    end subroutine Correction_Vertex_FPML_Veloc
    ! ###########################################################


    !###########################################################!
    !! Morceaux de code a rajouter dans define_arr.F90 :
    !###########################################################!

    if (Tdomain%specel(n)%FPML) then

        Tdomain%specel(n)%DumpSx(:,:,1) = (1.+Tdomain%sSubdomain(mat)%k) * Id(:,:) + 0.5 * Tdomain%sSubdomain(mat)%Dt * &
            wx(:,:) * Tdomain%sSubdomain(mat)%freq
        Tdomain%specel(n)%DumpSx (:,:,1) = 1./ Tdomain%specel(n)%DumpSx (:,:,1)
        Tdomain%specel(n)%DumpSx (:,:,0) = ((1+Tdomain%sSubdomain(mat)%k) * Id(:,:) - Tdomain%sSubdomain(mat)%Dt * 0.5 * &
            wx(:,:) * Tdomain%sSubdomain(mat)%freq) *  Tdomain%specel(n)%DumpSx(:,:,1)
        Tdomain%specel(n)%DumpMass(:,:,0) =  Tdomain%specel(n)%Density * Whei * Jac * wx * &
            ( Tdomain%sSubdomain(mat)%Dt * 0.5 * Tdomain%sSubdomain(mat)%freq + Tdomain%sSubdomain(mat)%k)
        Tdomain%specel(n)%DumpMass(:,:,1) =   Tdomain%specel(n)%Density * Whei *  Jac * wx * &
            (Tdomain%sSubdomain(mat)%k  - Tdomain%sSubdomain(mat)%Dt * 0.5 *  Tdomain%sSubdomain(mat)%freq)

        Tdomain%specel(n)%DumpSz(:,:,1) = (1+Tdomain%sSubdomain(mat)%k) * Id(:,:) + 0.5 * Tdomain%sSubdomain(mat)%Dt * &
            wz(:,:) * Tdomain%sSubdomain(mat)%freq
        Tdomain%specel(n)%DumpSz (:,:,1) = 1./ Tdomain%specel(n)%DumpSz (:,:,1)
        Tdomain%specel(n)%DumpSz (:,:,0) = ((1+Tdomain%sSubdomain(mat)%k) * Id(:,:) - Tdomain%sSubdomain(mat)%Dt * 0.5 *  &
            wz(:,:) * Tdomain%sSubdomain(mat)%freq) * Tdomain%specel(n)%DumpSz(:,:,1)
        Tdomain%specel(n)%DumpMass(:,:,2) =  Tdomain%specel(n)%Density * Whei *  Jac * wz * &
            ( Tdomain%sSubdomain(mat)%Dt * 0.5 * Tdomain%sSubdomain(mat)%freq + Tdomain%sSubdomain(mat)%k )
        Tdomain%specel(n)%DumpMass(:,:,3) = Tdomain%specel(n)%Density * Whei  * Jac * wz * &
            (Tdomain%sSubdomain(mat)%k  - Tdomain%sSubdomain(mat)%Dt * 0.5 *   Tdomain%sSubdomain(mat)%freq)

        Tdomain%specel(n)%Isx=Tdomain%sSubdomain(mat)%Dt*wx*Tdomain%sSubdomain(mat)%freq*Tdomain%specel(n)%DumpSx(:,:,1)
        Tdomain%specel(n)%Isz=Tdomain%sSubdomain(mat)%Dt*wz*Tdomain%sSubdomain(mat)%freq*Tdomain%specel(n)%DumpSz(:,:,1)

        Tdomain%specel(n)%Ivx=Tdomain%specel(n)%Density*Whei*Tdomain%sSubdomain(mat)%Dt*wx*Jac*Tdomain%sSubdomain(mat)%freq
        Tdomain%specel(n)%Ivz=Tdomain%specel(n)%Density*Whei*Tdomain%sSubdomain(mat)%Dt*wz*Jac*Tdomain%sSubdomain(mat)%freq
    endif

    ! Dans la boucle sur les faces
    if (Tdomain%sFace(nf)%FPML) call getIv_element2face (Tdomain,n_elem,nf,w_face,.true.)
    if (Tdomain%sFace(nf)%FPML ) call getIv_element2face (Tdomain,n_elem,nf,w_face,Tdomain%sFace(nf)%coherency)

    ! Boucles sur les vertexs :
    if (Tdomain%sVertex(nv_aus)%FPML) then
        Tdomain%sVertex(nv_aus)%Ivx(0) = Tdomain%sVertex(nv_aus)%Ivx(0)+ Tdomain%specel(n)%Ivx(0,0)
        Tdomain%sVertex(nv_aus)%Ivz(0) = Tdomain%sVertex(nv_aus)%Ivz(0) + Tdomain%specel(n)%Ivz(0,0)
    endif
    if (Tdomain%sVertex(nv_aus)%FPML) then
        Tdomain%sVertex(nv_aus)%Ivx(0) = Tdomain%sVertex(nv_aus)%Ivx(0)+ Tdomain%specel(n)%Ivx(ngllx-1,0)
        Tdomain%sVertex(nv_aus)%Ivz(0) = Tdomain%sVertex(nv_aus)%Ivz(0) + Tdomain%specel(n)%Ivz(ngllx-1,0)
    endif
    if (Tdomain%sVertex(nv_aus)%FPML) then
        Tdomain%sVertex(nv_aus)%Ivx(0) = Tdomain%sVertex(nv_aus)%Ivx(0)+ Tdomain%specel(n)%Ivx(ngllx-1,ngllz-1)
        Tdomain%sVertex(nv_aus)%Ivz(0) = Tdomain%sVertex(nv_aus)%Ivz(0) + Tdomain%specel(n)%Ivz(ngllx-1,ngllz-1)
    endif
    if (Tdomain%sVertex(nv_aus)%FPML) then
        Tdomain%sVertex(nv_aus)%Ivx(0) = Tdomain%sVertex(nv_aus)%Ivx(0)+ Tdomain%specel(n)%Ivx(0,ngllz-1)
        Tdomain%sVertex(nv_aus)%Ivz(0) = Tdomain%sVertex(nv_aus)%Ivz(0) + Tdomain%specel(n)%Ivz(0,ngllz-1)
    endif

    ! Invert Mass Matrix expression
    if (Tdomain%specel(n)%FPML) then
        LocMassMat(:,:) = Tdomain%specel(n)%MassMat(1:ngllx-2,1:ngllz-2)
        Tdomain%specel(n)%DumpVx (:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:ngllz-2,0)
        Tdomain%specel(n)%DumpVx (:,:,1) = 1./ Tdomain%specel(n)%DumpVx (:,:,1)
        Tdomain%specel(n)%DumpVx(:,:,0) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:ngllz-2,1)
        Tdomain%specel(n)%DumpVx (:,:,0) = Tdomain%specel(n)%DumpVx (:,:,0) * Tdomain%specel(n)%DumpVx (:,:,1)

        Tdomain%specel(n)%DumpVz (:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:ngllz-2,2)
        Tdomain%specel(n)%DumpVz (:,:,1) = 1./ Tdomain%specel(n)%DumpVz (:,:,1)
        Tdomain%specel(n)%DumpVz(:,:,0) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:ngllz-2,3)
        Tdomain%specel(n)%DumpVz (:,:,0) =    Tdomain%specel(n)%DumpVz (:,:,0) *   Tdomain%specel(n)%DumpVz (:,:,1)
        deallocate (Tdomain%specel(n)%MassMat) ; deallocate (Tdomain%specel(n)%DumpMass)
        allocate (Tdomain%specel(n)%MassMat(1:ngllx-2,1:ngllz-2) )
        Tdomain%specel(n)%MassMat =  LocMassMat
        LocMassMat = Tdomain%specel(n)%Ivx(1:ngllx-2,1:ngllz-2)
        deallocate (Tdomain%specel(n)%Ivx)
        allocate (Tdomain%specel(n)%Ivx(1:ngllx-2,1:ngllz-2) )
        Tdomain%specel(n)%Ivx = LocMassMat *  Tdomain%specel(n)%DumpVx (:,:,1)
        LocMassMat = Tdomain%specel(n)%Ivz(1:ngllx-2,1:ngllz-2)
        deallocate (Tdomain%specel(n)%Ivz)
        allocate (Tdomain%specel(n)%Ivz(1:ngllx-2,1:ngllz-2) )
        Tdomain%specel(n)%Ivz = LocMassMat * Tdomain%specel(n)%DumpVz (:,:,1)
    endif

    ! Dans boucle des faces
    if (Tdomain%sFace(nf)%FPML) then
        Tdomain%sFace(nf)%DumpVx (:,1) =  Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,0)
        Tdomain%sFace(nf)%DumpVx (:,1) = 1./Tdomain%sFace(nf)%DumpVx (:,1)
        Tdomain%sFace(nf)%DumpVx (:,0) =  Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,1)
        Tdomain%sFace(nf)%DumpVx (:,0) = Tdomain%sFace(nf)%DumpVx (:,0) *  Tdomain%sFace(nf)%DumpVx (:,1)

        Tdomain%sFace(nf)%DumpVz (:,1) =  Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,2)
        Tdomain%sFace(nf)%DumpVz (:,1) = 1./Tdomain%sFace(nf)%DumpVz (:,1)
        Tdomain%sFace(nf)%DumpVz (:,0) =  Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,3)
        Tdomain%sFace(nf)%DumpVz (:,0) = Tdomain%sFace(nf)%DumpVz (:,0) *  Tdomain%sFace(nf)%DumpVz (:,1)
        Tdomain%sFace(nf)%Ivx = Tdomain%sFace(nf)%Ivx * Tdomain%sFace(nf)%DumpVx(:,1)
        Tdomain%sFace(nf)%Ivz = Tdomain%sFace(nf)%Ivz * Tdomain%sFace(nf)%DumpVz(:,1)
        deallocate (Tdomain%sFace(nf)%DumpMass)
    endif

    ! Dans boucle des vertexs
    if (tdomain%sVertex(nv_aus)%FPML) then
        Tdomain%sVertex(nv_aus)%DumpVx (1) =  Tdomain%sVertex(nv_aus)%MassMat + Tdomain%sVertex(nv_aus)%DumpMass(0)
        Tdomain%sVertex(nv_aus)%DumpVx (1) = 1./Tdomain%sVertex(nv_aus)%DumpVx (1)
        Tdomain%sVertex(nv_aus)%DumpVx (0) =  Tdomain%sVertex(nv_aus)%MassMat + Tdomain%sVertex(nv_aus)%DumpMass(1)
        Tdomain%sVertex(nv_aus)%DumpVx (0) = Tdomain%sVertex(nv_aus)%DumpVx (0) * Tdomain%sVertex(nv_aus)%DumpVx(1)

        Tdomain%sVertex(nv_aus)%DumpVz (1) =  Tdomain%sVertex(nv_aus)%MassMat + Tdomain%sVertex(nv_aus)%DumpMass(2)
        Tdomain%sVertex(nv_aus)%DumpVz (1) = 1./Tdomain%sVertex(nv_aus)%DumpVz (1)
        Tdomain%sVertex(nv_aus)%DumpVz (0) =  Tdomain%sVertex(nv_aus)%MassMat + Tdomain%sVertex(nv_aus)%DumpMass(3)
        Tdomain%sVertex(nv_aus)%DumpVz (0) = Tdomain%sVertex(nv_aus)%DumpVz (0) * Tdomain%sVertex(nv_aus)%DumpVz(1)
        Tdomain%sVertex(nv_aus)%Ivx(0) = Tdomain%sVertex(nv_aus)%Ivx(0) * Tdomain%sVertex(nv_aus)%DumpVx(1)
        Tdomain%sVertex(nv_aus)%Ivz(0) = Tdomain%sVertex(nv_aus)%Ivz(0) * Tdomain%sVertex(nv_aus)%DumpVz(1)
        deallocate (Tdomain%sVertex(nv_aus)%DumpMass)
    endif


end module FPML
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!



!>
!!\file getiv_e2f.F90
!!\brief Contient la subroutine getIv_Element2Face.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) n_elem
!! \param integer, intent (IN) n_face
!! \param integer, intent (IN) w_face
!! \param logical logic
!<


subroutine getIv_Element2Face (Tdomain, n_elem, n_face, w_face, logic)
    use sdomain
    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n_elem, n_face, w_face
    logical :: logic

    ! local variables
    integer ::  ngll, ngllx, ngllz,i

    ! Modified by Gaetano Festa 01/06/2005

    ngll = Tdomain%sFace(n_face)%ngll
    ngllx = Tdomain%specel(n_elem)%ngllx
    ngllz = Tdomain%specel(n_elem)%ngllz
    if (logic) then
        if (w_face == 0 ) then
            Tdomain%sFace(n_face)%Ivx (1:ngll-2) = Tdomain%sFace(n_face)%Ivx (1:ngll-2) + &
                Tdomain%specel(n_elem)%Ivx (1:ngll-2, 0)
            Tdomain%sFace(n_face)%Ivz (1:ngll-2) = Tdomain%sFace(n_face)%Ivz (1:ngll-2) + &
                Tdomain%specel(n_elem)%Ivz (1:ngll-2, 0)
        else if (w_face == 1 ) then
            Tdomain%sFace(n_face)%Ivx (1:ngll-2) = Tdomain%sFace(n_face)%Ivx (1:ngll-2) + &
                Tdomain%specel(n_elem)%Ivx (ngllx-1,1:ngll-2)
            Tdomain%sFace(n_face)%Ivz (1:ngll-2) = Tdomain%sFace(n_face)%Ivz (1:ngll-2) + &
                Tdomain%specel(n_elem)%Ivz (ngllx-1,1:ngll-2)
        else if (w_face == 2 ) then
            Tdomain%sFace(n_face)%Ivx (1:ngll-2) = Tdomain%sFace(n_face)%Ivx (1:ngll-2) + &
                Tdomain%specel(n_elem)%Ivx (1:ngll-2, ngllz-1)
            Tdomain%sFace(n_face)%Ivz (1:ngll-2) = Tdomain%sFace(n_face)%Ivz (1:ngll-2) + &
                Tdomain%specel(n_elem)%Ivz (1:ngll-2, ngllz-1)
        else
            Tdomain%sFace(n_face)%Ivx (1:ngll-2) = Tdomain%sFace(n_face)%Ivx (1:ngll-2) + &
                Tdomain%specel(n_elem)%Ivx (0,1:ngll-2)
            Tdomain%sFace(n_face)%Ivz (1:ngll-2) = Tdomain%sFace(n_face)%Ivz (1:ngll-2) + &
                Tdomain%specel(n_elem)%Ivz (0,1:ngll-2)
        endif
    else
        if (w_face == 0 ) then
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%Ivx (i) = Tdomain%sFace(n_face)%Ivx (i) + &
                    Tdomain%specel(n_elem)%Ivx (ngll-1-i, 0)
                Tdomain%sFace(n_face)%Ivz (i) = Tdomain%sFace(n_face)%Ivz (i) + &
                    Tdomain%specel(n_elem)%Ivz (ngll-1-i, 0)
            enddo
        else if (w_face == 1 ) then
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%Ivx (i) = Tdomain%sFace(n_face)%Ivx (i) + &
                    Tdomain%specel(n_elem)%Ivx (ngllx-1,ngll-1-i)
                Tdomain%sFace(n_face)%Ivz (i) = Tdomain%sFace(n_face)%Ivz (i) + &
                    Tdomain%specel(n_elem)%Ivz (ngllx-1,ngll-1-i)
            enddo
        else if (w_face == 2 ) then
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%Ivx (i) = Tdomain%sFace(n_face)%Ivx (i) + &
                    Tdomain%specel(n_elem)%Ivx (ngll-1-i, ngllz-1)
                Tdomain%sFace(n_face)%Ivz (i) = Tdomain%sFace(n_face)%Ivz (i) + &
                    Tdomain%specel(n_elem)%Ivz (ngll-1-i, ngllz-1)
            enddo
        else
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%Ivx (i) = Tdomain%sFace(n_face)%Ivx (i) + &
                    Tdomain%specel(n_elem)%Ivx (0,ngll-1-i)
                Tdomain%sFace(n_face)%Ivz (i) = Tdomain%sFace(n_face)%Ivz (i) + &
                    Tdomain%specel(n_elem)%Ivz (0,ngll-1-i)
            enddo
        endif
    endif
    return
end subroutine  getIv_Element2Face

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
