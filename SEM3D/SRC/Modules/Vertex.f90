!>
!!\file Vertex.f90
!!\brief Assure la gestion des Vertex.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module svertices

    ! Modified by Gaetano 31/01/2005
    ! Modified by Paul 06/11/2005

    type :: vertex

       logical :: PML, Abs, FPML

       integer :: Iglobnum_Vertex

       real :: MassMat
       real, dimension(0:2) :: Forces, Displ, Veloc, Accel, V0
       real, dimension(0:2) :: Forces1, Forces2, Forces3
       real, dimension(:), pointer :: DumpMass
       real, dimension(:), pointer :: Veloc1, Veloc2, Veloc3
       real, dimension(:), pointer :: DumpVx, DumpVy, DumpVz
       real, dimension (:), pointer :: Iveloc1, Iveloc2,Iveloc3
       real, dimension (:), pointer :: Ivx, Ivy, Ivz

#ifdef MKA3D
       real, dimension (:), pointer :: ForcesMka
       !     integer, dimension (:), pointer :: FlagMka
       real, dimension (:), pointer :: tsurfsem
#endif

    end type vertex

contains

    ! ############################################################
    !>
    !! \fn subroutine Prediction_Vertex_Veloc (V,alpha,bega, gam1,dt)
    !! \brief
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !<


    !subroutine Prediction_Vertex_Veloc (V, alpha, bega, dt)
    subroutine Prediction_Vertex_Veloc (V)
        implicit none

        type (Vertex), intent (INOUT) :: V
        !real, intent (IN) :: alpha, bega, dt


        !  modification schema en temps mariotti pour couplage
        !     V%Forces = V%Displ + dt * V%Veloc + dt**2 * (0.5 - bega) * V%Accel
        !     V%V0 = V%Veloc
        !     V%Forces = alpha * V%Forces + (1-alpha) * V%Displ
        V%Forces = V%Displ
        V%V0 = V%Veloc

        return
    end subroutine Prediction_Vertex_Veloc

    ! ###########################################################
    !>
    !! \fn subroutine Correction_Vertex_Veloc (V, bega, gam1,dt)
    !! \brief
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !<


    !subroutine Correction_Vertex_Veloc (V, bega, gam1, dt)
    subroutine Correction_Vertex_Veloc (V, dt)
        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) :: dt
        integer :: i


        real :: xmas

#ifdef MKA3D
        xmas = 0.
        if (  V%tsurfsem(0) > 0. ) then
            xmas = 1.
            !    test volI2b
            !             xmas = 0.
        endif
#endif

        do i = 0,2
#ifdef MKA3D
            !   test inter1
            !       V%Forces(i) =  V%MassMat  * V%Forces(i)/(1. + V%FlagMka(i) )
            V%Forces(i) =  V%MassMat  * V%Forces(i)/(1. + xmas )
#else
            V%Forces(i) = V%MassMat  * V%Forces(i)
#endif
        enddo

        !  modification schema en temps mariotti pour couplage
        !V%Veloc = V%v0 + dt * V%Forces
        !V%Accel = V%Accel + gam1 / dt * (V%Veloc-V%V0)
        !V%Displ = V%Displ + bega * dt * (V%Veloc+V%V0)
        V%Veloc = V%v0 + dt * V%Forces
        V%Accel =  (V%Veloc-V%V0)/dt
        V%Displ = V%Displ +  dt * V%Veloc

        return
    end subroutine Correction_Vertex_Veloc

    ! ############################################################
    !>
    !! \fn subroutine Correction_Vertex_PML_Veloc (V, dt)
    !! \brief
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (IN) dt
    !<


    subroutine Correction_Vertex_PML_Veloc (V, dt)

        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) :: dt

        integer :: i


        do i = 0,2
            V%Veloc1(i) = V%DumpVx(0) * V%Veloc1(i) + dt * V%DumpVx(1) * V%Forces1(i)
            V%Veloc2(i) = V%DumpVy(0) * V%Veloc2(i) + dt * V%DumpVy(1) * V%Forces2(i)
            V%Veloc3(i) = V%DumpVz(0) * V%Veloc3(i) + dt * V%DumpVz(1) * V%Forces3(i)
        enddo

        V%Veloc = V%Veloc1 + V%Veloc2 + V%Veloc3

        if (V%Abs) then
            V%Veloc = 0
        endif

        return
    end subroutine Correction_Vertex_PML_Veloc

    ! ###########################################################
    subroutine Correction_Vertex_fPML_Veloc (V, dt,fil)

        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) :: dt,fil

        integer :: i
        real:: fil2, aus_v

        fil2 = fil**2


        do i = 0,2
            Aus_V = V%Veloc1(i)
            V%Veloc1(i) = V%DumpVx(0) * V%Veloc1(i) + dt * V%DumpVx(1) * V%Forces1(i) + V%Ivx(0) * V%Iveloc1(i)
            V%Iveloc1(i) = Fil2 * V%Iveloc1(i) + 0.5 * (1-Fil2) *  (Aus_V + V%Veloc1(i))

            Aus_V = V%Veloc2(i)
            V%Veloc2(i) = V%DumpVy(0) * V%Veloc2(i) + dt * V%DumpVy(1) * V%Forces2(i) + V%Ivy(0) * V%Iveloc2(i)
            V%Iveloc2(i) = Fil2 * V%Iveloc2(i) + 0.5 * (1-Fil2) *  (Aus_V + V%Veloc2(i))

            Aus_V = V%Veloc3(i)
            V%Veloc3(i) = V%DumpVz(0) * V%Veloc3(i) + dt * V%DumpVz(1) * V%Forces3(i) + V%Ivz(0) * V%Iveloc3(i)
            V%Iveloc3(i) = Fil2 * V%Iveloc3(i) + 0.5 * (1-Fil2) *  (Aus_V + V%Veloc3(i))
        enddo

        V%Veloc = V%Veloc1 + V%Veloc2 + V%Veloc3

        if (V%Abs) then
            V%Veloc = 0
        endif

        return
    end subroutine Correction_Vertex_FPML_Veloc
    ! ###########################################################

    subroutine get_vel_vertex(V,Vfree,dt,logic)
        implicit none

        type (Vertex), intent (IN) :: V
        real, dimension (0:2), intent (INOUT) :: Vfree
        logical, intent (IN) :: logic
        real, intent(IN) :: dt

        if (.not. V%PML) then

            if (logic) then
                Vfree(0:2) = Vfree(0:2) -  ( V%V0(0:2) + dt*V%MassMat*V%Forces(0:2) )
            else
                Vfree(0:2) =  V%V0(0:2) + dt*V%MassMat*V%Forces(0:2)
            endif

        else

            if (logic) then
                Vfree(0:2) = Vfree(0:2) - (V%DumpVz(0) * V%Veloc3(0:2) + dt * V%DumpVz(1) * V%Forces3(0:2) )
            else
                Vfree(0:2) =  V%DumpVz(0) * V%Veloc3(0:2) + dt * V%DumpVz(1) * V%Forces3(0:2)
            endif

        endif

        return
    end subroutine get_vel_vertex
    ! ###########################################################

end module svertices
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
