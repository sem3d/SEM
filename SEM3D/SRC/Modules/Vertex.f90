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

    type :: vertex_pml
       real, dimension(0:2) :: Forces1, Forces2, Forces3
       real, dimension(:), allocatable :: DumpMass
       real, dimension(:), allocatable :: Veloc1, Veloc2, Veloc3
       real, dimension(:), allocatable :: DumpVx, DumpVy, DumpVz
       real, dimension (:), allocatable :: Iveloc1, Iveloc2, Iveloc3
       real, dimension (:), allocatable :: Ivx, Ivy, Ivz
       real :: ForcesFl1, ForcesFl2, ForcesFl3, VelPhi1, VelPhi2, VelPhi3
    end type vertex_pml

    type :: vertex
       integer  :: mat_index
       logical :: PML, Abs, FPML
       integer :: Iglobnum_Vertex, global_numbering
       real :: MassMat
       real, dimension(0:2) :: Forces, Displ, Veloc, Accel, V0
       logical :: solid
       ! solid-fluid
       real :: ForcesFl, Phi, VelPhi, AccelPhi, VelPhi0

       type(vertex_pml), allocatable :: spml
#ifdef COUPLAGE
       real, dimension (:), allocatable :: ForcesMka
       real :: tsurfsem
#endif

    end type vertex

contains

    ! ############################################################
    !>
    !! \brief Predicteur pour les vertex
    !!
    !! \param type (Vertex), intent (INOUT) V
    !<
    subroutine Prediction_Vertex_Veloc (V, dt)
        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent(in) :: dt

        V%Forces = V%Displ
        V%V0 = V%Veloc
        return
    end subroutine Prediction_Vertex_Veloc

    ! ###########################################################
    !>
    !! \brief Integration
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (IN) dt
    !<
    subroutine Correction_Vertex_Veloc (V, dt)
        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) :: dt
        integer :: i

        do i = 0,2
            V%Forces(i) = V%MassMat  * V%Forces(i)
        enddo
        V%Veloc = V%V0 + dt * V%Forces
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
            V%spml%Veloc1(i) = V%spml%DumpVx(0) * V%spml%Veloc1(i) + dt * V%spml%DumpVx(1) * V%spml%Forces1(i)
            V%spml%Veloc2(i) = V%spml%DumpVy(0) * V%spml%Veloc2(i) + dt * V%spml%DumpVy(1) * V%spml%Forces2(i)
            V%spml%Veloc3(i) = V%spml%DumpVz(0) * V%spml%Veloc3(i) + dt * V%spml%DumpVz(1) * V%spml%Forces3(i)
        enddo

        V%Veloc = V%spml%Veloc1 + V%spml%Veloc2 + V%spml%Veloc3
        V%Displ = V%Displ +  dt * V%Veloc

        if (V%Abs) then
            V%Veloc = 0
        endif

        return
    end subroutine Correction_Vertex_PML_Veloc

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    subroutine Prediction_Vertex_VelPhi(V, dt)
        implicit none

        type(Vertex), intent(inout) :: V
        real, intent(in) :: dt

        V%VelPhi0 = V%VelPhi
        V%ForcesFl = V%Phi
        return
    end subroutine Prediction_Vertex_VelPhi

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    subroutine Correction_Vertex_VelPhi(V,dt)
        implicit none

        type(Vertex), intent(inout) :: V
        real, intent(in) :: dt


        V%ForcesFl = V%MassMat * V%ForcesFl
        V%VelPhi = V%VelPhi0 + dt * V%ForcesFl
        V%AccelPhi = (V%VelPhi-V%VelPhi0)/dt
        V%Phi = V%Phi + dt * V%VelPhi
        return
    end subroutine Correction_Vertex_VelPhi

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    subroutine Correction_Vertex_PML_VelPhi(V,dt)
        implicit none

        type(Vertex), intent(inout) :: V
        real, intent(in) :: dt

        V%spml%VelPhi1 = V%spml%DumpVx(0) * V%spml%VelPhi1 + dt * V%spml%DumpVx(1) * V%spml%ForcesFl1
        V%spml%VelPhi2 = V%spml%DumpVy(0) * V%spml%VelPhi2 + dt * V%spml%DumpVy(1) * V%spml%ForcesFl2
        V%spml%VelPhi3 = V%spml%DumpVz(0) * V%spml%VelPhi3 + dt * V%spml%DumpVz(1) * V%spml%ForcesFl3

        V%VelPhi = V%spml%VelPhi1 + V%spml%VelPhi2 + V%spml%VelPhi3

        if (V%Abs) then
            V%VelPhi = 0
        endif

        return
    end subroutine Correction_Vertex_PML_VelPhi

    ! ###########################################################
    subroutine Correction_Vertex_fPML_Veloc (V, dt, fil)

        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) :: dt,fil

        integer :: i
        real:: fil2, aus_v

        fil2 = fil**2

        do i = 0,2
            Aus_V = V%spml%Veloc1(i)
            V%spml%Veloc1(i) = V%spml%DumpVx(0) * V%spml%Veloc1(i) + dt * V%spml%DumpVx(1) * V%spml%Forces1(i) + V%spml%Ivx(0) * V%spml%Iveloc1(i)
            V%spml%Iveloc1(i) = Fil2 * V%spml%Iveloc1(i) + 0.5 * (1-Fil2) *  (Aus_V + V%spml%Veloc1(i))

            Aus_V = V%spml%Veloc2(i)
            V%spml%Veloc2(i) = V%spml%DumpVy(0) * V%spml%Veloc2(i) + dt * V%spml%DumpVy(1) * V%spml%Forces2(i) + V%spml%Ivy(0) * V%spml%Iveloc2(i)
            V%spml%Iveloc2(i) = Fil2 * V%spml%Iveloc2(i) + 0.5 * (1-Fil2) *  (Aus_V + V%spml%Veloc2(i))

            Aus_V = V%spml%Veloc3(i)
            V%spml%Veloc3(i) = V%spml%DumpVz(0) * V%spml%Veloc3(i) + dt * V%spml%DumpVz(1) * V%spml%Forces3(i) + V%spml%Ivz(0) * V%spml%Iveloc3(i)
            V%spml%Iveloc3(i) = Fil2 * V%spml%Iveloc3(i) + 0.5 * (1-Fil2) *  (Aus_V + V%spml%Veloc3(i))
        enddo

        V%Veloc = V%spml%Veloc1 + V%spml%Veloc2 + V%spml%Veloc3

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
                Vfree(0:2) = Vfree(0:2) - (V%spml%DumpVz(0) * V%spml%Veloc3(0:2) + dt * V%spml%DumpVz(1) * V%spml%Forces3(0:2) )
            else
                Vfree(0:2) =  V%spml%DumpVz(0) * V%spml%Veloc3(0:2) + dt * V%spml%DumpVz(1) * V%spml%Forces3(0:2)
            endif

        endif

        return
    end subroutine get_vel_vertex
    ! ###########################################################

    subroutine init_vertex(ve)
        type(Vertex), intent(inout) :: ve

        ve%PML = .false.
        ve%Abs = .false.
        ve%FPML = .false.
        ve%solid = .true.
        ve%global_numbering = -1
        ve%Iglobnum_Vertex = -1
        ve%MassMat = 0
        ve%Forces = 0.
        ve%Displ = 0.
        ve%Veloc = 0.
        ve%Accel = 0.
        ve%V0 = 0.
        ve%ForcesFl = 0.
        ve%Phi = 0.
        ve%VelPhi = 0.
        ve%AccelPhi = 0.
        ve%VelPhi0 = 0.
    end subroutine init_vertex

end module svertices
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
