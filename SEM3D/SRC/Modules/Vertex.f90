module svertices


    type :: vertex
       logical :: PML, Abs, FPML
       integer :: Iglobnum_Vertex, global_numbering,mat_index
       real :: MassMat
       real, dimension(0:2) :: Forces, Displ, Veloc, Accel, V0
       real, dimension(0:2) :: Forces1, Forces2, Forces3
       real, dimension(:), pointer :: DumpMass
       real, dimension(:), pointer :: Veloc1, Veloc2, Veloc3
       real, dimension(:), pointer :: DumpVx, DumpVy, DumpVz
       real, dimension(:), pointer :: Iveloc1, Iveloc2,Iveloc3
       real, dimension(:), pointer :: Ivx, Ivy, Ivz
       logical :: solid
       ! solid-fluid
       real :: ForcesFl, Phi, VelPhi, AccelPhi, VelPhi0
       real :: ForcesFl1, ForcesFl2, ForcesFl3, VelPhi1, VelPhi2, VelPhi3

    end type vertex

contains

    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    subroutine Prediction_Vertex_Veloc (V, alpha, bega, gam1, dt)

        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) :: alpha, bega, gam1, dt


        V%Forces = V%Displ + dt * V%Veloc + dt**2 * (0.5 - bega) * V%Accel
        V%V0 = V%Veloc
        V%Forces = alpha * V%Forces + (1-alpha) * V%Displ

        return
    end subroutine Prediction_Vertex_Veloc
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    subroutine Prediction_Vertex_VelPhi(V,alpha,bega,gam1,dt)

        implicit none

        type(Vertex), intent(inout) :: V
        real, intent(in) :: alpha, bega, gam1, dt


        V%ForcesFl = V%Phi + dt * V%VelPhi + dt**2 * (0.5d0 - bega) * V%AccelPhi
        V%VelPhi0 = V%VelPhi
        V%ForcesFl = alpha * V%ForcesFl + (1d0-alpha) * V%Phi

        return
    end subroutine Prediction_Vertex_VelPhi
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    subroutine Correction_Vertex_Veloc (V, bega, gam1, dt)

        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) :: bega, gam1, dt

        integer :: i


        do i = 0,2
            V%Forces(i) = V%MassMat * V%Forces(i)
        enddo
        V%Veloc = V%v0 + dt * V%Forces
        V%Accel = V%Accel + gam1 / dt * (V%Veloc-V%V0)
        V%Displ = V%Displ + bega * dt * (V%Veloc+V%V0)

        return
    end subroutine Correction_Vertex_Veloc
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    subroutine Correction_Vertex_VelPhi(V,bega,gam1,dt)
        implicit none

        type(Vertex), intent(inout) :: V
        real, intent(in) :: bega, gam1, dt


        V%ForcesFl = V%MassMat * V%ForcesFl

        V%VelPhi = V%VelPhi0 + dt * V%ForcesFl
        V%AccelPhi = V%AccelPhi + gam1 / dt * (V%VelPhi-V%VelPhi0)
        V%Phi = V%Phi + bega * dt * (V%VelPhi+V%VelPhi0)

        return
    end subroutine Correction_Vertex_VelPhi
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
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
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    subroutine Correction_Vertex_PML_VelPhi(V,dt)
        implicit none

        type(Vertex), intent(inout) :: V
        real, intent(in) :: dt

        V%VelPhi1 = V%DumpVx(0) * V%VelPhi1 + dt * V%DumpVx(1) * V%ForcesFl1
        V%VelPhi2 = V%DumpVy(0) * V%VelPhi2 + dt * V%DumpVy(1) * V%ForcesFl2
        V%VelPhi3 = V%DumpVz(0) * V%VelPhi3 + dt * V%DumpVz(1) * V%ForcesFl3

        V%VelPhi = V%VelPhi1 + V%VelPhi2 + V%VelPhi3

        if (V%Abs) then
            V%VelPhi = 0
        endif

        return
    end subroutine Correction_Vertex_PML_VelPhi
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
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
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
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
