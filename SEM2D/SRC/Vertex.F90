!>
!!\file Vertex.F90
!!\brief Assure la gestion des Vertex.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module svertices

    use smat_solver
    ! Modified by Gaetano 01/06/05

    type :: vertex

       logical :: PML, Abs, reflex, FPML, CPML, ADEPML, is_computed
       integer :: Glob_Numbering, mat_index, Type_DG
       real :: MassMat
       real, dimension (:), allocatable :: DumpMass, DumpVx, DumpVz, Forces1, Forces2, Veloc1, Veloc2
       real, dimension (:), allocatable :: Displ, Veloc, Forces, Accel, V0
       real, dimension (:), allocatable :: Ivx, Ivz,Iveloc1, Iveloc2
       real, dimension (:), allocatable :: Double_Value

#ifdef MKA3D
       real, dimension (:), allocatable :: ForcesMka
#endif
       ! DG
       integer :: valence
       real,   dimension (:), allocatable :: Vect_RK
       integer,dimension (:), allocatable :: Near_Face
       real, dimension (:,:), allocatable :: Kmat
       real,   dimension (:), allocatable :: SmbrLambda
       real,   dimension (:), allocatable :: Lambda, K_up

    end type vertex

contains

    ! ############################################################
    !>
    !! \brief
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) alpha
    !<


    !  subroutine Prediction_Vertex_Veloc (V,alpha,bega, dt)
    subroutine Prediction_Vertex_Veloc (V)
        implicit none

        type (Vertex), intent (INOUT) :: V
        !real, intent (IN) :: bega, dt, alpha

        !  test mariotti
        !     V%Forces = V%Displ + dt * V%Veloc + dt**2 * (0.5 - bega) * V%Accel
        !     V%V0 = V%Veloc
        !     V%Forces = alpha * V%Forces + (1-alpha) * V%Displ
        V%Forces = V%Displ
        V%V0 = V%Veloc

        return
    end subroutine Prediction_Vertex_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !<

    subroutine Correction_Vertex_Veloc (V, dt)
        implicit none

        type (Vertex), intent (INOUT) :: V
        !real, intent (IN) :: bega, gam1
        real, intent (IN) :: dt

        integer :: i

        do i = 0,1
            V%Forces(i) = V%MassMat * V%Forces(i)
        enddo

        if (V%Abs .or. V%reflex) V%Forces = 0
        !  test mariotti
        !      V%Veloc  = V%v0+ dt * V%Forces
        !      V%Accel  = V%Accel + gam1 /dt * (V%Veloc-V%V0)
        !      V%Displ  =  V%Displ + bega * dt * (V%Veloc+V%V0)
        V%Veloc  = V%v0+ dt * V%Forces
        V%Accel  =  (V%Veloc-V%V0)/dt
        V%Displ  =  V%Displ + dt * V%Veloc
        return
    end subroutine Correction_Vertex_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (IN) dt
    !<


    subroutine Correction_Vertex_PML_Veloc (V, dt)
        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) ::  dt

        integer :: i

        if  (V%Abs .or. V%reflex) then
            V%Veloc1 = 0; V%Veloc2 = 0; V%Veloc = 0
        else
            V%V0 = V%Veloc
            do i = 0,1
                V%Veloc1(i) = V%DumpVx(0) * V%Veloc1(i) + dt * V%DumpVx(1)*V%Forces1(i)
                V%Veloc2(i) =V%DumpVz(0) * V%Veloc2(i) + dt * V%DumpVz(1)*V%Forces2(i)
            enddo
            V%Veloc = V%Veloc1 + V%Veloc2
        endif

        return
    end subroutine Correction_Vertex_PML_Veloc

    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (IN) bega
    !! \param real, intent (IN) gam1
    !! \param real, intent (IN) dt
    !<

    subroutine Correction_Vertex_CPML_Veloc (V, dt)
        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) :: dt

        integer :: i

        V%V0 = V%Veloc
        if (V%Abs .or. V%reflex) then
            V%Veloc = 0
        else
            do i = 0,1
                V%Forces(i) = V%MassMat * V%Forces(i)
            enddo
            V%Veloc  = V%v0+ dt * V%Forces
            V%Accel  =  (V%Veloc-V%V0)/dt
            V%Displ  =  V%Displ + dt * V%Veloc
        endif
        return
    end subroutine Correction_Vertex_CPML_Veloc


    ! ###########################################################
    !>
    !! \brief
    !!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (IN) dt
    !! \param real, intent (IN) fil
    !<
    subroutine Correction_Vertex_FPML_Veloc (V, dt, fil)
        implicit none

        type (Vertex), intent (INOUT) :: V
        real, intent (IN) ::  dt, fil

        integer :: i
        real :: fil2
        real :: Ausiliar_Velocity

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

    !>
    !! \brief
    !!
    !! \param type (Vertex), intent (IN) V
    !! \param real, dimension (0:1), intent (INOUT) Vfree
    !! \param logical, intent (IN) logic
    !<


    subroutine get_vfree_vertex(V,Vfree,logic)
        implicit none

        type (Vertex), intent (IN) :: V
        real, dimension (0:1), intent (INOUT) :: Vfree
        logical, intent (IN) :: logic

        if (logic) then
            Vfree(0:1) = Vfree(0:1) - V%MassMat * V%Forces(0:1)
        else
            Vfree(0:1) =  V%MassMat * V%Forces(0:1)
        endif
        return

    end subroutine get_vfree_vertex

    ! ############################################################

    !>
    !! \brief This subroutine computes the kinetic energy on the vertex
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (INOUT) E_kin
    !<

    subroutine  compute_Kinetic_Energy_V (V, Dt, E_kin)
        implicit none

        type (Vertex), intent (IN) :: V
        real, intent (IN)    :: Dt
        real, intent (INOUT) :: E_kin
        real, dimension (0:1)  :: Vel_half

        Vel_half(:) = V%Veloc(:) + 0.5 * dt * V%Forces(:)
        E_kin       = 0.5 /V%MassMat * ( Vel_half(0)*Vel_half(0) &
                                      +Vel_half(1)*Vel_half(1))

    end subroutine compute_Kinetic_Energy_V


    ! ############################################################
    !>
    !! \brief This subroutine performs the cholesky factorisation of the
    !! matrix K defined at the vertex. So it is not really an inversion!
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (INOUT) E_kin
    !<
    subroutine  invert_K_vertex(V)
        implicit none

        type (Vertex), intent (INOUT) :: V
        integer                       :: n, INFO, i, imin, imax
        n = 2 * V%Valence

        ! On range dans K_up les colonnes successives de V%Kmat
        imin = 1
        do i=0,n-1
            imax = imin + i
            V%K_up(imin:imax) = V%Kmat(i,0:i)
            imin = imax +1
        enddo

        ! Appel de la routine Lapack pour factorisation Cholesky
        call SPPTRF( 'U', n, V%K_up, INFO)

        if (INFO .NE. 0) STOP "Factorisation of vertex matrix not successfull"

    end subroutine invert_K_vertex


    ! ############################################################
    !>
    !! \brief This subroutine solves a linear system K * Lambda = smbrLambda
    !! using a previously done cholesky factorisation K into U**T*U.
    !! The upper triangular matrix U is stored colomnwise into the vector K_up.
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (INOUT) E_kin
    !<
    subroutine  solve_lambda_vertex(V)
        implicit none

        type (Vertex), intent (INOUT)  :: V
        real, dimension(1:2*V%Valence) :: rhs
        integer                        :: n, INFO
        n = 2 * V%Valence

        ! Changing Right Hand Side numbering to use Lapack routine
        rhs(1:n) = V%smbrLambda(0:n-1)
        call SPPTRS( 'U', n, 1, V%K_up, rhs, n, INFO )

        ! If any error :
        if (INFO .NE. 0) STOP "Inversion of vertex matrix not successfull"

        ! Assigning the solution to V%Lambda
        V%Lambda(0:n-1) = rhs(1:n)
        V%SmbrLambda(:) = 0.

    end subroutine solve_lambda_vertex


    ! ############################################################

end module svertices
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
