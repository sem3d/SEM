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
    use constants
    ! Modified by Gaetano 01/06/05

    type :: vertex

       logical :: PML, Abs, reflex, CPML, ADEPML, is_computed
       integer :: Glob_Numbering, mat_index, Type_DG
       real :: MassMat
       real, dimension (:), allocatable :: DumpMass, DumpVx, DumpVz, Forces1, Forces2, Veloc1, Veloc2
       real, dimension (:), allocatable :: Displ, Veloc, Forces, Accel, V0
       real, dimension (:), allocatable :: Double_Value

#ifdef MKA3D
       real, dimension (:), allocatable :: ForcesMka
#endif
       ! DG
       integer :: valence
       real,   dimension (:), allocatable :: Vect_RK
       integer,dimension (:), allocatable :: Near_Face
       real, dimension (:,:), allocatable :: Kmat, Kmat_05dt
       real,   dimension (:), allocatable :: SmbrLambda
       real,   dimension (:), allocatable :: Lambda, K_up, K_up_05dt

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


    ! ############################################################
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
    subroutine  invert_K_vertex(V,nv)
        implicit none

        type (Vertex), intent (INOUT) :: V
        integer,       intent (IN)    :: nv
        integer                       :: n, INFO, i, imin, imax
        n = 2 * V%Valence

        ! On range dans K_up les colonnes successives de V%Kmat
        imin = 1
        do i=0,n-1
            imax = imin + i
            V%K_up(imin:imax) = V%Kmat(0:i,i)
            V%K_up_05dt(imin:imax) = V%Kmat_05dt(0:i,i)
            imin = imax +1
        enddo

        ! Appel de la routine Lapack pour factorisation Cholesky
        !call SPPTRF( 'U', n, V%K_up, INFO)
        !if (INFO .NE. 0) then
        !    WRITE(*,*) "Factorisation vertex matrix not successfull for vertex  ",nv," rank mat ",n
        !    STOP "Factorisation of vertex matrix not successfull. End of Program"
        !endif
        call SPPTRF( 'U', n, V%K_up_05dt, INFO)
        if (INFO .NE. 0) then
            WRITE(*,*) "Factorisation vertex matrix not successfull for vertex  ",nv," rank mat ",n
            STOP "Factorisation of vertex matrix not successfull. End of Program"
        endif

    end subroutine invert_K_vertex


    ! ############################################################
    !>
    !! \brief This subroutine solves a linear system K * Lambda = smbrLambda
    !! using a previously done cholesky factorisation K into U**T*U.
    !! The upper triangular matrix U is stored colomnwise into the vector K_up.
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (INOUT) E_kin
    !<
    subroutine  solve_lambda_vertex(V, nv, demi_dt)
        implicit none

        type (Vertex), intent (INOUT)  :: V
        integer,       intent (IN)     :: nv
        logical,       intent (IN)     :: demi_dt
        real, dimension(1:2*V%Valence) :: rhs
        integer                        :: n, INFO
        n = 2 * V%Valence

        ! Changing Right Hand Side numbering to use Lapack routine
        V%Lambda = 0.
        rhs(1:n) = V%smbrLambda(0:n-1)

        if (demi_dt .EQV. HALF_DT) then
            call SPPTRS( 'U', n, 1, V%K_up_05dt, rhs, n, INFO )
            V%Lambda(0:n-1) = rhs(1:n)
            call check_system_inversion(V, V%Kmat_05dt, nv)
        else
            call SPPTRS( 'U', n, 1, V%K_up, rhs, n, INFO )
            V%Lambda(0:n-1) = rhs(1:n)
            call check_system_inversion(V, V%Kmat, nv)
        endif
        ! If any error :
        if (INFO .NE. 0) STOP "Inversion of vertex matrix not successfull"

        V%SmbrLambda(:) = 0.

        ! Setting velocities to zero in case of reflective BC :
        if (V%reflex) V%Lambda(:) = 0.

    end subroutine solve_lambda_vertex


    ! ############################################################
    !>
    !! \brief This subroutine checks the inversion of linear system :
    !! K * Lambda = smbrLambda
    !! \param type (Vertex), intent (INOUT) V
    !! \param real, intent (INOUT) E_kin
    !<
    subroutine  check_system_inversion (V, Mat, nv)
        implicit none

        type (Vertex), intent (INOUT)  :: V
        integer,       intent (IN)     :: nv
        real, dimension(0:2*V%Valence-1,0:2*V%Valence-1), intent (IN) :: Mat
        integer                        :: n, i, j
        real                           :: tol, res, ratio

        n = 2 * V%Valence
        tol = 1.E-10

        do i=0,n-1
            res = 0.
            do j=0,n-1
                res = res + Mat(i,j)*V%Lambda(j)
            enddo
            if (abs(V%SmbrLambda(i)) .LT. 1.E-24) then
                if (abs(res) .GT. 1.E-24) &
                    write(*,*) "Problem-type 0 in solving system on vertex ",nv," on line :", i
            elseif (abs(res) .LT. 1.E-24) then
                if (abs(V%SmbrLambda(i)) .GT. 1.E-24) &
                    write(*,*) "Problem-type 1 in solving system on vertex ",nv," on line :", i
            else
                ratio = abs((res - V%SmbrLambda(i)) / V%SmbrLambda(i))
                if (ratio .GT. tol) then
                    write(*,*) "Problem in solving system on vertex ",nv," on line :", i," error :",ratio
                    !STOP "End of computation"
                endif
            endif
        enddo

    end subroutine check_system_inversion

    ! ############################################################


end module svertices
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
