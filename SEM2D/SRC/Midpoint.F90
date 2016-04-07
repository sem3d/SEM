!>
!!\file Newmark.F90
!!\brief Algorithmes de type midpoint avec une approche Galerkin Discontinu (HDG).
!!\version 1.0
!!\date 10/03/2009
!! La routine Midpoint assure l'intégration en temps via un des algorithmes Midpoint
!! qui peuvent etre iteratifs, ou en une seule passe. Les contraintes aux surfaces peuvent
!! etre prises en compte de facon implicites ou explicites, ce qui donne les differents
!! algorithmes et routines ci-dessous.
!<

module smidpoint
    use sdomain
    use scouplage
    use constants
    use mpi
    implicit none
contains



!>
!!\brief This Subroutine performs an iterative Implicit Midpoint time integration
!! The second member is computed at tn+1/2 using an iterative EXPLICIT approach.
!!\version 1.0
!!\date 20/11/2014
!! This subroutine is used only HDG elements
!! \param type (Domain), intent (INOUT) Tdomain
!! \param real         , intent (IN)    Dt
!<
subroutine Midpoint_impl_expl (Tdomain,Dt,n_it_max)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: dt
    integer, intent(in)   :: n_it_max

    ! local variables
    integer :: n, iter
    real :: timelocal

    ! Initialization Phase
    do n=0,Tdomain%n_elem-1
        ! Initialization of variables at tn :
        Tdomain%specel(n)%Strain0(:,:,:) = Tdomain%specel(n)%Strain(:,:,:)
        Tdomain%specel(n)%V0(:,:,:)      = Tdomain%specel(n)%Veloc (:,:,:)
        if (Tdomain%specel(n)%PML) call initialize_Psi(Tdomain%specel(n))
    enddo

    ! Initial Predictor phase.
    timelocal = Tdomain%TimeD%rtime
    call Forward_Euler_Resolution (Tdomain,timelocal,0.5*Dt,NOT_COMPUTE_VHAT)

    iter= 0
    timelocal = Tdomain%TimeD%rtime + 0.5*Dt
    do while (iter<n_it_max)
        ! Explicit resolution phase
        call Forward_Euler_Resolution (Tdomain,timelocal,0.5*Dt,COMPUTE_VHAT)
        iter = iter+1
    enddo

    ! Final Midpoint Evaluation using the values at tn+1/2 converged
    call Forward_Euler_Resolution(Tdomain,timelocal,Dt,COMPUTE_VHAT)

    ! Obtaining Vhat on faces at tn+1
    call get_vhat_from_current_state(Tdomain)

    return
end subroutine Midpoint_impl_expl


!>
!!\brief This Subroutine performs an iterative Implicit Midpoint time integration
!! The second member is computed at tn+1/2 using an iterative SEMI-IMPLICIT approach.
!!\version 1.0
!!\date 20/11/2014
!! This subroutine is used only HDG elements
!! \param type (Domain), intent (INOUT) Tdomain
!! \param real         , intent (IN)    Dt
!<
subroutine Midpoint_impl_semi_impl (Tdomain,Dt,n_it_max)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: dt
    integer, intent(in)   :: n_it_max

    ! local variables
    integer :: n, iter
    real    :: timelocal

    ! Initialization Phase
    do n=0,Tdomain%n_elem-1
        ! Initialization of variables at tn :
        Tdomain%specel(n)%Strain0(:,:,:) = Tdomain%specel(n)%Strain(:,:,:)
        Tdomain%specel(n)%V0(:,:,:)      = Tdomain%specel(n)%Veloc (:,:,:)
        if (Tdomain%specel(n)%PML) call initialize_Psi(Tdomain%specel(n))
    enddo
    do n=0,Tdomain%n_face-1
        ! Initialization of traces at tn :
        Tdomain%sFace(n)%V0(:,:) = Tdomain%sFace(n)%Veloc(:,:)
    enddo

    ! Initial Predictor phase.
    timelocal = Tdomain%TimeD%rtime
    call Forward_Euler_Resolution(Tdomain,timelocal,0.5*Dt,NOT_COMPUTE_VHAT)

    iter = 0
    timelocal = Tdomain%TimeD%rtime + 0.5*Dt
    do while (iter<n_it_max)
        ! Semi-Implicit resolution phase
        call Semi_Implicit_Resolution (Tdomain,timelocal,0.5*Dt)
        iter = iter+1
    enddo

    ! Final Midpoint Evaluation using the values at tn+1/2 converged
    call Forward_Euler_Resolution(Tdomain,timelocal,Dt,COMPUTE_VHAT)

    ! Obtaining Vhat on faces at tn+1
    call get_vhat_from_current_state(Tdomain)

    return
end subroutine Midpoint_impl_semi_impl


!>
!!\brief This Subroutine performs Predictor-Multicorrector (PMC) algorithm for
!! time integration using an SPLITTED semi-implicit approach.
!! Ideed, in this case, the traces terms are choosen at tn+1 while volumic fields
!! in RHS are choosen at tn+1/2. This method has been prooved to be not
!! so good (february 2015) precisely due to that splitting.
!!\version 1.0
!!\date 20/11/2014
!! This subroutine is used only HDG elements
!! \param type (Domain), intent (INOUT) Tdomain
!! \param real         , intent (IN)    Dt
!! \param integer      , intent (IN)    n_it_max
!<
subroutine Midpoint_Test(Tdomain,Dt,n_it_max)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: dt
    integer, intent(in)   :: n_it_max

    ! local variables
    integer :: n, iter
    real :: timelocal


    ! Initialization Phase
    do n=0,Tdomain%n_elem-1
        ! Computation of the  prediction :
        Tdomain%specel(n)%Strain0(:,:,:) = Tdomain%specel(n)%Strain(:,:,:)
        Tdomain%specel(n)%V0(:,:,:)      = Tdomain%specel(n)%Veloc (:,:,:)
        if (Tdomain%specel(n)%PML) call initialize_Psi(Tdomain%specel(n))
    enddo
    do n=0,Tdomain%n_face-1
        Tdomain%sFace(n)%V0 = Tdomain%sFace(n)%Veloc
    enddo

    ! Midpoint method :
    iter= 0
    timelocal = Tdomain%TimeD%rtime
    call Forward_Euler_Resolution(Tdomain,timelocal,0.5*Dt,NOT_COMPUTE_VHAT)

    timelocal = Tdomain%TimeD%rtime + 0.5*Dt

    do while (iter<n_it_max)

        ! Prediction Phase :
        if (iter > 0) then
            do n=0,Tdomain%n_elem-1
                Tdomain%specel(n)%Strain = 0.5 * (Tdomain%specel(n)%Strain0 + Tdomain%specel(n)%Strain)
                Tdomain%specel(n)%Veloc  = 0.5 * (Tdomain%specel(n)%V0      + Tdomain%specel(n)%Veloc )
                !if (Tdomain%specel(n)%ADEPML) call Prediction_Psi(Tdomain%specel(n))
            enddo
        endif

        ! Semi-Iplicit Resolution Phase
        call Semi_Implicit_Resolution_tnplus1 (Tdomain,timelocal,Dt)

        iter = iter+1
    enddo
    !call get_vhat_from_current_state(Tdomain)

    return
end subroutine Midpoint_Test



!>
!!\brief This Subroutine performs a midpoint method for advancing in time
!! using an explicit approach, contrary to the previous one.
!!\version 1.0
!!\date 20/11/2014
!! This subroutine is used only HDG elements
!! \param type (Domain), intent (INOUT) Tdomain
!! \param real         , intent (IN)    Dt
!<
subroutine Semi_Implicit_Resolution (Tdomain,timelocal,Dt)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: timelocal,dt

    ! local variables
    integer :: n, mat

    ! Building second members (= forces) of systems.
    do n=0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        call compute_InternalForces_HDG  (Tdomain%specel(n), &
            Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimex, &
            Tdomain%sSubDomain(mat)%hprimez,Tdomain%sSubDomain(mat)%hTprimez)
        !call add_tau_v (Tdomain%specel(n))
        call add_previous_state2forces (Tdomain%specel(n), Dt)
        if (Tdomain%specel(n)%PML) call update_Psi_ADEPML(Tdomain%specel(n), &
            Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimez, Dt)
    enddo

    ! External Forces computation
    call Compute_External_Forces(Tdomain,timelocal)

    ! Compute second member "R" of the Lagrange multiplicator system
    do n=0,Tdomain%n_elem-1
        ! Compute current element contribution to R
        call compute_smbr_R(Tdomain%specel(n),Dt)
        ! Assembles R on the faces and vertices
        call get_R_el2fv(Tdomain,n)
    enddo

    ! Solve linear systems on the vertices
    do n=0,Tdomain%n_vertex-1
        call solve_lambda_vertex(Tdomain%sVertex(n),n)
    enddo

    ! Constructing the Lambda (= velocities vhat) on the faces
    do n=0,Tdomain%n_face-1
        ! Computes lambda (= Vhat) on Face's inner nodes
        call compute_Vhat_face (Tdomain%sFace(n))
        ! Get lambda from near vertices
        call Get_lambda_v2f (Tdomain, n)
    enddo

!!!!!!!!!!!!!  MPI COMMUNICATIONS HERE  !!!!!!!!!!!!!!!!!!!

    ! Local Solvers at element level
    do n=0,Tdomain%n_elem-1
        ! Communication Lambda from faces to elements
        call get_Vhat_f2el(Tdomain,n)
        if(Tdomain%type_bc==DG_BC_REFL .OR. Tdomain%specel(n)%PML) call enforce_diriclet_BC(Tdomain,n)
        ! Local Solver
        call local_solver(Tdomain%specel(n),Dt)
    enddo

    return
end subroutine Semi_Implicit_Resolution


!>
!!\brief This Subroutine performs an explicit forward euler step fot time integration.
!! Therefore all the Right Hand Side is computed explicitly for advancing variables
!! with a Dt timestep.
!!\version 1.0
!!\date 20/11/2014
!! This subroutine is used only HDG elements
!! \param type (Domain), intent (INOUT) Tdomain
!! \param real         , intent (IN)    timelocal,Dt
!! \param logical      , intent (IN)    Vhat_computed
!<
subroutine Forward_Euler_Resolution (Tdomain,timelocal,Dt,computeVhat)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: timelocal,dt
    logical, intent(in)   :: computeVhat

    ! local variables
    integer :: n, mat, nface

    ! Building second members (= forces) of systems.
    do n=0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        call compute_InternalForces_HDG_Weak(Tdomain%specel(n), &
            Tdomain%sSubDomain(mat)%hprimex, &
            Tdomain%sSubDomain(mat)%hTprimez)
        if (computeVhat) &
            call compute_TracFace (Tdomain%specel(n))
        if (Tdomain%specel(n)%PML) call update_Psi_ADEPML(Tdomain%specel(n), &
            Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimez, Dt)
    enddo

    ! External Forces computation
    call Compute_External_Forces(Tdomain,timelocal)

    ! Calcul de la trace des vitesses Vhat
    if (computeVhat) then
        ! Envoi tractions sur faces
        do n = 0, Tdomain%n_elem-1
            call get_traction_el2f(Tdomain,n)
        enddo
        ! Calcul des Vhat sur les faces
        do nface = 0, Tdomain%n_face-1
            call Compute_Vhat_Face_Expl(Tdomain%sFace(nface))
        enddo
    endif

    ! Constructing the Lambda (= velocities vhat) on the faces
    do n=0,Tdomain%n_elem-1
        call get_Vhat_f2el(Tdomain,n)
        call Compute_Traces (Tdomain%specel(n),.true.)
        if(Tdomain%type_bc==DG_BC_REFL .OR. Tdomain%specel(n)%PML) call enforce_diriclet_BC(Tdomain,n)
        ! Add previous state to forces :
        call add_previous_state2forces (Tdomain%specel(n),Dt)
        Tdomain%specel(n)%Forces(:,:,2) = Tdomain%specel(n)%Forces(:,:,2) &
            - 1./Dt * Tdomain%specel(n)%Acoeff(:,:,12) * Tdomain%specel(n)%Strain0(:,:,2)
        if(Tdomain%type_bc==DG_BC_REFL .OR. Tdomain%specel(n)%PML) call enforce_diriclet_BC(Tdomain,n)
        call inversion_massmat(Tdomain%specel(n))
        Tdomain%specel(n)%Strain(:,:,:) = Dt * Tdomain%specel(n)%Forces(:,:,0:2)
        Tdomain%specel(n)%Veloc (:,:,:) = Dt * Tdomain%specel(n)%Forces(:,:,3:4)
    enddo

    return
end subroutine Forward_Euler_Resolution


!>
!!\brief This Subroutine performs a midpoint method for advancing in time
!! using an explicit approach, contrary to the previous one.
!!\version 1.0
!!\date 20/11/2014
!! This subroutine is used only HDG elements
!! \param type (Domain), intent (INOUT) Tdomain
!! \param real         , intent (IN)    Dt
!<
subroutine Semi_Implicit_Resolution_tnplus1 (Tdomain,timelocal,Dt)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: timelocal,dt

    ! local variables
    integer :: n, mat

    ! Building second members (= forces) of systems.
    do n=0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        call compute_InternalForces_HDG  (Tdomain%specel(n), &
            Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimex, &
            Tdomain%sSubDomain(mat)%hprimez,Tdomain%sSubDomain(mat)%hTprimez)
        call add_tau_v0 (Tdomain%specel(n))
        call add_previous_state2forces (Tdomain%specel(n), Dt)
        if (Tdomain%specel(n)%PML) call update_Psi_ADEPML(Tdomain%specel(n), &
            Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimez, Dt)
        call Get_V0_f2el (Tdomain, n) !!!! <-- MODIF !!!
        Tdomain%specel(n)%TracFace(:,:) = 0.       !!!! <-- MODIF !!!
        call compute_Traces (Tdomain%specel(n), .false.) !!!! <-- MODIF !!!
    enddo

    ! External Forces computation
    call Compute_External_Forces(Tdomain,timelocal)

    ! Compute second member "R" of the Lagrange multiplicator system
    do n=0,Tdomain%n_elem-1
        ! Compute current element contribution to R
        call compute_smbr_R(Tdomain%specel(n),Dt)
        ! Assembles R on the faces and vertices
        call get_R_el2fv(Tdomain,n)
    enddo

    ! Solve linear systems on the vertices
    do n=0,Tdomain%n_vertex-1
        call solve_lambda_vertex(Tdomain%sVertex(n),n)
    enddo

    ! Constructing the Lambda (= velocities vhat) on the faces
    do n=0,Tdomain%n_face-1
        ! Computes lambda (= Vhat) on Face's inner nodes
        call compute_Vhat_face (Tdomain%sFace(n))
        ! Get lambda from near vertices
        call Get_lambda_v2f (Tdomain, n)
    enddo

!!!!!!!!!!!!!  MPI COMMUNICATIONS HERE  !!!!!!!!!!!!!!!!!!!

    ! Local Solvers at element level
    do n=0,Tdomain%n_elem-1
        ! Communication Lambda from faces to elements
        call get_Vhat_f2el(Tdomain,n)
        Tdomain%specel(n)%Vhat = 0.5 * Tdomain%specel(n)%Vhat
        if(Tdomain%type_bc==DG_BC_REFL .OR. Tdomain%specel(n)%PML) call enforce_diriclet_BC(Tdomain,n)
        ! Local Solver
        call local_solver(Tdomain%specel(n),Dt)
    enddo

    return
end subroutine Semi_Implicit_Resolution_tnplus1


!>
!!\brief This Subroutine performs Vhat, the computation of the traces
!! velocities on the faces, from the current state of each elements.
!!\version 1.0
!!\date 10/05/2015
!! This subroutine is used only HDG elements
!! \param type (Domain), intent (INOUT) Tdomain
!<
subroutine get_vhat_from_current_state(Tdomain)

    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables
    integer :: n, nface

    do n=0,Tdomain%n_elem-1
        call compute_TracFace (Tdomain%specel(n))
        call get_traction_el2f(Tdomain,n)
    enddo
    ! Calcul des Vhat sur les faces
    do nface = 0, Tdomain%n_face-1
        call Compute_Vhat_Face_Expl(Tdomain%sFace(nface))
    enddo

end subroutine get_vhat_from_current_state

end module smidpoint
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
