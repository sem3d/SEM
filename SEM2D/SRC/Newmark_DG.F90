!>
!!\file Newmark.F90
!!\brief Algorithme de Newmark en cas de presence d'elements Galerkin Discontinu.
!!\version 1.0
!!\date 10/03/2009
!! La routine Newmark assure la résolution des équations via un algorithme de predicteur-multi-correcteur
!! des vitesses avec une formulation contrainte-vitesse décalée en temps dans les PML.
!<

module snewmark_pmc
    use sdomain
    use scouplage
    use mpi
    implicit none
contains

subroutine Newmark_PMC (Tdomain,Dt,n_it_max)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: dt
    integer, intent(in)   :: n_it_max

    ! local variables
    integer :: n, mat, iter
    real :: timelocal

    ! Predictor-MultiCorrector Newmark Veloci

    ! Initialization Phase
    do n=0,Tdomain%n_elem-1
        ! Computation of the  prediction :
        Tdomain%specel(n)%Strain0(:,:,:) = Tdomain%specel(n)%Strain(:,:,:)
        Tdomain%specel(n)%V0(:,:,:)      = Tdomain%specel(n)%Veloc (:,:,:)
        if (Tdomain%specel(n)%PML) call initialize_Psi(Tdomain%specel(n))
    enddo

    ! Midpoint method :
    iter= 0
    call Newmark_PMC_Explicit (Tdomain,Dt,1) ! <---- A SUPPRIMER
    iter = 1                                 ! <---- A SUPPRIMER

    do while (iter<n_it_max)

        ! Local time is tn for the first step, and tn+1/2 for the next ones
        if (iter==0) then
            timelocal = Tdomain%TimeD%rtime
        else
            timelocal = Tdomain%TimeD%rtime + 0.5*Dt
        endif

        ! Prediction Phase :
        do n=0,Tdomain%n_elem-1
            Tdomain%specel(n)%Strain = 0.5 * (Tdomain%specel(n)%Strain0 + Tdomain%specel(n)%Strain)
            Tdomain%specel(n)%Veloc  = 0.5 * (Tdomain%specel(n)%V0      + Tdomain%specel(n)%Veloc )
            !if (Tdomain%specel(n)%ADEPML) call Prediction_Psi(Tdomain%specel(n))
        enddo

        ! Building second members (= forces) of systems.
        do n=0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            call compute_InternalForces_HDG  (Tdomain%specel(n), &
              Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimex, &
              Tdomain%sSubDomain(mat)%hprimez,Tdomain%sSubDomain(mat)%hTprimez)
            call add_previous_state2forces (Tdomain%specel(n), Dt)
            if (Tdomain%specel(n)%PML) call update_Psi_ADEPML(Tdomain%specel(n), &
                            Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimez, Dt)
        enddo

        ! External Forces computation
        call Compute_External_Forces(Tdomain,timelocal)

        ! Compute second member "R" of the Lagrange multiplicator system
        do n=0,Tdomain%n_elem-1
            ! Compute current element contribution to R
            call compute_smbr_R(Tdomain%specel(n))
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
            ! Local Solver
            call local_solver(Tdomain%specel(n),Dt)
        enddo

        iter = iter+1
    enddo


    return
end subroutine Newmark_PMC



!>
!!\brief This Subroutine performs a midpoint method for advancing in time
!! using an explicit approach, contrary to the previous one.
!!\version 1.0
!!\date 20/11/2014
!! This subroutine is used only HDG elements
!! \param type (Domain), intent (INOUT) Tdomain
!! \param real         , intent (IN)    Dt
!<
subroutine Newmark_PMC_explicit (Tdomain,Dt,n_it_max)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: dt
    integer, intent(in)   :: n_it_max


    ! local variables
    integer :: n, mat, iter, nf, nface
    real :: timelocal

    ! Initialization Phase
    do n=0,Tdomain%n_elem-1
        ! Computation of the  prediction :
        Tdomain%specel(n)%Strain0(:,:,:) = Tdomain%specel(n)%Strain(:,:,:)
        Tdomain%specel(n)%V0(:,:,:)      = Tdomain%specel(n)%Veloc (:,:,:)
        if (Tdomain%specel(n)%PML) call initialize_Psi(Tdomain%specel(n))
    enddo

    ! Midpoint method :
    iter= 0

    do while (iter<n_it_max)

        ! Local time is tn for the first step, and tn+1/2 for the next ones
        if (iter==0) then
            timelocal = Tdomain%TimeD%rtime
        else
            timelocal = Tdomain%TimeD%rtime + 0.5*Dt
        endif

        ! Prediction Phase :
        do n=0,Tdomain%n_elem-1
            Tdomain%specel(n)%Strain = 0.5 * (Tdomain%specel(n)%Strain0 + Tdomain%specel(n)%Strain)
            Tdomain%specel(n)%Veloc  = 0.5 * (Tdomain%specel(n)%V0      + Tdomain%specel(n)%Veloc )
            !if (Tdomain%specel(n)%ADEPML) call Prediction_Psi(Tdomain%specel(n))
        enddo

        ! Building second members (= forces) of systems.
        do n=0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            call compute_InternalForces_DG_Weak(Tdomain%specel(n), &
                Tdomain%sSubDomain(mat)%hprimex, &
                Tdomain%sSubDomain(mat)%hTprimez)
            call compute_TracFace (Tdomain%specel(n))
            if (Tdomain%specel(n)%PML) call update_Psi_ADEPML(Tdomain%specel(n), &
                            Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimez, Dt)
        enddo

        ! External Forces computation
        call Compute_External_Forces(Tdomain,timelocal)

        ! Calcul des fluxs
        do n = 0, Tdomain%n_elem-1
            do nf = 0,3
                nface  = Tdomain%specel(n)%Near_Face(nf)
                call get_traction_el2f(Tdomain,n,nface,nf,timelocal)
            enddo
        enddo

        ! Constructing the Lambda (= velocities vhat) on the faces
        do n=0,Tdomain%n_elem-1
            do nf = 0,3        ! Computation of the Velocities Traces
                nface = Tdomain%specel(n)%Near_Face(nf)
                call Compute_Vhat(Tdomain%sFace(nface))
            enddo
            call get_Vhat_f2el(Tdomain,n)
            call Compute_Traces (Tdomain%specel(n),.true.)
            if(Tdomain%type_bc==DG_BC_REFL) call enforce_diriclet_BC(Tdomain,n)
            ! Add previous state to forces :
            call add_previous_state2forces (Tdomain%specel(n),Dt)
            Tdomain%specel(n)%Forces(:,:,2) = Tdomain%specel(n)%Forces(:,:,2) &
                                    - 1./Dt * Tdomain%specel(n)%Acoeff(:,:,12) * Tdomain%specel(n)%Strain0(:,:,2)
            if(Tdomain%type_bc==DG_BC_REFL) call enforce_diriclet_BC(Tdomain,n)
            call inversion_massmat(Tdomain%specel(n))
            Tdomain%specel(n)%Strain(:,:,:) = Dt * Tdomain%specel(n)%Forces(:,:,0:2)
            Tdomain%specel(n)%Veloc (:,:,:) = Dt * Tdomain%specel(n)%Forces(:,:,3:4)
        enddo

        iter = iter+1
    enddo


    return
end subroutine Newmark_PMC_explicit


end module snewmark_pmc
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
