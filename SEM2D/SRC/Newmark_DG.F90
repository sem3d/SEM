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

subroutine Newmark_PMC (Tdomain,Dt)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: dt

    ! local variables
    integer :: n, mat, nf, nface, iter, n_it_max
    real :: bega, gam1, alpha, timelocal
    !integer, dimension (MPI_STATUS_SIZE) :: status

    ! Predictor-MultiCorrector Newmark Velocity Scheme within a
    ! Time staggered Stress-Velocity formulation inside PML
    ! #################################################### !
    ! NO NO NO !!! SO FAR IT IS A SIMPLE MIDPOINT METHOD !!!
    ! #################################################### !

    alpha =Tdomain%TimeD%alpha
    bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
    gam1 = 1. / Tdomain%TimeD%gamma
    n_it_max = 4
    timelocal = Tdomain%TimeD%rtime + Dt


    ! Initialization Phase
    do n=0,Tdomain%n_elem-1
        ! Computation of the  prediction :
        Tdomain%specel(n)%Strain0(:,:,:) = Tdomain%specel(n)%Strain(:,:,:)
        Tdomain%specel(n)%V0(:,:,:)      = Tdomain%specel(n)%Veloc (:,:,:)
        !mat = Tdomain%specel(n)%mat_index
        !call Prediction_Elem_NPMC (Tdomain%specel(n), Tdomain%sSubDomain(mat)%hTprimex, &
        !                           Tdomain%sSubDomain(mat)%hprimez, Dt)
    enddo

    ! Midpoint method :
    iter= 0

    do while (iter<n_it_max)

        ! Prediction Phase :
        do n=0,Tdomain%n_elem-1
            if (iter == 0) then
                Tdomain%specel(n)%Strain = 0.5 * Tdomain%specel(n)%Strain0
                Tdomain%specel(n)%Veloc  = 0.5 * Tdomain%specel(n)%V0
            else
                Tdomain%specel(n)%Strain = 0.5 * (Tdomain%specel(n)%Strain0 + Tdomain%specel(n)%Strain)
                Tdomain%specel(n)%Veloc  = 0.5 * (Tdomain%specel(n)%V0      + Tdomain%specel(n)%Veloc )
            endif
        enddo

        ! Building second members (= forces) of systems.
        do n=0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            call compute_InternalForces_HDG  (Tdomain%specel(n), &
              Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimex, &
              Tdomain%sSubDomain(mat)%hprimez,Tdomain%sSubDomain(mat)%hTprimez)
            call add_previous_state2forces (Tdomain%specel(n), Dt)
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
            do nf = 0,3
                nface = Tdomain%specel(n)%Near_Face(nf)
                call get_Vhat_f2el(Tdomain,n,nface,nf)
            enddo
            ! Local Solver
            call local_solver(Tdomain%specel(n),Dt)
        enddo

        iter = iter+1
    enddo


    return
end subroutine Newmark_PMC

end module snewmark_pmc
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
