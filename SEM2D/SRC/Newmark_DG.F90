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
    integer :: i,j,n,np, ngllx, ngllz, mat, nelem,nf, w_face, nv_aus, nf_aus, nv, pos, iter
    integer :: n_face_pointed, tag_send, tag_receive, i_send, i_stock, ngll, ierr, i_proc
    integer, dimension (MPI_STATUS_SIZE) :: status
    real :: bega, gam1,alpha,dt

    real, dimension (:,:), allocatable :: Vxloc, Vzloc
    real, dimension (:,:), allocatable :: Smooth


    ! Predictor-MultiCorrector Newmark Velocity Scheme within a
    ! Time staggered Stress-Velocity formulation inside PML
    alpha =Tdomain%TimeD%alpha
    bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
    gam1 = 1. / Tdomain%TimeD%gamma
    n_it_max = 5

    ! Predictor Phase

    do n=0,Tdomain%n_elem-1
        ! Dans le tableau Vect_RK sont stockes les deformations et vitesses de l'instant tn
        Elem%S0(:,:,:) = Elem%Strain(:,:,:)
        Elem%V0(:,:,:) = Elem%Veloc(:,:,:)
        call Prediction_Elem_NPMC (Tdomain%specel(n))
    enddo

    ! Multicorrector phase :
    iter= 0

    do while (iter<n_it_max)

        do n=0,Tdomain%n_elem-1
            call compute_internal_forces ()
        enddo

        ! Send informations on vertices
        do n=0,Tdomain%n_face-1
            do i=0,1
                smbr = ...
                nv   = Tdomain%sFace(n)%Near_Vertex(i)
                pos  = Tdomain%sFace(n)%pos_in_VertMat(i)
                Tdomain%sVertex(nv)%smbrLambda((2*pos):(2*pos+1)) = smbr
            enddo
        enddo

        ! Solve linear systems on the vertices
        do n=0,Tdomain%n_vertex-1
            call solve_gradient_conjugue()
        enddo

        ! Constructing the Lambdas on the faces
        do n=0,Tdomain%n_face-1
            ngll = Tdomain%sFace(n)%ngll
            ! Faces ends (vertices)
            do i=0,1
                nv   = Tdomain%sFace(n)%Near_Vertex(i)
                pos  = Tdomain%sFace(n)%pos_in_VertMat(i)
                Tdomain%sFace(n)%Veloc(i*(ngll-1),:) = Tdomain%sVertex(nv)%Lambda((2*pos):(2*pos+1))
            enddo
            ! Faces inner nodes
            do i=1,ngll-1
                Tdomain%sFace(n)%Veloc(i) = .....
            enddo
        enddo

        ! Communication Lamdas from faces to elements
        do n=0,Tdomain%n_elem-1
            do nf = 0,3
                nface = Tdomain%specel(n)%Near_Face(nf)
                call get_Vhat_f2el(Tdomain,n,nface,nf)
            enddo
        enddo

        !!!!!!!!!!!!!  MPI COMMUNICATIONS HERE  !!!!!!!!!!!!!!!!!!!

        ! Solveurs locaux : gets Veloc and Strain from Lambda
        do n=0,Tdomain%n_elem-1
            ! Solveur Local
            ! Solveur Local
        enddo
        iter = iter+1
    enddo


    return
end subroutine Newmark_DG

end module snewmark_dg
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
