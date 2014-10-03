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
    integer :: ns, ncc,i,j,n,np, ngllx, ngllz, mat, nelem,nf, w_face, nv_aus, nf_aus, nv
    integer :: n_face_pointed, tag_send, tag_receive, i_send, i_stock, ngll, ierr, i_proc
    integer, dimension (MPI_STATUS_SIZE) :: status
    real :: bega, gam1,alpha,dt

    real, dimension (0:1) :: V_free_vertex
    real, dimension (:,:), allocatable :: Vxloc, Vzloc, V_free
    real, dimension (:,:), allocatable :: Smooth


    ! Predictor-MultiCorrector Newmark Velocity Scheme within a
    ! Time staggered Stress-Velocity formulation inside PML
    alpha =Tdomain%TimeD%alpha
    bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
    gam1 = 1. / Tdomain%TimeD%gamma


    ! Predictor Phase

    do n=0,Tdomain%n_elem-1
        ! Dans le tableau Vect_RK sont stockes les deformations et vitesses de l'instant tn
        Elem%S0(:,:,:) = Elem%Strain(:,:,:)
        Elem%V0(:,:,:) = Elem%Veloc(:,:,:)
        call Prediction_Elem_NPMC (Tdomain%specel(n))
    enddo

    ! Multicorrector phase :
    i = 0
    do while (resid >= something)
        do n=0,Tdomain%n_elem-1
            call compute_internal_forces ()
        enddo
        i = i+1
    enddo


    return
end subroutine Newmark_DG

end module snewmark_dg
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
