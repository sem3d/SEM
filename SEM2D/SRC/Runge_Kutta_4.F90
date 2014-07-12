!>
!!\file Runge_Kutta4.F90
!!\brief Algorithm for a Low Storage Explicit Runge Kutta (LSERK) scheme of order 4
!!\version 1.0
!!\date 15/11/2013
!! This algorithm comes from the following book :
!! Discontinuous Galerkin Methods, Algorithms, Analysis and Applications
!! by Jan S. Hesthaven & Tim Warburton, 2010, Springer-Verlag.
!! cf Algorithm pages 63-64
!<

module srungekutta
    use sdomain
    use scouplage
    use mpi
    use constants

    implicit none
contains

subroutine Runge_Kutta4 (Tdomain, dt)
    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real,    intent(in)   :: dt

    ! local variables
    integer :: i, n, mat, nf, nv, ngx, ngz
    integer :: tag_send, tag_receive, i_send, ierr, i_proc
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer               :: nface,  type_DG
    real                  :: timelocal
    real, dimension(3)    :: coeffs
    real, dimension (:,:), allocatable :: Vxloc, Vzloc
    logical :: acoustic


    ! Runge-Kutta Initialization
    do n = 0, Tdomain%n_elem-1
       if(Tdomain%specel(n)%Type_DG==GALERKIN_CONT) then
           ngx = Tdomain%specel(n)%ngllx
           ngz = Tdomain%specel(n)%ngllz
           Tdomain%specel(n)%Vect_RK(1:ngx-2,1:ngz-2,0:1) = Tdomain%specel(n)%Veloc(:,:,0:1)
           Tdomain%specel(n)%Vect_RK(:,:,2:4) = Tdomain%specel(n)%Stress(:,:,0:2)
          !Tdomain%specel(n)%Vect_RK(:,:,0:1) = Tdomain%specel(n)%Veloc(:,:,0:1)
          !Tdomain%specel(n)%Vect_RK(:,:,2:3) = Tdomain%specel(n)%Displ(:,:,0:1)
          do i=0,3
              nf = Tdomain%specel(n)%Near_Face(i)
              nv = Tdomain%specel(n)%Near_Vertex(i)
              Tdomain%sFace(nf)%Vect_RK(:,0:1) = Tdomain%sFace(nf)%Veloc(:,0:1)
              Tdomain%sVertex(nv)%Vect_RK(0:1) = Tdomain%sVertex(nv)%Veloc(0:1)
              !Tdomain%sFace(nf)%Vect_RK(:,0:1) = Tdomain%sFace(nf)%Veloc(:,0:1)
              !Tdomain%sFace(nf)%Vect_RK(:,2:3) = Tdomain%sFace(nf)%Displ(:,0:1)
              !Tdomain%sVertex(nv)%Vect_RK(0:1) = Tdomain%sVertex(nv)%Veloc(0:1)
              !Tdomain%sVertex(nv)%Vect_RK(2:3) = Tdomain%sVertex(nv)%Displ(0:1)
          enddo
       else
          Tdomain%specel(n)%Vect_RK(:,:,0:2) = Tdomain%specel(n)%Strain
          Tdomain%specel(n)%Vect_RK(:,:,3:4) = Tdomain%specel(n)%Veloc
       endif
       if(Tdomain%specel(n)%ADEPML) then
           Tdomain%specel(n)%Psi_RK(:,:,0) = Tdomain%specel(n)%PsiSxxx(:,:)
           Tdomain%specel(n)%Psi_RK(:,:,1) = Tdomain%specel(n)%PsiSzzz(:,:)
           Tdomain%specel(n)%Psi_RK(:,:,2) = Tdomain%specel(n)%PsiSxzx(:,:)
           Tdomain%specel(n)%Psi_RK(:,:,3) = Tdomain%specel(n)%PsiSxzz(:,:)
           Tdomain%specel(n)%Psi_RK(:,:,4) = Tdomain%specel(n)%PsiVxx(:,:)
           Tdomain%specel(n)%Psi_RK(:,:,5) = Tdomain%specel(n)%PsiVxz(:,:)
           Tdomain%specel(n)%Psi_RK(:,:,6) = Tdomain%specel(n)%PsiVzx(:,:)
           Tdomain%specel(n)%Psi_RK(:,:,7) = Tdomain%specel(n)%PsiVzz(:,:)
       endif
    enddo

    do i = 1,5
       ! Calcul des Coefficients RK low storage
       coeffs = Coeffs_LSERK(i)
       timelocal = Tdomain%TimeD%rtime + dt * coeffs(3)

       do n = 0, Tdomain%n_elem-1
          type_DG = Tdomain%specel(n)%Type_DG
          mat = Tdomain%specel(n)%mat_index
          select case (type_DG)
          case(GALERKIN_DG_STRONG) ! Discontinuous Galerkin Strong Formulation
             call compute_InternalForces_DG_Strong(Tdomain%specel(n), &
                                                   Tdomain%sSubDomain(mat)%hTprimex, &
                                                   Tdomain%sSubDomain(mat)%hprimez)
          case(GALERKIN_DG_WEAK) ! Discontinuous Galerkin Weak Formulation
             call compute_InternalForces_DG_Weak(Tdomain%specel(n), &
                                                 Tdomain%sSubDomain(mat)%hprimex, &
                                                 Tdomain%sSubDomain(mat)%hTprimez)
         case(GALERKIN_HDG_RP)   ! Hybridizable Discontinuous Galerkin
             call compute_InternalForces_DG_Weak(Tdomain%specel(n), &
                                                 Tdomain%sSubDomain(mat)%hprimex, &
                                                 Tdomain%sSubDomain(mat)%hTprimez)
             call compute_TracFace (Tdomain%specel(n))
          case(GALERKIN_CONT) ! Continuous Galerkin
             ngx = Tdomain%specel(n)%ngllx
             ngz = Tdomain%specel(n)%ngllz
             allocate (Vxloc(0:ngx-1, 0:ngz-1))
             allocate (Vzloc(0:ngx-1, 0:ngz-1))
             call get_PMLprediction_fv2el (Tdomain,n,Vxloc,vzloc,ngx,ngz,0.5,0.5,dt)
             call compute_2ndMember_Veloc_Stress(Tdomain%specel(n), Vxloc, Vzloc, &
                                                 Tdomain%sSubDomain(mat)%hprimex, &
                                                 Tdomain%sSubDomain(mat)%hTprimex, &
                                                 Tdomain%sSubDomain(mat)%hprimez, &
                                                 Tdomain%sSubDomain(mat)%hTprimez)
             deallocate (VxLoc, Vzloc)
             !call get_RealDispl_fv2el (Tdomain,n)
             !call compute_InternalForces_Elem(Tdomain%specel(n), &
             !                                 Tdomain%sSubDomain(mat)%hprimex, &
             !                                 Tdomain%sSubDomain(mat)%hTprimex, &
             !                                 Tdomain%sSubDomain(mat)%hprimez, &
             !                                 Tdomain%sSubDomain(mat)%hTprimez)
          end select
          ! Calcul des fluxs / Assemblage des forces
          do nf = 0,3
              nface  = Tdomain%specel(n)%Near_Face(nf)
              if(type_DG == GALERKIN_CONT) then
                  call Assemblage(Tdomain,n,nface,nf)
              elseif(type_DG == GALERKIN_HDG_RP) then
                  call get_traction_el2f(Tdomain,n,nface,nf)
              else  ! Usual Discontinuous Galerkin
                  call get_data_el2f(Tdomain,n,nface,nf)
              endif
          enddo
       enddo

       ! External Forces computation
       call Compute_External_Forces(Tdomain,timelocal)


       ! Communications MPI
       do i_proc = 0, Tdomain%n_communications - 1
          i_send = Tdomain%Communication_list (i_proc)
          call Create_send_data(TDomain,i_proc)
          tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank + 900
          tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 900
          call MPI_SEND (Tdomain%sWall(i_proc)%Send_data_2,2* Tdomain%sWall(i_proc)%n_points, &
               MPI_DOUBLE_PRECISION, i_send, tag_send, Tdomain%communicateur, ierr )
          call MPI_RECV (Tdomain%sWall(i_proc)%Receive_data_2,2*Tdomain%sWall(i_proc)%n_points,&
               MPI_DOUBLE_PRECISION, i_send, tag_receive, Tdomain%communicateur, status, ierr )
          call Assign_recv_data(TDomain,i_proc)
       enddo

       do n = 0, Tdomain%n_elem-1
          mat = Tdomain%specel(n)%mat_index
          type_DG = Tdomain%specel(n)%Type_DG
          if (type_DG==GALERKIN_CONT) then  ! Continuous Galerkin
             call inversion_massmat(Tdomain%specel(n))
             call update_Elem_RK4 (Tdomain%specel(n),coeffs(1),coeffs(2),dt)
          else                  ! Discontinuous Galerkin
             acoustic = Tdomain%specel(n)%Acoustic
             if (type_DG==GALERKIN_HDG_RP) then
                 do nf = 0,3        ! Computation of the Velocities Traces
                     nface = Tdomain%specel(n)%Near_Face(nf)
                     call Compute_Vhat(Tdomain%sFace(nface))
                     call get_Vhat_f2el(Tdomain,n,nface,nf)
                 enddo
                 !if(Tdomain%specel(n)%ADEPML) call enforce_diriclet_corners_vhat(Tdomain,n)
                 call Compute_Traces (Tdomain%specel(n))
             elseif (type_DG==GALERKIN_DG_WEAK .OR. type_DG==GALERKIN_DG_STRONG) then
                 do nf = 0,3        ! Computation of the fluxes
                     nface = Tdomain%specel(n)%Near_Face(nf)
                     call Compute_Flux(Tdomain%sFace(nface),n,type_DG,acoustic)
                     call get_flux_f2el(Tdomain,n,nface,nf)
                 enddo
             endif
             if(Tdomain%specel(n)%ADEPML) call add_Psi4PML(Tdomain%specel(n))
             if(Tdomain%type_bc==DG_BC_REFL) call enforce_diriclet_BC(Tdomain,n)
             call inversion_massmat(Tdomain%specel(n))
             if(Tdomain%specel(n)%ADEPML) call update_Psi_RK4(Tdomain%specel(n), &
                                                              Tdomain%sSubDomain(mat)%hprimex,  &
                                                              Tdomain%sSubDomain(mat)%hTprimex, &
                                                              Tdomain%sSubDomain(mat)%hprimez,  &
                                                              Tdomain%sSubDomain(mat)%hTprimez, &
                                                              coeffs(1),coeffs(2),dt)
             call update_Elem_RK4 (Tdomain%specel(n),coeffs(1),coeffs(2),dt)
          endif
       enddo

       ! Updates of faces and vertices if continuous elements
       if(Tdomain%type_elem==GALERKIN_CONT) call Update_FV_RK4 (Tdomain,coeffs(1),coeffs(2),dt)

    enddo ! End loop RK4

    return

  contains

    function Coeffs_LSERK(i)

      integer :: i
      real, dimension(3) :: Coeffs_LSERK

      select case (i)
      case (1)
         Coeffs_LSERK(1) = 0.
         Coeffs_LSERK(2) = 1432997174477./9575080441755.
         Coeffs_LSERK(3) = 0.
      case (2)
         Coeffs_LSERK(1) = -567301805773./1357537059087.
         Coeffs_LSERK(2) = 5161836677717./13612068292357.
         Coeffs_LSERK(3) = 1432997174477./9575080441755.
      case (3)
         Coeffs_LSERK(1) = -2404267990393./2016746695238.
         Coeffs_LSERK(2) = 1720146321549. /2090206949498.
         Coeffs_LSERK(3) = 2526269341429. /6820363962896.
      case (4)
         Coeffs_LSERK(1) = -3550918686646./2091501179385.
         Coeffs_LSERK(2) = 3134564353537. /4481467310338.
         Coeffs_LSERK(3) = 2006345519317. /3224310063776.
      case (5)
         Coeffs_LSERK(1) = -1275806237668./842570457699.
         Coeffs_LSERK(2) = 2277821191437. /14882151754819.
         Coeffs_LSERK(3) = 2802321613138. /2924317926251.
      end select

    end function Coeffs_LSERK

  end subroutine Runge_Kutta4


  ! ###########################################################
  !>
  !! \brief This subroutine updates the faces and the vertices of the domain
  !! in a Runge Kutta framework if the domain contains any Continuous Galerkin element.
  !!
  !! \param type (Domain), intent (INOUT) Tdomain
  !<
  subroutine Update_FV_RK4 (Tdomain,coeff1,coeff2,Dt)

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real, intent(IN) :: coeff1
    real, intent(IN) :: coeff2
    real, intent(IN) :: Dt
    ! local variables
    integer :: n, type_DG

    do n=0, Tdomain%n_face-1
       type_DG = Tdomain%sface(n)%Type_DG
       if (type_DG == GALERKIN_CONT) then
          ! Inversion of Mass Matrix to get 2nd member of velocity equation
          Tdomain%sface(n)%Forces(:,0)  = Tdomain%sface(n)%MassMat(:) * Tdomain%sface(n)%Forces(:,0)
          Tdomain%sface(n)%Forces(:,1)  = Tdomain%sface(n)%MassMat(:) * Tdomain%sface(n)%Forces(:,1)
          ! RK4 Updates of velocities and displacements
          Tdomain%sface(n)%Vect_RK(:,0:1) = coeff1 * Tdomain%sface(n)%Vect_RK(:,0:1) + Dt * Tdomain%sface(n)%Forces(:,0:1)
          Tdomain%sface(n)%Vect_RK(:,2:3) = coeff1 * Tdomain%sface(n)%Vect_RK(:,2:3) + Dt * Tdomain%sface(n)%Veloc (:,0:1)
          Tdomain%sface(n)%Veloc   = Tdomain%sface(n)%Veloc + coeff2 * Tdomain%sface(n)%Vect_RK(:,0:1)
          Tdomain%sface(n)%Displ   = Tdomain%sface(n)%Displ + coeff2 * Tdomain%sface(n)%Vect_RK(:,2:3)
          Tdomain%sface(n)%is_computed = .false.
       endif
    enddo
    do n=0, Tdomain%n_vertex-1
       type_DG = Tdomain%sVertex(n)%Type_DG
       if (type_DG == GALERKIN_CONT) then
          ! Inversion of Mass Matrix to get 2nd member of velocity equation
          Tdomain%sVertex(n)%Forces  = Tdomain%sVertex(n)%MassMat * Tdomain%sVertex(n)%Forces
          ! RK4 Updates of velocities and displacements
          Tdomain%sVertex(n)%Vect_RK(0:1) = coeff1 * Tdomain%sVertex(n)%Vect_RK(0:1) + Dt * Tdomain%sVertex(n)%Forces
          Tdomain%sVertex(n)%Vect_RK(2:3) = coeff1 * Tdomain%sVertex(n)%Vect_RK(2:3) + Dt * Tdomain%sVertex(n)%Veloc
          Tdomain%sVertex(n)%Veloc   = Tdomain%sVertex(n)%Veloc + coeff2 * Tdomain%sVertex(n)%Vect_RK(0:1)
          Tdomain%sVertex(n)%Displ   = Tdomain%sVertex(n)%Displ + coeff2 * Tdomain%sVertex(n)%Vect_RK(2:3)
          Tdomain%sVertex(n)%is_computed = .false.
       endif
    enddo

  end subroutine Update_FV_RK4

end module srungekutta
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
