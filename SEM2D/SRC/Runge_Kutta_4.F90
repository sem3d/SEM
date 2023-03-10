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
    use mpi
    use constants
    use smidpoint

    implicit none
contains

subroutine Runge_Kutta4 (Tdomain, dt)
    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real(fpp),    intent(in)   :: dt

    ! local variables
    integer :: i, n, mat
    integer :: tag_send, tag_receive, i_send, ierr, i_proc
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer               :: nface,  type_DG
    real(fpp)             :: timelocal
    real(fpp), dimension(3) :: coeffs


    do i = 1,5
       ! Calcul des Coefficients RK low storage
       coeffs = Coeffs_LSERK(i)
       timelocal = Tdomain%TimeD%rtime + dt * coeffs(3)

       do n = 0, Tdomain%n_elem-1
          type_DG = Tdomain%specel(n)%Type_DG
          mat = Tdomain%specel(n)%mat_index
          select case (type_DG)
          case(GALERKIN_CONT) ! Continuous Galerkin
             call get_RealDispl_fv2el (Tdomain,n)
             call compute_InternalForces_Elem (Tdomain%specel(n), &
                 Tdomain%sSubDomain(mat)%hprimex,  &
                 Tdomain%sSubDomain(mat)%hTprimex, &
                 Tdomain%sSubDomain(mat)%hprimez,  &
                 Tdomain%sSubDomain(mat)%hTprimez)
         case(GALERKIN_HDG_RP)   ! Hybridizable Discontinuous Galerkin
             call compute_InternalForces_HDG_Weak(Tdomain%specel(n), &
                                                 Tdomain%sSubDomain(mat)%hprimex, &
                                                 Tdomain%sSubDomain(mat)%hTprimez)
          case(GALERKIN_DG_STRONG) ! Discontinuous Galerkin Strong Formulation
             call compute_InternalForces_DG_Strong(Tdomain%specel(n), &
                                                   Tdomain%sSubDomain(mat)%hTprimex, &
                                                   Tdomain%sSubDomain(mat)%hprimez)
          case(GALERKIN_DG_WEAK) ! Discontinuous Galerkin Weak Formulation
             call compute_InternalForces_DG_Weak(Tdomain%specel(n), &
                                                 Tdomain%sSubDomain(mat)%hprimex, &
                                                 Tdomain%sSubDomain(mat)%hTprimez)
          end select
          if(Tdomain%specel(n)%ADEPML) call add_Psi4PML(Tdomain%specel(n))
       enddo

       ! External Forces computation
       call Compute_External_Forces(Tdomain,timelocal)

       ! Envoi des informations sur les faces
       do n = 0, Tdomain%n_elem-1
          type_DG = Tdomain%specel(n)%Type_DG
          select case (type_DG)
          case (GALERKIN_CONT)
              call Assemblage(Tdomain,n)
          case (GALERKIN_HDG_RP)
              if (i .NE. 1) call compute_TracFace (Tdomain%specel(n))
              if (i .NE. 1) call get_traction_el2f(Tdomain,n)
          case (GALERKIN_DG_STRONG)
              call get_data_el2f(Tdomain,n)
          case (GALERKIN_DG_WEAK)
              call get_data_el2f(Tdomain,n)
          end select
       enddo

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

       ! Face resolutions
       do nface = 0,Tdomain%n_face-1
           type_DG = Tdomain%sface(nface)%Type_DG
           select case (type_DG)
           case (GALERKIN_HDG_RP)
               if (i .NE. 1) call Compute_Vhat_Face_Expl(Tdomain%sFace(nface))
           case(GALERKIN_DG_STRONG)
               call Compute_Flux_DGstrong(Tdomain%sFace(nface))
           case(GALERKIN_DG_WEAK)
               call Compute_Flux_DGweak(Tdomain%sFace(nface))
           case(COUPLE_CG_HDG)
               call Compute_Flux_Coupling(Tdomain%sFace(nface))
           end select
       enddo

       ! Send Informations back from faces to elements
       do n = 0, Tdomain%n_elem-1
           type_DG = Tdomain%specel(n)%Type_DG
           select case (type_DG)
           case (GALERKIN_HDG_RP)
               call get_Vhat_f2el(Tdomain,n)
               call Compute_Traces (Tdomain%specel(n))
               if(Tdomain%type_bc==DG_BC_REFL .OR. Tdomain%specel(n)%PML) call enforce_diriclet_BC(Tdomain,n)
           case(GALERKIN_DG_STRONG)
               call get_flux_f2el_DGstrong(Tdomain,n)
           case(GALERKIN_DG_WEAK)
               call get_flux_f2el(Tdomain,n)
           end select
       enddo

       ! Update elements
       do n = 0, Tdomain%n_elem-1
          mat = Tdomain%specel(n)%mat_index
          type_DG = Tdomain%specel(n)%Type_DG
          call inversion_massmat(Tdomain%specel(n))
          if(Tdomain%specel(n)%ADEPML) call update_Psi_ADEPML_RK4(Tdomain%specel(n), &
                                                              Tdomain%sSubDomain(mat)%hTprimex, &
                                                              Tdomain%sSubDomain(mat)%hprimez,  &
                                                              coeffs(1),coeffs(2),dt)
          call update_Elem_RK4 (Tdomain%specel(n),coeffs(1),coeffs(2),dt)
       enddo

       ! Update Faces and Vertices (for CG only)
       call Update_FV_RK4 (Tdomain,coeffs(1),coeffs(2),dt)

    enddo ! End loop RK4

    ! Calcul Traces a la fin du pas de temps (HDG only)
    if (Tdomain%type_elem == GALERKIN_HDG_RP) &
        call get_vhat_from_current_state(Tdomain)


    return

  contains

    function Coeffs_LSERK(i)

      integer :: i
      real(fpp), dimension(3) :: Coeffs_LSERK

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
    real(fpp), intent(IN) :: coeff1
    real(fpp), intent(IN) :: coeff2
    real(fpp), intent(IN) :: Dt
    ! local variables
    integer :: n, type_DG, ngll


    do n=0, Tdomain%n_face-1
       ngll = Tdomain%sface(n)%ngll
       type_DG = Tdomain%sface(n)%Type_DG
       if (type_DG == GALERKIN_CONT .OR. type_DG == COUPLE_CG_HDG) then
          ! Reflexive Boundary Conditions (if any)
          if (Tdomain%sface(n)%reflex) Tdomain%sface(n)%Forces = 0.
          ! Sends the forces from face to vertex if the face is coupling CG-HDG
          if (type_DG == COUPLE_CG_HDG) call get_forces_f2v (Tdomain, n)
          ! Inversion of Mass Matrix to get 2nd member of velocity equation
          Tdomain%sface(n)%Forces(1:ngll-2,0) = Tdomain%sface(n)%MassMat(:) * Tdomain%sface(n)%Forces(1:ngll-2,0)
          Tdomain%sface(n)%Forces(1:ngll-2,1) = Tdomain%sface(n)%MassMat(:) * Tdomain%sface(n)%Forces(1:ngll-2,1)
          ! RK4 Updates of displacements and velocities (displacements first !)
          Tdomain%sface(n)%Vect_RK(:,2:3) = coeff1 * Tdomain%sface(n)%Vect_RK(:,2:3) &
                                          + Dt *Tdomain%sface(n)%Veloc(1:ngll-2,0:1)
          Tdomain%sface(n)%Vect_RK(:,0:1) = coeff1 * Tdomain%sface(n)%Vect_RK(:,0:1) &
                                          + Dt * Tdomain%sface(n)%Forces(1:ngll-2,0:1)
          Tdomain%sface(n)%Displ(1:ngll-2,0:1) = Tdomain%sface(n)%Displ(1:ngll-2,0:1) &
                                               + coeff2 * Tdomain%sface(n)%Vect_RK(:,2:3)
          Tdomain%sface(n)%Veloc(1:ngll-2,0:1) = Tdomain%sface(n)%Veloc(1:ngll-2,0:1) &
                                               + coeff2 * Tdomain%sface(n)%Vect_RK(:,0:1)
          Tdomain%sface(n)%Forces(:,:) = 0.
       endif
    enddo
    do n=0, Tdomain%n_vertex-1
       type_DG = Tdomain%sVertex(n)%Type_DG
       if (type_DG == GALERKIN_CONT) then
          ! Reflexive Boundary Conditions (if any)
          if (Tdomain%sVertex(n)%reflex) Tdomain%sVertex(n)%Forces = 0.
          ! Inversion of Mass Matrix to get 2nd member of velocity equation
          Tdomain%sVertex(n)%Forces  = Tdomain%sVertex(n)%MassMat * Tdomain%sVertex(n)%Forces
          ! RK4 Updates of displacements and velocities (displacements first !)
          Tdomain%sVertex(n)%Vect_RK(2:3) = coeff1 * Tdomain%sVertex(n)%Vect_RK(2:3) &
                                          + Dt * Tdomain%sVertex(n)%Veloc
          Tdomain%sVertex(n)%Vect_RK(0:1) = coeff1 * Tdomain%sVertex(n)%Vect_RK(0:1) &
                                          + Dt * Tdomain%sVertex(n)%Forces
          Tdomain%sVertex(n)%Displ = Tdomain%sVertex(n)%Displ + coeff2 * Tdomain%sVertex(n)%Vect_RK(2:3)
          Tdomain%sVertex(n)%Veloc = Tdomain%sVertex(n)%Veloc + coeff2 * Tdomain%sVertex(n)%Vect_RK(0:1)
          Tdomain%sVertex(n)%Forces  = 0.
       endif
    enddo

  end subroutine Update_FV_RK4

end module srungekutta
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
