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

subroutine Runge_Kutta4 (Tdomain, ntime, dt)
    use sdomain
    use scouplage
    use mpi

    implicit none
    type (domain), intent (INOUT) :: Tdomain
    integer, intent(in) :: ntime
    real,    intent(in)   :: dt

    ! local variables
    integer :: ns, i, j, k, n, np, mat, nelem, nf
    integer :: n_face_pointed, tag_send, tag_receive, i_send, i_stock, ngll, ierr, i_proc
    integer, dimension (MPI_STATUS_SIZE) :: status
    integer               :: nface,  type_DG
    logical               :: coherency
    real                  :: timelocal
    real, dimension(3)    :: coeffs


    ! Runge-Kutta Initialization
    do n = 0, Tdomain%n_elem-1
       if(Tdomain%specel(n)%Type_DG==2) then
          Tdomain%specel(n)%Vect_RK = Tdomain%specel(n)%Veloc
       else
          Tdomain%specel(n)%Vect_RK(:,:,0:2) = Tdomain%specel(n)%Strain
          Tdomain%specel(n)%Vect_RK(:,:,3:4) = Tdomain%specel(n)%Veloc
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
          case(0) ! Discontinuous Galerkin Strong Formulation
             call compute_InternalForces_DG_Strong(Tdomain%specel(n), &
                                                   Tdomain%sSubDomain(mat)%hTprimex, &
                                                   Tdomain%sSubDomain(mat)%hprimez)
          case(1) ! Discontinuous Galerkin Weak Formulation
             call compute_InternalForces_DG_Weak(Tdomain%specel(n), &
                                                 Tdomain%sSubDomain(mat)%hprimex, &
                                                 Tdomain%sSubDomain(mat)%hTprimez)
          case(2) ! Continuous Galerkin
             call get_Displ_fv2el (Tdomain,n)
             call compute_InternalForces_Elem(Tdomain%specel(n), &
                                              Tdomain%sSubDomain(mat)%hprimex, &
                                              Tdomain%sSubDomain(mat)%hTprimex, &
                                              Tdomain%sSubDomain(mat)%hprimez, &
                                              Tdomain%sSubDomain(mat)%hTprimez)
          end select
          ! Calcul des fluxs / Assemblage des forces
          do nf = 0,3
             nface  = Tdomain%specel(n)%Near_Face(nf)
             if(type_DG == 2) then
                call Assemblage(Tdomain,n,nface,nf)
             else
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
          call MPI_SEND (Tdomain%sWall(i_proc)%Send_data_2,2* Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION, i_send, &
               tag_send, Tdomain%communicateur, ierr )
          call MPI_RECV (Tdomain%sWall(i_proc)%Receive_data_2, 2* Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION, i_send, &
               tag_receive, Tdomain%communicateur, status, ierr )
          call Assign_recv_data(TDomain,i_proc)
       enddo

       do n = 0, Tdomain%n_elem-1
          type_DG = Tdomain%specel(n)%Type_DG
          if (type_DG==2) then  ! Continuous Galerkin
             call inversion_massmat(Tdomain%specel(n))
             Tdomain%specel(n)%Vect_RK = coeffs(1) * Tdomain%specel(n)%Vect_RK  + Tdomain%specel(n)%Forces * dt
             Tdomain%specel(n)%Veloc   = Tdomain%specel(n)%Veloc + coeffs(2) * Tdomain%specel(n)%Vect_RK
             Tdomain%specel(n)%Displ   = Tdomain%specel(n)%Displ + Tdomain%specel(n)%Veloc * dt ! ATTENTION : A REVOIR !!!!!
             Tdomain%specel(n)%Forces  = Tdomain%specel(n)%Displ
          else                  ! Discontinuous Galerkin
             do nf = 0,3        ! Computation of the fluxes
                nface = Tdomain%specel(n)%Near_Face(nf)
                call Compute_Flux(Tdomain%sFace(nface),n,type_DG)
                call get_flux_f2el(Tdomain,n,nface,nf)
             enddo
             call inversion_massmat(Tdomain%specel(n))
             Tdomain%specel(n)%Vect_RK = coeffs(1) * Tdomain%specel(n)%Vect_RK + Tdomain%specel(n)%Forces * dt
             Tdomain%specel(n)%Strain  = Tdomain%specel(n)%Strain + coeffs(2) * Tdomain%specel(n)%Vect_RK(:,:,0:2)
             Tdomain%specel(n)%Veloc   = Tdomain%specel(n)%Veloc  + coeffs(2) * Tdomain%specel(n)%Vect_RK(:,:,3:4)
          endif
       enddo

       ! Computing new values on the vertexes and the faces for Continuous Galerkin only
       do n=0, Tdomain%n_face-1
          type_DG = Tdomain%sface(n)%Type_DG
          if (type_DG == 2) then
             Tdomain%sface(n)%Forces(:,0)  = Tdomain%sface(n)%MassMat(:) * Tdomain%sface(n)%Forces(:,0)
             Tdomain%sface(n)%Forces(:,1)  = Tdomain%sface(n)%MassMat(:) * Tdomain%sface(n)%Forces(:,1)
             Tdomain%sface(n)%Vect_RK = coeffs(1) * Tdomain%sface(n)%Vect_RK + Tdomain%sface(n)%Forces * dt
             Tdomain%sface(n)%Veloc   = Tdomain%sface(n)%Veloc + coeffs(2) * Tdomain%sface(n)%Vect_RK
             Tdomain%sface(n)%Displ   = Tdomain%sface(n)%Displ + Tdomain%sface(n)%Veloc * dt ! ATTENTION : A REVOIR !!!!!
             Tdomain%sface(n)%Forces  = Tdomain%sface(n)%Displ
             Tdomain%sface(n)%is_computed = .false.
          endif
       enddo
       do n=0, Tdomain%n_vertex-1
          type_DG = Tdomain%sVertex(n)%Type_DG
          if (type_DG == 2) then
             Tdomain%sVertex(n)%Forces  = Tdomain%sVertex(n)%MassMat * Tdomain%sVertex(n)%Forces
             Tdomain%sVertex(n)%Vect_RK = coeffs(1) * Tdomain%sVertex(n)%Vect_RK + Tdomain%sVertex(n)%Forces * dt
             Tdomain%sVertex(n)%Veloc   = Tdomain%sVertex(n)%Veloc + coeffs(2) * Tdomain%sVertex(n)%Vect_RK
             Tdomain%sVertex(n)%Displ   = Tdomain%sVertex(n)%Displ + Tdomain%sVertex(n)%Veloc * dt ! ATTENTION : A REVOIR !!!!!
             Tdomain%sVertex(n)%Forces  = Tdomain%sVertex(n)%Displ
             Tdomain%sVertex(n)%is_computed = .false.
          endif
       enddo

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

