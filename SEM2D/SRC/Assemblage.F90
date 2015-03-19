!>
!!\file Assemblage.F90
!!\brief Contains subroutines for assembling forces and distributing
!! forces and displacements over the neighboring elements of a face
!!\version 2.0
!!\date 15/11/2013
!! Revision complete 10/03/2015
!! This algorithm come from a part of the previous Newmark routine
!! with several modifications.
!<
subroutine Assemblage (Tdomain, nelem)

  use sdomain
  use constants

  implicit none
  type (domain), intent (INOUT) :: Tdomain
  integer, intent(in)   :: nelem
  integer               :: nface, nf
  logical               :: coherency
  integer               :: nv, ngllx, ngllz

  ngllx =  Tdomain%specel(nelem)%ngllx
  ngllz =  Tdomain%specel(nelem)%ngllz

  do nf=0,3
      nface = Tdomain%specel(nelem)%Near_Face(nf)
      coherency = Tdomain%sFace(nface)%coherency
      ! Assemblage of the forces on the Face nface coming from element nelem
      if (nelem==Tdomain%sFace(nface)%Near_Element(0)) then
          call getInternalF_el2f (Tdomain,nelem,nface,nf,.true.)
      else
          call getInternalF_el2f (Tdomain,nelem,nface,nf,coherency)
      endif
      ! Couplage CG-HDG : Vitesses regroupees sur la face
      if (Tdomain%sFace(nface)%type_DG == COUPLE_CG_HDG) then
          call get_veloc_v2f (Tdomain, nface)
      endif
  enddo

  ! Assemblage of the forces of the first of the two vertexes of nface
  nv = Tdomain%specel(nelem)%Near_Vertex(0)
  Tdomain%sVertex(nv)%Forces = Tdomain%sVertex(nv)%Forces + Tdomain%specel(nelem)%Forces(0,0,0:1)
  nv = Tdomain%specel(nelem)%Near_Vertex(1)
  Tdomain%sVertex(nv)%Forces = Tdomain%sVertex(nv)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,0,0:1)
  nv = Tdomain%specel(nelem)%Near_Vertex(2)
  Tdomain%sVertex(nv)%Forces = Tdomain%sVertex(nv)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,ngllz-1,0:1)
  nv = Tdomain%specel(nelem)%Near_Vertex(3)
  Tdomain%sVertex(nv)%Forces = Tdomain%sVertex(nv)%Forces + Tdomain%specel(nelem)%Forces(0,ngllz-1,0:1)

end subroutine Assemblage


!>
!!\brief Subroutine that shares the data (velocities and strains) from an element
!! "nelem" to its face "nface" taking into account the problems of coherency.
!!\version 1.0
!!\date 18/11/2013
!! This subroutine is used only with DG elements
!<
subroutine Get_data_el2f (Tdomain, nelem)

  use sdomain
  implicit none

  type (domain), intent (INOUT) :: Tdomain
  integer, intent(in)   :: nelem

  ! local variables
  integer :: ngll, ngllx, ngllz, i, nface, w_face
  logical :: coherency

  ngllx = Tdomain%specel(nelem)%ngllx
  ngllz = Tdomain%specel(nelem)%ngllz

  do w_face=0,3
      nface = Tdomain%specel(nelem)%Near_Face(w_face)
      ngll  = Tdomain%sFace(nface)%ngll
      coherency = Tdomain%sFace(nface)%coherency
      if (Tdomain%sFace(nface)%Near_Element(0)==nelem) then
          if (w_face == 0 ) then
              Tdomain%sFace(nface)%Veloc_m  = Tdomain%specel(nelem)%Veloc (0:ngll-1,0,0:1)
              Tdomain%sFace(nface)%Strain_m = Tdomain%specel(nelem)%Strain(0:ngll-1,0,0:2)
          else if (w_face == 1 ) then
              Tdomain%sFace(nface)%Veloc_m  = Tdomain%specel(nelem)%Veloc (ngllx-1,0:ngll-1,0:1)
              Tdomain%sFace(nface)%Strain_m = Tdomain%specel(nelem)%Strain(ngllx-1,0:ngll-1,0:2)
          else if (w_face == 2 ) then
              Tdomain%sFace(nface)%Veloc_m  = Tdomain%specel(nelem)%Veloc (0:ngll-1,ngllz-1,0:1)
              Tdomain%sFace(nface)%Strain_m = Tdomain%specel(nelem)%Strain(0:ngll-1,ngllz-1,0:2)
          else
              Tdomain%sFace(nface)%Veloc_m  = Tdomain%specel(nelem)%Veloc (0,0:ngll-1,0:1)
              Tdomain%sFace(nface)%Strain_m = Tdomain%specel(nelem)%Strain(0,0:ngll-1,0:2)
          endif
      else if (coherency) then
          if (w_face == 0 ) then
              Tdomain%sFace(nface)%Veloc_p  = Tdomain%specel(nelem)%Veloc (0:ngll-1,0,0:1)
              Tdomain%sFace(nface)%Strain_p = Tdomain%specel(nelem)%Strain(0:ngll-1,0,0:2)
          else if (w_face == 1 ) then
              Tdomain%sFace(nface)%Veloc_p  = Tdomain%specel(nelem)%Veloc (ngllx-1,0:ngll-1,0:1)
              Tdomain%sFace(nface)%Strain_p = Tdomain%specel(nelem)%Strain(ngllx-1,0:ngll-1,0:2)
          else if (w_face == 2 ) then
              Tdomain%sFace(nface)%Veloc_p  = Tdomain%specel(nelem)%Veloc (0:ngll-1,ngllz-1,0:1)
              Tdomain%sFace(nface)%Strain_p = Tdomain%specel(nelem)%Strain(0:ngll-1,ngllz-1,0:2)
          else
              Tdomain%sFace(nface)%Veloc_p  = Tdomain%specel(nelem)%Veloc (0,0:ngll-1,0:1)
              Tdomain%sFace(nface)%Strain_p = Tdomain%specel(nelem)%Strain(0,0:ngll-1,0:2)
          end if
      else
          if (w_face == 0 ) then
              do i=0,ngll-1
                  Tdomain%sFace(nface)%Veloc_p(i,0:1)  = Tdomain%specel(nelem)%Veloc (ngll-1-i,0,0:1)
                  Tdomain%sFace(nface)%Strain_p(i,0:2) = Tdomain%specel(nelem)%Strain(ngll-1-i,0,0:2)
              end do
          else if (w_face == 1 ) then
              do i=0,ngll-1
                  Tdomain%sFace(nface)%Veloc_p(i,0:1)  = Tdomain%specel(nelem)%Veloc (ngllx-1,ngll-1-i,0:1)
                  Tdomain%sFace(nface)%Strain_p(i,0:2) = Tdomain%specel(nelem)%Strain(ngllx-1,ngll-1-i,0:2)
              end do
          else if (w_face == 2 ) then
              do i=0,ngll-1
                  Tdomain%sFace(nface)%Veloc_p(i,0:1)  = Tdomain%specel(nelem)%Veloc (ngll-1-i,ngllz-1,0:1)
                  Tdomain%sFace(nface)%Strain_p(i,0:2) = Tdomain%specel(nelem)%Strain(ngll-1-i,ngllz-1,0:2)
              end do
          else
              do i=0,ngll-1
                  Tdomain%sFace(nface)%Veloc_p(i,0:1)  = Tdomain%specel(nelem)%Veloc (0,ngll-1-i,0:1)
                  Tdomain%sFace(nface)%Strain_p(i,0:2) = Tdomain%specel(nelem)%Strain(0,ngll-1-i,0:2)
              end do
          endif
      end if
  enddo
  return

end subroutine Get_data_el2f


!>
!!\brief This Subroutine, for an element nelem, sends the numerical fluxes
!! computed on its neighboring faces to update the forces of the current element.
!! It takes into account the problems of coherency and weighting with the
!! right wheight for approximating the integral.
!!\version 1.0
!!\date 20/11/2013
!! This subroutine is used only with DG elements
!<
subroutine Get_flux_f2el (Tdomain, nelem)

  use sdomain
  implicit none

  type (domain), intent (INOUT) :: Tdomain
  integer, intent(in)   :: nelem

  ! local variables
  real, dimension (:,:), allocatable :: Flux
  integer :: ngll, ngllx, ngllz, i, imin, imax, nface, nf
  type(element), pointer :: Elem
  logical :: coherency

  Elem => Tdomain%specel(nelem)
  ngllx = Elem%ngllx
  ngllz = Elem%ngllz

  do nf = 0,3
      nface = Elem%Near_Face(nf)
      coherency = Tdomain%sFace(nface)%coherency
      ngll  = Tdomain%sFace(nface)%ngll
      allocate(Flux(0:ngll-1,0:4))
      call get_iminimax(Elem,nf,imin,imax)

      if (Tdomain%sFace(nface)%Near_Element(0)==nelem) then
          do i=0,ngll-1
              Flux(i,:) =  Tdomain%sFace(nface)%Flux(i,:)*Elem%Coeff_Integr_Faces(imin+i)
          enddo
      else ! nelem is Near_Element(1)
          if (coherency .AND. (.NOT. Tdomain%sFace(nface)%changing_media)) then
              do i=0,ngll-1
                  Flux(i,:) = -Tdomain%sFace(nface)%Flux(i,:)*Elem%Coeff_Integr_Faces(imin+i)
              enddo
          else if ((.NOT. coherency) .AND. (.NOT. Tdomain%sFace(nface)%changing_media)) then
              do i=0,ngll-1
                  Flux(i,:) = -Tdomain%sFace(nface)%Flux(ngll-1-i,:)*Elem%Coeff_Integr_Faces(imin+i)
              enddo
          else if (coherency .AND. Tdomain%sFace(nface)%changing_media) then
              do i=0,ngll-1
                  Flux(i,:) = Tdomain%sFace(nface)%Flux_p(i,:)*Elem%Coeff_Integr_Faces(imin+i)
              enddo
          else if ((.NOT. coherency) .AND. Tdomain%sFace(nface)%changing_media) then
              do i=0,ngll-1
                  Flux(i,:) = Tdomain%sFace(nface)%Flux_p(ngll-1-i,:)*Elem%Coeff_Integr_Faces(imin+i)
              enddo
          endif
      endif

      select case(nf)
      case(0) ! Bottom Face
          do i=0,ngll-1
              Tdomain%specel(nelem)%Forces(i,0,:) = Tdomain%specel(nelem)%Forces(i,0,:) - Flux(i,:)
          enddo
      case(1) ! Right Face
          do i=0,ngll-1
              Tdomain%specel(nelem)%Forces(ngllx-1,i,:) = Tdomain%specel(nelem)%Forces(ngllx-1,i,:) - Flux(i,:)
          enddo
      case(2) ! Top Face
          do i=0,ngll-1
              Tdomain%specel(nelem)%Forces(i,ngllz-1,:) = Tdomain%specel(nelem)%Forces(i,ngllz-1,:) - Flux(i,:)
          enddo
      case(3) ! Left Face
          do i=0,ngll-1
              Tdomain%specel(nelem)%Forces(0,i,:) = Tdomain%specel(nelem)%Forces(0,i,:) - Flux(i,:)
          enddo
      end select
  deallocate (Flux)
  enddo

  ! Dirichlet Boundary conditions (if any)
  if(Tdomain%type_bc==DG_BC_REFL) call enforce_diriclet_BC(Tdomain,nelem)

end subroutine Get_flux_f2el


!>
!!\brief Subroutine that shares the traction from an element
!! "nelem" to its face "nface" taking into account the problems of coherency.
!!\version 1.0
!!\date 18/11/2013
!! This subroutine is used only with HDG elements
!<
subroutine Get_traction_el2f (Tdomain, nelem)

    use sdomain
    use ssources
    use constants
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(in)   :: nelem

    ! local variables
    integer :: ngll, ngllx, ngllz, i, imin, imax, nface, nf
    logical :: coherency

    ngllx = Tdomain%specel(nelem)%ngllx
    ngllz = Tdomain%specel(nelem)%ngllz

    do nf = 0,3
        nface = Tdomain%specel(nelem)%Near_Face(nf)
        ngll  = Tdomain%sFace(nface)%ngll
        coherency  = Tdomain%sFace(nface)%coherency
        call get_iminimax(Tdomain%specel(nelem),nf,imin,imax)
        if (coherency .OR. (Tdomain%sFace(nface)%Near_Element(0)==nelem)) then
            Tdomain%sFace(nface)%Traction = Tdomain%sFace(nface)%Traction &
                                          + Tdomain%specel(nelem)%TracFace(imin:imax,0:1)
        else ! Case coherency = false
            do i=0,ngll-1
                Tdomain%sFace(nface)%Traction(i,0:1) = Tdomain%sFace(nface)%Traction(i,0:1) &
                                                     + Tdomain%specel(nelem)%TracFace(imax-i,0:1)
            end do
        end if
    enddo

    ! Ajout des sources surfaciques pour les problemes de Riemann autour des elements sources
    !if (Tdomain%specel(nelem)%is_source) then
    !    allocate(Fext(0:ngll-1,0:1))
    !    select case(w_face)
    !    case(0)
    !        Fext(0:ngll-1,0) =  CompSource(Tdomain%sSource(0),timelocal) &
    !                          * Tdomain%sSource(0)%Elem(0)%ExtForce(0:ngllx-1,0,0)
    !        Fext(0:ngll-1,1) =  CompSource(Tdomain%sSource(0),timelocal) &
    !                          * Tdomain%sSource(0)%Elem(0)%ExtForce(0:ngllx-1,0,1)
    !    case(1)
    !        Fext(0:ngll-1,0) =  CompSource(Tdomain%sSource(0),timelocal) &
    !                          * Tdomain%sSource(0)%Elem(0)%ExtForce(ngllx-1,0:ngllz-1,0)
    !        Fext(0:ngll-1,1) =  CompSource(Tdomain%sSource(0),timelocal) &
    !                          * Tdomain%sSource(0)%Elem(0)%ExtForce(ngllx-1,0:ngllz-1,1)
    !    case(2)
    !        Fext(0:ngll-1,0) =  CompSource(Tdomain%sSource(0),timelocal) &
    !                          * Tdomain%sSource(0)%Elem(0)%ExtForce(0:ngllx-1,ngllz-1,0)
    !        Fext(0:ngll-1,1) =  CompSource(Tdomain%sSource(0),timelocal) &
    !                          * Tdomain%sSource(0)%Elem(0)%ExtForce(0:ngllx-1,ngllz-1,1)
    !    case(3)
    !        Fext(0:ngll-1,0) =  CompSource(Tdomain%sSource(0),timelocal) &
    !                          * Tdomain%sSource(0)%Elem(0)%ExtForce(0,0:ngllz-1,0)
    !        Fext(0:ngll-1,1) =  CompSource(Tdomain%sSource(0),timelocal) &
    !                          * Tdomain%sSource(0)%Elem(0)%ExtForce(0,0:ngllz-1,1)
    !    end select
    !    Tdomain%sFace(nface)%Traction = Tdomain%sFace(nface)%Traction - 0.005 * Fext
    !endif
    return

end subroutine Get_traction_el2f

!>
!!\brief This Subroutine gets the Forces and Velocities from 2 neighboring vertices
!! This subroutine is used in a framework of coupling CG with HDG, and only
!! called for faces at the interface.
!!\version 1.0
!!\date 17/12/2014
!! This subroutine is used only with HDG elements
!<
subroutine get_veloc_v2f (Tdomain, nface)

    use sdomain
    use constants
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(IN)   :: nface
    integer               :: nv, ngll

    ngll = Tdomain%sFace(nface)%ngll

    ! Assemblage des Veloc et des Forces sur la Face
    nv = Tdomain%sFace(nface)%Near_Vertex(0)
    Tdomain%sFace(nface)%Veloc (0,:) = Tdomain%sVertex(nv)%Veloc(:)
    nv = Tdomain%sFace(nface)%Near_Vertex(1)
    Tdomain%sFace(nface)%Veloc (ngll-1,:) = Tdomain%sVertex(nv)%Veloc(:)

end subroutine get_veloc_v2f

!>
!!\brief This Subroutine sends the Forces from a face to its 2 neighboring vertices
!! This subroutine is used in a framework of coupling CG with HDG, and only
!! called for faces at the interface.
!!\version 1.0
!!\date 17/12/2014
!! This subroutine is used only with HDG elements
!<
subroutine get_forces_f2v (Tdomain, nface)

    use sdomain
    use constants
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(IN)   :: nface
    integer               :: nv, ngll

    ngll = Tdomain%sFace(nface)%ngll

    ! Assemblage des Veloc et des Forces sur la Face
    nv = Tdomain%sFace(nface)%Near_Vertex(0)
    Tdomain%sVertex(nv)%Forces(:) = Tdomain%sVertex(nv)%Forces(:) + Tdomain%sFace(nface)%Forces(0,:)
    nv = Tdomain%sFace(nface)%Near_Vertex(1)
    Tdomain%sVertex(nv)%Forces(:) = Tdomain%sVertex(nv)%Forces(:) + Tdomain%sFace(nface)%Forces(ngll-1,:)

end subroutine get_forces_f2v

!>
!!\brief Subroutine that sent the traces of velocities "Vhat" computed on a face
!! "nface" to the neighboring element "nelem" dealing with the problems of coherency.
!! Be carreful : here the wheight for integrals are not taked into account.
!!\version 1.0
!!\date 03/04/2014
!! This subroutine is used only with HDG elements
!<
subroutine Get_Vhat_f2el (Tdomain, nelem)

    use sdomain
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(in)   :: nelem
    integer               :: nface, w_face

    ! local variables
    integer :: ngll, ngllx, ngllz, i, imin, imax
    logical :: coherency

    ngllx = Tdomain%specel(nelem)%ngllx
    ngllz = Tdomain%specel(nelem)%ngllz

    do w_face=0,3
        nface = Tdomain%specel(nelem)%Near_Face(w_face)
        coherency = Tdomain%sFace(nface)%coherency
        ngll  = Tdomain%sFace(nface)%ngll
        call get_iminimax(Tdomain%specel(nelem),w_face,imin,imax)
        if (coherency .OR. (Tdomain%sFace(nface)%Near_Element(0)==nelem)) then
            Tdomain%specel(nelem)%Vhat(imin:imax,:)  = Tdomain%sFace(nface)%Veloc(:,:)
        else
            do i=0,ngll-1
                Tdomain%specel(nelem)%Vhat(imax-i,:) = Tdomain%sFace(nface)%Veloc(i,:)
            enddo
        endif
    enddo

end subroutine Get_Vhat_f2el


!>
!!\brief Subroutine in which the element elem sends its contribution to the
!! second member R (= tractions)of the system K*Lambda = R (for lagrange multiplicators)
!! to its neighbouring faces and vertices.
!!\version 1.0
!!\date 03/10/2014
!! This subroutine is used only with HDG elements in a semi-implicit framework
!<
subroutine Get_R_el2fv (Tdomain, nelem)

    use sdomain
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(IN)    :: nelem
    integer                :: imin, imax, i, nv, w_face, nface, ngll, pos1, pos2, n1, n2
    type(element), pointer :: Elem
    logical :: coherency
    Elem => Tdomain%specel(nelem)

    ! Send R from current element to its neghbouring faces
    do w_face=0,3
        call get_iminimax(Tdomain%specel(nelem),w_face,imin,imax)
        nface = Elem%Near_Face(w_face)
        ngll  = Tdomain%sFace(nface)%ngll
        coherency = Tdomain%sFace(nface)%coherency
        if (coherency .OR. (Tdomain%sFace(nface)%Near_Element(0)==nelem)) then
            Tdomain%sFace(nface)%smbr(:,0:1) = Tdomain%sFace(nface)%smbr(:,0:1) + Elem%TracFace(imin:imax,0:1)
        else
            do i=0,ngll-1
                Tdomain%sFace(nface)%smbr(i,0:1) = Tdomain%sFace(nface)%smbr(i,0:1) + Elem%TracFace(imax-i,0:1)
            enddo
        endif
    enddo

    ! Send R from current element to its neghbouring vertices
    do i=0,3
        nv = Elem%Near_Vertex(i)
        call get_gll_arround_corner(Elem,i,n1,n2)
        pos1 = Elem%pos_corner_in_VertMat(i,0)
        pos2 = Elem%pos_corner_in_VertMat(i,1)
        Tdomain%sVertex(nv)%smbrLambda(pos1:pos1+1) = Tdomain%sVertex(nv)%smbrLambda(pos1:pos1+1) &
                                                    + Elem%TracFace(n1,0:1)
        Tdomain%sVertex(nv)%smbrLambda(pos2:pos2+1) = Tdomain%sVertex(nv)%smbrLambda(pos2:pos2+1) &
                                                    + Elem%TracFace(n2,0:1)
    enddo

end subroutine Get_R_el2fv


!>
!!\brief Subroutine in which the face gets the velocities lambda (lagrange multiplicators)
!! from the vertices to the Face's ends.
!!\version 1.0
!!\date 03/10/2014
!! This subroutine is used only with HDG elements in a semi-implicit framework
!<
subroutine Get_lambda_v2f (Tdomain, nface)

    use sdomain
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(in)   :: nface
    integer               :: i, ngll, nv, pos

    ngll = Tdomain%sFace(nface)%ngll

    do i=0,1
        nv = Tdomain%sFace(nface)%Near_Vertex(i)
        pos = Tdomain%sFace(nface)%pos_in_VertMat(i)
        Tdomain%sFace(nface)%Veloc(i*(ngll-1),:) = Tdomain%sVertex(nv)%Lambda(pos:pos+1)
    enddo

end subroutine Get_lambda_v2f


!>
!!\brief This Subroutine enforces the Dirichlet Boundary Conditions (Veloc = 0)
!! in a DG or HDG context, setting the corresponding forces to zero.
!! Subroutine particularly usefull in a context using PML.
!!\version 1.0
!!\date 20/06/2014
!! This subroutine is used only with DG and HDG elements
!<
subroutine enforce_diriclet_BC (Tdomain, nelem)

    use sdomain
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(in)   :: nelem
    integer :: i, nface, ngllx, ngllz
    logical :: is_refl

    ngllx = Tdomain%specel(nelem)%ngllx
    ngllz = Tdomain%specel(nelem)%ngllz

    do i=0,3
        nface = Tdomain%specel(nelem)%Near_Face(i)
        is_refl = Tdomain%sFace(nface)%reflex
        if (is_refl) then
            select case (i)
            case(0)
                Tdomain%specel(nelem)%Forces(0:ngllx-1,0,3:4) = 0.
            case(1)
                Tdomain%specel(nelem)%Forces(ngllx-1,0:ngllz-1,3:4) = 0.
            case(2)
                Tdomain%specel(nelem)%Forces(0:ngllx-1,ngllz-1,3:4) = 0.
            case(3)
                Tdomain%specel(nelem)%Forces(0,0:ngllz-1,3:4) = 0.
            end select
        endif
    enddo

end subroutine enforce_diriclet_BC

!>
!!\brief This Subroutine enforces the Dirichlet Boundary Conditions (Veloc = 0)
!! in a DG or HDG context, setting the corresponding forces to zero.
!! Subroutine particularly usefull in a context using PML.
!!\version 1.0
!!\date 20/06/2014
!! This subroutine is used only with DG and HDG elements
!<
subroutine enforce_diriclet_corners_vhat (Tdomain, nelem)

    use sdomain
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(in)   :: nelem
    integer :: i, nface, ngllx, ngllz, nglltot, index1, index2
    logical :: is_refl

    ngllx = Tdomain%specel(nelem)%ngllx
    ngllz = Tdomain%specel(nelem)%ngllz
    nglltot = 2*(ngllx + ngllz)

    do i=0,3
        nface = Tdomain%specel(nelem)%Near_Face(i)
        is_refl = Tdomain%sFace(nface)%reflex
        if (is_refl) then
            select case (i)
                case(0)
                    index1 = 2*ngllx + ngllz
                    index2 = ngllx
                case(1)
                    index1 = ngllx-1
                    index2 = 2*ngllx + ngllz-1
                case(2)
                    index1 = 2*ngllx + 2*ngllz -1
                    index2 = ngllx + ngllz-1
                case(3)
                    index1 = 0
                    index2 = ngllx + ngllz
                end select
            Tdomain%specel(nelem)%Vhat(index1,:) = 0.
            Tdomain%specel(nelem)%Vhat(index2,:) = 0.
        endif
    enddo

end subroutine enforce_diriclet_corners_vhat

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
