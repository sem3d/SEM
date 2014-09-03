!>
!!\file Assemblage.F90
!!\brief Contains subroutines for assembling forces and distributing
!! forces and displacements over the neighboring elements of a face
!!\version 1.0
!!\date 15/11/2013
!! This algorithm come from a part of the previous Newmark routine
!! with several modifications.
!<
subroutine Assemblage (Tdomain, nelem, nface, w_face)

  use sdomain

  implicit none
  type (domain), intent (INOUT) :: Tdomain
  integer, intent(in)   :: nelem
  integer, intent(in)   :: nface
  integer, intent(in)   :: w_face
  logical               :: coherency
  integer               :: vertex0, vertex1, ngllx, ngllz

  ngllx =  Tdomain%specel(nelem)%ngllx
  ngllz =  Tdomain%specel(nelem)%ngllz
  coherency = Tdomain%sFace(nface)%coherency

  !if (.not. Tdomain%sFace(nface)%is_visited) then
  !   Tdomain%sFace(nface)%Forces = 0.
  !   Tdomain%sFace(nface)%is_visited = .true.
  !endif

  ! Assemblage of the forces on the Face nface coming from element nelem
  if (nelem==Tdomain%sFace(nface)%Near_Element(0)) then
     call getInternalF_el2f (Tdomain,nelem,nface,w_face,.true.)
  else
     call getInternalF_el2f (Tdomain,nelem,nface,w_face,coherency)
  endif
  ! Assemblage of the forces of the two vertexes of nface
  vertex0 = Tdomain%specel(nelem)%Near_Vertex(w_face)
  vertex1 = Tdomain%specel(nelem)%Near_Vertex(modulo(w_face+1,4))

  if (.not. Tdomain%sVertex(vertex0)%is_computed) then
     Tdomain%sVertex(vertex0)%Forces = 0.
     Tdomain%sVertex(vertex0)%is_computed = .true.
  endif
  if (.not. Tdomain%sVertex(vertex1)%is_computed) then
     Tdomain%sVertex(vertex1)%Forces = 0.
     Tdomain%sVertex(vertex1)%is_computed = .true.
  endif

  select case (w_face)
  case(0)
     Tdomain%sVertex(vertex0)%Forces = Tdomain%sVertex(vertex0)%Forces + Tdomain%specel(nelem)%Forces(0,0,0:1)
     Tdomain%sVertex(vertex1)%Forces = Tdomain%sVertex(vertex1)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,0,0:1)
  case(1)
     Tdomain%sVertex(vertex0)%Forces = Tdomain%sVertex(vertex0)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,0,0:1)
     Tdomain%sVertex(vertex1)%Forces = Tdomain%sVertex(vertex1)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,ngllz-1,0:1)
  case(2)
     Tdomain%sVertex(vertex0)%Forces = Tdomain%sVertex(vertex0)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,ngllz-1,0:1)
     Tdomain%sVertex(vertex1)%Forces = Tdomain%sVertex(vertex1)%Forces + Tdomain%specel(nelem)%Forces(0,ngllz-1,0:1)
  case(3)
     Tdomain%sVertex(vertex0)%Forces = Tdomain%sVertex(vertex0)%Forces + Tdomain%specel(nelem)%Forces(0,ngllz-1,0:1)
     Tdomain%sVertex(vertex1)%Forces = Tdomain%sVertex(vertex1)%Forces + Tdomain%specel(nelem)%Forces(0,0,0:1)
  end select

end subroutine Assemblage


!>
!!\brief Subroutine that shares the data (velocities and strains) from an element
!! "nelem" to its face "nface" taking into account the problems of coherency.
!!\version 1.0
!!\date 18/11/2013
!! This subroutine is used only with DG elements
!<
subroutine Get_data_el2f (Tdomain, nelem, nface, w_face)

  use sdomain
  implicit none

  type (domain), intent (INOUT) :: Tdomain
  integer, intent(in)   :: nelem
  integer, intent(in)   :: nface
  integer, intent(in)   :: w_face

  ! local variables
  integer :: ngll, ngllx, ngllz, i
  logical :: coherency

  ngll  = Tdomain%sFace(nface)%ngll
  ngllx = Tdomain%specel(nelem)%ngllx
  ngllz = Tdomain%specel(nelem)%ngllz
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
  return

end subroutine Get_data_el2f


!>
!!\brief Subroutine that sent the fluxes computed on a face "nface" to the neighboring
!! element "nelem" taking into account the problems of coherency and weighting with the
!! right wheight for approximating the integral.
!!\version 1.0
!!\date 20/11/2013
!! This subroutine is used only with DG elements
!<
subroutine Get_flux_f2el (Tdomain, nelem, nface, w_face)

  use sdomain
  implicit none

  type (domain), intent (INOUT) :: Tdomain
  integer, intent(in)   :: nelem
  integer, intent(in)   :: nface
  integer, intent(in)   :: w_face

  ! local variables
  integer :: ngll, ngllx, ngllz, i, imin, imax
  logical :: coherency

  ngll  = Tdomain%sFace(nface)%ngll
  ngllx = Tdomain%specel(nelem)%ngllx
  ngllz = Tdomain%specel(nelem)%ngllz
  coherency = Tdomain%sFace(nface)%coherency

  call get_iminimax(Tdomain%specel(nelem),w_face,imin,imax)

  if (coherency .OR. (Tdomain%sFace(nface)%Near_Element(0)==nelem)) then
     if (w_face == 0 ) then
        do i=0,ngll-1
           Tdomain%specel(nelem)%Forces(i,0,0:4) = Tdomain%specel(nelem)%Forces(i,0,0:4) - &
               Tdomain%sFace(nface)%Flux(i,0:4) * Tdomain%specel(nelem)%Coeff_Integr_Faces(imin+i)
        enddo
     else if (w_face == 1 ) then
         do i=0,ngll-1
             Tdomain%specel(nelem)%Forces(ngllx-1,i,0:4) = Tdomain%specel(nelem)%Forces(ngllx-1,i,0:4) - &
               Tdomain%sFace(nface)%Flux(i,0:4) * Tdomain%specel(nelem)%Coeff_Integr_Faces(imin+i)
         enddo
     else if (w_face == 2 ) then
         do i=0,ngll-1
             Tdomain%specel(nelem)%Forces(i,ngllz-1,0:4) = Tdomain%specel(nelem)%Forces(i,ngllz-1,0:4) - &
               Tdomain%sFace(nface)%Flux(i,0:4) * Tdomain%specel(nelem)%Coeff_Integr_Faces(imin+i)
         enddo
     else
         do i=0,ngll-1
             Tdomain%specel(nelem)%Forces(0,i,0:4) = Tdomain%specel(nelem)%Forces(0,i,0:4) - &
               Tdomain%sFace(nface)%Flux(i,0:4) * Tdomain%specel(nelem)%Coeff_Integr_Faces(imin+i)
         enddo
     endif
  else
     if (w_face == 0 ) then
        do i=0,ngll-1
           Tdomain%specel(nelem)%Forces(i,0,0:4) = Tdomain%specel(nelem)%Forces(i,0,0:4) - &
               Tdomain%sFace(nface)%Flux(ngll-1-i,0:4) * Tdomain%specel(nelem)%Coeff_Integr_Faces(imin+i)
        enddo
     else if (w_face == 1 ) then
        do i=0,ngll-1
           Tdomain%specel(nelem)%Forces(ngllx-1,i,0:4) = Tdomain%specel(nelem)%Forces(ngllx-1,i,0:4) - &
               Tdomain%sFace(nface)%Flux(ngll-1-i,0:4) * Tdomain%specel(nelem)%Coeff_Integr_Faces(imin+i)
        enddo
     else if (w_face == 2 ) then
        do i=0,ngll-1
           Tdomain%specel(nelem)%Forces(i,ngllz-1,0:4) = Tdomain%specel(nelem)%Forces(i,ngllz-1,0:4) - &
               Tdomain%sFace(nface)%Flux(ngll-1-i,0:4) * Tdomain%specel(nelem)%Coeff_Integr_Faces(imin+i)
        enddo
     else
        do i=0,ngll-1
           Tdomain%specel(nelem)%Forces(0,i,0:4) = Tdomain%specel(nelem)%Forces(0,i,0:4) - &
               Tdomain%sFace(nface)%Flux(ngll-1-i,0:4) * Tdomain%specel(nelem)%Coeff_Integr_Faces(imin+i)
        enddo
     endif
  endif


end subroutine Get_flux_f2el


!>
!!\brief Subroutine that shares the traction from an element
!! "nelem" to its face "nface" taking into account the problems of coherency.
!!\version 1.0
!!\date 18/11/2013
!! This subroutine is used only with HDG elements
!<
subroutine Get_traction_el2f (Tdomain, nelem, nface, w_face, timelocal)

    use sdomain
    use ssources
    use constants
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(in)   :: nelem
    integer, intent(in)   :: nface
    integer, intent(in)   :: w_face
    real,   intent (in)    :: timelocal

    ! local variables
    integer :: ngll, ngllx, ngllz, i, imin, imax, np
    logical :: coherency
    real, dimension(:,:), allocatable :: Fext

    ngll  = Tdomain%sFace(nface)%ngll
    ngllx = Tdomain%specel(nelem)%ngllx
    ngllz = Tdomain%specel(nelem)%ngllz
    coherency  = Tdomain%sFace(nface)%coherency

    call get_iminimax(Tdomain%specel(nelem),w_face,imin,imax)

    if (coherency .OR. (Tdomain%sFace(nface)%Near_Element(0)==nelem)) then
        Tdomain%sFace(nface)%Traction = Tdomain%sFace(nface)%Traction &
                                      + Tdomain%specel(nelem)%TracFace(imin:imax,0:1)
    else ! Case coherency = false
        do i=0,ngll-1
            Tdomain%sFace(nface)%Traction(i,0:1) = Tdomain%sFace(nface)%Traction(i,0:1) &
                                                 + Tdomain%specel(nelem)%TracFace(imax-i,0:1)
        end do
    end if

    ! Ajout des sources surfaciques pour les problemes de Riemann autour des elements sources
    if (Tdomain%specel(nelem)%is_source) then
        allocate(Fext(0:ngll-1,0:1))
        select case(w_face)
        case(0)
            Fext(0:ngll-1,0) =  CompSource(Tdomain%sSource(0),timelocal) &
                              * Tdomain%sSource(0)%Elem(0)%ExtForce(0:ngllx-1,0,0)
            Fext(0:ngll-1,1) =  CompSource(Tdomain%sSource(0),timelocal) &
                              * Tdomain%sSource(0)%Elem(0)%ExtForce(0:ngllx-1,0,1)
        case(1)
            Fext(0:ngll-1,0) =  CompSource(Tdomain%sSource(0),timelocal) &
                              * Tdomain%sSource(0)%Elem(0)%ExtForce(ngllx-1,0:ngllz-1,0)
            Fext(0:ngll-1,1) =  CompSource(Tdomain%sSource(0),timelocal) &
                              * Tdomain%sSource(0)%Elem(0)%ExtForce(ngllx-1,0:ngllz-1,1)
        case(2)
            Fext(0:ngll-1,0) =  CompSource(Tdomain%sSource(0),timelocal) &
                              * Tdomain%sSource(0)%Elem(0)%ExtForce(0:ngllx-1,ngllz-1,0)
            Fext(0:ngll-1,1) =  CompSource(Tdomain%sSource(0),timelocal) &
                              * Tdomain%sSource(0)%Elem(0)%ExtForce(0:ngllx-1,ngllz-1,1)
        case(3)
            Fext(0:ngll-1,0) =  CompSource(Tdomain%sSource(0),timelocal) &
                              * Tdomain%sSource(0)%Elem(0)%ExtForce(0,0:ngllz-1,0)
            Fext(0:ngll-1,1) =  CompSource(Tdomain%sSource(0),timelocal) &
                              * Tdomain%sSource(0)%Elem(0)%ExtForce(0,0:ngllz-1,1)
        end select
        Tdomain%sFace(nface)%Traction = Tdomain%sFace(nface)%Traction - 0.005 * Fext
    endif
    return

end subroutine Get_traction_el2f

!>
!!\brief Subroutine that sent the traces of velocities "Vhat" computed on a face
!! "nface" to the neighboring element "nelem" dealing with the problems of coherency.
!! Be carreful : here the wheight for integrals are not taked into account.
!!\version 1.0
!!\date 03/04/2014
!! This subroutine is used only with HDG elements
!<
subroutine Get_Vhat_f2el (Tdomain, nelem, nface, w_face)

    use sdomain
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(in)   :: nelem
    integer, intent(in)   :: nface
    integer, intent(in)   :: w_face

    ! local variables
    integer :: ngll, ngllx, ngllz, i, imin, imax
    logical :: coherency

    ngll  = Tdomain%sFace(nface)%ngll
    ngllx = Tdomain%specel(nelem)%ngllx
    ngllz = Tdomain%specel(nelem)%ngllz
    coherency = Tdomain%sFace(nface)%coherency

    call get_iminimax(Tdomain%specel(nelem),w_face,imin,imax)

    if (coherency .OR. (Tdomain%sFace(nface)%Near_Element(0)==nelem)) then
        Tdomain%specel(nelem)%Vhat(imin:imax,:)  = Tdomain%sFace(nface)%Veloc(:,:)
    else
        do i=0,ngll-1
            Tdomain%specel(nelem)%Vhat(imax-i,:) = Tdomain%sFace(nface)%Veloc(i,:)
        enddo
    endif

end subroutine Get_Vhat_f2el


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
    integer :: i, nface, ngllx, ngllz, nglltot, imin, imax, index1, index2
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
