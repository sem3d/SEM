!>
!!\file Assemblage.F90
!!\brief Contains subroutines for assembling forces and distributing
!! forces and displacements over the neighboring elements of a face
!!\version 1.0
!!\date 15/11/2013
!! This algorithm come from a part of the previous Newmark routine 
!! with several modifications.
!<

subroutine Assemblage (Tdomain, nelem, nface, nf)

  use sdomain

  implicit none
  type (domain), intent (INOUT) :: Tdomain
  integer, intent(in)   :: nelem
  integer, intent(in)   :: nface
  logical               :: coherency
  integer               :: vertex0, vertex1, ngllx, ngllz

  ngllx =  Tdomain%specel(nelem)%ngllx
  ngllz =  Tdomain%specel(nelem)%ngllz
  coherency = Tdomain%sFace(nface)%coherency

  ! Assemblage of the forces on the Face nface coming from element nelem 
  if (nelem==Tdomain%sFace(nface)%NearElement(0)) then
     call getInternalF_el2f (Tdomain,nelem,nface,nf,.true.)
  else
     call getInternalF_el2f (Tdomain,nelem,nface,nf,coherency)
  endif
  ! Assemblage of the forces of the two vertexes of nface
  vertex0 = Tdomain%specel(nelem)%NearVertex(nf)
  vertex1 = Tdomain%specel(nelem)%NearVertex(modulo(nf+1,4))
  select case (nf)
  case(0)
     Tdomain%sVertex(vertex0)%Forces = Tdomain%sVertex(vertex0)%Forces + Tdomain%specel(nelem)%Forces(0,0)
     Tdomain%sVertex(vertex1)%Forces = Tdomain%sVertex(vertex1)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,0)
  case(1)
     Tdomain%sVertex(vertex0)%Forces = Tdomain%sVertex(vertex0)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,0)
     Tdomain%sVertex(vertex1)%Forces = Tdomain%sVertex(vertex1)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,ngllz-1)
  case(2)
     Tdomain%sVertex(vertex0)%Forces = Tdomain%sVertex(vertex0)%Forces + Tdomain%specel(nelem)%Forces(ngllx-1,ngllz-1)
     Tdomain%sVertex(vertex1)%Forces = Tdomain%sVertex(vertex1)%Forces + Tdomain%specel(nelem)%Forces(0,ngllz-1)
  case(3)
     Tdomain%sVertex(vertex0)%Forces = Tdomain%sVertex(vertex0)%Forces + Tdomain%specel(nelem)%Forces(0,ngllz-1)
     Tdomain%sVertex(vertex1)%Forces = Tdomain%sVertex(vertex1)%Forces + Tdomain%specel(nelem)%Forces(0,0)
  end select
  
end subroutine Assemblage


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!>
!!\brief Set all the forces at the vertexes to zero
!!\version 1.0
!!\date 15/11/2013
!<
subroutine set_VerAndElem_Forces_Zero(Tdomain)

  use sdomain
  
  implicit none
  type (domain), intent (INOUT) :: Tdomain
  integer      :: n, ngll

  do n=0,Tdomain%n_vertex-1
     Tdomain%sVertex(n)%Forces = 0. 
  enddo
  do n=0,Tdomain%n_face-1
     ngll = Tdomain%sFace(n)%ngll
     Tdomain%sFace(n)%Forces(1:ngll-2,0:1) = 0.
  enddo
end subroutine set_VerAndElem_Forces_Zero
