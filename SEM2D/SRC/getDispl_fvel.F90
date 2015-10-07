!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file getDispl_fvel.F90
!!\brief Contient la subroutine get_Displ_fv2el.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) n
!<
subroutine get_Displ_fv2el(Tdomain,n)
    use sdomain
    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n

    type(element), pointer :: el
    type(face), pointer :: fc
    integer :: nx,nz,nv

    el => Tdomain%specel(n)
    nx = el%ngllx;  nz = el%ngllz

    fc => Tdomain%sFace(el%Near_Face(0))
    if ( fc%near_element(0) == n .or. fc%coherency) then
        el%Forces(1:nx-2,0,0:1) =  fc%Forces(:,0:1)
    else
        el%Forces(1:nx-2,0,0:1) =  fc%Forces(nx-2:1:-1,0:1)
    endif

    fc => Tdomain%sFace(el%Near_Face(1))
    if ( fc%near_element(0) == n.or. fc%coherency) then
        el%Forces(nx-1,1:nz-2,0:1) =  fc%Forces(:,0:1)
    else
        el%Forces(nx-1,1:nz-2,0:1) =  fc%Forces(nz-2:1:-1,0:1)
    endif

    fc => Tdomain%sFace(el%Near_Face(2))
    if ( fc%near_element(0) == n.or. fc%coherency) then
        el%Forces(1:nx-2,nz-1,0:1) =  fc%Forces(:,0:1)
    else
        el%Forces(1:nx-2,nz-1,0:1) =  fc%Forces(nx-2:1:-1,0:1)
    endif

    fc => Tdomain%sFace(el%Near_Face(3))
    if ( fc%near_element(0) == n.or. fc%coherency) then
        el%Forces(0,1:nz-2,0:1) =  fc%Forces(:,0:1)
    else
        el%Forces(0,1:nz-2,0:1) =  fc%Forces(nz-2:1:-1,0:1)
    endif

    nv = el%Near_Vertex(0)
    el%Forces (0,0,0:1) =Tdomain%sVertex(nv)%Forces(0:1)
    nv = el%Near_Vertex(1)
    el%Forces (nx-1,0,0:1) =Tdomain%sVertex(nv)%Forces(0:1)
    nv = el%Near_Vertex(2)
    el%Forces (nx-1,nz-1,0:1) =Tdomain%sVertex(nv)%Forces(0:1)
    nv = el%Near_Vertex(3)
    el%Forces (0,nz-1,0:1) =Tdomain%sVertex(nv)%Forces(0:1)

    return
end subroutine get_Displ_fv2el

!>
!! \brief This subroutine works like the precedent subroutine get_Displ_fv2el
!! but sending the real array Face%Displ and Vertex%Displ to the Element n.
!! This subroutine has been implemented because get_Displ_fv2el is confusing...
!! Indeed, get_Displ_fv2el assumes that Element%Forces, Face%Forces, and Vertex%Forces
!! are actually displacements, and not forces, which is true in the Newmark framework,
!! but false in others time schemes (such as Runge-Kutta).
!!
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) n
!<
subroutine get_RealDispl_fv2el(Tdomain,n)
    use sdomain
    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n

    type(element), pointer :: el
    type(face), pointer :: fc
    integer :: nx,nz,nv

    el => Tdomain%specel(n)
    nx = el%ngllx;  nz = el%ngllz
    el%Forces(1:nx-2,1:nz-2,0:1) = el%Displ(:,:,0:1)

    fc => Tdomain%sFace(el%Near_Face(0))
    if ( fc%near_element(0) == n .or. fc%coherency) then
        el%Forces(1:nx-2,0,0:1) =  fc%Displ(:,0:1)
    else
        el%Forces(1:nx-2,0,0:1) =  fc%Displ(nx-2:1:-1,0:1)
    endif

    fc => Tdomain%sFace(el%Near_Face(1))
    if ( fc%near_element(0) == n.or. fc%coherency) then
        el%Forces(nx-1,1:nz-2,0:1) =  fc%Displ(:,0:1)
    else
        el%Forces(nx-1,1:nz-2,0:1) =  fc%Displ(nz-2:1:-1,0:1)
    endif

    fc => Tdomain%sFace(el%Near_Face(2))
    if ( fc%near_element(0) == n.or. fc%coherency) then
        el%Forces(1:nx-2,nz-1,0:1) =  fc%Displ(:,0:1)
    else
        el%Forces(1:nx-2,nz-1,0:1) =  fc%Displ(nx-2:1:-1,0:1)
    endif

    fc => Tdomain%sFace(el%Near_Face(3))
    if ( fc%near_element(0) == n.or. fc%coherency) then
        el%Forces(0,1:nz-2,0:1) =  fc%Displ(:,0:1)
    else
        el%Forces(0,1:nz-2,0:1) =  fc%Displ(nz-2:1:-1,0:1)
    endif

    nv = el%Near_Vertex(0)
    el%Forces(0,0,0:1) =Tdomain%sVertex(nv)%Displ(0:1)
    nv = el%Near_Vertex(1)
    el%Forces(nx-1,0,0:1) =Tdomain%sVertex(nv)%Displ(0:1)
    nv = el%Near_Vertex(2)
    el%Forces(nx-1,nz-1,0:1) =Tdomain%sVertex(nv)%Displ(0:1)
    nv = el%Near_Vertex(3)
    el%Forces(0,nz-1,0:1) =Tdomain%sVertex(nv)%Displ(0:1)

    return
end subroutine get_RealDispl_fv2el

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
