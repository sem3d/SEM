!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file getif_el2f.F90
!!\brief Asure la sommation des forces pour évaluer la force exercée sur les faces
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) n_elem
!! \param integer, intent (IN) n_face
!! \param integer, intent (IN) w_face
!! \param logical logic
!<


subroutine getInternalF_el2f (Tdomain, n_elem, n_face, w_face, logic)
    use sdomain
    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n_elem, n_face, w_face
    logical :: logic

    ! local variables
    integer ::  ngll, ngllx, ngllz,i

    ! Modified by Gaetano Festa 01/06/2005

    ngll = Tdomain%sFace(n_face)%ngll
    ngllx = Tdomain%specel(n_elem)%ngllx
    ngllz = Tdomain%specel(n_elem)%ngllz
    if (logic) then
        if (w_face == 0 ) then
            Tdomain%sFace(n_face)%Forces (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces (1:ngll-2, 0,0:1)
        else if (w_face == 1 ) then
            Tdomain%sFace(n_face)%Forces (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces (ngllx-1,1:ngll-2,0:1)
        else if (w_face == 2 ) then
            Tdomain%sFace(n_face)%Forces (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces (1:ngll-2, ngllz-1,0:1)
        else
            Tdomain%sFace(n_face)%Forces (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces (0,1:ngll-2,0:1)
        endif
    else
        if (w_face == 0 ) then
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%Forces (i,0:1) = Tdomain%sFace(n_face)%Forces (i,0:1) + &
                    Tdomain%specel(n_elem)%Forces (ngll-1-i, 0,0:1)
            enddo
        else if (w_face == 1 ) then
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%Forces (i,0:1) = Tdomain%sFace(n_face)%Forces (i,0:1) + &
                    Tdomain%specel(n_elem)%Forces (ngllx-1,ngll-1-i,0:1)
            enddo
        else if (w_face == 2 ) then
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%Forces (i,0:1) = Tdomain%sFace(n_face)%Forces (i,0:1) + &
                    Tdomain%specel(n_elem)%Forces (ngll-1-i, ngllz-1,0:1)
            enddo
        else
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%Forces (i,0:1) = Tdomain%sFace(n_face)%Forces (i,0:1 ) + &
                    Tdomain%specel(n_elem)%Forces (0,ngll-1-i,0:1)
            enddo
        endif
    endif
    return
end subroutine  getInternalF_el2f

! ###################################################
!>
!! \brief
!!
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) n_elem
!! \param integer, intent (IN) n_face
!! \param integer, intent (IN) w_face
!! \param logical logic
!<


subroutine getInternalF_PML_el2f (Tdomain, n_elem, n_face, w_face, logic)
    use sdomain
    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n_elem, n_face, w_face
    logical :: logic

    ! local variables
    integer :: ngll, ngllx, ngllz,i

    ! Modified by Gaetano Festa 8/11/2004

    ngll = Tdomain%sFace(n_face)%ngll
    ngllx = Tdomain%specel(n_elem)%ngllx
    ngllz = Tdomain%specel(n_elem)%ngllz
    if (logic) then
        if (w_face == 0 ) then
            Tdomain%sFace(n_face)%Forces1 (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces1 (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces1 (1:ngll-2, 0,0:1)
            Tdomain%sFace(n_face)%Forces2 (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces2 (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces2 (1:ngll-2, 0,0:1)
        else if (w_face == 1 ) then
            Tdomain%sFace(n_face)%Forces1 (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces1 (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces1 (ngllx-1,1:ngll-2,0:1)
            Tdomain%sFace(n_face)%Forces2 (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces2 (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces2 (ngllx-1,1:ngll-2,0:1)
        else if (w_face == 2 ) then
            Tdomain%sFace(n_face)%Forces1 (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces1(1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces1 (1:ngll-2, ngllz-1,0:1)
            Tdomain%sFace(n_face)%Forces2 (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces2 (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces2 (1:ngll-2, ngllz-1,0:1)
        else
            Tdomain%sFace(n_face)%Forces1 (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces1 (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces1 (0,1:ngll-2,0:1)
            Tdomain%sFace(n_face)%Forces2 (1:ngll-2,0:1) = Tdomain%sFace(n_face)%Forces2 (1:ngll-2,0:1) + &
                Tdomain%specel(n_elem)%Forces2 (0,1:ngll-2,0:1)
        endif
    else
        if (w_face == 0 ) then
            do i = 0,ngll-1
                Tdomain%sFace(n_face)%Forces1 (i,0:1) = Tdomain%sFace(n_face)%Forces1 (i,0:1) + &
                    Tdomain%specel(n_elem)%Forces1 (ngll-1-i, 0,0:1)
                Tdomain%sFace(n_face)%Forces2 (i,0:1) = Tdomain%sFace(n_face)%Forces2 (i,0:1) + &
                    Tdomain%specel(n_elem)%Forces2 (ngll-1-i, 0,0:1)
            enddo
        else if (w_face == 1 ) then
            do i = 0,ngll-1
                Tdomain%sFace(n_face)%Forces1 (i,0:1) = Tdomain%sFace(n_face)%Forces1 (i,0:1) + &
                    Tdomain%specel(n_elem)%Forces1 (ngllx-1,ngll-1-i,0:1)
                Tdomain%sFace(n_face)%Forces2 (i,0:1) = Tdomain%sFace(n_face)%Forces2 (i,0:1) + &
                    Tdomain%specel(n_elem)%Forces2 (ngllx-1,ngll-1-i,0:1)
            enddo
        else if (w_face == 2 ) then
            do i = 0,ngll-1
                Tdomain%sFace(n_face)%Forces1 (i,0:1) = Tdomain%sFace(n_face)%Forces1 (i,0:1) + &
                    Tdomain%specel(n_elem)%Forces1 (ngll-1-i, ngllz-1,0:1)
                Tdomain%sFace(n_face)%Forces2 (i,0:1) = Tdomain%sFace(n_face)%Forces2 (i,0:1) + &
                    Tdomain%specel(n_elem)%Forces2 (ngll-1-i, ngllz-1,0:1)
            enddo
        else
            do i = 0,ngll-1
                Tdomain%sFace(n_face)%Forces1 (i,0:1) = Tdomain%sFace(n_face)%Forces1 (i,0:1 ) + &
                    Tdomain%specel(n_elem)%Forces1 (0,ngll-1-i,0:1)
                Tdomain%sFace(n_face)%Forces2 (i,0:1) = Tdomain%sFace(n_face)%Forces2 (i,0:1 ) + &
                    Tdomain%specel(n_elem)%Forces2 (0,ngll-1-i,0:1)
            enddo
        endif
    endif
    return
end subroutine  getInternalF_PML_el2f

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
