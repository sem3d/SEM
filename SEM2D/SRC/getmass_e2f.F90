!>
!!\file getmass_e2f.F90
!!\brief Contient la subroutine getMass_Element2Face.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Assure le calcul des masses associées à chaque face.
!!
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) n_elem
!! \param integer, intent (IN) n_face
!! \param integer, intent (IN) w_face
!! \param logical logic
!<


subroutine getMass_Element2Face (Tdomain, n_elem, n_face, w_face, logic)
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
            Tdomain%sFace(n_face)%MassMat (1:ngll-2) = Tdomain%sFace(n_face)%MassMat (1:ngll-2) + &
                Tdomain%specel(n_elem)%MassMat (1:ngll-2, 0)
        else if (w_face == 1 ) then
            Tdomain%sFace(n_face)%MassMat (1:ngll-2) = Tdomain%sFace(n_face)%MassMat (1:ngll-2) + &
                Tdomain%specel(n_elem)%MassMat (ngllx-1,1:ngll-2)
        else if (w_face == 2 ) then
            Tdomain%sFace(n_face)%MassMat (1:ngll-2) = Tdomain%sFace(n_face)%MassMat (1:ngll-2) + &
                Tdomain%specel(n_elem)%MassMat (1:ngll-2, ngllz-1)
        else
            Tdomain%sFace(n_face)%MassMat (1:ngll-2) = Tdomain%sFace(n_face)%MassMat (1:ngll-2) + &
                Tdomain%specel(n_elem)%MassMat (0,1:ngll-2)
        endif
    else
        if (w_face == 0 ) then
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%MassMat (i) = Tdomain%sFace(n_face)%MassMat (i) + &
                    Tdomain%specel(n_elem)%MassMat (ngll-1-i, 0)
            enddo
        else if (w_face == 1 ) then
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%MassMat (i) = Tdomain%sFace(n_face)%MassMat (i) + &
                    Tdomain%specel(n_elem)%MassMat (ngllx-1,ngll-1-i)
            enddo
        else if (w_face == 2 ) then
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%MassMat (i) = Tdomain%sFace(n_face)%MassMat (i) + &
                    Tdomain%specel(n_elem)%MassMat (ngll-1-i, ngllz-1)
            enddo
        else
            do i = 1,ngll-2
                Tdomain%sFace(n_face)%MassMat (i) = Tdomain%sFace(n_face)%MassMat (i ) + &
                    Tdomain%specel(n_elem)%MassMat (0,ngll-1-i)
            enddo
        endif
    endif
    return
end subroutine  getMass_Element2Face
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
