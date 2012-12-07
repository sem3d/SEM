!>
!!\file verify_in_quad.F90
!!\brief Contient principalement la subroutine verify_in_quad.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Vérifie si un point est dans un rectangle.
!!
!! \param real, dimension (0:3), intent (IN) xc
!! \param real, dimension (0:3), intent (IN) yc
!! \param real, intent (IN) x
!! \param real, intent (IN) y
!! \param logical, intent (OUT) inner
!<


subroutine verify_in_quad (xc,yc,x,y,inner)
    implicit none

    ! ##############################
    ! From mathforum
    ! Formula for a point in a rectangle within proof
    !
    ! Gaetano Festa, modified 16/10/2004
    ! #############################

    real, dimension (0:3), intent (IN) :: xc,yc
    real, intent (IN) :: x,y
    logical, intent (OUT) :: inner
    real :: signA,signA0,A

    inner = .false.
    call two_area(xc(0),yc(0),xc(1),yc(1),x,y,A)
    if (A /= 0)  then
        signA0 = A/abs(A)
        call two_area(xc(1),yc(1),xc(2),yc(2),x,y,A)
        if (A /= 0 ) then
            signA=signA0 * A/abs(A)
            if (signA < 0 ) return
        endif
    else
        call two_area(xc(1),yc(1),xc(2),yc(2),x,y,A)
        if (A==0) inner = .true.
        signA0 = A/abs(A)
    endif
    call two_area(xc(2),yc(2),xc(3),yc(3),x,y,A)
    if (A /= 0 ) then
        signA=signA0 * A/abs(A)
        if (signA < 0 ) return
    endif
    call two_area(xc(3),yc(3),xc(0),yc(0),x,y,A)
    if (A /= 0 ) then
        signA=signA0 * A/abs(A)
        if (signA < 0 ) return
    endif
    inner = .true.
    return
end subroutine verify_in_quad

! ################################################################
!>
!! \brief
!!
!! \param real, intent (IN) x1
!! \param real, intent (IN) y1
!! \param real, intent (IN) x2
!! \param real, intent (IN) y2
!! \param real, intent (IN) x3
!! \param real, intent (IN) y3
!! \param real, intent (OUT) A
!<


subroutine two_area (x1,y1,x2,y2,x3,y3,A)
    implicit none
    real, intent (IN) :: x1,y1,x2,y2,x3,y3
    real, intent (OUT) :: A

    A = x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2)
    return
end subroutine two_area
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
