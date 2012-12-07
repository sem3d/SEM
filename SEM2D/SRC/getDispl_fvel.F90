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


subroutine get_Displ_fv2el (Tdomain,n)

    ! Modified by Gaetano Festa 01/06/2005
    use sdomain
    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n

    integer :: nf,i,ngllx,ngllz

    ngllx = Tdomain%specel(n)%ngllx;  ngllz = Tdomain%specel(n)%ngllz
    nf = Tdomain%specel(n)%near_face(0)
    if ( Tdomain%sFace(nf)%near_element(0) == n) then
        Tdomain%specel(n)%Forces(1:ngllx-2,0,0) =  Tdomain%sFace(nf)%Forces(:,0)
        Tdomain%specel(n)%Forces(1:ngllx-2,0,1) =  Tdomain%sFace(nf)%Forces(:,1)
    else
        if ( Tdomain%sFace(nf)%coherency) then
            Tdomain%specel(n)%Forces(1:ngllx-2,0,0) =  Tdomain%sFace(nf)%Forces(:,0)
            ! 19/12/06 bug corrige par CS Tdomain%specel(n)%Forces(1:ngllx-1,0,1) =  Tdomain%sFace(nf)%Forces(:,1)
            Tdomain%specel(n)%Forces(1:ngllx-2,0,1) =  Tdomain%sFace(nf)%Forces(:,1)
        else
            do i = 1, ngllx-2
                Tdomain%specel(n)%Forces(i,0,0) =  Tdomain%sFace(nf)%Forces(ngllx-1-i,0)
                Tdomain%specel(n)%Forces(i,0,1) =  Tdomain%sFace(nf)%Forces(ngllx-1-i,1)
            enddo
        endif
    endif

    nf = Tdomain%specel(n)%near_face(1)
    if ( Tdomain%sFace(nf)%near_element(0) == n) then
        Tdomain%specel(n)%Forces(ngllx-1,1:ngllz-2,0) =  Tdomain%sFace(nf)%Forces(:,0)
        Tdomain%specel(n)%Forces(ngllx-1,1:ngllz-2,1) =  Tdomain%sFace(nf)%Forces(:,1)
    else
        if ( Tdomain%sFace(nf)%coherency) then
            Tdomain%specel(n)%Forces(ngllx-1,1:ngllz-2,0)=  Tdomain%sFace(nf)%Forces(:,0)
            Tdomain%specel(n)%Forces(ngllx-1,1:ngllz-2,1)=  Tdomain%sFace(nf)%Forces(:,1)
        else
            do i = 1, ngllz-2
                Tdomain%specel(n)%Forces(ngllx-1,i,0) =  Tdomain%sFace(nf)%Forces(ngllz-1-i,0)
                Tdomain%specel(n)%Forces(ngllx-1,i,1) =  Tdomain%sFace(nf)%Forces(ngllz-1-i,1)
            enddo
        endif
    endif
    nf = Tdomain%specel(n)%near_face(2)
    if ( Tdomain%sFace(nf)%near_element(0) == n) then
        Tdomain%specel(n)%Forces(1:ngllx-2,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(:,0)
        Tdomain%specel(n)%Forces(1:ngllx-2,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(:,1)
    else
        if ( Tdomain%sFace(nf)%coherency) then
            Tdomain%specel(n)%Forces(1:ngllx-2,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(:,0)
            Tdomain%specel(n)%Forces(1:ngllx-2,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(:,1)
        else
            do i = 1, ngllx-2
                Tdomain%specel(n)%Forces(i,ngllz-1,0) =  Tdomain%sFace(nf)%Forces(ngllx-1-i,0)
                Tdomain%specel(n)%Forces(i,ngllz-1,1) =  Tdomain%sFace(nf)%Forces(ngllx-1-i,1)
            enddo
        endif
    endif

    nf = Tdomain%specel(n)%near_face(3)
    if ( Tdomain%sFace(nf)%near_element(0) == n) then
        Tdomain%specel(n)%Forces(0,1:ngllz-2,0) =  Tdomain%sFace(nf)%Forces(:,0)
        Tdomain%specel(n)%Forces(0,1:ngllz-2,1) =  Tdomain%sFace(nf)%Forces(:,1)
    else
        if ( Tdomain%sFace(nf)%coherency) then
            Tdomain%specel(n)%Forces(0,1:ngllz-2,0)=  Tdomain%sFace(nf)%Forces(:,0)
            Tdomain%specel(n)%Forces(0,1:ngllz-2,1)=  Tdomain%sFace(nf)%Forces(:,1)
        else
            do i = 1, ngllz-2
                Tdomain%specel(n)%Forces(0,i,0) =  Tdomain%sFace(nf)%Forces(ngllz-1-i,0)
                Tdomain%specel(n)%Forces(0,i,1) =  Tdomain%sFace(nf)%Forces(ngllz-1-i,1)
            enddo
        endif
    endif

    nf = Tdomain%specel(n)%Near_Vertex(0)
    Tdomain%specel(n)%Forces (0,0,0:1) =Tdomain%sVertex(nf)%Forces(0:1)

    nf = Tdomain%specel(n)%Near_Vertex(1)
    Tdomain%specel(n)%Forces (ngllx-1,0,0:1) =Tdomain%sVertex(nf)%Forces(0:1)

    nf = Tdomain%specel(n)%Near_Vertex(2)
    Tdomain%specel(n)%Forces (ngllx-1,ngllz-1,0:1) =Tdomain%sVertex(nf)%Forces(0:1)

    nf = Tdomain%specel(n)%Near_Vertex(3)
    Tdomain%specel(n)%Forces (0,ngllz-1,0:1) =Tdomain%sVertex(nf)%Forces(0:1)

    return
end subroutine get_Displ_fv2el
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
