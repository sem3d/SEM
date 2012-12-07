!>
!! \file getPMLpred_e2el.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine get_PMLprediction_e2el (Tdomain, n, bega, dt)

    use sdomain

    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n
    real, intent(IN) :: dt, bega

    integer :: ne, ngllx, nglly, ngllz, i, ngll


    ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly;  ngllz = Tdomain%specel(n)%ngllz

    ne = Tdomain%specel(n)%near_edges(0)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(0) == 0 ) then
        Tdomain%specel(n)%Forces(1:ngllx-2,0,0,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(ngllx-1-i,0,0,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif

    ne = Tdomain%specel(n)%near_edges(1)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(1) == 0 ) then
        Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,0,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,0,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif

    ne = Tdomain%specel(n)%near_edges(2)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(2) == 0 ) then
        Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,0,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,0,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif


    ne = Tdomain%specel(n)%near_edges(3)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(3) == 0 ) then
        Tdomain%specel(n)%Forces(0,1:nglly-2,0,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(0,nglly-1-i,0,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif


    ne = Tdomain%specel(n)%near_edges(4)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(4) == 0 ) then
        Tdomain%specel(n)%Forces(ngllx-1,0,1:ngllz-2,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1-i,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif

    ne = Tdomain%specel(n)%near_edges(5)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(5) == 0 ) then
        Tdomain%specel(n)%Forces(1:ngllx-2,0,ngllz-1,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(ngllx-1-i,0,ngllz-1,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif

    ne = Tdomain%specel(n)%near_edges(6)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(6) == 0 ) then
        Tdomain%specel(n)%Forces(0,0,1:ngllz-2,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(0,0,ngllz-1-i,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif

    ne = Tdomain%specel(n)%near_edges(7)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(7) == 0 ) then
        Tdomain%specel(n)%Forces(ngllx-1,nglly-1,1:ngllz-2,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1-i,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif

    ne = Tdomain%specel(n)%near_edges(8)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(8) == 0 ) then
        Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,ngllz-1,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(ngllx-1,nglly-1-i,ngllz-1,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif

    ne = Tdomain%specel(n)%near_edges(9)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(9) == 0 ) then
        Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,ngllz-1,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(ngllx-1-i,nglly-1,ngllz-1,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif

    ne = Tdomain%specel(n)%near_edges(10)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(10) == 0 ) then
        Tdomain%specel(n)%Forces(0,nglly-1,1:ngllz-2,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1-i,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif

    ne = Tdomain%specel(n)%near_edges(11)
    ngll = Tdomain%sEdge(ne)%ngll
    if ( Tdomain%specel(n)%Orient_Edges(11) == 0 ) then
        Tdomain%specel(n)%Forces(0,1:nglly-2,ngllz-1,:) = Tdomain%sEdge(ne)%Veloc(:,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(:,:)
    else
        do i=1,ngll-2
            Tdomain%specel(n)%Forces(0,nglly-1-i,ngllz-1,:) = Tdomain%sEdge(ne)%Veloc(i,:) + dt * (0.5 - bega) *  Tdomain%sEdge(ne)%Accel(i,:)
        enddo
    endif


    return
end subroutine get_PMLprediction_e2el
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
