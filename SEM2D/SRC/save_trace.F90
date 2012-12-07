!>
!!\file save_trace.F90
!!\brief Contient la subroutine save_trace.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type (domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) it
!<


subroutine save_trace (Tdomain, it)

    !Modified by Gaetano Festa 16/06/2005
    use sdomain

    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: it

    integer :: ir, nr, nf, i,j, ngllx, ngllz, nsta
    real :: dum0, dum1
    real, dimension (:,:,:), allocatable :: Veloc


    nsta = Tdomain%n_receivers
    do ir = 0, nsta-1
        if (Tdomain%sReceiver(ir)%located_here) then
            nr = Tdomain%sReceiver(ir)%nr
            dum0 = 0; dum1 = 0
            ngllx = Tdomain%specel(nr)%ngllx
            ngllz = Tdomain%specel(nr)%ngllz

            allocate (Veloc(0:ngllx-1,0:ngllz-1,0:1))

            Veloc(1:ngllx-2,1:ngllz-2,0:1) = Tdomain%specel(nr)%Veloc(1:ngllx-2,1:ngllz-2,0:1)

            nf = Tdomain%specel(nr)%near_face(0)
            if ( Tdomain%sFace(nf)%near_element(0) == nr) then
                Veloc(1:ngllx-2,0,0) =  Tdomain%sFace(nf)%Veloc(:,0)
                Veloc(1:ngllx-2,0,1) =  Tdomain%sFace(nf)%Veloc(:,1)
            else
                if ( Tdomain%sFace(nf)%coherency) then
                    Veloc(1:ngllx-2,0,0) =  Tdomain%sFace(nf)%Veloc(:,0)
                    Veloc(1:ngllx-2,0,1) =  Tdomain%sFace(nf)%Veloc(:,1)
                else
                    do i = 1, ngllx-2
                        Veloc(i,0,0) =  Tdomain%sFace(nf)%Veloc(ngllx-1-i,0)
                        Veloc(i,0,1) =  Tdomain%sFace(nf)%Veloc(ngllx-1-i,1)
                    enddo
                endif
            endif

            nf = Tdomain%specel(nr)%near_face(1)
            if ( Tdomain%sFace(nf)%near_element(0) == nr) then
                Veloc(ngllx-1,1:ngllz-2,0) =  Tdomain%sFace(nf)%Veloc(:,0)
                Veloc(ngllx-1,1:ngllz-2,1) =  Tdomain%sFace(nf)%Veloc(:,1)
            else
                if ( Tdomain%sFace(nf)%coherency) then
                    Veloc(ngllx-1,1:ngllz-2,0)=  Tdomain%sFace(nf)%Veloc(:,0)
                    Veloc(ngllx-1,1:ngllz-2,1)=  Tdomain%sFace(nf)%Veloc(:,1)
                else
                    do i = 1, ngllz-2
                        Veloc(ngllx-1,i,0) =  Tdomain%sFace(nf)%Veloc(ngllz-1-i,0)
                        Veloc(ngllx-1,i,1) =  Tdomain%sFace(nf)%Veloc(ngllz-1-i,1)
                    enddo
                endif
            endif
            nf = Tdomain%specel(nr)%near_face(2)
            if ( Tdomain%sFace(nf)%near_element(0) == nr) then
                Veloc(1:ngllx-2,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(:,0)
                Veloc(1:ngllx-2,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(:,1)
            else
                if ( Tdomain%sFace(nf)%coherency) then
                    Veloc(1:ngllx-2,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(:,0)
                    Veloc(1:ngllx-2,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(:,1)
                else
                    do i = 1, ngllx-2
                        Veloc(i,ngllz-1,0) =  Tdomain%sFace(nf)%Veloc(ngllx-1-i,0)
                        Veloc(i,ngllz-1,1) =  Tdomain%sFace(nf)%Veloc(ngllx-1-i,1)
                    enddo
                endif
            endif

            nf = Tdomain%specel(nr)%near_face(3)
            if ( Tdomain%sFace(nf)%near_element(0) == nr) then
                Veloc(0,1:ngllz-2,0) =  Tdomain%sFace(nf)%Veloc(:,0)
                Veloc(0,1:ngllz-2,1) =  Tdomain%sFace(nf)%Veloc(:,1)
            else
                if ( Tdomain%sFace(nf)%coherency) then
                    Veloc(0,1:ngllz-2,0)=  Tdomain%sFace(nf)%Veloc(:,0)
                    Veloc(0,1:ngllz-2,1)=  Tdomain%sFace(nf)%Veloc(:,1)
                else
                    do i = 1, ngllz-2
                        Veloc(0,i,0) =  Tdomain%sFace(nf)%Veloc(ngllz-1-i,0)
                        Veloc(0,i,1) =  Tdomain%sFace(nf)%Veloc(ngllz-1-i,1)
                    enddo
                endif
            endif

            nf = Tdomain%specel(nr)%Near_Vertex(0)
            Veloc(0,0,0:1) = Tdomain%sVertex(nf)%Veloc(0:1)
            nf = Tdomain%specel(nr)%Near_Vertex(1)
            Veloc(ngllx-1,0,0:1) = Tdomain%sVertex(nf)%Veloc(0:1)
            nf = Tdomain%specel(nr)%Near_Vertex(2)
            Veloc(ngllx-1,ngllz-1,0:1) = Tdomain%sVertex(nf)%Veloc(0:1)
            nf = Tdomain%specel(nr)%Near_Vertex(3)
            Veloc(0,ngllz-1,0:1) = Tdomain%sVertex(nf)%Veloc(0:1)

            do j = 0,ngllz-1
                do i =0,ngllx -1
                    dum0 = dum0 + Tdomain%sReceiver(ir)%Interp_Coeff(i,j) * Veloc(i,j,0)
                    dum1 = dum1 + Tdomain%sReceiver(ir)%Interp_Coeff(i,j) * Veloc(i,j,1)
                enddo
            enddo
            Tdomain%Store_Trace(0,ir,it) = dum0;  Tdomain%Store_Trace(1,ir,it) = dum1;

            deallocate (Veloc)
        endif
    enddo

    return
end subroutine save_trace
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
