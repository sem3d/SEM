!>
!!\file getPMLpred_fv2el.F90
!!\brief Contient la subroutine get_PMLprediction_fv2el.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type (Domain), intent (IN) Tdomain
!! \param integer, intent (IN) n
!! \param integer, intent (IN) ngllx
!! \param integer, intent (IN) ngllz
!! \param real, intent(IN) dt
!! \param real, intent(IN) bega
!! \param real, intent(IN) alpha
!! \param real, dimension (0:ngllx-1,0:ngllz-1), intent(INOUT) Vxloc
!! \param real, dimension (0:ngllx-1,0:ngllz-1), intent(INOUT) Vzloc
!<


subroutine get_PMLprediction_fv2el (Tdomain,n,Vxloc,Vzloc,ngllx,ngllz,alpha, bega,dt)

    ! Modified by Gaetano Festa 01/06/2005
    use sdomain
    implicit none
    type (Domain), intent (IN) :: Tdomain
    integer, intent (IN) :: n, ngllx,ngllz
    real, intent(IN) :: dt, bega, alpha
    real, dimension (0:ngllx-1,0:ngllz-1), intent(INOUT) :: Vxloc,Vzloc

    integer :: nf,i

    nf = Tdomain%specel(n)%near_face(0)
    if ( Tdomain%sFace(nf)%near_element(0) == n) then
        Vxloc(1:ngllx-2,0) = (0.5+alpha)*  Tdomain%sFace(nf)%Veloc(:,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(:,0) + &
            (0.5-alpha) * Tdomain%sFace(nf)%V0(:,0)
        Vzloc(1:ngllx-2,0) =  (0.5+alpha)* Tdomain%sFace(nf)%Veloc(:,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(:,1) + &
            (0.5-alpha)*Tdomain%sFace(nf)%V0(:,1)
    else
        if ( Tdomain%sFace(nf)%coherency) then
            Vxloc(1:ngllx-2,0) =  (0.5+alpha)*  Tdomain%sFace(nf)%Veloc(:,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(:,0) + (0.5-alpha) * Tdomain%sFace(nf)%V0(:,0)
            Vzloc(1:ngllx-2,0) =   (0.5+alpha)* Tdomain%sFace(nf)%Veloc(:,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(:,1) + (0.5-alpha)*Tdomain%sFace(nf)%V0(:,1)
        else
            do i = 1, ngllx-2
                Vxloc(i,0) =  (0.5+alpha)*  Tdomain%sFace(nf)%Veloc(ngllx-1-i,0) + dt * (0.5 - bega) * &
                    Tdomain%sFace(nf)%Accel(ngllx-1-i,0) + (0.5-alpha)*Tdomain%sFace(nf)%V0(ngllx-1-i,0)
                Vzloc(i,0) =  (0.5+alpha)*  Tdomain%sFace(nf)%Veloc(ngllx-1-i,1) + dt * (0.5 - bega) *  &
                    Tdomain%sFace(nf)%Accel(ngllx-1-i,1)+ (0.5-alpha)*  Tdomain%sFace(nf)%V0(ngllx-1-i,1)
            enddo
        endif
    endif

    nf = Tdomain%specel(n)%near_face(1)
    if ( Tdomain%sFace(nf)%near_element(0) == n) then
        Vxloc(ngllx-1,1:ngllz-2) =   (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,0) + dt * (0.5 - bega) *  &
            Tdomain%sFace(nf)%Accel(:,0) + (0.5-alpha)*Tdomain%sFace(nf)%V0(:,0)
        Vzloc(ngllx-1,1:ngllz-2) =   (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,1) + dt * (0.5 - bega) *   &
            Tdomain%sFace(nf)%Accel(:,1) + (0.5-alpha)*Tdomain%sFace(nf)%V0(:,1)
    else
        if ( Tdomain%sFace(nf)%coherency) then
            Vxloc(ngllx-1,1:ngllz-2)=  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,0) + dt * (0.5 - bega) * &
                Tdomain%sFace(nf)%Accel(:,0) + (0.5-alpha)*Tdomain%sFace(nf)%V0(:,0)
            Vzloc(ngllx-1,1:ngllz-2)=  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,1) + dt * (0.5 - bega) *  &
                Tdomain%sFace(nf)%Accel(:,1) + (0.5-alpha)*Tdomain%sFace(nf)%V0(:,1)
        else
            do i = 1, ngllz-2
                Vxloc(ngllx-1,i) =   (0.5+alpha)*Tdomain%sFace(nf)%Veloc(ngllz-1-i,0) + dt * (0.5 - bega) *  &
                    Tdomain%sFace(nf)%Accel(ngllz-1-i,0) + (0.5-alpha)* Tdomain%sFace(nf)%V0(ngllz-1-i,0)
                Vzloc(ngllx-1,i) =   (0.5+alpha)*Tdomain%sFace(nf)%Veloc(ngllz-1-i,1) + dt * (0.5 - bega) *  &
                    Tdomain%sFace(nf)%Accel(ngllz-1-i,1) + (0.5-alpha)* Tdomain%sFace(nf)%V0(ngllz-1-i,1)
            enddo
        endif
    endif
    nf = Tdomain%specel(n)%near_face(2)
    if ( Tdomain%sFace(nf)%near_element(0) == n) then
        Vxloc(1:ngllx-2,ngllz-1) = (0.5+alpha)* Tdomain%sFace(nf)%Veloc(:,0) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(:,0) + &
            (0.5-alpha)* Tdomain%sFace(nf)%V0(:,0)
        Vzloc(1:ngllx-2,ngllz-1) =  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,1) + dt * (0.5 - bega) *  Tdomain%sFace(nf)%Accel(:,1) + &
            (0.5-alpha)*Tdomain%sFace(nf)%V0(:,1)
    else
        if ( Tdomain%sFace(nf)%coherency) then
            Vxloc(1:ngllx-2,ngllz-1) =  (0.5+alpha)* Tdomain%sFace(nf)%Veloc(:,0) + dt * (0.5 - bega) * &
                Tdomain%sFace(nf)%Accel(:,0) + (0.5-alpha)* Tdomain%sFace(nf)%V0(:,0)
            Vzloc(1:ngllx-2,ngllz-1) =  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,1) + dt * (0.5 - bega) *  &
                Tdomain%sFace(nf)%Accel(:,1) +  (0.5-alpha)*Tdomain%sFace(nf)%V0(:,1)
        else
            do i = 1, ngllx-2
                Vxloc(i,ngllz-1) =  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(ngllx-1-i,0) + dt * (0.5 - bega) * &
                    Tdomain%sFace(nf)%Accel(ngllx-1-i,0) + (0.5-alpha)*Tdomain%sFace(nf)%V0(ngllx-1-i,0)
                Vzloc(i,ngllz-1) =  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(ngllx-1-i,1) + dt * (0.5 - bega) * &
                    Tdomain%sFace(nf)%Accel(ngllx-1-i,1)+ (0.5-alpha)*Tdomain%sFace(nf)%V0(ngllx-1-i,1)
            enddo
        endif
    endif

    nf = Tdomain%specel(n)%near_face(3)
    if ( Tdomain%sFace(nf)%near_element(0) == n) then
        Vxloc(0,1:ngllz-2) =  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,0) + dt * (0.5 - bega) *  &
            Tdomain%sFace(nf)%Accel(:,0) + (0.5-alpha)*Tdomain%sFace(nf)%V0(:,0)
        Vzloc(0,1:ngllz-2) =  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,1) + dt * (0.5 - bega) *  &
            Tdomain%sFace(nf)%Accel(:,1)+ (0.5-alpha)*Tdomain%sFace(nf)%V0(:,1)
    else
        if ( Tdomain%sFace(nf)%coherency) then
            Vxloc(0,1:ngllz-2)=  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,0) + dt * (0.5 - bega) *  &
                Tdomain%sFace(nf)%Accel(:,0) + (0.5-alpha)*Tdomain%sFace(nf)%V0(:,0)
            Vzloc(0,1:ngllz-2)=  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(:,1) + dt * (0.5 - bega) *  &
                Tdomain%sFace(nf)%Accel(:,1)+ (0.5-alpha)*Tdomain%sFace(nf)%V0(:,1)
        else
            do i = 1, ngllz-2
                Vxloc(0,i) =  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(ngllz-1-i,0) + dt * (0.5 - bega) *  &
                    Tdomain%sFace(nf)%Accel(ngllz-1-i,0) + (0.5-alpha)*Tdomain%sFace(nf)%V0(ngllz-1-i,0)
                Vzloc(0,i) =  (0.5+alpha)*Tdomain%sFace(nf)%Veloc(ngllz-1-i,1) + dt * (0.5 - bega) *  &
                    Tdomain%sFace(nf)%Accel(ngllz-1-i,1) + (0.5-alpha)*Tdomain%sFace(nf)%V0(ngllz-1-i,1)
            enddo
        endif
    endif

    nf = Tdomain%specel(n)%Near_Vertex(0)
    Vxloc(0,0) =  (0.5+alpha)*Tdomain%sVertex(nf)%Veloc(0) + dt * (0.5 - bega) *  &
        Tdomain%sVertex(nf)%Accel(0) + (0.5-alpha)*Tdomain%sVertex(nf)%V0(0)
    Vzloc(0,0) =  (0.5+alpha)*Tdomain%sVertex(nf)%Veloc(1) + dt * (0.5 - bega) *  &
        Tdomain%sVertex(nf)%Accel(1) + (0.5-alpha)*Tdomain%sVertex(nf)%V0(1)

    nf = Tdomain%specel(n)%Near_Vertex(1)
    Vxloc(ngllx-1,0) =  (0.5+alpha)*Tdomain%sVertex(nf)%Veloc(0) + dt * (0.5 - bega) *  &
        Tdomain%sVertex(nf)%Accel(0) + (0.5-alpha)*Tdomain%sVertex(nf)%V0(0)
    Vzloc(ngllx-1,0) =  (0.5+alpha)*Tdomain%sVertex(nf)%Veloc(1) + dt * (0.5 - bega) *  &
        Tdomain%sVertex(nf)%Accel(1) + (0.5-alpha)*Tdomain%sVertex(nf)%V0(1)

    nf = Tdomain%specel(n)%Near_Vertex(2)
    Vxloc(ngllx-1,ngllz-1) =  (0.5+alpha)*Tdomain%sVertex(nf)%Veloc(0) + dt * (0.5 - bega) *  &
        Tdomain%sVertex(nf)%Accel(0) + (0.5-alpha)*Tdomain%sVertex(nf)%V0(0)
    Vzloc(ngllx-1,ngllz-1) =  (0.5+alpha)*Tdomain%sVertex(nf)%Veloc(1) + dt * (0.5 - bega) *  &
        Tdomain%sVertex(nf)%Accel(1) + (0.5-alpha)*Tdomain%sVertex(nf)%V0(1)

    nf = Tdomain%specel(n)%Near_Vertex(3)
    Vxloc(0,ngllz-1) =  (0.5+alpha)*Tdomain%sVertex(nf)%Veloc(0) + dt * (0.5 - bega) *  &
        Tdomain%sVertex(nf)%Accel(0) + (0.5-alpha)*Tdomain%sVertex(nf)%V0(0)
    Vzloc(0,ngllz-1) =  (0.5+alpha)*Tdomain%sVertex(nf)%Veloc(1) + dt * (0.5 - bega) *  &
        Tdomain%sVertex(nf)%Accel(1) + (0.5-alpha)*Tdomain%sVertex(nf)%V0(1)


    return
end subroutine get_PMLprediction_fv2el


subroutine get_simpler (Tdomain,n,Vxloc,Vzloc,ngllx,ngllz)

    ! Modified by S Terrana 09/09/2014
    use sdomain
    implicit none
    type (Domain), intent (IN) :: Tdomain
    integer, intent (IN) :: n, ngllx,ngllz
    real, dimension (0:ngllx-1,0:ngllz-1), intent(INOUT) :: Vxloc,Vzloc

    integer :: nf,nv

    nf = Tdomain%specel(n)%near_face(0)
    Vxloc(1:ngllx-2,0) = Tdomain%sFace(nf)%Veloc(:,0)
    Vzloc(1:ngllx-2,0) = Tdomain%sFace(nf)%Veloc(:,1)

    nf = Tdomain%specel(n)%near_face(1)
    Vxloc(ngllx-1,1:ngllz-2) = Tdomain%sFace(nf)%Veloc(:,0)
    Vzloc(ngllx-1,1:ngllz-2) = Tdomain%sFace(nf)%Veloc(:,1)

    nf = Tdomain%specel(n)%near_face(2)
    Vxloc(1:ngllx-2,ngllz-1) = Tdomain%sFace(nf)%Veloc(:,0)
    Vzloc(1:ngllx-2,ngllz-1) = Tdomain%sFace(nf)%Veloc(:,1)

    nf = Tdomain%specel(n)%near_face(3)
    Vxloc(0,1:ngllz-2) = Tdomain%sFace(nf)%Veloc(:,0)
    Vzloc(0,1:ngllz-2) = Tdomain%sFace(nf)%Veloc(:,1)

    nv = Tdomain%specel(n)%near_Vertex(0)
    Vxloc(0,0) = Tdomain%sVertex(nv)%Veloc(0)
    Vzloc(0,0) = Tdomain%sVertex(nv)%Veloc(1)

    nv = Tdomain%specel(n)%near_Vertex(1)
    Vxloc(ngllx-1,0) = Tdomain%sVertex(nv)%Veloc(0)
    Vzloc(ngllx-1,0) = Tdomain%sVertex(nv)%Veloc(1)

    nv = Tdomain%specel(n)%near_Vertex(2)
    Vxloc(ngllx-1,ngllz-1) = Tdomain%sVertex(nv)%Veloc(0)
    Vzloc(ngllx-1,ngllz-1) = Tdomain%sVertex(nv)%Veloc(1)

    nv = Tdomain%specel(n)%near_Vertex(3)
    Vxloc(0,ngllz-1) = Tdomain%sVertex(nv)%Veloc(0)
    Vzloc(0,ngllz-1) = Tdomain%sVertex(nv)%Veloc(1)

    return
end subroutine get_simpler
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
