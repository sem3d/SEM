!>
!! \file getPMLpred_v2el.f90
!! \brief
!!
!<

subroutine get_PMLprediction_v2el (Tdomain, n, bega, dt)

    use sdomain

    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n
    real, intent(IN) :: dt, bega

    integer :: nv, ngllx, nglly, ngllz, nnv

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nv = 0,7
        nnv = Tdomain%specel(n)%Near_Vertices(nv)
        ! now we call the general deassemblage routine
        call get_VectProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,                             &
            Tdomain%sVertex(nnv)%Veloc(0:2)+dt*(0.5-bega)*Tdomain%sVertex(nnv)%Accel(0:2),  &
            Tdomain%specel(n)%Forces(:,:,:,0:2))
    enddo

    return
end subroutine get_PMLprediction_v2el
!----------------------------------------------------------

!----------------------------------------------------------
subroutine get_PMLprediction_v2el_fl(Tdomain,n,bega,dt)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    real, intent(in)   :: bega,dt
    integer :: nv, ngllx,nglly,ngllz,nnv

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nv = 0,7
        nnv = Tdomain%specel(n)%Near_Vertices(nv)
        ! now we call the general deassemblage routine
        call get_ScalarProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,                     &
            Tdomain%sVertex(nnv)%VelPhi+dt*(0.5-bega)*Tdomain%sVertex(nnv)%AccelPhi,  &
            Tdomain%specel(n)%ForcesFl(:,:,:))
    enddo

    return
end subroutine get_PMLprediction_v2el_fl
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
