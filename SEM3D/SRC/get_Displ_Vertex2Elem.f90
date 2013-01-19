!>
!! \file get_Displ_Vertex2Elem.f90
!! \brief Recopie les Forces des vertex sur les elements
!!
!<

subroutine get_Displ_Vertex2Elem(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: nv, ngllx,nglly,ngllz,nnv

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nv = 0,7
        nnv = Tdomain%specel(n)%Near_Vertices(nv)
        ! now we call the general deassemblage routine
        call get_VectProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,  &
            Tdomain%sVertex(nnv)%Forces(0:2),Tdomain%specel(n)%Forces(:,:,:,0:2))
    enddo

    return
end subroutine get_Displ_Vertex2Elem

!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
subroutine get_Phi_Vertex2Elem(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: nv, ngllx,nglly,ngllz,nnv

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nv = 0,7
        nnv = Tdomain%specel(n)%Near_Vertices(nv)
        ! now we call the general deassemblage routine
        call get_ScalarProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,  &
            Tdomain%sVertex(nnv)%ForcesFl,Tdomain%specel(n)%ForcesFl(:,:,:))
    enddo

    return
end subroutine get_Phi_Vertex2Elem
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
