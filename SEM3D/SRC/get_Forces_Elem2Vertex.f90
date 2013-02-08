!>
!! \file get_Forces_Elem2Vertex.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine get_Forces_Elem2Vertex(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: ngllx,nglly,ngllz,nv,nnv


    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nv = 0,7
        nnv = Tdomain%specel(n)%Near_Vertices(nv)
        ! now we call the general assemblage routine
        call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
            Tdomain%sVertex(nnv)%Forces(0:2),Tdomain%specel(n)%Forces(:,:,:,0:2))
        if(Tdomain%sVertex(nnv)%PML)then
            call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%spml%Forces1(0:2),Tdomain%specel(n)%spml%Forces1(:,:,:,0:2))
            call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%spml%Forces2(0:2),Tdomain%specel(n)%spml%Forces2(:,:,:,0:2))
            call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%spml%Forces3(0:2),Tdomain%specel(n)%spml%Forces3(:,:,:,0:2))
        end if
    enddo

    return

end subroutine get_Forces_Elem2Vertex
!---------------------------------------------------------------------------------

!---------------------------------------------------------------------------------
subroutine get_ForcesFl_Elem2Vertex(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: ngllx,nglly,ngllz,nv,nnv


    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nv = 0,7
        nnv = Tdomain%specel(n)%Near_Vertices(nv)
        ! now we call the general assemblage routine
        call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
            Tdomain%sVertex(nnv)%ForcesFl,Tdomain%specel(n)%ForcesFl(:,:,:))
        if(Tdomain%sVertex(nnv)%PML)then
            call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%spml%ForcesFl1,Tdomain%specel(n)%spml%ForcesFl1(:,:,:))
            call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%spml%ForcesFl2,Tdomain%specel(n)%spml%ForcesFl2(:,:,:))
            call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%spml%ForcesFl3,Tdomain%specel(n)%spml%ForcesFl3(:,:,:))
        end if
    enddo

    return

end subroutine get_ForcesFl_Elem2Vertex
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
