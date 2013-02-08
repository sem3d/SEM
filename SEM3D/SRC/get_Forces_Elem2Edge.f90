!>
!! \file get_Forces_Elem2Edge.f90
!! \brief Ajout des forces elements sur les aretes
!!
!<
subroutine get_Forces_Elem2Edge(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: ngllx,nglly,ngllz,ngll,ne,nne,orient_e

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do ne = 0,11
        nne = Tdomain%specel(n)%Near_Edges(ne)
        orient_e = Tdomain%specel(n)%Orient_Edges(ne)
        ngll = Tdomain%sEdge(nne)%ngll
        ! now we call the general assemblage routine
        call get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
            Tdomain%sEdge(nne)%Forces(:,0:2),Tdomain%specel(n)%Forces(:,:,:,0:2))
        if(Tdomain%sEdge(nne)%PML)then
            call get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%spml%Forces1(:,0:2),Tdomain%specel(n)%spml%Forces1(:,:,:,0:2))
            call get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%spml%Forces2(:,0:2),Tdomain%specel(n)%spml%Forces2(:,:,:,0:2))
            call get_VectProperty_Elem2Edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%spml%Forces3(:,0:2),Tdomain%specel(n)%spml%Forces3(:,:,:,0:2))
        end if
    end do

    return

end subroutine get_Forces_Elem2Edge
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
subroutine get_ForcesFl_Elem2Edge(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: ngllx,nglly,ngllz,ngll,ne,nne,orient_e

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do ne = 0,11
        nne = Tdomain%specel(n)%Near_Edges(ne)
        orient_e = Tdomain%specel(n)%Orient_Edges(ne)
        ngll = Tdomain%sEdge(nne)%ngll
        ! now we call the general assemblage routine
        call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
            Tdomain%sEdge(nne)%ForcesFl(:),Tdomain%specel(n)%ForcesFl(:,:,:))
        if(Tdomain%sEdge(nne)%PML)then
            call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%spml%ForcesFl1(:),Tdomain%specel(n)%spml%ForcesFl1(:,:,:))
            call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%spml%ForcesFl2(:),Tdomain%specel(n)%spml%ForcesFl2(:,:,:))
            call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%spml%ForcesFl3(:),Tdomain%specel(n)%spml%ForcesFl3(:,:,:))
        end if
    end do

    return

end subroutine get_ForcesFl_Elem2Edge
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
