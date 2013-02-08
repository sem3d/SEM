!>
!! \file get_Mass_Elem2Edge.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine get_Mass_Elem2Edge(Tdomain,n)

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
            Tdomain%sEdge(nne)%MassMat,Tdomain%specel(n)%MassMat)
        if(Tdomain%sEdge(nne)%PML)then
            call get_VectProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                Tdomain%sEdge(nne)%spml%DumpMass(:,0:2),Tdomain%specel(n)%spml%DumpMass(:,:,:,0:2))
            if(Tdomain%sEdge(nne)%FPML)then
                call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%Ivx(:),Tdomain%specel(n)%spml%Ivx(:,:,:))
                call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%Ivy(:),Tdomain%specel(n)%spml%Ivy(:,:,:))
                call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,  &
                    Tdomain%sEdge(nne)%spml%Ivz(:),Tdomain%specel(n)%spml%Ivz(:,:,:))
            end if
        end if

    enddo

    return

end subroutine get_Mass_Elem2Edge
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
