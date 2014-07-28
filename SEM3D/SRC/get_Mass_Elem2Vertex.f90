!>
!! \file get_Mass_Elem2Vertex.f90
!! \brief
!!
!<
#if ! NEW_GLOBAL_METHOD
subroutine get_Mass_Elem2Vertex(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: ngllx,nglly,ngllz,i,nv,nnv


    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nv = 0,7
        nnv = Tdomain%specel(n)%Near_Vertices(nv)
        ! now we call the general assemblage routine
        call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
            Tdomain%sVertex(nnv)%MassMat,Tdomain%specel(n)%MassMat(:,:,:))
        if(Tdomain%sVertex(nnv)%PML)then
            call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                Tdomain%sVertex(nnv)%spml%DumpMass(0:2),Tdomain%specel(n)%spml%DumpMass(:,:,:,0:2))
            if(Tdomain%sVertex(nnv)%FPML)then
                call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%Ivx,Tdomain%specel(n)%spml%Ivx(:,:,:))
                call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%Ivy,Tdomain%specel(n)%spml%Ivy(:,:,:))
                call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,  &
                    Tdomain%sVertex(nnv)%spml%Ivz,Tdomain%specel(n)%spml%Ivz(:,:,:))
            end if
        end if

    enddo

    return

end subroutine get_Mass_Elem2Vertex
#endif
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
