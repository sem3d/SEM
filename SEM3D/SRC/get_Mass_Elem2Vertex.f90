subroutine get_Mass_Elem2Vertex(Tdomain,n,rank)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n,rank
    integer :: ngllx,nglly,ngllz,i,nv,nnv


    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nv = 0,7
        nnv = Tdomain%specel(n)%Near_Vertices(nv)
        ! now we call the general assemblage routine
        call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,rank,  &
            Tdomain%sVertex(nnv)%MassMat,Tdomain%specel(n)%MassMat(:,:,:))
        if(Tdomain%sVertex(nnv)%PML)then
            call get_VectProperty_Elem2vertex(nv,ngllx,nglly,ngllz,rank,  &
                Tdomain%sVertex(nnv)%DumpMass(0:2),Tdomain%specel(n)%DumpMass(:,:,:,0:2))
            if(Tdomain%sVertex(nnv)%FPML)then
                call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,rank,  &
                    Tdomain%sVertex(nnv)%Ivx,Tdomain%specel(n)%Ivx(:,:,:))
                call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,rank,  &
                    Tdomain%sVertex(nnv)%Ivy,Tdomain%specel(n)%Ivy(:,:,:))
                call get_ScalarProperty_Elem2vertex(nv,ngllx,nglly,ngllz,rank,  &
                    Tdomain%sVertex(nnv)%Ivz,Tdomain%specel(n)%Ivz(:,:,:))
            end if
        end if

    enddo

    return

end subroutine get_Mass_Elem2Vertex
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
