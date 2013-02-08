!>
!! \file get_Mass_Elem2Face.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine get_Mass_Elem2Face(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: ngllx,nglly,ngllz,ngll1,ngll2,i,j,k,nf,nnf,orient_f

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nf = 0,5
        nnf = Tdomain%specel(n)%Near_Faces(nf)
        orient_f = Tdomain%specel(n)%Orient_Faces(nf)
        ngll1 = Tdomain%sFace(nnf)%ngll1
        ngll2 = Tdomain%sFace(nnf)%ngll2

        ! now we call the general assemblage routine
        call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
            Tdomain%sFace(nnf)%MassMat(:,:),Tdomain%specel(n)%MassMat(:,:,:))
        if(Tdomain%sFace(nnf)%PML)then
            call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%spml%DumpMass(:,:,0:2),Tdomain%specel(n)%spml%DumpMass(:,:,:,0:2))
            if(Tdomain%sFace(nnf)%FPML)then
                call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%Ivx(:,:),Tdomain%specel(n)%spml%Ivx(:,:,:))
                call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%Ivy(:,:),Tdomain%specel(n)%spml%Ivy(:,:,:))
                call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                    Tdomain%sFace(nnf)%spml%Ivz(:,:),Tdomain%specel(n)%spml%Ivz(:,:,:))
            end if
        end if

    enddo

    return

end subroutine get_Mass_Elem2Face
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
