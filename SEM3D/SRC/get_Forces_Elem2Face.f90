!>
!! \file get_Forces_Elem2Face.f90
!! \brief Ajoute les Forces des Elements sur les faces
!!
!<

subroutine get_Forces_Elem2Face(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: ngllx,nglly,ngllz,ngll1,ngll2,nf,nnf,orient_f

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nf = 0,5
        nnf = Tdomain%specel(n)%Near_Faces(nf)
        orient_f = Tdomain%specel(n)%Orient_Faces(nf)
        ngll1 = Tdomain%sFace(nnf)%ngll1
        ngll2 = Tdomain%sFace(nnf)%ngll2
        ! now we call the general assemblage routine
        call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
            Tdomain%sFace(nnf)%Forces(:,:,0:2),Tdomain%specel(n)%Forces(:,:,:,0:2))
        if(Tdomain%sFace(nnf)%PML)then
            call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%spml%Forces1(:,:,0:2),Tdomain%specel(n)%spml%Forces1(:,:,:,0:2))
            call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%spml%Forces2(:,:,0:2),Tdomain%specel(n)%spml%Forces2(:,:,:,0:2))
            call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%spml%Forces3(:,:,0:2),Tdomain%specel(n)%spml%Forces3(:,:,:,0:2))
        end if
    enddo

    return

end subroutine get_Forces_Elem2Face
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
subroutine get_ForcesFl_Elem2Face(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: ngllx,nglly,ngllz,ngll1,ngll2,nf,nnf,orient_f

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
            Tdomain%sFace(nnf)%ForcesFl(:,:),Tdomain%specel(n)%ForcesFl(:,:,:))
        if(Tdomain%sFace(nnf)%PML)then
            call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%spml%ForcesFl1(:,:),Tdomain%specel(n)%spml%ForcesFl1(:,:,:))
            call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%spml%ForcesFl2(:,:),Tdomain%specel(n)%spml%ForcesFl2(:,:,:))
            call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
                Tdomain%sFace(nnf)%spml%ForcesFl3(:,:),Tdomain%specel(n)%spml%ForcesFl3(:,:,:))
        end if
    enddo

    return

end subroutine get_ForcesFl_Elem2Face
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
