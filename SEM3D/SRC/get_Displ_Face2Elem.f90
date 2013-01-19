!>
!! \file get_Displ_Face2Elem.f90
!! \brief Copie les Forces des faces sur les elements
!!
!<

subroutine get_Displ_Face2Elem(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: nf,i,j,ngllx,nglly,ngllz,ngll1,ngll2,nnf,orient_f

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nf = 0,5
        nnf = Tdomain%specel(n)%Near_Faces(nf)
        orient_f = Tdomain%specel(n)%Orient_Faces(nf)
        ngll1 = Tdomain%sFace(nnf)%ngll1
        ngll2 = Tdomain%sFace(nnf)%ngll2
        ! now we call the general deassemblage routine
        call get_VectProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
            Tdomain%sFace(nnf)%Forces(:,:,0:2),Tdomain%specel(n)%Forces(:,:,:,0:2))
    enddo

    return

end subroutine get_Displ_Face2Elem
!---------------------------------------------------------

!---------------------------------------------------------
subroutine get_Phi_Face2Elem(Tdomain,n)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: nf,i,j,ngllx,nglly,ngllz,ngll1,ngll2,nnf,orient_f

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    do nf = 0,5
        nnf = Tdomain%specel(n)%Near_Faces(nf)
        orient_f = Tdomain%specel(n)%Orient_Faces(nf)
        ngll1 = Tdomain%sFace(nnf)%ngll1
        ngll2 = Tdomain%sFace(nnf)%ngll2
        ! now we call the general deassemblage routine
        call get_ScalarProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,  &
            Tdomain%sFace(nnf)%ForcesFl(:,:),Tdomain%specel(n)%ForcesFl(:,:,:))
    enddo

    return

end subroutine get_Phi_Face2Elem
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
