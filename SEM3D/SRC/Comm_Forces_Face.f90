!>
!! \file Comm_Forces_Face.f90
!! \brief
!!
!<

subroutine Comm_Forces_Face(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngll,ngllPML,ngll_F,ngllPML_F
    integer :: ngll1,ngll2,nf,nnf,i,j,k,orient_f


    do nf = 0,Tdomain%sComm(n)%nb_faces-1
        nnf = Tdomain%sComm(n)%faces(nf)
        ngll1 = Tdomain%sFace(nnf)%ngll1
        ngll2 = Tdomain%sFace(nnf)%ngll2
        orient_f = Tdomain%sComm(n)%orient_faces(nf)
        if(Tdomain%sFace(nnf)%solid)then   ! solid part
            call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                      &
                Tdomain%sComm(n)%TakeForces(ngll:ngll+(ngll1-2)*(ngll2-2)-1,:),  &
                Tdomain%sFace(nnf)%Forces(:,:,:))
            ngll = ngll+(ngll1-2)*(ngll2-2)
            if(Tdomain%sFace(nnf)%PML)then
                call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                               &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+(ngll1-2)*(ngll2-2)-1,1,0:2),  &
                    Tdomain%sFace(nnf)%Forces1(:,:,:))
                call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                               &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+(ngll1-2)*(ngll2-2)-1,2,0:2),  &
                    Tdomain%sFace(nnf)%Forces2(:,:,:))
                call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                               &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+(ngll1-2)*(ngll2-2)-1,3,0:2),  &
                    Tdomain%sFace(nnf)%Forces3(:,:,:))
                ngllPML = ngllPML+(ngll1-2)*(ngll2-2)
            end if

        else   ! fluid part
            call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                          &
                Tdomain%sComm(n)%TakeForcesFl(ngll_F:ngll_F+(ngll1-2)*(ngll2-2)-1),  &
                Tdomain%sFace(nnf)%ForcesFl(:,:))
            ngll_F = ngll_F+(ngll1-2)*(ngll2-2)
            if(Tdomain%sFace(nnf)%PML)then
                call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+(ngll1-2)*(ngll2-2)-1,1),  &
                    Tdomain%sFace(nnf)%ForcesFl1(:,:))
                call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+(ngll1-2)*(ngll2-2)-1,2),  &
                    Tdomain%sFace(nnf)%ForcesFl2(:,:))
                call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+(ngll1-2)*(ngll2-2)-1,3),  &
                    Tdomain%sFace(nnf)%ForcesFl3(:,:))
                ngllPML_F = ngllPML_F+(ngll1-2)*(ngll2-2)
            end if

        end if
    end do

    return
end subroutine Comm_Forces_Face
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
