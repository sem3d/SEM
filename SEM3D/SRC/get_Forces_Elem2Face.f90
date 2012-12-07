subroutine get_Forces_Elem2Face(Tdomain,n,rank)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n,rank
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
    call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
             Tdomain%sFace(nnf)%Forces(:,:,0:2),Tdomain%specel(n)%Forces(:,:,:,0:2)) 
    if(Tdomain%sFace(nnf)%PML)then
        call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
                  Tdomain%sFace(nnf)%Forces1(:,:,0:2),Tdomain%specel(n)%Forces1(:,:,:,0:2)) 
        call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
                  Tdomain%sFace(nnf)%Forces2(:,:,0:2),Tdomain%specel(n)%Forces2(:,:,:,0:2)) 
        call get_VectProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
                  Tdomain%sFace(nnf)%Forces3(:,:,0:2),Tdomain%specel(n)%Forces3(:,:,:,0:2)) 
    end if
enddo

return

end subroutine get_Forces_Elem2Face
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine get_ForcesFl_Elem2Face(Tdomain,n,rank)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n,rank
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
    call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
             Tdomain%sFace(nnf)%ForcesFl(:,:),Tdomain%specel(n)%ForcesFl(:,:,:)) 
    if(Tdomain%sFace(nnf)%PML)then
        call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
                  Tdomain%sFace(nnf)%ForcesFl1(:,:),Tdomain%specel(n)%ForcesFl1(:,:,:)) 
        call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
                  Tdomain%sFace(nnf)%ForcesFl2(:,:),Tdomain%specel(n)%ForcesFl2(:,:,:)) 
        call get_ScalarProperty_Elem2face(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
                  Tdomain%sFace(nnf)%ForcesFl3(:,:),Tdomain%specel(n)%ForcesFl3(:,:,:)) 
    end if
enddo

return

end subroutine get_ForcesFl_Elem2Face
