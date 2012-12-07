subroutine get_PMLprediction_f2el(Tdomain,n,bega,dt,rank)

    use sdomain
    implicit none

    type(Domain), intent(inout) :: Tdomain
    integer, intent(in) :: n,rank
    real, intent(in) :: dt, bega
    integer :: nf,ngllx,nglly,ngllz,ngll1,ngll2,i,j,nnf,orient_f

ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

do nf = 0,5
    nnf = Tdomain%specel(n)%Near_Faces(nf)
    orient_f = Tdomain%specel(n)%Orient_Faces(nf)
    ngll1 = Tdomain%sFace(nnf)%ngll1
    ngll2 = Tdomain%sFace(nnf)%ngll2
 ! now we call the general deassemblage routine
    call get_VectProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,   &
                ngll2,rank,Tdomain%sFace(nnf)%Veloc(:,:,0:2)+              &
                dt*(0.5-bega)*Tdomain%sFace(nnf)%Accel(:,:,0:2),           &
                Tdomain%specel(n)%Forces(:,:,:,0:2)) 
enddo

return

end subroutine get_PMLprediction_f2el
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine get_PMLprediction_f2el_fl(Tdomain,n,bega,dt,rank)

    use sdomain
    implicit none

    type(Domain), intent(inout) :: Tdomain
    integer, intent(in) :: n,rank
    real, intent(in) :: dt, bega
    integer :: nf,ngllx,nglly,ngllz,ngll1,ngll2,i,j,nnf,orient_f

ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

do nf = 0,5
    nnf = Tdomain%specel(n)%Near_Faces(nf)
    orient_f = Tdomain%specel(n)%Orient_Faces(nf)
    ngll1 = Tdomain%sFace(nnf)%ngll1
    ngll2 = Tdomain%sFace(nnf)%ngll2
 ! now we call the general deassemblage routine
    call get_ScalarProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,   &
         Tdomain%sFace(nnf)%VelPhi(:,:)+dt*(0.5-bega)*Tdomain%sFace(nnf)%AccelPhi(:,:), &
         Tdomain%specel(n)%ForcesFl(:,:,:)) 
enddo

return

end subroutine get_PMLprediction_f2el_fl

