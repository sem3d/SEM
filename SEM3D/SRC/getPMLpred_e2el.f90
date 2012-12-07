subroutine get_PMLprediction_e2el(Tdomain,n,bega,dt,rank)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n,rank
    real, intent(in)  :: bega,dt
    integer :: i,ne,ngllx,nglly,ngllz,ngll,nne,orient_e
 
ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

do ne = 0,11
    nne = Tdomain%specel(n)%Near_Edges(ne)
    orient_e = Tdomain%specel(n)%Orient_Edges(ne)
    ngll = Tdomain%sEdge(nne)%ngll
 ! now we call the general deassemblage routine
    call get_VectProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,            &
         Tdomain%sEdge(nne)%Veloc(:,0:2)+dt*(0.5-bega)*Tdomain%sEdge(nne)%Accel(:,0:2), &
         Tdomain%specel(n)%Forces(:,:,:,0:2))
end do

return

end subroutine get_PMLprediction_e2el
!-----------------------------------------------------------------
!-----------------------------------------------------------------
subroutine get_PMLprediction_e2el_fl(Tdomain,n,bega,dt,rank)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n,rank
    real, intent(in)  :: bega,dt
    integer :: i,ne,ngllx,nglly,ngllz,ngll,nne,orient_e
 
ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

do ne = 0,11
    nne = Tdomain%specel(n)%Near_Edges(ne)
    orient_e = Tdomain%specel(n)%Orient_Edges(ne)
    ngll = Tdomain%sEdge(nne)%ngll
 ! now we call the general deassemblage routine
    call get_ScalarProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,          &
         Tdomain%sEdge(nne)%VelPhi(:)+dt*(0.5-bega)*Tdomain%sEdge(nne)%AccelPhi(:), &
         Tdomain%specel(n)%ForcesFl(:,:,:))
end do

return

end subroutine get_PMLprediction_e2el_fl
