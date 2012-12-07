subroutine get_Mass_Elem2Edge(Tdomain,n,rank)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n,rank
    integer :: ngllx,nglly,ngllz,ngll,i,j,ne,nne,orient_e
 
ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

do ne = 0,11
    nne = Tdomain%specel(n)%Near_Edges(ne)
    orient_e = Tdomain%specel(n)%Orient_Edges(ne)
    ngll = Tdomain%sEdge(nne)%ngll
  ! now we call the general assemblage routine
    call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,  &
               Tdomain%sEdge(nne)%MassMat(:),Tdomain%specel(n)%MassMat(:,:,:))
    if(Tdomain%sEdge(nne)%PML)then
        call get_VectProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,  &
                  Tdomain%sEdge(nne)%DumpMass(:,0:2),Tdomain%specel(n)%DumpMass(:,:,:,0:2)) 
      if(Tdomain%sEdge(nne)%FPML)then
        call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,  &
                  Tdomain%sEdge(nne)%Ivx(:),Tdomain%specel(n)%Ivx(:,:,:)) 
        call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,  &
                  Tdomain%sEdge(nne)%Ivy(:),Tdomain%specel(n)%Ivy(:,:,:)) 
        call get_ScalarProperty_Elem2edge(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,  &
                  Tdomain%sEdge(nne)%Ivz(:),Tdomain%specel(n)%Ivz(:,:,:)) 
      end if
    end if

enddo

return

end subroutine get_Mass_Elem2Edge
