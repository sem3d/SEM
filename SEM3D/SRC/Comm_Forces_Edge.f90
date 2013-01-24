subroutine Comm_Forces_Edge(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngll,ngllPML,ngll_F,ngllPML_F
    integer :: ngll1,i,ne,nne,orient_e


    do ne = 0,Tdomain%sComm(n)%nb_edges-1
        nne = Tdomain%sComm(n)%edges(ne)
        ngll1 = Tdomain%sEdge(nne)%ngll
        orient_e = Tdomain%sComm(n)%orient_edges(ne)
        if(Tdomain%sEdge(nne)%solid)then   ! solid part
            call Comm_Edge_VectorProperty(ngll1,orient_e,                &
                Tdomain%sComm(n)%TakeForces(ngll:ngll+ngll1-3,0:2),  &
                Tdomain%sEdge(nne)%Forces(:,:))
            ngll = ngll+ngll1-2
            if(Tdomain%sEdge(nne)%PML)then
                call Comm_Edge_VectorProperty(ngll1,orient_e,                         &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+ngll1-3,1,0:2),  &
                    Tdomain%sEdge(nne)%Forces1(:,:))
                call Comm_Edge_VectorProperty(ngll1,orient_e,                         &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+ngll1-3,2,0:2),  &
                    Tdomain%sEdge(nne)%Forces2(:,:))
                call Comm_Edge_VectorProperty(ngll1,orient_e,                         &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+ngll1-3,3,0:2),  &
                    Tdomain%sEdge(nne)%Forces3(:,:))
                ngllPML = ngllPML+ngll1-2
            end if

        else   ! fluid part
            call Comm_Edge_ScalarProperty(ngll1,orient_e,                    &
                Tdomain%sComm(n)%TakeForcesFl(ngll_F:ngll_F+ngll1-3),  &
                Tdomain%sEdge(nne)%ForcesFl(:))
            ngll_F = ngll_F+ngll1-2
            if(Tdomain%sEdge(nne)%PML)then
                call Comm_Edge_ScalarProperty(ngll1,orient_e,                             &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+ngll1-3,1),  &
                    Tdomain%sEdge(nne)%ForcesFl1(:))
                call Comm_Edge_ScalarProperty(ngll1,orient_e,                             &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+ngll1-3,2),  &
                    Tdomain%sEdge(nne)%ForcesFl2(:))
                call Comm_Edge_ScalarProperty(ngll1,orient_e,                             &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+ngll1-3,3),  &
                    Tdomain%sEdge(nne)%ForcesFl3(:))
                ngllPML_F = ngllPML_F+ngll1-2
            end if

        end if
    end do

    return
end subroutine Comm_Forces_Edge
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
