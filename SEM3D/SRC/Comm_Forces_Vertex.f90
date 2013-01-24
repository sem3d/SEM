subroutine Comm_Forces_Vertex(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngll,ngllPML,ngll_F,ngllPML_F
    integer :: i,nv


    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        if(Tdomain%sVertex(nv)%solid)then   ! solid part
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
            ngll = ngll + 1
            if (Tdomain%sVertex(nv)%PML) then
                Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                ngllPML = ngllPML + 1
            endif
        else   ! fluid part
            Tdomain%sVertex(nv)%ForcesFl = Tdomain%sVertex(nv)%ForcesFl + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
            ngll_F = ngll_F + 1
            if (Tdomain%sVertex(nv)%PML) then
                Tdomain%sVertex(nv)%ForcesFl1 = Tdomain%sVertex(nv)%ForcesFl1 + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
                Tdomain%sVertex(nv)%ForcesFl2 = Tdomain%sVertex(nv)%ForcesFl2 + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
                Tdomain%sVertex(nv)%ForcesFl3 = Tdomain%sVertex(nv)%ForcesFl3 + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
                ngllPML_F = ngllPML_F + 1
            endif

        end if
    enddo


    return
end subroutine Comm_Forces_Vertex
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
