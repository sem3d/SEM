!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine StoF_coupling(Tdomain,ngll_sol, ngll_flu, SF_ngll, SF_IGlobSol, SF_IGlobFlu, Veloc, BtN, ForcesFl)
    ! from solid to fluid: velocity (dot) normal
    use sdomain
    use scomm
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: ngll_sol, ngll_flu, SF_ngll
    integer, intent(in), dimension(0:SF_ngll-1) :: SF_IGlobSol, SF_IGlobFlu
    real, intent(in), dimension(0:ngll_sol-1,0:2) :: Veloc
    real, intent(in), dimension(0:SF_ngll-1,0:2) :: BtN
    real, intent(inout), dimension(0:ngll_flu-1) :: ForcesFl
    integer  :: i,j,k,n,idx
    real, dimension(0:SF_ngll-1) :: vn

    vn(:) = 0d0
    do i = 0,SF_ngll-1
        if (SF_IGlobSol(i) < 0) cycle ! solid face not on this proc
        do j = 0,2
            vn(i) = vn(i) + (BtN(i,j) * Veloc(SF_IGlobSol(i),j))
        enddo
    enddo

    if (Tdomain%Comm_SolFlu%ncomm > 0)then
        ! Give
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            k = 0
            do i=0,Tdomain%Comm_SolFlu%Data(n)%ndata-1
                idx = Tdomain%Comm_SolFlu%Data(n)%IGiveS(i)
                Tdomain%Comm_SolFlu%Data(n)%Give(k) = vn(idx)
                k=k+1
            end do
            Tdomain%Comm_SolFlu%Data(n)%nsend = k
        enddo

        ! Exchange interproc of vn
        call exchange_sem_var(Tdomain, 301, Tdomain%Comm_SolFlu)

        ! Take
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            k = 0
            do i=0,Tdomain%Comm_SolFlu%Data(n)%ndata-1
                idx = Tdomain%Comm_SolFlu%Data(n)%ITakeS(i)
                vn(idx) = vn(idx) + Tdomain%Comm_SolFlu%Data(n)%Take(k)
                k=k+1
            enddo
        enddo
    endif

    do i = 0,SF_ngll-1
        if (SF_IGlobFlu(i) < 0) cycle ! fluid face not on this proc
        ForcesFl(SF_IGlobFlu(i)) = ForcesFl(SF_IGlobFlu(i)) + vn(i)
    enddo

end subroutine StoF_coupling
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine FtoS_coupling(Tdomain,ngll_sol, ngll_flu, SF_ngll, SF_IGlobSol, SF_IGlobFlu, BtN, Save_forces, Save_depla, VelPhi, Forces, Depla)
    ! from fluid to solid: normal times pressure (= -rho . VelPhi)
    use sdomain
    use scomm
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: ngll_sol, ngll_flu, SF_ngll
    integer, intent(in), dimension(0:SF_ngll-1) :: SF_IGlobSol, SF_IGlobFlu
    real, intent(in), dimension(0:SF_ngll-1,0:2) :: Save_forces, Save_depla, BtN
    real, intent(in), dimension(0:ngll_flu-1) :: VelPhi
    real, intent(inout), dimension(0:ngll_sol-1,0:2) :: Forces, Depla

    integer  :: i,j,k,n,idx
    real, dimension(0:SF_ngll-1) :: pn

    pn(:) = 0d0
    do i = 0,SF_ngll-1
        do j = 0,2
            pn(i) = pn(i) + (BtN(i,j) * VelPhi(SF_IGlobFlu(i)))
        enddo
    enddo

    if (Tdomain%Comm_SolFlu%ncomm > 0)then
        ! Give
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            k = 0
            do i=0,Tdomain%Comm_SolFlu%Data(n)%ndata-1
                idx = Tdomain%Comm_SolFlu%Data(n)%IGiveS(i)
                Tdomain%Comm_SolFlu%Data(n)%Give(k) = pn(idx)
                k=k+1
            end do
            Tdomain%Comm_SolFlu%Data(n)%nsend = k
        enddo

        ! Exchange interproc of vn
        call exchange_sem_var(Tdomain, 302, Tdomain%Comm_SolFlu)

        ! Take
        do n = 0,Tdomain%Comm_SolFlu%ncomm-1
            k = 0
            do i=0,Tdomain%Comm_SolFlu%Data(n)%ndata-1
                idx = Tdomain%Comm_SolFlu%Data(n)%ITakeS(i)
                pn(idx) = pn(idx) + Tdomain%Comm_SolFlu%Data(n)%Take(k)
                k=k+1
            enddo
        enddo
    endif

    do j = 0,2
        do i = 0,SF_ngll-1
            Forces(SF_IGlobSol(i),j) = Save_forces(i,j) + pn(i)
            Depla(SF_IGlobSol(i),j) = Save_depla(i,j)
        enddo
    enddo

end subroutine FtoS_coupling
#if 0
!!TODO
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine Comm_Forces_Complete_StoF(n,Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: n
    integer  :: ngllSF,i,j,k,nf,ne,nv


    ngllSF = 0
    ! faces
    do i = 0,Tdomain%sComm(n)%SF_nf_shared-1
        nf = Tdomain%sComm(n)%SF_faces_shared(i)
        do j = 1,Tdomain%SF%SF_Face(nf)%ngll2-2
            do k = 1,Tdomain%SF%SF_Face(nf)%ngll1-2
                Tdomain%sComm(n)%GiveForcesSF_StoF(ngllSF) = Tdomain%SF%SF_Face(nf)%Vn(k,j)
                ngllSF = ngllSF + 1
            enddo
        enddo
    enddo
    ! edges
    do i = 0,Tdomain%sComm(n)%SF_ne_shared-1
        ne = Tdomain%sComm(n)%SF_edges_shared(i)
        do j = 1,Tdomain%SF%SF_Edge(ne)%ngll-2
            Tdomain%sComm(n)%GiveForcesSF_StoF(ngllSF) = Tdomain%SF%SF_Edge(ne)%Vn(j)
            ngllSF = ngllSF + 1
        enddo
    enddo
    ! vertices
    do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
        nv =  Tdomain%sComm(n)%SF_vertices_shared(i)
        Tdomain%sComm(n)%GiveForcesSF_StoF(ngllSF) = Tdomain%SF%SF_Vertex(nv)%Vn
        ngllSF = ngllSF + 1
    enddo

    if(ngllSF /= Tdomain%sComm(n)%ngllSF) stop "Bad counting of SF nodes."

    return
end subroutine Comm_Forces_Complete_StoF
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine Comm_Forces_Complete_StoF_PML(n,Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: n
    integer  :: ngllSF_PML,i,j,k,nf,ne,nv


    ngllSF_PML = 0
    ! faces
    do i = 0,Tdomain%sComm(n)%SF_nf_shared-1
        nf = Tdomain%sComm(n)%SF_faces_shared(i)
        if(Tdomain%SF%SF_Face(nf)%PML)then
            do j = 1,Tdomain%SF%SF_Face(nf)%ngll2-2
                do k = 1,Tdomain%SF%SF_Face(nf)%ngll1-2
                    Tdomain%sComm(n)%GiveForcesSF_StoF_PML(ngllSF_PML,1) = Tdomain%SF%SF_Face(nf)%Vn1(k,j)
                    Tdomain%sComm(n)%GiveForcesSF_StoF_PML(ngllSF_PML,2) = Tdomain%SF%SF_Face(nf)%Vn2(k,j)
                    Tdomain%sComm(n)%GiveForcesSF_StoF_PML(ngllSF_PML,3) = Tdomain%SF%SF_Face(nf)%Vn3(k,j)
                    ngllSF_PML = ngllSF_PML + 1
                enddo
            enddo
        end if
    enddo
    ! edges
    do i = 0,Tdomain%sComm(n)%SF_ne_shared-1
        ne = Tdomain%sComm(n)%SF_edges_shared(i)
        if(Tdomain%SF%SF_Edge(ne)%PML)then
            do j = 1,Tdomain%SF%SF_Edge(ne)%ngll-2
                Tdomain%sComm(n)%GiveForcesSF_StoF_PML(ngllSF_PML,1) = Tdomain%SF%SF_Edge(ne)%Vn1(j)
                Tdomain%sComm(n)%GiveForcesSF_StoF_PML(ngllSF_PML,2) = Tdomain%SF%SF_Edge(ne)%Vn2(j)
                Tdomain%sComm(n)%GiveForcesSF_StoF_PML(ngllSF_PML,3) = Tdomain%SF%SF_Edge(ne)%Vn3(j)
                ngllSF_PML = ngllSF_PML + 1
            enddo
        end if
    enddo
    ! vertices
    do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
        nv =  Tdomain%sComm(n)%SF_vertices_shared(i)
        if(Tdomain%SF%SF_Vertex(nv)%PML)then
            Tdomain%sComm(n)%GiveForcesSF_StoF_PML(ngllSF_PML,1) = Tdomain%SF%SF_Vertex(nv)%Vn1
            Tdomain%sComm(n)%GiveForcesSF_StoF_PML(ngllSF_PML,2) = Tdomain%SF%SF_Vertex(nv)%Vn2
            Tdomain%sComm(n)%GiveForcesSF_StoF_PML(ngllSF_PML,3) = Tdomain%SF%SF_Vertex(nv)%Vn3
            ngllSF_PML = ngllSF_PML + 1
        end if
    enddo

    if(ngllSF_PML /= Tdomain%sComm(n)%ngllSF_PML) stop "Bad counting of SF_PML nodes."

    return
end subroutine Comm_Forces_Complete_StoF_PML
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine Comm_Forces_Complete_FtoS(n,Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: n
    integer  :: ngllSF,i,j,k,nf,ne,nv

    ngllSF = 0
    ! faces
    do i = 0,Tdomain%sComm(n)%SF_nf_shared-1
        nf = Tdomain%sComm(n)%SF_faces_shared(i)
        do j = 1,Tdomain%SF%SF_Face(nf)%ngll2-2
            do k = 1,Tdomain%SF%SF_Face(nf)%ngll1-2
                Tdomain%sComm(n)%GiveForcesSF_FtoS(ngllSF,0:2) = Tdomain%SF%SF_Face(nf)%pn(k,j,0:2)
                ngllSF = ngllSF + 1
            enddo
        enddo
    enddo
    ! edges
    do i = 0,Tdomain%sComm(n)%SF_ne_shared-1
        ne = Tdomain%sComm(n)%SF_edges_shared(i)
        do j = 1,Tdomain%SF%SF_Edge(ne)%ngll-2
            Tdomain%sComm(n)%GiveForcesSF_FtoS(ngllSF,0:2) = Tdomain%SF%SF_Edge(ne)%pn(j,0:2)
            ngllSF = ngllSF + 1
        enddo
    enddo
    ! vertices
    do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
        nv =  Tdomain%sComm(n)%SF_vertices_shared(i)
        Tdomain%sComm(n)%GiveForcesSF_FtoS(ngllSF,0:2) = Tdomain%SF%SF_Vertex(nv)%pn(0:2)
        ngllSF = ngllSF + 1
    enddo

    if(ngllSF /= Tdomain%sComm(n)%ngllSF) stop "Bad counting of SF nodes."

    return
end subroutine Comm_Forces_Complete_FtoS
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine Comm_Forces_Complete_FtoS_PML(n,Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: n
    integer  :: ngllSF_PML,i,j,k,nf,ne,nv


    ngllSF_PML = 0
    ! faces
    do i = 0,Tdomain%sComm(n)%SF_nf_shared-1
        nf = Tdomain%sComm(n)%SF_faces_shared(i)
        if(Tdomain%SF%SF_Face(nf)%PML)then
            do j = 1,Tdomain%SF%SF_Face(nf)%ngll2-2
                do k = 1,Tdomain%SF%SF_Face(nf)%ngll1-2
                    Tdomain%sComm(n)%GiveForcesSF_FtoS_PML(ngllSF_PML,1,0:2) = Tdomain%SF%SF_Face(nf)%pn1(k,j,0:2)
                    Tdomain%sComm(n)%GiveForcesSF_FtoS_PML(ngllSF_PML,2,0:2) = Tdomain%SF%SF_Face(nf)%pn2(k,j,0:2)
                    Tdomain%sComm(n)%GiveForcesSF_FtoS_PML(ngllSF_PML,3,0:2) = Tdomain%SF%SF_Face(nf)%pn3(k,j,0:2)
                    ngllSF_PML = ngllSF_PML + 1
                enddo
            enddo
        end if
    enddo
    ! edges
    do i = 0,Tdomain%sComm(n)%SF_ne_shared-1
        ne = Tdomain%sComm(n)%SF_edges_shared(i)
        if(Tdomain%SF%SF_Edge(ne)%PML)then
            do j = 1,Tdomain%SF%SF_Edge(ne)%ngll-2
                Tdomain%sComm(n)%GiveForcesSF_FtoS_PML(ngllSF_PML,1,0:2) = Tdomain%SF%SF_Edge(ne)%pn1(j,0:2)
                Tdomain%sComm(n)%GiveForcesSF_FtoS_PML(ngllSF_PML,2,0:2) = Tdomain%SF%SF_Edge(ne)%pn2(j,0:2)
                Tdomain%sComm(n)%GiveForcesSF_FtoS_PML(ngllSF_PML,3,0:2) = Tdomain%SF%SF_Edge(ne)%pn3(j,0:2)
                ngllSF_PML = ngllSF_PML + 1
            enddo
        end if
    enddo
    ! vertices
    do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
        nv =  Tdomain%sComm(n)%SF_vertices_shared(i)
        if(Tdomain%SF%SF_Vertex(nv)%PML)then
            Tdomain%sComm(n)%GiveForcesSF_FtoS_PML(ngllSF_PML,1,0:2) = Tdomain%SF%SF_Vertex(nv)%pn1(0:2)
            Tdomain%sComm(n)%GiveForcesSF_FtoS_PML(ngllSF_PML,2,0:2) = Tdomain%SF%SF_Vertex(nv)%pn2(0:2)
            Tdomain%sComm(n)%GiveForcesSF_FtoS_PML(ngllSF_PML,3,0:2) = Tdomain%SF%SF_Vertex(nv)%pn3(0:2)
            ngllSF_PML = ngllSF_PML + 1
        end if
    enddo

    if(ngllSF_PML /= Tdomain%sComm(n)%ngllSF_PML) stop "Bad counting of SF_PML nodes."

    return
end subroutine Comm_Forces_Complete_FtoS_PML
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine Comm_Forces_FaceSF_StoF(Tdomain,n,ngllSF,ngllSF_PML)
    use sdomain
    use scommutils
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF,ngllSF_PML
    integer :: ngll1,ngll2,nf,nnf,orient_f


    do nf = 0,Tdomain%sComm(n)%SF_nf_shared-1
        nnf = Tdomain%sComm(n)%SF_faces_shared(nf)
        ngll1 = Tdomain%SF%SF_Face(nnf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nnf)%ngll2
        orient_f = Tdomain%SF%SF_Face(nnf)%Orient_Face
        call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                          &
            Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF:ngllSF+(ngll1-2)*(ngll2-2)-1),  &
            Tdomain%SF%SF_Face(nnf)%Vn(1:ngll1-2,1:ngll2-2))
        ngllSF = ngllSF+(ngll1-2)*(ngll2-2)
        if(Tdomain%SF%SF_Face(nnf)%PML)then
            call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                Tdomain%sComm(n)%TakeForcesSF_StoF_PML(ngllSF_PML:ngllSF_PML+(ngll1-2)*(ngll2-2)-1,1),  &
                Tdomain%SF%SF_Face(nnf)%Vn1(1:ngll1-2,1:ngll2-2))
            call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                Tdomain%sComm(n)%TakeForcesSF_StoF_PML(ngllSF_PML:ngllSF_PML+(ngll1-2)*(ngll2-2)-1,2),  &
                Tdomain%SF%SF_Face(nnf)%Vn2(1:ngll1-2,1:ngll2-2))
            call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                Tdomain%sComm(n)%TakeForcesSF_StoF_PML(ngllSF_PML:ngllSF_PML+(ngll1-2)*(ngll2-2)-1,3),  &
                Tdomain%SF%SF_Face(nnf)%Vn3(1:ngll1-2,1:ngll2-2))
            ngllSF_PML = ngllSF_PML+(ngll1-2)*(ngll2-2)
        end if

    end do

    return

end subroutine Comm_Forces_FaceSF_StoF
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine Comm_Forces_FaceSF_FtoS(Tdomain,n,ngllSF,ngllSF_PML)
    use sdomain
    use scommutils
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF,ngllSF_PML
    integer :: ngll1,ngll2,nf,nnf,orient_f


    do nf = 0,Tdomain%sComm(n)%SF_nf_shared-1
        nnf = Tdomain%sComm(n)%SF_faces_shared(nf)
        ngll1 = Tdomain%SF%SF_Face(nnf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nnf)%ngll2
        orient_f = Tdomain%SF%SF_Face(nnf)%Orient_Face
        call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                          &
            Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF:ngllSF+(ngll1-2)*(ngll2-2)-1,0:2),  &
            Tdomain%SF%SF_Face(nnf)%pn(1:ngll1-2,1:ngll2-2,0:2))
        ngllSF = ngllSF+(ngll1-2)*(ngll2-2)
        if(Tdomain%SF%SF_Face(nnf)%PML)then
            call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                                   &
                Tdomain%sComm(n)%TakeForcesSF_FtoS_PML(ngllSF_PML:ngllSF_PML+(ngll1-2)*(ngll2-2)-1,1,0:2),  &
                Tdomain%SF%SF_Face(nnf)%pn1(1:ngll1-2,1:ngll2-2,0:2))
            call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                                   &
                Tdomain%sComm(n)%TakeForcesSF_FtoS_PML(ngllSF_PML:ngllSF_PML+(ngll1-2)*(ngll2-2)-1,2,0:2),  &
                Tdomain%SF%SF_Face(nnf)%pn2(1:ngll1-2,1:ngll2-2,0:2))
            call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                                   &
                Tdomain%sComm(n)%TakeForcesSF_FtoS_PML(ngllSF_PML:ngllSF_PML+(ngll1-2)*(ngll2-2)-1,3,0:2),  &
                Tdomain%SF%SF_Face(nnf)%pn3(1:ngll1-2,1:ngll2-2,0:2))
            ngllSF_PML = ngllSF_PML+(ngll1-2)*(ngll2-2)
        end if

    end do

    return

end subroutine Comm_Forces_FaceSF_FtoS
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine Comm_Forces_EdgeSF_StoF(Tdomain,n,ngllSF,ngllSF_PML)
    use sdomain
    use scommutils
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF,ngllSF_PML
    integer :: ngll1,ne,nne,orient_e

    do ne = 0,Tdomain%sComm(n)%SF_ne_shared-1
        nne = Tdomain%sComm(n)%SF_edges_shared(ne)
        ngll1 = Tdomain%SF%SF_Edge(nne)%ngll
        orient_e = Tdomain%sComm(n)%SF_mapping_edges_shared(ne)
        call Comm_Edge_ScalarProperty(ngll1,orient_e,                    &
            Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF:ngllSF+ngll1-3),  &
            Tdomain%SF%SF_Edge(nne)%Vn(:))
        ngllSF = ngllSF+ngll1-2
        if(Tdomain%SF%SF_Edge(nne)%PML)then
            call Comm_Edge_ScalarProperty(ngll1,orient_e,                    &
                Tdomain%sComm(n)%TakeForcesSF_StoF_PML(ngllSF_PML:ngllSF_PML+ngll1-3,1),  &
                Tdomain%SF%SF_Edge(nne)%Vn1(:))
            call Comm_Edge_ScalarProperty(ngll1,orient_e,                    &
                Tdomain%sComm(n)%TakeForcesSF_StoF_PML(ngllSF_PML:ngllSF_PML+ngll1-3,2),  &
                Tdomain%SF%SF_Edge(nne)%Vn2(:))
            call Comm_Edge_ScalarProperty(ngll1,orient_e,                    &
                Tdomain%sComm(n)%TakeForcesSF_StoF_PML(ngllSF_PML:ngllSF_PML+ngll1-3,3),  &
                Tdomain%SF%SF_Edge(nne)%Vn3(:))
            ngllSF_PML = ngllSF_PML+ngll1-2
        end if
    end do

    return
end subroutine Comm_Forces_EdgeSF_StoF
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
subroutine Comm_Forces_EdgeSF_FtoS(Tdomain,n,ngllSF,ngllSF_PML)
    use sdomain
    use scommutils
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF,ngllSF_PML
    integer :: ngll1,ne,nne,orient_e

    do ne = 0,Tdomain%sComm(n)%SF_ne_shared-1
        nne = Tdomain%sComm(n)%SF_edges_shared(ne)
        ngll1 = Tdomain%SF%SF_Edge(nne)%ngll
        orient_e = Tdomain%sComm(n)%SF_mapping_edges_shared(ne)
        call Comm_Edge_VectorProperty(ngll1,orient_e,                    &
            Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF:ngllSF+ngll1-3,0:2),  &
            Tdomain%SF%SF_Edge(nne)%pn(:,0:2))
        ngllSF = ngllSF+ngll1-2
        if(Tdomain%SF%SF_Edge(nne)%PML)then
            call Comm_Edge_VectorProperty(ngll1,orient_e,                    &
                Tdomain%sComm(n)%TakeForcesSF_FtoS_PML(ngllSF_PML:ngllSF_PML+ngll1-3,1,0:2),  &
                Tdomain%SF%SF_Edge(nne)%pn1(:,0:2))
            call Comm_Edge_VectorProperty(ngll1,orient_e,                    &
                Tdomain%sComm(n)%TakeForcesSF_FtoS_PML(ngllSF_PML:ngllSF_PML+ngll1-3,2,0:2),  &
                Tdomain%SF%SF_Edge(nne)%pn2(:,0:2))
            call Comm_Edge_VectorProperty(ngll1,orient_e,                    &
                Tdomain%sComm(n)%TakeForcesSF_FtoS_PML(ngllSF_PML:ngllSF_PML+ngll1-3,3,0:2),  &
                Tdomain%SF%SF_Edge(nne)%pn3(:,0:2))
            ngllSF_PML = ngllSF_PML+ngll1-2
        end if
    end do

    return

end subroutine Comm_Forces_EdgeSF_FtoS
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine Comm_Forces_VertexSF_StoF(Tdomain,n,ngllSF,ngllSF_PML)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF,ngllSF_PML
    integer :: i,nv


    do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
        nv = Tdomain%sComm(n)%SF_vertices_shared(i)
        Tdomain%SF%SF_Vertex(nv)%Vn = Tdomain%SF%SF_Vertex(nv)%Vn + Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
        ngllSF = ngllSF + 1
        if(Tdomain%SF%SF_Vertex(nv)%PML)then
            Tdomain%SF%SF_Vertex(nv)%Vn1 = Tdomain%SF%SF_Vertex(nv)%Vn1 +   &
                Tdomain%sComm(n)%TakeForcesSF_StoF_PML(ngllSF_PML,1)
            Tdomain%SF%SF_Vertex(nv)%Vn2 = Tdomain%SF%SF_Vertex(nv)%Vn2 +   &
                Tdomain%sComm(n)%TakeForcesSF_StoF_PML(ngllSF_PML,2)
            Tdomain%SF%SF_Vertex(nv)%Vn3 = Tdomain%SF%SF_Vertex(nv)%Vn3 +   &
                Tdomain%sComm(n)%TakeForcesSF_StoF_PML(ngllSF_PML,3)
            ngllSF_PML = ngllSF_PML + 1
        end if
    enddo

    return
end subroutine Comm_Forces_VertexSF_StoF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine Comm_Forces_VertexSF_FtoS(Tdomain,n,ngllSF,ngllSF_PML)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF,ngllSF_PML
    integer :: i,nv


    do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
        nv =  Tdomain%sComm(n)%SF_vertices_shared(i)
        Tdomain%SF%SF_Vertex(nv)%pn(0:2) = Tdomain%SF%SF_Vertex(nv)%pn(0:2) + Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
        ngllSF = ngllSF + 1
        if(Tdomain%SF%SF_Vertex(nv)%PML)then
            Tdomain%SF%SF_Vertex(nv)%pn1(0:2) = Tdomain%SF%SF_Vertex(nv)%pn1(0:2) +   &
                Tdomain%sComm(n)%TakeForcesSF_FtoS_PML(ngllSF_PML,1,0:2)
            Tdomain%SF%SF_Vertex(nv)%pn2(0:2) = Tdomain%SF%SF_Vertex(nv)%pn2(0:2) +   &
                Tdomain%sComm(n)%TakeForcesSF_FtoS_PML(ngllSF_PML,2,0:2)
            Tdomain%SF%SF_Vertex(nv)%pn3(0:2) = Tdomain%SF%SF_Vertex(nv)%pn3(0:2) +   &
                Tdomain%sComm(n)%TakeForcesSF_FtoS_PML(ngllSF_PML,3,0:2)
            ngllSF_PML = ngllSF_PML + 1
        end if
    enddo

    return
end subroutine Comm_Forces_VertexSF_FtoS
#endif

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine SF_solid_values_saving(ngll_sol, SF_ngll, SF_IGlobSol, Forces, Depla, Save_forces, Save_depla)
    ! saves values of fields Forces, Displ on the solid side of a SF object -
    !    for they are flushed in the 1st correction phase
    implicit none

    integer, intent(in) :: ngll_sol, SF_ngll
    integer, dimension(0:SF_ngll-1), intent(in) :: SF_IGlobSol
    real, dimension(0:ngll_sol-1,0:2), intent(in) :: Forces, Depla
    real, dimension(0:SF_ngll-1,0:2), intent(inout) :: Save_forces, Save_depla
    integer  :: i,idx
    
    do i = 0, SF_ngll-1
        idx = SF_IGlobSol(i)
        if (idx /= -1) then
            Save_forces(i,:) = Forces(SF_IGlobSol(i),:)
            Save_depla(i,:) = Depla(SF_IGlobSol(i),:)
        else
            Save_forces(i,:) = 0.
            Save_depla(i,:) = 0.
        endif
    enddo

    return
end subroutine SF_solid_values_saving
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine Newmark_recorrect_solid(ngll_sol, SF_ngll, dt, SF_IGlobSol, MassMatSol, Forces, Veloc, Depla)
    implicit none

    integer, intent(in) :: ngll_sol, SF_ngll
    real, intent(in) :: dt
    integer, intent(in), dimension(0:SF_ngll-1) :: SF_IGlobSol
    real, intent(in), dimension(0:ngll_sol-1) :: MassMatSol
    real, intent(inout), dimension(0:ngll_sol-1,0:2) :: Forces, Veloc, Depla
    integer  :: i,j,n
    
    do j = 0,2
        do n = 0,SF_ngll-1
            i = SF_IGlobSol(n)
            Forces(i,j) = Forces(i,j) * MassMatSol(i)
            Veloc(i,j) = Veloc(i,j) + dt * Forces(i,j)
            Depla(i,j) = Depla(i,j) + dt * Veloc(i,j)
        enddo
    enddo

end subroutine Newmark_recorrect_solid
