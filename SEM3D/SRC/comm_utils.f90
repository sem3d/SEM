!>
!! \file Comm_Forces_Edge.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<
module scommutils

contains

subroutine Comm_Vertex_ScalarProperty(exch_val,prop_vertex)
    ! from vector of exchanged interproc values, to vertices - scalar case
    implicit none

    real, intent(in)  :: exch_val
    real, intent(inout) :: prop_vertex

    prop_vertex = prop_vertex + exch_val

    return

end subroutine Comm_Vertex_ScalarProperty

subroutine Comm_Vertex_VectorProperty(exch_val,prop_vertex)
    ! from vector of exchanged interproc values, to vertices - scalar case
    implicit none

    real, dimension(0:2), intent(in)  :: exch_val
    real, dimension(0:2), intent(inout) :: prop_vertex

    prop_vertex(:) = prop_vertex(:) + exch_val(:)

    return
end subroutine Comm_Vertex_VectorProperty

subroutine Comm_Edge_ScalarProperty(ngll1,orient_e,exch_val,prop_edge)
    ! from vector of exchanged interproc values, to edges - scalar case
    implicit none

    integer, intent(in)  :: ngll1,orient_e
    real, intent(in), dimension(ngll1-2)  :: exch_val
    real, dimension(1:ngll1-2), intent(inout) :: prop_edge
    integer  :: j

    select case(orient_e)
    case(0)
        do j = 1,ngll1-2
            prop_edge(j) = prop_edge(j) + exch_val(j)
        end do
    case(1)
        do j = 1,ngll1-2
            prop_edge(ngll1-1-j) = prop_edge(ngll1-1-j) + exch_val(j)
        end do
    end select

    return
end subroutine Comm_Edge_ScalarProperty

subroutine Comm_Edge_VectorProperty(ngll1,orient_e,exch_val,prop_edge)
    ! from vector of exchanged interproc values, to edges - vectorial case
    implicit none

    integer, intent(in)  :: ngll1,orient_e
    real, intent(in), dimension(ngll1-2,0:2)  :: exch_val
    real, dimension(1:ngll1-2,0:2), intent(inout) :: prop_edge
    integer  :: j

    select case(orient_e)
    case(0)
        do j = 1,ngll1-2
            prop_edge(j,0:2) = prop_edge(j,0:2) + exch_val(j,0:2)
        end do
    case(1)
        do j = 1,ngll1-2
            prop_edge(ngll1-1-j,0:2) = prop_edge(ngll1-1-j,0:2) + exch_val(j,0:2)
        end do
    end select
end subroutine Comm_Edge_VectorProperty

subroutine Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,exch_val,prop_face)
    ! from vector of exchanged interproc values, to faces - scalar case
    implicit none

    integer, intent(in)  :: ngll1,ngll2,orient_f
    real, intent(in), dimension((ngll1-2)*(ngll2-2))  :: exch_val
    real, dimension(1:ngll1-2,1:ngll2-2), intent(inout) :: prop_face
    integer  :: j,k,ngll

    ngll = 1
    select case(orient_f)
    case(0)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(k,j) = prop_face(k,j) + exch_val(ngll)
                ngll = ngll + 1
            end do
        end do
    case(1)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-k,j) = prop_face(ngll1-1-k,j) + exch_val(ngll)
                ngll = ngll + 1
            end do
        end do
    case(2)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(k,ngll2-1-j) = prop_face(k,ngll2-1-j) + exch_val(ngll)
                ngll = ngll + 1
            end do
        end do
    case(3)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-k,ngll2-1-j) = prop_face(ngll1-1-k,ngll2-1-j) + exch_val(ngll)
                ngll = ngll + 1
            end do
        end do
    case(4)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(j,k) = prop_face(j,k) + exch_val(ngll)
                ngll = ngll + 1
            end do
        end do
    case(5)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-j,k) = prop_face(ngll1-1-j,k) + exch_val(ngll)
                ngll = ngll + 1
            end do
        end do
    case(6)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(j,ngll2-1-k) = prop_face(j,ngll2-1-k) + exch_val(ngll)
                ngll = ngll + 1
            end do
        end do
    case(7)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-j,ngll2-1-k) = prop_face(ngll1-1-j,ngll2-1-k) + exch_val(ngll)
                ngll = ngll + 1
            end do
        end do
    end select

    if(ngll /= (ngll1-2)*(ngll2-2)+1) &
        stop "Pb in Comm_Face_ScalarProperty"

end subroutine Comm_Face_ScalarProperty

subroutine Comm_Face_VectorProperty(ngll1,ngll2,orient_f,exch_val,prop_face)
    ! from vector of exchanged interproc values, to faces - vectorial case
    implicit none

    integer, intent(in)  :: ngll1,ngll2,orient_f
    real, intent(in), dimension(1:(ngll1-2)*(ngll2-2),0:2)  :: exch_val
    real, dimension(1:ngll1-2,1:ngll2-2,0:2), intent(inout) :: prop_face
    integer  :: j,k,ngll

    ngll = 1
    select case(orient_f)
    case(0)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(k,j,0:2) = prop_face(k,j,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(1)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-k,j,0:2) = prop_face(ngll1-1-k,j,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(2)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(k,ngll2-1-j,0:2) = prop_face(k,ngll2-1-j,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(3)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-k,ngll2-1-j,0:2) = prop_face(ngll1-1-k,ngll2-1-j,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(4)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(j,k,0:2) = prop_face(j,k,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(5)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-j,k,0:2) = prop_face(ngll1-1-j,k,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(6)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(j,ngll2-1-k,0:2) = prop_face(j,ngll2-1-k,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    case(7)
        do j = 1,ngll2-2
            do k = 1,ngll1-2
                prop_face(ngll1-1-j,ngll2-1-k,0:2) = prop_face(ngll1-1-j,ngll2-1-k,0:2) + exch_val(ngll,0:2)
                ngll = ngll + 1
            end do
        end do
    end select

    if(ngll /= (ngll1-2)*(ngll2-2)+1) &
        stop "Pb in Comm_Face_VectorProperty"

end subroutine Comm_Face_VectorProperty

#if ! NEW_GLOBAL_METHOD
subroutine Comm_Forces_Edge (Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
    use sdomain

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngll,ngllPML,ngll_F,ngllPML_F
    integer :: ngll1,ne,nne,orient_e

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
                    Tdomain%sEdge(nne)%spml%Forces1(:,:))
                call Comm_Edge_VectorProperty(ngll1,orient_e,                         &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+ngll1-3,2,0:2),  &
                    Tdomain%sEdge(nne)%spml%Forces2(:,:))
                call Comm_Edge_VectorProperty(ngll1,orient_e,                         &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+ngll1-3,3,0:2),  &
                    Tdomain%sEdge(nne)%spml%Forces3(:,:))
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
                    Tdomain%sEdge(nne)%spml%ForcesFl1(:))
                call Comm_Edge_ScalarProperty(ngll1,orient_e,                             &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+ngll1-3,2),  &
                    Tdomain%sEdge(nne)%spml%ForcesFl2(:))
                call Comm_Edge_ScalarProperty(ngll1,orient_e,                             &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+ngll1-3,3),  &
                    Tdomain%sEdge(nne)%spml%ForcesFl3(:))
                ngllPML_F = ngllPML_F+ngll1-2
            end if

        end if
    end do

    return
end subroutine Comm_Forces_Edge
#endif

#if NEW_GLOBAL_METHOD
subroutine Comm_Forces_Assembl(Tdomain,n)
    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer :: ngll, ngll_F, ngllPML, ngllPML_F, i, j, k, nf, nnf
    integer :: ngll1, ngll2, idx, ne, nne, nv

    ngll = 0
    ngll_F = 0
    ngllPML = 0
    ngllPML_F = 0
    
    ! Faces
    do nf = 0,Tdomain%sComm(n)%nb_faces-1
        nnf = Tdomain%sComm(n)%faces(nf)
        ngll1 = Tdomain%sFace(nnf)%ngll1
        ngll2 = Tdomain%sFace(nnf)%ngll2
        if(Tdomain%sFace(nnf)%solid)then   ! solid part
            if(Tdomain%sFace(nnf)%PML)then
                do j = 1,Tdomain%sFace(nnf)%ngll2-2
                    do k = 1,Tdomain%sFace(nnf)%ngll1-2
                        idx = Tdomain%sFace(nnf)%Renum(k,j)
                        Tdomain%champs1%ForcesPML(idx,0:2) = Tdomain%champs1%ForcesPML(idx,0:2) + &
                                                             Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                        Tdomain%champs1%ForcesPML(idx+1,0:2) = Tdomain%champs1%ForcesPML(idx+1,0:2) + &
                                                               Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                        Tdomain%champs1%ForcesPML(idx+2,0:2) = Tdomain%champs1%ForcesPML(idx+2,0:2) + &
                                                               Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            else
                do j = 1,Tdomain%sFace(nnf)%ngll2-2
                    do k = 1,Tdomain%sFace(nnf)%ngll1-2
                        idx = Tdomain%sFace(nnf)%Renum(k,j)
                        Tdomain%champs1%Forces(idx,0:2) = Tdomain%champs1%Forces(idx,0:2) + &
                                                          Tdomain%sComm(n)%TakeForces(ngll,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
            endif
        else  ! fluid part
            if(Tdomain%sFace(nnf)%PML)then
                do j = 1,Tdomain%sFace(nnf)%ngll2-2
                    do k = 1,Tdomain%sFace(nnf)%ngll1-2
                        idx = Tdomain%sFace(nnf)%Renum(k,j)
!TODO
                        stop "comm_utils Comm_Forces_Assembl : PML fluide"
!                         Tdomain%champs1%ForcesFLPML(idx) = Tdomain%champs1%ForcesFLPML(idx) + &
!                                                            Tdomain%sComm(n)%TakeForcesPMLFL(ngllPML_F,1)
!                         Tdomain%champs1%ForcesFLPML(idx+1) = Tdomain%champs1%ForcesFLPML(idx+1) + &
!                                                              Tdomain%sComm(n)%TakeForcesPMLFL(ngllPML_F,2)
!                         Tdomain%champs1%ForcesFLPML(idx+2) = Tdomain%champs1%ForcesFLPML(idx+2) + &
!                                                              Tdomain%sComm(n)%TakeForcesPMLFL(ngllPML_F,3)
                        ngllPML_F = ngllPML_F + 1
                    enddo
                enddo
            else
                do j = 1,Tdomain%sFace(nnf)%ngll2-2
                    do k = 1,Tdomain%sFace(nnf)%ngll1-2
                        idx = Tdomain%sFace(nnf)%Renum(k,j)
                        Tdomain%champs1%ForcesFl(idx) = Tdomain%champs1%ForcesFl(idx) + &
                                                        Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
            endif
        endif
    enddo

    ! Edges
    do ne = 0,Tdomain%sComm(n)%nb_edges-1
        nne = Tdomain%sComm(n)%edges(ne)
        ngll1 = Tdomain%sEdge(nne)%ngll
        if(Tdomain%sEdge(nne)%solid)then   ! solid part
            if(Tdomain%sEdge(nne)%PML)then
                do k = 1,Tdomain%sEdge(nne)%ngll-2
                    idx = Tdomain%sEdge(nne)%Renum(k)
                    Tdomain%champs1%ForcesPML(idx,0:2) = Tdomain%champs1%ForcesPML(idx,0:2) + &
                                                         Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                    Tdomain%champs1%ForcesPML(idx+1,0:2) = Tdomain%champs1%ForcesPML(idx+1,0:2) + &
                                                           Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                    Tdomain%champs1%ForcesPML(idx+2,0:2) = Tdomain%champs1%ForcesPML(idx+2,0:2) + &
                                                           Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                    ngllPML = ngllPML + 1
                enddo
            else
                do k = 1,Tdomain%sEdge(nne)%ngll-2
                    idx = Tdomain%sEdge(nne)%Renum(k)
                    Tdomain%champs1%Forces(idx,0:2) = Tdomain%champs1%Forces(idx,0:2) + &
                                                      Tdomain%sComm(n)%TakeForces(ngll,0:2)
                    ngll = ngll + 1
                enddo
            endif
        else  ! fluid part
            if(Tdomain%sEdge(nne)%PML)then
                do k = 1,Tdomain%sEdge(nne)%ngll-2
                    idx = Tdomain%sEdge(nne)%Renum(k)
!TODO
                    stop "comm_utils Comm_Forces_Assembl : PML fluide"
!                     Tdomain%champs1%ForcesFLPML(idx) = Tdomain%champs1%ForcesFLPML(idx) + &
!                                                        Tdomain%sComm(n)%TakeForcesPMLFL(ngllPML_F,1)
!                     Tdomain%champs1%ForcesFLPML(idx+1) = Tdomain%champs1%ForcesFLPML(idx+1) + &
!                                                          Tdomain%sComm(n)%TakeForcesPMLFL(ngllPML_F,2)
!                     Tdomain%champs1%ForcesFLPML(idx+2) = Tdomain%champs1%ForcesFLPML(idx+2) + &
!                                                          Tdomain%sComm(n)%TakeForcesPMLFL(ngllPML_F,3)
                    ngllPML_F = ngllPML_F + 1
                enddo
            else
                do k = 1,Tdomain%sEdge(nne)%ngll-2
                    idx = Tdomain%sEdge(nne)%Renum(k)
                    Tdomain%champs1%ForcesFl(idx) = Tdomain%champs1%ForcesFl(idx) + &
                                                    Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                    ngll_F = ngll_F + 1
                enddo
            endif
        endif
    enddo

    ! Vertices
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        if(Tdomain%sVertex(nv)%solid)then   ! solid part
            if(Tdomain%sVertex(nv)%PML)then
                idx = Tdomain%sVertex(nv)%Renum
                Tdomain%champs1%ForcesPML(idx,0:2) = Tdomain%champs1%ForcesPML(idx,0:2) + &
                                                     Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                Tdomain%champs1%ForcesPML(idx+1,0:2) = Tdomain%champs1%ForcesPML(idx+1,0:2) + &
                                                       Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                Tdomain%champs1%ForcesPML(idx+2,0:2) = Tdomain%champs1%ForcesPML(idx+2,0:2) + &
                                                       Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                ngllPML = ngllPML + 1
            else
                idx = Tdomain%sVertex(nv)%Renum
                Tdomain%champs1%Forces(idx,0:2) = Tdomain%champs1%Forces(idx,0:2) + &
                                                  Tdomain%sComm(n)%TakeForces(ngll,0:2)
                ngll = ngll + 1
            endif
        else  ! fluid part
            if(Tdomain%sVertex(nv)%PML)then
                idx = Tdomain%sVertex(nv)%Renum
!TODO
                stop "comm_utils Comm_Forces_Assembl : PML fluide"
!                 Tdomain%champs1%ForcesFLPML(idx) = Tdomain%champs1%ForcesFLPML(idx) + &
!                                                    Tdomain%sComm(n)%TakeForcesPMLFL(ngllPML_F,1)
!                 Tdomain%champs1%ForcesFLPML(idx+1) = Tdomain%champs1%ForcesFLPML(idx+1) + &
!                                                      Tdomain%sComm(n)%TakeForcesPMLFL(ngllPML_F,2)
!                 Tdomain%champs1%ForcesFLPML(idx+2) = Tdomain%champs1%ForcesFLPML(idx+2) + &
!                                                      Tdomain%sComm(n)%TakeForcesPMLFL(ngllPML_F,3)
                ngllPML_F = ngllPML_F + 1
            else
                idx = Tdomain%sVertex(nv)%Renum
                Tdomain%champs1%ForcesFl(idx) = Tdomain%champs1%ForcesFl(idx) + &
                                                Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                ngll_F = ngll_F + 1
            endif
        endif
    enddo

    return
end subroutine Comm_Forces_Assembl
#endif

#if ! NEW_GLOBAL_METHOD
subroutine Comm_Forces_Face(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngll,ngllPML,ngll_F,ngllPML_F
    integer :: ngll1,ngll2,nf,nnf,orient_f

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
                    Tdomain%sFace(nnf)%spml%Forces1(:,:,:))
                call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                               &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+(ngll1-2)*(ngll2-2)-1,2,0:2),  &
                    Tdomain%sFace(nnf)%spml%Forces2(:,:,:))
                call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                               &
                    Tdomain%sComm(n)%TakeForcesPML(ngllPML:ngllPML+(ngll1-2)*(ngll2-2)-1,3,0:2),  &
                    Tdomain%sFace(nnf)%spml%Forces3(:,:,:))
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
                    Tdomain%sFace(nnf)%spml%ForcesFl1(:,:))
                call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+(ngll1-2)*(ngll2-2)-1,2),  &
                    Tdomain%sFace(nnf)%spml%ForcesFl2(:,:))
                call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                    Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F:ngllPML_F+(ngll1-2)*(ngll2-2)-1,3),  &
                    Tdomain%sFace(nnf)%spml%ForcesFl3(:,:))
                ngllPML_F = ngllPML_F+(ngll1-2)*(ngll2-2)
            end if

        end if
    end do

    return
end subroutine Comm_Forces_Face
#endif

#if ! NEW_GLOBAL_METHOD
subroutine Comm_Forces_Vertex (Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)

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
                Tdomain%sVertex(nv)%spml%Forces1(0:2) = Tdomain%sVertex(nv)%spml%Forces1(0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                Tdomain%sVertex(nv)%spml%Forces2(0:2) = Tdomain%sVertex(nv)%spml%Forces2(0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                Tdomain%sVertex(nv)%spml%Forces3(0:2) = Tdomain%sVertex(nv)%spml%Forces3(0:2) + Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                ngllPML = ngllPML + 1
            endif
        else   ! fluid part
            Tdomain%sVertex(nv)%ForcesFl = Tdomain%sVertex(nv)%ForcesFl + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
            ngll_F = ngll_F + 1
            if (Tdomain%sVertex(nv)%PML) then
                Tdomain%sVertex(nv)%spml%ForcesFl1 = Tdomain%sVertex(nv)%spml%ForcesFl1 + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
                Tdomain%sVertex(nv)%spml%ForcesFl2 = Tdomain%sVertex(nv)%spml%ForcesFl2 + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
                Tdomain%sVertex(nv)%spml%ForcesFl3 = Tdomain%sVertex(nv)%spml%ForcesFl3 + Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
                ngllPML_F = ngllPML_F + 1
            endif

        end if
    enddo

    return
end subroutine Comm_Forces_Vertex
#endif
#if ! NEW_GLOBAL_METHOD
subroutine Comm_Mass_Edge (Tdomain,n,ngll,ngllPML)

    use sdomain

    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: ngll1,i,j,ne


    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        ngll1 = Tdomain%sEdge(ne)%ngll

        if ( Tdomain%sComm(n)%orient_edges(i) == 0 ) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + Tdomain%sComm(n)%Take(ngll)
                ngll = ngll + 1
            enddo
            if (Tdomain%sEdge(ne)%PML) then
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    Tdomain%sEdge(ne)%spml%DumpMass(j,0:2) = Tdomain%sEdge(ne)%spml%DumpMass(j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                    if (Tdomain%any_FPML) then
                        Tdomain%sEdge(ne)%spml%Ivx(j) = Tdomain%sEdge(ne)%spml%Ivx(j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                        Tdomain%sEdge(ne)%spml%Ivy(j) = Tdomain%sEdge(ne)%spml%Ivy(j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                        Tdomain%sEdge(ne)%spml%Ivz(j) = Tdomain%sEdge(ne)%spml%Ivz(j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                    endif
                    ngllPML = ngllPML + 1
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_edges(i) == 1 ) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sEdge(ne)%MassMat(ngll1-1-j) = Tdomain%sEdge(ne)%MassMat(ngll1-1-j) + Tdomain%sComm(n)%Take(ngll)
                ngll = ngll + 1
            enddo
            if (Tdomain%sEdge(ne)%PML) then
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    Tdomain%sEdge(ne)%spml%DumpMass(ngll1-1-j,0:2) = Tdomain%sEdge(ne)%spml%DumpMass(ngll1-1-j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                    if (Tdomain%any_FPML) then
                        Tdomain%sEdge(ne)%spml%Ivx(ngll1-1-j) = Tdomain%sEdge(ne)%spml%Ivx(ngll1-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                        Tdomain%sEdge(ne)%spml%Ivy(ngll1-1-j) = Tdomain%sEdge(ne)%spml%Ivy(ngll1-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                        Tdomain%sEdge(ne)%spml%Ivz(ngll1-1-j) = Tdomain%sEdge(ne)%spml%Ivz(ngll1-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                    endif
                    ngllPML = ngllPML + 1
                enddo
            endif

        else
            print*,'Pb with coherency number for edges in Comm_Mass_Edge', i, Tdomain%sComm(n)%orient_edges(i)
            STOP 1
        endif

    enddo

    return
end subroutine Comm_Mass_Edge
#endif

#if NEW_GLOBAL_METHOD
subroutine Comm_Mass_Face_2 (Tdomain,n,ngll,ngllPML)
    use sdomain
    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: ngll1,ngll2,i,j,k,nf,idx

    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        if (Tdomain%sFace(nf)%solid) then ! solide
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,ngll2-2
                    do k = 1,ngll1-2
                        idx = Tdomain%sFace(nf)%Renum(k,j)
                        Tdomain%DumpMass(idx) = Tdomain%DumpMass(idx) + &
                                                Tdomain%sComm(n)%TakePML(ngllPML,0)
                        Tdomain%DumpMass(idx+1) = Tdomain%DumpMass(idx+1) + &
                                                  Tdomain%sComm(n)%TakePML(ngllPML,1)
                        Tdomain%DumpMass(idx+2) = Tdomain%DumpMass(idx+2) + &
                                                  Tdomain%sComm(n)%TakePML(ngllPML,2)
                        Tdomain%MassMatSolPml(idx) = Tdomain%MassMatSolPml(idx) + &
                                                     Tdomain%sComm(n)%TakePML(ngllPML,3)
                        Tdomain%MassMatSolPml(idx+1) = Tdomain%MassMatSolPml(idx+1) + &
                                                       Tdomain%sComm(n)%TakePML(ngllPML,3)
                        Tdomain%MassMatSolPml(idx+2) = Tdomain%MassMatSolPml(idx+2) + &
                                                       Tdomain%sComm(n)%TakePML(ngllPML,3)
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            else
                do j = 1,ngll2-2
                    do k = 1,ngll1-2
                        idx = Tdomain%sFace(nf)%Renum(k,j)
                        Tdomain%MassMatSol(idx) = Tdomain%MassMatSol(idx) + &
                                                  Tdomain%sComm(n)%Take(ngll)
                        ngll = ngll + 1
                    enddo
                enddo
            endif
        else ! fluide
            if (Tdomain%sFace(nf)%PML) then
                !TODO
                stop "comm_utils.f90 Comm_Mass_Face_2 PML fluide"
            else
                !TODO
                stop "comm_utils.f90 Comm_Mass_Face_2 fluide"
            endif
        endif
    enddo

    return
end subroutine Comm_Mass_Face_2

subroutine Comm_Mass_Edge_2(Tdomain,n,ngll,ngllPML)
    use sdomain
    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: ngll1,i,j,ne,idx

    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        ngll1 = Tdomain%sEdge(ne)%ngll
        if (Tdomain%sEdge(ne)%solid) then
            if (Tdomain%sEdge(ne)%PML) then
                do j = 1,ngll1-2
                    idx = Tdomain%sEdge(ne)%Renum(j)
                    Tdomain%DumpMass(idx) = Tdomain%DumpMass(idx) + &
                                            Tdomain%sComm(n)%TakePML(ngllPML,0)
                    Tdomain%DumpMass(idx+1) = Tdomain%DumpMass(idx+1) + &
                                              Tdomain%sComm(n)%TakePML(ngllPML,1)
                    Tdomain%DumpMass(idx+2) = Tdomain%DumpMass(idx+2) + &
                                              Tdomain%sComm(n)%TakePML(ngllPML,2)
                    Tdomain%MassMatSolPml(idx) = Tdomain%MassMatSolPml(idx) + &
                                                 Tdomain%sComm(n)%TakePML(ngllPML,3)
                    Tdomain%MassMatSolPml(idx+1) = Tdomain%MassMatSolPml(idx+1) + &
                                                   Tdomain%sComm(n)%TakePML(ngllPML,3)
                    Tdomain%MassMatSolPml(idx+2) = Tdomain%MassMatSolPml(idx+2) + &
                                                   Tdomain%sComm(n)%TakePML(ngllPML,3)
                    ngllPML = ngllPML + 1
                enddo
            else
                do j = 1,ngll1-2
                    idx = Tdomain%sEdge(ne)%Renum(j)
                    Tdomain%MassMatSol(idx) = Tdomain%MassMatSol(idx) + &
                                                Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            endif
        else
            if (Tdomain%sEdge(ne)%PML) then
                !TODO
                stop "comm_utils.f90 Comm_Mass_Edge_2 PML fluide"
            else
                !TODO
                stop "comm_utils.f90 Comm_Mass_Edge_2 fluide"
            endif
        endif
    enddo

    return
end subroutine Comm_Mass_Edge_2

subroutine Comm_Mass_Vertex_2(Tdomain,n,ngll,ngllPML)
    use sdomain
    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: i,nv,idx

    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        if (Tdomain%sVertex(nv)%solid) then
            if (Tdomain%sVertex(nv)%PML) then
                idx = Tdomain%sVertex(nv)%Renum
                Tdomain%DumpMass(idx) = Tdomain%DumpMass(idx) + &
                                        Tdomain%sComm(n)%TakePML(ngllPML,0)
                Tdomain%DumpMass(idx+1) = Tdomain%DumpMass(idx+1) + &
                                          Tdomain%sComm(n)%TakePML(ngllPML,1)
                Tdomain%DumpMass(idx+2) = Tdomain%DumpMass(idx+2) + &
                                          Tdomain%sComm(n)%TakePML(ngllPML,2)
                Tdomain%MassMatSolPml(idx) = Tdomain%MassMatSolPml(idx) + &
                                             Tdomain%sComm(n)%TakePML(ngllPML,3)
                Tdomain%MassMatSolPml(idx+1) = Tdomain%MassMatSolPml(idx+1) + &
                                               Tdomain%sComm(n)%TakePML(ngllPML,3)
                Tdomain%MassMatSolPml(idx+2) = Tdomain%MassMatSolPml(idx+2) + &
                                               Tdomain%sComm(n)%TakePML(ngllPML,3)
                ngllPML = ngllPML + 1
            else
                idx = Tdomain%sVertex(nv)%Renum
                Tdomain%MassMatSol(idx) = Tdomain%MassMatSol(idx) + &
                                          Tdomain%sComm(n)%Take(ngll)
                ngll = ngll + 1
            endif
        else
            if (Tdomain%sVertex(nv)%PML) then
                !TODO
                stop "comm_utils.f90 Comm_Mass_Vertex_2 PML fluide"
            else
                !TODO
                stop "comm_utils.f90 Comm_Mass_Vertex_2 fluide"
            endif
        endif
    enddo

    return
end subroutine Comm_Mass_Vertex_2
#endif

#if ! NEW_GLOBAL_METHOD
subroutine Comm_Mass_Face (Tdomain,n,ngll,ngllPML)
    use sdomain
    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: ngll1,ngll2,i,j,k,nf

    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2

        if ( Tdomain%sComm(n)%orient_faces(i) == 0 ) then

            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(k,j) = Tdomain%sFace(nf)%MassMat(k,j) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%spml%DumpMass(k,j,0:2) = Tdomain%sFace(nf)%spml%DumpMass(k,j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(k,j) = Tdomain%sFace(nf)%spml%Ivx(k,j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(k,j) = Tdomain%sFace(nf)%spml%Ivy(k,j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(k,j) = Tdomain%sFace(nf)%spml%Ivz(k,j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 1 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-k,j) = Tdomain%sFace(nf)%MassMat(ngll1-1-k,j) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-k,j,0:2) = Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-k,j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(ngll1-1-k,j) = Tdomain%sFace(nf)%spml%Ivx(ngll1-1-k,j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(ngll1-1-k,j) = Tdomain%sFace(nf)%spml%Ivy(ngll1-1-k,j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(ngll1-1-k,j) = Tdomain%sFace(nf)%spml%Ivz(ngll1-1-k,j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 2 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(k,ngll2-1-j) = Tdomain%sFace(nf)%MassMat(k,ngll2-1-j) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%spml%DumpMass(k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%spml%DumpMass(k,ngll2-1-j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivx(k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivy(k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivz(k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 3 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%MassMat(ngll1-1-k,ngll2-1-j) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-k,ngll2-1-j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivx(ngll1-1-k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivy(ngll1-1-k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivz(ngll1-1-k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 4 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%spml%DumpMass(j,k,0:2) = Tdomain%sFace(nf)%spml%DumpMass(j,k,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(j,k) = Tdomain%sFace(nf)%spml%Ivx(j,k) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(j,k) = Tdomain%sFace(nf)%spml%Ivy(j,k) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(j,k) = Tdomain%sFace(nf)%spml%Ivz(j,k) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 5 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-j,k) = Tdomain%sFace(nf)%MassMat(ngll1-1-j,k) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-j,k,0:2) = Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-j,k,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(ngll1-1-j,k) = Tdomain%sFace(nf)%spml%Ivx(ngll1-1-j,k) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(ngll1-1-j,k) = Tdomain%sFace(nf)%spml%Ivy(ngll1-1-j,k) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(ngll1-1-j,k) = Tdomain%sFace(nf)%spml%Ivz(ngll1-1-j,k) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 6 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(j,ngll2-1-k) = Tdomain%sFace(nf)%MassMat(j,ngll2-1-k) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%spml%DumpMass(j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%spml%DumpMass(j,ngll2-1-k,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivx(j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivy(j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivz(j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 7 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%MassMat(ngll1-1-j,ngll2-1-k) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-j,ngll2-1-k,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivx(ngll1-1-j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivy(ngll1-1-j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivz(ngll1-1-j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else
            print*,'Pb with coherency number for face'

        endif

    enddo

    return
end subroutine Comm_Mass_Face

subroutine Comm_Mass_Vertex (Tdomain,n,ngll,ngllPML)

    use sdomain

    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: i,nv


    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        Tdomain%svertex(nv)%MassMat = Tdomain%svertex(nv)%MassMat + Tdomain%sComm(n)%Take(ngll)
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%PML) then
            Tdomain%svertex(nv)%spml%DumpMass(0:2) = Tdomain%svertex(nv)%spml%DumpMass(0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
            if (Tdomain%any_FPML) then
                Tdomain%sVertex(nv)%spml%Ivx(0) = Tdomain%sVertex(nv)%spml%Ivx(0) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                Tdomain%sVertex(nv)%spml%Ivy(0) = Tdomain%sVertex(nv)%spml%Ivy(0) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                Tdomain%sVertex(nv)%spml%Ivz(0) = Tdomain%sVertex(nv)%spml%Ivz(0) + Tdomain%sComm(n)%TakePML(ngllPML,5)
            endif
            ngllPML = ngllPML + 1
        endif
    enddo

    return
end subroutine Comm_Mass_Vertex
#endif
end module scommutils
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
