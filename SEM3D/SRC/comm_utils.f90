!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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


subroutine Comm_Forces_Edge (Tdomain,io_comm,ngll,ngll_F,ngllPML,ngllPML_F)
    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    type(comm), intent(inout) :: io_comm
    integer, intent(inout) :: ngll,ngllPML,ngll_F,ngllPML_F
    integer :: ngll1,ne,nne,orient_e


    do ne = 0,io_comm%nb_edges-1
        nne = io_comm%edges(ne)
        ngll1 = Tdomain%sEdge(nne)%ngll
        orient_e = io_comm%orient_edges(ne)
        if(Tdomain%sEdge(nne)%solid)then   ! solid part
            call Comm_Edge_VectorProperty(ngll1,orient_e,                &
                io_comm%TakeForces(ngll:ngll+ngll1-3,0:2),  &
                Tdomain%sEdge(nne)%Forces(:,:))
            ngll = ngll+ngll1-2
            if(Tdomain%sEdge(nne)%PML)then
                call Comm_Edge_VectorProperty(ngll1,orient_e,                         &
                    io_comm%TakeForcesPML(ngllPML:ngllPML+ngll1-3,1,0:2),  &
                    Tdomain%sEdge(nne)%spml%Forces1(:,:))
                call Comm_Edge_VectorProperty(ngll1,orient_e,                         &
                    io_comm%TakeForcesPML(ngllPML:ngllPML+ngll1-3,2,0:2),  &
                    Tdomain%sEdge(nne)%spml%Forces2(:,:))
                call Comm_Edge_VectorProperty(ngll1,orient_e,                         &
                    io_comm%TakeForcesPML(ngllPML:ngllPML+ngll1-3,3,0:2),  &
                    Tdomain%sEdge(nne)%spml%Forces3(:,:))
                ngllPML = ngllPML+ngll1-2
            end if

        else   ! fluid part
            call Comm_Edge_ScalarProperty(ngll1,orient_e,                    &
                io_comm%TakeForcesFl(ngll_F:ngll_F+ngll1-3),  &
                Tdomain%sEdge(nne)%ForcesFl(:))
            ngll_F = ngll_F+ngll1-2
            if(Tdomain%sEdge(nne)%PML)then
                call Comm_Edge_ScalarProperty(ngll1,orient_e,                             &
                    io_comm%TakeForcesPMLFl(ngllPML_F:ngllPML_F+ngll1-3,1),  &
                    Tdomain%sEdge(nne)%spml%ForcesFl1(:))
                call Comm_Edge_ScalarProperty(ngll1,orient_e,                             &
                    io_comm%TakeForcesPMLFl(ngllPML_F:ngllPML_F+ngll1-3,2),  &
                    Tdomain%sEdge(nne)%spml%ForcesFl2(:))
                call Comm_Edge_ScalarProperty(ngll1,orient_e,                             &
                    io_comm%TakeForcesPMLFl(ngllPML_F:ngllPML_F+ngll1-3,3),  &
                    Tdomain%sEdge(nne)%spml%ForcesFl3(:))
                ngllPML_F = ngllPML_F+ngll1-2
            end if

        end if
    end do

    return
end subroutine Comm_Forces_Edge

subroutine Comm_Forces_Face(Tdomain,io_comm,ngll,ngll_F,ngllPML,ngllPML_F)
    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    type(comm), intent(inout) :: io_comm
    integer, intent(inout) :: ngll,ngllPML,ngll_F,ngllPML_F
    integer :: ngll1,ngll2,nf,nnf,orient_f


    do nf = 0,io_comm%nb_faces-1
        nnf = io_comm%faces(nf)
        ngll1 = Tdomain%sFace(nnf)%ngll1
        ngll2 = Tdomain%sFace(nnf)%ngll2
        orient_f = io_comm%orient_faces(nf)
        if(Tdomain%sFace(nnf)%solid)then   ! solid part
            call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                      &
                io_comm%TakeForces(ngll:ngll+(ngll1-2)*(ngll2-2)-1,:),  &
                Tdomain%sFace(nnf)%Forces(:,:,:))
            ngll = ngll+(ngll1-2)*(ngll2-2)
            if(Tdomain%sFace(nnf)%PML)then
                call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                               &
                    io_comm%TakeForcesPML(ngllPML:ngllPML+(ngll1-2)*(ngll2-2)-1,1,0:2),  &
                    Tdomain%sFace(nnf)%spml%Forces1(:,:,:))
                call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                               &
                    io_comm%TakeForcesPML(ngllPML:ngllPML+(ngll1-2)*(ngll2-2)-1,2,0:2),  &
                    Tdomain%sFace(nnf)%spml%Forces2(:,:,:))
                call Comm_Face_VectorProperty(ngll1,ngll2,orient_f,                               &
                    io_comm%TakeForcesPML(ngllPML:ngllPML+(ngll1-2)*(ngll2-2)-1,3,0:2),  &
                    Tdomain%sFace(nnf)%spml%Forces3(:,:,:))
                ngllPML = ngllPML+(ngll1-2)*(ngll2-2)
            end if

        else   ! fluid part
            call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                          &
                io_comm%TakeForcesFl(ngll_F:ngll_F+(ngll1-2)*(ngll2-2)-1),  &
                Tdomain%sFace(nnf)%ForcesFl(:,:))
            ngll_F = ngll_F+(ngll1-2)*(ngll2-2)
            if(Tdomain%sFace(nnf)%PML)then
                call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                    io_comm%TakeForcesPMLFl(ngllPML_F:ngllPML_F+(ngll1-2)*(ngll2-2)-1,1),  &
                    Tdomain%sFace(nnf)%spml%ForcesFl1(:,:))
                call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                    io_comm%TakeForcesPMLFl(ngllPML_F:ngllPML_F+(ngll1-2)*(ngll2-2)-1,2),  &
                    Tdomain%sFace(nnf)%spml%ForcesFl2(:,:))
                call Comm_Face_ScalarProperty(ngll1,ngll2,orient_f,                                   &
                    io_comm%TakeForcesPMLFl(ngllPML_F:ngllPML_F+(ngll1-2)*(ngll2-2)-1,3),  &
                    Tdomain%sFace(nnf)%spml%ForcesFl3(:,:))
                ngllPML_F = ngllPML_F+(ngll1-2)*(ngll2-2)
            end if

        end if
    end do

    return
end subroutine Comm_Forces_Face

subroutine Comm_Forces_Vertex (Tdomain,io_comm,ngll,ngll_F,ngllPML,ngllPML_F)
    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    type(comm), intent(inout) :: io_comm
    integer, intent(inout) :: ngll,ngllPML,ngll_F,ngllPML_F
    integer :: i,nv


    do i = 0,io_comm%nb_vertices-1
        nv =  io_comm%vertices(i)
        if(Tdomain%sVertex(nv)%solid)then   ! solid part
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + io_comm%TakeForces(ngll,0:2)
            ngll = ngll + 1
            if (Tdomain%sVertex(nv)%PML) then
                Tdomain%sVertex(nv)%spml%Forces1(0:2) = Tdomain%sVertex(nv)%spml%Forces1(0:2) + io_comm%TakeForcesPML(ngllPML,1,0:2)
                Tdomain%sVertex(nv)%spml%Forces2(0:2) = Tdomain%sVertex(nv)%spml%Forces2(0:2) + io_comm%TakeForcesPML(ngllPML,2,0:2)
                Tdomain%sVertex(nv)%spml%Forces3(0:2) = Tdomain%sVertex(nv)%spml%Forces3(0:2) + io_comm%TakeForcesPML(ngllPML,3,0:2)
                ngllPML = ngllPML + 1
            endif
        else   ! fluid part
            Tdomain%sVertex(nv)%ForcesFl = Tdomain%sVertex(nv)%ForcesFl + io_comm%TakeForcesFl(ngll_F)
            ngll_F = ngll_F + 1
            if (Tdomain%sVertex(nv)%PML) then
                Tdomain%sVertex(nv)%spml%ForcesFl1 = Tdomain%sVertex(nv)%spml%ForcesFl1 + io_comm%TakeForcesPMLFl(ngllPML_F,1)
                Tdomain%sVertex(nv)%spml%ForcesFl2 = Tdomain%sVertex(nv)%spml%ForcesFl2 + io_comm%TakeForcesPMLFl(ngllPML_F,2)
                Tdomain%sVertex(nv)%spml%ForcesFl3 = Tdomain%sVertex(nv)%spml%ForcesFl3 + io_comm%TakeForcesPMLFl(ngllPML_F,3)
                ngllPML_F = ngllPML_F + 1
            endif

        end if
    enddo


    return
end subroutine Comm_Forces_Vertex

subroutine Comm_Mass_Edge (Tdomain,io_comm,ngll,ngllPML)
    use sdomain
    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    type (comm), intent(inout) :: io_comm
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: ngll1,i,j,ne


    do i = 0,io_comm%nb_edges-1
        ne = io_comm%edges(i)
        ngll1 = Tdomain%sEdge(ne)%ngll

        if ( io_comm%orient_edges(i) == 0 ) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + io_comm%Take(ngll)
                ngll = ngll + 1
            enddo
            if (Tdomain%sEdge(ne)%PML) then
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    Tdomain%sEdge(ne)%spml%DumpMass(j,0:2) = Tdomain%sEdge(ne)%spml%DumpMass(j,0:2) + io_comm%TakePML(ngllPML,0:2)
                    if (Tdomain%any_FPML) then
                        Tdomain%sEdge(ne)%spml%Ivx(j) = Tdomain%sEdge(ne)%spml%Ivx(j) + io_comm%TakePML(ngllPML,3)
                        Tdomain%sEdge(ne)%spml%Ivy(j) = Tdomain%sEdge(ne)%spml%Ivy(j) + io_comm%TakePML(ngllPML,4)
                        Tdomain%sEdge(ne)%spml%Ivz(j) = Tdomain%sEdge(ne)%spml%Ivz(j) + io_comm%TakePML(ngllPML,5)
                    endif
                    ngllPML = ngllPML + 1
                enddo
            endif

        else if ( io_comm%orient_edges(i) == 1 ) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sEdge(ne)%MassMat(ngll1-1-j) = Tdomain%sEdge(ne)%MassMat(ngll1-1-j) + io_comm%Take(ngll)
                ngll = ngll + 1
            enddo
            if (Tdomain%sEdge(ne)%PML) then
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    Tdomain%sEdge(ne)%spml%DumpMass(ngll1-1-j,0:2) = Tdomain%sEdge(ne)%spml%DumpMass(ngll1-1-j,0:2) + io_comm%TakePML(ngllPML,0:2)
                    if (Tdomain%any_FPML) then
                        Tdomain%sEdge(ne)%spml%Ivx(ngll1-1-j) = Tdomain%sEdge(ne)%spml%Ivx(ngll1-1-j) + io_comm%TakePML(ngllPML,3)
                        Tdomain%sEdge(ne)%spml%Ivy(ngll1-1-j) = Tdomain%sEdge(ne)%spml%Ivy(ngll1-1-j) + io_comm%TakePML(ngllPML,4)
                        Tdomain%sEdge(ne)%spml%Ivz(ngll1-1-j) = Tdomain%sEdge(ne)%spml%Ivz(ngll1-1-j) + io_comm%TakePML(ngllPML,5)
                    endif
                    ngllPML = ngllPML + 1
                enddo
            endif

        else
            print*,'Pb with coherency number for edges in Comm_Mass_Edge', i, io_comm%orient_edges(i)
            STOP 1
        endif

    enddo

    return
end subroutine Comm_Mass_Edge

subroutine Comm_Mass_Face (Tdomain,io_comm,ngll,ngllPML)
    use sdomain
    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    type (comm), intent(inout) :: io_comm
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: ngll1,ngll2,i,j,k,nf


    do i = 0,io_comm%nb_faces-1
        nf = io_comm%faces(i)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2

        if ( io_comm%orient_faces(i) == 0 ) then

            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(k,j) = Tdomain%sFace(nf)%MassMat(k,j) + io_comm%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%spml%DumpMass(k,j,0:2) = Tdomain%sFace(nf)%spml%DumpMass(k,j,0:2) + io_comm%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(k,j) = Tdomain%sFace(nf)%spml%Ivx(k,j) + io_comm%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(k,j) = Tdomain%sFace(nf)%spml%Ivy(k,j) + io_comm%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(k,j) = Tdomain%sFace(nf)%spml%Ivz(k,j) + io_comm%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( io_comm%orient_faces(i) == 1 ) then
            nf = io_comm%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-k,j) = Tdomain%sFace(nf)%MassMat(ngll1-1-k,j) + io_comm%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-k,j,0:2) = Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-k,j,0:2) + io_comm%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(ngll1-1-k,j) = Tdomain%sFace(nf)%spml%Ivx(ngll1-1-k,j) + io_comm%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(ngll1-1-k,j) = Tdomain%sFace(nf)%spml%Ivy(ngll1-1-k,j) + io_comm%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(ngll1-1-k,j) = Tdomain%sFace(nf)%spml%Ivz(ngll1-1-k,j) + io_comm%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( io_comm%orient_faces(i) == 2 ) then
            nf = io_comm%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(k,ngll2-1-j) = Tdomain%sFace(nf)%MassMat(k,ngll2-1-j) + io_comm%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%spml%DumpMass(k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%spml%DumpMass(k,ngll2-1-j,0:2) + io_comm%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivx(k,ngll2-1-j) + io_comm%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivy(k,ngll2-1-j) + io_comm%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivz(k,ngll2-1-j) + io_comm%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( io_comm%orient_faces(i) == 3 ) then
            nf = io_comm%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%MassMat(ngll1-1-k,ngll2-1-j) + io_comm%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-k,ngll2-1-j,0:2) + io_comm%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivx(ngll1-1-k,ngll2-1-j) + io_comm%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivy(ngll1-1-k,ngll2-1-j) + io_comm%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%spml%Ivz(ngll1-1-k,ngll2-1-j) + io_comm%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( io_comm%orient_faces(i) == 4 ) then
            nf = io_comm%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + io_comm%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%spml%DumpMass(j,k,0:2) = Tdomain%sFace(nf)%spml%DumpMass(j,k,0:2) + io_comm%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(j,k) = Tdomain%sFace(nf)%spml%Ivx(j,k) + io_comm%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(j,k) = Tdomain%sFace(nf)%spml%Ivy(j,k) + io_comm%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(j,k) = Tdomain%sFace(nf)%spml%Ivz(j,k) + io_comm%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( io_comm%orient_faces(i) == 5 ) then
            nf = io_comm%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-j,k) = Tdomain%sFace(nf)%MassMat(ngll1-1-j,k) + io_comm%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-j,k,0:2) = Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-j,k,0:2) + io_comm%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(ngll1-1-j,k) = Tdomain%sFace(nf)%spml%Ivx(ngll1-1-j,k) + io_comm%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(ngll1-1-j,k) = Tdomain%sFace(nf)%spml%Ivy(ngll1-1-j,k) + io_comm%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(ngll1-1-j,k) = Tdomain%sFace(nf)%spml%Ivz(ngll1-1-j,k) + io_comm%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( io_comm%orient_faces(i) == 6 ) then
            nf = io_comm%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(j,ngll2-1-k) = Tdomain%sFace(nf)%MassMat(j,ngll2-1-k) + io_comm%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%spml%DumpMass(j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%spml%DumpMass(j,ngll2-1-k,0:2) + io_comm%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivx(j,ngll2-1-k) + io_comm%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivy(j,ngll2-1-k) + io_comm%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivz(j,ngll2-1-k) + io_comm%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( io_comm%orient_faces(i) == 7 ) then
            nf = io_comm%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%MassMat(ngll1-1-j,ngll2-1-k) + io_comm%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%spml%DumpMass(ngll1-1-j,ngll2-1-k,0:2) + io_comm%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%spml%Ivx(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivx(ngll1-1-j,ngll2-1-k) + io_comm%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%spml%Ivy(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivy(ngll1-1-j,ngll2-1-k) + io_comm%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%spml%Ivz(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%spml%Ivz(ngll1-1-j,ngll2-1-k) + io_comm%TakePML(ngllPML,5)
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

subroutine Comm_Mass_Vertex (Tdomain,io_comm,ngll,ngllPML)
    use sdomain
    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    type (comm), intent(inout) :: io_comm
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: i,nv


    do i = 0,io_comm%nb_vertices-1
        nv =  io_comm%vertices(i)
        Tdomain%svertex(nv)%MassMat = Tdomain%svertex(nv)%MassMat + io_comm%Take(ngll)
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%PML) then
            Tdomain%svertex(nv)%spml%DumpMass(0:2) = Tdomain%svertex(nv)%spml%DumpMass(0:2) + io_comm%TakePML(ngllPML,0:2)
            if (Tdomain%any_FPML) then
                Tdomain%sVertex(nv)%spml%Ivx(0) = Tdomain%sVertex(nv)%spml%Ivx(0) + io_comm%TakePML(ngllPML,3)
                Tdomain%sVertex(nv)%spml%Ivy(0) = Tdomain%sVertex(nv)%spml%Ivy(0) + io_comm%TakePML(ngllPML,4)
                Tdomain%sVertex(nv)%spml%Ivz(0) = Tdomain%sVertex(nv)%spml%Ivz(0) + io_comm%TakePML(ngllPML,5)
            endif
            ngllPML = ngllPML + 1
        endif
    enddo

    return
end subroutine Comm_Mass_Vertex

end module scommutils

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
