subroutine StoF_coupling(Tdomain)
    ! from solid to fluid: velocity (dot) normal
    use sdomain
    use mpi
    use scomm
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer  :: n,nf,ne,nv,nfs,ngll1,ngll2,ngll,i,j,nnf,nne,nnv,orient_e,    &
        code,ngllSF,ngllSF_PML
    real, dimension(:,:,:), allocatable  :: VelocFace

    ! init.
    do nf = 0,Tdomain%SF%SF_n_faces-1
        Tdomain%SF%SF_face(nf)%vn = 0d0
        if(Tdomain%SF%SF_Face(nf)%PML)then
            Tdomain%SF%SF_Face(nf)%vn1 = 0d0
            Tdomain%SF%SF_Face(nf)%vn2 = 0d0
            Tdomain%SF%SF_Face(nf)%vn3 = 0d0
        end if
    end do
    do ne = 0,Tdomain%SF%SF_n_edges-1
        Tdomain%SF%SF_edge(ne)%vn = 0d0
        if(Tdomain%SF%SF_Edge(ne)%PML)then
            Tdomain%SF%SF_Edge(ne)%vn1 = 0d0
            Tdomain%SF%SF_Edge(ne)%vn2 = 0d0
            Tdomain%SF%SF_Edge(ne)%vn3 = 0d0
        end if
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        Tdomain%SF%SF_vertex(nv)%vn = 0d0
        if(Tdomain%SF%SF_Vertex(nv)%PML)then
            Tdomain%SF%SF_Vertex(nv)%vn1 = 0d0
            Tdomain%SF%SF_Vertex(nv)%vn2 = 0d0
            Tdomain%SF%SF_Vertex(nv)%vn3 = 0d0
        end if
    end do

    do nf = 0,Tdomain%SF%SF_n_faces-1
        nfs = Tdomain%SF%SF_Face(nf)%Face(1)
        if(nfs < 0) cycle   ! solid face not on this proc
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        allocate(VelocFace(0:ngll1-1,0:ngll2-1,0:2))
        VelocFace = 0d0
        VelocFace(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nfs)%V0(:,:,:)
        ! deassemblage: particle velocity from edges and vertices to faces
        do i = 0,3
            ne = Tdomain%SF%SF_Face(nf)%Near_Edges(i)
            nne = Tdomain%SF%SF_Edge(ne)%Edge(1)
            ngll = Tdomain%SF%SF_Edge(ne)%ngll
            orient_e = Tdomain%SF%SF_face(nf)%Near_Edges_Orient(i)
            call Vector_Edge2Face(i,ngll1,ngll2,ngll,orient_e,         &
                VelocFace,Tdomain%sEdge(nne)%V0)
            nv = Tdomain%SF%SF_Face(nf)%Near_Vertices(i)
            nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(1)
            call Vector_Vertex2Face(i,ngll1,ngll2,    &
                VelocFace,Tdomain%sVertex(nnv)%V0)
        end do
        ! dot product on faces: normal.velocity
        do i = 0,2
            Tdomain%SF%SF_Face(nf)%vn(0:ngll1-1,0:ngll2-1) = Tdomain%SF%SF_Face(nf)%vn(0:ngll1-1,0:ngll2-1)+   &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,i)*VelocFace(0:ngll1-1,0:ngll2-1,i)
        end do
        if(Tdomain%SF%SF_Face(nf)%PML)then
            Tdomain%SF%SF_Face(nf)%vn1(0:ngll1-1,0:ngll2-1) =  &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,0)*VelocFace(0:ngll1-1,0:ngll2-1,0)
            Tdomain%SF%SF_Face(nf)%vn2(0:ngll1-1,0:ngll2-1) =  &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,1)*VelocFace(0:ngll1-1,0:ngll2-1,1)
            Tdomain%SF%SF_Face(nf)%vn3(0:ngll1-1,0:ngll2-1) =  &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,2)*VelocFace(0:ngll1-1,0:ngll2-1,2)
        end if

        deallocate(VelocFace)
        ! assemblage: (normal.velocity) from face to edges and vertices
        do i = 0,3
            ne = Tdomain%SF%SF_Face(nf)%Near_Edges(i)
            ngll = Tdomain%SF%SF_Edge(ne)%ngll
            orient_e = Tdomain%SF%SF_face(nf)%Near_Edges_Orient(i)
            nv = Tdomain%SF%SF_Face(nf)%Near_Vertices(i)
            if(Tdomain%SF%SF_Edge(ne)%PML)then
                call Scalar_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_Face(nf)%Vn1,Tdomain%SF%SF_Edge(ne)%Vn1)
                call Scalar_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_Face(nf)%Vn2,Tdomain%SF%SF_Edge(ne)%Vn2)
                call Scalar_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_Face(nf)%Vn3,Tdomain%SF%SF_Edge(ne)%Vn3)
            else
                call Scalar_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_Face(nf)%Vn,Tdomain%SF%SF_Edge(ne)%Vn)
            end if
            if(Tdomain%SF%SF_Vertex(nv)%PML)then
                call Scalar_Face2Vertex(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_Face(nf)%Vn1,Tdomain%SF%SF_Vertex(nv)%Vn1)
                call Scalar_Face2Vertex(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_Face(nf)%Vn2,Tdomain%SF%SF_Vertex(nv)%Vn2)
                call Scalar_Face2Vertex(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_Face(nf)%Vn3,Tdomain%SF%SF_Vertex(nv)%Vn3)
            else
                call Scalar_Face2Vertex(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_Face(nf)%Vn,Tdomain%SF%SF_Vertex(nv)%Vn)
            end if
        end do
    end do


    ! now we can exchange values for SF sides on different procs
    if(Tdomain%nb_procs > 1)then
        do n = 0,Tdomain%nb_procs-1
            call Comm_Forces_Complete_StoF(n,Tdomain)
            call Comm_Forces_Complete_StoF_PML(n,Tdomain)
        end do
        ! now we can exchange force values with other procs
        call exchange_sem_forces_StoF(Tdomain)

        ! assemblage on external GLLs
        do n = 0,Tdomain%nb_procs-1
            ngllSF = 0 ; ngllSF_PML = 0
            call Comm_Forces_FaceSF_StoF(Tdomain,n,ngllSF,ngllSF_PML)
            call Comm_Forces_EdgeSF_StoF(Tdomain,n,ngllSF,ngllSF_PML)
            call Comm_Forces_VertexSF_StoF(Tdomain,n,ngllSF,ngllSF_PML)
            if(ngllSF /= Tdomain%sComm(n)%ngllSF) stop "PB in counting SF GLL nodes."
        enddo
    end if


    ! now the coupling Solid -> Fluid can be done.
    do nf = 0,Tdomain%SF%SF_n_faces-1
        nnf = Tdomain%SF%SF_Face(nf)%Face(0)
        if(nnf < 0) cycle  ! fluid face not on this proc
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        if(Tdomain%SF%SF_Face(nf)%PML)then
            Tdomain%sFace(nnf)%spml%ForcesFl1(1:ngll1-2,1:ngll2-2) = Tdomain%sFace(nnf)%spml%ForcesFl1(1:ngll1-2,1:ngll2-2) +  &
                Tdomain%SF%SF_Face(nf)%Vn1(1:ngll1-2,1:ngll2-2)
            Tdomain%sFace(nnf)%spml%ForcesFl2(1:ngll1-2,1:ngll2-2) = Tdomain%sFace(nnf)%spml%ForcesFl2(1:ngll1-2,1:ngll2-2) +  &
                Tdomain%SF%SF_Face(nf)%Vn2(1:ngll1-2,1:ngll2-2)
            Tdomain%sFace(nnf)%spml%ForcesFl3(1:ngll1-2,1:ngll2-2) = Tdomain%sFace(nnf)%spml%ForcesFl3(1:ngll1-2,1:ngll2-2) +  &
                Tdomain%SF%SF_Face(nf)%Vn3(1:ngll1-2,1:ngll2-2)
        else
            Tdomain%sFace(nnf)%ForcesFl(1:ngll1-2,1:ngll2-2) = Tdomain%sFace(nnf)%ForcesFl(1:ngll1-2,1:ngll2-2) +  &
                Tdomain%SF%SF_Face(nf)%Vn(1:ngll1-2,1:ngll2-2)
        end if
    end do
    do ne = 0,Tdomain%SF%SF_n_edges-1
        nne = Tdomain%SF%SF_Edge(ne)%Edge(0)
        if(nne < 0) cycle  ! fluid edge not on this proc
        ngll = Tdomain%SF%SF_Edge(ne)%ngll
        if(Tdomain%SF%SF_Edge(ne)%PML)then
            Tdomain%sEdge(nne)%spml%ForcesFl1(1:ngll-2) = Tdomain%sEdge(nne)%spml%ForcesFl1(1:ngll-2) +  &
                Tdomain%SF%SF_Edge(ne)%Vn1(1:ngll-2)
            Tdomain%sEdge(nne)%spml%ForcesFl2(1:ngll-2) = Tdomain%sEdge(nne)%spml%ForcesFl2(1:ngll-2) +  &
                Tdomain%SF%SF_Edge(ne)%Vn2(1:ngll-2)
            Tdomain%sEdge(nne)%spml%ForcesFl3(1:ngll-2) = Tdomain%sEdge(nne)%spml%ForcesFl3(1:ngll-2) +  &
                Tdomain%SF%SF_Edge(ne)%Vn3(1:ngll-2)
        else
            Tdomain%sEdge(nne)%ForcesFl(1:ngll-2) = Tdomain%sEdge(nne)%ForcesFl(1:ngll-2) +  &
                Tdomain%SF%SF_Edge(ne)%Vn(1:ngll-2)
        end if
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(0)
        if(nnv < 0) cycle  ! fluid vertex not on this proc
        if(Tdomain%SF%SF_Vertex(nv)%PML)then
            Tdomain%sVertex(nnv)%spml%ForcesFl1 = Tdomain%sVertex(nnv)%spml%ForcesFl1 + Tdomain%SF%SF_Vertex(nv)%Vn1
            Tdomain%sVertex(nnv)%spml%ForcesFl2 = Tdomain%sVertex(nnv)%spml%ForcesFl2 + Tdomain%SF%SF_Vertex(nv)%Vn2
            Tdomain%sVertex(nnv)%spml%ForcesFl3 = Tdomain%sVertex(nnv)%spml%ForcesFl3 + Tdomain%SF%SF_Vertex(nv)%Vn3
        else
            Tdomain%sVertex(nnv)%ForcesFl = Tdomain%sVertex(nnv)%ForcesFl + Tdomain%SF%SF_Vertex(nv)%Vn
        end if
    end do


end subroutine StoF_coupling
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine FtoS_coupling(Tdomain)
    ! from fluid to solid: normal times pressure (= -rho . VelPhi)
    use sdomain
    use mpi
    use scomm
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer  :: n,nf,ne,nv,nff,ngll1,ngll2,ngll,i,j,k,nnf,nne,nnv,orient_e,    &
        n_rings,shift,I_give_to,I_take_from,code,which_elem,ngllSF,ngllSF_PML
    integer, dimension(MPI_STATUS_SIZE)   :: statut
    integer, parameter   :: etiquette = 100
    real, dimension(:,:), allocatable  :: PressFace

    ! init.
    do nf = 0,Tdomain%SF%SF_n_faces-1
        Tdomain%SF%SF_face(nf)%pn = 0d0
        if(Tdomain%SF%SF_Face(nf)%PML)then
            Tdomain%SF%SF_face(nf)%pn1 = 0d0
            Tdomain%SF%SF_face(nf)%pn2 = 0d0
            Tdomain%SF%SF_face(nf)%pn3 = 0d0
        end if
    end do
    do ne = 0,Tdomain%SF%SF_n_edges-1
        Tdomain%SF%SF_edge(ne)%pn = 0d0
        if(Tdomain%SF%SF_Edge(ne)%PML)then
            Tdomain%SF%SF_edge(ne)%pn1 = 0d0
            Tdomain%SF%SF_edge(ne)%pn2 = 0d0
            Tdomain%SF%SF_edge(ne)%pn3 = 0d0
        end if
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        Tdomain%SF%SF_vertex(nv)%pn = 0d0
        if(Tdomain%SF%SF_vertex(nv)%PML)then
            Tdomain%SF%SF_vertex(nv)%pn1 = 0d0
            Tdomain%SF%SF_vertex(nv)%pn2 = 0d0
            Tdomain%SF%SF_vertex(nv)%pn3 = 0d0
        end if
    end do


    do nf = 0,Tdomain%SF%SF_n_faces-1
        nff = Tdomain%SF%SF_Face(nf)%Face(0)
        if(nff < 0) cycle   ! fluid face not on this proc
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        allocate(PressFace(0:ngll1-1,0:ngll2-1))
        PressFace = 0d0
        PressFace(1:ngll1-2,1:ngll2-2) = Tdomain%sFace(nff)%VelPhi(:,:)
        ! deassemblage: potential velocity from edges and vertices to faces
        do i = 0,3
            ne = Tdomain%SF%SF_Face(nf)%Near_Edges(i)
            nne = Tdomain%SF%SF_Edge(ne)%Edge(0)
            ngll = Tdomain%SF%SF_Edge(ne)%ngll
            orient_e = Tdomain%SF%SF_face(nf)%Near_Edges_Orient(i)
            call Scalar_Edge2Face(i,ngll1,ngll2,ngll,orient_e,         &
                PressFace,Tdomain%sEdge(nne)%VelPhi)
            nv = Tdomain%SF%SF_Face(nf)%Near_Vertices(i)
            nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(0)
            call Scalar_Vertex2Face(i,ngll1,ngll2,    &
                PressFace,Tdomain%sVertex(nnv)%VelPhi)
        end do
        ! on faces: normal times pressure = -(rho * VelPhi) * normal
        do i = 0,2
            Tdomain%SF%SF_Face(nf)%pn(0:ngll1-1,0:ngll2-1,i) = -1d0*PressFace(0:ngll1-1,0:ngll2-1)*  &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,i)
        end do
        if(Tdomain%SF%SF_Face(nf)%PML)then
            Tdomain%SF%SF_Face(nf)%pn1(0:ngll1-1,0:ngll2-1,0) = -1d0*PressFace(0:ngll1-1,0:ngll2-1)*  &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,0)
            Tdomain%SF%SF_Face(nf)%pn2(0:ngll1-1,0:ngll2-1,1) = -1d0*PressFace(0:ngll1-1,0:ngll2-1)*  &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,1)
            Tdomain%SF%SF_Face(nf)%pn3(0:ngll1-1,0:ngll2-1,2) = -1d0*PressFace(0:ngll1-1,0:ngll2-1)*  &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,2)
        end if
        deallocate(PressFace)
        ! assemblage: (normal times pressure) from face to edges and vertices
        do i = 0,3
            ne = Tdomain%SF%SF_Face(nf)%Near_Edges(i)
            ngll = Tdomain%SF%SF_Edge(ne)%ngll
            orient_e = Tdomain%SF%SF_face(nf)%Near_Edges_Orient(i)
            nv = Tdomain%SF%SF_Face(nf)%Near_Vertices(i)
            if(Tdomain%SF%SF_Edge(ne)%PML)then
                call Vector_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_Face(nf)%pn1,Tdomain%SF%SF_Edge(ne)%pn1)
                call Vector_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_Face(nf)%pn2,Tdomain%SF%SF_Edge(ne)%pn2)
                call Vector_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_Face(nf)%pn3,Tdomain%SF%SF_Edge(ne)%pn3)
            else
                call Vector_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_Face(nf)%pn,Tdomain%SF%SF_Edge(ne)%pn)
            end if
            if(Tdomain%SF%SF_Vertex(nv)%PML)then
                call Vector_Face2Vertex(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_Face(nf)%pn1,Tdomain%SF%SF_Vertex(nv)%pn1)
                call Vector_Face2Vertex(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_Face(nf)%pn2,Tdomain%SF%SF_Vertex(nv)%pn2)
                call Vector_Face2Vertex(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_Face(nf)%pn3,Tdomain%SF%SF_Vertex(nv)%pn3)
            else
                call Vector_Face2Vertex(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_Face(nf)%pn,Tdomain%SF%SF_Vertex(nv)%pn)
            end if
        end do
    end do


    ! now we can exchange values for SF sides on different procs
    if(Tdomain%nb_procs > 1)then
        do n = 0,Tdomain%nb_procs-1
            call Comm_Forces_Complete_FtoS(n,Tdomain)
            call Comm_Forces_Complete_FtoS_PML(n,Tdomain)
        end do
        ! now we can exchange force values with proc n
        call exchange_sem_forces_FtoS(Tdomain)
        ! assemblage on external GLLs
        do n = 0,Tdomain%nb_procs-1
            ngllSF = 0 ; ngllSF_PML = 0
            call Comm_Forces_FaceSF_FtoS(Tdomain,n,ngllSF,ngllSF_PML)
            call Comm_Forces_EdgeSF_FtoS(Tdomain,n,ngllSF,ngllSF_PML)
            call Comm_Forces_VertexSF_FtoS(Tdomain,n,ngllSF,ngllSF_PML)
            if(ngllSF /= Tdomain%sComm(n)%ngllSF) stop "PB in counting SF GLL nodes."
        enddo

    end if

    ! now the coupling Fluid -> Solid can be done.
    do nf = 0,Tdomain%SF%SF_n_faces-1
        nnf = Tdomain%SF%SF_Face(nf)%Face(1)
        if(nnf < 0) cycle  ! solid face not on this proc
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        if(.not. Tdomain%sFace(nnf)%PML)then
            Tdomain%sFace(nnf)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%SF%SF_Face(nf)%save_forces(1:ngll1-2,1:ngll2-2,0:2) +  &
                Tdomain%SF%SF_Face(nf)%pn(1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nnf)%Displ(:,:,:) = Tdomain%SF%SF_Face(nf)%save_displ(:,:,:)
        else
            Tdomain%sFace(nnf)%spml%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nnf)%spml%Forces1(1:ngll1-2,1:ngll2-2,0:2) +  &
                Tdomain%SF%SF_Face(nf)%pn1(1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nnf)%spml%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nnf)%spml%Forces2(1:ngll1-2,1:ngll2-2,0:2) +  &
                Tdomain%SF%SF_Face(nf)%pn2(1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nnf)%spml%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nnf)%spml%Forces3(1:ngll1-2,1:ngll2-2,0:2) +  &
                Tdomain%SF%SF_Face(nf)%pn3(1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nnf)%spml%Veloc1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%SF%SF_Face(nf)%save_Veloc1(1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nnf)%spml%Veloc2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%SF%SF_Face(nf)%save_Veloc2(1:ngll1-2,1:ngll2-2,0:2)
            Tdomain%sFace(nnf)%spml%Veloc3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%SF%SF_Face(nf)%save_Veloc3(1:ngll1-2,1:ngll2-2,0:2)
        end if
    end do
    do ne = 0,Tdomain%SF%SF_n_edges-1
        nne = Tdomain%SF%SF_Edge(ne)%Edge(1)
        if(nne < 0) cycle  ! solid edge not on this proc
        ngll = Tdomain%SF%SF_Edge(ne)%ngll
        if(.not. Tdomain%sEdge(nne)%PML)then
            Tdomain%sEdge(nne)%Forces(1:ngll-2,0:2) = Tdomain%SF%SF_Edge(ne)%save_forces(1:ngll-2,0:2) +  &
                Tdomain%SF%SF_Edge(ne)%pn(1:ngll-2,0:2)
            Tdomain%sEdge(nne)%Displ(:,:) = Tdomain%SF%SF_Edge(ne)%save_displ(:,:)
        else
            Tdomain%sEdge(nne)%spml%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(nne)%spml%Forces1(1:ngll-2,0:2) +  &
                Tdomain%SF%SF_Edge(ne)%pn1(1:ngll-2,0:2)
            Tdomain%sEdge(nne)%spml%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(nne)%spml%Forces2(1:ngll-2,0:2)+  &
                Tdomain%SF%SF_Edge(ne)%pn2(1:ngll-2,0:2)
            Tdomain%sEdge(nne)%spml%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(nne)%spml%Forces3(1:ngll-2,0:2)+  &
                Tdomain%SF%SF_Edge(ne)%pn3(1:ngll-2,0:2)
            Tdomain%sEdge(nne)%spml%Veloc1(1:ngll-2,0:2) = Tdomain%SF%SF_Edge(ne)%save_Veloc1(1:ngll-2,0:2)
            Tdomain%sEdge(nne)%spml%Veloc2(1:ngll-2,0:2) = Tdomain%SF%SF_Edge(ne)%save_Veloc2(1:ngll-2,0:2)
            Tdomain%sEdge(nne)%spml%Veloc3(1:ngll-2,0:2) = Tdomain%SF%SF_Edge(ne)%save_Veloc3(1:ngll-2,0:2)
        end if
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(1)
        if(nnv < 0) cycle  ! solid vertex not on this proc
        if(.not. Tdomain%sVertex(nnv)%PML)then
            Tdomain%sVertex(nnv)%Forces(0:2) = Tdomain%SF%SF_Vertex(nv)%save_forces(0:2) + Tdomain%SF%SF_Vertex(nv)%pn(0:2)
            Tdomain%sVertex(nnv)%Displ(:) = Tdomain%SF%SF_Vertex(nv)%save_displ(:)
        else
            Tdomain%sVertex(nnv)%spml%Forces1(0:2) = Tdomain%sVertex(nnv)%spml%Forces1(0:2) + Tdomain%SF%SF_Vertex(nv)%pn1(0:2)
            Tdomain%sVertex(nnv)%spml%Forces2(0:2) = Tdomain%sVertex(nnv)%spml%Forces2(0:2) + Tdomain%SF%SF_Vertex(nv)%pn2(0:2)
            Tdomain%sVertex(nnv)%spml%Forces3(0:2) = Tdomain%sVertex(nnv)%spml%Forces3(0:2) + Tdomain%SF%SF_Vertex(nv)%pn3(0:2)
            Tdomain%sVertex(nnv)%spml%Veloc1(0:2) = Tdomain%SF%SF_Vertex(nv)%save_Veloc1(0:2)
            Tdomain%sVertex(nnv)%spml%Veloc2(0:2) = Tdomain%SF%SF_Vertex(nv)%save_Veloc2(0:2)
            Tdomain%sVertex(nnv)%spml%Veloc3(0:2) = Tdomain%SF%SF_Vertex(nv)%save_Veloc3(0:2)
        end if
    end do


end subroutine FtoS_coupling
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
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine SF_solid_values_saving(Tdomain)
    ! saves values of fields Forces, Displ on the solid side of a SF object -
    !    for they are flushed in the 1st correction phase
    use sdomain
    implicit none

    type(domain), intent(inout)   ::  Tdomain
    integer  :: nf,ne,nv,nnf,nne,nnv

    do nf = 0,Tdomain%SF%SF_n_faces-1
        nnf = Tdomain%SF%SF_Face(nf)%Face(1)
        if(nnf < 0) cycle  ! solid side not on this proc
        if(.not.Tdomain%sFace(nnf)%PML)then
            Tdomain%SF%SF_Face(nf)%save_forces(:,:,:) = Tdomain%sFace(nnf)%Forces(:,:,:)
            Tdomain%SF%SF_Face(nf)%save_displ(:,:,:) = Tdomain%sFace(nnf)%Displ(:,:,:)
        else
            Tdomain%SF%SF_Face(nf)%save_Veloc1(:,:,:) = Tdomain%sFace(nnf)%spml%Veloc1(:,:,:)
            Tdomain%SF%SF_Face(nf)%save_Veloc2(:,:,:) = Tdomain%sFace(nnf)%spml%Veloc2(:,:,:)
            Tdomain%SF%SF_Face(nf)%save_Veloc3(:,:,:) = Tdomain%sFace(nnf)%spml%Veloc3(:,:,:)
        end if
    end do
    do ne = 0,Tdomain%SF%SF_n_edges-1
        nne = Tdomain%SF%SF_Edge(ne)%Edge(1)
        if(nne < 0) cycle  ! solid side not on this proc
        if(.not.Tdomain%sEdge(nne)%PML)then
            Tdomain%SF%SF_Edge(ne)%save_forces(:,:) = Tdomain%sEdge(nne)%Forces(:,:)
            Tdomain%SF%SF_Edge(ne)%save_displ(:,:) = Tdomain%sEdge(nne)%Displ(:,:)
        else
            Tdomain%SF%SF_Edge(ne)%save_Veloc1(:,:) = Tdomain%sEdge(nne)%spml%Veloc1(:,:)
            Tdomain%SF%SF_Edge(ne)%save_Veloc2(:,:) = Tdomain%sEdge(nne)%spml%Veloc2(:,:)
            Tdomain%SF%SF_Edge(ne)%save_Veloc3(:,:) = Tdomain%sEdge(nne)%spml%Veloc3(:,:)
        end if
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(1)
        if(nnv < 0) cycle  ! solid side not on this proc
        if(.not.Tdomain%sVertex(nnv)%PML)then
            Tdomain%SF%SF_Vertex(nv)%save_forces(:) = Tdomain%sVertex(nnv)%Forces(:)
            Tdomain%SF%SF_Vertex(nv)%save_displ(:) = Tdomain%sVertex(nnv)%Displ(:)
        else
            Tdomain%SF%SF_Vertex(nv)%save_Veloc1(:) = Tdomain%sVertex(nnv)%spml%Veloc1(:)
            Tdomain%SF%SF_Vertex(nv)%save_Veloc2(:) = Tdomain%sVertex(nnv)%spml%Veloc2(:)
            Tdomain%SF%SF_Vertex(nv)%save_Veloc3(:) = Tdomain%sVertex(nnv)%spml%Veloc3(:)
        end if
    end do

    return
end subroutine SF_solid_values_saving
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine Newmark_recorrect_solid(Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer  :: nf,ne,nv,nnf,nne,nnv,mat
    real  :: bega,gam1,dt

    bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
    gam1 = 1. / Tdomain%TimeD%gamma

    do nf = 0,Tdomain%SF%SF_n_faces-1
        nnf = Tdomain%SF%SF_Face(nf)%Face(1)
        if(nnf >= 0)then
            mat = Tdomain%sFace(nnf)%mat_index
            dt = Tdomain%sSubdomain(mat)%dt
            if(.not.Tdomain%sface(nnf)%PML)then
                !   call Correction_Face_Veloc(Tdomain%sface(nnf),bega,gam1,dt)
                call Correction_Face_Veloc(Tdomain%sface(nnf),dt)
            else
                if(Tdomain%sface(nnf)%FPML)then
                    call Correction_Face_FPML_Veloc(Tdomain%sface(nnf),dt,Tdomain%sSubdomain(mat)%freq)
                else
                    call Correction_Face_PML_Veloc(Tdomain%sface(nnf),dt)
                endif
            endif
        end if
    end do

    do ne = 0,Tdomain%SF%SF_n_edges-1
        nne = Tdomain%SF%SF_Edge(ne)%Edge(1)
        if(nne >= 0)then
            mat = Tdomain%sEdge(nne)%mat_index
            dt = Tdomain%sSubdomain(mat)%dt
            if(.not.Tdomain%sEdge(nne)%PML)then
                !   call Correction_Edge_Veloc(Tdomain%sEdge(nne),bega,gam1,dt)
                call Correction_Edge_Veloc(Tdomain%sEdge(nne),dt)
            else
                if(Tdomain%sEdge(nne)%FPML)then
                    call Correction_Edge_FPML_Veloc(Tdomain%sEdge(nne),dt,Tdomain%sSubdomain(mat)%freq)
                else
                    call Correction_Edge_PML_Veloc(Tdomain%sEdge(nne),dt)
                endif
            endif
        end if
    end do

    do nv = 0,Tdomain%SF%SF_n_vertices-1
        nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(1)
        if(nnv >= 0)then
            mat = Tdomain%sVertex(nnv)%mat_index
            dt = Tdomain%sSubdomain(mat)%dt
            if(.not.Tdomain%sVertex(nnv)%PML)then
                !     call Correction_Vertex_Veloc(Tdomain%sVertex(nnv),bega,gam1,dt)
                call Correction_Vertex_Veloc(Tdomain%sVertex(nnv),dt)
            else
                if(Tdomain%sVertex(nnv)%FPML)then
                    call Correction_Vertex_FPML_Veloc(Tdomain%sVertex(nnv),dt,Tdomain%sSubdomain(mat)%freq)
                else
                    call Correction_Vertex_PML_Veloc(Tdomain%sVertex(nnv),dt)
                endif
            endif
        end if
    end do

end subroutine Newmark_recorrect_solid
