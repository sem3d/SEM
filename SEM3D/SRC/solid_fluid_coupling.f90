subroutine StoF_coupling(Tdomain,rg)
    ! from solid to fluid: velocity (dot) normal
    use sdomain
    use mpi
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer, intent(in)   :: rg
    integer  :: n,nf,ne,nv,nfs,ngll1,ngll2,ngll,i,j,k,nnf,nne,nnv,orient_e,    &
        n_rings,shift,I_give_to,I_take_from,code,ngllSF
    integer, dimension(MPI_STATUS_SIZE)   :: statut
    integer, parameter   :: etiquette = 100
    real, dimension(:,:,:), allocatable  :: VelocFace

    ! init.
    do ne = 0,Tdomain%SF%SF_n_edges-1
        Tdomain%SF%SF_edge(ne)%vn = 0d0
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        Tdomain%SF%SF_vertex(nv)%vn = 0d0
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
        Tdomain%SF%SF_Face(nf)%vn(0:ngll1-1,0:ngll2-1) = 0d0
        do i = 0,2
            Tdomain%SF%SF_Face(nf)%vn(0:ngll1-1,0:ngll2-1) = Tdomain%SF%SF_Face(nf)%vn(0:ngll1-1,0:ngll2-1)+   &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,i)*VelocFace(0:ngll1-1,0:ngll2-1,i)
        end do
        deallocate(VelocFace)
        ! assemblage: (normal.velocity) from face to edges and vertices
        do i = 0,3
            ne = Tdomain%SF%SF_Face(nf)%Near_Edges(i)
            ngll = Tdomain%SF%SF_Edge(ne)%ngll
            orient_e = Tdomain%SF%SF_face(nf)%Near_Edges_Orient(i)
            call Scalar_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                Tdomain%SF%SF_Face(nf)%Vn,Tdomain%SF%SF_Edge(ne)%Vn)
            nv = Tdomain%SF%SF_Face(nf)%Near_Vertices(i)
            call Scalar_Face2Vertex(i,ngll1,ngll2,    &
                Tdomain%SF%SF_Face(nf)%Vn,Tdomain%SF%SF_Vertex(nv)%Vn)
        end do
    end do


    ! now we can exchange values for SF sides on different procs
    if(Tdomain%n_proc > 1)then
        do n = 0,Tdomain%n_proc-1
            call Comm_Forces_Complete_StoF(n,Tdomain)
        end do
        ! now we can exchange force values with proc n
        n = Tdomain%n_proc
        do shift = 1,n-1
            if(rg == 1) print*,"SHIFTV",shift
            I_give_to = rg + shift
            if(I_give_to > n-1) I_give_to = I_give_to-n
            if(rg == 1) print*,"GIVEV",I_give_to
        end do


        do shift = 1,n-1
            if(rg == 1) print*,"SHIFT",shift
            I_give_to = rg + shift
            if (I_give_to > n-1)   I_give_to = I_give_to - n
            if(rg == 1) print*,"GIVE",I_give_to
            I_take_from = rg - shift
            if (I_take_from < 0)   I_take_from = I_take_from + n
            if(rg == 1) print*,"TAKE",I_take_from
            if (mod(n,shift)==0 .and. shift/=1) then
                n_rings = shift
            else if (mod(n,n-shift)==0 .and. shift/=n-1) then
                n_rings = n-shift
            else if (mod(n,2)==0 .and. mod(shift,2)==0) then
                n_rings = 2
            else
                n_rings = 1
            endif

            if(rg == 1) print*,"RINGS",n_rings
            do i = 0,n_rings-1
                if(rg == 1) print*,"VERIF",rg,shift,i,I_give_to,I_take_from
                if(rg == i)then
                    if (Tdomain%sComm(I_give_to)%ngllSF>0) then
                        call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForcesSF_StoF,Tdomain%sComm(I_give_to)%ngllSF, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                    if (Tdomain%sComm(I_take_from)%ngllSF > 0) then
                        call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForcesSF_StoF, Tdomain%sComm(I_take_from)%ngllSF, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                else
                    do j = 0,n/n_rings-1
                        if (rg == i + j*n_rings) then
                            if (Tdomain%sComm(I_take_from)%ngllSF>0) then
                                call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForcesSF_StoF,Tdomain%sComm(I_take_from)%ngllSF, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                            endif
                            if (Tdomain%sComm(I_give_to)%ngllSF>0) then
                                call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForcesSF_StoF,Tdomain%sComm(I_give_to)%ngllSF, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                            endif
                        end if
                    enddo
                endif
            enddo
            call MPI_BARRIER (MPI_COMM_WORLD, code)
        end do

        ! assemblage on external GLLs
        do n = 0,Tdomain%n_proc-1
            ngllSF = 0
            call Comm_Forces_FaceSF_StoF(Tdomain,n,ngllSF)
            call Comm_Forces_EdgeSF_StoF(Tdomain,n,ngllSF)
            call Comm_Forces_VertexSF_StoF(Tdomain,n,ngllSF)
        enddo

    end if

    ! now the coupling Solid -> Fluid can be done.
    do nf = 0,Tdomain%SF%SF_n_faces-1
        nnf = Tdomain%SF%SF_Face(nf)%Face(0)
        if(nnf < 0) cycle  ! fluid face not on this proc
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        Tdomain%sFace(nnf)%ForcesFl(1:ngll1-2,1:ngll2-2) = Tdomain%sFace(nnf)%ForcesFl(1:ngll1-2,1:ngll2-2) +  &
            Tdomain%SF%SF_Face(nf)%Vn(1:ngll1-2,1:ngll2-2)
    end do
    do ne = 0,Tdomain%SF%SF_n_edges-1
        nne = Tdomain%SF%SF_Edge(ne)%Edge(0)
        if(nne < 0) cycle  ! fluid edge not on this proc
        ngll = Tdomain%SF%SF_Edge(ne)%ngll
        Tdomain%sEdge(nne)%ForcesFl(1:ngll-2) = Tdomain%sEdge(nne)%ForcesFl(1:ngll-2) +  &
            Tdomain%SF%SF_Edge(ne)%Vn(1:ngll-2)
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(0)
        if(nnv < 0) cycle  ! fluid vertex not on this proc
        Tdomain%sVertex(nnv)%ForcesFl = Tdomain%sVertex(nnv)%ForcesFl + Tdomain%SF%SF_Vertex(nv)%Vn
    end do

end subroutine StoF_coupling
!----------------------------------------------------------------
!----------------------------------------------------------------
subroutine FtoS_coupling(Tdomain,rg)
    ! from fluid to solid: normal times pressure (= -rho . VelPhi)
    use sdomain
    use mpi
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer, intent(in)   :: rg
    integer  :: n,nf,ne,nv,nff,ngll1,ngll2,ngll,i,j,k,nnf,nne,nnv,orient_e,    &
        n_rings,shift,I_give_to,I_take_from,code,which_elem,ngllSF
    integer, dimension(MPI_STATUS_SIZE)   :: statut
    integer, parameter   :: etiquette = 100
    real, dimension(:,:), allocatable  :: PressFace

    ! init.
    do ne = 0,Tdomain%SF%SF_n_edges-1
        Tdomain%SF%SF_edge(ne)%pn = 0d0
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        Tdomain%SF%SF_vertex(nv)%pn = 0d0
    end do


    do nf = 0,Tdomain%SF%SF_n_faces-1
        nff = Tdomain%SF%SF_Face(nf)%Face(0)
        if(nff < 0) cycle   ! fluid face not on this proc
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        allocate(PressFace(0:ngll1-1,0:ngll2-1))
        PressFace = 0d0
        PressFace(1:ngll1-2,1:ngll2-2) = Tdomain%sFace(nff)%VelPhi(:,:)
        ! deassemblage: particle velocity from edges and vertices to faces
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
        ! on faces: normal times pressure = -rho * VelPhi * normal
        do i = 0,2
            Tdomain%SF%SF_Face(nf)%pn(0:ngll1-1,0:ngll2-1,i) = -1d0*PressFace(0:ngll1-1,0:ngll2-1)*  &
                Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,i)*       &
                Tdomain%SF%SF_Face(nf)%density(0:ngll1-1,0:ngll2-1)
        end do
        deallocate(PressFace)
        ! assemblage: (normal times pressure) from face to edges and vertices
        do i = 0,3
            ne = Tdomain%SF%SF_Face(nf)%Near_Edges(i)
            ngll = Tdomain%SF%SF_Edge(ne)%ngll
            orient_e = Tdomain%SF%SF_face(nf)%Near_Edges_Orient(i)
            call Vector_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                Tdomain%SF%SF_Face(nf)%pn,Tdomain%SF%SF_Edge(ne)%pn)
            nv = Tdomain%SF%SF_Face(nf)%Near_Vertices(i)
            call Vector_Face2Vertex(i,ngll1,ngll2,    &
                Tdomain%SF%SF_Face(nf)%pn,Tdomain%SF%SF_Vertex(nv)%pn)
        end do
    end do

    ! now we can exchange values for SF sides on different procs
    if(Tdomain%n_proc > 1)then
        do n = 0,Tdomain%n_proc-1
            call Comm_Forces_Complete_FtoS(n,Tdomain)
        end do
        ! now we can exchange force values with proc n
        n = Tdomain%n_proc
        do shift = 1,n-1
            I_give_to = rg + shift
            if (I_give_to > n-1)   I_give_to = I_give_to - n
            I_take_from = rg - shift
            if (I_take_from < 0)   I_take_from = I_take_from + n
            if (mod(n,shift)==0 .and. shift/=1) then
                n_rings = shift
            else if (mod(n,n-shift)==0 .and. shift/=n-1) then
                n_rings = n-shift
            else if (mod(n,2)==0 .and. mod(shift,2)==0) then
                n_rings = 2
            else
                n_rings = 1
            endif
            do i = 0,n_rings-1
                if(rg == i)then
                    if (Tdomain%sComm(I_give_to)%ngllSF>0) then
                        call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForcesSF_FtoS,Tdomain%sComm(I_give_to)%ngllSF, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                    if (Tdomain%sComm(I_take_from)%ngllSF > 0) then
                        call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForcesSF_FtoS, Tdomain%sComm(I_take_from)%ngllSF, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                else
                    do j = 0,n/n_rings-1
                        if (rg == i + j*n_rings) then
                            if (Tdomain%sComm(I_take_from)%ngllSF>0) then
                                call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForcesSF_FtoS,Tdomain%sComm(I_take_from)%ngllSF, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                            endif
                            if (Tdomain%sComm(I_give_to)%ngllSF>0) then
                                call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForcesSF_FtoS,Tdomain%sComm(I_give_to)%ngllSF, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                            endif
                        end if
                    enddo
                endif
            enddo
            call MPI_BARRIER (MPI_COMM_WORLD, code)
        end do

        ! assemblage on external GLLs
        do n = 0,Tdomain%n_proc-1
            ngllSF = 0
            call Comm_Forces_FaceSF_FtoS(Tdomain,n,ngllSF)
            call Comm_Forces_EdgeSF_FtoS(Tdomain,n,ngllSF)
            call Comm_Forces_VertexSF_FtoS(Tdomain,n,ngllSF)
        enddo

    end if

    ! now the coupling Fluid -> Solid can be done.
    do nf = 0,Tdomain%SF%SF_n_faces-1
        nnf = Tdomain%SF%SF_Face(nf)%Face(1)
        if(nnf < 0) cycle  ! solid face not on this proc
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        Tdomain%sFace(nnf)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%SF%SF_Face(nf)%save_forces(1:ngll1-2,1:ngll2-2,0:2) +  &
            Tdomain%SF%SF_Face(nf)%pn(1:ngll1-2,1:ngll2-2,0:2)
        if(.not. Tdomain%sFace(nnf)%PML)then
            Tdomain%sFace(nnf)%Displ(:,:,:) = Tdomain%SF%SF_Face(nf)%save_displ(:,:,:)
            Tdomain%sFace(nnf)%Accel(:,:,:) = Tdomain%SF%SF_Face(nf)%save_accel(:,:,:)
        end if
    end do
    do ne = 0,Tdomain%SF%SF_n_edges-1
        nne = Tdomain%SF%SF_Edge(ne)%Edge(1)
        if(nne < 0) cycle  ! solid edge not on this proc
        ngll = Tdomain%SF%SF_Edge(ne)%ngll
        Tdomain%sEdge(nne)%Forces(1:ngll-2,0:2) = Tdomain%SF%SF_Edge(ne)%save_forces(1:ngll-2,0:2) +  &
            Tdomain%SF%SF_Edge(ne)%pn(1:ngll-2,0:2)
        if(.not. Tdomain%sEdge(nne)%PML)then
            Tdomain%sEdge(nne)%Displ(:,:) = Tdomain%SF%SF_Edge(ne)%save_displ(:,:)
            Tdomain%sEdge(nne)%Accel(:,:) = Tdomain%SF%SF_Edge(ne)%save_accel(:,:)
        end if
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(1)
        if(nnv < 0) cycle  ! solid vertex not on this proc
        Tdomain%sVertex(nnv)%Forces(0:2) = Tdomain%SF%SF_Vertex(nv)%save_forces(0:2) + Tdomain%SF%SF_Vertex(nv)%pn(0:2)
        if(.not. Tdomain%sVertex(nnv)%PML)then
            Tdomain%sVertex(nnv)%Displ(:) = Tdomain%SF%SF_Vertex(nv)%save_displ(:)
            Tdomain%sVertex(nnv)%Accel(:) = Tdomain%SF%SF_Vertex(nv)%save_accel(:)
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


    print*,"YAALLLAH"

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

    print*,"YAALLLAH222"
    return
end subroutine Comm_Forces_Complete_StoF
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
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine Comm_Forces_FaceSF_StoF(Tdomain,n,ngllSF)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF

    integer :: ngll1,ngll2,i,j,k,nf


    do i = 0,Tdomain%sComm(n)%SF_nf_shared-1
        nf = Tdomain%sComm(n)%SF_faces_shared(i)
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2

        if(Tdomain%SF%SF_Face(nf)%Orient_Face == 0)then
            do j = 1,ngll2-2
                do k = 1,ngll1-2
                    Tdomain%SF%SF_Face(nf)%Vn(k,j) = Tdomain%SF%SF_Face(nf)%Vn(k,j) + Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 1)then
            do j = 1,ngll2-2
                do k = 1,ngll1-2
                    Tdomain%SF%SF_Face(nf)%Vn(ngll1-1-k,j) = Tdomain%SF%SF_Face(nf)%Vn(ngll1-1-k,j) + Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 2)then
            do j = 1,ngll2-2
                do k = 1,ngll1-2
                    Tdomain%SF%SF_Face(nf)%Vn(k,ngll2-1-j) = Tdomain%SF%SF_Face(nf)%Vn(k,ngll2-1-j) + Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 3)then
            do j = 1,ngll2-2
                do k = 1,ngll1-2
                    Tdomain%SF%SF_Face(nf)%Vn(ngll1-1-k,ngll2-1-j) = Tdomain%SF%SF_Face(nf)%Vn(ngll1-1-k,ngll2-1-j) +   &
                        Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 4)then
            do j = 1,ngll1-2
                do k = 1,ngll2-2
                    Tdomain%SF%SF_Face(nf)%Vn(j,k) = Tdomain%SF%SF_Face(nf)%Vn(j,k) + Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 5)then
            do j = 1,ngll1-2
                do k = 1,ngll2-2
                    Tdomain%SF%SF_Face(nf)%Vn(ngll1-1-j,k) = Tdomain%SF%SF_Face(nf)%Vn(ngll2-1-j,k) + Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 6)then
            do j = 1,ngll1-2
                do k = 1,ngll2-2
                    Tdomain%SF%SF_Face(nf)%Vn(j,ngll2-1-k) = Tdomain%SF%SF_Face(nf)%Vn(j,ngll2-1-k) + Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 7)then
            do j = 1,ngll1-2
                do k = 1,ngll2-2
                    Tdomain%SF%SF_Face(nf)%Vn(ngll1-1-j,ngll2-1-k) = Tdomain%SF%SF_Face(nf)%Vn(ngll1-1-j,ngll2-1-k) +   &
                        Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                    ngllSF = ngllSF + 1
                enddo
            enddo
        endif

    enddo


    return
end subroutine Comm_Forces_FaceSF_StoF
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine Comm_Forces_FaceSF_FtoS(Tdomain,n,ngllSF)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF

    integer :: ngll1,ngll2,i,j,k,nf


    do i = 0,Tdomain%sComm(n)%SF_nf_shared-1
        nf = Tdomain%sComm(n)%SF_faces_shared(i)
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2

        if(Tdomain%SF%SF_Face(nf)%Orient_Face == 0)then
            do j = 1,ngll2-2
                do k = 1,ngll1-2
                    Tdomain%SF%SF_Face(nf)%pn(k,j,0:2) = Tdomain%SF%SF_Face(nf)%pn(k,j,0:2) + Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 1)then
            do j = 1,ngll2-2
                do k = 1,ngll1-2
                    Tdomain%SF%SF_Face(nf)%pn(ngll1-1-k,j,0:2) = Tdomain%SF%SF_Face(nf)%pn(ngll1-1-k,j,0:2) +     &
                        Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 2)then
            do j = 1,ngll2-2
                do k = 1,ngll1-2
                    Tdomain%SF%SF_Face(nf)%pn(k,ngll2-1-j,0:2) = Tdomain%SF%SF_Face(nf)%pn(k,ngll2-1-j,0:2) +     &
                        Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 3)then
            do j = 1,ngll2-2
                do k = 1,ngll1-2
                    Tdomain%SF%SF_Face(nf)%pn(ngll1-1-k,ngll2-1-j,0:2) = Tdomain%SF%SF_Face(nf)%pn(ngll1-1-k,ngll2-1-j,0:2) +   &
                        Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 4)then
            do j = 1,ngll1-2
                do k = 1,ngll2-2
                    Tdomain%SF%SF_Face(nf)%pn(j,k,0:2) = Tdomain%SF%SF_Face(nf)%pn(j,k,0:2) + Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 5)then
            do j = 1,ngll1-2
                do k = 1,ngll2-2
                    Tdomain%SF%SF_Face(nf)%pn(ngll1-1-j,k,0:2) = Tdomain%SF%SF_Face(nf)%pn(ngll2-1-j,k,0:2) +   &
                        Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 6)then
            do j = 1,ngll1-2
                do k = 1,ngll2-2
                    Tdomain%SF%SF_Face(nf)%pn(j,ngll2-1-k,0:2) = Tdomain%SF%SF_Face(nf)%pn(j,ngll2-1-k,0:2) +   &
                        Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                    ngllSF = ngllSF + 1
                enddo
            enddo

        else if(Tdomain%SF%SF_Face(nf)%Orient_Face == 7)then
            do j = 1,ngll1-2
                do k = 1,ngll2-2
                    Tdomain%SF%SF_Face(nf)%pn(ngll1-1-j,ngll2-1-k,0:2) = Tdomain%SF%SF_Face(nf)%pn(ngll1-1-j,ngll2-1-k,0:2) +   &
                        Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                    ngllSF = ngllSF + 1
                enddo
            enddo
        endif

    enddo


    return
end subroutine Comm_Forces_FaceSF_FtoS
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
subroutine Comm_Forces_EdgeSF_StoF(Tdomain,n,ngllSF)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF
    integer :: ngll,i,j,ne


    do i = 0,Tdomain%sComm(n)%SF_ne_shared-1
        ne = Tdomain%sComm(n)%SF_edges_shared(i)
        ngll = Tdomain%SF%SF_Edge(ne)%ngll

        if(Tdomain%SF%SF_Edge(ne)%orient_edge == 0)then
            do j = 1,ngll-2
                Tdomain%SF%SF_Edge(ne)%Vn(j) = Tdomain%SF%SF_Edge(ne)%Vn(j) + Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                ngllSF = ngllSF + 1
            enddo

        else
            do j = 1,ngll-2
                Tdomain%SF%SF_Edge(ne)%Vn(ngll-1-j) = Tdomain%SF%SF_Edge(ne)%Vn(ngll-1-j) +    &
                    Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
                ngllSF = ngllSF + 1
            enddo

        endif

    enddo

    return
end subroutine Comm_Forces_EdgeSF_StoF
!------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------
subroutine Comm_Forces_EdgeSF_FtoS(Tdomain,n,ngllSF)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF
    integer :: ngll,i,j,ne


    do i = 0,Tdomain%sComm(n)%SF_ne_shared-1
        ne = Tdomain%sComm(n)%SF_edges_shared(i)
        ngll = Tdomain%SF%SF_Edge(ne)%ngll

        if(Tdomain%SF%SF_Edge(ne)%orient_edge == 0)then
            do j = 1,ngll-2
                Tdomain%SF%SF_Edge(ne)%pn(j,0:2) = Tdomain%SF%SF_Edge(ne)%pn(j,0:2) + Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                ngllSF = ngllSF + 1
            enddo

        else
            do j = 1,ngll-2
                Tdomain%SF%SF_Edge(ne)%pn(ngll-1-j,0:2) = Tdomain%SF%SF_Edge(ne)%pn(ngll-1-j,0:2) +    &
                    Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
                ngllSF = ngllSF + 1
            enddo

        endif

    enddo

    return
end subroutine Comm_Forces_EdgeSF_FtoS
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine Comm_Forces_VertexSF_StoF(Tdomain,n,ngllSF)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF
    integer :: i,nv


    do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
        nv =  Tdomain%sComm(n)%SF_vertices_shared(i)
        Tdomain%SF%SF_Vertex(nv)%Vn = Tdomain%SF%SF_Vertex(nv)%Vn + Tdomain%sComm(n)%TakeForcesSF_StoF(ngllSF)
        ngllSF = ngllSF + 1
    enddo


    return
end subroutine Comm_Forces_VertexSF_StoF
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
subroutine Comm_Forces_VertexSF_FtoS(Tdomain,n,ngllSF)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngllSF
    integer :: i,nv


    do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
        nv =  Tdomain%sComm(n)%SF_vertices_shared(i)
        Tdomain%SF%SF_Vertex(nv)%pn(0:2) = Tdomain%SF%SF_Vertex(nv)%pn(0:2) + Tdomain%sComm(n)%TakeForcesSF_FtoS(ngllSF,0:2)
        ngllSF = ngllSF + 1
    enddo


    return
end subroutine Comm_Forces_VertexSF_FtoS
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine SF_solid_values_saving(Tdomain)
    ! saves values of fields Forces, Displ and Accel on the solid side of a SF object -
    !    for they are flushed in the 1st correction phase
    use sdomain
    implicit none

    type(domain), intent(inout)   ::  Tdomain
    integer  :: nf,ne,nv,nnf,nne,nnv

    do nf = 0,Tdomain%SF%SF_n_faces-1
        nnf = Tdomain%SF%SF_Face(nf)%Face(1)
        if(nnf < 0) cycle  ! solid side not on this proc
        Tdomain%SF%SF_Face(nf)%save_forces(:,:,:) = Tdomain%sFace(nnf)%Forces(:,:,:)
        Tdomain%SF%SF_Face(nf)%save_displ(:,:,:) = Tdomain%sFace(nnf)%Displ(:,:,:)
        Tdomain%SF%SF_Face(nf)%save_accel(:,:,:) = Tdomain%sFace(nnf)%Accel(:,:,:)
    end do
    do ne = 0,Tdomain%SF%SF_n_edges-1
        nne = Tdomain%SF%SF_Edge(ne)%Edge(1)
        if(nne < 0) cycle  ! solid side not on this proc
        Tdomain%SF%SF_Edge(ne)%save_forces(:,:) = Tdomain%sEdge(nne)%Forces(:,:)
        Tdomain%SF%SF_Edge(ne)%save_displ(:,:) = Tdomain%sEdge(nne)%Displ(:,:)
        Tdomain%SF%SF_Edge(ne)%save_accel(:,:) = Tdomain%sEdge(nne)%Accel(:,:)
    end do
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(1)
        if(nnv < 0) cycle  ! solid side not on this proc
        Tdomain%SF%SF_Vertex(nv)%save_forces(:) = Tdomain%sVertex(nnv)%Forces(:)
        Tdomain%SF%SF_Vertex(nv)%save_displ(:) = Tdomain%sVertex(nnv)%Displ(:)
        Tdomain%SF%SF_Vertex(nv)%save_accel(:) = Tdomain%sVertex(nnv)%Accel(:)
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
                call Correction_Face_Veloc(Tdomain%sface(nnf),bega,gam1,dt)
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
                call Correction_Edge_Veloc(Tdomain%sEdge(nne),bega,gam1,dt)
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
                call Correction_Vertex_Veloc(Tdomain%sVertex(nnv),bega,gam1,dt)
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
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
