!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
subroutine normal_face_weighting(dir,ngllx,nglly,ngllz,ngll1,ngll2,normal,   &
    GLLwx,GLLwy,GLLwz,BtN)
    !- determination of the weighted normal term for a given (inter)face
    implicit none

    integer, intent(in)  :: dir,ngll1,ngll2,ngllx,nglly,ngllz
    real, dimension(0:ngllx-1), intent(in)  :: GLLwx
    real, dimension(0:nglly-1), intent(in)  :: GLLwy
    real, dimension(0:ngllz-1), intent(in)  :: GLLwz
    real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(in)  :: normal
    real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(out) :: BtN
    integer :: i,j

    if(dir == 0 .or. dir == 5)then
        do j = 0,ngll2-1
            do i = 0,ngll1-1
                Btn(i,j,0:2) = GLLwx(i)*GLLwy(j)*normal(i,j,0:2)
            enddo
        enddo
    else if(dir == 1 .or. dir == 3)then
        do j = 0,ngll2-1
            do i = 0,ngll1-1
                Btn(i,j,0:2) = GLLwx(i)*GLLwz(j)*normal(i,j,0:2)
            enddo
        enddo
    else
        do j = 0,ngll2-1
            do i = 0,ngll1-1
                Btn(i,j,0:2) = GLLwy(i)*GLLwz(j)*normal(i,j,0:2)
            enddo
        enddo
    endif

    return

end subroutine normal_face_weighting
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Vector_Face2Edge(ne,ngll1,ngll2,ngll,orient_e,BtNf,BtNe)
    !- distributes vector values BtN, from face to edges.
    implicit none

    integer, intent(in)   :: ne,ngll1,ngll2,ngll,orient_e
    real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(in)  :: BtNf
    real, dimension(1:ngll-2,0:2), intent(inout)  :: BtNe
    integer  :: j

    if(orient_e == 0)then
        select case(ne)
        case(0)
            BtNe(1:ngll-2,0:2) = BtNe(1:ngll-2,0:2)+BtNf(1:ngll1-2,0,0:2)
        case(1)
            BtNe(1:ngll-2,0:2) = BtNe(1:ngll-2,0:2)+BtNf(ngll1-1,1:ngll2-2,0:2)
        case(2)
            BtNe(1:ngll-2,0:2) = BtNe(1:ngll-2,0:2)+BtNf(1:ngll1-2,ngll2-1,0:2)
        case(3)
            BtNe(1:ngll-2,0:2) = BtNe(1:ngll-2,0:2)+BtNf(0,1:ngll2-2,0:2)
        end select
    else
        select case(ne)
        case(0)
            do j=1,ngll-2
                BtNe(j,0:2) = BtNe(j,0:2)+BtNf(ngll1-1-j,0,0:2)
            enddo
        case(1)
            do j=1,ngll-2
                BtNe(j,0:2) = BtNe(j,0:2)+BtNf(ngll1-1,ngll2-1-j,0:2)
            enddo
        case(2)
            do j=1,ngll-2
                BtNe(j,0:2) = BtNe(j,0:2)+BtNf(ngll1-1-j,ngll2-1,0:2)
            enddo
        case(3)
            do j=1,ngll-2
                BtNe(j,0:2) = BtNe(j,0:2)+BtNf(0,ngll2-1-j,0:2)
            enddo
        end select
    endif

end subroutine Vector_Face2Edge
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Vector_Face2Vertex(nv,ngll1,ngll2,BtNf,BtNv)
    !- distributes vector values (BtN), from face to vertices.
    implicit none

    integer, intent(in)   :: nv,ngll1,ngll2
    real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(in)  :: BtNf
    real, dimension(0:2), intent(inout)  :: BtNv

    select case(nv)
    case(0)
        BtNv(0:2) = BtNv(0:2)+BtNf(0,0,0:2)
    case(1)
        BtNv(0:2) = BtNv(0:2)+BtNf(ngll1-1,0,0:2)
    case(2)
        BtNv(0:2) = BtNv(0:2)+BtNf(ngll1-1,ngll2-1,0:2)
    case(3)
        BtNv(0:2) = BtNv(0:2)+BtNf(0,ngll2-1,0:2)
    end select


end subroutine Vector_Face2Vertex
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Vector_Edge2Face(ne,ngll1,ngll2,ngll,orient_e,BtNf,BtNe)
    !- assemblage vector values BtN, from edge to face.
    implicit none

    integer, intent(in)   :: ne,ngll1,ngll2,ngll,orient_e
    real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(inout)  :: BtNf
    real, dimension(1:ngll-2,0:2), intent(in)  :: BtNe
    integer  :: j

    if(orient_e == 0)then
        select case(ne)
        case(0)
            BtNf(1:ngll1-2,0,0:2) = BtNe(1:ngll-2,0:2)
        case(1)
            BtNf(ngll1-1,1:ngll2-2,0:2) = BtNe(1:ngll-2,0:2)
        case(2)
            BtNf(1:ngll1-2,ngll2-1,0:2) = BtNe(1:ngll-2,0:2)
        case(3)
            BtNf(0,1:ngll2-2,0:2) = BtNe(1:ngll-2,0:2)
        end select
    else
        select case(ne)
        case(0)
            do j=1,ngll-2
                BtNf(ngll1-1-j,0,0:2) = BtNe(j,0:2)
            enddo
        case(1)
            do j=1,ngll-2
                BtNf(ngll1-1,ngll2-1-j,0:2) = BtNe(j,0:2)
            enddo
        case(2)
            do j=1,ngll-2
                BtNf(ngll1-1-j,ngll2-1,0:2) = BtNe(j,0:2)
            enddo
        case(3)
            do j=1,ngll-2
                BtNf(0,ngll2-1-j,0:2) = BtNe(j,0:2)
            enddo
        end select
    endif

end subroutine Vector_Edge2Face
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Vector_Vertex2Face(nv,ngll1,ngll2,BtNf,BtNv)
    !- assemblage of vector values (BtN), from vertex to face.
    implicit none

    integer, intent(in)   :: nv,ngll1,ngll2
    real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(inout)  :: BtNf
    real, dimension(0:2), intent(in)  :: BtNv

    select case(nv)
    case(0)
        BtNf(0,0,0:2) = BtNv(0:2)
    case(1)
        BtNf(ngll1-1,0,0:2) = BtNv(0:2)
    case(2)
        BtNf(ngll1-1,ngll2-1,0:2) = BtNv(0:2)
    case(3)
        BtNf(0,ngll2-1,0:2) = BtNv(0:2)
    end select


end subroutine Vector_Vertex2Face
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Scalar_Face2Edge(ne,ngll1,ngll2,ngll,orient_e,BtNf,BtNe)
    !- distributes scalar values BtN, from face to edges.
    implicit none

    integer, intent(in)   :: ne,ngll1,ngll2,ngll,orient_e
    real, dimension(0:ngll1-1,0:ngll2-1), intent(in)  :: BtNf
    real, dimension(1:ngll-2), intent(inout)  :: BtNe
    integer  :: j

    if(orient_e == 0)then
        select case(ne)
        case(0)
            BtNe(1:ngll-2) = BtNe(1:ngll-2)+BtNf(1:ngll1-2,0)
        case(1)
            BtNe(1:ngll-2) = BtNe(1:ngll-2)+BtNf(ngll1-1,1:ngll2-2)
        case(2)
            BtNe(1:ngll-2) = BtNe(1:ngll-2)+BtNf(1:ngll1-2,ngll2-1)
        case(3)
            BtNe(1:ngll-2) = BtNe(1:ngll-2)+BtNf(0,1:ngll2-2)
        end select
    else
        select case(ne)
        case(0)
            do j=1,ngll-2
                BtNe(j) = BtNe(j)+BtNf(ngll1-1-j,0)
            enddo
        case(1)
            do j=1,ngll-2
                BtNe(j) = BtNe(j)+BtNf(ngll1-1,ngll2-1-j)
            enddo
        case(2)
            do j=1,ngll-2
                BtNe(j) = BtNe(j)+BtNf(ngll1-1-j,ngll2-1)
            enddo
        case(3)
            do j=1,ngll-2
                BtNe(j) = BtNe(j)+BtNf(0,ngll2-1-j)
            enddo
        end select
    endif

end subroutine Scalar_Face2Edge
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Scalar_Face2Vertex(nv,ngll1,ngll2,BtNf,BtNv)
    !- distributes scalar values (BtN), from face to vertices.
    implicit none

    integer, intent(in)   :: nv,ngll1,ngll2
    real, dimension(0:ngll1-1,0:ngll2-1), intent(in)  :: BtNf
    real, intent(inout)  :: BtNv

    select case(nv)
    case(0)
        BtNv = BtNv+BtNf(0,0)
    case(1)
        BtNv = BtNv+BtNf(ngll1-1,0)
    case(2)
        BtNv = BtNv+BtNf(ngll1-1,ngll2-1)
    case(3)
        BtNv = BtNv+BtNf(0,ngll2-1)
    end select


end subroutine Scalar_Face2Vertex
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Scalar_Edge2Face(ne,ngll1,ngll2,ngll,orient_e,BtNf,BtNe)
    !- assemblage scalar values BtN, from edge to face.
    implicit none

    integer, intent(in)   :: ne,ngll1,ngll2,ngll,orient_e
    real, dimension(0:ngll1-1,0:ngll2-1), intent(inout)  :: BtNf
    real, dimension(1:ngll-2), intent(in)  :: BtNe
    integer  :: j

    if(orient_e == 0)then
        select case(ne)
        case(0)
            BtNf(1:ngll1-2,0) = BtNe(1:ngll-2)
        case(1)
            BtNf(ngll1-1,1:ngll2-2) = BtNe(1:ngll-2)
        case(2)
            BtNf(1:ngll1-2,ngll2-1) = BtNe(1:ngll-2)
        case(3)
            BtNf(0,1:ngll2-2) = BtNe(1:ngll-2)
        end select
    else
        select case(ne)
        case(0)
            do j=1,ngll-2
                BtNf(ngll1-1-j,0) = BtNe(j)
            enddo
        case(1)
            do j=1,ngll-2
                BtNf(ngll1-1,ngll2-1-j) = BtNe(j)
            enddo
        case(2)
            do j=1,ngll-2
                BtNf(ngll1-1-j,ngll2-1) = BtNe(j)
            enddo
        case(3)
            do j=1,ngll-2
                BtNf(0,ngll2-1-j) = BtNe(j)
            enddo
        end select
    endif

end subroutine Scalar_Edge2Face
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Number_Edge2Face(ne,ngll1,ngll2,ngll,orient_e,BtNf,BtNe)
    !- assemblage scalar values BtN, from edge to face.
    implicit none

    integer, intent(in)   :: ne,ngll1,ngll2,ngll,orient_e
    integer, dimension(0:ngll1-1,0:ngll2-1), intent(inout)  :: BtNf
    integer, dimension(1:ngll-2), intent(in)  :: BtNe
    integer  :: j

    if(orient_e == 0)then
        select case(ne)
        case(0)
            BtNf(1:ngll1-2,0) = BtNe(1:ngll-2)
        case(1)
            BtNf(ngll1-1,1:ngll2-2) = BtNe(1:ngll-2)
        case(2)
            BtNf(1:ngll1-2,ngll2-1) = BtNe(1:ngll-2)
        case(3)
            BtNf(0,1:ngll2-2) = BtNe(1:ngll-2)
        end select
    else
        select case(ne)
        case(0)
            do j=1,ngll-2
                BtNf(ngll1-1-j,0) = BtNe(j)
            enddo
        case(1)
            do j=1,ngll-2
                BtNf(ngll1-1,ngll2-1-j) = BtNe(j)
            enddo
        case(2)
            do j=1,ngll-2
                BtNf(ngll1-1-j,ngll2-1) = BtNe(j)
            enddo
        case(3)
            do j=1,ngll-2
                BtNf(0,ngll2-1-j) = BtNe(j)
            enddo
        end select
    endif

end subroutine Number_Edge2Face
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Scalar_Vertex2Face(nv,ngll1,ngll2,BtNf,BtNv)
    !- assemblage of scalar values (BtN), from vertex to face.
    implicit none

    integer, intent(in)   :: nv,ngll1,ngll2
    real, dimension(0:ngll1-1,0:ngll2-1), intent(inout)  :: BtNf
    real, intent(in)  :: BtNv

    select case(nv)
    case(0)
        BtNf(0,0) = BtNv
    case(1)
        BtNf(ngll1-1,0) = BtNv
    case(2)
        BtNf(ngll1-1,ngll2-1) = BtNv
    case(3)
        BtNf(0,ngll2-1) = BtNv
    end select


end subroutine Scalar_Vertex2Face
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Number_Vertex2Face(nv,ngll1,ngll2,BtNf,BtNv)
    !- assemblage of scalar values (BtN), from vertex to face.
    implicit none

    integer, intent(in)   :: nv,ngll1,ngll2
    integer, dimension(0:ngll1-1,0:ngll2-1), intent(inout)  :: BtNf
    integer, intent(in)  :: BtNv

    select case(nv)
    case(0)
        BtNf(0,0) = BtNv
    case(1)
        BtNf(ngll1-1,0) = BtNv
    case(2)
        BtNf(ngll1-1,ngll2-1) = BtNv
    case(3)
        BtNf(0,ngll2-1) = BtNv
    end select


end subroutine Number_Vertex2Face
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine BtN_Face2Edge(ne,ngll1,ngll2,ngll,orient_e,BtNf,BtNe)
    !- distributes BtN (ponderated normal) values, from face to edges.
    implicit none

    integer, intent(in)   :: ne,ngll1,ngll2,ngll,orient_e
    real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(in)  :: BtNf
    real, dimension(1:ngll-2,0:2), intent(inout)  :: BtNe
    integer  :: j

    if(orient_e == 0)then
        select case(ne)
        case(0)
            BtNe(1:ngll-2,0:2) = BtNe(1:ngll-2,0:2)+BtNf(1:ngll1-2,0,0:2)
        case(1)
            BtNe(1:ngll-2,0:2) = BtNe(1:ngll-2,0:2)+BtNf(ngll1-1,1:ngll2-2,0:2)
        case(2)
            BtNe(1:ngll-2,0:2) = BtNe(1:ngll-2,0:2)+BtNf(1:ngll1-2,ngll2-1,0:2)
        case(3)
            BtNe(1:ngll-2,0:2) = BtNe(1:ngll-2,0:2)+BtNf(0,1:ngll2-2,0:2)
        end select
    else
        select case(ne)
        case(0)
            do j=1,ngll-2
                BtNe(j,0:2) = BtNe(j,0:2)+BtNf(ngll1-1-j,0,0:2)
            enddo
        case(1)
            do j=1,ngll-2
                BtNe(j,0:2) = BtNe(j,0:2)+BtNf(ngll1-1,ngll2-1-j,0:2)
            enddo
        case(2)
            do j=1,ngll-2
                BtNe(j,0:2) = BtNe(j,0:2)+BtNf(ngll1-1-j,ngll2-1,0:2)
            enddo
        case(3)
            do j=1,ngll-2
                BtNe(j,0:2) = BtNe(j,0:2)+BtNf(0,ngll2-1-j,0:2)
            enddo
        end select
    endif

end subroutine BtN_Face2Edge
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine BtN_Face2Vertex(nv,ngll1,ngll2,BtNf,BtNv)
    !- distributes BtN (ponderated normal) values, from face to vertices.
    implicit none

    integer, intent(in)   :: nv,ngll1,ngll2
    real, dimension(0:ngll1-1,0:ngll2-1,0:2), intent(in)  :: BtNf
    real, dimension(0:2), intent(inout)  :: BtNv

    select case(nv)
    case(0)
        BtNv(0:2) = BtNv(0:2)+BtNf(0,0,0:2)
    case(1)
        BtNv(0:2) = BtNv(0:2)+BtNf(ngll1-1,0,0:2)
    case(2)
        BtNv(0:2) = BtNv(0:2)+BtNf(ngll1-1,ngll2-1,0:2)
    case(3)
        BtNv(0:2) = BtNv(0:2)+BtNf(0,ngll2-1,0:2)
    end select


end subroutine BtN_Face2Vertex
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine define_FEV_Neumann(Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer :: nf,ne,nv,ngllx,nglly,ngllz,ngll1,ngll2,ngll,mat,i,orient_e
    real, dimension(:,:,:), allocatable :: Store_Btn

    do nf = 0,Tdomain%Neumann%Neu_n_faces-1
        ngll1 = Tdomain%Neumann%Neu_Face(nf)%ngll1 ; ngll2 = Tdomain%Neumann%neu_Face(nf)%ngll2
        mat = Tdomain%Neumann%Neu_Param%mat_index
        ngllx = Tdomain%sSubdomain(mat)%ngllx
        nglly = Tdomain%sSubdomain(mat)%nglly
        ngllz = Tdomain%sSubdomain(mat)%ngllz
        call normal_face_weighting(Tdomain%Neumann%Neu_face(nf)%dir,ngllx,nglly,ngllz,     &
            ngll1,ngll2,Tdomain%Neumann%Neu_face(nf)%normal,Tdomain%sSubdomain(mat)%GLLwx, &
            Tdomain%sSubdomain(mat)%GLLwy,Tdomain%sSubdomain(mat)%GLLwz,                   &
            Tdomain%Neumann%Neu_Face(nf)%Btn)
        !- internal communication of Btn: from Neumann faces to Neu edges and vertices
        do i = 0,3
            ne = Tdomain%Neumann%Neu_Face(nf)%Near_Edges(i)
            ngll = Tdomain%Neumann%Neu_Edge(ne)%ngll
            orient_e = Tdomain%Neumann%Neu_face(nf)%Near_Edges_Orient(i)
            call BtN_Face2Edge(i,ngll1,ngll2,ngll,orient_e,         &
                Tdomain%Neumann%Neu_Face(nf)%BtN,Tdomain%Neumann%Neu_Edge(ne)%BtN)
            nv = Tdomain%Neumann%Neu_Face(nf)%Near_Vertices(i)
            call BtN_Face2Vertex(i,ngll1,ngll2,    &
                Tdomain%Neumann%Neu_Face(nf)%BtN,Tdomain%Neumann%Neu_Vertex(nv)%BtN)
        end do

        !- changing size of face arrays.
        allocate(Store_Btn(1:ngll1-2,1:ngll2-2,0:2))
        Store_Btn(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%Neumann%Neu_Face(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2)
        deallocate(Tdomain%Neumann%Neu_Face(nf)%Btn)
        allocate(Tdomain%Neumann%Neu_Face(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2))
        Tdomain%Neumann%Neu_Face(nf)%Btn(1:ngll1-2,1:ngll2-2,0:2) = Store_Btn(1:ngll1-2,1:ngll2-2,0:2)
        deallocate(Store_Btn)
        deallocate(Tdomain%Neumann%Neu_Face(nf)%normal)

        !- end on the loop on a Neumann face
    enddo

end subroutine define_FEV_Neumann
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine define_Face_SF(Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer :: nf,ngllx,nglly,ngllz,ngll1,ngll2,mat,w_elem,dir

    integer :: i, j, kb
    Tdomain%SF%SF_BtN(:,:) = 0.

    kb = 0
    do nf = 0,Tdomain%SF%SF_n_faces-1
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1 ; ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        if(Tdomain%SF%SF_face(nf)%Face(0) > 0)then
            w_elem = Tdomain%sFace(Tdomain%SF%SF_face(nf)%Face(0))%which_elem
        else
            w_elem = Tdomain%sFace(Tdomain%SF%SF_face(nf)%Face(1))%which_elem
        end if
        mat = Tdomain%specel(w_elem)%mat_index
        ngllx = Tdomain%sSubdomain(mat)%ngllx
        nglly = Tdomain%sSubdomain(mat)%nglly
        ngllz = Tdomain%sSubdomain(mat)%ngllz
        dir = Tdomain%SF%SF_face(nf)%dir
        call normal_face_weighting(dir,ngllx,nglly,ngllz,     &
            ngll1,ngll2,Tdomain%SF%SF_face(nf)%normal,Tdomain%sSubdomain(mat)%GLLwx, &
            Tdomain%sSubdomain(mat)%GLLwy,Tdomain%sSubdomain(mat)%GLLwz,             &
            Tdomain%SF%SF_Face(nf)%Btn)


        ! On calcul le tableau de BtN global
        do j = 0,ngll2-1
           do i = 0,ngll1-1
               kb = Tdomain%SF%SF_Face(nf)%I_sf(i,j)
               Tdomain%SF%SF_BtN(kb,0:2) = Tdomain%SF%SF_BtN(kb,0:2) + Tdomain%SF%SF_Face(nf)%Btn(i,j,0:2)
           enddo
        enddo
    enddo

end subroutine define_Face_SF
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
subroutine Number_Face2Face(ngll1,ngll2,orient_e,BtNf,BtNe)
    !- assemblage scalar values BtN, from edge to face.
    implicit none

    integer, intent(in)   :: ngll1,ngll2,orient_e
    integer, dimension(0:ngll1-1,0:ngll2-1), intent(inout)  :: BtNf
    integer, dimension(1:ngll1-2,1:ngll2-2), intent(in)  :: BtNe
    integer  :: j


    select case(orient_e)
        case(0)
            BtNf(1:ngll1-2,1:ngll2-2) = BtNe
        case(1) ! TODO
            BtNf(1:ngll1-2,1:ngll2-2) = BtNe
        case(2) ! TODO
            BtNf(1:ngll1-2,1:ngll2-2) = BtNe
        case(3) ! TODO
            BtNf(1:ngll1-2,1:ngll2-2) = BtNe
        case(4) ! TODO
            BtNf(1:ngll1-2,1:ngll2-2) = BtNe
        case(5) ! TODO
            BtNf(1:ngll1-2,1:ngll2-2) = BtNe
        case(6) ! TODO
            BtNf(1:ngll1-2,1:ngll2-2) = BtNe
        case(7) ! TODO
            BtNf(1:ngll1-2,1:ngll2-2) = BtNe
    end select

end subroutine Number_Face2Face

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
