!>
!!\file global_numbering.f90
!!\brief Assure la correspondance entre les différentes numérotations.
!!
!<


!>
!! Définition de Iglobnum et  renvoi du nombre total de ddl: elements, faces, aretes, sommets
!!
!<
subroutine global_numbering(Tdomain,rank)

    ! routine different from the 2D case. Everything is independently numbered, here (inner
    !      points in elements, on faces, edges and vertices). And then associated in
    !      the "Iglobnum" field = global number of each GLL point.


    use sdomain
    use mindex
    implicit none
    integer, intent(in)  :: rank
    type(domain), intent (inout) :: Tdomain
    integer :: n, icount, i, j, k, ngllx, nglly, ngllz, nf, nnf, ne, nne, nv, ngll1, ngll2,   &
        orient_f, orient_e, ngll, nnv, idxS, idxF, idxSpml, idxFpml
    integer, dimension(0:6)  :: index_elem_f
    integer, dimension(0:4)  :: index_elem_e
    integer, dimension(0:2)  :: index_elem_v
    integer, dimension (:), allocatable :: renumS, renumF, renumSpml, renumFpml

    interface
        subroutine renum_solide(Tdomain, n, ngllx, nglly, ngllz, idxS, renumS)
            use sdomain
            implicit none
            type(domain), intent(inout) :: Tdomain
            integer, intent (inout) :: idxS
            integer, dimension (:), allocatable, intent (inout) :: renumS
            integer, intent (in) :: n, ngllx, nglly, ngllz
            integer :: i, j, k, num
        end subroutine renum_solide

        subroutine renum_solide_pml(Tdomain, n, ngllx, nglly, ngllz, idxS, renumS, idxSpml, renumSpml)
            use sdomain
            implicit none
            type(domain), intent(inout) :: Tdomain
            integer, intent (inout) :: idxS, idxSpml
            integer, dimension (:), allocatable, intent (inout) :: renumS, renumSpml
            integer, intent (in) :: n, ngllx, nglly, ngllz
            integer :: i, j, k, num
        end subroutine renum_solide_pml

        subroutine renum_fluide(Tdomain, n, ngllx, nglly, ngllz, idxF, renumF)
            use sdomain
            implicit none
            type(domain), intent(inout) :: Tdomain
            integer, intent (inout) :: idxF
            integer, dimension (:), allocatable, intent (inout) :: renumF
            integer, intent (in) :: n, ngllx, nglly, ngllz
            integer :: i, j, k, num
        end subroutine renum_fluide

        subroutine renum_fluide_pml(Tdomain, n, ngllx, nglly, ngllz, idxF, renumF, idxFpml, renumFpml)
            use sdomain
            implicit none
            type(domain), intent(inout) :: Tdomain
            integer, intent (inout) :: idxF, idxFpml
            integer, dimension (:), allocatable, intent (inout) :: renumF, renumFpml
            integer, intent (in) :: n, ngllx, nglly, ngllz
            integer :: i, j, k, num
        end subroutine renum_fluide_pml
    end interface

    icount = 0

    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        allocate(Tdomain%specel(n)%Iglobnum(0:ngllx-1,0:nglly-1,0:ngllz-1))
        Tdomain%specel(n)%Iglobnum = -1
        do k = 1,ngllz-2
            do j = 1,nglly-2
                do i = 1,ngllx-2
                    Tdomain%specel(n)%Iglobnum(i,j,k) = icount
                    icount = icount + 1
                enddo
            enddo
        enddo
    enddo


    do n = 0,Tdomain%n_face-1
        ngllx = Tdomain%sFace(n)%ngll1
        nglly = Tdomain%sFace(n)%ngll2
        allocate(Tdomain%sFace(n)%Iglobnum_Face(1:ngllx-2,1:nglly-2))
        Tdomain%sFace(n)%Iglobnum_Face = -1
        do j = 1,nglly-2
            do i = 1,ngllx-2
                Tdomain%sFace(n)%Iglobnum_Face(i,j) = icount
                icount = icount + 1
            enddo
        enddo
    enddo

    do n = 0,Tdomain%n_edge-1
        ngllx = Tdomain%sEdge(n)%ngll
        allocate(Tdomain%sEdge(n)%Iglobnum_Edge(1:ngllx-2))
        Tdomain%sEdge(n)%Iglobnum_Edge = -1
        do i = 1,ngllx-2
            Tdomain%sEdge(n)%Iglobnum_Edge(i) = icount
            icount = icount + 1
        enddo
    enddo

    do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%Iglobnum_Vertex = icount
        icount = icount + 1
    enddo

    ! total number of GLL points (= degrees of freedom)
    Tdomain%n_glob_points = icount


    ! recollecting at the element level, from faces, edges and vertices.

    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        ! taking information from faces
        do nf = 0,5
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            call ind_elem_face(nf,orient_f,ngllx,nglly,ngllz,index_elem_f)
            select case(orient_f)
            case(0,1,2,3)
                if(nf == 2 .or. nf == 4)then
                    Tdomain%specel(n)%Iglobnum(index_elem_f(0),                       &
                        index_elem_f(1):index_elem_f(2):index_elem_f(3),    &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6)) =  &
                        Tdomain%sFace(nnf)%Iglobnum_Face(1:ngll1-2,1:ngll2-2)
                else if(nf == 1 .or. nf == 3)then
                    Tdomain%specel(n)%Iglobnum(index_elem_f(1):index_elem_f(2):index_elem_f(3), &
                        index_elem_f(0),                                              &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6)) =            &
                        Tdomain%sFace(nnf)%Iglobnum_Face(1:ngll1-2,1:ngll2-2)
                else if(nf == 0 .or. nf == 5)then
                    Tdomain%specel(n)%Iglobnum(index_elem_f(1):index_elem_f(2):index_elem_f(3), &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6),    &
                        index_elem_f(0)) =                                  &
                        Tdomain%sFace(nnf)%Iglobnum_Face(1:ngll1-2,1:ngll2-2)

                end if
            case(4,5,6,7)
                if(nf == 2 .or. nf == 4)then
                    Tdomain%specel(n)%Iglobnum(index_elem_f(0),                                 &
                        index_elem_f(1):index_elem_f(2):index_elem_f(3),               &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6)) =             &
                        TRANSPOSE(Tdomain%sFace(nnf)%Iglobnum_Face(1:ngll1-2,1:ngll2-2))
                else if(nf == 1 .or. nf == 3)then
                    Tdomain%specel(n)%Iglobnum(index_elem_f(1):index_elem_f(2):index_elem_f(3), &
                        index_elem_f(0),                                               &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6)) =             &
                        TRANSPOSE(Tdomain%sFace(nnf)%Iglobnum_Face(1:ngll1-2,1:ngll2-2))
                else if(nf == 0 .or. nf == 5)then
                    Tdomain%specel(n)%Iglobnum(index_elem_f(1):index_elem_f(2):index_elem_f(3), &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6),               &
                        index_elem_f(0)) =                                             &
                        TRANSPOSE(Tdomain%sFace(nnf)%Iglobnum_Face(1:ngll1-2,1:ngll2-2))

                end if

            end select
        end do


        ! taking information from edges
        do ne = 0,11
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            nne = Tdomain%specel(n)%Near_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            call ind_elem_edge(ne,orient_e,ngllx,nglly,ngllz,index_elem_e)
            select case(ne)
            case(1,3,8,11)  ! only y-coordinate does vary
                Tdomain%specel(n)%Iglobnum(index_elem_e(0),                                 &
                    index_elem_e(2):index_elem_e(3):index_elem_e(4),                &
                    index_elem_e(1)) = Tdomain%sEdge(nne)%Iglobnum_Edge(1:ngll-2)
            case(0,2,5,9)   ! only x-coordinate does vary
                Tdomain%specel(n)%Iglobnum(index_elem_e(2):index_elem_e(3):index_elem_e(4), &
                    index_elem_e(0),index_elem_e(1)) =                              &
                    Tdomain%sEdge(nne)%Iglobnum_Edge(1:ngll-2)
            case(4,6,7,10)   ! only z-coordinate does vary
                Tdomain%specel(n)%Iglobnum(index_elem_e(0),index_elem_e(1),                 &
                    index_elem_e(2):index_elem_e(3):index_elem_e(4)) =              &
                    Tdomain%sEdge(nne)%Iglobnum_Edge(1:ngll-2)
            end select
        end do

        ! taking information from vertices
        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            call ind_elem_vertex(nv,ngllx,nglly,ngllz,index_elem_v)
            Tdomain%specel(n)%Iglobnum(index_elem_v(0),index_elem_v(1),index_elem_v(2)) =   &
                Tdomain%sVertex(nnv)%Iglobnum_Vertex
        end do

    enddo    ! end of the loop onto elements


    ! Renumerotation
    allocate(renumS(0:Tdomain%n_glob_points-1))
    allocate(renumF(0:Tdomain%n_glob_points-1))
    allocate(renumSpml(0:Tdomain%n_glob_points-1))
    allocate(renumFpml(0:Tdomain%n_glob_points-1))
    renumS(:) = -1
    renumF(:) = -1
    renumSpml(:) = -1
    renumFpml(:) = -1
    idxS = 0
    idxF = 0
    idxSpml = 0
    idxFpml = 0
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        if (Tdomain%specel(n)%solid)then
            ! Element solide avec ou sans PML
            allocate(Tdomain%specel(n)%ISol(0:ngllx-1,0:nglly-1,0:ngllz-1))
            Tdomain%specel(n)%ISol(:,:,:) = -1

            if (Tdomain%specel(n)%PML) then
                call renum_solide_pml(Tdomain, n, ngllx, nglly, ngllz, idxS, renumS, idxSpml, renumSpml)
            else
                call renum_solide(Tdomain, n, ngllx, nglly, ngllz, idxS, renumS)
            endif ! fin test pml on non

        else
            ! Element fluide avec ou sans PML
            allocate(Tdomain%specel(n)%IFlu(0:ngllx-1,0:nglly-1,0:ngllz-1))
            Tdomain%specel(n)%IFlu(:,:,:) = -1

            if (Tdomain%specel(n)%PML) then
                ! Element fluide PML
                call renum_fluide_pml(Tdomain, n, ngllx, nglly, ngllz, idxF, renumF, idxFpml, renumFpml)
            else
                ! Element fluide
                call renum_fluide(Tdomain, n, ngllx, nglly, ngllz, idxF, renumF)
            endif ! fin test pml ou non
        endif ! fin test solide ou liquide
    enddo ! fin boucle sur les elements

    Tdomain%ngll_s = idxS
    Tdomain%ngll_f = idxF
    Tdomain%ngll_pmls = idxSpml
    Tdomain%ngll_pmlf = idxFpml

    return
end subroutine global_numbering


!! Renumerotation des points de gauss d'element solide
subroutine renum_solide(Tdomain, n, ngllx, nglly, ngllz, idxS, renumS)
    use sdomain

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent (inout) :: idxS
    integer, dimension (:), allocatable, intent (inout) :: renumS
    integer, intent (in) :: n, ngllx, nglly, ngllz
    integer :: i, j, k, num

    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1
                num = renumS(Tdomain%specel(n)%Iglobnum(i,j,k))
                if (num == -1) then
                    renumS(Tdomain%specel(n)%Iglobnum(i,j,k)) = idxS
                    num = idxS
                    idxS = idxS + 1                                
                endif
                Tdomain%specel(n)%ISol(i,j,k) = num
            enddo
        enddo
    enddo

    return
end subroutine renum_solide


!! Renumerotation des points de gauss d'element solide et pml
subroutine renum_solide_pml(Tdomain, n, ngllx, nglly, ngllz, idxS, renumS, idxSpml, renumSpml)
    use sdomain

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent (inout) :: idxS, idxSpml
    integer, dimension (:), allocatable, intent (inout) :: renumS, renumSpml
    integer, intent (in) :: n, ngllx, nglly, ngllz
    integer :: i, j, k, num

    allocate(Tdomain%specel(n)%slpml%IPml(0:ngllx-1,0:nglly-1,0:ngllz-1))
    Tdomain%specel(n)%slpml%IPml(:,:,:) = -1
    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1
                num = renumS(Tdomain%specel(n)%Iglobnum(i,j,k))
                if (num == -1) then
                    renumS(Tdomain%specel(n)%Iglobnum(i,j,k)) = idxS
                    num = idxS
                    idxS = idxS + 1
                endif
                Tdomain%specel(n)%ISol(i,j,k) = num

                num = renumSpml(Tdomain%specel(n)%Iglobnum(i,j,k))
                if (num == -1) then
                    renumSpml(Tdomain%specel(n)%Iglobnum(i,j,k)) = idxSpml
                    num = idxSpml
                    idxSpml = idxSpml + 1
                endif
                Tdomain%specel(n)%slpml%IPml(i,j,k) = num
            enddo
        enddo
    enddo

    return
end subroutine renum_solide_pml

!! Renumerotation des points de gauss d'element fluide
subroutine renum_fluide(Tdomain, n, ngllx, nglly, ngllz, idxF, renumF)
    use sdomain

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent (inout) :: idxF
    integer, dimension (:), allocatable, intent (inout) :: renumF
    integer, intent (in) :: n, ngllx, nglly, ngllz
    integer :: i, j, k, num

    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1
                num = renumF(Tdomain%specel(n)%Iglobnum(i,j,k))
                if (num == -1) then
                    renumF(Tdomain%specel(n)%Iglobnum(i,j,k)) = idxF
                    num = idxF
                    idxF = idxF + 1
                endif
                Tdomain%specel(n)%IFlu(i,j,k) = num
            enddo
        enddo
    enddo

    return
end subroutine renum_fluide


!! Renumerotation des points de gauss d'element fluide et pml
subroutine renum_fluide_pml(Tdomain, n, ngllx, nglly, ngllz, idxF, renumF, idxFpml, renumFpml)
    use sdomain

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent (inout) :: idxF, idxFpml
    integer, dimension (:), allocatable, intent (inout) :: renumF, renumFpml
    integer, intent (in) :: n, ngllx, nglly, ngllz
    integer :: i, j, k, num

    allocate(Tdomain%specel(n)%flpml%IPml(0:ngllx-1,0:nglly-1,0:ngllz-1))
    Tdomain%specel(n)%flpml%IPml(:,:,:) = -1
    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1
                num = renumF(Tdomain%specel(n)%Iglobnum(i,j,k))
                if (num == -1) then
                    renumF(Tdomain%specel(n)%Iglobnum(i,j,k)) = idxF
                    num = idxF
                    idxF = idxF + 1
                endif
                Tdomain%specel(n)%IFlu(i,j,k) = num

                num = renumFpml(Tdomain%specel(n)%Iglobnum(i,j,k))
                if (num == -1) then
                    renumFpml(Tdomain%specel(n)%Iglobnum(i,j,k)) = idxFpml
                    num = idxFpml
                    idxFpml = idxFpml + 1
                endif
                Tdomain%specel(n)%flpml%IPml(i,j,k) = num
            enddo
        enddo
    enddo

    return
end subroutine renum_fluide_pml

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
