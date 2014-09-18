!>
!!\file global_numbering.f90
!!\brief Assure la correspondance entre les différentes numérotations.
!!
!<


!>
!! Définition de Iglobnum et  renvoi du nombre total de ddl: elements, faces, aretes, sommets
!!
!<
module mrenumber
  public :: global_numbering
  private :: populate_index_SF,renum_element, renum_face, renum_edge,renum_vertex
contains
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
        orient_f, orient_e, ngll, nnv
    integer, dimension(0:6)  :: index_elem_f
    integer, dimension(0:4)  :: index_elem_e
    integer, dimension(0:2)  :: index_elem_v

#if NEW_GLOBAL_METHOD
    integer :: idx, idxS, idxF, idxSpml, idxFpml, dir, ks, kl
    integer :: nbPtInterfSolPml, abscount, nproc
    integer, dimension (:), allocatable :: renumS, renumF, renumSpml, renumFpml
    integer :: nsol, nsolpml, nflu, nflupml

#endif


    icount = 0
#if NEW_GLOBAL_METHOD
    abscount = 0
#endif

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
#if NEW_GLOBAL_METHOD
        if (Tdomain%sFace(n)%Abs .and. Tdomain%sFace(n)%PML) then
            abscount = abscount + (ngllx-2) * (nglly-2)
        endif
#endif
    enddo

    do n = 0,Tdomain%n_edge-1
        ngllx = Tdomain%sEdge(n)%ngll
        allocate(Tdomain%sEdge(n)%Iglobnum_Edge(1:ngllx-2))
        Tdomain%sEdge(n)%Iglobnum_Edge = -1
        do i = 1,ngllx-2
            Tdomain%sEdge(n)%Iglobnum_Edge(i) = icount
            icount = icount + 1
        enddo
#if NEW_GLOBAL_METHOD
        if (Tdomain%sEdge(n)%Abs .and. Tdomain%sEdge(n)%PML) then
            abscount = abscount + ngllx - 2
        endif
#endif
    enddo

    do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%Iglobnum_Vertex = icount
        icount = icount + 1
#if NEW_GLOBAL_METHOD
        if (Tdomain%sVertex(n)%Abs .and. Tdomain%sVertex(n)%PML) then
            abscount = abscount + 1
        endif
#endif
    enddo

#if NEW_GLOBAL_METHOD
    Tdomain%nbOuterPMLNodes = abscount
    allocate(Tdomain%OuterPMLNodes(0:Tdomain%nbOuterPMLNodes-1))
#endif

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

#if NEW_GLOBAL_METHOD
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
                allocate(Tdomain%specel(n)%slpml)
                allocate(Tdomain%specel(n)%slpml%ISolPml(0:ngllx-1,0:nglly-1,0:ngllz-1))
                Tdomain%specel(n)%slpml%ISolPml(:,:,:) = -1
                call renum_element(ngllx, nglly, ngllz, Tdomain%specel(n)%Iglobnum, &
                                   Tdomain%specel(n)%slpml%ISolPml, idxSpml, renumSpml, .true.)
            else
                call renum_element(ngllx, nglly, ngllz, Tdomain%specel(n)%Iglobnum, &
                                   Tdomain%specel(n)%ISol, idxS, renumS, .false.)
            endif ! fin test pml on non

        else
            ! Element fluide avec ou sans PML
            allocate(Tdomain%specel(n)%IFlu(0:ngllx-1,0:nglly-1,0:ngllz-1))
            Tdomain%specel(n)%IFlu(:,:,:) = -1

            if (Tdomain%specel(n)%PML) then
                ! Element fluide PML
                allocate(Tdomain%specel(n)%flpml)
                allocate(Tdomain%specel(n)%flpml%IFluPml(0:ngllx-1,0:nglly-1,0:ngllz-1))
                Tdomain%specel(n)%flpml%IFluPml(:,:,:) = -1
                call renum_element(ngllx, nglly, ngllz, Tdomain%specel(n)%Iglobnum, &
                                   Tdomain%specel(n)%flpml%IFluPml, idxFpml, renumFpml, .true.)
            else
                ! Element fluide
                call renum_element(ngllx, nglly, ngllz, Tdomain%specel(n)%Iglobnum, &
                                   Tdomain%specel(n)%IFlu, idxF, renumF, .false.)
            endif ! fin test pml ou non
        endif ! fin test solide ou liquide
    enddo ! fin boucle sur les elements

    ! Faces
    do n = 0,Tdomain%n_face-1
        ngllx = Tdomain%sFace(n)%ngll1
        nglly = Tdomain%sFace(n)%ngll2
        if (Tdomain%sFace(n)%solid)then
            if (Tdomain%sFace(n)%PML) then
                call renum_face(ngllx, nglly, Tdomain%sFace(n)%Iglobnum_Face, &
                                idxSpml, renumSpml, .true.)
            else
                call renum_face(ngllx, nglly, Tdomain%sFace(n)%Iglobnum_Face, &
                                idxS, renumS, .false.)
            endif
        else
            if (Tdomain%sFace(n)%PML) then
                call renum_face(ngllx, nglly, Tdomain%sFace(n)%Iglobnum_Face, &
                                idxFpml, renumFpml, .true.)
            else
                call renum_face(ngllx, nglly, Tdomain%sFace(n)%Iglobnum_Face, &
                                idxF, renumF, .false.)
            endif
        endif
    enddo

    ! Edges
    do n = 0,Tdomain%n_edge-1
        ngllx = Tdomain%sEdge(n)%ngll
        if (Tdomain%sEdge(n)%solid)then
            if (Tdomain%sEdge(n)%PML) then
                call renum_edge(ngllx, Tdomain%sEdge(n)%Iglobnum_Edge, &
                                idxSpml, renumSpml, .true.)
            else
                call renum_edge(ngllx, Tdomain%sEdge(n)%Iglobnum_Edge, &
                                idxS, renumS, .false.)
            endif
        else
            if (Tdomain%sEdge(n)%PML) then
                call renum_edge(ngllx, Tdomain%sEdge(n)%Iglobnum_Edge, &
                                idxFpml, renumFpml, .true.)
            else
                call renum_edge(ngllx, Tdomain%sEdge(n)%Iglobnum_Edge, &
                                idxF, renumF, .false.)
            endif
        endif
    enddo

    ! Vertices
    do n = 0,Tdomain%n_vertex-1
        if (Tdomain%sVertex(n)%solid)then
            if (Tdomain%sVertex(n)%PML) then
                call renum_vertex(Tdomain%sVertex(n)%Iglobnum_Vertex, &
                                  idxSpml, renumSpml, .true.)
            else
                call renum_vertex(Tdomain%sVertex(n)%Iglobnum_Vertex, &
                                  idxS, renumS, .false.)
            endif
        else
            if (Tdomain%sVertex(n)%PML) then
                call renum_vertex(Tdomain%sVertex(n)%Iglobnum_Vertex, &
                                  idxFpml, renumFpml, .true.)
            else
                call renum_vertex(Tdomain%sVertex(n)%Iglobnum_Vertex, &
                                  idxF, renumF, .false.)
            endif
        endif
    enddo

    if (Tdomain%any_PML) then
        nbPtInterfSolPml = 0
        ! On cré l'interface de couplage Sol/PML
        do n = 0,Tdomain%n_glob_points-1
            if (renumS(n) > -1 .and. renumSpml(n) > -1) then
                ! On se trouve à l'interface Sol/PML : on compte le nombre de points à l'interface
                nbPtInterfSolPml = nbPtInterfSolPml + 1
            endif
        enddo
        Tdomain%nbInterfSolPml = nbPtInterfSolPml
        allocate(Tdomain%InterfSolPml(0:nbPtInterfSolPml-1,0:1))
        i = 0
        do n = 0,Tdomain%n_glob_points-1
            if (renumS(n) > -1 .and. renumSpml(n) > -1) then
                ! On se trouve à l'interface Sol/PML :
                Tdomain%InterfSolPml(i,0) = renumS(n)
                Tdomain%InterfSolPml(i,1) = renumSpml(n)
                !write(*,*) "n,idx sol, idx pml : ", n, renumS(n), renumSpml(n)
                i = i + 1
            endif
        enddo
    endif

    Tdomain%ngll_s = idxS
    Tdomain%ngll_f = idxF
    Tdomain%ngll_pmls = idxSpml
    Tdomain%ngll_pmlf = idxFpml

    !! Couplage solide/fluide
    !! On liste les numeraux globaux des points de gauss concernant les faces de couplage S/F
    if(Tdomain%logicD%SF_local_present)then
        ! On compte le nombre de points de gauss pour les faces de couplages
        k = 0
        do nf = 0,Tdomain%SF%SF_n_faces-1
            nnf = Tdomain%SF%SF_Face(nf)%Face(1)
            if(nnf < 0) cycle
            n = Tdomain%sFace(nnf)%which_elem
            dir = Tdomain%sFace(nnf)%dir
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            if(dir == 0 .OR. dir == 5)then
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        k = k + 1
                    enddo
                enddo
            else if(dir == 1 .OR. dir == 3)then
                do j = 0,ngllz-1
                    do i = 0,ngllx-1
                        k = k + 1
                    enddo
                enddo
            else
                do j = 0,ngllz-1
                    do i = 0,nglly-1
                        k = k + 1
                    enddo
                enddo
            endif
        enddo

        ! On alloue et on initialise les tableaux d'indice des points de couplage
        Tdomain%SF%ngll = k
        allocate(Tdomain%SF%SF_IGlobSol(0:k-1))
        allocate(Tdomain%SF%SF_IGlobFlu(0:k-1))
        Tdomain%SF%SF_IGlobSol = -1
        Tdomain%SF%SF_IGlobFlu = -1

        ! On peuple le tableau d'indice des ngll de couplage en fonction des direcions des faces
        ks = 0
        kl = 0
        k = 0

        do nf = 0,Tdomain%SF%SF_n_faces-1
            ! Partie fluide
            nnf = Tdomain%SF%SF_Face(nf)%Face(0)
            if(nnf > -1) then
                n = Tdomain%sFace(nnf)%which_elem
                dir = Tdomain%sFace(nnf)%dir
                ngllx = Tdomain%specel(n)%ngllx
                nglly = Tdomain%specel(n)%nglly
                ngllz = Tdomain%specel(n)%ngllz

                call populate_index_SF(Tdomain%SF%ngll, dir, ngllx, nglly, ngllz, &
                                       Tdomain%specel(n)%IFlu, kl, Tdomain%SF%SF_IGlobFlu)
            end if

            ! Partie Solide
            nnf = Tdomain%SF%SF_Face(nf)%Face(1)
            if(nnf > -1) then
                n = Tdomain%sFace(nnf)%which_elem
                dir = Tdomain%sFace(nnf)%dir
                ngllx = Tdomain%specel(n)%ngllx
                nglly = Tdomain%specel(n)%nglly
                ngllz = Tdomain%specel(n)%ngllz

                call populate_index_SF(Tdomain%SF%ngll, dir, ngllx, nglly, ngllz, &
                                       Tdomain%specel(n)%ISol, ks, Tdomain%SF%SF_IGlobSol)
            end if
        enddo

    endif ! fin traitement couplage Solide / Fluide

    abscount = 0
    do n = 0,Tdomain%n_face-1
        ngllx = Tdomain%sFace(n)%ngll1
        nglly = Tdomain%sFace(n)%ngll2
        if (Tdomain%sFace(n)%Abs .and. Tdomain%sFace(n)%PML) then
            do j = 1,nglly-2
                do i = 1,ngllx-2
                    idxSpml = renumSpml(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                    if (idxSpml .EQ. -1) stop "Unexpected non pml outer Face !"
                    Tdomain%OuterPMLNodes(abscount) = idxSpml
                    abscount = abscount + 1
                enddo
            enddo
        endif

        allocate(Tdomain%sFace(n)%Renum(1:ngllx-2,1:nglly-2))
        if (Tdomain%sFace(n)%solid) then
            if (Tdomain%sFace(n)%PML) then
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        Tdomain%sFace(n)%Renum(i,j) = renumSpml(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (Tdomain%sFace(n)%Renum(i,j) < 0) stop "Error renumbering SOLID PML faces"
                    enddo
                enddo
            else
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        Tdomain%sFace(n)%Renum(i,j) = renumS(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (Tdomain%sFace(n)%Renum(i,j) < 0) stop "Error renumbering SOLID faces"
                    enddo
                enddo
            endif
        else
            if (Tdomain%sFace(n)%PML) then
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        Tdomain%sFace(n)%Renum(i,j) = renumFpml(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (Tdomain%sFace(n)%Renum(i,j) < 0) stop "Error renumbering FLUID PML faces"
                    enddo
                enddo
            else
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        Tdomain%sFace(n)%Renum(i,j) = renumF(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (Tdomain%sFace(n)%Renum(i,j) < 0) stop "Error renumbering FLUID faces"
                    enddo
                enddo
            endif
        endif
    enddo

    do n = 0,Tdomain%n_edge-1
        ngllx = Tdomain%sEdge(n)%ngll
        if (Tdomain%sEdge(n)%Abs .and. Tdomain%sEdge(n)%PML) then
            do i = 1,ngllx-2
                idxSpml = renumSpml(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                if (idxSpml .EQ. -1) stop "Unexpected non pml outer Edge !"
                Tdomain%OuterPMLNodes(abscount) = idxSpml
                abscount = abscount + 1
            enddo
        endif

        allocate(Tdomain%sEdge(n)%Renum(1:ngllx-2))
        Tdomain%sEdge(n)%Renum(:) = -1
        if (Tdomain%sEdge(n)%solid) then
            if (Tdomain%sEdge(n)%PML) then
                do i = 1,ngllx-2
                    Tdomain%sEdge(n)%Renum(i) = renumSpml(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                enddo
            else
                do i = 1,ngllx-2
                    Tdomain%sEdge(n)%Renum(i) = renumS(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                enddo
            endif
        else
            if (Tdomain%sEdge(n)%PML) then
                do i = 1,ngllx-2
                    Tdomain%sEdge(n)%Renum(i) = renumFpml(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                enddo
            else
                do i = 1,ngllx-2
                    Tdomain%sEdge(n)%Renum(i) = renumF(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                enddo
            endif
        endif
    enddo

    do n = 0,Tdomain%n_vertex-1
        if (Tdomain%sVertex(n)%Abs .and. Tdomain%sVertex(n)%PML) then
            idxSpml = renumSpml(Tdomain%sVertex(n)%Iglobnum_Vertex)
            if (idxSpml .EQ. -1) stop "Unexpected non pml outer Vertex !"
            Tdomain%OuterPMLNodes(abscount) = idxSpml
            abscount = abscount + 1
        endif

        if (Tdomain%sVertex(n)%solid) then
            if (Tdomain%sVertex(n)%PML) then
                Tdomain%sVertex(n)%Renum = renumSpml(Tdomain%sVertex(n)%Iglobnum_Vertex)
            else
                Tdomain%sVertex(n)%Renum = renumS(Tdomain%sVertex(n)%Iglobnum_Vertex)
            endif
        else
            if (Tdomain%sVertex(n)%PML) then
                Tdomain%sVertex(n)%Renum = renumFpml(Tdomain%sVertex(n)%Iglobnum_Vertex)
            else
                Tdomain%sVertex(n)%Renum = renumF(Tdomain%sVertex(n)%Iglobnum_Vertex)
            endif
        endif
    enddo


    call prepare_comm_vector(Tdomain,rank,Tdomain%Comm_mass, 1, 4, 1, 1)
    call prepare_comm_vector(Tdomain,rank,Tdomain%Comm_forces, 3, 9, 1, 3)
!     call debug_comm_vector(Tdomain, rank, 0, 1, Tdomain%Comm_mass)
!     call debug_comm_vector(Tdomain, rank, 0, 2, Tdomain%Comm_mass)
!     call debug_comm_vector(Tdomain, rank, 0, 3, Tdomain%Comm_mass)
!     call debug_comm_vector(Tdomain, rank, 3, 0, Tdomain%Comm_mass)

    deallocate(renumS)
    deallocate(renumSpml)
    deallocate(renumF)
    deallocate(renumFpml)
#endif

    return
end subroutine global_numbering

#if NEW_GLOBAL_METHOD

subroutine prepare_comm_vector(Tdomain,rank,comm_data, nddlsol, nddlsolpml, nddlfluid, nddlfluidpml)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain
    type(comm_vector), intent(inout) :: comm_data
    integer, intent(in) :: nddlsol, nddlsolpml, nddlfluid, nddlfluidpml
    integer, intent(in) :: rank

    integer :: n,nproc,nsol,nsolpml,nflu,nflupml
    integer :: i,j,k,nf,ne,nv,idx,ngll1,ngll2

    ! Remplissage des IGive et ITake
    if(Tdomain%n_proc > 1)then
        call allocate_comm_vector(Tdomain,rank, comm_data, nddlsol, nddlsolpml, &
                                  nddlfluid, nddlfluidpml)

        do n = 0,Comm_data%ncomm-1
            nproc = Comm_data%Data(n)%dest
            nsol = 0
            nsolpml = 0
            nflu = 0
            nflupml = 0

            ! Remplissage des Igive
            ! Faces
            do i = 0,Tdomain%sComm(nproc)%nb_faces-1
                nf = Tdomain%sComm(nproc)%faces(i)
                if (Tdomain%sFace(nf)%solid) then
                    if (Tdomain%sFace(nf)%PML) then
                        do j = 1,Tdomain%sFace(nf)%ngll2-2
                            do k = 1,Tdomain%sFace(nf)%ngll1-2
                                idx = Tdomain%sFace(nf)%Renum(k,j)
                                Comm_data%Data(n)%IGiveSPML(nsolpml) = idx
                                nsolpml = nsolpml + 1
                            end do
                        end do
                    else
                        do j = 1,Tdomain%sFace(nf)%ngll2-2
                            do k = 1,Tdomain%sFace(nf)%ngll1-2
                                idx = Tdomain%sFace(nf)%Renum(k,j)
                                Comm_data%Data(n)%IGiveS(nsol) = idx
                                nsol = nsol + 1
                            end do
                        end do
                    endif
                else ! Fluid
                    if (Tdomain%sFace(nf)%PML) then
                        do j = 1,Tdomain%sFace(nf)%ngll2-2
                            do k = 1,Tdomain%sFace(nf)%ngll1-2
                                idx = Tdomain%sFace(nf)%Renum(k,j)
                                Comm_data%Data(n)%IGiveFPML(nflupml) = idx
                                nflupml = nflupml + 1
                            end do
                        end do
                    else
                        do j = 1,Tdomain%sFace(nf)%ngll2-2
                            do k = 1,Tdomain%sFace(nf)%ngll1-2
                                idx = Tdomain%sFace(nf)%Renum(k,j)
                                Comm_data%Data(n)%IGiveF(nflu) = idx
                                nflu = nflu + 1
                            end do
                        end do
                    endif
                endif
            enddo

            ! Edges
            do i = 0,Tdomain%sComm(nproc)%nb_edges-1
                ne = Tdomain%sComm(nproc)%edges(i)
                if (Tdomain%sEdge(ne)%solid) then
                    if (Tdomain%sEdge(ne)%PML) then
                        do j = 1,Tdomain%sEdge(ne)%ngll-2
                            idx = Tdomain%sEdge(ne)%Renum(j)
                            Comm_data%Data(n)%IGiveSPML(nsolpml) = idx
                            nsolpml = nsolpml + 1
                        enddo
                    else
                        do j = 1,Tdomain%sEdge(ne)%ngll-2
                            idx = Tdomain%sEdge(ne)%Renum(j)
                            Comm_data%Data(n)%IGiveS(nsol) = idx
                            nsol = nsol + 1
                        enddo
                    endif
                else ! Fluid
                    if (Tdomain%sEdge(ne)%PML) then
                        do j = 1,Tdomain%sEdge(ne)%ngll-2
                            idx = Tdomain%sEdge(ne)%Renum(j)
                            Comm_data%Data(n)%IGiveFPML(nflupml) = idx
                            nflupml = nflupml + 1
                        enddo
                    else
                        do j = 1,Tdomain%sEdge(ne)%ngll-2
                            idx = Tdomain%sEdge(ne)%Renum(j)
                            Comm_data%Data(n)%IGiveF(nflu) = idx
                            nflu = nflu + 1
                        enddo
                    endif
                endif
            enddo

            ! Vertices
            do i = 0,Tdomain%sComm(nproc)%nb_vertices-1
                nv =  Tdomain%sComm(nproc)%vertices(i)
                idx = Tdomain%sVertex(nv)%Renum
                if (Tdomain%svertex(nv)%solid) then
                    if (Tdomain%svertex(nv)%PML) then
                        Comm_data%Data(n)%IGiveSPML(nsolpml) = idx
                        nsolpml = nsolpml + 1
                    else
                        Comm_data%Data(n)%IGiveS(nsol) = idx
                        nsol = nsol + 1
                    endif
                else ! Fluid
                    if (Tdomain%svertex(nv)%PML) then
                        Comm_data%Data(n)%IGiveFPML(nflupml) = idx
                        nflupml = nflupml + 1
                    else
                        Comm_data%Data(n)%IGiveF(nflu) = idx
                        nflu = nflu + 1
                    endif
                endif
            enddo


            nsol = 0
            nsolpml = 0
            nflu = 0
            nflupml = 0
            ! Remplissage des ITake
            ! Faces
            do i = 0,Tdomain%sComm(nproc)%nb_faces-1
                nf = Tdomain%sComm(nproc)%faces(i)
                ngll1 = Tdomain%sFace(nf)%ngll1
                ngll2 = Tdomain%sFace(nf)%ngll2

                if ( Tdomain%sComm(nproc)%orient_faces(i) == 0 ) then
                    if (Tdomain%sFace(nf)%solid) then ! solide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(k,j)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(k,j)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                enddo
                            enddo
                        endif
                    else ! fluide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(k,j)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(k,j)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                enddo
                            enddo
                        endif
                    endif

                else if ( Tdomain%sComm(nproc)%orient_faces(i) == 1 ) then
                    if (Tdomain%sFace(nf)%solid) then ! solide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,j)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,j)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                enddo
                            enddo
                        endif
                    else ! fluide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,j)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,j)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                enddo
                            enddo
                        endif
                    endif

                else if ( Tdomain%sComm(nproc)%orient_faces(i) == 2 ) then
                    if (Tdomain%sFace(nf)%solid) then ! solide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                enddo
                            enddo
                        endif
                    else ! fluide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                enddo
                            enddo
                        endif
                    endif

                else if ( Tdomain%sComm(nproc)%orient_faces(i) == 3 ) then
                    if (Tdomain%sFace(nf)%solid) then ! solide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                enddo
                            enddo
                        endif
                    else ! fluide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll2-2
                                do k = 1,ngll1-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                enddo
                            enddo
                        endif
                    endif

                else if ( Tdomain%sComm(nproc)%orient_faces(i) == 4 ) then
                    if (Tdomain%sFace(nf)%solid) then ! solide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(j,k)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(j,k)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                enddo
                            enddo
                        endif
                    else ! fluide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(j,k)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(j,k)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                enddo
                            enddo
                        endif
                    endif

                else if ( Tdomain%sComm(nproc)%orient_faces(i) == 5 ) then
                    if (Tdomain%sFace(nf)%solid) then ! solide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,k)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,k)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                enddo
                            enddo
                        endif
                    else ! fluide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,k)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,k)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                enddo
                            enddo
                        endif
                    endif

                else if ( Tdomain%sComm(nproc)%orient_faces(i) == 6 ) then
                    if (Tdomain%sFace(nf)%solid) then ! solide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                enddo
                            enddo
                        endif
                    else ! fluide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                enddo
                            enddo
                        endif
                    endif

                else if ( Tdomain%sComm(nproc)%orient_faces(i) == 7 ) then
                    if (Tdomain%sFace(nf)%solid) then ! solide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                enddo
                            enddo
                        endif
                    else ! fluide
                        if (Tdomain%sFace(nf)%PML) then
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                enddo
                            enddo
                        else
                            do j = 1,ngll1-2
                                do k = 1,ngll2-2
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                enddo
                            enddo
                        endif
                    endif
                endif ! orient_f
            enddo            

            ! Edges
            do i = 0,Tdomain%sComm(nproc)%nb_edges-1
                ne = Tdomain%sComm(nproc)%edges(i)
                ngll1 = Tdomain%sEdge(ne)%ngll

                if ( Tdomain%sComm(nproc)%orient_edges(i) == 0 ) then
                    if (Tdomain%sEdge(ne)%solid) then
                        if (Tdomain%sEdge(ne)%PML) then
                            do j = 1,ngll1-2
                                idx = Tdomain%sEdge(ne)%Renum(j)
                                Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                nsolpml = nsolpml + 1
                            enddo
                        else
                            do j = 1,ngll1-2
                                idx = Tdomain%sEdge(ne)%Renum(j)
                                Comm_data%Data(n)%ITakeS(nsol) = idx
                                nsol = nsol + 1
                            enddo
                        endif
                    else
                        if (Tdomain%sEdge(ne)%PML) then
                            do j = 1,ngll1-2
                                idx = Tdomain%sEdge(ne)%Renum(j)
                                Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                nflupml = nflupml + 1
                            enddo
                        else
                            do j = 1,ngll1-2
                                idx = Tdomain%sEdge(ne)%Renum(j)
                                Comm_data%Data(n)%ITakeF(nflu) = idx
                                nflu = nflu + 1
                            enddo
                        endif
                    endif

                else if ( Tdomain%sComm(nproc)%orient_edges(i) == 1 ) then
                    if (Tdomain%sEdge(ne)%solid) then
                        if (Tdomain%sEdge(ne)%PML) then
                            do j = 1,ngll1-2
                                idx = Tdomain%sEdge(ne)%Renum(ngll1-1-j)
                                Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                nsolpml = nsolpml + 1
                            enddo
                        else
                            do j = 1,ngll1-2
                                idx = Tdomain%sEdge(ne)%Renum(ngll1-1-j)
                                Comm_data%Data(n)%ITakeS(nsol) = idx
                                nsol = nsol + 1
                            enddo
                        endif
                    else
                        if (Tdomain%sEdge(ne)%PML) then
                            do j = 1,ngll1-2
                                idx = Tdomain%sEdge(ne)%Renum(ngll1-1-j)
                                Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                nflupml = nflupml + 1
                            enddo
                        else
                            do j = 1,ngll1-2
                                idx = Tdomain%sEdge(ne)%Renum(ngll1-1-j)
                                Comm_data%Data(n)%ITakeF(nflu) = idx
                                nflu = nflu + 1
                            enddo
                        endif
                    endif
                endif ! orient_edges
            enddo

            ! Vertices
            do i = 0,Tdomain%sComm(nproc)%nb_vertices-1
                nv =  Tdomain%sComm(nproc)%vertices(i)
                idx = Tdomain%sVertex(nv)%Renum
                if (Tdomain%sVertex(nv)%solid) then
                    if (Tdomain%sVertex(nv)%PML) then
                        Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                        nsolpml = nsolpml + 1
                    else
                        Comm_data%Data(n)%ITakeS(nsol) = idx
                        nsol = nsol + 1
                    endif
                else
                    if (Tdomain%sVertex(nv)%PML) then
                        Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                        nflupml = nflupml + 1
                    else
                        Comm_data%Data(n)%ITakeF(nflu) = idx
                        nflu = nflu + 1
                    endif
                endif
            enddo
        enddo 
    endif

end subroutine prepare_comm_vector


subroutine allocate_comm_vector(Tdomain,rank, comm_data, nddlsol, nddlsolpml, nddlfluid, nddlfluidpml)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain
    type(comm_vector), intent(inout) :: comm_data
    integer, intent(in) :: nddlsol, nddlsolpml, nddlfluid, nddlfluidpml
    integer, intent(in) :: rank
    integer :: n_data, n_comm, nsol, nsolpml, nflu, nflupml
    integer :: n, nf, ne, nv, i, temp

    n_comm = 0
    do n = 0,Tdomain%n_proc-1
        if (Tdomain%sComm(n)%nb_faces > 0 .OR. &
            Tdomain%sComm(n)%nb_edges > 0 .OR. &
            Tdomain%sComm(n)%nb_vertices > 0) then
            n_comm = n_comm + 1
        endif
    enddo

    allocate(Comm_data%Data(0:n_comm-1))
    allocate(Comm_data%send_reqs(0:n_comm-1))
    allocate(Comm_data%recv_reqs(0:n_comm-1))
    Comm_data%ncomm = n_comm

    n_comm = 0
    do n = 0,Tdomain%n_proc-1
        if (Tdomain%sComm(n)%nb_faces < 1 .AND. &
            Tdomain%sComm(n)%nb_edges < 1 .AND. &
            Tdomain%sComm(n)%nb_vertices < 1) cycle

        nsolpml = 0
        nsol = 0
        nflupml = 0
        nflu = 0

        ! Faces
        do i = 0,Tdomain%sComm(n)%nb_faces-1
            nf = Tdomain%sComm(n)%faces(i)
            temp = (Tdomain%sFace(nf)%ngll1-2) * (Tdomain%sFace(nf)%ngll2-2)
            if (Tdomain%sFace(nf)%solid) then
                if (Tdomain%sFace(nf)%PML) then
                    nsolpml = nsolpml + temp
                else
                    nsol = nsol + temp
                endif
            else
                if (Tdomain%sFace(nf)%PML) then
                    nflupml = nflupml + temp
                else
                    nflu = nflu + temp
                endif
            endif
        enddo
        ! Edges
        do i = 0,Tdomain%sComm(n)%nb_edges-1
            ne = Tdomain%sComm(n)%edges(i)
            temp = Tdomain%sEdge(ne)%ngll-2
            if (Tdomain%sEdge(ne)%solid) then
                if (Tdomain%sEdge(ne)%PML) then
                    nsolpml = nsolpml + temp
                else
                    nsol = nsol + temp
                endif
            else
                if (Tdomain%sEdge(ne)%PML) then
                    nflupml = nflupml + temp
                else
                    nflu = nflu + temp
                endif
            endif
        enddo
        ! Vertices
        do i = 0,Tdomain%sComm(n)%nb_vertices-1
            nv = Tdomain%sComm(n)%vertices(i)
            if (Tdomain%svertex(nv)%solid) then
                if (Tdomain%svertex(nv)%PML) then
                    nsolpml = nsolpml + 1
                else
                    nsol = nsol + 1
                endif
            else
                if (Tdomain%svertex(nv)%PML) then
                    nflupml = nflupml + 1
                else
                    nflu = nflu + 1
                endif
            endif
        enddo

        n_data = nddlsol*nsol+nddlsolpml*nsolpml+nddlfluid*nflu+nddlfluidpml*nflupml
        ! Initialisation et allocation de Comm_vector_DumpMassAndMMSP
        Comm_data%Data(n_comm)%src = rank
        Comm_data%Data(n_comm)%dest = n
        Comm_data%Data(n_comm)%ndata = n_data
        Comm_data%Data(n_comm)%nsol = nsol
        Comm_data%Data(n_comm)%nsolpml = nsolpml
        Comm_data%Data(n_comm)%nflu = nflu
        Comm_data%Data(n_comm)%nflupml = nflupml
        allocate(Comm_data%Data(n_comm)%Give(0:n_data-1))
        allocate(Comm_data%Data(n_comm)%Take(0:n_data-1))
        allocate(Comm_data%Data(n_comm)%IGiveS(0:nsol-1))
        allocate(Comm_data%Data(n_comm)%ITakeS(0:nsol-1))
        allocate(Comm_data%Data(n_comm)%IGiveSPML(0:nsolpml-1))
        allocate(Comm_data%Data(n_comm)%ITakeSPML(0:nsolpml-1))
        allocate(Comm_data%Data(n_comm)%IGiveF(0:nflu-1))
        allocate(Comm_data%Data(n_comm)%ITakeF(0:nflu-1))
        allocate(Comm_data%Data(n_comm)%IGiveFPML(0:nflupml-1))
        allocate(Comm_data%Data(n_comm)%ITakeFPML(0:nflupml-1))

        n_comm = n_comm + 1
    enddo

    return
end subroutine allocate_comm_vector

subroutine populate_index_SF(ngll,dir,ngllx,nglly,ngllz,ISolFlu,k,SF_IGlob)
    implicit none

    integer, intent(in) :: ngll, dir, ngllx, nglly, ngllz
    integer, intent(in), dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: ISolFlu
    integer, intent(inout) :: k
    integer, intent(inout), dimension(0:ngll-1) :: SF_IGlob
    integer :: i,j

    if(dir == 0)then
        do j = 0,nglly-1
            do i = 0,ngllx-1
                SF_IGlob(k) = ISolFlu(i,j,0)
                k = k + 1
            enddo
        enddo
    else if(dir == 1)then
        do j = 0,ngllz-1
            do i = 0,ngllx-1
                SF_IGlob(k) = ISolFlu(i,0,j)
                k = k + 1
            enddo
        enddo
    else if(dir == 2)then
        do j = 0,ngllz-1
            do i = 0,nglly-1
                SF_IGlob(k) = ISolFlu(ngllx-1,i,j)
                k = k + 1
            enddo
        enddo
    else if(dir == 3)then
        do j = 0,ngllz-1
            do i = 0,ngllx-1
                SF_IGlob(k) = ISolFlu(i,nglly-1,j)
                k = k + 1
            enddo
        enddo
    else if(dir == 4)then
        do j = 0,ngllz-1
            do i = 0,nglly-1
                SF_IGlob(k) = ISolFlu(0,i,j)
                k = k + 1
            enddo
        enddo
    else if(dir == 5)then
        do j = 0,nglly-1
            do i = 0,ngllx-1
                SF_IGlob(k) = ISolFlu(i,j,ngllz-1)
                k = k + 1
            enddo
        enddo
    end if

    return
end subroutine populate_index_SF


subroutine renum_element(ngllx, nglly, ngllz, Iglobnum, Irenum, idx, renum, isPML)
    implicit none

    integer, intent(in) :: ngllx, nglly, ngllz
    integer, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Iglobnum
    integer, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: Irenum
    integer, intent (inout) :: idx
    integer, dimension (:), allocatable, intent (inout) :: renum
    logical, intent(in) :: isPML

    integer :: i, j, k, num, ind, ddl

    if (isPML) then
        ddl = 3
    else
        ddl = 1
    endif

    do k = 0,ngllz-1
        do j = 0,nglly-1
            do i = 0,ngllx-1
                ind = Iglobnum(i,j,k)
                num = renum(ind)
                if (num == -1) then
                    renum(ind) = idx
                    num = idx
                    idx = idx + ddl
                endif
                Irenum(i,j,k) = num
            enddo
        enddo
    enddo

    return
end subroutine renum_element

subroutine renum_face(ngllx, nglly, Iglobnum, idx, renum, isPML)
    implicit none

    integer, intent(in) :: ngllx, nglly
    integer, dimension(1:ngllx-2,1:nglly-2), intent(in) :: Iglobnum
    integer, intent (inout) :: idx
    integer, dimension (:), allocatable, intent (inout) :: renum
    logical, intent(in) :: isPML

    integer :: i, j, ind, ddl

    if (isPML) then
        ddl = 3
    else
        ddl = 1
    endif

    do j = 1,nglly-2
        do i = 1,ngllx-2
            ind = Iglobnum(i,j)
            if (renum(ind) == -1) then
                renum(ind) = idx
                idx = idx + ddl
            endif
        enddo
    enddo

    return
end subroutine renum_face

subroutine renum_edge(ngllx, Iglobnum, idx, renum, isPML)
    implicit none

    integer, intent(in) :: ngllx
    integer, dimension(1:ngllx-2), intent(in) :: Iglobnum
    integer, intent (inout) :: idx
    integer, dimension (:), allocatable, intent (inout) :: renum
    logical, intent(in) :: isPML

    integer :: i, ind, ddl

    if (isPML) then
        ddl = 3
    else
        ddl = 1
    endif

    do i = 1,ngllx-2
        ind = Iglobnum(i)
        if (renum(ind) == -1) then
            renum(ind) = idx
            idx = idx + ddl
        endif
    enddo

    return
end subroutine renum_edge

subroutine renum_vertex(Iglobnum, idx, renum, isPML)
    implicit none

    integer, intent(in) :: Iglobnum
    integer, intent (inout) :: idx
    integer, dimension (:), allocatable, intent (inout) :: renum
    logical, intent(in) :: isPML

    integer :: ddl

    if (isPML) then
        ddl = 3
    else
        ddl = 1
    endif

    if (renum(Iglobnum) == -1) then
        renum(Iglobnum) = idx
        idx = idx + ddl
    endif

    return
end subroutine renum_vertex

subroutine debug_comm_vector(Tdomain, rank, src, dest, comm)
    use sdomain
    implicit none
    type(domain), intent(in) :: Tdomain
    integer, intent(in) :: rank, src, dest
    type(comm_vector), intent(in) :: comm
    !
    integer :: i,k

    do i=0, comm%ncomm-1
        if (comm%Data(i)%src/=src .or. comm%Data(i)%dest/=dest) cycle

        write(*,*) rank, "COMM:", src, dest
        write(*,*) rank, "NSOL ", comm%Data(i)%nsol
        write(*,*) rank, "NSPML", comm%Data(i)%nsolpml
        write(*,*) rank, "NFLU ", comm%Data(i)%nflu
        write(*,*) rank, "NFPML", comm%Data(i)%nflupml

        write(*,*) rank, "ISOL>", (comm%Data(i)%IGiveS(k), k=0,10)
        write(*,*) rank, "ISOL<", (comm%Data(i)%ITakeS(k), k=0,10)
    end do
end subroutine debug_comm_vector


#endif
end module mrenumber
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
