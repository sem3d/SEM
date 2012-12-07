!>
!!\file global_numbering.f90
!!\brief Assure la correspondance entre les différentes numérotations.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<


!>
!! Définition de Iglobnum et  renvoi du nombre total de ddl: elements, faces, aretes, sommets
!!
!<

subroutine global_numbering (Tdomain,rg)

    ! Modified by Paul Cupillard 08/12/2005
    ! Commentaire GS:
    ! On parcourt tous les elements, on parcourt pour k fixe,j fixe, i fixe, en excluant les valeurs 0 et n-1, on definit Iglobnum en incrementant de 1 a chaque fois.
    ! Pour Iglobnum_face, on demarre de la derniere valeur de Iglobnum, on parcourt toutes les faces, on parcourt pour j fixe, i fixe, en excluant les valeurs 0 et n-1, on definit Iglobnum_face en incrementant de 1 a chaque fois.
    ! Pour Iglobnum_edge, on demarre de la derniere valeur de Iglobnum_face, on parcourt toutes les aretes, on parcourt pour i fixe, en excluant les valeurs 0 et n-1, on definit Iglobnum_edge en incrementant de 1 a chaque fois.
    ! Pour Iglobnum_vertex, on demarre de la derniere valeur de Iglobnum_edge, on parcourt tous les vertices, on definit Iglobnum_vertex en incrementant de 1 a chaque fois.
    ! on obtient finalement le nb total de ddl: elements, faces, aretes, sommets
    !
    use sdomain

    implicit none

    type(domain), target, intent (INOUT) :: Tdomain
    integer, intent(IN) :: rg

    integer :: n, icount, i,j,k, ngllx,nglly,ngllz, nf,ne,nv, ngll1,ngll2,nf_aus,ne_aus
    logical, dimension(:), allocatable :: L_Face, L_Edge
    integer, dimension(:), allocatable :: FaceNum, EdgeNum


    icount = 0

    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        allocate (Tdomain%specel(n)%Iglobnum(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
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
        allocate (Tdomain%sFace(n)%Iglobnum_Face(1:ngllx-2,1:nglly-2))
        Tdomain%sFace(n)%Iglobnum_Face = -1
        do j = 1,nglly-2
            do i = 1,ngllx-2
                Tdomain%sFace(n)%Iglobnum_Face(i,j) = icount
                icount = icount + 1 !icount non remis a zero
            enddo
        enddo
    enddo

    do n = 0,Tdomain%n_edge-1
        ngllx = Tdomain%sEdge(n)%ngll
        allocate (Tdomain%sEdge(n)%Iglobnum_Edge(1:ngllx-2))
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

    Tdomain%n_glob_points = icount
    print*,'proc Nb GLL points',rg,Tdomain%n_glob_points
    do n = 0,Tdomain%n_elem - 1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        ! Taking information from faces
        nf = Tdomain%specel(n)%near_faces(0)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        if ( Tdomain%specel(n)%Orient_Faces(0) == 0 ) then
            Tdomain%specel(n)%Iglobnum (1:ngll1-2, 1:ngll2-2, 0) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
        else if ( Tdomain%specel(n)%Orient_Faces(0) == 1 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll1-1-i, j, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(0) == 2 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (i, ngll2-1-j, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(0) == 3 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll1-1-i, ngll2-1-j, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(0) == 4 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (j, i, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(0) == 5 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll2-1-j, i, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(0) == 6 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (j, ngll1-1-i, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(0) == 7 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll2-1-j, ngll1-1-i, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        endif

        nf = Tdomain%specel(n)%near_faces(1)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        if ( Tdomain%specel(n)%Orient_Faces(1) == 0 ) then
            Tdomain%specel(n)%Iglobnum (1:ngll1-2, 0, 1:ngll2-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
        else if ( Tdomain%specel(n)%Orient_Faces(1) == 1 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll1-1-i, 0, j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(1) == 2 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (i, 0, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(1) == 3 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll1-1-i, 0, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(1) == 4 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (j, 0, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(1) == 5 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll2-1-j, 0, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(1) == 6 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (j, 0, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(1) == 7 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll2-1-j, 0, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        endif

        nf = Tdomain%specel(n)%near_faces(2)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        if ( Tdomain%specel(n)%Orient_Faces(2) == 0 ) then
            Tdomain%specel(n)%Iglobnum (ngllx-1, 1:ngll1-2, 1:ngll2-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
        else if ( Tdomain%specel(n)%Orient_Faces(2) == 1 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngllx-1, ngll1-1-i, j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(2) == 2 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngllx-1, i, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(2) == 3 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngllx-1, ngll1-1-i, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(2) == 4 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngllx-1, j, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(2) == 5 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngllx-1, ngll2-1-j, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(2) == 6 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngllx-1, j, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(2) == 7 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngllx-1, ngll2-1-j, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        endif

        nf = Tdomain%specel(n)%near_faces(3)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        if ( Tdomain%specel(n)%Orient_Faces(3) == 0 ) then
            Tdomain%specel(n)%Iglobnum (1:ngll1-2, nglly-1, 1:ngll2-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
        else if ( Tdomain%specel(n)%Orient_Faces(3) == 1 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll1-1-i, nglly-1, j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(3) == 2 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (i, nglly-1, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(3) == 3 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll1-1-i, nglly-1, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(3) == 4 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (j, nglly-1, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(3) == 5 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll2-1-j,nglly-1,i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(3) == 6 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (j, nglly-1, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(3) == 7 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll2-1-j, nglly-1, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        endif

        nf = Tdomain%specel(n)%near_faces(4)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        if ( Tdomain%specel(n)%Orient_Faces(4) == 0 ) then
            Tdomain%specel(n)%Iglobnum (0, 1:ngll1-2, 1:ngll2-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
        else if ( Tdomain%specel(n)%Orient_Faces(4) == 1 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (0, ngll1-1-i, j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(4) == 2 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (0, i, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(4) == 3 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (0, ngll1-1-i, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(4) == 4 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (0, j, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(4) == 5 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (0, ngll2-1-j, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(4) == 6 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (0, j, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(4) == 7 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (0, ngll2-1-j, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        endif

        nf = Tdomain%specel(n)%near_faces(5)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        if ( Tdomain%specel(n)%Orient_Faces(5) == 0 ) then
            Tdomain%specel(n)%Iglobnum (1:ngll1-2, 1:ngll2-2, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
        else if ( Tdomain%specel(n)%Orient_Faces(5) == 1 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll1-1-i, j, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(5) == 2 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (i, ngll2-1-j, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(5) == 3 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll1-1-i, ngll2-1-j, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(5) == 4 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (j, i, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(5) == 5 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll2-1-j, i, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(5) == 6 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (j, ngll1-1-i, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        else if ( Tdomain%specel(n)%Orient_Faces(5) == 7 ) then
            do i=1,ngll1-2
                do j=1,ngll2-2
                    Tdomain%specel(n)%Iglobnum (ngll2-1-j, ngll1-1-i, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
                enddo
            enddo
        endif


        ! Taking information from edges
        ne = Tdomain%specel(n)%Near_Edges(0)
        if ( Tdomain%specel(n)%Orient_Edges(0) == 0 ) then
            Tdomain%specel(n)%Iglobnum(1:ngllx-2,0,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
        else
            do i=1,ngllx-2
                Tdomain%specel(n)%Iglobnum(i,0,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllx-1-i)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(1)
        if ( Tdomain%specel(n)%Orient_Edges(1) == 0 ) then
            Tdomain%specel(n)%Iglobnum(ngllx-1,1:nglly-2,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
        else
            do j=1,nglly-2
                Tdomain%specel(n)%Iglobnum(ngllx-1,j,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(nglly-1-j)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(2)
        if ( Tdomain%specel(n)%Orient_Edges(2) == 0 ) then
            Tdomain%specel(n)%Iglobnum(1:ngllx-2,nglly-1,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
        else
            do i=1,ngllx-2
                Tdomain%specel(n)%Iglobnum(i,nglly-1,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllx-1-i)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(3)
        if ( Tdomain%specel(n)%Orient_Edges(3) == 0 ) then
            Tdomain%specel(n)%Iglobnum(0,1:nglly-2,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
        else
            do j=1,nglly-2
                Tdomain%specel(n)%Iglobnum(0,j,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(nglly-1-j)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(4)
        if ( Tdomain%specel(n)%Orient_Edges(4) == 0 ) then
            Tdomain%specel(n)%Iglobnum(ngllx-1,0,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
        else
            do k=1,ngllz-2
                Tdomain%specel(n)%Iglobnum(ngllx-1,0,k) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllz-1-k)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(5)
        if ( Tdomain%specel(n)%Orient_Edges(5) == 0 ) then
            Tdomain%specel(n)%Iglobnum(1:ngllx-2,0,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
        else
            do i=1,ngllx-2
                Tdomain%specel(n)%Iglobnum(i,0,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllx-1-i)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(6)
        if ( Tdomain%specel(n)%Orient_Edges(6) == 0 ) then
            Tdomain%specel(n)%Iglobnum(0,0,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
        else
            do k=1,ngllz-2
                Tdomain%specel(n)%Iglobnum(0,0,k) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllz-1-k)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(7)
        if ( Tdomain%specel(n)%Orient_Edges(7) == 0 ) then
            Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
        else
            do k=1,ngllz-2
                Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,k) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllz-1-k)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(8)
        if ( Tdomain%specel(n)%Orient_Edges(8) == 0 ) then
            Tdomain%specel(n)%Iglobnum(ngllx-1,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
        else
            do j=1,nglly-2
                Tdomain%specel(n)%Iglobnum(ngllx-1,j,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(nglly-1-j)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(9)
        if ( Tdomain%specel(n)%Orient_Edges(9) == 0 ) then
            Tdomain%specel(n)%Iglobnum(1:ngllx-2,nglly-1,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
        else
            do i=1,ngllx-2
                Tdomain%specel(n)%Iglobnum(i,nglly-1,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllx-1-i)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(10)
        if ( Tdomain%specel(n)%Orient_Edges(10) == 0 ) then
            Tdomain%specel(n)%Iglobnum(0,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
        else
            do k=1,ngllz-2
                Tdomain%specel(n)%Iglobnum(0,nglly-1,k) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllz-1-k)
            enddo
        endif
        ne = Tdomain%specel(n)%Near_Edges(11)
        if ( Tdomain%specel(n)%Orient_Edges(11) == 0 ) then
            Tdomain%specel(n)%Iglobnum(0,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
        else
            do j=1,nglly-2
                Tdomain%specel(n)%Iglobnum(0,j,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(nglly-1-j)
            enddo
        endif


        ! Taking information from vertices
        nv = Tdomain%specel(n)%Near_Vertices(0)
        Tdomain%specel(n)%Iglobnum(0,0,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
        nv = Tdomain%specel(n)%Near_Vertices(1)
        Tdomain%specel(n)%Iglobnum(ngllx-1,0,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
        nv = Tdomain%specel(n)%Near_Vertices(2)
        Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
        nv = Tdomain%specel(n)%Near_Vertices(3)
        Tdomain%specel(n)%Iglobnum(0,nglly-1,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
        nv = Tdomain%specel(n)%Near_Vertices(4)
        Tdomain%specel(n)%Iglobnum(0,0,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
        nv = Tdomain%specel(n)%Near_Vertices(5)
        Tdomain%specel(n)%Iglobnum(ngllx-1,0,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
        nv = Tdomain%specel(n)%Near_Vertices(6)
        Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
        nv = Tdomain%specel(n)%Near_Vertices(7)
        Tdomain%specel(n)%Iglobnum(0,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    enddo



    allocate (L_Face(0:Tdomain%n_face-1))
    L_Face = .true.
    allocate (L_Edge(0:Tdomain%n_edge-1))
    L_Edge = .true.
    do n = 0,Tdomain%n_elem-1

        ! Faces
        do i = 0,5
            nf = Tdomain%specel(n)%Near_Faces(i)
            if (L_Face(nf)) then
                if ( Tdomain%specel(n)%Orient_Faces(i) == 0 ) then
                    L_Face(nf) = .false.
                    allocate(Tdomain%sFace(nf)%FaceNum(0:3))
                    if (i==0) then
                        Tdomain%sFace(nf)%FaceNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                        Tdomain%sFace(nf)%FaceNum(1) = Tdomain%specel(n)%Control_Nodes(1)
                        Tdomain%sFace(nf)%FaceNum(2) = Tdomain%specel(n)%Control_Nodes(2)
                        Tdomain%sFace(nf)%FaceNum(3) = Tdomain%specel(n)%Control_Nodes(3)
                    else if (i==1) then
                        Tdomain%sFace(nf)%FaceNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                        Tdomain%sFace(nf)%FaceNum(1) = Tdomain%specel(n)%Control_Nodes(1)
                        Tdomain%sFace(nf)%FaceNum(2) = Tdomain%specel(n)%Control_Nodes(5)
                        Tdomain%sFace(nf)%FaceNum(3) = Tdomain%specel(n)%Control_Nodes(4)
                    else if (i==2) then
                        Tdomain%sFace(nf)%FaceNum(0) = Tdomain%specel(n)%Control_Nodes(1)
                        Tdomain%sFace(nf)%FaceNum(1) = Tdomain%specel(n)%Control_Nodes(2)
                        Tdomain%sFace(nf)%FaceNum(2) = Tdomain%specel(n)%Control_Nodes(6)
                        Tdomain%sFace(nf)%FaceNum(3) = Tdomain%specel(n)%Control_Nodes(5)
                    else if (i==3) then
                        Tdomain%sFace(nf)%FaceNum(0) = Tdomain%specel(n)%Control_Nodes(3)
                        Tdomain%sFace(nf)%FaceNum(1) = Tdomain%specel(n)%Control_Nodes(2)
                        Tdomain%sFace(nf)%FaceNum(2) = Tdomain%specel(n)%Control_Nodes(6)
                        Tdomain%sFace(nf)%FaceNum(3) = Tdomain%specel(n)%Control_Nodes(7)
                    else if (i==4) then
                        Tdomain%sFace(nf)%FaceNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                        Tdomain%sFace(nf)%FaceNum(1) = Tdomain%specel(n)%Control_Nodes(3)
                        Tdomain%sFace(nf)%FaceNum(2) = Tdomain%specel(n)%Control_Nodes(7)
                        Tdomain%sFace(nf)%FaceNum(3) = Tdomain%specel(n)%Control_Nodes(4)
                    else if (i==5) then
                        Tdomain%sFace(nf)%FaceNum(0) = Tdomain%specel(n)%Control_Nodes(4)
                        Tdomain%sFace(nf)%FaceNum(1) = Tdomain%specel(n)%Control_Nodes(5)
                        Tdomain%sFace(nf)%FaceNum(2) = Tdomain%specel(n)%Control_Nodes(6)
                        Tdomain%sFace(nf)%FaceNum(3) = Tdomain%specel(n)%Control_Nodes(7)
                    endif
                endif
            endif
        enddo

        ! Edges
        do i = 0,11
            ne = Tdomain%specel(n)%Near_Edges(i)
            if (L_Edge(ne)) then
                if ( Tdomain%specel(n)%Orient_Edges(i) == 0 ) then
                    L_Edge(ne) = .false.
                    allocate (Tdomain%sEdge(ne)%EdgeNum(0:1))
                    if (i==0) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(1)
                    else if (i==1) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(1)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(2)
                    else if (i==2) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(3)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(2)
                    else if (i==3) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(3)
                    else if (i==4) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(1)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(5)
                    else if (i==5) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(4)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(5)
                    else if (i==6) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(4)
                    else if (i==7) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(2)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(6)
                    else if (i==8) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(5)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(6)
                    else if (i==9) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(7)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(6)
                    else if (i==10) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(3)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(7)
                    else if (i==11) then
                        Tdomain%sEdge(ne)%EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(4)
                        Tdomain%sEdge(ne)%EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(7)
                    endif
                endif
            endif
        enddo
    enddo
    deallocate(L_Face,L_Edge)


    ! Facenum ne sert pas a grandchose
    do n = 0,Tdomain%n_elem-1

        ! Faces
        do i = 0,5
            allocate(FaceNum(0:3))
            FaceNum = -1
            nf = Tdomain%specel(n)%Near_Faces(i)

            if (i==0) then
                FaceNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                FaceNum(1) = Tdomain%specel(n)%Control_Nodes(1)
                FaceNum(2) = Tdomain%specel(n)%Control_Nodes(2)
                FaceNum(3) = Tdomain%specel(n)%Control_Nodes(3)
            else if (i==1) then
                FaceNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                FaceNum(1) = Tdomain%specel(n)%Control_Nodes(1)
                FaceNum(2) = Tdomain%specel(n)%Control_Nodes(5)
                FaceNum(3) = Tdomain%specel(n)%Control_Nodes(4)
            else if (i==2) then
                FaceNum(0) = Tdomain%specel(n)%Control_Nodes(1)
                FaceNum(1) = Tdomain%specel(n)%Control_Nodes(2)
                FaceNum(2) = Tdomain%specel(n)%Control_Nodes(6)
                FaceNum(3) = Tdomain%specel(n)%Control_Nodes(5)
            else if (i==3) then
                FaceNum(0) = Tdomain%specel(n)%Control_Nodes(3)
                FaceNum(1) = Tdomain%specel(n)%Control_Nodes(2)
                FaceNum(2) = Tdomain%specel(n)%Control_Nodes(6)
                FaceNum(3) = Tdomain%specel(n)%Control_Nodes(7)
            else if (i==4) then
                FaceNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                FaceNum(1) = Tdomain%specel(n)%Control_Nodes(3)
                FaceNum(2) = Tdomain%specel(n)%Control_Nodes(7)
                FaceNum(3) = Tdomain%specel(n)%Control_Nodes(4)
            else if (i==5) then
                FaceNum(0) = Tdomain%specel(n)%Control_Nodes(4)
                FaceNum(1) = Tdomain%specel(n)%Control_Nodes(5)
                FaceNum(2) = Tdomain%specel(n)%Control_Nodes(6)
                FaceNum(3) = Tdomain%specel(n)%Control_Nodes(7)
            endif

            if ( Tdomain%specel(n)%Orient_Faces(i) == 0 ) then
                if (.not. (  (FaceNum(0) == Tdomain%sFace(nf)%FaceNum(0)) .and. &
                    (FaceNum(1) == Tdomain%sFace(nf)%FaceNum(1)) .and. &
                    (FaceNum(2) == Tdomain%sFace(nf)%FaceNum(2)) .and. &
                    (FaceNum(3) == Tdomain%sFace(nf)%FaceNum(3)) ) ) then
                    print*,'TFace0',rg,n,i,nf
                endif
            else if ( Tdomain%specel(n)%Orient_Faces(i) == 1 ) then
                if (.not. (  (FaceNum(0) == Tdomain%sFace(nf)%FaceNum(1)) .and. &
                    (FaceNum(1) == Tdomain%sFace(nf)%FaceNum(0)) .and. &
                    (FaceNum(2) == Tdomain%sFace(nf)%FaceNum(3)) .and. &
                    (FaceNum(3) == Tdomain%sFace(nf)%FaceNum(2)) ) ) then
                    print*,'TFace1',rg,n,i,nf
                endif
            else if ( Tdomain%specel(n)%Orient_Faces(i) == 2 ) then
                if (.not. (  (FaceNum(0) == Tdomain%sFace(nf)%FaceNum(3)) .and. &
                    (FaceNum(1) == Tdomain%sFace(nf)%FaceNum(2)) .and. &
                    (FaceNum(2) == Tdomain%sFace(nf)%FaceNum(1)) .and. &
                    (FaceNum(3) == Tdomain%sFace(nf)%FaceNum(0)) ) ) then
                    print*,'TFace2',rg,n,i,nf
                endif
            else if ( Tdomain%specel(n)%Orient_Faces(i) == 3 ) then
                if (.not. (  (FaceNum(0) == Tdomain%sFace(nf)%FaceNum(2)) .and. &
                    (FaceNum(1) == Tdomain%sFace(nf)%FaceNum(3)) .and. &
                    (FaceNum(2) == Tdomain%sFace(nf)%FaceNum(0)) .and. &
                    (FaceNum(3) == Tdomain%sFace(nf)%FaceNum(1)) ) ) then
                    print*,'TFace3',rg,n,i,nf
                endif
            else if ( Tdomain%specel(n)%Orient_Faces(i) == 4 ) then
                if (.not. (  (FaceNum(0) == Tdomain%sFace(nf)%FaceNum(0)) .and. &
                    (FaceNum(1) == Tdomain%sFace(nf)%FaceNum(3)) .and. &
                    (FaceNum(2) == Tdomain%sFace(nf)%FaceNum(2)) .and. &
                    (FaceNum(3) == Tdomain%sFace(nf)%FaceNum(1)) ) ) then
                    print*,'TFace4',rg,n,i,nf
                endif
            else if ( Tdomain%specel(n)%Orient_Faces(i) == 5 ) then
                if (.not. (  (FaceNum(0) == Tdomain%sFace(nf)%FaceNum(3)) .and. &
                    (FaceNum(1) == Tdomain%sFace(nf)%FaceNum(0)) .and. &
                    (FaceNum(2) == Tdomain%sFace(nf)%FaceNum(1)) .and. &
                    (FaceNum(3) == Tdomain%sFace(nf)%FaceNum(2)) ) ) then
                    print*,'TFace5',rg,n,i,nf
                    print*,'TFace5',FaceNum(0),FaceNum(1),FaceNum(2),FaceNum(3)
                    print*,'TFace5',Tdomain%sFace(nf)%FaceNum(0),Tdomain%sFace(nf)%FaceNum(1),Tdomain%sFace(nf)%FaceNum(2),Tdomain%sFace(nf)%FaceNum(3)
                    print*,'&&&&&&&&&&&&&&'
                endif
            else if ( Tdomain%specel(n)%Orient_Faces(i) == 6 ) then
                if (.not. (  (FaceNum(0) == Tdomain%sFace(nf)%FaceNum(1)) .and. &
                    (FaceNum(1) == Tdomain%sFace(nf)%FaceNum(2)) .and. &
                    (FaceNum(2) == Tdomain%sFace(nf)%FaceNum(3)) .and. &
                    (FaceNum(3) == Tdomain%sFace(nf)%FaceNum(0)) ) ) then
                    print*,'TFace6',rg,n,i,nf
                endif
            else if ( Tdomain%specel(n)%Orient_Faces(i) == 7 ) then
                if (.not. (  (FaceNum(0) == Tdomain%sFace(nf)%FaceNum(2)) .and. &
                    (FaceNum(1) == Tdomain%sFace(nf)%FaceNum(1)) .and. &
                    (FaceNum(2) == Tdomain%sFace(nf)%FaceNum(0)) .and. &
                    (FaceNum(3) == Tdomain%sFace(nf)%FaceNum(3)) ) ) then
                    print*,'TFace7',rg,n,i,nf
                endif
            endif
            deallocate(FaceNum)
        enddo



        !Edges
        allocate(EdgeNum(0:3))
        do i = 0,11
            ne = Tdomain%specel(n)%Near_Edges(i)

            if (i==0) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(1)
            else if (i==1) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(1)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(2)
            else if (i==2) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(3)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(2)
            else if (i==3) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(3)
            else if (i==4) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(1)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(5)
            else if (i==5) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(4)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(5)
            else if (i==6) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(0)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(4)
            else if (i==7) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(2)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(6)
            else if (i==8) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(5)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(6)
            else if (i==9) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(7)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(6)
            else if (i==10) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(3)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(7)
            else if (i==11) then
                EdgeNum(0) = Tdomain%specel(n)%Control_Nodes(4)
                EdgeNum(1) = Tdomain%specel(n)%Control_Nodes(7)
            endif

            if ( Tdomain%specel(n)%Orient_Edges(i) == 0 ) then
                if (.not. (  (EdgeNum(0) == Tdomain%sEdge(ne)%EdgeNum(0)) .and. &
                    (EdgeNum(1) == Tdomain%sEdge(ne)%EdgeNum(1)) ) ) then
                    print*,'TEdge0',n,i,ne
                endif
            else
                if (.not. (  (EdgeNum(0) == Tdomain%sEdge(ne)%EdgeNum(1)) .and. &
                    (EdgeNum(1) == Tdomain%sEdge(ne)%EdgeNum(0)) ) ) then
                    print*,'TEdge0',n,i,ne
                endif
            endif
        enddo
        deallocate(EdgeNum)

    enddo
    !!

    !purge fuites memoire
    do nf=0,Tdomain%n_face-1
        deallocate (Tdomain%sFace(nf)%FaceNum)
    enddo

    do ne=0,Tdomain%n_edge-1
        deallocate (Tdomain%sEdge(ne)%EdgeNum)
    enddo



    !! Test pour verifier que l'orientation des edges par rapports aux faces est la meme dans le SO PW que dans le milieu UP et aussi pour Neu!!
    !
    if ( Tdomain%logicD%super_object_local_present ) then
        !print*,'Nb Edges PW',Tdomain%sPlaneW%n_edges
        do nf = 0, Tdomain%sPlaneW%n_faces-1

            nf_aus = Tdomain%sPlaneW%pFace(nf)%Face_UP
            do i = 0,3
                ne = Tdomain%sPlaneW%pFace(nf)%Near_Edges(i)
                ne_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_UP

                do icount = 0,10
                    if ( Tdomain%sEdge(ne_aus)%Which_Elem(icount) == Tdomain%sFace(nf_aus)%Which_Elem ) then
                        if ( .not.( Tdomain%sPlaneW%pFace(nf)%Orient_Edges(i) == Tdomain%specel(Tdomain%sEdge(ne_aus)%Which_Elem(icount))%Orient_Edges(Tdomain%sEdge(ne_aus)%Which_EdgeinElem(icount)) ) ) then
                            print*,'OrieEdgePW',rg,Tdomain%sEdge(ne_aus)%Which_Elem(icount),ne,ne_aus,Tdomain%sEdge(ne_aus)%Which_EdgeinElem(icount),&
                                Tdomain%specel(Tdomain%sEdge(ne_aus)%Which_Elem(icount))%Orient_Edges(Tdomain%sEdge(ne_aus)%Which_EdgeinElem(icount))
                        endif
                    endif
                enddo

            enddo

        enddo
    endif

    if ( Tdomain%logicD%Neumann_local_present ) then
        if(rg==0) print*,'Nb Faces Neu',Tdomain%sNeu%n_faces
        do nf = 0, Tdomain%sNeu%n_faces-1

            nf_aus = Tdomain%sNeu%nFace(nf)%Face
            do i = 0,3
                ne = Tdomain%sNeu%nFace(nf)%Near_Edges(i)
                ne_aus = Tdomain%sNeu%nEdge(ne)%Edge

                do icount = 0,10
                    if ( Tdomain%sEdge(ne_aus)%Which_Elem(icount) == Tdomain%sFace(nf_aus)%Which_Elem ) then
                        if ( .not.( Tdomain%sNeu%nFace(nf)%Orient_Edges(i) == Tdomain%specel(Tdomain%sEdge(ne_aus)%Which_Elem(icount))%Orient_Edges(Tdomain%sEdge(ne_aus)%Which_EdgeinElem(icount)) ) ) then
                            print*,'OrieEdgeNeu',rg,Tdomain%sEdge(ne_aus)%Which_Elem(icount),ne,ne_aus,Tdomain%sEdge(ne_aus)%Which_EdgeinElem(icount),&
                                Tdomain%specel(Tdomain%sEdge(ne_aus)%Which_Elem(icount))%Orient_Edges(Tdomain%sEdge(ne_aus)%Which_EdgeinElem(icount))
                        endif
                    endif
                enddo

            enddo

        enddo
    endif
    !
    !!
    !purge fuites memoire
    do ne=0,Tdomain%n_edge-1
        deallocate (Tdomain%sEdge(ne)%Which_EdgeinElem)
        deallocate (Tdomain%sEdge(ne)%Which_Elem)
    enddo

    return
end subroutine global_numbering
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
