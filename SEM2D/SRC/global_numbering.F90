!>
!!\file global_numbering.F90
!!\brief Assure la correspondance entre les différentes numérotations.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! Effectue la numérotation globale de tous les points de Gauss.
!! Donne la numérotation des points de Gauss d'une maille
!! Definition de Iglobnum.
!<
subroutine global_numbering ( Tdomain )
    use sdomain
    implicit none

    type(domain), target, intent (INOUT) :: Tdomain

    integer :: icount, n, ngllx, ngllz, i, k, ngll, v0, v1, nf

    logical :: goDie ! Debug
    goDie=.false.

    ! Corner GLL points
    icount = 0
    do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%Iglobnum_Vertex = icount
        icount = icount + 1
    enddo

    ! Faces Inner GLL points
    do n = 0,Tdomain%n_face-1
        ngll = Tdomain%sFace(n)%ngll
        allocate(Tdomain%sFace(n)%Iglobnum_Face(0:ngll-1))
        Tdomain%sFace(n)%Iglobnum_Face = -1
        do i = 1,ngll-2
            Tdomain%sFace(n)%Iglobnum_Face(i) = icount
            icount = icount + 1
        enddo
        ! Vertex 0
        v0 = Tdomain%sFace(n)%Near_Vertex(0)
        Tdomain%sFace(n)%Iglobnum_Face(0) = Tdomain%sVertex(v0)%Iglobnum_Vertex
        ! Vertex 1
        v1 = Tdomain%sFace(n)%Near_Vertex(1)
        Tdomain%sFace(n)%Iglobnum_Face(ngll-1) = Tdomain%sVertex(v1)%Iglobnum_Vertex
    enddo

    ! Elements Inner GLL points
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        allocate(Tdomain%specel(n)%Iglobnum(0:ngllx-1,0:ngllz-1))
        Tdomain%specel(n)%Iglobnum = -1
        do k = 1,ngllz-2
            do i = 1,ngllx-2
                Tdomain%specel(n)%Iglobnum(i,k) = icount
                icount = icount + 1
            enddo
        enddo
        ! Face 0 : k = 0
        nf = Tdomain%specel(n)%Near_Face(0)
        if ( Tdomain%sFace(nf)%Near_Element(1) == n .and. .not. Tdomain%sFace(nf)%coherency ) then
            do i = 0, ngllx-1
                Tdomain%specel(n)%Iglobnum(i,0) = Tdomain%sFace(nf)%Iglobnum_Face(ngllx-1-i)
            end do
        else
            do i = 0, ngllx-1
                Tdomain%specel(n)%Iglobnum(i,0) = Tdomain%sFace(nf)%Iglobnum_Face(i)
            end do
        endif
        ! Face 1 : i = ngllx-1
        nf = Tdomain%specel(n)%Near_Face(1)
        if ( Tdomain%sFace(nf)%Near_Element(1) == n .and. .not. Tdomain%sFace(nf)%coherency ) then
            if ( Tdomain%specel(n)%Iglobnum(ngllx-1,0) /= Tdomain%sFace(nf)%Iglobnum_Face(ngllz-1) ) stop "global_numbering KO"
            do k = 0, ngllz-1
                Tdomain%specel(n)%Iglobnum(ngllx-1,k) = Tdomain%sFace(nf)%Iglobnum_Face(ngllz-1-k)
            end do
        else
            if ( Tdomain%specel(n)%Iglobnum(ngllx-1,0) /= Tdomain%sFace(nf)%Iglobnum_Face(0) ) then
                do k = 0,ngllz-1
                    do i = 0,ngllx-1
                        write(*,'(i6)', advance='no') Tdomain%specel(n)%Iglobnum(i,k)
                    enddo
                    write(*,*) " "
                enddo
                write(*,*) "---------"
                write(*,*) Tdomain%sFace(nf)%Iglobnum_Face
                write(*,*) "---------"
                write(*,*) n, " : ", Tdomain%sFace(nf)%Near_Element , Tdomain%sFace(nf)%coherency
                write(*,*) "---------"
                write(*,*) Tdomain%sFace(nf)%which_face

                stop "global_numbering KO"
                goDie=.true.
            endif
            do k = 0, ngllz-1
                Tdomain%specel(n)%Iglobnum(ngllx-1,k) = Tdomain%sFace(nf)%Iglobnum_Face(k)
            end do
            if ( goDie ) then
                do k = 0,ngllz-1
                    do i = 0,ngllx-1
                        write(*,'(i6)', advance='no') Tdomain%specel(n)%Iglobnum(i,k)
                    enddo
                    write(*,*) " "
                enddo
                write(*,*) "---------"
                write(*,*) Tdomain%sFace(nf)%Iglobnum_Face
                write(*,*) "---------"
                write(*,*) n, " : ", Tdomain%sFace(nf)%Near_Element , Tdomain%sFace(nf)%coherency
                write(*,*) "---------"
                write(*,*) Tdomain%sFace(nf)%which_face

                stop "global_numbering KO"
            endif
        endif
        ! Face 2 : k = ngllz-1
        nf = Tdomain%specel(n)%Near_Face(2)
        if ( Tdomain%sFace(nf)%Near_Element(1) == n .and. .not. Tdomain%sFace(nf)%coherency ) then
            if ( Tdomain%specel(n)%Iglobnum(ngllx-1,ngllz-1) /= Tdomain%sFace(nf)%Iglobnum_Face(0) ) stop "global_numbering KO"
            do i = 0, ngllx-1
                Tdomain%specel(n)%Iglobnum(i,ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(ngllx-1-i)
            end do
        else
            if ( Tdomain%specel(n)%Iglobnum(ngllx-1,ngllz-1) /= Tdomain%sFace(nf)%Iglobnum_Face(ngllx-1) ) stop "global_numbering KO"
            do i = 0, ngllx-1
                Tdomain%specel(n)%Iglobnum(i,ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i)
            end do
        endif
        ! Face 3 : i = 0
        nf = Tdomain%specel(n)%Near_Face(3)
        if ( Tdomain%sFace(nf)%Near_Element(1) == n .and. .not. Tdomain%sFace(nf)%coherency ) then
            if ( Tdomain%specel(n)%Iglobnum(0,0) /= Tdomain%sFace(nf)%Iglobnum_Face(ngllz-1) ) stop "global_numbering KO"
            if ( Tdomain%specel(n)%Iglobnum(0,ngllz-1) /= Tdomain%sFace(nf)%Iglobnum_Face(0) ) stop "global_numbering KO"
            do k = 0, ngllz-1
                Tdomain%specel(n)%Iglobnum(0,k) = Tdomain%sFace(nf)%Iglobnum_Face(ngllz-1-k)
            end do
        else
            if ( Tdomain%specel(n)%Iglobnum(0,0) /= Tdomain%sFace(nf)%Iglobnum_Face(0) ) stop "global_numbering KO"
            if ( Tdomain%specel(n)%Iglobnum(0,ngllz-1) /= Tdomain%sFace(nf)%Iglobnum_Face(ngllz-1) ) stop "global_numbering KO"
            do k = 0, ngllz-1
                Tdomain%specel(n)%Iglobnum(0,k) = Tdomain%sFace(nf)%Iglobnum_Face(k)
            end do
        endif
    enddo
    Tdomain%n_glob_points = icount

    ! Debug
    do n = 0,Tdomain%n_elem-1
        write(*,*) " ================= ", n
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        do k = 0,ngllz-1
            do i = 0,ngllx-1
                write(*,'(i6)', advance='no') Tdomain%specel(n)%Iglobnum(i,k)
            enddo
            write(*,*) " "
        enddo
    enddo
end subroutine global_numbering

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
