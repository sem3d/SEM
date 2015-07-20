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

    ! Local variables
    integer :: n,nel,nv,icount,i,j,ii,jj,nrx,nrz,ngllx,ngllz
    integer ::nface_0, nface_1, nface_2, nface_3
    integer ::nvertex_0, nvertex_1, nvertex_2, nvertex_3
    integer :: nf0,nf1,nf2,nf3,wf0,wf1,wf2,wf3
    integer :: ne0,ne1,ne2,ne3,wv0,wv1,wv2,wv3
    integer, dimension (:), allocatable :: NumElement, WhichVertex

    logical :: change_ij,change_0ngll,lface_0,lface_1,lface_2,lface_3
    logical :: lvertex_0, lvertex_1, lvertex_2, lvertex_3
    logical, dimension (:), allocatable :: LogPointer

    !----------------------------------------------------------------------------
    ! This routine  defines a global numering for the points and an
    ! application between local and global numbering
    !
    ! For convention the indexing of GLL points is related to the
    ! edges of any element in this way
    ! j = 0                        => face (0)
    ! i = ngllx-1                => face (1)
    ! j = ngllz-1                => face (2)
    ! i = 0                        => face (3)
    !
    ! The strategy of this algorithm is the following
    ! 1) If a point is internal to the element, it has to be counted
    ! 2) If the point is on a edge, we search for the neighbor element
    ! 3) If its associated number is smaller, the point was already
    !      numbered, in that a case we search for its number, if the
    !      numer of the near element is greater, it has to be counted
    !
    ! Gaetano Festa, 3 October 2002 / Last Modification 18 October 2004
    ! -----------------------------------------------------------------------------

    nel = Tdomain%n_elem
    nv = Tdomain%n_vertex
    allocate(LogPointer(0:nv-1))
    allocate(NumElement(0:nv-1))
    allocate(WhichVertex(0:nv-1))
    LogPointer = .true.; NumElement = -100; WhichVertex = -100

    icount = 0
    do n = 0,nel-1

        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        allocate (Tdomain%specel(n)%Iglobnum(0:ngllx-1,0:ngllz-1))

        nface_0 = Tdomain%specel(n)%Near_Face(0)
        nface_1 = Tdomain%specel(n)%Near_Face(1)
        nface_2 = Tdomain%specel(n)%Near_Face(2)
        nface_3 = Tdomain%specel(n)%Near_Face(3)

        lface_0 =.false. ; lface_1 =.false.; lface_2 =.false.; lface_3 =.false.

        if  (Tdomain%sFace(nface_0)%Near_Element(0) == n ) then
            lface_0 = .true.
        else
            nf0 = Tdomain%sFace(nface_0)%Near_Element(0)
            wf0 = Tdomain%sFace(nface_0)%Which_Face(0)
        endif

        if  (Tdomain%sFace(nface_1)%Near_Element(0) == n ) then
            lface_1 = .true.
        else
            nf1 = Tdomain%sFace(nface_1)%Near_Element(0)
            wf1 = Tdomain%sFace(nface_1)%Which_Face(0)
        endif
        if  (Tdomain%sFace(nface_2)%Near_Element(0) == n ) then
            lface_2 = .true.
        else
            nf2 = Tdomain%sFace(nface_2)%Near_Element(0)
            wf2 = Tdomain%sFace(nface_2)%Which_Face(0)
        endif

        if  (Tdomain%sFace(nface_3)%Near_Element(0) == n ) then
            lface_3 = .true.
        else
            nf3 = Tdomain%sFace(nface_3)%Near_Element(0)
            wf3 = Tdomain%sFace(nface_3)%Which_Face(0)
        endif

        nvertex_0 = Tdomain%specel(n)%near_vertex(0)
        nvertex_1 = Tdomain%specel(n)%near_vertex(1)
        nvertex_2 = Tdomain%specel(n)%near_vertex(2)
        nvertex_3 = Tdomain%specel(n)%near_vertex(3)
        if (nvertex_0 < 0 .or. nvertex_0 >= nv .or. &
            nvertex_1 < 0 .or. nvertex_1 >= nv .or. &
            nvertex_2 < 0 .or. nvertex_2 >= nv .or. &
            nvertex_3 < 0 .or. nvertex_3 >= nv ) then
            write(*,*) "nvertex_0 ", nvertex_0, "nvertex_1 ", nvertex_1, "nvertex_2 ", nvertex_2, "nvertex_3 ", nvertex_3, ", nv ", nv
            stop "ERROR - global_numbering : near vertex index KO (bad value)"
        end if

        lvertex_0 = .false. ;    lvertex_1 = .false.;    lvertex_2 = .false.;    lvertex_3 = .false.

        if (LogPointer(nvertex_0)) then
            lvertex_0 = .true.
            LogPointer(nvertex_0) = .false.
            NumElement (nvertex_0) = n
            WhichVertex(nvertex_0) = 0
        else
            ne0 = NumElement (nvertex_0)
            wv0 = WhichVertex(nvertex_0)
        endif
        if (LogPointer(nvertex_1)) then
            lvertex_1 = .true.
            LogPointer(nvertex_1) = .false.
            NumElement (nvertex_1) = n
            WhichVertex(nvertex_1) = 1
        else
            ne1 = NumElement (nvertex_1)
            wv1 = WhichVertex(nvertex_1)
        endif
        if (LogPointer(nvertex_2)) then
            lvertex_2 = .true.
            LogPointer(nvertex_2) = .false.
            NumElement (nvertex_2) = n
            WhichVertex(nvertex_2) = 2
        else
            ne2 = NumElement (nvertex_2)
            wv2 = WhichVertex(nvertex_2)
        endif
        if (LogPointer(nvertex_3)) then
            lvertex_3 = .true.
            LogPointer(nvertex_3) = .false.
            NumElement (nvertex_3) = n
            WhichVertex(nvertex_3) = 3
        else
            ne3 = NumElement (nvertex_3)
            wv3 = WhichVertex(nvertex_3)
        endif

        do j = 0,ngllz-1
            do i = 0,ngllx-1

                if (j == 0)  then
                    if (i==0) then
                        if (lvertex_0) then
                            Tdomain%specel(n)%Iglobnum(i,j) = icount
                            icount = icount + 1
                        else
                            if (ne0 == -100 .or. ne0 < 0 .or. ne0 >= nel) then
                                write(*,*) "ne0 ", ne0
                                stop "ERROR - global_numbering : numelem KO (bad init)"
                            end if
                            if (wv0 == -100 .or. wv0 < 0 .or. wv0 >= nv) then
                                write(*,*) "wv0 ", wv0, ", nv ", nv
                                stop "ERROR - global_numbering : whichvertex KO (bad init)"
                            end if

                            call  Vertex_indexing_ij(ii,jj,wv0,ne0,Tdomain)
                            Tdomain%specel(n)%Iglobnum(i,j) = Tdomain%specel(ne0)%Iglobnum(ii,jj)
                        endif
                    else if (i == ngllx-1) then
                        if (lvertex_1) then
                            Tdomain%specel(n)%Iglobnum(i,j) = icount
                            icount = icount + 1
                        else
                            if (ne1 == -100 .or. ne1 < 0 .or. ne1 >= nel) then
                                write(*,*) "ne1 ", ne1
                                stop "ERROR - global_numbering : numelem KO (bad init)"
                            end if
                            if (wv1 == -100 .or. wv1 < 0 .or. wv1 >= nv) then
                                write(*,*) "wv1 ", wv1, ", nv ", nv
                                stop "ERROR - global_numbering : whichvertex KO (bad init)"
                            end if

                            call  Vertex_indexing_ij(ii,jj,wv1,ne1,Tdomain)
                            Tdomain%specel(n)%Iglobnum(i,j) = Tdomain%specel(ne1)%Iglobnum(ii,jj)
                        endif
                    else
                        if (lface_0) then
                            Tdomain%specel(n)%Iglobnum(i,j) = icount
                            icount = icount + 1
                        else
                            call whatface (0,wf0,change_ij,change_0ngll)
                            nrx = Tdomain%specel(nf0)%ngllx
                            nrz = Tdomain%specel(nf0)%ngllz
                            if (change_ij .and. change_0ngll) then
                                if (Tdomain%sFace(nface_0)%coherency) then
                                    Tdomain%specel(n)%Iglobnum(i,0)=Tdomain%specel(nf0)%Iglobnum(nrx-1,i)
                                else
                                    Tdomain%specel(n)%Iglobnum(i,0)=Tdomain%specel(nf0)%Iglobnum(nrx-1,nrz-1-i)
                                endif
                            else if (change_ij .and.  (.not. change_0ngll)) then
                                if (Tdomain%sFace(nface_0)%coherency) then
                                    Tdomain%specel(n)%Iglobnum(i,0)= Tdomain%specel(nf0)%Iglobnum(0,i)
                                else
                                    Tdomain%specel(n)%Iglobnum(i,0)= Tdomain%specel(nf0)%Iglobnum(0,nrz-1-i)
                                endif
                            else if ((.not.change_ij) .and.  change_0ngll) then
                                if (Tdomain%sFace(nface_0)%coherency) then
                                    Tdomain%specel(n)%Iglobnum(i,0)= Tdomain%specel(nf0)%Iglobnum(i,nrz-1)
                                else
                                    Tdomain%specel(n)%Iglobnum(i,0)= Tdomain%specel(nf0)%Iglobnum(nrx-1-i,nrz-1)
                                endif
                            else if ((.not.change_ij) .and.  (.not.change_0ngll)) then
                                if (Tdomain%sFace(nface_0)%coherency) then
                                    Tdomain%specel(n)%Iglobnum(i,0)= Tdomain%specel(nf0)%Iglobnum(i,0)
                                else
                                    Tdomain%specel(n)%Iglobnum(i,0)= Tdomain%specel(nf0)%Iglobnum(nrx-1-i,0)
                                endif
                            endif
                        endif
                    endif
                else if (j==ngllz-1) then
                    if (i==0) then
                        if (lvertex_3) then
                            Tdomain%specel(n)%Iglobnum(i,j) = icount
                            icount = icount + 1
                        else
                            if (ne3 == -100 .or. ne3 < 0 .or. ne3 >= nel) then
                                write(*,*) "ne3 ", ne3
                                stop "ERROR - global_numbering : numelem KO (bad init)"
                            end if
                            if (wv3 == -100 .or. wv3 < 0 .or. wv3 >= nv) then
                                write(*,*) "wv3 ", wv3, ", nv ", nv
                                stop "ERROR - global_numbering : whichvertex KO (bad init)"
                            end if

                            call  Vertex_indexing_ij(ii,jj,wv3,ne3,Tdomain)
                            Tdomain%specel(n)%Iglobnum(i,j) = Tdomain%specel(ne3)%Iglobnum(ii,jj)
                        endif
                    else if (i == ngllx-1) then
                        if (lvertex_2) then
                            Tdomain%specel(n)%Iglobnum(i,j) = icount
                            icount = icount + 1
                        else
                            if (ne2 == -100 .or. ne2 < 0 .or. ne2 >= nel) then
                                write(*,*) "ne2 ", ne2
                                stop "ERROR - global_numbering : numelem KO (bad init)"
                            end if
                            if (wv2 == -100 .or. wv2 < 0 .or. wv2 >= nv) then
                                write(*,*) "wv2 ", wv2, ", nv ", nv
                                stop "ERROR - global_numbering : whichvertex KO (bad init)"
                            end if

                            call  Vertex_indexing_ij(ii,jj,wv2,ne2,Tdomain)
                            Tdomain%specel(n)%Iglobnum(i,j) = Tdomain%specel(ne2)%Iglobnum(ii,jj)
                        endif
                    else
                        if (lface_2) then
                            Tdomain%specel(n)%Iglobnum(i,j) = icount
                            icount = icount + 1
                        else
                            call whatface (2,wf2,change_ij,change_0ngll)
                            nrx = Tdomain%specel(nf2)%ngllx
                            nrz = Tdomain%specel(nf2)%ngllz
                            if (change_ij .and. change_0ngll) then
                                if (Tdomain%sFace(nface_2)%coherency) then
                                    Tdomain%specel(n)%Iglobnum(i,ngllz-1)=Tdomain%specel(nf2)%Iglobnum(0,i)
                                else
                                    Tdomain%specel(n)%Iglobnum(i,ngllz-1)=Tdomain%specel(nf2)%Iglobnum(0,nrz-1-i)
                                endif
                            else if (change_ij .and.  (.not. change_0ngll)) then
                                if (Tdomain%sFace(nface_2)%coherency) then
                                    Tdomain%specel(n)%Iglobnum(i,ngllz-1)= Tdomain%specel(nf2)%Iglobnum(nrx-1,i)
                                else
                                    Tdomain%specel(n)%Iglobnum(i,ngllz-1)= Tdomain%specel(nf2)%Iglobnum(nrx-1,nrz-1-i)
                                endif
                            else if ((.not.change_ij) .and.  change_0ngll) then
                                if (Tdomain%sFace(nface_2)%coherency) then
                                    Tdomain%specel(n)%Iglobnum(i,ngllz-1)= Tdomain%specel(nf2)%Iglobnum(i,0)
                                else
                                    Tdomain%specel(n)%Iglobnum(i,ngllz-1)= Tdomain%specel(nf2)%Iglobnum(nrx-1-i,0)
                                endif
                            else if ((.not.change_ij) .and.  (.not.change_0ngll)) then
                                if (Tdomain%sFace(nface_2)%coherency) then
                                    Tdomain%specel(n)%Iglobnum(i,ngllz-1)= Tdomain%specel(nf2)%Iglobnum(i,nrz-1)
                                else
                                    Tdomain%specel(n)%Iglobnum(i,ngllz-1)= Tdomain%specel(nf2)%Iglobnum(nrx-1-i,nrz-1)
                                endif
                            endif
                        endif
                    endif

                else if (i == ngllx-1)  then
                    if (lface_1) then
                        Tdomain%specel(n)%Iglobnum(i,j) = icount
                        icount = icount + 1
                    else
                        call whatface (1,wf1,change_ij,change_0ngll)
                        nrx = Tdomain%specel(nf1)%ngllx
                        nrz = Tdomain%specel(nf1)%ngllz
                        if (change_ij .and. change_0ngll) then
                            if (Tdomain%sFace(nface_1)%coherency) then
                                Tdomain%specel(n)%Iglobnum(ngllx-1,j)=Tdomain%specel(nf1)%Iglobnum(j,0)
                            else
                                Tdomain%specel(n)%Iglobnum(ngllx-1,j)=Tdomain%specel(nf1)%Iglobnum(nrx-1-j,0)
                            endif
                        else if (change_ij .and.  (.not. change_0ngll)) then
                            if (Tdomain%sFace(nface_1)%coherency) then
                                Tdomain%specel(n)%Iglobnum(ngllx-1,j)= Tdomain%specel(nf1)%Iglobnum(j,nrz-1)
                            else
                                Tdomain%specel(n)%Iglobnum(ngllx-1,j)= Tdomain%specel(nf1)%Iglobnum(nrx-1-j,nrz-1)
                            endif
                        else if ((.not.change_ij) .and.  change_0ngll) then
                            if (Tdomain%sFace(nface_1)%coherency) then
                                Tdomain%specel(n)%Iglobnum(ngllx-1,j)= Tdomain%specel(nf1)%Iglobnum(0,j)
                            else
                                Tdomain%specel(n)%Iglobnum(ngllx-1,j)= Tdomain%specel(nf1)%Iglobnum(0,nrz-1-j)
                            endif
                        else if ((.not.change_ij) .and.  (.not.change_0ngll)) then
                            if (Tdomain%sFace(nface_1)%coherency) then
                                Tdomain%specel(n)%Iglobnum(ngllx-1,j)= Tdomain%specel(nf1)%Iglobnum(nrx-1,j)
                            else
                                Tdomain%specel(n)%Iglobnum(ngllx-1,j)= Tdomain%specel(nf1)%Iglobnum(nrx-1,nrz-1-j)
                            endif
                        endif
                    endif

                else if (i ==0 )  then
                    if (lface_3) then
                        Tdomain%specel(n)%Iglobnum(i,j) = icount
                        icount = icount + 1
                    else
                        call whatface (3,wf3,change_ij,change_0ngll)
                        nrx = Tdomain%specel(nf3)%ngllx
                        nrz = Tdomain%specel(nf3)%ngllz
                        if (change_ij .and. change_0ngll) then
                            if (Tdomain%sFace(nface_3)%coherency) then
                                Tdomain%specel(n)%Iglobnum(0,j)=Tdomain%specel(nf3)%Iglobnum(j,nrz-1)
                            else
                                Tdomain%specel(n)%Iglobnum(0,j)=Tdomain%specel(nf3)%Iglobnum(nrx-1-j,nrz-1)
                            endif
                        else if (change_ij .and.  (.not. change_0ngll)) then
                            if (Tdomain%sFace(nface_3)%coherency) then
                                Tdomain%specel(n)%Iglobnum(0,j)= Tdomain%specel(nf3)%Iglobnum(j,0)
                            else
                                Tdomain%specel(n)%Iglobnum(0,j)= Tdomain%specel(nf3)%Iglobnum(nrx-1-j,0)
                            endif
                        else if ((.not.change_ij) .and.  change_0ngll) then
                            if (Tdomain%sFace(nface_3)%coherency) then
                                Tdomain%specel(n)%Iglobnum(0,j)= Tdomain%specel(nf3)%Iglobnum(nrx-1,j)
                            else
                                Tdomain%specel(n)%Iglobnum(0,j)= Tdomain%specel(nf3)%Iglobnum(nrx-1,nrz-1-j)
                            endif
                        else if ((.not.change_ij) .and.  (.not.change_0ngll)) then
                            if (Tdomain%sFace(nface_3)%coherency) then
                                Tdomain%specel(n)%Iglobnum(0,j)= Tdomain%specel(nf3)%Iglobnum(0,j)
                            else
                                Tdomain%specel(n)%Iglobnum(0,j)= Tdomain%specel(nf3)%Iglobnum(0,nrz-1-j)
                            endif
                        endif
                    endif

                else
                    Tdomain%specel(n)%Iglobnum(i,j) = icount
                    icount = icount + 1
                endif

            enddo
        enddo
    enddo

    !
    Tdomain%n_glob_points = icount
    deallocate (LogPointer); deallocate (NumElement); deallocate (WhichVertex)
    return
end subroutine global_numbering

! #################################################
!>
!! \brief
!!
!! \param integer, intent (OUT) i
!! \param integer, intent (OUT) j
!! \param integer, intent (IN) wv
!! \param integer, intent (IN) ne
!! \param type (domain), intent (IN) Tdomain
!<


subroutine Vertex_indexing_ij(i,j,wv,ne,Tdomain)
    use sdomain
    implicit none

    integer, intent (OUT) :: i,j
    integer, intent (IN) :: wv,ne
    type (domain), intent (IN) :: Tdomain

    if (wv == 0) then
        i = 0; j= 0
    else if ( wv == 1) then
        i = Tdomain%specel(ne)%ngllx-1; j = 0
    else if ( wv == 2) then
        i = Tdomain%specel(ne)%ngllx-1; j = Tdomain%specel(ne)%ngllz-1
    else if ( wv == 3) then
        i = 0; j = Tdomain%specel(ne)%ngllz-1
    endif
    return
end subroutine  Vertex_indexing_ij

! #####################################################
!>
!! \brief
!!
!! \param integer, intent (IN) n1
!! \param integer, intent (IN) n2
!! \param logical, intent (OUT) invert_ij
!! \param logical, intent (OUT) invert_0ngll
!<


subroutine whatface (n1,n2,invert_ij,invert_0ngll)
    implicit none
    integer, intent (IN) :: n1,n2
    logical, intent (OUT) :: invert_ij,invert_0ngll

    ! Local variables
    integer :: nabs
    !-------------------------------------------------------------
    nabs = abs (n1-n2)
    if (nabs == 0) then
        invert_ij = .false.
        invert_0ngll = .false.
    else if ( nabs == 1 ) then
        invert_ij = .true.
        select case (n1)
        case (0)
            invert_0ngll  = .true.
        case (1)
            invert_0ngll= .true.
            if (n2 == 2 ) invert_0ngll=.false.
        case (2)
            invert_0ngll  = .true.
            if (n2 == 1 ) invert_0ngll=.false.
        case (3)
            invert_0ngll  = .true.
        end select
    else if (nabs == 2) then
        invert_ij = .false.
        invert_0ngll = .true.
    else if ( nabs == 3 ) then
        invert_ij = .true.
        invert_0ngll = .false.
    endif
    return
end subroutine whatface
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
