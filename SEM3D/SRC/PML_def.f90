!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file PML_def.f90
!!\brief La routine PML_definition applique les conditions initiales déclarées dans le fichier d'entrée.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<
subroutine PML_definition (Tdomain)

    use sdomain

    implicit none

    type (Domain), intent (INOUT) :: Tdomain

    integer :: n, i, nf, ne, nv, mat, nff, nfs,nef,nes,nvf,nvs
    logical :: logicPMLf,logicPMLs

    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        Tdomain%specel(n)%PML = .false.
        Tdomain%specel(n)%FPML = .false.
        if(Tdomain%sSubDomain(mat)%material_type == "P" .or.    &
            Tdomain%sSubDomain(mat)%material_type == "L") Tdomain%specel(n)%PML = .true.
        if(Tdomain%specel(n)%PML .and. Tdomain%sSubDomain(mat)%Filtering ) Tdomain%specel(n)%FPML =.true.

    enddo

    do n = 0,Tdomain%n_face-1
        Tdomain%sFace(n)%PML = .true.
        Tdomain%sFace(n)%FPML = .false.
        Tdomain%sFace(n)%Abs = .false.
    enddo
    do n = 0,Tdomain%n_edge-1
        Tdomain%sEdge(n)%PML = .true.
        Tdomain%sEdge(n)%FPML = .false.
        Tdomain%sEdge(n)%Abs = .false.
    enddo
    do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%PML = .true.
        Tdomain%sVertex(n)%FPML = .false.
        Tdomain%sVertex(n)%Abs = .false.
    enddo

    do n = 0,Tdomain%n_elem-1

        if (Tdomain%specel(n)%PML) then
            mat = Tdomain%specel(n)%mat_index
            if (Tdomain%sSubdomain(mat)%Px) then
                if (Tdomain%sSubdomain(mat)%Py) then
                    if (Tdomain%sSubdomain(mat)%Pz) then
                        if (Tdomain%sSubdomain(mat)%Left) then
                            if (Tdomain%sSubdomain(mat)%Forward) then
                                if (Tdomain%sSubdomain(mat)%Down) then
                                    nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                                else
                                    nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                                endif
                            else
                                if (Tdomain%sSubdomain(mat)%Down) then
                                    nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                                else
                                    nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                                endif
                            endif
                        else
                            if (Tdomain%sSubdomain(mat)%Forward) then
                                if (Tdomain%sSubdomain(mat)%Down) then
                                    nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                                else
                                    nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                                endif
                            else
                                if (Tdomain%sSubdomain(mat)%Down) then
                                    nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                                else
                                    nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                                endif
                            endif
                        endif
                    else
                        if (Tdomain%sSubdomain(mat)%Left) then
                            if (Tdomain%sSubdomain(mat)%Forward) then
                                ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                            else
                                ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                            endif
                        else
                            if (Tdomain%sSubdomain(mat)%Forward) then
                                ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                            else
                                ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                            endif
                        endif
                    endif
                else
                    if (Tdomain%sSubdomain(mat)%Pz) then
                        if (Tdomain%sSubdomain(mat)%Left) then
                            if (Tdomain%sSubdomain(mat)%Down) then
                                ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                            else
                                ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                            endif
                        else
                            if (Tdomain%sSubdomain(mat)%Down) then
                                ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                            else
                                ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                            endif
                        endif
                    else
                        if (Tdomain%sSubdomain(mat)%Left) then
                            nf = Tdomain%specel(n)%Near_Faces(2);   Tdomain%sFace(nf)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                        else
                            nf = Tdomain%specel(n)%Near_Faces(4);   Tdomain%sFace(nf)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                        endif
                    endif
                endif
            else
                if (Tdomain%sSubdomain(mat)%Py) then
                    if (Tdomain%sSubdomain(mat)%Pz) then
                        if (Tdomain%sSubdomain(mat)%Forward) then
                            if (Tdomain%sSubdomain(mat)%Down) then
                                ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                            else
                                ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                            endif
                        else
                            if (Tdomain%sSubdomain(mat)%Down) then
                                ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                            else
                                ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                            endif
                        endif
                    else
                        if (Tdomain%sSubdomain(mat)%Forward) then
                            nf = Tdomain%specel(n)%Near_Faces(3);   Tdomain%sFace(nf)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                        else
                            nf = Tdomain%specel(n)%Near_Faces(1);   Tdomain%sFace(nf)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                        endif
                    endif
                else
                    if (Tdomain%sSubdomain(mat)%Pz) then
                        if (Tdomain%sSubdomain(mat)%Down) then
                            nf = Tdomain%specel(n)%Near_Faces(5);   Tdomain%sFace(nf)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%PML = .false.
                        else
                            nf = Tdomain%specel(n)%Near_Faces(0);   Tdomain%sFace(nf)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%PML = .false.
                            ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%PML = .false.
                            nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%PML = .false.
                        endif
                    else
                        print *,"Pb: There is a PML element which is neither Px nor Py nor Pz !?!"
                    endif
                endif
            endif
        else
            do i = 0,5
                nf = Tdomain%specel(n)%Near_Faces(i)
                Tdomain%sFace(nf)%PML = .false.
            enddo
            do i = 0,11
                ne = Tdomain%specel(n)%Near_Edges(i)
                Tdomain%sEdge(ne)%PML = .false.
            enddo
            do i = 0,7
                nv = Tdomain%specel(n)%Near_Vertices(i)
                Tdomain%sVertex(nv)%PML = .false.
            enddo
        endif
    enddo


    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index

        if (.not. Tdomain%specel(n)%PML) cycle

        if (Tdomain%sSubdomain(mat)%Px) then
            if (Tdomain%sSubdomain(mat)%Left) then
                nf = Tdomain%specel(n)%Near_Faces(4);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%Abs = .true.
            else
                nf = Tdomain%specel(n)%Near_Faces(2);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%Abs = .true.
            endif
        endif
        if (Tdomain%sSubdomain(mat)%Py) then
            if (Tdomain%sSubdomain(mat)%Forward) then
                nf = Tdomain%specel(n)%Near_Faces(1);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(4);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(6);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%Abs = .true.
            else
                nf = Tdomain%specel(n)%Near_Faces(3);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(7);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(10);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%Abs = .true.
            endif
        endif
        if (Tdomain%sSubdomain(mat)%Pz) then
            if (Tdomain%sSubdomain(mat)%Down) then
                nf = Tdomain%specel(n)%Near_Faces(0);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(0);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(1);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(2);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(3);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(0);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(1);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(2);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(3);   Tdomain%sVertex(nv)%Abs = .true.
            else
                nf = Tdomain%specel(n)%Near_Faces(5);   Tdomain%sFace(nf)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(5);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(8);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(9);   Tdomain%sEdge(ne)%Abs = .true.
                ne = Tdomain%specel(n)%Near_Edges(11);   Tdomain%sEdge(ne)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(4);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(5);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(6);   Tdomain%sVertex(nv)%Abs = .true.
                nv = Tdomain%specel(n)%Near_Vertices(7);   Tdomain%sVertex(nv)%Abs = .true.
            endif
        endif
    enddo

    do n = 0, Tdomain%n_elem-1
        if  (Tdomain%specel(n)%FPML ) then
            do i = 0, 5
                nf = Tdomain%specel(n)%Near_Faces(i)
                if (Tdomain%sFace(nf)%PML)  Tdomain%sFace(nf)%FPML = .true.
            enddo
            do i = 0, 11
                ne = Tdomain%specel(n)%Near_Edges(i)
                if (Tdomain%sEdge(ne)%PML)  Tdomain%sEdge(ne)%FPML = .true.
            enddo
            do i = 0, 7
                nv = Tdomain%specel(n)%Near_Vertices(i)
                if (Tdomain%sVertex(nv)%PML)  Tdomain%sVertex(nv)%FPML = .true.
            enddo
        endif
    enddo

    ! eventual solid/fluid interfaces
    if(Tdomain%logicD%SF_local_present)then
        do nf = 0,Tdomain%SF%SF_n_faces-1
            logicPMLf = .false. ; logicPMLs = .false.
            nff = Tdomain%SF%SF_face(nf)%Face(0)
            nfs = Tdomain%SF%SF_face(nf)%Face(1)
            if(nff > -1)then
                if(Tdomain%sFace(nff)%PML) logicPMLf = .true.
            end if
            if(nfs > -1)then
                if(Tdomain%sFace(nfs)%PML) logicPMLs = .true.
            end if
            if(logicPMLf .or. logicPMLs) Tdomain%SF%SF_Face(nf)%PML = .true.
        end do
        do ne = 0,Tdomain%SF%SF_n_edges-1
            logicPMLf = .false. ; logicPMLs = .false.
            nef = Tdomain%SF%SF_edge(ne)%Edge(0)
            nes = Tdomain%SF%SF_edge(ne)%Edge(1)
            if(nef > -1)then
                if(Tdomain%sEdge(nef)%PML) logicPMLf = .true.
            end if
            if(nes > -1)then
                if(Tdomain%sEdge(nes)%PML) logicPMLs = .true.
            end if
            if(logicPMLf .or. logicPMLs) Tdomain%SF%SF_Edge(ne)%PML = .true.
        end do
        do nv = 0,Tdomain%SF%SF_n_vertices-1
            logicPMLf = .false. ; logicPMLs = .false.
            nvf = Tdomain%SF%SF_Vertex(nv)%Vertex(0)
            nvs = Tdomain%SF%SF_Vertex(nv)%Vertex(1)
            if(nvf > -1)then
                if(Tdomain%sVertex(nvf)%PML) logicPMLf = .true.
            end if
            if(nvs > -1)then
                if(Tdomain%sVertex(nvs)%PML) logicPMLs = .true.
            end if
            if(logicPMLf .or. logicPMLs) Tdomain%SF%SF_Vertex(nv)%PML = .true.
        end do

    end if


    return
end subroutine PML_definition

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
