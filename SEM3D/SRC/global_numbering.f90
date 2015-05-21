!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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
  private :: renumber_domains, renumber_sf_coupling, renumber_pml_domain
  private :: prepare_comm_vector, prepare_comm_vector_SF
contains

subroutine global_numbering(Tdomain)

    ! routine different from the 2D case. Everything is independently numbered, here (inner
    !      points in elements, on faces, edges and vertices). And then associated in
    !      the "Iglobnum" field = global number of each GLL point.


    use sdomain
    implicit none
    type(domain), intent (inout) :: Tdomain
    integer, dimension (:), allocatable :: renumS, renumF, renumSpml, renumFpml
    integer, dimension (:,:), allocatable :: renumSF

    ! Create a unique gll number per node location
    call renumber_global_gll_nodes(Tdomain)

    allocate(renumS(0:Tdomain%n_glob_points-1))
    allocate(renumF(0:Tdomain%n_glob_points-1))
    allocate(renumSpml(0:Tdomain%n_glob_points-1))
    allocate(renumFpml(0:Tdomain%n_glob_points-1))
    allocate(renumSF(0:Tdomain%n_glob_points-1,0:2))
    renumS = -1
    renumF = -1
    renumSpml = -1
    renumFpml = -1
    renumSF = -1

    ! Create individual numbering for each sub(physical) domain
    call renumber_domains(Tdomain, renumS, renumF, renumSpml, renumFpml)
    ! Create numbering for Solid-Fluid interface
    call renumber_sf_coupling(Tdomain, renumSF)
    !
    call renumber_pml_domain(Tdomain, renumS, renumF, renumSpml, renumFpml)
    ! Compute index of glls that participate in communications
    call prepare_comm_vector(Tdomain,Tdomain%Comm_data)
    call prepare_comm_vector_SF(Tdomain,Tdomain%n_glob_points,renumSF,Tdomain%Comm_SolFlu)
!     call debug_comm_vector(Tdomain, rank, 0, 1, Tdomain%Comm_data)
!     call debug_comm_vector(Tdomain, rank, 0, 2, Tdomain%Comm_data)
!     call debug_comm_vector(Tdomain, rank, 0, 3, Tdomain%Comm_data)
!     call debug_comm_vector(Tdomain, rank, 3, 0, Tdomain%Comm_data)

    deallocate(renumS)
    deallocate(renumSpml)
    deallocate(renumF)
    deallocate(renumFpml)

end subroutine global_numbering

subroutine renumber_global_gll_nodes(Tdomain)
    use sdomain
    use mindex
    implicit none
    type(domain), intent (inout) :: Tdomain
    integer :: n, icount, i, j, k, ngllx, nglly, ngllz, nf, nnf, ne, nne, nv, ngll1, ngll2,   &
        orient_f, orient_e, ngll, nnv
    integer, dimension(0:6)  :: index_elem_f
    integer, dimension(0:4)  :: index_elem_e
    integer, dimension(0:2)  :: index_elem_v
    ! Counts PML glls with abs flag in solid or fluid
    integer :: solid_abs_count, fluid_abs_count

	!Elements Inner GLL points
    icount = 0
    solid_abs_count = 0
    fluid_abs_count = 0

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

	!Faces Inner GLL points
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
        if (Tdomain%sFace(n)%Abs .and. Tdomain%sFace(n)%PML) then
            if (Tdomain%sFace(n)%Solid) then
                solid_abs_count = solid_abs_count + (ngllx-2) * (nglly-2)
            else
                fluid_abs_count = fluid_abs_count + (ngllx-2) * (nglly-2)
            endif
        endif
    enddo

	!Edges Inner GLL points
    do n = 0,Tdomain%n_edge-1
        ngllx = Tdomain%sEdge(n)%ngll
        allocate(Tdomain%sEdge(n)%Iglobnum_Edge(1:ngllx-2))
        Tdomain%sEdge(n)%Iglobnum_Edge = -1
        do i = 1,ngllx-2
            Tdomain%sEdge(n)%Iglobnum_Edge(i) = icount
            icount = icount + 1
        enddo
        if (Tdomain%sEdge(n)%Abs .and. Tdomain%sEdge(n)%PML) then
            if (Tdomain%sEdge(n)%Solid) then
                solid_abs_count = solid_abs_count + (ngllx-2)
            else
                fluid_abs_count = fluid_abs_count + (ngllx-2)
            endif
        endif
    enddo

	!Corner GLL points
    do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%Iglobnum_Vertex = icount
        icount = icount + 1
        if (Tdomain%sVertex(n)%Abs .and. Tdomain%sVertex(n)%PML) then
            if (Tdomain%sVertex(n)%Solid) then
                solid_abs_count = solid_abs_count + 1
            else
                fluid_abs_count = fluid_abs_count + 1
            endif
        endif
    enddo

    Tdomain%nbOuterSPMLNodes = solid_abs_count
    Tdomain%nbOuterFPMLNodes = fluid_abs_count
    allocate(Tdomain%OuterSPMLNodes(0:Tdomain%nbOuterSPMLNodes-1))
    allocate(Tdomain%OuterFPMLNodes(0:Tdomain%nbOuterFPMLNodes-1))
    ! total number of GLL points (= degrees of freedom)
    Tdomain%n_glob_points = icount

    !Recollecting at the element level, from faces, edges and vertices.
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        !Taking information from faces
        do nf = 0,5
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            call ind_elem_face(nf,orient_f,ngllx,nglly,ngllz,index_elem_f)
            select case(orient_f)
            case(0,1,2,3)
                if(nf == 2 .or. nf == 4)then
                    Tdomain%specel(n)%Iglobnum(                             &
                        index_elem_f(0),                                    &
                        index_elem_f(1):index_elem_f(2):index_elem_f(3),    &
                        index_elem_f(4):index_elem_f(5):index_elem_f(6))    &
                        = Tdomain%sFace(nnf)%Iglobnum_Face(1:ngll1-2,1:ngll2-2)
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


        !Taking information from edges
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

        !Taking information from vertices
        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            call ind_elem_vertex(nv,ngllx,nglly,ngllz,index_elem_v)
            Tdomain%specel(n)%Iglobnum(index_elem_v(0),index_elem_v(1),index_elem_v(2)) =   &
                Tdomain%sVertex(nnv)%Iglobnum_Vertex
        end do

    enddo    ! end of the loop onto elements

end subroutine renumber_global_gll_nodes


subroutine renumber_domains(Tdomain, renumS, renumF, renumSpml, renumFpml)
    use sdomain
    use mindex
    implicit none
    type(domain), intent (inout) :: Tdomain
    integer, dimension(0:), intent(inout) :: renumS, renumF, renumSpml, renumFpml

    integer :: n,  i, j
    integer :: ngllx, nglly, ngllz
    integer :: idxS, idxF, idxSpml, idxFpml, idxSF
    integer :: nbPtInterfSolPml, nbPtInterfFluPml
    !
    ! Renumerotation
    idxS = 0
    idxF = 0
    idxSpml = 0
    idxFpml = 0
    idxSF = 0
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        if (Tdomain%specel(n)%Solid)then
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

    Tdomain%nbInterfSolPml = 0
    Tdomain%nbInterfFluPml = 0
    if (Tdomain%any_PML) then
        nbPtInterfSolPml = 0
        nbPtInterfFluPml = 0
        ! On créé l'interface de couplage Sol/PML
        do n = 0,Tdomain%n_glob_points-1
            if (renumS(n) > -1 .and. renumSpml(n) > -1) then
                ! On se trouve ầ l'interface Sol/PML : on compte le nombre de points à l'interface
                nbPtInterfSolPml = nbPtInterfSolPml + 1
            endif
            if (renumF(n) > -1 .and. renumFpml(n) > -1) then
                ! On se trouve ầ l'interface Flu/PML : on compte le nombre de points à l'interface
                nbPtInterfFluPml = nbPtInterfFluPml + 1
            endif
            ! On Elimine les cas tordus (maille pml fluide en contact avec solide et vice-versa)
            if ((renumF(n) > -1) .and. (renumSpml(n) > -1) .and. (renumFpml(n)==-1)) then
                stop "Erreur maille pml solide en face d'une maille fluide"
            endif
            if ((renumS(n) > -1) .and. (renumFpml(n) > -1) .and. (renumSpml(n)==-1)) then
                stop "Erreur maille pml fluide en face d'une maille solide"
            endif
        enddo
        Tdomain%nbInterfSolPml = nbPtInterfSolPml
        Tdomain%nbInterfFluPml = nbPtInterfFluPml
        allocate(Tdomain%InterfSolPml(0:nbPtInterfSolPml-1,0:1))
        allocate(Tdomain%InterfFluPml(0:nbPtInterfFluPml-1,0:1))
        i = 0
        j = 0
        do n = 0,Tdomain%n_glob_points-1
            if (renumS(n) > -1 .and. renumSpml(n) > -1) then
                ! On se trouve à l'interface Sol/PML :
                Tdomain%InterfSolPml(i,0) = renumS(n)
                Tdomain%InterfSolPml(i,1) = renumSpml(n)
                i = i + 1
            endif
            if (renumF(n) > -1 .and. renumFpml(n) > -1) then
                ! On se trouve à l'interface Sol/PML :
                Tdomain%InterfFluPml(j,0) = renumF(n)
                Tdomain%InterfFluPml(j,1) = renumFpml(n)
                j = j + 1
            endif
        enddo
    endif

    Tdomain%ngll_s = idxS
    Tdomain%ngll_f = idxF
    Tdomain%ngll_pmls = idxSpml
    Tdomain%ngll_pmlf = idxFpml

end subroutine renumber_domains


subroutine renumber_sf_coupling(Tdomain, renumSF)
    use sdomain
    use mindex
    implicit none
    type(domain), intent (inout) :: Tdomain
    integer, dimension (0:,0:), intent(inout) :: renumSF
    !
    integer :: n, i, j, k, ngllx, nglly, ngllz, nf, nnf,  ngll1, ngll2
    integer :: idx,   dir, ks, kl
    integer :: sol, flu
    logical :: saut
    !! Couplage solide/fluide
    !! On liste les numeraux globaux des points de gauss concernant les faces de couplage S/F
    if(Tdomain%logicD%SF_local_present)then
        call renumSolFlu(Tdomain)

        k = 0
        allocate(Tdomain%SF%SF_IGlob(0:Tdomain%SF%ngll-1,0:2))
        Tdomain%SF%SF_IGlob = -2
        do nf = 0,Tdomain%SF%SF_n_faces-1
            ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
            ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
            do j = 0,Tdomain%sFace(nf)%ngll2-1
                do i = 0,Tdomain%sFace(nf)%ngll1-1
                    idx = Tdomain%SF%SF_Face(nf)%IG(i,j,0)
                    if (idx == -1) then
                        stop "Error (global_numbering): numerotation incorrect interface Sol/Flu"
                    endif
                    flu = Tdomain%SF%SF_Face(nf)%IG(i,j,1)
                    sol = Tdomain%SF%SF_Face(nf)%IG(i,j,2)
                    Tdomain%SF%SF_IGlob(idx,0) = flu
                    Tdomain%SF%SF_IGlob(idx,1) = sol
                    if (flu == -1 .or. sol == -1) then
                        ! Index of interface Sol/Flu
                        Tdomain%SF%SF_IGlob(idx,2) = k
                        !Comm_data%Data(n)%IGiveS(k)=idx
                        k = k + 1
                    endif
                enddo
            enddo
        enddo






        ! On peuple le tableau d'indice des ngll de couplage en fonction des direcions des faces
        k = 0
        ks = 0
        kl = 0
        saut = .false.
        do nf = 0,Tdomain%SF%SF_n_faces-1
            ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
            ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
            allocate(Tdomain%SF%SF_Face(nf)%I_sf(0:ngll1-1,0:ngll2-1))
            Tdomain%SF%SF_Face(nf)%I_sf = -1
            ! Partie fluide
            nnf = Tdomain%SF%SF_Face(nf)%Face(0)
            if(nnf > -1) then
                n = Tdomain%sFace(nnf)%which_elem
                dir = Tdomain%sFace(nnf)%dir
                ngllx = Tdomain%specel(n)%ngllx
                nglly = Tdomain%specel(n)%nglly
                ngllz = Tdomain%specel(n)%ngllz

                !write(*,*) "DEBUG: rank, which_elem fluide", rank, n, dir
                call populate_index_SF(1, dir, ngllx, nglly, ngllz, &
                                       Tdomain%specel(n)%IFlu,  &
                                       Tdomain%specel(n)%IGlobnum, kl, &
                                       Tdomain%n_glob_points, renumSF, Tdomain%SF%SF_Face(nf)%I_sf)
            else
                write(*,*) "ZZZZZZ Face fluide sur autre proc", Tdomain%rank, nf, Tdomain%SF%SF_Face(nf)%Face(1) 
                saut = .true.
            end if

            ! Partie Solide
            nnf = Tdomain%SF%SF_Face(nf)%Face(1)
            if(nnf > -1) then
                n = Tdomain%sFace(nnf)%which_elem
                dir = Tdomain%sFace(nnf)%dir
                ngllx = Tdomain%specel(n)%ngllx
                nglly = Tdomain%specel(n)%nglly
                ngllz = Tdomain%specel(n)%ngllz

                !write(*,*) "DEBUG: rank, which_elem solide", rank, n, dir
                call populate_index_SF(2, dir, ngllx, nglly, ngllz, &
                                       Tdomain%specel(n)%ISol,  &
                                       Tdomain%specel(n)%IGlobnum, ks, &
                                       Tdomain%n_glob_points, renumSF, Tdomain%SF%SF_Face(nf)%I_sf)
            else
                write(*,*) "ZZZZZZ Face solide sur autre proc", Tdomain%rank, nf, Tdomain%SF%SF_Face(nf)%Face(0)
                saut = .true.
            end if

            if (saut) then
                k = max(ks,kl)
                ks = k
                kl = k
                saut = .false.
            endif
        enddo

        ! On alloue et on initialise les tableaux d'indice des points de couplage
        k = max(ks,kl)
        Tdomain%SF%ngll = k
        allocate(Tdomain%SF%SF_IGlobSol(0:k-1))
        allocate(Tdomain%SF%SF_IGlobFlu(0:k-1))
        Tdomain%SF%SF_IGlobSol = -1
        Tdomain%SF%SF_IGlobFlu = -1
        do k=0,Tdomain%SF%ngll-1
            Tdomain%SF%SF_IGlobFlu(k) = renumSF(k,1)
            Tdomain%SF%SF_IGlobSol(k) = renumSF(k,2)
        end do
    endif ! fin traitement couplage Solide / Fluide

end subroutine renumber_sf_coupling


subroutine renumber_pml_domain(Tdomain, renumS, renumF, renumSpml, renumFpml)
    use sdomain
    use mindex
    implicit none
    type(domain), intent (inout) :: Tdomain
    integer, dimension (0:), intent(in) :: renumS, renumF, renumSpml, renumFpml
    integer :: n, i, j, ngll1, ngll2

    integer :: idxSpml, idxFpml
    integer :: solid_abs_count, fluid_abs_count

    solid_abs_count = 0
    fluid_abs_count = 0
    do n = 0,Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1
        ngll2 = Tdomain%sFace(n)%ngll2
        if (Tdomain%sFace(n)%Abs .and. Tdomain%sFace(n)%PML) then
            do j = 1,ngll2-2
                do i = 1,ngll1-2
                    if (Tdomain%sFace(n)%Solid) then
                        idxSpml = renumSpml(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (idxSpml .EQ. -1) stop "Unexpected non pml outer Face !"
                        Tdomain%OuterSPMLNodes(solid_abs_count) = idxSpml
                        solid_abs_count = solid_abs_count + 1
                    else
                        idxFpml = renumFpml(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (idxFpml .EQ. -1) stop "Unexpected non pml outer Face !"
                        Tdomain%OuterFPMLNodes(fluid_abs_count) = idxFpml
                        fluid_abs_count = fluid_abs_count + 1
                    endif
                enddo
            enddo
        endif

        allocate(Tdomain%sFace(n)%Renum(1:ngll1-2,1:ngll2-2))
        do j = 1,ngll2-2
            do i = 1,ngll1-2
                if (Tdomain%sFace(n)%solid) then
                    if (Tdomain%sFace(n)%PML) then
                        Tdomain%sFace(n)%Renum(i,j) = renumSpml(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (Tdomain%sFace(n)%Renum(i,j) < 0) stop "Error renumbering SOLID PML faces"
                    else
                        Tdomain%sFace(n)%Renum(i,j) = renumS(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (Tdomain%sFace(n)%Renum(i,j) < 0) stop "Error renumbering SOLID faces"
                    endif
                else
                    if (Tdomain%sFace(n)%PML) then
                        Tdomain%sFace(n)%Renum(i,j) = renumFpml(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (Tdomain%sFace(n)%Renum(i,j) < 0) stop "Error renumbering FLUID PML faces"
                    else
                        Tdomain%sFace(n)%Renum(i,j) = renumF(Tdomain%sFace(n)%Iglobnum_Face(i,j))
                        if (Tdomain%sFace(n)%Renum(i,j) < 0) stop "Error renumbering FLUID faces"
                    endif
                endif
            enddo
        enddo
    enddo

    do n = 0,Tdomain%n_edge-1
        ngll1 = Tdomain%sEdge(n)%ngll
        if (Tdomain%sEdge(n)%Abs .and. Tdomain%sEdge(n)%PML) then
            do i = 1,ngll1-2
                if (Tdomain%sEdge(n)%Solid) then
                    idxSpml = renumSpml(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                    if (idxSpml .EQ. -1) stop "Unexpected non pml outer Edge !"
                    Tdomain%OuterSPMLNodes(solid_abs_count) = idxSpml
                    solid_abs_count = solid_abs_count + 1
                else
                    idxFpml = renumFpml(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                    if (idxFpml .EQ. -1) stop "Unexpected non pml outer Edge !"
                    Tdomain%OuterFPMLNodes(fluid_abs_count) = idxFpml
                    fluid_abs_count = fluid_abs_count + 1
                endif
            enddo
        endif

        allocate(Tdomain%sEdge(n)%Renum(1:ngll1-2))
        Tdomain%sEdge(n)%Renum(:) = -1
        do i = 1,ngll1-2
            if (Tdomain%sEdge(n)%solid) then
                if (Tdomain%sEdge(n)%PML) then
                    Tdomain%sEdge(n)%Renum(i) = renumSpml(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                else
                    Tdomain%sEdge(n)%Renum(i) = renumS(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                endif
            else
                if (Tdomain%sEdge(n)%PML) then
                    Tdomain%sEdge(n)%Renum(i) = renumFpml(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                else
                    Tdomain%sEdge(n)%Renum(i) = renumF(Tdomain%sEdge(n)%Iglobnum_Edge(i))
                endif
            endif
        enddo
    enddo

    do n = 0,Tdomain%n_vertex-1
        if (Tdomain%sVertex(n)%Abs .and. Tdomain%sVertex(n)%PML) then
            if (Tdomain%sVertex(n)%solid) then
                idxSpml = renumSpml(Tdomain%sVertex(n)%Iglobnum_Vertex)
                if (idxSpml .EQ. -1) stop "Unexpected non pml outer Vertex !"
                Tdomain%OuterSPMLNodes(solid_abs_count) = idxSpml
                solid_abs_count = solid_abs_count + 1
            else
                idxFpml = renumFpml(Tdomain%sVertex(n)%Iglobnum_Vertex)
                if (idxFpml .EQ. -1) stop "Unexpected non pml outer Vertex !"
                Tdomain%OuterFPMLNodes(fluid_abs_count) = idxFpml
                fluid_abs_count = fluid_abs_count + 1
            endif
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

    return
end subroutine renumber_pml_domain

subroutine prepare_comm_vector(Tdomain,comm_data)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain
    type(comm_vector), intent(inout) :: comm_data

    integer :: n,ncomm,nsol,nsolpml,nflu,nflupml
    integer :: i,j,k,nf,ne,nv,idx,ngll1,ngll2

    ! Remplissage des IGive et ITake
    if(Tdomain%nb_procs < 2)then
        Comm_data%ncomm = 0
        return
    endif

    call allocate_comm_vector(Tdomain,comm_data)

    do n = 0,Comm_data%ncomm-1
            ncomm = Comm_data%Data(n)%ncomm
            nsol = 0
            nsolpml = 0
            nflu = 0
            nflupml = 0

            ! Remplissage des Igive
            ! Faces
            do i = 0,Tdomain%sComm(ncomm)%nb_faces-1
                nf = Tdomain%sComm(ncomm)%faces(i)
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        idx = Tdomain%sFace(nf)%Renum(k,j)
                        if (Tdomain%sFace(nf)%solid) then
                            if (Tdomain%sFace(nf)%PML) then
                                Comm_data%Data(n)%IGiveSPML(nsolpml) = idx
                                nsolpml = nsolpml + 1
                            else
                                Comm_data%Data(n)%IGiveS(nsol) = idx
                                nsol = nsol + 1
                            endif
                        else ! Fluid
                            if (Tdomain%sFace(nf)%PML) then
                                Comm_data%Data(n)%IGiveFPML(nflupml) = idx
                                nflupml = nflupml + 1
                            else
                                Comm_data%Data(n)%IGiveF(nflu) = idx
                                nflu = nflu + 1
                            endif
                        endif
                    end do
                end do
            enddo

            ! Edges
            do i = 0,Tdomain%sComm(ncomm)%nb_edges-1
                ne = Tdomain%sComm(ncomm)%edges(i)
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    idx = Tdomain%sEdge(ne)%Renum(j)
                    if (Tdomain%sEdge(ne)%solid) then
                        if (Tdomain%sEdge(ne)%PML) then
                            Comm_data%Data(n)%IGiveSPML(nsolpml) = idx
                            nsolpml = nsolpml + 1
                        else
                            Comm_data%Data(n)%IGiveS(nsol) = idx
                            nsol = nsol + 1
                        endif
                    else ! Fluid
                        if (Tdomain%sEdge(ne)%PML) then
                            Comm_data%Data(n)%IGiveFPML(nflupml) = idx
                            nflupml = nflupml + 1
                        else
                            Comm_data%Data(n)%IGiveF(nflu) = idx
                            nflu = nflu + 1
                        endif
                    endif
                enddo
            enddo

            ! Vertices
            do i = 0,Tdomain%sComm(ncomm)%nb_vertices-1
                nv =  Tdomain%sComm(ncomm)%vertices(i)
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
            do i = 0,Tdomain%sComm(ncomm)%nb_faces-1
                nf = Tdomain%sComm(ncomm)%faces(i)
                ngll1 = Tdomain%sFace(nf)%ngll1
                ngll2 = Tdomain%sFace(nf)%ngll2

                if ( Tdomain%sComm(ncomm)%orient_faces(i) == 0 ) then
                    do j = 1,ngll2-2
                        do k = 1,ngll1-2
                            if (Tdomain%sFace(nf)%solid) then ! solide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(k,j)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(k,j)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                endif
                            else ! fluide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(k,j)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(k,j)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                endif
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 1 ) then
                    do j = 1,ngll2-2
                        do k = 1,ngll1-2
                            if (Tdomain%sFace(nf)%solid) then ! solide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,j)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,j)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                endif
                            else ! fluide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,j)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,j)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                endif
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 2 ) then
                    do j = 1,ngll2-2
                        do k = 1,ngll1-2
                            if (Tdomain%sFace(nf)%solid) then ! solide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                endif
                            else ! fluide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                endif
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 3 ) then
                    do j = 1,ngll2-2
                        do k = 1,ngll1-2
                            if (Tdomain%sFace(nf)%solid) then ! solide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                endif
                            else ! fluide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-k,ngll2-1-j)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                endif
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 4 ) then
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            if (Tdomain%sFace(nf)%solid) then ! solide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(j,k)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(j,k)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                endif
                            else ! fluide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(j,k)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(j,k)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                endif
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 5 ) then
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            if (Tdomain%sFace(nf)%solid) then ! solide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,k)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,k)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                endif
                            else ! fluide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,k)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,k)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                endif
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 6 ) then
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            if (Tdomain%sFace(nf)%solid) then ! solide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                endif
                            else ! fluide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                endif
                            endif
                        enddo
                    enddo
                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 7 ) then
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            if (Tdomain%sFace(nf)%solid) then ! solide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                    nsolpml = nsolpml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeS(nsol) = idx
                                    nsol = nsol + 1
                                endif
                            else ! fluide
                                if (Tdomain%sFace(nf)%PML) then
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                    nflupml = nflupml + 1
                                else
                                    idx = Tdomain%sFace(nf)%Renum(ngll1-1-j,ngll2-1-k)
                                    Comm_data%Data(n)%ITakeF(nflu) = idx
                                    nflu = nflu + 1
                                endif
                            endif
                        enddo
                    enddo
                endif ! orient_f
            enddo

            ! Edges
            do i = 0,Tdomain%sComm(ncomm)%nb_edges-1
                ne = Tdomain%sComm(ncomm)%edges(i)
                ngll1 = Tdomain%sEdge(ne)%ngll

                if ( Tdomain%sComm(ncomm)%orient_edges(i) == 0 ) then
                    do j = 1,ngll1-2
                        if (Tdomain%sEdge(ne)%solid) then
                            if (Tdomain%sEdge(ne)%PML) then
                                idx = Tdomain%sEdge(ne)%Renum(j)
                                Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                nsolpml = nsolpml + 1
                            else
                                idx = Tdomain%sEdge(ne)%Renum(j)
                                Comm_data%Data(n)%ITakeS(nsol) = idx
                                nsol = nsol + 1
                            endif
                        else
                            if (Tdomain%sEdge(ne)%PML) then
                                idx = Tdomain%sEdge(ne)%Renum(j)
                                Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                nflupml = nflupml + 1
                            else
                                idx = Tdomain%sEdge(ne)%Renum(j)
                                Comm_data%Data(n)%ITakeF(nflu) = idx
                                nflu = nflu + 1
                            endif
                        endif
                    enddo
                else if ( Tdomain%sComm(ncomm)%orient_edges(i) == 1 ) then
                    do j = 1,ngll1-2
                        if (Tdomain%sEdge(ne)%solid) then
                            if (Tdomain%sEdge(ne)%PML) then
                                idx = Tdomain%sEdge(ne)%Renum(ngll1-1-j)
                                Comm_data%Data(n)%ITakeSPML(nsolpml) = idx
                                nsolpml = nsolpml + 1
                            else
                                idx = Tdomain%sEdge(ne)%Renum(ngll1-1-j)
                                Comm_data%Data(n)%ITakeS(nsol) = idx
                                nsol = nsol + 1
                            endif
                        else
                            if (Tdomain%sEdge(ne)%PML) then
                                idx = Tdomain%sEdge(ne)%Renum(ngll1-1-j)
                                Comm_data%Data(n)%ITakeFPML(nflupml) = idx
                                nflupml = nflupml + 1
                            else
                                idx = Tdomain%sEdge(ne)%Renum(ngll1-1-j)
                                Comm_data%Data(n)%ITakeF(nflu) = idx
                                nflu = nflu + 1
                            endif
                        endif
                    enddo
                endif ! orient_edges
            enddo

            ! Vertices
            do i = 0,Tdomain%sComm(ncomm)%nb_vertices-1
                nv =  Tdomain%sComm(ncomm)%vertices(i)
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

end subroutine prepare_comm_vector


subroutine allocate_comm_vector(Tdomain,comm_data)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain
    type(comm_vector), intent(inout) :: comm_data
    integer :: n_data, n_comm, nsol, nsolpml, nflu, nflupml
    integer :: n, nf, ne, nv, i, temp

    n_comm = 0
    do n = 0,Tdomain%tot_comm_proc-1
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
    do n = 0,Tdomain%tot_comm_proc-1
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

        ! the size of data items (nddlxxx) is for communication from comm_forces
        ! the amount of data exchanged from define_arrays is different but lower for now:
        ! eg: DumpMass and MassMatSolPml acount for 6 real, compared to 9 for forcesPml
        ! for fluid we need only 4 during computation but 6 for mass exchange (but only
        ! to simplify code since 2 are really needed)
        n_data = 3*nsol+9*nsolpml+1*nflu+6*nflupml
        ! Initialisation et allocation de Comm_vector_DumpMassAndMMSP
        Comm_data%Data(n_comm)%src = Tdomain%rank
        Comm_data%Data(n_comm)%dest = Tdomain%sComm(n)%dest
        Comm_data%Data(n_comm)%ncomm = n
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

subroutine prepare_comm_vector_SF(Tdomain,ngll,renumSF,comm_data)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain
    type(comm_vector), intent(inout) :: comm_data
    integer, intent(in) :: ngll
    integer, intent(in), dimension(0:ngll-1,0:2) :: renumSF

    integer :: n,ncomm,nSF
    integer :: i,j,k,nf,ne,nv,idx,ngll1,ngll2, nb

    ! Remplissage des IGive et ITake
    if(Tdomain%nb_procs < 2)then
        comm_data%ncomm = 0
        return
    endif

    call allocate_comm_vector_SF(Tdomain, comm_data)

        do n = 0,Comm_data%ncomm-1
            ncomm = Comm_data%Data(n)%ncomm
            nSF = 0
            nb = 0

            ! Remplissage des Igive
            ! Faces
            do i = 0,Tdomain%sComm(ncomm)%SF_nf_shared-1
                nf = Tdomain%sComm(ncomm)%SF_faces_shared(i)
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Comm_data%Data(n)%IGiveS(nSF) = Tdomain%SF%SF_Face(i)%IG(k,j,2) ! on donne le solide
                        nSF = nSF + 1
                    end do
                end do
            enddo

            ! Edges
            do i = 0,Tdomain%sComm(ncomm)%SF_ne_shared-1
                ne = Tdomain%sComm(ncomm)%SF_edges_shared(i)
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    Comm_data%Data(n)%IGiveS(nSF) = Tdomain%SF%SF_Edge(i)%IG(j,2) ! on donne le solide
                    nSF = nSF + 1
                enddo
            enddo

            ! Vertices
            do i = 0,Tdomain%sComm(ncomm)%SF_nv_shared-1
                Comm_data%Data(n)%IGiveS(nSF) = Tdomain%SF%SF_Vertex(i)%IG(2) ! on donne le solide
                nSF = nSF + 1
            enddo

            write(*,*) "GGGGG rank,n,nb", Tdomain%rank, n, nb

            nSF = 0
            ! Remplissage des ITake
            ! Faces
            do i = 0,Tdomain%sComm(ncomm)%SF_nf_shared-1
                nf = Tdomain%sComm(ncomm)%SF_faces_shared(i)
                ngll1 = Tdomain%sFace(nf)%ngll1
                ngll2 = Tdomain%sFace(nf)%ngll2

                if ( Tdomain%sComm(ncomm)%orient_faces(i) == 0 ) then
                    do j = 1,ngll2-2
                        do k = 1,ngll1-2
                            idx = Tdomain%sFace(nf)%Iglobnum_Face(k,j)
                            Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                            nSF = nSF + 1
                        enddo
                    enddo

                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 1 ) then
                    do j = 1,ngll2-2
                        do k = 1,ngll1-2
                            idx = Tdomain%sFace(nf)%Iglobnum_Face(ngll1-1-k,j)
                            Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                            nSF = nSF + 1
                        enddo
                    enddo

                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 2 ) then
                   do j = 1,ngll2-2
                        do k = 1,ngll1-2
                            idx = Tdomain%sFace(nf)%Iglobnum_Face(k,ngll2-1-j)
                            Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                            nSF = nSF + 1
                        enddo
                    enddo

                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 3 ) then
                    do j = 1,ngll2-2
                        do k = 1,ngll1-2
                            idx = Tdomain%sFace(nf)%Iglobnum_Face(ngll1-1-k,ngll2-1-j)
                            Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                            nSF = nSF + 1
                        enddo
                    enddo

                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 4 ) then
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            idx = Tdomain%sFace(nf)%Iglobnum_Face(j,k)
                            Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                            nSF = nSF + 1
                        enddo
                    enddo

                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 5 ) then
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            idx = Tdomain%sFace(nf)%Iglobnum_Face(ngll1-1-j,k)
                            Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                            nSF = nSF + 1
                        enddo
                    enddo

                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 6 ) then
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            idx = Tdomain%sFace(nf)%Iglobnum_Face(j,ngll2-1-k)
                            Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                            nSF = nSF + 1
                        enddo
                    enddo

                else if ( Tdomain%sComm(ncomm)%orient_faces(i) == 7 ) then
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            idx = Tdomain%sFace(nf)%Iglobnum_Face(ngll1-1-j,ngll2-1-k)
                            Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                            nSF = nSF + 1
                        enddo
                    enddo
                endif ! orient_f
            enddo

            ! Edges
            do i = 0,Tdomain%sComm(ncomm)%SF_ne_shared-1
                ne = Tdomain%sComm(ncomm)%SF_edges_shared(i)
                ngll1 = Tdomain%sEdge(ne)%ngll

                if ( Tdomain%sComm(ncomm)%orient_edges(i) == 0 ) then
                    do j = 1,ngll1-2
                        idx = Tdomain%sEdge(ne)%Iglobnum_Edge(j)
                        Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                        nSF = nSF + 1
                    enddo

                else if ( Tdomain%sComm(ncomm)%orient_edges(i) == 1 ) then
                    do j = 1,ngll1-2
                        idx = Tdomain%sEdge(ne)%Iglobnum_Edge(ngll1-1-j)
                        Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                        nSF = nSF + 1
                    enddo
                endif ! orient_edges
            enddo

            ! Vertices
            do i = 0,Tdomain%sComm(ncomm)%SF_nv_shared-1
                nv =  Tdomain%sComm(ncomm)%SF_vertices_shared(i)
                idx = Tdomain%sVertex(nv)%Iglobnum_Vertex
                Comm_data%Data(n)%ITakeS(nSF) = renumSF(idx,0)
                nSF = nSF + 1
            enddo
        enddo

end subroutine prepare_comm_vector_SF

subroutine allocate_comm_vector_SF(Tdomain, comm_data)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain
    type(comm_vector), intent(inout) :: comm_data
    integer :: n_data, n_comm, nSF
    integer :: n, nf, ne, nv, i

    n_comm = 0
    do n = 0,Tdomain%tot_comm_proc-1
        if (Tdomain%sComm(n)%SF_nf_shared > 0 .OR. &
            Tdomain%sComm(n)%SF_ne_shared > 0 .OR. &
            Tdomain%sComm(n)%SF_nv_shared > 0) then
            n_comm = n_comm + 1
        endif
    enddo

    allocate(Comm_data%Data(0:n_comm-1))
    allocate(Comm_data%send_reqs(0:n_comm-1))
    allocate(Comm_data%recv_reqs(0:n_comm-1))
    Comm_data%ncomm = n_comm

    n_comm = 0
    do n = 0,Tdomain%tot_comm_proc-1
        if (Tdomain%sComm(n)%SF_nf_shared < 1 .AND. &
            Tdomain%sComm(n)%SF_ne_shared < 1 .AND. &
            Tdomain%sComm(n)%SF_nv_shared < 1) cycle

        nSF = 0

        ! Faces
        do i = 0,Tdomain%sComm(n)%SF_nf_shared-1
            nf = Tdomain%sComm(n)%SF_faces_shared(i)
            nSF = nSF + ((Tdomain%sFace(nf)%ngll1-2) * (Tdomain%sFace(nf)%ngll2-2))
        enddo
        ! Edges
        do i = 0,Tdomain%sComm(n)%SF_ne_shared-1
            ne = Tdomain%sComm(n)%SF_edges_shared(i)
            nSF = nSF + Tdomain%sEdge(ne)%ngll-2
        enddo
        ! Vertices
        do i = 0,Tdomain%sComm(n)%SF_nv_shared-1
            nv = Tdomain%sComm(n)%SF_vertices_shared(i)
            nSF = nSF + 1
        enddo

        n_data = nSF

        ! Initialisation et allocation de Comm_vector
        Comm_data%Data(n_comm)%src = Tdomain%rank
        Comm_data%Data(n_comm)%dest = Tdomain%sComm(n)%dest
        Comm_data%Data(n_comm)%ncomm = n
        Comm_data%Data(n_comm)%ndata = n_data
        allocate(Comm_data%Data(n_comm)%Give(0:n_data-1))
        allocate(Comm_data%Data(n_comm)%Take(0:n_data-1))
        allocate(Comm_data%Data(n_comm)%IGiveS(0:n_data-1))
        allocate(Comm_data%Data(n_comm)%ITakeS(0:n_data-1))

        n_comm = n_comm + 1
    enddo

    return
end subroutine allocate_comm_vector_SF

subroutine populate_index_SF(sf, dir,ngllx,nglly,ngllz, &
                             ISolFlu, IGlobnum,k,nglobpoints,renumSF, Isf)
    implicit none

    integer, intent(in) :: sf, dir, ngllx, nglly, ngllz, nglobpoints
    integer, intent(in), dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: ISolFlu
    integer, intent(in), dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: IGlobnum
    integer, intent(inout) :: k
    integer, intent(inout), dimension(0:nglobpoints-1, 0:2) :: renumSF
    integer, intent(inout), dimension(:,:), allocatable :: Isf

    integer :: i,j

    if(dir == 0)then
        do j = 0,nglly-1
            do i = 0,ngllx-1
                renumSF(k, 0) = Iglobnum(i,j,0)
                renumSF(k, sf) = ISolFlu(i,j,0)
                if (Isf(i,j)/=-1 .and. Isf(i,j)/=k) then
                    write(*,*) "Warning :...."
                else
                    Isf(i,j) = k
                endif
                k = k + 1
            enddo
        enddo
    else if(dir == 1)then
        do j = 0,ngllz-1
            do i = 0,ngllx-1
                renumSF(k, 0) = Iglobnum(i,0,j)
                renumSF(k, sf) = ISolFlu(i,0,j)
                if (Isf(i,j)/=-1 .and. Isf(i,j)/=k) then
                    write(*,*) "Warning :...."
                else
                    Isf(i,j) = k
                endif
                k = k + 1
            enddo
        enddo
    else if(dir == 2)then
        do j = 0,ngllz-1
            do i = 0,nglly-1
                renumSF(k, 0) = Iglobnum(ngllx-1,i,j)
                renumSF(k, sf) = ISolFlu(ngllx-1,i,j)
                if (Isf(i,j)/=-1 .and. Isf(i,j)/=k) then
                    write(*,*) "Warning :...."
                else
                    Isf(i,j) = k
                endif
                k = k + 1
            enddo
        enddo
    else if(dir == 3)then
        do j = 0,ngllz-1
            do i = 0,ngllx-1
                renumSF(k, 0) = Iglobnum(i,nglly-1,j)
                renumSF(k, sf) = ISolFlu(i,nglly-1,j)
                if (Isf(i,j)/=-1 .and. Isf(i,j)/=k) then
                    write(*,*) "Warning :...."
                else
                    Isf(i,j) = k
                endif
                k = k + 1
            enddo
        enddo
    else if(dir == 4)then
        do j = 0,ngllz-1
            do i = 0,nglly-1
                renumSF(k, 0) = Iglobnum(0,i,j)
                renumSF(k, sf) = ISolFlu(0,i,j)
                if (Isf(i,j)/=-1 .and. Isf(i,j)/=k) then
                    write(*,*) "Warning :...."
                else
                    Isf(i,j) = k
                endif
                k = k + 1
            enddo
        enddo
    else if(dir == 5)then
        do j = 0,nglly-1
            do i = 0,ngllx-1
                renumSF(k, 0) = Iglobnum(i,j,ngllz-1)
                renumSF(k, sf) = ISolFlu(i,j,ngllz-1)
                if (Isf(i,j)/=-1 .and. Isf(i,j)/=k) then
                    write(*,*) "Warning :...."
                else
                    Isf(i,j) = k
                endif
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
    integer, dimension(0:), intent (inout) :: renum
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
    integer, dimension(0:), intent (inout) :: renum
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
    integer, dimension(0:), intent (inout) :: renum
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
    integer, dimension(0:), intent (inout) :: renum
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

subroutine debug_comm_vector(Tdomain, src, dest, commvec)
    use sdomain
    implicit none
    type(domain), intent(in) :: Tdomain
    integer, intent(in) :: src, dest
    type(comm_vector), intent(in) :: commvec
    !
    integer :: i,k,rank

    rank = Tdomain%rank

    do i=0, commvec%ncomm-1
        if (commvec%Data(i)%src/=src .or. commvec%Data(i)%dest/=dest) cycle

        write(*,*) rank, "COMM:", src, dest
        write(*,*) rank, "NSOL ", commvec%Data(i)%nsol
        write(*,*) rank, "NSPML", commvec%Data(i)%nsolpml
        write(*,*) rank, "NFLU ", commvec%Data(i)%nflu
        write(*,*) rank, "NFPML", commvec%Data(i)%nflupml

        write(*,*) rank, "ISOL>", (commvec%Data(i)%IGiveS(k), k=0,10)
        write(*,*) rank, "ISOL<", (commvec%Data(i)%ITakeS(k), k=0,10)
    end do
end subroutine debug_comm_vector

integer function num_gll_face(Tdomain, nnf)
    use sdomain
    implicit none
    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: nnf
    !
    integer :: ngllx, nglly, ngllz, dir, n
    !
    if(nnf < 0) then
        num_gll_face = 0
    else
        n = Tdomain%sFace(nnf)%which_elem
        dir = Tdomain%sFace(nnf)%dir
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        if(dir == 0 .OR. dir == 5)then
            num_gll_face = (ngllx * nglly)
        else if(dir == 1 .OR. dir == 3)then
            num_gll_face = (ngllx * ngllz)
        else
            num_gll_face = (nglly * ngllz)
        endif
    end if
    return
end function num_gll_face

subroutine renumSolFlu(Tdomain)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain

    integer :: i,j,k,nf,ne,nv,ngll,ngll1,ngll2,nff,nfs,nne,nnv
    integer :: orient_e, orient_f

    k = 0
    do nv = 0,Tdomain%SF%SF_n_vertices-1
        Tdomain%SF%SF_Vertex(nv)%IG(0) = k
        k = k + 1
    end do

    do ne = 0,Tdomain%SF%SF_n_edges-1
        ngll = Tdomain%SF%SF_Edge(ne)%ngll
        allocate(Tdomain%SF%SF_Edge(ne)%IG(1:ngll-2,0:2))
        Tdomain%SF%SF_Edge(ne)%IG = -1
        do i=1,ngll-2
            Tdomain%SF%SF_Edge(ne)%IG(i,0) = k
            k = k +1
        end do
    end do

    do nf = 0,Tdomain%SF%SF_n_faces-1
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        ! gll to Index for : 0 sf_index, 1 : Iglobnum gll Fluide, 2: Iglobnum gll Solid
        allocate(Tdomain%SF%SF_Face(nf)%IG(0:ngll1-1,0:ngll2-1,0:2))
        Tdomain%SF%SF_Face(nf)%IG = -1
        do i=1,ngll1-2
            do j=1,ngll2-2
                Tdomain%SF%SF_Face(nf)%IG(i,j,0) = k
                k = k +1
            end do
        end do
        do i = 0,3
            ne = Tdomain%SF%SF_Face(nf)%Near_Edges(i)
            ! The fluid edge exists on this proc, Since the fluid element exists
            ngll = Tdomain%SF%SF_Edge(ne)%ngll
            orient_e = Tdomain%SF%SF_face(nf)%Near_Edges_Orient(i)
            call Number_Edge2Face(i,ngll1,ngll2,ngll,orient_e,         &
                Tdomain%SF%SF_face(nf)%IG(:,:,0),Tdomain%SF%SF_Edge(ne)%IG)
            nv = Tdomain%SF%SF_Face(nf)%Near_Vertices(i)
            call Number_Vertex2Face(i,ngll1,ngll2,    &
                Tdomain%SF%SF_face(nf)%IG(:,:,0),Tdomain%SF%SF_Vertex(nv)%IG)
        end do
    end do
    !! Here we have Tdomain%SF%SF_Face(nf)%IG(:,:,0) that contains unique shared numbering of each SF nodes

    Tdomain%SF%ngll = k

    !! Fluide
    do nf = 0,Tdomain%SF%SF_n_faces-1
        nff = Tdomain%SF%SF_Face(nf)%Face(0)
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        if (nff>=0) then
            orient_f = 0
            call Number_Face2Face(ngll1,ngll2, orient_f, &
                                  Tdomain%SF%SF_face(nf)%IG(:,:,1),Tdomain%sFace(nff)%Iglobnum_Face)
        end if

        do i = 0,3
            ne = Tdomain%SF%SF_Face(nf)%Near_Edges(i)
            nne = Tdomain%SF%SF_Edge(ne)%Edge(0)
            if (nne>=0) then
                ! The fluid edge exists on this proc, Since the fluid element exists
                ngll = Tdomain%SF%SF_Edge(ne)%ngll
                !!! XXX orient_e = ref car fluide present
                orient_e = Tdomain%SF%SF_face(nf)%Near_Edges_Orient(i)
                call Number_Edge2Face(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_face(nf)%IG(:,:,1),Tdomain%sEdge(nne)%Iglobnum_Edge)
            end if
            nv = Tdomain%SF%SF_Face(nf)%Near_Vertices(i)
            nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(0)
            if (nnv>=0) then
                call Number_Vertex2Face(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_face(nf)%IG(:,:,1),Tdomain%sVertex(nnv)%Iglobnum_Vertex)
            end if
        end do
    end do

    !! Solid
    do nf = 0,Tdomain%SF%SF_n_faces-1
        nfs = Tdomain%SF%SF_Face(nf)%Face(1)
        ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
        ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
        if (nfs>=0) then
            orient_f = Tdomain%SF%SF_Face(nf)%Orient_Face
            call Number_Face2Face(ngll1,ngll2, orient_f, &
                                  Tdomain%SF%SF_face(nf)%IG(:,:,2),Tdomain%sFace(nfs)%Iglobnum_Face)
            
        end if

        do i = 0,3
            ne = Tdomain%SF%SF_Face(nf)%Near_Edges(i)
            nne = Tdomain%SF%SF_Edge(ne)%Edge(1)
            if (nne>=0) then
                ! The solid edge exists on this proc, Since the solid element exists
                ngll = Tdomain%SF%SF_Edge(ne)%ngll
                !!! XXX si fluide present ok sinon orient_e = ref ? 
                orient_e = Tdomain%SF%SF_face(nf)%Near_Edges_Orient(i)
                call Number_Edge2Face(i,ngll1,ngll2,ngll,orient_e,         &
                    Tdomain%SF%SF_face(nf)%IG(:,:,2),Tdomain%sEdge(nne)%Iglobnum_Edge)
            endif
            nv = Tdomain%SF%SF_Face(nf)%Near_Vertices(i)
            nnv = Tdomain%SF%SF_Vertex(nv)%Vertex(1)
            if (nnv>=0) then
                call Number_Vertex2Face(i,ngll1,ngll2,    &
                    Tdomain%SF%SF_face(nf)%IG(:,:,2),Tdomain%sVertex(nnv)%Iglobnum_Vertex)
            endif
        end do
    end do

    return
end subroutine renumSolFlu

end module mrenumber
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
