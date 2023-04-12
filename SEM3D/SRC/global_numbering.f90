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
  implicit none
  public :: global_numbering, get_surface_numbering
  private :: allocate_comm_vector, prepare_comm_vector


contains

subroutine global_numbering(Tdomain)

    ! routine different from the 2D case. Everything is independently numbered, here (inner
    !      points in elements, on faces, edges and vertices). And then associated in
    !      the "Iglobnum" field = global number of each GLL point.


    use sdomain
    implicit none
    type(domain), intent (inout) :: Tdomain
    !
    integer :: k

    ! Create a unique gll number per node location
    call renumber_global_gll_nodes(Tdomain)

    ! Reorder domain nodes to improve locality in memory (avoid cache misses)
    ! buggy
    !call reorder_domains(Tdomain)
    call renumber_interface(Tdomain, Tdomain%intSolPml, DM_SOLID_CG, DM_SOLID_CG_PML)
    call renumber_interface(Tdomain, Tdomain%intFluPml, DM_FLUID_CG, DM_FLUID_CG_PML)
    if (Tdomain%logicD%SF_local_present) then
        call renumber_interface(Tdomain, Tdomain%SF%intSolFlu, DM_SOLID_CG, DM_FLUID_CG)
        call renumber_interface(Tdomain, Tdomain%SF%intSolFluPml, DM_SOLID_CG_PML, DM_FLUID_CG_PML)
    end if
    do k=0,size(Tdomain%sSurfaces)-1
        call renumber_surface(Tdomain, Tdomain%sSurfaces(k)%surf_sl, DM_SOLID_CG)
        call renumber_surface(Tdomain, Tdomain%sSurfaces(k)%surf_sldg, DM_SOLID_DG)
        call renumber_surface(Tdomain, Tdomain%sSurfaces(k)%surf_fl, DM_FLUID_CG)
        call renumber_surface(Tdomain, Tdomain%sSurfaces(k)%surf_spml, DM_SOLID_CG_PML)
        call renumber_surface(Tdomain, Tdomain%sSurfaces(k)%surf_fpml, DM_FLUID_CG_PML)
    end do

    call prepare_comm_vector(Tdomain, Tdomain%Comm_data)

    call prepare_comm_surface(Tdomain, Tdomain%Comm_SolFlu)

end subroutine global_numbering

subroutine renumber_global_gll_nodes(Tdomain)
    use sdomain
    use mindex
    implicit none
    type(domain), intent (inout) :: Tdomain
    integer :: n, i, j, k
    integer :: nf, nnf
    integer :: ne, nne
    integer :: nv, nnv
    ! Counts PML glls with abs flag in solid or fluid
    integer :: solid_abs_count, fluid_abs_count
    ! A counter for each domain (0: global)
    integer, dimension(0:DM_MAX) :: icount ! nombre de glls    par domaine
    integer, dimension(0:DM_MAX) :: ecount ! nombre d'elements par domaine
    integer, dimension(0:DM_MAX) :: fcount ! nombre de faces par domaine
    integer :: ngll, n3, n0, bnum, ee, lnum, vchunk0
    integer, dimension(0:3) :: elface
    integer, dimension(0:1) :: eledge
    integer :: dom
    integer :: idxi, idxj, idxk, ind
    integer, dimension(0:2) :: i0, di, dj
    logical :: fail

    icount = 0
    ecount = 0
    fcount = 0
    solid_abs_count = 0
    fluid_abs_count = 0

    !Elements Inner GLL points
    do n = 0,Tdomain%n_elem-1
        ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
        dom = Tdomain%specel(n)%domain
        Tdomain%specel(n)%lnum = ecount(dom)
        ecount(dom) = ecount(dom)+1
        allocate(Tdomain%specel(n)%Iglobnum(0:ngll-1,0:ngll-1,0:ngll-1))
        allocate(Tdomain%specel(n)%Idom    (0:ngll-1,0:ngll-1,0:ngll-1))
        Tdomain%specel(n)%Iglobnum = -1
        Tdomain%specel(n)%Idom = -1
        ! assign a global number first
        do k = 1,ngll-2
            do j = 1,ngll-2
                do i = 1,ngll-2
                    Tdomain%specel(n)%Iglobnum(i,j,k) = icount(0)
                    icount(0) = icount(0) + 1
!                    Tdomain%specel(n)%Idom(i,j,k) = icount(dom)
!                    icount(dom) = icount(dom) + 1
                enddo
            enddo
        enddo
    end do
    do n = 0,Tdomain%n_elem-1
        ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
        dom = Tdomain%specel(n)%domain
        ! assign local numbering for inside gll nodes
        n3 = (ngll-2)*(ngll-2)*(ngll-2)
        lnum = Tdomain%specel(n)%lnum
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        n0 = bnum*VCHUNK*n3+ee
        vchunk0 = VCHUNK
        ! make sure we have no 'holes' in the numbering
        if ((bnum+1)*VCHUNK>ecount(dom)) then
            ! last block of elems may not be complete
            vchunk0 = mod(ecount(dom),VCHUNK)
        endif
        do k = 1,ngll-2
            do j = 1,ngll-2
                do i = 1,ngll-2
                    Tdomain%specel(n)%Idom(i,j,k) = n0
                    n0 = n0 + vchunk0
                enddo
            enddo
        enddo
        icount(dom) = icount(dom) + n3
    enddo
    Tdomain%sdom   %nbelem = ecount(DM_SOLID_CG)
    Tdomain%sdomdg %nbelem = ecount(DM_SOLID_DG)
    Tdomain%fdom   %nbelem = ecount(DM_FLUID_CG)
    Tdomain%spmldom%nbelem = ecount(DM_SOLID_CG_PML)
    Tdomain%fpmldom%nbelem = ecount(DM_FLUID_CG_PML)

    !Faces Inner GLL points
    do n = 0,Tdomain%n_face-1
        dom   = Tdomain%sFace(n)%domain
        ngll = domain_ngll(Tdomain,dom)
        Tdomain%sFace(n)%lnum = fcount(dom)
        fcount(dom) = fcount(dom)+1
        allocate(Tdomain%sFace(n)%Iglobnum_Face(0:ngll-1,0:ngll-1))
        allocate(Tdomain%sFace(n)%Idom         (0:ngll-1,0:ngll-1))
        Tdomain%sFace(n)%Iglobnum_Face = -1
        Tdomain%sFace(n)%Idom = -1
        do j = 1,ngll-2
            do i = 1,ngll-2
                Tdomain%sFace(n)%Iglobnum_Face(i,j) = icount(0)
                icount(0) = icount(0) + 1
                Tdomain%sFace(n)%Idom(i,j) = icount(dom)
                icount(dom) = icount(dom) + 1
            enddo
        enddo
    enddo
    Tdomain%sdom   %nbface = fcount(DM_SOLID_CG)
    Tdomain%sdomdg %nbface = fcount(DM_SOLID_DG)
    Tdomain%fdom   %nbface = fcount(DM_FLUID_CG)
    Tdomain%spmldom%nbface = fcount(DM_SOLID_CG_PML)
    Tdomain%fpmldom%nbface = fcount(DM_FLUID_CG_PML)

    !Edges Inner GLL points
    do n = 0,Tdomain%n_edge-1
        dom   = Tdomain%sEdge(n)%domain
        ngll = domain_ngll(Tdomain, dom)
        allocate(Tdomain%sEdge(n)%Iglobnum_Edge(1:ngll-2))
        allocate(Tdomain%sEdge(n)%Idom         (1:ngll-2))
        Tdomain%sEdge(n)%Iglobnum_Edge = -1
        Tdomain%sEdge(n)%Idom = -1
        do i = 1,ngll-2
            Tdomain%sEdge(n)%Iglobnum_Edge(i) = icount(0)
            icount(0) = icount(0) + 1
            Tdomain%sEdge(n)%Idom(i) = icount(dom)
            icount(dom) = icount(dom) + 1
        enddo
    enddo

    !Corner GLL points
    do n = 0,Tdomain%n_vertex-1
        dom = Tdomain%sVertex(n)%domain
        Tdomain%sVertex(n)%Iglobnum_Vertex = icount(0)
        icount(0) = icount(0) + 1
        Tdomain%sVertex(n)%Idom = icount(dom)
        icount(dom) = icount(dom) + 1
    enddo

    ! total number of GLL points (= degrees of freedom)
    Tdomain%n_glob_points   = icount(0)
    Tdomain%sdom%nglltot    = icount(DM_SOLID_CG)
    Tdomain%sdomdg%nglltot  = icount(DM_SOLID_DG)
    Tdomain%fdom%nglltot    = icount(DM_FLUID_CG)
    Tdomain%spmldom%nglltot = icount(DM_SOLID_CG_PML)
    Tdomain%fpmldom%nglltot = icount(DM_FLUID_CG_PML)

    !Recollecting at the element level, from faces, edges and vertices.
    do n = 0,Tdomain%n_elem-1
        dom = Tdomain%specel(n)%domain
        ngll = domain_ngll(Tdomain, dom)
        !Taking information from faces
        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            do k=0,3
                elface(k) = Tdomain%specel(n)%Control_nodes(face_def(k,nf))
            end do
            call ind_elem_face(ngll, nf, Tdomain%sFace(nnf)%inodes, elface, i0, di, dj)

            do i=1,ngll-2
                do j=1,ngll-2
                    idxi = i0(0)+i*di(0)+j*dj(0)
                    idxj = i0(1)+i*di(1)+j*dj(1)
                    idxk = i0(2)+i*di(2)+j*dj(2)
                    ind = Tdomain%sFace(nnf)%Iglobnum_Face(i,j)
                    Tdomain%specel(n)%Iglobnum(idxi,idxj,idxk) = ind
                    ind = Tdomain%sFace(nnf)%Idom(i,j)
                    Tdomain%specel(n)%Idom(idxi,idxj,idxk) = ind
                end do
            end do
        end do

        !Taking information from edges
        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            do k=0,1
                eledge(k) = Tdomain%specel(n)%Control_nodes(edge_def(k,ne))
            end do
            call ind_elem_edge(ngll, ne, Tdomain%sEdge(nne)%inodes, eledge, i0, di)
            !write(*,*) ne, i0, di, eledge
            do i=1,ngll-2
                idxi = i0(0)+i*di(0)
                idxj = i0(1)+i*di(1)
                idxk = i0(2)+i*di(2)
                Tdomain%specel(n)%Iglobnum(idxi,idxj,idxk) = Tdomain%sEdge(nne)%Iglobnum_Edge(i)
                Tdomain%specel(n)%Idom(idxi,idxj,idxk) = Tdomain%sEdge(nne)%Idom(i)
            end do
        end do

        !Taking information from vertices
        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            idxi = vertex_def(0,nv)*(ngll-1)
            idxj = vertex_def(1,nv)*(ngll-1)
            idxk = vertex_def(2,nv)*(ngll-1)

            Tdomain%specel(n)%Iglobnum(idxi, idxj, idxk) = Tdomain%sVertex(nnv)%Iglobnum_Vertex
            Tdomain%specel(n)%Idom(idxi, idxj, idxk) = Tdomain%sVertex(nnv)%Idom
        end do

    enddo    ! end of the loop onto elements

    ! One last step, copy the borders (edges, vertices) from the elements to the faces
    do n = 0,Tdomain%n_elem-1
        dom = Tdomain%specel(n)%domain
        ngll = domain_ngll(Tdomain, dom)
        !Taking information from faces
        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            do k=0,3
                elface(k) = Tdomain%specel(n)%Control_nodes(face_def(k,nf))
            end do
            call ind_elem_face(ngll, nf, Tdomain%sFace(nnf)%inodes, elface, i0, di, dj)

            do i=0,ngll-1
                j = 0
                idxi = i0(0)+i*di(0)+j*dj(0)
                idxj = i0(1)+i*di(1)+j*dj(1)
                idxk = i0(2)+i*di(2)+j*dj(2)
                ind = Tdomain%specel(n)%Iglobnum(idxi,idxj,idxk)
                Tdomain%sFace(nnf)%Iglobnum_Face(i,j) = ind
                ind = Tdomain%specel(n)%Idom(idxi,idxj,idxk)
                Tdomain%sFace(nnf)%Idom(i,j) = ind
                !
                j = ngll-1
                idxi = i0(0)+i*di(0)+j*dj(0)
                idxj = i0(1)+i*di(1)+j*dj(1)
                idxk = i0(2)+i*di(2)+j*dj(2)
                ind = Tdomain%specel(n)%Iglobnum(idxi,idxj,idxk)
                Tdomain%sFace(nnf)%Iglobnum_Face(i,j) = ind
                ind = Tdomain%specel(n)%Idom(idxi,idxj,idxk)
                Tdomain%sFace(nnf)%Idom(i,j) = ind
            end do
            do j=0,ngll-1
                i = 0
                idxi = i0(0)+i*di(0)+j*dj(0)
                idxj = i0(1)+i*di(1)+j*dj(1)
                idxk = i0(2)+i*di(2)+j*dj(2)
                ind = Tdomain%specel(n)%Iglobnum(idxi,idxj,idxk)
                Tdomain%sFace(nnf)%Iglobnum_Face(i,j) = ind
                ind = Tdomain%specel(n)%Idom(idxi,idxj,idxk)
                Tdomain%sFace(nnf)%Idom(i,j) = ind
                !
                i = ngll-1
                idxi = i0(0)+i*di(0)+j*dj(0)
                idxj = i0(1)+i*di(1)+j*dj(1)
                idxk = i0(2)+i*di(2)+j*dj(2)
                ind = Tdomain%specel(n)%Iglobnum(idxi,idxj,idxk)
                Tdomain%sFace(nnf)%Iglobnum_Face(i,j) = ind
                ind = Tdomain%specel(n)%Idom(idxi,idxj,idxk)
                Tdomain%sFace(nnf)%Idom(i,j) = ind
            end do
        end do
    end do
#if 1
    ! Check that all glls have an Iglobnum
    do n = 0,Tdomain%n_elem-1
        ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
        dom  = Tdomain%specel(n)%domain
        fail = .false.
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    ind = Tdomain%specel(n)%Iglobnum(i,j,k)
                    !write(*,*) i,j,k,ind
                    if (ind<0) fail=.true.
                    ind = Tdomain%specel(n)%Idom(i,j,k)
                    if (ind<0) fail=.true.
                enddo
            enddo
        enddo
        if (fail) stop 1
    enddo
#endif
end subroutine renumber_global_gll_nodes

! Compute the mapping between face%Idom(i,j) to surf%map so that
! map(renum(face%Idom(i,j))) == face%Idom(i,j) for each face of surf
subroutine get_surface_numbering(Tdomain, surf, dom, renum)
    use sdomain
    use sinterface
    implicit none
    type(domain), intent(inout) :: Tdomain
    type(surf_num), intent(in) :: surf
    integer, intent(in) :: dom
    integer, dimension(:), allocatable, intent(out) :: renum
    !
    integer :: nglltot, i

    nglltot = domain_nglltot(Tdomain, dom)
    allocate(renum(0:nglltot-1))
    renum = -1
    do i=0,surf%nbtot-1
        renum(surf%map(i)) = i
    end do
end subroutine get_surface_numbering

subroutine renumber_surface(Tdomain, surf, dom0)
    use sdomain
    use mindex
    implicit none
    type(domain), intent (inout) :: Tdomain
    type(surf_num), intent(inout) :: surf
    integer, intent(in) :: dom0
    !
    integer :: nf, nfs
    integer :: ne, nes
    integer :: nv, nvs
    integer :: ngll_if
    integer :: i, j
    integer :: ngll
    integer :: doms
    integer :: idx0
    ! Count number of gll on the interface
    ngll_if = 0
    ! The order in which we store the gll in map is important because it is reused
    ! elsewhere (ie computation of the face normals)
    do nf=0,surf%n_faces-1
        nfs = surf%if_faces(nf)
        doms = Tdomain%sFace(nfs)%domain
        ngll = domain_ngll(Tdomain, doms)
        ngll_if = ngll_if + (ngll-2)*(ngll-2)
        ! While we're at it, check coherency...
        if (doms/=dom0) then
            write(*,*) 'Incoherent Solid-PML interface', doms, dom0
            stop 1
        end if
    end do
    do ne=0,surf%n_edges-1
        nes = surf%if_edges(ne)
        doms = Tdomain%sEdge(nes)%domain
        ngll = domain_ngll(Tdomain, doms)
        ngll_if = ngll_if + ngll-2
        ! While we're at it, check coherency...
        if (doms/=dom0) then
            write(*,*) 'Incoherent Solid-PML interface', doms, dom0
            stop 1
        end if
    end do
    ! VERTICES
    do nv=0,surf%n_vertices-1
        nvs = surf%if_vertices(nv)
        doms = Tdomain%sVertex(nvs)%domain
        ngll_if = ngll_if + 1
        if (doms/=dom0) then
            write(*,*) 'Incoherent Solid-PML interface', doms, dom0
            stop 1
        end if
    end do
    surf%nbtot = ngll_if
    allocate(surf%map(0:ngll_if-1))
    ! copy gll numbers on the interface
    ngll_if = 0
    ! FACES
    do nf=0,surf%n_faces-1
        nfs = surf%if_faces(nf)
        doms = Tdomain%sFace(nfs)%domain
        ngll = domain_ngll(Tdomain, doms)
        do j=1,ngll-2
            do i=1,ngll-2
                surf%map(ngll_if) = Tdomain%sFace(nfs)%Idom(i,j)
                ngll_if = ngll_if + 1
            end do
        end do
    end do
    ! EDGES
    do ne=0,surf%n_edges-1
        nes = surf%if_edges(ne)
        doms = Tdomain%sEdge(nes)%domain
        ngll = domain_ngll(Tdomain, doms)
        do i=1,ngll-2
            idx0 = Tdomain%sEdge(nes)%Idom(i)
            surf%map(ngll_if) = idx0
            ngll_if = ngll_if + 1
        end do
    end do
    ! VERTICES
    do nv=0,surf%n_vertices-1
        nvs = surf%if_vertices(nv)
        surf%map(ngll_if) = Tdomain%sVertex(nvs)%Idom
        ngll_if = ngll_if + 1
    end do
    ! Check
    if (ngll_if/=surf%nbtot) then
        write(*,*) "Incoherent interface face+edge+vert != face"
        stop 1
    end if
end subroutine renumber_surface

subroutine renumber_interface(Tdomain, inter, dom0, dom1)
    use sdomain
    use mindex
    implicit none
    type(domain), intent (inout) :: Tdomain
    type(inter_num), intent(inout) :: inter
    integer, intent(in) :: dom0, dom1
    ! Count number of gll on the interface
    call renumber_surface(Tdomain, inter%surf0, dom0)
    call renumber_surface(Tdomain, inter%surf1, dom1)
    ! We need to make sure both faces of an interface are numbered.
    ! Orphan faces (ie whose elements belong to another processor
    ! will not have Idom or Iglobnum filled correctly at this point
    call apply_numbering_coherency(Tdomain, inter, dom0, dom1)

end subroutine renumber_interface

subroutine apply_numbering_coherency(Tdomain, inter, dom0, dom1)
    use sdomain
    use mindex
    implicit none
    type(domain), intent (inout) :: Tdomain
    type(inter_num), intent(inout) :: inter
    integer, intent(in) :: dom0, dom1
    !
    integer, dimension(:), allocatable :: renum0, renum1
    !
    integer :: nglltot0, nglltot1, i, j
    integer :: idx, imap, idom
    integer :: ngll, ngll1
    integer :: nf, nfs0, nfs1
    logical :: face0_orphan, face1_orphan
    nglltot0 = domain_nglltot(Tdomain, dom0)
    nglltot1 = domain_nglltot(Tdomain, dom1)
    ngll     = domain_ngll(Tdomain, dom0)
    ngll1    = domain_ngll(Tdomain, dom1)

    allocate(renum0(0:nglltot0-1))
    allocate(renum1(0:nglltot1-1))

    if (ngll1/=ngll .and. inter%surf0%nbtot>0) then
        write(*,*) "Error incompatible domains ", dom0, ngll, dom1, ngll1
        stop 1
    end if

    do i=0,inter%surf0%nbtot-1
        idx = inter%surf0%map(i)
        if (idx>=0) then
            renum0(idx) = i
        else
            write(*,*) "Error dom=", dom0, "map:", i
        end if
    end do
    do i=0,inter%surf1%nbtot-1
        idx = inter%surf1%map(i)
        if (idx>=0) then
            renum1(idx) = i
        else
            write(*,*) "Error dom=", dom1, "map:", i
        end if
    end do
    ! FACES
    do nf=0,inter%surf0%n_faces-1
        nfs0 = inter%surf0%if_faces(nf)
        nfs1 = inter%surf1%if_faces(nf)
        face0_orphan = .false.
        face1_orphan = .false.

        do j=0,ngll-1
            do i=0,ngll-1
                if (Tdomain%sFace(nfs0)%Idom(i,j)==-1) then
                    face0_orphan = .true.
                end if
                if (Tdomain%sFace(nfs1)%Idom(i,j)==-1) then
                    face1_orphan = .true.
                end if
            end do
        end do
        if (face0_orphan .and. face1_orphan) then
            write(*,*) "Coherency problem"
            stop 1
        end if

        if (face0_orphan) then
            do j=0,ngll-1
                do i=0,ngll-1
                    ! Unconditionnally copy Iglobnum from 'good' face
                    Tdomain%sFace(nfs0)%Iglobnum_Face(i,j) = Tdomain%sFace(nfs1)%Iglobnum_Face(i,j)
                    if (Tdomain%sFace(nfs0)%Idom(i,j)==-1) then
                        ! This
                        idx = Tdomain%sFace(nfs1)%Idom(i,j) ! can't be -1
                        imap = renum1(idx)
                        idom = inter%surf0%map(imap)
                        Tdomain%sFace(nfs0)%Idom(i,j) = idom
                    end if
                end do
            end do
        end if
        if (face1_orphan) then
            do j=0,ngll-1
                do i=0,ngll-1
                    ! Unconditionnally copy Iglobnum from 'good' face
                    Tdomain%sFace(nfs1)%Iglobnum_Face(i,j) = Tdomain%sFace(nfs0)%Iglobnum_Face(i,j)
                    if (Tdomain%sFace(nfs1)%Idom(i,j)==-1) then
                        ! This
                        idx = Tdomain%sFace(nfs0)%Idom(i,j) ! can't be -1
                        imap = renum0(idx)
                        idom = inter%surf1%map(imap)
                        Tdomain%sFace(nfs1)%Idom(i,j) = idom
                    end if
                end do
            end do
        end if
    end do

    deallocate(renum0, renum1)
end subroutine apply_numbering_coherency

subroutine prepare_comm_vector(Tdomain,comm_data)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain
    type(comm_vector), intent(inout) :: comm_data

    integer :: n,ncomm,nsol,nsolpml,nflu,nflupml,nsoldg
    integer :: i,j,k,nf,ne,nv,idx
    integer :: dom, ngll

    ! Remplissage des IGive et ITake
    if(Tdomain%nb_procs < 2)then
        Comm_data%ncomm = 0
        return
    endif

    call allocate_comm_vector(Tdomain,comm_data)

    do n = 0,Comm_data%ncomm-1
        ncomm = Comm_data%Data(n)%ncomm
        nsol = 0
        nsoldg = 0
        nsolpml = 0
        nflu = 0
        nflupml = 0

        ! Remplissage des Igive
        ! Faces
        do i = 0,Tdomain%sComm(ncomm)%nb_faces-1
            nf = Tdomain%sComm(ncomm)%faces(i)
            dom = Tdomain%sFace(nf)%domain
            ngll = domain_ngll(Tdomain, dom)
            do j = 1,ngll-2
                do k = 1,ngll-2
                    idx = Tdomain%sFace(nf)%Idom(k,j)
                    select case(dom)
                    case (DM_SOLID_CG)
                        Comm_data%Data(n)%IGiveS(nsol) = idx
                        nsol = nsol + 1
                    case (DM_SOLID_DG)
                        Comm_data%Data(n)%IGiveSDG(nsoldg) = idx
                        nsoldg = nsoldg + 1
                    case (DM_FLUID_CG)
                        Comm_data%Data(n)%IGiveF(nflu) = idx
                        nflu = nflu + 1
                    case (DM_SOLID_CG_PML)
                        Comm_data%Data(n)%IGiveSPML(nsolpml) = idx
                        nsolpml = nsolpml + 1
                    case (DM_FLUID_CG_PML)
                        Comm_data%Data(n)%IGiveFPML(nflupml) = idx
                        nflupml = nflupml + 1
                    case default
                        stop "unknown domain"
                    end select
                end do
            end do
        enddo

        ! Edges
        do i = 0,Tdomain%sComm(ncomm)%nb_edges-1
            ne = Tdomain%sComm(ncomm)%edges(i)
            dom = Tdomain%sEdge(ne)%domain
            ngll = domain_ngll(Tdomain, dom)
            do j = 1,ngll-2
                idx = Tdomain%sEdge(ne)%Idom(j)
                select case(dom)
                case (DM_SOLID_CG)
                    Comm_data%Data(n)%IGiveS(nsol) = idx
                    nsol = nsol + 1
                case (DM_SOLID_DG)
                    Comm_data%Data(n)%IGiveSDG(nsoldg) = idx
                    nsoldg = nsoldg + 1
                case (DM_FLUID_CG)
                    Comm_data%Data(n)%IGiveF(nflu) = idx
                    nflu = nflu + 1
                case (DM_SOLID_CG_PML)
                    Comm_data%Data(n)%IGiveSPML(nsolpml) = idx
                    nsolpml = nsolpml + 1
                case (DM_FLUID_CG_PML)
                    Comm_data%Data(n)%IGiveFPML(nflupml) = idx
                    nflupml = nflupml + 1
                case default
                    stop "unknown domain"
                end select
            enddo
        enddo

        ! Vertices
        do i = 0,Tdomain%sComm(ncomm)%nb_vertices-1
            nv =  Tdomain%sComm(ncomm)%vertices(i)
            idx = Tdomain%sVertex(nv)%Idom
            select case(Tdomain%sVertex(nv)%domain)
            case (DM_SOLID_CG)
                Comm_data%Data(n)%IGiveS(nsol) = idx
                nsol = nsol + 1
            case (DM_SOLID_DG)
                Comm_data%Data(n)%IGiveSDG(nsoldg) = idx
                nsoldg = nsoldg + 1
            case (DM_FLUID_CG)
                Comm_data%Data(n)%IGiveF(nflu) = idx
                nflu = nflu + 1
            case (DM_SOLID_CG_PML)
                Comm_data%Data(n)%IGiveSPML(nsolpml) = idx
                nsolpml = nsolpml + 1
            case (DM_FLUID_CG_PML)
                Comm_data%Data(n)%IGiveFPML(nflupml) = idx
                nflupml = nflupml + 1
            case default
                stop "unknown domain"
            end select
        enddo
    end do
end subroutine prepare_comm_vector


subroutine allocate_comm_vector(Tdomain,comm_data)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain
    type(comm_vector), intent(inout) :: comm_data
    integer :: n_data, n_comm, nsol, nsolpml, nflu, nflupml, nsoldg
    integer :: n, nf, ne, nv, i, temp
    integer :: dom, ngll

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
        nsoldg = 0
        nflupml = 0
        nflu = 0

        ! Faces
        do i = 0,Tdomain%sComm(n)%nb_faces-1
            nf = Tdomain%sComm(n)%faces(i)
            dom = Tdomain%sFace(nf)%domain
            ngll = domain_ngll(Tdomain, dom)
            temp = (ngll-2) * (ngll-2)
            select case (dom)
            case (DM_SOLID_CG_PML)
                nsolpml = nsolpml + temp
            case (DM_SOLID_CG)
                nsol = nsol + temp
            case (DM_SOLID_DG)
                nsoldg = nsoldg + temp
            case (DM_FLUID_CG_PML)
                nflupml = nflupml + temp
            case (DM_FLUID_CG)
                nflu = nflu + temp
            case default
                stop "unknown domain"
            end select
        enddo
        ! Edges
        do i = 0,Tdomain%sComm(n)%nb_edges-1
            ne = Tdomain%sComm(n)%edges(i)
            dom = Tdomain%sEdge(ne)%domain
            ngll = domain_ngll(Tdomain, dom)
            temp = ngll-2
            select case (dom)
            case (DM_SOLID_CG_PML)
                nsolpml = nsolpml + temp
            case (DM_SOLID_CG)
                nsol = nsol + temp
            case (DM_SOLID_DG)
                nsoldg = nsoldg + temp
            case (DM_FLUID_CG_PML)
                nflupml = nflupml + temp
            case (DM_FLUID_CG)
                nflu = nflu + temp
            case default
                stop "unknown domain"
            end select
        enddo
        ! Vertices
        do i = 0,Tdomain%sComm(n)%nb_vertices-1
            nv = Tdomain%sComm(n)%vertices(i)
            select case (Tdomain%sVertex(nv)%domain)
            case (DM_SOLID_CG_PML)
                nsolpml = nsolpml + 1
            case (DM_SOLID_CG)
                nsol = nsol + 1
            case (DM_SOLID_DG)
                nsoldg = nsoldg + 1
            case (DM_FLUID_CG_PML)
                nflupml = nflupml + 1
            case (DM_FLUID_CG)
                nflu = nflu + 1
            case default
                stop "unknown domain"
            end select
        enddo

        ! the size of data items (nddlxxx) is for communication from comm_forces
        ! the amount of data exchanged from define_arrays is different but lower for now:
        ! eg: DumpMass and MassMatSolPml acount for 6 real, compared to 9 for forcesPml
        ! for fluid we need only 4 during computation but 6 for mass exchange (but only
        ! to simplify code since 2 are really needed)
        n_data = 3*nsol+3*nsoldg + 9*nsolpml+1*nflu+6*nflupml
        ! Initialisation et allocation de Comm_vector_DumpMassAndMMSP
        Comm_data%Data(n_comm)%src = Tdomain%rank
        Comm_data%Data(n_comm)%dest = Tdomain%sComm(n)%dest
        Comm_data%Data(n_comm)%ncomm = n
        Comm_data%Data(n_comm)%ndata = n_data
        Comm_data%Data(n_comm)%nsol = nsol
        Comm_data%Data(n_comm)%nsoldg = nsoldg
        Comm_data%Data(n_comm)%nsolpml = nsolpml
        Comm_data%Data(n_comm)%nflu = nflu
        Comm_data%Data(n_comm)%nflupml = nflupml
        allocate(Comm_data%Data(n_comm)%Give(0:n_data-1))
        allocate(Comm_data%Data(n_comm)%Take(0:n_data-1))
        allocate(Comm_data%Data(n_comm)%IGiveS(0:nsol-1))
        allocate(Comm_data%Data(n_comm)%IGiveSDG(0:nsoldg-1))
        allocate(Comm_data%Data(n_comm)%IGiveSPML(0:nsolpml-1))
        allocate(Comm_data%Data(n_comm)%IGiveF(0:nflu-1))
        allocate(Comm_data%Data(n_comm)%IGiveFPML(0:nflupml-1))

!        write(*,*) "COMM:", Tdomain%rank, "->", Comm_Data%Data(n_comm)%dest, ": NGLLS ", nsol
!        write(*,*) "COMM:", Tdomain%rank, "->", Comm_Data%Data(n_comm)%dest, ": NGLLF ", nflu
!        write(*,*) "COMM:", Tdomain%rank, "->", Comm_Data%Data(n_comm)%dest, ": NGLLSP", nsolpml
!        write(*,*) "COMM:", Tdomain%rank, "->", Comm_Data%Data(n_comm)%dest, ": NGLLFP", nflupml
        n_comm = n_comm + 1
    enddo

    return
end subroutine allocate_comm_vector


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
!        write(*,*) rank, "ISOL<", (commvec%Data(i)%ITakeS(k), k=0,10)
    end do
end subroutine debug_comm_vector

subroutine reorder_domains(Tdomain)
    use sdomain
    implicit none
    type(domain), intent (inout) :: Tdomain
    !
    integer, dimension(:,:), allocatable :: reorder
    integer :: n, i, j, k
    integer :: ind, dom, new_ind
    integer :: ngll
    integer, dimension(0:5) :: oind
    ! Initialize new order : set new order to -1

    allocate(reorder(0:5,0:Tdomain%n_glob_points-1))
    reorder(:,:)   =-1
    oind(:) = 0 ! Ordered index

    ! Reorder over elements
    ngll = 0
    new_ind = 0
    do n = 0,Tdomain%n_elem-1
        ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    ind = Tdomain%specel(n)%Idom(i,j,k)
                    dom = Tdomain%specel(n)%domain
                    new_ind = reorder(dom, ind)
                    if (new_ind==-1) then
                        new_ind = oind(dom)
                        oind(dom) = oind(dom) + 1
                        reorder(dom, ind) = new_ind
                    end if
                    Tdomain%specel(n)%Idom(i,j,k) = new_ind
                end do
            end do
        end do
    end do

    ! This is still usefull for ghost faces/edges/vertices
    ! from other cpus
    ! Reorder over faces
    do n = 0, Tdomain%n_face-1
        dom = Tdomain%sFace(n)%domain
        ngll = domain_ngll(Tdomain, dom)
        do j = 0,ngll-1
            do i = 0,ngll-1
                ind = Tdomain%sFace(n)%Idom(i,j)
                new_ind = reorder(dom, ind)
                if (new_ind==-1) then
                    new_ind = oind(dom)
                    oind(dom) = oind(dom) + 1
                    reorder(dom, ind) = new_ind
                end if
                Tdomain%sFace(n)%Idom(i,j) = new_ind
            enddo
        enddo
    end do

    ! Reorder over edges
    do n = 0, Tdomain%n_edge-1
        dom = Tdomain%sEdge(n)%domain
        ngll = domain_ngll(Tdomain, dom)
        do i = 0,ngll-1
            ind = Tdomain%sEdge(n)%Idom(i)
            new_ind = reorder(dom, ind)
            if (new_ind==-1) then
                new_ind = oind(dom)
                oind(dom) = oind(dom) + 1
                reorder(dom, ind) = new_ind
            end if
            Tdomain%sEdge(n)%Idom(i) = new_ind
        enddo
    end do

    ! Reorder over vertices
    do n = 0, Tdomain%n_vertex-1
        ind = Tdomain%sVertex(n)%Idom
        dom = Tdomain%sVertex(n)%domain
        if (new_ind==-1) then
            new_ind = oind(dom)
            oind(dom) = oind(dom) + 1
            reorder(dom, ind) = new_ind
        end if
        Tdomain%sVertex(n)%Idom = new_ind
    end do

    deallocate(reorder)
end subroutine reorder_domains

subroutine build_comms_surface(Tdomain, comm_data, surface, dom)
    use sdomain
    use mintset
    implicit none
    type(domain), intent (inout) :: Tdomain
    type(comm_vector), intent(inout) :: comm_data
    type(surf_num), intent(in) :: surface
    integer, intent(in) :: dom
    !
    integer :: n, src, dst, n1, n2
    integer :: count, idx, i
    integer :: ngll
    integer, dimension(:), allocatable :: igive
    integer, allocatable, dimension(:) :: renum
    !
    src = Tdomain%rank
!    write(*,*) "Surface:", surface%n_vertices, surface%n_edges
    call get_surface_numbering(Tdomain, surface, dom, renum)

    ngll = domain_ngll(Tdomain, dom)
    do n = 0,Tdomain%Comm_data%ncomm-1
        dst = Tdomain%Comm_data%Data(n)%dest
        count = 0
        n2 = surface%nbtot
        select case(dom)
        case (DM_FLUID_CG)
            n1 = Tdomain%Comm_data%Data(n)%nflu
            if (n1>0.and.n2>0) then
                count = intersect_arrays(n1, Tdomain%Comm_data%Data(n)%IGiveF, n2, surface%map, igive)
            else
                count = 0
                allocate(igive(0:-1))
            endif
        case (DM_FLUID_CG_PML)
            n1 = Tdomain%Comm_data%Data(n)%nflupml
            if (n1>0.and.n2>0) then
                count = intersect_arrays(n1, Tdomain%Comm_data%Data(n)%IGiveFPML, n2, surface%map, igive)
            else
                count = 0
                allocate(igive(0:-1))
            endif
        end select
        do i=0,count-1
            idx = igive(i)
            idx = renum(idx)
            igive(i) = idx
        enddo
        select case(dom)
        case (DM_SOLID_CG)
            allocate(comm_data%Data(n)%IGiveS(0:count-1))
            comm_data%Data(n)%IGiveS = igive
            comm_data%Data(n)%nsol = count
        case (DM_SOLID_DG)
            allocate(comm_data%Data(n)%IGiveSDG(0:count-1))
            comm_data%Data(n)%IGiveSDG = igive
            comm_data%Data(n)%nsoldg = count
        case (DM_FLUID_CG)
            allocate(comm_data%Data(n)%IGiveF(0:count-1))
            comm_data%Data(n)%IGiveF = igive
            comm_data%Data(n)%nflu = count
        case (DM_SOLID_CG_PML)
            allocate(comm_data%Data(n)%IGiveSPML(0:count-1))
            comm_data%Data(n)%IGiveSPML = igive
            comm_data%Data(n)%nsolpml = count
        case (DM_FLUID_CG_PML)
            allocate(comm_data%Data(n)%IGiveFPML(0:count-1))
            comm_data%Data(n)%IGiveFPML = igive
            comm_data%Data(n)%nflupml = count
        end select
        comm_data%Data(n)%src = src
        comm_data%Data(n)%dest = dst

!        write(*,*) "COMM:", src, "->", dst, ":", count, "(",dom,")"
        deallocate(igive)
    end do
end subroutine build_comms_surface

subroutine prepare_comm_surface(Tdomain, comm_data)
    use sdomain
    implicit none

    type(domain), intent (inout) :: Tdomain
!    type(comm_vector), intent(in) :: comm_src
    type(comm_vector), intent(inout) :: comm_data
!    type(surf_num), intent(in) :: surface
    integer :: ncomm, sz, n

    ncomm = Tdomain%tot_comm_proc
    allocate(Comm_data%Data(0:ncomm-1))
    allocate(Comm_data%send_reqs(0:ncomm-1))
    allocate(Comm_data%recv_reqs(0:ncomm-1))
    Comm_data%ncomm = ncomm

    ! Compte le nb de points GLL commun entre Tdomain%SF%intSolFlu%surf1 et les points de Tdomain
    call build_comms_surface(Tdomain, comm_data, Tdomain%SF%intSolFlu%surf1, DM_FLUID_CG)
    call build_comms_surface(Tdomain, comm_data, Tdomain%SF%intSolFluPml%surf1, DM_FLUID_CG_PML)

    do n = 0, comm_data%ncomm-1
        sz = comm_data%Data(n)%nflu + comm_data%Data(n)%nflupml
        comm_data%Data(n)%ndata = 3*sz

        allocate(comm_data%Data(n)%Give(0:3*sz-1))
        allocate(comm_data%Data(n)%Take(0:3*sz-1))
    end do

end subroutine prepare_comm_surface

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
