!! ============================================================================
!! Solid <-> POTENTIAL-fluid coupling for SEM2D (faithful port of SEM3D).
!!
!! STATUS: WORKING (validated 2026-07-09 on t4/t5: stable, arrivals match the native
!! iso reference within 0.2us, fluid amplitude obeys p = rho.c.v within ~6%).
!!
!! SCOPE: only solid <-> POTENTIAL-fluid (Fluid_Aniso / Cstar-fluid) interfaces. The
!! ISO fluid is a DISPLACEMENT formulation (elastic kernel, mu=0) and couples natively
!! through shared CG DOFs -- it is deliberately NOT split (sf_needs_split).
!!
!! ARCHITECTURE (approach B): at setup, each interface face AND its vertices are
!! DUPLICATED so the solid and fluid sides own separate one-sided objects updated by
!! the STANDARD face/vertex machinery (predictor/assembly/corrector). The coupling is
!! then the exact SEM3D non-CPML sequence (Newmark.F90):
!!   1. sf_stof (after assembly, before correctors):
!!        fluid Forces += BtN . Veloc_solid(old)          [natural BC of the psi eq.]
!!   2. fluid interface DOFs corrected normally in the main loop (VelPhi becomes NEW);
!!      solid interface DOFs are SKIPPED (deferred).
!!   3. sf_ftos: solid Forces -= BtN . VelPhi(new)        [pressure load, p = -VelPhi]
!!      (mirror of SEM3D FtoS: champs(f1)%Veloc there is the FORCE accumulator).
!!   4. deferred standard correction of the solid interface faces/vertices
!!      (M^-1 and dt applied by the correctors themselves).
!! BtN = unit normal OUTWARD FROM THE FLUID * line-Jacobian * 1D GLL weight (mirror of
!! SEM3D SF_BtN orientation), covering face-interior nodes AND endpoint vertices.
!!
!! MPI: bit-exact on any number of ranks. The mesher keeps the S-F interface off partition
!! cuts (heavy solid<->fluid edge weight), so only interface VERTICES at cut crossings are
!! cross-rank. (B.5) sf_augment_walls puts the fluid vertex copy in the comm walls (mass +
!! Forces summed by the standard exchanges); sf_exchange_ftos sums the post-exchange FtoS
!! force delta at shared interface solid vertices. Validated np1==np4 to ~1e-14.
!!
!! Limits: CG non-PML interfaces (the interface faces themselves are not in a PML layer).
!! ============================================================================
module solid_fluid_coupling_2d

    use constants
    use sdomain
    use selement
    use mpi
    implicit none


    ! (B.4) split-DOF coupling: paired solid/fluid faces + endpoint vertices, with BtN.
    type sfpair_t
        integer :: fsol = -1, fflu = -1, ngll = 0
        integer :: vsol(0:1) = -1, vflu(0:1) = -1
        real(fpp), allocatable :: btn(:,:)   ! (0:1, 0:ngll-1) = normal*lineJac*GLLw at face nodes
    end type sfpair_t
    type(sfpair_t), allocatable :: sfp(:)
    integer :: n_sfp = 0

    ! (B.5) cross-rank coupling. The S-F interface is never a partition cut (the mesher gives
    ! solid<->fluid dual-graph edges a heavy weight), so interface FACES stay on one rank; only
    ! interface VERTICES where a partition cut crosses the interface line are shared. For those:
    !   (i)  the fluid-side copy vflu is appended to the comm walls (sf_augment_walls), so the
    !        standard mass + per-step Forces exchanges sum its contributions across ranks;
    !   (ii) the FtoS force delta -- applied AFTER the main exchange -- is summed with a small
    !        dedicated exchange (sf_exchange_ftos) over the interface SOLID vertices.
    type sfwall_t
        integer :: n = 0
        integer, allocatable :: pos(:)   ! positions of interface-solid vertices in sWall%Vertex_List
    end type sfwall_t
    type(sfwall_t), allocatable :: sfw(:)               ! (0:n_communications-1)
    real(fpp), allocatable :: ftos_acc(:,:)             ! (0:1, 0:n_vertex-1) local FtoS delta, per step

contains

    !-----------------------------------------------------------------------
    !> True when the (e0,e1) pair is a solid <-> POTENTIAL-fluid interface: exactly one
    !! side acoustic AND that side's subdomain is a potential formulation (Fluid_Aniso /
    !! Cstar fluid). ISO fluid (displacement, mu=0 elastic) returns false: it couples
    !! natively through shared CG DOFs.
    logical function sf_needs_split(Tdomain, e0, e1) result(res)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: e0, e1
        logical :: a0, a1
        a0 = Tdomain%specel(e0)%acoustic; a1 = Tdomain%specel(e1)%acoustic
        res = (a0 .neqv. a1)
    end function sf_needs_split

    !-----------------------------------------------------------------------
    !> (B.1) Duplicate each solid-fluid interface face so the solid and fluid elements own
    !! SEPARATE one-sided face objects (Near_Element(1)=-1), instead of sharing one Face.
    !! MUST run after read_mesh (acoustic/changing_media set) and BEFORE global_numbering.
    !! No-op when there is no S-F interface (other models untouched). Sets is_sf_iface so
    !! PML_definition does NOT flag the split faces as free-surface.
    subroutine split_sf_interface_faces(Tdomain)
        type(domain), intent(inout) :: Tdomain
        type(face), dimension(:), pointer :: newf
        integer :: nf, e0, e1, esol, eflu, nnew, idx, oldn, ed, wf0, wf1
        logical :: a0, a1

        ! Only solid <-> POTENTIAL-fluid interfaces are split. The ISO fluid is a
        ! displacement formulation (elastic kernel, mu=0): its interface couples natively
        ! through the shared CG DOFs and must NOT be split (see 2026-07-09 solution plan).
        nnew = 0
        do nf = 0, Tdomain%n_face-1
            if (.not. Tdomain%sFace(nf)%changing_media) cycle
            e0 = Tdomain%sFace(nf)%Near_Element(0); e1 = Tdomain%sFace(nf)%Near_Element(1)
            if (e0 < 0 .or. e1 < 0) cycle
            if (.not. sf_needs_split(Tdomain, e0, e1)) cycle
            nnew = nnew + 1
        end do
        if (nnew == 0) return
        allocate(sfp(nnew)); n_sfp = 0   ! record solid/fluid face pairs while splitting

        oldn = Tdomain%n_face
        allocate(newf(0:oldn+nnew-1))
        newf(0:oldn-1) = Tdomain%sFace(0:oldn-1)   ! deep copy (Face allocatable comps not yet allocated)

        idx = oldn - 1
        do nf = 0, oldn-1
            if (.not. newf(nf)%changing_media) cycle
            e0 = newf(nf)%Near_Element(0); e1 = newf(nf)%Near_Element(1)
            if (e0 < 0 .or. e1 < 0) cycle
            if (.not. sf_needs_split(Tdomain, e0, e1)) cycle
            a0 = Tdomain%specel(e0)%acoustic
            if (a0) then; eflu = e0; esol = e1; else; eflu = e1; esol = e0; end if
            wf0 = newf(nf)%Which_face(0); wf1 = newf(nf)%Which_face(1)
            ! mark the two interface-endpoint vertices for the (B.2) vertex split
            Tdomain%sVertex(newf(nf)%Near_Vertex(0))%is_sf_vertex = .true.
            Tdomain%sVertex(newf(nf)%Near_Vertex(1))%is_sf_vertex = .true.

            idx = idx + 1
            newf(idx) = newf(nf)                       ! duplicate for the FLUID side

            ! SOLID side keeps face nf (one-sided)
            newf(nf)%Near_Element(0) = esol; newf(nf)%Near_Element(1) = -1
            newf(nf)%Which_face(0) = merge(wf0, wf1, esol == e0); newf(nf)%Which_face(1) = -1
            newf(nf)%mat_index = Tdomain%specel(esol)%mat_index
            newf(nf)%is_sf_iface = .true.; newf(nf)%changing_media = .false.

            ! FLUID side -> new face idx (one-sided)
            newf(idx)%Near_Element(0) = eflu; newf(idx)%Near_Element(1) = -1
            newf(idx)%Which_face(0) = merge(wf0, wf1, eflu == e0); newf(idx)%Which_face(1) = -1
            newf(idx)%mat_index = Tdomain%specel(eflu)%mat_index
            newf(idx)%is_sf_iface = .true.; newf(idx)%changing_media = .false.

            ! repoint the fluid element's Near_Face from nf -> idx
            do ed = 0, 3
                if (Tdomain%specel(eflu)%Near_Face(ed) == nf) then
                    Tdomain%specel(eflu)%Near_Face(ed) = idx
                    exit
                end if
            end do

            ! record the pair (BtN + endpoint vertices filled later)
            n_sfp = n_sfp + 1
            sfp(n_sfp)%fsol = nf; sfp(n_sfp)%fflu = idx
            sfp(n_sfp)%ngll = newf(nf)%ngll
        end do

        deallocate(Tdomain%sFace)
        Tdomain%sFace => newf
        Tdomain%n_face = oldn + nnew
        if (Tdomain%Mpi_var%my_rank == 0) &
            write(*,'(a,i0,a,i0)') ' [sf] split ', nnew, ' interface faces -> n_face = ', Tdomain%n_face

        ! ---- (B.2) split the interface VERTICES so solid- and fluid-side entities own
        !      separate vertex objects (removes the solid-velocity + fluid-potential mixing).
        block
            type(vertex), dimension(:), pointer :: newv
            integer, allocatable :: vmap(:)
            integer :: nv, oldnv, nvsplit, vidx, kk, ov
            oldnv = Tdomain%n_vertex
            nvsplit = 0
            do nv = 0, oldnv-1
                if (Tdomain%sVertex(nv)%is_sf_vertex) nvsplit = nvsplit + 1
            end do
            if (nvsplit == 0) return
            allocate(newv(0:oldnv+nvsplit-1))
            newv(0:oldnv-1) = Tdomain%sVertex(0:oldnv-1)
            allocate(vmap(0:oldnv-1)); vmap = -1
            vidx = oldnv - 1
            do nv = 0, oldnv-1
                if (.not. newv(nv)%is_sf_vertex) cycle
                vidx = vidx + 1
                newv(vidx) = newv(nv)          ! fluid-side copy of the vertex
                newv(vidx)%is_sf_vertex = .false.   ! flag stays ONLY on the solid side
                vmap(nv) = vidx
            end do
            ! repoint FLUID (acoustic) elements' Near_Vertex to the fluid copy
            do kk = 0, Tdomain%n_elem-1
                if (.not. Tdomain%specel(kk)%acoustic) cycle
                do nv = 0, 3
                    ov = Tdomain%specel(kk)%Near_Vertex(nv)
                    if (vmap(ov) >= 0) then
                        newv(vmap(ov))%mat_index = Tdomain%specel(kk)%mat_index
                        Tdomain%specel(kk)%Near_Vertex(nv) = vmap(ov)
                    end if
                end do
            end do
            ! repoint FLUID faces' Near_Vertex (Near_Element(0) is the fluid element)
            do kk = 0, Tdomain%n_face-1
                if (Tdomain%sFace(kk)%Near_Element(0) < 0) cycle
                if (.not. Tdomain%specel(Tdomain%sFace(kk)%Near_Element(0))%acoustic) cycle
                do nv = 0, 1
                    ov = Tdomain%sFace(kk)%Near_Vertex(nv)
                    if (vmap(ov) >= 0) Tdomain%sFace(kk)%Near_Vertex(nv) = vmap(ov)
                end do
            end do
            deallocate(Tdomain%sVertex)
            Tdomain%sVertex => newv
            Tdomain%n_vertex = oldnv + nvsplit
            ! (i) append each fluid-side vertex copy vflu=vmap(nv) to every comm wall that lists
            ! the original interface vertex nv (now the solid-side vsol). Both ranks iterate their
            ! matched-order Vertex_List, so appended vflu stay position-matched across the pair.
            ! wall_transfer (after the split) sizes the exchange buffers from the new n_vertices.
            call sf_augment_walls(Tdomain, vmap)
            deallocate(vmap)
            ! fill the endpoint vertex pairs (solid face keeps Vsol, fluid face has Vflu)
            do kk = 1, n_sfp
                sfp(kk)%vsol(0:1) = Tdomain%sFace(sfp(kk)%fsol)%Near_Vertex(0:1)
                sfp(kk)%vflu(0:1) = Tdomain%sFace(sfp(kk)%fflu)%Near_Vertex(0:1)
            end do
            if (Tdomain%Mpi_var%my_rank == 0) &
                write(*,'(a,i0,a,i0)') ' [sf] split ', nvsplit, ' interface vertices -> n_vertex = ', Tdomain%n_vertex
        end block
    end subroutine split_sf_interface_faces


    !-----------------------------------------------------------------------
    !> (i) Append the fluid-side vertex copies to the comm walls. vmap(nv) is the new fluid
    !! vertex index for a split interface vertex nv (or <0). For every wall listing nv, append
    !! vmap(nv) at the end (kept in the order the wall's Vertex_List is scanned, so paired ranks
    !! agree). Original entries keep their positions -> other position-matched exchanges (PML
    !! handshake, regular vertex Forces) are undisturbed.
    subroutine sf_augment_walls(Tdomain, vmap)
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: vmap(0:)
        integer :: w, q, nw, cnt, add
        integer, allocatable :: newlist(:)
        do w = 0, Tdomain%n_communications-1
            nw = Tdomain%sWall(w)%n_vertices
            cnt = 0
            do q = 0, nw-1
                if (vmap(Tdomain%sWall(w)%Vertex_List(q)) >= 0) cnt = cnt + 1
            end do
            if (cnt == 0) cycle
            allocate(newlist(0:nw+cnt-1))
            newlist(0:nw-1) = Tdomain%sWall(w)%Vertex_List(0:nw-1)
            add = nw
            do q = 0, nw-1
                if (vmap(Tdomain%sWall(w)%Vertex_List(q)) >= 0) then
                    newlist(add) = vmap(Tdomain%sWall(w)%Vertex_List(q)); add = add + 1
                end if
            end do
            deallocate(Tdomain%sWall(w)%Vertex_List)
            allocate(Tdomain%sWall(w)%Vertex_List(0:nw+cnt-1))
            Tdomain%sWall(w)%Vertex_List(:) = newlist(:)
            Tdomain%sWall(w)%n_vertices = nw + cnt
            deallocate(newlist)
        end do
    end subroutine sf_augment_walls

    !-----------------------------------------------------------------------
    !> (ii) Record, per comm wall, the positions of interface SOLID vertices in Vertex_List, and
    !! allocate the FtoS delta accumulator. Call once, after the walls are final (post-split /
    !! wall_transfer). is_sf_vertex is set symmetrically by the split on both ranks that share an
    !! interface vertex, so the recorded positions are position-matched across the pair.
    subroutine sf_build_comm(Tdomain)
        type(domain), intent(inout) :: Tdomain
        integer :: w, q, nw, cnt
        if (Tdomain%n_communications <= 0 .or. n_sfp == 0) return
        allocate(sfw(0:Tdomain%n_communications-1))
        do w = 0, Tdomain%n_communications-1
            nw = Tdomain%sWall(w)%n_vertices
            cnt = 0
            do q = 0, nw-1
                if (Tdomain%sVertex(Tdomain%sWall(w)%Vertex_List(q))%is_sf_vertex) cnt = cnt + 1
            end do
            sfw(w)%n = cnt
            allocate(sfw(w)%pos(0:max(cnt-1,0)))
            cnt = 0
            do q = 0, nw-1
                if (Tdomain%sVertex(Tdomain%sWall(w)%Vertex_List(q))%is_sf_vertex) then
                    sfw(w)%pos(cnt) = q; cnt = cnt + 1
                end if
            end do
        end do
        allocate(ftos_acc(0:1, 0:Tdomain%n_vertex-1)); ftos_acc = 0._fpp
    end subroutine sf_build_comm

    !> (ii) Sum the FtoS force delta across ranks at shared interface solid vertices. Call right
    !! after sf_ftos, before the deferred solid-interface correction. Each rank applied only its
    !! own BtN.VelPhi; this adds the neighbour's, so the deferred corrector sees the full load.
    subroutine sf_exchange_ftos(Tdomain)
        type(domain), intent(inout) :: Tdomain
        integer :: w, k, nv, partner, ierr, cnt, st(MPI_STATUS_SIZE)
        real(fpp), allocatable :: sbuf(:,:), rbuf(:,:)
        if (.not. allocated(sfw)) return
        do w = 0, Tdomain%n_communications-1
            cnt = sfw(w)%n
            if (cnt == 0) cycle
            allocate(sbuf(0:1,0:cnt-1), rbuf(0:1,0:cnt-1)); rbuf = 0._fpp
            do k = 0, cnt-1
                nv = Tdomain%sWall(w)%Vertex_List(sfw(w)%pos(k))
                sbuf(0:1,k) = ftos_acc(0:1, nv)
            end do
            partner = Tdomain%Communication_list(w)
            call MPI_Sendrecv(sbuf, 2*cnt, MPI_DOUBLE_PRECISION, partner, 830, &
                              rbuf, 2*cnt, MPI_DOUBLE_PRECISION, partner, 830, &
                              Tdomain%communicateur, st, ierr)
            do k = 0, cnt-1
                nv = Tdomain%sWall(w)%Vertex_List(sfw(w)%pos(k))
                Tdomain%sVertex(nv)%Forces(0:1) = Tdomain%sVertex(nv)%Forces(0:1) + rbuf(0:1,k)
            end do
            deallocate(sbuf, rbuf)
        end do
    end subroutine sf_exchange_ftos

    !-----------------------------------------------------------------------
    !> (B.4) Compute BtN = normal(solid->fluid) * lineJac * GLLw at every interface-face node
    !! (interior + endpoints). Call after geometry is built (shape/define_arrays).
    subroutine build_sf_coupling(Tdomain)
        use shape_lin, only: compute_Jacobian_1D, n_from_vertices
        type(domain), intent(inout) :: Tdomain
        integer :: m, p, fsol, ngll, mat, wf
        real(fpp) :: fn(0:1), Jac1D
        if (n_sfp == 0) return
        do m = 1, n_sfp
            fsol = sfp(m)%fsol; ngll = sfp(m)%ngll
            wf = Tdomain%sFace(fsol)%Which_face(0); mat = Tdomain%sFace(fsol)%mat_index
            allocate(sfp(m)%btn(0:1, 0:ngll-1)); sfp(m)%btn = 0._fpp
            call n_from_vertices(Tdomain, fn, Tdomain%sFace(fsol)%Near_Vertex(0), &
                                 Tdomain%sFace(fsol)%Near_Vertex(1))
            if (wf >= 2) fn = -fn                       ! outward from the solid element
            fn = -fn   ! flip: outward from the FLUID (mirror SEM3D SF_BtN orientation),
                       ! so StoF/FtoS below can use SEM3D's formulas verbatim.
            call compute_Jacobian_1D(Tdomain, fsol, Jac1D)
            do p = 0, ngll-1
                if (wf == 0 .or. wf == 2) then
                    sfp(m)%btn(0,p) = fn(0)*Jac1D*Tdomain%sSubdomain(mat)%GLLwx(p)
                    sfp(m)%btn(1,p) = fn(1)*Jac1D*Tdomain%sSubdomain(mat)%GLLwx(p)
                else
                    sfp(m)%btn(0,p) = fn(0)*Jac1D*Tdomain%sSubdomain(mat)%GLLwz(p)
                    sfp(m)%btn(1,p) = fn(1)*Jac1D*Tdomain%sSubdomain(mat)%GLLwz(p)
                end if
            end do
        end do
        if (Tdomain%Mpi_var%my_rank == 0) &
            write(*,'(a,i0,a)') ' [sf] coupling built on ', n_sfp, ' face pairs'
        call sf_build_comm(Tdomain)   ! (ii) cross-rank FtoS exchange tables
    end subroutine build_sf_coupling

    !> StoF: inject the solid normal velocity into the fluid RHS. Call AFTER assembly,
    !! BEFORE the (standard) corrector loop. Handles face-interior nodes + endpoint vertices.
    subroutine sf_stof(Tdomain)
        type(domain), intent(inout) :: Tdomain
        integer :: m, p, ngll, fsol, fflu
        real(fpp) :: vn
        do m = 1, n_sfp
            ngll = sfp(m)%ngll; fsol = sfp(m)%fsol; fflu = sfp(m)%fflu
            do p = 1, ngll-2
                vn = sfp(m)%btn(0,p)*Tdomain%sFace(fsol)%Veloc(p,0) &
                   + sfp(m)%btn(1,p)*Tdomain%sFace(fsol)%Veloc(p,1)
                Tdomain%sFace(fflu)%Forces(p,0) = Tdomain%sFace(fflu)%Forces(p,0) + vn
            end do
            call stof_vertex(Tdomain, sfp(m)%vsol(0), sfp(m)%vflu(0), sfp(m)%btn(:,0))
            call stof_vertex(Tdomain, sfp(m)%vsol(1), sfp(m)%vflu(1), sfp(m)%btn(:,ngll-1))
        end do
    end subroutine sf_stof

    subroutine stof_vertex(Tdomain, vs, vf, btn)
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: vs, vf
        real(fpp), intent(in) :: btn(0:1)
        Tdomain%sVertex(vf)%Forces(0) = Tdomain%sVertex(vf)%Forces(0) &
            + btn(0)*Tdomain%sVertex(vs)%Veloc(0) + btn(1)*Tdomain%sVertex(vs)%Veloc(1)
    end subroutine stof_vertex

    !> FtoS: the fluid pressure (p = -VelPhi) loads the solid as a FORCE contribution
    !! (mirror SEM3D FtoS_coupling: champs(f1)%Veloc there is the FORCE accumulator --
    !! newmark_predictor_solid zeroes it, newmark_corrector_solid applies M^-1 and dt).
    !! Timing (SEM3D order): AFTER the fluid DOFs are corrected (VelPhi is NEW), BEFORE
    !! the solid interface DOFs are corrected -- so those corrections are DEFERRED in
    !! Newmark and run right after this call, via the standard correctors.
    subroutine sf_ftos(Tdomain)
        type(domain), intent(inout) :: Tdomain
        integer :: m, p, ngll, fsol, fflu
        real(fpp) :: vphi
        if (allocated(ftos_acc)) ftos_acc = 0._fpp   ! (ii) reset this step's local FtoS delta
        do m = 1, n_sfp
            ngll = sfp(m)%ngll; fsol = sfp(m)%fsol; fflu = sfp(m)%fflu
            do p = 1, ngll-2
                vphi = Tdomain%sFace(fflu)%Veloc(p,0)      ! NEW VelPhi (fluid corrected)
                Tdomain%sFace(fsol)%Forces(p,0) = Tdomain%sFace(fsol)%Forces(p,0) - sfp(m)%btn(0,p)*vphi
                Tdomain%sFace(fsol)%Forces(p,1) = Tdomain%sFace(fsol)%Forces(p,1) - sfp(m)%btn(1,p)*vphi
            end do
            call ftos_vertex(Tdomain, sfp(m)%vsol(0), sfp(m)%vflu(0), sfp(m)%btn(:,0))
            call ftos_vertex(Tdomain, sfp(m)%vsol(1), sfp(m)%vflu(1), sfp(m)%btn(:,ngll-1))
        end do
    end subroutine sf_ftos

    subroutine ftos_vertex(Tdomain, vs, vf, btn)
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: vs, vf
        real(fpp), intent(in) :: btn(0:1)
        real(fpp) :: vphi
        vphi = Tdomain%sVertex(vf)%Veloc(0)                ! NEW VelPhi at the fluid vertex
        Tdomain%sVertex(vs)%Forces(0) = Tdomain%sVertex(vs)%Forces(0) - btn(0)*vphi
        Tdomain%sVertex(vs)%Forces(1) = Tdomain%sVertex(vs)%Forces(1) - btn(1)*vphi
        ! (ii) record the local delta so shared interface solid vertices can sum it across ranks
        if (allocated(ftos_acc)) then
            ftos_acc(0,vs) = ftos_acc(0,vs) - btn(0)*vphi
            ftos_acc(1,vs) = ftos_acc(1,vs) - btn(1)*vphi
        end if
    end subroutine ftos_vertex

end module solid_fluid_coupling_2d
