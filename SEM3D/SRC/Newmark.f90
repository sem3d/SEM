!>
!!\file Newmark.f90
!!\brief Algorithme de Newmark
!! La routine Newmark assure la r√©solution des √©quations via un algorithme de predicteur-multi-correcteur
!! des vitesses avec une formulation contrainte-vitesse d√©cal√©e en temps dans les PML.
!<

subroutine Newmark(Tdomain,ntime)
    ! Predictor-MultiCorrector Newmark Velocity Scheme within a
    ! Time staggered Stress-Velocity formulation inside PML
    use sdomain
#ifdef COUPLAGE
    use sCouplage
#endif
    use mcapteur
    use mpi
    use scomm
    use scommutils
    use orientation
    use assembly
    use schamps

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: ntime
    integer :: n, mat
    integer :: nf, ne, nv
    integer :: nf_aus, ne_aus, nv_aus
    integer :: ngll, ngll1, ngll2, ngllPML, ngll_F, ngllPML_F
    integer, parameter :: etiquette = 100
#ifdef COUPLAGE
    integer :: i,j
#endif


    ! Predictor-MultiCorrector Newmark Velocity Scheme within a
    ! Time staggered Stress-Velocity formulation inside PML
    ! PML needs to be implemented
    if(.not. Tdomain%TimeD%velocity_scheme)   &
        stop "Newmark scheme implemented only in velocity form."


    !- Prediction Phase
    call Newmark_Predictor(Tdomain,Tdomain%champs1)

    !- Solution phase
    call internal_forces(Tdomain,Tdomain%champs1)


    ! External Forces
    if(Tdomain%logicD%any_source)then
        call external_forces(Tdomain,Tdomain%TimeD%rtime,ntime,Tdomain%champs1)
    end if

#ifdef COUPLAGE
#if 0
    if (ntime>0) then
        call calcul_couplage_force(Tdomain, ntime)
    endif

    !Gsa Ipsis (tout le passage)
    ! AJOUT DES FORCES de couplage

    ! Les forces de couplage correspondent a la contrainte multiplie par la surface
    ! du point de gauss correspondant sur un meme proc la somme est
    ! deja faite par contre lorsqu un vertex est partage sur plusieurs
    ! proc alors chaque proc n a qu une partie de la somme il faut
    ! donc lui ajouter les contributions des autres proc pour prendre
    ! en compte les forces imposee lors du couplage avec mka sur les
    ! points de gauss internes aux faces
    do nf = 0, Tdomain%n_face-1
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2
        do j=1,ngll2-2
            do i=1,ngll1-2
                Tdomain%sFace(nf)%Forces(i,j,0:2) = Tdomain%sFace(nf)%ForcesExt(i,j,0:2) + Tdomain%sFace(nf)%Forces(i,j,0:2)
            enddo
        enddo

    enddo

    ! aretes
    do nf = 0, Tdomain%n_edge-1
        ngll = Tdomain%sEdge(nf)%ngll
        do i=1,ngll-2
            Tdomain%sEdge(nf)%Forces(i,0:2) = Tdomain%sEdge(nf)%ForcesExt(i,0:2) + Tdomain%sEdge(nf)%Forces(i,0:2)
        enddo
    enddo

    ! pour prendre en compte les forces imposee lors du couplage avec mka sur les points de gauss des vertex
    do nv = 0, Tdomain%n_vertex-1
        Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%ForcesExt(0:2) + Tdomain%sVertex(nv)%Forces(0:2)
    enddo
#endif
#endif


    ! MPI communications
    call comm_forces(Tdomain)

    ! Neumann B.C.: associated forces
    if(Tdomain%logicD%neumann_local_present)then
#if 1
    !TODO
    stop "TODO: conditions de Neumann non prises en compte (Newmark)"
#else
        mat = Tdomain%Neumann%Neu_Param%mat_index
        do nf = 0,Tdomain%Neumann%Neu_n_faces-1
            ngll1 = Tdomain%Neumann%Neu_Face(nf)%ngll1
            ngll2 = Tdomain%Neumann%Neu_Face(nf)%ngll2
            call compute_Neu_forces_on_face(Tdomain%Neumann%Neu_Face(nf),     &
                Tdomain%Neumann%Neu_Param,Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime)
        enddo
        do ne = 0,Tdomain%Neumann%Neu_n_edges-1
            ngll = Tdomain%Neumann%Neu_Edge(ne)%ngll
            ne_aus = Tdomain%Neumann%Neu_Edge(ne)%Edge
            call compute_Neu_forces_on_edge(Tdomain%Neumann%Neu_Edge(ne),     &
                Tdomain%Neumann%Neu_Param,Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime)
        enddo
        do nv = 0,Tdomain%Neumann%Neu_n_vertices-1
            nv_aus = Tdomain%Neumann%Neu_Vertex(nv)%Vertex
            n = merge(0,1,nv == 4)
            call compute_Neu_forces_on_vertex(Tdomain%Neumann%Neu_Vertex(nv),n,  &
                Tdomain%Neumann%Neu_Param,Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime)
        enddo
        ! addition of Neumann forces
        do nf = 0, Tdomain%Neumann%Neu_n_faces-1
            nf_aus = Tdomain%Neumann%Neu_Face(nf)%Face
            if(.not.Tdomain%sface(nf_aus)%PML)then
                Tdomain%sFace(nf_aus)%Forces(:,:,0:2) = Tdomain%sFace(nf_aus)%Forces(:,:,0:2) - &
                    Tdomain%Neumann%Neu_Face(nf)%Forces(:,:,0:2)
            else
                Tdomain%sFace(nf_aus)%spml%Forces3(:,:,0:2) = Tdomain%sFace(nf_aus)%spml%Forces3(:,:,0:2) - &
                    Tdomain%Neumann%Neu_Face(nf)%Forces(:,:,0:2)
            endif
        enddo
        do ne = 0, Tdomain%Neumann%Neu_n_edges-1
            ne_aus = Tdomain%Neumann%Neu_Edge(ne)%Edge
            if(.not.Tdomain%sedge(ne_aus)%PML)then
                Tdomain%sEdge(ne_aus)%Forces(:,0:2) = Tdomain%sEdge(ne_aus)%Forces(:,0:2) - &
                    Tdomain%Neumann%Neu_Edge(ne)%Forces(:,0:2)
            else
                Tdomain%sEdge(ne_aus)%spml%Forces3(:,0:2) = Tdomain%sEdge(ne_aus)%spml%Forces3(:,0:2) - &
                    Tdomain%Neumann%Neu_Edge(ne)%Forces(:,0:2)
            endif
        enddo
        do nv = 0, Tdomain%Neumann%Neu_n_vertices-1
            nv_aus = Tdomain%Neumann%Neu_Vertex(nv)%Vertex
            if(.not.Tdomain%svertex(nv_aus)%PML)then
                Tdomain%sVertex(nv_aus)%Forces(0:2) = Tdomain%sVertex(nv_aus)%Forces(0:2) -  &
                    Tdomain%Neumann%Neu_Vertex(nv)%Forces(0:2)
            else
                Tdomain%sVertex(nv_aus)%spml%Forces3(0:2) = Tdomain%sVertex(nv_aus)%spml%Forces3(0:2) - &
                    Tdomain%Neumann%Neu_Vertex(nv)%Forces(0:2)
            endif
        enddo
#endif
    endif


    !- solid -> fluid coupling (normal dot velocity)
    if(Tdomain%logicD%SF_local_present)then
        call SF_solid_values_saving(Tdomain%ngll_s, Tdomain%SF%ngll, Tdomain%SF%SF_IGlobSol, &
                                      Tdomain%champs1%Forces, Tdomain%champs1%Depla, &
                                      Tdomain%champs0%Save_forces, &
                                      Tdomain%champs0%Save_depla)
        call StoF_coupling(Tdomain, Tdomain%ngll_s, Tdomain%ngll_f, Tdomain%SF%ngll, Tdomain%SF%SF_IGlobSol, &
            Tdomain%SF%SF_IGlobFlu, Tdomain%champs0%Veloc, &
            Tdomain%SF%SF_BtN, Tdomain%champs1%ForcesFl)
    end if


    !- correction phase
    call Newmark_Corrector(Tdomain,Tdomain%champs1)

    if(Tdomain%logicD%SF_local_present)then
        !- fluid -> solid coupling (pressure times velocity)
        call FtoS_coupling(Tdomain, Tdomain%ngll_s, Tdomain%ngll_f, Tdomain%SF%ngll, &
                             Tdomain%SF%SF_IGlobSol, Tdomain%SF%SF_IGlobFlu, &
                             Tdomain%SF%SF_BtN, Tdomain%champs0%Save_forces, &
                             Tdomain%champs0%Save_depla, Tdomain%champs0%VelPhi, &
                             Tdomain%champs0%Forces, Tdomain%champs0%Depla)
        !- recorrecting on solid faces
        call Newmark_recorrect_solid(Tdomain%ngll_s, Tdomain%SF%ngll, Tdomain%TimeD%dtmin, &
                                       Tdomain%SF%SF_IGlobSol, Tdomain%MassMatSol, &
                                       Tdomain%champs0%Forces, Tdomain%champs0%Veloc, &
                                       Tdomain%champs0%Depla)
    end if

    if (Tdomain%rank==0 .and. mod(ntime,20)==0) print *,' Iteration  =  ',ntime,'    temps  = ',Tdomain%TimeD%rtime

    return

end subroutine Newmark
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine comm_forces(Tdomain)
    use sdomain
    use scomm
    implicit none

    type(domain), intent(inout)   :: Tdomain

    integer :: n, i, j, k, m, idx

    if(Tdomain%Comm_data%ncomm > 0)then
        do n = 0,Tdomain%Comm_data%ncomm-1
            ! Domain SOLID
            k = 0
            do i=0,Tdomain%Comm_data%Data(n)%nsol-1
                idx = Tdomain%Comm_data%Data(n)%IGiveS(i)
                do j=0,2
                    Tdomain%Comm_data%Data(n)%Give(k) = Tdomain%champs1%Forces(idx,j)
                    k=k+1
                enddo
            end do

            ! Domain SOLID PML
            do i=0,Tdomain%Comm_data%Data(n)%nsolpml-1
                idx = Tdomain%Comm_data%Data(n)%IGiveSPML(i)
                do j=0,2
                    do m=0,2
                        Tdomain%Comm_data%Data(n)%Give(k) = Tdomain%champs1%ForcesPML(idx+j,m)
                        k=k+1
                    enddo
                enddo
            end do

            ! Domain FLUID
            do i=0,Tdomain%Comm_data%Data(n)%nflu-1
                idx = Tdomain%Comm_data%Data(n)%IGiveF(i)
                Tdomain%Comm_data%Data(n)%Give(k) = Tdomain%champs1%ForcesFl(idx)
                k=k+1
            end do

            ! Domain FLUID PML
            do i=0,Tdomain%Comm_data%Data(n)%nflupml-1
                idx = Tdomain%Comm_data%Data(n)%IGiveFPML(i)
                do j=0,2
                    Tdomain%Comm_data%Data(n)%Give(k) = Tdomain%champs1%fpml_Forces(idx+j)
                    k=k+1
                enddo
            end do

            Tdomain%Comm_data%Data(n)%nsend = k
        end do
        
        ! Exchange
        call exchange_sem_var(Tdomain, 104, Tdomain%Comm_data)

        ! Take
        do n = 0,Tdomain%Comm_data%ncomm-1
            ! Domain SOLID
            k = 0
            do i=0,Tdomain%Comm_data%Data(n)%nsol-1
                idx = Tdomain%Comm_data%Data(n)%ITakeS(i)
                do j=0,2
                    Tdomain%champs1%Forces(idx,j) = Tdomain%champs1%Forces(idx,j) + &
                                                    Tdomain%Comm_data%Data(n)%Take(k)
                    k = k + 1
                enddo
            end do

            ! Domain SOLID PML
            do i=0,Tdomain%Comm_data%Data(n)%nsolpml-1
                idx = Tdomain%Comm_data%Data(n)%ITakeSPML(i)
                do j=0,2
                    do m=0,2
                        Tdomain%champs1%ForcesPML(idx+j,m) = Tdomain%champs1%ForcesPML(idx+j,m) + &
                                                             Tdomain%Comm_data%Data(n)%Take(k)
                        k=k+1
                    enddo
                enddo
            end do

            ! Domain FLUID
            do i=0,Tdomain%Comm_data%Data(n)%nflu-1
                idx = Tdomain%Comm_data%Data(n)%ITakeF(i)
                Tdomain%champs1%ForcesFl(idx) = Tdomain%champs1%ForcesFl(idx) + &
                                                Tdomain%Comm_data%Data(n)%Take(k)
                k=k+1
            end do

            ! Domain FLUID PML
            do i=0,Tdomain%Comm_data%Data(n)%nflupml-1
                idx = Tdomain%Comm_data%Data(n)%ITakeFPML(i)
                do j=0,2
                    Tdomain%champs1%fpml_Forces(idx+j) = Tdomain%champs1%fpml_Forces(idx+j) + &
                        Tdomain%Comm_data%Data(n)%Take(k)
                    k=k+1
                 enddo
            end do
        end do
    endif

    return

end subroutine comm_forces

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine Newmark_Predictor(Tdomain,champs1)
    use schamps
    use sdomain
    implicit none

    type(domain), intent(inout)   :: Tdomain
    type(champs), intent(inout) :: champs1
    integer :: nel, mat, n, indsol, indpml, indflu
    real :: bega, dt

    bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
    dt = Tdomain%TimeD%dtmin

    ! Elements solide
    if (Tdomain%ngll_s /= 0) then
        champs1%Depla = Tdomain%champs0%Depla
        champs1%Veloc = Tdomain%champs0%Veloc
        champs1%Forces = 0d0
    endif

    ! Elements fluide
    if (Tdomain%ngll_f /= 0) then
        champs1%VelPhi = Tdomain%champs0%VelPhi
        champs1%Phi = Tdomain%champs0%Phi
        champs1%ForcesFl = 0d0
    endif

    ! Elements solide pml
    if (Tdomain%ngll_pmls /= 0) then
        champs1%ForcesPML = 0.
        do n = 0,Tdomain%nbInterfSolPml-1
            ! Couplage ‡ l'interface solide / PML
            indsol = Tdomain%InterfSolPml(n,0)
            indpml = Tdomain%InterfSolPml(n,1)
            Tdomain%champs0%VelocPML(indpml,:) = Tdomain%champs0%Veloc(indsol,:)
            Tdomain%champs0%VelocPML(indpml+1,:) = 0.
            Tdomain%champs0%VelocPML(indpml+2,:) = 0.
        enddo

        do n = 0, Tdomain%nbOuterSPMLNodes-1
            indpml = Tdomain%OuterSPMLNodes(n)
            Tdomain%champs0%VelocPML(indpml,:) = 0.
            Tdomain%champs0%VelocPML(indpml+1,:) = 0.
            Tdomain%champs0%VelocPML(indpml+2,:) = 0.
        enddo

        do nel = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(nel)%Solid) cycle
            if (.not. Tdomain%specel(nel)%PML) cycle

            mat = Tdomain%specel(nel)%mat_index
            if (Tdomain%curve) then
                ! TODO
            else
                if (Tdomain%specel(nel)%FPML) then
                    ! TODO
                else
                    call Prediction_Elem_PML_Veloc(Tdomain%specel(nel), bega, dt, &
                        Tdomain%sSubDomain(mat)%hPrimex, &
                        Tdomain%sSubDomain(mat)%hPrimey, &
                        Tdomain%sSubDomain(mat)%hprimez, &
                        Tdomain%ngll_pmls, &
                        Tdomain%champs0%VelocPML, &
                        Tdomain%champs1%ForcesPML)
                endif
            endif
        enddo
    endif

    ! Elements fluide pml
    if (Tdomain%ngll_pmlf /= 0) then
        champs1%fpml_Forces = 0.
        do n = 0,Tdomain%nbInterfFluPml-1
            ! Couplage ‡ l'interface fluide / PML
            indflu = Tdomain%InterfFluPml(n,0)
            indpml = Tdomain%InterfFluPml(n,1)
            Tdomain%champs0%fpml_VelPhi(indpml) = Tdomain%champs0%VelPhi(indflu)
            Tdomain%champs0%fpml_VelPhi(indpml+1) = 0.
            Tdomain%champs0%fpml_VelPhi(indpml+2) = 0.
        enddo

        do n = 0, Tdomain%nbOuterFPMLNodes-1
            indpml = Tdomain%OuterFPMLNodes(n)
            Tdomain%champs0%fpml_VelPhi(indpml) = 0.
            Tdomain%champs0%fpml_VelPhi(indpml+1) = 0.
            Tdomain%champs0%fpml_VelPhi(indpml+2) = 0.
        enddo

        do nel = 0,Tdomain%n_elem-1
            if (Tdomain%specel(nel)%Solid) cycle
            if (.not. Tdomain%specel(nel)%PML) cycle

            mat = Tdomain%specel(nel)%mat_index
            if (Tdomain%curve) then
                ! TODO
            else
                if (Tdomain%specel(nel)%FPML) then
                    ! TODO
                else
                    call Prediction_Elem_PML_VelPhi(Tdomain%specel(nel), bega, dt, &
                        Tdomain%sSubDomain(mat)%hPrimex, &
                        Tdomain%sSubDomain(mat)%hPrimey, &
                        Tdomain%sSubDomain(mat)%hprimez, &
                        Tdomain%ngll_pmlf, &
                        Tdomain%champs0%fpml_VelPhi, &
                        Tdomain%champs1%fpml_Forces)
                endif
            endif
        enddo
    endif

    return

end subroutine Newmark_Predictor
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine Newmark_Corrector(Tdomain,champs1)
    use sdomain
    use schamps
    implicit none

    type(domain), intent(inout)   :: Tdomain
    type(champs), intent(in)   :: champs1
    integer  :: n, i_dir
    double precision :: dt

    dt = Tdomain%TimeD%dtmin
    ! Si il existe des ÈlÈments PML solide
    if (Tdomain%ngll_pmls /= 0) then
        do i_dir = 0,2
            Tdomain%champs0%VelocPML(:,i_dir) = Tdomain%champs0%DumpV(:,0) * &
                                                Tdomain%champs0%VelocPML(:,i_dir) + &
                                                dt * &
                                                Tdomain%champs0%DumpV(:,1) * &
                                                champs1%ForcesPML(:,i_dir)
        enddo
        !TODO Eventuellement : DeplaPML(:,:) = DeplaPML(:,:) + dt * VelocPML(:,:)
    endif

    ! Si il existe des ÈlÈments PML fluide
    if (Tdomain%ngll_pmlf /= 0) then
        Tdomain%champs0%fpml_VelPhi(:) = Tdomain%champs0%fpml_DumpV(:,0) * &
            Tdomain%champs0%fpml_VelPhi(:) + &
            dt * &
            Tdomain%champs0%fpml_DumpV(:,1) * &
            champs1%fpml_Forces(:)
    endif

    ! Si il existe des ÈlÈments solide
    if (Tdomain%ngll_s /= 0) then
        do i_dir = 0,2
            Tdomain%champs0%Forces(:,i_dir) = champs1%Forces(:,i_dir) * Tdomain%MassMatSol(:)
        enddo
        Tdomain%champs0%Veloc = Tdomain%champs0%Veloc + dt * Tdomain%champs0%Forces
        Tdomain%champs0%Depla = Tdomain%champs0%Depla + dt * Tdomain%champs0%Veloc
    endif

    ! Si il existe des ÈlÈments fluide
    if (Tdomain%ngll_f /= 0) then
        Tdomain%champs0%ForcesFl = champs1%ForcesFl * Tdomain%MassMatFlu
        Tdomain%champs0%VelPhi = (Tdomain%champs0%VelPhi + dt * Tdomain%champs0%ForcesFl) * Tdomain%champs0%Fluid_dirich
        Tdomain%champs0%Phi = Tdomain%champs0%Phi + dt * Tdomain%champs0%VelPhi
    endif

    return
end subroutine Newmark_Corrector
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine internal_forces(Tdomain,champs1)
    ! volume forces - depending on rheology
    use sdomain
    use schamps
    use forces_aniso
    implicit none

    type(domain), intent(inout)  :: Tdomain
    type(champs), intent(inout) :: champs1
    integer  :: n,mat, indsol, indflu, indpml

    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        if(.not. Tdomain%specel(n)%PML)then
            call forces_int(Tdomain%specel(n), Tdomain%sSubDomain(mat),           &
                Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                Tdomain%sSubDomain(mat)%hTprimez, Tdomain%n_sls,Tdomain%aniso,     &
                Tdomain%specel(n)%solid, champs1)
        else ! PML
            if (Tdomain%specel(n)%solid)then
                call compute_InternalForces_PML_Elem(Tdomain%specel(n),                &
                    Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimey, &
                    Tdomain%sSubDomain(mat)%hTprimez, Tdomain%ngll_pmls, champs1%ForcesPML)
            else
                call compute_InternalForces_PML_Elem_Fl(Tdomain%specel(n),                &
                    Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimey, &
                    Tdomain%sSubDomain(mat)%hTprimez, Tdomain%ngll_pmlf, champs1%fpml_Forces)
            endif
        endif
    enddo

    ! Couplage interface solide / PML
    if (Tdomain%ngll_pmls > 0) then
        do n = 0,Tdomain%nbInterfSolPml-1
            indsol = Tdomain%InterfSolPml(n,0)
            indpml = Tdomain%InterfSolPml(n,1)
            champs1%Forces(indsol,:) = champs1%Forces(indsol,:) + &
                                       champs1%ForcesPML(indpml,:) + &
                                       champs1%ForcesPML(indpml+1,:) + &
                                       champs1%ForcesPML(indpml+2,:)
        enddo
    endif
    ! Couplage interface fluid / PML
    if (Tdomain%ngll_pmlf > 0) then
        do n = 0,Tdomain%nbInterfFluPml-1
            indflu = Tdomain%InterfFluPml(n,0)
            indpml = Tdomain%InterfFluPml(n,1)
            champs1%ForcesFl(indflu) = champs1%ForcesFl(indflu) + &
                                       champs1%fpml_Forces(indpml) + &
                                       champs1%fpml_Forces(indpml+1) + &
                                       champs1%fpml_Forces(indpml+2)
        enddo
    endif

    return
end subroutine internal_forces
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine external_forces(Tdomain,timer,ntime,champs1)
    use sdomain
    use schamps
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)  :: ntime
    real, intent(in)  :: timer
    type(champs), intent(inout)  :: champs1
    integer  :: ns,nel,i_dir, i,j,k, idx
    real :: t, ft

    do ns = 0, Tdomain%n_source-1
        if(Tdomain%rank == Tdomain%sSource(ns)%proc)then
            nel = Tdomain%Ssource(ns)%elem

            !  vieille version:
            ! time : t_(n+1/2) for solid ; t_n for fluid
             ! t = merge(timer+Tdomain%TimeD%dtmin/2d0,timer,Tdomain%specel(nel)%solid)
            ! nouvelle version:
            ! le temps n'est plus decale pour les sources, pour un saute-mouton
            !   on rajoute le 1/2 pas de temps qui correspond au fait que la
            !    exterieure doive etre prise a t_(n+1/2)
            t = timer+Tdomain%TimeD%dtmin/2d0
            !
            ft = CompSource(Tdomain%sSource(ns), t, ntime)
            if(Tdomain%sSource(ns)%i_type_source == 1 .or. Tdomain%sSource(ns)%i_type_source == 2) then
                ! collocated force in solid
                !
                do i_dir = 0,2
                    do k = 0,Tdomain%specel(nel)%ngllz-1
                        do j = 0,Tdomain%specel(nel)%nglly-1
                            do i = 0,Tdomain%specel(nel)%ngllx-1
                                idx = Tdomain%specel(nel)%Isol(i,j,k)
                                champs1%Forces(idx, i_dir) = champs1%Forces(idx, i_dir) + &
                                    ft*Tdomain%sSource(ns)%ExtForce(i,j,k,i_dir)
                            enddo
                        enddo
                    enddo
                enddo
            else if(Tdomain%sSource(ns)%i_type_source == 3)then    ! pressure pulse in fluid
                do k = 0,Tdomain%specel(nel)%ngllz-1
                    do j = 0,Tdomain%specel(nel)%nglly-1
                        do i = 0,Tdomain%specel(nel)%ngllx-1
                            idx = Tdomain%specel(nel)%Iflu(i,j,k)
                            champs1%ForcesFl(idx) = champs1%ForcesFl(idx) +    &
                                ft*Tdomain%sSource(ns)%ExtForce(i,j,k,0)
                        enddo
                    enddo
                enddo
            end if
        endif
    enddo

    return
end subroutine external_forces
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
