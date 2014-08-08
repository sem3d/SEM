!>
!!\file Newmark.f90
!!\brief Algorithme de Newmark
!!\author
!!\version 1.0
!!\date 10/03/2009
!! La routine Newmark assure la rï¿½solution des ï¿½quations via un algorithme de predicteur-multi-correcteur
!! des vitesses avec une formulation contrainte-vitesse dï¿½calï¿½e en temps dans les PML.
!<

subroutine Newmark(Tdomain,rg,ntime)
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
    integer, intent(in) :: rg,ntime

    integer :: n, mat,code
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
#if NEW_GLOBAL_METHOD
    call Newmark_Predictor2(Tdomain,Tdomain%champs1)
#else
    call Newmark_Predictor(Tdomain,rg)
#endif

    !- Solution phase
#if NEW_GLOBAL_METHOD
    call internal_forces2(Tdomain,Tdomain%champs1)
#else
    call internal_forces(Tdomain,rg)
#endif


    ! External Forces
    if(Tdomain%logicD%any_source)then
#if NEW_GLOBAL_METHOD
        call external_forces2(Tdomain,Tdomain%TimeD%rtime,ntime,rg,Tdomain%champs1)
#else
        call external_forces(Tdomain,Tdomain%TimeD%rtime,ntime,rg)
#endif
    end if

    ! Communication of Forces within a single process
#if ! NEW_GLOBAL_METHOD
    call inside_proc_forces(Tdomain)
#endif



#ifdef COUPLAGE
    if (ntime>0) then
        call calcul_couplage_force(Tdomain, ntime, rg)
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
                Tdomain%sFace(nf)%Forces(i,j,0:2) = Tdomain%sFace(nf)%ForcesMka(i,j,0:2) + Tdomain%sFace(nf)%Forces(i,j,0:2)
            enddo
        enddo

    enddo

    ! aretes
    do nf = 0, Tdomain%n_edge-1
        ngll = Tdomain%sEdge(nf)%ngll
        do i=1,ngll-2
            Tdomain%sEdge(nf)%Forces(i,0:2) = Tdomain%sEdge(nf)%ForcesMka(i,0:2) + Tdomain%sEdge(nf)%Forces(i,0:2)
        enddo
    enddo

    ! pour prendre en compte les forces imposee lors du couplage avec mka sur les points de gauss des vertex
    do nv = 0, Tdomain%n_vertex-1
        Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%ForcesMka(0:2) + Tdomain%sVertex(nv)%Forces(0:2)
    enddo
#endif


    ! MPI communications
    if(Tdomain%n_proc > 1)then
        ! from external faces, edges and vertices to Communication global arrays
        do n = 0,Tdomain%n_proc-1
#if NEW_GLOBAL_METHOD
            call Comm_Forces_Complete_2(n,Tdomain)
#else
            call Comm_Forces_Complete(n,Tdomain)
            call Comm_Forces_PML_Complete(n,Tdomain)
#endif
        end do

        call exchange_sem_forces(Tdomain, rg)

        ! now: assemblage on external faces, edges and vertices
        do n = 0,Tdomain%n_proc-1
#if NEW_GLOBAL_METHOD
            call Comm_Forces_Assembl(Tdomain,n)
#else
            ngll = 0
            ngll_F = 0
            ngllPML = 0
            ngllPML_F = 0
            call Comm_Forces_Face  (Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
            call Comm_Forces_Edge  (Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
            call Comm_Forces_Vertex(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
#endif
        enddo

    endif  ! if nproc > 1

    ! Neumann B.C.: associated forces
    if(Tdomain%logicD%neumann_local_present)then
#if NEW_GLOBAL_METHOD
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
#if NEW_GLOBAL_METHOD
        call SF_solid_values_saving_2(Tdomain%ngll_s, Tdomain%SF%ngll, Tdomain%SF%SF_IGlobSol, &
                                      Tdomain%champs1%Forces, Tdomain%champs1%Depla, &
                                      Tdomain%champs0%Save_forces, &
                                      Tdomain%champs0%Save_depla)
        call StoF_coupling_2(Tdomain%ngll_s, Tdomain%ngll_f, Tdomain%SF%ngll, Tdomain%SF%SF_IGlobSol, &
                             Tdomain%SF%SF_IGlobFlu, Tdomain%champs0%Veloc, &
                             Tdomain%SF%SF_BtN, Tdomain%champs1%ForcesFl)
#else
        call SF_solid_values_saving(Tdomain)
        call StoF_coupling(Tdomain,rg)
#endif
    end if


    !- correction phase
#if NEW_GLOBAL_METHOD
    call Newmark_Corrector2(Tdomain,Tdomain%champs1)
#else
    call Newmark_Corrector(Tdomain,rg)
#endif

    if(Tdomain%logicD%SF_local_present)then
#if NEW_GLOBAL_METHOD
        !- fluid -> solid coupling (pressure times velocity)
        call FtoS_coupling_2(Tdomain%ngll_s, Tdomain%ngll_f, Tdomain%SF%ngll, &
                             Tdomain%SF%SF_IGlobSol, Tdomain%SF%SF_IGlobFlu, &
                             Tdomain%SF%SF_BtN, Tdomain%champs0%Save_forces, &
                             Tdomain%champs0%Save_depla, Tdomain%champs0%VelPhi, &
                             Tdomain%champs0%Forces, Tdomain%champs0%Depla)
        !- recorrecting on solid faces
        call Newmark_recorrect_solid_2(Tdomain%ngll_s, Tdomain%SF%ngll, Tdomain%TimeD%dtmin, &
                                       Tdomain%SF%SF_IGlobSol, Tdomain%MassMatSol, &
                                       Tdomain%champs0%Forces, Tdomain%champs0%Veloc, &
                                       Tdomain%champs0%Depla)
#else
        !- fluid -> solid coupling (pressure times velocity)
        call FtoS_coupling(Tdomain,rg)
        !- recorrecting on solid faces, edges and vertices
        call Newmark_recorrect_solid(Tdomain)
#endif
    end if

    if (rg==0 .and. mod(ntime,20)==0) print *,' Iteration  =  ',ntime,'    temps  = ',Tdomain%TimeD%rtime

    return

end subroutine Newmark
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
#if ! NEW_GLOBAL_METHOD
subroutine Newmark_Predictor(Tdomain,rg)
    use sdomain
    use assembly
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer, intent(in)  :: rg
    real  :: alpha,bega,gam1,dt
    integer  :: i,n,mat,nf,ne,nv
    logical, dimension(:), allocatable  :: L_Face,L_Edge,L_Vertex

    alpha = Tdomain%TimeD%alpha
    bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
    gam1 = 1. / Tdomain%TimeD%gamma


    do n = 0,Tdomain%n_elem-1
        ! attention au pas de temps a envoyer dans la loi
        ! il faut le pas de temps global
        ! calcule dans courant.f90
        mat = Tdomain%specel(n)%mat_index
        Dt = Tdomain%sSubDomain(mat)%dt
        ! SOLID PART
        if(Tdomain%specel(n)%solid)then
            if (.not. Tdomain%specel(n)%PML) then  ! physical part
                call Prediction_Elem_Veloc(Tdomain%specel(n))
            else    ! PML part
                call get_PMLprediction_v2el(Tdomain,n,bega,dt)
                call get_PMLprediction_e2el(Tdomain,n,bega,dt)
                call get_PMLprediction_f2el(Tdomain,n,bega,dt)
                if (Tdomain%curve) then
                    call Prediction_Elem_PML_Veloc_curve (Tdomain%specel(n),bega, dt, Tdomain%sSubDomain(mat)%hTPrimex, &
                        Tdomain%sSubDomain(mat)%hPrimey, Tdomain%sSubDomain(mat)%hprimez)
                else
                    if (Tdomain%specel(n)%FPML) then
                        call Prediction_Elem_FPML_Veloc(Tdomain%specel(n),bega,dt, &
                            Tdomain%sSubDomain(mat)%hTPrimex,Tdomain%sSubDomain(mat)%hPrimey, &
                            Tdomain%sSubDomain(mat)%hprimez,rg,n,Tdomain%sSubDomain(mat)%freq)
                    else
                        call Prediction_Elem_PML_Veloc(Tdomain%specel(n),bega,dt, &
                            Tdomain%sSubDomain(mat)%hTPrimex,Tdomain%sSubDomain(mat)%hPrimey, &
                            Tdomain%sSubDomain(mat)%hprimez,rg,n)
                    endif
                    !   fin test sur curve
                endif
            endif
            ! FLUID PART
        else
            if(.not. Tdomain%specel(n)%PML)then  ! physical part
                call Prediction_Elem_VelPhi(Tdomain%specel(n),Dt)
            else    ! PML part
                call get_PMLprediction_v2el_fl(Tdomain,n,bega,dt)
                call get_PMLprediction_e2el_fl(Tdomain,n,bega,dt)
                call get_PMLprediction_f2el_fl(Tdomain,n,bega,dt)
                call Prediction_Elem_PML_VelPhi(Tdomain%specel(n),bega,dt,             &
                    Tdomain%sSubDomain(mat)%hTPrimex,Tdomain%sSubDomain(mat)%hPrimey, &
                    Tdomain%sSubDomain(mat)%hprimez)
            endif
        end if
    enddo

    allocate(L_Face(0:Tdomain%n_face-1))
    L_Face = .true.
    allocate(L_Edge(0:Tdomain%n_edge-1))
    L_Edge = .true.
    allocate(L_Vertex(0:Tdomain%n_vertex-1))
    L_Vertex = .true.
    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        Dt = Tdomain%sSubDomain(mat)%dt
        do i = 0,5
            nf = Tdomain%specel(n)%Near_Faces(i)
            if(L_Face(nf))then
                L_Face(nf) = .false.
                if(.not.Tdomain%sface(nf)%PML)then
                    if(Tdomain%sFace(nf)%solid)then
                        call Prediction_Face_Veloc(Tdomain%sface(nf),dt)
                    else
                        call Prediction_Face_VelPhi(Tdomain%sface(nf),dt)
                    end if
                end if
            endif
        enddo
        do i = 0,11
            ne = Tdomain%specel(n)%Near_Edges(i)
            if(L_Edge(ne))then
                L_Edge(ne) = .false.
                if(.not.Tdomain%sedge(ne)%PML)then
                    if(Tdomain%sEdge(ne)%solid)then
                        call Prediction_Edge_Veloc(Tdomain%sedge(ne),dt)
                    else
                        call Prediction_Edge_VelPhi(Tdomain%sedge(ne),dt)
                    end if
                end if
            endif
        enddo
        do i = 0,7
            nv = Tdomain%specel(n)%Near_Vertices(i)
            if(L_Vertex(nv))then
                L_Vertex(nv) = .false.
                if(.not.Tdomain%svertex(nv)%PML)then
                    if(Tdomain%sVertex(nv)%solid)then
                        call Prediction_Vertex_Veloc(Tdomain%svertex(nv),dt)
                    else
                        call Prediction_Vertex_VelPhi(Tdomain%svertex(nv),dt)
                    end if
                end if
            endif
        enddo
    enddo
    deallocate(L_Face,L_Edge,L_Vertex)

    return

end subroutine Newmark_Predictor
#endif
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
#if NEW_GLOBAL_METHOD
subroutine Newmark_Predictor2(Tdomain,champs1)
    use schamps
    use sdomain
    implicit none

    type(domain), intent(inout)   :: Tdomain
    type(champs), intent(inout) :: champs1
    integer :: nel, mat, n, indsol, indpml
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
            ! Couplage à l'interface solide / PML
            indsol = Tdomain%InterfSolPml(n,0)
            indpml = Tdomain%InterfSolPml(n,1)
            Tdomain%champs0%VelocPML(indpml,:) = Tdomain%champs0%Veloc(indsol,:)
            Tdomain%champs0%VelocPML(indpml+1,:) = 0.
            Tdomain%champs0%VelocPML(indpml+2,:) = 0.
        enddo

        do n = 0, Tdomain%nbOuterPMLNodes-1
            indpml = Tdomain%OuterPMLNodes(n)
            Tdomain%champs0%VelocPML(indpml,:) = 0.
            Tdomain%champs0%VelocPML(indpml+1,:) = 0.
            Tdomain%champs0%VelocPML(indpml+2,:) = 0.
        enddo

        do nel = 0,Tdomain%n_elem-1
            if (Tdomain%specel(nel)%PML) then
                mat = Tdomain%specel(nel)%mat_index
                if (Tdomain%curve) then
                    ! TODO
                else
                    if (Tdomain%specel(nel)%FPML) then
                        ! TODO
                    else
                        call Prediction_Elem_PML_Veloc_2(Tdomain%specel(nel), bega, dt, &
                                                         Tdomain%sSubDomain(mat)%hTPrimex, &
                                                         Tdomain%sSubDomain(mat)%hPrimey, &
                                                         Tdomain%sSubDomain(mat)%hprimez, &
                                                         Tdomain%ngll_pmls, &
                                                         Tdomain%champs0%VelocPML, &
                                                         Tdomain%champs1%ForcesPML)
                    endif
                endif
            endif
        enddo
    endif

    return

end subroutine Newmark_Predictor2
#endif
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
#if ! NEW_GLOBAL_METHOD
subroutine Newmark_Corrector(Tdomain,rg)
    use sdomain
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer, intent(in)  :: rg
    real  :: dt
    integer  :: i,n,mat,nf,ne,nv
    logical, dimension(:), allocatable  :: L_Face,L_Edge,L_Vertex

    allocate(L_Face(0:Tdomain%n_face-1))
    L_Face = .true.
    allocate(L_Edge(0:Tdomain%n_edge-1))
    L_Edge = .true.
    allocate(L_Vertex(0:Tdomain%n_vertex-1))
    L_Vertex = .true.


    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        dt = Tdomain%sSubDomain(mat)%Dt
        if(Tdomain%specel(n)%solid)then   ! solid part
            ! inside element
            if(.not. Tdomain%specel(n)%PML)then
                call Correction_Elem_Veloc(Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt)
            else
                if(Tdomain%specel(n)%FPML)then
                    call Correction_Elem_FPML_Veloc(Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt,Tdomain%sSubdomain(mat)%freq)
                else
                    call Correction_Elem_PML_Veloc(Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt)
                endif
            endif
            ! faces
            do i = 0,5
                nf = Tdomain%specel(n)%Near_Faces(i)
                if (L_Face(nf)) then
                    L_Face(nf) = .false.
                    if (.not.Tdomain%sface(nf)%PML) then
                        call Correction_Face_Veloc (Tdomain%sface(nf), dt)
                    else
                        if (Tdomain%sface(nf)%FPML) then
                            call Correction_Face_FPML_Veloc (Tdomain%sface(nf), dt, Tdomain%sSubdomain(mat)%freq)
                        else
                            call Correction_Face_PML_Veloc (Tdomain%sface(nf), dt)
                        endif
                    endif
                endif
            enddo
            ! edges
            do i = 0,11
                ne = Tdomain%specel(n)%Near_Edges(i)
                if (L_Edge(ne)) then
                    L_Edge(ne) = .false.
                    if (.not.Tdomain%sedge(ne)%PML) then
                        call Correction_Edge_Veloc(Tdomain%sedge(ne),dt)
                    else
                        if (Tdomain%sedge(ne)%FPML) then
                            call Correction_Edge_FPML_Veloc(Tdomain%sedge(ne),dt,Tdomain%sSubdomain(mat)%freq)
                        else
                            call Correction_Edge_PML_Veloc(Tdomain%sedge(ne),dt)
                        endif
                    endif
                endif
            enddo
            ! vertices
            do i = 0,7
                nv = Tdomain%specel(n)%Near_Vertices(i)
                if(L_Vertex(nv))then
                    L_Vertex(nv) = .false.
                    if(.not. Tdomain%svertex(nv)%PML)then
                        call Correction_Vertex_Veloc (Tdomain%svertex(nv), dt)
                    else
                        if(Tdomain%svertex(nv)%FPML)then
                            call Correction_Vertex_FPML_Veloc(Tdomain%svertex(nv),dt,Tdomain%sSubdomain(mat)%freq)
                        else
                            call Correction_Vertex_PML_Veloc(Tdomain%svertex(nv),dt)
                        endif
                    endif
                endif
            enddo
        else    ! fluid part
            ! inside element
            if(.not. Tdomain%specel(n)%PML)then
                call Correction_Elem_VelPhi(Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt)
            else
                call Correction_Elem_PML_VelPhi(Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt)
            endif
            ! faces
            do i = 0,5
                nf = Tdomain%specel(n)%Near_Faces(i)
                if(L_Face(nf))then
                    L_Face(nf) = .false.
                    if(.not.Tdomain%sface(nf)%PML)then
                        call Correction_Face_VelPhi(Tdomain%sface(nf),dt)
                    else
                        call Correction_Face_PML_VelPhi(Tdomain%sface(nf),dt)
                    endif
                endif
            enddo
            ! edges
            do i = 0,11
                ne = Tdomain%specel(n)%Near_Edges(i)
                if (L_Edge(ne)) then
                    L_Edge(ne) = .false.
                    if (.not.Tdomain%sedge(ne)%PML) then
                        call Correction_Edge_VelPhi(Tdomain%sedge(ne),dt)
                    else
                        call Correction_Edge_PML_VelPhi(Tdomain%sedge(ne),dt)
                    endif
                endif
            enddo
            ! vertices
            do i = 0,7
                nv = Tdomain%specel(n)%Near_Vertices(i)
                if(L_Vertex(nv))then
                    L_Vertex(nv) = .false.
                    if(.not. Tdomain%svertex(nv)%PML)then
                        call Correction_Vertex_VelPhi(Tdomain%svertex(nv),dt)
                    else
                        call Correction_Vertex_PML_VelPhi(Tdomain%svertex(nv),dt)
                    endif
                endif
            enddo

        end if
    enddo
    deallocate(L_Face,L_Edge,L_Vertex)

    return
end subroutine Newmark_Corrector
#endif
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
#if NEW_GLOBAL_METHOD
subroutine Newmark_Corrector2(Tdomain,champs1)
    use sdomain
    use schamps
    implicit none

    type(domain), intent(inout)   :: Tdomain
    type(champs), intent(in)   :: champs1
    integer  :: n, i_dir

    ! Si il existe des éléments PML solide
    if (Tdomain%ngll_pmls /= 0) then
        do i_dir = 0,2
            Tdomain%champs0%VelocPML(:,i_dir) = Tdomain%champs0%DumpV(:,0) * &
                                                Tdomain%champs0%VelocPML(:,i_dir) + &
                                                Tdomain%TimeD%dtmin * &
                                                Tdomain%champs0%DumpV(:,1) * &
                                                champs1%ForcesPML(:,i_dir)
        enddo
        !TODO Eventuellement : DeplaPML(:,:) = DeplaPML(:,:) + dt * VelocPML(:,:)
    endif

    ! Si il existe des éléments solide
    if (Tdomain%ngll_s /= 0) then
        do n = 0,Tdomain%ngll_s-1
            do i_dir = 0,2
                Tdomain%champs0%Forces(n,i_dir) = champs1%Forces(n,i_dir) * Tdomain%MassMatSol(n)
            enddo
        enddo
        Tdomain%champs0%Veloc = Tdomain%champs0%Veloc + Tdomain%TimeD%dtmin * Tdomain%champs0%Forces
        Tdomain%champs0%Depla = Tdomain%champs0%Depla + Tdomain%TimeD%dtmin * Tdomain%champs0%Veloc
    endif

    ! Si il existe des éléments fluide
    if (Tdomain%ngll_f /= 0) then
        Tdomain%champs0%ForcesFl = champs1%ForcesFl * Tdomain%MassMatFlu
        Tdomain%champs0%VelPhi = (Tdomain%champs0%VelPhi + Tdomain%TimeD%dtmin * Tdomain%champs0%ForcesFl) * Tdomain%champs0%Fluid_dirich
        Tdomain%champs0%Phi = Tdomain%champs0%Phi + Tdomain%TimeD%dtmin * Tdomain%champs0%VelPhi
    endif

    return
end subroutine Newmark_Corrector2
#endif
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
#if ! NEW_GLOBAL_METHOD
subroutine internal_forces(Tdomain,rank)
    ! volume forces - depending on rheology
    use sdomain
    use forces_aniso
    use assembly
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: rank
    integer  :: n,mat


    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        if(.not. Tdomain%specel(n)%PML)then
            if(Tdomain%specel(n)%solid)then
                call get_Displ_Face2Elem(Tdomain,n)
                call get_Displ_Edge2Elem(Tdomain,n)
                call get_Displ_Vertex2Elem(Tdomain,n)
            else
                call get_Phi_Face2Elem(Tdomain,n)
                call get_Phi_Edge2Elem(Tdomain,n)
                call get_Phi_Vertex2Elem(Tdomain,n)
            end if
            call forces_int(Tdomain%specel(n), Tdomain%sSubDomain(mat),            &
                Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                Tdomain%sSubDomain(mat)%hTprimez, Tdomain%n_sls,Tdomain%aniso,     &
                Tdomain%specel(n)%solid)
           
        else   ! PML part
            if(Tdomain%specel(n)%solid)then
                call compute_InternalForces_PML_Elem(Tdomain%specel(n),                &
                    Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimey, &
                    Tdomain%sSubDomain(mat)%hTprimez)
            else
                call compute_InternalForces_PML_Elem_Fl(Tdomain%specel(n),            &
                    Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimey, &
                    Tdomain%sSubDomain(mat)%hTprimez)
            end if
        end if
    enddo

    return
end subroutine internal_forces
#endif
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
#if NEW_GLOBAL_METHOD
subroutine internal_forces2(Tdomain,champs1)
    ! volume forces - depending on rheology
    use sdomain
    use schamps
    use forces_aniso
    implicit none

    type(domain), intent(inout)  :: Tdomain
    type(champs), intent(inout) :: champs1
    integer  :: n,mat, indsol, indpml

    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        if(.not. Tdomain%specel(n)%PML)then
            call forces_int2(Tdomain%specel(n), Tdomain%sSubDomain(mat),           &
                Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                Tdomain%sSubDomain(mat)%hTprimez, Tdomain%n_sls,Tdomain%aniso,     &
                Tdomain%specel(n)%solid, champs1)
        else ! PML
            if (Tdomain%specel(n)%solid)then
                call compute_InternalForces_PML_Elem_2(Tdomain%specel(n),                &
                    Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimey, &
                    Tdomain%sSubDomain(mat)%hTprimez, Tdomain%ngll_pmls, champs1%ForcesPML)
            else
                ! TODO
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

    return
end subroutine internal_forces2
#endif
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
#if ! NEW_GLOBAL_METHOD
subroutine external_forces(Tdomain,timer,ntime,rank)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)  :: rank, ntime
    real, intent(in)  :: timer
    integer  :: ns,nel,i_dir
    real :: t, ft

    do ns = 0, Tdomain%n_source-1
        if(rank == Tdomain%sSource(ns)%proc)then
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
                    Tdomain%specel(nel)%sl%Forces(:,:,:,i_dir) = Tdomain%specel(nel)%sl%Forces(:,:,:,i_dir)+ &
                        ft*Tdomain%sSource(ns)%ExtForce(:,:,:,i_dir)
                end do
            else if(Tdomain%sSource(ns)%i_type_source == 3)then    ! pressure pulse in fluid
                Tdomain%specel(nel)%fl%ForcesFl(:,:,:) = Tdomain%specel(nel)%fl%ForcesFl(:,:,:)+    &
                    ft*Tdomain%sSource(ns)%ExtForce(:,:,:,0)
            end if
        endif
    enddo

    return
end subroutine external_forces
#endif
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
#if NEW_GLOBAL_METHOD
subroutine external_forces2(Tdomain,timer,ntime,rank,champs1)
    use sdomain
    use schamps
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)  :: rank, ntime
    real, intent(in)  :: timer
    type(champs), intent(inout)  :: champs1
    integer  :: ns,nel,i_dir, i,j,k
    real :: t, ft

    do ns = 0, Tdomain%n_source-1
        if(rank == Tdomain%sSource(ns)%proc)then
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
                                champs1%Forces(Tdomain%specel(nel)%Isol(i,j,k),i_dir) = champs1%Forces(Tdomain%specel(nel)%Isol(i,j,k),i_dir) + &
                                    ft*Tdomain%sSource(ns)%ExtForce(i,j,k,i_dir)
                            enddo
                        enddo
                    enddo
                enddo
            else if(Tdomain%sSource(ns)%i_type_source == 3)then    ! pressure pulse in fluid
                do k = 0,Tdomain%specel(nel)%ngllz-1
                    do j = 0,Tdomain%specel(nel)%nglly-1
                        do i = 0,Tdomain%specel(nel)%ngllx-1
                            champs1%ForcesFl(Tdomain%specel(nel)%Iflu(i,j,k)) = champs1%ForcesFl(Tdomain%specel(nel)%Iflu(i,j,k)) +    &
                                ft*Tdomain%sSource(ns)%ExtForce(i,j,k,0)
                        enddo
                    enddo
                enddo
            end if
        endif
    enddo

    return
end subroutine external_forces2
#endif
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
#if ! NEW_GLOBAL_METHOD
subroutine inside_proc_forces(Tdomain)
    use sdomain
    use assembly
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer  :: n,nf,ne,nv

    ! Communication of Forces

    do nf = 0,Tdomain%n_face-1
        if(Tdomain%sFace(nf)%solid)then
            Tdomain%sFace(nf)%Forces = 0
            if(Tdomain%sFace(nf)%PML)then
                Tdomain%sFace(nf)%spml%Forces1 = 0
                Tdomain%sFace(nf)%spml%Forces2 = 0
                Tdomain%sFace(nf)%spml%Forces3 = 0
            endif
        else
            Tdomain%sFace(nf)%ForcesFl = 0
            if(Tdomain%sFace(nf)%PML)then
                Tdomain%sFace(nf)%spml%ForcesFl1 = 0
                Tdomain%sFace(nf)%spml%ForcesFl2 = 0
                Tdomain%sFace(nf)%spml%ForcesFl3 = 0
            endif
        end if
    enddo
    do ne = 0,Tdomain%n_edge-1
        if(Tdomain%sEdge(ne)%solid)then
            Tdomain%sEdge(ne)%Forces = 0
            if(Tdomain%sEdge(ne)%PML)then
                Tdomain%sEdge(ne)%spml%Forces1 = 0
                Tdomain%sEdge(ne)%spml%Forces2 = 0
                Tdomain%sEdge(ne)%spml%Forces3 = 0
            endif
        else
            Tdomain%sEdge(ne)%ForcesFl = 0
            if(Tdomain%sEdge(ne)%PML)then
                Tdomain%sEdge(ne)%spml%ForcesFl1 = 0
                Tdomain%sEdge(ne)%spml%ForcesFl2 = 0
                Tdomain%sEdge(ne)%spml%ForcesFl3 = 0
            end if
        end if
    enddo
    do nv = 0,Tdomain%n_vertex-1
        if(Tdomain%sVertex(nv)%solid)then
            Tdomain%sVertex(nv)%Forces = 0
            if(Tdomain%sVertex(nv)%PML)then
                Tdomain%sVertex(nv)%spml%Forces1 = 0
                Tdomain%sVertex(nv)%spml%Forces2 = 0
                Tdomain%sVertex(nv)%spml%Forces3 = 0
            endif
        else
            Tdomain%sVertex(nv)%ForcesFl = 0
            if(Tdomain%sVertex(nv)%PML)then
                Tdomain%sVertex(nv)%spml%ForcesFl1 = 0
                Tdomain%sVertex(nv)%spml%ForcesFl2 = 0
                Tdomain%sVertex(nv)%spml%ForcesFl3 = 0
            endif

        end if
    enddo


    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%solid)then
            call get_Forces_Elem2Face(Tdomain,n)
            call get_Forces_Elem2Edge(Tdomain,n)
            call get_Forces_Elem2Vertex(Tdomain,n)
        else
            call get_ForcesFl_Elem2Face(Tdomain,n)
            call get_ForcesFl_Elem2Edge(Tdomain,n)
            call get_ForcesFl_Elem2Vertex(Tdomain,n)
        end if
    enddo

    return
end subroutine inside_proc_forces
#endif
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
#if NEW_GLOBAL_METHOD
subroutine Comm_Forces_Complete_2(n,Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: n
    integer :: ngll, ngll_F, idx, i, j, k, nf, ne, nv, ngllPML, ngllPML_F

    ngll = 0 ; ngll_F = 0
    ngllPML = 0 ; ngllPML_F = 0

    ! Faces
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        if(Tdomain%sFace(nf)%solid)then
            if(Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        idx = Tdomain%sFace(nf)%Renum(k,j)
                        Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%champs1%ForcesPML(idx,0:2)
                        Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%champs1%ForcesPML(idx+1,0:2)
                        Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%champs1%ForcesPML(idx+2,0:2)
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            else
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        idx = Tdomain%sFace(nf)%Renum(k,j)
                        Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%champs1%Forces(idx,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
            endif
        else
            if(Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        idx = Tdomain%sFace(nf)%Renum(k,j)
! TODO
                        stop "Newmark Comm_Forces_PML_Complete : PML fluide"
!                         Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,1) = Tdomain%champs1%ForcesFLPML(idx)
!                         Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,2) = Tdomain%champs1%ForcesFLPML(idx+1)
!                         Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,3) = Tdomain%champs1%ForcesFLPML(idx+2)
                        ngllPML_F = ngllPML_F + 1                       
                    enddo                                               
                enddo
            else
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        idx = Tdomain%sFace(nf)%Renum(k,j)
                        Tdomain%sComm(n)%GiveForcesFl(ngll_F) = Tdomain%champs1%ForcesFl(idx)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
            endif
        end if
    enddo
    ! Edges
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        if(Tdomain%sEdge(ne)%solid)then
            if(Tdomain%sEdge(ne)%PML)then
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    idx = Tdomain%sEdge(ne)%Renum(j)
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%champs1%ForcesPML(idx,0:2)
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%champs1%ForcesPML(idx+1,0:2)
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%champs1%ForcesPML(idx+2,0:2)
                    ngllPML = ngllPML + 1
                enddo
            else
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    idx = Tdomain%sEdge(ne)%Renum(j)
                    Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%champs1%Forces(idx,0:2)
                    ngll = ngll + 1
                enddo
            endif
        else
            if(Tdomain%sEdge(ne)%PML)then
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    idx = Tdomain%sEdge(ne)%Renum(j)
! TODO
                    stop "Newmark Comm_Forces_PML_Complete : PML fluide"
!                     Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,1) = Tdomain%champs1%ForcesFLPML(idx)
!                     Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,2) = Tdomain%champs1%ForcesFLPML(idx+1)
!                     Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,3) = Tdomain%champs1%ForcesFLPML(idx+2)
                    ngllPML_F = ngllPML_F + 1
                enddo
            else
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    idx = Tdomain%sEdge(ne)%Renum(j)
                    Tdomain%sComm(n)%GiveForcesFl(ngll_F) = Tdomain%champs1%ForcesFl(idx)
                    ngll_F = ngll_F + 1
                enddo
            endif
        end if
    enddo
    ! Vertices
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        if(Tdomain%sVertex(nv)%solid)then
            if(Tdomain%sVertex(nv)%PML) then
                idx = Tdomain%sVertex(nv)%Renum
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%champs1%ForcesPML(idx,0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%champs1%ForcesPML(idx+1,0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%champs1%ForcesPML(idx+2,0:2)
                ngllPML = ngllPML + 1
            else
                idx = Tdomain%sVertex(nv)%Renum
                Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%champs1%Forces(idx,0:2)
                ngll = ngll + 1
            endif
        else
            if(Tdomain%sVertex(nv)%PML) then
                 idx = Tdomain%sVertex(nv)%Renum
! TODO
                    stop "Newmark Comm_Forces_PML_Complete : PML fluide"
!                 Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,1) = Tdomain%champs1%ForcesFLPML(idx)
!                 Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,2) = Tdomain%champs1%ForcesFLPML(idx+1)
!                 Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,3) = Tdomain%champs1%ForcesFLPML(idx+2)
                ngllPML_F = ngllPML_F + 1
            else
                idx = Tdomain%sVertex(nv)%Renum
                Tdomain%sComm(n)%GiveForcesFl(ngll_F) = Tdomain%champs1%ForcesFl(idx)
                ngll_F = ngll_F + 1
            endif
        end if
    enddo
end subroutine Comm_Forces_Complete_2
#endif
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
#if ! NEW_GLOBAL_METHOD
subroutine Comm_Forces_Complete(n,Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: n
    integer  :: ngll,ngll_F,i,j,k,nf,ne,nv

    ngll = 0 ; ngll_F = 0
    ! faces
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        if(Tdomain%sFace(nf)%solid)then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sFace(nf)%Forces(k,j,0:2)
                    ngll = ngll + 1
                enddo
            enddo
        else
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GiveForcesFl(ngll_F) = Tdomain%sFace(nf)%ForcesFl(k,j)
                    ngll_F = ngll_F + 1
                enddo
            enddo
        end if
    enddo
    ! edges
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        if(Tdomain%sEdge(ne)%solid)then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2)
                ngll = ngll + 1
            enddo
        else
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GiveForcesFl(ngll_F) = Tdomain%sEdge(ne)%ForcesFl(j)
                ngll_F = ngll_F + 1
            enddo
        end if
    enddo
    ! vertices
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        if(Tdomain%sVertex(nv)%solid)then
            Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sVertex(nv)%Forces(0:2)
            ngll = ngll + 1
        else
            Tdomain%sComm(n)%GiveForcesFl(ngll_F) = Tdomain%sVertex(nv)%ForcesFl
            ngll_F = ngll_F + 1
        end if
    enddo

    return
end subroutine Comm_Forces_Complete
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine Comm_Forces_PML_Complete(n,Tdomain)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: n
    integer  :: ngllPML,ngllPML_F,i,j,k,nf,ne,nv

    ngllPML = 0 ; ngllPML_F = 0

    ! faces
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        if(Tdomain%sFace(nf)%PML)then
            if(Tdomain%sFace(nf)%solid)then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sFace(nf)%spml%Forces1(k,j,0:2)
                        Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sFace(nf)%spml%Forces2(k,j,0:2)
                        Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sFace(nf)%spml%Forces3(k,j,0:2)
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            else
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,1) = Tdomain%sFace(nf)%spml%ForcesFl1(k,j)
                        Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,2) = Tdomain%sFace(nf)%spml%ForcesFl2(k,j)
                        Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,3) = Tdomain%sFace(nf)%spml%ForcesFl3(k,j)
                        ngllPML_F = ngllPML_F + 1
                    enddo
                enddo
            end if
        endif
    enddo
    ! edges
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        if(Tdomain%sEdge(ne)%PML)then
            if(Tdomain%sEdge(ne)%solid)then
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sEdge(ne)%spml%Forces1(j,0:2)
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sEdge(ne)%spml%Forces2(j,0:2)
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sEdge(ne)%spml%Forces3(j,0:2)
                    ngllPML = ngllPML + 1
                enddo
            else
                do j = 1,Tdomain%sEdge(ne)%ngll-2
                    Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,1) = Tdomain%sEdge(ne)%spml%ForcesFl1(j)
                    Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,2) = Tdomain%sEdge(ne)%spml%ForcesFl2(j)
                    Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,3) = Tdomain%sEdge(ne)%spml%ForcesFl3(j)
                    ngllPML_F = ngllPML_F + 1
                enddo
            end if
        endif
    enddo
    ! vertices
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        if(Tdomain%sVertex(nv)%PML)then
            if(Tdomain%sVertex(nv)%solid)then
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sVertex(nv)%spml%Forces1(0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sVertex(nv)%spml%Forces2(0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sVertex(nv)%spml%Forces3(0:2)
                ngllPML = ngllPML + 1
            else
                Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,1) = Tdomain%sVertex(nv)%spml%ForcesFl1
                Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,2) = Tdomain%sVertex(nv)%spml%ForcesFl2
                Tdomain%sComm(n)%GiveForcesPMLFl(ngllPML_F,3) = Tdomain%sVertex(nv)%spml%ForcesFl3
                ngllPML_F = ngllPML_F + 1
            end if
        endif
    enddo

    return
end subroutine Comm_Forces_PML_Complete
#endif
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!! subroutine GetProp_Elem(Tdomain,n,nr,rank)
!!     use sdomain
!!     implicit none
!! 
!!     type(domain), intent(inout)  :: Tdomain
!!     integer, intent(in)  :: n,nr,rank
!!     integer  :: nf,ne,nv,orient_f,orient_e,ngllx,nglly,ngllz,   &
!!         ngll1,ngll2,ngll,nnf,nne,nnv,i,j,k
!! 
!!     ngllx = Tdomain%specel(n)%ngllx
!!     nglly = Tdomain%specel(n)%nglly
!!     ngllz = Tdomain%specel(n)%ngllz
!! 
!!     if(Tdomain%specel(n)%solid)then
!!         do k =1,ngllz-2
!!             do j =1,nglly-2
!!                 do i =1,ngllx-2
!!                     Tdomain%sReceiver(nr)%coeff(i,j,k,:) = Tdomain%specel(n)%Veloc(i,j,k,:)
!!                 end do
!!             end do
!!         end do
!!         do nf = 0,5
!!             nnf = Tdomain%specel(n)%Near_Faces(nf)
!!             orient_f = Tdomain%specel(n)%Orient_Faces(nf)
!!             ngll1 = Tdomain%sFace(nnf)%ngll1
!!             ngll2 = Tdomain%sFace(nnf)%ngll2
!!             call get_VectProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
!!                 Tdomain%sFace(nnf)%Veloc(:,:,:),Tdomain%sReceiver(nr)%coeff(:,:,:,:))
!!         enddo
!!         do ne = 0,11
!!             nne = Tdomain%specel(n)%Near_Edges(ne)
!!             orient_e = Tdomain%specel(n)%Orient_Edges(ne)
!!             ngll = Tdomain%sEdge(nne)%ngll
!!             call get_VectProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,  &
!!                 Tdomain%sEdge(nne)%Veloc(:,:),Tdomain%sReceiver(nr)%coeff(:,:,:,:))
!!         end do
!!         do nv = 0,7
!!             nnv = Tdomain%specel(n)%Near_Vertices(nv)
!!             call get_VectProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,rank,  &
!!                 Tdomain%sVertex(nnv)%Veloc(:),Tdomain%sReceiver(nr)%coeff(:,:,:,:))
!!         enddo
!! 
!!     else  ! liquid
!!         Tdomain%sReceiver(nr)%coeff_fl(1:ngllx-2,1:nglly-2,1:ngllz-2) =    &
!!             Tdomain%specel(n)%VelPhi(:,:,:)
!!         do nf = 0,5
!!             nnf = Tdomain%specel(n)%Near_Faces(nf)
!!             orient_f = Tdomain%specel(n)%Orient_Faces(nf)
!!             ngll1 = Tdomain%sFace(nnf)%ngll1
!!             ngll2 = Tdomain%sFace(nnf)%ngll2
!!             call get_ScalarProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
!!                 Tdomain%sFace(nnf)%VelPhi(:,:),Tdomain%sReceiver(nr)%coeff_fl(:,:,:))
!!         enddo
!!         do ne = 0,11
!!             nne = Tdomain%specel(n)%Near_Edges(ne)
!!             orient_e = Tdomain%specel(n)%Orient_Edges(ne)
!!             ngll = Tdomain%sEdge(nne)%ngll
!!             call get_ScalarProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,  &
!!                 Tdomain%sEdge(nne)%VelPhi(:),Tdomain%sReceiver(nr)%coeff_fl(:,:,:))
!!         end do
!!         do nv = 0,7
!!             nnv = Tdomain%specel(n)%Near_Vertices(nv)
!!             call get_ScalarProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,rank,  &
!!                 Tdomain%sVertex(nnv)%VelPhi,Tdomain%sReceiver(nr)%coeff_fl(:,:,:))
!!         enddo
!! 
!!     end if
!! 
!!     return
!! end subroutine GetProp_Elem
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
