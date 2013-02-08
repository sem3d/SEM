!>
!!\file Newmark.f90
!!\brief Algorithme de Newmark
!!\author
!!\version 1.0
!!\date 10/03/2009
!! La routine Newmark assure la r�solution des �quations via un algorithme de predicteur-multi-correcteur
!! des vitesses avec une formulation contrainte-vitesse d�cal�e en temps dans les PML.
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

    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: rg,ntime

    logical, dimension(:), allocatable :: L_Face, L_Edge, L_Vertex
    integer ::  n, mat,  ngll1,ngll2,ngll3,    &
        nf,ne,nv, nf_aus,ne_aus,nv_aus, ngll, code, &
        i,j, x,y,z,ngllPML,shift, I_give_to, I_take_from,     &
        n_rings, ntimetrace,ngll_F,ngllPML_F
    integer, parameter :: etiquette = 100
    integer, dimension(mpi_status_size) :: statut
    real :: Dt


    ! Predictor-MultiCorrector Newmark Velocity Scheme within a
    ! Time staggered Stress-Velocity formulation inside PML
    ! PML needs to be implemented
    if(.not. Tdomain%TimeD%velocity_scheme)   &
        stop "Newmark scheme implemented only in velocity form."

    !- Prediction Phase
    call Newmark_Predictor(Tdomain,rg)

    !- Solution phase
    call internal_forces(Tdomain,rg)


    ! External Forces
    if(Tdomain%logicD%any_source)then
        call external_forces(Tdomain,Tdomain%TimeD%rtime,ntime,rg)
    end if

    ! Communication of Forces within a single process
    call inside_proc_forces(Tdomain)



#ifdef COUPLAGE
    if (ntime>0) then
        call calcul_couplage_force(Tdomain, ntime, rg)
    endif

    !Gsa Ipsis (tout le passage)
    ! AJOUT DES FORCES MKA3D

    ! la ForcesMka corespond a la contrainte multiplie par la surface du point de gauss correspondant
    ! sur un meme proc la somme est deja faite
    ! par contre lorsqu un vertex est partage sur plusieurs proc alors chaque proc n a qu une partie de la somme
    ! il faut donc lui ajouter les contributions des autres proc
    ! pour prendre en compte les forces imposee lors du couplage avec mka sur les points de gauss internes aux faces
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
            call Comm_Forces_Complete(n,Tdomain)
            call Comm_Forces_PML_Complete(n,Tdomain)
        end do

        call exchange_sem_forces(Tdomain, rg)

        ! now: assemblage on external faces, edges and vertices
        do n = 0,Tdomain%n_proc-1
            ngll = 0
            ngll_F = 0
            ngllPML = 0
            ngllPML_F = 0
            call Comm_Forces_Face  (Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
            call Comm_Forces_Edge  (Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
            call Comm_Forces_Vertex(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)
        enddo

    endif  ! if nproc > 1


    ! Neumann B.C.: associated forces
    if(Tdomain%logicD%neumann_local_present)then
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
    endif

    !- solid -> fluid coupling (normal dot velocity)
    if(Tdomain%logicD%SF_local_present)then
        call SF_solid_values_saving(Tdomain)
        call StoF_coupling(Tdomain,rg)
        dt = Tdomain%sSubdomain(mat)%dt
    end if

    call Newmark_Corrector(Tdomain,rg)
    if(Tdomain%logicD%SF_local_present)then
        !- fluid -> solid coupling (pressure times velocity)
        call FtoS_coupling(Tdomain,rg)
        !- recorrecting on solid faces, edges and vertices
        call Newmark_recorrect_solid(Tdomain)
    end if
    !  modif mariotti fevrier 2007 cea capteur displ
    ! si on veut soritr la vitesse mettre Veloc a la place de Displ
    ! il faut en fait stocker les deux et pouvoir sortir l un ou l autre au choix
    ! dans le fichier de data
    if (Tdomain%logicD%save_trace) then
        call save_traces(Tdomain, ntime, rg)
    end if


    ! Save Trace
    !if(Tdomain%logicD%save_trace) call dumptrace(Tdomain,rg,ntime)

    if (rg==0 .and. mod(ntime,20)==0) print *,' Iteration  =  ',ntime,'    temps  = ',Tdomain%TimeD%rtime

    return

end subroutine Newmark
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine Newmark_Predictor(Tdomain,rg)

    use sdomain
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
                call get_PMLprediction_v2el(Tdomain,n,bega,dt,rg)
                call get_PMLprediction_e2el(Tdomain,n,bega,dt,rg)
                call get_PMLprediction_f2el(Tdomain,n,bega,dt,rg)
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
                call get_PMLprediction_v2el_fl(Tdomain,n,bega,dt,rg)
                call get_PMLprediction_e2el_fl(Tdomain,n,bega,dt,rg)
                call get_PMLprediction_f2el_fl(Tdomain,n,bega,dt,rg)
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
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
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine internal_forces(Tdomain,rank)
    ! volume forces - depending on rheology
    use sdomain
    use forces_aniso
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: rank
    integer  :: n,mat


    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        if(Tdomain%specel(n)%solid)then  ! solid part
            if(.not.Tdomain%specel(n)%PML)then
                call get_Displ_Face2Elem(Tdomain,n)
                call get_Displ_Edge2Elem(Tdomain,n)
                call get_Displ_Vertex2Elem(Tdomain,n)
                if (.true.) then
                    call forces_int (Tdomain%specel(n),&
                        Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                        Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                        Tdomain%sSubDomain(mat)%hTprimez, Tdomain%n_sls, Tdomain%aniso)
                else
                    ! !! Si on veut utiliser ce qui suit a la place de la routine forces_int, alors il faut
                    ! !! decommenter l'allocation de Acoeff dans "allocate_domain.f90" et la definition de
                    ! !! Acoeff dans "define_arrays.f90".
                    ! !! ATTENTION: la routine ci-dessous ne prend en compte ni l'attenuation ni l'anisotropie.
                    call compute_InternalForces_Elem(Tdomain%specel(n), &
                        Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimex, &
                        Tdomain%sSubDomain(mat)%hprimey,Tdomain%sSubDomain(mat)%hTprimey, &
                        Tdomain%sSubDomain(mat)%hprimez,Tdomain%sSubDomain(mat)%hTprimez)
                endif
            else
                call compute_InternalForces_PML_Elem(Tdomain%specel(n),                &
                    Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimey, &
                    Tdomain%sSubDomain(mat)%hTprimez)
            endif
        else   ! fluid part
            if(.not.Tdomain%specel(n)%PML)then
                call get_Phi_Face2Elem(Tdomain,n)
                call get_Phi_Edge2Elem(Tdomain,n)
                call get_Phi_Vertex2Elem(Tdomain,n)
                call compute_InternalForces_Elem_Fluid(Tdomain%specel(n),             &
                    Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimex, &
                    Tdomain%sSubDomain(mat)%hprimey,Tdomain%sSubDomain(mat)%hTprimey, &
                    Tdomain%sSubDomain(mat)%hprimez,Tdomain%sSubDomain(mat)%hTprimez)
            else
                call compute_InternalForces_PML_Elem_Fl(Tdomain%specel(n),            &
                    Tdomain%sSubDomain(mat)%hprimex,Tdomain%sSubDomain(mat)%hTprimey, &
                    Tdomain%sSubDomain(mat)%hTprimez)
            end if

        end if
    enddo

    return
end subroutine internal_forces
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
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

            ! time : t_(n+1/2) for solid ; t_n for fluid
            t = merge(timer+Tdomain%TimeD%dtmin/2d0,timer,Tdomain%specel(nel)%solid)
            !
            if(Tdomain%sSource(ns)%i_type_source == 1)then  ! collocated force in solid
                ft = CompSource(Tdomain%sSource(ns), t, ntime)
                !
                i_dir = Tdomain%Ssource(ns)%i_dir
                Tdomain%specel(nel)%Forces(:,:,:,i_dir) = Tdomain%specel(nel)%Forces(:,:,:,i_dir)+ &
                    ft*Tdomain%sSource(ns)%ExtForce(:,:,:)
            else if(Tdomain%sSource(ns)%i_type_source == 2)then  ! moment tensor source
                ft = CompSource(Tdomain%sSource(ns), t, ntime)
                !
                do i_dir = 0,2
                    Tdomain%specel(nel)%Forces(:,:,:,i_dir) = Tdomain%specel(nel)%Forces(:,:,:,i_dir)+ &
                        ft*Tdomain%sSource(ns)%MomForce(:,:,:,i_dir)
                end do
            else if(Tdomain%sSource(ns)%i_type_source == 3)then    ! pressure pulse in fluid
                Tdomain%specel(nel)%ForcesFl(:,:,:) = Tdomain%specel(nel)%ForcesFl(:,:,:)+    &
                    CompSource_Fl(Tdomain%sSource(ns),t)*Tdomain%sSource(ns)%ExtForce(:,:,:)
            end if
        endif
    enddo

    return
end subroutine external_forces
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine inside_proc_forces(Tdomain)
    use sdomain
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
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
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
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
subroutine dumptrace(Tdomain,rank,ntime)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)   :: rank,ntime
    integer  :: n,nr,ndt2,ngllx,nglly,ngllz,ntimetrace,i,j,k

    do nr = 0, Tdomain%n_receivers-1
        ndt2 = Tdomain%sReceiver(nr)%ndt
        n = Tdomain%sReceiver(nr)%elem
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz
        if(rank == Tdomain%sReceiver(nr)%proc)then
            if(mod(ntime,Tdomain%TimeD%ntrace) == 0)then
                if(Tdomain%sReceiver(nr)%flag == 1)then
                    if(Tdomain%specel(n)%solid)then
                        allocate(Tdomain%sReceiver(nr)%StoreTrace(0:Tdomain%TimeD%ntrace-1,0:2))
                    else
                        allocate(Tdomain%sReceiver(nr)%StoreTrace_Fl(0:Tdomain%TimeD%ntrace-1))
                    end if
                end if
                if(Tdomain%sReceiver(nr)%flag == 2)then
                    if(Tdomain%specel(n)%solid)then
                        allocate(Tdomain%sReceiver(nr)%StoreTrace(0:(Tdomain%TimeD%ntrace-1)/ndt2,0:2))
                    else
                        allocate(Tdomain%sReceiver(nr)%StoreTrace_Fl(0:(Tdomain%TimeD%ntrace-1)/ndt2))
                    end if
                end if
                if(Tdomain%specel(n)%solid) Tdomain%sReceiver(nr)%StoreTrace = 0.
                if(.not. Tdomain%specel(n)%solid) Tdomain%sReceiver(nr)%StoreTrace_Fl = 0.
            endif

            if(Tdomain%sReceiver(nr)%flag == 1 .or. ((Tdomain%sReceiver(nr)%flag == 2) .and.   &
                (mod(ntime+1,ndt2) == 0)))then
                if(Tdomain%sReceiver(nr)%flag == 1)      &
                    ntimetrace = mod(ntime,Tdomain%TimeD%ntrace)
                if(Tdomain%sReceiver(nr)%flag == 2)      &
                    ntimetrace = mod((ntime+1)/ndt2-1,Tdomain%TimeD%ntrace/ndt2)

                call getProp_Elem(Tdomain,n,nr,rank)

                if(Tdomain%specel(n)%solid)then
                    do k = 0,ngllz-1
                        do j = 0,nglly-1
                            do i = 0,ngllx-1
                                Tdomain%sReceiver(nr)%StoreTrace(ntimetrace,:) = Tdomain%sReceiver(nr)%StoreTrace(ntimetrace,:) + &
                                    Tdomain%sReceiver(nr)%coeff(i,j,k,:) * Tdomain%sReceiver(nr)%pol(i,j,k)
                            enddo
                        enddo
                    enddo
                else
                    do k = 0,ngllz-1
                        do j = 0,nglly-1
                            do i = 0,ngllx-1
                                Tdomain%sReceiver(nr)%StoreTrace_Fl(ntimetrace) = Tdomain%sReceiver(nr)%StoreTrace_Fl(ntimetrace) + &
                                    Tdomain%sReceiver(nr)%coeff_fl(i,j,k) * Tdomain%sReceiver(nr)%pol(i,j,k)
                            enddo
                        enddo
                    enddo

                end if
            end if

        endif
    end do

    return

end subroutine dumptrace
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine GetProp_Elem(Tdomain,n,nr,rank)
    use sdomain
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)  :: n,nr,rank
    integer  :: nf,ne,nv,orient_f,orient_e,ngllx,nglly,ngllz,   &
        ngll1,ngll2,ngll,nnf,nne,nnv,i,j,k

    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    if(Tdomain%specel(n)%solid)then
        do k =1,ngllz-2
            do j =1,nglly-2
                do i =1,ngllx-2
                    Tdomain%sReceiver(nr)%coeff(i,j,k,:) = Tdomain%specel(n)%Veloc(i,j,k,:)
                end do
            end do
        end do
        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            call get_VectProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
                Tdomain%sFace(nnf)%Veloc(:,:,:),Tdomain%sReceiver(nr)%coeff(:,:,:,:))
        enddo
        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            call get_VectProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,  &
                Tdomain%sEdge(nne)%Veloc(:,:),Tdomain%sReceiver(nr)%coeff(:,:,:,:))
        end do
        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            call get_VectProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,rank,  &
                Tdomain%sVertex(nnv)%Veloc(:),Tdomain%sReceiver(nr)%coeff(:,:,:,:))
        enddo

    else  ! liquid
        Tdomain%sReceiver(nr)%coeff_fl(1:ngllx-2,1:nglly-2,1:ngllz-2) =    &
            Tdomain%specel(n)%VelPhi(:,:,:)
        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            orient_f = Tdomain%specel(n)%Orient_Faces(nf)
            ngll1 = Tdomain%sFace(nnf)%ngll1
            ngll2 = Tdomain%sFace(nnf)%ngll2
            call get_ScalarProperty_Face2Elem(nf,orient_f,ngllx,nglly,ngllz,ngll1,ngll2,rank,  &
                Tdomain%sFace(nnf)%VelPhi(:,:),Tdomain%sReceiver(nr)%coeff_fl(:,:,:))
        enddo
        do ne = 0,11
            nne = Tdomain%specel(n)%Near_Edges(ne)
            orient_e = Tdomain%specel(n)%Orient_Edges(ne)
            ngll = Tdomain%sEdge(nne)%ngll
            call get_ScalarProperty_Edge2Elem(ne,orient_e,ngllx,nglly,ngllz,ngll,rank,  &
                Tdomain%sEdge(nne)%VelPhi(:),Tdomain%sReceiver(nr)%coeff_fl(:,:,:))
        end do
        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            call get_ScalarProperty_Vertex2Elem(nv,ngllx,nglly,ngllz,rank,  &
                Tdomain%sVertex(nnv)%VelPhi,Tdomain%sReceiver(nr)%coeff_fl(:,:,:))
        enddo

    end if

    return
end subroutine GetProp_Elem
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
