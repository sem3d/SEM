!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Newmark.f90
!!\brief Algorithme de Newmark
!! La routine Newmark assure la résolution des équations via un algorithme de predicteur-multi-correcteur
!! des vitesses avec une formulation contrainte-vitesse décalée en temps dans les PML.
!<

#include "index.h"

subroutine Newmark(Tdomain,ntime)
    ! Predictor-MultiCorrector Newmark Velocity Scheme within a
    ! Time staggered Stress-Velocity formulation inside PML
    use sdomain
#ifdef COUPLAGE
    use sCouplage
#endif
    use mcapteur
    use mpi
    use scomm, only : exchange_sem_var, comm_give_data, comm_take_data
    use scommutils
    use stat, only : stat_starttick, stat_stoptick
    use sf_coupling
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: ntime
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
    call Newmark_Predictor(Tdomain)

    !- Solution phase
    call internal_forces(Tdomain)


    ! External Forces
    if(Tdomain%logicD%any_source)then
        call stat_starttick()
        call external_forces(Tdomain,Tdomain%TimeD%rtime,ntime)
        call stat_stoptick(STAT_FEXT)
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
        call StoF_coupling(Tdomain)
    end if

    !- correction phase
    call Newmark_Corrector_Fluid(Tdomain)

    if(Tdomain%logicD%SF_local_present)then
        !- fluid -> solid coupling (pressure times velocity)
        call FtoS_coupling(Tdomain)
    end if
    call Newmark_Corrector_Solid(Tdomain)

    if (Tdomain%rank==0 .and. mod(ntime,20)==0) print *,' Iteration  =  ',ntime,'    temps  = ',Tdomain%TimeD%rtime

    return

end subroutine Newmark
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine comm_forces(Tdomain)
    use sdomain
    use scomm
    use stat, only : stat_starttick, stat_stoptick
    implicit none

    type(domain), intent(inout)   :: Tdomain

    integer :: n, k

    if(Tdomain%Comm_data%ncomm > 0)then
        call stat_starttick()
        do n = 0,Tdomain%Comm_data%ncomm-1
            ! Domain SOLID
            k = 0
            call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                Tdomain%Comm_data%Data(n)%IGiveS, Tdomain%sdom%champs1%Forces, k)

            ! Domain SOLID PML
            if (Tdomain%Comm_data%Data(n)%nsolpml>0) then
                call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                    Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%champs1%ForcesPML, k)
            end if

            ! Domain FLUID
            call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                Tdomain%Comm_data%Data(n)%IGiveF, Tdomain%fdom%champs1%ForcesFl, k)

            ! Domain FLUID PML
            if (Tdomain%Comm_data%Data(n)%nflupml>0) then
                call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                    Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%champs1%fpml_Forces, k)
            end if
            Tdomain%Comm_data%Data(n)%nsend = k
        end do
        call stat_stoptick(STAT_GIVE)

        ! Exchange
        call exchange_sem_var(Tdomain, 104, Tdomain%Comm_data)

        ! Take
        call stat_starttick()
        do n = 0,Tdomain%Comm_data%ncomm-1
            ! Domain SOLID
            k = 0
            call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                Tdomain%Comm_data%Data(n)%IGiveS, Tdomain%sdom%champs1%Forces, k)

            ! Domain SOLID PML
            if (Tdomain%Comm_data%Data(n)%nsolpml>0) then
                call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                    Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%champs1%ForcesPML, k)
            end if

            ! Domain FLUID
            call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                Tdomain%Comm_data%Data(n)%IGiveF, Tdomain%fdom%champs1%ForcesFl, k)

            ! Domain FLUID PML
            if (Tdomain%Comm_data%Data(n)%nflupml>0) then
                call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                    Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%champs1%fpml_Forces, k)
            end if
        end do
        call stat_stoptick(STAT_TAKE)
    endif

    return

end subroutine comm_forces

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine Newmark_Predictor(Tdomain)
    use sdomain
    use dom_fluidpml
    use stat, only : stat_starttick, stat_stoptick
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer :: n, indsol, indpml, indflu
    real :: bega, dt

    bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
    dt = Tdomain%TimeD%dtmin

    ! Elements solide
    if (Tdomain%sdom%nglltot /= 0) then
        call stat_starttick()
        Tdomain%sdom%champs1%Depla = Tdomain%sdom%champs0%Depla
        Tdomain%sdom%champs1%Veloc = Tdomain%sdom%champs0%Veloc
        Tdomain%sdom%champs1%Forces = 0d0
        call stat_stoptick(STAT_FSOL)
    endif

    ! Elements fluide
    if (Tdomain%fdom%nglltot /= 0) then
        call stat_starttick()
        Tdomain%fdom%champs1%VelPhi = Tdomain%fdom%champs0%VelPhi
        Tdomain%fdom%champs1%Phi    = Tdomain%fdom%champs0%Phi
        Tdomain%fdom%champs1%ForcesFl = 0d0
        call stat_stoptick(STAT_FFLU)
    endif

    ! Elements solide pml
    if (Tdomain%spmldom%nglltot /= 0) then
        call stat_starttick()
        Tdomain%spmldom%champs1%ForcesPML = 0.
        do n = 0,Tdomain%intSolPml%surf0%nbtot-1
            ! Couplage à l'interface solide / PML
            indsol = Tdomain%intSolPml%surf0%map(n)
            indpml = Tdomain%intSolPml%surf1%map(n)
            Tdomain%spmldom%champs0%VelocPML(indpml,:,0) = Tdomain%sdom%champs0%Veloc(indsol,:)
            Tdomain%spmldom%champs0%VelocPML(indpml,:,1) = 0.
            Tdomain%spmldom%champs0%VelocPML(indpml,:,2) = 0.
        enddo
        ! Prediction
        Tdomain%spmldom%champs1%VelocPML = Tdomain%spmldom%champs0%VelocPML + dt*(0.5-bega)*Tdomain%spmldom%champs1%ForcesPML
        call stat_stoptick(STAT_PSOL)
    endif

    ! Elements fluide pml
    if (Tdomain%fpmldom%nglltot /= 0) then
        call stat_starttick()
        Tdomain%fpmldom%champs1%fpml_Forces = 0.
        do n = 0,Tdomain%intFluPml%surf0%nbtot-1
            ! Couplage à l'interface fluide / PML
            indflu = Tdomain%intFluPml%surf0%map(n)
            indpml = Tdomain%intFluPml%surf1%map(n)
            Tdomain%fpmldom%champs0%fpml_VelPhi(indpml,0) = Tdomain%fdom%champs0%VelPhi(indflu)
            Tdomain%fpmldom%champs0%fpml_VelPhi(indpml,1) = 0.
            Tdomain%fpmldom%champs0%fpml_VelPhi(indpml,2) = 0.
        enddo
        ! Prediction
        Tdomain%fpmldom%champs1%fpml_Velphi = Tdomain%fpmldom%champs0%fpml_VelPhi + dt*(0.5-bega)*Tdomain%fpmldom%champs1%fpml_Forces
        call stat_stoptick(STAT_PFLU)
    endif

    return

end subroutine Newmark_Predictor
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine Newmark_Corrector_Fluid(Tdomain)
    use sdomain
    use stat, only : stat_starttick, stat_stoptick
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer  :: n,  indpml
    double precision :: dt

    dt = Tdomain%TimeD%dtmin
    ! Si il existe des éléments PML fluides
    if (Tdomain%fpmldom%nglltot /= 0) then
        call stat_starttick()
        Tdomain%fpmldom%champs0%fpml_VelPhi(:,:) = Tdomain%fpmldom%champs0%fpml_DumpV(:,0,:) * &
                                                   Tdomain%fpmldom%champs0%fpml_VelPhi(:,:) + &
                                                   dt * &
                                                   Tdomain%fpmldom%champs0%fpml_DumpV(:,1,:) * &
                                                   Tdomain%fpmldom%champs1%fpml_Forces(:,:)
        do n = 0, Tdomain%fpmldom%n_dirich-1
            indpml = Tdomain%fpmldom%dirich(n)
                     Tdomain%fpmldom%champs0%fpml_VelPhi(indpml,0) = 0.
                     Tdomain%fpmldom%champs0%fpml_VelPhi(indpml,1) = 0.
                     Tdomain%fpmldom%champs0%fpml_VelPhi(indpml,2) = 0.
        enddo


        Tdomain%fpmldom%champs0%fpml_Phi = Tdomain%fpmldom%champs0%fpml_Phi + &
                                           dt*Tdomain%fpmldom%champs0%fpml_VelPhi
        call stat_stoptick(STAT_PFLU)
    endif
    ! Si il existe des éléments fluides
    if (Tdomain%fdom%nglltot /= 0) then
        call stat_starttick()
        Tdomain%fdom%champs0%ForcesFl = Tdomain%fdom%champs1%ForcesFl * Tdomain%fdom%MassMat
        Tdomain%fdom%champs0%VelPhi = (Tdomain%fdom%champs0%VelPhi + dt * Tdomain%fdom%champs0%ForcesFl)
        do n = 0, Tdomain%fdom%n_dirich-1
            indpml = Tdomain%fdom%dirich(n)
            Tdomain%fdom%champs0%VelPhi(indpml) = 0.
        enddo
        Tdomain%fdom%champs0%Phi = Tdomain%fdom%champs0%Phi + dt * Tdomain%fdom%champs0%VelPhi
        call stat_stoptick(STAT_FFLU)
    endif

    return
end subroutine Newmark_Corrector_Fluid
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine Newmark_Corrector_Solid(Tdomain)
    use sdomain
    use stat, only : stat_starttick, stat_stoptick
    implicit none

    type(domain), intent(inout)   :: Tdomain
    integer  :: n, i_dir, indpml
    double precision :: dt

    dt = Tdomain%TimeD%dtmin
    ! Si il existe des éléments PML solides
    if (Tdomain%spmldom%nglltot /= 0) then
        call stat_starttick()
        do i_dir = 0,2
            Tdomain%spmldom%champs0%VelocPML(:,i_dir,:) = Tdomain%spmldom%champs0%DumpV(:,0,:) * &
                                                Tdomain%spmldom%champs0%VelocPML(:,i_dir,:) + &
                                                dt * &
                                                Tdomain%spmldom%champs0%DumpV(:,1,:) * &
                                                Tdomain%spmldom%champs1%ForcesPML(:,i_dir,:)
        enddo
        !TODO Eventuellement : DeplaPML(:,:) = DeplaPML(:,:) + dt * VelocPML(:,:)
        do n = 0, Tdomain%spmldom%n_dirich-1
            indpml = Tdomain%spmldom%dirich(n)
            Tdomain%spmldom%champs0%VelocPML(indpml,:,:) = 0.
        enddo
        call stat_stoptick(STAT_PSOL)
    endif

    ! Si il existe des éléments solides
    if (Tdomain%sdom%nglltot /= 0) then
        call stat_starttick()
        do i_dir = 0,2
            Tdomain%sdom%champs0%Forces(:,i_dir) = Tdomain%sdom%champs1%Forces(:,i_dir) * Tdomain%sdom%MassMat(:)
        enddo
        Tdomain%sdom%champs0%Veloc = Tdomain%sdom%champs0%Veloc + dt * Tdomain%sdom%champs0%Forces
        do n = 0, Tdomain%sdom%n_dirich-1
            indpml = Tdomain%sdom%dirich(n)
            Tdomain%sdom%champs0%Veloc(indpml,:) = 0.
        enddo
        Tdomain%sdom%champs0%Depla = Tdomain%sdom%champs0%Depla + dt * Tdomain%sdom%champs0%Veloc
        call stat_stoptick(STAT_FSOL)
    endif
    return
end subroutine Newmark_Corrector_Solid
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine internal_forces(Tdomain)
    ! volume forces - depending on rheology
    use sdomain
    use dom_solid
    use dom_solidpml
    use dom_fluid
    use dom_fluidpml
    use stat, only : stat_starttick, stat_stoptick
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer  :: n, indsol, indflu, indpml

    if (Tdomain%fdom%nbelem>0) then
        call stat_starttick()
        do n = 0,Tdomain%fdom%nbelem-1,CHUNK
            call forces_int_fluid(Tdomain%fdom, Tdomain%fdom%champs1, n)
        end do
        call stat_stoptick(STAT_FFLU)
    end if
    if (Tdomain%fpmldom%nbelem>0) then
        call stat_starttick()
        do n = 0,Tdomain%fpmldom%nbelem-1,CHUNK
            call pred_flu_pml(Tdomain%fpmldom, Tdomain%TimeD%dtmin, Tdomain%fpmldom%champs1, n)
            call forces_int_flu_pml(Tdomain%fpmldom, Tdomain%fpmldom%champs1, n)
        end do
        call stat_stoptick(STAT_PFLU)
    end if
    if (Tdomain%sdom%nbelem>0) then
        call stat_starttick()
        do n = 0,Tdomain%sdom%nbelem-1,CHUNK
            call forces_int_solid(Tdomain%sdom, Tdomain%sdom%champs1, n, Tdomain%nl_flag)
        end do
        call stat_stoptick(STAT_FSOL)
    end if
    if (Tdomain%spmldom%nbelem>0) then
        call stat_starttick()
        do n = 0,Tdomain%spmldom%nbelem-1,CHUNK
            call pred_sol_pml(Tdomain%spmldom, Tdomain%TimeD%dtmin, Tdomain%spmldom%champs1, n)
            call forces_int_sol_pml(Tdomain%spmldom, Tdomain%spmldom%champs1, n)
        end do
        call stat_stoptick(STAT_PSOL)
    end if

    ! Couplage interface solide / PML
    if (Tdomain%spmldom%nglltot > 0) then
        do n = 0,Tdomain%intSolPml%surf0%nbtot-1
            indsol = Tdomain%intSolPml%surf0%map(n)
            indpml = Tdomain%intSolPml%surf1%map(n)
            Tdomain%sdom%champs1%Forces(indsol,:) = Tdomain%sdom%champs1%Forces(indsol,:) + &
                                                    Tdomain%spmldom%champs1%ForcesPML(indpml,:,0) + &
                                                    Tdomain%spmldom%champs1%ForcesPML(indpml,:,1) + &
                                                    Tdomain%spmldom%champs1%ForcesPML(indpml,:,2)
        enddo
    endif
    ! Couplage interface fluid / PML
    if (Tdomain%fpmldom%nglltot > 0) then
        do n = 0,Tdomain%intFluPml%surf0%nbtot-1
            indflu = Tdomain%intFluPml%surf0%map(n)
            indpml = Tdomain%intFluPml%surf1%map(n)
            Tdomain%fdom%champs1%ForcesFl(indflu) = Tdomain%fdom%champs1%ForcesFl(indflu) + &
                                                    Tdomain%fpmldom%champs1%fpml_Forces(indpml,0) + &
                                                    Tdomain%fpmldom%champs1%fpml_Forces(indpml,1) + &
                                                    Tdomain%fpmldom%champs1%fpml_Forces(indpml,2)
        enddo
    endif

    return
end subroutine internal_forces
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
subroutine external_forces(Tdomain,timer,ntime)
    use sdomain
    implicit none
#include "index.h"

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)  :: ntime
    real, intent(in)  :: timer
    integer  :: ns,nel,i_dir, i,j,k, idx, lnum,ngll
    real :: t, ft

    do ns = 0, Tdomain%n_source-1
        if(Tdomain%rank == Tdomain%sSource(ns)%proc)then
            nel = Tdomain%Ssource(ns)%elem
            lnum = Tdomain%specel(nel)%lnum
            ngll = 0
            select case (Tdomain%specel(nel)%domain)
                 case (DM_SOLID)
                     ngll = Tdomain%sdom%ngll
                 case (DM_FLUID)
                     ngll = Tdomain%fdom%ngll
                 case (DM_SOLID_PML)
                     ngll = Tdomain%spmldom%ngll
                 case (DM_FLUID_PML)
                     ngll = Tdomain%fpmldom%ngll
            end select

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
                    do k = 0,ngll-1
                        do j = 0,ngll-1
                            do i = 0,ngll-1
                                idx = Tdomain%sdom%Idom_(i,j,k,lnum)
                                Tdomain%sdom%champs1%Forces(idx, i_dir) = Tdomain%sdom%champs1%Forces(idx, i_dir) + &
                                    ft*Tdomain%sSource(ns)%ExtForce(i,j,k,i_dir)
                            enddo
                        enddo
                    enddo
                enddo
            else if(Tdomain%sSource(ns)%i_type_source == 3)then    ! pressure pulse in fluid
                do k = 0,ngll-1
                    do j = 0,ngll-1
                        do i = 0,ngll-1
                            idx = Tdomain%fdom%Idom_(i,j,k,lnum)
                            Tdomain%fdom%champs1%ForcesFl(idx) = Tdomain%fdom%champs1%ForcesFl(idx) +    &
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
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
