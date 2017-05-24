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

module mtimestep
    use constants
    use sdomain
    implicit none

    real(fpp), dimension(6), parameter :: LDDRK_beta = (/ &
        0.0_fpp, &
        -0.737101392796_fpp, &
        -1.634740794341_fpp, &
        -0.744739003780_fpp, &
        -1.469897351522_fpp, &
        -2.813971388035_fpp /)
    real(fpp), dimension(6), parameter :: LDDRK_gamma = (/ &
        0.032918605146_fpp, &
        0.823256998200_fpp, &
        0.381530948900_fpp, &
        0.200092213184_fpp, &
        1.718581042715_fpp, &
        0.27_fpp /)
    real(fpp), dimension(6), parameter :: LDDRK_c = (/ &
        0.0_fpp, &
        0.032918605146_fpp, &
        0.249351723343_fpp, &
        0.466911705055_fpp, &
        0.582030414044_fpp, &
        0.847252983783_fpp /)
contains

    subroutine Timestep_LDDRK(Tdomain,ntime)
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: ntime
        !
        real(fpp) :: cb, cg, cc, t0, t, dt
        integer :: i0, i1, tmp, i
        i0 = 0
        i1 = 1
        t0 = Tdomain%TimeD%rtime
        dt = Tdomain%TimeD%dtmin
        Tdomain%sdom%champs(i1)%Veloc = Tdomain%sdom%champs(i0)%Veloc
        Tdomain%sdom%champs(i1)%Depla = Tdomain%sdom%champs(i0)%Depla
        Tdomain%sdom%champs(i0)%Forces = 0d0
        do i=1,6
            cb = LDDRK_beta(i)
            cg = LDDRK_gamma(i)
            cc = LDDRK_c(i)
            t = t0 + cc*dt


            call internal_forces(Tdomain, i0, i1, ntime)
            call external_forces(Tdomain,t,ntime, i1)
            call comm_forces(Tdomain,i1)


            ! Only solid for starters
            Tdomain%sdom%champs(i1)%Veloc = cb*Tdomain%sdom%champs(i0)%Veloc + dt*Tdomain%sdom%champs(i1)%Forces
            Tdomain%sdom%champs(i1)%Depla = cb*Tdomain%sdom%champs(i0)%Depla + dt*Tdomain%sdom%champs(i1)%Veloc

            Tdomain%sdom%champs(i1)%Veloc = Tdomain%sdom%champs(i0)%Veloc+Tdomain%sdom%champs(i0)%Veloc
            Tdomain%sdom%champs(i1)%Depla = Tdomain%sdom%champs(i0)%Depla

            ! SWAP i0,i1 for storage
            tmp = i1
            i1 = i0
            i0 = tmp
        end do
    end subroutine Timestep_LDDRK
    !
    subroutine Newmark(Tdomain,ntime)
        ! Predictor-MultiCorrector Newmark Velocity Scheme within a
        ! Time staggered Stress-Velocity formulation inside PML
#ifdef COUPLAGE
        use sCouplage
#endif
        use mcapteur
        use mpi
        use scomm, only : exchange_sem_var, comm_give_data, comm_take_data
        use scommutils
        use stat, only : stat_starttick, stat_stoptick, STAT_FEXT
        use sf_coupling
        use surface_load

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
        call internal_forces(Tdomain, 0, 1, ntime)

        ! External Forces
        if(Tdomain%logicD%any_source)then
            call stat_starttick()
            call external_forces(Tdomain,Tdomain%TimeD%rtime,ntime, 1)
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
        call comm_forces(Tdomain,1)
        ! Neumann B.C.: associated forces
        !    if(Tdomain%logicD%neumann_local_present)then
        !#if 0
        !    !TODO
        !    stop "TODO: conditions de Neumann non prises en compte (Newmark)"
        !#else
        !    ! call add_Newman_forces(Tdomain)
        !#endif
        !    endif

        if (Tdomain%logicD%surfBC) then
            call add_surface_force(Tdomain,0,1)
        endif
        !- solid -> fluid coupling (normal dot velocity)
        call Newmark_Corrector_S(Tdomain)
        if(Tdomain%logicD%SF_local_present)then
            call StoF_coupling(Tdomain,0,1)
        end if

        !- correction phase
        call Newmark_Corrector_F(Tdomain)
        if(Tdomain%logicD%SF_local_present)then
            !- fluid -> solid coupling (pressure times velocity)
            call FtoS_coupling(Tdomain,0,1)
        end if

        if (Tdomain%rank==0 .and. mod(ntime,20)==0) print *,' Iteration  =  ',ntime,'    temps  = ',Tdomain%TimeD%rtime

        return
    end subroutine Newmark
    !---------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------
    subroutine comm_forces(Tdomain, f1)
        use scomm
        use stat, only : stat_starttick, stat_stoptick, STAT_GIVE, STAT_TAKE
        implicit none

        type(domain), intent(inout)   :: Tdomain
        integer, intent(in) :: f1

        integer :: n, k

        if(Tdomain%Comm_data%ncomm > 0)then
            call stat_starttick()
            do n = 0,Tdomain%Comm_data%ncomm-1
                ! Domain SOLID
                k = 0
                call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                    Tdomain%Comm_data%Data(n)%IGiveS, Tdomain%sdom%champs(f1)%Forces, k)

                ! Domain SOLID PML
                if (Tdomain%Comm_data%Data(n)%nsolpml>0) then
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
#ifdef CPML
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%champs(f1)%Forces, k)
#else
                    Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%champs(f1)%ForcesPML, k)
#endif
                end if

                ! Domain FLUID
                call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                    Tdomain%Comm_data%Data(n)%IGiveF, Tdomain%fdom%champs(f1)%ForcesFl, k)

                ! Domain FLUID PML
                if (Tdomain%Comm_data%Data(n)%nflupml>0) then
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
#ifdef CPML
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%champs(f1)%ForcesFl, k)
#else
                    Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%champs(f1)%fpml_Forces, k)
#endif
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
                    Tdomain%Comm_data%Data(n)%IGiveS, Tdomain%sdom%champs(f1)%Forces, k)

                ! Domain SOLID PML
                if (Tdomain%Comm_data%Data(n)%nsolpml>0) then
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
#ifdef CPML
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%champs(f1)%Forces, k)
#else
                    Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%champs(f1)%ForcesPML, k)
#endif
                end if

                ! Domain FLUID
                call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                    Tdomain%Comm_data%Data(n)%IGiveF, Tdomain%fdom%champs(f1)%ForcesFl, k)

                ! Domain FLUID PML
                if (Tdomain%Comm_data%Data(n)%nflupml>0) then
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
#ifdef CPML
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%champs(f1)%ForcesFl, k)
#else
                    Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%champs(f1)%fpml_Forces, k)
#endif
                end if
            end do
            call stat_stoptick(STAT_TAKE)
        endif

        return

    end subroutine comm_forces

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine Newmark_Predictor(Tdomain)
        use dom_fluid
        use dom_fluidpml
        use dom_solid
        use dom_solidpml
        use stat, only : stat_starttick, stat_stoptick, STAT_FSOL, STAT_FFLU, STAT_PSOL, STAT_PFLU
        implicit none

        type(domain), intent(inout)   :: Tdomain

        ! Elements solide
        if (Tdomain%sdom%nglltot /= 0) then
            call stat_starttick()
            call newmark_predictor_solid(Tdomain%sdom,0,1)
            call stat_stoptick(STAT_FSOL)
        endif

        ! Elements fluide
        if (Tdomain%fdom%nglltot /= 0) then
            call stat_starttick()
            call newmark_predictor_fluid(Tdomain%fdom,0,1)
            call stat_stoptick(STAT_FFLU)
        endif

        ! Elements solide pml
        if (Tdomain%spmldom%nglltot /= 0) then
            call stat_starttick()
            call newmark_predictor_solidpml(Tdomain%spmldom, Tdomain,0,1)
            call stat_stoptick(STAT_PSOL)
        endif

        ! Elements fluide pml
        if (Tdomain%fpmldom%nglltot /= 0) then
            call stat_starttick()
            call newmark_predictor_fluidpml(Tdomain%fpmldom, Tdomain,0,1)
            call stat_stoptick(STAT_PFLU)
        endif

        return

    end subroutine Newmark_Predictor
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine Newmark_Corrector_F(Tdomain)
        use dom_fluid
        use dom_fluidpml
        use stat, only : stat_starttick, stat_stoptick, STAT_PFLU, STAT_FFLU
        implicit none

        type(domain), intent(inout)   :: Tdomain
        real(fpp) :: dt

        dt = Tdomain%TimeD%dtmin
        ! Si il existe des éléments PML fluides
        if (Tdomain%fpmldom%nglltot /= 0) then
            call stat_starttick()
            call newmark_corrector_fluidpml(Tdomain%fpmldom, dt, 0, 1)
            call stat_stoptick(STAT_PFLU)
        endif
        ! Si il existe des éléments fluides
        if (Tdomain%fdom%nglltot /= 0) then
            call stat_starttick()
            call newmark_corrector_fluid(Tdomain%fdom, dt, 0, 1)
            call stat_stoptick(STAT_FFLU)
        endif

        return
    end subroutine Newmark_Corrector_F
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    subroutine Newmark_Corrector_S(Tdomain)
        use dom_solid
        use dom_solidpml
        use stat, only : stat_starttick, stat_stoptick, STAT_PSOL, STAT_FSOL
        implicit none

        type(domain), intent(inout)   :: Tdomain
        real(fpp) :: dt, t

        dt = Tdomain%TimeD%dtmin
        t = Tdomain%TimeD%rtime
        ! Si il existe des éléments PML solides
        if (Tdomain%spmldom%nglltot /= 0) then
            call stat_starttick()
            call newmark_corrector_solidpml(Tdomain%spmldom, dt, t, 0, 1)
            call stat_stoptick(STAT_PSOL)
        endif

        ! Si il existe des éléments solides
        if (Tdomain%sdom%nglltot /= 0) then
            call stat_starttick()
            call newmark_corrector_solid(Tdomain%sdom, dt, 0, 1)
            call stat_stoptick(STAT_FSOL)
        endif
        return
    end subroutine Newmark_Corrector_S
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    subroutine internal_forces(Tdomain, i0, i1, ntime)
        ! volume forces - depending on rheology
        use dom_solid
        use dom_solidpml
        use dom_fluid
        use dom_fluidpml
        use smirror
        use stat, only : stat_starttick, stat_stoptick, STAT_FFLU, STAT_PFLU, STAT_FSOL, STAT_PSOL
        implicit none

        type(domain), intent(inout)  :: Tdomain
        integer, intent(in) :: i0, i1, ntime
        integer  :: n, indsol, indflu, indpml

        if (Tdomain%fdom%nbelem>0) then
            call stat_starttick()
            if (Tdomain%fdom%use_mirror.and.Tdomain%fdom%mirror_fl%n_glltot>0.and.Tdomain%fdom%mirror_type>0) then
                call load_mirror_fl(Tdomain%fdom, ntime)
            endif
            do n = 0,Tdomain%fdom%nblocks-1
                call forces_int_fluid(Tdomain%fdom, Tdomain%fdom%champs(i1), n)
            end do
            if (Tdomain%fdom%use_mirror.and.Tdomain%fdom%mirror_fl%n_glltot>0.and.Tdomain%fdom%mirror_type==0) then
                call dump_mirror_fl(Tdomain%fdom, ntime)
            endif
            call stat_stoptick(STAT_FFLU)
        end if
        if (Tdomain%fpmldom%nbelem>0) then
            call stat_starttick()
            do n = 0,Tdomain%fpmldom%nblocks-1
                call pred_flu_pml(Tdomain%fpmldom, Tdomain%TimeD%dtmin, Tdomain%fpmldom%champs(i1), n)
                call forces_int_flu_pml(Tdomain%fpmldom, Tdomain%fpmldom%champs(i1), n, Tdomain)
            end do
            call stat_stoptick(STAT_PFLU)
        end if
        if (Tdomain%sdom%nbelem>0) then
            call stat_starttick()
            if (Tdomain%sdom%use_mirror.and.Tdomain%sdom%mirror_sl%n_glltot>0.and.Tdomain%sdom%mirror_type>0) then
                call load_mirror_sl(Tdomain%sdom, ntime)
            endif
            do n = 0,Tdomain%sdom%nblocks-1
                call forces_int_solid(Tdomain%sdom, Tdomain%sdom%champs(i1), n, Tdomain%nl_flag)
            end do
            if (Tdomain%sdom%use_mirror.and.Tdomain%sdom%mirror_sl%n_glltot>0.and.Tdomain%sdom%mirror_type==0) then
                call dump_mirror_sl(Tdomain%sdom, ntime)
            endif
            call stat_stoptick(STAT_FSOL)
        end if
        if (Tdomain%spmldom%nbelem>0) then
            call stat_starttick()
            do n = 0,Tdomain%spmldom%nblocks-1
                call pred_sol_pml(Tdomain%spmldom, Tdomain%TimeD%dtmin, Tdomain%spmldom%champs(i1), n)
                call forces_int_sol_pml(Tdomain%spmldom, Tdomain%spmldom%champs(i1), n, Tdomain)
            end do
            call stat_stoptick(STAT_PSOL)
        end if

        ! Couplage interface solide / PML
        if (Tdomain%spmldom%nglltot > 0) then
            do n = 0,Tdomain%intSolPml%surf0%nbtot-1
                indsol = Tdomain%intSolPml%surf0%map(n)
                indpml = Tdomain%intSolPml%surf1%map(n)
                Tdomain%sdom%champs(i1)%Forces(indsol,:) = Tdomain%sdom%champs(i1)%Forces(indsol,:) + &
#ifdef CPML
                    Tdomain%spmldom%champs(i1)%Forces(indpml,:) &
                    - Tdomain%spmldom%DumpMat(indpml)*Tdomain%spmldom%champs(i0)%Veloc(indpml,:) &
                    - Tdomain%spmldom%MasUMat(indpml)*Tdomain%spmldom%champs(i0)%Depla(indpml,:)
#else
                Tdomain%spmldom%champs(i1)%ForcesPML(indpml,:,0) + &
                    Tdomain%spmldom%champs(i1)%ForcesPML(indpml,:,1) + &
                    Tdomain%spmldom%champs(i1)%ForcesPML(indpml,:,2)
#endif
            enddo
        endif
        ! Couplage interface fluid / PML
        if (Tdomain%fpmldom%nglltot > 0) then
            do n = 0,Tdomain%intFluPml%surf0%nbtot-1
                indflu = Tdomain%intFluPml%surf0%map(n)
                indpml = Tdomain%intFluPml%surf1%map(n)
                Tdomain%fdom%champs(i1)%ForcesFl(indflu) = Tdomain%fdom%champs(i1)%ForcesFl(indflu) + &
#ifdef CPML
                    Tdomain%fpmldom%champs(i1)%ForcesFl(indpml)
#else
                Tdomain%fpmldom%champs(i1)%fpml_Forces(indpml,0) + &
                    Tdomain%fpmldom%champs(i1)%fpml_Forces(indpml,1) + &
                    Tdomain%fpmldom%champs(i1)%fpml_Forces(indpml,2)
#endif
            enddo
        endif

        return
    end subroutine internal_forces
    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine external_forces(Tdomain,timer,ntime, i1)

        implicit none
#include "index.h"

        type(domain), intent(inout)  :: Tdomain
        integer, intent(in)  :: ntime
        real(kind=fpp), intent(in)  :: timer
        integer, intent(in) :: i1
        integer  :: ns,nel,i_dir, i,j,k, idx, lnum,ngll, bnum, ee
        real(kind=fpp) :: t, ft, val

        do ns = 0, Tdomain%n_source-1
            if(Tdomain%rank == Tdomain%sSource(ns)%proc)then
                nel = Tdomain%Ssource(ns)%elem
                lnum = Tdomain%specel(nel)%lnum
                ngll = domain_ngll(Tdomain, Tdomain%specel(nel)%domain)
                bnum = lnum/VCHUNK
                ee = mod(lnum,VCHUNK)

                !  vieille version:
                ! time : t_(n+1/2) for solid ; t_n for fluid
                ! t = merge(timer+Tdomain%TimeD%dtmin/2d0,timer,Tdomain%specel(nel)%solid)
                ! nouvelle version:
                ! le temps n'est plus decale pour les sources, pour un saute-mouton
                !   on rajoute le 1/2 pas de temps qui correspond au fait que la
                !    exterieure doive etre prise a t_(n+1/2)
                t = timer+Tdomain%TimeD%dtmin/2_fpp
                ! TAG_SOURCE
                ft = CompSource(Tdomain%sSource(ns), t, ntime)
                if(Tdomain%sSource(ns)%i_type_source == 1 .or. Tdomain%sSource(ns)%i_type_source == 2) then
                    ! collocated force in solid
                    !
                    do i_dir = 0,2
                        do k = 0,ngll-1
                            do j = 0,ngll-1
                                do i = 0,ngll-1
                                    idx = Tdomain%sdom%Idom_(i,j,k,bnum,ee)
                                    val = Tdomain%sdom%champs(i1)%Forces(idx, i_dir) + ft*Tdomain%sSource(ns)%ExtForce(i,j,k,i_dir)
                                    Tdomain%sdom%champs(i1)%Forces(idx, i_dir) = val
                                enddo
                            enddo
                        enddo
                    enddo
                else if(Tdomain%sSource(ns)%i_type_source == 3)then    ! pressure pulse in fluid
                    do k = 0,ngll-1
                        do j = 0,ngll-1
                            do i = 0,ngll-1
                                idx = Tdomain%fdom%Idom_(i,j,k,bnum,ee)
                                val = Tdomain%fdom%champs(i1)%ForcesFl(idx)
                                val = val + ft*Tdomain%sSource(ns)%ExtForce(i,j,k,0)
                                !write(*,*) ntime,nel,i,j,k,val
                                Tdomain%fdom%champs(i1)%ForcesFl(idx) = val
                            enddo
                        enddo
                    enddo
                end if
            endif
        enddo

        return
    end subroutine external_forces

end module mtimestep
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
