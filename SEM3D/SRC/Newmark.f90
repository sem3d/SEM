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
    ! From Berland, 2007 : doi:10.1016/j.compfluid.2005.04.003
    ! 
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
        use dom_solid
        use dom_solid_dg
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: ntime
        !
        real(fpp) :: cb, cg, cc, t0, t, dt
        integer :: i
        t0 = Tdomain%TimeD%rtime
        dt = Tdomain%TimeD%dtmin
        do i=1,6
            cb = LDDRK_beta(i)
            cg = LDDRK_gamma(i)
            cc = LDDRK_c(i)
            t = t0 + cc*dt

            call lddrk_init_solid(Tdomain%sdom, 0)
            call lddrk_init_solid_dg(Tdomain%sdomdg, 0)
            call internal_forces(Tdomain, 0, 0, ntime)
            call external_forces(Tdomain, t, ntime, 0)
            call comm_forces(Tdomain, 0)
            call lddrk_update_solid(Tdomain%sdom, 0, 1, dt, cb, cg)
            call lddrk_update_solid_dg(Tdomain%sdomdg, 0, 1, dt, cb, cg)
        end do
    end subroutine Timestep_LDDRK
    !
    subroutine Newmark(Tdomain,ntime)
        ! Predictor-MultiCorrector Newmark Velocity Scheme within a
        ! Time staggered Stress-Velocity formulation inside PML
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
            call stat_starttick(STAT_FEXT)
            call external_forces(Tdomain,Tdomain%TimeD%rtime,ntime, 1)
            call stat_stoptick(STAT_FEXT)
        end if

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
#ifdef CPML
        if(Tdomain%logicD%SF_local_present)then
            call StoF_coupling(Tdomain,0,1)
            !- fluid -> solid coupling (pressure times velocity)
            call FtoS_coupling(Tdomain,0,1)
        end if
        !- correction phase
        call Newmark_Corrector_F(Tdomain)
        !- solid -> fluid coupling (normal dot velocity)
        call Newmark_Corrector_S(Tdomain)
#else
        if(Tdomain%logicD%SF_local_present)then
            call StoF_coupling(Tdomain,0,1)
        end if
        !- correction phase
        call Newmark_Corrector_F(Tdomain)
        if(Tdomain%logicD%SF_local_present)then
            !- fluid -> solid coupling (pressure times velocity)
            call FtoS_coupling(Tdomain,0,1)
        end if
        !- solid -> fluid coupling (normal dot velocity)
        call Newmark_Corrector_S(Tdomain)
#endif

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
            call stat_starttick(STAT_GIVE)
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
            call stat_starttick(STAT_TAKE)
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
        use dom_solid_dg
        use dom_solidpml
        use stat, only : stat_starttick, stat_stoptick, STAT_FSOL, STAT_FFLU, STAT_PSOL, STAT_PFLU, STAT_FSOL_DG
        implicit none

        type(domain), intent(inout)   :: Tdomain

        ! Elements solide
        if (Tdomain%sdom%nglltot /= 0) then
            call stat_starttick(STAT_FSOL)
            call newmark_predictor_solid(Tdomain%sdom,0,1)
            call stat_stoptick(STAT_FSOL)
        endif

        ! Elements solide  DG
        if (Tdomain%sdomdg%nglltot /= 0) then
            call stat_starttick(STAT_FSOL_DG)
            call newmark_predictor_solid_dg(Tdomain%sdomdg,0,1)
            call stat_stoptick(STAT_FSOL_DG)
        endif

        ! Elements fluide
        if (Tdomain%fdom%nglltot /= 0) then
            call stat_starttick(STAT_FFLU)
            call newmark_predictor_fluid(Tdomain%fdom,0,1)
            call stat_stoptick(STAT_FFLU)
        endif

        ! Elements solide pml
        if (Tdomain%spmldom%nglltot /= 0) then
            call stat_starttick(STAT_PSOL)
            call newmark_predictor_solidpml(Tdomain%spmldom, Tdomain,0,1)
            call stat_stoptick(STAT_PSOL)
        endif

        ! Elements fluide pml
        if (Tdomain%fpmldom%nglltot /= 0) then
            call stat_starttick(STAT_PFLU)
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
            call stat_starttick(STAT_PFLU)
            call newmark_corrector_fluidpml(Tdomain%fpmldom, dt, 0, 1)
            call stat_stoptick(STAT_PFLU)
        endif
        ! Si il existe des éléments fluides
        if (Tdomain%fdom%nglltot /= 0) then
            call stat_starttick(STAT_FFLU)
            call newmark_corrector_fluid(Tdomain%fdom, dt, 0, 1)
            call stat_stoptick(STAT_FFLU)
        endif

        return
    end subroutine Newmark_Corrector_F
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    subroutine Newmark_Corrector_S(Tdomain)
        use dom_solid
        use dom_solid_dg
        use dom_solidpml
        use stat, only : stat_starttick, stat_stoptick, STAT_PSOL, STAT_FSOL, STAT_FSOL_DG
        implicit none

        type(domain), intent(inout)   :: Tdomain
        real(fpp) :: dt, t

        dt = Tdomain%TimeD%dtmin
        t = Tdomain%TimeD%rtime
        ! Si il existe des éléments PML solides
        if (Tdomain%spmldom%nglltot /= 0) then
            call stat_starttick(STAT_PSOL)
            call newmark_corrector_solidpml(Tdomain%spmldom, dt, t, 0, 1)
            call stat_stoptick(STAT_PSOL)
        endif

        ! Si il existe des éléments solides
        if (Tdomain%sdom%nglltot /= 0) then
            call stat_starttick(STAT_FSOL)
            call newmark_corrector_solid(Tdomain%sdom, dt, 0, 1)
            call stat_stoptick(STAT_FSOL)
        endif

        ! Si il existe des éléments solides DG
        if (Tdomain%sdomdg%nglltot /= 0) then
            call stat_starttick(STAT_FSOL_DG)
            call newmark_corrector_solid_dg(Tdomain%sdomdg, dt, 0, 1)
            call stat_stoptick(STAT_FSOL_DG)
        endif
        return
    end subroutine Newmark_Corrector_S
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    subroutine internal_forces(Tdomain, i0, i1, ntime)
        ! volume forces - depending on rheology
        use dom_solid
        use dom_solid_dg
        use dom_solidpml
        use dom_fluid
        use dom_fluidpml
        use smirror
        use stat, only : stat_starttick, stat_stoptick, STAT_FFLU, STAT_PFLU, STAT_FSOL, STAT_PSOL, STAT_FSOL_DG
        implicit none

        type(domain), intent(inout)  :: Tdomain
        integer, intent(in) :: i0, i1, ntime
        integer  :: n, indsol, indflu, indpml

        ! DOMAIN FLUID
        if (Tdomain%fdom%nbelem>0) then
            call stat_starttick(STAT_FFLU)
            ! Mirror Record
            if (Tdomain%mirror_type==0.and.Tdomain%fdom%mirror_fl%n_glltot>0) then
                do n = 0,Tdomain%fdom%nblocks-1
                    call forces_int_fluid_mirror_dump(Tdomain%fdom,Tdomain%fdom%champs(i1),n)
                enddo
                call dump_mirror_fl(Tdomain%fdom,ntime)
            ! Mirror Forward/Backward
            elseif (Tdomain%mirror_type>0.and.Tdomain%fdom%mirror_fl%n_glltot>0) then
                call load_mirror_fl(Tdomain%fdom,ntime)
                do n = 0,Tdomain%fdom%nblocks-1
                    call forces_int_fluid_mirror_load(Tdomain%fdom,Tdomain%fdom%champs(i1),n)
                enddo
            ! Without Mirror
            else
                do n = 0,Tdomain%fdom%nblocks-1
                    call forces_int_fluid(Tdomain%fdom,Tdomain%fdom%champs(i1),n)
                enddo
            endif
            call stat_stoptick(STAT_FFLU)
        endif
        ! DOMAIN PML FLUID
        if (Tdomain%fpmldom%nbelem>0) then
            call stat_starttick(STAT_PFLU)
            do n = 0,Tdomain%fpmldom%nblocks-1
                call pred_flu_pml(Tdomain%fpmldom, Tdomain%TimeD%dtmin, Tdomain%fpmldom%champs(i1), n)
                call forces_int_flu_pml(Tdomain%fpmldom, Tdomain%fpmldom%champs(i1), n, Tdomain)
            end do
            call stat_stoptick(STAT_PFLU)
        end if
        ! DOMAIN SOLID
        if (Tdomain%sdom%nbelem>0) then
            call stat_starttick(STAT_FSOL)
            ! Mirror Record
            if (Tdomain%mirror_type==0.and.Tdomain%sdom%mirror_sl%n_glltot>0) then
                if (Tdomain%mirror_expl) then
                    do n = 0,Tdomain%sdom%nblocks-1
                        call forces_int_solid_mirror_dump_expl(Tdomain%sdom,Tdomain%sdom%champs(i1),n)
                    enddo
                else
                    do n = 0,Tdomain%sdom%nblocks-1
                        call forces_int_solid_mirror_dump(Tdomain%sdom,Tdomain%sdom%champs(i1),n)
                    enddo
                endif
                call dump_mirror_sl(Tdomain%sdom,ntime)
            ! Mirror Forward/Backward
            elseif (Tdomain%mirror_type>0.and.Tdomain%sdom%mirror_sl%n_glltot>0) then
                call load_mirror_sl(Tdomain%sdom,ntime)
                if (Tdomain%mirror_recalc) then
                    do n = 0,Tdomain%sdom%nblocks-1
                        call forces_int_solid_mirror_load_recalc(Tdomain%sdom,Tdomain%sdom%champs(i1),n)
                    enddo
                else
                    if (Tdomain%mirror_expl) then
                        do n = 0,Tdomain%sdom%nblocks-1
                            call forces_int_solid_mirror_load_expl(Tdomain%sdom,Tdomain%sdom%champs(i1),n)
                        enddo
                    else
                        do n = 0,Tdomain%sdom%nblocks-1
                            call forces_int_solid_mirror_load(Tdomain%sdom,Tdomain%sdom%champs(i1),n)
                        enddo
                    endif
                endif
            ! Without Mirror
            else
                do n = 0,Tdomain%sdom%nblocks-1
                    call forces_int_solid(Tdomain%sdom,Tdomain%sdom%champs(i1),n,Tdomain%nl_flag)
                enddo
            endif
            call stat_stoptick(STAT_FSOL)
        endif
        !! DOMAIN SOLID DG
        if (Tdomain%sdomdg%nbelem>0) then
            call stat_starttick(STAT_FSOL_DG)
            do n = 0,Tdomain%sdomdg%nblocks-1
                call forces_int_solid_dg(Tdomain%sdomdg,Tdomain%sdomdg%champs(i1),n)
            enddo
            call stat_stoptick(STAT_FSOL_DG)
        endif
        ! DOMAIN PML SOLID
        if (Tdomain%spmldom%nbelem>0) then
            call stat_starttick(STAT_PSOL)
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
    subroutine external_forces(Tdomain,timer,ntime,i1)

        implicit none
#include "index.h"

        type(domain), intent(inout)  :: Tdomain
        integer, intent(in)  :: ntime
        real(kind=fpp), intent(in)  :: timer
        integer, intent(in) :: i1
        integer :: ns,nel,i_dir, i,j,k, idx, lnum,ngll, bnum, ee, ntimecur, dom
        real(kind=fpp) :: t, ft, val, timercur

        if (Tdomain%mirror_type==1) return
        !!!if (Tdomain%mirror_type>=1) return

        timercur = timer
        ntimecur = ntime
        if (Tdomain%mirror_type==2) then
            timercur = Tdomain%TimeD%Duration-(Tdomain%TimeD%dtmin*3_fpp/2_fpp)-timer
            ntimecur = Tdomain%TimeD%ntimeMax-1-ntime
        endif

        do ns = 0, Tdomain%n_source-1
            if(Tdomain%rank == Tdomain%sSource(ns)%proc)then
                nel = Tdomain%Ssource(ns)%elem
                lnum = Tdomain%specel(nel)%lnum
                dom = Tdomain%specel(nel)%domain
                ngll = domain_ngll(Tdomain, dom)
                bnum = lnum/VCHUNK
                ee = mod(lnum,VCHUNK)

                !  vieille version:
                ! time : t_(n+1/2) for solid ; t_n for fluid
                ! t = merge(timer+Tdomain%TimeD%dtmin/2d0,timer,Tdomain%specel(nel)%solid)
                ! nouvelle version:
                ! le temps n'est plus decale pour les sources, pour un saute-mouton
                !   on rajoute le 1/2 pas de temps qui correspond au fait que la
                !    exterieure doive etre prise a t_(n+1/2)
                t = timercur+Tdomain%TimeD%dtmin/2_fpp
                ! TAG_SOURCE
                ft = CompSource(Tdomain%sSource(ns), t, ntimecur)
                if(Tdomain%sSource(ns)%i_type_source == 1 .or. Tdomain%sSource(ns)%i_type_source == 2) then
                    ! collocated force in solid
                    !
                    do i_dir = 0,2
                        do k = 0,ngll-1
                            do j = 0,ngll-1
                                do i = 0,ngll-1
                                    if (dom==DM_SOLID_CG) then
                                        idx = Tdomain%sdom%Idom_(i,j,k,bnum,ee)
                                        val = Tdomain%sdom%champs(i1)%Forces(idx, i_dir) + ft*Tdomain%sSource(ns)%ExtForce(i,j,k,i_dir)
                                        Tdomain%sdom%champs(i1)%Forces(idx, i_dir) = val
                                    else if (dom==DM_SOLID_DG) then
                                        idx = Tdomain%sdomdg%Idom_(i,j,k,bnum,ee)
                                        val = Tdomain%sdomdg%champs(i1)%Forces(idx, i_dir) + ft*Tdomain%sSource(ns)%ExtForce(i,j,k,i_dir)
                                        Tdomain%sdomdg%champs(i1)%Forces(idx, i_dir) = val
                                    end if
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
                                !write(*,*) ntimecur,nel,i,j,k,val
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
