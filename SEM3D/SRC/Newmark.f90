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

        real(fpp) :: cb, cg, cc, t0, t, dt, bnum
        integer :: l, i, ngll, ii, j, k, ee, n

        t0 = Tdomain%TimeD%rtime
        dt = Tdomain%TimeD%dtmin
        ngll = Tdomain%sdomdg%ngll

        do l = 1,6

            cb = LDDRK_beta(l)
            cg = LDDRK_gamma(l)
            cc = LDDRK_c(l)

            t  = t0 + cc*dt

            call lddrk_init_solid(Tdomain%sdom, 2)
            call lddrk_init_solid_dg(Tdomain%sdomdg, 2)
            call internal_forces(Tdomain, 0, 2, ntime)
            call external_forces(Tdomain, t, ntime, 2)

            call assemble_forces_dg(Tdomain, Tdomain%sdomdg, 2)
            call comm_forces(Tdomain, 2)
            call deassemble_forces_dg(Tdomain, Tdomain%sdomdg,2)

            call lddrk_update_solid(Tdomain%sdom, 0, 1, 2, dt, cb, cg)
            call lddrk_update_solid_dg(Tdomain%sdomdg, 0, 1, 2, dt, cb, cg)

         end do
        !ngll = Tdomain%sdomdg%ngll
        !do n = 0,Tdomain%n_elem-1
        !    bnum = Tdomain%specel(n)%lnum/VCHUNK
        !    ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
        !    do k = 0,ngll-1
        !        do j = 0,ngll-1
        !            do ii = 0,ngll-1
        !                do ee = 0, VCHUNK-1
        !
        !                    write(*,*) "ee = ", ee, " i = ", ii, " j = ", j, " k = ", k
        !                    write(*,*) ii,j,k,"E11 = ", Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,0,:)
        !                    write(*,*) "E22 = ", Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,1,:)
        !                    write(*,*) "E33 = ", Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,2,:)
        !                    write(*,*) "E12 = ", Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,3,:)
        !                    write(*,*) "E13 = ", Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,4,:)
        !                    write(*,*) "E23 = ", Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,5,:)
        !                    write(*,*) ii,j,k,"Vx = ",  Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,6,:)
        !                    write(*,*) i,j,k,"Vy = ",  Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,7,:)
        !                    write(*,*) i,j,k,"Vz = " , Tdomain%sdomdg%champs(0)%Q(ee,ii,j,k,8,:)
        !
        !                enddo
        !            enddo
        !        enddo
        !    enddo
        !enddo

    end subroutine Timestep_LDDRK

    subroutine assemble_forces_dg(Tdomain, dom, f1)

        implicit none

        type(domain), intent(inout) :: Tdomain
        type(domain_solid_dg),intent(inout)       :: dom
        integer, intent(in)         :: f1
        integer                     :: nbelem, ngll
        integer                     :: n,i,j,k,ee,idx, bnum

        nbelem = dom%nbelem
        ngll   = dom%ngll

        dom%Qasm(:,:) = 0.d0

        do n = 0,Tdomain%n_elem-1
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee   = mod(Tdomain%specel(n)%lnum,VCHUNK)
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        dom%Qasm(idx,6:8) = dom%Qasm(idx,6:8) + dom%champs(f1)%Q(ee,i,j,k,6:8,bnum)
                    enddo
                enddo
            enddo
        enddo

    end subroutine assemble_forces_dg

    subroutine deassemble_forces_dg(Tdomain, dom, f1)

        implicit none

        type(domain), intent(inout) :: Tdomain
        type(domain_solid_dg),intent(inout)       :: dom
        integer, intent(in)         :: f1
        integer                     :: ngll
        integer                     :: n,i,j,k,ee, bnum, idx
        real(fpp), dimension(0:8)   :: val
        real(fpp)                   :: wh

        ngll   = dom%ngll
 
        do n = 0,Tdomain%n_elem-1
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        wh = dom%GLLw(i)*dom%GLLw(j)*dom%GLLw(k)
                        do ee = 0, VCHUNK-1
                            idx = dom%Idom_(i,j,k,bnum,ee)
                            dom%champs(f1)%Q(ee,i,j,k,0:5,bnum) = dom%champs(f1)%Q(ee,i,j,k,0:5,bnum) / (wh*dom%Jacob_(i,j,k,bnum,ee))
                            dom%champs(f1)%Q(ee,i,j,k,6:8,bnum) = dom%Qasm(idx,6:8) * dom%MassMat(idx)
                        enddo
                    enddo
                enddo
            enddo
        enddo

    end subroutine deassemble_forces_dg

    subroutine Newmark(Tdomain,ntime)
        ! Predictor-MultiCorrector Newmark Velocity Scheme within a
        ! Time staggered Stress-Velocity formulation inside PML
        use mcapteur
        use sem_mpi
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
        !$acc wait(1)

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
                if (Tdomain%Comm_data%Data(n)%nsol>0) then
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveS, Tdomain%sdom%champs(f1)%Veloc, k)
                end if
                if (Tdomain%Comm_data%Data(n)%nsoldg>0) then
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveSDG, Tdomain%sdomdg%Qasm, k)
                end if

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
                if (Tdomain%Comm_data%Data(n)%nflu>0) then
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveF, Tdomain%fdom%champs(f1)%ForcesFl, k)
                end if
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
                if (Tdomain%Comm_data%Data(n)%nsol>0) then
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveS, Tdomain%sdom%champs(f1)%Veloc, k)
                end if
                ! Domain SOLID DG
                if (Tdomain%Comm_data%Data(n)%nsoldg>0) then
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveSDG, Tdomain%sdomdg%Qasm, k)
                end if

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
                if (Tdomain%Comm_data%Data(n)%nflu>0) then
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveF, Tdomain%fdom%champs(f1)%ForcesFl, k)
                end if

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
            !call stat_starttick(STAT_FSOL_DG)
            !call newmark_predictor_solid_dg(Tdomain%sdomdg,0,1)
            !call stat_stoptick(STAT_FSOL_DG)
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
            !call stat_starttick(STAT_FSOL_DG)
            !call newmark_corrector_solid_dg(Tdomain%sdomdg, dt, 0, 1)
            !call stat_stoptick(STAT_FSOL_DG)
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
        logical :: m_dump, m_load, m_expl, m_recalc
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
            m_load = .false.
            m_dump = .false.
            m_recalc = .false.
            m_expl = .false.
            call stat_starttick(STAT_FSOL)
            ! Mirror Record
            if (Tdomain%mirror_type==0.and.Tdomain%sdom%mirror_sl%n_glltot>0) then
                m_dump = .true.
                if (Tdomain%mirror_expl) then
                    m_expl = .true.
                endif
            ! Mirror Forward/Backward
            elseif (Tdomain%mirror_type>0.and.Tdomain%sdom%mirror_sl%n_glltot>0) then
                m_load = .true.
                call load_mirror_sl(Tdomain%sdom,ntime)
                if (Tdomain%mirror_recalc) then
                    m_recalc = .true.
                else
                    if (Tdomain%mirror_expl) then
                        m_expl = .true.
                    endif
                endif
            ! Without Mirror
            end if
            call forces_int_solid_mainloop(Tdomain%sdom, i0, i1, Tdomain%nl_flag, m_dump, m_load, m_expl,m_recalc)
            if (m_dump) then
                call dump_mirror_sl(Tdomain%sdom,ntime)
            endif

            call stat_stoptick(STAT_FSOL)
        endif
        !! DOMAIN SOLID DG
        if (Tdomain%sdomdg%nbelem>0) then
            call stat_starttick(STAT_FSOL_DG)
            do n = 0,Tdomain%sdomdg%nblocks-1
                call forces_int_solid_dg(Tdomain%sdomdg, Tdomain%sdomdg%champs(i0), Tdomain%sdomdg%champs(i1),n)
            enddo

            call stat_stoptick(STAT_FSOL_DG)
        endif
        ! DOMAIN PML SOLID
        if (Tdomain%spmldom%nbelem>0) then
            call stat_starttick(STAT_PSOL)
            call forces_int_sol_pml_mainloop(Tdomain%spmldom, i1)
            call stat_stoptick(STAT_PSOL)
        end if

        ! Couplage interface solide / PML
        if (Tdomain%spmldom%nglltot > 0) then
            call couplage_pml_solid(Tdomain, Tdomain%sdom, Tdomain%spmldom, i1)
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
        use dom_solid
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
                    call apply_source_solid(Tdomain%sSource(ns), Tdomain%sdom, i1, ft, lnum)
                else if(Tdomain%sSource(ns)%i_type_source == 3)then    ! pressure pulse in fluid
                    do k = 0,ngll-1
                        do j = 0,ngll-1
                            do i = 0,ngll-1
                                idx = Tdomain%fdom%Idom_(i,j,k,bnum,ee)
                                val = Tdomain%fdom%champs(i1)%ForcesFl(idx)
                                val = val + ft*Tdomain%sSource(ns)%ExtForce(i,j,k,0)
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
