module scomm

contains
    subroutine shift_to_parameters(rank,n,shift,give_to,take_from,rings)

        implicit none

        integer, intent(in)  :: rank,n,shift
        integer, intent(out) :: give_to,take_from,rings

        give_to = rank+shift
        if(give_to > n-1) give_to = give_to-n
        take_from = rank-shift
        if(take_from < 0) take_from = take_from+n
        if(mod(n,shift) == 0 .and. shift /= 1)then
            rings = shift
        else if(mod(n,n-shift) == 0 .and. shift /= n-1)then
            rings = n-shift
        else if(mod(n,2) == 0 .and. mod(shift,2) == 0)then
            rings = 2
        else
            rings = 1
        endif

    end subroutine shift_to_parameters
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    subroutine ALGO_COMM(rank,n,rings,ngll_give,ngll_take,give_to,take_from,fact_mult,  &
        TabGiven,TabTaken)
        ! old algo to exchange properties between processors
        use mpi
        implicit none
        integer, intent(in) :: rings,rank,n,ngll_give,ngll_take,give_to,take_from,fact_mult
        real, dimension(0:fact_mult*ngll_give-1), intent(in) :: TabGiven
        real, dimension(0:fact_mult*ngll_take-1), intent(out) :: TabTaken
        integer  :: i,j,code,etiquette
        integer, dimension(MPI_STATUS_SIZE) :: statut

        etiquette = 100

        do i = 0,rings-1
            if(rank == i)then
                if(ngll_give > 0)then
                    call MPI_SEND(TabGiven,fact_mult*ngll_give,MPI_DOUBLE_PRECISION,  &
                        give_to,etiquette,MPI_COMM_WORLD,code)
                endif
                if(ngll_take > 0)then
                    call MPI_RECV(TabTaken,fact_mult*ngll_take,MPI_DOUBLE_PRECISION,  &
                        take_from,etiquette,MPI_COMM_WORLD,statut,code)
                endif
            else
                do j = 0,n/rings-1
                    if(rank == i + j*rings)then
                        if(ngll_take > 0)then
                            call MPI_RECV(TabTaken,fact_mult*ngll_take,MPI_DOUBLE_PRECISION, &
                                take_from,etiquette,MPI_COMM_WORLD,statut,code)
                        endif
                        if(ngll_give > 0)then
                            call MPI_SEND(TabGiven,fact_mult*ngll_give,MPI_DOUBLE_PRECISION, &
                                give_to,etiquette,MPI_COMM_WORLD,code)
                        endif
                    endif
                enddo
            endif
        enddo

        return

    end subroutine ALGO_COMM


    subroutine exchange_sem_1(Tdomain, rg)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer :: n, shift, I_give_to, I_take_from, n_rings, code

        write(*,*) "COMM 1"
        !- now we can exchange (communication global arrays)
        n = Tdomain%n_proc
        do shift = 1,n-1
            call shift_to_parameters(rg,n,shift,I_give_to,I_take_from,n_rings)
            ! physical mass comm. (solid + fluid)
            call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngll_tot,                  &
                Tdomain%sComm(I_take_from)%ngll_tot,I_give_to,I_take_from,1,     &
                Tdomain%sComm(I_give_to)%Give,Tdomain%sComm(I_take_from)%Take)
            call MPI_BARRIER(MPI_COMM_WORLD,code)
            ! PML mass comm. (solid + fluid)
            if(Tdomain%any_PML)then
                if(Tdomain%any_FPML)then
                    call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngllPML_tot,           &
                        Tdomain%sComm(I_take_from)%ngllPML_tot,I_give_to,I_take_from,6,       &
                        Tdomain%sComm(I_give_to)%GivePML,Tdomain%sComm(I_take_from)%TakePML)
                else
                    call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngllPML_tot,           &
                        Tdomain%sComm(I_take_from)%ngllPML_tot,I_give_to,I_take_from,3,       &
                        Tdomain%sComm(I_give_to)%GivePML,Tdomain%sComm(I_take_from)%TakePML)
                end if
            end if
            call MPI_BARRIER(MPI_COMM_WORLD,code)
            ! normal Neumann comm.
            if(Tdomain%logicD%Neumann_local_present)then
                call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngllNeu,           &
                    Tdomain%sComm(I_take_from)%ngllNeu,I_give_to,I_take_from,3,       &
                    Tdomain%sComm(I_give_to)%GiveNeu,Tdomain%sComm(I_take_from)%TakeNeu)
            end if
            call MPI_BARRIER(MPI_COMM_WORLD,code)

        enddo
        write(*,*) "COMM 2"

    end subroutine exchange_sem_1

    subroutine exchange_sem_forces_1(Tdomain, rg)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer :: n, shift, I_give_to, I_take_from, n_rings, code

        write(*,*) "FCOMM 1"
        ! now we can exchange force values with proc n
        n = Tdomain%n_proc
        do shift = 1,n-1
            call shift_to_parameters(rg,n,shift,I_give_to,I_take_from,n_rings)
            ! solid forces exchange
            call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngll,        &
                Tdomain%sComm(I_take_from)%ngll,I_give_to,I_take_from,3,   &
                Tdomain%sComm(I_give_to)%GiveForces,Tdomain%sComm(I_take_from)%TakeForces)
            call MPI_BARRIER(MPI_COMM_WORLD,code)
            ! fluid forces exchange
            call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngll_F,        &
                Tdomain%sComm(I_take_from)%ngll_F,I_give_to,I_take_from,1,   &
                Tdomain%sComm(I_give_to)%GiveForcesFl,Tdomain%sComm(I_take_from)%TakeForcesFl)
            call MPI_BARRIER(MPI_COMM_WORLD,code)
            ! solid PML forces exchange
            call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngllPML,        &
                Tdomain%sComm(I_take_from)%ngllPML,I_give_to,I_take_from,9,   &
                Tdomain%sComm(I_give_to)%GiveForcesPML,Tdomain%sComm(I_take_from)%TakeForcesPML)
            call MPI_BARRIER(MPI_COMM_WORLD,code)
            ! fluid PML forces exchange
            call ALGO_COMM(rg,n,n_rings,Tdomain%sComm(I_give_to)%ngllPML_F,        &
                Tdomain%sComm(I_take_from)%ngllPML_F,I_give_to,I_take_from,3,   &
                Tdomain%sComm(I_give_to)%GiveForcesPMLFl,Tdomain%sComm(I_take_from)%TakeForcesPMLFl)
            call MPI_BARRIER(MPI_COMM_WORLD,code)

        end do  ! do shift
        write(*,*) "FCOMM 2"

    end subroutine exchange_sem_forces_1


    subroutine exchange_sem_1d(Tdomain, rg, tag, vector)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer, intent(in) :: tag
        type(comm_vector_1d), intent(inout) :: vector
        integer :: other, ierr, n

        !- now we can exchange (communication global arrays)
        n = Tdomain%n_proc
        vector%send_reqs = MPI_REQUEST_NULL
        vector%recv_reqs = MPI_REQUEST_NULL

        do other = 0,n-1
            if (other==rg) continue

            if (vector%n(other)/=0) then
                call MPI_Isend(vector%Give(:,other), vector%n(other), &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, vector%send_reqs(other), ierr)
                call MPI_Irecv(vector%Take(:,other), vector%n(other), &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, vector%recv_reqs(other), ierr)
            endif
        enddo
    end subroutine exchange_sem_1d

    subroutine  exchange_sem_wait_1d(Tdomain, vector)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        type(comm_vector_1d), intent(inout) :: vector
        integer, dimension(MPI_STATUS_SIZE, Tdomain%n_proc) :: statuses
        integer :: ierr

        call MPI_Waitall(Tdomain%n_proc, vector%recv_reqs, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, vector%send_reqs, statuses, ierr)
    end subroutine exchange_sem_wait_1d


    subroutine exchange_sem(Tdomain, rg)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer :: n, k
        integer :: other, ierr
        integer, dimension(0:Tdomain%n_proc) :: send_req, recv_req, send_pml_req, recv_pml_req
        integer, parameter :: tag=101, tag_pml=102
        integer, dimension(MPI_STATUS_SIZE,Tdomain%n_proc) :: statuses
        !write(*,*) "COMM 1"
        !- now we can exchange (communication global arrays)
        n = Tdomain%n_proc
        send_req = MPI_REQUEST_NULL
        recv_req = MPI_REQUEST_NULL
        send_pml_req = MPI_REQUEST_NULL
        recv_pml_req = MPI_REQUEST_NULL

        do other = 0,n-1
            if (other==rg) continue

            if (Tdomain%sComm(other)%ngll_tot/=0) then
                call MPI_Isend(Tdomain%sComm(other)%Give, Tdomain%sComm(other)%ngll_tot, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, send_req(other), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%Take, Tdomain%sComm(other)%ngll_tot, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, recv_req(other), ierr)
            endif
            if(Tdomain%any_PML.or.Tdomain%any_FPML)then
                if (Tdomain%any_FPML) then
                    k = 6
                else
                    k = 3
                end if
                if (Tdomain%sComm(other)%ngllPML_tot/=0) then
                    call MPI_Isend(Tdomain%sComm(other)%GivePML, k*Tdomain%sComm(other)%ngllPML_tot, &
                        MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, &
                        send_pml_req(other), ierr)
                    call MPI_Irecv(Tdomain%sComm(other)%TakePML, k*Tdomain%sComm(other)%ngllPML_tot, &
                        MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, &
                        recv_pml_req(other), ierr)
                endif
            endif
        enddo

        !write(*,*) "COMM done 1"
        call MPI_Waitall(Tdomain%n_proc, recv_req, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, recv_pml_req, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, send_req, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, send_pml_req, statuses, ierr)

    end subroutine exchange_sem

    subroutine exchange_sem_forces(Tdomain, rg)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer :: n
        integer :: other, ierr
        integer, dimension(0:Tdomain%n_proc) :: req_s_f, req_r_f, req_s_pml, req_r_pml
        integer, dimension(0:Tdomain%n_proc) :: req_s_fl, req_r_fl, req_s_fpml, req_r_fpml
        integer, parameter :: tag_sl=101, tag_fl=102, tag_pml=103, tag_fpml=104
        integer, dimension(MPI_STATUS_SIZE,Tdomain%n_proc) :: statuses
        !write(*,*) "ENTER Exchange sem forces", rg
        !- now we can exchange (communication global arrays)
        n = Tdomain%n_proc
        req_s_f = MPI_REQUEST_NULL
        req_r_f = MPI_REQUEST_NULL
        req_s_fl = MPI_REQUEST_NULL
        req_r_fl = MPI_REQUEST_NULL
        req_s_pml = MPI_REQUEST_NULL
        req_r_pml = MPI_REQUEST_NULL
        req_s_fpml = MPI_REQUEST_NULL
        req_r_fpml = MPI_REQUEST_NULL

        do other = 0,n-1
            if (other==rg) continue

            if (Tdomain%sComm(other)%ngll_tot/=0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForces, 3*Tdomain%sComm(other)%ngll, &
                    MPI_DOUBLE_PRECISION, other, tag_sl, Tdomain%communicateur, req_s_f(other), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForces, 3*Tdomain%sComm(other)%ngll, &
                    MPI_DOUBLE_PRECISION, other, tag_sl, Tdomain%communicateur, req_r_f(other), ierr)
            endif
            if (Tdomain%sComm(other)%ngll_F/=0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForcesFl, 1*Tdomain%sComm(other)%ngll_F, &
                    MPI_DOUBLE_PRECISION, other, tag_fl, Tdomain%communicateur, req_s_fl(other), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForces, 1*Tdomain%sComm(other)%ngll, &
                    MPI_DOUBLE_PRECISION, other, tag_fl, Tdomain%communicateur, req_r_f(other), ierr)
            endif
            if (Tdomain%sComm(other)%ngllPML/=0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForcesPML, 9*Tdomain%sComm(other)%ngllPML, &
                    MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, req_s_pml(other), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForcesPML, 9*Tdomain%sComm(other)%ngllPML, &
                    MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, req_r_pml(other), ierr)
            endif
            if (Tdomain%sComm(other)%ngllPML_F/=0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForcesPMLFl, 3*Tdomain%sComm(other)%ngllPML_F, &
                    MPI_DOUBLE_PRECISION, other, tag_fpml, Tdomain%communicateur, req_s_fpml(other), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForcesPMLFl, 3*Tdomain%sComm(other)%ngllPML_F, &
                    MPI_DOUBLE_PRECISION, other, tag_fpml, Tdomain%communicateur, req_r_fpml(other), ierr)
            endif
        enddo
        !write(*,*) "WAIT Exchange sem forces", rg

        call MPI_Waitall(Tdomain%n_proc, req_s_f, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_r_f, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_s_fl, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_r_fl, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_s_pml, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_r_pml, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_s_fpml, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_r_fpml, statuses, ierr)
        !write(*,*) "END Exchange sem forces", rg
    end subroutine exchange_sem_forces

end module scomm
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
