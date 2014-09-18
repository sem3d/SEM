module scomm

contains

    subroutine exchange_sem_var(Tdomain, tag, vector)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: tag
        type(comm_vector), intent(inout) :: vector

        integer, dimension(MPI_STATUS_SIZE,vector%ncomm) :: statuses
        integer :: dest, src, ierr, n, i

        !- now we can exchange (communication global arrays)
        vector%send_reqs = MPI_REQUEST_NULL
        vector%recv_reqs = MPI_REQUEST_NULL

        do i = 0,vector%ncomm-1
            dest = vector%Data(i)%dest
            src = vector%Data(i)%src
            call MPI_Isend(vector%Data(i)%Give, vector%Data(i)%ndata, &
                           MPI_DOUBLE_PRECISION, dest, tag, Tdomain%communicateur, &
                           vector%send_reqs(i), ierr)

            call MPI_Irecv(vector%Data(i)%Take, vector%Data(i)%ndata, &
                           MPI_DOUBLE_PRECISION, dest, tag, Tdomain%communicateur, &
                           vector%recv_reqs(i), ierr)
        enddo

        call MPI_Waitall(vector%ncomm, vector%recv_reqs, statuses, ierr)
        call MPI_Waitall(vector%ncomm, vector%send_reqs, statuses, ierr)

    end subroutine exchange_sem_var

!    subroutine exchange_sem_1d(Tdomain, rg, tag, vector)
!        use sdomain
!        use mpi
!        implicit none
!        type(domain), intent(inout) :: Tdomain
!        integer, intent(in) :: rg
!        integer, intent(in) :: tag
!        type(comm_vector_1d), intent(inout) :: vector
!        integer :: other, ierr, n
!
!        !- now we can exchange (communication global arrays)
!        n = Tdomain%n_proc
!        vector%send_reqs = MPI_REQUEST_NULL
!        vector%recv_reqs = MPI_REQUEST_NULL
!
!        do other = 0,n-1
!            if (other==rg) cycle
!
!            if (vector%n(other)/=0) then
!                call MPI_Isend(vector%Give(:,other), vector%n(other), &
!                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, vector%send_reqs(other), ierr)
!                call MPI_Irecv(vector%Take(:,other), vector%n(other), &
!                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, vector%recv_reqs(other), ierr)
!            endif
!        enddo
!    end subroutine exchange_sem_1d
!
!    subroutine  exchange_sem_wait_1d(Tdomain, vector)
!        use sdomain
!        use mpi
!        implicit none
!        type(domain), intent(inout) :: Tdomain
!        type(comm_vector_1d), intent(inout) :: vector
!        integer, dimension(MPI_STATUS_SIZE, Tdomain%n_proc) :: statuses
!        integer :: ierr
!
!        call MPI_Waitall(Tdomain%n_proc, vector%recv_reqs, statuses, ierr)
!        call MPI_Waitall(Tdomain%n_proc, vector%send_reqs, statuses, ierr)
!    end subroutine exchange_sem_wait_1d

#if ! NEW_GLOBAL_METHOD
    subroutine exchange_sem(Tdomain, rg)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer :: n, k
        integer :: other, ierr
        integer, dimension(Tdomain%n_proc) :: send_req, recv_req, send_pml_req, recv_pml_req
        integer, parameter :: tag=101, tag_pml=102
        integer, dimension(MPI_STATUS_SIZE,Tdomain%n_proc) :: statuses

        !- now we can exchange (communication global arrays)
        n = Tdomain%n_proc
        send_req = MPI_REQUEST_NULL
        recv_req = MPI_REQUEST_NULL
        send_pml_req = MPI_REQUEST_NULL
        recv_pml_req = MPI_REQUEST_NULL

        do other = 0,n-1
            if (other==rg) cycle
            if (Tdomain%sComm(other)%ngll_tot/=0) then
                call MPI_Isend(Tdomain%sComm(other)%Give, Tdomain%sComm(other)%ngll_tot, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, send_req(other+1), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%Take, Tdomain%sComm(other)%ngll_tot, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, recv_req(other+1), ierr)
            endif
            if(Tdomain%any_PML.or.Tdomain%any_FPML)then
                if (Tdomain%any_FPML) then
                    k = 6
                else
                    k = 4
                end if
                if (Tdomain%sComm(other)%ngllPML_tot/=0) then
                    call MPI_Isend(Tdomain%sComm(other)%GivePML, k*Tdomain%sComm(other)%ngllPML_tot, &
                        MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, &
                        send_pml_req(other+1), ierr)
                    call MPI_Irecv(Tdomain%sComm(other)%TakePML, k*Tdomain%sComm(other)%ngllPML_tot, &
                        MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, &
                        recv_pml_req(other+1), ierr)
                endif
            endif
        enddo

        !write(*,*) "COMM done 1"
        call MPI_Waitall(Tdomain%n_proc, recv_req, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, send_req, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, send_pml_req, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, recv_pml_req, statuses, ierr)

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
            if (other==rg) cycle

            if (Tdomain%sComm(other)%ngll/=0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForces, 3*Tdomain%sComm(other)%ngll, &
                    MPI_DOUBLE_PRECISION, other, tag_sl, Tdomain%communicateur, req_s_f(other), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForces, 3*Tdomain%sComm(other)%ngll, &
                    MPI_DOUBLE_PRECISION, other, tag_sl, Tdomain%communicateur, req_r_f(other), ierr)
            endif
            if (Tdomain%sComm(other)%ngll_F/=0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForcesFl, 1*Tdomain%sComm(other)%ngll_F, &
                    MPI_DOUBLE_PRECISION, other, tag_fl, Tdomain%communicateur, req_s_fl(other), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForcesFL,1*Tdomain%sComm(other)%ngll_F, &
                    MPI_DOUBLE_PRECISION, other, tag_fl, Tdomain%communicateur, req_r_fl(other), ierr)
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
#endif
    subroutine exchange_sem_forces_StoF(Tdomain, rg)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer :: n
        integer :: other, ierr
        integer, dimension(Tdomain%n_proc) :: req_s, req_r, req_s_pml,req_r_pml
        integer, parameter :: tag_sf = 301, tag_sf_pml = 302
        integer, dimension(MPI_STATUS_SIZE,Tdomain%n_proc) :: statuses
        !- now we can exchange (communication global arrays)
        n = Tdomain%n_proc
        req_s = MPI_REQUEST_NULL
        req_r = MPI_REQUEST_NULL
        req_s_pml = MPI_REQUEST_NULL
        req_r_pml = MPI_REQUEST_NULL


        do other = 0,n-1
            if (other == rg) cycle
            if (Tdomain%sComm(other)%ngllSF > 0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForcesSF_StoF,Tdomain%sComm(other)%ngllSF, &
                     MPI_DOUBLE_PRECISION, other, tag_sf, Tdomain%communicateur,req_s(other+1), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForcesSF_StoF,Tdomain%sComm(other)%ngllSF, &
                     MPI_DOUBLE_PRECISION, other, tag_sf,Tdomain%communicateur,req_r(other+1), ierr)
            endif
            if (Tdomain%sComm(other)%ngllSF_PML > 0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForcesSF_StoF_PML,3*Tdomain%sComm(other)%ngllSF_PML, &
                     MPI_DOUBLE_PRECISION, other, tag_sf_pml,Tdomain%communicateur,req_s_pml(other+1), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForcesSF_StoF_PML,3*Tdomain%sComm(other)%ngllSF_PML, &
                     MPI_DOUBLE_PRECISION, other,tag_sf_pml,Tdomain%communicateur,req_r_pml(other+1), ierr)
            endif
        enddo

        call MPI_Waitall(Tdomain%n_proc, req_s, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_r, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_s_pml, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_r_pml, statuses, ierr)

    end subroutine exchange_sem_forces_StoF

    subroutine exchange_sem_forces_FtoS(Tdomain, rg)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer, intent(in) :: rg
        integer :: n
        integer :: other, ierr
        integer, dimension(0:Tdomain%n_proc) :: req_s, req_r, req_s_pml,req_r_pml
        integer, parameter :: tag = 101
        integer, dimension(MPI_STATUS_SIZE,Tdomain%n_proc) :: statuses
        !- now we can exchange (communication global arrays)
        n = Tdomain%n_proc
        req_s = MPI_REQUEST_NULL
        req_r = MPI_REQUEST_NULL
        req_s_pml = MPI_REQUEST_NULL
        req_r_pml = MPI_REQUEST_NULL

        do other = 0,n-1
            if (other == rg) cycle
            if (Tdomain%sComm(other)%ngllSF > 0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForcesSF_FtoS,3*Tdomain%sComm(other)%ngllSF, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur,req_s(other), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForcesSF_FtoS,3*Tdomain%sComm(other)%ngllSF, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur,req_r(other), ierr)
            endif
            if (Tdomain%sComm(other)%ngllSF_PML > 0) then
                call MPI_Isend(Tdomain%sComm(other)%GiveForcesSF_FtoS_PML,9*Tdomain%sComm(other)%ngllSF_PML, &
                     MPI_DOUBLE_PRECISION, other, tag,Tdomain%communicateur,req_s_pml(other), ierr)
                call MPI_Irecv(Tdomain%sComm(other)%TakeForcesSF_FtoS_PML,9*Tdomain%sComm(other)%ngllSF_PML, &
                     MPI_DOUBLE_PRECISION, other, tag,Tdomain%communicateur,req_r_pml(other), ierr)
            endif

        enddo

        call MPI_Waitall(Tdomain%n_proc, req_s, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_r, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_s_pml, statuses, ierr)
        call MPI_Waitall(Tdomain%n_proc, req_r_pml, statuses, ierr)

    end subroutine exchange_sem_forces_FtoS


end module scomm
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
