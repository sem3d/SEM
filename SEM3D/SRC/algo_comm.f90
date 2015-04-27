!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module scomm

contains

    subroutine exchange_sem(Tdomain)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: n, rg
        integer :: k
        integer :: other, ierr
        integer, dimension(Tdomain%tot_comm_proc) :: send_req, recv_req, send_pml_req, recv_pml_req
        integer, parameter :: tag=101, tag_pml=102
        integer, dimension(MPI_STATUS_SIZE,Tdomain%tot_comm_proc) :: statuses
        !write(*,*) "COMM 1"
        !- now we can exchange (communication global arrays)
        rg = Tdomain%rank
        send_req = MPI_REQUEST_NULL
        recv_req = MPI_REQUEST_NULL
        send_pml_req = MPI_REQUEST_NULL
        recv_pml_req = MPI_REQUEST_NULL

        do n = 0,Tdomain%tot_comm_proc-1
            other = Tdomain%sComm(n)%dest
            if (other==rg) cycle
            if (Tdomain%sComm(n)%ngll_tot/=0) then
                call MPI_Isend(Tdomain%sComm(n)%Give, Tdomain%sComm(n)%ngll_tot, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, send_req(n+1), ierr)
                call MPI_Irecv(Tdomain%sComm(n)%Take, Tdomain%sComm(n)%ngll_tot, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur, recv_req(n+1), ierr)
            endif
            if(Tdomain%any_PML.or.Tdomain%any_FPML)then
                if (Tdomain%any_FPML) then
                    k = 6
                else
                    k = 3
                end if
                if (Tdomain%sComm(n)%ngllPML_tot/=0) then
                    call MPI_Isend(Tdomain%sComm(n)%GivePML, k*Tdomain%sComm(n)%ngllPML_tot, &
                        MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, &
                        send_pml_req(n+1), ierr)
                    call MPI_Irecv(Tdomain%sComm(n)%TakePML, k*Tdomain%sComm(n)%ngllPML_tot, &
                        MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, &
                        recv_pml_req(n+1), ierr)
                endif
            endif
        enddo

        !write(*,*) "COMM done 1"
        call MPI_Waitall(Tdomain%tot_comm_proc, recv_req, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, send_req, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, send_pml_req, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, recv_pml_req, statuses, ierr)

    end subroutine exchange_sem

    subroutine exchange_sem_forces(Tdomain)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: rg
        integer :: n
        integer :: other, ierr
        integer, dimension(Tdomain%tot_comm_proc) :: req_s_f, req_r_f, req_s_pml, req_r_pml
        integer, dimension(Tdomain%tot_comm_proc) :: req_s_fl, req_r_fl, req_s_fpml, req_r_fpml
        integer, parameter :: tag_sl=101, tag_fl=102, tag_pml=103, tag_fpml=104
        integer, dimension(MPI_STATUS_SIZE,Tdomain%tot_comm_proc) :: statuses
        !write(*,*) "ENTER Exchange sem forces", rg
        !- now we can exchange (communication global arrays)
        rg = Tdomain%rank
        req_s_f = MPI_REQUEST_NULL
        req_r_f = MPI_REQUEST_NULL
        req_s_fl = MPI_REQUEST_NULL
        req_r_fl = MPI_REQUEST_NULL
        req_s_pml = MPI_REQUEST_NULL
        req_r_pml = MPI_REQUEST_NULL
        req_s_fpml = MPI_REQUEST_NULL
        req_r_fpml = MPI_REQUEST_NULL

        do n = 0,Tdomain%tot_comm_proc-1
            other = Tdomain%sComm(n)%dest
            if (other==rg) cycle

            if (Tdomain%sComm(n)%ngll/=0) then
                call MPI_Isend(Tdomain%sComm(n)%GiveForces, 3*Tdomain%sComm(n)%ngll, &
                    MPI_DOUBLE_PRECISION, other, tag_sl, Tdomain%communicateur, req_s_f(n+1), ierr)
                call MPI_Irecv(Tdomain%sComm(n)%TakeForces, 3*Tdomain%sComm(n)%ngll, &
                    MPI_DOUBLE_PRECISION, other, tag_sl, Tdomain%communicateur, req_r_f(n+1), ierr)
            endif
            if (Tdomain%sComm(n)%ngll_F/=0) then
                call MPI_Isend(Tdomain%sComm(n)%GiveForcesFl, 1*Tdomain%sComm(n)%ngll_F, &
                    MPI_DOUBLE_PRECISION, other, tag_fl, Tdomain%communicateur, req_s_fl(n+1), ierr)
                call MPI_Irecv(Tdomain%sComm(n)%TakeForcesFL,1*Tdomain%sComm(n)%ngll_F, &
                    MPI_DOUBLE_PRECISION, other, tag_fl, Tdomain%communicateur, req_r_fl(n+1), ierr)
            endif
            if (Tdomain%sComm(n)%ngllPML/=0) then
                call MPI_Isend(Tdomain%sComm(n)%GiveForcesPML, 9*Tdomain%sComm(n)%ngllPML, &
                    MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, req_s_pml(n+1), ierr)
                call MPI_Irecv(Tdomain%sComm(n)%TakeForcesPML, 9*Tdomain%sComm(n)%ngllPML, &
                    MPI_DOUBLE_PRECISION, other, tag_pml, Tdomain%communicateur, req_r_pml(n+1), ierr)
            endif
            if (Tdomain%sComm(n)%ngllPML_F/=0) then
                call MPI_Isend(Tdomain%sComm(n)%GiveForcesPMLFl, 3*Tdomain%sComm(n)%ngllPML_F, &
                    MPI_DOUBLE_PRECISION, other, tag_fpml, Tdomain%communicateur, req_s_fpml(n+1), ierr)
                call MPI_Irecv(Tdomain%sComm(n)%TakeForcesPMLFl, 3*Tdomain%sComm(n)%ngllPML_F, &
                    MPI_DOUBLE_PRECISION, other, tag_fpml, Tdomain%communicateur, req_r_fpml(n+1), ierr)
            endif
        enddo
        !write(*,*) "WAIT Exchange sem forces", rg

        call MPI_Waitall(Tdomain%tot_comm_proc, req_s_f, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_r_f, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_s_fl, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_r_fl, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_s_pml, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_r_pml, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_s_fpml, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_r_fpml, statuses, ierr)
        !write(*,*) "END Exchange sem forces", rg
    end subroutine exchange_sem_forces

    subroutine exchange_sem_forces_StoF(Tdomain)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: rg
        integer :: n
        integer :: other, ierr
        integer, dimension(Tdomain%tot_comm_proc) :: req_s, req_r, req_s_pml,req_r_pml
        integer, parameter :: tag_sf = 301, tag_sf_pml = 302
        integer, dimension(MPI_STATUS_SIZE,Tdomain%tot_comm_proc) :: statuses
        !- now we can exchange (communication global arrays)
        rg = Tdomain%rank
        req_s = MPI_REQUEST_NULL
        req_r = MPI_REQUEST_NULL
        req_s_pml = MPI_REQUEST_NULL
        req_r_pml = MPI_REQUEST_NULL


        do n = 0,Tdomain%tot_comm_proc-1
            other = Tdomain%sComm(n)%dest
            if (other == rg) cycle
            if (Tdomain%sComm(n)%ngllSF > 0) then
                call MPI_Isend(Tdomain%sComm(n)%GiveForcesSF_StoF,Tdomain%sComm(n)%ngllSF, &
                     MPI_DOUBLE_PRECISION, other, tag_sf, Tdomain%communicateur,req_s(n+1), ierr)
                call MPI_Irecv(Tdomain%sComm(n)%TakeForcesSF_StoF,Tdomain%sComm(n)%ngllSF, &
                     MPI_DOUBLE_PRECISION, other, tag_sf,Tdomain%communicateur,req_r(n+1), ierr)
            endif
            if (Tdomain%sComm(n)%ngllSF_PML > 0) then
                call MPI_Isend(Tdomain%sComm(n)%GiveForcesSF_StoF_PML,3*Tdomain%sComm(n)%ngllSF_PML, &
                     MPI_DOUBLE_PRECISION, other, tag_sf_pml,Tdomain%communicateur,req_s_pml(n+1), ierr)
                call MPI_Irecv(Tdomain%sComm(n)%TakeForcesSF_StoF_PML,3*Tdomain%sComm(n)%ngllSF_PML, &
                     MPI_DOUBLE_PRECISION, other,tag_sf_pml,Tdomain%communicateur,req_r_pml(n+1), ierr)
            endif
        enddo

        call MPI_Waitall(Tdomain%tot_comm_proc, req_s, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_r, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_s_pml, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_r_pml, statuses, ierr)

    end subroutine exchange_sem_forces_StoF

    subroutine exchange_sem_forces_FtoS(Tdomain)
        use sdomain
        use mpi
        implicit none
        type(domain), intent(inout) :: Tdomain
        integer :: rg
        integer :: n
        integer :: other, ierr
        integer, dimension(Tdomain%tot_comm_proc) :: req_s, req_r, req_s_pml,req_r_pml
        integer, parameter :: tag = 101
        integer, dimension(MPI_STATUS_SIZE,Tdomain%tot_comm_proc) :: statuses
        !- now we can exchange (communication global arrays)
        rg = Tdomain%rank
        req_s = MPI_REQUEST_NULL
        req_r = MPI_REQUEST_NULL
        req_s_pml = MPI_REQUEST_NULL
        req_r_pml = MPI_REQUEST_NULL

        do n = 0,Tdomain%tot_comm_proc-1
            other = Tdomain%sComm(n)%dest
            if (other == rg) cycle
            if (Tdomain%sComm(n)%ngllSF > 0) then
                call MPI_Isend(Tdomain%sComm(n)%GiveForcesSF_FtoS,3*Tdomain%sComm(n)%ngllSF, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur,req_s(n+1), ierr)
                call MPI_Irecv(Tdomain%sComm(n)%TakeForcesSF_FtoS,3*Tdomain%sComm(n)%ngllSF, &
                    MPI_DOUBLE_PRECISION, other, tag, Tdomain%communicateur,req_r(n+1), ierr)
            endif
            if (Tdomain%sComm(n)%ngllSF_PML > 0) then
                call MPI_Isend(Tdomain%sComm(n)%GiveForcesSF_FtoS_PML,9*Tdomain%sComm(n)%ngllSF_PML, &
                     MPI_DOUBLE_PRECISION, other, tag,Tdomain%communicateur,req_s_pml(n+1), ierr)
                call MPI_Irecv(Tdomain%sComm(n)%TakeForcesSF_FtoS_PML,9*Tdomain%sComm(n)%ngllSF_PML, &
                     MPI_DOUBLE_PRECISION, other, tag,Tdomain%communicateur,req_r_pml(n+1), ierr)
            endif

        enddo

        call MPI_Waitall(Tdomain%tot_comm_proc, req_s, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_r, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_s_pml, statuses, ierr)
        call MPI_Waitall(Tdomain%tot_comm_proc, req_r_pml, statuses, ierr)

    end subroutine exchange_sem_forces_FtoS


end module scomm

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
