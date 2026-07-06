!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file allocate_domain.f90
!!\brief Gére l'allocation des domaines.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Allocation des attributs de la structure domain
!!
!! \param type(domain), intent (INOUT) Tdomain
!<
module sdomain_alloc
    use sdomain
    use dom_solid
    use dom_solid_dg
    use dom_solidpml
    use dom_fluid
    use dom_fluidpml
    implicit none
contains

subroutine allocate_domain (Tdomain)

    type(domain), intent (INOUT) :: Tdomain
    integer :: n,ngll
    !integer :: mat, randSize, assocMat

    if(.false.) then
        write(*,*) "Tdomain%any_sdom = ", Tdomain%any_sdom
        write(*,*) "Tdomain%any_fdom = ", Tdomain%any_fdom
        write(*,*) "Tdomain%any_spml = ", Tdomain%any_spml
        write(*,*) "Tdomain%any_fpml = ", Tdomain%any_fpml
    end if
    allocate(Tdomain%sdom     %champs(0:Tdomain%TimeD%nsubsteps))
    allocate(Tdomain%sdomdg   %champs(0:Tdomain%TimeD%nsubsteps))
    allocate(Tdomain%spmldom  %champs(0:Tdomain%TimeD%nsubsteps))
    allocate(Tdomain%fdom      %champs(0:Tdomain%TimeD%nsubsteps))
    allocate(Tdomain%fpmldom   %champs(0:Tdomain%TimeD%nsubsteps))
    if(Tdomain%any_sdom)      call allocate_dom_solid      (Tdomain, Tdomain%sdom)
    if(Tdomain%any_sdomdg)    call allocate_dom_solid_dg   (Tdomain, Tdomain%sdomdg)
    if(Tdomain%any_fdom)      call allocate_dom_fluid      (Tdomain, Tdomain%fdom)
    if(Tdomain%any_spml)      call allocate_dom_solidpml   (Tdomain, Tdomain%spmldom)
    if(Tdomain%any_fpml)      call allocate_dom_fluidpml   (Tdomain, Tdomain%fpmldom)

    do n = 0,Tdomain%n_elem-1
        ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
        allocate(Tdomain%specel(n)%MassMat(0:ngll-1, 0:ngll-1, 0:ngll-1))
        allocate(Tdomain%specel(n)%En_S_avg(0:ngll-2, 0:ngll-2, 0:ngll-2))
        allocate(Tdomain%specel(n)%En_P_avg(0:ngll-2, 0:ngll-2, 0:ngll-2))

    enddo

    call report_domain_totals(Tdomain)
end subroutine allocate_domain

!! Reports, from rank 0, the model-wide (all-ranks) number of elements, GLL points and DOF
!! per domain type, plus a load-balance line (min/max of per-rank DOF). DOF/GLL differs by
!! domain (solid 3, fluid 1, PML/DG more) and reflects the relative computational weight.
subroutine report_domain_totals(Tdomain)
    use mpi
    type(domain), intent(in) :: Tdomain
    integer, parameter :: ND = 5
    ! order: 1 solid, 2 solid DG, 3 fluid (iso or aniso, see dom%aniso), 4 solid PML, 5 fluid PML
    integer(kind=8) :: loc(3,ND), tot(3,ND)              ! (nbelem, nglltot, ndof) per domain
    integer(kind=8) :: rk_ndof, sum_ndof, min_ndof, max_ndof
    integer         :: dof(ND), ierr, d, nprocs
    character(len=12) :: nm(ND)
    real(fpp)       :: avg
    dof = (/ 3, 9, 1, 9, 3 /)     ! DOF per GLL point (see champs_* field allocations)
    nm  = (/ 'solid       ', 'solid DG    ', 'fluid       ', &
             'solid PML   ', 'fluid PML   ' /)

    loc = 0
    if (Tdomain%any_sdom)      then; loc(1,1)=Tdomain%sdom%nbelem;      loc(2,1)=Tdomain%sdom%nglltot;      endif
    if (Tdomain%any_sdomdg)    then; loc(1,2)=Tdomain%sdomdg%nbelem;    loc(2,2)=Tdomain%sdomdg%nglltot;    endif
    if (Tdomain%any_fdom)      then; loc(1,3)=Tdomain%fdom%nbelem;      loc(2,3)=Tdomain%fdom%nglltot;      endif
    if (Tdomain%any_spml)      then; loc(1,4)=Tdomain%spmldom%nbelem;   loc(2,4)=Tdomain%spmldom%nglltot;   endif
    if (Tdomain%any_fpml)      then; loc(1,5)=Tdomain%fpmldom%nbelem;   loc(2,5)=Tdomain%fpmldom%nglltot;   endif
    do d = 1, ND
        loc(3,d) = loc(2,d) * int(dof(d), kind=8)
    end do
    rk_ndof = sum(loc(3,:))     ! this rank's total DOF (compute weight)

    call MPI_Reduce(loc,     tot,      3*ND, MPI_INTEGER8, MPI_SUM, 0, Tdomain%communicateur, ierr)
    call MPI_Reduce(rk_ndof, sum_ndof, 1,    MPI_INTEGER8, MPI_SUM, 0, Tdomain%communicateur, ierr)
    call MPI_Reduce(rk_ndof, min_ndof, 1,    MPI_INTEGER8, MPI_MIN, 0, Tdomain%communicateur, ierr)
    call MPI_Reduce(rk_ndof, max_ndof, 1,    MPI_INTEGER8, MPI_MAX, 0, Tdomain%communicateur, ierr)
    call MPI_Comm_size(Tdomain%communicateur, nprocs, ierr)

    if (Tdomain%rank == 0) then
        write(*,*)
        write(*,*) "===== Domain summary (all ranks) : elements / GLL points / DOF ====="
        do d = 1, ND
            if (tot(1,d) > 0) then
                write(*,'(A,A,I12,A,I14,A,I14,A,F6.2,A)') "  ", nm(d), &
                    tot(1,d), " elem", tot(2,d), " gll", tot(3,d), " dof (", &
                    100.0d0*real(tot(3,d),8)/max(sum_ndof,1_8), " % dof)"
            end if
        end do
        write(*,'(A,I12,A,I14,A,I14,A)') "  TOTAL       ", &
            sum(tot(1,:)), " elem", sum(tot(2,:)), " gll", sum_ndof, " dof"
        if (nprocs > 1) then
            avg = real(sum_ndof,8)/nprocs
            write(*,'(A,I14,A,I14,A,F14.1,A,F6.2)') "  DOF/rank  min=", min_ndof, &
                "  max=", max_ndof, "  avg=", avg, "  imbalance(max/avg)=", real(max_ndof,8)/max(avg,1.0d0)
        end if
        write(*,*)
    end if
end subroutine report_domain_totals

end module sdomain_alloc

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
