!>
!!\file Compute_external_forces.F90
!!\brief Computes the external Forces (i.e. the source of the seism)
!!\version 1.0
!!\date 20/11/2013
!! This algorithm contains the part of the previous Newmark routine
!! which was calculating the external forces.
!<
subroutine Compute_external_forces (Tdomain,timelocal)

    use sdomain
    use ssources
    use constants
    implicit none
    type (domain), intent (INOUT) :: Tdomain
    real(fpp), intent (INOUT)          :: timelocal
    real(fpp), dimension(0:1)          :: Fext
    real(fpp)                          :: srcval
    integer  :: n, ns, ncc, ngllx, ngllz, i, j, np, nDG, nmat0
    logical  :: use_time_integral

    do n = 0, Tdomain%n_source-1
        ! Time factor for this source at this step. Type 7 (pressure) injects the running
        ! time integral int(f dt) instead of f(t): p=-VelPhi and the source enters at the
        ! phi-acceleration level, so the integral makes the pressure equal f(t) (mirrors
        ! SEM3D Newmark.f90 type-7). Type 3 (fluidpulse) routed to the moment-tensor
        ! (displacement-vector) iso-acoustic path (source_excit_fluid) ALSO needs the
        ! integral: M0=kappa*deltaV requires the injected VOLUME (integral of the declared
        ! rate f(t)), not the rate itself, whereas the potential-fluid phi equation takes
        ! f(t) directly as a rate term -- see 2026-07-17 fluid-aniso velocity investigation.
        ! Accumulate ONCE per source per step, before the element loop.
        use_time_integral = (Tdomain%sSource(n)%i_type_source == 7)
        if (Tdomain%sSource(n)%i_type_source == 3 .and. Tdomain%sSource(n)%ine > 0) then
            nmat0 = Tdomain%specel(Tdomain%sSource(n)%Elem(0)%nr)%mat_index
            if (.not. (Tdomain%sSubdomain(nmat0)%deftype == MATDEF_FLUID_ANISO .or. &
                       Tdomain%sSubdomain(nmat0)%deftype == CSTAR_FLUID)) then
                use_time_integral = .true.
            endif
        endif
        if (use_time_integral) then
            Tdomain%sSource(n)%time_integral = Tdomain%sSource(n)%time_integral &
                + CompSource(Tdomain%sSource(n), timelocal) * Tdomain%TimeD%dtmin
            srcval = Tdomain%sSource(n)%time_integral
        else
            srcval = CompSource(Tdomain%sSource(n), timelocal)
        end if
        do ns =0, Tdomain%sSource(n)%ine-1
            ncc = Tdomain%sSource(n)%Elem(ns)%nr
            ngllx = Tdomain%specel(ncc)%ngllx; ngllz = Tdomain%specel(ncc)%ngllz
            if (Tdomain%specel(ncc)%Type_DG == GALERKIN_CONT) then
                nDG = 0
            else if (Tdomain%specel(ncc)%acoustic) then
                nDG = 1
            else! In DG case, the forces has to be put in Forces(:,:,3:4)
                nDG = 3
            endif
            if (Tdomain%sSource(n)%i_type_source  == 6) then ! Special sources on the strain equation.
                do j = 0,ngllz-1
                    do i = 0,ngllx-1
                        do np = 0,1
                            Fext(np) = srcval &
                                     * Tdomain%sSource(n)%Elem(ns)%ExtForce(i,j,np)
                            Tdomain%specel(ncc)%Forces(i,j,np) = Tdomain%specel(ncc)%Forces(i,j,np) + Fext(np)
                        enddo
                    enddo
                enddo
            else ! Usual sources added on the velocity equation.
                do j = 0,ngllz-1
                    do i = 0,ngllx-1
                        do np = 0,1
                            Fext(np) = srcval &
                                     * Tdomain%sSource(n)%Elem(ns)%ExtForce(i,j,np)
                        Tdomain%specel(ncc)%Forces(i,j,np+nDG) = Tdomain%specel(ncc)%Forces(i,j,np+nDG) + Fext(np)
                        enddo
                    enddo
                enddo
            endif
        enddo
    enddo

end subroutine Compute_external_forces
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
