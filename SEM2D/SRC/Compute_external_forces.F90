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
    real, intent (INOUT)          :: timelocal
    real, dimension(0:1)          :: Fext
    integer  :: n, ns, ncc, ngllx, ngllz, i, j, np, nDG

    do n = 0, Tdomain%n_source-1
        do ns =0, Tdomain%sSource(n)%ine-1
            ncc = Tdomain%sSource(n)%Elem(ns)%nr
            ngllx = Tdomain%specel(ncc)%ngllx; ngllz = Tdomain%specel(ncc)%ngllz
            if ( Tdomain%specel(ncc)%Type_DG == GALERKIN_CONT) then
                nDG = 0
            else ! In DG case, the forces has to be put in Forces(:,:,3:4)
                nDG = 3
            endif
            if (Tdomain%sSource(n)%i_type_source  == 6) then ! Special sources on the strain equation.
                do j = 0,ngllz-1
                    do i = 0,ngllx-1
                        do np = 0,2
                            Fext(np) = CompSource (Tdomain%sSource(n),timelocal) &
                                     * Tdomain%sSource(n)%Elem(ns)%ExtForce(i,j,np)
                            Tdomain%specel(ncc)%Forces(i,j,np) = Tdomain%specel(ncc)%Forces(i,j,np) + Fext(np)
                        enddo
                    enddo
                enddo
            else ! Usual sources added on the velocity equation.
                do j = 0,ngllz-1
                    do i = 0,ngllx-1
                        do np = 0,1
                            Fext(np) = CompSource (Tdomain%sSource(n),timelocal) &
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
