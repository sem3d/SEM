!>
!!\file check_inputs_and_mesh.F90
!!\brief Contains the subroutine check_inputs_and_mesh() which
!! checks in the inputs entered in the input file are not contradictory
!! and are already implemented. If not, the program terminates.
!!\author S. TERRANA
!!\version 1.0
!!\date 10/02/2014
!!
!<

subroutine check_inputs_and_mesh(Tdomain)
    use sdomain
    use constants

    implicit none

    type(domain), intent(inout)  :: Tdomain

    if (Tdomain%Implicitness .EQ. TIME_INTEG_IMPLICIT) then
        STOP "ERROR : Implicit methods not yet implemented. Please choose explicit|semi_implicit"
    endif

    if((Tdomain%Implicitness .EQ. TIME_INTEG_IMPLICIT) .AND. .NOT. &
      ((Tdomain%type_timeInteg .NE. TIME_INTEG_MIDPOINT_ITER) .OR. &
       (Tdomain%type_timeInteg .NE. TIME_INTEG_MIDPOINT)))  then
        STOP "ERROR : Semi-Implicit methods are implemented only for midpoint time integration."
    endif

    if ((Tdomain%type_timeInteg .NE. TIME_INTEG_NEWMARK)     .AND. &
        (Tdomain%type_timeInteg .NE. TIME_INTEG_RK4)         .AND. &
        (Tdomain%type_timeInteg .NE. TIME_INTEG_MIDPOINT) .AND. &
        (Tdomain%type_timeInteg .NE. TIME_INTEG_MIDPOINT_ITER)) then
        WRITE (*,*) "This choice for time integration does not exist !!!!!"
        STOP "Please choose RK4, Newmark, Midpoint, or Midpoint_iter"
    endif

    if ((Tdomain%type_timeInteg .EQ. TIME_INTEG_NEWMARK) .AND. &
        (Tdomain%type_elem .NE. GALERKIN_CONT)) then
        STOP "Error : For Newmark, only continuous Galerkin is available"
    endif

    if ((Tdomain%type_timeInteg .EQ. TIME_INTEG_MIDPOINT_ITER) .AND. &
        (Tdomain%type_elem .NE. GALERKIN_HDG_RP)) then
        WRITE (*,*) "This choice of element for iterative Midpoint is not available"
        !STOP "Error : Please choose dg_type = hdg_rp if you want to use Midpoint_iter"
    endif

    if ((Tdomain%logicD%post_proc) .AND. (Tdomain%type_elem .NE. GALERKIN_HDG_RP)) then
        WRITE (*,*) "Post-processing is available for HDG only !!!!!"
        STOP "Please choose dg_type = hdg_rp or disable the post_processing option."
    endif

    if ((Tdomain%logicD%post_proc) .AND. (Tdomain%capt_loc_type == CAPT_NEAREST_NODE)) then
        WRITE (*,*) "Cannot Post-process AND relocalize receiver... results will be wrong !!!!!"
        STOP "Please choose capt_loc_type = 0 or disable the post_processing option."
    endif

    if ((Tdomain%type_timeInteg .EQ. TIME_INTEG_MIDPOINT) .AND. &
        (Tdomain%type_elem .NE. GALERKIN_HDG_RP)) then
        WRITE (*,*) "This choice of element for the Midpoint time integration is not available"
        STOP "Error : Please choose dg_type = hdg_rp if you want to use Midpoint"
    endif

    if ((Tdomain%type_bc==DG_BC_FREE) .AND. (Tdomain%type_flux==FLUX_CENTERED)) then
        STOP "FREE Surface not implemented for Centered flux"
    endif

    if ((Tdomain%type_elem == GALERKIN_CONT) .AND. (Tdomain%type_flux .NE. FLUX_NONE)) then
        STOP "If you want to use Continuous Galerkin elements, please choose Flux NONE"
    endif

    if ((Tdomain%type_bc.NE.DG_BC_FREE) .AND. &
        (Tdomain%type_bc.NE.DG_BC_ABS) .AND. &
        (Tdomain%type_bc.NE.DG_BC_REFL))then
        STOP "This choice for boundary condition does not exist. Please choose free or absorbing or reflexing"
    endif

    if ((Tdomain%type_elem.NE.GALERKIN_CONT)     .AND. &
        (Tdomain%type_elem.NE.GALERKIN_DG_STRONG).AND. &
        (Tdomain%type_elem.NE.GALERKIN_DG_WEAK)  .AND. &
        (Tdomain%type_elem.NE.GALERKIN_HDG_RP)) then
        STOP "This choice for elem_type does not exist. Choose : continuous, dg-strong, dg-weak, or hdg_rp"
    endif

    ! Check the regularity of the mesh :
    call check_mesh(Tdomain)

end subroutine check_inputs_and_mesh


!>
!!\brief Subroutine check_mesh() checks the right orientation
!! of each element of the mesh and its convexity.
!! and are already implemented. If not, the program terminates.
!!\author S. TERRANA
!!\version 1.0
!!\date 10/02/2014
!!
!<
subroutine check_mesh(Tdomain)
    use sdomain
    use constants
    implicit none

    type(domain), intent(inout)  :: Tdomain
    real(fpp), dimension(0:1)    :: vect1, vect2, vect3
    integer, dimension (0:3)     :: Vertices
    integer  :: n, nv0, nv1, nv2, nv3
    real(fpp):: crossp1, crossp2

    do n=0,Tdomain%n_elem-1
        Vertices = Tdomain%specel(n)%Near_Vertex
        nv0 = Tdomain%sVertex(Vertices(0))%Glob_Numbering
        nv1 = Tdomain%sVertex(Vertices(1))%Glob_Numbering
        nv2 = Tdomain%sVertex(Vertices(2))%Glob_Numbering
        nv3 = Tdomain%sVertex(Vertices(3))%Glob_Numbering

        vect1(0) = Tdomain%coord_nodes(0,nv1) - Tdomain%coord_nodes(0,nv0)
        vect1(1) = Tdomain%coord_nodes(1,nv1) - Tdomain%coord_nodes(1,nv0)
        vect2(0) = Tdomain%coord_nodes(0,nv2) - Tdomain%coord_nodes(0,nv1)
        vect2(1) = Tdomain%coord_nodes(1,nv2) - Tdomain%coord_nodes(1,nv1)
        vect3(0) = Tdomain%coord_nodes(0,nv3) - Tdomain%coord_nodes(0,nv2)
        vect3(1) = Tdomain%coord_nodes(1,nv3) - Tdomain%coord_nodes(1,nv2)

        crossp1 = Vect1(0)*Vect2(1) - Vect1(1)*Vect2(0)
        crossp2 = Vect2(0)*Vect3(1) - Vect2(1)*Vect3(0)

        if (crossp1 == 0. .OR. crossp2 == 0.) then
            write(*,*) "Element number ", n, " is either a flat element has an edge with length 0"
            STOP "ERROR IN MESH : Element flat or degerated !"
        endif
        if ((crossp1 .LT. 0.) .OR. (crossp2 .LT. 0.)) then
            if ((crossp1 .LT. 0.) .AND. (crossp2 .LT. 0.)) then
                write(*,*) "Element number ", n, " is not properly oriented (bad node local numbering)"
                STOP "ERROR IN MESH : Element not properly oriented !"
            else
                write(*,*) "Element number ", n, " is a concave-quad or a cross-quad"
                STOP "ERROR IN MESH : Element concave or crossed!"
            endif
        endif
    enddo

end subroutine check_mesh

!>
!!\brief Subroutine check_modified_newmark() stops the run if
!! newmark_modified is requested on a mesh/config it does not support.
!! The modified-equation scheme (NewmarkModified.F90) is a linear,
!! matrix-free eigenvalue-based extension of EXPLICIT central-difference
!! Newmark (beta=0, gamma=0.5): its k=1 term is only mathematically
!! equivalent to one classic-Newmark step under that specific
!! parametrization (see NewmarkModified.F90's header), so any other
!! alpha/beta/gamma silently gives wrong physics rather than an error at
!! the numerics level - must be caught here. Also incompatible with PML
!! memory variables (real, time-persistent convolution state, not a
!! stateless linear operator), fault/friction nonlinearity, and
!! solid-fluid interfaces (ApplyA_global does not include the sf_stof/
!! sf_ftos coupling terms). Must run after PML_definition (which sets
!! specel%PML), so it is called separately from check_inputs_and_mesh
!! (which runs before the mesh's PML flags are known).
!<
subroutine check_modified_newmark(Tdomain)
    use sdomain
    implicit none

    type(domain), intent(in) :: Tdomain
    integer :: n
    logical :: has_pml, has_sf_iface

    if (.not. Tdomain%TimeD%modified) return

    ! Regional selection (select_modified_region, irons_dtcrit.F90) flags shared
    ! faces/vertices per-rank from local element data only; a partition-boundary
    ! face could then get inconsistent %modified on each side, desyncing that
    ! shared DOF across ranks. Not handled yet -- restrict to single-rank runs.
    if (Tdomain%Mpi_var%n_proc > 1) then
        STOP "ERROR : newmark_modified=true (regional selection) does not support MPI (n_proc>1) yet -- run with a single process."
    endif

    if (Tdomain%TimeD%beta /= 0._fpp .or. Tdomain%TimeD%gamma /= 0.5_fpp) then
        STOP "ERROR : newmark_modified=true requires explicit central-difference Newmark (beta=0, gamma=0.5) in time_scheme - its k=1 term is only equivalent to classic Newmark under that parametrization."
    endif

    has_pml = .false.
    do n = 0, Tdomain%n_elem - 1
        if (Tdomain%specel(n)%PML) has_pml = .true.
    enddo

    if (has_pml) then
        STOP "ERROR : newmark_modified=true does not support PML/CPML elements. Remove PML from the mesh or disable newmark_modified."
    endif

    if (Tdomain%n_fault > 0) then
        STOP "ERROR : newmark_modified=true does not support fault elements. Remove faults from the mesh or disable newmark_modified."
    endif

    has_sf_iface = .false.
    do n = 0, Tdomain%n_face - 1
        if (Tdomain%sFace(n)%is_sf_iface) has_sf_iface = .true.
    enddo

    if (has_sf_iface) then
        STOP "ERROR : newmark_modified=true does not support solid-fluid interfaces. Use a single-medium mesh or disable newmark_modified."
    endif
end subroutine check_modified_newmark

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
