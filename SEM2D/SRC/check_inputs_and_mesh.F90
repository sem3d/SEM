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

    if ((Tdomain%type_timeInteg .NE. TIME_INTEG_NEWMARK) .AND. &
        (Tdomain%type_timeInteg .NE. TIME_INTEG_RK4)) then
        STOP "This choice for time_integ does not exist. Please choose RK4 or Newmark"
    endif

    if ((Tdomain%type_timeInteg .EQ. TIME_INTEG_NEWMARK) .AND. &
        (Tdomain%type_elem .NE. GALERKIN_CONT)) then
        STOP "Error : For Newmark, only continuous Galerkin is available"
    endif

    if ((Tdomain%type_timeInteg .EQ. TIME_INTEG_RK4) .AND. &
        (Tdomain%type_elem .EQ. GALERKIN_CONT)) then
        STOP "Error : For RK4, only discontinuous Galerkin methods are available"
    endif

    if ((Tdomain%type_bc==DG_BC_FREE) .AND. (Tdomain%type_flux==FLUX_CENTERED)) then
        STOP "FREE Surface not implemented for Centered flux"
    endif

    if ((Tdomain%type_bc.NE.DG_BC_FREE) .AND. (Tdomain%type_bc.NE.DG_BC_ABS)) then
        STOP "This choice for boundary condition does not exist. Please choose free or absorbing"
    endif

    if ((Tdomain%type_elem.NE.GALERKIN_CONT) .AND. &
        (Tdomain%type_elem.NE.GALERKIN_DG_STRONG) .AND. &
        (Tdomain%type_elem.NE.GALERKIN_DG_WEAK)) then
        STOP "This choice for elem_type does not exist. Choose : continuous, dg-strong, or dg-weak"
    endif

end subroutine check_inputs_and_mesh
