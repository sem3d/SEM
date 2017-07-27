!>
!!\file invariants.f90
!!\brief Contains subroutines to compute stress invariants
!!
!<

module invariants

    use constants

    implicit none

contains

    !****************************************************************************
    ! FIRST INVARIANT
    !****************************************************************************
    
    subroutine first_tensor_invariant(tensor,I1)
        ! intent IN
        real(fpp), dimension(0:5), intent(in) :: tensor
        ! intent INOUT
        real(fpp), intent(inout)              :: I1
        !
        I1 = dot_product(tensor,Mvector)
        !
        return
        !
    end subroutine first_tensor_invariant 

    !****************************************************************************
    ! DEVIATORIC TENSOR
    !****************************************************************************
    
    subroutine tensor_deviator(tensor,deviator)
        ! intent IN
        real(fpp), dimension(0:5), intent(in)           :: tensor
        ! intent INOUT
        real(fpp), dimension(0:5), intent(inout)        :: deviator
        !
        real(fpp) :: I1
        !
        call first_tensor_invariant(tensor,I1)
        I1 = I1/three
        !
        deviator(0:5) = tensor(0:5) - I1*Mvector(0:5)
        !
        return
        !
    end subroutine tensor_deviator 

end module invariants

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!

