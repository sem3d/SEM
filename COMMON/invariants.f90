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
    ! FIRST INVARIANT I1
    !****************************************************************************
    
    subroutine first_tensor_invariant(tensor,I1)
        ! intent IN
        real(fpp), intent(in), dimension(0:5) :: tensor
        ! intent INOUT
        real(fpp), intent(inout)              :: I1
        !
        I1 = dot_product(tensor,Mvector)
        !
        return
        !
    end subroutine first_tensor_invariant 

    !****************************************************************************
    ! DEVIATORIC TENSOR Sij
    !****************************************************************************
    
    subroutine tensor_deviator(tensor,deviator)
        ! intent IN
        real(fpp), intent(in), dimension(0:5)    :: tensor
        ! intent INOUT
        real(fpp), intent(inout), dimension(0:5) :: deviator
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

    !****************************************************************************
    ! SECOND TENSOR INVARIANT I2
    !****************************************************************************

    subroutine second_tensor_invariant(tensor,I2)
        ! intent IN
        real(fpp), intent(in), dimension(0:5) :: tensor
        ! intent INOUT
        real(fpp), intent(inout)              :: I2
        !
        I2 = tensor(0)*tensor(1)+tensor(1)*tensor(2)+tensor(2)*tensor(0)-&
            (tensor(3)**2+tensor(4)**2+tensor(5)**2)
        !
        return
        !
    end subroutine second_tensor_invariant

    !****************************************************************************
    ! THIRD TENSOR INVARIANT I3
    !****************************************************************************

    subroutine third_tensor_invariant(tensor,I3)
        ! intent IN
        real(fpp), intent(in), dimension(0:5) :: tensor
        ! intent INOUT
        real(fpp), intent(inout)              :: I3
        !
        I3 = tensor(0)*tensor(1)*tensor(2)+&
            2*tensor(3)*tensor(4)*tensor(5)-&
            tensor(2)*tensor(3)**2-&
            tensor(1)*tensor(4)**2-&
            tensor(0)*tensor(5)**2
        !
        return
        !
    end subroutine third_tensor_invariant

    !****************************************************************************
    ! SECOND DEVIATOR INVARIANT J2
    !****************************************************************************

    subroutine second_deviator_invariant(deviator,J2)
        ! input IN
        real(fpp), intent(in), dimension(0:5) :: deviator
        ! intent INOUT
        real(fpp), intent(inout)              :: J2
        !
        real(fpp), dimension(0:5) :: temp
        !
        temp = half*Avector*deviator
        J2 = dot_product(deviator,temp)
        !
        return
        !
    end subroutine second_deviator_invariant

    !****************************************************************************
    ! THIRD DEVIATOR INVARIANT J3
    !****************************************************************************

    subroutine third_deviator_invariant(tensor,J3)
        ! input IN
        real(fpp), intent(in), dimension(0:5) :: tensor
        ! intent INOUT
        real(fpp), intent(inout)              :: J3
        !
        real(fpp) :: I1,I2,I3
        !
        call first_tensor_invariant(tensor,I1)
        call second_tensor_invariant(tensor,I2)
        call third_tensor_invariant(tensor,I3)
        !
        J3 = two/27.0d0*I1**3-one/three*I1*I2+I3
        !
        return
        !
    end subroutine third_deviator_invariant

end module invariants

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=4 et tw=80 smartindent : !!
