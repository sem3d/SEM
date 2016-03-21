!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>


module mIntset
    type :: IntSet
        integer :: val
        type(IntSet), pointer :: next
    end type IntSet

contains

    subroutine free_intset(set)
        type(IntSet), pointer, intent(inout) :: set
        !
        type(IntSet), pointer :: temp
        do while(associated(set))
            temp => set%next
            deallocate(set)
            set => temp
        end do
    end subroutine free_intset

    function intersect_arrays(n1,v1,n2,v2,res) result(count)
        integer, dimension(n1), intent(in) :: v1
        integer, dimension(n2), intent(in) :: v2
        integer, dimension(:), intent(out), allocatable :: res
        !
        type(IntSet), pointer :: set, temp
        integer :: count, i, j

        nullify(set)
        count = 0
        do i=1,n1
            do j=1,n2
                if (v1(i)==v2(j)) then
                    allocate(temp)
                    temp%val = v1(i)
                    temp%next=>set
                    set=>temp
                    count = count +1
                    exit
                endif
            end do
        end do
        allocate(res(0:count-1))
        count = 0
        do while(associated(set))
            res(count) = set%val
            temp => set
            set=>set%next
            count = count+1
            deallocate(temp)
        end do
    end function intersect_arrays

end module mIntset


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
