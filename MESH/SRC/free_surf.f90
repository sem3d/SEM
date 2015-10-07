!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module free_surf_fluid

    !- finding elements with faces, edges, vertices at the free surface in fluids:
    !- will allow to chose the desired boundary condition there:
    !   1- default = wall (Neumann boundary condition with normal velocity = 0)
    !   2- pressure equals 0 = Dirichlet condition
    !   3- Neumann boundary condition: normal velocity imposed

    ! ATM: Any fluid element with a free surface is marked as having condition 2

    !- these fluid elements are then positively flagged => attribution for face, edges,
    !    vertices in the SEM run.

    private  :: face2corner

contains

    subroutine find_fluid_elem_freesurf(n_elem,dxadj,dxadjncy,Ipointer,elem_solid,elem_fluid_free)

        implicit none
        integer, intent(in)  :: n_elem
        integer, dimension(0:), intent(in)  :: dxadj,dxadjncy
        logical, dimension(0:), intent(in)  :: elem_solid
        logical, dimension(0:), intent(out)  :: elem_fluid_free
        integer, dimension(0:,0:), intent(in) :: Ipointer
        integer :: i,j,k,n,nel,num,ok,neighbor,ind_elem
        integer, dimension(0:3)  ::  corner

        find_el: do nel = 0,n_elem-1
            elem_fluid_free(nel) = .false.
            if(elem_solid(nel) .or. (dxadj(nel+1)- dxadj(nel) == 6))then
                cycle find_el   ! solid element, or each face sees another element
            end if
            ! now we look if the face #5 is free
            call face2corner(Ipointer(0:,nel),5,corner)
            find0 : do i = dxadj(nel),dxadj(nel+1)-1
                neighbor = dxadjncy(i)
                find1 : do j = 0,3   ! does the face belong to this neighbour?
                    num = corner(j)
                    ok = 0
                    find2 : do k = 0,7
                        if(Ipointer(k,neighbor) == num)then
                            ok = 1
                            exit find2
                        endif
                    enddo find2
                    if(ok == 0) exit find1   ! face not shared with this neighbor
                    if(j == 3) then
                        cycle find_el  ! a neighbor found: not a free surface
                    end if
                end do find1
            end do find0
            ! no neighbor for the face #5: it'a a free face, Dirichlet condition
            elem_fluid_free(nel) = .true.
        end do find_el

    end subroutine find_fluid_elem_freesurf

    !-------------------------------------
    subroutine face2corner(Ipoint,nfa,corn)
        implicit none
        integer, intent(in)   :: Ipoint(0:7),nfa
        integer, intent(out)  :: corn(0:3)

        select case(nfa)
        case(0)
            corn(0) = Ipoint(0)
            corn(1) = Ipoint(1)
            corn(2) = Ipoint(2)
            corn(3) = Ipoint(3)
        case(1)
            corn(0) = Ipoint(0)
            corn(1) = Ipoint(1)
            corn(2) = Ipoint(5)
            corn(3) = Ipoint(4)
        case(2)
            corn(0) = Ipoint(1)
            corn(1) = Ipoint(2)
            corn(2) = Ipoint(6)
            corn(3) = Ipoint(5)
        case(3)
            corn(0) = Ipoint(3)
            corn(1) = Ipoint(2)
            corn(2) = Ipoint(6)
            corn(3) = Ipoint(7)
        case(4)
            corn(0) = Ipoint(0)
            corn(1) = Ipoint(3)
            corn(2) = Ipoint(7)
            corn(3) = Ipoint(4)
        case(5)
            corn(0) = Ipoint(4)
            corn(1) = Ipoint(5)
            corn(2) = Ipoint(6)
            corn(3) = Ipoint(7)
        end select

    end subroutine face2corner
    !---------------------------------------


end module free_surf_fluid

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
