module free_surf_fluid

  !- finding elements with faces, edges, vertices at the free surface in fluids:
  !    pressure equals 0 = Dirichlet condition 

  !- in general: face #5 in a element, is oriented vertically outward; TO DO: check this (eventual
  !    problem when deformed fluid elements near SF interface)

  !- these fluid elements are then positively flagged => attribution for face, edges, 
  !    vertices in the SEM run.

    use mesh_properties

    contains

    subroutine find_fluid_elem_dirich(n_elem,dxadj,dxadjncy,Ipointer,elem_solid,elem_fluid_dirich)
    
        implicit none
        integer, intent(in)  :: n_elem
        integer, dimension(0:), intent(in)  :: dxadj,dxadjncy
        logical, dimension(0:), intent(in)  :: elem_solid
        logical, dimension(0:), intent(out)  :: elem_fluid_dirich
        integer, dimension(0:,0:), intent(in) :: Ipointer
        integer :: i,j,k,n,nel,num,ok,neighbor,ind_elem
        integer, dimension(0:3)  ::  corner


   find_el: do nel = 0,n_elem-1
            elem_fluid_dirich(nel) = .false.
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
              elem_fluid_dirich(nel) = .true.
          end do find_el

    end subroutine find_fluid_elem_dirich

end module free_surf_fluid

