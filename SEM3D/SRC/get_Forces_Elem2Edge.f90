subroutine get_Forces_Elem2Edge (Tdomain,n)

use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n

integer :: ngllx,nglly,ngllz,ngll,i,j,ne
 

ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

do i = 0,11
    ne = Tdomain%specel(n)%Near_Edges(i)
    ngll = Tdomain%sEdge(ne)%ngll

    if ( Tdomain%specel(n)%Orient_Edges(i) == 0 ) then

        select case (i)
         case (0)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(1:ngllx-2,0,0,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(1:ngllx-2,0,0,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(1:ngllx-2,0,0,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(1:ngllx-2,0,0,0:2)
            endif
          case (1)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,0,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(ngllx-1,1:nglly-2,0,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(ngllx-1,1:nglly-2,0,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(ngllx-1,1:nglly-2,0,0:2)
            endif
         case (2)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,0,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(1:ngllx-2,nglly-1,0,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(1:ngllx-2,nglly-1,0,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(1:ngllx-2,nglly-1,0,0:2)
            endif
         case (3)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(0,1:nglly-2,0,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(0,1:nglly-2,0,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(0,1:nglly-2,0,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(0,1:nglly-2,0,0:2)
            endif
         case (4)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(ngllx-1,0,1:ngllz-2,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(ngllx-1,0,1:ngllz-2,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(ngllx-1,0,1:ngllz-2,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(ngllx-1,0,1:ngllz-2,0:2)
            endif
         case (5)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(1:ngllx-2,0,ngllz-1,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(1:ngllx-2,0,ngllz-1,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(1:ngllx-2,0,ngllz-1,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(1:ngllx-2,0,ngllz-1,0:2)
            endif
         case (6)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(0,0,1:ngllz-2,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(0,0,1:ngllz-2,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(0,0,1:ngllz-2,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(0,0,1:ngllz-2,0:2)
            endif
         case (7)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(ngllx-1,nglly-1,1:ngllz-2,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(ngllx-1,nglly-1,1:ngllz-2,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(ngllx-1,nglly-1,1:ngllz-2,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(ngllx-1,nglly-1,1:ngllz-2,0:2)
            endif
         case (8)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,ngllz-1,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(ngllx-1,1:nglly-2,ngllz-1,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(ngllx-1,1:nglly-2,ngllz-1,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(ngllx-1,1:nglly-2,ngllz-1,0:2)
            endif
         case (9)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,ngllz-1,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(1:ngllx-2,nglly-1,ngllz-1,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(1:ngllx-2,nglly-1,ngllz-1,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(1:ngllx-2,nglly-1,ngllz-1,0:2)
            endif
         case (10)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(0,nglly-1,1:ngllz-2,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(0,nglly-1,1:ngllz-2,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(0,nglly-1,1:ngllz-2,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(0,nglly-1,1:ngllz-2,0:2)
            endif
         case (11)
            Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces(0,1:nglly-2,ngllz-1,0:2)
            if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces1(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces1(0,1:nglly-2,ngllz-1,0:2)
             Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces2(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces2(0,1:nglly-2,ngllz-1,0:2)
             Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) = Tdomain%sEdge(ne)%Forces3(1:ngll-2,0:2) + &
                                                     Tdomain%specel(n)%Forces3(0,1:nglly-2,ngllz-1,0:2)
            endif
        end select
    
    else 

    select case (i)
     case (0)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(ngllx-1-j,0,0,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(ngllx-1-j,0,0,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(ngllx-1-j,0,0,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(ngllx-1-j,0,0,0:2)
        enddo
      endif
     case (1)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,0,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(ngllx-1,nglly-1-j,0,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(ngllx-1,nglly-1-j,0,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(ngllx-1,nglly-1-j,0,0:2)
        enddo
      endif
     case (2)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,0,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(ngllx-1-j,nglly-1,0,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(ngllx-1-j,nglly-1,0,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(ngllx-1-j,nglly-1,0,0:2)
        enddo
      endif
     case (3)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(0,nglly-1-j,0,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(0,nglly-1-j,0,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(0,nglly-1-j,0,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(0,nglly-1-j,0,0:2)
        enddo
      endif
     case (4)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1-j,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(ngllx-1,0,ngllz-1-j,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(ngllx-1,0,ngllz-1-j,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(ngllx-1,0,ngllz-1-j,0:2)
        enddo
      endif
     case (5)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(ngllx-1-j,0,ngllz-1,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(ngllx-1-j,0,ngllz-1,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(ngllx-1-j,0,ngllz-1,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(ngllx-1-j,0,ngllz-1,0:2)
        enddo
      endif
     case (6)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(0,0,ngllz-1-j,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(0,0,ngllz-1-j,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(0,0,ngllz-1-j,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(0,0,ngllz-1-j,0:2)
        enddo
      endif
     case (7)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1-j,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(ngllx-1,nglly-1,ngllz-1-j,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(ngllx-1,nglly-1,ngllz-1-j,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(ngllx-1,nglly-1,ngllz-1-j,0:2)
        enddo
      endif
     case (8)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,ngllz-1,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(ngllx-1,nglly-1-j,ngllz-1,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(ngllx-1,nglly-1-j,ngllz-1,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(ngllx-1,nglly-1-j,ngllz-1,0:2)
        enddo
      endif
     case (9)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,ngllz-1,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(ngllx-1-j,nglly-1,ngllz-1,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(ngllx-1-j,nglly-1,ngllz-1,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(ngllx-1-j,nglly-1,ngllz-1,0:2)
        enddo
      endif
     case (10)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1-j,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(0,nglly-1,ngllz-1-j,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(0,nglly-1,ngllz-1-j,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(0,nglly-1,ngllz-1-j,0:2)
        enddo
      endif
     case (11)
      do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces(j,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2) + &
                                       Tdomain%specel(n)%Forces(0,nglly-1-j,ngllz-1,0:2)
      enddo
      if (Tdomain%sEdge(ne)%PML .and. Tdomain%specel(n)%PML) then
        do j = 1,ngll-2
        Tdomain%sEdge(ne)%Forces1(j,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2) + &
                                       Tdomain%specel(n)%Forces1(0,nglly-1-j,ngllz-1,0:2)
        Tdomain%sEdge(ne)%Forces2(j,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2) + &
                                       Tdomain%specel(n)%Forces2(0,nglly-1-j,ngllz-1,0:2)
        Tdomain%sEdge(ne)%Forces3(j,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2) + &
                                       Tdomain%specel(n)%Forces3(0,nglly-1-j,ngllz-1,0:2)
        enddo
      endif
    end select

    endif

enddo


return
end subroutine get_Forces_Elem2Edge
