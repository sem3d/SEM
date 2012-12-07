subroutine get_Forces_Elem2Vertex (Tdomain,n)

use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n

integer :: ngllx,nglly,ngllz,i,nv
 

ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

    do i = 0,7
        nv = Tdomain%specel(n)%Near_Vertices(i)
        select case (i)
         case (0)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(0,0,0,0:2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(0,0,0,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(0,0,0,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(0,0,0,0:2)
            endif
         case (1)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(ngllx-1,0,0,0:2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(ngllx-1,0,0,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(ngllx-1,0,0,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(ngllx-1,0,0,0:2)
            endif
         case (2)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(ngllx-1,nglly-1,0,0:2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(ngllx-1,nglly-1,0,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(ngllx-1,nglly-1,0,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(ngllx-1,nglly-1,0,0:2)
            endif
         case (3)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(0,nglly-1,0,0:2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(0,nglly-1,0,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(0,nglly-1,0,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(0,nglly-1,0,0:2)
            endif
         case (4)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(0,0,ngllz-1,0:2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(0,0,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(0,0,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(0,0,ngllz-1,0:2)
            endif
         case (5)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(ngllx-1,0,ngllz-1,0:2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(ngllx-1,0,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(ngllx-1,0,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(ngllx-1,0,ngllz-1,0:2)
            endif
         case (6)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(ngllx-1,nglly-1,ngllz-1,0:2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(ngllx-1,nglly-1,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(ngllx-1,nglly-1,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(ngllx-1,nglly-1,ngllz-1,0:2)
            endif
         case (7)
            Tdomain%sVertex(nv)%Forces(0:2) = Tdomain%sVertex(nv)%Forces(0:2) + &
                                              Tdomain%specel(n)%Forces(0,nglly-1,ngllz-1,0:2)
            if (Tdomain%sVertex(nv)%PML .and. Tdomain%specel(n)%PML) then
             Tdomain%sVertex(nv)%Forces1(0:2) = Tdomain%sVertex(nv)%Forces1(0:2) + &
                                              Tdomain%specel(n)%Forces1(0,nglly-1,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces2(0:2) = Tdomain%sVertex(nv)%Forces2(0:2) + &
                                              Tdomain%specel(n)%Forces2(0,nglly-1,ngllz-1,0:2)
             Tdomain%sVertex(nv)%Forces3(0:2) = Tdomain%sVertex(nv)%Forces3(0:2) + &
                                              Tdomain%specel(n)%Forces3(0,nglly-1,ngllz-1,0:2)
            endif
        end select
    enddo


return
end subroutine get_Forces_Elem2Vertex
