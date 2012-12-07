subroutine get_Mass_Elem2Vertex (Tdomain,n)

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
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                      Tdomain%specel(n)%MassMat(0,0,0)
     case (1)
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                      Tdomain%specel(n)%MassMat(ngllx-1,0,0)
     case (2)
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                      Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,0)
     case (3)
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                      Tdomain%specel(n)%MassMat(0,nglly-1,0)
     case (4)
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                      Tdomain%specel(n)%MassMat(0,0,ngllz-1)
     case (5)
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                      Tdomain%specel(n)%MassMat(ngllx-1,0,ngllz-1)
     case (6)
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                      Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,ngllz-1)
     case (7)
        Tdomain%sVertex(nv)%MassMat = Tdomain%sVertex(nv)%MassMat + &
                                      Tdomain%specel(n)%MassMat(0,nglly-1,ngllz-1)
    end select

    if (Tdomain%sVertex(nv)%PML) then

        select case (i)
         case (0)
            Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,0,0,:)
         if (Tdomain%sVertex(nv)%FPML) then
                 Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) + &
                                          Tdomain%specel(n)%Ivx(0,0,0)
                 Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) + &
                                          Tdomain%specel(n)%Ivy(0,0,0)
                 Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) + &
                                          Tdomain%specel(n)%Ivz(0,0,0)
        endif
         case (1)
            Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,0,0,:)
         if (Tdomain%sVertex(nv)%FPML) then
                 Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) + &
                                          Tdomain%specel(n)%Ivx(ngllx-1,0,0)
                 Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) + &
                                          Tdomain%specel(n)%Ivy(ngllx-1,0,0)
                 Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) + &
                                          Tdomain%specel(n)%Ivz(ngllx-1,0,0)
        endif
         case (2)
            Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,0,:)
         if (Tdomain%sVertex(nv)%FPML) then
                 Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) + &
                                          Tdomain%specel(n)%Ivx(ngllx-1,nglly-1,0)
                 Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) + &
                                          Tdomain%specel(n)%Ivy(ngllx-1,nglly-1,0)
                 Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) + &
                                          Tdomain%specel(n)%Ivz(ngllx-1,nglly-1,0)
        endif
         case (3)
            Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,nglly-1,0,:)
         if (Tdomain%sVertex(nv)%FPML) then
                 Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) + &
                                          Tdomain%specel(n)%Ivx(0,nglly-1,0)
                 Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) + &
                                          Tdomain%specel(n)%Ivy(0,nglly-1,0)
                 Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) + &
                                          Tdomain%specel(n)%Ivz(0,nglly-1,0)
        endif
         case (4)
            Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,0,ngllz-1,:)
         if (Tdomain%sVertex(nv)%FPML) then
                 Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) + &
                                          Tdomain%specel(n)%Ivx(0,0,ngllz-1)
                 Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) + &
                                          Tdomain%specel(n)%Ivy(0,0,ngllz-1)
                 Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) + &
                                          Tdomain%specel(n)%Ivz(0,0,ngllz-1)
        endif
         case (5)
            Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,0,ngllz-1,:)
          if (Tdomain%sVertex(nv)%FPML) then
                 Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) + &
                                          Tdomain%specel(n)%Ivx(ngllx-1,0,ngllz-1)
                 Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) + &
                                          Tdomain%specel(n)%Ivy(ngllx-1,0,ngllz-1)
                 Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) + &
                                          Tdomain%specel(n)%Ivz(ngllx-1,0,ngllz-1)
        endif
         case (6)
            Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,ngllz-1,:)
         if (Tdomain%sVertex(nv)%FPML) then
                 Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) + &
                                          Tdomain%specel(n)%Ivx(ngllx-1,nglly-1,ngllz-1)
                 Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) + &
                                          Tdomain%specel(n)%Ivy(ngllx-1,nglly-1,ngllz-1)
                 Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) + &
                                          Tdomain%specel(n)%Ivz(ngllx-1,nglly-1,ngllz-1)
        endif
         case (7)
            Tdomain%sVertex(nv)%DumpMass(:) = Tdomain%sVertex(nv)%DumpMass(:) + &
                                              Tdomain%specel(n)%DumpMass(0,nglly-1,ngllz-1,:)
         if (Tdomain%sVertex(nv)%FPML) then
                 Tdomain%sVertex(nv)%Ivx(0) = Tdomain%sVertex(nv)%Ivx(0) + &
                                          Tdomain%specel(n)%Ivx(0,nglly-1,ngllz-1)
                 Tdomain%sVertex(nv)%Ivy(0) = Tdomain%sVertex(nv)%Ivy(0) + &
                                          Tdomain%specel(n)%Ivy(0,nglly-1,ngllz-1)
                 Tdomain%sVertex(nv)%Ivz(0) = Tdomain%sVertex(nv)%Ivz(0) + &
                                          Tdomain%specel(n)%Ivz(0,nglly-1,ngllz-1)
        endif
        end select
    endif

enddo


return
end subroutine get_Mass_Elem2Vertex
