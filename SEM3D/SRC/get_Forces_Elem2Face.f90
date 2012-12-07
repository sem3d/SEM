subroutine get_Forces_Elem2Face (Tdomain,n)

use sdomain

implicit none

type (Domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: n

integer :: ngllx,nglly,ngllz,ngll1,ngll2,i,j,k,nf
 

ngllx = Tdomain%specel(n)%ngllx
nglly = Tdomain%specel(n)%nglly
ngllz = Tdomain%specel(n)%ngllz

do i = 0,5
    nf = Tdomain%specel(n)%Near_Faces(i)
    ngll1 = Tdomain%sFace(nf)%ngll1
    ngll2 = Tdomain%sFace(nf)%ngll2
  
    if ( Tdomain%specel(n)%Orient_Faces(i) == 0 ) then

        select case (i)
         case (0)
            Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                Tdomain%specel(n)%Forces(1:ngllx-2,1:nglly-2,0,0:2)
            if (Tdomain%sFace(nf)%PML) then
             Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                  Tdomain%specel(n)%Forces1(1:ngllx-2,1:nglly-2,0,0:2)
             Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                  Tdomain%specel(n)%Forces2(1:ngllx-2,1:nglly-2,0,0:2)
             Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                  Tdomain%specel(n)%Forces3(1:ngllx-2,1:nglly-2,0,0:2)
            endif
         case (1)
            Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                Tdomain%specel(n)%Forces(1:ngllx-2,0,1:ngllz-2,0:2)
            if (Tdomain%sFace(nf)%PML) then
             Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                  Tdomain%specel(n)%Forces1(1:ngllx-2,0,1:ngllz-2,0:2)
             Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                  Tdomain%specel(n)%Forces2(1:ngllx-2,0,1:ngllz-2,0:2)
             Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                  Tdomain%specel(n)%Forces3(1:ngllx-2,0,1:ngllz-2,0:2)
            endif
         case (2)
            Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                            Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,1:ngllz-2,0:2)
            if (Tdomain%sFace(nf)%PML) then
             Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Tdomain%specel(n)%Forces1(ngllx-1,1:nglly-2,1:ngllz-2,0:2)
             Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Tdomain%specel(n)%Forces2(ngllx-1,1:nglly-2,1:ngllz-2,0:2)
             Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Tdomain%specel(n)%Forces3(ngllx-1,1:nglly-2,1:ngllz-2,0:2)
            endif
         case (3)
            Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                            Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,1:ngllz-2,0:2)
            if (Tdomain%sFace(nf)%PML) then
             Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Tdomain%specel(n)%Forces1(1:ngllx-2,nglly-1,1:ngllz-2,0:2)
             Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Tdomain%specel(n)%Forces2(1:ngllx-2,nglly-1,1:ngllz-2,0:2)
             Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Tdomain%specel(n)%Forces3(1:ngllx-2,nglly-1,1:ngllz-2,0:2)
            endif
         case (4)
            Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                Tdomain%specel(n)%Forces(0,1:nglly-2,1:ngllz-2,0:2)
            if (Tdomain%sFace(nf)%PML) then
             Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                  Tdomain%specel(n)%Forces1(0,1:nglly-2,1:ngllz-2,0:2)
             Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                  Tdomain%specel(n)%Forces2(0,1:nglly-2,1:ngllz-2,0:2)
             Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                                  Tdomain%specel(n)%Forces3(0,1:nglly-2,1:ngllz-2,0:2)
            endif
         case (5)
            Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2) + &
                                                            Tdomain%specel(n)%Forces(1:ngllx-2,1:nglly-2,ngllz-1,0:2)
            if (Tdomain%sFace(nf)%PML) then
             Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces1(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Tdomain%specel(n)%Forces1(1:ngllx-2,1:nglly-2,ngllz-1,0:2)
             Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces2(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Tdomain%specel(n)%Forces2(1:ngllx-2,1:nglly-2,ngllz-1,0:2)
             Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) = Tdomain%sFace(nf)%Forces3(1:ngll1-2,1:ngll2-2,0:2) + &
                                                              Tdomain%specel(n)%Forces3(1:ngllx-2,1:nglly-2,ngllz-1,0:2)
            endif
         end select

    else if ( Tdomain%specel(n)%Orient_Faces(i) == 1 ) then
 
    select case (i)
      case (0)
        do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces(j,:,0:2) = Tdomain%sFace(nf)%Forces(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-j,1:nglly-2,0,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
            Tdomain%sFace(nf)%Forces1(j,:,0:2) = Tdomain%sFace(nf)%Forces1(j,:,0:2) + &
                                                Tdomain%specel(n)%Forces1(ngllx-1-j,1:nglly-2,0,0:2)
            Tdomain%sFace(nf)%Forces2(j,:,0:2) = Tdomain%sFace(nf)%Forces2(j,:,0:2) + &
                                                Tdomain%specel(n)%Forces2(ngllx-1-j,1:nglly-2,0,0:2)
            Tdomain%sFace(nf)%Forces3(j,:,0:2) = Tdomain%sFace(nf)%Forces3(j,:,0:2) + &
                                                Tdomain%specel(n)%Forces3(ngllx-1-j,1:nglly-2,0,0:2)
          enddo
        endif
      case (1)
        do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces(j,:,0:2) = Tdomain%sFace(nf)%Forces(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-j,0,1:ngllz-2,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces1(j,:,0:2) = Tdomain%sFace(nf)%Forces1(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-j,0,1:ngllz-2,0:2)
         Tdomain%sFace(nf)%Forces2(j,:,0:2) = Tdomain%sFace(nf)%Forces2(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-j,0,1:ngllz-2,0:2)
         Tdomain%sFace(nf)%Forces3(j,:,0:2) = Tdomain%sFace(nf)%Forces3(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-j,0,1:ngllz-2,0:2)
          enddo
        endif
      case (2)
        do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces(j,:,0:2) = Tdomain%sFace(nf)%Forces(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,1:ngllz-2,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces1(j,:,0:2) = Tdomain%sFace(nf)%Forces1(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1,nglly-1-j,1:ngllz-2,0:2)
         Tdomain%sFace(nf)%Forces2(j,:,0:2) = Tdomain%sFace(nf)%Forces2(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1,nglly-1-j,1:ngllz-2,0:2)
         Tdomain%sFace(nf)%Forces3(j,:,0:2) = Tdomain%sFace(nf)%Forces3(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1,nglly-1-j,1:ngllz-2,0:2)
          enddo
        endif
      case (3)
        do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces(j,:,0:2) = Tdomain%sFace(nf)%Forces(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,1:ngllz-2,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces1(j,:,0:2) = Tdomain%sFace(nf)%Forces1(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-j,nglly-1,1:ngllz-2,0:2)
         Tdomain%sFace(nf)%Forces2(j,:,0:2) = Tdomain%sFace(nf)%Forces2(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-j,nglly-1,1:ngllz-2,0:2)
         Tdomain%sFace(nf)%Forces3(j,:,0:2) = Tdomain%sFace(nf)%Forces3(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-j,nglly-1,1:ngllz-2,0:2)
          enddo
        endif
      case (4)
        do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces(j,:,0:2) = Tdomain%sFace(nf)%Forces(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces(0,nglly-1-j,1:ngllz-2,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces1(j,:,0:2) = Tdomain%sFace(nf)%Forces1(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces1(0,nglly-1-j,1:ngllz-2,0:2)
         Tdomain%sFace(nf)%Forces2(j,:,0:2) = Tdomain%sFace(nf)%Forces2(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces2(0,nglly-1-j,1:ngllz-2,0:2)
         Tdomain%sFace(nf)%Forces3(j,:,0:2) = Tdomain%sFace(nf)%Forces3(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces3(0,nglly-1-j,1:ngllz-2,0:2)
          enddo
        endif
      case (5)
        do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces(j,:,0:2) = Tdomain%sFace(nf)%Forces(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-j,1:nglly-2,ngllz-1,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
         Tdomain%sFace(nf)%Forces1(j,:,0:2) = Tdomain%sFace(nf)%Forces1(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-j,1:nglly-2,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces2(j,:,0:2) = Tdomain%sFace(nf)%Forces2(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-j,1:nglly-2,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces3(j,:,0:2) = Tdomain%sFace(nf)%Forces3(j,:,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-j,1:nglly-2,ngllz-1,0:2)
          enddo
        endif
    end select

    else if ( Tdomain%specel(n)%Orient_Faces(i) == 2 ) then
 
    select case (i)
      case (0)
        do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces(:,j,0:2) = Tdomain%sFace(nf)%Forces(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1-j,0,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces1(:,j,0:2) = Tdomain%sFace(nf)%Forces1(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces1(1:ngllx-2,nglly-1-j,0,0:2)
         Tdomain%sFace(nf)%Forces2(:,j,0:2) = Tdomain%sFace(nf)%Forces2(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces2(1:ngllx-2,nglly-1-j,0,0:2)
         Tdomain%sFace(nf)%Forces3(:,j,0:2) = Tdomain%sFace(nf)%Forces3(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces3(1:ngllx-2,nglly-1-j,0,0:2)
          enddo
        endif
      case (1)
        do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces(:,j,0:2) = Tdomain%sFace(nf)%Forces(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces(1:ngllx-2,0,ngllz-1-j,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces1(:,j,0:2) = Tdomain%sFace(nf)%Forces1(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces1(1:ngllx-2,0,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(:,j,0:2) = Tdomain%sFace(nf)%Forces2(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces2(1:ngllx-2,0,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(:,j,0:2) = Tdomain%sFace(nf)%Forces3(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces3(1:ngllx-2,0,ngllz-1-j,0:2)
          enddo
        endif
      case (2)
        do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces(:,j,0:2) = Tdomain%sFace(nf)%Forces(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1,1:nglly-2,ngllz-1-j,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces1(:,j,0:2) = Tdomain%sFace(nf)%Forces1(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1,1:nglly-2,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(:,j,0:2) = Tdomain%sFace(nf)%Forces2(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1,1:nglly-2,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(:,j,0:2) = Tdomain%sFace(nf)%Forces3(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1,1:nglly-2,ngllz-1-j,0:2)
          enddo
        endif
      case (3)
        do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces(:,j,0:2) = Tdomain%sFace(nf)%Forces(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1,ngllz-1-j,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces1(:,j,0:2) = Tdomain%sFace(nf)%Forces1(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces1(1:ngllx-2,nglly-1,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(:,j,0:2) = Tdomain%sFace(nf)%Forces2(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces2(1:ngllx-2,nglly-1,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(:,j,0:2) = Tdomain%sFace(nf)%Forces3(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces3(1:ngllx-2,nglly-1,ngllz-1-j,0:2)
          enddo
        endif
      case (4)
        do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces(:,j,0:2) = Tdomain%sFace(nf)%Forces(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces(0,1:nglly-2,ngllz-1-j,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces1(:,j,0:2) = Tdomain%sFace(nf)%Forces1(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces1(0,1:nglly-2,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(:,j,0:2) = Tdomain%sFace(nf)%Forces2(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces2(0,1:nglly-2,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(:,j,0:2) = Tdomain%sFace(nf)%Forces3(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces3(0,1:nglly-2,ngllz-1-j,0:2)
          enddo
        endif
      case (5)
        do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces(:,j,0:2) = Tdomain%sFace(nf)%Forces(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces(1:ngllx-2,nglly-1-j,ngllz-1,0:2)
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll2-2 
         Tdomain%sFace(nf)%Forces1(:,j,0:2) = Tdomain%sFace(nf)%Forces1(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces1(1:ngllx-2,nglly-1-j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces2(:,j,0:2) = Tdomain%sFace(nf)%Forces2(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces2(1:ngllx-2,nglly-1-j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces3(:,j,0:2) = Tdomain%sFace(nf)%Forces3(:,j,0:2) + &
                                          Tdomain%specel(n)%Forces3(1:ngllx-2,nglly-1-j,ngllz-1,0:2)
          enddo
        endif
    end select

    else if ( Tdomain%specel(n)%Orient_Faces(i) == 3 ) then
 
    select case (i)
      case (0)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-k,0,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-j,nglly-1-k,0,0:2)
           Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-j,nglly-1-k,0,0:2)
           Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-j,nglly-1-k,0,0:2)
           enddo
          enddo
        endif
      case (1)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-j,0,ngllz-1-k,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-j,0,ngllz-1-k,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-j,0,ngllz-1-k,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-j,0,ngllz-1-k,0:2)
           enddo
          enddo
        endif
      case (2)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1,nglly-1-j,ngllz-1-k,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1,nglly-1-j,ngllz-1-k,0:2) 
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1,nglly-1-j,ngllz-1-k,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1,nglly-1-j,ngllz-1-k,0:2)
           enddo
          enddo
        endif
      case (3)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1,ngllz-1-k,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-j,nglly-1,ngllz-1-k,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-j,nglly-1,ngllz-1-k,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-j,nglly-1,ngllz-1-k,0:2)
           enddo
          enddo
        endif
      case (4)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(0,nglly-1-j,ngllz-1-k,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(0,nglly-1-j,ngllz-1-k,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(0,nglly-1-j,ngllz-1-k,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(0,nglly-1-j,ngllz-1-k,0:2)
           enddo
          enddo
        endif
      case (5)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-j,nglly-1-k,ngllz-1,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-j,nglly-1-k,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-j,nglly-1-k,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-j,nglly-1-k,ngllz-1,0:2)
           enddo
          enddo
        endif
    end select

    else if ( Tdomain%specel(n)%Orient_Faces(i) == 4 ) then
 
    select case (i)
      case (0)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(k,j,0,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(k,j,0,0:2)
           Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(k,j,0,0:2)
           Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(k,j,0,0:2)
           enddo
          enddo
        endif
      case (1)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(k,0,j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(k,0,j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(k,0,j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(k,0,j,0:2)
           enddo
          enddo
        endif
      case (2)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1,k,j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1,k,j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1,k,j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1,k,j,0:2)
           enddo
          enddo
        endif
      case (3)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(k,nglly-1,j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(k,nglly-1,j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(k,nglly-1,j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(k,nglly-1,j,0:2)
           enddo
          enddo
        endif
      case (4)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(0,k,j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(0,k,j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(0,k,j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(0,k,j,0:2)
           enddo
          enddo
        endif
      case (5)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(k,j,ngllz-1,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(k,j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(k,j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(k,j,ngllz-1,0:2)
           enddo
          enddo
        endif
    end select

    else if ( Tdomain%specel(n)%Orient_Faces(i) == 5 ) then
 
    select case (i)
      case (0)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-k,j,0,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-k,j,0,0:2)
           Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-k,j,0,0:2)
           Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-k,j,0,0:2)
           enddo
          enddo
        endif
      case (1)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-k,0,j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-k,0,j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-k,0,j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-k,0,j,0:2)
           enddo
          enddo
        endif
      case (2)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1,nglly-1-k,j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1,nglly-1-k,j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1,nglly-1-k,j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1,nglly-1-k,j,0:2)
           enddo
          enddo
        endif
      case (3)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-k,nglly-1,j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-k,nglly-1,j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-k,nglly-1,j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-k,nglly-1,j,0:2)
           enddo
          enddo
        endif
      case (4)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(0,nglly-1-k,j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(0,nglly-1-k,j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(0,nglly-1-k,j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(0,nglly-1-k,j,0:2)
           enddo
          enddo
        endif
      case (5)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-k,j,ngllz-1,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-k,j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-k,j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-k,j,ngllz-1,0:2)
           enddo
          enddo
        endif
    end select

    else if ( Tdomain%specel(n)%Orient_Faces(i) == 6 ) then
 
    select case (i)
      case (0)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(k,nglly-1-j,0,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(k,nglly-1-j,0,0:2)
           Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(k,nglly-1-j,0,0:2)
           Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(k,nglly-1-j,0,0:2)
           enddo
          enddo
        endif
      case (1)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(k,0,ngllz-1-j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(k,0,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(k,0,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(k,0,ngllz-1-j,0:2)
           enddo
          enddo
        endif
      case (2)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1,k,ngllz-1-j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1,k,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1,k,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1,k,ngllz-1-j,0:2)
           enddo
          enddo
        endif
      case (3)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(k,nglly-1,ngllz-1-j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(k,nglly-1,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(k,nglly-1,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(k,nglly-1,ngllz-1-j,0:2)
           enddo
          enddo
        endif
      case (4)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(0,k,ngllz-1-j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(0,k,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(0,k,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(0,k,ngllz-1-j,0:2)
           enddo
          enddo
        endif
      case (5)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(k,nglly-1-j,ngllz-1,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(k,nglly-1-j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(k,nglly-1-j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(k,nglly-1-j,ngllz-1,0:2)
           enddo
          enddo
        endif
    end select

    else if ( Tdomain%specel(n)%Orient_Faces(i) == 7 ) then
 
    select case (i)
      case (0)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-k,nglly-1-j,0,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
           Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-k,nglly-1-j,0,0:2)
           Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-k,nglly-1-j,0,0:2)
           Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-k,nglly-1-j,0,0:2)
           enddo
          enddo
        endif
      case (1)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-k,0,ngllz-1-j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-k,0,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-k,0,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-k,0,ngllz-1-j,0:2)
           enddo
          enddo
        endif
      case (2)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1,nglly-1-k,ngllz-1-j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1,nglly-1-k,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1,nglly-1-k,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1,nglly-1-k,ngllz-1-j,0:2)
           enddo
          enddo
        endif
      case (3)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-k,nglly-1,ngllz-1-j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-k,nglly-1,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-k,nglly-1,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-k,nglly-1,ngllz-1-j,0:2)
           enddo
          enddo
        endif
      case (4)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(0,nglly-1-k,ngllz-1-j,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(0,nglly-1-k,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(0,nglly-1-k,ngllz-1-j,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(0,nglly-1-k,ngllz-1-j,0:2)
           enddo
          enddo
        endif
      case (5)
        do j = 1,ngll1-2 
         do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces(ngllx-1-k,nglly-1-j,ngllz-1,0:2)
         enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
          do j = 1,ngll1-2 
           do k = 1,ngll2-2
         Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces1(ngllx-1-k,nglly-1-j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces2(ngllx-1-k,nglly-1-j,ngllz-1,0:2)
         Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                          Tdomain%specel(n)%Forces3(ngllx-1-k,nglly-1-j,ngllz-1,0:2)
           enddo
          enddo
        endif
    end select

  endif
enddo

return
end subroutine get_Forces_Elem2Face
