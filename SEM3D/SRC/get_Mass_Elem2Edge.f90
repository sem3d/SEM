!>
!! \file get_Mass_Elem2Edge.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine get_Mass_Elem2Edge (Tdomain,n)

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
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(1:ngllx-2,0,0)
            case (1)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(ngllx-1,1:ngllz-2,0)
            case (2)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(1:ngllx-2,nglly-1,0)
            case (3)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(0,1:nglly-2,0)
            case (4)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(ngllx-1,0,1:ngllz-2)
            case (5)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(1:ngllx-2,0,ngllz-1)
            case (6)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(0,0,1:ngllz-2)
            case (7)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,1:ngllz-2)
            case (8)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(ngllx-1,1:nglly-2,ngllz-1)
            case (9)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(1:ngllx-2,nglly-1,ngllz-1)
            case (10)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(0,nglly-1,1:ngllz-2)
            case (11)
                Tdomain%sEdge(ne)%MassMat(:) = Tdomain%sEdge(ne)%MassMat(:) + &
                    Tdomain%specel(n)%MassMat(0,1:nglly-2,ngllz-1)
            end select

            if (Tdomain%sEdge(ne)%PML) then
                select case (i)
                case (0)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(1:ngllx-2,0,0,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(1:ngllx-2,0,0)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(1:ngllx-2,0,0)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(1:ngllx-2,0,0)
                    endif
                case (1)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(ngllx-1,1:nglly-2,0,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(ngllx-1,1:nglly-2,0)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(ngllx-1,1:nglly-2,0)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(ngllx-1,1:nglly-2,0)
                    endif
                case (2)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(1:ngllx-2,nglly-1,0,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(1:ngllx-2,nglly-1,0)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(1:ngllx-2,nglly-1,0)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(1:ngllx-2,nglly-1,0)
                    endif
                case (3)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(0,1:nglly-2,0,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(0,1:nglly-2,0)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(0,1:nglly-2,0)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(0,1:nglly-2,0)
                    endif
                case (4)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(ngllx-1,0,1:ngllz-2,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(ngllx-1,0,1:ngllz-2)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(ngllx-1,0,1:ngllz-2)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(ngllx-1,0,1:ngllz-2)
                    endif
                case (5)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(1:ngllx-2,0,ngllz-1,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(1:ngllx-2,0,ngllz-1)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(1:ngllx-2,0,ngllz-1)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(1:ngllx-2,0,ngllz-1)
                    endif
                case (6)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(0,0,1:ngllz-2,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(0,0,1:ngllz-2)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(0,0,1:ngllz-2)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(0,0,1:ngllz-2)
                    endif
                case (7)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,1:ngllz-2,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(ngllx-1,nglly-1,1:ngllz-2)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(ngllx-1,nglly-1,1:ngllz-2)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(ngllx-1,nglly-1,1:ngllz-2)
                    endif
                case (8)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(ngllx-1,1:nglly-2,ngllz-1,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(ngllx-1,1:nglly-2,ngllz-1)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(ngllx-1,1:nglly-2,ngllz-1)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(ngllx-1,1:nglly-2,ngllz-1)
                    endif
                case (9)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(1:ngllx-2,nglly-1,ngllz-1,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(1:ngllx-2,nglly-1,ngllz-1)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(1:ngllx-2,nglly-1,ngllz-1)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(1:ngllx-2,nglly-1,ngllz-1)
                    endif
                case (10)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(0,nglly-1,1:ngllz-2,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(0,nglly-1,1:ngllz-2)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(0,nglly-1,1:ngllz-2)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(0,nglly-1,1:ngllz-2)
                    endif
                case (11)
                    Tdomain%sEdge(ne)%DumpMass(:,:) = Tdomain%sEdge(ne)%DumpMass(:,:) + &
                        Tdomain%specel(n)%DumpMass(0,1:nglly-2,ngllz-1,:)
                    if (Tdomain%sEdge(ne)%FPML) then
                        Tdomain%sEdge(ne)%Ivx(:) = Tdomain%sEdge(ne)%Ivx(:) + &
                            Tdomain%specel(n)%Ivx(0,1:nglly-2,ngllz-1)
                        Tdomain%sEdge(ne)%Ivy(:) = Tdomain%sEdge(ne)%Ivy(:) + &
                            Tdomain%specel(n)%Ivy(0,1:nglly-2,ngllz-1)
                        Tdomain%sEdge(ne)%Ivz(:) = Tdomain%sEdge(ne)%Ivz(:) + &
                            Tdomain%specel(n)%Ivz(0,1:nglly-2,ngllz-1)
                    endif
                end select
            endif

        else

            select case (i)
            case (0)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(ngllx-1-j,0,0)
                enddo
            case (1)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(ngllx-1,nglly-1-j,0)
                enddo
            case (2)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(ngllx-1-j,nglly-1,0)
                enddo
            case (3)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(0,nglly-1-j,0)
                enddo
            case (4)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(ngllx-1,0,ngllz-1-j)
                enddo
            case (5)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(ngllx-1-j,0,ngllz-1)
                enddo
            case (6)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(0,0,ngllz-1-j)
                enddo
            case (7)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(ngllx-1,nglly-1,ngllz-1-j)
                enddo
            case (8)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(ngllx-1,nglly-1-j,ngllz-1)
                enddo
            case (9)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(ngllx-1-j,nglly-1,ngllz-1)
                enddo
            case (10)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(0,nglly-1,ngllz-1-j)
                enddo
            case (11)
                do j = 1,ngll-2
                    Tdomain%sEdge(ne)%MassMat(j) = Tdomain%sEdge(ne)%MassMat(j) + &
                        Tdomain%specel(n)%MassMat(0,nglly-1-j,ngllz-1)
                enddo
            end select

            if (Tdomain%sEdge(ne)%PML) then
                select case (i)
                case (0)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1-j,0,0,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(ngllx-1-j,0,0)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(ngllx-1-j,0,0)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(ngllx-1-j,0,0)
                        endif
                    enddo
                case (1)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1-j,0,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(ngllx-1,nglly-1-j,0)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(ngllx-1,nglly-1-j,0)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(ngllx-1,nglly-1-j,0)
                        endif
                    enddo
                case (2)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1-j,nglly-1,0,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(ngllx-1-j,nglly-1,0)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(ngllx-1-j,nglly-1,0)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(ngllx-1-j,nglly-1,0)
                        endif
                    enddo
                case (3)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(0,nglly-1-j,0,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(0,nglly-1-j,0)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(0,nglly-1-j,0)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(0,nglly-1-j,0)
                        endif
                    enddo
                case (4)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1,0,ngllz-1-j,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(ngllx-1,0,ngllz-1-j)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(ngllx-1,0,ngllz-1-j)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(ngllx-1,0,ngllz-1-j)
                        endif
                    enddo
                case (5)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1-j,0,ngllz-1,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(ngllx-1-j,0,ngllz-1)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(ngllx-1-j,0,ngllz-1)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(ngllx-1-j,0,ngllz-1)
                        endif
                    enddo
                case (6)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(0,0,ngllz-1-j,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(0,0,ngllz-1-j)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(0,0,ngllz-1-j)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(0,0,ngllz-1-j)
                        endif
                    enddo
                case (7)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1,ngllz-1-j,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(ngllx-1,nglly-1,ngllz-1-j)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(ngllx-1,nglly-1,ngllz-1-j)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(ngllx-1,nglly-1,ngllz-1-j)
                        endif
                    enddo
                case (8)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1-j,ngllz-1,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(ngllx-1,nglly-1-j,ngllz-1)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(ngllx-1,nglly-1-j,ngllz-1)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(ngllx-1,nglly-1-j,ngllz-1)
                        endif
                    enddo
                case (9)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1-j,nglly-1,ngllz-1,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(ngllx-1-j,nglly-1,ngllz-1)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(ngllx-1-j,nglly-1,ngllz-1)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(ngllx-1-j,nglly-1,ngllz-1)
                        endif
                    enddo
                case (10)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(0,nglly-1,ngllz-1-j,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(0,nglly-1,ngllz-1-j)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(0,nglly-1,ngllz-1-j)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(0,nglly-1,ngllz-1-j)
                        endif
                    enddo
                case (11)
                    do j = 1,ngll-2
                        Tdomain%sEdge(ne)%DumpMass(j,:) = Tdomain%sEdge(ne)%DumpMass(j,:) + &
                            Tdomain%specel(n)%DumpMass(0,nglly-1-j,ngllz-1,:)
                        if (Tdomain%sEdge(ne)%FPML) then
                            Tdomain%sEdge(ne)%Ivx(j) = Tdomain%sEdge(ne)%Ivx(j) + &
                                Tdomain%specel(n)%Ivx(0,nglly-1-j,ngllz-1)
                            Tdomain%sEdge(ne)%Ivy(j) = Tdomain%sEdge(ne)%Ivy(j) + &
                                Tdomain%specel(n)%Ivy(0,nglly-1-j,ngllz-1)
                            Tdomain%sEdge(ne)%Ivz(j) = Tdomain%sEdge(ne)%Ivz(j) + &
                                Tdomain%specel(n)%Ivz(0,nglly-1-j,ngllz-1)
                        endif
                    enddo
                end select
            endif

        endif

    enddo


    return
end subroutine get_Mass_Elem2Edge
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
