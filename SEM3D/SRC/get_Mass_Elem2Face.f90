!>
!! \file get_Mass_Elem2Face.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine get_Mass_Elem2Face (Tdomain,n)

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
                Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                    Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,0)
            case (1)
                Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                    Tdomain%specel(n)%MassMat(1:ngllx-2,0,1:ngllz-2)
            case (2)
                Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                    Tdomain%specel(n)%MassMat(ngllx-1,1:nglly-2,1:ngllz-2)
            case (3)
                Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                    Tdomain%specel(n)%MassMat(1:ngllx-2,nglly-1,1:ngllz-2)
            case (4)
                Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                    Tdomain%specel(n)%MassMat(0,1:nglly-2,1:ngllz-2)
            case (5)
                Tdomain%sFace(nf)%MassMat(:,:) = Tdomain%sFace(nf)%MassMat(:,:) + &
                    Tdomain%specel(n)%MassMat(1:ngllx-2,1:nglly-2,ngllz-1)
            end select
            if (Tdomain%sFace(nf)%PML) then
                select case (i)
                case (0)
                    Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                        Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,0,:)
                    if (Tdomain%sFace(nf)%FPML) then
                        Tdomain%sFace(nf)%Ivx(:,:) = Tdomain%sFace(nf)%Ivx(:,:) + &
                            Tdomain%specel(n)%Ivx(1:ngllx-2,1:nglly-2,0)
                        Tdomain%sFace(nf)%Ivy(:,:) = Tdomain%sFace(nf)%Ivy(:,:) + &
                            Tdomain%specel(n)%Ivy(1:ngllx-2,1:nglly-2,0)
                        Tdomain%sFace(nf)%Ivz(:,:) = Tdomain%sFace(nf)%Ivz(:,:) + &
                            Tdomain%specel(n)%Ivz(1:ngllx-2,1:nglly-2,0)
                    endif
                case (1)
                    Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                        Tdomain%specel(n)%DumpMass(1:ngllx-2,0,1:ngllz-2,:)
                    if (Tdomain%sFace(nf)%FPML) then
                        Tdomain%sFace(nf)%Ivx(:,:) = Tdomain%sFace(nf)%Ivx(:,:) + &
                            Tdomain%specel(n)%Ivx(1:ngllx-2,0,1:ngllz-2)
                        Tdomain%sFace(nf)%Ivy(:,:) = Tdomain%sFace(nf)%Ivy(:,:) + &
                            Tdomain%specel(n)%Ivy(1:ngllx-2,0,1:ngllz-2)
                        Tdomain%sFace(nf)%Ivz(:,:) = Tdomain%sFace(nf)%Ivz(:,:) + &
                            Tdomain%specel(n)%Ivz(1:ngllx-2,0,1:ngllz-2)
                    endif
                case (2)
                    Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                        Tdomain%specel(n)%DumpMass(ngllx-1,1:nglly-2,1:ngllz-2,:)
                    if (Tdomain%sFace(nf)%FPML) then
                        Tdomain%sFace(nf)%Ivx(:,:) = Tdomain%sFace(nf)%Ivx(:,:) + &
                            Tdomain%specel(n)%Ivx(ngllx-1,1:nglly-2,1:ngllz-2)
                        Tdomain%sFace(nf)%Ivy(:,:) = Tdomain%sFace(nf)%Ivy(:,:) + &
                            Tdomain%specel(n)%Ivy(ngllx-1,1:nglly-2,1:ngllz-2)
                        Tdomain%sFace(nf)%Ivz(:,:) = Tdomain%sFace(nf)%Ivz(:,:) + &
                            Tdomain%specel(n)%Ivz(ngllx-1,1:nglly-2,1:ngllz-2)
                    endif
                case (3)
                    Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                        Tdomain%specel(n)%DumpMass(1:ngllx-2,nglly-1,1:ngllz-2,:)
                    if (Tdomain%sFace(nf)%FPML) then
                        Tdomain%sFace(nf)%Ivx(:,:) = Tdomain%sFace(nf)%Ivx(:,:) + &
                            Tdomain%specel(n)%Ivx(1:ngllx-2,nglly-1,1:ngllz-2)
                        Tdomain%sFace(nf)%Ivy(:,:) = Tdomain%sFace(nf)%Ivy(:,:) + &
                            Tdomain%specel(n)%Ivy(1:ngllx-2,nglly-1,1:ngllz-2)
                        Tdomain%sFace(nf)%Ivz(:,:) = Tdomain%sFace(nf)%Ivz(:,:) + &
                            Tdomain%specel(n)%Ivz(1:ngllx-2,nglly-1,1:ngllz-2)
                    endif
                case (4)
                    Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                        Tdomain%specel(n)%DumpMass(0,1:nglly-2,1:ngllz-2,:)
                    if (Tdomain%sFace(nf)%FPML) then
                        Tdomain%sFace(nf)%Ivx(:,:) = Tdomain%sFace(nf)%Ivx(:,:) + &
                            Tdomain%specel(n)%Ivx(0,1:nglly-2,1:ngllz-2)
                        Tdomain%sFace(nf)%Ivy(:,:) = Tdomain%sFace(nf)%Ivy(:,:) + &
                            Tdomain%specel(n)%Ivy(0,1:nglly-2,1:ngllz-2)
                        Tdomain%sFace(nf)%Ivz(:,:) = Tdomain%sFace(nf)%Ivz(:,:) + &
                            Tdomain%specel(n)%Ivz(0,1:nglly-2,1:ngllz-2)
                    endif
                case (5)
                    Tdomain%sFace(nf)%DumpMass(:,:,:) = Tdomain%sFace(nf)%DumpMass(:,:,:) + &
                        Tdomain%specel(n)%DumpMass(1:ngllx-2,1:nglly-2,ngllz-1,:)
                    if (Tdomain%sFace(nf)%FPML) then
                        Tdomain%sFace(nf)%Ivx(:,:) = Tdomain%sFace(nf)%Ivx(:,:) + &
                            Tdomain%specel(n)%Ivx(1:ngllx-2,1:nglly-2,ngllz-1)
                        Tdomain%sFace(nf)%Ivy(:,:) = Tdomain%sFace(nf)%Ivy(:,:) + &
                            Tdomain%specel(n)%Ivy(1:ngllx-2,1:nglly-2,ngllz-1)
                        Tdomain%sFace(nf)%Ivz(:,:) = Tdomain%sFace(nf)%Ivz(:,:) + &
                            Tdomain%specel(n)%Ivz(1:ngllx-2,1:nglly-2,ngllz-1)
                    endif
                end select
            endif

        else if ( Tdomain%specel(n)%Orient_Faces(i) == 1 ) then

            select case (i)
            case (0)
                do j = 1,ngll1-2
                    Tdomain%sFace(nf)%MassMat(j,:) = Tdomain%sFace(nf)%MassMat(j,:) + &
                        Tdomain%specel(n)%MassMat(ngllx-1-j,1:nglly-2,0)
                enddo
            case (1)
                do j = 1,ngll1-2
                    Tdomain%sFace(nf)%MassMat(j,:) = Tdomain%sFace(nf)%MassMat(j,:) + &
                        Tdomain%specel(n)%MassMat(ngllx-1-j,0,1:ngllz-2)
                enddo
            case (2)
                do j = 1,ngll1-2
                    Tdomain%sFace(nf)%MassMat(j,:) = Tdomain%sFace(nf)%MassMat(j,:) + &
                        Tdomain%specel(n)%MassMat(ngllx-1,nglly-1-j,1:ngllz-2)
                enddo
            case (3)
                do j = 1,ngll1-2
                    Tdomain%sFace(nf)%MassMat(j,:) = Tdomain%sFace(nf)%MassMat(j,:) + &
                        Tdomain%specel(n)%MassMat(ngllx-1-j,nglly-1,1:ngllz-2)
                enddo
            case (4)
                do j = 1,ngll1-2
                    Tdomain%sFace(nf)%MassMat(j,:) = Tdomain%sFace(nf)%MassMat(j,:) + &
                        Tdomain%specel(n)%MassMat(0,nglly-1-j,1:ngllz-2)
                enddo
            case (5)
                do j = 1,ngll1-2
                    Tdomain%sFace(nf)%MassMat(j,:) = Tdomain%sFace(nf)%MassMat(j,:) + &
                        Tdomain%specel(n)%MassMat(ngllx-1-j,1:nglly-2,ngllz-1)
                enddo
            end select

            if (Tdomain%sFace(nf)%PML) then
                select case (i)
                case (0)
                    do j = 1,ngll1-2
                        Tdomain%sFace(nf)%DumpMass(j,:,:) = Tdomain%sFace(nf)%DumpMass(j,:,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1-j,1:nglly-2,0,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(j,:) = Tdomain%sFace(nf)%Ivx(j,:) + &
                                Tdomain%specel(n)%Ivx(ngllx-1-j,1:nglly-2,0)
                            Tdomain%sFace(nf)%Ivy(j,:) = Tdomain%sFace(nf)%Ivy(j,:) + &
                                Tdomain%specel(n)%Ivy(ngllx-1-j,1:nglly-2,0)
                            Tdomain%sFace(nf)%Ivz(j,:) = Tdomain%sFace(nf)%Ivz(j,:) + &
                                Tdomain%specel(n)%Ivz(ngllx-1-j,1:nglly-2,0)
                        endif
                    enddo
                case (1)
                    do j = 1,ngll1-2
                        Tdomain%sFace(nf)%DumpMass(j,:,:) = Tdomain%sFace(nf)%DumpMass(j,:,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1-j,0,1:ngllz-2,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(j,:) = Tdomain%sFace(nf)%Ivx(j,:) + &
                                Tdomain%specel(n)%Ivx(ngllx-1-j,0,1:ngllz-2)
                            Tdomain%sFace(nf)%Ivy(j,:) = Tdomain%sFace(nf)%Ivy(j,:) + &
                                Tdomain%specel(n)%Ivy(ngllx-1-j,0,1:ngllz-2)
                            Tdomain%sFace(nf)%Ivz(j,:) = Tdomain%sFace(nf)%Ivz(j,:) + &
                                Tdomain%specel(n)%Ivz(ngllx-1-j,0,1:ngllz-2)
                        endif
                    enddo
                case (2)
                    do j = 1,ngll1-2
                        Tdomain%sFace(nf)%DumpMass(j,:,:) = Tdomain%sFace(nf)%DumpMass(j,:,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1-j,1:ngllz-2,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(j,:) = Tdomain%sFace(nf)%Ivx(j,:) + &
                                Tdomain%specel(n)%Ivx(ngllx-1,nglly-1-j,1:ngllz-2)
                            Tdomain%sFace(nf)%Ivy(j,:) = Tdomain%sFace(nf)%Ivy(j,:) + &
                                Tdomain%specel(n)%Ivy(ngllx-1,nglly-1-j,1:ngllz-2)
                            Tdomain%sFace(nf)%Ivz(j,:) = Tdomain%sFace(nf)%Ivz(j,:) + &
                                Tdomain%specel(n)%Ivz(ngllx-1,nglly-1-j,1:ngllz-2)
                        endif
                    enddo
                case (3)
                    do j = 1,ngll1-2
                        Tdomain%sFace(nf)%DumpMass(j,:,:) = Tdomain%sFace(nf)%DumpMass(j,:,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1-j,nglly-1,1:ngllz-2,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(j,:) = Tdomain%sFace(nf)%Ivx(j,:) + &
                                Tdomain%specel(n)%Ivx(ngllx-1-j,nglly-1,1:ngllz-2)
                            Tdomain%sFace(nf)%Ivy(j,:) = Tdomain%sFace(nf)%Ivy(j,:) + &
                                Tdomain%specel(n)%Ivy(ngllx-1-j,nglly-1,1:ngllz-2)
                            Tdomain%sFace(nf)%Ivz(j,:) = Tdomain%sFace(nf)%Ivz(j,:) + &
                                Tdomain%specel(n)%Ivz(ngllx-1-j,nglly-1,1:ngllz-2)
                        endif
                    enddo
                case (4)
                    do j = 1,ngll1-2
                        Tdomain%sFace(nf)%DumpMass(j,:,:) = Tdomain%sFace(nf)%DumpMass(j,:,:) + &
                            Tdomain%specel(n)%DumpMass(0,nglly-1-j,1:ngllz-2,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(j,:) = Tdomain%sFace(nf)%Ivx(j,:) + &
                                Tdomain%specel(n)%Ivx(0,nglly-1-j,1:ngllz-2)
                            Tdomain%sFace(nf)%Ivy(j,:) = Tdomain%sFace(nf)%Ivy(j,:) + &
                                Tdomain%specel(n)%Ivy(0,nglly-1-j,1:ngllz-2)
                            Tdomain%sFace(nf)%Ivz(j,:) = Tdomain%sFace(nf)%Ivz(j,:) + &
                                Tdomain%specel(n)%Ivz(0,nglly-1-j,1:ngllz-2)
                        endif
                    enddo
                case (5)
                    do j = 1,ngll1-2
                        Tdomain%sFace(nf)%DumpMass(j,:,:) = Tdomain%sFace(nf)%DumpMass(j,:,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1-j,1:nglly-2,ngllz-1,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(j,:) = Tdomain%sFace(nf)%Ivx(j,:) + &
                                Tdomain%specel(n)%Ivx(ngllx-1-j,1:nglly-2,ngllz-1)
                            Tdomain%sFace(nf)%Ivy(j,:) = Tdomain%sFace(nf)%Ivy(j,:) + &
                                Tdomain%specel(n)%Ivy(ngllx-1-j,1:nglly-2,ngllz-1)
                            Tdomain%sFace(nf)%Ivz(j,:) = Tdomain%sFace(nf)%Ivz(j,:) + &
                                Tdomain%specel(n)%Ivz(ngllx-1-j,1:nglly-2,ngllz-1)
                        endif
                    enddo
                end select
            endif

        else if ( Tdomain%specel(n)%Orient_Faces(i) == 2 ) then

            select case (i)
            case (0)
                do j = 1,ngll2-2
                    Tdomain%sFace(nf)%MassMat(:,j) = Tdomain%sFace(nf)%MassMat(:,j) + &
                        Tdomain%specel(n)%MassMat(1:ngllx-2,nglly-1-j,0)
                enddo
            case (1)
                do j = 1,ngll2-2
                    Tdomain%sFace(nf)%MassMat(:,j) = Tdomain%sFace(nf)%MassMat(:,j) + &
                        Tdomain%specel(n)%MassMat(1:ngllx-2,0,ngllz-1-j)
                enddo
            case (2)
                do j = 1,ngll2-2
                    Tdomain%sFace(nf)%MassMat(:,j) = Tdomain%sFace(nf)%MassMat(:,j) + &
                        Tdomain%specel(n)%MassMat(ngllx-1,1:nglly-2,ngllz-1-j)
                enddo
            case (3)
                do j = 1,ngll2-2
                    Tdomain%sFace(nf)%MassMat(:,j) = Tdomain%sFace(nf)%MassMat(:,j) + &
                        Tdomain%specel(n)%MassMat(1:ngllx-2,nglly-1,ngllz-1-j)
                enddo
            case (4)
                do j = 1,ngll2-2
                    Tdomain%sFace(nf)%MassMat(:,j) = Tdomain%sFace(nf)%MassMat(:,j) + &
                        Tdomain%specel(n)%MassMat(0,1:nglly-2,ngllz-1-j)
                enddo
            case (5)
                do j = 1,ngll2-2
                    Tdomain%sFace(nf)%MassMat(:,j) = Tdomain%sFace(nf)%MassMat(:,j) + &
                        Tdomain%specel(n)%MassMat(1:ngllx-2,nglly-1-j,ngllz-1)
                enddo
            end select

            if (Tdomain%sFace(nf)%PML) then
                select case (i)
                case (0)
                    do j = 1,ngll2-2
                        Tdomain%sFace(nf)%DumpMass(:,j,:) = Tdomain%sFace(nf)%DumpMass(:,j,:) + &
                            Tdomain%specel(n)%DumpMass(1:ngllx-2,nglly-1-j,0,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(:,j) = Tdomain%sFace(nf)%Ivx(:,j) + &
                                Tdomain%specel(n)%Ivx(1:ngllx-2,nglly-1-j,0)
                            Tdomain%sFace(nf)%Ivy(:,j) = Tdomain%sFace(nf)%Ivy(:,j) + &
                                Tdomain%specel(n)%Ivy(1:ngllx-2,nglly-1-j,0)
                            Tdomain%sFace(nf)%Ivz(:,j) = Tdomain%sFace(nf)%Ivz(:,j) + &
                                Tdomain%specel(n)%Ivz(1:ngllx-2,nglly-1-j,0)
                        endif
                    enddo
                case (1)
                    do j = 1,ngll2-2
                        Tdomain%sFace(nf)%DumpMass(:,j,:) = Tdomain%sFace(nf)%DumpMass(:,j,:) + &
                            Tdomain%specel(n)%DumpMass(1:ngllx-2,0,ngllz-1-j,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(:,j) = Tdomain%sFace(nf)%Ivx(:,j) + &
                                Tdomain%specel(n)%Ivx(1:ngllx-2,0,ngllz-1-j)
                            Tdomain%sFace(nf)%Ivy(:,j) = Tdomain%sFace(nf)%Ivy(:,j) + &
                                Tdomain%specel(n)%Ivy(1:ngllx-2,0,ngllz-1-j)
                            Tdomain%sFace(nf)%Ivz(:,j) = Tdomain%sFace(nf)%Ivz(:,j) + &
                                Tdomain%specel(n)%Ivz(1:ngllx-2,0,ngllz-1-j)
                        endif
                    enddo
                case (2)
                    do j = 1,ngll2-2
                        Tdomain%sFace(nf)%DumpMass(:,j,:) = Tdomain%sFace(nf)%DumpMass(:,j,:) + &
                            Tdomain%specel(n)%DumpMass(ngllx-1,1:nglly-2,ngllz-1-j,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(:,j) = Tdomain%sFace(nf)%Ivx(:,j) + &
                                Tdomain%specel(n)%Ivx(ngllx-1,1:nglly-2,ngllz-1-j)
                            Tdomain%sFace(nf)%Ivy(:,j) = Tdomain%sFace(nf)%Ivy(:,j) + &
                                Tdomain%specel(n)%Ivy(ngllx-1,1:nglly-2,ngllz-1-j)
                            Tdomain%sFace(nf)%Ivz(:,j) = Tdomain%sFace(nf)%Ivz(:,j) + &
                                Tdomain%specel(n)%Ivz(ngllx-1,1:nglly-2,ngllz-1-j)
                        endif
                    enddo
                case (3)
                    do j = 1,ngll2-2
                        Tdomain%sFace(nf)%DumpMass(:,j,:) = Tdomain%sFace(nf)%DumpMass(:,j,:) + &
                            Tdomain%specel(n)%DumpMass(1:ngllx-2,nglly-1,ngllz-1-j,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(:,j) = Tdomain%sFace(nf)%Ivx(:,j) + &
                                Tdomain%specel(n)%Ivx(1:ngllx-2,nglly-1,ngllz-1-j)
                            Tdomain%sFace(nf)%Ivy(:,j) = Tdomain%sFace(nf)%Ivy(:,j) + &
                                Tdomain%specel(n)%Ivy(1:ngllx-2,nglly-1,ngllz-1-j)
                            Tdomain%sFace(nf)%Ivz(:,j) = Tdomain%sFace(nf)%Ivz(:,j) + &
                                Tdomain%specel(n)%Ivz(1:ngllx-2,nglly-1,ngllz-1-j)
                        endif
                    enddo
                case (4)
                    do j = 1,ngll2-2
                        Tdomain%sFace(nf)%DumpMass(:,j,:) = Tdomain%sFace(nf)%DumpMass(:,j,:) + &
                            Tdomain%specel(n)%DumpMass(0,1:nglly-2,ngllz-1-j,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(:,j) = Tdomain%sFace(nf)%Ivx(:,j) + &
                                Tdomain%specel(n)%Ivx(0,1:nglly-2,ngllz-1-j)
                            Tdomain%sFace(nf)%Ivy(:,j) = Tdomain%sFace(nf)%Ivy(:,j) + &
                                Tdomain%specel(n)%Ivy(0,1:nglly-2,ngllz-1-j)
                            Tdomain%sFace(nf)%Ivz(:,j) = Tdomain%sFace(nf)%Ivz(:,j) + &
                                Tdomain%specel(n)%Ivz(0,1:nglly-2,ngllz-1-j)
                        endif
                    enddo
                case (5)
                    do j = 1,ngll2-2
                        Tdomain%sFace(nf)%DumpMass(:,j,:) = Tdomain%sFace(nf)%DumpMass(:,j,:) + &
                            Tdomain%specel(n)%DumpMass(1:ngllx-2,nglly-1-j,ngllz-1,:)
                        if (Tdomain%sFace(nf)%FPML) then
                            Tdomain%sFace(nf)%Ivx(:,j) = Tdomain%sFace(nf)%Ivx(:,j) + &
                                Tdomain%specel(n)%Ivx(1:ngllx-2,nglly-1-j,ngllz-1)
                            Tdomain%sFace(nf)%Ivy(:,j) = Tdomain%sFace(nf)%Ivy(:,j) + &
                                Tdomain%specel(n)%Ivy(1:ngllx-2,nglly-1-j,ngllz-1)
                            Tdomain%sFace(nf)%Ivz(:,j) = Tdomain%sFace(nf)%Ivz(:,j) + &
                                Tdomain%specel(n)%Ivz(1:ngllx-2,nglly-1-j,ngllz-1)
                        endif
                    enddo
                end select
            endif

        else if ( Tdomain%specel(n)%Orient_Faces(i) == 3 ) then

            select case (i)
            case (0)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-j,nglly-1-k,0)
                    enddo
                enddo
            case (1)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-j,0,ngllz-1-k)
                    enddo
                enddo
            case (2)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1,nglly-1-j,ngllz-1-k)
                    enddo
                enddo
            case (3)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-j,nglly-1,ngllz-1-k)
                    enddo
                enddo
            case (4)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(0,nglly-1-j,ngllz-1-k)
                    enddo
                enddo
            case (5)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-j,nglly-1-k,ngllz-1)
                    enddo
                enddo
            end select

            if (Tdomain%sFace(nf)%PML) then
                select case (i)
                case (0)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-j,nglly-1-k,0,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-j,nglly-1-k,0)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-j,nglly-1-k,0)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-j,nglly-1-k,0)
                            endif
                        enddo
                    enddo
                case (1)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-j,0,ngllz-1-k,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-j,0,ngllz-1-k)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-j,0,ngllz-1-k)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-j,0,ngllz-1-k)
                            endif
                        enddo
                    enddo
                case (2)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1-j,ngllz-1-k,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1,nglly-1-j,ngllz-1-k)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1,nglly-1-j,ngllz-1-k)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1,nglly-1-j,ngllz-1-k)
                            endif
                        enddo
                    enddo
                case (3)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-j,nglly-1,ngllz-1-k,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-j,nglly-1,ngllz-1-k)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-j,nglly-1,ngllz-1-k)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-j,nglly-1,ngllz-1-k)
                            endif
                        enddo
                    enddo
                case (4)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(0,nglly-1-j,ngllz-1-k,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(0,nglly-1-j,ngllz-1-k)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(0,nglly-1-j,ngllz-1-k)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(0,nglly-1-j,ngllz-1-k)
                            endif
                        enddo
                    enddo
                case (5)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-j,nglly-1-k,ngllz-1,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-j,nglly-1-k,ngllz-1)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-j,nglly-1-k,ngllz-1)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-j,nglly-1-k,ngllz-1)
                            endif
                        enddo
                    enddo
                end select
            endif

        else if ( Tdomain%specel(n)%Orient_Faces(i) == 4 ) then

            select case (i)
            case (0)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(k,j,0)
                    enddo
                enddo
            case (1)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(k,0,j)
                    enddo
                enddo
            case (2)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1,k,j)
                    enddo
                enddo
            case (3)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(k,nglly-1,j)
                    enddo
                enddo
            case (4)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(0,k,j)
                    enddo
                enddo
            case (5)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(k,j,ngllz-1)
                    enddo
                enddo
            end select

            if (Tdomain%sFace(nf)%PML) then
                select case (i)
                case (0)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(k,j,0,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(k,j,0)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(k,j,0)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(k,j,0)
                            endif
                        enddo
                    enddo
                case (1)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(k,0,j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(k,0,j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(k,0,j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(k,0,j)
                            endif
                        enddo
                    enddo
                case (2)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1,k,j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1,k,j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1,k,j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1,k,j)
                            endif
                        enddo
                    enddo
                case (3)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(k,nglly-1,j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(k,nglly-1,j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(k,nglly-1,j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(k,nglly-1,j)
                            endif
                        enddo
                    enddo
                case (4)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(0,k,j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(0,k,j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(0,k,j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(0,k,j)
                            endif
                        enddo
                    enddo
                case (5)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(k,j,ngllz-1,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(k,j,ngllz-1)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(k,j,ngllz-1)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(k,j,ngllz-1)
                            endif
                        enddo
                    enddo
                end select
            endif

        else if ( Tdomain%specel(n)%Orient_Faces(i) == 5 ) then

            select case (i)
            case (0)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-k,j,0)
                    enddo
                enddo
            case (1)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-k,0,j)
                    enddo
                enddo
            case (2)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1,nglly-1-k,j)
                    enddo
                enddo
            case (3)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-k,nglly-1,j)
                    enddo
                enddo
            case (4)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(0,nglly-1-k,j)
                    enddo
                enddo
            case (5)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-k,j,ngllz-1)
                    enddo
                enddo
            end select

            if (Tdomain%sFace(nf)%PML) then
                select case (i)
                case (0)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-k,j,0,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-k,j,0)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-k,j,0)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-k,j,0)
                            endif
                        enddo
                    enddo
                case (1)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-k,0,j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-k,0,j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-k,0,j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-k,0,j)
                            endif
                        enddo
                    enddo
                case (2)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1-k,j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1,nglly-1-k,j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1,nglly-1-k,j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1,nglly-1-k,j)
                            endif
                        enddo
                    enddo
                case (3)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-k,nglly-1,j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-k,nglly-1,j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-k,nglly-1,j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-k,nglly-1,j)
                            endif
                        enddo
                    enddo
                case (4)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(0,nglly-1-k,j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(0,nglly-1-k,j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(0,nglly-1-k,j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(0,nglly-1-k,j)
                            endif
                        enddo
                    enddo
                case (5)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-k,j,ngllz-1,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-k,j,ngllz-1)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-k,j,ngllz-1)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-k,j,ngllz-1)
                            endif
                        enddo
                    enddo
                end select
            endif

        else if ( Tdomain%specel(n)%Orient_Faces(i) == 6 ) then

            select case (i)
            case (0)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(k,nglly-1-j,0)
                    enddo
                enddo
            case (1)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(k,0,ngllz-1-j)
                    enddo
                enddo
            case (2)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1,k,ngllz-1-j)
                    enddo
                enddo
            case (3)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(k,nglly-1,ngllz-1-j)
                    enddo
                enddo
            case (4)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(0,k,ngllz-1-j)
                    enddo
                enddo
            case (5)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(k,nglly-1-j,ngllz-1)
                    enddo
                enddo
            end select

            if (Tdomain%sFace(nf)%PML) then
                select case (i)
                case (0)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(k,nglly-1-j,0,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(k,nglly-1-j,0)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(k,nglly-1-j,0)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(k,nglly-1-j,0)
                            endif
                        enddo
                    enddo
                case (1)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(k,0,ngllz-1-j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(k,0,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(k,0,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(k,0,ngllz-1-j)
                            endif
                        enddo
                    enddo
                case (2)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1,k,ngllz-1-j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1,k,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1,k,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1,k,ngllz-1-j)
                            endif
                        enddo
                    enddo
                case (3)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(k,nglly-1,ngllz-1-j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(k,nglly-1,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(k,nglly-1,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(k,nglly-1,ngllz-1-j)
                            endif
                        enddo
                    enddo
                case (4)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(0,k,ngllz-1-j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(0,k,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(0,k,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(0,k,ngllz-1-j)
                            endif
                        enddo
                    enddo
                case (5)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(k,nglly-1-j,ngllz-1,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(k,nglly-1-j,ngllz-1)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(k,nglly-1-j,ngllz-1)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(k,nglly-1-j,ngllz-1)
                            endif
                        enddo
                    enddo
                end select
            endif

        else if ( Tdomain%specel(n)%Orient_Faces(i) == 7 ) then

            select case (i)
            case (0)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-k,nglly-1-j,0)
                    enddo
                enddo
            case (1)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-k,0,ngllz-1-j)
                    enddo
                enddo
            case (2)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1,nglly-1-k,ngllz-1-j)
                    enddo
                enddo
            case (3)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-k,nglly-1,ngllz-1-j)
                    enddo
                enddo
            case (4)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(0,nglly-1-k,ngllz-1-j)
                    enddo
                enddo
            case (5)
                do j = 1,ngll1-2
                    do k = 1,ngll2-2
                        Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + &
                            Tdomain%specel(n)%MassMat(ngllx-1-k,nglly-1-j,ngllz-1)
                    enddo
                enddo
            end select

            if (Tdomain%sFace(nf)%PML) then
                select case (i)
                case (0)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-k,nglly-1-j,0,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-k,nglly-1-j,0)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-k,nglly-1-j,0)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-k,nglly-1-j,0)
                            endif
                        enddo
                    enddo
                case (1)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-k,0,ngllz-1-j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-k,0,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-k,0,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-k,0,ngllz-1-j)
                            endif
                        enddo
                    enddo
                case (2)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1,nglly-1-k,ngllz-1-j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1,nglly-1-k,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1,nglly-1-k,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1,nglly-1-k,ngllz-1-j)
                            endif
                        enddo
                    enddo
                case (3)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-k,nglly-1,ngllz-1-j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-k,nglly-1,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-k,nglly-1,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-k,nglly-1,ngllz-1-j)
                            endif
                        enddo
                    enddo
                case (4)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(0,nglly-1-k,ngllz-1-j,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(0,nglly-1-k,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(0,nglly-1-k,ngllz-1-j)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(0,nglly-1-k,ngllz-1-j)
                            endif
                        enddo
                    enddo
                case (5)
                    do j = 1,ngll1-2
                        do k = 1,ngll2-2
                            Tdomain%sFace(nf)%DumpMass(j,k,:) = Tdomain%sFace(nf)%DumpMass(j,k,:) + &
                                Tdomain%specel(n)%DumpMass(ngllx-1-k,nglly-1-j,ngllz-1,:)
                            if (Tdomain%sFace(nf)%FPML) then
                                Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + &
                                    Tdomain%specel(n)%Ivx(ngllx-1-k,nglly-1-j,ngllz-1)
                                Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + &
                                    Tdomain%specel(n)%Ivy(ngllx-1-k,nglly-1-j,ngllz-1)
                                Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + &
                                    Tdomain%specel(n)%Ivz(ngllx-1-k,nglly-1-j,ngllz-1)
                            endif
                        enddo
                    enddo
                end select
            endif

        endif

    enddo


    return
end subroutine get_Mass_Elem2Face
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
