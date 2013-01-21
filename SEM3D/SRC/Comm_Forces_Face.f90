!>
!! \file Comm_Forces_Face.f90
!! \brief
!!
!<

subroutine Comm_Forces_Face(Tdomain,n,ngll,ngll_F,ngllPML,ngllPML_F)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in) :: n
    integer, intent(inout) :: ngll,ngllPML,ngll_F,ngllPML_F

    integer :: ngll1,ngll2,i,j,k,nf


    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2

        if(Tdomain%sComm(n)%orient_faces(i) == 0)then
            if(Tdomain%sFace(nf)%solid)then   ! solid part
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%Forces(k,j,0:2) = Tdomain%sFace(nf)%Forces(k,j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll2-2
                        do k = 1,Tdomain%sFace(nf)%ngll1-2
                            Tdomain%sFace(nf)%Forces1(k,j,0:2) = Tdomain%sFace(nf)%Forces1(k,j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                            Tdomain%sFace(nf)%Forces2(k,j,0:2) = Tdomain%sFace(nf)%Forces2(k,j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                            Tdomain%sFace(nf)%Forces3(k,j,0:2) = Tdomain%sFace(nf)%Forces3(k,j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                            ngllPML = ngllPML + 1
                        enddo
                    enddo
                endif
            else   ! fluid part
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%ForcesFl(k,j) = Tdomain%sFace(nf)%ForcesFl(k,j) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll2-2
                        do k = 1,Tdomain%sFace(nf)%ngll1-2
                            Tdomain%sFace(nf)%ForcesFl1(k,j) = Tdomain%sFace(nf)%ForcesFl1(k,j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
                            Tdomain%sFace(nf)%ForcesFl2(k,j) = Tdomain%sFace(nf)%ForcesFl2(k,j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
                            Tdomain%sFace(nf)%ForcesFl3(k,j) = Tdomain%sFace(nf)%ForcesFl3(k,j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
                            ngllPML_F = ngllPML_F + 1
                        enddo
                    enddo
                endif
            end if

        else if(Tdomain%sComm(n)%orient_faces(i) == 1)then
            if(Tdomain%sFace(nf)%solid)then   ! solid part
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%Forces(ngll1-1-k,j,0:2) = Tdomain%sFace(nf)%Forces(ngll1-1-k,j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll2-2
                        do k = 1,Tdomain%sFace(nf)%ngll1-2
                            Tdomain%sFace(nf)%Forces1(ngll1-1-k,j,0:2) = Tdomain%sFace(nf)%Forces1(ngll1-1-k,j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                            Tdomain%sFace(nf)%Forces2(ngll1-1-k,j,0:2) = Tdomain%sFace(nf)%Forces2(ngll1-1-k,j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                            Tdomain%sFace(nf)%Forces3(ngll1-1-k,j,0:2) = Tdomain%sFace(nf)%Forces3(ngll1-1-k,j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                            ngllPML = ngllPML + 1
                        enddo
                    enddo
                endif
            else   ! fluid part
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%ForcesFl(ngll1-1-k,j) = Tdomain%sFace(nf)%ForcesFl(ngll1-1-k,j) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll2-2
                        do k = 1,Tdomain%sFace(nf)%ngll1-2
                            Tdomain%sFace(nf)%ForcesFl1(ngll1-1-k,j) = Tdomain%sFace(nf)%ForcesFl1(ngll1-1-k,j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
                            Tdomain%sFace(nf)%ForcesFl2(ngll1-1-k,j) = Tdomain%sFace(nf)%ForcesFl2(ngll1-1-k,j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
                            Tdomain%sFace(nf)%ForcesFl3(ngll1-1-k,j) = Tdomain%sFace(nf)%ForcesFl3(ngll1-1-k,j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
                            ngllPML_F = ngllPML_F + 1
                        enddo
                    enddo
                endif
            end if

        else if(Tdomain%sComm(n)%orient_faces(i) == 2)then
            if(Tdomain%sFace(nf)%solid)then  ! solid part
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%Forces(k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces(k,ngll2-1-j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll2-2
                        do k = 1,Tdomain%sFace(nf)%ngll1-2
                            Tdomain%sFace(nf)%Forces1(k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces1(k,ngll2-1-j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                            Tdomain%sFace(nf)%Forces2(k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces2(k,ngll2-1-j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                            Tdomain%sFace(nf)%Forces3(k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces3(k,ngll2-1-j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                            ngllPML = ngllPML + 1
                        enddo
                    enddo
                endif
            else  ! fluid part
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%ForcesFl(k,ngll2-1-j) = Tdomain%sFace(nf)%ForcesFl(k,ngll2-1-j) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll2-2
                        do k = 1,Tdomain%sFace(nf)%ngll1-2
                            Tdomain%sFace(nf)%ForcesFl1(k,ngll2-1-j) = Tdomain%sFace(nf)%ForcesFl1(k,ngll2-1-j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
                            Tdomain%sFace(nf)%ForcesFl2(k,ngll2-1-j) = Tdomain%sFace(nf)%ForcesFl2(k,ngll2-1-j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
                            Tdomain%sFace(nf)%ForcesFl3(k,ngll2-1-j) = Tdomain%sFace(nf)%ForcesFl3(k,ngll2-1-j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
                            ngllPML_F = ngllPML_F + 1
                        enddo
                    enddo
                endif

            end if

        else if(Tdomain%sComm(n)%orient_faces(i) == 3)then
            if(Tdomain%sFace(nf)%solid)then   ! solid part
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%Forces(ngll1-1-k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces(ngll1-1-k,ngll2-1-j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll2-2
                        do k = 1,Tdomain%sFace(nf)%ngll1-2
                            Tdomain%sFace(nf)%Forces1(ngll1-1-k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces1(ngll1-1-k,ngll2-1-j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                            Tdomain%sFace(nf)%Forces2(ngll1-1-k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces2(ngll1-1-k,ngll2-1-j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                            Tdomain%sFace(nf)%Forces3(ngll1-1-k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces3(ngll1-1-k,ngll2-1-j,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                            ngllPML = ngllPML + 1
                        enddo
                    enddo
                endif
            else    ! fluid part
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%ForcesFl(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%ForcesFl(ngll1-1-k,ngll2-1-j) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll2-2
                        do k = 1,Tdomain%sFace(nf)%ngll1-2
                            Tdomain%sFace(nf)%ForcesFl1(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%ForcesFl1(ngll1-1-k,ngll2-1-j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
                            Tdomain%sFace(nf)%ForcesFl2(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%ForcesFl2(ngll1-1-k,ngll2-1-j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
                            Tdomain%sFace(nf)%ForcesFl3(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%ForcesFl3(ngll1-1-k,ngll2-1-j) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
                            ngllPML_F = ngllPML_F + 1
                        enddo
                    enddo
                endif
            end if

        else if(Tdomain%sComm(n)%orient_faces(i) == 4)then
            if(Tdomain%sFace(nf)%solid)then   ! solid part
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll1-2
                        do k = 1,Tdomain%sFace(nf)%ngll2-2
                            Tdomain%sFace(nf)%Forces1(j,k,0:2) = Tdomain%sFace(nf)%Forces1(j,k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                            Tdomain%sFace(nf)%Forces2(j,k,0:2) = Tdomain%sFace(nf)%Forces2(j,k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                            Tdomain%sFace(nf)%Forces3(j,k,0:2) = Tdomain%sFace(nf)%Forces3(j,k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                            ngllPML = ngllPML + 1
                        enddo
                    enddo
                endif
            else   ! fluid part
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%ForcesFl(j,k) = Tdomain%sFace(nf)%ForcesFl(j,k) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll1-2
                        do k = 1,Tdomain%sFace(nf)%ngll2-2
                            Tdomain%sFace(nf)%ForcesFl1(j,k) = Tdomain%sFace(nf)%ForcesFl1(j,k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
                            Tdomain%sFace(nf)%ForcesFl2(j,k) = Tdomain%sFace(nf)%ForcesFl2(j,k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
                            Tdomain%sFace(nf)%ForcesFl3(j,k) = Tdomain%sFace(nf)%ForcesFl3(j,k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
                            ngllPML_F = ngllPML_F + 1
                        enddo
                    enddo
                endif
            end if

        else if(Tdomain%sComm(n)%orient_faces(i) == 5)then
            if(Tdomain%sFace(nf)%solid)then  ! solid part
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%Forces(ngll1-1-j,k,0:2) = Tdomain%sFace(nf)%Forces(ngll2-1-j,k,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll1-2
                        do k = 1,Tdomain%sFace(nf)%ngll2-2
                            Tdomain%sFace(nf)%Forces1(ngll1-1-j,k,0:2) = Tdomain%sFace(nf)%Forces1(ngll1-1-j,k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                            Tdomain%sFace(nf)%Forces2(ngll1-1-j,k,0:2) = Tdomain%sFace(nf)%Forces2(ngll1-1-j,k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                            Tdomain%sFace(nf)%Forces3(ngll1-1-j,k,0:2) = Tdomain%sFace(nf)%Forces3(ngll1-1-j,k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                            ngllPML = ngllPML + 1
                        enddo
                    enddo
                endif
            else  ! fluid part
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%ForcesFl(ngll1-1-j,k) = Tdomain%sFace(nf)%ForcesFl(ngll2-1-j,k) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll1-2
                        do k = 1,Tdomain%sFace(nf)%ngll2-2
                            Tdomain%sFace(nf)%ForcesFl1(ngll1-1-j,k) = Tdomain%sFace(nf)%ForcesFl1(ngll1-1-j,k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
                            Tdomain%sFace(nf)%ForcesFl2(ngll1-1-j,k) = Tdomain%sFace(nf)%ForcesFl2(ngll1-1-j,k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
                            Tdomain%sFace(nf)%ForcesFl3(ngll1-1-j,k) = Tdomain%sFace(nf)%ForcesFl3(ngll1-1-j,k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
                            ngllPML_F = ngllPML_F + 1
                        enddo
                    enddo
                endif

            end if

        else if(Tdomain%sComm(n)%orient_faces(i) == 6)then
            if(Tdomain%sFace(nf)%solid)then  ! solid part
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%Forces(j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces(j,ngll2-1-k,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll1-2
                        do k = 1,Tdomain%sFace(nf)%ngll2-2
                            Tdomain%sFace(nf)%Forces1(j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces1(j,ngll2-1-k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                            Tdomain%sFace(nf)%Forces2(j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces2(j,ngll2-1-k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                            Tdomain%sFace(nf)%Forces3(j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces3(j,ngll2-1-k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                            ngllPML = ngllPML + 1
                        enddo
                    enddo
                endif
            else   ! fluid part
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%ForcesFl(j,ngll2-1-k) = Tdomain%sFace(nf)%ForcesFl(j,ngll2-1-k) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll1-2
                        do k = 1,Tdomain%sFace(nf)%ngll2-2
                            Tdomain%sFace(nf)%ForcesFl1(j,ngll2-1-k) = Tdomain%sFace(nf)%ForcesFl1(j,ngll2-1-k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML,1)
                            Tdomain%sFace(nf)%ForcesFl2(j,ngll2-1-k) = Tdomain%sFace(nf)%ForcesFl2(j,ngll2-1-k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML,2)
                            Tdomain%sFace(nf)%ForcesFl3(j,ngll2-1-k) = Tdomain%sFace(nf)%ForcesFl3(j,ngll2-1-k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML,3)
                            ngllPML_F = ngllPML_F + 1
                        enddo
                    enddo
                endif

            end if

        else if(Tdomain%sComm(n)%orient_faces(i) == 7)then
            if(Tdomain%sFace(nf)%solid)then   ! solid part
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%Forces(ngll1-1-j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces(ngll1-1-j,ngll2-1-k,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                        ngll = ngll + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll1-2
                        do k = 1,Tdomain%sFace(nf)%ngll2-2
                            Tdomain%sFace(nf)%Forces1(ngll1-1-j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces1(ngll1-1-j,ngll2-1-k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,1,0:2)
                            Tdomain%sFace(nf)%Forces2(ngll1-1-j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces2(ngll1-1-j,ngll2-1-k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,2,0:2)
                            Tdomain%sFace(nf)%Forces3(ngll1-1-j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces3(ngll1-1-j,ngll2-1-k,0:2) + &
                                Tdomain%sComm(n)%TakeForcesPML(ngllPML,3,0:2)
                            ngllPML = ngllPML + 1
                        enddo
                    enddo
                endif
            else   ! fluid part
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%ForcesFl(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%ForcesFl(ngll1-1-j,ngll2-1-k) + Tdomain%sComm(n)%TakeForcesFl(ngll_F)
                        ngll_F = ngll_F + 1
                    enddo
                enddo
                if(Tdomain%sFace(nf)%PML)then
                    do j = 1,Tdomain%sFace(nf)%ngll1-2
                        do k = 1,Tdomain%sFace(nf)%ngll2-2
                            Tdomain%sFace(nf)%ForcesFl1(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%ForcesFl1(ngll1-1-j,ngll2-1-k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,1)
                            Tdomain%sFace(nf)%ForcesFl2(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%ForcesFl2(ngll1-1-j,ngll2-1-k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,2)
                            Tdomain%sFace(nf)%ForcesFl3(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%ForcesFl3(ngll1-1-j,ngll2-1-k) + &
                                Tdomain%sComm(n)%TakeForcesPMLFl(ngllPML_F,3)
                            ngllPML_F = ngllPML_F + 1
                        enddo
                    enddo
                endif

            end if
        endif

    enddo

    return
end subroutine Comm_Forces_Face
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
