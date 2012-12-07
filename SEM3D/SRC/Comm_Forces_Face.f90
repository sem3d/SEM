!>
!! \file Comm_Forces_Face.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine Comm_Forces_Face (Tdomain,n,ngll,ngllPML)

    use sdomain

    implicit none

    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: n
    integer, intent (INOUT) :: ngll,ngllPML

    integer :: ngll1,ngll2,i,j,k,nf


    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        ngll1 = Tdomain%sFace(nf)%ngll1
        ngll2 = Tdomain%sFace(nf)%ngll2

        if ( Tdomain%sComm(n)%orient_faces(i) == 0 ) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%Forces(k,j,0:2) = Tdomain%sFace(nf)%Forces(k,j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
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

        else if ( Tdomain%sComm(n)%orient_faces(i) == 1 ) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%Forces(ngll1-1-k,j,0:2) = Tdomain%sFace(nf)%Forces(ngll1-1-k,j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
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

        else if ( Tdomain%sComm(n)%orient_faces(i) == 2 ) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%Forces(k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces(k,ngll2-1-j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
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

        else if ( Tdomain%sComm(n)%orient_faces(i) == 3 ) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%Forces(ngll1-1-k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%Forces(ngll1-1-k,ngll2-1-j,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
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

        else if ( Tdomain%sComm(n)%orient_faces(i) == 4 ) then
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%Forces(j,k,0:2) = Tdomain%sFace(nf)%Forces(j,k,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
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

        else if ( Tdomain%sComm(n)%orient_faces(i) == 5 ) then
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%Forces(ngll1-1-j,k,0:2) = Tdomain%sFace(nf)%Forces(ngll2-1-j,k,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
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

        else if ( Tdomain%sComm(n)%orient_faces(i) == 6 ) then
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%Forces(j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces(j,ngll2-1-k,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
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

        else if ( Tdomain%sComm(n)%orient_faces(i) == 7 ) then
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%Forces(ngll1-1-j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%Forces(ngll1-1-j,ngll2-1-k,0:2) + Tdomain%sComm(n)%TakeForces(ngll,0:2)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
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

        endif

    enddo


    return
end subroutine Comm_Forces_Face
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
