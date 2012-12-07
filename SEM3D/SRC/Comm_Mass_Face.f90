subroutine Comm_Mass_Face (Tdomain,n,ngll,ngllPML)

    ! Modified by Elise Delavaud 08/02/2006


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
                    Tdomain%sFace(nf)%MassMat(k,j) = Tdomain%sFace(nf)%MassMat(k,j) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%DumpMass(k,j,0:2) = Tdomain%sFace(nf)%DumpMass(k,j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%Ivx(k,j) = Tdomain%sFace(nf)%Ivx(k,j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%Ivy(k,j) = Tdomain%sFace(nf)%Ivy(k,j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%Ivz(k,j) = Tdomain%sFace(nf)%Ivz(k,j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 1 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-k,j) = Tdomain%sFace(nf)%MassMat(ngll1-1-k,j) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%DumpMass(ngll1-1-k,j,0:2) = Tdomain%sFace(nf)%DumpMass(ngll1-1-k,j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%Ivx(ngll1-1-k,j) = Tdomain%sFace(nf)%Ivx(ngll1-1-k,j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%Ivy(ngll1-1-k,j) = Tdomain%sFace(nf)%Ivy(ngll1-1-k,j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%Ivz(ngll1-1-k,j) = Tdomain%sFace(nf)%Ivz(ngll1-1-k,j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 2 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(k,ngll2-1-j) = Tdomain%sFace(nf)%MassMat(k,ngll2-1-j) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%DumpMass(k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%DumpMass(k,ngll2-1-j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%Ivx(k,ngll2-1-j) = Tdomain%sFace(nf)%Ivx(k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%Ivy(k,ngll2-1-j) = Tdomain%sFace(nf)%Ivy(k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%Ivz(k,ngll2-1-j) = Tdomain%sFace(nf)%Ivz(k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 3 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%MassMat(ngll1-1-k,ngll2-1-j) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll2-2
                    do k = 1,Tdomain%sFace(nf)%ngll1-2
                        Tdomain%sFace(nf)%DumpMass(ngll1-1-k,ngll2-1-j,0:2) = Tdomain%sFace(nf)%DumpMass(ngll1-1-k,ngll2-1-j,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%Ivx(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%Ivx(ngll1-1-k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%Ivy(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%Ivy(ngll1-1-k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%Ivz(ngll1-1-k,ngll2-1-j) = Tdomain%sFace(nf)%Ivz(ngll1-1-k,ngll2-1-j) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 4 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(j,k) = Tdomain%sFace(nf)%MassMat(j,k) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%DumpMass(j,k,0:2) = Tdomain%sFace(nf)%DumpMass(j,k,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%Ivx(j,k) = Tdomain%sFace(nf)%Ivx(j,k) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%Ivy(j,k) = Tdomain%sFace(nf)%Ivy(j,k) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%Ivz(j,k) = Tdomain%sFace(nf)%Ivz(j,k) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 5 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-j,k) = Tdomain%sFace(nf)%MassMat(ngll1-1-j,k) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%DumpMass(ngll1-1-j,k,0:2) = Tdomain%sFace(nf)%DumpMass(ngll1-1-j,k,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%Ivx(ngll1-1-j,k) = Tdomain%sFace(nf)%Ivx(ngll1-1-j,k) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%Ivy(ngll1-1-j,k) = Tdomain%sFace(nf)%Ivy(ngll1-1-j,k) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%Ivz(ngll1-1-j,k) = Tdomain%sFace(nf)%Ivz(ngll1-1-j,k) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 6 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(j,ngll2-1-k) = Tdomain%sFace(nf)%MassMat(j,ngll2-1-k) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%DumpMass(j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%DumpMass(j,ngll2-1-k,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%Ivx(j,ngll2-1-k) = Tdomain%sFace(nf)%Ivx(j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%Ivy(j,ngll2-1-k) = Tdomain%sFace(nf)%Ivy(j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%Ivz(j,ngll2-1-k) = Tdomain%sFace(nf)%Ivz(j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else if ( Tdomain%sComm(n)%orient_faces(i) == 7 ) then
            nf = Tdomain%sComm(n)%faces(i)
            do j = 1,Tdomain%sFace(nf)%ngll1-2
                do k = 1,Tdomain%sFace(nf)%ngll2-2
                    Tdomain%sFace(nf)%MassMat(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%MassMat(ngll1-1-j,ngll2-1-k) + Tdomain%sComm(n)%Take(ngll)
                    ngll = ngll + 1
                enddo
            enddo
            if (Tdomain%sFace(nf)%PML) then
                do j = 1,Tdomain%sFace(nf)%ngll1-2
                    do k = 1,Tdomain%sFace(nf)%ngll2-2
                        Tdomain%sFace(nf)%DumpMass(ngll1-1-j,ngll2-1-k,0:2) = Tdomain%sFace(nf)%DumpMass(ngll1-1-j,ngll2-1-k,0:2) + Tdomain%sComm(n)%TakePML(ngllPML,0:2)
                        if (Tdomain%any_FPML) then
                            Tdomain%sFace(nf)%Ivx(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%Ivx(ngll1-1-j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,3)
                            Tdomain%sFace(nf)%Ivy(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%Ivy(ngll1-1-j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,4)
                            Tdomain%sFace(nf)%Ivz(ngll1-1-j,ngll2-1-k) = Tdomain%sFace(nf)%Ivz(ngll1-1-j,ngll2-1-k) + Tdomain%sComm(n)%TakePML(ngllPML,5)
                        endif
                        ngllPML = ngllPML + 1
                    enddo
                enddo
            endif

        else
            print*,'Pb with coherency number for face'

        endif

    enddo



    return
end subroutine Comm_Mass_Face
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
