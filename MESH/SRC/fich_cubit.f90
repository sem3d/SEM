module fich_cubit


contains

    !-----------------------------------------
    subroutine lec1_cubit(indunit,ln,n_mat)

        implicit none

        character(len=70)    :: ligne
        integer, intent(in)  :: indunit
        integer, intent(out) :: ln, n_mat
        integer              :: i


        i = 0

        ! counting of global nodes
        do
            i = i+1
            read(indunit,"(a)",end=100) ligne
            if(ligne(1:2) == "*N")then
                i = 0
            end if
            if(ligne(1:2) == "*E") then
                ln = i-1
                exit
            endif
        end do


        ! Number of materials
100     n_mat = 0
        rewind(indunit)
        do
            read(indunit,"(a)",end=101) ligne
            if(ligne(1:2) == "*E") n_mat = n_mat+1
        end do

101     return

    end subroutine lec1_cubit
    !--------------------------------------------
    subroutine lec2_cubit(indunit,n_elem,elem,ind_mat,n_neu,n_PW)

        implicit none

        integer, intent(in) :: indunit
        character(len=70)   :: ligne
        character           :: c
        integer   :: l,i,n_mat,nn
        logical  :: bound_cond
        integer, intent(out) :: n_elem, elem(:), ind_mat(:)
        integer, intent(out) :: n_neu,n_PW


        n_elem = 0
        bound_cond = .false.
        n_neu = 0 ; n_PW = 0

        ! Number of elements per block
        elem(:) = 0
        nn = 0
        rewind(indunit)
        do
            read(indunit,"(a)",end=102) ligne
            if(ligne(1:2) == "*E") nn = nn+1  ! material increment
            if(nn < 1) cycle
            ! Neumann (traction) or Dirichlet boundary conditions or plane wave interface
            if(ligne(1:2) == "*T" .or. ligne(1:2) == "*D" .or. ligne(1:2) == "*P")then
                bound_cond = .true.
                exit
            end if
            elem(nn) = elem(nn)+1
        end do

102     elem(:) = elem(:)-1

        ! total number of elements
        n_elem = SUM(elem)


        ! Indices of materials
        ind_mat(:) = 0
        rewind(indunit)
        nn = 0
        do
            read(indunit,"(a)",end=121) ligne
            if(ligne(1:2) == "*E")then
                nn = nn+1
                l = len_trim(ligne)
                c = ligne(l-1:l-1)
                i = iachar(c)-iachar('0')
                if(i>=1 .and. i<=9) ind_mat(nn) = 10*i
                c = ligne(l:l)
                i = iachar(c)-iachar('0')
                if(i>=1 .and. i<=9) ind_mat(nn) = ind_mat(nn)+i
            end if
        end do


        ! eventual boundary conditions
121     if(bound_cond)then
            ! Neumann
            rewind(indunit)
            do
                read(indunit,"(a)",end=123) ligne
                if(ligne(1:2) == "*T")then
                    print*,"  --> Neumann B.C. imposed."
                    n_neu = 0
                    cycle
                end if
                if(ligne(1:2) == "*D" .or. ligne(1:2) == "*P") exit
                n_neu = n_neu+1
            end do
            ! Plane wave
123         rewind(indunit)
            do
                read(indunit,"(a)",end=122) ligne
                if(ligne(1:2) == "*P")then
                    print*,"  --> Plane wave imposed on an interface."
                    n_PW = 0
                    cycle
                end if
                if(ligne(1:2) == "*D" .or. (ligne(1:2) == "*T" .and. n_neu == 0)) exit
                n_PW = n_PW+1
            end do

        end if

        ! now we can rewind to read node and element properties
122     rewind(indunit)

        return

    end subroutine lec2_cubit
    !---------------------------------------------------------
    !---------------------------------------------------------
    subroutine lec_neu_cubit(indunit,n_neu,Faces_Neu)
        ! reading of nodes for each face on which Neumann BC are applied
        implicit none
        integer, intent(in)    :: indunit,n_neu
        integer, dimension(0:3,0:n_neu-1), intent(out)  :: Faces_Neu
        integer   :: i,j,k
        character(len=70)   :: ligne

        rewind(indunit)
        do
            read(indunit,"(a)") ligne
            if(ligne(1:2) == "*T") exit
        end do
        do i = 0,n_neu-1
            read(indunit,*) k,(Faces_Neu(j,i),j=0,3)
            ! changing index
            Faces_Neu(0:,i) = Faces_Neu(0:,i)-1
        end do

    end subroutine lec_neu_cubit
    !---------------------------------------------------------
end module fich_cubit
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
