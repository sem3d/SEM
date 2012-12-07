!>
!!\file plot_grid.F90
!!\brief Assure la sortie des coordonnées des noeuds.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type(domain), target, intent (IN) Tdomain
!<


subroutine plot_grid (Tdomain)
    use sdomain
    use semdatafiles

    implicit none
    type(domain), target, intent (IN) :: Tdomain
    integer :: n,i,i_aus,j, i_proc
    real :: x0,z0
    character (len=MAX_FILE_SIZE) :: fnamef

    ! Plotting the grid

    i_proc = Tdomain%Mpi_var%my_rank

    call semname_plot_grid_grid(i_proc,fnamef)

    open (22,file=fnamef,form="formatted",status="unknown")
    write (22,"(a2)") ">>"
    do n = 0,Tdomain%n_elem-1
        do i = 0,3
            i_aus = Tdomain%specel(n)%Control_Nodes(i)
            x0 = Tdomain%Coord_Nodes(0,i_aus)
            z0 = Tdomain%Coord_Nodes(1,i_aus)
            write (22,*) x0, z0
        enddo
        write (22,"(a2)") ">>"
    enddo
    close (22)

    if (Tdomain%logicD%super_object_local_present) then

        call semname_plot_grid_surface(i_proc,fnamef)

        open (23,file=fnamef, form="formatted",status="unknown")
        write (23,"(a2)") ">>"
        do n = 0, Tdomain%n_fault-1
            do i = 0, Tdomain%sFault(n)%n_face-1
                do j = 0,1
                    x0 = Tdomain%sFault(n)%Fface(i)%X_Vertex(j)
                    z0 = Tdomain%sFault(n)%Fface(i)%Z_Vertex(j)
                    write (23,*) x0, z0
                enddo
                write (23,"(a2)") ">>"
            enddo
        enddo
    endif
    close (23)
    return
end subroutine plot_grid
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
