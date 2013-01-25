program mesher_mat

    !- Building of mesh files compatible with SEM 3D

    use mesh2spec
    implicit none

    character, dimension(:), allocatable  ::  matarray
    real, allocatable, dimension(:)       :: rho, Pspeed, Sspeed, Qpression, Qmu
    integer                   :: n_parts, pml_top, pml_bool, pml_bottom, n_mat, n_mat_tot, i_err,   &
        nn_mat, NGLLx, NGLLy, NGLLz, NGLLxP, NGLLyP, NGLLzP, npow, mat_b, i
    real                      :: freq, Apow, dt
    logical                   :: filter
    logical, parameter        :: VRAI = .true., FAUX = .false.
    integer, parameter        :: NMAX_PROCS = 1024

    !- PML present or not?
    pml_bool = 0 ; pml_top = 0 ; pml_bottom = 0

    call system("clear")

    write(*,*) "***********************************************************************"
    write(*,*) "----------  Construction of input files for SEM runs (3D)  -----------"
    write(*,*) "***********************************************************************"
    write(*,*)
    write(*,*) " -> How many procs for the run?"
    read(*,*) n_parts
    if(n_parts < 1 .or. n_parts > NMAX_PROCS)   &
        stop "In mesher_mat: sure, for the proc number?"
    write(*,*)
    write(*,*) " -> Construction of the material model and mesh file? (0: no, 1: yes) "
    read(*,*), mat_b
    if(mat_b /= 0 .and. mat_b /= 1) stop "Wrong character.."

    if(mat_b == 0)then  ! material model already constructed
        write(*,*) "    --> Material file 'material.input' already exists."
        write(*,*) "      --> construction of the array for material nature (fluid, solid, PML,..)."
        open(10,file="material.input",action="read",status="old",iostat=i_err)
        if(i_err > 0) stop "File 'material.input' not found in the working directory"
        read(10,*) ; read(10,*)
        read(10,*) n_mat
        close(10)
        allocate(matarray(0:n_mat-1))
        call mat_table_construct(matarray)
        !- SEM mesh file generated
        call gen_mesh(matarray,n_parts)
    else   ! we have to construct the material model
        write(*,*)
        write(*,*) "**********************************************"
        write(*,*) " --> Construction of the material model: the file 'mat.dat' must be in the working directory."
        write(*,*)
        open(10,file="mat.dat",action="read",status="old")
        !- number of materials
        read(10,*) n_mat
        n_mat_tot = n_mat
        !- PMLs
        write(*,*) "  --> Number of non-PML materials:",n_mat
        read(10,*) pml_bool
        if(pml_bool /= 0 .and. pml_bool /= 1) stop "In mesh2spec: PML or not?"
        if(pml_bool == 1)then   ! PMLs added
            write(*,*) "  --> PMLs added."
            n_mat_tot = n_mat_tot+8
            read(10,*) pml_top,pml_bottom
            if(pml_top /= 0 .and. pml_top /= 1) stop "In mesh2spec: PML on top or not?"
            if(pml_bottom /= 0 .and. pml_bottom /= 1) stop "In mesh2spec: PML at the bottom or not?"
            if(pml_top == 1)then
                write(*,*) "    --> PMLs on the top."
                n_mat_tot = n_mat_tot+9
            end if
            if(pml_bottom == 1)then
                write(*,*) "    --> PMLs at the bottom."
                n_mat_tot = n_mat_tot+9
            end if
        else
            read(10,*)
            write(*,*)
            write(*,*) "  --> Run without absorbing boundaries: you may observe some reflected phases."
            write(*,*) "       Ok? (type enter if..)"
            read*
        end if
        !- construction of the material table
        allocate(matarray(0:n_mat_tot-1),rho(0:n_mat_tot-1),Pspeed(0:n_mat_tot-1),Sspeed(0:n_mat_tot-1),    &
            Qpression(0:n_mat_tot-1),Qmu(0:n_mat_tot-1))
        write(*,*)
        write(*,*) "************************************************"
        write(*,*) "  --> Construction of the material table: OK."
        do i = 0,n_mat-1
            read(10,"(a1)") matarray(i)
        end do
        close(10)

        !- SEM mesh file generated
        call gen_mesh(matarray,n_parts,pml_bool,pml_top,pml_bottom,n_mat)

        deallocate(rho,Pspeed,Sspeed,Qpression,Qmu)

    end if

    deallocate(matarray)

!!!!!!!!!!!!!!!!!!!!!!
contains
    !-------------------------------------
    subroutine mat_table_construct(mattab)
        !- obtention of the material table array: important for the fluid/solid interfaces in mesh2spec
        character, intent(out)  :: mattab(0:)
        integer                 :: i

        open(10,file="material.input",action="read",status="old")
        read(10,*) ; read(10,*) ; read(10,*)
        do i = 0,size(mattab)-1
            read(10,"(a1)") mattab(i)
        end do
        close(10)

    end subroutine mat_table_construct
    !-------------------------------------
end program mesher_mat
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
