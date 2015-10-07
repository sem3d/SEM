!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module mesh_properties

    use fich_cubit
    use fich_unv
    use sets, only : sort
    use mesh_earthchunk

    implicit none

contains


    !---------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------

    subroutine mesh_init_3D(n_nods,n_points,n_elem,n_blocks,          &
        xco,yco,zco,Ipointer,Material,tabmat,n_neu,n_PW,     &
        Faces_Neumann, Faces_PW)
        implicit none
        integer, intent(out)   :: n_nods, n_points, n_elem, n_blocks, n_neu, &
            n_PW
        integer, allocatable, dimension(:), intent(out)   :: Material
        integer, allocatable, dimension(:,:), intent(out) :: Ipointer,    &
            Faces_Neumann, Faces_PW
        character, dimension(:), allocatable, intent(out)              :: tabmat
        real, allocatable, dimension(:), intent(out)   :: xco,yco,zco
        real, allocatable, dimension(:)   :: thick,z_dep
        integer   :: choice, i_ex, idummy,icount,iunit, nfile, i, j, n,    &
            mesh_type,pml_b, pml_t, pml_bottom, n_mat_tot,nmatref
        integer, allocatable, dimension(:)  :: ind_mat, n_elem_mat, nlay
        real      :: xmin,xmax,ymin,ymax,zmin,zmax,step_x,step_y,   &
            xminref,xmaxref,yminref,ymaxref,zminref,zmaxref
        integer :: n_mate
        character(len=30)    :: cubitmesh
        character(len=60), dimension(:), allocatable  :: unv_files


        type(EarthChunk_t) :: earthchunk

        ! no Neumann or Plane Wave faces for the time being.
        n_neu = 0 ; n_PW = 0

        ! menu
        choice = -1
        do while(choice < 0 .or. choice > 5)
            write(*,*)
            write(*,*) "****************************************"
            write(*,*) "****************************************"
            write(*,*) " WHICH INITIAL MESH?"
            write(*,*) "     1- On the fly"
            write(*,*) "     2- Abaqus from Cubit"
            write(*,*) "     3- Ideas (.unv) files"
            write(*,*) "     4- HDF5 Hex8 files"
            write(*,*) "     5- Earth Chunk"
            read(*,*) choice
            write(*,*)
            write(*,*) "****************************************"
        end do

        select case(choice)
        case(1)  ! mesh created on the fly
            ! general information input file 
            !- construction of the material table
            call  create_model_ondafly(n_mat_tot,tabmat,pml_b,pml_t,pml_bottom,   &
                xmin,xmax,step_x,ymin,ymax,step_y,zmax,mesh_type,nlay,thick)
            write(*,*)
            write(*,*) "************************************************"
            write(*,*) "  --> Construction of the material table: OK."

            call init_ondafly(xmin,xmax,ymin,ymax,zmax,xminref,xmaxref,   &
                yminref,ymaxref,zminref,zmaxref, step_x,step_y,  &
                mesh_type,n_nods,n_points,n_elem,pml_b,   &
                pml_t,pml_bottom,nlay,thick,z_dep)
            allocate(xco(0:n_points-1),yco(0:n_points-1),zco(0:n_points-1))
            allocate(Ipointer(0:n_nods-1,0:n_elem-1))
            allocate(Material(0:n_elem-1))
            n_blocks = size(tabmat)

            !- just for Solid/Fluid trials..
            !     Material(0:) = 0
            !     Material(5) = 1
            !     Material(6) = 1
            !     Material(4) = 1
            !     Material(7) = 1

            call mesh_on_the_fly(xmin,xmax,ymin,ymax,zmin,zmax,step_x,step_y,  &
                z_dep,mesh_type,xco,yco,zco,Ipointer)

            !- PML materials added
            call nature_elem(Ipointer,xco,yco,zco,size(tabmat),Material,xminref,    &
                xmaxref,yminref,ymaxref,zminref,zmaxref,pml_b,pml_t,pml_bottom,thick)

            !- Cubit file
        case(2)
            write(*,*)
            write(*,*) "------------------------------------"
            write(*,*) "  --> CUBIT file to be analyzed."
            write(*,*) "------------------------------------"
            write(*,*)
            write(*,*) "    --> Warning: there must be no empty line in the Cubit file (depending on the version)"
            write(*,*) "    --> Warning 2: file mater.in must be in the directory to construct material table."
            write(*,*)
            ! construction of the material table
            call mat_table_construct(tabmat)
            write(*,*) "  --> Name of the Cubit file:"
            read(*,*) cubitmesh
            iunit = 10
            open(iunit,file=cubitmesh,status="old",iostat=i_ex)
            if(i_ex > 0) stop " Cubit file does not exist"
            write(*,*)
            write(*,*) "  --> Number of control nodes per element: "
            read(*,*) n_nods
            if(n_nods /= 8) stop "In mesh_properties:    &
                & Only elements with 8 nodes for the time being."
            write(*,*)
            ! some changes, relative to Elise Delavaud's version: here we do not add
            !   informations, in the Cubit file, but get them by analysing the
            !   Cubit file structure.
            call lec1_cubit(iunit,n_points,n_blocks)
            if(n_blocks /= size(tabmat)) stop "In mesh_properties : &
                & Incoherency between mesh and material files."
            allocate(ind_mat(n_blocks),n_elem_mat(n_blocks))
            call lec2_cubit(iunit,n_elem,n_elem_mat,ind_mat,n_neu,n_PW)
            read(10,*)   ! Header
            read(10,*)   ! Name of the mesh
            read(10,*)   ! *NODE
            write(*,*) "  --> Number of control points, elements, materials: "
            write(*,*) "    ",n_points, n_elem, n_blocks
            write(*,*)
            !- co-ordinates of control points
            allocate(xco(0:n_points-1),yco(0:n_points-1),zco(0:n_points-1))
            do i = 0,n_points-1
                read(10,*) idummy,xco(i),yco(i),zco(i)
            enddo

            !- general index for each control point of an element
            allocate(Ipointer(0:n_nods-1,0:n_elem-1))
            allocate(Material(0:n_elem-1))

            do i = 0,n_blocks-1
                read(10,*)    ! *ELEMENTS
                do n = 0,n_elem_mat(i+1)-1
                    read(10,*) icount,(Ipointer(j,icount-1),j=0,n_nods-1)
                    Material(icount-1) = ind_mat(i+1) - 1
                end do
            end do
            !- eventual Neumann BC
            if(n_neu > 0)then
                allocate(Faces_Neumann(0:3,0:n_neu-1))
                call lec_neu_cubit(iunit,n_neu,Faces_Neumann)
            end if
            !- end of reading Cubit file
            close(10)
            deallocate(ind_mat,n_elem_mat)
            !- transfo: initial mesh with a given orientation/element => SEM orientation
            write(*,*) "  -->  Transformation -> SEM element nodes' orientation"
            call transfo_8(n_elem,n_points,Ipointer,xco,yco,zco)
            write(*,*)
            write(*,*) "  - END of Cubit-style file READING -"
            write(*,*) "****************************************"
            write(*,*)

            ! UNV files
        case(3,4)
            write(*,*) "****************************************"
            write(*,*) "  --> files to be read"
            write(*,*) "    --> How many files ?"
            read*, nfile
            n_blocks = nfile   ! number of materials
            allocate(unv_files(0:nfile-1))
            call mat_table_construct_unv(tabmat)
            call lec_init_unv(unv_files)

            n_mate = size(tabmat,1)
            n_nods = 8
            if (choice==3) then
                call lec_unv(unv_files,n_points,n_elem,Material,Ipointer,xco,yco,zco, n_blocks, n_mate)
            else
                call lec_hdf5(unv_files,n_points,n_elem,Material,Ipointer,xco,yco,zco, n_blocks)
            endif

            write(*,*) "***********************************"
            write(*,*) "  - END of files READING -"
            write(*,*) "***********************************"
            deallocate(unv_files)

        case(5)

            call init_earthchunk(earthchunk)

            n_points = earthchunk%total_pt
            n_elem = earthchunk%total_elem
            n_nods = earthchunk%nods
            allocate(xco(0:n_points-1),yco(0:n_points-1),zco(0:n_points-1))
            allocate(Ipointer(0:n_nods-1,0:n_elem-1))
            allocate(Material(0:n_elem-1))

            pml_b = 1  ! PMLs in the domain: lateral and bottom
            nmatref = 1   ! only one physical material
            n_blocks = 1 + 9+8 ! physical and PMLs
            allocate(tabmat(0:17))
            tabmat(0) = 'S' ; tabmat(1:17) = 'P'

            call create_earthchunk(earthchunk, pml_b, nmatref, xco,yco,zco, Ipointer, Material)

            call clean_earthchunk(earthchunk)


        case default
            stop " Please make a choice for the initial meshing."
        end select

    end subroutine mesh_init_3D
    !-------------------------------------
    !-------------------------------------
    subroutine transfo_8(nelem,npts,Ipointer,xco,yco,zco)
        !- modifies the indexing in each element: from Cubit to SEM orientation and convention
        implicit none
        integer, intent(in)   :: nelem,npts
        integer, intent(out)  :: Ipointer(0:7,0:nelem-1)
        real, intent(in)      :: xco(0:npts-1),yco(0:npts-1),zco(0:npts-1)
        integer  :: icount, IpointerN(0:7)

        write(*,*) "        -- Transformation in progress  --"
        do icount = 1,nelem
            if((xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ))     then
                IpointerN(0) = Ipointer(0,icount-1)
                IpointerN(1) = Ipointer(1,icount-1)
                IpointerN(2) = Ipointer(5,icount-1)
                IpointerN(3) = Ipointer(4,icount-1)
                IpointerN(4) = Ipointer(3,icount-1)
                IpointerN(5) = Ipointer(2,icount-1)
                IpointerN(6) = Ipointer(6,icount-1)
                IpointerN(7) = Ipointer(7,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if((xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(3,icount-1)
                IpointerN(1) = Ipointer(0,icount-1)
                IpointerN(2) = Ipointer(1,icount-1)
                IpointerN(3) = Ipointer(2,icount-1)
                IpointerN(4) = Ipointer(7,icount-1)
                IpointerN(5) = Ipointer(4,icount-1)
                IpointerN(6) = Ipointer(5,icount-1)
                IpointerN(7) = Ipointer(6,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(4,icount-1)
                IpointerN(1) = Ipointer(0,icount-1)
                IpointerN(2) = Ipointer(1,icount-1)
                IpointerN(3) = Ipointer(5,icount-1)
                IpointerN(4) = Ipointer(7,icount-1)
                IpointerN(5) = Ipointer(3,icount-1)
                IpointerN(6) = Ipointer(2,icount-1)
                IpointerN(7) = Ipointer(6,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(5,icount-1)
                IpointerN(1) = Ipointer(4,icount-1)
                IpointerN(2) = Ipointer(0,icount-1)
                IpointerN(3) = Ipointer(1,icount-1)
                IpointerN(4) = Ipointer(6,icount-1)
                IpointerN(5) = Ipointer(7,icount-1)
                IpointerN(6) = Ipointer(3,icount-1)
                IpointerN(7) = Ipointer(2,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(2,icount-1)
                IpointerN(1) = Ipointer(3,icount-1)
                IpointerN(2) = Ipointer(0,icount-1)
                IpointerN(3) = Ipointer(1,icount-1)
                IpointerN(4) = Ipointer(6,icount-1)
                IpointerN(5) = Ipointer(7,icount-1)
                IpointerN(6) = Ipointer(4,icount-1)
                IpointerN(7) = Ipointer(5,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(1,icount-1)
                IpointerN(1) = Ipointer(2,icount-1)
                IpointerN(2) = Ipointer(3,icount-1)
                IpointerN(3) = Ipointer(0,icount-1)
                IpointerN(4) = Ipointer(5,icount-1)
                IpointerN(5) = Ipointer(6,icount-1)
                IpointerN(6) = Ipointer(7,icount-1)
                IpointerN(7) = Ipointer(4,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(1,icount-1)
                IpointerN(1) = Ipointer(5,icount-1)
                IpointerN(2) = Ipointer(4,icount-1)
                IpointerN(3) = Ipointer(0,icount-1)
                IpointerN(4) = Ipointer(2,icount-1)
                IpointerN(5) = Ipointer(6,icount-1)
                IpointerN(6) = Ipointer(7,icount-1)
                IpointerN(7) = Ipointer(3,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(1,icount-1)
                IpointerN(1) = Ipointer(0,icount-1)
                IpointerN(2) = Ipointer(3,icount-1)
                IpointerN(3) = Ipointer(2,icount-1)
                IpointerN(4) = Ipointer(5,icount-1)
                IpointerN(5) = Ipointer(4,icount-1)
                IpointerN(6) = Ipointer(7,icount-1)
                IpointerN(7) = Ipointer(6,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(1,icount-1)
                IpointerN(1) = Ipointer(0,icount-1)
                IpointerN(2) = Ipointer(4,icount-1)
                IpointerN(3) = Ipointer(5,icount-1)
                IpointerN(4) = Ipointer(2,icount-1)
                IpointerN(5) = Ipointer(3,icount-1)
                IpointerN(6) = Ipointer(7,icount-1)
                IpointerN(7) = Ipointer(6,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(2,icount-1)
                IpointerN(1) = Ipointer(1,icount-1)
                IpointerN(2) = Ipointer(0,icount-1)
                IpointerN(3) = Ipointer(3,icount-1)
                IpointerN(4) = Ipointer(6,icount-1)
                IpointerN(5) = Ipointer(5,icount-1)
                IpointerN(6) = Ipointer(4,icount-1)
                IpointerN(7) = Ipointer(7,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(5,icount-1)
                IpointerN(1) = Ipointer(1,icount-1)
                IpointerN(2) = Ipointer(0,icount-1)
                IpointerN(3) = Ipointer(4,icount-1)
                IpointerN(4) = Ipointer(6,icount-1)
                IpointerN(5) = Ipointer(2,icount-1)
                IpointerN(6) = Ipointer(3,icount-1)
                IpointerN(7) = Ipointer(7,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(3,icount-1)
                IpointerN(1) = Ipointer(2,icount-1)
                IpointerN(2) = Ipointer(1,icount-1)
                IpointerN(3) = Ipointer(0,icount-1)
                IpointerN(4) = Ipointer(7,icount-1)
                IpointerN(5) = Ipointer(6,icount-1)
                IpointerN(6) = Ipointer(5,icount-1)
                IpointerN(7) = Ipointer(4,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(4,icount-1)
                IpointerN(1) = Ipointer(5,icount-1)
                IpointerN(2) = Ipointer(1,icount-1)
                IpointerN(3) = Ipointer(0,icount-1)
                IpointerN(4) = Ipointer(7,icount-1)
                IpointerN(5) = Ipointer(6,icount-1)
                IpointerN(6) = Ipointer(2,icount-1)
                IpointerN(7) = Ipointer(3,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(0,icount-1)
                IpointerN(1) = Ipointer(3,icount-1)
                IpointerN(2) = Ipointer(2,icount-1)
                IpointerN(3) = Ipointer(1,icount-1)
                IpointerN(4) = Ipointer(4,icount-1)
                IpointerN(5) = Ipointer(7,icount-1)
                IpointerN(6) = Ipointer(6,icount-1)
                IpointerN(7) = Ipointer(5,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(0,icount-1)
                IpointerN(1) = Ipointer(4,icount-1)
                IpointerN(2) = Ipointer(5,icount-1)
                IpointerN(3) = Ipointer(1,icount-1)
                IpointerN(4) = Ipointer(3,icount-1)
                IpointerN(5) = Ipointer(7,icount-1)
                IpointerN(6) = Ipointer(6,icount-1)
                IpointerN(7) = Ipointer(2,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(4,icount-1)
                IpointerN(1) = Ipointer(5,icount-1)
                IpointerN(2) = Ipointer(6,icount-1)
                IpointerN(3) = Ipointer(7,icount-1)
                IpointerN(4) = Ipointer(0,icount-1)
                IpointerN(5) = Ipointer(1,icount-1)
                IpointerN(6) = Ipointer(2,icount-1)
                IpointerN(7) = Ipointer(3,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(3,icount-1)
                IpointerN(1) = Ipointer(2,icount-1)
                IpointerN(2) = Ipointer(6,icount-1)
                IpointerN(3) = Ipointer(7,icount-1)
                IpointerN(4) = Ipointer(0,icount-1)
                IpointerN(5) = Ipointer(1,icount-1)
                IpointerN(6) = Ipointer(5,icount-1)
                IpointerN(7) = Ipointer(4,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(7,icount-1)
                IpointerN(1) = Ipointer(4,icount-1)
                IpointerN(2) = Ipointer(5,icount-1)
                IpointerN(3) = Ipointer(6,icount-1)
                IpointerN(4) = Ipointer(3,icount-1)
                IpointerN(5) = Ipointer(0,icount-1)
                IpointerN(6) = Ipointer(1,icount-1)
                IpointerN(7) = Ipointer(2,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(7,icount-1)
                IpointerN(1) = Ipointer(3,icount-1)
                IpointerN(2) = Ipointer(2,icount-1)
                IpointerN(3) = Ipointer(6,icount-1)
                IpointerN(4) = Ipointer(4,icount-1)
                IpointerN(5) = Ipointer(0,icount-1)
                IpointerN(6) = Ipointer(1,icount-1)
                IpointerN(7) = Ipointer(5,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(6,icount-1)
                IpointerN(1) = Ipointer(7,icount-1)
                IpointerN(2) = Ipointer(4,icount-1)
                IpointerN(3) = Ipointer(5,icount-1)
                IpointerN(4) = Ipointer(2,icount-1)
                IpointerN(5) = Ipointer(3,icount-1)
                IpointerN(6) = Ipointer(0,icount-1)
                IpointerN(7) = Ipointer(1,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(6,icount-1)
                IpointerN(1) = Ipointer(7,icount-1)
                IpointerN(2) = Ipointer(3,icount-1)
                IpointerN(3) = Ipointer(2,icount-1)
                IpointerN(4) = Ipointer(5,icount-1)
                IpointerN(5) = Ipointer(4,icount-1)
                IpointerN(6) = Ipointer(0,icount-1)
                IpointerN(7) = Ipointer(1,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(2,icount-1)
                IpointerN(1) = Ipointer(6,icount-1)
                IpointerN(2) = Ipointer(7,icount-1)
                IpointerN(3) = Ipointer(3,icount-1)
                IpointerN(4) = Ipointer(1,icount-1)
                IpointerN(5) = Ipointer(5,icount-1)
                IpointerN(6) = Ipointer(4,icount-1)
                IpointerN(7) = Ipointer(0,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(5,icount-1)
                IpointerN(1) = Ipointer(6,icount-1)
                IpointerN(2) = Ipointer(7,icount-1)
                IpointerN(3) = Ipointer(4,icount-1)
                IpointerN(4) = Ipointer(1,icount-1)
                IpointerN(5) = Ipointer(2,icount-1)
                IpointerN(6) = Ipointer(3,icount-1)
                IpointerN(7) = Ipointer(0,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(5,icount-1)
                IpointerN(1) = Ipointer(4,icount-1)
                IpointerN(2) = Ipointer(7,icount-1)
                IpointerN(3) = Ipointer(6,icount-1)
                IpointerN(4) = Ipointer(1,icount-1)
                IpointerN(5) = Ipointer(0,icount-1)
                IpointerN(6) = Ipointer(3,icount-1)
                IpointerN(7) = Ipointer(2,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(4,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(2,icount-1)
                IpointerN(1) = Ipointer(3,icount-1)
                IpointerN(2) = Ipointer(7,icount-1)
                IpointerN(3) = Ipointer(6,icount-1)
                IpointerN(4) = Ipointer(1,icount-1)
                IpointerN(5) = Ipointer(0,icount-1)
                IpointerN(6) = Ipointer(4,icount-1)
                IpointerN(7) = Ipointer(5,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(6,icount-1)
                IpointerN(1) = Ipointer(5,icount-1)
                IpointerN(2) = Ipointer(4,icount-1)
                IpointerN(3) = Ipointer(7,icount-1)
                IpointerN(4) = Ipointer(2,icount-1)
                IpointerN(5) = Ipointer(1,icount-1)
                IpointerN(6) = Ipointer(0,icount-1)
                IpointerN(7) = Ipointer(3,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(6,icount-1)
                IpointerN(1) = Ipointer(2,icount-1)
                IpointerN(2) = Ipointer(3,icount-1)
                IpointerN(3) = Ipointer(7,icount-1)
                IpointerN(4) = Ipointer(5,icount-1)
                IpointerN(5) = Ipointer(1,icount-1)
                IpointerN(6) = Ipointer(0,icount-1)
                IpointerN(7) = Ipointer(4,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(7,icount-1)
                IpointerN(1) = Ipointer(6,icount-1)
                IpointerN(2) = Ipointer(5,icount-1)
                IpointerN(3) = Ipointer(4,icount-1)
                IpointerN(4) = Ipointer(3,icount-1)
                IpointerN(5) = Ipointer(2,icount-1)
                IpointerN(6) = Ipointer(1,icount-1)
                IpointerN(7) = Ipointer(0,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(7,icount-1)
                IpointerN(1) = Ipointer(6,icount-1)
                IpointerN(2) = Ipointer(2,icount-1)
                IpointerN(3) = Ipointer(3,icount-1)
                IpointerN(4) = Ipointer(4,icount-1)
                IpointerN(5) = Ipointer(5,icount-1)
                IpointerN(6) = Ipointer(1,icount-1)
                IpointerN(7) = Ipointer(0,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(4,icount-1)
                IpointerN(1) = Ipointer(7,icount-1)
                IpointerN(2) = Ipointer(6,icount-1)
                IpointerN(3) = Ipointer(5,icount-1)
                IpointerN(4) = Ipointer(0,icount-1)
                IpointerN(5) = Ipointer(3,icount-1)
                IpointerN(6) = Ipointer(2,icount-1)
                IpointerN(7) = Ipointer(1,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(3,icount-1)
                IpointerN(1) = Ipointer(7,icount-1)
                IpointerN(2) = Ipointer(6,icount-1)
                IpointerN(3) = Ipointer(2,icount-1)
                IpointerN(4) = Ipointer(0,icount-1)
                IpointerN(5) = Ipointer(4,icount-1)
                IpointerN(6) = Ipointer(5,icount-1)
                IpointerN(7) = Ipointer(1,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(0,icount-1)
                IpointerN(1) = Ipointer(3,icount-1)
                IpointerN(2) = Ipointer(7,icount-1)
                IpointerN(3) = Ipointer(4,icount-1)
                IpointerN(4) = Ipointer(1,icount-1)
                IpointerN(5) = Ipointer(2,icount-1)
                IpointerN(6) = Ipointer(6,icount-1)
                IpointerN(7) = Ipointer(5,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(0,icount-1)
                IpointerN(1) = Ipointer(4,icount-1)
                IpointerN(2) = Ipointer(7,icount-1)
                IpointerN(3) = Ipointer(3,icount-1)
                IpointerN(4) = Ipointer(1,icount-1)
                IpointerN(5) = Ipointer(5,icount-1)
                IpointerN(6) = Ipointer(6,icount-1)
                IpointerN(7) = Ipointer(2,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(4,icount-1)
                IpointerN(1) = Ipointer(0,icount-1)
                IpointerN(2) = Ipointer(3,icount-1)
                IpointerN(3) = Ipointer(7,icount-1)
                IpointerN(4) = Ipointer(5,icount-1)
                IpointerN(5) = Ipointer(1,icount-1)
                IpointerN(6) = Ipointer(2,icount-1)
                IpointerN(7) = Ipointer(6,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(3,icount-1)
                IpointerN(1) = Ipointer(0,icount-1)
                IpointerN(2) = Ipointer(4,icount-1)
                IpointerN(3) = Ipointer(7,icount-1)
                IpointerN(4) = Ipointer(2,icount-1)
                IpointerN(5) = Ipointer(1,icount-1)
                IpointerN(6) = Ipointer(5,icount-1)
                IpointerN(7) = Ipointer(6,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(7,icount-1)
                IpointerN(1) = Ipointer(4,icount-1)
                IpointerN(2) = Ipointer(0,icount-1)
                IpointerN(3) = Ipointer(3,icount-1)
                IpointerN(4) = Ipointer(6,icount-1)
                IpointerN(5) = Ipointer(5,icount-1)
                IpointerN(6) = Ipointer(1,icount-1)
                IpointerN(7) = Ipointer(2,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(7,icount-1)
                IpointerN(1) = Ipointer(3,icount-1)
                IpointerN(2) = Ipointer(0,icount-1)
                IpointerN(3) = Ipointer(4,icount-1)
                IpointerN(4) = Ipointer(6,icount-1)
                IpointerN(5) = Ipointer(2,icount-1)
                IpointerN(6) = Ipointer(1,icount-1)
                IpointerN(7) = Ipointer(5,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(3,icount-1)
                IpointerN(1) = Ipointer(7,icount-1)
                IpointerN(2) = Ipointer(4,icount-1)
                IpointerN(3) = Ipointer(0,icount-1)
                IpointerN(4) = Ipointer(2,icount-1)
                IpointerN(5) = Ipointer(6,icount-1)
                IpointerN(6) = Ipointer(5,icount-1)
                IpointerN(7) = Ipointer(1,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(1,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(2,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(7,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(4,icount-1)
                IpointerN(1) = Ipointer(7,icount-1)
                IpointerN(2) = Ipointer(3,icount-1)
                IpointerN(3) = Ipointer(0,icount-1)
                IpointerN(4) = Ipointer(5,icount-1)
                IpointerN(5) = Ipointer(6,icount-1)
                IpointerN(6) = Ipointer(2,icount-1)
                IpointerN(7) = Ipointer(1,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(1,icount-1)
                IpointerN(1) = Ipointer(2,icount-1)
                IpointerN(2) = Ipointer(6,icount-1)
                IpointerN(3) = Ipointer(5,icount-1)
                IpointerN(4) = Ipointer(0,icount-1)
                IpointerN(5) = Ipointer(3,icount-1)
                IpointerN(6) = Ipointer(7,icount-1)
                IpointerN(7) = Ipointer(4,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(1,icount-1)
                IpointerN(1) = Ipointer(5,icount-1)
                IpointerN(2) = Ipointer(6,icount-1)
                IpointerN(3) = Ipointer(2,icount-1)
                IpointerN(4) = Ipointer(0,icount-1)
                IpointerN(5) = Ipointer(4,icount-1)
                IpointerN(6) = Ipointer(7,icount-1)
                IpointerN(7) = Ipointer(3,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(5,icount-1)
                IpointerN(1) = Ipointer(1,icount-1)
                IpointerN(2) = Ipointer(2,icount-1)
                IpointerN(3) = Ipointer(6,icount-1)
                IpointerN(4) = Ipointer(4,icount-1)
                IpointerN(5) = Ipointer(0,icount-1)
                IpointerN(6) = Ipointer(3,icount-1)
                IpointerN(7) = Ipointer(7,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(2,icount-1)
                IpointerN(1) = Ipointer(1,icount-1)
                IpointerN(2) = Ipointer(5,icount-1)
                IpointerN(3) = Ipointer(6,icount-1)
                IpointerN(4) = Ipointer(3,icount-1)
                IpointerN(5) = Ipointer(0,icount-1)
                IpointerN(6) = Ipointer(4,icount-1)
                IpointerN(7) = Ipointer(7,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(2,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(6,icount-1)
                IpointerN(1) = Ipointer(2,icount-1)
                IpointerN(2) = Ipointer(1,icount-1)
                IpointerN(3) = Ipointer(5,icount-1)
                IpointerN(4) = Ipointer(7,icount-1)
                IpointerN(5) = Ipointer(3,icount-1)
                IpointerN(6) = Ipointer(0,icount-1)
                IpointerN(7) = Ipointer(4,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(5,icount-1)-1)-xco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(1,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(0,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(6,icount-1)
                IpointerN(1) = Ipointer(5,icount-1)
                IpointerN(2) = Ipointer(1,icount-1)
                IpointerN(3) = Ipointer(2,icount-1)
                IpointerN(4) = Ipointer(7,icount-1)
                IpointerN(5) = Ipointer(4,icount-1)
                IpointerN(6) = Ipointer(0,icount-1)
                IpointerN(7) = Ipointer(3,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(6,icount-1)-1)-xco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(4,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(5,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(4,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(2,icount-1)
                IpointerN(1) = Ipointer(6,icount-1)
                IpointerN(2) = Ipointer(5,icount-1)
                IpointerN(3) = Ipointer(1,icount-1)
                IpointerN(4) = Ipointer(3,icount-1)
                IpointerN(5) = Ipointer(7,icount-1)
                IpointerN(6) = Ipointer(4,icount-1)
                IpointerN(7) = Ipointer(0,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            else if ((xco(Ipointer(6,icount-1)-1)-xco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(7,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(3,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(6,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(1,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(0,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(0,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(3,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(6,icount-1)-1)>0 ) )     then
                IpointerN(0) = Ipointer(5,icount-1)
                IpointerN(1) = Ipointer(6,icount-1)
                IpointerN(2) = Ipointer(2,icount-1)
                IpointerN(3) = Ipointer(1,icount-1)
                IpointerN(4) = Ipointer(4,icount-1)
                IpointerN(5) = Ipointer(7,icount-1)
                IpointerN(6) = Ipointer(3,icount-1)
                IpointerN(7) = Ipointer(0,icount-1)
                Ipointer(0:7,icount-1) = IpointerN(0:7)
            endif
            ! checking
            if (.not.((xco(Ipointer(1,icount-1)-1)-xco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(2,icount-1)-1)-xco(Ipointer(3,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(5,icount-1)-1)-xco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (xco(Ipointer(6,icount-1)-1)-xco(Ipointer(7,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(3,icount-1)-1)-yco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(2,icount-1)-1)-yco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(6,icount-1)-1)-yco(Ipointer(5,icount-1)-1)>0 ) .and. &
                (yco(Ipointer(7,icount-1)-1)-yco(Ipointer(4,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(4,icount-1)-1)-zco(Ipointer(0,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(5,icount-1)-1)-zco(Ipointer(1,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(6,icount-1)-1)-zco(Ipointer(2,icount-1)-1)>0 ) .and. &
                (zco(Ipointer(7,icount-1)-1)-zco(Ipointer(3,icount-1)-1)>0 ) ) ) then
                write(*,*) "!!! Problem with element: ",icount
            endif
        end do
        !- Fortran -> C notation
        do icount = 0, nelem-1
            Ipointer(:,icount) = Ipointer(:,icount)-1
        end do
    end subroutine transfo_8

    !-------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------
    subroutine create_model_ondafly(n_mat_tot,mattab,pml_bool,pml_top,pml_bottom,   &
        xmin,xmax,step_x,ymin,ymax,step_y,zmax,mesh_type,nlay,thick)
        ! reads 'mat.dat' with info related to construction of a simple model on
        ! the fly 
        implicit none
        integer, intent(out)  :: n_mat_tot,pml_bool,pml_top,pml_bottom,mesh_type
        character, allocatable, dimension(:),intent(out) :: mattab
        real, intent(out)  :: step_x,step_y,xmin,xmax,ymin,ymax,zmax
        integer, dimension(:), allocatable, intent(out) :: nlay
        real, dimension(:), allocatable, intent(out) :: thick
        integer    :: i_err,n_mat,i,n_lay,ngllPML
        real       :: tr
        real, dimension(:), allocatable :: rho,pvel,svel,qp,qs
        integer, dimension(:), allocatable :: ngll
        character, allocatable, dimension(:):: matarray
        logical  :: strat_bool

        write(*,*)
        write(*,*) "-------------------------------------------------"
        write(*,*) "  -->  Mesh to be created on the fly: input data."
        write(*,*) "-------------------------------------------------"
        write(*,*)
        write(*,*) "  --> Construction of the model:"
        write(*,*) "        files to be found in the directory where 'mesher' is launched:" 
        write(*,*) "      1- general file 'mat.dat'" 
        write(*,*) "      2- the material file 'mater.in'"
        write(*,*)


        open(10,file="mater.in",action="read",status="old",iostat=i_err)
        if(i_err > 0) stop "File 'mater.in' not found in the working directory"
        read(10,*) n_mat
        allocate(rho(0:n_mat-1),pvel(0:n_mat-1),svel(0:n_mat-1),qp(0:n_mat-1),  &
            qs(0:n_mat-1),ngll(0:n_mat-1),matarray(0:n_mat-1))
        do i = 0,n_mat-1
            read(10,*) matarray(i),pvel(i),svel(i),rho(i),ngll(i),tr,tr,tr,qp(i),qs(i)
        end do
        close(10)


        open(10,file="mat.dat",action="read",status="old",iostat=i_err)
        if(i_err > 0) stop "File 'mat.dat' not found in the working directory"
        read(10,*) xmin ; read(10,*) xmax ; read(10,*) step_x
        read(10,*) ymin ; read(10,*) ymax ; read(10,*) step_y
        read(10,*) zmax
        read(10,*) n_lay  !- number of materials
        if(n_lay /= n_mat) stop "  --> Pb in construction of the model: number of materials diff. betw. mater.in and mat.dat"
        allocate(thick(0:n_lay-1),nlay(0:n_lay-1))
        do i = 0,n_lay-1
            read(10,*) thick(i),nlay(i)
        end do
        n_mat_tot = n_mat
        write(*,*) "  --> Number of non-PML materials:",n_mat
        strat_bool = .false.
        if(n_lay > 1)then
            print*,"  --> Stratified medium."
            strat_bool = .true.
        end if
        !  PMLs
        read(10,*) pml_bool
        if(pml_bool /= 0 .and. pml_bool /= 1) stop "In mesh2spec: PML or not?"
        if(pml_bool == 1)then   ! PMLs added
            write(*,*) "  --> PMLs added."
            n_mat_tot = n_mat_tot + 8*n_lay
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
            read(10,*) ngllPML
        else
            read(10,*)
            read(10,*)
            write(*,*)
            write(*,*) "  --> Run without absorbing boundaries: you may observe some reflected phases."
            write(*,*) "       Ok? (type enter if..)"
            read*
        end if
        ! number of control nodes
        read(10,*) mesh_type
        close(10)

        ! material file created
        allocate(mattab(0:n_mat_tot-1))

        call mat_table_construct_ondafly(n_mat,n_mat_tot,matarray,mattab,pml_bool,pml_top,pml_bottom)

        if(pml_bool == 1)then
            call write_mater_ondafly(n_mat,mattab,pvel,svel,rho,ngll,qp,qs,pml_bool,   &
                pml_top,pml_bottom,ngllPML)
        else
            call write_mater_ondafly(n_mat,mattab,pvel,svel,rho,ngll,qp,qs)
        end if

        deallocate(matarray,rho,pvel,svel,qp,qs,ngll)

    end subroutine create_model_ondafly
    !-------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------
    subroutine mat_table_construct(mattab)
        !- obtention of the material table array: important for the fluid/solid interfaces in mesh2spec
        character, dimension(:), allocatable, intent(inout)  :: mattab
        integer                 :: i,nblock

        open(10,file="mater.in",action="read",status="old")
        read(10,*) nblock
        allocate(mattab(0:nblock-1))
        do i = 0,size(mattab)-1
            read(10,"(a1)") mattab(i)
        end do
        close(10)

    end subroutine mat_table_construct
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------
    !-------------------------------------------------------------------------------------------------
    subroutine mat_table_construct_unv(mattab)
        !- obtention of the material table array: important for the fluid/solid interfaces in mesh2spec
        character, dimension(:), allocatable, intent(inout)  :: mattab
        integer                 :: i,nblock

        open(10,file="material.input",action="read",status="old")
        read(10,*) nblock
        allocate(mattab(0:nblock-1))
        do i = 0,size(mattab)-1
            read(10,"(a1)") mattab(i)
        end do
        close(10)

    end subroutine mat_table_construct_unv
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    subroutine mat_table_construct_ondafly(nmat,nmat_tot,matarray,mattab,pml_bool,pml_t,pml_b)
        !- obtention of the material table array: important for the fluid/solid interfaces in mesh2spec
        !     (case where mesh built on the fly)
        !- construction of the "material.input" file directly used in SEM3D
        character, dimension(0:), intent(in)   :: matarray
        character, dimension(0:), intent(out)  :: mattab
        integer, intent(in)     :: nmat,nmat_tot,pml_bool,pml_t,pml_b
        integer                 :: i,icount

        icount = nmat

        do i = 0,nmat-1
            mattab(i) = matarray(i)
        end do
        if(size(mattab) == nmat) return
        do i = 0,nmat-1
            if(mattab(i) == 'F') mattab(icount:icount+7) = 'L' 
            if(mattab(i) == 'S' .or. mattab(i) == 'R') mattab(icount:icount+7) = 'P'
            icount = icount + 8
        end do

        if(pml_b == 1)then
            if(mattab(nmat-1) == 'F') mattab(icount:icount+8) = 'L' 
            if(mattab(nmat-1) == 'S' .or. mattab(nmat-1) == 'R') mattab(icount:icount+8) = 'P'
            icount = icount + 9
        end if

        if(pml_t == 1)then
            if(mattab(0) == 'F') mattab(icount:icount+8) = 'L' 
            if(mattab(0) == 'S' .or. mattab(0) == 'R') mattab(icount:icount+8) = 'P'
            icount = icount + 9
        end if


    end subroutine mat_table_construct_ondafly
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------

    subroutine init_ondafly(xmin,xmax,ymin,ymax,zmax,xminref,xmaxref,yminref,ymaxref,zminref,zmaxref,    &
        step_x,step_y,mesh_type,nnods,npts,nelem,pml_b,pml_t,pml_bott,nlay,thick,z_dep)

        implicit none
        integer, intent(in)  :: pml_b,pml_t,pml_bott
        integer, dimension(0:), intent(in)  :: nlay
        real, dimension(0:), intent(in)  :: thick
        real, intent(inout)     :: xmin,xmax,ymin,ymax
        real, intent(in)     :: zmax,step_x,step_y
        real, intent(out)    :: xminref,xmaxref,yminref,ymaxref,zminref,zmaxref
        real, dimension(:), allocatable, intent(out)    ::  z_dep
        integer, intent(in)  :: mesh_type
        integer, intent(out) :: npts, nelem, nnods
        integer              :: i,j,n_z, nelemx,nelemy,nelemz,npx,npy,npz
        real   ::  step_max,step_z

        write(*,*) "*****************************************"
        write(*,*) "  --> Xmin: ", xmin
        write(*,*) "  --> Xmax: ",  xmax
        write(*,*) "  --> Ymin: ",  ymin
        write(*,*) "  --> Ymax: ",  ymax
        write(*,*) "  --> Zmax: ",  zmax
        write(*,*) "  --> Step in x: ", step_x
        write(*,*) "  --> Step in y: ", step_y
        write(*,*) "  --> 8 or 27 nodes (1 or 2)?: ", mesh_type

        if(mesh_type /= 1 .and. mesh_type /= 2) &
            stop "In init_ondafly: we said 1 or 2.."
        if(xmin > xmax .or. ymin > ymax .or. step_x <= 0d0 .or. step_y <= 0d0) &
            stop "Pb in init_ondafly at the outset."

        !- number of nodes
        nnods = merge(8,27,mesh_type == 1)


        !- construction of the z-dependence
        step_max = max(step_x,step_y)
        n_z = 1
        do i = 0,size(nlay)-1
            n_z = n_z+nlay(i)
        end do
        if(pml_t == 1) n_z = n_z+1
        if(pml_bott == 1) n_z = n_z+1

        allocate(z_dep(0:n_z-1))
        if(pml_t == 1)then
            z_dep(n_z-1) = zmax+step_max
            n_z = n_z-1
        end if
        z_dep(n_z-1) = zmax ; n_z = n_z-1 
        do i = 0,size(nlay)-1
            step_z = thick(i)/nlay(i)
            do j = 1,nlay(i)
                z_dep(n_z-1) = z_dep(n_z)-step_z
                n_z = n_z-1 
            end do
        end do
        if(pml_bott == 1)then
            z_dep(0) = z_dep(1)-step_max
        end if


        !- reference boundaries
        xminref = xmin ; xmaxref = xmax
        yminref = ymin ; ymaxref = ymax
        zmaxref = zmax
        if(pml_bott == 1)then
            zminref = z_dep(1)
        else
            zminref = z_dep(0)
        end if

        !- changes in lengths if PMLs added:
        if(pml_b == 1)then
            xmin = xmin-step_x ; xmax = xmax+step_x
            ymin = ymin-step_y ; ymax = ymax+step_y
        end if

        !- number of elements
        nelemx = ceiling((xmax-xmin)/step_x)
        nelemy = ceiling((ymax-ymin)/step_y)
        nelemz = size(z_dep)-1
        nelem = nelemx*nelemy*nelemz
        !- number of points
        npx = merge(nelemx+1,2*nelemx+1,mesh_type==1)
        npy = merge(nelemy+1,2*nelemy+1,mesh_type==1)
        npz = merge(nelemz+1,2*nelemz+1,mesh_type==1)
        npts = npx*npy*npz
        return

    end subroutine init_ondafly
    !----------------------------------
    subroutine mesh_on_the_fly(xmin,xmax,ymin,ymax,zmin,zmax,xstep,ystep,z_dep,meshtype,xp,yp,zp,Ipoint)

        implicit none
        real, intent(in)     :: xmin,xmax,ymin,ymax,zmin,zmax,xstep,ystep
        real, dimension(0:), intent(in)   :: z_dep
        integer, intent(in)  :: meshtype
        real, intent(out)    :: xp(0:),yp(0:),zp(0:)
        integer, intent(out) :: Ipoint(0:,0:)
        integer              :: nelemx,nelemy,nelemz,nx,ny,nz, indelem, aux_pt

        !- number of elements
        nelemx = ceiling((xmax-xmin)/xstep)
        nelemy = ceiling((ymax-ymin)/ystep)
        nelemz = size(z_dep)-1

        indelem = 0

        select case(meshtype)
        case(1)      ! 8 nodes
            write(*,*) "  --> 8 nodes per element."
            do nz = 0,nelemz-1
                do ny = 0,nelemy-1
                    do nx = 0, nelemx-1
                        Ipoint(0,indelem) = nx+ny*(nelemx+1)+nz*(nelemx+1)*(nelemy+1)
                        Ipoint(1,indelem) = nx+ny*(nelemx+1)+nz*(nelemx+1)*(nelemy+1)+1
                        Ipoint(2,indelem) = nx+(ny+1)*(nelemx+1)+nz*(nelemx+1)*(nelemy+1)+1
                        Ipoint(3,indelem) = nx+(ny+1)*(nelemx+1)+nz*(nelemx+1)*(nelemy+1)
                        Ipoint(4,indelem) = nx+ny*(nelemx+1)+(nz+1)*(nelemx+1)*(nelemy+1)
                        Ipoint(5,indelem) = nx+ny*(nelemx+1)+(nz+1)*(nelemx+1)*(nelemy+1)+1
                        Ipoint(6,indelem) = nx+(ny+1)*(nelemx+1)+(nz+1)*(nelemx+1)*(nelemy+1)+1
                        Ipoint(7,indelem) = nx+(ny+1)*(nelemx+1)+(nz+1)*(nelemx+1)*(nelemy+1)
                        indelem = indelem+1
                    end do
                end do
            end do

            !- co-ordinates
            aux_pt = 0
            do nz = 0,nelemz
                do ny = 0, nelemy
                    do nx = 0, nelemx
                        xp(aux_pt) = xmin+nx*xstep
                        yp(aux_pt) = ymin+ny*ystep
                        zp(aux_pt) = z_dep(nz)
                        aux_pt = aux_pt+1
                    end do
                end do
            end do

        case(2)      ! 27 nodes
            write(*,*) "  --> 27 nodes per element."
            !-
            do nz = 0,nelemz-1
                do ny = 0,nelemy-1
                    do nx = 0, nelemx-1
                        Ipoint(0,indelem) = 2*nx+2*ny*(2*nelemx+1)+2*nz*(2*nelemx+1)*(2*nelemy+1)
                        Ipoint(1,indelem) = 2*nx+2*ny*(2*nelemx+1)+2*nz*(2*nelemx+1)*(2*nelemy+1)+2
                        Ipoint(2,indelem) = 2*nx+2*(ny+1)*(2*nelemx+1)+2*nz*(2*nelemx+1)*(2*nelemy+1)+2
                        Ipoint(3,indelem) = 2*nx+2*(ny+1)*(2*nelemx+1)+2*nz*(2*nelemx+1)*(2*nelemy+1)
                        Ipoint(4,indelem) = 2*nx+2*ny*(2*nelemx+1)+2*(nz+1)*(2*nelemx+1)*(2*nelemy+1)
                        Ipoint(5,indelem) = 2*nx+2*ny*(2*nelemx+1)+2*(nz+1)*(2*nelemx+1)*(2*nelemy+1)+2
                        Ipoint(6,indelem) = 2*nx+2*(ny+1)*(2*nelemx+1)+2*(nz+1)*(2*nelemx+1)*(2*nelemy+1)+2
                        Ipoint(7,indelem) = 2*nx+2*(ny+1)*(2*nelemx+1)+2*(nz+1)*(2*nelemx+1)*(2*nelemy+1)

                        Ipoint(8,indelem) = 2*nx+2*ny*(2*nelemx+1)+2*nz*(2*nelemx+1)*(2*nelemy+1)+1
                        Ipoint(9,indelem) = 2*nx+(2*ny+1)*(2*nelemx+1)+2*nz*(2*nelemx+1)*(2*nelemy+1)+2
                        Ipoint(10,indelem) = 2*nx+2*(ny+1)*(2*nelemx+1)+2*nz*(2*nelemx+1)*(2*nelemy+1)+1
                        Ipoint(11,indelem) = 2*nx+(2*ny+1)*(2*nelemx+1)+2*nz*(2*nelemx+1)*(2*nelemy+1)
                        Ipoint(12,indelem) = 2*nx+2*ny*(2*nelemx+1)+(2*nz+1)*(2*nelemx+1)*(2*nelemy+1)
                        Ipoint(13,indelem) = 2*nx+2*ny*(2*nelemx+1)+(2*nz+1)*(2*nelemx+1)*(2*nelemy+1)+2
                        Ipoint(14,indelem) = 2*nx+2*(ny+1)*(2*nelemx+1)+(2*nz+1)*(2*nelemx+1)*(2*nelemy+1)+2
                        Ipoint(15,indelem) = 2*nx+2*(ny+1)*(2*nelemx+1)+(2*nz+1)*(2*nelemx+1)*(2*nelemy+1)
                        Ipoint(16,indelem) = 2*nx+2*ny*(2*nelemx+1)+2*(nz+1)*(2*nelemx+1)*(2*nelemy+1)+1
                        Ipoint(17,indelem) = 2*nx+(2*ny+1)*(2*nelemx+1)+2*(nz+1)*(2*nelemx+1)*(2*nelemy+1)+2
                        Ipoint(18,indelem) = 2*nx+2*(ny+1)*(2*nelemx+1)+2*(nz+1)*(2*nelemx+1)*(2*nelemy+1)+1
                        Ipoint(19,indelem) = 2*nx+(2*ny+1)*(2*nelemx+1)+2*(nz+1)*(2*nelemx+1)*(2*nelemy+1)
                        Ipoint(20,indelem) = 2*nx+(2*ny+1)*(2*nelemx+1)+2*nz*(2*nelemx+1)*(2*nelemy+1)+1
                        Ipoint(21,indelem) = 2*nx+2*ny*(2*nelemx+1)+(2*nz+1)*(2*nelemx+1)*(2*nelemy+1)+1
                        Ipoint(22,indelem) = 2*nx+(2*ny+1)*(2*nelemx+1)+(2*nz+1)*(2*nelemx+1)*(2*nelemy+1)+2
                        Ipoint(23,indelem) = 2*nx+2*(ny+1)*(2*nelemx+1)+(2*nz+1)*(2*nelemx+1)*(2*nelemy+1)+1
                        Ipoint(24,indelem) = 2*nx+(2*ny+1)*(2*nelemx+1)+(2*nz+1)*(2*nelemx+1)*(2*nelemy+1)
                        Ipoint(25,indelem) = 2*nx+(2*ny+1)*(2*nelemx+1)+2*(nz+1)*(2*nelemx+1)*(2*nelemy+1)+1
                        Ipoint(26,indelem) = 2*nx+(2*ny+1)*(2*nelemx+1)+(2*nz+1)*(2*nelemx+1)*(2*nelemy+1)+1
                        indelem = indelem+1
                    end do
                end do
            end do
            !- co-ordinates
            aux_pt = 0
            do nz = 0,2*nelemz
                do ny = 0, 2*nelemy
                    do nx = 0, 2*nelemx
                        xp(aux_pt) = xmin+nx*xstep/2d0
                        yp(aux_pt) = ymin+ny*ystep/2d0
                        if(mod(nz,2) == 0)then
                            zp(aux_pt) = z_dep(nz/2)
                        else
                            zp(aux_pt) = (z_dep(nz/2)+z_dep(nz/2+1))/2d0
                        end if
                        aux_pt = aux_pt+1
                    end do
                end do
            end do
        end select

        return

    end subroutine mesh_on_the_fly
    !---------------------
    subroutine nature_elem(Ipoint,xp,yp,zp,nmat,mat,xminref,xmaxref,       &
        yminref,ymaxref,zminref,zmaxref,pml_bool,pml_t,pml_b,thick)
        !- when on the fly construction: allows to determine the PML layers
        integer, intent(in)   :: Ipoint(0:,0:),nmat,pml_bool,pml_t,pml_b
        real, intent(in)      :: xp(0:),yp(0:),zp(0:), xminref, xmaxref,   &
            yminref,ymaxref,zminref,zmaxref
        real, dimension(0:),intent(in)   :: thick
        integer, intent(inout)  :: mat(0:)
        integer               :: i,j,n,nelem,ind_pml
        real                  :: coord(0:7,0:2), bary(0:2), z(0:size(thick))

        nelem = size(Ipoint,2)

        ind_pml = 0
        if(pml_t == 1) ind_pml = 9

        z(0) = zmaxref
        do j = 1,size(thick)
            z(j) = z(j-1)-thick(j-1)
        end do

        do n = 0,nelem-1
            do i = 0,7
                coord(i,0) = xp(Ipoint(i,n))
                coord(i,1) = yp(Ipoint(i,n))
                coord(i,2) = zp(Ipoint(i,n))
            end do
            call barycentre(coord,bary)

            ! eventual bottom and top PMLs 
            if(bary(2) > zmaxref) then   ! top PMLs
                if(bary(1) < yminref)then
                    if(bary(0) < xminref)then
                        Mat(n) = nmat-9
                    else if(bary(0) > xmaxref)then
                        Mat(n) = nmat-8
                    else
                        Mat(n) = nmat-5
                    end if
                else if(bary(1) > ymaxref)then
                    if(bary(0) < xminref)then
                        Mat(n) = nmat-6
                    else if(bary(0) > xmaxref)then
                        Mat(n) = nmat-7
                    else
                        Mat(n) = nmat-3
                    end if
                else
                    if(bary(0) < xminref)then
                        Mat(n) = nmat-2
                    else if(bary(0) > xmaxref)then
                        Mat(n) = nmat-4
                    else
                        Mat(n) = nmat-1
                    end if
                end if
            end if
            if(bary(2) < zminref) then   ! bottom PMLs
                if(bary(1) < yminref)then
                    if(bary(0) < xminref)then
                        Mat(n) = nmat-9-ind_pml
                    else if(bary(0) > xmaxref)then
                        Mat(n) = nmat-8-ind_pml
                    else
                        Mat(n) = nmat-5-ind_pml
                    end if
                else if(bary(1) > ymaxref)then
                    if(bary(0) < xminref)then
                        Mat(n) = nmat-6-ind_pml
                    else if(bary(0) > xmaxref)then
                        Mat(n) = nmat-7-ind_pml
                    else
                        Mat(n) = nmat-3-ind_pml
                    end if
                else
                    if(bary(0) < xminref)then
                        Mat(n) = nmat-2-ind_pml
                    else if(bary(0) > xmaxref)then
                        Mat(n) = nmat-4-ind_pml
                    else
                        Mat(n) = nmat-1-ind_pml
                    end if
                end if
            end if

            if(bary(2) < zmaxref .and. bary(2) > zminref)then
                do j = 0,size(z)-1
                    if(bary(2) < z(j) .and. bary(2) > z(j+1)) exit 
                end do

                if(bary(1) < yminref)then
                    if(bary(0) < xminref)then
                        Mat(n) = size(thick)+8*j
                    else if(bary(0) > xmaxref)then
                        Mat(n) = size(thick)+8*j+1
                    else
                        Mat(n) = size(thick)+8*j+4
                    end if
                else if(bary(1) > ymaxref)then
                    if(bary(0) < xminref)then
                        Mat(n) = size(thick)+8*j+3
                    else if(bary(0) > xmaxref)then
                        Mat(n) = size(thick)+8*j+2
                    else
                        Mat(n) = size(thick)+8*j+6
                    end if
                else
                    if(bary(0) < xminref)then
                        Mat(n) = size(thick)+8*j+7
                    else if(bary(0) > xmaxref)then
                        Mat(n) = size(thick)+8*j+5
                    else
                        Mat(n) = j
                    end if
                end if
            end if

        end do

    end subroutine nature_elem
    !--------------------------------------------------------
    !--------------------------------------------------------

    subroutine barycentre(coord,bary)
        real, intent(in), dimension(0:7,0:2)   :: coord
        real, intent(out), dimension(0:2)  :: bary
        integer :: i

        bary(0:) = 0d0
        do i = 0,7
            bary(0:) = bary(0:)+coord(i,0:)
        end do
        bary = bary/8

    end subroutine barycentre
    !--------------------------------------------------------------------
    !--------------------------------------------------------------------
    subroutine write_mater_ondafly(nmat,mattab,vp,vs,rho,ngll,qp,qs,pml_bool,pml_t,pml_b,ngll_PML)
        !- construction of the "material.input" file directly used in SEM3D
        integer, intent(in)   :: nmat
        character, dimension(0:), intent(in)   :: mattab
        real, dimension(0:nmat-1), intent(in)  :: vp,vs,rho,qp,qs
        integer, dimension(0:nmat-1), intent(in)  :: ngll
        integer, intent(in),optional     :: pml_bool,pml_t,pml_b,ngll_PML
        integer                 :: i,icount
        real  :: tr
        character(len=*), parameter  :: FMT1="(a1,2x,f9.3,2x,f9.3,2x,f9.3,2x,3(i2,2x),f8.5,2x,f8.3,2x,f8.3)"
        character(len=*), parameter  :: FMT2="(l1,1x,i2,1x,f5.2,1x,6(l1,1x),f5.2,i3)"
        real, parameter  :: a = 10d0, k =0d0
        integer, parameter  :: n = 2
        logical, parameter   :: VRAI = .true. , FAUX = .false.

        integer :: nblock,j,jj,count_matR
        integer :: seed
        character(len=10) :: corr, marga,margb,margc
        real(kind=8) :: lx,ly,lz,sigma2a,sigma2b,sigma2c

        icount = 0
        tr = 0.0001d0

        open(10,file="material.input",action="write",status="replace")
        write(10,*) size(mattab)
        !  physical layers
        do i = 0,nmat-1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll(i),ngll(i),   &
                ngll(i),tr,qp(i),qs(i)
            icount = icount+1
        end do

        ! lateral PMLs
        if(present(pml_bool) .and. pml_bool == 1)then
            do i = 0,nmat-1
                write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                    ngll(i),tr,qp(i),qs(i)
                icount = icount+1
                write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                    ngll(i),tr,qp(i),qs(i)
                icount = icount+1
                write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                    ngll(i),tr,qp(i),qs(i)
                icount = icount+1
                write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                    ngll(i),tr,qp(i),qs(i)
                icount = icount+1
                write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll(i),ngll_PML,   &
                    ngll(i),tr,qp(i),qs(i)
                icount = icount+1
                write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll(i),   &
                    ngll(i),tr,qp(i),qs(i)
                icount = icount+1
                write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll(i),ngll_PML,   &
                    ngll(i),tr,qp(i),qs(i)
                icount = icount+1
                write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll(i),   &
                    ngll(i),tr,qp(i),qs(i)
                icount = icount+1

            end do
        end if
        ! bottom PMLs
        if(present(pml_b) .and. pml_b == 1)then
            i = nmat-1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll(i),ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll(i),   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll(i),ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll(i),   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll(i),ngll(i),   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
        end if
        ! top PMLs
        if(present(pml_t) .and. pml_t == 1)then
            i = 0
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll(i),ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll(i),   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll(i),ngll_PML,   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll_PML,ngll(i),   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
            write(10,FMT=FMT1) mattab(icount),vp(i),vs(i),rho(i),ngll(i),ngll(i),   &
                ngll_PML,tr,qp(i),qs(i)
            icount = icount+1
        end if

        !! lines dedicated to PMLs
        write(10,*) ; write(10,*)
        ! lateral PMLs
        if(present(pml_bool) .and. pml_bool == 1)then
            do i = 0,nmat-1
                write(10,FMT=FMT2) FAUX,n,a,VRAI,VRAI,VRAI,VRAI,FAUX,FAUX,k,i
                write(10,FMT=FMT2) FAUX,n,a,VRAI,FAUX,VRAI,VRAI,FAUX,FAUX,k,i
                write(10,FMT=FMT2) FAUX,n,a,VRAI,FAUX,VRAI,FAUX,FAUX,FAUX,k,i
                write(10,FMT=FMT2) FAUX,n,a,VRAI,VRAI,VRAI,FAUX,FAUX,FAUX,k,i
                write(10,FMT=FMT2) FAUX,n,a,FAUX,VRAI,VRAI,VRAI,FAUX,FAUX,k,i
                write(10,FMT=FMT2) FAUX,n,a,VRAI,FAUX,FAUX,VRAI,FAUX,FAUX,k,i
                write(10,FMT=FMT2) FAUX,n,a,FAUX,VRAI,VRAI,FAUX,FAUX,FAUX,k,i
                write(10,FMT=FMT2) FAUX,n,a,VRAI,VRAI,FAUX,VRAI,FAUX,FAUX,k,i
            end do
        end if
        ! bottom PMLs
        if(present(pml_b) .and. pml_b == 1)then
            write(10,FMT=FMT2) FAUX,n,a,VRAI,VRAI,VRAI,VRAI,VRAI,VRAI,k,nmat-1
            write(10,FMT=FMT2) FAUX,n,a,VRAI,FAUX,VRAI,VRAI,VRAI,VRAI,k,nmat-1
            write(10,FMT=FMT2) FAUX,n,a,VRAI,FAUX,VRAI,FAUX,VRAI,VRAI,k,nmat-1
            write(10,FMT=FMT2) FAUX,n,a,VRAI,VRAI,VRAI,FAUX,VRAI,VRAI,k,nmat-1
            write(10,FMT=FMT2) FAUX,n,a,FAUX,VRAI,VRAI,VRAI,VRAI,VRAI,k,nmat-1
            write(10,FMT=FMT2) FAUX,n,a,VRAI,FAUX,FAUX,VRAI,VRAI,VRAI,k,nmat-1
            write(10,FMT=FMT2) FAUX,n,a,FAUX,VRAI,VRAI,FAUX,VRAI,VRAI,k,nmat-1
            write(10,FMT=FMT2) FAUX,n,a,VRAI,VRAI,FAUX,VRAI,VRAI,VRAI,k,nmat-1
            write(10,FMT=FMT2) FAUX,n,a,FAUX,VRAI,FAUX,VRAI,VRAI,VRAI,k,nmat-1
        end if
        ! top PMLs
        if(present(pml_t) .and. pml_t == 1)then
            write(10,FMT=FMT2) FAUX,n,a,VRAI,VRAI,VRAI,VRAI,VRAI,FAUX,k,0
            write(10,FMT=FMT2) FAUX,n,a,VRAI,FAUX,VRAI,VRAI,VRAI,FAUX,k,0
            write(10,FMT=FMT2) FAUX,n,a,VRAI,FAUX,VRAI,FAUX,VRAI,FAUX,k,0
            write(10,FMT=FMT2) FAUX,n,a,VRAI,VRAI,VRAI,FAUX,VRAI,FAUX,k,0
            write(10,FMT=FMT2) FAUX,n,a,FAUX,VRAI,VRAI,VRAI,VRAI,FAUX,k,0
            write(10,FMT=FMT2) FAUX,n,a,VRAI,FAUX,FAUX,VRAI,VRAI,FAUX,k,0
            write(10,FMT=FMT2) FAUX,n,a,FAUX,VRAI,VRAI,FAUX,VRAI,FAUX,k,0
            write(10,FMT=FMT2) FAUX,n,a,VRAI,VRAI,FAUX,VRAI,VRAI,FAUX,k,0
            write(10,FMT=FMT2) FAUX,n,a,FAUX,VRAI,FAUX,VRAI,VRAI,FAUX,k,0
        end if

        !! Lines dedicated to Random Field
        write(10,*) ; write(10,*)
        count_matR = 0
        do j = 0,nmat-1
            if(mattab(j)=='R') then
                open(20,file="mater.in",action="read",status="old")
                read(20,*) nblock
                do jj = 0,(nblock+1+count_matR)
                    read(20,*)
                end do
                read(20,*)  corr, lx, ly, lz, marga ,sigma2a, margb ,sigma2b, margc ,sigma2c,seed
                write(10,*) corr, lx, ly, lz, marga ,sigma2a, margb ,sigma2b, margc ,sigma2c,seed
                count_matR = count_matR + 1
                close(20)
            end if
        end do

        close(10)



    end subroutine write_mater_ondafly

    !--------------------------------------------------------
    !--------------------------------------------------------
    integer function vertices2edge(vert)

        implicit none
        integer, intent(in)    :: vert(0:1)
        integer                :: corner(0:1), edge

        corner = vert
        call sort(corner,2)
        if(corner(0)==0) then
            if(corner(1)==1) edge = 0
            if(corner(1)==3) edge = 3
            if(corner(1)==4) edge = 6
        else if(corner(0)==1) then
            if(corner(1)==2) edge = 1
            if(corner(1)==5) edge = 4
        else if(corner(0)==2) then
            if(corner(1)==3) edge = 2
            if(corner(1)==6) edge = 7
        else if(corner(0)==3) then
            edge = 10
        else if(corner(0)==4) then
            if(corner(1)==5) edge = 5
            if(corner(1)==7) edge = 11
        else if(corner(0)==5) then
            edge = 8
        else if(corner(0)==6) then
            edge = 9
        endif

        vertices2edge = edge

    end function vertices2edge
    !------------------------------------------------
    integer function edge2vertex(edge)
        ! returns the index of the 1st of the 2 vertices for a given edge
        integer, intent(in)  :: edge

        if(edge == 0 .or. edge == 3 .or. edge == 6) edge2vertex = 0
        if(edge == 1 .or. edge == 4) edge2vertex = 1
        if(edge == 7) edge2vertex = 2
        if(edge == 2 .or. edge == 10) edge2vertex = 3
        if(edge == 5 .or. edge == 11) edge2vertex = 4
        if(edge == 8) edge2vertex = 5
        if(edge == 9) edge2vertex = 7

    end function edge2vertex
    !------------------------------------------------
    subroutine edge_orientation(corner,length,ne,orient)

        implicit none

        integer, intent(IN) :: length, ne
        integer, dimension(0:length-1), intent(IN) :: corner
        integer, intent(INOUT) :: orient

        orient = 0
        if (ne==2 .or. ne==9) then
            if (corner(0) < corner(1)) orient = 1
        else
            if (corner(0) > corner(1)) orient = 1
        endif

        return
    end subroutine edge_orientation
    !-------------------------------------
    subroutine face2corner(Ipoint,nfa,corn)
        implicit none
        integer, intent(in)   :: Ipoint(0:7),nfa
        integer, intent(out)  :: corn(0:3)

        select case(nfa)
        case(0)
            corn(0) = Ipoint(0)
            corn(1) = Ipoint(1)
            corn(2) = Ipoint(2)
            corn(3) = Ipoint(3)
        case(1)
            corn(0) = Ipoint(0)
            corn(1) = Ipoint(1)
            corn(2) = Ipoint(5)
            corn(3) = Ipoint(4)
        case(2)
            corn(0) = Ipoint(1)
            corn(1) = Ipoint(2)
            corn(2) = Ipoint(6)
            corn(3) = Ipoint(5)
        case(3)
            corn(0) = Ipoint(3)
            corn(1) = Ipoint(2)
            corn(2) = Ipoint(6)
            corn(3) = Ipoint(7)
        case(4)
            corn(0) = Ipoint(0)
            corn(1) = Ipoint(3)
            corn(2) = Ipoint(7)
            corn(3) = Ipoint(4)
        case(5)
            corn(0) = Ipoint(4)
            corn(1) = Ipoint(5)
            corn(2) = Ipoint(6)
            corn(3) = Ipoint(7)
        end select

    end subroutine face2corner
    !---------------------------------------
    subroutine edge2corner(Ipoint,ne,e_corn)

        implicit none
        integer, intent(in)  :: Ipoint(0:),ne
        integer, intent(out) :: e_corn(0:1)

        select case (ne)
        case (0)
            e_corn(0) = Ipoint(0)
            e_corn(1) = Ipoint(1)
        case (1)
            e_corn(0) = Ipoint(1)
            e_corn(1) = Ipoint(2)
        case (2)
            e_corn(0) = Ipoint(3)
            e_corn(1) = Ipoint(2)
        case (3)
            e_corn(0) = Ipoint(0)
            e_corn(1) = Ipoint(3)
        case (4)
            e_corn(0) = Ipoint(1)
            e_corn(1) = Ipoint(5)
        case (5)
            e_corn(0) = Ipoint(4)
            e_corn(1) = Ipoint(5)
        case (6)
            e_corn(0) = Ipoint(0)
            e_corn(1) = Ipoint(4)
        case (7)
            e_corn(0) = Ipoint(2)
            e_corn(1) = Ipoint(6)
        case (8)
            e_corn(0) = Ipoint(5)
            e_corn(1) = Ipoint(6)
        case (9)
            e_corn(0) = Ipoint(7)
            e_corn(1) = Ipoint(6)
        case (10)
            e_corn(0) = Ipoint(3)
            e_corn(1) = Ipoint(7)
        case (11)
            e_corn(0) = Ipoint(4)
            e_corn(1) = Ipoint(7)
        end select

    end subroutine edge2corner
    !-------------------------------------------
    integer function neighb_face(neighb_corn)
        implicit none
        integer, dimension(0:3), intent(in)  :: neighb_corn


        if(neighb_corn(0)==0) then   ! So which face of the neighbor is it ?
            if(neighb_corn(3)==3) neighb_face = 0
            if(neighb_corn(3)==5) neighb_face = 1
            if(neighb_corn(3)==7) neighb_face = 4
        else if(neighb_corn(0)==1) then
            neighb_face = 2
        else if(neighb_corn(0)==2) then
            neighb_face = 3
        else if(neighb_corn(0)==4) then
            neighb_face = 5
        else
            stop "Coherency Pb between faces and nodes of an element"
        endif

    end function neighb_face
    !--------------------------------------------
    integer function neighb_edge(e_neighbor_corn)

        implicit none
        integer, intent(in), dimension(0:1)  :: e_neighbor_corn

        if (e_neighbor_corn(0)==0) then
            if (e_neighbor_corn(1)==1) neighb_edge = 0
            if (e_neighbor_corn(1)==3) neighb_edge = 3
            if (e_neighbor_corn(1)==4) neighb_edge = 6
        else if (e_neighbor_corn(0)==1) then
            if (e_neighbor_corn(1)==2) neighb_edge = 1
            if (e_neighbor_corn(1)==5) neighb_edge = 4
        else if (e_neighbor_corn(0)==2) then
            if (e_neighbor_corn(1)==3) neighb_edge = 2
            if (e_neighbor_corn(1)==6) neighb_edge = 7
        else if (e_neighbor_corn(0)==3) then
            neighb_edge = 10
        else if (e_neighbor_corn(0)==4) then
            if (e_neighbor_corn(1)==5) neighb_edge = 5
            if (e_neighbor_corn(1)==7) neighb_edge = 11
        else if (e_neighbor_corn(0)==5) then
            neighb_edge = 8
        else if (e_neighbor_corn(0)==6) then
            neighb_edge = 9
        endif

    end function neighb_edge
    !----------------------------------------
    subroutine caract_face(nf,corner_d,edge_d)

        implicit none

        integer, intent(IN) :: nf
        integer, dimension(0:3), intent(INOUT) :: corner_d,edge_d


        select case (nf)
        case (0)
            corner_d(0) = 0
            corner_d(1) = 1
            corner_d(2) = 2
            corner_d(3) = 3
            edge_d(0) = 0
            edge_d(1) = 1
            edge_d(2) = 2
            edge_d(3) = 3
        case (1)
            corner_d(0) = 0
            corner_d(1) = 1
            corner_d(2) = 5
            corner_d(3) = 4
            edge_d(0) = 0
            edge_d(1) = 4
            edge_d(2) = 5
            edge_d(3) = 6
        case (2)
            corner_d(0) = 1
            corner_d(1) = 2
            corner_d(2) = 6
            corner_d(3) = 5
            edge_d(0) = 1
            edge_d(1) = 7
            edge_d(2) = 8
            edge_d(3) = 4
        case (3)
            corner_d(0) = 3
            corner_d(1) = 2
            corner_d(2) = 6
            corner_d(3) = 7
            edge_d(0) = 2
            edge_d(1) = 7
            edge_d(2) = 9
            edge_d(3) = 10
        case (4)
            corner_d(0) = 0
            corner_d(1) = 3
            corner_d(2) = 7
            corner_d(3) = 4
            edge_d(0) = 3
            edge_d(1) = 10
            edge_d(2) = 11
            edge_d(3) = 6
        case (5)
            corner_d(0) = 4
            corner_d(1) = 5
            corner_d(2) = 6
            corner_d(3) = 7
            edge_d(0) = 5
            edge_d(1) = 8
            edge_d(2) = 9
            edge_d(3) = 11
        end select

        return

    end subroutine caract_face
    !--------------------------------------------------------
    subroutine face_orientation(corner,length,nf,orient)

        implicit none

        integer, intent(IN) :: length, nf
        integer, dimension(0:length-1), intent(IN) :: corner
        integer, intent(INOUT) :: orient

        integer :: n1,n2,n3

        n1 = corner(0)
        n2 = corner(1)
        n3 = corner(2)
        if (nf==0) then
            if ( n1==0 .and. n2==1 .and. n3==2 ) orient = 0
            if ( n1==1 .and. n2==2 .and. n3==3 ) orient = 6
            if ( n1==2 .and. n2==3 .and. n3==0 ) orient = 3
            if ( n1==3 .and. n2==0 .and. n3==1 ) orient = 5
            if ( n1==0 .and. n2==3 .and. n3==2 ) orient = 4
            if ( n1==3 .and. n2==2 .and. n3==1 ) orient = 2
            if ( n1==2 .and. n2==1 .and. n3==0 ) orient = 7
            if ( n1==1 .and. n2==0 .and. n3==3 ) orient = 1
        else if (nf==1) then
            if ( n1==0 .and. n2==1 .and. n3==5 ) orient = 0
            if ( n1==1 .and. n2==5 .and. n3==4 ) orient = 6
            if ( n1==5 .and. n2==4 .and. n3==0 ) orient = 3
            if ( n1==4 .and. n2==0 .and. n3==1 ) orient = 5
            if ( n1==0 .and. n2==4 .and. n3==5 ) orient = 4
            if ( n1==4 .and. n2==5 .and. n3==1 ) orient = 2
            if ( n1==5 .and. n2==1 .and. n3==0 ) orient = 7
            if ( n1==1 .and. n2==0 .and. n3==4 ) orient = 1
        else if (nf==2) then
            if ( n1==1 .and. n2==2 .and. n3==6 ) orient = 0
            if ( n1==2 .and. n2==6 .and. n3==5 ) orient = 6
            if ( n1==6 .and. n2==5 .and. n3==1 ) orient = 3
            if ( n1==5 .and. n2==1 .and. n3==2 ) orient = 5
            if ( n1==1 .and. n2==5 .and. n3==6 ) orient = 4
            if ( n1==5 .and. n2==6 .and. n3==2 ) orient = 2
            if ( n1==6 .and. n2==2 .and. n3==1 ) orient = 7
            if ( n1==2 .and. n2==1 .and. n3==5 ) orient = 1
        else if (nf==3) then
            if ( n1==3 .and. n2==2 .and. n3==6 ) orient = 0
            if ( n1==2 .and. n2==6 .and. n3==7 ) orient = 6
            if ( n1==6 .and. n2==7 .and. n3==3 ) orient = 3
            if ( n1==7 .and. n2==3 .and. n3==2 ) orient = 5
            if ( n1==3 .and. n2==7 .and. n3==6 ) orient = 4
            if ( n1==7 .and. n2==6 .and. n3==2 ) orient = 2
            if ( n1==6 .and. n2==2 .and. n3==3 ) orient = 7
            if ( n1==2 .and. n2==3 .and. n3==7 ) orient = 1
        else if (nf==4) then
            if ( n1==0 .and. n2==3 .and. n3==7 ) orient = 0
            if ( n1==3 .and. n2==7 .and. n3==4 ) orient = 6
            if ( n1==7 .and. n2==4 .and. n3==0 ) orient = 3
            if ( n1==4 .and. n2==0 .and. n3==3 ) orient = 5
            if ( n1==0 .and. n2==4 .and. n3==7 ) orient = 4
            if ( n1==4 .and. n2==7 .and. n3==3 ) orient = 2
            if ( n1==7 .and. n2==3 .and. n3==0 ) orient = 7
            if ( n1==3 .and. n2==0 .and. n3==4 ) orient = 1
        else if (nf==5) then
            if ( n1==4 .and. n2==5 .and. n3==6 ) orient = 0
            if ( n1==5 .and. n2==6 .and. n3==7 ) orient = 6
            if ( n1==6 .and. n2==7 .and. n3==4 ) orient = 3
            if ( n1==7 .and. n2==4 .and. n3==5 ) orient = 5
            if ( n1==4 .and. n2==7 .and. n3==6 ) orient = 4
            if ( n1==7 .and. n2==6 .and. n3==5 ) orient = 2
            if ( n1==6 .and. n2==5 .and. n3==4 ) orient = 7
            if ( n1==5 .and. n2==4 .and. n3==7 ) orient = 1
        endif

    end subroutine face_orientation
    !-----------------------------------------------------
    subroutine which_edge(corner,edge)

        implicit none

        integer, dimension(0:3), intent(IN) :: corner
        integer, dimension(0:3), intent(INOUT) :: edge

        integer :: i,j1,j2


        do i = 0,3
            j1 = corner(i)
            if (i==3) then
                j2 = corner(0)
            else
                j2 = corner(i+1)
            endif
            if ( (j1==0 .and. j2==1) .or. (j1==1 .and. j2==0) ) then
                edge(i) = 0
            else if ( (j1==1 .and. j2==2) .or. (j1==2 .and. j2==1) ) then
                edge(i) = 1
            else if ( (j1==2 .and. j2==3) .or. (j1==3 .and. j2==2) ) then
                edge(i) = 2
            else if ( (j1==0 .and. j2==3) .or. (j1==3 .and. j2==0) ) then
                edge(i) = 3
            else if ( (j1==1 .and. j2==5) .or. (j1==5 .and. j2==1) ) then
                edge(i) = 4
            else if ( (j1==4 .and. j2==5) .or. (j1==5 .and. j2==4) ) then
                edge(i) = 5
            else if ( (j1==0 .and. j2==4) .or. (j1==4 .and. j2==0) ) then
                edge(i) = 6
            else if ( (j1==2 .and. j2==6) .or. (j1==6 .and. j2==2) ) then
                edge(i) = 7
            else if ( (j1==5 .and. j2==6) .or. (j1==6 .and. j2==5) ) then
                edge(i) = 8
            else if ( (j1==6 .and. j2==7) .or. (j1==7 .and. j2==6) ) then
                edge(i) = 9
            else if ( (j1==3 .and. j2==7) .or. (j1==7 .and. j2==3) ) then
                edge(i) = 10
            else if ( (j1==4 .and. j2==7) .or. (j1==7 .and. j2==4) ) then
                edge(i) = 11
            endif
        enddo

    end subroutine which_edge
    !-----------------------------------------------------
    !--------------------------------------------------------

end module mesh_properties


!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
