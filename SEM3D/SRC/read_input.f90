!>
!!\file read_input.f90
!!\brief Contient la subroutine read_input().
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Assure la lecture des fichiers de donn�es en entr�e.
!!
!! \param type (domain), intent (INOUT) Tdomain
!<

module semconfig
    use semdatafiles
    use iso_c_binding
    integer, parameter :: NAME_MAX=200
    ! Ce type doit correspondre au type sem_config_t du module read_input.c **a l'ordre pres**
    type, bind(c) :: sem_config
       type(C_PTR)    :: run_name

       !! Integration
       integer(C_INT) :: accel_scheme
       integer(C_INT) :: veloc_scheme
       real(C_DOUBLE) :: sim_time
       real(C_DOUBLE) :: alpha
       real(C_DOUBLE) :: beta
       real(C_DOUBLE) :: gamma
       real(C_DOUBLE) :: courant

       !! Modele, maillage
       type(C_PTR)    :: mesh_file
       integer(C_INT) :: model
       integer(C_INT) :: anisotropy
       type(C_PTR)    :: mat_file
       integer(C_INT) :: nsources
       type(C_PTR)    :: source

       !! Capteurs
       integer(C_INT) :: save_traces
       integer(C_INT) :: traces_interval
       type(C_PTR)    :: station_file

       !! Snapshots
       integer(C_INT) :: save_snap
       real(C_DOUBLE) :: snap_interval

       !! Protection reprise
       integer(C_INT) :: prorep
       integer(C_INT) :: prorep_iter
       integer(C_INT) :: prorep_restart_iter

       integer(C_INT) :: verbose_level
       real(C_DOUBLE) :: mpml

       !! Amortissement
       integer(C_INT) :: nsolids
       real(C_DOUBLE), dimension(2) :: atn_band
       real(C_DOUBLE) :: atn_period

       !! Neumann
       integer(C_INT) :: neu_present
       integer(C_INT) :: neu_type
       integer(C_INT) :: neu_mat
       real(C_DOUBLE), dimension(3) :: neu_L
       real(C_DOUBLE), dimension(3) :: neu_C
       real(C_DOUBLE) :: neu_f0
    end type sem_config


    type, bind(c) :: sem_source
       type(C_PTR) :: next
       real(C_DOUBLE), dimension(3) :: coords;
       integer(C_INT) :: type;
       integer(C_INT) :: dir;
       integer(C_INT) :: func;
       real(C_DOUBLE), dimension(6) :: moments;
       real(C_DOUBLE) :: tau;
       real(C_DOUBLE) :: freq;
       real(C_DOUBLE), dimension(4) :: band;
       real(C_DOUBLE) :: ts
       real(C_DOUBLE) :: gamma
       type(C_PTR) :: time_file
    end type sem_source

    interface
       subroutine read_sem_config(config, spec, err) bind(c)
           use iso_c_binding
           import :: sem_config
           type(sem_config), intent(in) :: config
           character(C_CHAR), dimension(*) :: spec
           integer(C_INT), intent(out) :: err
       end subroutine read_sem_config

       subroutine dump_config(config) bind(c)
           use iso_c_binding
           import :: sem_config
           type(sem_config), intent(in) :: config
       end subroutine dump_config

       function strlen(s) bind(c)
           use iso_c_binding
           type(c_ptr), intent(in), value :: s
           integer(C_INT) :: strlen
       end function strlen
    end interface

contains

    function fromcstr(cstr)
        use iso_c_binding
        type(c_ptr), intent(in) :: cstr
        character(Len=MAX_FILE_SIZE) :: fromcstr
        character,pointer,dimension(:) :: ctemp
        integer, dimension(1) :: clen
        integer :: i
        clen(1) = strlen(cstr)
        call c_f_pointer(cstr, ctemp, clen)
        fromcstr=''
        do i=1,clen(1)
            fromcstr(i:i) = ctemp(i)
        end do
    end function fromcstr

subroutine read_gradient_desc(Tdomain, rg)
    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    integer :: i,j

    open(12,file=Tdomain%file_bassin,action="read",status="old")
    !   x_type 0 : on impose un materiau homogene dans chaque sous domaine
    !   x_type 1 : on impose un gradient de proprietes en fonction de z et ceci pour chaque colonne
    read(12,*) Tdomain%sBassin%x_type
    if ( Tdomain%sBassin%x_type == 2 ) then
        read(12,*) Tdomain%sBassin%zmin
        read(12,*) Tdomain%sBassin%zmax
        if ( rg==0) then
            write(*,*) ' gradient uniquement pour z entre ', Tdomain%sBassin%zmin,'  et ', Tdomain%sBassin%zmax
        endif
    endif
    !   nombre de colonnes
    read(12,*) Tdomain%sBassin%n_colonne
    allocate(Tdomain%sBassin%x_coord(0:Tdomain%sBassin%n_colonne))
    !  lecture des coordonnees en x des colonnes
    do i = 0,Tdomain%sBassin%n_colonne
        read(12,*) Tdomain%sBassin%x_coord(i)
        if ( rg==0) then
            write(*,*) ' colonne ',Tdomain%sBassin%x_coord(i)
        endif
    enddo
    !   nombre de couches
    read(12,*) Tdomain%sBassin%n_layer
    if ( rg==0) then
        write(*,*) ' n_layer ', Tdomain%sBassin%n_layer
    endif
    allocate(Tdomain%sBassin%z_layer(0:Tdomain%sBassin%n_colonne,0:Tdomain%sBassin%n_layer))
    allocate(  Tdomain%sBassin%z_rho(0:Tdomain%sBassin%n_colonne,0:Tdomain%sBassin%n_layer))
    allocate(   Tdomain%sBassin%z_Cp(0:Tdomain%sBassin%n_colonne,0:Tdomain%sBassin%n_layer))
    allocate(   Tdomain%sBassin%z_Cs(0:Tdomain%sBassin%n_colonne,0:Tdomain%sBassin%n_layer))
    ! on lit pour chaque colonne
    !  on lit chaque profondeur et les proprietes mecaniques associees
    do i = 0,Tdomain%sBassin%n_colonne
        do j = 0,Tdomain%sBassin%n_layer
            read(12,*) Tdomain%sBassin%z_layer(i,j),Tdomain%sBassin%z_rho(i,j),Tdomain%sBassin%z_Cp(i,j),Tdomain%sBassin%z_Cs(i,j)
            if ( rg==0) then
                write(*,*) ' gradient ',i,j
                write(*,*) Tdomain%sBassin%z_layer(i,j),Tdomain%sBassin%z_rho(i,j),Tdomain%sBassin%z_Cp(i,j),Tdomain%sBassin%z_Cs(i,j)
            endif
        enddo
    enddo
    if ( rg==0) then
        write(*,*) ' fin lecture gradient '
    endif
    close(12)
end subroutine read_gradient_desc

subroutine echo_input_params(Tdomain, rg)
    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    character(Len=MAX_FILE_SIZE) :: fnamef
    integer :: i

    call semname_read_input_spec(fnamef)

    !    debut ajout ecriture par un seul proc
    if ( rg == 0 ) then
        open(91,file=fnamef, form="formatted", status="unknown")

        write (91,*) Tdomain%Title_simulation
        write (91,*) Tdomain%TimeD%acceleration_scheme
        write (91,*) Tdomain%TimeD%velocity_scheme
        write (91,*) Tdomain%TimeD%duration
        write (91,*) Tdomain%TimeD%alpha
        write (91,*) Tdomain%TimeD%beta
        write (91,*) Tdomain%TimeD%gamma
        write (91,*) Tdomain%mesh_file
        write (91,*) Tdomain%material_file
        write (91,*) Tdomain%logicD%save_trace
        write (91,*) Tdomain%logicD%save_snapshots
        write (91,*) Tdomain%logicD%save_energy
        write (91,*) Tdomain%logicD%plot_grid
        write (91,*) Tdomain%logicD%run_exec
        write (91,*) Tdomain%logicD%run_debug
        write (91,*) Tdomain%logicD%run_echo

        if (Tdomain%logicD%save_trace) then
            write (91,*) Tdomain%station_file
        else
            write (91,*) " No parameter ned here"
        endif

        if (Tdomain%logicD%save_snapshots) then
            write (91,*) Tdomain%TimeD%time_snapshots
        else
            write (91,*) " No parameter ned here"
        endif

        write(11,*) Tdomain%logicD%Neumann, "  Neumann B.C.?"

        if (Tdomain%logicD%any_source) then
            write (91,*) Tdomain%n_source
            do i = 0, Tdomain%n_source - 1
                write (91,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Ysource, Tdomain%Ssource(i)%Zsource
                write (91,*) Tdomain%Ssource(i)%i_type_source
                if (Tdomain%Ssource(i)%i_type_source == 1 ) then
                    write (91,*) Tdomain%Ssource(i)%i_dir
                else
                    write (91,*) "No parameter need here"
                endif
                write (91,*) Tdomain%Ssource(i)%i_time_function
                write (91,*) Tdomain%Ssource(i)%tau_b
                write (91,*) Tdomain%Ssource(i)%cutoff_freq
            enddo
        else
            write (91,*)  "No available sources "
        endif
        write (91,*) "All right, runner ?"
        close (91)
    endif
    !   fin  ajout ecriture par un seul proc
end subroutine echo_input_params


subroutine read_mesh_file(Tdomain, rg)
    use sdomain
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    integer :: i,j, code, ok
    logical :: neumann_log

    !-- Reading mesh properties
    open(12, file=Tdomain%mesh_file, iostat=ok, status="old", form="formatted")
    if (ok/=0) then
        write (*,*) "Process ",rg, " can't open his mesh_file ", Tdomain%mesh_file
        stop
    endif
    read (12,*) Tdomain%n_dime
    if(Tdomain%n_dime /=3)   &
        stop "No general code for the time being: only 3D propagation"
    read(12,*) Tdomain%logicD%solid_fluid
    read(12,*) Tdomain%logicD%all_fluid
    if(rg == 0)then
        if(Tdomain%logicD%solid_fluid)then
            write(*,*) "  --> Propagation in solid-fluid media."
        else if(Tdomain%logicD%all_fluid)then
            write(*,*) "  --> Propagation in fluid media."
        else
            write(*,*) "  --> Propagation in solid media."
        end if
    end if
    call MPI_BARRIER(Tdomain%communicateur, code)
    read(12,*) neumann_log
    if(neumann_log .neqv. Tdomain%logicD%Neumann)  &
        stop "Introduction of Neumann B.C.: mesh and input files not in coincidence."
    read(12,*)   ! Global nodes for each proc.
    read (12,*) Tdomain%n_glob_nodes
    write(6,*) rg, ': nb noeuds',Tdomain%n_glob_nodes
    read (12,*) Tdomain%curve
    allocate (Tdomain%Coord_nodes(0:Tdomain%n_dime-1,0:Tdomain%n_glob_nodes-1))
    do i = 0,Tdomain%n_glob_nodes-1
        read(12,*) (Tdomain%Coord_nodes(j,i), j=0,Tdomain%n_dime-1)
    enddo
    read(12,*)  ! Elements
    read(12,*) Tdomain%n_elem
    write(6,*) rg, ': nb elts  ',Tdomain%n_elem
    allocate(Tdomain%specel(0:Tdomain%n_elem-1))
    do i=0,Tdomain%n_elem-1
        call init_element(Tdomain%specel(i))
    enddo
    read(12,*)  ! Materials
    read (12,*) Tdomain%n_mat
    do i = 0, Tdomain%n_elem-1
        read(12,*) Tdomain%specel(i)%mat_index, Tdomain%specel(i)%solid
    enddo
    read(12,*) ! Index of nodes for elements
    read(12,*) Tdomain%n_nodes
    do i = 0, Tdomain%n_elem-1
        allocate(Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1))
        read(12,*) Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1)
    enddo
    read(12,*)  ! Faces and elements properties related to faces
    read(12,*) Tdomain%n_face
    write(*,*) "***** NFACES=", Tdomain%n_face
    allocate(Tdomain%sFace(0:Tdomain%n_face-1))
    do i=0,Tdomain%n_face-1
        call init_face(Tdomain%sFace(i))
    enddo
    do i = 0, Tdomain%n_elem-1
        read(12,*) Tdomain%specel(i)%Near_Faces(0:5)
        read(12,*) Tdomain%specel(i)%Orient_Faces(0:5)
    enddo
    read(12,*)  ! Edges
    read(12,*) Tdomain%n_edge
    allocate(Tdomain%sEdge(0:Tdomain%n_edge-1))
    do i=0,Tdomain%n_edge-1
        call init_edge(Tdomain%sEdge(i))
    enddo
    do i = 0, Tdomain%n_elem-1
        read(12,*) Tdomain%specel(i)%Near_Edges(0:11)
        read(12,*) Tdomain%specel(i)%Orient_Edges(0:11)
    enddo
    read(12,*)  ! Vertices
    read(12,*) Tdomain%n_vertex
    allocate(Tdomain%sVertex(0:Tdomain%n_vertex-1))
    do i=0,Tdomain%n_vertex-1
        call init_vertex(Tdomain%sVertex(i))
    enddo
    do i = 0,Tdomain%n_elem-1
        read(12,*) Tdomain%specel(i)%Near_Vertices(0:7)
    enddo
    read(12,*)  ! relationship vertex <-> global node
    do i = 0,Tdomain%n_vertex-1
        read(12,*) Tdomain%sVertex(i)%global_numbering
    end do
    ! Solid-fluid properties, eventually
    !   AJOUTER Une routine de verif dans le fichier materiaux
    Tdomain%logicD%SF_local_present = .false.
    if(Tdomain%logicD%solid_fluid)then
        read(12,*)
        read(12,*) Tdomain%logicD%SF_local_present
        if(Tdomain%logicD%SF_local_present)then
            read(12,*) ! Solid-fluid properties
            read(12,*) ! SF faces
            read(12,*) Tdomain%SF%SF_n_faces
            allocate(Tdomain%SF%SF_face(0:Tdomain%SF%SF_n_faces-1))
            read(12,*) ! Edges and their orientation, for each SF face
            do i = 0, Tdomain%SF%SF_n_faces-1
                read(12,*) Tdomain%SF%SF_face(i)%Near_Edges(0:3)
                read(12,*) Tdomain%SF%SF_face(i)%Near_Edges_Orient(0:3)
            end do
            read(12,*) ! Vertices for each SF face
            do i = 0, Tdomain%SF%SF_n_faces-1
                read(12,*) Tdomain%SF%SF_face(i)%Near_Vertices(0:3)
            end do
            read(12,*) ! associated fluid (0) and solid (1) faces; orientation Solid/Fluid
            do i = 0, Tdomain%SF%SF_n_faces-1
                read(12,*) Tdomain%SF%SF_face(i)%Face(0), Tdomain%SF%SF_face(i)%Face(1),  &
                    Tdomain%SF%SF_face(i)%Orient_Face
            end do
            read(12,*) ! SF edges
            read(12,*) Tdomain%SF%SF_n_edges
            allocate(Tdomain%SF%SF_edge(0:Tdomain%SF%SF_n_edges-1))
            do i = 0, Tdomain%SF%SF_n_edges-1
                read(12,*) Tdomain%SF%SF_edge(i)%Edge(0), Tdomain%SF%SF_edge(i)%Edge(1),  &
                    Tdomain%SF%SF_edge(i)%Orient_Edge
            end do
            read(12,*) ! SF vertices
            read(12,*) Tdomain%SF%SF_n_vertices
            allocate(Tdomain%SF%SF_vertex(0:Tdomain%SF%SF_n_vertices-1))
            do i = 0, Tdomain%SF%SF_n_vertices-1
                read(12,*) Tdomain%SF%SF_vertex(i)%Vertex(0), Tdomain%SF%SF_vertex(i)%Vertex(1)
            end do
        end if
    end if

    ! Neumann B.C. properties, eventually
    Tdomain%logicD%Neumann_local_present = .false.
    if(Tdomain%logicD%Neumann)then
        read(12,*)
        read(12,*) Tdomain%logicD%Neumann_local_present
        if(Tdomain%logicD%Neumann_local_present)then
            read(12,*) ! Neumann properties
            read(12,*) ! Neumann faces
            read(12,*) Tdomain%Neumann%Neu_n_faces
            allocate(Tdomain%Neumann%Neu_face(0:Tdomain%Neumann%Neu_n_faces-1))
            read(12,*) ! Edges and their orientation, for each Neumann face
            do i = 0, Tdomain%Neumann%Neu_n_faces-1
                read(12,*) Tdomain%Neumann%Neu_face(i)%Near_Edges(0:3)
                read(12,*) Tdomain%Neumann%Neu_face(i)%Near_Edges_Orient(0:3)
            end do
            read(12,*) ! Vertices for each Neumann face
            do i = 0, Tdomain%Neumann%Neu_n_faces-1
                read(12,*) Tdomain%Neumann%Neu_face(i)%Near_Vertices(0:3)
            end do
            read(12,*) ! associated face
            do i = 0, Tdomain%Neumann%Neu_n_faces-1
                read(12,*) Tdomain%Neumann%Neu_face(i)%Face
            end do
            read(12,*) ! Neumann edges
            read(12,*) Tdomain%Neumann%Neu_n_edges
            allocate(Tdomain%Neumann%Neu_edge(0:Tdomain%Neumann%Neu_n_edges-1))
            do i = 0, Tdomain%Neumann%Neu_n_edges-1
                read(12,*) Tdomain%Neumann%Neu_edge(i)%Edge
            end do
            read(12,*) ! Neumann vertices
            read(12,*) Tdomain%Neumann%Neu_n_vertices
            allocate(Tdomain%Neumann%Neu_vertex(0:Tdomain%Neumann%Neu_n_vertices-1))
            do i = 0, Tdomain%Neumann%Neu_n_vertices-1
                read(12,*) Tdomain%Neumann%Neu_vertex(i)%Vertex
            end do
        end if
    end if

    read(12,*)
    read(12,*)  ! Interproc communications
    read(12,*) Tdomain%n_proc
    allocate (Tdomain%sComm(0:Tdomain%n_proc-1))
    do i = 0,Tdomain%n_proc-1
        read(12,*) Tdomain%sComm(i)%nb_faces, Tdomain%sComm(i)%nb_edges, Tdomain%sComm(i)%nb_vertices
        if(Tdomain%logicD%SF_local_present)then
            read(12,*) Tdomain%sComm(i)%SF_nf_shared, Tdomain%sComm(i)%SF_ne_shared, Tdomain%sComm(i)%SF_nv_shared
        else
            Tdomain%sComm(i)%SF_nf_shared = 0
            Tdomain%sComm(i)%SF_ne_shared = 0
            Tdomain%sComm(i)%SF_nv_shared = 0
        end if
        if(Tdomain%logicD%Neumann_local_present)then
            read(12,*) Tdomain%sComm(i)%Neu_ne_shared, Tdomain%sComm(i)%Neu_nv_shared
        else
            Tdomain%sComm(i)%Neu_ne_shared = 0
            Tdomain%sComm(i)%Neu_nv_shared = 0
        end if
        if(Tdomain%sComm(i)%nb_faces > 0)then
            allocate(Tdomain%sComm(i)%faces(0:Tdomain%sComm(i)%nb_faces-1))
            allocate(Tdomain%sComm(i)%orient_faces(0:Tdomain%sComm(i)%nb_faces-1))
            do j = 0,Tdomain%sComm(i)%nb_faces-1
                read(12,*) Tdomain%sComm(i)%faces(j),Tdomain%sComm(i)%orient_faces(j)
            enddo
        else
            nullify(Tdomain%sComm(i)%faces)
            nullify(Tdomain%sComm(i)%orient_faces)
        endif
        if(Tdomain%sComm(i)%nb_edges > 0)then
            allocate(Tdomain%sComm(i)%edges(0:Tdomain%sComm(i)%nb_edges-1))
            allocate(Tdomain%sComm(i)%orient_edges(0:Tdomain%sComm(i)%nb_edges-1))
            do j = 0,Tdomain%sComm(i)%nb_edges-1
                read(12,*) Tdomain%sComm(i)%edges(j),Tdomain%sComm(i)%orient_edges(j)
            enddo
        else
            nullify(Tdomain%sComm(i)%edges)
            nullify(Tdomain%sComm(i)%orient_edges)
        endif
        if(Tdomain%sComm(i)%nb_vertices > 0)then
            allocate(Tdomain%sComm(i)%vertices(0:Tdomain%sComm(i)%nb_vertices-1))
            do j = 0,Tdomain%sComm(i)%nb_vertices-1
                read(12,*) Tdomain%sComm(i)%vertices(j)
            enddo
        endif
        if(Tdomain%logicD%SF_local_present)then
            if(Tdomain%sComm(i)%SF_nf_shared > 0)then
                allocate(Tdomain%sComm(i)%SF_faces_shared(0:Tdomain%sComm(i)%SF_nf_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_nf_shared-1
                    read(12,*) Tdomain%sComm(i)%SF_faces_shared(j)
                enddo
            else
                nullify(Tdomain%sComm(i)%SF_faces_shared)
            endif
            if(Tdomain%sComm(i)%SF_ne_shared > 0)then
                allocate(Tdomain%sComm(i)%SF_edges_shared(0:Tdomain%sComm(i)%SF_ne_shared-1))
                allocate(Tdomain%sComm(i)%SF_mapping_edges_shared(0:Tdomain%sComm(i)%SF_ne_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_ne_shared-1
                    read(12,*) Tdomain%sComm(i)%SF_edges_shared(j),Tdomain%sComm(i)%SF_mapping_edges_shared(j)
                enddo
            else
                nullify(Tdomain%sComm(i)%SF_edges_shared)
                nullify(Tdomain%sComm(i)%SF_mapping_edges_shared)
            endif
            if(Tdomain%sComm(i)%SF_nv_shared > 0)then
                allocate(Tdomain%sComm(i)%SF_vertices_shared(0:Tdomain%sComm(i)%SF_nv_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_nv_shared-1
                    read(12,*) Tdomain%sComm(i)%SF_vertices_shared(j)
                enddo
            else
                nullify(Tdomain%sComm(i)%SF_vertices_shared)
            endif
        end if
        if(Tdomain%logicD%Neumann_local_present)then
            if(Tdomain%sComm(i)%Neu_ne_shared > 0)then
                allocate(Tdomain%sComm(i)%Neu_edges_shared(0:Tdomain%sComm(i)%Neu_ne_shared-1))
                allocate(Tdomain%sComm(i)%Neu_mapping_edges_shared(0:Tdomain%sComm(i)%Neu_ne_shared-1))
                do j = 0,Tdomain%sComm(i)%Neu_ne_shared-1
                    read(12,*) Tdomain%sComm(i)%Neu_edges_shared(j),Tdomain%sComm(i)%Neu_mapping_edges_shared(j)
                enddo
            else
                nullify(Tdomain%sComm(i)%Neu_edges_shared)
                nullify(Tdomain%sComm(i)%Neu_mapping_edges_shared)
            endif
            if(Tdomain%sComm(i)%Neu_nv_shared > 0)then
                allocate(Tdomain%sComm(i)%Neu_vertices_shared(0:Tdomain%sComm(i)%Neu_nv_shared-1))
                do j = 0,Tdomain%sComm(i)%Neu_nv_shared-1
                    read(12,*) Tdomain%sComm(i)%Neu_vertices_shared(j)
                enddo
            else
                nullify(Tdomain%sComm(i)%Neu_vertices_shared)
            endif
        else
            nullify(Tdomain%sComm(i)%Neu_edges_shared)
            nullify(Tdomain%sComm(i)%Neu_mapping_edges_shared)
            nullify(Tdomain%sComm(i)%Neu_vertices_shared)
        end if

    end do
    close(12)

    write(*,*) "Mesh read correctly for proc #", rg
end subroutine read_mesh_file

subroutine finalize_mesh_connectivity(Tdomain, rg)
    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    integer :: i, j, k, n, nf, nnf, mat, ne, nv, nne, nnv
    integer :: icount, n_aus


    ! faces and edges => which element?
    do i = 0,Tdomain%n_face-1
        do j = 0, Tdomain%n_elem-1
            do k = 0,5
                if(Tdomain%specel(j)%Near_Faces(k) == i)then
                    Tdomain%sFace(i)%Which_Elem = j
                endif
            enddo
        enddo
    enddo

    do i = 0,Tdomain%n_edge-1
        allocate(Tdomain%sEdge(i)%Which_Elem(0:11))
        allocate(Tdomain%sEdge(i)%Which_EdgeinElem(0:11))
        Tdomain%sEdge(i)%Which_Elem = -1
        Tdomain%sEdge(i)%Which_EdgeinElem = -1
        icount = 0
        do j = 0, Tdomain%n_elem - 1
            do k=0,11
                if ( Tdomain%specel(j)%Near_Edges(k) == i )  then
                    Tdomain%sEdge(i)%Which_Elem(icount) = j
                    Tdomain%sEdge(i)%Which_EdgeinElem(icount) = k
                    icount = icount+1
                endif
            enddo
        enddo
    enddo

    ! material => time steps ; solid/liquid attribution
    do n = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(n)%mat_index
        do nf = 0,5
            nnf = Tdomain%specel(n)%Near_Faces(nf)
            !Tdomain%sFace(nnf)%mat_index = mat
            Tdomain%sFace(nnf)%solid = Tdomain%specel(n)%solid
        end do
        do ne = 0,11
            nne = Tdomain%specel(n)%Near_edges(ne)
            !Tdomain%sEdge(nne)%mat_index = mat
            Tdomain%sEdge(nne)%solid = Tdomain%specel(n)%solid
        end do
        do nv = 0,7
            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            !Tdomain%sVertex(nnv)%mat_index = mat
            Tdomain%sVertex(nnv)%solid = Tdomain%specel(n)%solid
        end do
    end do

    !- Neumann local properties
    if(Tdomain%logicD%neumann_local_present)then
        do nf = 0, Tdomain%Neumann%Neu_n_faces-1
            n_aus = Tdomain%Neumann%Neu_Face(nf)%Face
            Tdomain%Neumann%Neu_Face(nf)%ngll1 = Tdomain%sFace(n_aus)%ngll1
            Tdomain%Neumann%Neu_Face(nf)%ngll2 = Tdomain%sFace(n_aus)%ngll2
            Tdomain%Neumann%Neu_Face(nf)%dir = Tdomain%sFace(n_aus)%dir
        enddo
        do ne = 0, Tdomain%Neumann%Neu_n_edges-1
            n_aus = Tdomain%Neumann%Neu_Edge(ne)%Edge
            Tdomain%Neumann%Neu_Edge(ne)%ngll = Tdomain%sEdge(n_aus)%ngll
        enddo
    endif

    !- Solid/fluid interfaces local properties
    if(Tdomain%logicD%SF_local_present)then
        do nf = 0, Tdomain%SF%SF_n_faces-1
            n_aus = Tdomain%SF%SF_Face(nf)%Face(0)
            if(n_aus < 0) n_aus = Tdomain%SF%SF_Face(nf)%Face(1)
            Tdomain%SF%SF_Face(nf)%ngll1 = Tdomain%sFace(n_aus)%ngll1
            Tdomain%SF%SF_Face(nf)%ngll2 = Tdomain%sFace(n_aus)%ngll2
            Tdomain%SF%SF_Face(nf)%dir = Tdomain%sFace(n_aus)%dir
        enddo
        do ne = 0, Tdomain%SF%SF_n_edges-1
            n_aus = Tdomain%SF%SF_Edge(ne)%Edge(0)
            if(n_aus < 0) n_aus = Tdomain%SF%SF_Edge(ne)%Edge(1)
            Tdomain%SF%SF_Edge(ne)%ngll = Tdomain%sEdge(n_aus)%ngll
        enddo
    endif

end subroutine finalize_mesh_connectivity

subroutine read_material_file(Tdomain, rg)
    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    character(Len=MAX_FILE_SIZE) :: fnamef
    integer :: i, j, n_aus, npml, mat, ne, nf
    logical, dimension(:), allocatable :: L_Face, L_Edge
    real :: dtmin

    npml = 0
    allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))

    call semname_read_inputmesh_parametrage(Tdomain%material_file,fnamef)
    open (13, file=fnamef, status="old", form="formatted")

    read(13,*) n_aus

    if(n_aus /= Tdomain%n_mat)   &
        stop "Incompatibility between the mesh file and the material file for n_mat"



    if (Tdomain%aniso) then
        print *,"The code can't put anisotropy in a homogeneous media"
        stop
    endif
    do i = 0,Tdomain%n_mat-1
        read(13,*) Tdomain%sSubDomain(i)%material_type, Tdomain%sSubDomain(i)%Pspeed, &
            Tdomain%sSubDomain(i)%Sspeed, Tdomain%sSubDomain(i)%dDensity,      &
            Tdomain%sSubDomain(i)%NGLLx, Tdomain%sSubDomain(i)%NGLLy, Tdomain%sSubDomain(i)%NGLLz, &
            Tdomain%sSubDomain(i)%Dt, Tdomain%sSubDomain(i)%Qpression,  Tdomain%sSubDomain(i)%Qmu

        if(rg==0) then
            write (*,*) 'Material :', i
            write (*,*) 'type:', Tdomain%sSubDomain(i)%material_type
            write (*,*) 'Pspeed:', Tdomain%sSubDomain(i)%Pspeed
            write (*,*) 'Sspeed:', Tdomain%sSubDomain(i)%Sspeed
            write (*,*) 'Density:', Tdomain%sSubDomain(i)%dDensity
            write (*,*) 'NGLL:', Tdomain%sSubDomain(i)%NGLLx, Tdomain%sSubDomain(i)%NGLLy, Tdomain%sSubDomain(i)%NGLLz
            write (*,*) 'Dt:', Tdomain%sSubDomain(i)%Dt
            write (*,*) 'Qp:', Tdomain%sSubDomain(i)%Qpression
            write (*,*) 'Qmu:', Tdomain%sSubDomain(i)%Qmu
        endif

        call Lame_coefficients (Tdomain%sSubDomain(i))
        if(rg==0) &
            print*,' lame ',Tdomain%sSubDomain(i)%DMu,Tdomain%sSubDomain(i)%DLambda ,Tdomain%sSubDomain(i)%DKappa
        if (Tdomain%sSubDomain(i)%material_type == "P" .or. Tdomain%sSubDomain(i)%material_type == "L")  then
            Tdomain%sSubDomain(i)%wpml = npml
            npml = npml + 1
        endif
    enddo

    Tdomain%any_PML = .false.
    Tdomain%any_FPML = .false.
    if(npml > 0) then
        Tdomain%any_PML = .true.
        read(13,*); read(13,*)
        do i = 0,Tdomain%n_mat-1
            if(Tdomain%sSubdomain(i)%material_type == "P" .or.     &
                Tdomain%sSubDomain(i)%material_type == "L") then
                read(13,*) Tdomain%sSubdomain(i)%Filtering, Tdomain%sSubdomain(i)%npow,    &
                    Tdomain%sSubdomain(i)%Apow, Tdomain%sSubdomain(i)%Px,                  &
                    Tdomain%sSubdomain(i)%Left, Tdomain%sSubdomain(i)%Py,                  &
                    Tdomain%sSubdomain(i)%Forward, Tdomain%sSubdomain(i)%Pz,               &
                    Tdomain%sSubdomain(i)%Down, Tdomain%sSubdomain(i)%freq
                if(Tdomain%sSubdomain(i)%Filtering) Tdomain%any_FPML = .true.
            endif
        enddo
    endif
    close(13)

    !- GLL properties in elements, on faces, edges.
    allocate(L_Face(0:Tdomain%n_face-1))
    L_Face = .true.
    allocate(L_Edge(0:Tdomain%n_edge-1))
    L_Edge = .true.
    do i = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(i)%mat_index
        Tdomain%specel(i)%ngllx = Tdomain%sSubDomain(mat)%NGLLx
        Tdomain%specel(i)%nglly = Tdomain%sSubDomain(mat)%NGLLy
        Tdomain%specel(i)%ngllz = Tdomain%sSubDomain(mat)%NGLLz
        do j = 0,5
            nf = Tdomain%specel(i)%Near_Faces(j)
            if(L_Face(nf) .and. Tdomain%specel(i)%Orient_Faces(j) == 0)then
                L_Face(nf) = .false. !L_Face est pour eviter que l'on modifie plusieurs fois la meme face
                if(j == 0 .or. j == 5)then
                    Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%ngllx
                    Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%nglly
                else if(j == 1 .or. j == 3)then
                    Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%ngllx
                    Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%ngllz
                else
                    Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%nglly
                    Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%ngllz
                endif
                Tdomain%sFace(nf)%dir = j
            endif
        enddo
        do j = 0,11
            ne = Tdomain%specel(i)%Near_Edges(j)
            if(L_Edge(ne) .and. Tdomain%specel(i)%Orient_Edges(j) == 0)then
                L_Edge(ne) = .false.
                if (j == 0 .or. j == 2 .or. j == 5 .or. j == 9)then
                    Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%ngllx
                else if (j == 1 .or. j == 3 .or. j == 8 .or. j == 11)then
                    Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%nglly
                else
                    Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%ngllz
                endif
            endif
        enddo
    enddo
    deallocate(L_Face,L_Edge)

    ! MODIF to be done here.: Gaetano's formulae are wrong for filtering PMLs
    do i = 0, Tdomain%n_mat-1
        if((Tdomain%sSubdomain(i)%material_type == "P" .or.   &
            Tdomain%sSubdomain(i)%material_type == "L") .and. &
            Tdomain%sSubdomain(i)%Filtering)                  &
            Tdomain%sSubdomain(i)%freq = exp(-Tdomain%sSubdomain(i)%freq*Tdomain%sSubdomain(i)%dt/2)
    enddo

    dtmin = 1e20
    do i = 0,Tdomain%n_mat-1
        if(Tdomain%sSubDomain(i)%Dt < dtmin) dtmin = Tdomain%sSubDomain(i)%Dt
    enddo
    Tdomain%TimeD%dtmin = dtmin
    if(dtmin > 0)then
        Tdomain%TimeD%ntimeMax = int(Tdomain%TimeD%Duration/dtmin)
    else
        stop "Your dt min is zero : verify it"
    endif
    if(rg==0) &
        print *,'ntimemax',Tdomain%TimeD%ntimeMax,Tdomain%TimeD%Duration,dtmin

end subroutine read_material_file


subroutine create_sem_sources(Tdomain, config)
    use sdomain
    implicit none
    type(domain), intent(inout)  :: Tdomain
    type(sem_config), intent(in) :: config
    type(sem_source), pointer :: src
    integer :: nsrc

    Tdomain%n_source = config%nsources
    allocate (Tdomain%Ssource(0:Tdomain%n_source-1))

    call c_f_pointer(config%source, src)
    nsrc = 0
    do while(associated(src))
        Tdomain%Ssource(nsrc)%Xsource = src%coords(1)
        Tdomain%Ssource(nsrc)%Ysource = src%coords(2)
        Tdomain%Ssource(nsrc)%Zsource = src%coords(3)
        Tdomain%Ssource(nsrc)%i_type_source = src%type
        ! Comportement temporel
        Tdomain%Ssource(nsrc)%i_time_function = src%func
        Tdomain%Ssource(nsrc)%cutoff_freq = src%freq ! func=2,4
        Tdomain%Ssource(nsrc)%tau_b = src%tau ! func=1,2,3,4,5
        Tdomain%Ssource(nsrc)%fh = src%band  ! func=3
        Tdomain%Ssource(nsrc)%gamma = src%gamma ! func=4
        Tdomain%Ssource(nsrc)%ts = src%ts   ! func=4
        ! Comportement Spacial
        ! i_type_source==1
        Tdomain%Ssource(nsrc)%i_dir = src%dir
        ! i_type_source==2
        Tdomain%Ssource(nsrc)%Moment(0,0) = src%moments(1)
        Tdomain%Ssource(nsrc)%Moment(1,1) = src%moments(2)
        Tdomain%Ssource(nsrc)%Moment(2,1) = src%moments(3)
        Tdomain%Ssource(nsrc)%Moment(0,1) = src%moments(4)
        Tdomain%Ssource(nsrc)%Moment(1,0) = src%moments(4)
        Tdomain%Ssource(nsrc)%Moment(0,2) = src%moments(5)
        Tdomain%Ssource(nsrc)%Moment(2,0) = src%moments(5)
        Tdomain%Ssource(nsrc)%Moment(1,2) = src%moments(6)
        Tdomain%Ssource(nsrc)%Moment(2,1) = src%moments(6)

        nsrc = nsrc + 1
        Tdomain%logicD%any_source = .true.
        call c_f_pointer(src%next, src)
    end do


end subroutine create_sem_sources

subroutine read_input (Tdomain, rg, code)
    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout)  :: Tdomain
    integer, intent(in)          :: rg
    integer, intent(out)         :: code
    character(Len=MAX_FILE_SIZE) :: fnamef
    type(sem_config)             :: config
    logical                      :: logic_scheme

    call semname_read_input_input(fnamef)

    call read_sem_config(config, trim(fnamef)//C_NULL_CHAR, code)

    if (code/=1) then
        stop 1
    endif
    call dump_config(config)

    ! On copie les parametres renvoyes dans la structure C
    Tdomain%Title_simulation = fromcstr(config%run_name)
    Tdomain%TimeD%acceleration_scheme = config%accel_scheme .ne. 0
    Tdomain%TimeD%velocity_scheme = config%veloc_scheme .ne. 0
    Tdomain%TimeD%duration = config%sim_time
    Tdomain%TimeD%alpha = config%alpha
    Tdomain%TimeD%beta = config%beta
    Tdomain%TimeD%gamma = config%gamma
    Tdomain%TimeD%courant = config%courant

    Tdomain%mesh_file = fromcstr(config%mesh_file)
    call semname_read_input_meshfile(rg,Tdomain%mesh_file,fnamef)
    Tdomain%mesh_file = fnamef
    Tdomain%aniso = config%anisotropy .ne. 0
    Tdomain%material_file = fromcstr(config%mat_file)
    Tdomain%logicD%save_trace = config%save_traces .ne. 0
    Tdomain%logicD%save_snapshots = config%save_snap .ne. 0
    ! MODIF ICI: energie? deformation?..
    !Tdomain%logicD%save_energy = !?
    Tdomain%logicD%save_restart = config%prorep .ne. 0
    !Tdomain%logicD%plot_grid
    !Tdomain%logicD%run_exec
    !Tdomain%logicD%run_debug
    !Tdomain%logicD%run_echo
    Tdomain%logicD%run_restart = config%prorep .ne. 0
    Tdomain%TimeD%iter_reprise = config%prorep_iter
    Tdomain%TimeD%ncheck = config%prorep_iter ! frequence de sauvegarde
    Tdomain%station_file = fromcstr(config%station_file)
    Tdomain%TimeD%ntrace = config%traces_interval ! XXX
    Tdomain%TimeD%time_snapshots = config%snap_interval
    logic_scheme = Tdomain%TimeD%acceleration_scheme .neqv. Tdomain%TimeD%velocity_scheme
    if(.not. logic_scheme) then
        stop "Both acceleration and velocity schemes: no compatibility, chose only one."
    end if

    ! Amortissement
    Tdomain%n_sls = config%nsolids
    Tdomain%T1_att = config%atn_band(1)
    Tdomain%T2_att = config%atn_band(2)
    Tdomain%T0_modele = config%atn_period
    write(*,*) "SLS=", Tdomain%n_sls


    ! Neumann boundary conditions? If yes: geometrical properties read in the mesh files.
    Tdomain%logicD%Neumann = config%neu_present /= 0
    Tdomain%Neumann%Neu_Param%what_bc = 'S'
    Tdomain%Neumann%Neu_Param%mat_index = config%neu_mat
    if (config%neu_type==1) then
        Tdomain%Neumann%Neu_Param%wtype = 'P'
    else
        Tdomain%Neumann%Neu_Param%wtype = 'S'
    end if
    Tdomain%Neumann%Neu_Param%lx = config%neu_L(1)
    Tdomain%Neumann%Neu_Param%ly = config%neu_L(2)
    Tdomain%Neumann%Neu_Param%lz = config%neu_L(3)
    Tdomain%Neumann%Neu_Param%xs = config%neu_C(1)
    Tdomain%Neumann%Neu_Param%ys = config%neu_C(2)
    Tdomain%Neumann%Neu_Param%zs = config%neu_C(3)
    Tdomain%Neumann%Neu_Param%f0 = config%neu_f0


    ! Create sources from C structures
    call create_sem_sources(Tdomain, config)

    !call echo_input_params(Tdomain, rg)

    !- Parametrage super object desactive
    Tdomain%logicD%super_object_local_present = .false.

    call read_mesh_file(Tdomain, rg)

    !---   Properties of materials.
    call read_material_file(Tdomain, rg)

    call finalize_mesh_connectivity(Tdomain, rg)

end subroutine read_input

end module semconfig
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
