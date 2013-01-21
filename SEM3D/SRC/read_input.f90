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


subroutine read_attn_desc(Tdomain, rg)
    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    character(Len=MAX_FILE_SIZE) :: fnamef
    integer :: unit_att
    logical :: trouve

    unit_att=22
    call semname_read_input_amortissement(fnamef)
    inquire(file=fnamef, exist=trouve)
    if( .not. trouve) then
        Tdomain%n_sls=0
        if(rg==0) write(*,*) 'No attenuation model: file not found: ', trim(adjustl(fnamef))
        return
    endif

    Tdomain%n_sls=1
    open (unit=unit_att,file=fnamef,form="formatted",status="old")

    !  modif mariotti fevrier 2007 cea
    !   n_sls est le nombre de sous decoupage en pas de frequence pour avoir un
    !   facteur de qualite a peu pres constant dans l intervalle
    !   voir article de Liu Kanamori du GJRA de 1976
    read (unit_att,*) Tdomain%n_sls

    if (Tdomain%n_sls>0) then
        read (unit_att,*) Tdomain%T1_att, Tdomain%T2_att
        read (unit_att,*) Tdomain%T0_modele
        if(rg==0) &
            write(*,*) ' modele attenuation sur temps ', Tdomain%T1_att, Tdomain%T2_att,Tdomain%T0_modele
    else
        if(rg==0) &
            write(*,*) 'pas de  modele attenuation sur temps '
    endif

    close(unit_att)
end subroutine read_attn_desc

subroutine read_source_desc(Tdomain, rg)
    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    character(Len=MAX_FILE_SIZE) :: fnamef
    integer :: unit_src, i
    logical :: trouve

    unit_src=21
    Tdomain%logicD%any_source=.true.
    call semname_read_input_source(fnamef)
    inquire(file=fnamef, exist=trouve)
    if( .not. trouve) then
        Tdomain%logicD%any_source=.false.
        if(rg==0) &
            write(*,*) 'No source for Sem: could not find file:', trim(adjustl(fnamef))
        return
    endif

    open (unit=unit_src,file=fnamef,form="formatted",status="old")

    read (unit_src,*) Tdomain%n_source

    if(rg==0) write (*,*) 'nombre de sources ', Tdomain%n_source
    allocate (Tdomain%Ssource(0:Tdomain%n_source-1))

    do i = 0, Tdomain%n_source - 1
        read (unit_src,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Ysource, Tdomain%Ssource(i)%Zsource
        read (unit_src,*) Tdomain%Ssource(i)%i_type_source
        if (rg==0) then
            write (*,*) 'SOURCE:', i
            write (*,*) '  type spatial :',Tdomain%Ssource(i)%i_type_source
        end if

        select case (Tdomain%Ssource(i)%i_type_source)
        case (1)

            read (unit_src,*) Tdomain%Ssource(i)%i_dir
            if (rg==0) then
                write (*,*) '  direction ',Tdomain%Ssource(i)%i_dir
            end if
        case (2)
            read (unit_src,*) Tdomain%Ssource(i)%Moment(0,0), Tdomain%Ssource(i)%Moment(1,1), Tdomain%Ssource(i)%Moment(2,2)
            read (unit_src,*) Tdomain%Ssource(i)%Moment(0,1), Tdomain%Ssource(i)%Moment(0,2), Tdomain%Ssource(i)%Moment(1,2)
            Tdomain%Ssource(i)%Moment(1,0) = Tdomain%Ssource(i)%Moment(0,1)
            Tdomain%Ssource(i)%Moment(2,0) = Tdomain%Ssource(i)%Moment(0,2)
            Tdomain%Ssource(i)%Moment(2,1) = Tdomain%Ssource(i)%Moment(1,2)

        end select

        read (unit_src,*) Tdomain%Ssource(i)%i_time_function
        read (unit_src,*) Tdomain%Ssource(i)%tau_b
        if(rg==0) then
            write (*,*) '  type temporel :', Tdomain%Ssource(i)%i_time_function,  Tdomain%Ssource(i)%tau_b
        endif

        ! Les parametres suivants dependent du type de source
        select case (Tdomain%Ssource(i)%i_time_function)
        case (2)
            read (unit_src,*) Tdomain%Ssource(i)%cutoff_freq
            if (rg==0) then
                write (*,*) '  Cutoff:', Tdomain%Ssource(i)%cutoff_freq
            end if
        case (3)
            read (unit_src,*) Tdomain%Ssource(i)%fh(0:3)
        case (4)
            !   Gabor signal
            read (unit_src,*) Tdomain%Ssource(i)%cutoff_freq
            read (unit_src,*) Tdomain%Ssource(i)%gama
            read (unit_src,*) Tdomain%Ssource(i)%ts
            if (rg==0) then
                write(*,*)  '   Gabor signal '
                write (*,*) '  frequencep:', Tdomain%Ssource(i)%cutoff_freq
                write (*,*) '  gamma:', Tdomain%Ssource(i)%gama
                write (*,*) '  ts:', Tdomain%Ssource(i)%ts
            end if
        case (5)
            read (unit_src,*) Tdomain%Ssource(i)%time_file
        end select
    enddo

    close(unit_src)

end subroutine read_source_desc


subroutine read_gradient_desc(Tdomain, rg)
    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    character(Len=MAX_FILE_SIZE) :: fnamef
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
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    character(Len=MAX_FILE_SIZE) :: fnamef
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
        end if
        if(Tdomain%logicD%Neumann_local_present)then
            read(12,*) Tdomain%sComm(i)%Neu_ne_shared, Tdomain%sComm(i)%Neu_nv_shared
        end if
        if(Tdomain%sComm(i)%nb_faces > 0)then
            allocate(Tdomain%sComm(i)%faces(0:Tdomain%sComm(i)%nb_faces-1))
            allocate(Tdomain%sComm(i)%orient_faces(0:Tdomain%sComm(i)%nb_faces-1))
            do j = 0,Tdomain%sComm(i)%nb_faces-1
                read(12,*) Tdomain%sComm(i)%faces(j),Tdomain%sComm(i)%orient_faces(j)
            enddo
        endif
        if(Tdomain%sComm(i)%nb_edges > 0)then
            allocate(Tdomain%sComm(i)%edges(0:Tdomain%sComm(i)%nb_edges-1))
            allocate(Tdomain%sComm(i)%orient_edges(0:Tdomain%sComm(i)%nb_edges-1))
            do j = 0,Tdomain%sComm(i)%nb_edges-1
                read(12,*) Tdomain%sComm(i)%edges(j),Tdomain%sComm(i)%orient_edges(j)
            enddo
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
            endif
            if(Tdomain%sComm(i)%SF_ne_shared > 0)then
                allocate(Tdomain%sComm(i)%SF_edges_shared(0:Tdomain%sComm(i)%SF_ne_shared-1))
                allocate(Tdomain%sComm(i)%SF_mapping_edges_shared(0:Tdomain%sComm(i)%SF_ne_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_ne_shared-1
                    read(12,*) Tdomain%sComm(i)%SF_edges_shared(j),Tdomain%sComm(i)%SF_mapping_edges_shared(j)
                enddo
            endif
            if(Tdomain%sComm(i)%SF_nv_shared > 0)then
                allocate(Tdomain%sComm(i)%SF_vertices_shared(0:Tdomain%sComm(i)%SF_nv_shared-1))
                do j = 0,Tdomain%sComm(i)%SF_nv_shared-1
                    read(12,*) Tdomain%sComm(i)%SF_vertices_shared(j)
                enddo
            endif
        end if
        if(Tdomain%logicD%Neumann_local_present)then
            if(Tdomain%sComm(i)%Neu_ne_shared > 0)then
                allocate(Tdomain%sComm(i)%Neu_edges_shared(0:Tdomain%sComm(i)%Neu_ne_shared-1))
                allocate(Tdomain%sComm(i)%Neu_mapping_edges_shared(0:Tdomain%sComm(i)%Neu_ne_shared-1))
                do j = 0,Tdomain%sComm(i)%Neu_ne_shared-1
                    read(12,*) Tdomain%sComm(i)%Neu_edges_shared(j),Tdomain%sComm(i)%Neu_mapping_edges_shared(j)
                enddo
            endif
            if(Tdomain%sComm(i)%Neu_nv_shared > 0)then
                allocate(Tdomain%sComm(i)%Neu_vertices_shared(0:Tdomain%sComm(i)%Neu_nv_shared-1))
                do j = 0,Tdomain%sComm(i)%Neu_nv_shared-1
                    read(12,*) Tdomain%sComm(i)%Neu_vertices_shared(j)
                enddo
            endif
        end if

    end do
    close(12)

    write(*,*) "Mesh read correctly for proc #", rg
end subroutine read_mesh_file

subroutine read_material_file(Tdomain, rg)
    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg
    character(Len=MAX_FILE_SIZE) :: fnamef
    integer :: i, n_aus, npml

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
end subroutine read_material_file


subroutine read_input (Tdomain, rg, code)

    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg

    logical :: logic_scheme, sortie, neumann_log
    logical, dimension(:), allocatable :: L_Face, L_Edge
    integer :: length,i,j,n_aus,mat,ok,nf,ne,nv,k,icount,n,i_aus,  &
        ipoint,code,nnf,nne,nnv
    real :: dtmin, x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7
    character(Len=MAX_FILE_SIZE) :: fnamef

    call semname_read_input_input(fnamef)

    open (11,file=fnamef,form="formatted",status="old")
    read(11,*) Tdomain%Title_simulation
    if (rg==0) write (*,*) 'title:', Tdomain%Title_simulation
    read(11,*) Tdomain%TimeD%acceleration_scheme
    if (rg==0) write (*,*) 'accel:', Tdomain%TimeD%acceleration_scheme
    read(11,*) Tdomain%TimeD%velocity_scheme
    if (rg==0) write (*,*) 'veloc:', Tdomain%TimeD%velocity_scheme
    read(11,*) Tdomain%TimeD%duration
    if (rg==0) write (*,*) 'duration:', Tdomain%TimeD%duration
    read(11,*) Tdomain%TimeD%alpha
    if (rg==0) write (*,*) 'alpha:', Tdomain%TimeD%alpha
    read(11,*) Tdomain%TimeD%beta
    if (rg==0) write (*,*) 'beta:', Tdomain%TimeD%beta
    read(11,*) Tdomain%TimeD%gamma
    if (rg==0) write (*,*) 'gamma:', Tdomain%TimeD%gamma
    read(11,*) Tdomain%mesh_file
    if (rg==0) write (*,*) 'mesh_file:', Tdomain%mesh_file
    call semname_read_input_meshfile(rg,Tdomain%mesh_file,fnamef)
    Tdomain%mesh_file = fnamef
    write (*,*)  rg,':',' file=', Tdomain%mesh_file
    call MPI_BARRIER(Tdomain%communicateur, code)

    Tdomain%aniso = .FALSE.

    read(11,*) Tdomain%material_file
    if (rg==0) write (*,*) 'material_file:', Tdomain%material_file
    read(11,*) Tdomain%logicD%save_trace
    if (rg==0) write (*,*) 'save_trace:', Tdomain%logicD%save_trace
    read(11,*) Tdomain%logicD%save_snapshots
    if (rg==0) write (*,*) 'save_snapshots:', Tdomain%logicD%save_snapshots
    ! MODIF ICI: energie? deformation?..
    read(11,*) Tdomain%logicD%save_energy
    if (rg==0) write (*,*) 'save_energy:', Tdomain%logicD%save_energy
    read(11,*) Tdomain%logicD%save_restart
    if (rg==0) write (*,*) 'save_restart:', Tdomain%logicD%save_restart
    read(11,*) Tdomain%logicD%plot_grid
    if (rg==0) write (*,*) 'plot_grid:', Tdomain%logicD%plot_grid
    read(11,*) Tdomain%logicD%run_exec
    if (rg==0) write (*,*) 'run_exec:', Tdomain%logicD%run_exec
    read(11,*) Tdomain%logicD%run_debug
    if (rg==0) write (*,*) 'run_debug:', Tdomain%logicD%run_debug
    read(11,*) Tdomain%logicD%run_echo
    if (rg==0) write (*,*) 'run_echo:', Tdomain%logicD%run_echo
    read(11,*) Tdomain%logicD%run_restart
    if (rg==0) write (*,*) 'run_restart:', Tdomain%logicD%run_restart
    if(rg==0) &
        write (6,*) 'run_restart',Tdomain%logicD%run_restart
    ! numero de l iteration de reprise
    read (11,*) Tdomain%TimeD%iter_reprise
    if (rg==0) write (*,*) 'restart from iter:', Tdomain%TimeD%iter_reprise

    ! creation de fichiers de reprise
    read (11,*) Tdomain%TimeD%ncheck ! frequence de sauvegarde
    if (rg==0) write (*,*) 'restart snap freq:', Tdomain%TimeD%ncheck

    read(11,*) Tdomain%station_file
    if (rg==0) write (*,*) 'Fichier capteurs:', Tdomain%station_file
    read(11,*) Tdomain%TimeD%ntrace
    if (rg==0) write (*,*) 'Freq capteurs:', Tdomain%TimeD%ntrace
    read(11,*) Tdomain%TimeD%time_snapshots
    if (rg==0) write (*,*) 'Freq snapshots:', Tdomain%TimeD%time_snapshots

    logic_scheme = Tdomain%TimeD%acceleration_scheme .neqv. Tdomain%TimeD%velocity_scheme
    if(.not. logic_scheme) then
        stop "Both acceleration and velocity schemes: no compatibility, chose only one."
    end if

    ! Neumann boundary conditions? If yes: geometrical properties read in the mesh files.
    read(11,*) Tdomain%logicD%Neumann
    if (rg==0) write (*,*) 'Neumann BC:', Tdomain%logicD%Neumann
    read(11,*) Tdomain%neumann_file
    if (rg==0) write (*,*) 'Neumann file:', Tdomain%neumann_file


    call read_source_desc(Tdomain, rg)


    call read_attn_desc(Tdomain, rg)

    ! conversion d'un maillage unv en maillage sem
    read(11,*,ERR=101) Tdomain%bMailUnv
    if (rg==0) write (*,*) 'Maillage au format unv:', Tdomain%bMailUnv

    ! prise en compte du fichier des capteurs
    read (11,*,ERR=102) Tdomain%bCapteur
    if (rg==0) write (*,*) 'Utilisation des capteurs:', Tdomain%bCapteur

    !   ajout lecture du flag pour les MPML
    !   (Guillot et Mariotti fevrier 2012)
    read(11,*) Tdomain%logicD%MPML
    if (rg==0) write(*,*) 'Prise en compte MPML:', Tdomain%logicD%MPML
    read(11,*) Tdomain%MPML_coeff
    if (rg==0) write(*,*) 'Coefficient MPML:', Tdomain%MPML_coeff

    !  debut ajout Gradient sur proprietes
    read (11,*) Tdomain%logicD%grad_bassin
    if (rg==0) write(*,*) 'Prise en compte Gradient:', Tdomain%logicD%grad_bassin
    read(11,*) Tdomain%file_bassin
    if (rg==0) write(*,*) 'Fichier desc Gradient:', Tdomain%file_bassin

    if ( Tdomain%logicD%grad_bassin ) then
        call read_gradient_desc(Tdomain, rg)
    endif
    !  fin ajout Gradient sur proprietes

    close (11)

    ! If echo modality write the read parameter in a file
    if (Tdomain%logicD%run_echo) then

        call echo_input_params(Tdomain, rg)

    endif

    !print *,'ap echo',tDomain%bMailUnv  !Gsa ipsis
    ! si le format du maillage initial est UNV,il faut le convertir en maillage sem2D
    if (tDomain%bMailUnv) call convertirUnv(tDomain, rg)
    if (tDomain%bMailUnv) then
        write(6,*) 'End of the conversion unv -> SEM'
        write(6,'(A,1X,A)') 'The name of the mesh file is now:', Tdomain%mesh_file
    endif


    call read_mesh_file(Tdomain, rg)



    !-- Checking mesh properties
    !Tdomain%check_mesh_file = trim(Tdomain%mesh_file) // "_chk"
    !print*,Tdomain%check_mesh_file
    !open(12, file=Tdomain%check_mesh_file, iostat=ok, status="replace", form="formatted",action="write")
    !if(ok /= 0)then
    !    write(*,*) "Process ",rg, " can't open its check mesh_file."
    !    stop
    !end if
    !write(12,*) Tdomain%n_dime
    !write(12,*) Tdomain%logicD%solid_fluid
    !write(12,*) Tdomain%logicD%all_fluid
    !write(12,*) Tdomain%logicD%Neumann
    !write(12,*)  "Global nodes for each proc."
    !write(12,*) Tdomain%n_glob_nodes
    !write(12,*) Tdomain%curve
    !do i = 0,Tdomain%n_glob_nodes-1
    !   write(12,*) (Tdomain%Coord_nodes(j,i), j=0,Tdomain%n_dime-1)
    !enddo
    !write(12,*) "Elements"
    !write(12,*) Tdomain%n_elem
    !write(12,*) "Materials"
    !write(12,*) Tdomain%n_mat
    !do i = 0, Tdomain%n_elem-1
    !   write(12,*) Tdomain%specel(i)%mat_index, Tdomain%specel(i)%solid
    !enddo
    !write(12,*) "Index of nodes for elements"
    !write(12,*) Tdomain%n_nodes
    !do i = 0, Tdomain%n_elem-1
    !    write(12,*) Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1)
    !enddo
    !write(12,*) "Faces and elements properties related to faces"
    !write(12,*) Tdomain%n_face
    !do i = 0, Tdomain%n_elem-1
    !    write(12,*) Tdomain%specel(i)%Near_Faces(0:5)
    !    write(12,*) Tdomain%specel(i)%Orient_Faces(0:5)
    !enddo
    !write(12,*) "Edges"
    !write(12,*) Tdomain%n_edge
    !do i = 0, Tdomain%n_elem-1
    !    write(12,*) Tdomain%specel(i)%Near_Edges(0:11)
    !    write(12,*) Tdomain%specel(i)%Orient_Edges(0:11)
    !enddo
    !write(12,*) "Vertices"
    !write(12,*) Tdomain%n_vertex
    !do i = 0,Tdomain%n_elem-1
    !    write(12,*) Tdomain%specel(i)%Near_Vertices(0:7)
    !enddo
    !write(12,*) "vertex <-> global node"
    !do i = 0,Tdomain%n_vertex-1
    !    write(12,*) Tdomain%sVertex(i)%global_numbering
    !end do

    !if(Tdomain%logicD%solid_fluid)then
    !    write(12,*)
    !    write(12,*) Tdomain%logicD%SF_local_present
    !    if(Tdomain%logicD%SF_local_present)then
    !      write(12,*) "Solid-fluid properties"
    !      write(12,*) "SF faces"
    !      write(12,*) Tdomain%SF%SF_n_faces
    !      write(12,*) "Edges and their orientation, for each SF face"
    !      do i = 0, Tdomain%SF%SF_n_faces-1
    !        write(12,*) Tdomain%SF%SF_face(i)%Near_Edges(0:3)
    !        write(12,*) Tdomain%SF%SF_face(i)%Near_Edges_Orient(0:3)
    !      end do
    !      write(12,*) "Vertices for each SF face"
    !      do i = 0, Tdomain%SF%SF_n_faces-1
    !        write(12,*) Tdomain%SF%SF_face(i)%Near_Vertices(0:3)
    !      end do
    !      write(12,*) "associated fluid (0) and solid (1) faces; orientation Solid/Fluid"
    !      do i = 0, Tdomain%SF%SF_n_faces-1
    !        write(12,*) Tdomain%SF%SF_face(i)%Face(0), Tdomain%SF%SF_face(i)%Face(1),  &
    !                   Tdomain%SF%SF_face(i)%Orient_Face
    !      end do
    !      write(12,*) "SF edges"
    !      write(12,*) Tdomain%SF%SF_n_edges
    !      do i = 0, Tdomain%SF%SF_n_edges-1
    !        write(12,*) Tdomain%SF%SF_edge(i)%Edge(0), Tdomain%SF%SF_edge(i)%Edge(1),  &
    !                   Tdomain%SF%SF_edge(i)%Orient_Edge
    !      end do
    !      write(12,*) "SF vertices"
    !      write(12,*) Tdomain%SF%SF_n_vertices
    !      do i = 0, Tdomain%SF%SF_n_vertices-1
    !        write(12,*) Tdomain%SF%SF_vertex(i)%vertex(0), Tdomain%SF%SF_vertex(i)%vertex(1)
    !      end do
    !    end if
    !end if

    ! Neumann B.C. properties, eventually
    !if(Tdomain%logicD%Neumann)then
    !    write(12,*)
    !    write(12,*) Tdomain%logicD%Neumann_local_present
    !    if(Tdomain%logicD%Neumann_local_present)then
    !      write(12,*) "Neumann properties"
    !      write(12,*) "Neumann faces"
    !      write(12,*) Tdomain%Neumann%Neu_n_faces
    !      write(12,*) "Edges and their orientation, for each Neumann face"
    !      do i = 0, Tdomain%Neumann%Neu_n_faces-1
    !        write(12,*) Tdomain%Neumann%Neu_face(i)%Near_Edges(0:3)
    !        write(12,*) Tdomain%Neumann%Neu_face(i)%Near_Edges_Orient(0:3)
    !      end do
    !      write(12,*) "Vertices for each Neumann face"
    !      do i = 0, Tdomain%Neumann%Neu_n_faces-1
    !        write(12,*) Tdomain%Neumann%Neu_face(i)%Near_Vertices(0:3)
    !      end do
    !      write(12,*) "associated face"
    !      do i = 0, Tdomain%Neumann%Neu_n_faces-1
    !        write(12,*) Tdomain%Neumann%Neu_face(i)%Face
    !      end do
    !      write(12,*) "Neumann edges"
    !      write(12,*) Tdomain%Neumann%Neu_n_edges
    !      do i = 0, Tdomain%Neumann%Neu_n_edges-1
    !        write(12,*) Tdomain%Neumann%Neu_edge(i)%Edge
    !      end do
    !      write(12,*) "Neumann vertices"
    !      write(12,*) Tdomain%Neumann%Neu_n_vertices
    !      do i = 0, Tdomain%Neumann%Neu_n_vertices-1
    !        write(12,*) Tdomain%Neumann%Neu_vertex(i)%vertex
    !      end do
    !    end if
    !end if

    !write(12,*)
    !write(12,*) "Interproc communications"
    !write(12,*) Tdomain%n_proc
    !do i = 0,Tdomain%n_proc-1
    !    write(12,*) Tdomain%sComm(i)%nb_faces, Tdomain%sComm(i)%nb_edges, Tdomain%sComm(i)%nb_vertices
    !    if(Tdomain%logicD%SF_local_present)then
    !      write(12,*) Tdomain%sComm(i)%SF_nf_shared, Tdomain%sComm(i)%SF_ne_shared, Tdomain%sComm(i)%SF_nv_shared
    !    end if
    !    if(Tdomain%logicD%Neumann_local_present)then
    !      write(12,*) Tdomain%sComm(i)%Neu_ne_shared, Tdomain%sComm(i)%Neu_nv_shared
    !    end if
    !    if(Tdomain%sComm(i)%nb_faces > 0)then
    !        do j = 0,Tdomain%sComm(i)%nb_faces-1
    !            write(12,*) Tdomain%sComm(i)%faces(j),Tdomain%sComm(i)%orient_faces(j)
    !        enddo
    !    endif
    !    if(Tdomain%sComm(i)%nb_edges > 0)then
    !        do j = 0,Tdomain%sComm(i)%nb_edges-1
    !            write(12,*) Tdomain%sComm(i)%edges(j),Tdomain%sComm(i)%orient_edges(j)
    !        enddo
    !    endif
    !    if(Tdomain%sComm(i)%nb_vertices > 0)then
    !        do j = 0,Tdomain%sComm(i)%nb_vertices-1
    !            write(12,*) Tdomain%sComm(i)%vertices(j)
    !        enddo
    !    endif
    !    if(Tdomain%logicD%SF_local_present)then
    !      if(Tdomain%sComm(i)%SF_nf_shared > 0)then
    !        do j = 0,Tdomain%sComm(i)%SF_nf_shared-1
    !            write(12,*) Tdomain%sComm(i)%SF_faces_shared(j)
    !        enddo
    !      endif
    !      if(Tdomain%sComm(i)%SF_ne_shared > 0)then
    !        do j = 0,Tdomain%sComm(i)%SF_ne_shared-1
    !            write(12,*) Tdomain%sComm(i)%SF_edges_shared(j),Tdomain%sComm(i)%SF_mapping_edges_shared(j)
    !        enddo
    !      endif
    !      if(Tdomain%sComm(i)%SF_nv_shared > 0)then
    !        do j = 0,Tdomain%sComm(i)%SF_nv_shared-1
    !            write(12,*) Tdomain%sComm(i)%SF_vertices_shared(j)
    !        enddo
    !      endif
    !    end if
    !    if(Tdomain%logicD%Neumann_local_present)then
    !      if(Tdomain%sComm(i)%Neu_ne_shared > 0)then
    !        do j = 0,Tdomain%sComm(i)%Neu_ne_shared-1
    !            write(12,*) Tdomain%sComm(i)%Neu_edges_shared(j),Tdomain%sComm(i)%Neu_mapping_edges_shared(j)
    !        enddo
    !      endif
    !      if(Tdomain%sComm(i)%Neu_nv_shared > 0)then
    !        do j = 0,Tdomain%sComm(i)%Neu_nv_shared-1
    !            write(12,*) Tdomain%sComm(i)%Neu_vertices_shared(j)
    !        enddo
    !      endif
    !    end if

    !end do
    !close(12)


    !---   Properties of materials.
    call read_material_file(Tdomain, rg)

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

    !- receivers'properties
    if(Tdomain%logicD%save_trace)then

        call semname_read_inputmesh_parametrage(Tdomain%station_file,fnamef)
        open (14, file=fnamef, status="old")
        read(14,*) Tdomain%n_receivers
        allocate(Tdomain%sReceiver(0:Tdomain%n_receivers-1))
        do i = 0, Tdomain%n_receivers-1
            read(14,*) Tdomain%sReceiver(i)%xRec, Tdomain%sReceiver(i)%yRec,    &
                Tdomain%sReceiver(i)%zRec, Tdomain%sReceiver(i)%flag,    &
                Tdomain%sReceiver(i)%ndt
        enddo
        close(14)
    endif

    if(Tdomain%logicD%plot_grid .or. Tdomain%logicD%save_snapshots) then
        ! PostProcess Gid
        write(fnamef,"(a,I2.2,a)") "Proc",rg,".flavia.msh"
        open(22, file=fnamef, iostat=ok, status="unknown", form="formatted")
        write(22,*) Tdomain%n_glob_nodes,Tdomain%n_elem
        do i = 0,Tdomain%n_glob_nodes-1
            write(22,"(I4.4,3G12.4)") i+1,(Tdomain%Coord_nodes(j,i), j=0,Tdomain%n_dime-1)
        enddo
        write(22,*)
        do i = 0, Tdomain%n_elem - 1
            write (22,"(I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I4.4,a,I2.2)") i+1,' ',&
                Tdomain%specel(i)%Control_Nodes(0)+1,' ',Tdomain%specel(i)%Control_Nodes(1)+1,' ',&
                Tdomain%specel(i)%Control_Nodes(2)+1,' ',Tdomain%specel(i)%Control_Nodes(3)+1,' ',&
                Tdomain%specel(i)%Control_Nodes(4)+1,' ',Tdomain%specel(i)%Control_Nodes(5)+1,' ',&
                Tdomain%specel(i)%Control_Nodes(6)+1,' ',Tdomain%specel(i)%Control_Nodes(7)+1,' ',&
                Tdomain%specel(i)%mat_index+1
        enddo
        close(22)
    endif

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

    return
101 write(*,*) 'Mauvaise lecture de bmailunv (conversion unv->sem). Arret'
    stop

102 write(*,*) 'Mauvaise lecture de bCapteur (gestion des capteurs). Arret'
    stop

end subroutine read_input
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
