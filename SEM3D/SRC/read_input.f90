!>
!!\file read_input.f90
!!\brief Contient la subroutine read_input().
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Assure la lecture des fichiers de données en entrée.
!!
!! \param type (domain), intent (INOUT) Tdomain
!<


subroutine read_input (Tdomain, rg, code)

    use sdomain
    use semdatafiles
    use mpi
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent(IN) :: rg

    logical :: logic_scheme
    logical, dimension(:), allocatable :: L_Face, L_Edge
    integer :: i, j, npml, n_aus, mat, ok, nf, ne, nv, k, icount, code
#ifndef MKA3D
    integer length
#endif
    real :: dtmin
    character(Len=MAX_FILE_SIZE) :: fnamef
    integer :: unit_src, unit_att
    logical :: trouve

    call semname_read_input_input(fnamef)

    open (11,file=fnamef,form="formatted",status="old")
    read (11,*) Tdomain%Title_simulation
    read (11,*) Tdomain%TimeD%acceleration_scheme
    read (11,*) Tdomain%TimeD%velocity_scheme
    read (11,*) Tdomain%TimeD%duration
    read (11,*) Tdomain%TimeD%alpha
    read (11,*) Tdomain%TimeD%beta
    read (11,*) Tdomain%TimeD%gamma
    read (11,*) Tdomain%mesh_file
    if(rg==0) print*,'MESH_FILE ', Tdomain%mesh_file
    print*,'MESH_FILE processeur sem  ',rg, Tdomain%mesh_file
#ifndef MKA3D
    length = len_trim(Tdomain%mesh_file) + 1
    write (Tdomain%mesh_file(length:length+3),'(i4.4)') rg
#endif
    Tdomain%aniso = .FALSE.

    read (11,*) Tdomain%material_file
    read (11,*) Tdomain%logicD%save_trace
    read (11,*) Tdomain%logicD%save_snapshots
    read (11,*) Tdomain%logicD%save_energy
    read (11,*) Tdomain%logicD%save_restart
    read (11,*) Tdomain%logicD%plot_grid
    read (11,*) Tdomain%logicD%run_exec
    read (11,*) Tdomain%logicD%run_debug
    read (11,*) Tdomain%logicD%run_echo
    read (11,*) Tdomain%logicD%run_restart
    if(rg==0) &
        write (6,*) 'run_restart',Tdomain%logicD%run_restart
    ! numero de l iteration de reprise
    read (11,*) Tdomain%TimeD%iter_reprise

    ! creation de fichiers de reprise
    if (Tdomain%logicD%save_restart) then
        read (11,*) Tdomain%TimeD%ncheck ! frequence de sauvegarde
    else
        read( 11,*)
    endif

    if (Tdomain%logicD%save_trace) then
        read (11,*) Tdomain%station_file
        read (11,*) Tdomain%TimeD%ntrace
        if(rg==0) &
            write (*,*) "Sauvegarde demandee sur processeur ",rg
    else
        read (11,*)
        read (11,*)
    endif

    if (Tdomain%logicD%save_snapshots) then
        read (11,*) Tdomain%TimeD%time_snapshots
    else
        read (11,*)
    endif

    if(.NOT. (Tdomain%TimeD%acceleration_scheme .AND. Tdomain%TimeD%velocity_scheme)) then
        logic_scheme = Tdomain%TimeD%acceleration_scheme .or. Tdomain%TimeD%velocity_scheme
    else
        logic_scheme = .false.
    endif

    if (.not. logic_scheme) then
        write (*,*) "No compatible acceleration and velocity schemes"
        stop
    endif

    read (11,*) Tdomain%logicD%super_object
    if (Tdomain%logicD%super_object) then
        read (11,*) Tdomain%Super_object_type, Tdomain%super_object_file
    else
        read (11,*)
    endif
    read (11,*) Tdomain%logicD%Neumann
    if ( Tdomain%logicD%Neumann ) then
        read (11,*) Tdomain%neumann_file
    else
        read (11,*)
    endif


#ifdef MKA3D
    unit_src=21
    unit_att=22
    Tdomain%logicD%any_source=.true.
    call semname_read_input_source(fnamef)
    inquire(file=fnamef, exist=trouve)
    if( .not. trouve) then
        Tdomain%logicD%any_source=.false.
        if(rg==0) &
            write(*,*) 'No source for Sem'
    else
        open (unit=unit_src,file=fnamef,form="formatted",status="old")
    endif

    call semname_read_input_amortissement(fnamef)
    inquire(file=fnamef, exist=trouve)
    if( .not. trouve) then
        Tdomain%n_sls=0
        if(rg==0) &
            write(*,*) 'No attenuation model'
    else
        Tdomain%n_sls=1
        open (unit=unit_att,file=fnamef,form="formatted",status="old")
    endif
#else
    unit_src=11
    unit_att=11
    read(unit_src,*) Tdomain%logicD%any_source
#endif

    if (Tdomain%logicD%any_source) then
        read (unit_src,*) Tdomain%n_source

        if(rg==0) &
            write (*,*) 'nombre de sources ', Tdomain%n_source
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
    endif

#ifdef MKA3D
    close(unit_src)
#endif


    !  modif mariotti fevrier 2007 cea
    !   n_sls est le nombre de sous decoupage en pas de frequence pour avoir un
    !   facteur de qualite a peu pres constant dans l intervalle
    !   voir article de Liu Kanamori du GJRA de 1976
#ifdef MKA3D
    if (Tdomain%n_sls>0) then
        read (unit_att,*) Tdomain%n_sls
#else
        read (unit_att,*) Tdomain%n_sls
#endif

        if (Tdomain%n_sls>0) then
            read (unit_att,*) Tdomain%T1_att, Tdomain%T2_att
            read (unit_att,*) Tdomain%T0_modele
            if(rg==0) &
                write(*,*) ' modele attenuation sur temps ', Tdomain%T1_att, Tdomain%T2_att,Tdomain%T0_modele
        else
            if(rg==0) &
                write(*,*) 'pas de  modele attenuation sur temps '
        endif

#ifdef MKA3D
    endif
    close(unit_att)
#endif


    ! conversion d'un maillage unv en maillage sem
    read(11,*,ERR=101) Tdomain%bMailUnv

    ! prise en compte du fichier des capteurs
    read (11,*,ERR=102) Tdomain%bCapteur
    !   ajout lecture du flag pour les MPML
    !   (Guillot et Mariotti fevrier 2012)
    read (11,*) Tdomain%logicD%MPML
    if ( Tdomain%logicD%MPML ) then
        read(11,*) Tdomain%MPML_coeff
        if (rg==0) print*,' prise en compte des MPML avec coefficient de ',Tdomain%MPML_coeff
    endif

    !  debut ajout Gradient sur proprietes
    read (11,*) Tdomain%logicD%grad_bassin
    if ( Tdomain%logicD%grad_bassin ) then
        read(11,*) Tdomain%file_bassin
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
    endif
    !  fin ajout Gradient sur proprietes

    close (11)

    ! If echo modality write the read parameter in a file
    if (Tdomain%logicD%run_echo) then

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

            if (Tdomain%logicD%super_object) then
                write (91,*) Tdomain%Super_object_type,"   ", Tdomain%super_object_file
            else
                write (91,*) "no super objects present"
            endif

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
    endif

    !print *,'ap echo',tDomain%bMailUnv  !Gsa ipsis
    ! si le format du maillage initial est UNV,il faut le convertir en maillage sem2D
    if (tDomain%bMailUnv) call convertirUnv(tDomain, rg)
    if (tDomain%bMailUnv) then
        write(6,*) 'End of the conversion unv -> SEM'
        write(6,'(A,1X,A)') 'The name of the mesh file is now:', Tdomain%mesh_file
    endif

    ! Read Mesh properties
#ifdef MKA3D
    call semname_read_input_meshfile(rg,Tdomain%mesh_file,fnamef)
    Tdomain%mesh_file = fnamef
    if(rg==0) print*,'MESH_FILE', Tdomain%mesh_file
    print*,'rg  MESH_FILE',rg, Tdomain%mesh_file
#endif
    open (12, file=Tdomain%mesh_file, iostat=ok, status="old", form="formatted")
    print*,'ok rg  MESH_FILE',ok, rg, Tdomain%mesh_file
    if (ok/=0) then
        write (*,*) "Process ",rg, " can't open his mesh_file ", Tdomain%mesh_file
        stop
        call mpi_finalize(code)
    endif
    read (12,*) Tdomain%n_dime
    read (12,*) Tdomain%n_glob_nodes
    write(6,*) 'nb noeuds',Tdomain%n_glob_nodes,'proc ',rg
    read (12,*) Tdomain%curve
    allocate (Tdomain%Coord_nodes(0:Tdomain%n_dime-1,0:Tdomain%n_glob_nodes-1))
    do i = 0,Tdomain%n_glob_nodes-1
        read (12,*) (Tdomain%Coord_nodes(j,i), j=0,Tdomain%n_dime-1)
    enddo
    read (12,*) Tdomain%n_elem
    write(6,*) 'nb elts',Tdomain%n_elem,'proc ',rg
    allocate (Tdomain%specel(0:Tdomain%n_elem-1))
    do i=0,Tdomain%n_elem-1
        Tdomain%specel(i)%PML = .FALSE.
    enddo
    read (12,*) Tdomain%n_mat

    do i = 0, Tdomain%n_elem - 1
        read(12,*) Tdomain%specel(i)%mat_index
    enddo
    if (Tdomain%n_dime == 3) then
        read (12,*) Tdomain%n_nodes
        do i = 0, Tdomain%n_elem - 1
            allocate (Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1))
            read(12,*) Tdomain%specel(i)%Control_Nodes(0:Tdomain%n_nodes-1)
        enddo
        read (12,*) Tdomain%n_face
        allocate (Tdomain%sFace(0:Tdomain%n_face-1))
        do i=0,Tdomain%n_face-1
            Tdomain%sFace(i)%PML = .FALSE.
        enddo
        do i = 0, Tdomain%n_elem - 1
            read(12,*) Tdomain%specel(i)%Near_Faces(0:5)
            read(12,*) Tdomain%specel(i)%Orient_Faces(0:5)
        enddo
        read (12,*) Tdomain%n_edge
        allocate (Tdomain%sEdge(0:Tdomain%n_edge-1))
        do i=0,Tdomain%n_edge-1
            Tdomain%sEdge(i)%PML = .FALSE.
        enddo
        do i = 0, Tdomain%n_elem - 1
            read(12,*) Tdomain%specel(i)%Near_Edges(0:11)
            read(12,*) Tdomain%specel(i)%Orient_Edges(0:11)
        enddo
        read (12,*) Tdomain%n_vertex
        allocate (Tdomain%sVertex(0:Tdomain%n_vertex-1))
        do i=0,Tdomain%n_vertex-1
            Tdomain%sVertex(i)%PML = .FALSE.
        enddo
        do i = 0, Tdomain%n_elem - 1
            read(12,*) Tdomain%specel(i)%Near_Vertices(0:7)
        enddo
        read(12,*)
        ! Information about super-object
        read(12,*) ! Super Object
        read(12,*) Tdomain%logicD%super_object_local_present
        if (Tdomain%logicD%super_object) then
            if ( Tdomain%logicD%super_object_local_present ) then
                read(12,*) ! Faces in Super Object
                read(12,*) Tdomain%sPlaneW%n_faces
                allocate (Tdomain%sPlaneW%pFace(0:Tdomain%sPlaneW%n_faces-1))
                read(12,*) ! 4 Edges for each face
                do nf = 0, Tdomain%sPlaneW%n_faces-1
                    read(12,*) Tdomain%sPlaneW%pFace(nf)%Near_Edges(0:3)
                    read(12,*) Tdomain%sPlaneW%pFace(nf)%Orient_Edges(0:3)
                enddo
                read(12,*) ! 4 Vertices for each face
                do nf = 0, Tdomain%sPlaneW%n_faces-1
                    read(12,*) Tdomain%sPlaneW%pFace(nf)%Near_Vertices(0:3)
                enddo
                read(12,*) ! Glob number of up and down faces and orientation of down compared to up
                read(12,*)
                do nf = 0, Tdomain%sPlaneW%n_faces-1
                    read(12,*) Tdomain%sPlaneW%pFace(nf)%Face_UP, Tdomain%sPlaneW%pFace(nf)%Face_DOWN, Tdomain%sPlaneW%pFace(nf)%Orient
                enddo
                read(12,*) ! Glob number of up and down edges and orientation of down compared to up
                read(12,*) Tdomain%sPlaneW%n_edges
                allocate (Tdomain%sPlaneW%pEdge(0:Tdomain%sPlaneW%n_edges-1))
                do ne = 0, Tdomain%sPlaneW%n_edges-1
                    read(12,*) Tdomain%sPlaneW%pEdge(ne)%Edge_UP, Tdomain%sPlaneW%pEdge(ne)%Edge_DOWN, Tdomain%sPlaneW%pEdge(ne)%Orient
                enddo
                read(12,*) ! Glob number of up and down vertices and orientation of down compared to up
                read(12,*) Tdomain%sPlaneW%n_vertices
                allocate (Tdomain%sPlaneW%pVertex(0:Tdomain%sPlaneW%n_vertices-1))
                do nv = 0, Tdomain%sPlaneW%n_vertices-1
                    read(12,*) Tdomain%sPlaneW%pVertex(nv)%Vertex_UP, Tdomain%sPlaneW%pVertex(nv)%Vertex_DOWN
                enddo
            endif
        endif
        read(12,*)
        read(12,*) ! Neumann
        read(12,*) Tdomain%logicD%Neumann_local_present
        if ( Tdomain%logicD%Neumann_local_present ) then
            read(12,*) ! Faces in Neumann
            read(12,*) Tdomain%sNeu%n_faces
            allocate (Tdomain%sNeu%nFace(0:Tdomain%sNeu%n_faces-1))
            read(12,*) ! 4 Edges for each face
            do nf = 0, Tdomain%sNeu%n_faces-1
                read(12,*) Tdomain%sNeu%nFace(nf)%Near_Edges(0:3)
                read(12,*) Tdomain%sNeu%nFace(nf)%Orient_Edges(0:3)
            enddo
            read(12,*) ! 4 Vertices for each face
            do nf = 0, Tdomain%sNeu%n_faces-1
                read(12,*) Tdomain%sNeu%nFace(nf)%Near_Vertices(0:3)
            enddo
            read(12,*) ! Glob number of faces
            read(12,*)
            do nf = 0, Tdomain%sNeu%n_faces-1
                read(12,*) Tdomain%sNeu%nFace(nf)%Face
            enddo
            read(12,*) ! Glob number of edges
            read(12,*) Tdomain%sNeu%n_edges
            allocate (Tdomain%sNeu%nEdge(0:Tdomain%sNeu%n_edges-1))
            do ne = 0, Tdomain%sNeu%n_edges-1
                read(12,*) Tdomain%sNeu%nEdge(ne)%Edge
            enddo
            read(12,*) ! Glob number vertices
            read(12,*) Tdomain%sNeu%n_vertices
            allocate (Tdomain%sNeu%nVertex(0:Tdomain%sNeu%n_vertices-1))
            do nv = 0, Tdomain%sNeu%n_vertices-1
                read(12,*) Tdomain%sNeu%nVertex(nv)%Vertex
            enddo
        endif
        read(12,*) ! Communication between procs
        read (12,*) Tdomain%n_proc
        allocate (Tdomain%sComm(0:Tdomain%n_proc-1))
        do i = 0,Tdomain%n_proc-1
            read(12,*) Tdomain%sComm(i)%nb_faces, Tdomain%sComm(i)%nb_edges, Tdomain%sComm(i)%nb_vertices, &
                Tdomain%sComm(i)%nb_edges_so, Tdomain%sComm(i)%nb_vertices_so, Tdomain%sComm(i)%nb_edges_neu, Tdomain%sComm(i)%nb_vertices_neu
            if (Tdomain%sComm(i)%nb_faces>0) then
                allocate (Tdomain%sComm(i)%faces(0:Tdomain%sComm(i)%nb_faces-1))
                allocate (Tdomain%sComm(i)%orient_faces(0:Tdomain%sComm(i)%nb_faces-1))
                do j = 0,Tdomain%sComm(i)%nb_faces-1
                    read(12,*) Tdomain%sComm(i)%faces(j),Tdomain%sComm(i)%orient_faces(j)
                enddo
            endif
            if (Tdomain%sComm(i)%nb_edges>0) then
                allocate (Tdomain%sComm(i)%edges(0:Tdomain%sComm(i)%nb_edges-1))
                allocate (Tdomain%sComm(i)%orient_edges(0:Tdomain%sComm(i)%nb_edges-1))
                do j = 0,Tdomain%sComm(i)%nb_edges-1
                    read(12,*) Tdomain%sComm(i)%edges(j),Tdomain%sComm(i)%orient_edges(j)
                enddo
            endif
            if (Tdomain%sComm(i)%nb_vertices>0) then
                allocate (Tdomain%sComm(i)%vertices(0:Tdomain%sComm(i)%nb_vertices-1))
                do j = 0,Tdomain%sComm(i)%nb_vertices-1
                    read(12,*) Tdomain%sComm(i)%vertices(j)
                enddo
            endif
            if (Tdomain%sComm(i)%nb_edges_so>0) then
                allocate (Tdomain%sComm(i)%edges_SO(0:Tdomain%sComm(i)%nb_edges_so-1))
                allocate (Tdomain%sComm(i)%orient_edges_SO(0:Tdomain%sComm(i)%nb_edges_so-1))
                do j = 0,Tdomain%sComm(i)%nb_edges_so-1
                    read(12,*) Tdomain%sComm(i)%edges_SO(j),Tdomain%sComm(i)%orient_edges_SO(j)
                enddo
            endif
            if (Tdomain%sComm(i)%nb_vertices_so>0) then
                allocate (Tdomain%sComm(i)%vertices_SO(0:Tdomain%sComm(i)%nb_vertices_so-1))
                do j = 0,Tdomain%sComm(i)%nb_vertices_so-1
                    read(12,*) Tdomain%sComm(i)%vertices_SO(j)
                enddo
            endif
            if (Tdomain%sComm(i)%nb_edges_neu>0) then
                allocate (Tdomain%sComm(i)%edges_Neu(0:Tdomain%sComm(i)%nb_edges_neu-1))
                allocate (Tdomain%sComm(i)%orient_edges_Neu(0:Tdomain%sComm(i)%nb_edges_neu-1))
                do j = 0,Tdomain%sComm(i)%nb_edges_neu-1
                    read(12,*) Tdomain%sComm(i)%edges_Neu(j),Tdomain%sComm(i)%orient_edges_Neu(j)
                enddo
            endif
            if (Tdomain%sComm(i)%nb_vertices_neu>0) then
                allocate (Tdomain%sComm(i)%vertices_Neu(0:Tdomain%sComm(i)%nb_vertices_neu-1))
                do j = 0,Tdomain%sComm(i)%nb_vertices_neu-1
                    read(12,*) Tdomain%sComm(i)%vertices_Neu(j)
                enddo
            endif
        enddo
        read(12,*)
        read(12,*) ! Glob number vertices
        read(12,*) Tdomain%logicD%Save_Surface
        if ( Tdomain%logicD%Save_Surface ) then
            read(12,*) Tdomain%sSurf%n_vertices
            allocate (Tdomain%sSurf%nVertex(0:Tdomain%sSurf%n_vertices-1))
            do nv = 0, Tdomain%sSurf%n_vertices-1
                read(12,*) Tdomain%sSurf%nVertex(nv)%Vertex
            enddo
        endif
    else
        write (*,*) "A dimension different from 3 is not yet taken into account"
        stop
    endif
    close (12)


    npml = 0
    allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))

    call semname_read_inputmesh_parametrage(Tdomain%material_file,fnamef)
    open (13, file=fnamef, status="old", form="formatted")

    read (13,*) n_aus

    if (n_aus /= Tdomain%n_mat) then
        write (*,*) "Incompatibility between the mesh file and the material file for n_mat "
        stop
    endif


    if (Tdomain%aniso) then
        print *,"The code can't put anisotropy in a homogeneous media"
        stop
    endif
    do i = 0,Tdomain%n_mat-1
        read (13,*) Tdomain%sSubDomain(i)%material_type, Tdomain%sSubDomain(i)%Pspeed, &
            Tdomain%sSubDomain(i)%Sspeed, Tdomain%sSubDomain(i)%dDensity, &
            Tdomain%sSubDomain(i)%NGLLx, Tdomain%sSubDomain(i)%NGLLy, Tdomain%sSubDomain(i)%NGLLz, &
            Tdomain%sSubDomain(i)%Dt, Tdomain%sSubDomain(i)%Qpression,  Tdomain%sSubDomain(i)%Qmu

        if(rg==0) &
            write (*,*) Tdomain%sSubDomain(i)%material_type, Tdomain%sSubDomain(i)%Pspeed, &
            Tdomain%sSubDomain(i)%Sspeed, Tdomain%sSubDomain(i)%dDensity, &
            Tdomain%sSubDomain(i)%NGLLx, Tdomain%sSubDomain(i)%NGLLy, Tdomain%sSubDomain(i)%NGLLz, &
            Tdomain%sSubDomain(i)%Dt, Tdomain%sSubDomain(i)%Qpression,  Tdomain%sSubDomain(i)%Qmu

        call Lame_coefficients (Tdomain%sSubDomain(i))
        if(rg==0) &
            print*,' lame ',Tdomain%sSubDomain(i)%DMu,Tdomain%sSubDomain(i)%DLambda ,Tdomain%sSubDomain(i)%DKappa
        if (Tdomain%sSubDomain(i)%material_type == "P")  then
            Tdomain%sSubDomain(i)%wpml = npml
            npml = npml + 1
        endif
    enddo

    Tdomain%any_PML = .false.
    Tdomain%any_FPML = .false.
    if (npml > 0) then
        Tdomain%any_PML = .true.
        read(13,*); read(13,*)
        do i = 0,Tdomain%n_mat-1
            if (Tdomain%sSubdomain(i)%material_type == "P") then
                read (13,*) Tdomain%sSubdomain(i)%Filtering, Tdomain%sSubdomain(i)%npow, &
                    Tdomain%sSubdomain(i)%Apow, Tdomain%sSubdomain(i)%Px, &
                    Tdomain%sSubdomain(i)%Left, Tdomain%sSubdomain(i)%Py, &
                    Tdomain%sSubdomain(i)%Forward, Tdomain%sSubdomain(i)%Pz, &
                    Tdomain%sSubdomain(i)%Down, Tdomain%sSubdomain(i)%freq
                if (Tdomain%sSubdomain(i)%Filtering) Tdomain%any_FPML = .true.
            endif
        enddo
    endif
    close (13)


    allocate (L_Face(0:Tdomain%n_face-1))
    L_Face = .true.
    allocate (L_Edge(0:Tdomain%n_edge-1))
    L_Edge = .true.
    do i = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(i)%mat_index
        Tdomain%specel(i)%ngllx = Tdomain%sSubDomain(mat)%NGLLx
        Tdomain%specel(i)%nglly = Tdomain%sSubDomain(mat)%NGLLy
        Tdomain%specel(i)%ngllz = Tdomain%sSubDomain(mat)%NGLLz
        do j = 0,5
            nf = Tdomain%specel(i)%Near_Faces(j)
            if (L_Face(nf) .and. Tdomain%specel(i)%Orient_Faces(j)==0) then
                L_Face(nf) = .false. !L_Face est pour eviter que l'on modifie plusieurs fois la meme face
                if (j==0 .or. j==5) then
                    Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%ngllx
                    Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%nglly
                    Tdomain%sFace(nf)%dir = j
                else if (j==1 .or. j==3) then
                    Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%ngllx
                    Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%ngllz
                    Tdomain%sFace(nf)%dir = j
                else
                    Tdomain%sFace(nf)%ngll1 = Tdomain%specel(i)%nglly
                    Tdomain%sFace(nf)%ngll2 = Tdomain%specel(i)%ngllz
                    Tdomain%sFace(nf)%dir = j
                endif
            endif
        enddo
        do j = 0,11
            ne = Tdomain%specel(i)%Near_Edges(j)
            if (L_Edge(ne) .and. Tdomain%specel(i)%Orient_Edges(j)==0) then
                L_Edge(ne) = .false.
                if (j==0 .or. j==2 .or. j==5 .or. j==9) then
                    Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%ngllx
                else if (j==1 .or. j==3 .or. j==8 .or. j==11) then
                    Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%nglly
                else
                    Tdomain%sEdge(ne)%ngll = Tdomain%specel(i)%ngllz
                endif
            endif
        enddo
    enddo
    deallocate (L_Face,L_Edge)

    if (Tdomain%logicD%super_object_local_present) then
        do nf = 0, Tdomain%sPlaneW%n_faces-1
            n_aus = Tdomain%sPlaneW%pFace(nf)%Face_UP
            Tdomain%sPlaneW%pFace(nf)%ngll1 = Tdomain%sFace(n_aus)%ngll1
            Tdomain%sPlaneW%pFace(nf)%ngll2 = Tdomain%sFace(n_aus)%ngll2
            Tdomain%sPlaneW%pFace(nf)%dir = Tdomain%sFace(n_aus)%dir
        enddo
        do ne = 0, Tdomain%sPlaneW%n_edges-1
            n_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_UP
            Tdomain%sPlaneW%pEdge(ne)%ngll = Tdomain%sEdge(n_aus)%ngll
        enddo
    endif
    if (Tdomain%logicD%neumann_local_present) then
        do nf = 0, Tdomain%sNeu%n_faces-1
            n_aus = Tdomain%sNeu%nFace(nf)%Face
            Tdomain%sNeu%nFace(nf)%ngll1 = Tdomain%sFace(n_aus)%ngll1
            Tdomain%sNeu%nFace(nf)%ngll2 = Tdomain%sFace(n_aus)%ngll2
            Tdomain%sNeu%nFace(nf)%dir = Tdomain%sFace(n_aus)%dir
        enddo
        do ne = 0, Tdomain%sNeu%n_edges-1
            n_aus = Tdomain%sNeu%nEdge(ne)%Edge
            Tdomain%sNeu%nEdge(ne)%ngll = Tdomain%sEdge(n_aus)%ngll
        enddo
    endif

    do i = 0, Tdomain%n_mat-1
        if (Tdomain%sSubdomain(i)%material_type == "P" .and. Tdomain%sSubdomain(i)%Filtering ) &
            Tdomain%sSubdomain(i)%freq = exp (-Tdomain%sSubdomain(i)%freq*Tdomain%sSubdomain(i)%dt/2)
    enddo

    dtmin = 1e20
    do i = 0,Tdomain%n_mat-1
        if (Tdomain%sSubDomain(i)%Dt < dtmin) dtmin = Tdomain%sSubDomain(i)%Dt
    enddo
    Tdomain%TimeD%dtmin = dtmin
    if (dtmin > 0) then
        Tdomain%TimeD%ntimeMax = int (Tdomain%TimeD%Duration/dtmin)
    else
        write (*,*) "Your dt min is zero : verify it"
        stop
    endif
    if(rg==0) &
        print *,'ntimemax',Tdomain%TimeD%ntimeMax,Tdomain%TimeD%Duration,dtmin


    if (Tdomain%logicD%save_trace) then

        call semname_read_inputmesh_parametrage(Tdomain%station_file,fnamef)
        open (14, file=fnamef, status="old")
        read (14,*) Tdomain%n_receivers
        allocate (Tdomain%sReceiver(0:Tdomain%n_receivers-1))
        do i = 0, Tdomain%n_receivers-1
            read(14,*) Tdomain%sReceiver(i)%xRec, Tdomain%sReceiver(i)%yRec, Tdomain%sReceiver(i)%zRec
        enddo
        close (14)
    endif


    do i = 0,Tdomain%n_face-1
        do j = 0, Tdomain%n_elem - 1
            do k=0,5
                if ( Tdomain%specel(j)%Near_Faces(k) == i )  then
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

    return

101 write(*,*) 'Mauvaise lecture de bmailunv (conversion unv->sem). Arret'
    stop

102 write(*,*) 'Mauvaise lecture de bCapteur (gestion des capteurs). Arret'
    stop

end subroutine read_Input
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
