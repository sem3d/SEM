subroutine read_input(Tdomain,rg)

    use sdomain
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: rg

    logical :: logic_scheme, sortie, neumann_log
    logical, dimension(:), allocatable :: L_Face, L_Edge
    integer :: length,i,j,npml,n_aus,mat,ok,nf,ne,nv,k,icount,n,i_aus,  &
               ipoint,code,nnf,nne,nnv
    real :: dtmin, x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7
    character*20 :: fnamef


open(11,file="input.spec",form="formatted",status="old")
read(11,*) Tdomain%Title_simulation
read(11,*) Tdomain%TimeD%acceleration_scheme
read(11,*) Tdomain%TimeD%velocity_scheme
read(11,*) Tdomain%TimeD%duration
read(11,*) Tdomain%TimeD%alpha
read(11,*) Tdomain%TimeD%beta
read(11,*) Tdomain%TimeD%gamma
read(11,*) Tdomain%mesh_file
length = len_trim(Tdomain%mesh_file)+1
write(Tdomain%mesh_file(length:length+2),'(i3.3)') rg
read(11,*) Tdomain%material_file
read(11,*) Tdomain%logicD%save_trace
read(11,*) Tdomain%logicD%save_snapshots
! MODIF ICI: energie? deformation?..
read(11,*) Tdomain%logicD%save_energy
read(11,*) Tdomain%logicD%save_restart
read(11,*) Tdomain%logicD%plot_grid
read(11,*) Tdomain%logicD%run_exec
read(11,*) Tdomain%logicD%run_debug
read(11,*) Tdomain%logicD%run_echo
read(11,*) Tdomain%logicD%run_restart

if(Tdomain%logicD%save_restart) then
    read(11,*) Tdomain%TimeD%ncheck
else 
    read(11,*)
    if(rg == 0)then
       write(*,*) "Sure that you wish a run without any backup? (Press enter)" ; read*
    end if
    call MPI_BARRIER(MPI_COMM_WORLD,code)
endif

if(Tdomain%logicD%save_trace)then
    read(11,*) Tdomain%station_file
    read(11,*) Tdomain%TimeD%ntrace
else 
    read(11,*)
    read(11,*)
endif

if(Tdomain%logicD%save_snapshots)then
    read(11,*) Tdomain%TimeD%time_snapshots
else 
    read(11,*)
endif

logic_scheme = Tdomain%TimeD%acceleration_scheme .xor. Tdomain%TimeD%velocity_scheme
if(.not. logic_scheme) then
    stop "Both acceleration and velocity schemes: no compatibility, chose only one."
end if

! Neumann boundary conditions? If yes: geometrical properties read in the mesh files.
read(11,*) Tdomain%logicD%Neumann
if(Tdomain%logicD%Neumann)then
    read(11,*) Tdomain%neumann_file
else
    read(11,*)
end if

read(11,*) Tdomain%logicD%any_source
if(Tdomain%logicD%any_source) then
    read(11,*) Tdomain%n_source 
    allocate(Tdomain%sSource(0:Tdomain%n_source-1))
    do i = 0, Tdomain%n_source-1
        read(11,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Ysource,   &
                   Tdomain%Ssource(i)%Zsource
        read(11,*) Tdomain%Ssource(i)%i_type_source
   ! type_source : 1 = collocated pulse, 2 = moment tensor , 3 = pulse in pressure (fluid)
        if(Tdomain%sSource(i)%i_type_source == 1) then   ! pulse in one space direction
            read(11,*) Tdomain%sSource(i)%i_dir
            read(11,*)
        else if(Tdomain%sSource(i)%i_type_source == 2)then
            read(11,*) Tdomain%sSource(i)%Moment(0,0), Tdomain%sSource(i)%Moment(1,1), Tdomain%sSource(i)%Moment(2,2)
            read(11,*) Tdomain%sSource(i)%Moment(0,1), Tdomain%sSource(i)%Moment(0,2), Tdomain%sSource(i)%Moment(1,2)
            Tdomain%sSource(i)%Moment(1,0) = Tdomain%sSource(i)%Moment(0,1)
            Tdomain%sSource(i)%Moment(2,0) = Tdomain%sSource(i)%Moment(0,2)
            Tdomain%sSource(i)%Moment(2,1) = Tdomain%sSource(i)%Moment(1,2)
        endif
        read(11,*) Tdomain%Ssource(i)%i_time_function
        read(11,*) Tdomain%Ssource(i)%tau_b
        read(11,*) Tdomain%Ssource(i)%cutoff_freq
    end do
end if
! MPML?
read(11,*) Tdomain%logicD%MPML
if(Tdomain%logicD%MPML) read(11,*) Tdomain%MPML_coeff

close(11)

! If echo modality, write the read parameter in a file
if(Tdomain%logicD%run_echo) then
    open(91,file="input_spec_echo", form="formatted", status="unknown")
    write(91,*) Tdomain%Title_simulation
    write(91,*) Tdomain%TimeD%acceleration_scheme
    write(91,*) Tdomain%TimeD%velocity_scheme
    write(91,*) Tdomain%TimeD%duration
    write(91,*) Tdomain%TimeD%alpha
    write(91,*) Tdomain%TimeD%beta
    write(91,*) Tdomain%TimeD%gamma
    write(91,*) Tdomain%mesh_file
    write(91,*) Tdomain%material_file
    write(91,*) Tdomain%logicD%save_trace
    write(91,*) Tdomain%logicD%save_snapshots
    write(91,*) Tdomain%logicD%save_energy
    write(91,*) Tdomain%logicD%plot_grid
    write(91,*) Tdomain%logicD%run_exec
    write(91,*) Tdomain%logicD%run_debug
    write(91,*) Tdomain%logicD%run_echo

    if(Tdomain%logicD%save_trace) then
        write(91,*) Tdomain%station_file
    else 
        write(91,*) " No parameter given here"
    endif

    if(Tdomain%logicD%save_snapshots) then
        write(91,*) Tdomain%TimeD%nsnap
    else 
        write(91,*) " No parameter given here"
    endif

    write(11,*) Tdomain%logicD%Neumann, "  Neumann B.C.?"

    if(Tdomain%logicD%any_source) then
        write(91,*) Tdomain%n_source 
        do i = 0, Tdomain%n_source-1
            write(91,*) Tdomain%Ssource(i)%Xsource, Tdomain%Ssource(i)%Ysource, Tdomain%Ssource(i)%Zsource
            write(91,*) Tdomain%Ssource(i)%i_type_source
            if(Tdomain%Ssource(i)%i_type_source == 1) then
                write(91,*) Tdomain%Ssource(i)%i_dir
            else
                write(91,*) "No parameter need here"
            endif
            write(91,*) Tdomain%Ssource(i)%i_time_function
            write(91,*) Tdomain%Ssource(i)%tau_b
            write(91,*) Tdomain%Ssource(i)%cutoff_freq
        enddo
    else
        write(91,*)  "No available sources "
    end if
    write(91,*) "All right, runner?"
    close(91)
end if


!-- Reading mesh properties
open(12, file=Tdomain%mesh_file, iostat=ok, status="old", form="formatted")
if(ok /= 0)then
    write(*,*) "Process ",rg, " can't open its mesh_file."
    stop
end if
read(12,*) Tdomain%n_dime
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
call MPI_BARRIER(MPI_COMM_WORLD, code)
read(12,*) neumann_log
if(neumann_log .neqv. Tdomain%logicD%Neumann)  &
    stop "Introduction of Neumann B.C.: mesh and input files not in coincidence." 
read(12,*)   ! Global nodes for each proc.
read(12,*) Tdomain%n_glob_nodes
read(12,*) Tdomain%curve
allocate(Tdomain%Coord_nodes(0:Tdomain%n_dime-1,0:Tdomain%n_glob_nodes-1))
do i = 0,Tdomain%n_glob_nodes-1
   read(12,*) (Tdomain%Coord_nodes(j,i), j=0,Tdomain%n_dime-1)
enddo
read(12,*)  ! Elements
read(12,*) Tdomain%n_elem
allocate(Tdomain%specel(0:Tdomain%n_elem-1))
read(12,*) ! Materials
read(12,*) Tdomain%n_mat
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
do i = 0, Tdomain%n_elem-1
    read(12,*) Tdomain%specel(i)%Near_Faces(0:5)
    read(12,*) Tdomain%specel(i)%Orient_Faces(0:5)
enddo
read(12,*)  ! Edges
read(12,*) Tdomain%n_edge
allocate(Tdomain%sEdge(0:Tdomain%n_edge-1))
do i = 0, Tdomain%n_elem-1
    read(12,*) Tdomain%specel(i)%Near_Edges(0:11)
    read(12,*) Tdomain%specel(i)%Orient_Edges(0:11)
enddo
read(12,*)  ! Vertices
read(12,*) Tdomain%n_vertex
allocate(Tdomain%sVertex(0:Tdomain%n_vertex-1))
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
npml = 0
allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))
open(13,file=Tdomain%material_file,status="old",action="read")
read(13,*) n_aus

if(n_aus /= Tdomain%n_mat)   &
    stop "Incompatibility between the mesh file and the material file for n_mat"

do i = 0,Tdomain%n_mat-1
    read(13,*) Tdomain%sSubDomain(i)%material_type, Tdomain%sSubDomain(i)%Pspeed, &
               Tdomain%sSubDomain(i)%Sspeed, Tdomain%sSubDomain(i)%dDensity,      &
               Tdomain%sSubDomain(i)%Dt,Tdomain%sSubDomain(i)%NGLLx,              &
               Tdomain%sSubDomain(i)%NGLLy, Tdomain%sSubDomain(i)%NGLLz
    if(Tdomain%sSubDomain(i)%material_type == "P" .or. Tdomain%sSubDomain(i)%material_type == "L")  then
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
            L_Face(nf) = .false.
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

! Lame's coefficient in case of isotropic solid or fluid
do i = 0, Tdomain%n_mat-1
    call Lame_coefficients(Tdomain%sSubDomain(i))
enddo

!- receivers'properties
if(Tdomain%logicD%save_trace)then
    open(14,file=Tdomain%station_file,status="old",action="read")
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
   write(22,"(I4.4,3f)") i+1,(Tdomain%Coord_nodes(j,i), j=0,Tdomain%n_dime-1)
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
   allocate(Tdomain%sEdge(i)%Which_Elem(0:10))
   allocate(Tdomain%sEdge(i)%Which_EdgeinElem(0:10))
   Tdomain%sEdge(i)%Which_Elem = -1
   Tdomain%sEdge(i)%Which_EdgeinElem = -1
   icount = 0
   do j = 0, Tdomain%n_elem - 1
      do k=0,11
        if(Tdomain%specel(j)%Near_Edges(k) == i)then
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
        Tdomain%sFace(nnf)%mat_index = mat
        Tdomain%sFace(nnf)%solid = Tdomain%specel(n)%solid
    end do
    do ne = 0,11
        nne = Tdomain%specel(n)%Near_edges(ne)
        Tdomain%sEdge(nne)%mat_index = mat
        Tdomain%sEdge(nne)%solid = Tdomain%specel(n)%solid
    end do
    do nv = 0,7
        nnv = Tdomain%specel(n)%Near_Vertices(nv)
        Tdomain%sVertex(nnv)%mat_index = mat
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
end subroutine read_input
