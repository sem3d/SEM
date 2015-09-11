!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
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
    use mesh3d
    use iso_c_binding
    use sem_c_config

    implicit none

    public :: read_input

contains


    subroutine read_gradient_desc(Tdomain)
        use sdomain
        use semdatafiles
        use mpi
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer :: rg
        integer :: i,j

        rg = Tdomain%rank

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
                        write (91,*) Tdomain%Ssource(i)%dir
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



    subroutine finalize_mesh_connectivity(Tdomain)
        use sdomain
        use semdatafiles
        use mpi
        implicit none

        type(domain), intent(inout) :: Tdomain
        integer :: i, j, k, n, nf, nnf, mat, ne, nv, nne, nnv
        integer :: n_aus


        ! faces and edges => which element?
        !do i = 0,Tdomain%n_face-1
        !    do j = 0, Tdomain%n_elem-1
        !        do k = 0,5
        !            if(Tdomain%specel(j)%Near_Faces(k) == i)then
        !                Tdomain%sFace(i)%Which_Elem = j
        !            endif
        !        enddo
        !    enddo
        !enddo
        do j = 0, Tdomain%n_elem-1
            do k = 0,5
                i = Tdomain%specel(j)%Near_Faces(k)
                Tdomain%sFace(i)%Which_Elem = j
            enddo
        enddo

        !do i = 0,Tdomain%n_edge-1
        !    allocate(Tdomain%sEdge(i)%Which_Elem(0:11))
        !    allocate(Tdomain%sEdge(i)%Which_EdgeinElem(0:11))
        !    Tdomain%sEdge(i)%Which_Elem = -1
        !    Tdomain%sEdge(i)%Which_EdgeinElem = -1
        !    icount = 0
        !    do j = 0, Tdomain%n_elem - 1
        !        do k=0,11
        !            if ( Tdomain%specel(j)%Near_Edges(k) == i )  then
        !                Tdomain%sEdge(i)%Which_Elem(icount) = j
        !                Tdomain%sEdge(i)%Which_EdgeinElem(icount) = k
        !                icount = icount+1
        !            endif
        !        enddo
        !    enddo
        !enddo

        ! material => time steps ; solid/liquid attributions
        do n = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            do nf = 0,5
                nnf = Tdomain%specel(n)%Near_Faces(nf)
                Tdomain%sFace(nnf)%mat_index = mat
                Tdomain%sFace(nnf)%solid = Tdomain%specel(n)%solid
                Tdomain%sFace(nnf)%fluid_dirich = .false.
              ! if(Tdomain%specel(n)%fluid_dirich .and. (nf == 5))then
              !     Tdomain%sFace(nnf)%fluid_dirich = .true.
              !  end if
            end do
            do ne = 0,11
                nne = Tdomain%specel(n)%Near_edges(ne)
                Tdomain%sEdge(nne)%mat_index = mat
                Tdomain%sEdge(nne)%solid = Tdomain%specel(n)%solid
                Tdomain%sEdge(nne)%fluid_dirich = .false.
               ! if(Tdomain%specel(n)%fluid_dirich .and.   &
               !    ((ne == 5).or.(ne == 8).or.(ne == 9).or.(ne == 11)))then
               !    Tdomain%sEdge(nne)%fluid_dirich = .true.
               ! end if
            end do
            do nv = 0,7
                nnv = Tdomain%specel(n)%Near_Vertices(nv)
                Tdomain%sVertex(nnv)%mat_index = mat
                Tdomain%sVertex(nnv)%solid = Tdomain%specel(n)%solid
                Tdomain%sVertex(nnv)%fluid_dirich = .false.
               ! if(Tdomain%specel(n)%fluid_dirich .and.   &
               !    ((nv == 4).or.(nv == 5).or.(nv == 6).or.(nv == 7)))then
               !   Tdomain%sVertex(nnv)%fluid_dirich = .true.
               ! end if
            end do

            !        !Faces
            !        do nf = 0, 5
            !            nnf = Tdomain%specel(n)%Near_Faces(nf)
            !            if(nnf == 2 .or. nnf == 3 .or. nnf == 5) then
            !                Tdomain%sFace(nnf)%mat_list(0) = mat
            !            else(nnf == 0 .or. nnf == 1 .or. nnf == 4) then
            !                Tdomain%sFace(nnf)%mat_list(1) = mat
            !            end if
            !        end do
            !        !Edges
            !        do ne = 0, 11
            !            nne = Tdomain%specel(n)%Near_edges(ne)
            !            if(nne == 0 .or. nne == 3 .or. nne == 6) then
            !                Tdomain%sEdge(nne)%mat_list(0) = mat
            !            else if (nne == 1 .or. nne == 2 .or. nne == 4) then
            !                Tdomain%sEdge(nne)%mat_list(1) = mat
            !            else if (nne == 7 .or. nne == 8 .or. nne == 9) then
            !                Tdomain%sEdge(nne)%mat_list(2) = mat
            !            else if (nne == 5 .or. nne == 10 .or. nne == 11) then
            !                Tdomain%sEdge(nne)%mat_list(3) = mat
            !            end if
            !        end do
            !        !Vertex
            !        do nv = 0, 7
            !            nnv = Tdomain%specel(n)%Near_Vertices(nv)
            !            Tdomain%sVertex(nnv)%mat_list(nnv) = mat
            !        end do

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

    subroutine read_material_file(Tdomain)
        use sdomain
        use semdatafiles
        use mpi
        implicit none

        type(domain), intent(inout) :: Tdomain
        character(Len=MAX_FILE_SIZE) :: fnamef
        integer :: i, j, n_aus, npml, mat, ne, nf, nRandom
        logical, dimension(:), allocatable :: L_Face, L_Edge
        real :: dtmin
        integer :: rg

        rg = Tdomain%rank
        npml = 0
        nRandom = 0
        !allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1)) !Changed to mesh3d.f90 line 337

        call semname_read_inputmesh_parametrage(Tdomain%material_file,fnamef)
        open (13, file=fnamef, status="old", form="formatted")

        read(13,*) n_aus

        if(n_aus /= Tdomain%n_mat) then
            write(*,*) trim(fnamef), n_aus, Tdomain%n_mat
            stop "Incompatibility between the mesh file and the material file for n_mat"
        endif

        if (Tdomain%aniso) then
            print *,"The code can't put anisotropy in a homogeneous media"
            stop
        endif

        allocate(Tdomain%not_PML_List(0:Tdomain%n_mat-1))
        Tdomain%any_PML      = .false.
        Tdomain%any_FPML     = .false.
        Tdomain%any_Random   = .false.
        Tdomain%not_PML_List = .true.

        do i = 0,Tdomain%n_mat-1
            read(13,*) Tdomain%sSubDomain(i)%material_type, &
                Tdomain%sSubDomain(i)%Pspeed,        &
                Tdomain%sSubDomain(i)%Sspeed,        &
                Tdomain%sSubDomain(i)%dDensity,      &
                Tdomain%sSubDomain(i)%NGLLx,         &
                Tdomain%sSubDomain(i)%NGLLy,         &
                Tdomain%sSubDomain(i)%NGLLz,         &
                Tdomain%sSubDomain(i)%Dt,            &
                Tdomain%sSubDomain(i)%Qpression,     &
                Tdomain%sSubDomain(i)%Qmu,           &
                Tdomain%sSubDomain(i)%nl_prop%LMC_prop%sigma_yld &
                Tdomain%sSubDomain(i)%nl_prop%LMC_prop%C_kin     &
                Tdomain%sSubDomain(i)%nl_prop%LMC_prop%kapa_kin  &
                Tdomain%sSubDomain(i)%nl_prop%LMC_prop%b_iso     &
                Tdomain%sSubDomain(i)%nl_prop%LMC_prop%Rinf_iso

            Tdomain%sSubdomain(i)%Filtering = .false.

            if(rg==0 .and. .false.) then
                write (*,*) 'Material :', i
                write (*,*) 'type     :', Tdomain%sSubDomain(i)%material_type
                write (*,*) 'Pspeed   :', Tdomain%sSubDomain(i)%Pspeed
                write (*,*) 'Sspeed   :', Tdomain%sSubDomain(i)%Sspeed
                write (*,*) 'Density  :', Tdomain%sSubDomain(i)%dDensity
                write (*,*) 'NGLL     :', Tdomain%sSubDomain(i)%NGLLx, Tdomain%sSubDomain(i)%NGLLy, Tdomain%sSubDomain(i)%NGLLz
                write (*,*) 'Dt       :', Tdomain%sSubDomain(i)%Dt
                write (*,*) 'Qp       :', Tdomain%sSubDomain(i)%Qpression
                write (*,*) 'Qmu      :', Tdomain%sSubDomain(i)%Qmu
                if (Tdomain%nl_flag) then
                    write (*,*) 'NL_LMC - sigma_yld:', Tdomain%sSubDomain(i)%nl_prop%LMC_prop%sigma_yld
                    write (*,*) 'NL_LMC - C_kin:',     Tdomain%sSubDomain(i)%nl_prop%LMC_prop%C_kin
                    write (*,*) 'NL_LMC - kapa_kin:',  Tdomain%sSubDomain(i)%nl_prop%LMC_prop%kapa_kin
                    write (*,*) 'NL_LMC - b_iso:',     Tdomain%sSubDomain(i)%nl_prop%LMC_prop%b_iso
                    write (*,*) 'NL_LMC - Rinf:',      Tdomain%sSubDomain(i)%nl_prop%LMC_prop%Rinf_iso
                end if
            endif

            Tdomain%sSubdomain(i)%assocMat = i

            call Lame_coefficients (Tdomain%sSubDomain(i))

            !        if(rg==0) &
            !            print*,' lame ',Tdomain%sSubDomain(i)%DMu,Tdomain%sSubDomain(i)%DLambda ,Tdomain%sSubDomain(i)%DKappa
            if (Tdomain%sSubDomain(i)%material_type == "P" .or. Tdomain%sSubDomain(i)%material_type == "L")  then
                npml = npml + 1
                Tdomain%not_PML_List(i) = .false.
            else
            endif

            if (Tdomain%sSubDomain(i)%material_type == "R") then
                nRandom = nRandom + 1
            end if
        enddo

        if(npml > 0) then
            if(rg==0) write (*,*) "!!WARNING change on 'material.input', put associated material after PML existing definition', "
            Tdomain%any_PML = .true.
            read(13,*); read(13,*)
            do i = 0,Tdomain%n_mat-1
                if(.not. Tdomain%not_PML_List(i)) then
                    read(13,*) Tdomain%sSubdomain(i)%Filtering, &
                        Tdomain%sSubdomain(i)%npow,      &
                        Tdomain%sSubdomain(i)%Apow,      &
                        Tdomain%sSubdomain(i)%Px,        &
                        Tdomain%sSubdomain(i)%Left,      &
                        Tdomain%sSubdomain(i)%Py,        &
                        Tdomain%sSubdomain(i)%Forward,   &
                        Tdomain%sSubdomain(i)%Pz,        &
                        Tdomain%sSubdomain(i)%Down,      &
                        Tdomain%sSubdomain(i)%freq,      &
                        Tdomain%sSubdomain(i)%assocMat
                    if(Tdomain%sSubdomain(i)%Filtering) Tdomain%any_FPML = .true.
                endif
            enddo
        endif

        if(nRandom > 0) then
            !Building element list in each subdomain
      !      do i=0,Tdomain%n_mat-1
      !          allocate (Tdomain%sSubdomain(i)%elemList(0:Tdomain%sSubdomain(i)%nElem-1))
      !          Tdomain%sSubdomain(i)%elemList(:) = -1 !-1 to detect errors
      !          Tdomain%sSubdomain(i)%nElem       = 0 !Using to count the elements in the next loop
      !      enddo
      !      do i=0,Tdomain%n_elem-1
      !          mat = Tdomain%specel(i)%mat_index
      !          Tdomain%sSubdomain(mat)%elemList(Tdomain%sSubdomain(mat)%nElem) = i
      !          Tdomain%sSubdomain(mat)%nElem = Tdomain%sSubdomain(mat)%nElem + 1
      !      enddo

      !      !Defining existing subdomain list in each domain
      !      allocate (Tdomain%subD_exist(0:Tdomain%n_mat-1))
      !      allocate (Tdomain%subDComm(0:Tdomain%n_mat - 1))
      !      Tdomain%subD_exist(:) = .true.
      !      do mat=0,Tdomain%n_mat-1
      !          if(Tdomain%sSubdomain(mat)%nElem == 0) Tdomain%subD_exist(mat) = .false.
      !      enddo

            Tdomain%any_Random = .true.
            read(13,*); read(13,*)
            do i = 0,Tdomain%n_mat-1
                !write(*,*) "Reading Random Material (", i, ")"
                if(Tdomain%sSubdomain(i)%material_type == "R") then
                    allocate(Tdomain%sSubdomain(i)%corrL(0:2))
                    allocate(Tdomain%sSubdomain(i)%varProp(0:2))
                    allocate(Tdomain%sSubdomain(i)%margiFirst(0:2))
                    !write(*,*) "Allocation finished"
                    !corrMod, corrLX, corrLY, corrLZ, margiFirstDens, varDens, margiFirstLambda, varLambda, margiFirstMu, varMu
                    read(13,*) Tdomain%sSubdomain(i)%corrMod,       &
                        Tdomain%sSubdomain(i)%corrL(0),      &
                        Tdomain%sSubdomain(i)%corrL(1),      &
                        Tdomain%sSubdomain(i)%corrL(2),      &
                        Tdomain%sSubdomain(i)%margiFirst(0), &
                        Tdomain%sSubdomain(i)%varProp(0),    &
                        Tdomain%sSubdomain(i)%margiFirst(1), &
                        Tdomain%sSubdomain(i)%varProp(1),    &
                        Tdomain%sSubdomain(i)%margiFirst(2), &
                        Tdomain%sSubdomain(i)%varProp(2),    &
                        Tdomain%sSubdomain(i)%seedStart
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
        real :: ndir

        Tdomain%n_source = config%nsources
        allocate (Tdomain%Ssource(0:Tdomain%n_source-1))

        call c_f_pointer(config%source, src)
        nsrc = 0
        do while(associated(src))
            Tdomain%Ssource(nsrc)%Xsource = src%coords(1)
            Tdomain%Ssource(nsrc)%Ysource = src%coords(2)
            Tdomain%Ssource(nsrc)%Zsource = src%coords(3)
            Tdomain%Ssource(nsrc)%i_type_source = src%type
            Tdomain%Ssource(nsrc)%amplitude_factor = src%amplitude
            if (src%func .eq. 5) then
                Tdomain%Ssource(nsrc)%time_file = fromcstr(src%time_file)
            end if
            ! Comportement temporel
            Tdomain%Ssource(nsrc)%i_time_function = src%func
            Tdomain%Ssource(nsrc)%cutoff_freq = src%freq ! func=2,4
            Tdomain%Ssource(nsrc)%tau_b = src%tau ! func=1,2,3,4,5
            Tdomain%Ssource(nsrc)%fh = src%band  ! func=3
            Tdomain%Ssource(nsrc)%gamma = src%gamma ! func=4
            Tdomain%Ssource(nsrc)%ts = src%ts   ! func=4
            Tdomain%Ssource(nsrc)%Q = src%Q
            Tdomain%Ssource(nsrc)%Y = src%Y
            Tdomain%Ssource(nsrc)%X = src%X
            Tdomain%Ssource(nsrc)%L = src%L
            Tdomain%Ssource(nsrc)%v = src%v
            Tdomain%Ssource(nsrc)%a = src%a
            Tdomain%Ssource(nsrc)%d = src%d

            ! Comportement Spacial
            ! i_type_source==1
            ndir = sqrt(src%dir(1)**2 + src%dir(2)**2 + src%dir(3)**2)
            if (ndir>0) Tdomain%Ssource(nsrc)%dir(0:2) = src%dir(1:3)/ndir

            ! i_type_source==2
            Tdomain%Ssource(nsrc)%Moment(0,0) = src%moments(1)
            Tdomain%Ssource(nsrc)%Moment(1,1) = src%moments(2)
            Tdomain%Ssource(nsrc)%Moment(2,2) = src%moments(3)
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


    function is_in_box(pos, box)
        real, dimension(3), intent(in) :: pos
        real, dimension(6), intent(in) :: box
        logical :: is_in_box
        !
        integer :: i

        is_in_box = .true.
        do i=1,3
            if (pos(i)<box(i)) is_in_box = .false.
            if (pos(i)>box(3+i)) is_in_box = .false.
        end do
    end function is_in_box
    !>
    ! Selectionne les elements pour les inclure ou non dans les snapshots
    !<
    subroutine select_output_elements(Tdomain, config)
        use sdomain
        implicit none
        type(domain), intent(inout)  :: Tdomain
        type(sem_config), intent(in) :: config

        type(sem_snapshot_cond), pointer :: selection
        integer :: n, i, ipoint
        real, dimension(3) :: pos
        logical :: sel

        do n = 0, Tdomain%n_elem-1
            call c_f_pointer(config%snapshot_selection, selection)
            do while(associated(selection))
                if (selection%include==1) then
                    sel = .true.
                else
                    sel = .false.
                end if
                pos = 0.
                do i=0,7
                    ipoint = Tdomain%specel(n)%Control_Nodes(i)
                    pos = pos + Tdomain%Coord_Nodes(:,ipoint)
                end do
                pos = pos/8
                select case (selection%type)
                case (1)
                    ! All
                    Tdomain%specel(n)%output = sel
                case (2)
                    ! Material
                    if (Tdomain%specel(n)%mat_index == selection%material) Tdomain%specel(n)%output = sel
                case (3)
                    ! Box
                    if (is_in_box(pos, selection%box)) Tdomain%specel(n)%output = sel
                end select

                call c_f_pointer(selection%next, selection)
            end do
        end do

    end subroutine select_output_elements

    subroutine read_input (Tdomain, code)
        use sdomain
        use semdatafiles
        use mpi
        use constants
        use mcapteur
        implicit none

        type(domain), intent(inout)  :: Tdomain
        integer, intent(out)         :: code
        character(Len=MAX_FILE_SIZE) :: fnamef
        logical                      :: logic_scheme
        integer                      :: imat
        integer                      :: rg

        rg = Tdomain%rank

        call semname_file_input_spec(fnamef)

        call read_sem_config(Tdomain%config, trim(fnamef)//C_NULL_CHAR, code)

        if (code/=1) then
            stop 1
        endif
        if (rg==0) call dump_config(Tdomain%config) !Print of configuration on the screen

        ! On copie les parametres renvoyes dans la structure C
        Tdomain%Title_simulation          = fromcstr(Tdomain%config%run_name)
        Tdomain%TimeD%acceleration_scheme = Tdomain%config%accel_scheme .ne. 0
        Tdomain%TimeD%velocity_scheme     = Tdomain%config%veloc_scheme .ne. 0
        Tdomain%TimeD%duration            = Tdomain%config%sim_time
        Tdomain%TimeD%alpha               = Tdomain%config%alpha
        Tdomain%TimeD%beta                = Tdomain%config%beta
        Tdomain%TimeD%gamma               = Tdomain%config%gamma
        Tdomain%nl_flag                   = Tdomain%config%nl_flag
        if (rg==0) then
            if (Tdomain%TimeD%alpha /= 0.5 .or. Tdomain%TimeD%beta /= 0.5 .or. Tdomain%TimeD%gamma /= 1.) then
                write(*,*) "***WARNING*** : Les parametres alpha,beta,gamma sont ignores dans cette version"
                write(*,*) "***WARNING*** : on prend: alpha=0.5, beta=0.5, gamma=1"
            endif
        end if
        Tdomain%TimeD%alpha = 0.5
        Tdomain%TimeD%beta = 0.5
        Tdomain%TimeD%gamma = 1.
        ! OUTPUT FIELDS
        Tdomain%out_variables(0:8) = Tdomain%config%out_variables
        Tdomain%nReqOut = sum(Tdomain%out_variables(0:3)) + 3*sum(Tdomain%out_variables(4:6)) + 6*sum(Tdomain%out_variables(7:8))

        Tdomain%TimeD%courant             = Tdomain%config%courant
        Tdomain%mesh_file                 = fromcstr(Tdomain%config%mesh_file)
        call semname_read_input_meshfile(rg,Tdomain%mesh_file,fnamef) !indicates the path to the mesh file for this proc"
        Tdomain%mesh_file             = fnamef
        Tdomain%aniso                 = Tdomain%config%anisotropy .ne. 0
        Tdomain%material_file         = fromcstr(Tdomain%config%mat_file)
        Tdomain%logicD%save_trace     = Tdomain%config%save_traces .ne. 0
        Tdomain%logicD%save_snapshots = Tdomain%config%save_snap .ne. 0
        Tdomain%ngroup                = Tdomain%config%n_group_outputs
        Tdomain%logicD%save_restart   = Tdomain%config%prorep_iter .ne. 0
        Tdomain%logicD%MPML           = .false.
        Tdomain%MPML_coeff            = Tdomain%config%mpml
        if (Tdomain%config%mpml/=0) then
            Tdomain%logicD%MPML = .true.
        end if
        Tdomain%logicD%grad_bassin = .false.
        Tdomain%logicD%run_restart = Tdomain%config%prorep .ne. 0
        Tdomain%TimeD%iter_reprise = Tdomain%config%prorep_restart_iter
        Tdomain%TimeD%ncheck       = Tdomain%config%prorep_iter ! frequence de sauvegarde

        Tdomain%TimeD%ntrace        = Tdomain%config%traces_interval ! XXX
        Tdomain%traces_format       = Tdomain%config%traces_format
        Tdomain%TimeD%time_snapshots = Tdomain%config%snap_interval
        logic_scheme                 = Tdomain%TimeD%acceleration_scheme .neqv. Tdomain%TimeD%velocity_scheme
        if(.not. logic_scheme) then
            stop "Both acceleration and velocity schemes: no compatibility, chose only one."
        end if

        ! Amortissement
        Tdomain%n_sls     = Tdomain%config%nsolids
        Tdomain%T1_att    = Tdomain%config%atn_band(1)
        Tdomain%T2_att    = Tdomain%config%atn_band(2)
        Tdomain%T0_modele = Tdomain%config%atn_period
        if (rg==0) then
            write(*,*) "Attenuation SLS =", Tdomain%n_sls
            write(*,*) "         period =", Tdomain%T0_modele
            write(*,*) "           band =", Tdomain%T1_att, Tdomain%T2_att
        end if


        ! Neumann boundary conditions? If yes: geometrical properties read in the mesh files.
        Tdomain%logicD%Neumann = Tdomain%config%neu_present /= 0
        Tdomain%Neumann%Neu_Param%what_bc = 'S'
        Tdomain%Neumann%Neu_Param%mat_index = Tdomain%config%neu_mat
        if (Tdomain%config%neu_type==1) then
            Tdomain%Neumann%Neu_Param%wtype = 'P'
        else
            Tdomain%Neumann%Neu_Param%wtype = 'S'
        end if
        Tdomain%Neumann%Neu_Param%lx = Tdomain%config%neu_L(1)
        Tdomain%Neumann%Neu_Param%ly = Tdomain%config%neu_L(2)
        Tdomain%Neumann%Neu_Param%lz = Tdomain%config%neu_L(3)
        Tdomain%Neumann%Neu_Param%xs = Tdomain%config%neu_C(1)
        Tdomain%Neumann%Neu_Param%ys = Tdomain%config%neu_C(2)
        Tdomain%Neumann%Neu_Param%zs = Tdomain%config%neu_C(3)
        Tdomain%Neumann%Neu_Param%f0 = Tdomain%config%neu_f0

        ! Create sources from C structures
        call create_sem_sources(Tdomain, Tdomain%config)

        !- Parametrage super object desactive
        Tdomain%logicD%super_object_local_present = .false.

        !---   Reading mesh file
        call read_mesh_file_h5(Tdomain)

        !---   Properties of materials.
        call read_material_file(Tdomain)

        ! Material Earthchunk

        Tdomain%earthchunk_isInit=0
        if( Tdomain%config%material_type == MATERIAL_EARTHCHUNK) then
            Tdomain%earthchunk_isInit=1

        endif


        if( Tdomain%config%material_present == 1) then

            select case (Tdomain%config%material_type)

            case (MATERIAL_PREM)
                Tdomain%aniso=.true.
            case (MATERIAL_EARTHCHUNK)
                Tdomain%earthchunk_isInit=1
                Tdomain%aniso=.true.
                Tdomain%earthchunk_file = fromcstr(Tdomain%config%model_file)
                Tdomain%earthchunk_delta_lon = Tdomain%config%delta_lon
                Tdomain%earthchunk_delta_lat = Tdomain%config%delta_lat

            end select

            do imat=0,Tdomain%n_mat-1
                Tdomain%sSubDomain(imat)%material_definition = Tdomain%config%material_type
            enddo
        else
            do imat=0,Tdomain%n_mat-1
                Tdomain%sSubDomain(imat)%material_definition = MATERIAL_CONSTANT
            enddo
        endif

        !write(*,*) rg, "Reading materials done"

        call finalize_mesh_connectivity(Tdomain)
        call select_output_elements(Tdomain, Tdomain%config)

    end subroutine read_input

end module semconfig

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
