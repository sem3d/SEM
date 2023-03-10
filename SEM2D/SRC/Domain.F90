!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Domain.F90
!!\brief Contient le definition du type domain
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module sdomain

    ! Modified by Gaetano 12/10/05

    use selement
    use sfaces
    use svertices
    use smortars
    use ssources
    use stimeparam
    use logical_input
    use ssubdomains
    use sreceivers
    use sfaults
    use mpi_list
    use communication_object
    use semdatafiles
    use constants

    type :: domain
       ! Communicateur pour tous les processus SEM
       integer :: communicateur
       ! Communicateurs utilises uniquement pour le couplage
       integer :: communicateur_global
       integer :: master_superviseur
       !
       integer :: n_elem, n_face, n_vertex, n_source, n_glob_nodes, n_line ,n_receivers, n_mortar
       integer :: n_nodes, n_mat,n_glob_points, n_super_object, n_fault, n_communications
       integer :: type_timeInteg, type_elem, type_flux, type_bc, pml_type, capt_loc_type, Implicitness

       integer, dimension (:), pointer :: Line_index, Communication_list
       integer :: n_quad ! Total number of quad elements to output (including subelements)
       real(fpp), dimension (:,:), pointer :: Coord_nodes, GlobCoord
       real(fpp), dimension(:,:,:), pointer :: Store_Trace


       real(fpp), dimension(:,:), pointer :: GrandeurVitesse     ! tableau conservatoire des grandeurs
       real(fpp), dimension(:,:), pointer :: GrandeurDepla     ! tableau conservatoire des grandeurs
       real(fpp), dimension(:), pointer :: GrandeurDeformation ! pour tous les points de gauss a une iteration donnee (utilise pour les sorties capteurs)
       ! Additions S. Terrana Janv 2014
       integer :: nCapt
       integer, dimension (:,:), pointer :: elems_capteurs, faces_capteurs
       integer, dimension (:),   pointer :: type_capteurs


       character (len=MAX_FILE_SIZE) :: Title_simulation, mesh_file, station_file, material_file
       character (len=1), dimension (:), pointer ::  Name_Line, Super_object_type
       character (len=MAX_FILE_SIZE), dimension (:), pointer :: Super_object_file

       logical :: any_PML,bMailUnv,bCapteur,openfilescapt

       type (time) :: TimeD
       type (logical_array) :: logicD
       type (mpi_objects) :: Mpi_var
       type(element), dimension(:), pointer :: specel
       type(face), dimension (:), pointer :: sface
       type (vertex), dimension (:), pointer :: svertex
       type (mortar), dimension (:), pointer :: smortar
       type (source), dimension (:), pointer :: sSource
       type (subdomain), dimension (:), pointer :: sSubDomain
       type (receiver), dimension (:), pointer :: sReceiver
       type (fault), dimension (:), pointer :: sFault
       type (communicating_wall), dimension (:), pointer :: sWall



    end type domain

contains

    subroutine dist_max_elem(Tdomain)

        type (Domain), intent (INOUT) :: Tdomain
        integer ipoint, jpoint
        integer n
        real(fpp) :: coor_i(0:1), coor_j(0:1)
        real(fpp) :: dist_max


        do n = 0,Tdomain%n_elem-1
            dist_max = 0.

            do i=0,Tdomain%n_nodes-1
                ipoint = Tdomain%specel(n)%Control_Nodes(i)
                coor_i = Tdomain%Coord_nodes(0:1,ipoint)
                do j=i+1,Tdomain%n_nodes-1
                    jpoint = Tdomain%specel(n)%Control_Nodes(j)
                    coor_j = Tdomain%Coord_nodes(0:1,jpoint)
                    dist_max = max(dist_max, sqrt((coor_i(0)-coor_j(0))**2+(coor_i(1)-coor_j(1))**2))
                enddo
            enddo
            Tdomain%specel(n)%dist_max = dist_max
        enddo

    end subroutine dist_max_elem


subroutine read_material_file(Tdomain)
    use semdatafiles
    use mpi
    implicit none

    type(domain), intent(inout) :: Tdomain
    character(Len=MAX_FILE_SIZE) :: fnamef
    !
    real(fpp) :: dtmin
    integer   :: i, j, mat, mat2, npml, nmortar, ngll1, ngll2
    integer   :: w_face, w_face2, n_aus, k_aus
    real(fpp) :: Qp, Qs

    ! Read material properties
    npml = 0
    allocate(Tdomain%sSubdomain(0:Tdomain%n_mat-1))

    call semname_read_inputmesh_parametrage(Tdomain%material_file,fnamef)
    open (13,file=fnamef, status="old", form="formatted")

    read (13,*) n_aus  ! Number of material in file
    if (n_aus<Tdomain%n_mat) then
        write(*,*) "Material file doesn't contain enough definitions"
        stop 1
    end if
    do i = 0, Tdomain%n_mat-1
        read (13,*) Tdomain%sSubDomain(i)%material_type, Tdomain%sSubDomain(i)%Pspeed, &
            Tdomain%sSubDomain(i)%Sspeed, Tdomain%sSubDomain(i)%dDensity, &
            Tdomain%sSubDomain(i)%NGLLx, n_aus, Tdomain%sSubDomain(i)%NGLLz, Tdomain%sSubDomain(i)%Dt, &
            Qp, Qs
        Tdomain%sSubDomain(i)%n_loc_dim = 2
        if ( Tdomain%sSubDomain(i)%NGLLx == Tdomain%sSubDomain(i)%NGLLz)  Tdomain%sSubDomain(i)%n_loc_dim = 1
        if (Tdomain%sSubDomain(i)%material_type == "P" )  then
            npml = npml + 1
        endif
        ! Pour l'instant, on a un seul type de Flux et d'Elements pour TOUT le domaine
        Tdomain%sSubDomain(i)%type_DG   = Tdomain%type_Elem
        Tdomain%sSubDomain(i)%type_Flux = Tdomain%type_Flux
        ! POUR COUPLAGE CG-DG
        !Tdomain%sSubDomain(1)%type_DG   = 0 ! <---- A SUPPRIMER !!!!!!
        !Tdomain%sSubDomain(1)%type_Flux = 0 ! <---- A SUPPRIMER !!!!!!
        !Tdomain%sSubDomain(0)%type_DG   = 3 ! <---- A SUPPRIMER !!!!!!
        !Tdomain%sSubDomain(0)%type_Flux = 4 ! <---- A SUPPRIMER !!!!!!
    enddo
    Tdomain%any_PML  = .false.
    if (npml > 0 ) then
        Tdomain%any_PML = .true.
        read(13,*); read(13,*)
        do i = 0,Tdomain%n_mat-1
            if (Tdomain%sSubdomain(i)%material_type == "P" ) then
                read (13,*) Tdomain%sSubdomain(i)%Filtering,  Tdomain%sSubdomain(i)%npow, Tdomain%sSubdomain(i)%Apow, &
                    Tdomain%sSubdomain(i)%Px, Tdomain%sSubdomain(i)%Left, Tdomain%sSubdomain(i)%Pz,  &
                    Tdomain%sSubdomain(i)%Down, Tdomain%sSubdomain(i)%freq, Tdomain%sSubdomain(i)%k
                ! Warning : The variable "Filtering" is no longer used : the kind of PML|FPML|CPML|ADEPML is
                ! assigned directly in the file : input.spec :
                Tdomain%sSubdomain(i)%pml_type = Tdomain%pml_type
            endif
        enddo
    endif

    close (13)

    do i = 0,Tdomain%n_elem-1
        mat = Tdomain%specel(i)%mat_index
        Tdomain%specel(i)%ngllx = Tdomain%sSubDomain(mat)%NGLLx
        Tdomain%specel(i)%ngllz = Tdomain%sSubDomain(mat)%NGLLz
        Tdomain%specel(i)%type_DG = Tdomain%sSubDomain(mat)%type_DG
        ! Checks if the current element is on acoustic part of the domain
        if (Tdomain%sSubdomain(mat)%material_type .EQ. "F") then
           Tdomain%specel(i)%acoustic = .true.
        elseif ((Tdomain%sSubdomain(mat)%material_type .EQ. "P") .AND. &
                (Tdomain%sSubdomain(mat)%Sspeed .EQ. 0.)) then
            Tdomain%specel(i)%acoustic = .true.
        else
           Tdomain%specel(i)%acoustic = .false.
        endif
    enddo

    do i = 0, Tdomain%n_face-1
        n_aus = Tdomain%sFace(i)%Near_Element(0)
        w_face = Tdomain%sFace(i)%Which_face(0)
        mat = Tdomain%specel(n_aus)%mat_index
        if (w_face == 0 .or. w_face==2) then
            Tdomain%sFace(i)%ngll = Tdomain%sSubDomain(mat)%ngllx
        else
            Tdomain%sFace(i)%ngll = Tdomain%sSubDomain(mat)%ngllz
        endif
        Tdomain%sFace(i)%mat_index = mat
        n_aus = Tdomain%sFace(i)%Near_Vertex(0)
        Tdomain%sVertex(n_aus)%mat_index = mat
        n_aus = Tdomain%sFace(i)%Near_Vertex(1)
        Tdomain%sVertex(n_aus)%mat_index = mat
        ! Sebadditions for DG - Sept 2013
        Tdomain%sFace(i)%type_Flux = Tdomain%sSubDomain(mat)%type_Flux
        Tdomain%sFace(i)%type_DG   = Tdomain%sSubDomain(mat)%type_DG
        ! Check if the Face is an elastic-acoustic interface :
        n_aus = Tdomain%sFace(i)%Near_Element(0)
        k_aus = Tdomain%sFace(i)%Near_Element(1)
        ! Set the Acoustic Face's Flag to be the same as Near_Element(0) :
        Tdomain%sFace(i)%acoustic = Tdomain%specel(n_aus)%acoustic
        if (k_aus .NE.-1) then
           if (Tdomain%specel(n_aus)%acoustic .EQV. Tdomain%specel(k_aus)%acoustic) then
              Tdomain%sFace(i)%changing_media = .false.
           else
              Tdomain%sFace(i)%changing_media = .true.
              Tdomain%sFace(i)%acoustic = .false.
           endif
        else
           Tdomain%sFace(i)%changing_media = .false.
        endif
    enddo

    ! CREATION DES MORTARS
    nmortar = 0
    do i = 0, Tdomain%n_face-1
        Tdomain%sFace(i)%mortar = .false.
        n_aus = Tdomain%sFace(i)%Near_Element(0)
        k_aus = Tdomain%sFace(i)%Near_Element(1)
        mat = Tdomain%specel(n_aus)%mat_index
        if (k_aus .NE.-1) then
            mat2 = Tdomain%specel(k_aus)%mat_index
            if (mat .NE. mat2) then
                w_face = Tdomain%sFace(i)%Which_face(0)
                w_face2= Tdomain%sFace(i)%Which_face(1)
                if (w_face == 0 .or. w_face==2) then
                    ngll1 = Tdomain%sSubDomain(mat)%ngllx
                else
                    ngll1 = Tdomain%sSubDomain(mat)%ngllz
                endif
                if (w_face2 == 0 .or. w_face2==2) then
                    ngll2 = Tdomain%sSubDomain(mat2)%ngllx
                else
                    ngll2 = Tdomain%sSubDomain(mat2)%ngllz
                endif
                if (ngll1.NE.ngll2) then
                    Tdomain%sFace(i)%ngll = max(ngll1,ngll2)
                    nmortar = nmortar + 1
                    Tdomain%sFace(i)%mortar = .true.
                endif
            endif
        endif
    enddo

    ! Traitement des mortars :
    Tdomain%n_mortar = nmortar
    write(*,*) "Nombre de faces mortars : ", nmortar
    allocate(Tdomain%sMortar(0:nmortar-1))
    nmortar = 0
    do i = 0, Tdomain%n_face-1
        if (Tdomain%sFace(i)%mortar) then
            Tdomain%sMortar(nmortar)%Near_Face(0) = i
            Tdomain%sFace(i)%mortarId = nmortar
            n_aus = Tdomain%sFace(i)%Near_Element(0)
            k_aus = Tdomain%sFace(i)%Near_Element(1)
            w_face = Tdomain%sFace(i)%Which_face(0)
            w_face2= Tdomain%sFace(i)%Which_face(1)
            mat = Tdomain%specel(n_aus)%mat_index
            mat2 = Tdomain%specel(k_aus)%mat_index
            if (w_face == 0 .or. w_face==2) then
                ngll1 = Tdomain%sSubDomain(mat)%ngllx
            else
                ngll1 = Tdomain%sSubDomain(mat)%ngllz
            endif
            if (w_face2 == 0 .or. w_face2==2) then
                ngll2 = Tdomain%sSubDomain(mat2)%ngllx
            else
                ngll2 = Tdomain%sSubDomain(mat2)%ngllz
            endif
            Tdomain%sMortar(nmortar)%ngllmin = min(ngll1,ngll2)
            Tdomain%sMortar(nmortar)%ngllmax = max(ngll1,ngll2)
            if (ngll1>ngll2) then ! Near_Element(0) doit etre le "slave"
                Tdomain%sMortar(nmortar)%Near_Element(0) = k_aus
                Tdomain%sMortar(nmortar)%Near_Element(1) = n_aus
            else
                Tdomain%sMortar(nmortar)%Near_Element(0) = n_aus
                Tdomain%sMortar(nmortar)%Near_Element(1) = k_aus
            endif
            nmortar = nmortar + 1
        endif
    enddo

    ! Modifs pour pr??server DG_STRONG et DG_WEAK avec la nouvelle formulation HDG-acoustique
    if (Tdomain%type_elem==GALERKIN_DG_STRONG .OR. Tdomain%type_elem==GALERKIN_DG_WEAK) then
        do i = 0, Tdomain%n_face-1
            if (Tdomain%sFace(i)%changing_media) then
                n_aus = Tdomain%sFace(i)%Near_Element(0)
                if (Tdomain%specel(n_aus)%acoustic) then
                    Tdomain%sFace(i)%acoustic = .true.
                else
                    Tdomain%sFace(i)%acoustic = .false.
                endif
            endif
        enddo
        do i = 0,Tdomain%n_elem-1
            Tdomain%specel(i)%acoustic = .false.
        enddo
    endif

    do i = 0, Tdomain%n_vertex-1
        mat = Tdomain%sVertex(i)%mat_index
        Tdomain%sVertex(i)%Type_DG = Tdomain%sSubDomain(mat)%Type_DG
    enddo

    if (Tdomain%logicD%super_object_local_present) then
        do i = 0, Tdomain%n_fault-1
            do j = 0, Tdomain%sFault(i)%n_face-1
                n_aus = Tdomain%sFault(i)%fFace(j)%Face_UP
                Tdomain%sFault(i)%fFace(j)%ngll = Tdomain%sFace(n_aus)%ngll
                w_face = Tdomain%sFace(n_aus)%Which_face(0)
                if (w_face == 0 .or. w_face ==2) then
                    Tdomain%sFault(i)%fFace(j)%mat_dir = 1
                else
                    Tdomain%sFault(i)%fFace(j)%mat_dir = 0
                endif
                Tdomain%sFault(i)%fFace(j)%mat_index = Tdomain%sFace(n_aus)%mat_index
                n_aus = Tdomain%sFault(i)%fFace(j)%Face_to_Vertex(0)
                Tdomain%sFault(i)%fVertex(n_aus)%mat_index =  Tdomain%sFault(i)%fFace(j)%mat_index
                n_aus = Tdomain%sFault(i)%fFace(j)%Face_to_Vertex(1)
                Tdomain%sFault(i)%fVertex(n_aus)%mat_index =  Tdomain%sFault(i)%fFace(j)%mat_index
            enddo
            do j = 0, Tdomain%sFault(i)%n_vertex-1
                Tdomain%sFault(i)%fVertex(j)%Termination = .false.
                if (Tdomain%sFault(i)%fVertex(j)%Vertex_UP == -3) &
                    Tdomain%sFault(i)%fVertex(j)%Termination = .true.
            enddo
        enddo
    endif

    do i = 0, Tdomain%n_mat-1
        if (Tdomain%sSubdomain(i)%material_type == "P" &
            .and. Tdomain%sSubdomain(i)%Filtering &
            .and. (Tdomain%pml_type == 2) ) then
            Tdomain%sSubdomain(i)%freq = exp (-Tdomain%sSubdomain(i)%freq*Tdomain%sSubdomain(i)%dt/2)
        endif
    enddo

    dtmin =1e20
    do i = 0,Tdomain%n_mat-1
        if (Tdomain%sSubDomain(i)%Dt < dtmin ) dtmin = Tdomain%sSubDomain(i)%Dt
    enddo
    Tdomain%TimeD%dtmin = dtmin
    if (dtmin > 0) then
        Tdomain%TimeD%ntimeMax = int (Tdomain%TimeD%Duration/dtmin)
    else
        write (*,*) "Your dt min is zero : verify it"
        stop
    endif

    do i = 0, Tdomain%n_mat-1
        call Lame_coefficients (Tdomain%sSubDomain(i))
    enddo

end subroutine read_material_file


end module sdomain

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
