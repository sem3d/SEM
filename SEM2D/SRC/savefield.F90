!>
!!\file savefield.F90
!!\brief Contient la subroutine savefield.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Sauve les champs : Field.
!!
!! \param type (domain), intent (INOUT) TDomain
!! \param integer, intent (IN) it
!! \param integer, intent (IN) i_snap
!! \param logical sortie_capteur_vitesse
!<

subroutine savefield (Tdomain,it,sortie_capteur,i_snap,nom_grandeur,nom_dir_sorties)
    
    use sdomain
    use semdatafiles

    implicit none

    type (domain), intent (INOUT):: TDomain
    integer, intent (IN) :: it,i_snap
    logical :: sortie_capteur
    character(len=10), intent(INOUT) :: nom_grandeur

    ! local variables
    integer :: ngll, ngllx,ngllz,ipoint, i,j,n_face,n, n_vertex
    real, dimension (0:1,0:Tdomain%n_glob_points-1) :: Field
    character (len=MAX_FILE_SIZE) :: fnamef
    character(len=24), intent(IN) :: nom_dir_sorties
    real, dimension(:,:,:,:), pointer :: buf_specel
    real, dimension (:,:,:), pointer :: buf_face
    real, dimension (:,:), pointer :: buf_vertex
    integer ngllmax

    ngllmax = 1
    do n = 0,Tdomain%n_elem-1
        ngllmax = max(ngllmax,Tdomain%specel(n)%ngllx)
        ngllmax = max(ngllmax,Tdomain%specel(n)%ngllz)
    enddo

    allocate (buf_specel(0:Tdomain%n_elem-1,0:ngllmax-1,0:ngllmax-1,0:1))
    allocate (buf_face(0:Tdomain%n_face-1,0:ngllmax-1,0:1))
    allocate (buf_vertex(0:Tdomain%n_vertex-1,0:1))

    if(trim(adjustl(nom_grandeur))=="vel") then
        do n = 0,Tdomain%n_elem-1
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            buf_specel(n,1:ngllx-2,1:ngllz-2,:) = Tdomain%specel(n)%Veloc(1:ngllx-2,1:ngllz-2,:)
        enddo

        do n = 0,Tdomain%n_face-1
            ngll = Tdomain%sFace(n)%ngll
            buf_face(n,1:ngll-2,:) = Tdomain%sFace(n)%Veloc(1:ngll-2,:)
        enddo

        do n = 0,Tdomain%n_vertex-1
            buf_vertex(n,:) = Tdomain%svertex(n)%Veloc(:)
        enddo
    else if (trim(adjustl(nom_grandeur))=="depla") then
        do n = 0,Tdomain%n_elem-1
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            if (Tdomain%specel(n)%PML ) then
                buf_specel(n,1:ngllx-2,1:ngllz-2,:) = 0.
            else
                buf_specel(n,1:ngllx-2,1:ngllz-2,:) = Tdomain%specel(n)%Displ(1:ngllx-2,1:ngllz-2,:)
            endif
        enddo

        do n = 0,Tdomain%n_face-1
            ngll = Tdomain%sFace(n)%ngll
            if (Tdomain%sFace(n)%PML ) then
                buf_face(n,1:ngll-2,:) = 0.
            else
                buf_face(n,1:ngll-2,:) = Tdomain%sFace(n)%Displ(1:ngll-2,:)
            endif
        enddo

        do n = 0,Tdomain%n_vertex-1
            if (Tdomain%sVertex(n)%PML) then
                buf_vertex(n,:) = 0.
            else
                buf_vertex(n,:) = Tdomain%svertex(n)%Displ(:)
            endif
        enddo
    else
        write(6,*) 'Erreur variable non reconnue',nom_grandeur
    endif


    Field = -2e6
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx;  ngllz = Tdomain%specel(n)%ngllz

        do j = 1,ngllz-2
            do i = 1,ngllx-2
                ipoint = Tdomain%specel(n)%Iglobnum(i,j)
                Field (0,ipoint) = buf_specel(n,i,j,0)
                Field (1,ipoint) = buf_specel(n,i,j,1)

            enddo
        enddo

        ! Face 0
        do i = 1, ngllx-2
            n_face = Tdomain%specel(n)%Near_face(0)
            if (Tdomain%sface(n_face)%Near_Element(0) == n) then
                ipoint = Tdomain%specel(n)%Iglobnum(i,0)
                Field (0,ipoint) = buf_face(n_face,i,0)
                Field (1,ipoint) = buf_face(n_face,i,1)

            endif
            n_face = Tdomain%specel(n)%Near_face(2)
            if (Tdomain%sface(n_face)%Near_Element(0) == n) then
                ipoint = Tdomain%specel(n)%Iglobnum(i,ngllz-1)
                Field (0,ipoint) = buf_face(n_face,i,0)
                Field (1,ipoint) = buf_face(n_face,i,1)

            endif
        enddo
        do j = 1, ngllz-2
            n_face = Tdomain%specel(n)%Near_face(1)
            if (Tdomain%sface(n_face)%Near_Element(0) == n) then
                ipoint = Tdomain%specel(n)%Iglobnum(ngllx-1,j)
                Field (0,ipoint) = buf_face(n_face,j,0)
                Field (1,ipoint) = buf_face(n_face,j,1)

            endif
            n_face = Tdomain%specel(n)%Near_face(3)
            if (Tdomain%sface(n_face)%Near_Element(0) == n) then
                ipoint = Tdomain%specel(n)%Iglobnum(0,j)
                Field (0,ipoint) = buf_face(n_face,j,0)
                Field (1,ipoint) = buf_face(n_face,j,1)

            endif
        enddo

        n_vertex = Tdomain%specel(n)%Near_Vertex(0)
        ipoint = Tdomain%specel(n)%Iglobnum(0,0)

        Field (0,ipoint) = buf_vertex(n_vertex,0)
        Field (1,ipoint) = buf_vertex(n_vertex,1)

        n_vertex = Tdomain%specel(n)%Near_Vertex(1)
        ipoint = Tdomain%specel(n)%Iglobnum(ngllx-1,0)
        Field (0,ipoint) = buf_vertex(n_vertex,0)
        Field (1,ipoint) = buf_vertex(n_vertex,1)

        n_vertex = Tdomain%specel(n)%Near_Vertex(2)
        ipoint = Tdomain%specel(n)%Iglobnum(ngllx-1,ngllz-1)
        Field (0,ipoint) = buf_vertex(n_vertex,0)
        Field (1,ipoint) = buf_vertex(n_vertex,1)

        n_vertex = Tdomain%specel(n)%Near_Vertex(3)
        ipoint = Tdomain%specel(n)%Iglobnum(0,ngllz-1)
        Field (0,ipoint) = buf_vertex(n_vertex,0)
        Field (1,ipoint) = buf_vertex(n_vertex,1)
    enddo

    do n = 0, Tdomain%n_glob_points-1
        if (abs(Field(0,n)) < 1e-40) Field(0,n) = 0.
        if (abs(Field(1,n)) < 1e-40) Field(1,n) = 0.
    enddo

    !    write(*,*) 'Ecriture des fichiers de sortie ok'

    deallocate(buf_specel)
    deallocate(buf_face)
    deallocate(buf_vertex)

    if (sortie_capteur) then
        if(trim(adjustl(nom_grandeur))=="vel") then
            do ipoint=0,Tdomain%n_glob_points-1
                Tdomain%GrandeurVitesse(0,ipoint) = Field (0,ipoint)
                Tdomain%GrandeurVitesse(1,ipoint) = Field (1,ipoint)
            enddo
        else if(trim(adjustl(nom_grandeur))=="depla") then
            do ipoint=0,Tdomain%n_glob_points-1
                Tdomain%GrandeurDepla(0,ipoint) = Field (0,ipoint)
                Tdomain%GrandeurDepla(1,ipoint) = Field (1,ipoint)
            enddo
        endif
    endif


    if (Tdomain%logicD%save_snapshots.and.i_snap == 0) then

#ifdef MKA3D
        !! les fichiers de sortie de deplacements pour Sem sont appeles displ
        !! alors que les fichiers de capteurs s'appellent *depla*
        if(trim(adjustl(nom_grandeur))=="depla") then
            nom_grandeur="displ"
        endif
        call semname_savefield_dat(nom_dir_sorties,nom_grandeur,Tdomain%mpi_var%my_rank,it,fnamef)
        open (61,file=fnamef,status="unknown",form="formatted")
#else
        call semname_savefield_fieldt(Tdomain%mpi_var%my_rank,it,fnamef)
        open (61,file=fnamef,status="unknown",form="formatted")
#endif


        do ipoint=0,Tdomain%n_glob_points-1

#ifdef MKA3D
            write (61,"(I8,3G17.8)") ipoint, Field (0,ipoint),Field (1,ipoint),0.
#else
            if (.not. Field(0,ipoint) < -1e6 ) &
                write (61,"(3G17.8)") Tdomain%GlobCoord (0,ipoint), Tdomain%GlobCoord (1,ipoint), Field (0,ipoint)**2+Field (1,ipoint)**2
#endif
        enddo
        call flush(61)
        close(61)

    endif


    return
end subroutine savefield
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
