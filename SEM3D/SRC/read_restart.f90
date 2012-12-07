!>
!!\file read_restart.f90
!!\brief Contient la subroutine read_restart().
!!
!! Gère la reprise de Sem3d

!>
!! \brief Assure la reprise par la lecture des fichiers de protection.
!!
!! \param type (domain), intent (inout) Tdomain
!! \param integer,intent (inout) isort
!<

#ifdef MKA3D
subroutine init_restart(Tdomain, rg, isort, fnamef)
    use sdomain
    use semdatafiles
    implicit none
    type (domain), intent (IN)       :: Tdomain
    integer, intent (in)             :: isort
    integer, intent (IN)             :: rg
    character (len=MAX_FILE_SIZE), intent(out) :: fnamef
    character (len=MAX_FILE_SIZE)    :: fnamer, fnamec
    character (len=MAX_FILE_SIZE)    :: commande


    call semname_couplage_iter(Tdomain%TimeD%iter_reprise,rg,fnamef)
    call semname_couplage_iterr(Tdomain%TimeD%iter_reprise,fnamer)
    if (rg == 0) then
        ! copie du fichier temps.dat dans le rep de Resultat
        call semname_couplage_commandecpt(fnamer, fnamec)
        commande="cp "//trim(adjustl(fnamec)) !!modif 09/11
        call system(commande)

        ! copie du repertoire des sorties capteurs sem dans le rep de resultats
        call semname_couplage_commanderm(fnamec)
        commande="rm -Rf "//trim(adjustl(fnamec))
        call system(commande)
        call semname_couplage_commandecp(fnamer,fnamec)
        commande="cp -r "//trim(adjustl(fnamec))
        call system(commande)
    endif
end subroutine init_restart

#else
subroutine init_restart(Tdomain, rg, isort, fnamef)
    use sdomain
    use semdatafiles
    implicit none
    type (domain), intent (IN)       :: Tdomain
    integer, intent (in)             :: isort
    integer, intent (IN)             :: rg
    character (len=MAX_FILE_SIZE), intent(out) :: fnamef

    call semname_read_restart_save_checkpoint_rank(rg,fnamef)
end subroutine init_restart
#endif

subroutine clean_prot(Tdomain, rg)
    use sdomain
    use semdatafiles
    implicit none
    type (domain), intent (IN)       :: Tdomain
    integer, intent (IN)             :: rg
    character (len=MAX_FILE_SIZE)    :: commande
    character (len=6) :: sit

    if (rg == 0) then
        write(sit,'(I6)') Tdomain%TimeD%prot_m0
        commande='find ./ProRep/sem -name "Prot*" ! -name "Prot*'//trim(adjustl(sit))//'*" -exec rm -fr  {} \; '
        ! on supprime les fichiers et repertoire de protection autres que celui contenant prot_m0
        !write(6,*) 'commande',commande
        call system(commande)
    endif
end subroutine clean_prot

#ifdef USE_HDF5
subroutine read_Veloc(Tdomain, elem_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    real*8, allocatable, dimension(:) :: veloc, veloc1, veloc2, veloc3, displ
    integer :: idx1, idx2, idx3, ngllx, nglly, ngllz
    integer :: n, i, j, k


    call read_dset_1d_real(elem_id, "Veloc", veloc)
    call read_dset_1d_real(elem_id, "Veloc1", veloc1)
    call read_dset_1d_real(elem_id, "Veloc2", veloc2)
    call read_dset_1d_real(elem_id, "Veloc3", veloc3)
    call read_dset_1d_real(elem_id, "Displ", displ)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 1,ngllz-2
            do j = 1,nglly-2
                do i = 1,ngllx-2
                    Tdomain%specel(n)%Veloc(i,j,k,0) = veloc(idx1+0)
                    Tdomain%specel(n)%Veloc(i,j,k,1) = veloc(idx1+1)
                    Tdomain%specel(n)%Veloc(i,j,k,2) = veloc(idx1+2)
                    idx1 = idx1 + 3
                    if ( .not. Tdomain%specel(n)%PML ) then
                        Tdomain%specel(n)%Displ(i,j,k,0) = displ(idx2+0)
                        Tdomain%specel(n)%Displ(i,j,k,1) = displ(idx2+1)
                        Tdomain%specel(n)%Displ(i,j,k,2) = displ(idx2+2)
                        idx2 = idx2 + 3
                    else
                        Tdomain%specel(n)%Veloc1(i,j,k,0) = veloc1(idx3+0)
                        Tdomain%specel(n)%Veloc1(i,j,k,1) = veloc1(idx3+1)
                        Tdomain%specel(n)%Veloc1(i,j,k,2) = veloc1(idx3+2)
                        Tdomain%specel(n)%Veloc2(i,j,k,0) = veloc2(idx3+0)
                        Tdomain%specel(n)%Veloc2(i,j,k,1) = veloc2(idx3+1)
                        Tdomain%specel(n)%Veloc2(i,j,k,2) = veloc2(idx3+2)
                        Tdomain%specel(n)%Veloc3(i,j,k,0) = veloc3(idx3+0)
                        Tdomain%specel(n)%Veloc3(i,j,k,1) = veloc3(idx3+1)
                        Tdomain%specel(n)%Veloc3(i,j,k,2) = veloc3(idx3+2)
                        idx3 = idx3 + 3
                    end if
                end do
            end do
        end do
    end do
    deallocate(veloc)
    deallocate(veloc1)
    deallocate(veloc2)
    deallocate(veloc3)
    deallocate(displ)
end subroutine read_Veloc

subroutine read_EpsilonVol(Tdomain, elem_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    real*8, allocatable, dimension(:) :: epsilonvol, rvol, rxx, ryy, rxy, rxz, ryz
    integer :: idx1, idx2, idx3, ngllx, nglly, ngllz
    integer :: n, i, j, k, i_sls
    integer :: n_solid

    call read_dset_1d_real(elem_id, "EpsilonVol", epsilonVol)
    call read_dset_1d_real(elem_id, "Rvol", rvol)
    call read_dset_1d_real(elem_id, "R_xx", rxx)
    call read_dset_1d_real(elem_id, "R_yy", ryy)
    call read_dset_1d_real(elem_id, "R_xy", rxy)
    call read_dset_1d_real(elem_id, "R_xz", rxz)
    call read_dset_1d_real(elem_id, "R_yz", ryz)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    n_solid = Tdomain%n_sls
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if ( .not. Tdomain%specel(n)%PML ) then
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        if (n_solid>0) then
                            if (Tdomain%aniso) then
                            else
                                Tdomain%specel(n)%epsilonvol_(i,j,k) = epsilonVol(idx1)
                                idx1 = idx1 + 1
                                do i_sls = 0, n_solid-1
                                    Tdomain%specel(n)%R_vol_(i_sls,i,j,k) = rvol(idx2+i_sls)
                                end do
                                idx2 = idx2 + n_solid
                            end if
                            do i_sls = 0, n_solid-1
                                Tdomain%specel(n)%R_xx_(i_sls,i,j,k) = rxx(idx3+i_sls)
                                Tdomain%specel(n)%R_yy_(i_sls,i,j,k) = ryy(idx3+i_sls)
                                Tdomain%specel(n)%R_xy_(i_sls,i,j,k) = rxy(idx3+i_sls)
                                Tdomain%specel(n)%R_xz_(i_sls,i,j,k) = rxz(idx3+i_sls)
                                Tdomain%specel(n)%R_yz_(i_sls,i,j,k) = ryz(idx3+i_sls)
                            end do
                            idx3 = idx3 + n_solid
                        end if
                    end do
                end do
            end do
        end if
    end do
    deallocate(epsilonvol)
    deallocate(rvol)
    deallocate(rxx)
    deallocate(ryy)
    deallocate(rxy)
    deallocate(rxz)
    deallocate(ryz)
end subroutine read_EpsilonVol

subroutine read_EpsilonDev(Tdomain, elem_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    real*8, allocatable, dimension(:) :: eps_dev_xx, eps_dev_yy, eps_dev_xy, eps_dev_xz, eps_dev_yz
    integer :: idx, ngllx, nglly, ngllz
    integer :: n, i, j, k, i_sls
    integer :: n_solid

    call read_dset_1d_real(elem_id, "EpsilonDev_xx", eps_dev_xx)
    call read_dset_1d_real(elem_id, "EpsilonDev_yy", eps_dev_yy)
    call read_dset_1d_real(elem_id, "EpsilonDev_xy", eps_dev_xy)
    call read_dset_1d_real(elem_id, "EpsilonDev_xz", eps_dev_xz)
    call read_dset_1d_real(elem_id, "EpsilonDev_yz", eps_dev_yz)
    idx = 1
    n_solid = Tdomain%n_sls
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if ( .not. Tdomain%specel(n)%PML ) then
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        if (n_solid>0) then
                            Tdomain%specel(n)%epsilondev_xx_(i,j,k) = eps_dev_xx(idx)
                            Tdomain%specel(n)%epsilondev_yy_(i,j,k) = eps_dev_yy(idx)
                            Tdomain%specel(n)%epsilondev_xy_(i,j,k) = eps_dev_xy(idx)
                            Tdomain%specel(n)%epsilondev_xz_(i,j,k) = eps_dev_xz(idx)
                            Tdomain%specel(n)%epsilondev_yz_(i,j,k) = eps_dev_yz(idx)
                            idx = idx + 1
                        end if
                    end do
                end do
            end do
        end if
    end do
    deallocate(eps_dev_xx)
    deallocate(eps_dev_yy)
    deallocate(eps_dev_xy)
    deallocate(eps_dev_xz)
    deallocate(eps_dev_yz)
end subroutine read_EpsilonDev

subroutine read_Stress(Tdomain, elem_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    real*8, allocatable, dimension(:) :: stress
    integer :: idx, ngllx, nglly, ngllz
    integer :: n, i, j, k, i_sls
    integer :: n_solid

    call read_dset_1d_real(elem_id, "Stress", stress)
    idx = 1
    n_solid = Tdomain%n_sls
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if (Tdomain%specel(n)%PML ) then
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        Tdomain%specel(n)%Diagonal_Stress1(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%Diagonal_Stress1(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%Diagonal_Stress1(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%Diagonal_Stress2(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%Diagonal_Stress2(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%Diagonal_Stress2(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%Diagonal_Stress3(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%Diagonal_Stress3(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%Diagonal_Stress3(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%Residual_Stress1(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%Residual_Stress1(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%Residual_Stress1(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%Residual_Stress2(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%Residual_Stress2(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%Residual_Stress2(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                    end do
                end do
            end do
        end if
    end do
    deallocate(stress)
end subroutine read_Stress

subroutine read_Faces(Tdomain, face_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: face_id
    real*8, allocatable, dimension(:) :: veloc, veloc1, veloc2, veloc3, displ
    integer :: idx1, idx2, idx3, ngll1, ngll2
    integer :: n, i, j


    call read_dset_1d_real(face_id, "Veloc", veloc)
    call read_dset_1d_real(face_id, "Veloc1", veloc1)
    call read_dset_1d_real(face_id, "Veloc2", veloc2)
    call read_dset_1d_real(face_id, "Veloc3", veloc3)
    call read_dset_1d_real(face_id, "Displ", displ)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    do n = 0,Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1
        ngll2 = Tdomain%sFace(n)%ngll2

        do j = 1,ngll2-2
            do i = 1,ngll1-2
                Tdomain%sFace(n)%Veloc(i,j,0) = veloc(idx1+0)
                Tdomain%sFace(n)%Veloc(i,j,1) = veloc(idx1+1)
                Tdomain%sFace(n)%Veloc(i,j,2) = veloc(idx1+2)
                idx1 = idx1 + 3
                if ( .not. Tdomain%sFace(n)%PML ) then
                    Tdomain%sFace(n)%Displ(i,j,0) = displ(idx2+0)
                    Tdomain%sFace(n)%Displ(i,j,1) = displ(idx2+1)
                    Tdomain%sFace(n)%Displ(i,j,2) = displ(idx2+2)
                    idx2 = idx2 + 3
                else
                    Tdomain%sFace(n)%Veloc1(i,j,0) = veloc1(idx3+0)
                    Tdomain%sFace(n)%Veloc1(i,j,1) = veloc1(idx3+1)
                    Tdomain%sFace(n)%Veloc1(i,j,2) = veloc1(idx3+2)
                    Tdomain%sFace(n)%Veloc2(i,j,0) = veloc2(idx3+0)
                    Tdomain%sFace(n)%Veloc2(i,j,1) = veloc2(idx3+1)
                    Tdomain%sFace(n)%Veloc2(i,j,2) = veloc2(idx3+2)
                    Tdomain%sFace(n)%Veloc3(i,j,0) = veloc3(idx3+0)
                    Tdomain%sFace(n)%Veloc3(i,j,1) = veloc3(idx3+1)
                    Tdomain%sFace(n)%Veloc3(i,j,2) = veloc3(idx3+2)
                    idx3 = idx3 + 3
                end if
            end do
        end do
    end do
    deallocate(veloc)
    deallocate(veloc1)
    deallocate(veloc2)
    deallocate(veloc3)
    deallocate(displ)
end subroutine read_Faces

subroutine read_Edges(Tdomain, edge_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: edge_id
    real*8, allocatable, dimension(:) :: veloc, veloc1, veloc2, veloc3, displ
    integer :: idx1, idx2, idx3, ngll
    integer :: n, i


    call read_dset_1d_real(edge_id, "Veloc", veloc)
    call read_dset_1d_real(edge_id, "Veloc1", veloc1)
    call read_dset_1d_real(edge_id, "Veloc2", veloc2)
    call read_dset_1d_real(edge_id, "Veloc3", veloc3)
    call read_dset_1d_real(edge_id, "Displ", displ)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll

        do i = 1,ngll-2
            Tdomain%sEdge(n)%Veloc(i,0) = veloc(idx1+0)
            Tdomain%sEdge(n)%Veloc(i,1) = veloc(idx1+1)
            Tdomain%sEdge(n)%Veloc(i,2) = veloc(idx1+2)
            idx1 = idx1 + 3
            if ( .not. Tdomain%sEdge(n)%PML ) then
                Tdomain%sEdge(n)%Displ(i,0) = displ(idx2+0)
                Tdomain%sEdge(n)%Displ(i,1) = displ(idx2+1)
                Tdomain%sEdge(n)%Displ(i,2) = displ(idx2+2)
                idx2 = idx2 + 3
            else
                Tdomain%sEdge(n)%Veloc1(i,0) = veloc1(idx3+0)
                Tdomain%sEdge(n)%Veloc1(i,1) = veloc1(idx3+1)
                Tdomain%sEdge(n)%Veloc1(i,2) = veloc1(idx3+2)
                Tdomain%sEdge(n)%Veloc2(i,0) = veloc2(idx3+0)
                Tdomain%sEdge(n)%Veloc2(i,1) = veloc2(idx3+1)
                Tdomain%sEdge(n)%Veloc2(i,2) = veloc2(idx3+2)
                Tdomain%sEdge(n)%Veloc3(i,0) = veloc3(idx3+0)
                Tdomain%sEdge(n)%Veloc3(i,1) = veloc3(idx3+1)
                Tdomain%sEdge(n)%Veloc3(i,2) = veloc3(idx3+2)
                idx3 = idx3 + 3
            end if
        end do
    end do
    deallocate(veloc)
    deallocate(veloc1)
    deallocate(veloc2)
    deallocate(veloc3)
    deallocate(displ)
end subroutine read_Edges

subroutine read_Vertices(Tdomain, vertex_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: vertex_id
    real*8, allocatable, dimension(:) :: veloc, veloc1, veloc2, veloc3, displ
    integer :: idx1, idx2, idx3
    integer :: n


    call read_dset_1d_real(vertex_id, "Veloc", veloc)
    call read_dset_1d_real(vertex_id, "Veloc1", veloc1)
    call read_dset_1d_real(vertex_id, "Veloc2", veloc2)
    call read_dset_1d_real(vertex_id, "Veloc3", veloc3)
    call read_dset_1d_real(vertex_id, "Displ", displ)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    do n = 0,Tdomain%n_vertex-1
        Tdomain%sVertex(n)%Veloc(0) = veloc(idx1+0)
        Tdomain%sVertex(n)%Veloc(1) = veloc(idx1+1)
        Tdomain%sVertex(n)%Veloc(2) = veloc(idx1+2)
        idx1 = idx1 + 3
        if ( .not. Tdomain%sVertex(n)%PML ) then
            Tdomain%sVertex(n)%Displ(0) = displ(idx2+0)
            Tdomain%sVertex(n)%Displ(1) = displ(idx2+1)
            Tdomain%sVertex(n)%Displ(2) = displ(idx2+2)
            idx2 = idx2 + 3
        else
            Tdomain%sVertex(n)%Veloc1(0) = veloc1(idx3+0)
            Tdomain%sVertex(n)%Veloc1(1) = veloc1(idx3+1)
            Tdomain%sVertex(n)%Veloc1(2) = veloc1(idx3+2)
            Tdomain%sVertex(n)%Veloc2(0) = veloc2(idx3+0)
            Tdomain%sVertex(n)%Veloc2(1) = veloc2(idx3+1)
            Tdomain%sVertex(n)%Veloc2(2) = veloc2(idx3+2)
            Tdomain%sVertex(n)%Veloc3(0) = veloc3(idx3+0)
            Tdomain%sVertex(n)%Veloc3(1) = veloc3(idx3+1)
            Tdomain%sVertex(n)%Veloc3(2) = veloc3(idx3+2)
            idx3 = idx3 + 3
        end if
    end do
    deallocate(veloc)
    deallocate(veloc1)
    deallocate(veloc2)
    deallocate(veloc3)
    deallocate(displ)
end subroutine read_Vertices

subroutine read_restart (Tdomain,rg, isort)
    use HDF5
    use sem_hdf5
    use sdomain
    use semdatafiles
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer,intent (inout)::isort
    integer, intent (IN) :: rg
    character (len=MAX_FILE_SIZE) :: fnamef

    ! Vars sem
    real :: rtime, dtmin
    integer :: it

    ! HDF5 Variables
    integer(HID_T) :: fid, dset_id, elem_id, face_id, edge_id, vertex_id
    integer :: hdferr

    call init_hdf5()
    call init_restart(Tdomain, rg, isort, fnamef)


    write(*,*) "OPENING RESTART FILE:", trim(fnamef)
    call h5fopen_f(fnamef, H5F_ACC_RDONLY_F, fid, hdferr)

    call h5gopen_f(fid, 'Elements', elem_id, hdferr)
    call h5gopen_f(fid, 'Faces', face_id, hdferr)
    call h5gopen_f(fid, 'Edges', edge_id, hdferr)
    call h5gopen_f(fid, 'Vertices', vertex_id, hdferr)

    call read_attr_real(fid, "rtime", rtime)
    call read_attr_real(fid, "dtmin", dtmin)
    call read_attr_int(fid, "iteration", it)
    call read_attr_int(fid, "isort", isort)
    Tdomain%TimeD%rtime = rtime
    Tdomain%TimeD%dtmin = dtmin
    Tdomain%TimeD%NtimeMin = it

    Tdomain%TimeD%prot_m0 = Tdomain%TimeD%NtimeMin !!init des iterations de protection
    Tdomain%TimeD%prot_m1 = Tdomain%TimeD%NtimeMin

    !preparation pour le pas de temps suivant (absent de la version initiale de Sem3d)
    Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin
    Tdomain%TimeD%NtimeMin=Tdomain%TimeD%NtimeMin+1

    if (rg == 0) then
        write (*,'(A40,I8,A1,f10.6)') "SEM : REPRISE a iteration et tps:", it," ",rtime
    endif

    call read_Veloc(Tdomain, elem_id)
    call read_EpsilonVol(Tdomain, elem_id)
    call read_EpsilonDev(Tdomain, elem_id)
    call read_Stress(Tdomain, elem_id)
    call read_Faces(Tdomain, face_id)
    call read_Edges(Tdomain, edge_id)
    call read_Vertices(Tdomain, vertex_id)

    call h5gclose_f(elem_id, hdferr)
    call h5gclose_f(face_id, hdferr)
    call h5gclose_f(edge_id, hdferr)
    call h5gclose_f(vertex_id, hdferr)
    call h5fclose_f(fid, hdferr)

    call clean_prot(Tdomain, rg)
    return
end subroutine read_restart

#else

subroutine read_restart (Tdomain,rg, isort)
    use sdomain
    use semdatafiles
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer,intent (inout)::isort
    integer, intent (IN) :: rg
    character (len=MAX_FILE_SIZE) :: fnamef

    !  complement de sauvegarde pour le partie facteur de qualite Qp et Qs
    integer :: n_solid, i_sls

    ! local variables
    integer :: n, ngllx, nglly, ngllz, i, j, k, ngll, ngll1, ngll2
    call init_restart(Tdomain, rg, isort, fnamef)

    open (61, file=fnamef, status="unknown", form="formatted")
    read(61,*) Tdomain%TimeD%rtime, Tdomain%TimeD%dtmin
    read(61,*) Tdomain%TimeD%NtimeMin,isort !la version initiale de Sem3d comportait slt rtime, Ntimemin

    Tdomain%TimeD%prot_m0 = Tdomain%TimeD%NtimeMin !!GSa 11/2009 !!init des iterations de protection
    Tdomain%TimeD%prot_m1 = Tdomain%TimeD%NtimeMin !!GSa 11/2009

    !preparation pour le pas de temps suivant (absent de la version initiale de Sem3d)
    Tdomain%TimeD%rtime = Tdomain%TimeD%rtime + Tdomain%TimeD%dtmin
    Tdomain%TimeD%NtimeMin=Tdomain%TimeD%NtimeMin+1

    ! Save Fields for Elements
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz

        n_solid = Tdomain%n_sls

        if ( .not. Tdomain%specel(n)%PML ) then
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        read(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
                        read(61,*) Tdomain%specel(n)%Displ(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Displ(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Displ(i,j,k,2)
                        if ( n_solid > 0 ) then
                            if (Tdomain%aniso) then
                            else
                                read(61,*)  Tdomain%specel(n)%epsilonvol_ (i,j,k)
                                do i_sls = 0,n_solid-1
                                    read(61,*)  Tdomain%specel(n)%R_vol_ (i_sls,i,j,k)
                                enddo
                            endif
                            do i_sls = 0,n_solid-1
                                read(61,*)  Tdomain%specel(n)%R_xx_ (i_sls,i,j,k)
                                read(61,*)  Tdomain%specel(n)%R_yy_ (i_sls,i,j,k)
                                read(61,*)  Tdomain%specel(n)%R_xy_ (i_sls,i,j,k)
                                read(61,*)  Tdomain%specel(n)%R_xz_ (i_sls,i,j,k)
                                read(61,*)  Tdomain%specel(n)%R_yz_ (i_sls,i,j,k)
                            enddo
                            read(61,*) Tdomain%specel(n)%epsilondev_xx_ (i,j,k)
                            read(61,*) Tdomain%specel(n)%epsilondev_yy_ (i,j,k)
                            read(61,*) Tdomain%specel(n)%epsilondev_xy_ (i,j,k)
                            read(61,*) Tdomain%specel(n)%epsilondev_xz_ (i,j,k)
                            read(61,*) Tdomain%specel(n)%epsilondev_yz_ (i,j,k)
                        endif
                    enddo
                enddo
            enddo
        else
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        read(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
                        read(61,*) Tdomain%specel(n)%Veloc1(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Veloc1(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Veloc1(i,j,k,2)
                        read(61,*) Tdomain%specel(n)%Veloc2(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Veloc2(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Veloc2(i,j,k,2)
                        read(61,*) Tdomain%specel(n)%Veloc3(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Veloc3(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Veloc3(i,j,k,2)
                    enddo
                enddo
            enddo
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        read(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,2)
                        read(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,2)
                        read(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,2)
                        read(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,2)
                        read(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,0)
                        read(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,1)
                        read(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,2)
                    enddo
                enddo
            enddo
        endif
    enddo

    ! Save Fields for Faces
    do n = 0,Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
        if (.not. Tdomain%sFace(n)%PML ) then
            do j = 1,ngll2-2
                do i = 1,ngll1-2
                    read(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
                    read(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
                    read(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
                    read(61,*) Tdomain%sFace(n)%Displ(i,j,0)
                    read(61,*) Tdomain%sFace(n)%Displ(i,j,1)
                    read(61,*) Tdomain%sFace(n)%Displ(i,j,2)
                enddo
            enddo
        else
            do j = 1,ngll2-2
                do i = 1,ngll1-2
                    read(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
                    read(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
                    read(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
                    read(61,*) Tdomain%sFace(n)%Veloc1(i,j,0)
                    read(61,*) Tdomain%sFace(n)%Veloc1(i,j,1)
                    read(61,*) Tdomain%sFace(n)%Veloc1(i,j,2)
                    read(61,*) Tdomain%sFace(n)%Veloc2(i,j,0)
                    read(61,*) Tdomain%sFace(n)%Veloc2(i,j,1)
                    read(61,*) Tdomain%sFace(n)%Veloc2(i,j,2)
                    read(61,*) Tdomain%sFace(n)%Veloc3(i,j,0)
                    read(61,*) Tdomain%sFace(n)%Veloc3(i,j,1)
                    read(61,*) Tdomain%sFace(n)%Veloc3(i,j,2)
                enddo
            enddo
        endif
    enddo

    ! Save Fields for Edges
    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll
        if (.not. Tdomain%sEdge(n)%PML ) then
            do i = 1,ngll-2
                read(61,*) Tdomain%sEdge(n)%Veloc(i,0)
                read(61,*) Tdomain%sEdge(n)%Veloc(i,1)
                read(61,*) Tdomain%sEdge(n)%Veloc(i,2)
                read(61,*) Tdomain%sEdge(n)%Displ(i,0)
                read(61,*) Tdomain%sEdge(n)%Displ(i,1)
                read(61,*) Tdomain%sEdge(n)%Displ(i,2)
            enddo
        else
            do i = 1,ngll-2
                read(61,*) Tdomain%sEdge(n)%Veloc(i,0)
                read(61,*) Tdomain%sEdge(n)%Veloc(i,1)
                read(61,*) Tdomain%sEdge(n)%Veloc(i,2)
                read(61,*) Tdomain%sEdge(n)%Veloc1(i,0)
                read(61,*) Tdomain%sEdge(n)%Veloc1(i,1)
                read(61,*) Tdomain%sEdge(n)%Veloc1(i,2)
                read(61,*) Tdomain%sEdge(n)%Veloc2(i,0)
                read(61,*) Tdomain%sEdge(n)%Veloc2(i,1)
                read(61,*) Tdomain%sEdge(n)%Veloc2(i,2)
                read(61,*) Tdomain%sEdge(n)%Veloc3(i,0)
                read(61,*) Tdomain%sEdge(n)%Veloc3(i,1)
                read(61,*) Tdomain%sEdge(n)%Veloc3(i,2)
            enddo
        endif
    enddo

    ! Save Fields for Vertices
    do n = 0,Tdomain%n_vertex-1
        if (.not. Tdomain%sVertex(n)%PML ) then
            do i = 0,2
                read(61,*) Tdomain%sVertex(n)%Veloc(i)
                read(61,*) Tdomain%sVertex(n)%Displ(i)
            enddo
        else
            do i = 0,2
                read(61,*) Tdomain%sVertex(n)%Veloc(i)
                read(61,*) Tdomain%sVertex(n)%Veloc1(i)
                read(61,*) Tdomain%sVertex(n)%Veloc2(i)
                read(61,*) Tdomain%sVertex(n)%Veloc3(i)
            enddo
        endif
    enddo
    close(61)

    call clean_prot(Tdomain, rg)

    return
end subroutine read_restart
#endif
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
