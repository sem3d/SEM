!>
!! \file save_checkpoint.f90
!! \brief Gère la protection de Sem3d
!! \author
!! \version 1.0
!! \date
!!
!<

#ifdef MKA3D
subroutine init_protection(Tdomain, it, rg, fnamef)
    use sdomain
    use semdatafiles
    use mpi

    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: it, rg
    character (len=MAX_FILE_SIZE), INTENT(OUT) :: fnamef
    character (len=MAX_FILE_SIZE) :: fnamer, fnamec
    character (len=MAX_FILE_SIZE) :: commande
    character (len=MAX_FILE_SIZE) :: commande_lg
    integer :: ierr

    call semname_couplage_iter(it,rg,fnamef)
    call semname_couplage_iterr(it,fnamer)

    ! creation du repertoire data/sem/Protection_<it> (par tous les procs)
    commande="mkdir -p "//trim(fnamer)

    ! recherche et destruction au fur et a mesure des anciennes prots
    if (rg == 0) then

        call system(trim(commande))

        Tdomain%TimeD%prot_m2 = Tdomain%TimeD%prot_m1
        Tdomain%TimeD%prot_m1 = Tdomain%TimeD%prot_m0
        Tdomain%TimeD%prot_m0 = it

        if (Tdomain%TimeD%prot_m2>0) then
            call semname_couplage_iterr(Tdomain%TimeD%prot_m2, fnamec)
            commande="rm -Rf "//trim(adjustl(fnamec))
            !print *,'commande suppression ',commande
            call system(trim(commande))        ! suppression de l avant avant derniere protection
        endif


        ! copie du fichier temps.dat dans le rep de protection
        !!  commande="cp ./Resultats/sem/temps.dat "//trim(fnamer) !!initial
        call semname_save_checkpoint_cp3(fnamer,fnamec)
        commande_lg="mkdir -p "//trim(adjustl(fnamec))
        !print *,'commande cp ',commande_lg
        call system(trim(commande_lg))

        ! copie du repertoire des sorties capteurs sem dans le rep de protection
        call semname_save_checkpoint_cp2(fnamer,fnamec)
        commande_lg="mkdir -p "//trim(adjustl(fnamec))
        !print *,'commande cp2 ',commande_lg
        call system(trim(commande_lg))
    endif
    call MPI_Barrier(Tdomain%communicateur, ierr)

end subroutine init_protection
#else
subroutine init_protection(Tdomain, it, rg, fnamef)
    use sdomain
    use semdatafiles

    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: it, rg
    character (len=MAX_FILE_SIZE), INTENT(OUT) :: fnamef

    call semname_save_checkpoint_rank(rg,fnamef)
end subroutine init_protection
#endif

#ifdef USE_HDF5

subroutine compute_save_offsets(Tdomain, offset, offset_f, offset_e, offset_v)
    use sdomain
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(kind=4), intent(out), dimension (8) :: offset
    integer(kind=4), intent(out), dimension (3) :: offset_f
    integer(kind=4), intent(out), dimension (3) :: offset_e
    integer(kind=4), intent(out), dimension (3) :: offset_v
    integer :: n,ngllx,nglly,ngllz,ngll,ngll1,ngll2,n_solid

    ! Calcul des offsets de positions dans les tableaux
    n_solid = Tdomain%n_sls
    offset(1:8)=0
    offset_f(1:3)=0
    offset_e(1:3)=0
    offset_v(1:3)=0
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz
        ngll = (ngllx-2)*(nglly-2)*(ngllz-2)
        ngll2 = ngllx*nglly*ngllz

        ! pour Veloc : 1
        offset(1) = offset(1) + ngll*3
        if ( .not. Tdomain%specel(n)%PML ) then
            ! Veloc1, Veloc2, Veloc3 : 2
            offset(2) = offset(2) + 0
            ! pour Displ : 3
            offset(3) = offset(3) + ngll*3

            if ( n_solid > 0 ) then
                if (Tdomain%aniso) then
                    ! pour EpsilonVol : 4
                    offset(4) = offset(4) + 0
                    ! pour R_vol : 5
                    offset(5) = offset(5) + 0
                else
                    ! pour EpsilonVol : 4
                    offset(4) = offset(4) + ngll
                    ! pour R_vol : 5
                    offset(5) = offset(5) + ngll*n_solid
                endif
                ! pour R_xx, R_yy, R_xy, R_xz, R_yz : 6
                offset(6) = offset(6) + ngll*n_solid
                ! pour epsilondev_xx, yy, xy, xz, yz : 7
                offset(7) = offset(7) + ngll
            else
                ! pour EpsilonVol : 4
                offset(4) = offset(4) + 0
                ! pour R_vol : 5
                offset(5) = offset(5) + 0
                ! pour R_xx, R_yy, R_xy, R_xz, R_yz : 6
                offset(6) = offset(6) + 0
                ! pour epsilondev_xx, yy, xy, xz, yz : 7
                offset(7) = offset(7) + 0
            end if
            ! pour Stress : 8
            offset(8) = offset(8) + 0
        else
            ! Veloc1, Veloc2, Veloc3 : 2
            offset(2) = offset(2) + ngll*3
            ! pour Displ : 3
            offset(3) = offset(3) + 0
            ! pour EpsilonVol : 4
            offset(4) = offset(4) + 0
            ! pour R_vol : 5
            offset(5) = offset(5) + 0
            ! pour R_xx, R_yy, R_xy, R_xz, R_yz : 6
            offset(6) = offset(6) + 0
            ! pour epsilondev_xx, yy, xy, xz, yz : 7
            offset(7) = offset(7) + 0
            ! pour Stress : 8
            offset(8) = offset(8) + ngll2*15
        end if
    end do
    ! Save Fields for Faces
    do n = 0,Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
        ngll = (ngll1-2)*(ngll2-2)
        if (.not. Tdomain%sFace(n)%PML ) then
            offset_f(1) = offset_f(1) + 3*ngll
            offset_f(2) = offset_f(2) + 3*ngll
            offset_f(3) = offset_f(3) + 0
        else
            offset_f(1) = offset_f(1) + 3*ngll
            offset_f(2) = offset_f(2) + 0
            offset_f(3) = offset_f(3) + 3*ngll
        end if
    end do
    ! Save Fields for Edges
    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll-2
        if (.not. Tdomain%sEdge(n)%PML ) then
            offset_e(1) = offset_e(1) + 3*ngll
            offset_e(2) = offset_e(2) + 3*ngll
            offset_e(3) = offset_e(3) + 0
        else
            offset_e(1) = offset_e(1) + 3*ngll
            offset_e(2) = offset_e(2) + 0
            offset_e(3) = offset_e(3) + 3*ngll
        end if
    enddo

    ! Save Fields for Vertices
    do n = 0,Tdomain%n_vertex-1
        if (.not. Tdomain%sVertex(n)%PML ) then
            offset_v(1) = offset_v(1) + 3
            offset_v(2) = offset_v(2) + 3
            offset_v(3) = offset_v(3) + 0
        else
            offset_v(1) = offset_v(1) + 3
            offset_v(2) = offset_v(2) + 0
            offset_v(3) = offset_v(3) + 3
        end if
    enddo

    write(*,*) "Offsets (C):", offset
    write(*,*) "Offsets (F):", offset_f
    write(*,*) "Offsets (E):", offset_e
    write(*,*) "Offsets (V):", offset_v
end subroutine compute_save_offsets


subroutine write_Veloc(Tdomain, nmax, elem_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id
    integer(kind=4), intent(IN) :: nmax

    integer :: n,ngllx,nglly,ngllz,idx,i,j,k,hdferr
    real(kind=8), dimension(1:nmax) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Veloc", H5T_IEEE_F64LE, nmax, dset_id)
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 1,ngllz-2
            do j = 1,nglly-2
                do i = 1,ngllx-2
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
                    data(idx+0) = Tdomain%specel(n)%Veloc(i,j,k,0)
                    data(idx+1) = Tdomain%specel(n)%Veloc(i,j,k,1)
                    data(idx+2) = Tdomain%specel(n)%Veloc(i,j,k,2)
                    idx = idx + 3
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Veloc

subroutine write_Veloc123(Tdomain, nmax, elem_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id1, dset_id2, dset_id3
    integer(kind=4), intent(IN) :: nmax

    integer :: n,ngllx,nglly,ngllz,idx,i,j,k,hdferr
    real(kind=8), dimension(1:nmax) :: data1, data2, data3
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Veloc1", H5T_IEEE_F64LE, nmax, dset_id1)
    call create_dset(elem_id, "Veloc2", H5T_IEEE_F64LE, nmax, dset_id2)
    call create_dset(elem_id, "Veloc3", H5T_IEEE_F64LE, nmax, dset_id3)

    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if (Tdomain%specel(n)%PML ) then
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        if (idx.gt.nmax) then
                            write(*,*) "Erreur fatale sauvegarde des protections"
                            stop 1
                        end if
                        data1(idx+0) = Tdomain%specel(n)%Veloc1(i,j,k,0)
                        data1(idx+1) = Tdomain%specel(n)%Veloc1(i,j,k,1)
                        data1(idx+2) = Tdomain%specel(n)%Veloc1(i,j,k,2)
                        data2(idx+0) = Tdomain%specel(n)%Veloc2(i,j,k,0)
                        data2(idx+1) = Tdomain%specel(n)%Veloc2(i,j,k,1)
                        data2(idx+2) = Tdomain%specel(n)%Veloc2(i,j,k,2)
                        data3(idx+0) = Tdomain%specel(n)%Veloc3(i,j,k,0)
                        data3(idx+1) = Tdomain%specel(n)%Veloc3(i,j,k,1)
                        data3(idx+2) = Tdomain%specel(n)%Veloc3(i,j,k,2)
                        idx = idx + 3
                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id1, H5T_NATIVE_DOUBLE, data1, dims, hdferr)
    call h5dwrite_f(dset_id2, H5T_NATIVE_DOUBLE, data2, dims, hdferr)
    call h5dwrite_f(dset_id3, H5T_NATIVE_DOUBLE, data3, dims, hdferr)
    call h5dclose_f(dset_id1, hdferr)
    call h5dclose_f(dset_id2, hdferr)
    call h5dclose_f(dset_id3, hdferr)
end subroutine write_Veloc123

subroutine write_Disp(Tdomain, nmax, elem_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id
    integer(kind=4), intent(IN) :: nmax
    integer :: n,ngllx,nglly,ngllz,idx,i,j,k,hdferr
    real(kind=8), dimension(1:nmax) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Displ", H5T_IEEE_F64LE, nmax, dset_id)
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if ( .not. Tdomain%specel(n)%PML ) then
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        if (idx.gt.nmax) then
                            write(*,*) "Erreur fatale sauvegarde des protections"
                            stop 1
                        end if
                        data(idx+0) = Tdomain%specel(n)%Displ(i,j,k,0)
                        data(idx+1) = Tdomain%specel(n)%Displ(i,j,k,1)
                        data(idx+2) = Tdomain%specel(n)%Displ(i,j,k,2)
                        idx = idx + 3
                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Disp

subroutine write_EpsilonVol(Tdomain, nmax, elem_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id
    integer(kind=4), intent(IN) :: nmax
    integer :: n,ngllx,nglly,ngllz,idx,i,j,k,hdferr
    real(kind=8), dimension(1:nmax) :: data
    integer(HSIZE_T), dimension(1) :: dims
    integer :: n_solid

    call create_dset(elem_id, "EpsilonVol", H5T_IEEE_F64LE, nmax, dset_id)
    dims(1) = nmax
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
                            if (Tdomain%aniso) then
                            else
                                if (idx.gt.nmax) then
                                    write(*,*) "Erreur fatale sauvegarde des protections"
                                    stop 1
                                end if
                                data(idx+0) = Tdomain%specel(n)%epsilonvol_(i,j,k)
                                idx = idx + 1
                            end if
                        end if
                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_EpsilonVol

subroutine write_Rvol(Tdomain, nmax, elem_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id
    integer(kind=4), intent(IN) :: nmax
    integer :: n,ngllx,nglly,ngllz,idx,i,j,k,hdferr,i_sls
    real(kind=8), dimension(1:nmax) :: data
    integer(HSIZE_T), dimension(1) :: dims
    integer :: n_solid

    call create_dset(elem_id, "Rvol", H5T_IEEE_F64LE, nmax, dset_id)
    dims(1) = nmax
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
                            if (Tdomain%aniso) then
                            else
                                if (idx.gt.nmax) then
                                    write(*,*) "Erreur fatale sauvegarde des protections"
                                    stop 1
                                end if
                                do i_sls = 0, n_solid-1
                                    data(idx+i_sls) = Tdomain%specel(n)%R_vol_(i_sls,i,j,k)
                                end do
                                idx = idx + n_solid
                            end if
                        end if
                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Rvol

subroutine write_Rxyz(Tdomain, nmax, elem_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id_xx, dset_id_yy, dset_id_xy, dset_id_xz, dset_id_yz
    integer(kind=4), intent(IN) :: nmax
    integer :: n,ngllx,nglly,ngllz,idx,i,j,k,hdferr,i_sls
    real(kind=8), dimension(1:nmax) :: data_xx, data_yy, data_xy, data_xz, data_yz
    integer(HSIZE_T), dimension(1) :: dims
    integer :: n_solid

    call create_dset(elem_id, "R_xx", H5T_IEEE_F64LE, nmax, dset_id_xx)
    call create_dset(elem_id, "R_yy", H5T_IEEE_F64LE, nmax, dset_id_yy)
    call create_dset(elem_id, "R_xy", H5T_IEEE_F64LE, nmax, dset_id_xy)
    call create_dset(elem_id, "R_xz", H5T_IEEE_F64LE, nmax, dset_id_xz)
    call create_dset(elem_id, "R_yz", H5T_IEEE_F64LE, nmax, dset_id_yz)
    dims(1) = nmax
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
                            if (idx.gt.nmax) then
                                write(*,*) "Erreur fatale sauvegarde des protections"
                                stop 1
                            end if
                            do i_sls=0, n_solid-1
                                data_xx(idx+i_sls) = Tdomain%specel(n)%R_xx_(i_sls,i,j,k)
                                data_yy(idx+i_sls) = Tdomain%specel(n)%R_yy_(i_sls,i,j,k)
                                data_xy(idx+i_sls) = Tdomain%specel(n)%R_xy_(i_sls,i,j,k)
                                data_xz(idx+i_sls) = Tdomain%specel(n)%R_xz_(i_sls,i,j,k)
                                data_yz(idx+i_sls) = Tdomain%specel(n)%R_yz_(i_sls,i,j,k)
                            end do
                            idx = idx + n_solid
                        end if
                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id_xx, H5T_NATIVE_DOUBLE, data_xx, dims, hdferr)
    call h5dwrite_f(dset_id_yy, H5T_NATIVE_DOUBLE, data_yy, dims, hdferr)
    call h5dwrite_f(dset_id_xy, H5T_NATIVE_DOUBLE, data_xy, dims, hdferr)
    call h5dwrite_f(dset_id_xz, H5T_NATIVE_DOUBLE, data_xz, dims, hdferr)
    call h5dwrite_f(dset_id_yz, H5T_NATIVE_DOUBLE, data_yz, dims, hdferr)
    call h5dclose_f(dset_id_xx, hdferr)
    call h5dclose_f(dset_id_yy, hdferr)
    call h5dclose_f(dset_id_xy, hdferr)
    call h5dclose_f(dset_id_xz, hdferr)
    call h5dclose_f(dset_id_yz, hdferr)
end subroutine write_Rxyz

subroutine write_EpsilonDev(Tdomain, nmax, elem_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id_xx, dset_id_yy, dset_id_xy, dset_id_xz, dset_id_yz
    integer(kind=4), intent(IN) :: nmax
    integer :: n,ngllx,nglly,ngllz,idx,i,j,k,hdferr
    real(kind=8), dimension(1:nmax) :: data_xx, data_yy, data_xy, data_xz, data_yz
    integer(HSIZE_T), dimension(1) :: dims
    integer :: n_solid

    call create_dset(elem_id, "EpsilonDev_xx", H5T_IEEE_F64LE, nmax, dset_id_xx)
    call create_dset(elem_id, "EpsilonDev_yy", H5T_IEEE_F64LE, nmax, dset_id_yy)
    call create_dset(elem_id, "EpsilonDev_xy", H5T_IEEE_F64LE, nmax, dset_id_xy)
    call create_dset(elem_id, "EpsilonDev_xz", H5T_IEEE_F64LE, nmax, dset_id_xz)
    call create_dset(elem_id, "EpsilonDev_yz", H5T_IEEE_F64LE, nmax, dset_id_yz)
    dims(1) = nmax
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
                            if (idx.gt.nmax) then
                                write(*,*) "Erreur fatale sauvegarde des protections"
                                stop 1
                            end if
                            data_xx(idx) = Tdomain%specel(n)%epsilondev_xx_(i,j,k)
                            data_yy(idx) = Tdomain%specel(n)%epsilondev_yy_(i,j,k)
                            data_xy(idx) = Tdomain%specel(n)%epsilondev_xy_(i,j,k)
                            data_xz(idx) = Tdomain%specel(n)%epsilondev_xz_(i,j,k)
                            data_yz(idx) = Tdomain%specel(n)%epsilondev_yz_(i,j,k)
                            idx = idx + 1
                        end if
                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id_xx, H5T_NATIVE_DOUBLE, data_xx, dims, hdferr)
    call h5dwrite_f(dset_id_yy, H5T_NATIVE_DOUBLE, data_yy, dims, hdferr)
    call h5dwrite_f(dset_id_xy, H5T_NATIVE_DOUBLE, data_xy, dims, hdferr)
    call h5dwrite_f(dset_id_xz, H5T_NATIVE_DOUBLE, data_xz, dims, hdferr)
    call h5dwrite_f(dset_id_yz, H5T_NATIVE_DOUBLE, data_yz, dims, hdferr)
    call h5dclose_f(dset_id_xx, hdferr)
    call h5dclose_f(dset_id_yy, hdferr)
    call h5dclose_f(dset_id_xy, hdferr)
    call h5dclose_f(dset_id_xz, hdferr)
    call h5dclose_f(dset_id_yz, hdferr)
end subroutine write_EpsilonDev

subroutine write_Stress(Tdomain, nmax, elem_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id
    integer(kind=4), intent(IN) :: nmax

    integer :: n,ngllx,nglly,ngllz,idx,i,j,k,hdferr
    real(kind=8), dimension(1:nmax) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Stress", H5T_IEEE_F64LE, nmax, dset_id)

    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if (Tdomain%specel(n)%PML ) then
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        if (idx.gt.nmax) then
                            write(*,*) "Erreur fatale sauvegarde des protections"
                            stop 1
                        end if
                        data(idx+ 0) = Tdomain%specel(n)%Diagonal_Stress1(i,j,k,0)
                        data(idx+ 1) = Tdomain%specel(n)%Diagonal_Stress1(i,j,k,1)
                        data(idx+ 2) = Tdomain%specel(n)%Diagonal_Stress1(i,j,k,2)
                        idx = idx + 3
                        data(idx+ 0) = Tdomain%specel(n)%Diagonal_Stress2(i,j,k,0)
                        data(idx+ 1) = Tdomain%specel(n)%Diagonal_Stress2(i,j,k,1)
                        data(idx+ 2) = Tdomain%specel(n)%Diagonal_Stress2(i,j,k,2)
                        idx = idx + 3
                        data(idx+ 0) = Tdomain%specel(n)%Diagonal_Stress3(i,j,k,0)
                        data(idx+ 1) = Tdomain%specel(n)%Diagonal_Stress3(i,j,k,1)
                        data(idx+ 2) = Tdomain%specel(n)%Diagonal_Stress3(i,j,k,2)
                        idx = idx + 3
                        data(idx+ 0) = Tdomain%specel(n)%Residual_Stress1(i,j,k,0)
                        data(idx+ 1) = Tdomain%specel(n)%Residual_Stress1(i,j,k,1)
                        data(idx+ 2) = Tdomain%specel(n)%Residual_Stress1(i,j,k,2)
                        idx = idx + 3
                        data(idx+ 0) = Tdomain%specel(n)%Residual_Stress2(i,j,k,0)
                        data(idx+ 1) = Tdomain%specel(n)%Residual_Stress2(i,j,k,1)
                        data(idx+ 2) = Tdomain%specel(n)%Residual_Stress2(i,j,k,2)
                        idx = idx + 3
                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Stress

subroutine write_Faces(Tdomain, offset_f, face_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: face_id
    integer(kind=4), dimension(3), intent(IN) :: offset_f

    integer(HID_T) :: veloc_id, displ_id, veloc1_id, veloc2_id, veloc3_id
    integer :: n,ngll1,ngll2,idx1,idx2,idx3,i,j,hdferr
    real(kind=8), dimension(1:offset_f(1)) :: veloc
    real(kind=8), dimension(1:offset_f(2)) :: displ
    real(kind=8), dimension(1:offset_f(3)) :: veloc1
    real(kind=8), dimension(1:offset_f(3)) :: veloc2
    real(kind=8), dimension(1:offset_f(3)) :: veloc3
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(face_id, "Veloc", H5T_IEEE_F64LE, offset_f(1), veloc_id)
    call create_dset(face_id, "Displ", H5T_IEEE_F64LE, offset_f(2), displ_id)
    call create_dset(face_id, "Veloc1", H5T_IEEE_F64LE, offset_f(3), veloc1_id)
    call create_dset(face_id, "Veloc2", H5T_IEEE_F64LE, offset_f(3), veloc2_id)
    call create_dset(face_id, "Veloc3", H5T_IEEE_F64LE, offset_f(3), veloc3_id)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    do n = 0,Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1
        ngll2 = Tdomain%sFace(n)%ngll2

        do j = 1,ngll2-2
            do i = 1,ngll1-2
                if (idx1.gt.offset_f(1)) then
                    write(*,*) "Erreur fatale sauvegarde des protections"
                    stop 1
                end if
                veloc(idx1+0) = Tdomain%sFace(n)%Veloc(i,j,0)
                veloc(idx1+1) = Tdomain%sFace(n)%Veloc(i,j,1)
                veloc(idx1+2) = Tdomain%sFace(n)%Veloc(i,j,2)
                idx1 = idx1 + 3
                if (.not. Tdomain%sFace(n)%PML ) then
                    displ(idx2+0) = Tdomain%sFace(n)%Displ(i,j,0)
                    displ(idx2+1) = Tdomain%sFace(n)%Displ(i,j,1)
                    displ(idx2+2) = Tdomain%sFace(n)%Displ(i,j,2)
                    idx2 = idx2 + 3
                else
                    veloc1(idx3+0) = Tdomain%sFace(n)%Veloc1(i,j,0)
                    veloc1(idx3+1) = Tdomain%sFace(n)%Veloc1(i,j,1)
                    veloc1(idx3+2) = Tdomain%sFace(n)%Veloc1(i,j,2)
                    veloc2(idx3+0) = Tdomain%sFace(n)%Veloc2(i,j,0)
                    veloc2(idx3+1) = Tdomain%sFace(n)%Veloc2(i,j,1)
                    veloc2(idx3+2) = Tdomain%sFace(n)%Veloc2(i,j,2)
                    veloc3(idx3+0) = Tdomain%sFace(n)%Veloc3(i,j,0)
                    veloc3(idx3+1) = Tdomain%sFace(n)%Veloc3(i,j,1)
                    veloc3(idx3+2) = Tdomain%sFace(n)%Veloc3(i,j,2)
                    idx3 = idx3 + 3
                end if
            enddo
        enddo
    enddo
    dims(1) = offset_f(1)
    call h5dwrite_f(veloc_id, H5T_NATIVE_DOUBLE, veloc, dims, hdferr)
    dims(1) = offset_f(2)
    call h5dwrite_f(displ_id, H5T_NATIVE_DOUBLE, displ, dims, hdferr)
    dims(1) = offset_f(3)
    call h5dwrite_f(veloc1_id, H5T_NATIVE_DOUBLE, veloc1, dims, hdferr)
    call h5dwrite_f(veloc2_id, H5T_NATIVE_DOUBLE, veloc2, dims, hdferr)
    call h5dwrite_f(veloc3_id, H5T_NATIVE_DOUBLE, veloc3, dims, hdferr)

    call h5dclose_f(veloc_id, hdferr)
    call h5dclose_f(displ_id, hdferr)
    call h5dclose_f(veloc1_id, hdferr)
    call h5dclose_f(veloc2_id, hdferr)
    call h5dclose_f(veloc3_id, hdferr)

end subroutine write_Faces

subroutine write_Edges(Tdomain, offset_e, edge_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: edge_id
    integer(kind=4), dimension(3), intent(IN) :: offset_e

    integer(HID_T) :: veloc_id, displ_id, veloc1_id, veloc2_id, veloc3_id
    integer :: n,ngll,idx1,idx2,idx3,i,j,hdferr
    real(kind=8), dimension(1:offset_e(1)) :: veloc
    real(kind=8), dimension(1:offset_e(2)) :: displ
    real(kind=8), dimension(1:offset_e(3)) :: veloc1
    real(kind=8), dimension(1:offset_e(3)) :: veloc2
    real(kind=8), dimension(1:offset_e(3)) :: veloc3
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(edge_id, "Veloc", H5T_IEEE_F64LE, offset_e(1), veloc_id)
    call create_dset(edge_id, "Displ", H5T_IEEE_F64LE, offset_e(2), displ_id)
    call create_dset(edge_id, "Veloc1", H5T_IEEE_F64LE, offset_e(3), veloc1_id)
    call create_dset(edge_id, "Veloc2", H5T_IEEE_F64LE, offset_e(3), veloc2_id)
    call create_dset(edge_id, "Veloc3", H5T_IEEE_F64LE, offset_e(3), veloc3_id)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll
        do i = 1,ngll-2
            if (idx1.gt.offset_e(1)) then
                write(*,*) "Erreur fatale sauvegarde des protections"
                stop 1
            end if
            if (.not. Tdomain%sEdge(n)%PML ) then
                veloc(idx1+0) = Tdomain%sEdge(n)%Veloc(i,0)
                veloc(idx1+1) = Tdomain%sEdge(n)%Veloc(i,1)
                veloc(idx1+2) = Tdomain%sEdge(n)%Veloc(i,2)
                idx1 = idx1 + 3
                displ(idx2+0) = Tdomain%sEdge(n)%Displ(i,0)
                displ(idx2+1) = Tdomain%sEdge(n)%Displ(i,1)
                displ(idx2+2) = Tdomain%sEdge(n)%Displ(i,2)
                idx2 = idx2 + 3
            else
                veloc(idx1+0) = Tdomain%sEdge(n)%Veloc(i,0)
                veloc(idx1+1) = Tdomain%sEdge(n)%Veloc(i,1)
                veloc(idx1+2) = Tdomain%sEdge(n)%Veloc(i,2)
                idx1 = idx1 + 3
                veloc1(idx3+0) = Tdomain%sEdge(n)%Veloc1(i,0)
                veloc1(idx3+1) = Tdomain%sEdge(n)%Veloc1(i,1)
                veloc1(idx3+2) = Tdomain%sEdge(n)%Veloc1(i,2)
                veloc2(idx3+0) = Tdomain%sEdge(n)%Veloc2(i,0)
                veloc2(idx3+1) = Tdomain%sEdge(n)%Veloc2(i,1)
                veloc2(idx3+2) = Tdomain%sEdge(n)%Veloc2(i,2)
                veloc3(idx3+0) = Tdomain%sEdge(n)%Veloc3(i,0)
                veloc3(idx3+1) = Tdomain%sEdge(n)%Veloc3(i,1)
                veloc3(idx3+2) = Tdomain%sEdge(n)%Veloc3(i,2)
                idx3 = idx3 + 3
            end if
        enddo
    enddo
    dims(1) = offset_e(1)
    call h5dwrite_f(veloc_id, H5T_NATIVE_DOUBLE, veloc, dims, hdferr)
    dims(1) = offset_e(2)
    call h5dwrite_f(displ_id, H5T_NATIVE_DOUBLE, displ, dims, hdferr)
    dims(1) = offset_e(3)
    call h5dwrite_f(veloc1_id, H5T_NATIVE_DOUBLE, veloc1, dims, hdferr)
    call h5dwrite_f(veloc2_id, H5T_NATIVE_DOUBLE, veloc2, dims, hdferr)
    call h5dwrite_f(veloc3_id, H5T_NATIVE_DOUBLE, veloc3, dims, hdferr)

    call h5dclose_f(veloc_id, hdferr)
    call h5dclose_f(displ_id, hdferr)
    call h5dclose_f(veloc1_id, hdferr)
    call h5dclose_f(veloc2_id, hdferr)
    call h5dclose_f(veloc3_id, hdferr)
end subroutine write_Edges

subroutine write_Vertices(Tdomain, offset_v, vertex_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: vertex_id
    integer(kind=4), dimension(3), intent(IN) :: offset_v

    integer(HID_T) :: veloc_id, displ_id, veloc1_id, veloc2_id, veloc3_id
    integer :: n,idx1,idx2,idx3,i,j,hdferr
    real(kind=8), dimension(1:offset_v(1)) :: veloc
    real(kind=8), dimension(1:offset_v(2)) :: displ
    real(kind=8), dimension(1:offset_v(3)) :: veloc1
    real(kind=8), dimension(1:offset_v(3)) :: veloc2
    real(kind=8), dimension(1:offset_v(3)) :: veloc3
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(vertex_id, "Veloc", H5T_IEEE_F64LE, offset_v(1), veloc_id)
    call create_dset(vertex_id, "Displ", H5T_IEEE_F64LE, offset_v(2), displ_id)
    call create_dset(vertex_id, "Veloc1", H5T_IEEE_F64LE, offset_v(3), veloc1_id)
    call create_dset(vertex_id, "Veloc2", H5T_IEEE_F64LE, offset_v(3), veloc2_id)
    call create_dset(vertex_id, "Veloc3", H5T_IEEE_F64LE, offset_v(3), veloc3_id)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    do n = 0,Tdomain%n_vertex-1
        if (idx1.gt.offset_v(1)) then
            write(*,*) "Erreur fatale sauvegarde des protections"
            stop 1
        end if
        if (.not. Tdomain%sVertex(n)%PML ) then
            veloc(idx1+0) = Tdomain%sVertex(n)%Veloc(0)
            veloc(idx1+1) = Tdomain%sVertex(n)%Veloc(1)
            veloc(idx1+2) = Tdomain%sVertex(n)%Veloc(2)
            idx1 = idx1 + 3
            displ(idx2+0) = Tdomain%sVertex(n)%Displ(0)
            displ(idx2+1) = Tdomain%sVertex(n)%Displ(1)
            displ(idx2+2) = Tdomain%sVertex(n)%Displ(2)
            idx2 = idx2 + 3
        else
            veloc(idx1+0) = Tdomain%sVertex(n)%Veloc(0)
            veloc(idx1+1) = Tdomain%sVertex(n)%Veloc(1)
            veloc(idx1+2) = Tdomain%sVertex(n)%Veloc(2)
            idx1 = idx1 + 3
            veloc1(idx3+0) = Tdomain%sVertex(n)%Veloc1(0)
            veloc1(idx3+1) = Tdomain%sVertex(n)%Veloc1(1)
            veloc1(idx3+2) = Tdomain%sVertex(n)%Veloc1(2)
            veloc2(idx3+0) = Tdomain%sVertex(n)%Veloc2(0)
            veloc2(idx3+1) = Tdomain%sVertex(n)%Veloc2(1)
            veloc2(idx3+2) = Tdomain%sVertex(n)%Veloc2(2)
            veloc3(idx3+0) = Tdomain%sVertex(n)%Veloc3(0)
            veloc3(idx3+1) = Tdomain%sVertex(n)%Veloc3(1)
            veloc3(idx3+2) = Tdomain%sVertex(n)%Veloc3(2)
            idx3 = idx3 + 3
        end if
    enddo
    dims(1) = offset_v(1)
    call h5dwrite_f(veloc_id, H5T_NATIVE_DOUBLE, veloc, dims, hdferr)
    dims(1) = offset_v(2)
    call h5dwrite_f(displ_id, H5T_NATIVE_DOUBLE, displ, dims, hdferr)
    dims(1) = offset_v(3)
    call h5dwrite_f(veloc1_id, H5T_NATIVE_DOUBLE, veloc1, dims, hdferr)
    call h5dwrite_f(veloc2_id, H5T_NATIVE_DOUBLE, veloc2, dims, hdferr)
    call h5dwrite_f(veloc3_id, H5T_NATIVE_DOUBLE, veloc3, dims, hdferr)

    call h5dclose_f(veloc_id, hdferr)
    call h5dclose_f(displ_id, hdferr)
    call h5dclose_f(veloc1_id, hdferr)
    call h5dclose_f(veloc2_id, hdferr)
    call h5dclose_f(veloc3_id, hdferr)

end subroutine write_Vertices

subroutine save_checkpoint (Tdomain, rtime, it, rg, dtmin, isort)
    use sdomain
    use HDF5
    use sem_hdf5
    use semdatafiles
    implicit none

    type (domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: it, rg, isort
    real, intent (IN) :: rtime, dtmin
    character (len=MAX_FILE_SIZE) :: fnamef, fnamer
    !  complement de sauvegarde pour le partie facteur de qualite Qp et Qs
    integer :: n_solid , i_sls
    ! HDF5 stuff
    integer :: hdferr
    integer(HID_T) :: fid, dset_id, elem_id, face_id, edge_id, vertex_id
    logical :: avail
    integer :: nelem, noffset
    integer :: size_vec, size_eps, size_epsaniso
    integer(kind=4), dimension (8) :: offset
    integer(kind=4), dimension (3) :: offset_f ! Veloc / Displ / (Veloc1,Veloc2,Veloc3)
    integer(kind=4), dimension (3) :: offset_e ! Veloc / Displ / (Veloc1,Veloc2,Veloc3)
    integer(kind=4), dimension (3) :: offset_v ! Veloc / Displ / (Veloc1,Veloc2,Veloc3)
    integer(HSIZE_T), dimension(2) :: off_dims

    if (rg == 0) then
        write (*,'(A40,I8,A1,f10.6)') "SEM : protection a iteration et tps:",it," ",rtime
    endif
    call init_protection(Tdomain, it, rg, fnamef)

    nelem = Tdomain%n_elem

    off_dims(1) = nelem+1
    off_dims(2) = 8

    call init_hdf5()

    call compute_save_offsets(Tdomain, offset, offset_f, offset_e, offset_v)

    call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

    call h5gcreate_f(fid, 'Elements', elem_id, hdferr, 0)
    call h5gcreate_f(fid, 'Faces', face_id, hdferr, 0)
    call h5gcreate_f(fid, 'Edges', edge_id, hdferr, 0)
    call h5gcreate_f(fid, 'Vertices', vertex_id, hdferr, 0)

    !  call create_dset_2d(fid, "Offsets", H5T_STD_I32LE, nelem+1, 8, dset_id)
    !  call h5dwrite_f(dset_id, H5T_STD_I32LE, offset, off_dims, hdferr)
    !  call h5dclose_f(dset_id, hdferr)

    call write_attr_real(fid, "rtime", rtime)
    call write_attr_real(fid, "dtmin", dtmin)
    call write_attr_int(fid, "iteration", it)
    call write_attr_int(fid, "isort", isort)
    call write_Veloc(Tdomain, offset(1), elem_id)
    call write_Veloc123(Tdomain, offset(2), elem_id)
    call write_Disp(Tdomain, offset(3), elem_id)
    call write_EpsilonVol(Tdomain, offset(4), elem_id)
    call write_Rvol(Tdomain, offset(5), elem_id)
    call write_Rxyz(Tdomain, offset(6), elem_id)
    call write_EpsilonDev(Tdomain, offset(7), elem_id)
    call write_Stress(Tdomain, offset(8), elem_id)

    call write_Faces(Tdomain, offset_f, face_id)
    call write_Edges(Tdomain, offset_e, edge_id)
    call write_Vertices(Tdomain, offset_v, vertex_id)

    call h5gclose_f(elem_id, hdferr)
    call h5gclose_f(face_id, hdferr)
    call h5gclose_f(edge_id, hdferr)
    call h5gclose_f(vertex_id, hdferr)
    call h5fclose_f(fid, hdferr)
end subroutine save_checkpoint

#else
subroutine save_checkpoint (Tdomain, rtime, it, rg, dtmin, icountc)
    use sdomain
    use semdatafiles
    implicit none

    type (domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: it, rg, icountc
    real, intent (IN) :: rtime, dtmin
    ! local variables
    integer :: n,ngllx,nglly,ngllz,i,j,k,ngll,ngll1,ngll2
    character (len=MAX_FILE_SIZE) :: fnamef, fnamer

    !  complement de sauvegarde pour le partie facteur de qualite Qp et Qs
    integer :: n_solid , i_sls

    if (rg == 0) then
        write (*,'(A40,I8,A1,f10.6)') "SEM : protection a iteration et tps:",it," ",rtime
        !write(6,*) 'save_checkpoint isort vaut ',icountc !!verif Gsa
    endif

    call init_protection(Tdomain, it, rg, fnamef)

    open (61, file=fnamef, status="unknown", form="formatted",err=101)
    write(61,*) rtime, dtmin !initialement Sem3d ecrivait ici uniquement rtime et i
    write(61,*) it, icountc !initialement Sem3d ecrivait ici uniquement rtime et it
101 continue
    ! Save Fields for Elements
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz

        n_solid = Tdomain%n_sls

        if ( .not. Tdomain%specel(n)%PML ) then
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        write(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
                        write(61,*) Tdomain%specel(n)%Displ(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Displ(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Displ(i,j,k,2)
                        if ( n_solid > 0 ) then
                            if (Tdomain%aniso) then
                            else
                                write(61,*)  Tdomain%specel(n)%epsilonvol_ (i,j,k)
                                do i_sls = 0,n_solid-1
                                    write(61,*)  Tdomain%specel(n)%R_vol_ (i_sls,i,j,k)
                                enddo
                            endif
                            do i_sls = 0,n_solid-1
                                write(61,*)  Tdomain%specel(n)%R_xx_ (i_sls,i,j,k)
                                write(61,*)  Tdomain%specel(n)%R_yy_ (i_sls,i,j,k)
                                write(61,*)  Tdomain%specel(n)%R_xy_ (i_sls,i,j,k)
                                write(61,*)  Tdomain%specel(n)%R_xz_ (i_sls,i,j,k)
                                write(61,*)  Tdomain%specel(n)%R_yz_ (i_sls,i,j,k)
                            enddo
                            write(61,*) Tdomain%specel(n)%epsilondev_xx_ (i,j,k)
                            write(61,*) Tdomain%specel(n)%epsilondev_yy_ (i,j,k)
                            write(61,*) Tdomain%specel(n)%epsilondev_xy_ (i,j,k)
                            write(61,*) Tdomain%specel(n)%epsilondev_xz_ (i,j,k)
                            write(61,*) Tdomain%specel(n)%epsilondev_yz_ (i,j,k)
                        endif
                    enddo
                enddo
            enddo
        else
            do k = 1,ngllz-2
                do j = 1,nglly-2
                    do i = 1,ngllx-2
                        write(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
                        write(61,*) Tdomain%specel(n)%Veloc1(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Veloc1(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Veloc1(i,j,k,2)
                        write(61,*) Tdomain%specel(n)%Veloc2(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Veloc2(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Veloc2(i,j,k,2)
                        write(61,*) Tdomain%specel(n)%Veloc3(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Veloc3(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Veloc3(i,j,k,2)
                    enddo
                enddo
            enddo
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        write(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,2)
                        write(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,2)
                        write(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,2)
                        write(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,2)
                        write(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,0)
                        write(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,1)
                        write(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,2)
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
                    write(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
                    write(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
                    write(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
                    write(61,*) Tdomain%sFace(n)%Displ(i,j,0)
                    write(61,*) Tdomain%sFace(n)%Displ(i,j,1)
                    write(61,*) Tdomain%sFace(n)%Displ(i,j,2)
                enddo
            enddo
        else
            do j = 1,ngll2-2
                do i = 1,ngll1-2
                    write(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
                    write(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
                    write(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
                    write(61,*) Tdomain%sFace(n)%Veloc1(i,j,0)
                    write(61,*) Tdomain%sFace(n)%Veloc1(i,j,1)
                    write(61,*) Tdomain%sFace(n)%Veloc1(i,j,2)
                    write(61,*) Tdomain%sFace(n)%Veloc2(i,j,0)
                    write(61,*) Tdomain%sFace(n)%Veloc2(i,j,1)
                    write(61,*) Tdomain%sFace(n)%Veloc2(i,j,2)
                    write(61,*) Tdomain%sFace(n)%Veloc3(i,j,0)
                    write(61,*) Tdomain%sFace(n)%Veloc3(i,j,1)
                    write(61,*) Tdomain%sFace(n)%Veloc3(i,j,2)
                enddo
            enddo
        endif
    enddo

    ! Save Fields for Edges
    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll
        if (.not. Tdomain%sEdge(n)%PML ) then
            do i = 1,ngll-2
                write(61,*) Tdomain%sEdge(n)%Veloc(i,0)
                write(61,*) Tdomain%sEdge(n)%Veloc(i,1)
                write(61,*) Tdomain%sEdge(n)%Veloc(i,2)
                write(61,*) Tdomain%sEdge(n)%Displ(i,0)
                write(61,*) Tdomain%sEdge(n)%Displ(i,1)
                write(61,*) Tdomain%sEdge(n)%Displ(i,2)
            enddo
        else
            do i = 1,ngll-2
                write(61,*) Tdomain%sEdge(n)%Veloc(i,0)
                write(61,*) Tdomain%sEdge(n)%Veloc(i,1)
                write(61,*) Tdomain%sEdge(n)%Veloc(i,2)
                write(61,*) Tdomain%sEdge(n)%Veloc1(i,0)
                write(61,*) Tdomain%sEdge(n)%Veloc1(i,1)
                write(61,*) Tdomain%sEdge(n)%Veloc1(i,2)
                write(61,*) Tdomain%sEdge(n)%Veloc2(i,0)
                write(61,*) Tdomain%sEdge(n)%Veloc2(i,1)
                write(61,*) Tdomain%sEdge(n)%Veloc2(i,2)
                write(61,*) Tdomain%sEdge(n)%Veloc3(i,0)
                write(61,*) Tdomain%sEdge(n)%Veloc3(i,1)
                write(61,*) Tdomain%sEdge(n)%Veloc3(i,2)
            enddo
        endif
    enddo

    ! Save Fields for Vertices
    do n = 0,Tdomain%n_vertex-1
        if (.not. Tdomain%sVertex(n)%PML ) then
            do i = 0,2
                write(61,*) Tdomain%sVertex(n)%Veloc(i)
                write(61,*) Tdomain%sVertex(n)%Displ(i)
            enddo
        else
            do i = 0,2
                write(61,*) Tdomain%sVertex(n)%Veloc(i)
                write(61,*) Tdomain%sVertex(n)%Veloc1(i)
                write(61,*) Tdomain%sVertex(n)%Veloc2(i)
                write(61,*) Tdomain%sVertex(n)%Veloc3(i)
            enddo
        endif
    enddo
    close(61)



    return
end subroutine save_checkpoint
#endif
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
