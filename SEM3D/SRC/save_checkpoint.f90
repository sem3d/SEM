!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file save_checkpoint.f90
!! \brief Gère la protection de Sem3d
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine init_protection(Tdomain, it, prot_file)
    use sdomain
    use semdatafiles
    use mpi
    use sem_c_config, only : sem_mkdir
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: it
    integer :: rg
    character (len=MAX_FILE_SIZE), INTENT(OUT) :: prot_file
    character (len=MAX_FILE_SIZE) :: dir_prot, dir_prot_prev, times_file
    character (len=MAX_FILE_SIZE) :: dir_traces, dir_prot_traces
    character (len=MAX_FILE_SIZE) :: commande
    integer :: ierr

    rg = Tdomain%rank

    call semname_protection_iter_rank_file(it,rg,prot_file)

    call MPI_Barrier(Tdomain%communicateur, ierr)

    ! recherche et destruction au fur et a mesure des anciennes prots
    if (rg == 0) then

        call semname_protection_iter_dir(it,dir_prot)
        ! creation du repertoire data/sem/Protection_<it> (par tous les procs)
        ierr = sem_mkdir(trim(adjustl(dir_prot)))

        Tdomain%TimeD%prot_m2 = Tdomain%TimeD%prot_m1
        Tdomain%TimeD%prot_m1 = Tdomain%TimeD%prot_m0
        Tdomain%TimeD%prot_m0 = it

        if (Tdomain%TimeD%prot_m2>0) then
            call semname_protection_iter_dir(Tdomain%TimeD%prot_m2, dir_prot_prev)
            commande="rm -Rf "//trim(adjustl(dir_prot_prev))
            call system(trim(commande))        ! suppression de l avant avant derniere protection
        endif

        ! copie du fichier temps.dat dans le rep de protection
        call semname_results_temps_sem(times_file)
        commande="cp "//trim(adjustl(times_file))//" "//trim(adjustl(dir_prot))
        call system(trim(commande))

        ! copie du repertoire des sorties capteurs sem dans le rep de protection
        call semname_dir_capteurs(dir_traces)
        call semname_protection_iter_dir_capteurs(it,dir_prot_traces)
        commande="cp -R "//trim(adjustl(dir_traces))//" "//dir_prot_traces
        call system(trim(commande))
    endif
    call MPI_Barrier(Tdomain%communicateur, ierr)

end subroutine init_protection

subroutine compute_save_offsets(Tdomain, offset)
    use constants
    use sdomain
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(kind=4), intent(out), dimension (12) :: offset
    integer :: n,ngllx,nglly,ngllz,ngll,ngll2,n_solid

    ! Calcul des offsets de positions dans les tableaux
    n_solid = Tdomain%n_sls
    offset(1:12)=0
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz
        ngll = (ngllx-2)*(nglly-2)*(ngllz-2)
        ngll2 = ngllx*nglly*ngllz

        select case (Tdomain%specel(n)%domain)
        case (DM_SOLID)
            ! pour Veloc : 1
            offset(1) = offset(1) + ngll*3
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
                    offset(4) = offset(4) + ngll2
                    ! pour R_vol : 5
                    offset(5) = offset(5) + ngll2*n_solid
                endif
                ! pour R_xx, R_yy, R_xy, R_xz, R_yz : 6
                offset(6) = offset(6) + ngll2*n_solid
                ! pour epsilondev_xx, yy, xy, xz, yz : 7
                offset(7) = offset(7) + ngll2
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
        case (DM_SOLID_PML)
            ! pour Veloc : 1
            offset(1) = offset(1) + ngll*3
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
            offset(8) = offset(8) + ngll2*18
        case (DM_FLUID)
            ! pour VelPhi : 9
            offset(9) = offset(9) + ngll
            ! VelPhi1, VelPhi2, VelPhi3 : 10
            offset(10) = offset(10) + 0
            ! pour Phi : 11
            offset(11) = offset(11) + ngll
            ! pour Veloc : 12
            offset(12) = offset(12) + 0
        case (DM_FLUID_PML)
            ! pour VelPhi : 9
            offset(9) = offset(9) + ngll
            ! VelPhi1, VelPhi2, VelPhi3 : 10
            offset(10) = offset(10) + ngll
            ! pour Phi : 11
            offset(11) = offset(11) + 0
            ! pour Veloc : 12
            offset(12) = offset(12) + ngll2*9
        end select
    end do
end subroutine compute_save_offsets

subroutine write_EpsilonVol(Tdomain, nmax, elem_id)
    use constants, only : DM_SOLID
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
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 0,ngllz-1
            do j = 0,nglly-1
                do i = 0,ngllx-1
                    if (n_solid>0) then
                        if (Tdomain%aniso) then
                        else
                            if (idx.gt.nmax) then
                                write(*,*) "Erreur fatale sauvegarde des protections"
                                stop 1
                            end if
                            data(idx+0) = Tdomain%specel(n)%sl%epsilonvol_(i,j,k)
                            idx = idx + 1
                        end if
                    end if
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_EpsilonVol

subroutine write_Rvol(Tdomain, nmax, elem_id)
    use constants, only : DM_SOLID
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
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if


    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 0,ngllz-1
            do j = 0,nglly-1
                do i = 0,ngllx-1
                    if (n_solid>0) then
                        if (Tdomain%aniso) then
                        else
                            if (idx.gt.nmax) then
                                write(*,*) "Erreur fatale sauvegarde des protections"
                                stop 1
                            end if
                            do i_sls = 0, n_solid-1
                                data(idx+i_sls) = Tdomain%specel(n)%sl%R_vol_(i_sls,i,j,k)
                            end do
                            idx = idx + n_solid
                        end if
                    end if
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Rvol

subroutine write_Rxyz(Tdomain, nmax, elem_id)
    use constants, only : DM_SOLID
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
    if(nmax <= 0)then
        call h5dclose_f(dset_id_xx, hdferr)
        call h5dclose_f(dset_id_yy, hdferr)
        call h5dclose_f(dset_id_xy, hdferr)
        call h5dclose_f(dset_id_xz, hdferr)
        call h5dclose_f(dset_id_yz, hdferr)
        return
    end if


    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 0,ngllz-1
            do j = 0,nglly-1
                do i = 0,ngllx-1
                    if (n_solid>0) then
                        if (idx.gt.nmax) then
                            write(*,*) "Erreur fatale sauvegarde des protections"
                            stop 1
                        end if
                        do i_sls=0, n_solid-1
                            data_xx(idx+i_sls) = Tdomain%specel(n)%sl%R_xx_(i_sls,i,j,k)
                            data_yy(idx+i_sls) = Tdomain%specel(n)%sl%R_yy_(i_sls,i,j,k)
                            data_xy(idx+i_sls) = Tdomain%specel(n)%sl%R_xy_(i_sls,i,j,k)
                            data_xz(idx+i_sls) = Tdomain%specel(n)%sl%R_xz_(i_sls,i,j,k)
                            data_yz(idx+i_sls) = Tdomain%specel(n)%sl%R_yz_(i_sls,i,j,k)
                        end do
                        idx = idx + n_solid
                    end if
                enddo
            enddo
        enddo
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
    use constants, only : DM_SOLID
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
    if(nmax <= 0)then
        call h5dclose_f(dset_id_xx, hdferr)
        call h5dclose_f(dset_id_yy, hdferr)
        call h5dclose_f(dset_id_xy, hdferr)
        call h5dclose_f(dset_id_xz, hdferr)
        call h5dclose_f(dset_id_yz, hdferr)
        return
    end if
    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 0,ngllz-1
            do j = 0,nglly-1
                do i = 0,ngllx-1
                    if (n_solid>0) then
                        if (idx.gt.nmax) then
                            write(*,*) "Erreur fatale sauvegarde des protections"
                            stop 1
                        end if
                        data_xx(idx) = Tdomain%specel(n)%sl%epsilondev_xx_(i,j,k)
                        data_yy(idx) = Tdomain%specel(n)%sl%epsilondev_yy_(i,j,k)
                        data_xy(idx) = Tdomain%specel(n)%sl%epsilondev_xy_(i,j,k)
                        data_xz(idx) = Tdomain%specel(n)%sl%epsilondev_xz_(i,j,k)
                        data_yz(idx) = Tdomain%specel(n)%sl%epsilondev_yz_(i,j,k)
                        idx = idx + 1
                    end if
                enddo
            enddo
        enddo
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
    use constants, only : DM_SOLID_PML
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
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if

    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID_PML) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 0,ngllz-1
            do j = 0,nglly-1
                do i = 0,ngllx-1
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
                    data(idx+ 0) = Tdomain%specel(n)%slpml%Diagonal_Stress1(i,j,k,0)
                    data(idx+ 1) = Tdomain%specel(n)%slpml%Diagonal_Stress1(i,j,k,1)
                    data(idx+ 2) = Tdomain%specel(n)%slpml%Diagonal_Stress1(i,j,k,2)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%specel(n)%slpml%Diagonal_Stress2(i,j,k,0)
                    data(idx+ 1) = Tdomain%specel(n)%slpml%Diagonal_Stress2(i,j,k,1)
                    data(idx+ 2) = Tdomain%specel(n)%slpml%Diagonal_Stress2(i,j,k,2)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%specel(n)%slpml%Diagonal_Stress3(i,j,k,0)
                    data(idx+ 1) = Tdomain%specel(n)%slpml%Diagonal_Stress3(i,j,k,1)
                    data(idx+ 2) = Tdomain%specel(n)%slpml%Diagonal_Stress3(i,j,k,2)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%specel(n)%slpml%Residual_Stress1(i,j,k,0)
                    data(idx+ 1) = Tdomain%specel(n)%slpml%Residual_Stress1(i,j,k,1)
                    data(idx+ 2) = Tdomain%specel(n)%slpml%Residual_Stress1(i,j,k,2)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%specel(n)%slpml%Residual_Stress2(i,j,k,0)
                    data(idx+ 1) = Tdomain%specel(n)%slpml%Residual_Stress2(i,j,k,1)
                    data(idx+ 2) = Tdomain%specel(n)%slpml%Residual_Stress2(i,j,k,2)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%specel(n)%slpml%Residual_Stress3(i,j,k,0)
                    data(idx+ 1) = Tdomain%specel(n)%slpml%Residual_Stress3(i,j,k,1)
                    data(idx+ 2) = Tdomain%specel(n)%slpml%Residual_Stress3(i,j,k,2)
                    idx = idx + 3
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Stress

subroutine write_Veloc_Fluid_PML(Tdomain, nmax, elem_id)
    use constants, only : DM_FLUID_PML
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

    call create_dset(elem_id, "Veloc_Fl_PML", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if

    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_FLUID_PML) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 0,ngllz-1
            do j = 0,nglly-1
                do i = 0,ngllx-1
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
                    data(idx+ 0) = Tdomain%specel(n)%flpml%Veloc1(i,j,k,0)
                    data(idx+ 1) = Tdomain%specel(n)%flpml%Veloc1(i,j,k,1)
                    data(idx+ 2) = Tdomain%specel(n)%flpml%Veloc1(i,j,k,2)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%specel(n)%flpml%Veloc2(i,j,k,0)
                    data(idx+ 1) = Tdomain%specel(n)%flpml%Veloc2(i,j,k,1)
                    data(idx+ 2) = Tdomain%specel(n)%flpml%Veloc2(i,j,k,2)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%specel(n)%flpml%Veloc3(i,j,k,0)
                    data(idx+ 1) = Tdomain%specel(n)%flpml%Veloc3(i,j,k,1)
                    data(idx+ 2) = Tdomain%specel(n)%flpml%Veloc3(i,j,k,2)
                    idx = idx + 3
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Veloc_Fluid_PML


subroutine save_checkpoint (Tdomain, rtime, it, dtmin, isort)
    use sdomain
    use HDF5
    use sem_hdf5
    use semdatafiles
    implicit none

    type (domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: it, isort
    real, intent (IN) :: rtime, dtmin
    character (len=MAX_FILE_SIZE) :: fnamef
    ! HDF5 stuff
    integer :: hdferr
    integer(HID_T) :: fid, elem_id
    integer :: nelem
    integer(kind=4), dimension (12) :: offset
    integer(HSIZE_T), dimension(2) :: off_dims
    integer :: rg

    rg = Tdomain%rank

    if (rg == 0) then
        write (*,'(A44,I8,A1,f10.6)') "--> SEM : protection at iteration and time :",it," ",rtime
    endif
    call init_protection(Tdomain, it, fnamef)

    nelem = Tdomain%n_elem

    off_dims(1) = nelem+1
    off_dims(2) = 12

    call init_hdf5()
    call compute_save_offsets(Tdomain, offset)

    call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

    if (hdferr/=0) then
        write(*,*) "Detected HDF error at open:", fnamef, hdferr
        stop "Error writing HDF file"
    end if

    ! ifort doesn't care, but gfortran complains that the last integer should be 8 bytes
    call h5gcreate_f(fid, 'Elements', elem_id, hdferr, 0_SIZE_T)

    call write_attr_real(fid, "rtime", rtime)
    call write_attr_real(fid, "dtmin", dtmin)
    call write_attr_int(fid, "iteration", it)
    call write_attr_int(fid, "isort", isort)
    if (Tdomain%ngll_s.gt.0) then
        call write_dataset(elem_id, "sl_Veloc", Tdomain%champs0%Veloc)
        call write_dataset(elem_id, "sl_Displ", Tdomain%champs0%Depla)
    end if
    if (Tdomain%ngll_f.gt.0) then
        call write_dataset(elem_id, "fl_VelPhi", Tdomain%champs0%VelPhi)
        call write_dataset(elem_id, "fl_Phi", Tdomain%champs0%Phi)
    end if
    if (Tdomain%ngll_pmls.gt.0) then
        call write_dataset(elem_id, "spml_Veloc", Tdomain%champs0%VelocPML)
    end if
    if (Tdomain%ngll_pmlf.gt.0) then
        call write_dataset(elem_id, "fpml_VelPhi", Tdomain%champs0%fpml_VelPhi)
    end if
    call write_EpsilonVol(Tdomain, offset(4), elem_id)
    call write_EpsilonDev(Tdomain, offset(7), elem_id)
    call write_Stress(Tdomain, offset(8), elem_id)
    call write_Veloc_Fluid_PML(Tdomain, offset(12), elem_id)
    call write_Rvol(Tdomain, offset(5), elem_id)
    call write_Rxyz(Tdomain, offset(6), elem_id)

    call h5gclose_f(elem_id, hdferr)
    call h5fclose_f(fid, hdferr)
end subroutine save_checkpoint

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
