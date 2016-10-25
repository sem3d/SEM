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
    integer(kind=4), intent(out), dimension (13) :: offset
    integer :: n,ngllx,nglly,ngllz,ngll,ngll2,n_solid

    ! Calcul des offsets de positions dans les tableaux
    n_solid = Tdomain%n_sls
    offset(1:13)=0
    do n = 0,Tdomain%n_elem-1
        ngllx = 0; nglly = 0; ngllz = 0;
        select case (Tdomain%specel(n)%domain)
             case (DM_SOLID)
                 ngllx = Tdomain%sdom%ngll
                 nglly = Tdomain%sdom%ngll
                 ngllz = Tdomain%sdom%ngll
             case (DM_FLUID)
                 ngllx = Tdomain%fdom%ngll
                 nglly = Tdomain%fdom%ngll
                 ngllz = Tdomain%fdom%ngll
             case (DM_SOLID_PML)
                 ngllx = Tdomain%spmldom%ngll
                 nglly = Tdomain%spmldom%ngll
                 ngllz = Tdomain%spmldom%ngll
             case (DM_FLUID_PML)
                 ngllx = Tdomain%fpmldom%ngll
                 nglly = Tdomain%fpmldom%ngll
                 ngllz = Tdomain%fpmldom%ngll
             case default
                 stop "unknown domain"
        end select
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

            if (Tdomain%nl_flag) then
                ! for stress computed via nonlinear algorithm
                offset(9) = offset(9)+ngll2*6
                ! for strain computed via nonlinear algorithm
                offset(10) = offset(10)+ngll2*6
                ! for internal variables computed via nonlinear algorithm
                offset(11) = offset(11)+ngll2*6 ! back stress
                offset(12) = offset(12)+ngll2   ! radius
                ! for plastic strain
                offset(13) = offset(13)+ngll2*6
            else
                ! for stress computed via nonlinear algorithm
                offset(9) = offset(9)+0
                ! for strain computed via nonlinear algorithm
                offset(10) = offset(10)+0
                ! for internal variables computed via nonlinear algorithm
                offset(11) = offset(11)+0
                offset(12) = offset(12)+0
                ! for plastic strain
                offset(13) = offset(13)+0
            endif
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
        case default
            stop "unknown domain"
        end select
    end do
end subroutine compute_save_offsets

subroutine write_EpsilonVol(Tdomain, nmax, elem_id)
    use constants, only : DM_SOLID
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
#include "index.h"
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id
    integer(kind=4), intent(IN) :: nmax
    integer :: n,ngll,idx,i,j,k,hdferr,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data
    integer(HSIZE_T), dimension(1) :: dims
    integer :: n_solid

    call create_dset(elem_id, "EpsilonVol", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    allocate(data(1:nmax))
    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngll = Tdomain%sdom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (n_solid>0) then
                        if (Tdomain%aniso) then
                        else
                            if (idx.gt.nmax) then
                                write(*,*) "Erreur fatale sauvegarde des protections"
                                stop 1
                            end if
                            data(idx+0) = Tdomain%sdom%epsilonvol_(i,j,k,bnum,ee)
                            idx = idx + 1
                        end if
                    end if
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    deallocate(data)
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
    integer :: n,ngll,idx,i,j,k,hdferr,i_sls,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data
    integer(HSIZE_T), dimension(1) :: dims
    integer :: n_solid

    call create_dset(elem_id, "Rvol", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    allocate(data(1:nmax))

    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngll = Tdomain%sdom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (n_solid>0) then
                        if (Tdomain%aniso) then
                        else
                            if (idx.gt.nmax) then
                                write(*,*) "Erreur fatale sauvegarde des protections"
                                stop 1
                            end if
                            do i_sls = 0, n_solid-1
                                data(idx+i_sls) = Tdomain%sdom%R_vol_(i_sls,i,j,k,bnum,ee)
                            end do
                            idx = idx + n_solid
                        end if
                    end if
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    deallocate(data)
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
    integer :: n,ngll,idx,i,j,k,hdferr,i_sls,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data_xx, data_yy, data_xy, data_xz, data_yz
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
    allocate(data_xx(1:nmax))
    allocate(data_yy(1:nmax))
    allocate(data_xy(1:nmax))
    allocate(data_xz(1:nmax))
    allocate(data_yz(1:nmax))


    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngll = Tdomain%sdom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (n_solid>0) then
                        if (idx.gt.nmax) then
                            write(*,*) "Erreur fatale sauvegarde des protections"
                            stop 1
                        end if
                        do i_sls=0, n_solid-1
                            data_xx(idx+i_sls) = Tdomain%sdom%R_xx_(i_sls,i,j,k,bnum,ee)
                            data_yy(idx+i_sls) = Tdomain%sdom%R_yy_(i_sls,i,j,k,bnum,ee)
                            data_xy(idx+i_sls) = Tdomain%sdom%R_xy_(i_sls,i,j,k,bnum,ee)
                            data_xz(idx+i_sls) = Tdomain%sdom%R_xz_(i_sls,i,j,k,bnum,ee)
                            data_yz(idx+i_sls) = Tdomain%sdom%R_yz_(i_sls,i,j,k,bnum,ee)
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
    deallocate(data_xx,data_yy,data_xy,data_xz,data_yz)
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
    integer :: n,ngll,idx,i,j,k,hdferr,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data_xx, data_yy, data_xy, data_xz, data_yz
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
    allocate(data_xx(1:nmax))
    allocate(data_yy(1:nmax))
    allocate(data_xy(1:nmax))
    allocate(data_xz(1:nmax))
    allocate(data_yz(1:nmax))

    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngll = Tdomain%sdom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (n_solid>0) then
                        if (idx.gt.nmax) then
                            write(*,*) "Erreur fatale sauvegarde des protections"
                            stop 1
                        end if
                        data_xx(idx) = Tdomain%sdom%epsilondev_xx_(i,j,k,bnum,ee)
                        data_yy(idx) = Tdomain%sdom%epsilondev_yy_(i,j,k,bnum,ee)
                        data_xy(idx) = Tdomain%sdom%epsilondev_xy_(i,j,k,bnum,ee)
                        data_xz(idx) = Tdomain%sdom%epsilondev_xz_(i,j,k,bnum,ee)
                        data_yz(idx) = Tdomain%sdom%epsilondev_yz_(i,j,k,bnum,ee)
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
    deallocate(data_xx,data_yy,data_xy,data_xz,data_yz)
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
#include "index.h"
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(kind=4), intent(IN) :: nmax

#ifndef CPML
    integer(HID_T) :: dset_id
    integer :: n,ngll,idx,i,j,k,hdferr,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Stress", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    allocate(data(1:nmax))
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID_PML) cycle
        ngll = Tdomain%spmldom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
                    data(idx+ 0) = Tdomain%spmldom%Diagonal_Stress1_(i,j,k,0,bnum,ee)
                    data(idx+ 1) = Tdomain%spmldom%Diagonal_Stress1_(i,j,k,1,bnum,ee)
                    data(idx+ 2) = Tdomain%spmldom%Diagonal_Stress1_(i,j,k,2,bnum,ee)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%spmldom%Diagonal_Stress2_(i,j,k,0,bnum,ee)
                    data(idx+ 1) = Tdomain%spmldom%Diagonal_Stress2_(i,j,k,1,bnum,ee)
                    data(idx+ 2) = Tdomain%spmldom%Diagonal_Stress2_(i,j,k,2,bnum,ee)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%spmldom%Diagonal_Stress3_(i,j,k,0,bnum,ee)
                    data(idx+ 1) = Tdomain%spmldom%Diagonal_Stress3_(i,j,k,1,bnum,ee)
                    data(idx+ 2) = Tdomain%spmldom%Diagonal_Stress3_(i,j,k,2,bnum,ee)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%spmldom%Residual_Stress1_(i,j,k,0,bnum,ee)
                    data(idx+ 1) = Tdomain%spmldom%Residual_Stress1_(i,j,k,1,bnum,ee)
                    data(idx+ 2) = Tdomain%spmldom%Residual_Stress1_(i,j,k,2,bnum,ee)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%spmldom%Residual_Stress2_(i,j,k,0,bnum,ee)
                    data(idx+ 1) = Tdomain%spmldom%Residual_Stress2_(i,j,k,1,bnum,ee)
                    data(idx+ 2) = Tdomain%spmldom%Residual_Stress2_(i,j,k,2,bnum,ee)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%spmldom%Residual_Stress3_(i,j,k,0,bnum,ee)
                    data(idx+ 1) = Tdomain%spmldom%Residual_Stress3_(i,j,k,1,bnum,ee)
                    data(idx+ 2) = Tdomain%spmldom%Residual_Stress3_(i,j,k,2,bnum,ee)
                    idx = idx + 3
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    deallocate(data)
    call h5dclose_f(dset_id, hdferr)
#endif
end subroutine write_Stress

subroutine write_Stress_nl(Tdomain, nmax, elem_id)
    use constants, only : DM_SOLID
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
#include "index.h"
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(kind=4), intent(IN) :: nmax

#ifndef CPML
    integer(HID_T) :: dset_id
    integer :: n,ngll,idx,i,j,k,hdferr,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Stress_nl", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    allocate(data(1:nmax))
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngll = Tdomain%sdom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
                    data(idx+ 0) = Tdomain%sdom%stress_(0,i,j,k,bnum,ee)
                    data(idx+ 1) = Tdomain%sdom%stress_(1,i,j,k,bnum,ee)
                    data(idx+ 2) = Tdomain%sdom%stress_(2,i,j,k,bnum,ee)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%sdom%stress_(3,i,j,k,bnum,ee)
                    data(idx+ 1) = Tdomain%sdom%stress_(4,i,j,k,bnum,ee)
                    data(idx+ 2) = Tdomain%sdom%stress_(5,i,j,k,bnum,ee)
                    idx = idx + 3
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    deallocate(data)
    call h5dclose_f(dset_id, hdferr)
#endif
end subroutine write_Stress_nl
!
subroutine write_Strain_nl(Tdomain, nmax, elem_id)
    use constants, only : DM_SOLID
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
#include "index.h"
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(kind=4), intent(IN) :: nmax

#ifndef CPML
    integer(HID_T) :: dset_id
    integer :: n,ngll,idx,i,j,k,hdferr,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Strain_nl", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    allocate(data(1:nmax))
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngll = Tdomain%sdom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
                    data(idx+ 0) = Tdomain%sdom%strain_(0,i,j,k,bnum,ee)
                    data(idx+ 1) = Tdomain%sdom%strain_(1,i,j,k,bnum,ee)
                    data(idx+ 2) = Tdomain%sdom%strain_(2,i,j,k,bnum,ee)
                    idx = idx + 3
                    data(idx+ 0) = Tdomain%sdom%strain_(3,i,j,k,bnum,ee)
                    data(idx+ 1) = Tdomain%sdom%strain_(4,i,j,k,bnum,ee)
                    data(idx+ 2) = Tdomain%sdom%strain_(5,i,j,k,bnum,ee)
                    idx = idx + 3
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    deallocate(data)
    call h5dclose_f(dset_id, hdferr)
#endif
end subroutine write_Strain_nl
!
subroutine write_Pstrain_nl(Tdomain, nmax, elem_id)
    use constants, only : DM_SOLID
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
#include "index.h"
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(kind=4), intent(IN) :: nmax

#ifndef CPML
    integer(HID_T) :: dset_id
    integer :: n,ngll,idx,i,j,k,hdferr,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Pstrain_nl", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    allocate(data(1:nmax))
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngll = Tdomain%sdom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
                    data(idx+ 0) = Tdomain%sdom%plstrain_(0,i,j,k,bnum,ee)
                    data(idx+ 1) = Tdomain%sdom%plstrain_(1,i,j,k,bnum,ee)
                    data(idx+ 2) = Tdomain%sdom%plstrain_(2,i,j,k,bnum,ee)
                    idx = idx + 3                                 
                    data(idx+ 0) = Tdomain%sdom%plstrain_(3,i,j,k,bnum,ee)
                    data(idx+ 1) = Tdomain%sdom%plstrain_(4,i,j,k,bnum,ee)
                    data(idx+ 2) = Tdomain%sdom%plstrain_(5,i,j,k,bnum,ee)
                    idx = idx + 3
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    deallocate(data)
    call h5dclose_f(dset_id, hdferr)
#endif
end subroutine write_Pstrain_nl
!
subroutine write_Center_nl(Tdomain, nmax, elem_id)
    use constants, only : DM_SOLID
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
#include "index.h"
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(kind=4), intent(IN) :: nmax

#ifndef CPML
    integer(HID_T) :: dset_id
    integer :: n,ngll,idx,i,j,k,hdferr,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Center_nl", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    allocate(data(1:nmax))
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngll = Tdomain%sdom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
                    data(idx+ 0) = Tdomain%sdom%center_(0,i,j,k,bnum,ee)
                    data(idx+ 1) = Tdomain%sdom%center_(1,i,j,k,bnum,ee)
                    data(idx+ 2) = Tdomain%sdom%center_(2,i,j,k,bnum,ee)
                    idx = idx + 3                               
                    data(idx+ 0) = Tdomain%sdom%center_(3,i,j,k,bnum,ee)
                    data(idx+ 1) = Tdomain%sdom%center_(4,i,j,k,bnum,ee)
                    data(idx+ 2) = Tdomain%sdom%center_(5,i,j,k,bnum,ee)
                    idx = idx + 3
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    deallocate(data)
    call h5dclose_f(dset_id, hdferr)
#endif
end subroutine write_Center_nl
!
subroutine write_Radius_nl(Tdomain, nmax, elem_id)
    use constants, only : DM_SOLID
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
#include "index.h"
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(kind=4), intent(IN) :: nmax

#ifndef CPML
    integer(HID_T) :: dset_id
    integer :: n,ngll,idx,i,j,k,hdferr,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Radius_nl", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    allocate(data(1:nmax))
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_SOLID) cycle
        ngll = Tdomain%sdom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
                    data(idx+ 0) = Tdomain%sdom%radius_(i,j,k,bnum,ee)
                    idx = idx +1
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    deallocate(data)
    call h5dclose_f(dset_id, hdferr)
#endif
end subroutine write_Radius_nl
!
subroutine write_Veloc_Fluid_PML(Tdomain, nmax, elem_id)
    use constants, only : DM_FLUID_PML
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
#include "index.h"
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    integer(HID_T) :: dset_id
    integer(kind=4), intent(IN) :: nmax

    integer :: n,ngll,idx,i,j,k,hdferr,bnum,ee
    real(kind=8), allocatable, dimension(:) :: data
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(elem_id, "Veloc_Fl_PML", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    allocate(data(1:nmax))

    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%domain/=DM_FLUID_PML) cycle
        ngll = Tdomain%fpmldom%ngll
        bnum = Tdomain%specel(n)%lnum/VCHUNK
        ee = mod(Tdomain%specel(n)%lnum,VCHUNK)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (idx.gt.nmax) then
                        write(*,*) "Erreur fatale sauvegarde des protections"
                        stop 1
                    end if
#ifdef CPML
                    ! TODO
#else
                    data(idx+ 0) = Tdomain%fpmldom%PMLVeloc_(i,j,k,0,bnum,ee)
                    data(idx+ 1) = Tdomain%fpmldom%PMLVeloc_(i,j,k,1,bnum,ee)
                    data(idx+ 2) = Tdomain%fpmldom%PMLVeloc_(i,j,k,2,bnum,ee)
                    idx = idx + 3
#endif
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    deallocate(data)
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
    integer(kind=4), dimension (13) :: offset
    integer(HSIZE_T), dimension(2) :: off_dims
    integer :: rg

    rg = Tdomain%rank

    if (rg == 0) then
        write (*,'(A44,I8,A1,f10.6)') "--> SEM : protection at iteration and time :",it," ",rtime
    endif
    call init_protection(Tdomain, it, fnamef)

    nelem = Tdomain%n_elem

    off_dims(1) = nelem+1
    off_dims(2) = 13

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

    if (Tdomain%sdom%nglltot.gt.0) then
        call write_dataset(elem_id, "sl_Veloc", Tdomain%sdom%champs0%Veloc)
        call write_dataset(elem_id, "sl_Displ", Tdomain%sdom%champs0%Depla)
        ! nonlinear
        if (Tdomain%nl_flag) then
            call write_Stress_nl(Tdomain,offset(9),elem_id)
            call write_Strain_nl(Tdomain,offset(10),elem_id)
            call write_Center_nl(Tdomain,offset(11),elem_id)
            call write_Radius_nl(Tdomain,offset(12),elem_id)
            call write_Pstrain_nl(Tdomain,offset(13),elem_id)
        endif
    end if
    if (Tdomain%fdom%nglltot.gt.0) then
        call write_dataset(elem_id, "fl_VelPhi", Tdomain%fdom%champs0%VelPhi)
        call write_dataset(elem_id, "fl_Phi",    Tdomain%fdom%champs0%Phi)
    end if
    if (Tdomain%spmldom%nglltot.gt.0) then
#ifdef CPML
        call write_dataset(elem_id, "spml_Veloc", Tdomain%spmldom%champs0%Veloc)
#else
        call write_dataset(elem_id, "spml_Veloc", Tdomain%spmldom%champs0%VelocPML)
#endif
    end if
    if (Tdomain%fpmldom%nglltot.gt.0) then
#ifdef CPML
        call write_dataset(elem_id, "fpml_VelPhi", Tdomain%fpmldom%champs0%VelPhi)
#else
        call write_dataset(elem_id, "fpml_VelPhi", Tdomain%fpmldom%champs0%fpml_VelPhi)
#endif
        call write_Veloc_Fluid_PML(Tdomain, offset(12), elem_id)
    end if

    call write_EpsilonVol(Tdomain, offset(4), elem_id)
    call write_EpsilonDev(Tdomain, offset(7), elem_id)
    call write_Stress(Tdomain, offset(8), elem_id)

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
