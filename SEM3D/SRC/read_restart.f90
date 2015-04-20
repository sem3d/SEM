!>
!!\file read_restart.f90
!!\brief Contient la subroutine read_restart().
!!
!! Gère la reprise de Sem3d

subroutine read_EpsilonVol(Tdomain, elem_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    double precision, allocatable, dimension(:) :: epsilonvol, rvol, rxx, ryy, rxy, rxz, ryz
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

        if ( Tdomain%specel(n)%solid .and. ( .not. Tdomain%specel(n)%PML ) ) then
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        if (n_solid>0) then
                            if (Tdomain%aniso) then
                            else
                                Tdomain%specel(n)%sl%epsilonvol_(i,j,k) = epsilonVol(idx1)
                                idx1 = idx1 + 1
                                do i_sls = 0, n_solid-1
                                    Tdomain%specel(n)%sl%R_vol_(i_sls,i,j,k) = rvol(idx2+i_sls)
                                end do
                                idx2 = idx2 + n_solid
                            end if
                            do i_sls = 0, n_solid-1
                                Tdomain%specel(n)%sl%R_xx_(i_sls,i,j,k) = rxx(idx3+i_sls)
                                Tdomain%specel(n)%sl%R_yy_(i_sls,i,j,k) = ryy(idx3+i_sls)
                                Tdomain%specel(n)%sl%R_xy_(i_sls,i,j,k) = rxy(idx3+i_sls)
                                Tdomain%specel(n)%sl%R_xz_(i_sls,i,j,k) = rxz(idx3+i_sls)
                                Tdomain%specel(n)%sl%R_yz_(i_sls,i,j,k) = ryz(idx3+i_sls)
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
    double precision, allocatable, dimension(:) :: eps_dev_xx, eps_dev_yy, eps_dev_xy, eps_dev_xz, eps_dev_yz
    integer :: idx, ngllx, nglly, ngllz
    integer :: n, i, j, k
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

        if ( Tdomain%specel(n)%solid .and. ( .not. Tdomain%specel(n)%PML ) ) then
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        if (n_solid>0) then
                            Tdomain%specel(n)%sl%epsilondev_xx_(i,j,k) = eps_dev_xx(idx)
                            Tdomain%specel(n)%sl%epsilondev_yy_(i,j,k) = eps_dev_yy(idx)
                            Tdomain%specel(n)%sl%epsilondev_xy_(i,j,k) = eps_dev_xy(idx)
                            Tdomain%specel(n)%sl%epsilondev_xz_(i,j,k) = eps_dev_xz(idx)
                            Tdomain%specel(n)%sl%epsilondev_yz_(i,j,k) = eps_dev_yz(idx)
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
    double precision, allocatable, dimension(:) :: stress
    integer :: idx, ngllx, nglly, ngllz
    integer :: n, i, j, k
    integer :: n_solid

    call read_dset_1d_real(elem_id, "Stress", stress)
    idx = 1
    n_solid = Tdomain%n_sls
    do n = 0,Tdomain%n_elem-1
        if(.not. Tdomain%specel(n)%solid) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if (Tdomain%specel(n)%PML ) then
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        Tdomain%specel(n)%slpml%Diagonal_Stress1(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%slpml%Diagonal_Stress1(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%slpml%Diagonal_Stress1(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%slpml%Diagonal_Stress2(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%slpml%Diagonal_Stress2(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%slpml%Diagonal_Stress2(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%slpml%Diagonal_Stress3(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%slpml%Diagonal_Stress3(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%slpml%Diagonal_Stress3(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%slpml%Residual_Stress1(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%slpml%Residual_Stress1(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%slpml%Residual_Stress1(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%slpml%Residual_Stress2(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%slpml%Residual_Stress2(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%slpml%Residual_Stress2(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%slpml%Residual_Stress3(i,j,k,0) = stress(idx+0)
                        Tdomain%specel(n)%slpml%Residual_Stress3(i,j,k,1) = stress(idx+1)
                        Tdomain%specel(n)%slpml%Residual_Stress3(i,j,k,2) = stress(idx+2)
                        idx = idx + 3
                    end do
                end do
            end do
        end if
    end do
    deallocate(stress)
end subroutine read_Stress

subroutine read_Veloc_Fluid_PML(Tdomain, elem_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    double precision, allocatable, dimension(:) :: Veloc
    integer :: idx, ngllx, nglly, ngllz
    integer :: n, i, j, k
    integer :: id

    call read_dset_1d_real(elem_id, "Veloc_Fl_PML", veloc)
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%solid) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if (Tdomain%specel(n)%PML ) then
            do k = 0,ngllz-1
                do j = 0,nglly-1
                    do i = 0,ngllx-1
                        id = Tdomain%specel(n)%slpml%ISolPML(i,j,k)
                        Tdomain%specel(n)%flpml%Veloc1(i,j,k,0) = veloc(idx+0)
                        Tdomain%specel(n)%flpml%Veloc1(i,j,k,1) = veloc(idx+1)
                        Tdomain%specel(n)%flpml%Veloc1(i,j,k,2) = veloc(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%flpml%Veloc2(i,j,k,0) = veloc(idx+0)
                        Tdomain%specel(n)%flpml%Veloc2(i,j,k,1) = veloc(idx+1)
                        Tdomain%specel(n)%flpml%Veloc2(i,j,k,2) = veloc(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%flpml%Veloc3(i,j,k,0) = veloc(idx+0)
                        Tdomain%specel(n)%flpml%Veloc3(i,j,k,1) = veloc(idx+1)
                        Tdomain%specel(n)%flpml%Veloc3(i,j,k,2) = veloc(idx+2)
                        idx = idx + 3
                    end do
                end do
            end do
        end if
    end do
    deallocate(veloc)
end subroutine read_Veloc_Fluid_PML


subroutine read_restart (Tdomain,rg, isort)
    use HDF5
    use sem_hdf5
    use sdomain
    use semdatafiles
    use protrep
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer,intent (inout)::isort
    integer, intent (IN) :: rg
    character (len=MAX_FILE_SIZE) :: fnamef

    ! Vars sem
    real :: rtime, dtmin
    integer :: it

    ! HDF5 Variables
    integer(HID_T) :: fid, elem_id
    integer :: hdferr

    call init_hdf5()
    call init_restart(Tdomain%communicateur, rg, Tdomain%TimeD%iter_reprise, fnamef)

    call h5fopen_f(fnamef, H5F_ACC_RDONLY_F, fid, hdferr)

    call h5gopen_f(fid, 'Elements', elem_id, hdferr)
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

    if (Tdomain%ngll_s.gt.0) then
        call read_dataset(elem_id, "sl_Veloc", Tdomain%champs0%Veloc, ibase=0)
        call read_dataset(elem_id, "sl_Displ", Tdomain%champs0%Depla, ibase=0)
    end if
    if (Tdomain%ngll_f.gt.0) then
        call read_dataset(elem_id, "fl_VelPhi", Tdomain%champs0%VelPhi, ibase=0)
        call read_dataset(elem_id, "fl_Phi", Tdomain%champs0%Phi, ibase=0)
    end if
    if (Tdomain%ngll_pmls.gt.0) then
        call read_dataset(elem_id, "spml_Veloc", Tdomain%champs0%VelocPML, ibase=0)
    end if
    if (Tdomain%ngll_pmlf.gt.0) then
        call read_dataset(elem_id, "fpml_VelPhi", Tdomain%champs0%fpml_VelPhi, ibase=0)
    end if
    call read_EpsilonVol(Tdomain, elem_id)
    call read_EpsilonDev(Tdomain, elem_id)
    call read_Stress(Tdomain, elem_id)
    call read_Veloc_Fluid_PML(Tdomain, elem_id)
    call h5gclose_f(elem_id, hdferr)
    call h5fclose_f(fid, hdferr)

    call clean_prot(Tdomain%TimeD%prot_m0, rg)
    return
end subroutine read_restart

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
