!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module msnapshots
    use sdomain
    use hdf5
    use sem_hdf5
    use semdatafiles
    use mpi
    use constants
    use orientation
    use sem_c_bindings
    implicit none
contains


    subroutine compute_saved_elements(Tdomain, irenum, nnodes)
        type (domain), intent (INOUT):: Tdomain
        integer, allocatable, dimension(:), intent(out) :: irenum ! maps Iglobnum to file node number
        integer, intent(out) :: nnodes
        integer :: n, i, k, ngllx, ngllz, ig, gn, ne

        allocate(irenum(0:Tdomain%n_glob_points-1))

        irenum = -1
        ig = 0
        ne = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            ne = ne + 1
            do k = 0,ngllz - 1
                do i = 0,ngllx - 1
                    gn = Tdomain%specel(n)%Iglobnum(i,k)
                    if (irenum(gn) == -1) then
                        irenum(gn) = ig
                        ig = ig + 1
                    end if
                end do
            end do
        end do
        nnodes = ig
    end subroutine compute_saved_elements

    subroutine create_dir_sorties(Tdomain, rg, isort)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer,intent(in) :: isort, rg
        character(Len=MAX_FILE_SIZE) :: temp
        integer :: code

        if (rg==0) then
            call semname_snap_result_dir(isort, temp)
            code = sem_mkdir(temp)
        end if
        call mpi_barrier(Tdomain%communicateur, code)
    end subroutine create_dir_sorties

    !>
    !! Ecrit la geometrie pour les sorties dans un fichier HDF5
    !! Le format est destine a etre relu par paraview / xdmf
    !<
    subroutine write_snapshot_geom(Tdomain, rg)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: rg
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid
        integer :: hdferr, code
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes

        call init_hdf5()
        if (rg==0) then
            code = sem_mkdir(path_results)
        end if
        call mpi_barrier(Tdomain%communicateur, code)
        call semname_snap_geom_file(rg, fnamef)

        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

        call compute_saved_elements(Tdomain, irenum, nnodes)

        call write_global_nodes(Tdomain, fid, irenum, nnodes)

        call write_elem_connectivity(Tdomain, fid, irenum)

        call write_constant_fields(Tdomain, fid, irenum, nnodes)

        call write_material_fields(Tdomain, fid, irenum, nnodes)

        call h5fclose_f(fid, hdferr)

        if (rg==0) call write_master_xdmf(Tdomain)
    end subroutine write_snapshot_geom


    subroutine write_global_nodes(Tdomain, fid, irenum, nnodes)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, intent(in) :: nnodes
        integer, dimension(:), intent(in), allocatable :: irenum
        !
        real(fpp), dimension(:,:), allocatable :: nodes
        integer(HSIZE_T), dimension(2) :: dims
        integer :: n
        integer(HID_T) :: nodes_id
        integer :: hdferr

        allocate(nodes(0:2,0:nnodes-1))
        do n = 0, Tdomain%n_glob_points-1
            if (irenum(n)>=0) then
                nodes(0,irenum(n)) = Tdomain%GlobCoord(0,n)
                nodes(1,irenum(n)) = Tdomain%GlobCoord(1,n)
                nodes(2,irenum(n)) = 0.
            end if
        end do

        dims(1) = 3
        dims(2) = nnodes
        call create_dset_2d(fid, "Nodes", H5T_IEEE_F64LE, dims(1), dims(2), nodes_id)
        call h5dwrite_f(nodes_id, H5T_NATIVE_DOUBLE, nodes, dims, hdferr)
        call h5dclose_f(nodes_id, hdferr)
        deallocate(nodes)
    end subroutine write_global_nodes

    subroutine write_elem_connectivity(Tdomain, fid, irenum)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, dimension(:), intent(in), allocatable :: irenum
        !
        integer(HID_T) :: elem_id, mat_id, ngll_id, globnum_id, elem_num_id
        integer :: ngllx, ngllz
        integer(HSIZE_T), dimension(2) :: dims
        integer, dimension(:,:), allocatable :: data
        integer, dimension(:), allocatable :: mat, elem_num, iglobnum
        integer, dimension(3,0:Tdomain%n_elem-1) :: ngll
        integer :: count, ig, nglobnum
        integer :: i, k, n, nb_elem
        integer :: hdferr

        ! First we count the number of hexaedrons
        count = 0
        nglobnum = 0
        k = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngll(1,k) = Tdomain%specel(n)%ngllx
            ngll(2,k) = Tdomain%specel(n)%ngllz
            ngll(3,k) = 0
            ! Max number of global points (can count elem vertices twice)
            nglobnum = nglobnum + ngll(1,k)*ngll(2,k)
            ! Number of subelements
            count = count+(ngll(1,k)-1)*(ngll(2,k)-1)
            k = k + 1
        enddo
        nb_elem = k
        !! Nombre de points de gauss par element
        dims(1) = 3
        dims(2) = nb_elem
        call create_dset_2d(fid, "NGLL", H5T_STD_I16LE, dims(1), dims(2), ngll_id)
        call h5dwrite_f(ngll_id, H5T_NATIVE_INTEGER, ngll, dims, hdferr)
        call h5dclose_f(ngll_id, hdferr)

        Tdomain%n_quad = count
        allocate( data(1:4,0:count-1))
        allocate( mat(0:count-1))
        allocate( elem_num(0:count-1))

        call create_dset_2d(fid, "Elements", H5T_STD_I32LE, 4, count, elem_id)
        call create_dset(fid, "Material", H5T_STD_I32LE, count, mat_id)
        call create_dset(fid, "ElemID", H5T_STD_I32LE, count, elem_num_id)
        !call create_dset(fid, "ElemID", H5T_STD_I32LE, Tdomain%n_elem, elem_num_id)
        call create_dset(fid, "Iglobnum", H5T_STD_I32LE, nglobnum, globnum_id)

        allocate (iglobnum(nglobnum))
        dims(1) = 4
        dims(2) = count
        count = 0
        ig = 1
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            do k = 0,ngllz-2
                do i = 0,ngllx-2
                    data(1,count) = irenum(Tdomain%specel(n)%Iglobnum(i+0,k+0))
                    data(2,count) = irenum(Tdomain%specel(n)%Iglobnum(i+1,k+0))
                    data(3,count) = irenum(Tdomain%specel(n)%Iglobnum(i+1,k+1))
                    data(4,count) = irenum(Tdomain%specel(n)%Iglobnum(i+0,k+1))
                    mat(count) = Tdomain%specel(n)%mat_index
                    elem_num(count) = n
                    count=count+1
                end do
            end do
            do k = 0,ngllz - 1
                do i = 0,ngllx - 1
                    iglobnum(ig) = Tdomain%specel(n)%Iglobnum(i,k)
                    ig = ig + 1
                end do
            end do
        end do
        dims(1) = nglobnum
        dims(2) = 0
        call h5dwrite_f(globnum_id, H5T_NATIVE_INTEGER, iglobnum, dims, hdferr)
        deallocate(iglobnum)
        call h5dclose_f(globnum_id, hdferr)
        call h5dwrite_f(elem_id, H5T_NATIVE_INTEGER, data, dims, hdferr)
        call h5dclose_f(elem_id, hdferr)
        dims(1) = count
        dims(2) = 0
        call h5dwrite_f(mat_id, H5T_NATIVE_INTEGER, mat, dims, hdferr)
        call h5dclose_f(mat_id, hdferr)
        call h5dwrite_f(elem_num_id, H5T_NATIVE_INTEGER, elem_num, dims, hdferr)
        call h5dclose_f(elem_num_id, hdferr)
        deallocate(data)
        !deallocate(elem_num)
    end subroutine write_elem_connectivity

    subroutine save_field_h5(Tdomain, rg, isort)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer, intent(in) :: rg, isort
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer(HID_T) :: fid, displ_id, veloc_id, press_id, accel_id, rotat_id
        integer(HID_T) :: P_energy_id, K_energy_id, eps_vol_id
        integer(HID_T) :: eps_dev_xx_id, eps_dev_zz_id, eps_dev_xz_id
        integer(HID_T) :: sig_dev_xx_id, sig_dev_zz_id, sig_dev_xz_id
        integer(HID_T) :: dUxdx_id, dUxdz_id, dUzdx_id, dUzdz_id
        integer(HID_T) :: L_energy_id, S_energy_id, R_energy_id

        integer(HSIZE_T), dimension(2) :: dims
        integer(HSIZE_T), dimension(1) :: dimr
        real(fpp), dimension(:,:),allocatable :: displ, veloc, accel, field_rotat
        real(fpp), dimension(:), allocatable :: press, rotat
        real(fpp), dimension(:), allocatable :: P_energy, K_energy, eps_vol
        real(fpp), dimension(:), allocatable :: eps_dev_xx, eps_dev_zz, eps_dev_xz
        real(fpp), dimension(:), allocatable :: sig_dev_xx, sig_dev_zz, sig_dev_xz
        real(fpp), dimension(:), allocatable :: dUxdx, dUxdz, dUzdx, dUzdz
        real(fpp), dimension(:), allocatable :: L_energy, S_energy, R_energy

        real(fpp), dimension(:,:,:),allocatable :: field_displ, field_veloc, field_accel
        integer, dimension(:), allocatable :: valence
        integer :: hdferr
        integer :: ngllx, ngllz, idx, mat
        integer :: i, k, n
        integer, allocatable, dimension(:) :: irenum ! maps Iglobnum to file node number
        integer :: nnodes
        real(fpp), dimension(0:1,0:1) :: invgrad_ij
        real(fpp) :: dUx_dxi, dUx_deta, dUz_dxi, dUz_deta
        real(fpp) :: DXX, DXZ, DZX, DZZ
        real(fpp) :: eps_xx, eps_zz, eps_xz, eps_v
        real(fpp) :: sig_xx, sig_zz, sig_xz, sig_mean
        real(fpp), dimension(3) :: strain_v, stress_v

        call create_dir_sorties(Tdomain, rg, isort)
        call semname_snap_result_file(rg, isort, fnamef)

        call compute_saved_elements(Tdomain, irenum, nnodes)


        call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

        ! Dimensions pour champs vectoriels et scalaires
        dims(1) = 3
        dims(2) = nnodes
        dimr(1) = nnodes

        call create_dset_2d(fid, "displ", H5T_IEEE_F64LE, 3, nnodes, displ_id)
        call create_dset_2d(fid, "veloc", H5T_IEEE_F64LE, 3, nnodes, veloc_id)
        call create_dset_2d(fid, "accel", H5T_IEEE_F64LE, 3, nnodes, accel_id)
        call create_dset(fid, "rotat", H5T_IEEE_F64LE, nnodes, rotat_id)

        if (Tdomain%out_var_snap(OUT_PRESSION) == 1) then
            call create_dset(fid, "pressure", H5T_IEEE_F64LE, nnodes, press_id)
            allocate(press(0:nnodes-1)); press = 0.0_fpp
        endif
        if (Tdomain%out_var_snap(OUT_ENERGYP) == 1) then
            call create_dset(fid, "P_energy", H5T_IEEE_F64LE, nnodes, P_energy_id)
            allocate(P_energy(0:nnodes-1)); P_energy = 0.0_fpp
        endif
        if (Tdomain%out_var_snap(OUT_ENERGYK) == 1) then
            call create_dset(fid, "K_energy", H5T_IEEE_F64LE, nnodes, K_energy_id)
            allocate(K_energy(0:nnodes-1)); K_energy = 0.0_fpp
        endif
        if (Tdomain%out_var_snap(OUT_EPS_VOL) == 1) then
            call create_dset(fid, "eps_vol", H5T_IEEE_F64LE, nnodes, eps_vol_id)
            allocate(eps_vol(0:nnodes-1)); eps_vol = 0.0_fpp
        endif
        if (Tdomain%out_var_snap(OUT_EPS_DEV) == 1) then
            call create_dset(fid, "eps_dev_xx", H5T_IEEE_F64LE, nnodes, eps_dev_xx_id)
            call create_dset(fid, "eps_dev_zz", H5T_IEEE_F64LE, nnodes, eps_dev_zz_id)
            call create_dset(fid, "eps_dev_xz", H5T_IEEE_F64LE, nnodes, eps_dev_xz_id)
            allocate(eps_dev_xx(0:nnodes-1)); eps_dev_xx = 0.0_fpp
            allocate(eps_dev_zz(0:nnodes-1)); eps_dev_zz = 0.0_fpp
            allocate(eps_dev_xz(0:nnodes-1)); eps_dev_xz = 0.0_fpp
        endif
        if (Tdomain%out_var_snap(OUT_STRESS_DEV) == 1) then
            call create_dset(fid, "sig_dev_xx", H5T_IEEE_F64LE, nnodes, sig_dev_xx_id)
            call create_dset(fid, "sig_dev_zz", H5T_IEEE_F64LE, nnodes, sig_dev_zz_id)
            call create_dset(fid, "sig_dev_xz", H5T_IEEE_F64LE, nnodes, sig_dev_xz_id)
            allocate(sig_dev_xx(0:nnodes-1)); sig_dev_xx = 0.0_fpp
            allocate(sig_dev_zz(0:nnodes-1)); sig_dev_zz = 0.0_fpp
            allocate(sig_dev_xz(0:nnodes-1)); sig_dev_xz = 0.0_fpp
        endif
        if (Tdomain%out_var_snap(OUT_DUDX) == 1) then
            call create_dset(fid, "dUxdx", H5T_IEEE_F64LE, nnodes, dUxdx_id)
            call create_dset(fid, "dUxdz", H5T_IEEE_F64LE, nnodes, dUxdz_id)
            call create_dset(fid, "dUzdx", H5T_IEEE_F64LE, nnodes, dUzdx_id)
            call create_dset(fid, "dUzdz", H5T_IEEE_F64LE, nnodes, dUzdz_id)
            allocate(dUxdx(0:nnodes-1)); dUxdx = 0.0_fpp
            allocate(dUxdz(0:nnodes-1)); dUxdz = 0.0_fpp
            allocate(dUzdx(0:nnodes-1)); dUzdx = 0.0_fpp
            allocate(dUzdz(0:nnodes-1)); dUzdz = 0.0_fpp
        endif
        if (Tdomain%out_var_snap(OUT_ENERGYD) == 1) then
            call create_dset(fid, "L_energy", H5T_IEEE_F64LE, nnodes, L_energy_id)
            call create_dset(fid, "S_energy", H5T_IEEE_F64LE, nnodes, S_energy_id)
            call create_dset(fid, "R_energy", H5T_IEEE_F64LE, nnodes, R_energy_id)
            allocate(L_energy(0:nnodes-1)); L_energy = 0.0_fpp
            allocate(S_energy(0:nnodes-1)); S_energy = 0.0_fpp
            allocate(R_energy(0:nnodes-1)); R_energy = 0.0_fpp
        endif

        allocate(displ(0:2,0:nnodes-1))
        allocate(veloc(0:2,0:nnodes-1))
        allocate(accel(0:2,0:nnodes-1))
        allocate(rotat(0:nnodes-1))
        allocate(valence(0:nnodes-1))
        allocate(field_rotat(0:0,0:0))

        ngllx = 0
        ngllz = 0
        valence(:) = 0
        displ(:,:) = 0
        veloc(:,:) = 0
        accel(:,:) = 0
        rotat(:) = 0
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            if (ngllx /= Tdomain%specel(n)%ngllx .or. &
                ngllz /= Tdomain%specel(n)%ngllz) then
                ngllx = Tdomain%specel(n)%ngllx
                ngllz = Tdomain%specel(n)%ngllz
                if (allocated(field_displ)) deallocate(field_displ)
                if (allocated(field_veloc)) deallocate(field_veloc)
                if (allocated(field_accel)) deallocate(field_accel)
                if (allocated(field_rotat)) deallocate(field_rotat)
                allocate(field_displ(0:ngllx-1,0:ngllz-1,2))
                allocate(field_veloc(0:ngllx-1,0:ngllz-1,2))
                allocate(field_accel(0:ngllx-1,0:ngllz-1,2))
                allocate(field_rotat(0:ngllx-1,0:ngllz-1))
            endif

            call gather_elem_displ(Tdomain, n, field_displ)
            call gather_elem_veloc(Tdomain, n, field_veloc, .true.)
            call gather_elem_accel(Tdomain, n, field_accel)
            call compute_rotational(Tdomain,n,ngllx,ngllz,field_veloc,field_rotat)

            do k = 0,ngllz-1
                do i = 0,ngllx-1
                    idx = irenum(Tdomain%specel(n)%Iglobnum(i,k))
                    valence(idx) = valence(idx)+1
                    displ(0:1,idx) = field_displ(i,k,:)
                    veloc(0:1,idx) = veloc(0:1,idx)+field_veloc(i,k,:)
                    accel(0:1,idx) = accel(0:1,idx)+field_accel(i,k,:)
                    rotat(idx) = rotat(idx)+field_rotat(i,k)

                    mat = Tdomain%specel(n)%mat_index
                    ! Kinetic energy
                    if (allocated(K_energy)) then
                        K_energy(idx) = K_energy(idx) + 0.5_fpp * Tdomain%specel(n)%Density(i,k) * (field_veloc(i,k,0)**2 + field_veloc(i,k,1)**2)
                    endif

                    if ((.not. Tdomain%specel(n)%acoustic) .and. (.not. Tdomain%specel(n)%PML)) then
                        ! Physical derivatives
                        dUx_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * field_displ(:,k,0))
                        dUx_deta = sum(field_displ(i,:,0) * Tdomain%sSubdomain(mat)%hprimez(:,k))
                        dUz_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * field_displ(:,k,1))
                        dUz_deta = sum(field_displ(i,:,1) * Tdomain%sSubdomain(mat)%hprimez(:,k))

                        invgrad_ij = Tdomain%specel(n)%InvGrad(i,k,:,:)
                        DXX = invgrad_ij(0,0)*dUx_dxi + invgrad_ij(0,1)*dUx_deta
                        DXZ = invgrad_ij(1,0)*dUx_dxi + invgrad_ij(1,1)*dUx_deta
                        DZX = invgrad_ij(0,0)*dUz_dxi + invgrad_ij(0,1)*dUz_deta
                        DZZ = invgrad_ij(1,0)*dUz_dxi + invgrad_ij(1,1)*dUz_deta

                        if (allocated(dUxdx)) then
                            dUxdx(idx) = dUxdx(idx) + DXX
                            dUxdz(idx) = dUxdz(idx) + DXZ
                            dUzdx(idx) = dUzdx(idx) + DZX
                            dUzdz(idx) = dUzdz(idx) + DZZ
                        endif

                        eps_xx = DXX
                        eps_zz = DZZ
                        eps_xz = 0.5_fpp * (DXZ + DZX)
                        eps_v = eps_xx + eps_zz

                        if (allocated(eps_vol)) then
                            eps_vol(idx) = eps_vol(idx) + eps_v
                        endif

                        if (allocated(eps_dev_xx)) then
                            eps_dev_xx(idx) = eps_dev_xx(idx) + (eps_xx - 0.5_fpp * eps_v)
                            eps_dev_zz(idx) = eps_dev_zz(idx) + (eps_zz - 0.5_fpp * eps_v)
                            eps_dev_xz(idx) = eps_dev_xz(idx) + eps_xz
                        endif

                        if (allocated(Tdomain%specel(n)%Cij2d)) then
                            strain_v = (/ eps_xx, eps_zz, 2.0_fpp * eps_xz /)
                            stress_v = matmul(Tdomain%specel(n)%Cij2d(:,:,i,k), strain_v)
                            sig_xx = stress_v(1)
                            sig_zz = stress_v(2)
                            sig_xz = stress_v(3)
                        else
                            sig_xx = Tdomain%specel(n)%Lambda(i,k) * eps_v + 2.0_fpp * Tdomain%specel(n)%Mu(i,k) * eps_xx
                            sig_zz = Tdomain%specel(n)%Lambda(i,k) * eps_v + 2.0_fpp * Tdomain%specel(n)%Mu(i,k) * eps_zz
                            sig_xz = 2.0_fpp * Tdomain%specel(n)%Mu(i,k) * eps_xz
                        endif
                        sig_mean = 0.5_fpp * (sig_xx + sig_zz)

                        if (allocated(press)) then
                            press(idx) = press(idx) - sig_mean
                        endif

                        if (allocated(sig_dev_xx)) then
                            sig_dev_xx(idx) = sig_dev_xx(idx) + (sig_xx - sig_mean)
                            sig_dev_zz(idx) = sig_dev_zz(idx) + (sig_zz - sig_mean)
                            sig_dev_xz(idx) = sig_dev_xz(idx) + sig_xz
                        endif

                        if (allocated(P_energy)) then
                            P_energy(idx) = P_energy(idx) + 0.5_fpp * (sig_xx * eps_xx + sig_zz * eps_zz + 2.0_fpp * sig_xz * eps_xz)
                        endif

                        if (allocated(L_energy)) then
                            L_energy(idx) = L_energy(idx) + Tdomain%specel(n)%Mu(i,k)/2.0_fpp * (DXZ - DZX)**2
                            S_energy(idx) = S_energy(idx) + (0.5_fpp * Tdomain%specel(n)%Lambda(i,k) + Tdomain%specel(n)%Mu(i,k)) * eps_v**2
                            R_energy(idx) = R_energy(idx) + 2.0_fpp * Tdomain%specel(n)%Mu(i,k) * DXZ * DZX - 2.0_fpp * Tdomain%specel(n)%Mu(i,k) * eps_xx * eps_zz
                        endif
                    else if (Tdomain%specel(n)%acoustic) then
                        if (allocated(press)) then
                            if (allocated(Tdomain%specel(n)%IDensTensor2d)) then
                                sig_mean = field_veloc(i,k,0) ! VelPhi = -p
                            else
                                dUx_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * field_displ(:,k,0))
                                dUx_deta = sum(field_displ(i,:,0) * Tdomain%sSubdomain(mat)%hprimez(:,k))
                                dUz_dxi = sum(Tdomain%sSubdomain(mat)%hTprimex(:,i) * field_displ(:,k,1))
                                dUz_deta = sum(field_displ(i,:,1) * Tdomain%sSubdomain(mat)%hprimez(:,k))

                                invgrad_ij = Tdomain%specel(n)%InvGrad(i,k,:,:)
                                DXX = invgrad_ij(0,0)*dUx_dxi + invgrad_ij(0,1)*dUx_deta
                                DZZ = invgrad_ij(1,0)*dUz_dxi + invgrad_ij(1,1)*dUz_deta
                                eps_v = DXX + DZZ
                                sig_mean = -Tdomain%specel(n)%Lambda(i,k) * eps_v
                            endif
                            press(idx) = press(idx) - sig_mean
                        endif
                    endif
                end do
            end do
        end do
        ! normalization
        do i = 0,nnodes-1
            if (valence(i)/=0) then
                veloc(:,i) = veloc(:,i)/valence(i)
                accel(:,i) = accel(:,i)/valence(i)
                rotat(i) = rotat(i)/valence(i)
                if (allocated(press)) press(i) = press(i)/valence(i)
                if (allocated(P_energy)) P_energy(i) = P_energy(i)/valence(i)
                if (allocated(K_energy)) K_energy(i) = K_energy(i)/valence(i)
                if (allocated(eps_vol)) eps_vol(i) = eps_vol(i)/valence(i)
                if (allocated(eps_dev_xx)) then
                    eps_dev_xx(i) = eps_dev_xx(i)/valence(i)
                    eps_dev_zz(i) = eps_dev_zz(i)/valence(i)
                    eps_dev_xz(i) = eps_dev_xz(i)/valence(i)
                endif
                if (allocated(sig_dev_xx)) then
                    sig_dev_xx(i) = sig_dev_xx(i)/valence(i)
                    sig_dev_zz(i) = sig_dev_zz(i)/valence(i)
                    sig_dev_xz(i) = sig_dev_xz(i)/valence(i)
                endif
                if (allocated(dUxdx)) then
                    dUxdx(i) = dUxdx(i)/valence(i)
                    dUxdz(i) = dUxdz(i)/valence(i)
                    dUzdx(i) = dUzdx(i)/valence(i)
                    dUzdz(i) = dUzdz(i)/valence(i)
                endif
                if (allocated(L_energy)) then
                    L_energy(i) = L_energy(i)/valence(i)
                    S_energy(i) = S_energy(i)/valence(i)
                    R_energy(i) = R_energy(i)/valence(i)
                endif
            else
                write(*,*) "Elem",i," non traite"
            end if
        end do
        call h5dwrite_f(displ_id, H5T_NATIVE_DOUBLE, displ, dims, hdferr)
        call h5dwrite_f(veloc_id, H5T_NATIVE_DOUBLE, veloc, dims, hdferr)
        call h5dwrite_f(accel_id, H5T_NATIVE_DOUBLE, accel, dims, hdferr)
        call h5dwrite_f(rotat_id, H5T_NATIVE_DOUBLE, rotat, dimr, hdferr)

        if (allocated(press)) call h5dwrite_f(press_id, H5T_NATIVE_DOUBLE, press, dimr, hdferr)
        if (allocated(P_energy)) call h5dwrite_f(P_energy_id, H5T_NATIVE_DOUBLE, P_energy, dimr, hdferr)
        if (allocated(K_energy)) call h5dwrite_f(K_energy_id, H5T_NATIVE_DOUBLE, K_energy, dimr, hdferr)
        if (allocated(eps_vol)) call h5dwrite_f(eps_vol_id, H5T_NATIVE_DOUBLE, eps_vol, dimr, hdferr)
        if (allocated(eps_dev_xx)) then
            call h5dwrite_f(eps_dev_xx_id, H5T_NATIVE_DOUBLE, eps_dev_xx, dimr, hdferr)
            call h5dwrite_f(eps_dev_zz_id, H5T_NATIVE_DOUBLE, eps_dev_zz, dimr, hdferr)
            call h5dwrite_f(eps_dev_xz_id, H5T_NATIVE_DOUBLE, eps_dev_xz, dimr, hdferr)
        endif
        if (allocated(sig_dev_xx)) then
            call h5dwrite_f(sig_dev_xx_id, H5T_NATIVE_DOUBLE, sig_dev_xx, dimr, hdferr)
            call h5dwrite_f(sig_dev_zz_id, H5T_NATIVE_DOUBLE, sig_dev_zz, dimr, hdferr)
            call h5dwrite_f(sig_dev_xz_id, H5T_NATIVE_DOUBLE, sig_dev_xz, dimr, hdferr)
        endif
        if (allocated(dUxdx)) then
            call h5dwrite_f(dUxdx_id, H5T_NATIVE_DOUBLE, dUxdx, dimr, hdferr)
            call h5dwrite_f(dUxdz_id, H5T_NATIVE_DOUBLE, dUxdz, dimr, hdferr)
            call h5dwrite_f(dUzdx_id, H5T_NATIVE_DOUBLE, dUzdx, dimr, hdferr)
            call h5dwrite_f(dUzdz_id, H5T_NATIVE_DOUBLE, dUzdz, dimr, hdferr)
        endif
        if (allocated(L_energy)) then
            call h5dwrite_f(L_energy_id, H5T_NATIVE_DOUBLE, L_energy, dimr, hdferr)
            call h5dwrite_f(S_energy_id, H5T_NATIVE_DOUBLE, S_energy, dimr, hdferr)
            call h5dwrite_f(R_energy_id, H5T_NATIVE_DOUBLE, R_energy, dimr, hdferr)
        endif

        call h5dclose_f(displ_id, hdferr)
        call h5dclose_f(veloc_id, hdferr)
        call h5dclose_f(accel_id, hdferr)
        if (Tdomain%out_var_snap(OUT_PRESSION) == 1) call h5dclose_f(press_id, hdferr)
        if (Tdomain%out_var_snap(OUT_ENERGYP) == 1) call h5dclose_f(P_energy_id, hdferr)
        if (Tdomain%out_var_snap(OUT_ENERGYK) == 1) call h5dclose_f(K_energy_id, hdferr)
        if (Tdomain%out_var_snap(OUT_EPS_VOL) == 1) call h5dclose_f(eps_vol_id, hdferr)
        if (Tdomain%out_var_snap(OUT_EPS_DEV) == 1) then
            call h5dclose_f(eps_dev_xx_id, hdferr)
            call h5dclose_f(eps_dev_zz_id, hdferr)
            call h5dclose_f(eps_dev_xz_id, hdferr)
        endif
        if (Tdomain%out_var_snap(OUT_STRESS_DEV) == 1) then
            call h5dclose_f(sig_dev_xx_id, hdferr)
            call h5dclose_f(sig_dev_zz_id, hdferr)
            call h5dclose_f(sig_dev_xz_id, hdferr)
        endif
        if (Tdomain%out_var_snap(OUT_DUDX) == 1) then
            call h5dclose_f(dUxdx_id, hdferr)
            call h5dclose_f(dUxdz_id, hdferr)
            call h5dclose_f(dUzdx_id, hdferr)
            call h5dclose_f(dUzdz_id, hdferr)
        endif
        if (Tdomain%out_var_snap(OUT_ENERGYD) == 1) then
            call h5dclose_f(L_energy_id, hdferr)
            call h5dclose_f(S_energy_id, hdferr)
            call h5dclose_f(R_energy_id, hdferr)
        endif
        call h5dclose_f(rotat_id, hdferr)
        call h5fclose_f(fid, hdferr)

        if (allocated(press)) deallocate(press)
        if (allocated(P_energy)) deallocate(P_energy)
        if (allocated(K_energy)) deallocate(K_energy)
        if (allocated(eps_vol)) deallocate(eps_vol)
        if (allocated(eps_dev_xx)) deallocate(eps_dev_xx, eps_dev_zz, eps_dev_xz)
        if (allocated(sig_dev_xx)) deallocate(sig_dev_xx, sig_dev_zz, sig_dev_xz)
        if (allocated(dUxdx)) deallocate(dUxdx, dUxdz, dUzdx, dUzdz)
        if (allocated(L_energy)) deallocate(L_energy, S_energy, R_energy)
        deallocate(displ,veloc,valence,rotat)


        call write_xdmf(Tdomain, rg, isort, nnodes)
    end subroutine save_field_h5

    subroutine write_master_xdmf(Tdomain)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer :: n_procs, nelem
        character (len=MAX_FILE_SIZE) :: fnamef
        integer :: rg
        n_procs = Tdomain%mpi_var%n_proc
        nelem = Tdomain%n_elem
        call semname_xdmf_master(fnamef)

        open (61,file=fnamef,status="unknown",form="formatted")
        write(61,"(a)") '<?xml version="1.0" ?>'
        write(61,"(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write(61,"(a)") '<Xdmf Version="2.0" xmlns:xi="http://www.w3.org/2001/XInclude">'
        write(61,"(a)") '<Domain>'
        write(61,"(a)") '<Grid CollectionType="Spatial" GridType="Collection">'
        !!! XXX: recuperer le nom par semname_*
        do rg=0,n_procs-1
            write(61,"(a,I4.4,a)") '<xi:include href="mesh.',rg,'.xmf" xpointer="xpointer(//Xdmf/Domain/Grid)"/>'
        end do
        write(61,"(a)") '</Grid>'
        write(61,"(a)") '</Domain>'
        write(61,"(a)") '</Xdmf>'
        close(61)

    end subroutine write_master_xdmf

    subroutine write_xdmf(Tdomain, rg, isort, nnodes)
        implicit none
        type (domain), intent (IN):: Tdomain
        integer, intent(in) :: rg, isort, nnodes
        !
        character (len=MAX_FILE_SIZE) :: fnamef
        integer   :: i, nn, ne
        real(fpp) :: time
        logical   :: sa, fa
        call semname_xdmf(rg, fnamef)

        nn = nnodes
        ne = Tdomain%n_quad
        sa = has_aniso_solid_2d(Tdomain)   ! this rank has solid-aniso (Cij) output?
        fa = has_aniso_fluid_2d(Tdomain)   ! this rank has fluid-aniso (rho_ij/kappa) output?
        open (61,file=fnamef,status="unknown",form="formatted")
        write(61,"(a)") '<?xml version="1.0" ?>'
        write(61,"(a)") '<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd">'
        write(61,"(a)") '<Xdmf Version="2.0">'
        write(61,"(a)") '<Domain>'
        write(61,"(a)") '<Grid CollectionType="Temporal" GridType="Collection">'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Mat" Format="HDF" Datatype="Int"  Dimensions="',ne, &
            '">geometry',rg,'.h5:/Material</DataItem>'

        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Mass" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',rg,'.h5:/Mass</DataItem>'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="ID" Format="HDF" Datatype="Int"  Dimensions="',ne, &
            '">geometry',rg,'.h5:/ElemID</DataItem>'
        write(61,"(a,I8,a,I4.4,a)") '<DataItem Name="Jac" Format="HDF" Datatype="Int"  Dimensions="',nn, &
            '">geometry',rg,'.h5:/Jac</DataItem>'
        time = 0
        do i=1,isort
            write(61,"(a,I4.4,a,I4.4,a)") '<Grid Name="mesh.',i,'.',rg,'">'
            write(61,"(a,F20.10,a)") '<Time Value="', time,'"/>'
            write(61,"(a,I8,a)") '<Topology Type="Quadrilateral" NumberOfElements="',ne,'">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Int" Dimensions="',ne,' 4">'
            write(61,"(a,I4.4,a)") 'geometry',rg,'.h5:/Elements'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Topology>'
            write(61,"(a)") '<Geometry Type="XYZ">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a)") 'geometry',rg,'.h5:/Nodes'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Geometry>'

            write(61,"(a,I8,a)") '<Attribute Name="Displ" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/displ'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'

            write(61,"(a,I8,a)") '<Attribute Name="Veloc" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/veloc'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'

            write(61,"(a,I8,a)") '<Attribute Name="Accel" Center="Node" AttributeType="Vector" Dimensions="',nn,' 3">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,' 3">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/accel'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'

            write(61,"(a,I4.4,a)") '<Attribute Name="Rotat" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
            write(61,"(a,I4.4,a,I4.4,a)") 'Rsem',i,'/sem_field.',rg,'.h5:/rotat'
            write(61,"(a)") '</DataItem>'
            write(61,"(a)") '</Attribute>'

            if (Tdomain%out_var_snap(OUT_PRESSION) == 1) then
                call write_xdmf_res_scalar(61, "Pressure", i, rg, nn, "pressure")
            endif
            if (Tdomain%out_var_snap(OUT_ENERGYP) == 1) then
                call write_xdmf_res_scalar(61, "P_energy", i, rg, nn, "P_energy")
            endif
            if (Tdomain%out_var_snap(OUT_ENERGYK) == 1) then
                call write_xdmf_res_scalar(61, "K_energy", i, rg, nn, "K_energy")
            endif
            if (Tdomain%out_var_snap(OUT_EPS_VOL) == 1) then
                call write_xdmf_res_scalar(61, "eps_vol", i, rg, nn, "eps_vol")
            endif
            if (Tdomain%out_var_snap(OUT_EPS_DEV) == 1) then
                call write_xdmf_res_scalar(61, "eps_dev_xx", i, rg, nn, "eps_dev_xx")
                call write_xdmf_res_scalar(61, "eps_dev_zz", i, rg, nn, "eps_dev_zz")
                call write_xdmf_res_scalar(61, "eps_dev_xz", i, rg, nn, "eps_dev_xz")
            endif
            if (Tdomain%out_var_snap(OUT_STRESS_DEV) == 1) then
                call write_xdmf_res_scalar(61, "sig_dev_xx", i, rg, nn, "sig_dev_xx")
                call write_xdmf_res_scalar(61, "sig_dev_zz", i, rg, nn, "sig_dev_zz")
                call write_xdmf_res_scalar(61, "sig_dev_xz", i, rg, nn, "sig_dev_xz")
            endif
            if (Tdomain%out_var_snap(OUT_DUDX) == 1) then
                call write_xdmf_res_scalar(61, "dUxdx", i, rg, nn, "dUxdx")
                call write_xdmf_res_scalar(61, "dUxdz", i, rg, nn, "dUxdz")
                call write_xdmf_res_scalar(61, "dUzdx", i, rg, nn, "dUzdx")
                call write_xdmf_res_scalar(61, "dUzdz", i, rg, nn, "dUzdz")
            endif
            if (Tdomain%out_var_snap(OUT_ENERGYD) == 1) then
                call write_xdmf_res_scalar(61, "L_energy", i, rg, nn, "L_energy")
                call write_xdmf_res_scalar(61, "S_energy", i, rg, nn, "S_energy")
                call write_xdmf_res_scalar(61, "R_energy", i, rg, nn, "R_energy")
            endif


            write(61,"(a)") '<Attribute Name="Domain" Center="Grid" AttributeType="Scalar" Dimensions="1">'
            write(61,"(a,I4,a)") '<DataItem Format="XML" Datatype="Int"  Dimensions="1">',rg,'</DataItem>'
            write(61,"(a)") '</Attribute>'

            write(61,"(a,I8,a)") '<Attribute Name="Mat" Center="Cell" AttributeType="Scalar" Dimensions="',ne,'">'
            !write(61,"(a,I8,a,I4.4,a)") '<DataItem Format="HDF" Datatype="Int"  Dimensions="',ne, &
            !    '">geometry',rg,'.h5:/Material</DataItem>'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="Mat"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I8,a)") '<Attribute Name="Mass" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="Mass"]</DataItem>'
!            write(61,"(a,I8,a,I4.4,a)") '<DataItem Format="HDF" Datatype="Int"  Dimensions="',nn, &
!                '">geometry',rg,'.h5:/Mass</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I8,a)") '<Attribute Name="ID" Center="Cell" AttributeType="Scalar" Dimensions="',ne,'">'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="ID"]</DataItem>'
            write(61,"(a)") '</Attribute>'
            write(61,"(a,I8,a)") '<Attribute Name="Jac" Center="Node" AttributeType="Scalar" Dimensions="',nn,'">'
            write(61,"(a,I4,a)") '<DataItem Reference="XML">/Xdmf/Domain/Grid/Grid[',rg+1, &
                ']/DataItem[@Name="Jac"]</DataItem>'
!            write(61,"(a,I8,a,I4.4,a)") '<DataItem Format="HDF" Datatype="Int"  Dimensions="',nn, &
!                '">geometry',rg,'.h5:/Jac</DataItem>'
            write(61,"(a)") '</Attribute>'
            ! Material properties (time-independent, from the geometry h5)
            call write_xdmf_mat_attr(61, "Lambda",  rg, nn)
            call write_xdmf_mat_attr(61, "Mu",      rg, nn)
            call write_xdmf_mat_attr(61, "Density", rg, nn)
            if (sa) then
                call write_xdmf_mat_attr(61, "C11", rg, nn)
                call write_xdmf_mat_attr(61, "C22", rg, nn)
                call write_xdmf_mat_attr(61, "C33", rg, nn)
                call write_xdmf_mat_attr(61, "C12", rg, nn)
                call write_xdmf_mat_attr(61, "C13", rg, nn)
                call write_xdmf_mat_attr(61, "C23", rg, nn)
            end if
            if (fa) then
                call write_xdmf_mat_attr(61, "rho11", rg, nn)
                call write_xdmf_mat_attr(61, "rho22", rg, nn)
                call write_xdmf_mat_attr(61, "rho12", rg, nn)
                call write_xdmf_mat_attr(61, "Kappa", rg, nn)
            end if
            write(61,"(a)") '</Grid>'
            ! XXX inexact pour l'instant
            time = time+Tdomain%TimeD%time_snapshots
        end do
        write(61,"(a)") '</Grid>'
        write(61,"(a)") '</Domain>'
        write(61,"(a)") '</Xdmf>'
        close(61)
    end subroutine write_xdmf

    !! \brief Emit one XDMF node-scalar Attribute referencing geometry<rg>.h5:/<name>.
    subroutine write_xdmf_mat_attr(unit, name, rg, nn)
        implicit none
        integer, intent(in) :: unit, rg, nn
        character(len=*), intent(in) :: name
        write(unit,"(a)") '<Attribute Name="'//trim(name)//'" Center="Node" AttributeType="Scalar">'
        write(unit,"(a,I8,a,I4.4,a)") '<DataItem Format="HDF" NumberType="Float" Precision="8" Dimensions="', &
            nn,'">geometry',rg,'.h5:/'//trim(name)//'</DataItem>'
        write(unit,"(a)") '</Attribute>'
    end subroutine write_xdmf_mat_attr

    subroutine write_constant_fields(Tdomain, fid, irenum, nnodes)
        implicit none
        type (domain), intent (INOUT):: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, dimension(:), intent(in), allocatable :: irenum
        integer, intent(in) :: nnodes
        !
        integer(HID_T) :: mass_id, jac_id
        integer(HSIZE_T), dimension(1) :: dims
        real(fpp), dimension(:),allocatable :: mass, jac
        real(fpp), dimension(:,:), allocatable :: locmass
        integer :: hdferr
        integer :: ngllx, ngllz, idx
        integer :: i, k, n
        
        call create_dset(fid, "Mass", H5T_IEEE_F64LE, nnodes, mass_id)
        call create_dset(fid, "Jac",  H5T_IEEE_F64LE, nnodes, jac_id)

        dims(1) = Tdomain%n_glob_points
        allocate(mass(0:nnodes-1))
        allocate(jac(0:nnodes-1))
        ! mass
        ngllx=-1
        ngllz=-1
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            if (ngllx /= Tdomain%specel(n)%ngllx .or. &
                ngllz /= Tdomain%specel(n)%ngllz) then
                ngllx = Tdomain%specel(n)%ngllx
                ngllz = Tdomain%specel(n)%ngllz
                if (allocated(locmass)) deallocate(locmass)
                allocate(locmass(0:ngllx-1,0:ngllz-1))
            endif

            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            call gather_elem_mass(Tdomain, n, locmass)
            do k = 0,ngllz-1
                do i = 0,ngllx-1
                    idx = irenum(Tdomain%specel(n)%Iglobnum(i,k))
                    mass(idx) = locmass(i,k)
                end do
            end do
        end do
        ! jac
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            do k = 0,ngllz-1
                do i = 0,ngllx-1
                    idx = irenum(Tdomain%specel(n)%Iglobnum(i,k))
                    jac(idx) = Tdomain%specel(n)%Jacob(i,k)
                end do
            end do
        end do

        call h5dwrite_f(mass_id, H5T_NATIVE_DOUBLE, mass, dims, hdferr)
        call h5dclose_f(mass_id, hdferr)
        call h5dwrite_f(jac_id, H5T_NATIVE_DOUBLE, jac, dims, hdferr)
        call h5dclose_f(jac_id, hdferr)
        deallocate(mass,jac,locmass)

    end subroutine write_constant_fields

    !! \brief Write one time-independent per-node scalar field into the geometry h5.
    !! Mirrors the Mass/Jac pattern in write_constant_fields (create + write + close).
    subroutine write_node_scalar(Tdomain, fid, name, nnodes, values)
        implicit none
        type (domain), intent(in) :: Tdomain
        integer(HID_T), intent(in) :: fid
        character(len=*), intent(in) :: name
        integer, intent(in) :: nnodes
        real(fpp), dimension(0:nnodes-1), intent(in) :: values
        integer(HID_T) :: dset_id
        integer(HSIZE_T), dimension(1) :: dims
        integer :: hdferr
        call create_dset(fid, trim(name), H5T_IEEE_F64LE, nnodes, dset_id)
        dims(1) = Tdomain%n_glob_points
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, values, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
    end subroutine write_node_scalar

    !! \brief Write material properties as time-independent per-node scalar fields.
    !! Always: Lambda, Mu, Density (isotropic moduli + rho; for aniso elements these are
    !! the isotropic projection filled in define_arr). Solid aniso: Cij (6 indep comps,
    !! Voigt 2D 1=xx,2=zz,3=xz). Fluid aniso: rho11/rho22/rho12 (inverse-density tensor)
    !! + Kappa (=1/invKappa2d). Mirrors SEM3D snapshots.F90 material output.
    subroutine write_material_fields(Tdomain, fid, irenum, nnodes)
        implicit none
        type (domain), intent(inout) :: Tdomain
        integer(HID_T), intent(in) :: fid
        integer, dimension(:), intent(in), allocatable :: irenum
        integer, intent(in) :: nnodes
        !
        real(fpp), dimension(:), allocatable :: lamb, mu, dens
        real(fpp), dimension(:), allocatable :: c11,c22,c33,c12,c13,c23
        real(fpp), dimension(:), allocatable :: rho11,rho22,rho12,kappa
        logical :: has_solid_aniso, has_fluid_aniso
        integer :: n, i, k, idx, ngllx, ngllz

        allocate(lamb(0:nnodes-1), mu(0:nnodes-1), dens(0:nnodes-1))
        lamb=0._fpp; mu=0._fpp; dens=0._fpp

        has_solid_aniso = has_aniso_solid_2d(Tdomain)
        has_fluid_aniso = has_aniso_fluid_2d(Tdomain)
        if (has_solid_aniso) then
            allocate(c11(0:nnodes-1),c22(0:nnodes-1),c33(0:nnodes-1), &
                     c12(0:nnodes-1),c13(0:nnodes-1),c23(0:nnodes-1))
            c11=0._fpp;c22=0._fpp;c33=0._fpp;c12=0._fpp;c13=0._fpp;c23=0._fpp
        end if
        if (has_fluid_aniso) then
            allocate(rho11(0:nnodes-1),rho22(0:nnodes-1),rho12(0:nnodes-1),kappa(0:nnodes-1))
            rho11=0._fpp;rho22=0._fpp;rho12=0._fpp;kappa=0._fpp
        end if

        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            do k = 0,ngllz-1
                do i = 0,ngllx-1
                    idx = irenum(Tdomain%specel(n)%Iglobnum(i,k))
                    lamb(idx) = Tdomain%specel(n)%Lambda(i,k)
                    mu(idx)   = Tdomain%specel(n)%Mu(i,k)
                    dens(idx) = Tdomain%specel(n)%Density(i,k)
                    if (has_solid_aniso .and. allocated(Tdomain%specel(n)%Cij2d)) then
                        c11(idx)=Tdomain%specel(n)%Cij2d(1,1,i,k)
                        c22(idx)=Tdomain%specel(n)%Cij2d(2,2,i,k)
                        c33(idx)=Tdomain%specel(n)%Cij2d(3,3,i,k)
                        c12(idx)=Tdomain%specel(n)%Cij2d(1,2,i,k)
                        c13(idx)=Tdomain%specel(n)%Cij2d(1,3,i,k)
                        c23(idx)=Tdomain%specel(n)%Cij2d(2,3,i,k)
                    end if
                    if (has_fluid_aniso .and. allocated(Tdomain%specel(n)%IDensTensor2d)) then
                        rho11(idx)=Tdomain%specel(n)%IDensTensor2d(1,1,i,k)
                        rho22(idx)=Tdomain%specel(n)%IDensTensor2d(2,2,i,k)
                        rho12(idx)=Tdomain%specel(n)%IDensTensor2d(1,2,i,k)
                        if (Tdomain%specel(n)%invKappa2d(i,k) /= 0._fpp) &
                            kappa(idx)=1._fpp/Tdomain%specel(n)%invKappa2d(i,k)
                    end if
                end do
            end do
        end do

        call write_node_scalar(Tdomain, fid, "Lambda",  nnodes, lamb)
        call write_node_scalar(Tdomain, fid, "Mu",      nnodes, mu)
        call write_node_scalar(Tdomain, fid, "Density", nnodes, dens)
        deallocate(lamb, mu, dens)
        if (has_solid_aniso) then
            call write_node_scalar(Tdomain, fid, "C11", nnodes, c11)
            call write_node_scalar(Tdomain, fid, "C22", nnodes, c22)
            call write_node_scalar(Tdomain, fid, "C33", nnodes, c33)
            call write_node_scalar(Tdomain, fid, "C12", nnodes, c12)
            call write_node_scalar(Tdomain, fid, "C13", nnodes, c13)
            call write_node_scalar(Tdomain, fid, "C23", nnodes, c23)
            deallocate(c11,c22,c33,c12,c13,c23)
        end if
        if (has_fluid_aniso) then
            call write_node_scalar(Tdomain, fid, "rho11", nnodes, rho11)
            call write_node_scalar(Tdomain, fid, "rho22", nnodes, rho22)
            call write_node_scalar(Tdomain, fid, "rho12", nnodes, rho12)
            call write_node_scalar(Tdomain, fid, "Kappa", nnodes, kappa)
            deallocate(rho11,rho22,rho12,kappa)
        end if
    end subroutine write_material_fields

    !! Does this rank have any solid-aniso (Cij2d) / fluid-aniso (IDensTensor2d) output element?
    logical function has_aniso_solid_2d(Tdomain)
        type (domain), intent(in) :: Tdomain
        integer :: n
        has_aniso_solid_2d = .false.
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            if (allocated(Tdomain%specel(n)%Cij2d)) then
                has_aniso_solid_2d = .true.; return
            end if
        end do
    end function has_aniso_solid_2d

    logical function has_aniso_fluid_2d(Tdomain)
        type (domain), intent(in) :: Tdomain
        integer :: n
        has_aniso_fluid_2d = .false.
        do n = 0,Tdomain%n_elem-1
            if (.not. Tdomain%specel(n)%OUTPUT) cycle
            if (allocated(Tdomain%specel(n)%IDensTensor2d)) then
                has_aniso_fluid_2d = .true.; return
            end if
        end do
    end function has_aniso_fluid_2d


    !! \brief subroutine calculant le rotationel d'un champ de vitesse
    !! pour un element.
    !!
    !! \author Sebastien Terrana
    !! \date 09/03/2015
    !! \param type (Domain), intent (IN) Tdomain
    !! \param integer, intent (IN) nelem
    !<
    subroutine compute_rotational(Tdomain,nel,ngllx,ngllz,field_veloc,field_rotat)
        implicit none
        type (domain), intent (IN) :: Tdomain
        integer, intent (IN) :: nel, ngllx, ngllz
        real(fpp), dimension(0:ngllx-1,0:ngllz-1,0:1), intent (IN):: field_veloc
        real(fpp), dimension(0:ngllx-1,0:ngllz-1), intent (INOUT) :: field_rotat
        integer :: mat

        mat = Tdomain%specel(nel)%mat_index

        field_rotat = Tdomain%specel(nel)%InvGrad(:,:,0,0) * &
                      MATMUL (Tdomain%Ssubdomain(mat)%hTprimex, field_veloc(:,:,1)) &
                    + Tdomain%specel(nel)%InvGrad(:,:,0,1) * &
                      MATMUL (field_veloc(:,:,1), Tdomain%Ssubdomain(mat)%hprimez) &
                    - Tdomain%specel(nel)%InvGrad(:,:,1,0) * &
                      MATMUL (Tdomain%Ssubdomain(mat)%hTprimex, field_veloc(:,:,0)) &
                    - Tdomain%specel(nel)%InvGrad(:,:,1,1) * &
                      MATMUL (field_veloc(:,:,0), Tdomain%Ssubdomain(mat)%hprimez)

        return
    end subroutine compute_rotational

    subroutine write_xdmf_res_scalar(unit, name, isort, rg, nn, dataset_name)
        implicit none
        integer, intent(in) :: unit, isort, rg, nn
        character(len=*), intent(in) :: name, dataset_name
        write(unit,"(a)") '<Attribute Name="'//trim(name)//'" Center="Node" AttributeType="Scalar">'
        write(unit,"(a,I8,a)") '<DataItem Format="HDF" Datatype="Float" Precision="8" Dimensions="',nn,'">'
        write(unit,"(a,I4.4,a,I4.4,a)") 'Rsem',isort,'/sem_field.',rg,'.h5:/'//trim(dataset_name)
        write(unit,"(a)") '</DataItem>'
        write(unit,"(a)") '</Attribute>'
    end subroutine write_xdmf_res_scalar


end module msnapshots

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
