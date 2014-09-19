!>
!! \file save_checkpoint.f90
!! \brief Gère la protection de Sem3d
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine init_protection(Tdomain, it, rg, prot_file)
    use sdomain
    use semdatafiles
    use mpi
    use sem_c_config, only : sem_mkdir
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: it, rg
    character (len=MAX_FILE_SIZE), INTENT(OUT) :: prot_file
    character (len=MAX_FILE_SIZE) :: dir_prot, dir_prot_prev, times_file
    character (len=MAX_FILE_SIZE) :: dir_traces, dir_prot_traces
    character (len=MAX_FILE_SIZE) :: commande
    integer :: ierr

    call semname_protection_iter_rank_file(it,rg,prot_file)


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
#if ! NEW_GLOBAL_METHOD
subroutine compute_save_offsets(Tdomain, offset, offset_f, offset_e, offset_v)
    use sdomain
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(kind=4), intent(out), dimension (12) :: offset
    integer(kind=4), intent(out), dimension (6) :: offset_f
    integer(kind=4), intent(out), dimension (6) :: offset_e
    integer(kind=4), intent(out), dimension (6) :: offset_v
    integer :: n,ngllx,nglly,ngllz,ngll,ngll1,ngll2,n_solid

    ! Calcul des offsets de positions dans les tableaux
    n_solid = Tdomain%n_sls
    offset(1:12)=0
    offset_f(1:6)=0
    offset_e(1:6)=0
    offset_v(1:6)=0
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz
        ngll = (ngllx-2)*(nglly-2)*(ngllz-2)
        ngll2 = ngllx*nglly*ngllz
        if(Tdomain%specel(n)%solid)then  ! solid part
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
            offset(8) = offset(8) + ngll2*18
        end if
        else   ! fluid part
        ! pour VelPhi : 9
        offset(9) = offset(9) + ngll
        if ( .not. Tdomain%specel(n)%PML ) then
            ! VelPhi1, VelPhi2, VelPhi3 : 10
            offset(10) = offset(10) + 0
            ! pour Phi : 11
            offset(11) = offset(11) + ngll
            ! pour Veloc : 12
            offset(12) = offset(12) + 0
        else  ! PML
            ! VelPhi1, VelPhi2, VelPhi3 : 10
            offset(10) = offset(10) + ngll
            ! pour Phi : 11
            offset(11) = offset(11) + 0
            ! pour Veloc : 12
            offset(12) = offset(12) + ngll2*9
        end if

        end if
    end do
    ! Save Fields for Faces
    do n = 0,Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
        ngll = (ngll1-2)*(ngll2-2)
        if(Tdomain%sFace(n)%solid)then
        if (.not. Tdomain%sFace(n)%PML ) then
            offset_f(1) = offset_f(1) + 3*ngll
            offset_f(2) = offset_f(2) + 3*ngll
            offset_f(3) = offset_f(3) + 0
        else
            offset_f(1) = offset_f(1) + 3*ngll
            offset_f(2) = offset_f(2) + 0
            offset_f(3) = offset_f(3) + 3*ngll
        end if
        else   ! fluid case
        if (.not. Tdomain%sFace(n)%PML ) then
            offset_f(4) = offset_f(4) + ngll
            offset_f(5) = offset_f(5) + ngll
            offset_f(6) = offset_f(6) + 0
        else
            offset_f(4) = offset_f(4) + ngll
            offset_f(5) = offset_f(5) + 0
            offset_f(6) = offset_f(6) + ngll
        end if
        end if
    end do
    ! Save Fields for Edges
    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll-2
        if(Tdomain%sEdge(n)%solid)then
        if (.not. Tdomain%sEdge(n)%PML ) then
            offset_e(1) = offset_e(1) + 3*ngll
            offset_e(2) = offset_e(2) + 3*ngll
            offset_e(3) = offset_e(3) + 0
        else
            offset_e(1) = offset_e(1) + 3*ngll
            offset_e(2) = offset_e(2) + 0
            offset_e(3) = offset_e(3) + 3*ngll
        end if
        else
        if (.not. Tdomain%sEdge(n)%PML ) then
            offset_e(4) = offset_e(4) + ngll
            offset_e(5) = offset_e(5) + ngll
            offset_e(6) = offset_e(6) + 0
        else
            offset_e(4) = offset_e(4) + ngll
            offset_e(5) = offset_e(5) + 0
            offset_e(6) = offset_e(6) + ngll
        end if
        end if
    enddo

    ! Save Fields for Vertices
    do n = 0,Tdomain%n_vertex-1
        if(Tdomain%sVertex(n)%solid)then
        if (.not. Tdomain%sVertex(n)%PML ) then
            offset_v(1) = offset_v(1) + 3
            offset_v(2) = offset_v(2) + 3
            offset_v(3) = offset_v(3) + 0
        else
            offset_v(1) = offset_v(1) + 3
            offset_v(2) = offset_v(2) + 0
            offset_v(3) = offset_v(3) + 3
        end if
        else
        if (.not. Tdomain%sVertex(n)%PML ) then
            offset_v(4) = offset_v(4) + 1
            offset_v(5) = offset_v(5) + 1
            offset_v(6) = offset_v(6) + 0
        else
            offset_v(4) = offset_v(4) + 1
            offset_v(5) = offset_v(5) + 0
            offset_v(6) = offset_v(6) + 1
        end if

        end if
    enddo

    !write(*,*) "Offsets (C):", offset
    !write(*,*) "Offsets (F):", offset_f
    !write(*,*) "Offsets (E):", offset_e
    !write(*,*) "Offsets (V):", offset_v
end subroutine compute_save_offsets
#else
subroutine compute_save_offsets_2(Tdomain, offset)
    use sdomain
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(kind=4), intent(out), dimension (12) :: offset
    integer :: n,ngllx,nglly,ngllz,ngll,ngll1,ngll2,n_solid

    ! Calcul des offsets de positions dans les tableaux
    n_solid = Tdomain%n_sls
    offset(1:12)=0
    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz
        ngll = (ngllx-2)*(nglly-2)*(ngllz-2)
        ngll2 = ngllx*nglly*ngllz
        if(Tdomain%specel(n)%solid)then  ! solid part
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
            offset(8) = offset(8) + ngll2*18
        end if
        else   ! fluid part
        ! pour VelPhi : 9
        offset(9) = offset(9) + ngll
        if ( .not. Tdomain%specel(n)%PML ) then
            ! VelPhi1, VelPhi2, VelPhi3 : 10
            offset(10) = offset(10) + 0
            ! pour Phi : 11
            offset(11) = offset(11) + ngll
            ! pour Veloc : 12
            offset(12) = offset(12) + 0
        else  ! PML
            ! VelPhi1, VelPhi2, VelPhi3 : 10
            offset(10) = offset(10) + ngll
            ! pour Phi : 11
            offset(11) = offset(11) + 0
            ! pour Veloc : 12
            offset(12) = offset(12) + ngll2*9
        end if

        end if
    end do
end subroutine compute_save_offsets_2
#endif

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
#if NEW_GLOBAL_METHOD
    integer :: id
#endif


    call create_dset(elem_id, "Veloc", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then  ! all fluid simulation
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(.not. Tdomain%specel(n)%solid) cycle
        if(Tdomain%specel(n)%PML) cycle
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
#if ! NEW_GLOBAL_METHOD
                    data(idx+0) = Tdomain%specel(n)%sl%Veloc(i,j,k,0)
                    data(idx+1) = Tdomain%specel(n)%sl%Veloc(i,j,k,1)
                    data(idx+2) = Tdomain%specel(n)%sl%Veloc(i,j,k,2)
#else
                    id = Tdomain%specel(n)%ISol(i,j,k)
                    data(idx+0) = Tdomain%champs0%Veloc(id,0)
                    data(idx+1) = Tdomain%champs0%Veloc(id,1)
                    data(idx+2) = Tdomain%champs0%Veloc(id,2)
#endif
                    idx = idx + 3
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Veloc

subroutine write_VelPhi(Tdomain, nmax, elem_id)
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


    call create_dset(elem_id, "VelPhi", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then   ! all solid simulation
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%solid) cycle
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
#if ! NEW_GLOBAL_METHOD
                    data(idx) = Tdomain%specel(n)%fl%VelPhi(i,j,k)
#else
                    data(idx) = Tdomain%champs0%VelPhi(Tdomain%specel(n)%IFlu(i,j,k))
#endif
                    idx = idx + 1
                enddo
            enddo
        enddo
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_VelPhi



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
#if NEW_GLOBAL_METHOD
    integer :: id
#endif
    real(kind=8), dimension(1:nmax) :: data1, data2, data3
    integer(HSIZE_T), dimension(1) :: dims



    call create_dset(elem_id, "Veloc1", H5T_IEEE_F64LE, nmax, dset_id1)
    call create_dset(elem_id, "Veloc2", H5T_IEEE_F64LE, nmax, dset_id2)
    call create_dset(elem_id, "Veloc3", H5T_IEEE_F64LE, nmax, dset_id3)
    if(nmax <= 0)then
        call h5dclose_f(dset_id1, hdferr)
        call h5dclose_f(dset_id2, hdferr)
        call h5dclose_f(dset_id3, hdferr)
        return
    end if

    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(.not. Tdomain%specel(n)%solid) cycle
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
#if ! NEW_GLOBAL_METHOD
                        data1(idx+0) = Tdomain%specel(n)%slpml%Veloc1(i,j,k,0)
                        data1(idx+1) = Tdomain%specel(n)%slpml%Veloc1(i,j,k,1)
                        data1(idx+2) = Tdomain%specel(n)%slpml%Veloc1(i,j,k,2)
                        data2(idx+0) = Tdomain%specel(n)%slpml%Veloc2(i,j,k,0)
                        data2(idx+1) = Tdomain%specel(n)%slpml%Veloc2(i,j,k,1)
                        data2(idx+2) = Tdomain%specel(n)%slpml%Veloc2(i,j,k,2)
                        data3(idx+0) = Tdomain%specel(n)%slpml%Veloc3(i,j,k,0)
                        data3(idx+1) = Tdomain%specel(n)%slpml%Veloc3(i,j,k,1)
                        data3(idx+2) = Tdomain%specel(n)%slpml%Veloc3(i,j,k,2)
#else
                        id = Tdomain%specel(n)%slpml%ISolPML(i,j,k)
                        data1(idx+0) = Tdomain%champs0%VelocPML(id,0)
                        data1(idx+1) = Tdomain%champs0%VelocPML(id,1)
                        data1(idx+2) = Tdomain%champs0%VelocPML(id,2)
                        data2(idx+0) = Tdomain%champs0%VelocPML(id+1,0)
                        data2(idx+1) = Tdomain%champs0%VelocPML(id+1,1)
                        data2(idx+2) = Tdomain%champs0%VelocPML(id+1,2)
                        data3(idx+0) = Tdomain%champs0%VelocPML(id+2,0)
                        data3(idx+1) = Tdomain%champs0%VelocPML(id+2,1)
                        data3(idx+2) = Tdomain%champs0%VelocPML(id+2,2)
#endif
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

subroutine write_VelPhi123(Tdomain, nmax, elem_id)
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

    call create_dset(elem_id, "VelPhi1", H5T_IEEE_F64LE, nmax, dset_id1)
    call create_dset(elem_id, "VelPhi2", H5T_IEEE_F64LE, nmax, dset_id2)
    call create_dset(elem_id, "VelPhi3", H5T_IEEE_F64LE, nmax, dset_id3)
    if(nmax <= 0)then
        call h5dclose_f(dset_id1, hdferr)
        call h5dclose_f(dset_id2, hdferr)
        call h5dclose_f(dset_id3, hdferr)
        return
    end if
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%solid) cycle
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
                        data1(idx) = Tdomain%specel(n)%flpml%VelPhi1(i,j,k)
                        data2(idx) = Tdomain%specel(n)%flpml%VelPhi2(i,j,k)
                        data3(idx) = Tdomain%specel(n)%flpml%VelPhi3(i,j,k)
                        idx = idx + 1
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
end subroutine write_VelPhi123



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
#if NEW_GLOBAL_METHOD
    integer :: id
#endif

    call create_dset(elem_id, "Displ", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(.not. Tdomain%specel(n)%solid) cycle
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
#if ! NEW_GLOBAL_METHOD
                        data(idx+0) = Tdomain%specel(n)%sl%Displ(i,j,k,0)
                        data(idx+1) = Tdomain%specel(n)%sl%Displ(i,j,k,1)
                        data(idx+2) = Tdomain%specel(n)%sl%Displ(i,j,k,2)
#else
                        id = Tdomain%specel(n)%ISol(i,j,k)
                        data(idx+0) = Tdomain%champs0%Depla(id,0)
                        data(idx+1) = Tdomain%champs0%Depla(id,1)
                        data(idx+2) = Tdomain%champs0%Depla(id,2)
#endif
                        idx = idx + 3
                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Disp


subroutine write_Phi(Tdomain, nmax, elem_id)
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

    call create_dset(elem_id, "Phi", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%solid) cycle
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
#if ! NEW_GLOBAL_METHOD
                        data(idx) = Tdomain%specel(n)%fl%Phi(i,j,k)
#else
                        data(idx) = Tdomain%champs0%Phi(Tdomain%specel(n)%IFlu(i,j,k))
#endif
                        idx = idx + 1
                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Phi


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
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if
    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(.not. Tdomain%specel(n)%solid) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if ( .not. Tdomain%specel(n)%PML ) then
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
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if


    dims(1) = nmax
    idx = 1
    n_solid = Tdomain%n_sls

    do n = 0,Tdomain%n_elem-1
        if(.not. Tdomain%specel(n)%solid) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if ( .not. Tdomain%specel(n)%PML ) then
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
        if(.not. Tdomain%specel(n)%solid) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if ( .not. Tdomain%specel(n)%PML ) then
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
        if(.not. Tdomain%specel(n)%solid) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if ( .not. Tdomain%specel(n)%PML ) then
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
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if

    dims(1) = nmax
    idx = 1
    do n = 0,Tdomain%n_elem-1
        if(.not. Tdomain%specel(n)%solid) cycle
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
        end if
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Stress

subroutine write_Veloc_Fluid_PML(Tdomain, nmax, elem_id)
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
#if NEW_GLOBAL_METHOD
    integer :: id
#endif

    call create_dset(elem_id, "Veloc_Fl", H5T_IEEE_F64LE, nmax, dset_id)
    if(nmax <= 0)then
        call h5dclose_f(dset_id, hdferr)
        return
    end if

    dims(1) = nmax
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
                        if (idx.gt.nmax) then
                            write(*,*) "Erreur fatale sauvegarde des protections"
                            stop 1
                        end if
#if ! NEW_GLOBAL_METHOD
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
#else
                        id = Tdomain%specel(n)%slpml%ISolPML(i,j,k)
                        data(idx+ 0) = Tdomain%champs0%VelocPML(id,0)
                        data(idx+ 1) = Tdomain%champs0%VelocPML(id,1)
                        data(idx+ 2) = Tdomain%champs0%VelocPML(id,2)
                        idx = idx + 3
                        data(idx+ 0) = Tdomain%champs0%VelocPML(id+1,0)
                        data(idx+ 1) = Tdomain%champs0%VelocPML(id+1,1)
                        data(idx+ 2) = Tdomain%champs0%VelocPML(id+1,2)
                        idx = idx + 3
                        data(idx+ 0) = Tdomain%champs0%VelocPML(id+2,0)
                        data(idx+ 1) = Tdomain%champs0%VelocPML(id+2,1)
                        data(idx+ 2) = Tdomain%champs0%VelocPML(id+2,2)
                        idx = idx + 3
#endif

                    enddo
                enddo
            enddo
        end if
    enddo
    call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, hdferr)
    call h5dclose_f(dset_id, hdferr)
end subroutine write_Veloc_Fluid_PML


#if ! NEW_GLOBAL_METHOD
subroutine write_Faces(Tdomain, offset_f, face_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: face_id
    integer(kind=4), dimension(6), intent(IN) :: offset_f

    integer(HID_T) :: veloc_id, displ_id, veloc1_id, veloc2_id, veloc3_id,  &
                      velphi_id, phi_id, velphi1_id, velphi2_id, velphi3_id
    integer :: n,ngll1,ngll2,idx1,idx2,idx3,idx4,idx5,idx6,i,j,hdferr
    real(kind=8), dimension(1:offset_f(1)) :: veloc
    real(kind=8), dimension(1:offset_f(2)) :: displ
    real(kind=8), dimension(1:offset_f(3)) :: veloc1
    real(kind=8), dimension(1:offset_f(3)) :: veloc2
    real(kind=8), dimension(1:offset_f(3)) :: veloc3
    real(kind=8), dimension(1:offset_f(4)) :: velphi
    real(kind=8), dimension(1:offset_f(5)) :: phi
    real(kind=8), dimension(1:offset_f(6)) :: velphi1
    real(kind=8), dimension(1:offset_f(6)) :: velphi2
    real(kind=8), dimension(1:offset_f(6)) :: velphi3

    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(face_id, "Veloc", H5T_IEEE_F64LE, offset_f(1), veloc_id)
    call create_dset(face_id, "Displ", H5T_IEEE_F64LE, offset_f(2), displ_id)
    call create_dset(face_id, "Veloc1", H5T_IEEE_F64LE, offset_f(3), veloc1_id)
    call create_dset(face_id, "Veloc2", H5T_IEEE_F64LE, offset_f(3), veloc2_id)
    call create_dset(face_id, "Veloc3", H5T_IEEE_F64LE, offset_f(3), veloc3_id)
    call create_dset(face_id, "VelPhi", H5T_IEEE_F64LE, offset_f(4), velphi_id)
    call create_dset(face_id, "Phi", H5T_IEEE_F64LE, offset_f(5), phi_id)
    call create_dset(face_id, "VelPhi1", H5T_IEEE_F64LE, offset_f(6), velphi1_id)
    call create_dset(face_id, "VelPhi2", H5T_IEEE_F64LE, offset_f(6), velphi2_id)
    call create_dset(face_id, "VelPhi3", H5T_IEEE_F64LE, offset_f(6), velphi3_id)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    idx4 = 1
    idx5 = 1
    idx6 = 1
    do n = 0,Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1
        ngll2 = Tdomain%sFace(n)%ngll2
        if(Tdomain%sFace(n)%solid)then   ! solid face
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
                    veloc1(idx3+0) = Tdomain%sFace(n)%spml%Veloc1(i,j,0)
                    veloc1(idx3+1) = Tdomain%sFace(n)%spml%Veloc1(i,j,1)
                    veloc1(idx3+2) = Tdomain%sFace(n)%spml%Veloc1(i,j,2)
                    veloc2(idx3+0) = Tdomain%sFace(n)%spml%Veloc2(i,j,0)
                    veloc2(idx3+1) = Tdomain%sFace(n)%spml%Veloc2(i,j,1)
                    veloc2(idx3+2) = Tdomain%sFace(n)%spml%Veloc2(i,j,2)
                    veloc3(idx3+0) = Tdomain%sFace(n)%spml%Veloc3(i,j,0)
                    veloc3(idx3+1) = Tdomain%sFace(n)%spml%Veloc3(i,j,1)
                    veloc3(idx3+2) = Tdomain%sFace(n)%spml%Veloc3(i,j,2)
                    idx3 = idx3 + 3
                end if
            enddo
        enddo
        else    ! fluid face
        do j = 1,ngll2-2
            do i = 1,ngll1-2
                if (idx4.gt.offset_f(4)) then
                    write(*,*) "Erreur fatale sauvegarde des protections"
                    stop 1
                end if
                velphi(idx4) = Tdomain%sFace(n)%VelPhi(i,j)
                idx4 = idx4 + 1
                if (.not. Tdomain%sFace(n)%PML ) then
                    phi(idx5) = Tdomain%sFace(n)%Phi(i,j)
                    idx5 = idx5 + 1
                else
                    velphi1(idx6) = Tdomain%sFace(n)%spml%VelPhi1(i,j)
                    velphi2(idx6) = Tdomain%sFace(n)%spml%VelPhi2(i,j)
                    velphi3(idx6) = Tdomain%sFace(n)%spml%VelPhi3(i,j)
                    idx6 = idx6 + 1
                end if
            enddo
        enddo

        end if
    enddo
    dims(1) = offset_f(1)
    if(dims(1) > 0) call h5dwrite_f(veloc_id, H5T_NATIVE_DOUBLE, veloc, dims, hdferr)
    dims(1) = offset_f(2)
    if(dims(1) > 0) call h5dwrite_f(displ_id, H5T_NATIVE_DOUBLE, displ, dims, hdferr)
    dims(1) = offset_f(3)
    if(dims(1) > 0)then
        call h5dwrite_f(veloc1_id, H5T_NATIVE_DOUBLE, veloc1, dims, hdferr)
        call h5dwrite_f(veloc2_id, H5T_NATIVE_DOUBLE, veloc2, dims, hdferr)
        call h5dwrite_f(veloc3_id, H5T_NATIVE_DOUBLE, veloc3, dims, hdferr)
    end if
    dims(1) = offset_f(4)
    if(dims(1) > 0) call h5dwrite_f(velphi_id, H5T_NATIVE_DOUBLE, velphi, dims, hdferr)
    dims(1) = offset_f(5)
    if(dims(1) > 0) call h5dwrite_f(phi_id, H5T_NATIVE_DOUBLE, phi, dims, hdferr)
    dims(1) = offset_f(6)
    if(dims(1) > 0)then
        call h5dwrite_f(velphi1_id, H5T_NATIVE_DOUBLE, velphi1, dims, hdferr)
        call h5dwrite_f(velphi2_id, H5T_NATIVE_DOUBLE, velphi2, dims, hdferr)
        call h5dwrite_f(velphi3_id, H5T_NATIVE_DOUBLE, velphi3, dims, hdferr)
    end if


    call h5dclose_f(veloc_id, hdferr)
    call h5dclose_f(displ_id, hdferr)
    call h5dclose_f(veloc1_id, hdferr)
    call h5dclose_f(veloc2_id, hdferr)
    call h5dclose_f(veloc3_id, hdferr)
    call h5dclose_f(velphi_id, hdferr)
    call h5dclose_f(phi_id, hdferr)
    call h5dclose_f(velphi1_id, hdferr)
    call h5dclose_f(velphi2_id, hdferr)
    call h5dclose_f(velphi3_id, hdferr)

end subroutine write_Faces

subroutine write_Edges(Tdomain, offset_e, edge_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: edge_id
    integer(kind=4), dimension(6), intent(IN) :: offset_e

    integer(HID_T) :: veloc_id, displ_id, veloc1_id, veloc2_id, veloc3_id,   &
                      velphi_id, phi_id, velphi1_id, velphi2_id, velphi3_id
    integer :: n,ngll,idx1,idx2,idx3,idx4,idx5,idx6,i,j,hdferr
    real(kind=8), dimension(1:offset_e(1)) :: veloc
    real(kind=8), dimension(1:offset_e(2)) :: displ
    real(kind=8), dimension(1:offset_e(3)) :: veloc1
    real(kind=8), dimension(1:offset_e(3)) :: veloc2
    real(kind=8), dimension(1:offset_e(3)) :: veloc3
    real(kind=8), dimension(1:offset_e(4)) :: velphi
    real(kind=8), dimension(1:offset_e(5)) :: phi
    real(kind=8), dimension(1:offset_e(6)) :: velphi1
    real(kind=8), dimension(1:offset_e(6)) :: velphi2
    real(kind=8), dimension(1:offset_e(6)) :: velphi3
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(edge_id, "Veloc", H5T_IEEE_F64LE, offset_e(1), veloc_id)
    call create_dset(edge_id, "Displ", H5T_IEEE_F64LE, offset_e(2), displ_id)
    call create_dset(edge_id, "Veloc1", H5T_IEEE_F64LE, offset_e(3), veloc1_id)
    call create_dset(edge_id, "Veloc2", H5T_IEEE_F64LE, offset_e(3), veloc2_id)
    call create_dset(edge_id, "Veloc3", H5T_IEEE_F64LE, offset_e(3), veloc3_id)
    call create_dset(edge_id, "VelPhi", H5T_IEEE_F64LE, offset_e(4), velphi_id)
    call create_dset(edge_id, "Phi", H5T_IEEE_F64LE, offset_e(5), phi_id)
    call create_dset(edge_id, "VelPhi1", H5T_IEEE_F64LE, offset_e(6), velphi1_id)
    call create_dset(edge_id, "VelPhi2", H5T_IEEE_F64LE, offset_e(6), velphi2_id)
    call create_dset(edge_id, "VelPhi3", H5T_IEEE_F64LE, offset_e(6), velphi3_id)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    idx4 = 1
    idx5 = 1
    idx6 = 1
    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll
        if(Tdomain%sEdge(n)%solid)then
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
                veloc1(idx3+0) = Tdomain%sEdge(n)%spml%Veloc1(i,0)
                veloc1(idx3+1) = Tdomain%sEdge(n)%spml%Veloc1(i,1)
                veloc1(idx3+2) = Tdomain%sEdge(n)%spml%Veloc1(i,2)
                veloc2(idx3+0) = Tdomain%sEdge(n)%spml%Veloc2(i,0)
                veloc2(idx3+1) = Tdomain%sEdge(n)%spml%Veloc2(i,1)
                veloc2(idx3+2) = Tdomain%sEdge(n)%spml%Veloc2(i,2)
                veloc3(idx3+0) = Tdomain%sEdge(n)%spml%Veloc3(i,0)
                veloc3(idx3+1) = Tdomain%sEdge(n)%spml%Veloc3(i,1)
                veloc3(idx3+2) = Tdomain%sEdge(n)%spml%Veloc3(i,2)
                idx3 = idx3 + 3
            end if
        enddo
        else    ! fluid
        do i = 1,ngll-2
            if (idx4.gt.offset_e(4)) then
                write(*,*) "Erreur fatale sauvegarde des protections"
                stop 1
            end if
            if (.not. Tdomain%sEdge(n)%PML ) then
                velPhi(idx4) = Tdomain%sEdge(n)%VelPhi(i)
                idx4 = idx4 + 1
                phi(idx5) = Tdomain%sEdge(n)%Phi(i)
                idx5 = idx5 + 1
            else
                velphi(idx4) = Tdomain%sEdge(n)%VelPhi(i)
                idx4 = idx4 + 1
                velphi1(idx6) = Tdomain%sEdge(n)%spml%VelPhi1(i)
                velphi2(idx6) = Tdomain%sEdge(n)%spml%VelPhi2(i)
                velphi3(idx6) = Tdomain%sEdge(n)%spml%VelPhi3(i)
                idx6 = idx6 + 1
            end if
        enddo

        end if
    enddo
    dims(1) = offset_e(1)
    if(dims(1) > 0) call h5dwrite_f(veloc_id, H5T_NATIVE_DOUBLE, veloc, dims, hdferr)
    dims(1) = offset_e(2)
    if(dims(1) > 0) call h5dwrite_f(displ_id, H5T_NATIVE_DOUBLE, displ, dims, hdferr)
    dims(1) = offset_e(3)
    if(dims(1) > 0)then
        call h5dwrite_f(veloc1_id, H5T_NATIVE_DOUBLE, veloc1, dims, hdferr)
        call h5dwrite_f(veloc2_id, H5T_NATIVE_DOUBLE, veloc2, dims, hdferr)
        call h5dwrite_f(veloc3_id, H5T_NATIVE_DOUBLE, veloc3, dims, hdferr)
    end if
    dims(1) = offset_e(4)
    if(dims(1) > 0) call h5dwrite_f(velphi_id, H5T_NATIVE_DOUBLE, velphi, dims, hdferr)
    dims(1) = offset_e(5)
    if(dims(1) > 0) call h5dwrite_f(phi_id, H5T_NATIVE_DOUBLE, phi, dims, hdferr)
    dims(1) = offset_e(6)
    if(dims(1) > 0)then
        call h5dwrite_f(velphi1_id, H5T_NATIVE_DOUBLE, velphi1, dims, hdferr)
        call h5dwrite_f(velphi2_id, H5T_NATIVE_DOUBLE, velphi2, dims, hdferr)
        call h5dwrite_f(velphi3_id, H5T_NATIVE_DOUBLE, velphi3, dims, hdferr)
    end if

    call h5dclose_f(veloc_id, hdferr)
    call h5dclose_f(displ_id, hdferr)
    call h5dclose_f(veloc1_id, hdferr)
    call h5dclose_f(veloc2_id, hdferr)
    call h5dclose_f(veloc3_id, hdferr)
    call h5dclose_f(velphi_id, hdferr)
    call h5dclose_f(phi_id, hdferr)
    call h5dclose_f(velphi1_id, hdferr)
    call h5dclose_f(velphi2_id, hdferr)
    call h5dclose_f(velphi3_id, hdferr)

end subroutine write_Edges

subroutine write_Vertices(Tdomain, offset_v, vertex_id)
    use sdomain, only : domain
    use sem_hdf5, only : create_dset
    use HDF5
    implicit none
    type (domain), intent (IN):: Tdomain
    integer(HID_T), intent(IN) :: vertex_id
    integer(kind=4), dimension(6), intent(IN) :: offset_v

    integer(HID_T) :: veloc_id, displ_id, veloc1_id, veloc2_id, veloc3_id,  &
                      velphi_id, phi_id, velphi1_id, velphi2_id, velphi3_id
    integer :: n,idx1,idx2,idx3,idx4,idx5,idx6,i,j,hdferr
    real(kind=8), dimension(1:offset_v(1)) :: veloc
    real(kind=8), dimension(1:offset_v(2)) :: displ
    real(kind=8), dimension(1:offset_v(3)) :: veloc1
    real(kind=8), dimension(1:offset_v(3)) :: veloc2
    real(kind=8), dimension(1:offset_v(3)) :: veloc3
    real(kind=8), dimension(1:offset_v(4)) :: velphi
    real(kind=8), dimension(1:offset_v(5)) :: phi
    real(kind=8), dimension(1:offset_v(6)) :: velphi1
    real(kind=8), dimension(1:offset_v(6)) :: velphi2
    real(kind=8), dimension(1:offset_v(6)) :: velphi3
    integer(HSIZE_T), dimension(1) :: dims

    call create_dset(vertex_id, "Veloc", H5T_IEEE_F64LE, offset_v(1), veloc_id)
    call create_dset(vertex_id, "Displ", H5T_IEEE_F64LE, offset_v(2), displ_id)
    call create_dset(vertex_id, "Veloc1", H5T_IEEE_F64LE, offset_v(3), veloc1_id)
    call create_dset(vertex_id, "Veloc2", H5T_IEEE_F64LE, offset_v(3), veloc2_id)
    call create_dset(vertex_id, "Veloc3", H5T_IEEE_F64LE, offset_v(3), veloc3_id)
    call create_dset(vertex_id, "VelPhi", H5T_IEEE_F64LE, offset_v(4), velphi_id)
    call create_dset(vertex_id, "Phi", H5T_IEEE_F64LE, offset_v(5), phi_id)
    call create_dset(vertex_id, "VelPhi1", H5T_IEEE_F64LE, offset_v(6), velphi1_id)
    call create_dset(vertex_id, "VelPhi2", H5T_IEEE_F64LE, offset_v(6), velphi2_id)
    call create_dset(vertex_id, "VelPhi3", H5T_IEEE_F64LE, offset_v(6), velphi3_id)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    idx4 = 1
    idx5 = 1
    idx6 = 1
    do n = 0,Tdomain%n_vertex-1
        if(Tdomain%sVertex(n)%solid)then
        if (idx1 > offset_v(1)) then
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
            veloc1(idx3+0) = Tdomain%sVertex(n)%spml%Veloc1(0)
            veloc1(idx3+1) = Tdomain%sVertex(n)%spml%Veloc1(1)
            veloc1(idx3+2) = Tdomain%sVertex(n)%spml%Veloc1(2)
            veloc2(idx3+0) = Tdomain%sVertex(n)%spml%Veloc2(0)
            veloc2(idx3+1) = Tdomain%sVertex(n)%spml%Veloc2(1)
            veloc2(idx3+2) = Tdomain%sVertex(n)%spml%Veloc2(2)
            veloc3(idx3+0) = Tdomain%sVertex(n)%spml%Veloc3(0)
            veloc3(idx3+1) = Tdomain%sVertex(n)%spml%Veloc3(1)
            veloc3(idx3+2) = Tdomain%sVertex(n)%spml%Veloc3(2)
            idx3 = idx3 + 3
        end if
        else   ! fluid
        if (idx4 > offset_v(4)) then
            write(*,*) "Erreur fatale sauvegarde des protections"
            stop 1
        end if
        if(.not. Tdomain%sVertex(n)%PML)then
            velphi(idx4) = Tdomain%sVertex(n)%VelPhi
            idx4 = idx4 + 1
            phi(idx5) = Tdomain%sVertex(n)%Phi
            idx5 = idx5 + 1
        else
            velphi(idx4) = Tdomain%sVertex(n)%VelPhi
            idx4 = idx4 + 1
            velphi1(idx6) = Tdomain%sVertex(n)%spml%VelPhi1
            velphi2(idx6) = Tdomain%sVertex(n)%spml%Velphi2
            velphi3(idx6) = Tdomain%sVertex(n)%spml%VelPhi3
            idx6 = idx6 + 1
        end if

        end if
    enddo
    dims(1) = offset_v(1)
    if(dims(1) > 0) call h5dwrite_f(veloc_id, H5T_NATIVE_DOUBLE, veloc, dims, hdferr)
    dims(1) = offset_v(2)
    if(dims(1) > 0) call h5dwrite_f(displ_id, H5T_NATIVE_DOUBLE, displ, dims, hdferr)
    dims(1) = offset_v(3)
    if(dims(1) > 0)then
        call h5dwrite_f(veloc1_id, H5T_NATIVE_DOUBLE, veloc1, dims, hdferr)
        call h5dwrite_f(veloc2_id, H5T_NATIVE_DOUBLE, veloc2, dims, hdferr)
        call h5dwrite_f(veloc3_id, H5T_NATIVE_DOUBLE, veloc3, dims, hdferr)
    end if
    dims(1) = offset_v(4)
    if(dims(1) > 0) call h5dwrite_f(velphi_id, H5T_NATIVE_DOUBLE, velphi, dims, hdferr)
    dims(1) = offset_v(5)
    if(dims(1) > 0) call h5dwrite_f(phi_id, H5T_NATIVE_DOUBLE, phi, dims, hdferr)
    dims(1) = offset_v(6)
    if(dims(1) > 0)then
        call h5dwrite_f(velphi1_id, H5T_NATIVE_DOUBLE, velphi1, dims, hdferr)
        call h5dwrite_f(velphi2_id, H5T_NATIVE_DOUBLE, velphi2, dims, hdferr)
        call h5dwrite_f(velphi3_id, H5T_NATIVE_DOUBLE, velphi3, dims, hdferr)
    end if

    call h5dclose_f(veloc_id, hdferr)
    call h5dclose_f(displ_id, hdferr)
    call h5dclose_f(veloc1_id, hdferr)
    call h5dclose_f(veloc2_id, hdferr)
    call h5dclose_f(veloc3_id, hdferr)
    call h5dclose_f(velphi_id, hdferr)
    call h5dclose_f(phi_id, hdferr)
    call h5dclose_f(velphi1_id, hdferr)
    call h5dclose_f(velphi2_id, hdferr)
    call h5dclose_f(velphi3_id, hdferr)

end subroutine write_Vertices
#endif
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
    integer(HID_T) :: fid, dset_id, elem_id
#if ! NEW_GLOBAL_METHOD
    integer(HID_T) :: face_id, edge_id, vertex_id
#endif
    logical :: avail
    integer :: nelem, noffset
    integer :: size_vec, size_eps, size_epsaniso
    integer(kind=4), dimension (12) :: offset
#if ! NEW_GLOBAL_METHOD
    integer(kind=4), dimension (6) :: offset_f ! Veloc / Displ / (Veloc1,Veloc2,Veloc3)
    integer(kind=4), dimension (6) :: offset_e ! Veloc / Displ / (Veloc1,Veloc2,Veloc3)
    integer(kind=4), dimension (6) :: offset_v ! Veloc / Displ / (Veloc1,Veloc2,Veloc3)
#endif
    integer(HSIZE_T), dimension(2) :: off_dims

    if (rg == 0) then
        write (*,'(A44,I8,A1,f10.6)') "--> SEM : protection at iteration and time :",it," ",rtime
    endif
    call init_protection(Tdomain, it, rg, fnamef)

    nelem = Tdomain%n_elem

    off_dims(1) = nelem+1
    off_dims(2) = 12

    call init_hdf5()
#if ! NEW_GLOBAL_METHOD
    call compute_save_offsets(Tdomain, offset, offset_f, offset_e, offset_v)
#else
    call compute_save_offsets_2(Tdomain, offset)
#endif

    call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)

    if (hdferr/=0) then
        write(*,*) "Detected HDF error at open:", fnamef, hdferr
        stop "Error writing HDF file"
    end if

    !write(*,*) "DBG: Create groups"
    ! ifort doesn't care, but gfortran complains that the last integer should be 8 bytes
    call h5gcreate_f(fid, 'Elements', elem_id, hdferr, 0_SIZE_T)
#if ! NEW_GLOBAL_METHOD
    call h5gcreate_f(fid, 'Faces', face_id, hdferr, 0_SIZE_T)
    call h5gcreate_f(fid, 'Edges', edge_id, hdferr, 0_SIZE_T)
    call h5gcreate_f(fid, 'Vertices', vertex_id, hdferr, 0_SIZE_T)
#endif

    !  call create_dset_2d(fid, "Offsets", H5T_STD_I32LE, nelem+1, 8, dset_id)
    !  call h5dwrite_f(dset_id, H5T_STD_I32LE, offset, off_dims, hdferr)
    !  call h5dclose_f(dset_id, hdferr)

    !write(*,*) "DBG: Create attrs"
    call write_attr_real(fid, "rtime", rtime)
    call write_attr_real(fid, "dtmin", dtmin)
    call write_attr_int(fid, "iteration", it)
    call write_attr_int(fid, "isort", isort)
    !write(*,*) "DBG: Write veloc"
    call write_Veloc(Tdomain, offset(1), elem_id)
    !write(*,*) "DBG: Write veloc PML"
    call write_Veloc123(Tdomain, offset(2), elem_id)
    !write(*,*) "DBG: Write disp"
    call write_Disp(Tdomain, offset(3), elem_id)
    !write(*,*) "DBG: Write espvol"
    call write_EpsilonVol(Tdomain, offset(4), elem_id)
    !write(*,*) "DBG: Write rvol"
    call write_Rvol(Tdomain, offset(5), elem_id)
    call write_Rxyz(Tdomain, offset(6), elem_id)
    call write_EpsilonDev(Tdomain, offset(7), elem_id)
    call write_Stress(Tdomain, offset(8), elem_id)
    !write(*,*) "DBG: Write velphi"
    call write_VelPhi(Tdomain, offset(9), elem_id)
    !write(*,*) "DBG: Write velphi123"
    call write_VelPhi123(Tdomain, offset(10), elem_id)
    !write(*,*) "DBG: Write phi"
    call write_Phi(Tdomain, offset(11), elem_id)
    !write(*,*) "DBG: Write fluid pml"
    call write_Veloc_Fluid_PML(Tdomain, offset(12), elem_id)
#if ! NEW_GLOBAL_METHOD
    !write(*,*) "DBG: Write faces"
    call write_Faces(Tdomain, offset_f, face_id)
    !write(*,*) "DBG: Write edges"
    call write_Edges(Tdomain, offset_e, edge_id)
    !write(*,*) "DBG: Write vertices"
    call write_Vertices(Tdomain, offset_v, vertex_id)
#endif

    call h5gclose_f(elem_id, hdferr)
#if ! NEW_GLOBAL_METHOD
    call h5gclose_f(face_id, hdferr)
    call h5gclose_f(edge_id, hdferr)
    call h5gclose_f(vertex_id, hdferr)
#endif
    call h5fclose_f(fid, hdferr)
end subroutine save_checkpoint
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
