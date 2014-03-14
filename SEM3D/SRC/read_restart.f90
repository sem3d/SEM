!>
!!\file read_restart.f90
!!\brief Contient la subroutine read_restart().
!!
!! Gère la reprise de Sem3d



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
        if(.not. Tdomain%specel(n)%solid) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 1,ngllz-2
            do j = 1,nglly-2
                do i = 1,ngllx-2
                    Tdomain%specel(n)%sl%Veloc(i,j,k,0) = veloc(idx1+0)
                    Tdomain%specel(n)%sl%Veloc(i,j,k,1) = veloc(idx1+1)
                    Tdomain%specel(n)%sl%Veloc(i,j,k,2) = veloc(idx1+2)
                    idx1 = idx1 + 3
                    if ( .not. Tdomain%specel(n)%PML ) then
                        Tdomain%specel(n)%sl%Displ(i,j,k,0) = displ(idx2+0)
                        Tdomain%specel(n)%sl%Displ(i,j,k,1) = displ(idx2+1)
                        Tdomain%specel(n)%sl%Displ(i,j,k,2) = displ(idx2+2)
                        idx2 = idx2 + 3
                    else
                        Tdomain%specel(n)%slpml%Veloc1(i,j,k,0) = veloc1(idx3+0)
                        Tdomain%specel(n)%slpml%Veloc1(i,j,k,1) = veloc1(idx3+1)
                        Tdomain%specel(n)%slpml%Veloc1(i,j,k,2) = veloc1(idx3+2)
                        Tdomain%specel(n)%slpml%Veloc2(i,j,k,0) = veloc2(idx3+0)
                        Tdomain%specel(n)%slpml%Veloc2(i,j,k,1) = veloc2(idx3+1)
                        Tdomain%specel(n)%slpml%Veloc2(i,j,k,2) = veloc2(idx3+2)
                        Tdomain%specel(n)%slpml%Veloc3(i,j,k,0) = veloc3(idx3+0)
                        Tdomain%specel(n)%slpml%Veloc3(i,j,k,1) = veloc3(idx3+1)
                        Tdomain%specel(n)%slpml%Veloc3(i,j,k,2) = veloc3(idx3+2)
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

subroutine read_VelPhi(Tdomain, elem_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: elem_id
    real*8, allocatable, dimension(:) :: velphi, velphi1, velphi2, velphi3, phi
    integer :: idx1, idx2, idx3, ngllx, nglly, ngllz
    integer :: n, i, j, k


    call read_dset_1d_real(elem_id, "VelPhi", velphi)
    call read_dset_1d_real(elem_id, "VelPhi1", velphi1)
    call read_dset_1d_real(elem_id, "VelPhi2", velphi2)
    call read_dset_1d_real(elem_id, "VelPhi3", velphi3)
    call read_dset_1d_real(elem_id, "Phi", phi)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    do n = 0,Tdomain%n_elem-1
        if(Tdomain%specel(n)%solid) cycle
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        do k = 1,ngllz-2
            do j = 1,nglly-2
                do i = 1,ngllx-2
                    Tdomain%specel(n)%fl%VelPhi(i,j,k) = velphi(idx1)
                    idx1 = idx1 + 1
                    if ( .not. Tdomain%specel(n)%PML ) then
                        Tdomain%specel(n)%fl%Phi(i,j,k) = phi(idx2)
                        idx2 = idx2 + 1
                    else
                        Tdomain%specel(n)%flpml%VelPhi1(i,j,k) = velphi1(idx3)
                        Tdomain%specel(n)%flpml%VelPhi2(i,j,k) = velphi2(idx3)
                        Tdomain%specel(n)%flpml%VelPhi3(i,j,k) = velphi3(idx3)
                        idx3 = idx3 + 1
                    end if
                end do
            end do
        end do
    end do
    deallocate(velphi)
    deallocate(velphi1)
    deallocate(velphi2)
    deallocate(velphi3)
    deallocate(phi)
end subroutine read_VelPhi



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
    real*8, allocatable, dimension(:) :: eps_dev_xx, eps_dev_yy, eps_dev_xy, eps_dev_xz, eps_dev_yz
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

        if ( .not. Tdomain%specel(n)%PML ) then
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
    real*8, allocatable, dimension(:) :: stress
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
    real*8, allocatable, dimension(:) :: Veloc
    integer :: idx, ngllx, nglly, ngllz
    integer :: n, i, j, k

    call read_dset_1d_real(elem_id, "Veloc_Fl", veloc)
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
                        Tdomain%specel(n)%slpml%Veloc1(i,j,k,0) = veloc(idx+0)
                        Tdomain%specel(n)%slpml%Veloc1(i,j,k,1) = veloc(idx+1)
                        Tdomain%specel(n)%slpml%Veloc1(i,j,k,2) = veloc(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%slpml%Veloc2(i,j,k,0) = veloc(idx+0)
                        Tdomain%specel(n)%slpml%Veloc2(i,j,k,1) = veloc(idx+1)
                        Tdomain%specel(n)%slpml%Veloc2(i,j,k,2) = veloc(idx+2)
                        idx = idx + 3
                        Tdomain%specel(n)%slpml%Veloc3(i,j,k,0) = veloc(idx+0)
                        Tdomain%specel(n)%slpml%Veloc3(i,j,k,1) = veloc(idx+1)
                        Tdomain%specel(n)%slpml%Veloc3(i,j,k,2) = veloc(idx+2)
                        idx = idx + 3
                    end do
                end do
            end do
        end if
    end do
    deallocate(veloc)
end subroutine read_Veloc_Fluid_PML



subroutine read_Faces(Tdomain, face_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: face_id
    real*8, allocatable, dimension(:) :: veloc, veloc1, veloc2, veloc3, displ, &
                                         velphi, velphi1, velphi2, velphi3, phi
    integer :: idx1, idx2, idx3,idx4,idx5,idx6, ngll1, ngll2
    integer :: n, i, j


    call read_dset_1d_real(face_id, "Veloc", veloc)
    call read_dset_1d_real(face_id, "Veloc1", veloc1)
    call read_dset_1d_real(face_id, "Veloc2", veloc2)
    call read_dset_1d_real(face_id, "Veloc3", veloc3)
    call read_dset_1d_real(face_id, "Displ", displ)
    call read_dset_1d_real(face_id, "VelPhi", velphi)
    call read_dset_1d_real(face_id, "VelPhi1", velphi1)
    call read_dset_1d_real(face_id, "VelPhi2", velphi2)
    call read_dset_1d_real(face_id, "VelPhi3", velphi3)
    call read_dset_1d_real(face_id, "Phi", phi)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    idx4 = 1
    idx5 = 1
    idx6 = 1
    do n = 0,Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1
        ngll2 = Tdomain%sFace(n)%ngll2
        if(Tdomain%sFace(n)%solid)then
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
                    Tdomain%sFace(n)%spml%Veloc1(i,j,0) = veloc1(idx3+0)
                    Tdomain%sFace(n)%spml%Veloc1(i,j,1) = veloc1(idx3+1)
                    Tdomain%sFace(n)%spml%Veloc1(i,j,2) = veloc1(idx3+2)
                    Tdomain%sFace(n)%spml%Veloc2(i,j,0) = veloc2(idx3+0)
                    Tdomain%sFace(n)%spml%Veloc2(i,j,1) = veloc2(idx3+1)
                    Tdomain%sFace(n)%spml%Veloc2(i,j,2) = veloc2(idx3+2)
                    Tdomain%sFace(n)%spml%Veloc3(i,j,0) = veloc3(idx3+0)
                    Tdomain%sFace(n)%spml%Veloc3(i,j,1) = veloc3(idx3+1)
                    Tdomain%sFace(n)%spml%Veloc3(i,j,2) = veloc3(idx3+2)
                    idx3 = idx3 + 3
                end if
            end do
        end do
        else    ! fluid part
        do j = 1,ngll2-2
            do i = 1,ngll1-2
                Tdomain%sFace(n)%VelPhi(i,j) = velphi(idx4)
                idx4 = idx4 + 1
                if ( .not. Tdomain%sFace(n)%PML ) then
                    Tdomain%sFace(n)%Phi(i,j) = phi(idx5)
                    idx5 = idx5 + 1
                else
                    Tdomain%sFace(n)%spml%VelPhi1(i,j) = velphi1(idx6)
                    Tdomain%sFace(n)%spml%VelPhi2(i,j) = velphi2(idx6)
                    Tdomain%sFace(n)%spml%VelPhi3(i,j) = velphi3(idx6)
                    idx6 = idx6 + 1
                end if
            end do
        end do

        end if
    end do
    deallocate(veloc)
    deallocate(veloc1)
    deallocate(veloc2)
    deallocate(veloc3)
    deallocate(displ)
    deallocate(velphi)
    deallocate(velphi1)
    deallocate(velphi2)
    deallocate(velphi3)
    deallocate(phi)
end subroutine read_Faces

subroutine read_Edges(Tdomain, edge_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: edge_id
    real*8, allocatable, dimension(:) :: veloc, veloc1, veloc2, veloc3, displ, &
                                         velphi, velphi1, velphi2, velphi3, phi
    integer :: idx1, idx2, idx3, idx4,idx5,idx6,ngll
    integer :: n, i


    call read_dset_1d_real(edge_id, "Veloc", veloc)
    call read_dset_1d_real(edge_id, "Veloc1", veloc1)
    call read_dset_1d_real(edge_id, "Veloc2", veloc2)
    call read_dset_1d_real(edge_id, "Veloc3", veloc3)
    call read_dset_1d_real(edge_id, "Displ", displ)
    call read_dset_1d_real(edge_id, "VelPhi", velphi)
    call read_dset_1d_real(edge_id, "VelPhi1", velphi1)
    call read_dset_1d_real(edge_id, "VelPhi2", velphi2)
    call read_dset_1d_real(edge_id, "VelPhi3", velphi3)
    call read_dset_1d_real(edge_id, "Phi", phi)
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
                Tdomain%sEdge(n)%spml%Veloc1(i,0) = veloc1(idx3+0)
                Tdomain%sEdge(n)%spml%Veloc1(i,1) = veloc1(idx3+1)
                Tdomain%sEdge(n)%spml%Veloc1(i,2) = veloc1(idx3+2)
                Tdomain%sEdge(n)%spml%Veloc2(i,0) = veloc2(idx3+0)
                Tdomain%sEdge(n)%spml%Veloc2(i,1) = veloc2(idx3+1)
                Tdomain%sEdge(n)%spml%Veloc2(i,2) = veloc2(idx3+2)
                Tdomain%sEdge(n)%spml%Veloc3(i,0) = veloc3(idx3+0)
                Tdomain%sEdge(n)%spml%Veloc3(i,1) = veloc3(idx3+1)
                Tdomain%sEdge(n)%spml%Veloc3(i,2) = veloc3(idx3+2)
                idx3 = idx3 + 3
            end if
        end do
        else   ! fluid part
        do i = 1,ngll-2
            Tdomain%sEdge(n)%VelPhi(i) = velphi(idx4)
            idx4 = idx4 + 1
            if ( .not. Tdomain%sEdge(n)%PML ) then
                Tdomain%sEdge(n)%Phi(i) = phi(idx5)
                idx5 = idx5 + 1
            else
                Tdomain%sEdge(n)%spml%VelPhi1(i) = velphi1(idx6)
                Tdomain%sEdge(n)%spml%VelPhi2(i) = velphi2(idx6)
                Tdomain%sEdge(n)%spml%VelPhi3(i) = velphi3(idx6)
                idx6 = idx6 + 1
            end if
        end do

        end if
    end do
    deallocate(veloc)
    deallocate(veloc1)
    deallocate(veloc2)
    deallocate(veloc3)
    deallocate(displ)
    deallocate(velphi)
    deallocate(velphi1)
    deallocate(velphi2)
    deallocate(velphi3)
    deallocate(phi)
end subroutine read_Edges

subroutine read_Vertices(Tdomain, vertex_id)
    use sdomain
    use HDF5
    use sem_hdf5, only : read_dset_1d_real
    implicit none
    type (domain), intent (INOUT):: Tdomain
    integer(HID_T), intent(IN) :: vertex_id
    real*8, allocatable, dimension(:) :: veloc, veloc1, veloc2, veloc3, displ, &
                                         velphi, velphi1, velphi2, velphi3, phi
    integer :: idx1, idx2, idx3,idx4,idx5,idx6
    integer :: n


    call read_dset_1d_real(vertex_id, "Veloc", veloc)
    call read_dset_1d_real(vertex_id, "Veloc1", veloc1)
    call read_dset_1d_real(vertex_id, "Veloc2", veloc2)
    call read_dset_1d_real(vertex_id, "Veloc3", veloc3)
    call read_dset_1d_real(vertex_id, "Displ", displ)
    call read_dset_1d_real(vertex_id, "VelPhi", velphi)
    call read_dset_1d_real(vertex_id, "VelPhi1", velphi1)
    call read_dset_1d_real(vertex_id, "VelPhi2", velphi2)
    call read_dset_1d_real(vertex_id, "VelPhi3", velphi3)
    call read_dset_1d_real(vertex_id, "Phi", phi)
    idx1 = 1
    idx2 = 1
    idx3 = 1
    idx4 = 1
    idx5 = 1
    idx6 = 1
    do n = 0,Tdomain%n_vertex-1
        if(Tdomain%sVertex(n)%solid)then
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
            Tdomain%sVertex(n)%spml%Veloc1(0) = veloc1(idx3+0)
            Tdomain%sVertex(n)%spml%Veloc1(1) = veloc1(idx3+1)
            Tdomain%sVertex(n)%spml%Veloc1(2) = veloc1(idx3+2)
            Tdomain%sVertex(n)%spml%Veloc2(0) = veloc2(idx3+0)
            Tdomain%sVertex(n)%spml%Veloc2(1) = veloc2(idx3+1)
            Tdomain%sVertex(n)%spml%Veloc2(2) = veloc2(idx3+2)
            Tdomain%sVertex(n)%spml%Veloc3(0) = veloc3(idx3+0)
            Tdomain%sVertex(n)%spml%Veloc3(1) = veloc3(idx3+1)
            Tdomain%sVertex(n)%spml%Veloc3(2) = veloc3(idx3+2)
            idx3 = idx3 + 3
        end if
        else   ! fluid part
        Tdomain%sVertex(n)%VelPhi = velphi(idx4)
        idx4 = idx4 + 1
        if ( .not. Tdomain%sVertex(n)%PML ) then
            Tdomain%sVertex(n)%Phi = phi(idx5)
            idx5 = idx5 + 1
        else
            Tdomain%sVertex(n)%spml%VelPhi1 = velphi1(idx6)
            Tdomain%sVertex(n)%spml%VelPhi2 = velphi2(idx6)
            Tdomain%sVertex(n)%spml%VelPhi3 = velphi3(idx6)
            idx6 = idx6 + 1
        end if

        end if
    end do
    deallocate(veloc)
    deallocate(veloc1)
    deallocate(veloc2)
    deallocate(veloc3)
    deallocate(displ)
    deallocate(velphi)
    deallocate(velphi1)
    deallocate(velphi2)
    deallocate(velphi3)
    deallocate(phi)
end subroutine read_Vertices

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
    integer(HID_T) :: fid, elem_id, face_id, edge_id, vertex_id
    integer :: hdferr

    call init_hdf5()
    call init_restart(Tdomain%communicateur, rg, Tdomain%TimeD%iter_reprise, fnamef)

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
    call read_VelPhi(Tdomain, elem_id)
    call read_EpsilonVol(Tdomain, elem_id)
    call read_EpsilonDev(Tdomain, elem_id)
    call read_Stress(Tdomain, elem_id)
    call read_Veloc_Fluid_PML(Tdomain, elem_id)
    call read_Faces(Tdomain, face_id)
    call read_Edges(Tdomain, edge_id)
    call read_Vertices(Tdomain, vertex_id)

    call h5gclose_f(elem_id, hdferr)
    call h5gclose_f(face_id, hdferr)
    call h5gclose_f(edge_id, hdferr)
    call h5gclose_f(vertex_id, hdferr)
    call h5fclose_f(fid, hdferr)

    call clean_prot(Tdomain%TimeD%prot_m0, rg)
    return
end subroutine read_restart

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
