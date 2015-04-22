module mfields
    use sdomain
    use orientation
    use deriv3d
contains

    subroutine fluid_velocity(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,InvGrad,    &
        density,phi,veloc)
        ! gives the physical particle velocity in the fluid = 1/dens grad(dens.Phi)
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: density,phi
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(out) :: Veloc
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: dphi_dx,dphi_dy,dphi_dz

        ! physical gradient
        call physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
            phi,dphi_dx,dphi_dy,dphi_dz)

        !
        Veloc(:,:,:,0) = dphi_dx(:,:,:)/density(:,:,:) 
        Veloc(:,:,:,1) = dphi_dy(:,:,:)/density(:,:,:) 
        Veloc(:,:,:,2) = dphi_dz(:,:,:)/density(:,:,:) 

    end subroutine fluid_velocity

    !---------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------

    subroutine grad_displ_solid(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,InvGrad,displ,    &
        dUx_dx,dUx_dy,dUx_dz,dUy_dx,dUy_dy,dUy_dz,dUz_dx,dUz_dy,dUz_dz)
        ! gives the displacement gradient in the solid part
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(in) :: displ 
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: dUx_dx,dUx_dy,dUx_dz,  &
            dUy_dx,dUy_dy,dUy_dz,dUz_dx,dUz_dy,dUz_dz

        call physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
            displ(:,:,:,0),dUx_dx,dUx_dy,dUx_dz)
        call physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
            displ(:,:,:,1),dUy_dx,dUy_dy,dUy_dz)
        call physical_part_deriv(ngllx,nglly,ngllz,hTprimex,hprimey,hprimez,InvGrad,   &
            displ(:,:,:,2),dUz_dx,dUz_dy,dUz_dz)


    end subroutine grad_displ_solid

    !---------------------------------------------------------------------------------------
    !---------------------------------------------------------------------------------------

    subroutine pressure_solid(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,InvGrad,displ,    &
        lambda,mu,pressure)
        ! gives the isotropic term of the stress tensor = pressure = -1/3 trace(stress tensor)
        implicit none
        integer, intent(in)  :: ngllx,nglly,ngllz
        real, dimension(0:ngllx-1,0:ngllx-1), intent(in) :: hTprimex
        real, dimension(0:nglly-1,0:nglly-1), intent(in) :: hprimey
        real, dimension(0:ngllz-1,0:ngllz-1), intent(in) :: hprimez
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:2), intent(in) :: InvGrad
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2), intent(in) :: displ 
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: mu,lambda
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: pressure
        real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1) :: dUx_dx,dUx_dy,dUx_dz,  &
            dUy_dx,dUy_dy,dUy_dz,dUz_dx,dUz_dy,dUz_dz


        call grad_displ_solid(ngllx,nglly,ngllz,htprimex,hprimey,hprimez,InvGrad,displ,    &
            dUx_dx,dUx_dy,dUx_dz,dUy_dx,dUy_dy,dUy_dz,dUz_dx,dUz_dy,dUz_dz)

        pressure(:,:,:) = -(lambda(:,:,:)+2d0/3d0*mu(:,:,:))*(dUx_dx(:,:,:)+dUy_dy(:,:,:)+dUz_dz(:,:,:))

    end subroutine pressure_solid


    subroutine check_field(nel, field, nx, ny, nz)
        integer, intent(in) :: nel, nx, ny, nz
        real, dimension(0:nx-1,0:ny-1,0:nz-1), intent(in) :: field
        !
        integer :: i,j,k
        do i=0,nx-1
            do j=0,ny-1
                do k=0,nz-1
                    if (field(i,j,k)>1e2) then
                        write(*,*) "Elem:",nel, "(",i,j,k,")"
                    end if
                end do
            end do
        end do
    end subroutine check_field

    subroutine gather_field(el, field, src_field, inum)
        type(element), intent(in), pointer :: el
        real, dimension(0:,0:,0:,0:), intent(out) :: field
        real, dimension(0:,0:), intent(in) :: src_field
        integer, dimension(0:,0:,0:), intent(in) :: inum
        !
        integer :: i,j,k

        do k=0,el%ngllz-1
            do j=0,el%nglly-1
                do i=0,el%ngllx-1
                    ind = inum(i,j,k)
                    field(i,j,k,:) = src_field(ind,:) 
                enddo
            enddo
        enddo
    end subroutine gather_field

    subroutine gather_field_pml(el, field, src_field, inum)
        type(element), intent(in), pointer :: el
        real, dimension(0:,0:,0:,0:), intent(out) :: field
        real, dimension(0:,0:), intent(in) :: src_field
        integer, dimension(0:,0:,0:), intent(in) :: inum
        !
        integer :: i,j,k

        do k=0,el%ngllz-1
            do j=0,el%nglly-1
                do i=0,el%ngllx-1
                    ind = inum(i,j,k)
                    field(i,j,k,:) = src_field(ind,:) + src_field(ind+1,:) + src_field(ind+2,:)
                enddo
            enddo
        enddo
    end subroutine gather_field_pml

    subroutine gather_field_fpml(el, field, src_field, inum)
        type(element), intent(in), pointer :: el
        real, dimension(0:,0:,0:), intent(out) :: field
        real, dimension(0:), intent(in) :: src_field
        integer, dimension(0:,0:,0:), intent(in) :: inum
        !
        integer :: i,j,k

        do k=0,el%ngllz-1
            do j=0,el%nglly-1
                do i=0,el%ngllx-1
                    ind = inum(i,j,k)
                    field(i,j,k) = src_field(ind) + src_field(ind+1) + src_field(ind+2)
                enddo
            enddo
        enddo
    end subroutine gather_field_fpml

    subroutine gather_field_fluid(el, field, src_field, inum)
        type(element), intent(in), pointer :: el
        real, dimension(0:,0:,0:), intent(out) :: field
        real, dimension(0:), intent(in) :: src_field
        integer, dimension(0:,0:,0:), intent(in) :: inum
        !
        integer :: i,j,k

        do k=0,el%ngllz-1
            do j=0,el%nglly-1
                do i=0,el%ngllx-1
                    ind = inum(i,j,k)
                    field(i,j,k) = src_field(ind) 
                enddo
            enddo
        enddo
    end subroutine gather_field_fluid

    subroutine gather_elem_veloc(Tdomain, nel, field)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:,0:,0:), intent(out) :: field
        real, dimension(:,:,:), allocatable :: phi
        type(element), pointer :: el
        integer :: nx, ny, nz, i, j, k, mat, ind
        nx = Tdomain%specel(nel)%ngllx
        ny = Tdomain%specel(nel)%nglly
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)
        if (el%solid) then
            if (el%PML) then
                do k=0,nz-1
                    do j=0,ny-1
                        do i=0,nx-1
                            ind = el%slpml%ISolPML(i,j,k)
                            field(i,j,k,:) = Tdomain%champs0%VelocPml(ind,:) + &
                                             Tdomain%champs0%VelocPml(ind+1,:) + &
                                             Tdomain%champs0%VelocPml(ind+2,:)
                        enddo
                    enddo
                enddo
            else
                do k=0,nz-1
                    do j=0,ny-1
                        do i=0,nx-1
                            ind = el%ISol(i,j,k)
                            field(i,j,k,:) = Tdomain%champs0%Veloc(ind,:)
                        enddo
                    enddo
                enddo
            endif
        else ! liquid
            if (el%PML) then
                field = 0d0
            else
                allocate(phi(0:nx-1,0:ny-1,0:nz-1))
                do k=0,nz-1
                    do j=0,ny-1
                        do i=0,nx-1
                            ind = el%IFlu(i,j,k)
                            phi(i,j,k) = Tdomain%champs0%Phi(ind)
                        enddo
                    enddo
                enddo
                mat = el%mat_index
                call fluid_velocity(nx,ny,nz,Tdomain%sSubdomain(mat)%htprimex,              &
                            Tdomain%sSubdomain(mat)%hprimey,Tdomain%sSubdomain(mat)%hprimez, &
                            el%InvGrad,el%density,phi,field)
                deallocate(phi)
            endif
        end if

    end subroutine gather_elem_veloc

    subroutine gather_elem_displ(Tdomain, nel, field)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:,0:,0:), intent(out) :: field
        type(element), pointer :: el
        integer :: nx, ny, nz, i, j, k, ind
        nx = Tdomain%specel(nel)%ngllx
        ny = Tdomain%specel(nel)%nglly
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)

        if (el%solid) then
            if (el%PML) then
                field = 0d0
            else
                do k=0,nz-1
                    do j=0,ny-1
                        do i=0,nx-1
                            ind = el%ISol(i,j,k)
                            field(i,j,k,:) = Tdomain%champs0%Depla(ind,:)
                        enddo
                    enddo
                enddo
            endif
        else ! liquid
            field = 0d0
        end if

    end subroutine gather_elem_displ

    subroutine gather_elem_accel(Tdomain, nel, field)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:,0:,0:), intent(out) :: field
        real, dimension(:,:,:), allocatable :: vphi
        type(element), pointer :: el
        integer :: nx, ny, nz, i, j, k, ind
        nx = Tdomain%specel(nel)%ngllx
        ny = Tdomain%specel(nel)%nglly
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)
        if (el%solid) then
            if (el%PML) then
                do k=0,nz-1
                    do j=0,ny-1
                        do i=0,nx-1
                            ind = el%slpml%ISolPML(i,j,k)
                            field(i,j,k,:) = Tdomain%champs1%ForcesPML(ind,:) + &
                                             Tdomain%champs1%ForcesPML(ind+1,:) + &
                                             Tdomain%champs1%ForcesPML(ind+2,:)
                        enddo
                    enddo
                enddo
            else
                do k=0,nz-1
                    do j=0,ny-1
                        do i=0,nx-1
                            ind = el%ISol(i,j,k)
                            field(i,j,k,:) = Tdomain%champs0%Forces(ind,:)
                        enddo
                    enddo
                enddo
            endif
        else ! liquid
            if (el%PML) then
                field = 0d0
            else
                allocate(vphi(0:nx-1,0:ny-1,0:nz-1))
                do k=0,nz-1
                    do j=0,ny-1
                        do i=0,nx-1
                            ind = el%IFlu(i,j,k)
                            vphi(i,j,k) = Tdomain%champs0%VelPhi(ind)
                        enddo
                    enddo
                enddo
                mat = el%mat_index
                call fluid_velocity(nx,ny,nz,Tdomain%sSubdomain(mat)%htprimex,              &
                            Tdomain%sSubdomain(mat)%hprimey,Tdomain%sSubdomain(mat)%hprimez, &
                            el%InvGrad,el%density,vphi,field)
                deallocate(vphi)
            endif
        end if
    end subroutine gather_elem_accel


    subroutine gather_elem_press(Tdomain, nel, field)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:,0:), intent(out) :: field
        real, dimension(:,:,:,:), allocatable :: displ
        type(element), pointer :: el
        integer :: nx, ny, nz, i, mat
        nx = Tdomain%specel(nel)%ngllx
        ny = Tdomain%specel(nel)%nglly
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)
        if (el%solid) then
            allocate(displ(0:nx-1,0:ny-1,0:nz-1,0:2))
            mat = el%mat_index
            call gather_elem_displ(Tdomain, nel, displ)
            call pressure_solid(nx,ny,nz,Tdomain%sSubdomain(mat)%htprimex,              &
                 Tdomain%sSubdomain(mat)%hprimey,Tdomain%sSubdomain(mat)%hprimez, &
                 el%InvGrad, displ, el%Lambda, el%Mu,field)
            deallocate(displ)
        else ! liquid
            do k=0,nz-1
                do j=0,ny-1
                    do i=0,nx-1
                        ind = el%IFlu(i,j,k)
                        field(i,j,k) = -Tdomain%champs0%VelPhi(ind)
                    enddo
                enddo
            enddo
        end if

    end subroutine gather_elem_press

end module mfields
