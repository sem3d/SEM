module mfields
    use sdomain

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
            density*phi,dphi_dx,dphi_dy,dphi_dz)

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

    !! Recopie dans field le champs de deplacement reparti sur les Elem, Face, Edge, Vertex
    subroutine gather_elem_displ(Tdomain, nel, field)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:,0:,0:), intent(out) :: field
        type(element), pointer :: el
        type(face), pointer :: fc
        type(edge), pointer :: ed
        type(vertex), pointer :: vx
        integer :: nx, ny, nz, i
        nx = Tdomain%specel(nel)%ngllx
        ny = Tdomain%specel(nel)%nglly
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)

        if (el%solid) then
            field(1:nx-2,1:ny-2,1:nz-2,0:2) = el%Displ(:,:,:,:)
            do i=0,5
                fc => Tdomain%sFace(el%Near_Faces(i))
                call get_VectProperty_Face2Elem(i,el%Orient_Faces(i), nx, ny, nz, fc%ngll1, fc%ngll2, &
                    fc%Displ, field)
            end do

            do i=0,11
                ed => Tdomain%sEdge(el%Near_Edges(i))
                call get_VectProperty_Edge2Elem(i,el%Orient_Edges(i), nx, ny, nz, ed%ngll, &
                    ed%Displ, field)
            end do
            do i=0,7
                vx => Tdomain%sVertex(el%Near_Vertices(i))
                call get_VectProperty_Vertex2Elem(i, nx, ny, nz, vx%Displ, field)
            end do
        else ! liquid
            field = 0d0
        end if

    end subroutine gather_elem_displ

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

    subroutine gather_elem_veloc(Tdomain, nel, field)
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:,0:,0:), intent(out) :: field
        real, dimension(:,:,:), allocatable :: phi
        type(element), pointer :: el
        type(face), pointer :: fc
        type(edge), pointer :: ed
        type(vertex), pointer :: vx
        integer :: nx, ny, nz, i, mat
        nx = Tdomain%specel(nel)%ngllx
        ny = Tdomain%specel(nel)%nglly
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)
        if (el%solid) then
            field(1:nx-2,1:ny-2,1:nz-2,0:2) = el%Veloc(:,:,:,:)
            do i=0,5
                fc => Tdomain%sFace(el%Near_Faces(i))
                call get_VectProperty_Face2Elem(i,el%Orient_Faces(i), nx, ny, nz, fc%ngll1, fc%ngll2, &
                    fc%Veloc, field)
            end do

            do i=0,11
                ed => Tdomain%sEdge(el%Near_Edges(i))
                call get_VectProperty_Edge2Elem(i,el%Orient_Edges(i), nx, ny, nz, ed%ngll, &
                    ed%Veloc, field)
            end do
            do i=0,7
                vx => Tdomain%sVertex(el%Near_Vertices(i))
                call get_VectProperty_Vertex2Elem(i, nx, ny, nz, vx%Veloc, field)
            end do
        else ! liquid
            allocate(phi(0:nx-1,0:ny-1,0:nz-1))
            phi(1:nx-2,1:ny-2,1:nz-2) = el%Phi(:,:,:)
            do i=0,5
                fc => Tdomain%sFace(el%Near_Faces(i))
                call get_ScalarProperty_Face2Elem(i,el%Orient_Faces(i), nx, ny, nz, fc%ngll1, fc%ngll2, &
                    fc%Phi, phi)
            end do

            do i=0,11
                ed => Tdomain%sEdge(el%Near_Edges(i))
                call get_ScalarProperty_Edge2Elem(i,el%Orient_Edges(i), nx, ny, nz, ed%ngll, &
                    ed%Phi, phi)
            end do
            do i=0,7
                vx => Tdomain%sVertex(el%Near_Vertices(i))
                call get_ScalarProperty_Vertex2Elem(i, nx, ny, nz, vx%Phi, phi)
            end do
            mat = el%mat_index
            call fluid_velocity(nx,ny,nz,Tdomain%sSubdomain(mat)%htprimex,              &
                          Tdomain%sSubdomain(mat)%hprimey,Tdomain%sSubdomain(mat)%hprimez, &
                          el%InvGrad,el%density,phi,field)
            deallocate(phi)
        end if

    end subroutine gather_elem_veloc



end module mfields
