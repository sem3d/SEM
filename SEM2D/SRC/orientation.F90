module orientation
    use sdomain
contains

    subroutine get_VectProperty_Face2Elem(nf, orient, ngllx, ngllz, ngll, Sfield, Dfield)
        integer, intent(in)  :: nf       ! numero local de la face
        logical, intent(in)  :: orient   ! orientation
        integer, intent(in)  :: ngllx    ! ngllx de l'element cible
        integer, intent(in)  :: ngllz    ! ngllz de l'element cible
        integer, intent(in)  :: ngll     ! ngll de la face
        real, dimension(0:ngllx-1,0:ngllz-1,0:1), intent(inout) :: Dfield ! Valeur E/S du champ de l'elem
        real, dimension(1:ngll-2,0:1), intent(in) :: Sfield ! Valeur E du champ de la face
        select case(nf)
        case (0)
            if (orient) then
                Dfield(1:ngllx-2,0,:) =  Sfield(1:ngllx-2,:)
            else
                Dfield(1:ngllx-2,0,:) =  Sfield(ngllx-2:1:-1,:)
            endif
        case(1)
            if (orient) then
                Dfield(ngllx-1,1:ngllz-2,:) =  Sfield(1:ngllz-2,:)
            else
                Dfield(ngllx-1,1:ngllz-2,:) =  Sfield(ngllz-2:1:-1,:)
            endif
        case(2)
            if (orient) then
                Dfield(1:ngllx-2,ngllz-1,:) =  Sfield(1:ngllx-2,:)
            else
                Dfield(1:ngllx-2,ngllz-1,:) =  Sfield(ngllx-2:1:-1,:)
            endif
         case(3)
             if (orient) then
                 Dfield(0,1:ngllz-2,:) =  Sfield(1:ngllz-2,:)
             else
                 Dfield(0,1:ngllz-2,:) =  Sfield(ngllz-2:1:-1,:)
             endif
        end select
    end subroutine get_VectProperty_Face2Elem

    subroutine get_ScalProperty_Face2Elem(nf, orient, ngllx, ngllz, ngll, Sfield, Dfield)
        integer, intent(in)  :: nf       ! numero local de la face
        logical, intent(in)  :: orient   ! orientation
        integer, intent(in)  :: ngllx    ! ngllx de l'element cible
        integer, intent(in)  :: ngllz    ! ngllz de l'element cible
        integer, intent(in)  :: ngll     ! ngll de la face
        real, dimension(0:ngllx-1,0:ngllz-1), intent(inout) :: Dfield ! Valeur E/S du champ de l'elem
        real, dimension(1:ngll-2), intent(in) :: Sfield ! Valeur E du champ de la face
        select case(nf)
        case (0)
            if (orient) then
                Dfield(1:ngllx-2,0) =  Sfield(1:ngllx-2)
            else
                Dfield(1:ngllx-2,0) =  Sfield(ngllx-2:1:-1)
            endif
        case(1)
            if (orient) then
                Dfield(ngllx-1,1:ngllz-2) =  Sfield(1:ngllz-2)
            else
                Dfield(ngllx-1,1:ngllz-2) =  Sfield(ngllz-2:1:-1)
            endif
        case(2)
            if (orient) then
                Dfield(1:ngllx-2,ngllz-1) =  Sfield(1:ngllx-2)
            else
                Dfield(1:ngllx-2,ngllz-1) =  Sfield(ngllx-2:1:-1)
            endif
         case(3)
             if (orient) then
                 Dfield(0,1:ngllz-2) =  Sfield(1:ngllz-2)
             else
                 Dfield(0,1:ngllz-2) =  Sfield(ngllz-2:1:-1)
             endif
        end select
    end subroutine get_ScalProperty_Face2Elem


   !! Recopie dans field le champs de deplacement reparti sur les Elem, Face, Edge, Vertex
    subroutine gather_elem_displ(Tdomain, nel, field)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:,0:), intent(out) :: field
        type(element), pointer :: el
        type(face), pointer :: fc
        type(vertex), pointer :: vx
        integer :: nx, nz, i
        logical :: orient
        nx = Tdomain%specel(nel)%ngllx
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)

        field(1:nx-2,1:nz-2,0:1) = el%Displ(:,:,:)
        do i=0,3
            fc => Tdomain%sFace(el%Near_Face(i))
            orient = fc%coherency .or. fc%near_element(0) == nel
!            write(*,*) "EL", nel, "fc=",i, "orient:",orient, " nf=", el%Near_Face(i)
            call get_VectProperty_Face2Elem(i, orient, nx, nz, fc%ngll, fc%Displ, field)
        end do


        vx => Tdomain%sVertex(el%Near_Vertex(0))
        field(0,0,:) = vx%Displ
        vx => Tdomain%sVertex(el%Near_Vertex(1))
        field(nx-1,0,:) = vx%Displ
        vx => Tdomain%sVertex(el%Near_Vertex(2))
        field(nx-1,nz-1,:) = vx%Displ
        vx => Tdomain%sVertex(el%Near_Vertex(3))
        field(0,nz-1,:) = vx%Displ

    end subroutine gather_elem_displ

    subroutine gather_elem_veloc(Tdomain, nel, field)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:,0:), intent(out) :: field
        type(element), pointer :: el
        type(face), pointer :: fc
        type(vertex), pointer :: vx
        integer :: nx, nz, i, nx0, nx1, nz0, nz1
        logical :: orient
        nx = Tdomain%specel(nel)%ngllx
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)

        nx0 = 1
        nx1 = nx-2
        nz0 = 1
        nz1 = nz-2

        field(nx0:nx1,nz0:nz1,0:1) = el%Veloc(:,:,:)

        do i=0,3
            fc => Tdomain%sFace(el%Near_Face(i))
            orient = fc%coherency .or. fc%near_element(0) == nel
            call get_VectProperty_Face2Elem(i, orient, nx, nz, fc%ngll, fc%Veloc, field)
        end do

        vx => Tdomain%sVertex(el%Near_Vertex(0))
        field(0,0,:) = vx%Veloc
        vx => Tdomain%sVertex(el%Near_Vertex(1))
        field(nx-1,0,:) = vx%Veloc
        vx => Tdomain%sVertex(el%Near_Vertex(2))
        field(nx-1,nz-1,:) = vx%Veloc
        vx => Tdomain%sVertex(el%Near_Vertex(3))
        field(0,nz-1,:) = vx%Veloc

    end subroutine gather_elem_veloc

    subroutine gather_elem_accel(Tdomain, nel, field)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:,0:), intent(out) :: field
        type(element), pointer :: el
        type(face), pointer :: fc
        type(vertex), pointer :: vx
        integer :: nx, nz, i, nx0, nx1, nz0, nz1
        logical :: orient
        nx = Tdomain%specel(nel)%ngllx
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)

        nx0 = 1
        nx1 = nx-2
        nz0 = 1
        nz1 = nz-2

        field(nx0:nx1,nz0:nz1,0:1) = el%Accel(:,:,:)

        do i=0,3
            fc => Tdomain%sFace(el%Near_Face(i))
            orient = fc%coherency .or. fc%near_element(0) == nel
            call get_VectProperty_Face2Elem(i, orient, nx, nz, fc%ngll, fc%Accel, field)
        end do

        vx => Tdomain%sVertex(el%Near_Vertex(0))
        field(0,0,:) = vx%Accel
        vx => Tdomain%sVertex(el%Near_Vertex(1))
        field(nx-1,0,:) = vx%Accel
        vx => Tdomain%sVertex(el%Near_Vertex(2))
        field(nx-1,nz-1,:) = vx%Accel
        vx => Tdomain%sVertex(el%Near_Vertex(3))
        field(0,nz-1,:) = vx%Accel

    end subroutine gather_elem_accel

    subroutine gather_elem_mass(Tdomain, nel, field)
        implicit none
        type(domain), intent(in) :: Tdomain
        integer, intent(in) :: nel
        real, dimension(0:,0:), intent(out) :: field
        type(element), pointer :: el
        type(face), pointer :: fc
        type(vertex), pointer :: vx
        integer :: nx, nz, i, nx0, nx1, nz0, nz1
        logical :: orient
        nx = Tdomain%specel(nel)%ngllx
        nz = Tdomain%specel(nel)%ngllz
        el => Tdomain%specel(nel)

        nx0 = 1
        nx1 = nx-2
        nz0 = 1
        nz1 = nz-2

        field(nx0:nx1,nz0:nz1) = el%MassMat(:,:)

        do i=0,3
            fc => Tdomain%sFace(el%Near_Face(i))
            orient = fc%coherency .or. fc%near_element(0) == nel
            call get_ScalProperty_Face2Elem(i, orient, nx, nz, fc%ngll, fc%MassMat, field)
        end do

        vx => Tdomain%sVertex(el%Near_Vertex(0))
        field(0,0) = vx%MassMat
        vx => Tdomain%sVertex(el%Near_Vertex(1))
        field(nx-1,0) = vx%MassMat
        vx => Tdomain%sVertex(el%Near_Vertex(2))
        field(nx-1,nz-1) = vx%MassMat
        vx => Tdomain%sVertex(el%Near_Vertex(3))
        field(0,nz-1) = vx%MassMat

    end subroutine gather_elem_mass

end module orientation
