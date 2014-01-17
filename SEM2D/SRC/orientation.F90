module orientation

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

end module orientation
