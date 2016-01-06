!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_fluid
    use constants
    use sdomain
    implicit none

contains

    subroutine forces_int_fluid(Elem, mat, htprimex, hprimey, htprimey, hprimez, htprimez,  &
        champs1)

        type (Element), intent (INOUT) :: Elem
        type (subdomain), intent(IN) :: mat
        real, dimension (0:Elem%ngllx-1, 0:Elem%ngllx-1), intent (IN) :: htprimex
        real, dimension (0:Elem%nglly-1, 0:Elem%nglly-1), intent (IN) :: hprimey, htprimey
        real, dimension (0:Elem%ngllz-1, 0:Elem%ngllz-1), intent (IN) :: hprimez, hTprimez
        type(champs), intent(inout) :: champs1

        integer :: m1,m2,m3, i,j,k
        real, dimension (0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dPhiX,dPhiY,dPhiZ,Fo_Fl

        real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: Phi


        m1 = Elem%ngllx;   m2 = Elem%nglly;   m3 = Elem%ngllz

        ! d(rho*Phi)_dX
        ! d(rho*Phi)_dY
        ! d(rho*Phi)_dZ
        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    Phi(i,j,k) = champs1%Phi(Elem%Idom(i,j,k))
                enddo
            enddo
        enddo
        call physical_part_deriv(m1,m2,m3,htprimex,hprimey,hprimez,Elem%InvGrad, Phi, &
            dPhiX, dPhiY, dPhiZ)

        ! internal forces
        call calcul_forces_fluid(Fo_Fl,                &
            Elem%Invgrad, &
            htprimex,htprimey,htprimez, &
            Elem%Jacob,mat%GLLwx,mat%GLLwy,mat%GLLwz, &
            dPhiX,dPhiY,dPhiZ,       &
            Elem%Density,            &
            m1,m2,m3)

        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    champs1%ForcesFl(Elem%Idom(i,j,k)) = champs1%ForcesFl(Elem%Idom(i,j,k))-Fo_Fl(i,j,k)
                enddo
            enddo
        enddo


        return
    end subroutine forces_int_fluid

end module dom_fluid

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
