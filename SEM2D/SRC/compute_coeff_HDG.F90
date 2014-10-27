!>
!!\file compute_coeff_HDG.F90
!!\brief Contains subroutines for computation of coefficients needed for HDG
!! methods for either explicit or implicit schemes. These coefficients are usually
!! defined on inter-element interfaces.
!!\version 1.0
!!\date 27/10/2014
!<

module scompute_coeff_HDG
    use selement
    use mpi
    use constants
    implicit none
contains

    ! ###########################################################
    !>
    !! \brief This subroutine computes the Penalization Matrix MatPen
    !! of an element on each external node of the element.
    !! It suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_MatPen (Elem)

        type (Element), intent (INOUT)   :: Elem
        real, dimension (0:Elem%ngllx-1) :: Zp_x, Zs_x
        real, dimension (0:Elem%ngllz-1) :: Zp_z, Zs_z
        integer    :: imin, imax, ngllx, ngllz

        ngllx = Elem%ngllx ; ngllz = Elem%ngllz

        ! Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        Zp_x(:) = sqrt(Elem%Density(0:ngllx-1,0) *(Elem%Lambda(0:ngllx-1,0)+2.*Elem%Mu(0:ngllx-1,0)))
        Zs_x(:) = sqrt(Elem%Density(0:ngllx-1,0) * Elem%Mu(0:ngllx-1,0))
        Elem%MatPen(imin:imax,0) = Zp_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zs_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,1) = Zs_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zp_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,2) =(Zp_x(:)-Zs_x(:))*Elem%Normal_Nodes(imin:imax,0) *Elem%Normal_Nodes(imin:imax,1)

        ! Right Face :
        call get_iminimax(Elem,1,imin,imax)
        Zp_z(:) = sqrt(Elem%Density(ngllx-1,0:ngllz-1) * (Elem%Lambda(ngllx-1,0:ngllz-1) &
                                                         + 2.*Elem%Mu(ngllx-1,0:ngllz-1)))
        Zs_z(:) = sqrt(Elem%Density(ngllx-1,0:ngllz-1) *  Elem%Mu(ngllx-1,0:ngllz-1))
        Elem%MatPen(imin:imax,0) = Zp_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zs_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,1) = Zs_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zp_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,2) =(Zp_z(:)-Zs_z(:))*Elem%Normal_Nodes(imin:imax,0) *Elem%Normal_Nodes(imin:imax,1)

        ! Top Face :
        call get_iminimax(Elem,2,imin,imax)
        Zp_x(:) = sqrt(Elem%Density(0:ngllx-1,ngllz-1) * (Elem%Lambda(0:ngllx-1,ngllz-1) &
                                                         + 2.*Elem%Mu(0:ngllx-1,ngllz-1)))
        Zs_x(:) = sqrt(Elem%Density(0:ngllx-1,ngllz-1) *  Elem%Mu(0:ngllx-1,ngllz-1))
        Elem%MatPen(imin:imax,0) = Zp_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zs_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,1) = Zs_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zp_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,2) =(Zp_x(:)-Zs_x(:))*Elem%Normal_Nodes(imin:imax,0) *Elem%Normal_Nodes(imin:imax,1)

        ! Left Face :
        call get_iminimax(Elem,3,imin,imax)
        Zp_z(:) = sqrt(Elem%Density(0,0:ngllz-1) *(Elem%Lambda(0,0:ngllz-1)+2.*Elem%Mu(0,0:ngllz-1)))
        Zs_z(:) = sqrt(Elem%Density(0,0:ngllz-1) * Elem%Mu(0,0:ngllz-1))
        Elem%MatPen(imin:imax,0) = Zp_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zs_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,1) = Zs_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                 + Zp_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
        Elem%MatPen(imin:imax,2) =(Zp_z(:)-Zs_z(:))*Elem%Normal_Nodes(imin:imax,0) *Elem%Normal_Nodes(imin:imax,1)

    end subroutine compute_MatPen

    ! ###########################################################
    !>
    !! \brief This subroutine computes the local matrix C.A^-1 which
    !! will be used to build the system on Lagrange multiplicators K * Lambda = R
    !! Indeed C.A^-1 is usefull for computing both K and R.
    !! This matrix is defined for each node lying on the element border.
    !! It suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !<
    subroutine  compute_CAinv (Elem)

        type (Element), intent (INOUT)   :: Elem

        

    end subroutine compute_CAinv

    ! ###########################################################
    !>
    !! \brief This subroutine computes the local matrix E.D^-1 which
    !! will be used to build the system on Lagrange multiplicators K * Lambda = R
    !! Indeed E.D^-1 is usefull for computing both K and R.
    !! This matrix is defined for each node lying on the element border.
    !! It suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (Element), intent (INOUT) Elem
    !<
    subroutine  compute_EDinv (Elem)

        type (Element), intent (INOUT)   :: Elem

        

    end subroutine compute_EDinv


end module scompute_coeff_HDG
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
