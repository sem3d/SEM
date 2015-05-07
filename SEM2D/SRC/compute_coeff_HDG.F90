!>
!!\file compute_coeff_HDG.F90
!!\brief Contains subroutines for computation of coefficients needed for HDG
!! methods for either explicit or implicit schemes. These coefficients are usually
!! defined on inter-element interfaces.
!!\version 1.0
!!\date 27/10/2014
!<

module scompute_coeff_HDG
    use sdomain
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
    !! IMPORTANT NOTE : TO OBTAIN THE WILCOX'S METHOD'S RESULTS USING HDG,
    !! THE Elem%MatPen SHOULD NEVER BE COMPUTED USING THE ACOUSTIC FORMULA.
    !! \param type (Element), intent (INOUT) Elem
    !!
    !<
    subroutine  compute_MatPen (Tdomain,nelem)

        type (Domain), intent (INOUT) :: Tdomain
        integer,       intent (IN)    :: nelem
        type(element), pointer :: Elem
        real, dimension (0:Tdomain%specel(nelem)%ngllx-1) :: Zp_x, Zs_x
        real, dimension (0:Tdomain%specel(nelem)%ngllz-1) :: Zp_z, Zs_z
        integer    :: imin, imax, ngllx, ngllz, nf

        Elem => Tdomain%specel(nelem)
        ngllx = Elem%ngllx ; ngllz = Elem%ngllz

        ! Bottom Face :
        call get_iminimax(Elem,0,imin,imax)
        nf = Elem%Near_Face(0)
        Zp_x(:) = sqrt(Elem%Density(0:ngllx-1,0) *(Elem%Lambda(0:ngllx-1,0)+2.*Elem%Mu(0:ngllx-1,0)))
        Zs_x(:) = sqrt(Elem%Density(0:ngllx-1,0) * Elem%Mu(0:ngllx-1,0))
        if (Elem%Acoustic .AND. (.NOT. Tdomain%sFace(nf)%changing_media)) then
            Elem%MatPen(imin:imax,0) = Zp_x(:)
            Elem%MatPen(imin:imax,1) = Zp_x(:)
            Elem%MatPen(imin:imax,2) = 0.
        else
            Elem%MatPen(imin:imax,0) = Zp_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                     + Zs_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
            Elem%MatPen(imin:imax,1) = Zs_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                     + Zp_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
            Elem%MatPen(imin:imax,2) =(Zp_x(:)-Zs_x(:))*Elem%Normal_Nodes(imin:imax,0) &
                                     * Elem%Normal_Nodes(imin:imax,1)
        endif

        ! Right Face :
        call get_iminimax(Elem,1,imin,imax)
        nf = Elem%Near_Face(1)
        Zp_z(:) = sqrt(Elem%Density(ngllx-1,0:ngllz-1) * (Elem%Lambda(ngllx-1,0:ngllz-1) &
                                                         + 2.*Elem%Mu(ngllx-1,0:ngllz-1)))
        Zs_z(:) = sqrt(Elem%Density(ngllx-1,0:ngllz-1) *  Elem%Mu(ngllx-1,0:ngllz-1))
        if (Elem%Acoustic .AND. (.NOT. Tdomain%sFace(nf)%changing_media)) then
            Elem%MatPen(imin:imax,0) = Zp_z(:)
            Elem%MatPen(imin:imax,1) = Zp_z(:)
            Elem%MatPen(imin:imax,2) = 0.
        else
            Elem%MatPen(imin:imax,0) = Zp_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                     + Zs_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
            Elem%MatPen(imin:imax,1) = Zs_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                     + Zp_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
            Elem%MatPen(imin:imax,2) =(Zp_z(:)-Zs_z(:))*Elem%Normal_Nodes(imin:imax,0) &
                                     * Elem%Normal_Nodes(imin:imax,1)
        endif

        ! Top Face :
        call get_iminimax(Elem,2,imin,imax)
        nf = Elem%Near_Face(2)
        Zp_x(:) = sqrt(Elem%Density(0:ngllx-1,ngllz-1) * (Elem%Lambda(0:ngllx-1,ngllz-1) &
                                                         + 2.*Elem%Mu(0:ngllx-1,ngllz-1)))
        Zs_x(:) = sqrt(Elem%Density(0:ngllx-1,ngllz-1) *  Elem%Mu(0:ngllx-1,ngllz-1))
        if (Elem%Acoustic .AND. (.NOT. Tdomain%sFace(nf)%changing_media)) then
            Elem%MatPen(imin:imax,0) = Zp_x(:)
            Elem%MatPen(imin:imax,1) = Zp_x(:)
            Elem%MatPen(imin:imax,2) = 0.
        else
            Elem%MatPen(imin:imax,0) = Zp_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                     + Zs_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
            Elem%MatPen(imin:imax,1) = Zs_x(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                     + Zp_x(:)*Elem%Normal_Nodes(imin:imax,1)**2
            Elem%MatPen(imin:imax,2) =(Zp_x(:)-Zs_x(:))*Elem%Normal_Nodes(imin:imax,0) &
                                     * Elem%Normal_Nodes(imin:imax,1)
        endif

        ! Left Face :
        call get_iminimax(Elem,3,imin,imax)
        nf = Elem%Near_Face(3)
        Zp_z(:) = sqrt(Elem%Density(0,0:ngllz-1) *(Elem%Lambda(0,0:ngllz-1)+2.*Elem%Mu(0,0:ngllz-1)))
        Zs_z(:) = sqrt(Elem%Density(0,0:ngllz-1) * Elem%Mu(0,0:ngllz-1))
        if (Elem%Acoustic .AND. (.NOT. Tdomain%sFace(nf)%changing_media)) then
            Elem%MatPen(imin:imax,0) = Zp_z(:)
            Elem%MatPen(imin:imax,1) = Zp_z(:)
            Elem%MatPen(imin:imax,2) = 0.
        else
            Elem%MatPen(imin:imax,0) = Zp_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                     + Zs_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
            Elem%MatPen(imin:imax,1) = Zs_z(:)*Elem%Normal_Nodes(imin:imax,0)**2 &
                                     + Zp_z(:)*Elem%Normal_Nodes(imin:imax,1)**2
            Elem%MatPen(imin:imax,2) =(Zp_z(:)-Zs_z(:))*Elem%Normal_Nodes(imin:imax,0) &
                                     * Elem%Normal_Nodes(imin:imax,1)
        endif

    end subroutine compute_MatPen


    ! ###########################################################
    !>
    !! \brief This subroutine computes the "G-part" of matrix K on a given face
    !! (for the system on Lagrange multiplicators K * Lambda = R)
    !! on the interiors nodes of the face only. It computes matrix G at
    !! elementary level, and then dend its contributiuon to neighbour faces.
    !! It suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (domain), intent (INOUT) Tdomain
    !! \param integer, intent (IN) nelem
    !<
    subroutine build_K_expl(Tdomain, nelem)

        implicit none
        type (domain), intent (INOUT) :: Tdomain
        integer, intent(IN) :: nelem
        real, dimension(0:2*(Tdomain%specel(nelem)%ngllx+Tdomain%specel(nelem)%ngllz)-1,0:2) :: G
        type(element), pointer :: Elem
        integer :: nf, nface, i, imin, imax
        logical :: coherency

        Elem => Tdomain%specel(nelem)

        ! Calcul du terme G qui contribue e la matrice K
        G(:,0) = Elem%Coeff_integr_Faces(:) * Elem%MatPen(:,0)
        G(:,1) = Elem%Coeff_integr_Faces(:) * Elem%MatPen(:,1)
        G(:,2) = Elem%Coeff_integr_Faces(:) * Elem%MatPen(:,2)

        ! Envoi des matrices sur les faces :
        do nf=0,3
            nface = Elem%Near_Face(nf)
            call get_iminimax(Elem,nf,imin,imax)
            coherency  = Tdomain%sFace(nface)%coherency
            if (coherency .OR. (Tdomain%sFace(nface)%Near_Element(0)==nelem)) then
                Tdomain%sFace(nface)%KinvExpl(:,:) = Tdomain%sFace(nface)%KinvExpl(:,:) + G(imin:imax,:)
            else
                do i=0,Tdomain%sFace(nface)%ngll-1
                    Tdomain%sFace(nface)%KinvExpl(i,:) = Tdomain%sFace(nface)%KinvExpl(i,:) + G(imax-i,:)
                end do
            endif
        enddo

    end subroutine build_K_expl

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

        type (Element), intent (INOUT) :: Elem
        integer :: imin, imax, ngx, ngz

        ngx = Elem%ngllx ; ngz = Elem%ngllz

        Elem%CAinv(:,0,0) = Elem%Coeff_Integr_Faces(:) * Elem%Normal_Nodes(:,0)
        Elem%CAinv(:,0,1) = Elem%Coeff_Integr_Faces(:) * Elem%Normal_Nodes(:,0)
        Elem%CAinv(:,0,2) = Elem%Coeff_Integr_Faces(:) * Elem%Normal_Nodes(:,1)
        Elem%CAinv(:,1,0) = Elem%Coeff_Integr_Faces(:) * Elem%Normal_Nodes(:,1)
        Elem%CAinv(:,1,1) = Elem%Coeff_Integr_Faces(:) * Elem%Normal_Nodes(:,1)
        Elem%CAinv(:,1,2) = Elem%Coeff_Integr_Faces(:) * Elem%Normal_Nodes(:,0)

        call get_iminimax(Elem,0,imin,imax)
        Elem%CAinv(imin:imax,0,0) = (Elem%Lambda(0:ngx-1,0) + 2*Elem%Mu(0:ngx-1,0)) &
                                   / Elem%Acoeff(0:ngx-1,0,12) * Elem%CAinv(imin:imax,0,0)
        Elem%CAinv(imin:imax,0,1) = (Elem%Lambda(0:ngx-1,0)) &
                                   / Elem%Acoeff(0:ngx-1,0,12) * Elem%CAinv(imin:imax,0,1)
        Elem%CAinv(imin:imax,0,2) = (Elem%Mu(0:ngx-1,0)) &
                                   / Elem%Acoeff(0:ngx-1,0,12) * Elem%CAinv(imin:imax,0,2)
        Elem%CAinv(imin:imax,1,0) = (Elem%Lambda(0:ngx-1,0)) &
                                   / Elem%Acoeff(0:ngx-1,0,12) * Elem%CAinv(imin:imax,1,0)
        Elem%CAinv(imin:imax,1,1) = (Elem%Lambda(0:ngx-1,0) + 2*Elem%Mu(0:ngx-1,0)) &
                                   / Elem%Acoeff(0:ngx-1,0,12) * Elem%CAinv(imin:imax,1,1)
        Elem%CAinv(imin:imax,1,2) = (Elem%Mu(0:ngx-1,0)) &
                                   / Elem%Acoeff(0:ngx-1,0,12) * Elem%CAinv(imin:imax,1,2)
        call get_iminimax(Elem,1,imin,imax)
        Elem%CAinv(imin:imax,0,0) = (Elem%Lambda(ngx-1,0:ngz-1) + 2*Elem%Mu(ngx-1,0:ngz-1)) &
                                   / Elem%Acoeff(ngx-1,0:ngz-1,12) * Elem%CAinv(imin:imax,0,0)
        Elem%CAinv(imin:imax,0,1) = (Elem%Lambda(ngx-1,0:ngz-1)) &
                                   / Elem%Acoeff(ngx-1,0:ngz-1,12) * Elem%CAinv(imin:imax,0,1)
        Elem%CAinv(imin:imax,0,2) = (Elem%Mu(ngx-1,0:ngz-1)) &
                                   / Elem%Acoeff(ngx-1,0:ngz-1,12) * Elem%CAinv(imin:imax,0,2)
        Elem%CAinv(imin:imax,1,0) = (Elem%Lambda(ngx-1,0:ngz-1)) &
                                   / Elem%Acoeff(ngx-1,0:ngz-1,12) * Elem%CAinv(imin:imax,1,0)
        Elem%CAinv(imin:imax,1,1) = (Elem%Lambda(ngx-1,0:ngz-1) + 2*Elem%Mu(ngx-1,0:ngz-1)) &
                                   / Elem%Acoeff(ngx-1,0:ngz-1,12) * Elem%CAinv(imin:imax,1,1)
        Elem%CAinv(imin:imax,1,2) = (Elem%Mu(ngx-1,0:ngz-1)) &
                                   / Elem%Acoeff(ngx-1,0:ngz-1,12) * Elem%CAinv(imin:imax,1,2)
        call get_iminimax(Elem,2,imin,imax)
        Elem%CAinv(imin:imax,0,0) = (Elem%Lambda(0:ngx-1,ngz-1) + 2*Elem%Mu(0:ngx-1,ngz-1)) &
                                   / Elem%Acoeff(0:ngx-1,ngz-1,12) * Elem%CAinv(imin:imax,0,0)
        Elem%CAinv(imin:imax,0,1) = (Elem%Lambda(0:ngx-1,ngz-1)) &
                                   / Elem%Acoeff(0:ngx-1,ngz-1,12) * Elem%CAinv(imin:imax,0,1)
        Elem%CAinv(imin:imax,0,2) = (Elem%Mu(0:ngx-1,ngz-1)) &
                                   / Elem%Acoeff(0:ngx-1,ngz-1,12) * Elem%CAinv(imin:imax,0,2)
        Elem%CAinv(imin:imax,1,0) = (Elem%Lambda(0:ngx-1,ngz-1)) &
                                   / Elem%Acoeff(0:ngx-1,ngz-1,12) * Elem%CAinv(imin:imax,1,0)
        Elem%CAinv(imin:imax,1,1) = (Elem%Lambda(0:ngx-1,ngz-1) + 2*Elem%Mu(0:ngx-1,ngz-1)) &
                                   / Elem%Acoeff(0:ngx-1,ngz-1,12) * Elem%CAinv(imin:imax,1,1)
        Elem%CAinv(imin:imax,1,2) = (Elem%Mu(0:ngx-1,ngz-1)) &
                                   / Elem%Acoeff(0:ngx-1,ngz-1,12) * Elem%CAinv(imin:imax,1,2)
        call get_iminimax(Elem,3,imin,imax)
        Elem%CAinv(imin:imax,0,0) = (Elem%Lambda(0,0:ngz-1) + 2*Elem%Mu(0,0:ngz-1)) &
                                   / Elem%Acoeff(0,0:ngz-1,12) * Elem%CAinv(imin:imax,0,0)
        Elem%CAinv(imin:imax,0,1) = (Elem%Lambda(0,0:ngz-1)) &
                                   / Elem%Acoeff(0,0:ngz-1,12) * Elem%CAinv(imin:imax,0,1)
        Elem%CAinv(imin:imax,0,2) = (Elem%Mu(0,0:ngz-1)) &
                                   / Elem%Acoeff(0,0:ngz-1,12) * Elem%CAinv(imin:imax,0,2)
        Elem%CAinv(imin:imax,1,0) = (Elem%Lambda(0,0:ngz-1)) &
                                   / Elem%Acoeff(0,0:ngz-1,12) * Elem%CAinv(imin:imax,1,0)
        Elem%CAinv(imin:imax,1,1) = (Elem%Lambda(0,0:ngz-1) + 2*Elem%Mu(0,0:ngz-1)) &
                                   / Elem%Acoeff(0,0:ngz-1,12) * Elem%CAinv(imin:imax,1,1)
        Elem%CAinv(imin:imax,1,2) = (Elem%Mu(0,0:ngz-1)) &
                                   / Elem%Acoeff(0,0:ngz-1,12) * Elem%CAinv(imin:imax,1,2)

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
    subroutine  compute_EDinv (Elem, Dt)

        type (Element), intent (INOUT)   :: Elem
        real, intent(IN) :: Dt
        integer :: imin, imax, ngx, ngz, nc, n1, n2
        real, dimension(0:2*(Elem%ngllx+Elem%ngllz)-1,0:1,0:1) :: matD
        real, dimension(0:2*(Elem%ngllx+Elem%ngllz)-1) :: det, tmp

        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! Calcul de la matrice D :
        ! Attention la matrice D ne devrait pas etre dedoublee aux noeuds des coin
        ! Ici cependant elle est dedoublee car il faudra la multiplier avec E qui, elle,
        ! est en effet dedoublee. Du coup les 2 entres de D pour un coin ont la meme valeur.

        ! Termes provenant du terme tau=matpen
        matD(:,0,0) = 0.5*Dt * Elem%Coeff_Integr_Faces(:) * Elem%MatPen(:,0)
        matD(:,1,1) = 0.5*Dt * Elem%Coeff_Integr_Faces(:) * Elem%MatPen(:,1)
        matD(:,0,1) = 0.5*Dt * Elem%Coeff_Integr_Faces(:) * Elem%MatPen(:,2)
        matD(:,1,0) = 0.5*Dt * Elem%Coeff_Integr_Faces(:) * Elem%MatPen(:,2)
        ! Couplage aux coins :
        do nc=0,3
            call get_gll_arround_corner(Elem,nc,n1,n2)
            matD(n1,:,:) = matD(n1,:,:) + matD(n2,:,:)
            matD(n2,:,:) = matD(n1,:,:)
        enddo
        !matD = 0.

        ! Termes provenant de la matrice de masse :
        ! Bottom face :
        call get_iminimax(Elem,0,imin,imax)
        matD(imin:imax,0,0) = matD(imin:imax,0,0) + Elem%Acoeff(0:ngx-1,0,12)*Elem%Density(0:ngx-1,0)
        matD(imin:imax,1,1) = matD(imin:imax,1,1) + Elem%Acoeff(0:ngx-1,0,12)*Elem%Density(0:ngx-1,0)
        ! Right face :
        call get_iminimax(Elem,1,imin,imax)
        matD(imin:imax,0,0) = matD(imin:imax,0,0) + Elem%Acoeff(ngx-1,0:ngz-1,12)*Elem%Density(ngx-1,0:ngz-1)
        matD(imin:imax,1,1) = matD(imin:imax,1,1) + Elem%Acoeff(ngx-1,0:ngz-1,12)*Elem%Density(ngx-1,0:ngz-1)
        ! Top Face :
        call get_iminimax(Elem,2,imin,imax)
        matD(imin:imax,0,0) = matD(imin:imax,0,0) + Elem%Acoeff(0:ngx-1,ngz-1,12)*Elem%Density(0:ngx-1,ngz-1)
        matD(imin:imax,1,1) = matD(imin:imax,1,1) + Elem%Acoeff(0:ngx-1,ngz-1,12)*Elem%Density(0:ngx-1,ngz-1)
        ! Left Face :
        call get_iminimax(Elem,3,imin,imax)
        matD(imin:imax,0,0) = matD(imin:imax,0,0) + Elem%Acoeff(0,0:ngz-1,12)*Elem%Density(0,0:ngz-1)
        matD(imin:imax,1,1) = matD(imin:imax,1,1) + Elem%Acoeff(0,0:ngz-1,12)*Elem%Density(0,0:ngz-1)

        ! Inversion de la matrice D sur tous les noeuds de bord :
        det(:) = matD(:,0,0) * matD(:,1,1) - matD(:,0,1) * matD(:,1,0)
        tmp(:) = matD(:,0,0)
        matD(:,0,0) = 1./det(:)  * matD(:,1,1)
        matD(:,1,1) = 1./det(:)  * tmp(:)
        matD(:,0,1) = -1./det(:) * matD(:,0,1)
        matD(:,1,0) = -1./det(:) * matD(:,1,0)

        ! Stockage de la matrice D^-1 :
        Elem%Dinv(:,0) = matD(:,0,0)
        Elem%Dinv(:,1) = matD(:,1,1)
        Elem%Dinv(:,2) = matD(:,0,1)

        ! Calcul du produit matriciel E * D^-1 :
        Elem%EDinv(:,0,0) = Elem%Coeff_Integr_Faces(:) * (Elem%MatPen(:,0)*matD(:,0,0)+Elem%MatPen(:,2)*matD(:,1,0))
        Elem%EDinv(:,0,1) = Elem%Coeff_Integr_Faces(:) * (Elem%MatPen(:,0)*matD(:,0,1)+Elem%MatPen(:,2)*matD(:,1,1))
        Elem%EDinv(:,1,0) = Elem%Coeff_Integr_Faces(:) * (Elem%MatPen(:,2)*matD(:,0,0)+Elem%MatPen(:,1)*matD(:,1,0))
        Elem%EDinv(:,1,1) = Elem%Coeff_Integr_Faces(:) * (Elem%MatPen(:,2)*matD(:,0,1)+Elem%MatPen(:,1)*matD(:,1,1))

    end subroutine compute_EDinv

    ! ###########################################################
    !>
    !! \brief This subroutine computes the local matrix E.J which
    !! will be used to build the system on Lagrange multiplicators K * Lambda = R
    !! Indeed E.J is usefull for computing both K and R.
    !! Note that J = D^-1 . C . A^-1 and therefore, E.J = EDinv.CAinv.
    !! Moreover, E.J will stored such as CAinv = CAinv + E;J.
    !! This matrix is defined for each node lying on the element border.
    !! It suitable for Hybridizable Discontinuous Galerkin elements only,
    !! with a semi-implicit time scheme strategy.
    !! \param type (Element), intent (INOUT) Elem
    !<
    subroutine compute_EJ (Elem)

        type (Element), intent (INOUT)   :: Elem
        integer :: i, ngx, ngz, nc, n1, n2
        real, dimension(0:2*(Elem%ngllx+Elem%ngllz)-1,0:1,0:2) :: CAinv_tmp
        real, dimension(0:1,0:2) :: EJ

        ngx = Elem%ngllx ; ngz = Elem%ngllz

        ! Calcul de la matrice  CAinv_tmp :
        ! Attention la matrice CAinv_tmp ne devrait pas etre dedoublee aux noeuds des coin
        ! Ici cependant elle est dedoublee car il faudra la multiplier avec E qui, elle,
        ! est en effet dedoublee. Du coup les 2 entres de CAinv_tmp pour un coin ont la meme valeur,
        ! qui est la somme des contributions des deux normales au coin.

        ! Termes provenant du terme tau=matpen
        CAinv_tmp = Elem%CAinv
        ! Couplage aux coins :
        do nc=0,3
            call get_gll_arround_corner(Elem,nc,n1,n2)
            CAinv_tmp(n1,:,:) = CAinv_tmp(n1,:,:) + CAinv_tmp(n2,:,:)
            CAinv_tmp(n2,:,:) = CAinv_tmp(n1,:,:)
        enddo

        do i=0,(2*(ngx+ngz)-1)
            EJ(:,:) = MATMUL(Elem%EDinv(i,:,:),CAinv_tmp(i,:,:))
            Elem%CAinv(i,:,:) = Elem%CAinv(i,:,:) - EJ(:,:)
        enddo

    end subroutine compute_EJ

    ! ###########################################################
    !>
    !! \brief This subroutine computes the matrix K on a given face
    !! (for the system on Lagrange multiplicators K * Lambda = R)
    !! on the interiors nodes of the face only.
    !! This matrix is defined for each node lying on the element border.
    !! It suitable for Hybridizable Discontinuous Galerkin elements only.
    !! \param type (domain), intent (INOUT) Tdomain
    !<
    subroutine build_K_on_face(Tdomain, nelem, Dt)

        implicit none
        type (domain), intent (INOUT) :: Tdomain
        integer, intent(IN) :: nelem
        real,    intent(IN) :: Dt
        real, dimension(0:2*(Tdomain%specel(nelem)%ngllx+Tdomain%specel(nelem)%ngllz)-1,0:1,0:1) :: CtAC, EtDE
        real, dimension(0:2*(Tdomain%specel(nelem)%ngllx+Tdomain%specel(nelem)%ngllz)-1,0:2) :: K
        type(element), pointer :: Elem
        integer :: nf, nface, i, imin, imax, n1, n2, nv, pos1, pos2
        logical :: coherency

        Elem => Tdomain%specel(nelem)

        ! Calcul du terme Ct*Ainv*C qui contribue a la matrice K
        CtAC(:,0,0) = Elem%Coeff_integr_Faces(:) * (Elem%CAinv(:,0,0)*Elem%Normal_Nodes(:,0) &
                                                   +Elem%CAinv(:,0,2)*Elem%Normal_Nodes(:,1))
        CtAC(:,0,1) = Elem%Coeff_integr_Faces(:) * (Elem%CAinv(:,0,1)*Elem%Normal_Nodes(:,1) &
                                                   +Elem%CAinv(:,0,2)*Elem%Normal_Nodes(:,0))
        CtAC(:,1,0) = Elem%Coeff_integr_Faces(:) * (Elem%CAinv(:,1,0)*Elem%Normal_Nodes(:,0) &
                                                   +Elem%CAinv(:,1,2)*Elem%Normal_Nodes(:,1))
        CtAC(:,1,1) = Elem%Coeff_integr_Faces(:) * (Elem%CAinv(:,1,1)*Elem%Normal_Nodes(:,1) &
                                                   +Elem%CAinv(:,1,2)*Elem%Normal_Nodes(:,0))

        ! Calcul du terme Et*Dinv*E qui contribue a la matrice K
        EtDE(:,0,0) = Elem%Coeff_integr_Faces(:) * (Elem%EDinv(:,0,0)*Elem%MatPen(:,0) &
                                                   +Elem%EDinv(:,0,1)*Elem%MatPen(:,2))
        EtDE(:,0,1) = Elem%Coeff_integr_Faces(:) * (Elem%EDinv(:,0,0)*Elem%MatPen(:,2) &
                                                   +Elem%EDinv(:,0,1)*Elem%MatPen(:,1))
        EtDE(:,1,0) = Elem%Coeff_integr_Faces(:) * (Elem%EDinv(:,1,0)*Elem%MatPen(:,0) &
                                                   +Elem%EDinv(:,1,1)*Elem%MatPen(:,2))
        EtDE(:,1,1) = Elem%Coeff_integr_Faces(:) * (Elem%EDinv(:,1,0)*Elem%MatPen(:,2) &
                                                   +Elem%EDinv(:,1,1)*Elem%MatPen(:,1))

        ! Addition de toutes les contributions pour la matrice K : (etant sym, elle n'a que 3 composantes)
        K(:,0) = 0.5*Dt * (CtAC(:,0,0) - EtDE(:,0,0)) + Elem%Coeff_integr_Faces(:) * Elem%MatPen(:,0)
        K(:,1) = 0.5*Dt * (CtAC(:,1,1) - EtDE(:,1,1)) + Elem%Coeff_integr_Faces(:) * Elem%MatPen(:,1)
        K(:,2) = 0.5*Dt * (CtAC(:,0,1) - EtDE(:,0,1)) + Elem%Coeff_integr_Faces(:) * Elem%MatPen(:,2)

        ! Envoi des matrices sur les faces :
        do nf=0,3
            nface = Elem%Near_Face(nf)
            call get_iminimax(Elem,nf,imin,imax)
            coherency  = Tdomain%sFace(nface)%coherency
            if (coherency .OR. (Tdomain%sFace(nface)%Near_Element(0)==nelem)) then
                Tdomain%sFace(nface)%Kinv(:,:) = Tdomain%sFace(nface)%Kinv(:,:) + K(imin:imax,:)
            else
                do i=0,Tdomain%sFace(nface)%ngll-1
                    Tdomain%sFace(nface)%Kinv(i,:) = Tdomain%sFace(nface)%Kinv(i,:) + K(imax-i,:)
                end do
            endif
        enddo

        ! Termes diagonaux seulement intervenant dans les systemes aux vertexs :
        do i=0,3 ! ieme coin
            call get_gll_arround_corner(Elem,i,n1,n2)
            nv = Elem%Near_Vertex(i)
            ! Position dans les matrices des vertexs :
            pos1 = Elem%pos_corner_in_VertMat(i,0)
            pos2 = Elem%pos_corner_in_VertMat(i,1)
            ! Termes diagonaux des matrices sur les vertexs :
            Tdomain%sVertex(nv)%Kmat(pos1,pos1)     = Tdomain%sVertex(nv)%Kmat(pos1,pos1)    + K(n1,0)
            Tdomain%sVertex(nv)%Kmat(pos1,pos1+1)   = Tdomain%sVertex(nv)%Kmat(pos1,pos1+1)  + K(n1,2)
            Tdomain%sVertex(nv)%Kmat(pos1+1,pos1)   = Tdomain%sVertex(nv)%Kmat(pos1+1,pos1)  + K(n1,2)
            Tdomain%sVertex(nv)%Kmat(pos1+1,pos1+1) = Tdomain%sVertex(nv)%Kmat(pos1+1,pos1+1)+ K(n1,1)
            Tdomain%sVertex(nv)%Kmat(pos2,pos2)     = Tdomain%sVertex(nv)%Kmat(pos2,pos2)    + K(n2,0)
            Tdomain%sVertex(nv)%Kmat(pos2,pos2+1)   = Tdomain%sVertex(nv)%Kmat(pos2,pos2+1)  + K(n2,2)
            Tdomain%sVertex(nv)%Kmat(pos2+1,pos2)   = Tdomain%sVertex(nv)%Kmat(pos2+1,pos2)  + K(n2,2)
            Tdomain%sVertex(nv)%Kmat(pos2+1,pos2+1) = Tdomain%sVertex(nv)%Kmat(pos2+1,pos2+1)+ K(n2,1)
            !Tdomain%sVertex(n1)%Kmat(pos1:pos1+1,pos1:pos1+1) = K(imin,0:1,0:1)
            !Tdomain%sVertex(n2)%Kmat(pos2:pos2+1,pos2:pos2+1) = K(imax,0:1,0:1)
        enddo

    end subroutine build_K_on_face


    ! ###########################################################
    !>
    !! \brief This subroutine sends the contributions of a given element
    !! Elem to the K-matrices of its neighbouring vertices. It actually
    !! computes the non-diagonal terms only. Diagonal terms have been already
    !! dealt with subroutine build_K_on_face.
    !! (K refers to the system on Lagrange multiplicators K * Lambda = R)
    !! It suitable for Hybridizable Discontinuous Galerkin elements only,
    !! and only for a semi-implicit time scheme.
    !! \param type (domain), intent (INOUT) Tdomain
    !<
    subroutine build_K_on_vertex(Tdomain, nelem, Dt)

        implicit none
        type (domain), intent (INOUT) :: Tdomain
        integer, intent(IN) :: nelem
        real,    intent(IN) :: Dt
        real, dimension(0:2,0:1) :: C1, C2
        real, dimension(0:1,0:1) :: E1, E2, K12, K21
        type(element), pointer :: Elem
        integer :: i, nv, n1, n2, pos1, pos2

        Elem => Tdomain%specel(nelem)

        ! Termes extra-diagonaux correspondant aux coins :
        C1 = 0. ; C2 = 0. ; E1 = 0. ; E2 = 0.
        do i=0,3 ! i_eme coin
            nv = Elem%Near_Vertex(i)
            call get_gll_arround_corner(Elem,i,n1,n2)
            ! Computation of C matrices at the 2 sides of the corner
            C1(0,0) = Elem%Coeff_integr_Faces(n1)*Elem%normal_nodes(n1,0)
            C1(1,1) = Elem%Coeff_integr_Faces(n1)*Elem%normal_nodes(n1,1)
            C1(2,0) = Elem%Coeff_integr_Faces(n1)*Elem%normal_nodes(n1,1)
            C1(2,1) = Elem%Coeff_integr_Faces(n1)*Elem%normal_nodes(n1,0)
            C2(0,0) = Elem%Coeff_integr_Faces(n2)*Elem%normal_nodes(n2,0)
            C2(1,1) = Elem%Coeff_integr_Faces(n2)*Elem%normal_nodes(n2,1)
            C2(2,0) = Elem%Coeff_integr_Faces(n2)*Elem%normal_nodes(n2,1)
            C2(2,1) = Elem%Coeff_integr_Faces(n2)*Elem%normal_nodes(n2,0)
            ! Computation of C matrices at the 2 sides of the corner
            E1(0,0) = Elem%Coeff_integr_Faces(n1)*Elem%MatPen(n1,0)
            E1(0,1) = Elem%Coeff_integr_Faces(n1)*Elem%MatPen(n1,2)
            E1(1,0) = Elem%Coeff_integr_Faces(n1)*Elem%MatPen(n1,2)
            E1(1,1) = Elem%Coeff_integr_Faces(n1)*Elem%MatPen(n1,1)
            E2(0,0) = Elem%Coeff_integr_Faces(n2)*Elem%MatPen(n2,0)
            E2(0,1) = Elem%Coeff_integr_Faces(n2)*Elem%MatPen(n2,2)
            E2(1,0) = Elem%Coeff_integr_Faces(n2)*Elem%MatPen(n2,2)
            E2(1,1) = Elem%Coeff_integr_Faces(n2)*Elem%MatPen(n2,1)
            ! Computation of the 2 cross-terms K12 and K21
            K12(:,:) = MATMUL(Elem%CAinv(n1,:,:),C2(:,:))
            K21(:,:) = MATMUL(Elem%CAinv(n2,:,:),C1(:,:))
            K12(:,:) = K12(:,:) - MATMUL(Elem%EDinv(n1,:,:),E2(:,:))
            K21(:,:) = K21(:,:) - MATMUL(Elem%EDinv(n2,:,:),E1(:,:))
            ! Envoi des matrices K12 et K21 sur la matrice du vertex
            pos1 = Elem%pos_corner_in_VertMat(i,0) ; pos2 = Elem%pos_corner_in_VertMat(i,1)
            Tdomain%sVertex(nv)%Kmat(pos1:pos1+1,pos2:pos2+1) = 0.5*Dt * K12(:,:)
            Tdomain%sVertex(nv)%Kmat(pos2:pos2+1,pos1:pos1+1) = 0.5*Dt * K21(:,:)
        enddo

    end subroutine build_K_on_vertex

    ! ############################################################


end module scompute_coeff_HDG
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
