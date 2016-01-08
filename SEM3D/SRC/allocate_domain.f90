!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file allocate_domain.f90
!!\brief G�re l'allocation des domaines.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Allocation des attributs de la structure domain
!!
!! \param type(domain), intent (INOUT) Tdomain
!<
module sdomain_alloc
    use sdomain
    implicit none
contains
subroutine allocate_domain (Tdomain)

    type(domain), intent (INOUT) :: Tdomain
    integer :: n,ngllx,nglly,ngllz
    integer :: n_solid
    integer :: mat, randSize, assocMat
    logical :: ispml, issolid

    do mat = 0,Tdomain%n_mat-1
        assocMat = Tdomain%sSubdomain(mat)%assocMat
        if(Tdomain%sSubDomain(assocMat)%material_type == "R") then
            if (Tdomain%subD_exist(mat)) then
                call random_seed(size = randSize)

                if(.not.(allocated(Tdomain%sSubdomain(mat)%chosenSeed))) then
                    allocate(Tdomain%sSubdomain(mat)%chosenSeed(randSize))
                end if
            end if
        end if
    end do

    do n = 0,Tdomain%n_elem-1

        n_solid = Tdomain%n_sls

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        ispml = Tdomain%specel(n)%domain==DM_SOLID_PML .or. Tdomain%specel(n)%domain==DM_FLUID_PML
        issolid = Tdomain%specel(n)%domain==DM_SOLID_PML .or. Tdomain%specel(n)%domain==DM_SOLID

        allocate(Tdomain%specel(n)%Density(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%MassMat(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%Lambda (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%Mu     (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%Kappa  (0:ngllx-1, 0:nglly-1, 0:ngllz-1))

        if(ispml)then ! PML Common parts
            allocate(Tdomain%specel(n)%xpml)
            allocate(Tdomain%specel(n)%xpml%DumpSx(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
            allocate(Tdomain%specel(n)%xpml%DumpSy(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
            allocate(Tdomain%specel(n)%xpml%DumpSz(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
            allocate(Tdomain%specel(n)%xpml%DumpMass(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
        end if

        if(issolid)then  ! SOLID PART
            allocate(Tdomain%specel(n)%sl)
            if(Tdomain%TimeD%velocity_scheme)then
                if(ispml)then
                    allocate(Tdomain%specel(n)%slpml)
                    allocate(Tdomain%specel(n)%slpml%Diagonal_Stress(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Diagonal_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Diagonal_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Diagonal_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Residual_Stress(0:ngllx-1, 0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Residual_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Residual_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Residual_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    Tdomain%specel(n)%slpml%Diagonal_Stress = 0d0
                    Tdomain%specel(n)%slpml%Diagonal_Stress1 = 0d0
                    Tdomain%specel(n)%slpml%Diagonal_Stress2 = 0d0
                    Tdomain%specel(n)%slpml%Diagonal_Stress3 = 0d0
                    Tdomain%specel(n)%slpml%Residual_Stress = 0d0
                    Tdomain%specel(n)%slpml%Residual_Stress1 = 0d0
                    Tdomain%specel(n)%slpml%Residual_Stress2 = 0d0
                    Tdomain%specel(n)%slpml%Residual_Stress3 = 0d0

                else ! PML
                    if (Tdomain%aniso) then
                        allocate (Tdomain%specel(n)%sl%Cij (0:20, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                    endif
                    if (n_solid>0) then
                        if (Tdomain%aniso) then
                            allocate (Tdomain%specel(n)%sl%Q (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        else
                            !              allocate (Tdomain%specel(n)%Kappa (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate (Tdomain%specel(n)%sl%Qs (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate (Tdomain%specel(n)%sl%Qp (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate (Tdomain%specel(n)%sl%onemPbeta (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate (Tdomain%specel(n)%sl%epsilonvol_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            Tdomain%specel(n)%sl%epsilonvol_ = 0
                            allocate (Tdomain%specel(n)%sl%factor_common_P (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate (Tdomain%specel(n)%sl%alphaval_P (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate (Tdomain%specel(n)%sl%betaval_P (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate (Tdomain%specel(n)%sl%gammaval_P (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate (Tdomain%specel(n)%sl%R_vol_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            Tdomain%specel(n)%sl%R_vol_ = 0
                        endif
                        allocate (Tdomain%specel(n)%sl%onemSbeta (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%epsilondev_xx_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%epsilondev_yy_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%epsilondev_xy_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%epsilondev_xz_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%epsilondev_yz_ (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        Tdomain%specel(n)%sl%epsilondev_xx_ = 0
                        Tdomain%specel(n)%sl%epsilondev_yy_ = 0
                        Tdomain%specel(n)%sl%epsilondev_xy_ = 0
                        Tdomain%specel(n)%sl%epsilondev_xz_ = 0
                        Tdomain%specel(n)%sl%epsilondev_yz_ = 0
                        allocate (Tdomain%specel(n)%sl%factor_common_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%alphaval_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%betaval_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%gammaval_3 (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%R_xx_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%R_yy_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%R_xy_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%R_xz_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        allocate (Tdomain%specel(n)%sl%R_yz_ (0:n_solid-1, 0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        Tdomain%specel(n)%sl%R_xx_ = 0
                        Tdomain%specel(n)%sl%R_yy_ = 0
                        Tdomain%specel(n)%sl%R_xy_ = 0
                        Tdomain%specel(n)%sl%R_xz_ = 0
                        Tdomain%specel(n)%sl%R_yz_ = 0
                    endif ! n_solid
                endif !PML
            endif
        else   ! FLUID PART
            allocate(Tdomain%specel(n)%fl)
            if(Tdomain%TimeD%velocity_scheme)then
                if(ispml)then
                    allocate(Tdomain%specel(n)%flpml)
                    allocate(Tdomain%specel(n)%flpml%Veloc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    Tdomain%specel(n)%flpml%Veloc = 0d0
                endif
            endif

        end if
    enddo

    ! Allocation et initialisation de Tdomain%champs0 et champs1 pour les solides
    if (Tdomain%sdom%ngll /= 0) then
        allocate(Tdomain%sdom%champs0%Forces(0:Tdomain%sdom%ngll-1,0:2))
        allocate(Tdomain%sdom%champs0%Depla (0:Tdomain%sdom%ngll-1,0:2))
        allocate(Tdomain%sdom%champs0%Veloc (0:Tdomain%sdom%ngll-1,0:2))
        allocate(Tdomain%sdom%champs1%Forces(0:Tdomain%sdom%ngll-1,0:2))
        allocate(Tdomain%sdom%champs1%Depla (0:Tdomain%sdom%ngll-1,0:2))
        allocate(Tdomain%sdom%champs1%Veloc (0:Tdomain%sdom%ngll-1,0:2))

        Tdomain%sdom%champs0%Forces = 0d0
        Tdomain%sdom%champs0%Depla = 0d0
        Tdomain%sdom%champs0%Veloc = 0d0

        ! Allocation de Tdomain%MassMat pour les solides
        allocate(Tdomain%sdom%MassMat(0:Tdomain%sdom%ngll-1))
        Tdomain%sdom%MassMat = 0d0
    endif

    ! Allocation et initialisation de Tdomain%champs0 pour les PML solides
    if (Tdomain%spmldom%ngll /= 0) then
        allocate(Tdomain%spmldom%champs1%ForcesPML(0:Tdomain%spmldom%ngll-1,0:2,0:2))
        allocate(Tdomain%spmldom%champs0%VelocPML (0:Tdomain%spmldom%ngll-1,0:2,0:2))
        allocate(Tdomain%spmldom%champs1%VelocPML (0:Tdomain%spmldom%ngll-1,0:2,0:2))
        allocate(Tdomain%spmldom%champs0%DumpV    (0:Tdomain%spmldom%ngll-1,0:1,0:2))
        Tdomain%spmldom%champs1%ForcesPML = 0d0
        Tdomain%spmldom%champs0%VelocPML = 0d0
        Tdomain%spmldom%champs0%DumpV = 0d0

        ! Allocation de Tdomain%MassMat pour les PML solides
        allocate(Tdomain%spmldom%MassMat(0:Tdomain%spmldom%ngll-1))
        Tdomain%spmldom%MassMat = 0d0

        allocate(Tdomain%spmldom%DumpMass(0:Tdomain%spmldom%ngll-1,0:2))
        Tdomain%spmldom%DumpMass = 0d0
    endif

    ! Allocation et initialisation de Tdomain%champs0 et champs1 pour les fluides
    if (Tdomain%fdom%ngll /= 0) then
        allocate(Tdomain%fdom%champs0%ForcesFl(0:Tdomain%fdom%ngll-1))
        allocate(Tdomain%fdom%champs0%Phi     (0:Tdomain%fdom%ngll-1))
        allocate(Tdomain%fdom%champs0%VelPhi  (0:Tdomain%fdom%ngll-1))
        allocate(Tdomain%fdom%champs1%ForcesFl(0:Tdomain%fdom%ngll-1))
        allocate(Tdomain%fdom%champs1%Phi     (0:Tdomain%fdom%ngll-1))
        allocate(Tdomain%fdom%champs1%VelPhi  (0:Tdomain%fdom%ngll-1))

        Tdomain%fdom%champs0%ForcesFl = 0d0
        Tdomain%fdom%champs0%Phi = 0d0
        Tdomain%fdom%champs0%VelPhi = 0d0

        ! Allocation de Tdomain%MassMat pour les fluides
        allocate(Tdomain%fdom%MassMat(0:Tdomain%fdom%ngll-1))
        Tdomain%fdom%MassMat = 0d0
    endif

    ! Allocation et initialisation de Tdomain%champs0 pour les PML fluides
    if (Tdomain%fpmldom%ngll /= 0) then
        allocate(Tdomain%fpmldom%champs1%fpml_Forces(0:Tdomain%fpmldom%ngll-1,0:2))
        allocate(Tdomain%fpmldom%champs0%fpml_VelPhi(0:Tdomain%fpmldom%ngll-1,0:2))
        allocate(Tdomain%fpmldom%champs0%fpml_Phi   (0:Tdomain%fpmldom%ngll-1,0:2))
        allocate(Tdomain%fpmldom%champs1%fpml_VelPhi(0:Tdomain%fpmldom%ngll-1,0:2))
        allocate(Tdomain%fpmldom%champs0%fpml_DumpV (0:Tdomain%fpmldom%ngll-1,0:1,0:2))
        Tdomain%fpmldom%champs1%fpml_Forces = 0d0
        Tdomain%fpmldom%champs0%fpml_VelPhi = 0d0
        Tdomain%fpmldom%champs0%fpml_Phi = 0d0
        Tdomain%fpmldom%champs0%fpml_DumpV = 0d0

        ! Allocation de Tdomain%MassMat pour les PML fluides
        allocate(Tdomain%fpmldom%MassMat(0:Tdomain%fpmldom%ngll-1))
        Tdomain%fpmldom%MassMat = 0d0

        allocate(Tdomain%fpmldom%DumpMass(0:Tdomain%fpmldom%ngll-1,0:2))
        Tdomain%fpmldom%DumpMass = 0d0

    endif

    return
end subroutine allocate_domain


end module sdomain_alloc

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
