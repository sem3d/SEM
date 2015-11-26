!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file allocate_domain.f90
!!\brief Gére l'allocation des domaines.
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
                    else
                        if (Tdomain%nl_flag == 1) then ! NL variables
                            allocate(Tdomain%specel(n)%sl%nl_param_el)
                            allocate(Tdomain%specel(n)%sl%nl_param_el%lmc_param_el)
                            allocate(Tdomain%specel(n)%sl%nl_param_el%lmc_param_el%sigma_yld (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate(Tdomain%specel(n)%sl%nl_param_el%lmc_param_el%b_iso     (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate(Tdomain%specel(n)%sl%nl_param_el%lmc_param_el%Rinf_iso  (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate(Tdomain%specel(n)%sl%nl_param_el%lmc_param_el%C_kin     (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                            allocate(Tdomain%specel(n)%sl%nl_param_el%lmc_param_el%kapa_kin  (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
                        end if
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
                        end if
                    end if
                end if
            else   ! FLUID PART
                allocate(Tdomain%specel(n)%fl)
                if(Tdomain%TimeD%velocity_scheme)then
                    if(ispml)then
                        allocate(Tdomain%specel(n)%flpml)
                        allocate(Tdomain%specel(n)%flpml%Veloc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                        Tdomain%specel(n)%flpml%Veloc = 0d0
                    end if
                end if
            end if
        enddo

        ! Allocation et initialisation de Tdomain%champs0 et champs1 pour les solides
        if (Tdomain%ngll_s /= 0) then
            allocate(Tdomain%champs0%Forces(0:Tdomain%ngll_s-1,0:2))
            allocate(Tdomain%champs0%Depla(0:Tdomain%ngll_s-1,0:2))
            allocate(Tdomain%champs0%Veloc(0:Tdomain%ngll_s-1,0:2))
            allocate(Tdomain%champs1%Forces(0:Tdomain%ngll_s-1,0:2))
            allocate(Tdomain%champs1%Depla(0:Tdomain%ngll_s-1,0:2))
            allocate(Tdomain%champs1%Veloc(0:Tdomain%ngll_s-1,0:2))
            Tdomain%champs0%Forces = 0d0
            Tdomain%champs0%Depla  = 0d0
            Tdomain%champs0%Veloc  = 0d0

            if (Tdomain%nl_flag == 1) then
                allocate(Tdomain%champs0%Epsilon_pl (0:Tdomain%ngll_s-1,0:5))
                allocate(Tdomain%champs0%Stress     (0:Tdomain%ngll_s-1,0:5))
                allocate(Tdomain%champs0%Xkin       (0:Tdomain%ngll_s-1,0:5))
                allocate(Tdomain%champs0%Riso       (0:Tdomain%ngll_s-1))
                allocate(Tdomain%champs1%Epsilon_pl (0:Tdomain%ngll_s-1,0:5))
                allocate(Tdomain%champs1%Stress     (0:Tdomain%ngll_s-1,0:5))
                allocate(Tdomain%champs1%Xkin       (0:Tdomain%ngll_s-1,0:5))
                allocate(Tdomain%champs1%Riso       (0:Tdomain%ngll_s-1))
                Tdomain%champs0%Epsilon_pl = 0d0
                Tdomain%champs0%Stress     = 0d0
                Tdomain%champs0%Xkin       = 0d0
                Tdomain%champs0%Riso       = 0d0
            end if

            ! Allocation de Tdomain%MassMatSol pour les solides
            allocate(Tdomain%MassMatSol(0:Tdomain%ngll_s-1))
            Tdomain%MassMatSol = 0d0
        endif

        ! Allocation et initialisation de Tdomain%champs0 pour les PML solides
        if (Tdomain%ngll_pmls /= 0) then
            allocate(Tdomain%champs1%ForcesPML(0:Tdomain%ngll_pmls-1,0:2,0:2))
            allocate(Tdomain%champs0%VelocPML(0:Tdomain%ngll_pmls-1,0:2,0:2))
            allocate(Tdomain%champs1%VelocPML(0:Tdomain%ngll_pmls-1,0:2,0:2))
            allocate(Tdomain%champs0%DumpV(0:Tdomain%ngll_pmls-1,0:1,0:2))
            Tdomain%champs1%ForcesPML = 0d0
            Tdomain%champs0%VelocPML = 0d0
            Tdomain%champs0%DumpV = 0d0

            ! Allocation de Tdomain%MassMatSolPml pour les PML solides
            allocate(Tdomain%MassMatSolPml(0:Tdomain%ngll_pmls-1))
            Tdomain%MassMatSolPml = 0d0

            allocate(Tdomain%DumpMass(0:Tdomain%ngll_pmls-1,0:2))
            Tdomain%DumpMass = 0d0
        endif

        ! Allocation et initialisation de Tdomain%champs0 et champs1 pour les fluides
        if (Tdomain%ngll_f /= 0) then
            allocate(Tdomain%champs0%ForcesFl(0:Tdomain%ngll_f-1))
            allocate(Tdomain%champs0%Phi(0:Tdomain%ngll_f-1))
            allocate(Tdomain%champs0%VelPhi(0:Tdomain%ngll_f-1))
            allocate(Tdomain%champs1%ForcesFl(0:Tdomain%ngll_f-1))
            allocate(Tdomain%champs1%Phi(0:Tdomain%ngll_f-1))
            allocate(Tdomain%champs1%VelPhi(0:Tdomain%ngll_f-1))

            Tdomain%champs0%ForcesFl = 0d0
            Tdomain%champs0%Phi = 0d0
            Tdomain%champs0%VelPhi = 0d0

            ! Allocation de Tdomain%MassMatFlu pour les fluides
            allocate(Tdomain%MassMatFlu(0:Tdomain%ngll_f-1))
            Tdomain%MassMatFlu = 0d0
        endif

        ! Allocation et initialisation de Tdomain%champs0 pour les PML fluides
        if (Tdomain%ngll_pmlf /= 0) then
            allocate(Tdomain%champs1%fpml_Forces(0:Tdomain%ngll_pmlf-1,0:2))
            allocate(Tdomain%champs0%fpml_VelPhi(0:Tdomain%ngll_pmlf-1,0:2))
            allocate(Tdomain%champs0%fpml_Phi(0:Tdomain%ngll_pmlf-1,0:2))
            allocate(Tdomain%champs1%fpml_VelPhi(0:Tdomain%ngll_pmlf-1,0:2))
            allocate(Tdomain%champs0%fpml_DumpV(0:Tdomain%ngll_pmlf-1,0:1,0:2))
            Tdomain%champs1%fpml_Forces = 0d0
            Tdomain%champs0%fpml_VelPhi = 0d0
            Tdomain%champs0%fpml_Phi = 0d0
            Tdomain%champs0%fpml_DumpV = 0d0

            ! Allocation de Tdomain%MassMatSolPml pour les PML solides
            allocate(Tdomain%MassMatFluPml(0:Tdomain%ngll_pmlf-1))
            Tdomain%MassMatFluPml = 0d0

            allocate(Tdomain%fpml_DumpMass(0:Tdomain%ngll_pmlf-1,0:2))
            Tdomain%fpml_DumpMass = 0d0

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
