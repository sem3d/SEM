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
    integer :: n,nf,ne,nv,ngllx,nglly,ngllz,ngll1,ngll2,ngll
    integer :: n_solid
    integer :: mat, randSize, assocMat

    do mat = 0,Tdomain%n_mat-1
        assocMat = Tdomain%sSubdomain(mat)%assocMat
        if(Tdomain%sSubDomain(assocMat)%material_type == "R") then
            if (Tdomain%subD_exist(mat)) then
                call random_seed(size = randSize)
                allocate(Tdomain%sSubdomain(mat)%chosenSeed(randSize))
                allocate(Tdomain%sSubDomain(mat)%MinBound(0:2))
                allocate(Tdomain%sSubDomain(mat)%MaxBound(0:2))
            end if
        end if
    end do

    do n = 0,Tdomain%n_elem-1

        n_solid = Tdomain%n_sls

        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        allocate(Tdomain%specel(n)%Density(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%MassMat(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%Lambda (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%Mu     (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%Kappa  (0:ngllx-1, 0:nglly-1, 0:ngllz-1))

        if(Tdomain%specel(n)%PML)then ! PML Common parts
            allocate(Tdomain%specel(n)%xpml)
            allocate(Tdomain%specel(n)%xpml%DumpSx(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
            allocate(Tdomain%specel(n)%xpml%DumpSy(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
            allocate(Tdomain%specel(n)%xpml%DumpSz(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1))
            allocate(Tdomain%specel(n)%xpml%DumpMass(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
            allocate(Tdomain%specel(n)%xpml%DumpVx(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1))
            allocate(Tdomain%specel(n)%xpml%DumpVy(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1))
            allocate(Tdomain%specel(n)%xpml%DumpVz(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1))
        end if

        if(Tdomain%specel(n)%solid)then  ! SOLID PART
            allocate(Tdomain%specel(n)%sl)
            if(Tdomain%TimeD%velocity_scheme)then
                allocate(Tdomain%specel(n)%sl%Displ (1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                allocate(Tdomain%specel(n)%sl%Veloc (1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                allocate(Tdomain%specel(n)%sl%Accel (1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                allocate(Tdomain%specel(n)%sl%V0    (1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                allocate(Tdomain%specel(n)%sl%Forces(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                Tdomain%specel(n)%sl%Displ = 0
                Tdomain%specel(n)%sl%Veloc = 0d0
                Tdomain%specel(n)%sl%Accel = 0d0
                Tdomain%specel(n)%sl%V0 = 0d0
                Tdomain%specel(n)%sl%Forces = 0d0
                if(Tdomain%specel(n)%PML)then
                    allocate(Tdomain%specel(n)%slpml)
                    allocate(Tdomain%specel(n)%sl%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:35))
                    allocate(Tdomain%specel(n)%slpml%Diagonal_Stress(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Diagonal_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Diagonal_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Diagonal_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Residual_Stress(0:ngllx-1, 0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Residual_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Residual_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Residual_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Veloc1(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                    allocate(Tdomain%specel(n)%slpml%Veloc2(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                    allocate(Tdomain%specel(n)%slpml%Veloc3(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                    allocate(Tdomain%specel(n)%slpml%Forces1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Forces2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%slpml%Forces3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    Tdomain%specel(n)%slpml%Diagonal_Stress = 0d0
                    Tdomain%specel(n)%slpml%Diagonal_Stress1 = 0d0
                    Tdomain%specel(n)%slpml%Diagonal_Stress2 = 0d0
                    Tdomain%specel(n)%slpml%Diagonal_Stress3 = 0d0
                    Tdomain%specel(n)%slpml%Residual_Stress = 0d0
                    Tdomain%specel(n)%slpml%Residual_Stress1 = 0d0
                    Tdomain%specel(n)%slpml%Residual_Stress2 = 0d0
                    Tdomain%specel(n)%slpml%Residual_Stress3 = 0d0
                    Tdomain%specel(n)%slpml%Veloc1 = 0d0
                    Tdomain%specel(n)%slpml%Veloc2 = 0d0
                    Tdomain%specel(n)%slpml%Veloc3 = 0d0

                    if(Tdomain%specel(n)%FPML)then
                        allocate(Tdomain%specel(n)%slpml%Isx(0:ngllx-1,0:nglly-1,0:ngllz-1))
                        allocate(Tdomain%specel(n)%slpml%Isy(0:ngllx-1,0:nglly-1,0:ngllz-1))
                        allocate(Tdomain%specel(n)%slpml%Isz(0:ngllx-1,0:nglly-1,0:ngllz-1))
                        allocate(Tdomain%specel(n)%slpml%Ivx(0:ngllx-1,0:nglly-1,0:ngllz-1))
                        allocate(Tdomain%specel(n)%slpml%Ivy(0:ngllx-1,0:nglly-1,0:ngllz-1))
                        allocate(Tdomain%specel(n)%slpml%Ivz(0:ngllx-1,0:nglly-1,0:ngllz-1))
                        allocate(Tdomain%specel(n)%slpml%I_Diagonal_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                        allocate(Tdomain%specel(n)%slpml%I_Diagonal_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                        allocate(Tdomain%specel(n)%slpml%I_Diagonal_Stress3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                        allocate(Tdomain%specel(n)%slpml%I_Residual_Stress1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                        allocate(Tdomain%specel(n)%slpml%I_Residual_Stress2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                        allocate(Tdomain%specel(n)%slpml%IVeloc1(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                        allocate(Tdomain%specel(n)%slpml%IVeloc2(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                        allocate(Tdomain%specel(n)%slpml%IVeloc3(1:ngllx-2,1:nglly-2,1:ngllz-2,0:2))
                        Tdomain%specel(n)%slpml%I_Diagonal_Stress1 = 0.
                        Tdomain%specel(n)%slpml%I_Diagonal_Stress2 = 0.
                        Tdomain%specel(n)%slpml%I_Diagonal_Stress3 = 0.
                        Tdomain%specel(n)%slpml%I_Residual_Stress1 = 0.
                        Tdomain%specel(n)%slpml%I_Residual_Stress2 = 0.
                        Tdomain%specel(n)%slpml%IVeloc1 = 0.
                        Tdomain%specel(n)%slpml%IVeloc2 = 0.
                        Tdomain%specel(n)%slpml%IVeloc3 = 0.
                    endif ! FPML


                    if (Tdomain%curve) then
                        allocate (Tdomain%specel(n)%slpml%Normales (0:2, 0:2))
                        allocate (Tdomain%specel(n)%slpml%Inv_Normales (0:2, 0:2))
                    endif
                else ! PML
                    !  modif mariotti fevrier 2007 cea
                    !          allocate (Tdomain%specel(n)%Displ (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
                    !          Tdomain%specel(n)%Displ = 0
                    !  modif mariotti fevrier 2007 cea capteur displ
                    !          allocate (Tdomain%specel(n)%Acoeff (0:ngllx-1, 0:nglly-1, 0:ngllz-1, 0:44))
                    !          allocate (Tdomain%specel(n)%Displ (1:ngllx-2, 1:nglly-2, 1:ngllz-2, 0:2))
                    !          Tdomain%specel(n)%Displ = 0
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
                allocate(Tdomain%specel(n)%fl%VelPhi(1:ngllx-2,1:nglly-2,1:ngllz-2))
                allocate(Tdomain%specel(n)%fl%AccelPhi(1:ngllx-2,1:nglly-2,1:ngllz-2))
                allocate(Tdomain%specel(n)%fl%VelPhi0(1:ngllx-2,1:nglly-2,1:ngllz-2))
                allocate(Tdomain%specel(n)%fl%ForcesFl(0:ngllx-1,0:nglly-1,0:ngllz-1))
                allocate(Tdomain%specel(n)%fl%Phi(1:ngllx-2,1:nglly-2,1:ngllz-2))
                Tdomain%specel(n)%fl%Phi = 0
                Tdomain%specel(n)%fl%VelPhi = 0d0
                Tdomain%specel(n)%fl%AccelPhi = 0d0
                Tdomain%specel(n)%fl%VelPhi0 = 0d0
                Tdomain%specel(n)%fl%ForcesFl = 0d0
                if(Tdomain%specel(n)%PML)then
                    allocate(Tdomain%specel(n)%flpml)
                    allocate(Tdomain%specel(n)%fl%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:17))
                    allocate(Tdomain%specel(n)%flpml%Veloc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%flpml%Veloc1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%flpml%Veloc2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%flpml%Veloc3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%flpml%VelPhi1(1:ngllx-2,1:nglly-2,1:ngllz-2))
                    allocate(Tdomain%specel(n)%flpml%VelPhi2(1:ngllx-2,1:nglly-2,1:ngllz-2))
                    allocate(Tdomain%specel(n)%flpml%VelPhi3(1:ngllx-2,1:nglly-2,1:ngllz-2))
                    allocate(Tdomain%specel(n)%flpml%ForcesFl1(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(Tdomain%specel(n)%flpml%ForcesFl2(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    allocate(Tdomain%specel(n)%flpml%ForcesFl3(0:ngllx-1,0:nglly-1,0:ngllz-1))
                    Tdomain%specel(n)%flpml%Veloc = 0d0
                    Tdomain%specel(n)%flpml%Veloc1 = 0d0
                    Tdomain%specel(n)%flpml%Veloc2 = 0d0
                    Tdomain%specel(n)%flpml%Veloc3 = 0d0
                    Tdomain%specel(n)%flpml%VelPhi1 = 0d0
                    Tdomain%specel(n)%flpml%VelPhi2 = 0d0
                    Tdomain%specel(n)%flpml%VelPhi3 = 0d0
                else
                    allocate(Tdomain%specel(n)%fl%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:5))
                endif
            endif

        end if
    enddo


    do n = 0, Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
        if(Tdomain%sFace(n)%solid)then    ! SOLID PART
            allocate(Tdomain%sFace(n)%MassMat(1:ngll1-2,1:ngll2-2))
            allocate(Tdomain%sFace(n)%Veloc(1:ngll1-2,1:ngll2-2,0:2))

            !  modif mariotti fevrier 2007 cea capteur displ
            allocate (Tdomain%sFace(n)%Displ (1:ngll1-2, 1:ngll2-2, 0:2))
            Tdomain%sFace(n)%Displ = 0

            allocate(Tdomain%sFace(n)%Forces(1:ngll1-2,1:ngll2-2,0:2))

#ifdef COUPLAGE
            allocate (Tdomain%sFace(n)%ForcesExt(1:ngll1-2,1:ngll2-2,0:2 ) )
            Tdomain%sFace(n)%ForcesExt = 0.
            !        allocate (Tdomain%sFace(n)%FlagMka(1:ngll1-2,1:ngll2-2,0:2 ) )
            !        Tdomain%sFace(n)%FlagMka = 0
            allocate (Tdomain%sFace(n)%tsurfsem(1:ngll1-2,1:ngll2-2 ) )
            Tdomain%sFace(n)%tsurfsem = 0.
#endif

            allocate(Tdomain%sFace(n)%Accel(1:ngll1-2,1:ngll2-2,0:2))
            allocate(Tdomain%sFace(n)%V0(1:ngll1-2,1:ngll2-2,0:2))
            Tdomain%sFace(n)%MassMat = 0d0
            Tdomain%sFace(n)%Veloc = 0d0
            Tdomain%sFace(n)%Accel = 0d0
            Tdomain%sFace(n)%V0 = 0d0
            Tdomain%sFace(n)%Forces = 0d0
            if(Tdomain%sFace(n)%PML)then
                allocate(Tdomain%sFace(n)%spml)
                allocate(Tdomain%sFace(n)%spml%Forces1(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%sFace(n)%spml%Forces2(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%sFace(n)%spml%Forces3(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%sFace(n)%spml%DumpMass(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%sFace(n)%spml%Veloc1(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%sFace(n)%spml%Veloc2(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%sFace(n)%spml%Veloc3(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%sFace(n)%spml%DumpVx(1:ngll1-2,1:ngll2-2,0:1))
                allocate(Tdomain%sFace(n)%spml%DumpVy(1:ngll1-2,1:ngll2-2,0:1))
                allocate(Tdomain%sFace(n)%spml%DumpVz(1:ngll1-2,1:ngll2-2,0:1))
                Tdomain%sFace(n)%spml%DumpMass = 0.
                Tdomain%sFace(n)%spml%Veloc1 = 0.
                Tdomain%sFace(n)%spml%Veloc2 = 0.
                Tdomain%sFace(n)%spml%Veloc3 = 0.
                if(Tdomain%sFace(n)%FPML)then
                    allocate(Tdomain%sFace(n)%spml%Ivx(1:ngll1-2,1:ngll2-1))
                    allocate(Tdomain%sFace(n)%spml%Ivy(1:ngll1-2,1:ngll2-1))
                    allocate(Tdomain%sFace(n)%spml%Ivz(1:ngll1-2,1:ngll2-1))
                    allocate(Tdomain%sFace(n)%spml%IVeloc1(1:ngll1-2,1:ngll2-1,0:2))
                    allocate(Tdomain%sFace(n)%spml%IVeloc2(1:ngll1-2,1:ngll2-1,0:2))
                    allocate(Tdomain%sFace(n)%spml%IVeloc3(1:ngll1-2,1:ngll2-1,0:2))
                    Tdomain%sFace(n)%spml%IVeloc1 = 0.
                    Tdomain%sFace(n)%spml%IVeloc2 = 0.
                    Tdomain%sFace(n)%spml%IVeloc3 = 0.
                    Tdomain%sFace(n)%spml%Ivx = 0.
                    Tdomain%sFace(n)%spml%Ivy = 0.
                    Tdomain%sFace(n)%spml%Ivz = 0.
                endif


            else
                !  modif mariotti fevrier 2007 cea capteur displ
                !       allocate (Tdomain%sFace(n)%Displ (1:ngll1-2, 1:ngll2-2, 0:2))
                !       Tdomain%sFace(n)%Displ = 0
            endif
        else   ! FLUID PART
            allocate(Tdomain%sFace(n)%MassMat(1:ngll1-2,1:ngll2-2))
            allocate(Tdomain%sFace(n)%VelPhi(1:ngll1-2,1:ngll2-2))
            allocate(Tdomain%sFace(n)%ForcesFl(1:ngll1-2,1:ngll2-2))
            allocate(Tdomain%sFace(n)%AccelPhi(1:ngll1-2,1:ngll2-2))
            allocate(Tdomain%sFace(n)%VelPhi0(1:ngll1-2,1:ngll2-2))
            allocate(Tdomain%sFace(n)%Phi(1:ngll1-2,1:ngll2-2))
            Tdomain%sFace(n)%Phi = 0d0
            Tdomain%sFace(n)%MassMat = 0d0
            Tdomain%sFace(n)%VelPhi = 0d0
            Tdomain%sFace(n)%AccelPhi = 0d0
            Tdomain%sFace(n)%VelPhi0 = 0d0
            Tdomain%sFace(n)%ForcesFl = 0d0
            if(Tdomain%sFace(n)%PML)then
                allocate(Tdomain%sFace(n)%spml)
                allocate(Tdomain%sFace(n)%spml%ForcesFl1(1:ngll1-2,1:ngll2-2))
                allocate(Tdomain%sFace(n)%spml%ForcesFl2(1:ngll1-2,1:ngll2-2))
                allocate(Tdomain%sFace(n)%spml%ForcesFl3(1:ngll1-2,1:ngll2-2))
                allocate(Tdomain%sFace(n)%spml%DumpMass(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%sFace(n)%spml%VelPhi1(1:ngll1-2,1:ngll2-2))
                allocate(Tdomain%sFace(n)%spml%VelPhi2(1:ngll1-2,1:ngll2-2))
                allocate(Tdomain%sFace(n)%spml%VelPhi3(1:ngll1-2,1:ngll2-2))
                allocate(Tdomain%sFace(n)%spml%DumpVx(1:ngll1-2,1:ngll2-2,0:1))
                allocate(Tdomain%sFace(n)%spml%DumpVy(1:ngll1-2,1:ngll2-2,0:1))
                allocate(Tdomain%sFace(n)%spml%DumpVz(1:ngll1-2,1:ngll2-2,0:1))
                Tdomain%sFace(n)%spml%DumpMass = 0.
                Tdomain%sFace(n)%spml%VelPhi1 = 0.
                Tdomain%sFace(n)%spml%VelPhi2 = 0.
                Tdomain%sFace(n)%spml%VelPhi3 = 0.
            endif
        endif
    enddo



    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll
        if(Tdomain%sEdge(n)%solid)then   ! SOLID PART
            allocate(Tdomain%sEdge(n)%MassMat(1:ngll-2))
            allocate(Tdomain%sEdge(n)%Veloc(1:ngll-2,0:2))

            !  modif mariotti fevrier 2007 cea capteur displ
            allocate (Tdomain%sEdge(n)%Displ (1:ngll-2, 0:2))
            Tdomain%sEdge(n)%Displ = 0

            allocate(Tdomain%sEdge(n)%Forces(1:ngll-2,0:2))
            allocate(Tdomain%sEdge(n)%Accel(1:ngll-2,0:2))
            allocate(Tdomain%sEdge(n)%V0(1:ngll-2,0:2))
            Tdomain%sEdge(n)%MassMat = 0d0
            Tdomain%sEdge(n)%Veloc = 0d0
            Tdomain%sEdge(n)%Accel = 0d0
            Tdomain%sEdge(n)%V0 = 0d0
            Tdomain%sEdge(n)%Forces = 0d0
            if(Tdomain%sEdge(n)%PML)then
                allocate(Tdomain%sEdge(n)%spml)
                allocate(Tdomain%sEdge(n)%spml%Forces1(1:ngll-2,0:2))
                allocate(Tdomain%sEdge(n)%spml%Forces2(1:ngll-2,0:2))
                allocate(Tdomain%sEdge(n)%spml%Forces3(1:ngll-2,0:2))
                allocate(Tdomain%sEdge(n)%spml%DumpMass(1:ngll-2, 0:2))
                allocate(Tdomain%sEdge(n)%spml%Veloc1(1:ngll-2,0:2))
                allocate(Tdomain%sEdge(n)%spml%Veloc2(1:ngll-2,0:2))
                allocate(Tdomain%sEdge(n)%spml%Veloc3(1:ngll-2,0:2))
                allocate(Tdomain%sEdge(n)%spml%DumpVx(1:ngll-2,0:1))
                allocate(Tdomain%sEdge(n)%spml%DumpVy(1:ngll-2,0:1))
                allocate(Tdomain%sEdge(n)%spml%DumpVz(1:ngll-2,0:1))
                Tdomain%sEdge(n)%spml%DumpMass = 0d0
                Tdomain%sEdge(n)%spml%Veloc1 = 0d0
                Tdomain%sEdge(n)%spml%Veloc2 = 0d0
                Tdomain%sEdge(n)%spml%Veloc3 = 0d0
                if(Tdomain%sEdge(n)%FPML)then
                    allocate(Tdomain%sEdge(n)%spml%Ivx(1:ngll-2))
                    allocate(Tdomain%sEdge(n)%spml%Ivy(1:ngll-2))
                    allocate(Tdomain%sEdge(n)%spml%Ivz(1:ngll-2))
                    allocate(Tdomain%sEdge(n)%spml%IVeloc1(1:ngll-2,0:2))
                    allocate(Tdomain%sEdge(n)%spml%IVeloc2(1:ngll-2,0:2))
                    allocate(Tdomain%sEdge(n)%spml%IVeloc3(1:ngll-2,0:2))
                    Tdomain%sEdge(n)%spml%IVeloc1 = 0.
                    Tdomain%sEdge(n)%spml%IVeloc2 = 0.
                    Tdomain%sEdge(n)%spml%IVeloc3 = 0.
                    Tdomain%sEdge(n)%spml%Ivx = 0.
                    Tdomain%sEdge(n)%spml%Ivy = 0.
                    Tdomain%sEdge(n)%spml%Ivz = 0.
                endif



            else
                !        allocate (Tdomain%sEdge(n)%Displ (1:ngll-2, 0:2))
                !        Tdomain%sEdge(n)%Displ = 0
            endif
        else   ! FLUID PART
            allocate(Tdomain%sEdge(n)%MassMat(1:ngll-2))
            allocate(Tdomain%sEdge(n)%VelPhi(1:ngll-2))
            allocate(Tdomain%sEdge(n)%ForcesFl(1:ngll-2))
            allocate(Tdomain%sEdge(n)%AccelPhi(1:ngll-2))
            allocate(Tdomain%sEdge(n)%VelPhi0(1:ngll-2))
            allocate(Tdomain%sEdge(n)%Phi(1:ngll-2))
            Tdomain%sEdge(n)%Phi = 0d0
            Tdomain%sEdge(n)%MassMat = 0d0
            Tdomain%sEdge(n)%VelPhi = 0d0
            Tdomain%sEdge(n)%AccelPhi = 0d0
            Tdomain%sEdge(n)%VelPhi0 = 0d0
            Tdomain%sEdge(n)%ForcesFl = 0d0
            if(Tdomain%sEdge(n)%PML)then
                allocate(Tdomain%sEdge(n)%spml)
                allocate(Tdomain%sEdge(n)%spml%ForcesFl1(1:ngll-2))
                allocate(Tdomain%sEdge(n)%spml%ForcesFl2(1:ngll-2))
                allocate(Tdomain%sEdge(n)%spml%ForcesFl3(1:ngll-2))
                allocate(Tdomain%sEdge(n)%spml%DumpMass(1:ngll-2,0:2))
                allocate(Tdomain%sEdge(n)%spml%VelPhi1(1:ngll-2))
                allocate(Tdomain%sEdge(n)%spml%VelPhi2(1:ngll-2))
                allocate(Tdomain%sEdge(n)%spml%VelPhi3(1:ngll-2))
                allocate(Tdomain%sEdge(n)%spml%DumpVx(1:ngll-2,0:1))
                allocate(Tdomain%sEdge(n)%spml%DumpVy(1:ngll-2,0:1))
                allocate(Tdomain%sEdge(n)%spml%DumpVz(1:ngll-2,0:1))
                Tdomain%sEdge(n)%spml%DumpMass = 0d0
                Tdomain%sEdge(n)%spml%VelPhi1 = 0d0
                Tdomain%sEdge(n)%spml%VelPhi2 = 0d0
                Tdomain%sEdge(n)%spml%VelPhi3 = 0d0
            endif

        endif

        !Gsa Ipsis 3D
#ifdef COUPLAGE

        allocate (Tdomain%sEdge(n)%ForcesExt(1:ngll-2,0:2) )
        Tdomain%sEdge(n)%ForcesExt = 0
        !        allocate (Tdomain%sEdge(n)%FlagMka(1:ngll-2,0:2) )
        !        Tdomain%sEdge(n)%FlagMka = 0
        allocate (Tdomain%sEdge(n)%tsurfsem(1:ngll-2) )
        Tdomain%sEdge(n)%tsurfsem = 0.
#endif

    enddo


    do n = 0,Tdomain%n_vertex-1
        if(Tdomain%sVertex(n)%solid)then  ! SOLID PART
            Tdomain%sVertex(n)%MassMat = 0d0
            Tdomain%sVertex(n)%Veloc = 0d0
            Tdomain%sVertex(n)%Accel = 0d0
            Tdomain%sVertex(n)%V0 = 0d0
            Tdomain%sVertex(n)%Forces = 0d0
#ifdef COUPLAGE
            allocate (Tdomain%sVertex(n)%ForcesExt(0:2) )
            Tdomain%sVertex(n)%ForcesExt = 0
            Tdomain%sVertex(n)%tsurfsem = 0.
#endif

            if(Tdomain%sVertex(n)%PML)then
                allocate(Tdomain%sVertex(n)%spml)
                allocate(Tdomain%sVertex(n)%spml%Veloc1(0:2))
                allocate(Tdomain%sVertex(n)%spml%Veloc2(0:2))
                allocate(Tdomain%sVertex(n)%spml%Veloc3(0:2))
                allocate(Tdomain%sVertex(n)%spml%DumpVx(0:2))
                allocate(Tdomain%sVertex(n)%spml%DumpVy(0:2))
                allocate(Tdomain%sVertex(n)%spml%DumpVz(0:2))
                allocate(Tdomain%sVertex(n)%spml%DumpMass(0:2))
                Tdomain%sVertex(n)%spml%Veloc1 = 0d0
                Tdomain%sVertex(n)%spml%Veloc2 = 0d0
                Tdomain%sVertex(n)%spml%Veloc3 = 0d0
                Tdomain%sVertex(n)%spml%DumpMass = 0d0
                if(Tdomain%sVertex(n)%FPML)then
                    allocate(Tdomain%sVertex(n)%spml%IVeloc1(0:2))
                    allocate(Tdomain%sVertex(n)%spml%IVeloc2(0:2))
                    allocate(Tdomain%sVertex(n)%spml%IVeloc3(0:2))
                    allocate(Tdomain%sVertex(n)%spml%Ivx(0:0))
                    allocate(Tdomain%sVertex(n)%spml%Ivy(0:0))
                    allocate(Tdomain%sVertex(n)%spml%Ivz(0:0))
                    Tdomain%sVertex(n)%spml%IVeloc1 = 0d0
                    Tdomain%sVertex(n)%spml%IVeloc2 = 0d0
                    Tdomain%sVertex(n)%spml%IVeloc3 = 0d0
                    Tdomain%sVertex(n)%spml%Ivx = 0d0
                    Tdomain%sVertex(n)%spml%Ivy = 0d0
                    Tdomain%sVertex(n)%spml%Ivz = 0d0
                endif
            else
                Tdomain%sVertex(n)%Displ = 0d0
            endif
        else   ! FLUID PART
            Tdomain%sVertex(n)%MassMat = 0d0
            Tdomain%sVertex(n)%VelPhi = 0d0
            Tdomain%sVertex(n)%AccelPhi = 0d0
            Tdomain%sVertex(n)%VelPhi0 = 0d0
            Tdomain%sVertex(n)%ForcesFl = 0d0
            Tdomain%sVertex(n)%Phi = 0d0
            if(Tdomain%sVertex(n)%PML)then
                allocate(Tdomain%sVertex(n)%spml)
                allocate(Tdomain%sVertex(n)%spml%DumpVx(0:2))
                allocate(Tdomain%sVertex(n)%spml%DumpVy(0:2))
                allocate(Tdomain%sVertex(n)%spml%DumpVz(0:2))
                allocate(Tdomain%sVertex(n)%spml%DumpMass(0:2))
                Tdomain%sVertex(n)%spml%VelPhi1 = 0d0
                Tdomain%sVertex(n)%spml%VelPhi2 = 0d0
                Tdomain%sVertex(n)%spml%VelPhi3 = 0d0
                Tdomain%sVertex(n)%spml%DumpMass = 0d0
            endif
        endif
    enddo

    ! solid-fluid related terms
    if(Tdomain%logicD%SF_local_present)then
        do nf = 0,Tdomain%SF%SF_n_faces-1
            ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
            ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
            allocate(Tdomain%SF%SF_Face(nf)%BtN(0:ngll1-1,0:ngll2-1,0:2))
            allocate(Tdomain%SF%SF_Face(nf)%density(0:ngll1-1,0:ngll2-1))
            Tdomain%SF%SF_Face(nf)%BtN = 0d0
            Tdomain%SF%SF_Face(nf)%density = 0d0
            allocate(Tdomain%SF%SF_Face(nf)%pn(0:ngll1-1,0:ngll2-1,0:2))
            allocate(Tdomain%SF%SF_Face(nf)%vn(0:ngll1-1,0:ngll2-1))
            allocate(Tdomain%SF%SF_Face(nf)%save_forces(1:ngll1-2,1:ngll2-2,0:2))
            Tdomain%SF%SF_Face(nf)%pn = 0d0
            Tdomain%SF%SF_Face(nf)%Vn = 0d0
            Tdomain%SF%SF_Face(nf)%save_forces = 0d0
            if(Tdomain%SF%SF_Face(nf)%PML)then
                allocate(Tdomain%SF%SF_Face(nf)%pn1(0:ngll1-1,0:ngll2-1,0:2))
                allocate(Tdomain%SF%SF_Face(nf)%pn2(0:ngll1-1,0:ngll2-1,0:2))
                allocate(Tdomain%SF%SF_Face(nf)%pn3(0:ngll1-1,0:ngll2-1,0:2))
                allocate(Tdomain%SF%SF_Face(nf)%vn1(0:ngll1-1,0:ngll2-1))
                allocate(Tdomain%SF%SF_Face(nf)%vn2(0:ngll1-1,0:ngll2-1))
                allocate(Tdomain%SF%SF_Face(nf)%vn3(0:ngll1-1,0:ngll2-1))
                allocate(Tdomain%SF%SF_Face(nf)%save_veloc1(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%SF%SF_Face(nf)%save_veloc2(1:ngll1-2,1:ngll2-2,0:2))
                allocate(Tdomain%SF%SF_Face(nf)%save_veloc3(1:ngll1-2,1:ngll2-2,0:2))
                Tdomain%SF%SF_Face(nf)%pn1 = 0d0
                Tdomain%SF%SF_Face(nf)%pn2 = 0d0
                Tdomain%SF%SF_Face(nf)%pn3 = 0d0
                Tdomain%SF%SF_Face(nf)%Vn1 = 0d0
                Tdomain%SF%SF_Face(nf)%Vn2 = 0d0
                Tdomain%SF%SF_Face(nf)%Vn3 = 0d0
                Tdomain%SF%SF_Face(nf)%save_veloc1 = 0d0
                Tdomain%SF%SF_Face(nf)%save_veloc2 = 0d0
                Tdomain%SF%SF_Face(nf)%save_veloc3 = 0d0
            else
                allocate(Tdomain%SF%SF_Face(nf)%save_displ(1:ngll1-2,1:ngll2-2,0:2))
                Tdomain%SF%SF_Face(nf)%save_displ = 0d0
            end if
        end do

        do ne = 0,Tdomain%SF%SF_n_edges-1
            ngll = Tdomain%SF%SF_Edge(ne)%ngll
            allocate(Tdomain%SF%SF_Edge(ne)%pn(1:ngll-2,0:2))
            allocate(Tdomain%SF%SF_Edge(ne)%save_forces(1:ngll-2,0:2))
            allocate(Tdomain%SF%SF_Edge(ne)%Vn(1:ngll-2))
            Tdomain%SF%SF_Edge(ne)%pn = 0d0
            Tdomain%SF%SF_Edge(ne)%save_forces = 0d0
            Tdomain%SF%SF_Edge(ne)%Vn = 0d0
            if(Tdomain%SF%SF_Edge(ne)%PML)then
                allocate(Tdomain%SF%SF_Edge(ne)%pn1(1:ngll-2,0:2))
                allocate(Tdomain%SF%SF_Edge(ne)%pn2(1:ngll-2,0:2))
                allocate(Tdomain%SF%SF_Edge(ne)%pn3(1:ngll-2,0:2))
                allocate(Tdomain%SF%SF_Edge(ne)%save_veloc1(1:ngll-2,0:2))
                allocate(Tdomain%SF%SF_Edge(ne)%save_veloc2(1:ngll-2,0:2))
                allocate(Tdomain%SF%SF_Edge(ne)%save_veloc3(1:ngll-2,0:2))
                allocate(Tdomain%SF%SF_Edge(ne)%Vn1(1:ngll-2))
                allocate(Tdomain%SF%SF_Edge(ne)%Vn2(1:ngll-2))
                allocate(Tdomain%SF%SF_Edge(ne)%Vn3(1:ngll-2))
                Tdomain%SF%SF_Edge(ne)%pn1 = 0d0
                Tdomain%SF%SF_Edge(ne)%pn2 = 0d0
                Tdomain%SF%SF_Edge(ne)%pn3 = 0d0
                Tdomain%SF%SF_Edge(ne)%save_veloc1 = 0d0
                Tdomain%SF%SF_Edge(ne)%save_veloc2 = 0d0
                Tdomain%SF%SF_Edge(ne)%save_veloc3 = 0d0
                Tdomain%SF%SF_Edge(ne)%Vn1 = 0d0
                Tdomain%SF%SF_Edge(ne)%Vn2 = 0d0
                Tdomain%SF%SF_Edge(ne)%Vn3 = 0d0
            else
                allocate(Tdomain%SF%SF_Edge(ne)%save_displ(1:ngll-2,0:2))
                Tdomain%SF%SF_Edge(ne)%save_displ = 0d0
            end if
        enddo

        do nv = 0,Tdomain%SF%SF_n_vertices-1
            Tdomain%SF%SF_Vertex(nv)%pn = 0d0
            Tdomain%SF%SF_Vertex(nv)%Vn = 0d0
            Tdomain%SF%SF_Vertex(nv)%save_forces = 0d0
            if(Tdomain%SF%SF_Vertex(nv)%PML)then
                Tdomain%SF%SF_Vertex(nv)%pn1 = 0d0
                Tdomain%SF%SF_Vertex(nv)%pn2 = 0d0
                Tdomain%SF%SF_Vertex(nv)%pn3 = 0d0
                Tdomain%SF%SF_Vertex(nv)%Vn1 = 0d0
                Tdomain%SF%SF_Vertex(nv)%Vn2 = 0d0
                Tdomain%SF%SF_Vertex(nv)%Vn3 = 0d0
                Tdomain%SF%SF_Vertex(nv)%save_veloc1 = 0d0
                Tdomain%SF%SF_Vertex(nv)%save_veloc2 = 0d0
                Tdomain%SF%SF_Vertex(nv)%save_veloc3 = 0d0
            else
                Tdomain%SF%SF_Vertex(nv)%save_displ = 0d0
            end if
        enddo

    end if


    !------------------------
    ! Inter-proc communications
    !------------------------
    if(Tdomain%tot_comm_proc > 0)then
        do n = 0,Tdomain%tot_comm_proc-1
            call allocate_comm_proc (Tdomain, Tdomain%sComm(n))
        enddo
    else
!        Tdomain%sComm(0)%ngll = 0
!        Tdomain%sComm(0)%ngllPML = 0
!        Tdomain%sComm(0)%ngllSO = 0
!        Tdomain%sComm(0)%ngllPML_F = 0
    endif

    return
end subroutine allocate_domain

    subroutine allocate_comm_proc(Tdomain, io_comm)
        type(domain), intent (INOUT) :: Tdomain
        type(comm), intent(inout) :: io_comm
        integer :: nf,ne,nv,i,ngll1,ngll2,   &
            ngll,ngllPML,ngllSO,ngllNeu,ngllSF,ngllSF_PML,ngll_F,ngllPML_F

            ngll = 0
            ngll_F = 0
            ngllPML = 0
            ngllSO = 0    ! XXX
            ngllPML_F = 0 ! XXX
            ngllNeu = 0
            ngllSF = 0
            ngllSF_PML = 0
            do i = 0,io_comm%nb_faces-1
                nf = io_comm%faces(i)
                ngll1 = Tdomain%sFace(nf)%ngll1 ; ngll2 = Tdomain%sFace(nf)%ngll2
                if(Tdomain%sFace(nf)%solid)then  ! solid face
                    ngll = ngll + (ngll1-2)*(ngll2-2)
                    if(Tdomain%sFace(nf)%PML)then
                        ngllPML = ngllPML + (ngll1-2)*(ngll2-2)
                    endif
                else   ! fluid face
                    ngll_F = ngll_F + (ngll1-2)*(ngll2-2)
                    if(Tdomain%sFace(nf)%PML)then
                        ngllPML_F = ngllPML_F + (ngll1-2)*(ngll2-2)
                    endif
                end if
            enddo
            do i = 0,io_comm%nb_edges-1
                ne = io_comm%edges(i)
                if(Tdomain%sEdge(ne)%solid)then  ! solid edge
                    ngll = ngll + Tdomain%sEdge(ne)%ngll-2
                    if(Tdomain%sEdge(ne)%PML) then
                        ngllPML = ngllPML + Tdomain%sEdge(ne)%ngll-2
                    endif
                else   ! fluid edge
                    ngll_F = ngll_F + Tdomain%sEdge(ne)%ngll-2
                    if(Tdomain%sEdge(ne)%PML) then
                        ngllPML_F = ngllPML_F + Tdomain%sEdge(ne)%ngll-2
                    endif
                end if
            enddo
            do i = 0,io_comm%nb_vertices-1
                nv = io_comm%vertices(i)
                if(Tdomain%sVertex(nv)%solid)then  ! solid vertex
                    ngll = ngll + 1
                    if(Tdomain%sVertex(nv)%PML) then
                        ngllPML = ngllPML + 1
                    endif
                else    ! fluid vertex
                    ngll_F = ngll_F + 1
                    if(Tdomain%sVertex(nv)%PML) then
                        ngllPML_F = ngllPML_F + 1
                    endif
                end if
            enddo
            ! Neumann
            do i = 0,io_comm%Neu_ne_shared-1
                ne = io_comm%Neu_edges_shared(i)
                ngllNeu = ngllNeu + Tdomain%Neumann%Neu_Edge(ne)%ngll-2
            enddo
            do i = 0,io_comm%Neu_nv_shared-1
                ngllNeu = ngllNeu + 1
            enddo
            ! Solid/fluid
            do i = 0,io_comm%SF_nf_shared-1
                nf = io_comm%SF_faces_shared(i)
                ngll1 = Tdomain%SF%SF_Face(nf)%ngll1; ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
                ngllSF = ngllSF + (ngll1-2)*(ngll2-2)
                if(Tdomain%SF%SF_Face(nf)%PML) ngllSF_PML =     &
                    ngllSF_PML+(ngll1-2)*(ngll2-2)
            enddo
            do i = 0,io_comm%SF_ne_shared-1
                ne = io_comm%SF_edges_shared(i)
                ngllSF = ngllSF + Tdomain%SF%SF_Edge(ne)%ngll-2
                if(Tdomain%SF%SF_Edge(ne)%PML) ngllSF_PML =     &
                    ngllSF_PML+Tdomain%SF%SF_Edge(ne)%ngll-2
            enddo
            do i = 0,io_comm%SF_nv_shared-1
                nv = io_comm%SF_vertices_shared(i)
                ngllSF = ngllSF + 1
                if(Tdomain%SF%SF_Vertex(nv)%PML) ngllSF_PML = ngllSF_PML+1
            enddo

            ! allocations for inter-proc communications
            ! for mass communications
            if(ngll+ngll_F > 0)then
                allocate(io_comm%Give(0:ngll+ngll_F-1))
                allocate(io_comm%Take(0:ngll+ngll_F-1))
            end if
            if(ngllPML+ngllPML_F > 0)then
                if(Tdomain%any_FPML)then
                    allocate(io_comm%GivePML(0:ngllPML+ngllPML_F-1,0:5))
                    allocate(io_comm%TakePML(0:ngllPML+ngllPML_F-1,0:5))
                else
                    allocate(io_comm%GivePML(0:ngllPML+ngllPML_F-1,0:2))
                    allocate(io_comm%TakePML(0:ngllPML+ngllPML_F-1,0:2))
                end if
            end if
            ! solid
            if(ngll > 0)then
                allocate(io_comm%GiveForces(0:ngll-1,0:2))
                allocate(io_comm%TakeForces(0:ngll-1,0:2))
            endif
            if(ngllPML > 0)then
                allocate(io_comm%GiveForcesPML(0:ngllPML-1,1:3,0:2))
                allocate(io_comm%TakeForcesPML(0:ngllPML-1,1:3,0:2))
            endif
            ! fluid
            if(ngll_F > 0)then
                allocate(io_comm%GiveForcesFl(0:ngll_F-1))
                allocate(io_comm%TakeForcesFl(0:ngll_F-1))
            endif
            if(ngllPML_F > 0)then
                allocate(io_comm%GiveForcesPMLFl(0:ngllPML_F-1,1:3))
                allocate(io_comm%TakeForcesPMLFl(0:ngllPML_F-1,1:3))
            endif
            ! Neumann
            if(ngllNeu > 0)then
                allocate(io_comm%GiveNeu(0:ngllNeu-1,0:2))
                allocate(io_comm%TakeNeu(0:ngllNeu-1,0:2))
            endif
            ! solid-fluid interface
            if(ngllSF > 0)then
                allocate(io_comm%GiveForcesSF_FtoS(0:ngllSF-1,0:2))
                allocate(io_comm%TakeForcesSF_FtoS(0:ngllSF-1,0:2))
                allocate(io_comm%GiveForcesSF_StoF(0:ngllSF-1))
                allocate(io_comm%TakeForcesSF_StoF(0:ngllSF-1))
            endif
            if(ngllSF_PML > 0)then
                allocate(io_comm%GiveForcesSF_FtoS_PML(0:ngllSF_PML-1,1:3,0:2))
                allocate(io_comm%TakeForcesSF_FtoS_PML(0:ngllSF_PML-1,1:3,0:2))
                allocate(io_comm%GiveForcesSF_StoF_PML(0:ngllSF_PML-1,1:3))
                allocate(io_comm%TakeForcesSF_StoF_PML(0:ngllSF_PML-1,1:3))
            endif

            io_comm%ngll = ngll
            io_comm%ngll_F = ngll_F
            io_comm%ngll_tot = ngll+ngll_F
            io_comm%ngllPML = ngllPML
            io_comm%ngllPML_F = ngllPML_F
            io_comm%ngllPML_tot = ngllPML+ngllPML_F
            io_comm%ngllNeu = ngllNeu
            io_comm%ngllSF = ngllSF
            io_comm%ngllSF_PML = ngllSF_PML

        end subroutine allocate_comm_proc
end module sdomain_alloc
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
