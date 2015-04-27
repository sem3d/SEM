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
    integer :: n,nf,ne,nv,i,j,k,ngllx,nglly,ngllz,ngll1,ngll2,ngll
    integer :: n_solid, idx
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
        end if

        if(Tdomain%specel(n)%solid)then  ! SOLID PART
            allocate(Tdomain%specel(n)%sl)
            if(Tdomain%TimeD%velocity_scheme)then
                if(Tdomain%specel(n)%PML)then
                    allocate(Tdomain%specel(n)%sl%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:35))
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
                if(Tdomain%specel(n)%PML)then
                    allocate(Tdomain%specel(n)%fl%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:17))
                    allocate(Tdomain%specel(n)%flpml%Veloc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%flpml%Veloc1(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%flpml%Veloc2(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    allocate(Tdomain%specel(n)%flpml%Veloc3(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2))
                    Tdomain%specel(n)%flpml%Veloc = 0d0
                    Tdomain%specel(n)%flpml%Veloc1 = 0d0
                    Tdomain%specel(n)%flpml%Veloc2 = 0d0
                    Tdomain%specel(n)%flpml%Veloc3 = 0d0
                else
                    allocate(Tdomain%specel(n)%fl%Acoeff(0:ngllx-1,0:nglly-1,0:ngllz-1,0:5))
                endif
            endif

        end if
    enddo


    do n = 0, Tdomain%n_face-1
        ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
        if(Tdomain%sFace(n)%solid)then    ! SOLID PART
#ifdef COUPLAGE
            allocate (Tdomain%sFace(n)%ForcesExt(1:ngll1-2,1:ngll2-2,0:2 ) )
            Tdomain%sFace(n)%ForcesExt = 0.
            !        allocate (Tdomain%sFace(n)%FlagMka(1:ngll1-2,1:ngll2-2,0:2 ) )
            !        Tdomain%sFace(n)%FlagMka = 0
            allocate (Tdomain%sFace(n)%tsurfsem(1:ngll1-2,1:ngll2-2 ) )
            Tdomain%sFace(n)%tsurfsem = 0.
#endif
            if(Tdomain%sFace(n)%PML)then
                allocate(Tdomain%sFace(n)%spml)
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
            if(Tdomain%sFace(n)%PML)then
                allocate(Tdomain%sFace(n)%spml)
            endif
        endif
    enddo



    do n = 0,Tdomain%n_edge-1
        ngll = Tdomain%sEdge(n)%ngll
        if(Tdomain%sEdge(n)%solid)then   ! SOLID PART
            if(Tdomain%sEdge(n)%PML)then
                allocate(Tdomain%sEdge(n)%spml)
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

            endif
        else   ! FLUID PART
            if(Tdomain%sEdge(n)%PML)then
                allocate(Tdomain%sEdge(n)%spml)
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
#ifdef COUPLAGE
            allocate (Tdomain%sVertex(n)%ForcesExt(0:2) )
            Tdomain%sVertex(n)%ForcesExt = 0
            Tdomain%sVertex(n)%tsurfsem = 0.
#endif

            if(Tdomain%sVertex(n)%PML)then
                allocate(Tdomain%sVertex(n)%spml)
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
            endif
        else   ! FLUID PART
            if(Tdomain%sVertex(n)%PML)then
                allocate(Tdomain%sVertex(n)%spml)
            endif
        endif
    enddo

    ! solid-fluid related terms
    k = 0
    if(Tdomain%logicD%SF_local_present)then
        ! On alloue le tableau de BtN
        allocate(Tdomain%SF%SF_BtN(0:Tdomain%SF%ngll-1,0:2))
        Tdomain%SF%SF_BtN(:,:) = 0d0
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

    ! Allocation et initialisation de Tdomain%champs0 et champs1 pour les solides
    if (Tdomain%ngll_s /= 0) then
        allocate(Tdomain%champs0%Forces(0:Tdomain%ngll_s-1,0:2))
        allocate(Tdomain%champs0%Depla(0:Tdomain%ngll_s-1,0:2))
        allocate(Tdomain%champs0%Veloc(0:Tdomain%ngll_s-1,0:2))
        allocate(Tdomain%champs1%Forces(0:Tdomain%ngll_s-1,0:2))
        allocate(Tdomain%champs1%Depla(0:Tdomain%ngll_s-1,0:2))
        allocate(Tdomain%champs1%Veloc(0:Tdomain%ngll_s-1,0:2))

        Tdomain%champs0%Forces = 0d0
        Tdomain%champs0%Depla = 0d0
        Tdomain%champs0%Veloc = 0d0

        ! Allocation de Tdomain%MassMatSol pour les solides
        allocate(Tdomain%MassMatSol(0:Tdomain%ngll_s-1))
        Tdomain%MassMatSol = 0d0
    endif

    ! Allocation et initialisation de Tdomain%champs0 pour les PML solides
    if (Tdomain%ngll_pmls /= 0) then
        allocate(Tdomain%champs1%ForcesPML(0:Tdomain%ngll_pmls-1,0:2))
        allocate(Tdomain%champs0%VelocPML(0:Tdomain%ngll_pmls-1,0:2))
        allocate(Tdomain%champs0%DumpV(0:Tdomain%ngll_pmls-1,0:1))
        Tdomain%champs1%ForcesPML = 0d0
        Tdomain%champs0%VelocPML = 0d0
        Tdomain%champs0%DumpV = 0d0

        ! Allocation de Tdomain%MassMatSolPml pour les PML solides
        allocate(Tdomain%MassMatSolPml(0:Tdomain%ngll_pmls-1))
        Tdomain%MassMatSolPml = 0d0

        allocate(Tdomain%DumpMass(0:Tdomain%ngll_pmls-1))
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

        ! Allocation et Initialisation de Fluid_dirich
        ! permet d'annuler VelPhi sur les faces libres dans Newmark_corrector
        allocate(Tdomain%fl_dirich(0:Tdomain%ngll_f-1))
        Tdomain%fl_dirich = 1.0
        do n = 0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%Solid) cycle
            if (Tdomain%specel(n)%PML) cycle
            if (.not. Tdomain%specel(n)%fluid_dirich) cycle
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            do j = 0,nglly-1
                do i = 0,ngllx-1
                    idx = Tdomain%specel(n)%IFlu(i,j,ngllz-1)
                    if (idx==-1) stop "Error"
                    Tdomain%fl_dirich(idx) = 0.
                enddo
            enddo
        enddo

        ! Allocation de Tdomain%MassMatFlu pour les fluides
        allocate(Tdomain%MassMatFlu(0:Tdomain%ngll_f-1))
        Tdomain%MassMatFlu = 0d0
    endif

    ! Allocation et initialisation de Tdomain%champs0 pour les PML fluides
    if (Tdomain%ngll_pmlf /= 0) then
        allocate(Tdomain%champs1%fpml_Forces(0:Tdomain%ngll_pmlf-1))
        allocate(Tdomain%champs0%fpml_VelPhi(0:Tdomain%ngll_pmlf-1))
        allocate(Tdomain%champs0%fpml_Phi(0:Tdomain%ngll_pmlf-1))
        allocate(Tdomain%champs0%fpml_DumpV(0:Tdomain%ngll_pmlf-1,0:1))
        Tdomain%champs1%fpml_Forces = 0d0
        Tdomain%champs0%fpml_VelPhi = 0d0
        Tdomain%champs0%fpml_Phi = 0d0
        Tdomain%champs0%fpml_DumpV = 0d0

        ! Allocation de Tdomain%MassMatSolPml pour les PML solides
        allocate(Tdomain%MassMatFluPml(0:Tdomain%ngll_pmlf-1))
        Tdomain%MassMatFluPml = 0d0

        allocate(Tdomain%fpml_DumpMass(0:Tdomain%ngll_pmlf-1))
        Tdomain%fpml_DumpMass = 0d0

        ! Allocation et Initialisation de Fluid_dirich
        ! permet d'annuler VelPhi sur les faces libres dans Newmark_corrector
        allocate(Tdomain%fpml_dirich(0:Tdomain%ngll_pmlf-1))
        Tdomain%fpml_dirich = 1.0
        do n = 0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%Solid) cycle
            if (.not. Tdomain%specel(n)%PML) cycle
            if (.not. Tdomain%specel(n)%fluid_dirich) cycle
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            do j = 0,nglly-1
                do i = 0,ngllx-1
                    idx = Tdomain%specel(n)%flpml%IFluPml(i,j,ngllz-1)
                    if (idx==-1) stop "Error"
                    Tdomain%fpml_dirich(idx+0) = 0.
                    Tdomain%fpml_dirich(idx+1) = 0.
                    Tdomain%fpml_dirich(idx+2) = 0.
                enddo
            enddo
        enddo
    endif

    ! Allocation et initialisation des champs pour le couplage solide / fluide
    if (Tdomain%logicD%SF_local_present)then
        allocate(Tdomain%champs0%Save_forces(0:Tdomain%SF%ngll-1,0:2))
        allocate(Tdomain%champs0%Save_depla(0:Tdomain%SF%ngll-1,0:2))
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
