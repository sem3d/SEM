!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file define_arrays.f90
!! \brief Contient la routine Define_Arrays()
!! \author
!! \version 1.0
!! \date
!!
!<
module mdefinitions
    use sdomain
    use mpi
    use scomm
    use constants
    use dom_solid
    use dom_fluid
    use dom_solidpml
    use dom_fluidpml
    implicit none
#include "index.h"

    public :: define_arrays
    private :: define_PML_DumpEnd

contains

    subroutine Define_Arrays(Tdomain)
        use constants
        use model_earthchunk
        use model_prem
        use build_prop_files
        implicit none

        type (domain), intent (INOUT), target :: Tdomain
        integer :: n, mat, rg
        rg = Tdomain%rank

        if( Tdomain%earthchunk_isInit/=0) then
            call load_earthchunk(Tdomain%earthchunk_file, Tdomain%earthchunk_delta_lon, Tdomain%earthchunk_delta_lat)
        endif

        !Applying properties that were written on files
        if (rg == 0) write (*,*) "--> APPLYING PROPERTIES FILES "
        !- applying properties files
        call apply_prop_files (Tdomain, rg)

        !write (*,*) "--> after apply_prop_files "

        do n = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            if ( mat < 0 .or. mat >= Tdomain%n_mat ) then
                print*, "ERROR : inconsistent material index = ", mat
                stop
            end if
! Attribute elastic properties from material !!!
            ! Sets Lambda, Mu, Qmu, ... from mat
            call init_material_properties(Tdomain, Tdomain%specel(n), Tdomain%sSubdomain(mat))
            ! Compute MassMat
            call init_local_mass(Tdomain, Tdomain%specel(n))
            ! Computes DumpS, DumpMass (local),and for FPML :  Iv and Is
            if (Tdomain%specel(n)%domain==DM_SOLID_PML) then
                call init_solidpml_properties(Tdomain, Tdomain%specel(n), Tdomain%sSubdomain(mat))
            end if
            if (Tdomain%specel(n)%domain==DM_FLUID_PML) then
                call init_fluidpml_properties(Tdomain, Tdomain%specel(n), Tdomain%sSubdomain(mat))
            end if
        enddo
        ! Here we have local mass matrix (not assembled) on elements and
        ! each of faces, edges, vertices containing assembled (on local processor only) mass matrix
        if( Tdomain%earthchunk_isInit/=0) then
            call clean_earthchunk()
        endif
        !- defining Neumann properties (Btn: the complete normal term, ponderated
        !      by Gaussian weights)
#if 0
        if(Tdomain%logicD%neumann_local_present)then
            call define_FEV_Neumann(Tdomain)
        endif
#endif

        call assemble_mass_matrices(Tdomain)
        call finalize_pml_properties(Tdomain)
        call inverse_mass_mat(Tdomain)
        do n = 0,size(Tdomain%sSurfaces)-1
            if (trim(Tdomain%sSurfaces(n)%name)=="dirichlet") then
                call init_dirichlet_surface(Tdomain, Tdomain%sSurfaces(n))
                exit
            end if
        end do

        ! Copy Idom from element to domain_XXX
        do n = 0,Tdomain%n_elem-1
            if      (Tdomain%specel(n)%domain==DM_SOLID    ) then
                Tdomain%sdom%Idom_(:,:,:,Tdomain%specel(n)%lnum)    = Tdomain%specel(n)%Idom
            else if (Tdomain%specel(n)%domain==DM_SOLID_PML) then
                Tdomain%spmldom%Idom_(:,:,:,Tdomain%specel(n)%lnum) = Tdomain%specel(n)%Idom
            else if (Tdomain%specel(n)%domain==DM_FLUID    ) then
                Tdomain%fdom%Idom_(:,:,:,Tdomain%specel(n)%lnum)    = Tdomain%specel(n)%Idom
            else if (Tdomain%specel(n)%domain==DM_FLUID_PML) then
                Tdomain%fpmldom%Idom_(:,:,:,Tdomain%specel(n)%lnum) = Tdomain%specel(n)%Idom
            else
                stop "unknown domain"
            end if
            deallocate(Tdomain%specel(n)%Idom)
        end do

    end subroutine define_arrays


    subroutine init_dirichlet_surface(Tdomain, surf)
        type (domain), intent (INOUT) :: Tdomain
        type (SurfaceT), intent(INOUT) :: surf
        !
        Tdomain%sdom%n_dirich = surf%surf_sl%nbtot
        if (Tdomain%sdom%n_dirich/=0) then
            allocate(Tdomain%sdom%dirich(0:surf%surf_sl%nbtot-1))
            Tdomain%sdom%dirich = surf%surf_sl%map
        end if
        Tdomain%fdom%n_dirich = surf%surf_fl%nbtot
        if (Tdomain%fdom%n_dirich/=0) then
            allocate(Tdomain%fdom%dirich(0:surf%surf_fl%nbtot-1))
            Tdomain%fdom%dirich = surf%surf_fl%map
        end if
        Tdomain%spmldom%n_dirich = surf%surf_spml%nbtot
        if (Tdomain%spmldom%n_dirich/=0) then
            allocate(Tdomain%spmldom%dirich(0:surf%surf_spml%nbtot-1))
            Tdomain%spmldom%dirich = surf%surf_spml%map
        end if
        Tdomain%fpmldom%n_dirich = surf%surf_fpml%nbtot
        if (Tdomain%fpmldom%n_dirich/=0) then
            allocate(Tdomain%fpmldom%dirich(0:surf%surf_fpml%nbtot-1))
            Tdomain%fpmldom%dirich = surf%surf_fpml%map
        end if

    end subroutine init_dirichlet_surface

    subroutine assemble_mass_matrices(Tdomain)
        implicit none
        type (domain), intent (INOUT), target :: Tdomain

        integer :: n, indsol, indflu, indpml
        integer :: k
        real :: Mass

        ! Couplage à l'interface solide / PML
        do n = 0,Tdomain%intSolPml%surf0%nbtot-1
            indsol = Tdomain%intSolPml%surf0%map(n)
            indpml = Tdomain%intSolPml%surf1%map(n)
            if (indsol<0 .or. indsol>=Tdomain%sdom%nglltot) then
                write(*,*) "Pb indexation Sol"
                stop 1
            end if
            if (indpml<0 .or. indpml>=Tdomain%spmldom%nglltot) then
                write(*,*) "Pb indexation Sol-pml"
                stop 1
            end if
            Mass = Tdomain%sdom%MassMat(indsol) + Tdomain%spmldom%MassMat(indpml)
            Tdomain%sdom%MassMat(indsol) = Mass
            Tdomain%spmldom%MassMat(indpml) = Mass
        enddo

        ! Couplage à l'interface fluid / PML
        do n = 0,Tdomain%intFluPml%surf0%nbtot-1
            indflu = Tdomain%intFluPml%surf0%map(n)
            indpml = Tdomain%intFluPml%surf1%map(n)
            Mass = Tdomain%fdom%MassMat(indflu) + Tdomain%fpmldom%MassMat(indpml)
            Tdomain%fdom%MassMat(indflu) = Mass
            Tdomain%fpmldom%MassMat(indpml) = Mass
        enddo

        if(Tdomain%Comm_data%ncomm > 0)then
            do n = 0,Tdomain%Comm_data%ncomm-1
                ! Domain SOLID
                k = 0
                call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                    Tdomain%Comm_data%Data(n)%IGiveS, Tdomain%sdom%MassMat, k)

                ! Domain SOLID PML
                if (Tdomain%Comm_data%Data(n)%nsolpml>0) then
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%DumpMass, k)
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%MassMat, k)
                end if

                ! Domain FLUID
                call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                    Tdomain%Comm_data%Data(n)%IGiveF, Tdomain%fdom%MassMat, k)

                ! Domain FLUID PML
                if (Tdomain%Comm_data%Data(n)%nflupml>0) then
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%DumpMass, k)
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%MassMat, k)
                end if

                Tdomain%Comm_data%Data(n)%nsend = k
            end do

            ! Exchange
            call exchange_sem_var(Tdomain, 103, Tdomain%Comm_data)

            ! Take
            do n = 0,Tdomain%Comm_data%ncomm-1
                ! Domain SOLID
                k = 0
                call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                    Tdomain%Comm_data%Data(n)%IGiveS, Tdomain%sdom%MassMat, k)

                ! Domain SOLID PML
                if (Tdomain%Comm_data%Data(n)%nsolpml>0) then
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%DumpMass, k)
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%MassMat, k)
                end if

                ! Domain FLUID
                call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                    Tdomain%Comm_data%Data(n)%IGiveF,  Tdomain%fdom%MassMat, k)

                ! Domain FLUID PML
                if (Tdomain%Comm_data%Data(n)%nflupml>0) then
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%DumpMass, k)
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%MassMat, k)
                end if
            end do
        end if

!!! XXXXXXXXX QUID echange PML->SOLIDE inter CPU... ?
        return
    end subroutine assemble_mass_matrices
    !---------------------------------------------------------------------------------------
    subroutine inverse_mass_mat(Tdomain)
        type (domain), intent (INOUT), target :: Tdomain
        !

        if (Tdomain%sdom%nglltot    /= 0)    Tdomain%sdom%MassMat(:) = 1d0/Tdomain%sdom%MassMat(:)
        if (Tdomain%fdom%nglltot    /= 0)    Tdomain%fdom%MassMat(:) = 1d0/Tdomain%fdom%MassMat(:)
        if (Tdomain%spmldom%nglltot /= 0) Tdomain%spmldom%MassMat(:) = 1d0/Tdomain%spmldom%MassMat(:)
        if (Tdomain%fpmldom%nglltot /= 0) Tdomain%fpmldom%MassMat(:) = 1d0/Tdomain%fpmldom%MassMat(:)

    end subroutine inverse_mass_mat

    subroutine init_material_properties(Tdomain, specel, mat)
        type (domain), intent (INOUT), target :: Tdomain
        type (element), intent(inout) :: specel
        type (subdomain), intent(in) :: mat
        !

        ! integration de la prise en compte du gradient de proprietes

        select case(mat%material_definition)
            case(MATERIAL_CONSTANT)
                !    on copie toujours le materiau de base
                select case (specel%domain)
                    case (DM_SOLID)
                        call init_material_properties_solid(Tdomain%sdom,specel%lnum,-1,-1,-1,&
                             mat%DDensity,mat%DLambda,mat%DMu,mat%DKappa, mat)
                    case (DM_FLUID)
                        call init_material_properties_fluid(Tdomain%fdom,specel%lnum,-1,-1,-1,&
                             mat%DDensity,mat%DLambda)
                    case (DM_SOLID_PML)
                        call init_material_properties_solidpml(Tdomain%spmldom,specel%lnum,-1,-1,-1,&
                             mat%DDensity,mat%DLambda,mat%DMu)
                    case (DM_FLUID_PML)
                        call init_material_properties_fluidpml(Tdomain%fpmldom,specel%lnum,-1,-1,-1,&
                             mat%DDensity,mat%DLambda)
                    case default
                        stop "unknown domain"
                end select
                !    si le flag gradient est actif alors on peut changer les proprietes
            case( MATERIAL_EARTHCHUNK )
                call initialize_material_earthchunk(Tdomain, specel, Tdomain%GlobCoord, size(Tdomain%GlobCoord,2))
            case( MATERIAL_PREM )
                call initialize_material_prem(Tdomain, specel, Tdomain%GlobCoord, size(Tdomain%GlobCoord,2))
            case( MATERIAL_GRADIENT )
                !    on copie toujours le materiau de base
                select case (specel%domain)
                    case (DM_SOLID)
                        call init_material_properties_solid(Tdomain%sdom,specel%lnum,-1,-1,-1,&
                             mat%DDensity,mat%DLambda,mat%DMu,mat%DKappa,mat)
                    case (DM_FLUID)
                        call init_material_properties_fluid(Tdomain%fdom,specel%lnum,-1,-1,-1,&
                             mat%DDensity,mat%DLambda)
                    case (DM_SOLID_PML)
                        call init_material_properties_solidpml(Tdomain%spmldom,specel%lnum,-1,-1,-1,&
                             mat%DDensity,mat%DLambda,mat%DMu)
                    case (DM_FLUID_PML)
                        call init_material_properties_fluidpml(Tdomain%fpmldom,specel%lnum,-1,-1,-1,&
                             mat%DDensity,mat%DLambda)
                    case default
                        stop "unknown domain"
                end select
                !    si le flag gradient est actif alors on peut changer les proprietes
                if ( Tdomain%logicD%grad_bassin ) then
                    call initialize_material_gradient(Tdomain, specel)
                endif
            case( MATERIAL_RANDOM )
                !Don`t do anything, the basic properties were initialized by file
                if(materialIsConstant(Tdomain, mat)) then
                    select case (specel%domain)
                        case (DM_SOLID)
                            call init_material_properties_solid(Tdomain%sdom,specel%lnum,-1,-1,-1,&
                                 mat%DDensity,mat%DLambda,mat%DMu,mat%DKappa,mat)
                        case (DM_FLUID)
                            call init_material_properties_fluid(Tdomain%fdom,specel%lnum,-1,-1,-1,&
                                 mat%DDensity,mat%DLambda)
                        case (DM_SOLID_PML)
                            call init_material_properties_solidpml(Tdomain%spmldom,specel%lnum,-1,-1,-1,&
                                 mat%DDensity,mat%DLambda,mat%DMu)
                        case (DM_FLUID_PML)
                            call init_material_properties_fluidpml(Tdomain%fpmldom,specel%lnum,-1,-1,-1,&
                                 mat%DDensity,mat%DLambda)
                        case default
                            stop "unknown domain"
                    end select
                end if
        end select
    end subroutine init_material_properties

    subroutine finalize_pml_properties(Tdomain)
        type (domain), intent (INOUT), target :: Tdomain
        !
        if (Tdomain%spmldom%nbelem>0) call finalize_solidpml_properties(Tdomain%spmldom)
        if (Tdomain%fpmldom%nbelem>0) call finalize_fluidpml_properties(Tdomain%fpmldom)
    end subroutine finalize_pml_properties

    subroutine init_local_mass(Tdomain, specel)
        type (domain), intent (INOUT), target :: Tdomain
        type (element), intent(inout) :: specel
        !
        integer :: i, j, k, ind
        real(fpp) :: Whei
        real, dimension(:), allocatable :: GLLw

        integer ngll

        ngll = domain_ngll(Tdomain, specel%domain)
        call domain_gllw(Tdomain, specel%domain, GLLw)

        !- general (element) weighting: tensorial property..
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    Whei = GLLw(i)*GLLw(j)*GLLw(k)
                    ind = specel%Idom(i,j,k)
                    select case (specel%domain)
                        case (DM_SOLID)
                            call init_local_mass_solid(Tdomain%sdom,specel,i,j,k,ind,Whei)
                        case (DM_SOLID_PML)
                            call init_local_mass_solidpml(Tdomain%spmldom,specel,i,j,k,ind,Whei)
                        case (DM_FLUID)
                            call init_local_mass_fluid(Tdomain%fdom,specel,i,j,k,ind,Whei)
                        case (DM_FLUID_PML)
                            call init_local_mass_fluidpml(Tdomain%fpmldom,specel,i,j,k,ind,Whei)
                        case default
                            stop "unknown domain"
                    end select
                enddo
            enddo
        enddo
        deallocate(GLLw)
    end subroutine init_local_mass

    subroutine initialize_material_gradient(Tdomain, specel)
        type (domain), intent (INOUT), target :: Tdomain
        type (element), intent(inout) :: specel
        !
        integer :: i,j,k,ipoint,iflag
        real :: Mu,Lambda,Kappa
        integer :: imx,imy,imz
        real :: xp,yp,zp,xfact
        real :: zg1,zd1,zg2,zd2,zz1,zz2,zfact
        real :: xd1,xg1
        real :: zrho,zrho1,zrho2,zCp,zCp1,zCp2,zCs,zCs1,zCs2
        integer :: ngll
        integer :: icolonne,jlayer
        !
        ngll = domain_ngll(Tdomain, specel%domain)

        !    debut modification des proprietes des couches de materiaux
        !    bassin    voir programme Surface.f90

        !     n_layer nombre de couches
        !     n_colonne nombre de colonnes en x ici uniquement
        !     x_type == 0 on remet des materiaux  homogenes dans chaque bloc
        !     x_type == 1 on met des gradients pour chaque colonne en interpolant
        !     suivant z
        !       integer  :: n_colonne, n_layer, x_type
        !    x_coord correspond aux abscisses des colonnes
        !       real, pointer, dimension(:) :: x_coord
        !      z_layer profondeur de  linterface pour chaque x de colonne
        !      on definit egalement le materiaux par rho, Cp , Cs
        !       real, pointer, dimension(:,:) :: z_layer, z_rho, z_Cp, z_Cs


        !     on cherche tout d abord a localiser la maille a partir d un
        !     point de Gauss interne milieux (imx,imy,imz)
        imx = 1+(ngll-1)/2
        imy = 1+(ngll-1)/2
        imz = 1+(ngll-1)/2
        !     on impose qu une maille appartienne a un seul groupe de gradient de
        !     proprietes
        ipoint = specel%Iglobnum(imx,imy,imz)
        xp = Tdomain%GlobCoord(0,ipoint)
        yp = Tdomain%GlobCoord(1,ipoint)
        zp = Tdomain%GlobCoord(2,ipoint)
        iflag = 0
        if ( Tdomain%sBassin%x_type .eq. 2 ) then
            if ( zp .gt. Tdomain%sBassin%zmax) then
                iflag = 1
            endif
            if ( zp .lt. Tdomain%sBassin%zmin) then
                iflag = 1
            endif
        endif
        !  si iflag nul on peut faire les modifications  pour toute la maille
        if ( iflag .eq. 0 ) then
            icolonne = 0
            xfact = 0.D0
            do i = 1, Tdomain%sBassin%n_colonne
                if ( xp .ge. Tdomain%sBassin%x_coord(i-1) .and.  xp .lt. Tdomain%sBassin%x_coord(i) ) then
                    icolonne = i-1
                    xfact = (xp - Tdomain%sBassin%x_coord(i-1))/(Tdomain%sBassin%x_coord(i)-Tdomain%sBassin%x_coord(i-1))
                endif
            enddo

            jlayer = 0
            zfact = 0.D0
            do j = 1,Tdomain%sBassin%n_layer
                zg1 = Tdomain%sBassin%z_layer(icolonne,j-1)
                zd1 = Tdomain%sBassin%z_layer(icolonne+1,j-1)
                zz1 = zg1 + xfact*(zd1-zg1)
                zg2 = Tdomain%sBassin%z_layer(icolonne,j)
                zd2 = Tdomain%sBassin%z_layer(icolonne+1,j)
                zz2 = zg2 + xfact*(zd2-zg2)
                if ( zp .ge. zz1 .and. zp .lt. zz2 ) then
                    jlayer = j-1
                    zfact = ( zp -zz1)/(zz2-zz1)
                endif
            enddo
            !        limite du sous-domaine de gradient
            xg1 = Tdomain%sBassin%x_coord(icolonne)
            xd1 = Tdomain%sBassin%x_coord(icolonne+1)
            zg1 = Tdomain%sBassin%z_layer(icolonne,jlayer)
            zd1 = Tdomain%sBassin%z_layer(icolonne+1,jlayer)
            zg2 = Tdomain%sBassin%z_layer(icolonne,jlayer+1)
            zd2 = Tdomain%sBassin%z_layer(icolonne+1,jlayer+1)
            !
            zrho1 = Tdomain%sBassin%z_rho(icolonne,jlayer)
            zrho2 = Tdomain%sBassin%z_rho(icolonne,jlayer+1)
            zCp1 = Tdomain%sBassin%z_Cp(icolonne,jlayer)
            zCp2 = Tdomain%sBassin%z_Cp(icolonne,jlayer+1)
            zCs1 = Tdomain%sBassin%z_Cs(icolonne,jlayer)
            zCs2 = Tdomain%sBassin%z_Cs(icolonne,jlayer+1)

            !     boucle sur les points de Gauss de la maille
            !     xp, yp, zp coordonnees du point de Gauss
            do k = 0, ngll -1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        ipoint = specel%Iglobnum(i,j,k)
                        xp = Tdomain%GlobCoord(0,ipoint)
                        yp = Tdomain%GlobCoord(1,ipoint)
                        zp = Tdomain%GlobCoord(2,ipoint)
                        if ( Tdomain%sBassin%x_type .ge. 1 ) then
                            !    interpolations  pour le calcul du gradient
                            xfact = ( xp - xg1)/(xd1-xg1)
                            zz1 = zg1 + xfact*(zd1-zg1)
                            zz2 = zg2 + xfact*(zd2-zg2)
                            zfact = ( zp - zz1)/(zz2-zz1)
                        else
                            !   on met les memes proprietes dans toute la maille
                            zfact = 0.D0
                        endif
                        zrho   = zrho1 + zfact*(zrho2-zrho1)
                        zCp   = zCp1 + zfact*(zCp2-zCp1)
                        zCs   = zCs1 + zfact*(zCs2-zCs1)
                        !     calcul des coeffcients elastiques
                        Mu     = zrho*zCs*zCs
                        Lambda = zrho*(zCp*zCp - zCs*zCs)
                        Kappa  = Lambda + 2.D0*Mu/3.D0

                        select case (specel%domain)
                            case (DM_SOLID)
                                call init_material_properties_solid(Tdomain%sdom,specel%lnum,i,j,k,&
                                     zrho,Lambda,Mu,Kappa)
                            case (DM_FLUID)
                                call init_material_properties_fluid(Tdomain%fdom,specel%lnum,i,j,k,&
                                     zrho,Lambda)
                            case (DM_SOLID_PML)
                                call init_material_properties_solidpml(Tdomain%spmldom,specel%lnum,i,j,k,&
                                     zrho,Lambda,Mu)
                            case (DM_FLUID_PML)
                                call init_material_properties_fluidpml(Tdomain%fpmldom,specel%lnum,i,j,k,&
                                     zrho,Lambda)
                            case default
                                stop "unknown domain"
                        end select
                    enddo
                enddo
            enddo
            !    fin test iflag nul
        endif
        !    fin modification des proprietes des couches de materiaux
    end subroutine initialize_material_gradient

    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    function materialIsConstant(Tdomain, mat) result(authorization)

        !INPUTS
        type (domain), intent (in), target :: Tdomain
        type (subdomain), intent(in) :: mat

        !OUTPUT
        logical :: authorization

        !LOCAL
        integer :: assocMat

        assocMat = mat%assocMat
        authorization = .false.


        if(Tdomain%sSubDomain(assocMat)%material_type == "S" .or. &
            Tdomain%sSubDomain(assocMat)%material_type == "P" .or. &
            Tdomain%sSubDomain(assocMat)%material_type == "F" .or. &
            Tdomain%sSubDomain(assocMat)%material_type == "L" .or. &
            Tdomain%sSubDomain(assocMat)%material_type == "T") then
            authorization = .true.
        end if

    end function materialIsConstant

end module mdefinitions
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------

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
