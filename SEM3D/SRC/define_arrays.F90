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
    use scomm
    use constants
    use dom_solid
    use dom_solid_dg
    use dom_fluid
    use dom_solidpml
    use dom_fluidpml
    implicit none
#include "index.h"

    public :: define_arrays

contains

    subroutine define_arrays(Tdomain)
        use constants
        use model_earthchunk
        use model_prem
        use mCourant
        use smirror
        implicit none

        type (domain), intent (INOUT), target :: Tdomain
        integer :: n, mat, rg, bnum, ee

        ! Copy Idom from element to domain_XXX
        ! Note : this must be done first as CPLM domain needs dom%Idom_ to build M, K, ...
        ! Note : the non-CPML domains (Solid, Fluid, ...) use specel%Idom instead of dom%Idom_ to build M, K...
        ! TODO : for non-CPML domains, use dom%Idom_ instead of specel%Idom (better memory/cache locality)
        do n = 0,Tdomain%n_elem-1
            bnum = Tdomain%specel(n)%lnum/VCHUNK
            ee = mod(Tdomain%specel(n)%lnum,VCHUNK)

            if      (Tdomain%specel(n)%domain==DM_SOLID_CG    ) then
                Tdomain%sdom%Idom_(:,:,:,bnum,ee)    = Tdomain%specel(n)%Idom
            else if (Tdomain%specel(n)%domain==DM_SOLID_CG_PML) then
                Tdomain%spmldom%Idom_(:,:,:,bnum,ee) = Tdomain%specel(n)%Idom
            else if (Tdomain%specel(n)%domain==DM_FLUID_CG    ) then
                Tdomain%fdom%Idom_(:,:,:,bnum,ee)    = Tdomain%specel(n)%Idom
            else if (Tdomain%specel(n)%domain==DM_FLUID_CG_PML) then
                Tdomain%fpmldom%Idom_(:,:,:,bnum,ee) = Tdomain%specel(n)%Idom
            else if (Tdomain%specel(n)%domain==DM_SOLID_DG) then
                Tdomain%sdomdg%Idom_(:,:,:,bnum,ee) = Tdomain%specel(n)%Idom
            else
                stop "unknown domain"
            end if
            !deallocate(Tdomain%specel(n)%Idom) ! TODO : to uncomment when non-CPML domains use dom%Idom_ instead of specel%Idom
        end do

        rg = Tdomain%rank

        if( Tdomain%earthchunk_isInit/=0) then
            call load_earthchunk(Tdomain%earthchunk_file, Tdomain%earthchunk_delta_lon, Tdomain%earthchunk_delta_lat)
        endif


        call init_materials(Tdomain)
        ! We need the timestep to continue with PML defs...
        call Compute_Courant(Tdomain,rg)
!
        call init_domains(Tdomain)
        do n = 0,Tdomain%n_elem-1
            mat = Tdomain%specel(n)%mat_index
            ! Compute MassMat
            call init_local_mass(Tdomain, Tdomain%specel(n))

            ! Computes DumpS, DumpMass (local),and for FPML :  Iv and Is
            if (Tdomain%specel(n)%domain==DM_SOLID_CG_PML) then
                call init_solidpml_properties(Tdomain, Tdomain%specel(n), Tdomain%sSubdomain(mat))
            end if
            if (Tdomain%specel(n)%domain==DM_FLUID_CG_PML) then
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

        ! Copy Idom from element to domain_XXX

        if ((Tdomain%n_DIRIC /= 0).and.(Tdomain%logicD%surfBC)) then
           call init_dirichlet_unstructuredMesh(Tdomain)
        else
             do n = 0,size(Tdomain%sSurfaces)-1
                if (trim(Tdomain%sSurfaces(n)%name)=="dirichlet") then
                   call init_dirichlet_surface(Tdomain, Tdomain%sSurfaces(n))
                   exit
                end if
             enddo
        endif

        ! Mirror
        if (Tdomain%use_mirror) call init_mirror(Tdomain)

        do n = 0,Tdomain%n_elem-1
            if(allocated(Tdomain%specel(n)%Idom)) &
                deallocate(Tdomain%specel(n)%Idom) ! TODO : delete when non-CPML domains use dom%Idom_ instead of specel%Idom
        end do
    end subroutine define_arrays

    subroutine init_dirichlet_unstructuredMesh(Tdomain)

         implicit none
         type (domain), intent (INOUT)  :: Tdomain
         integer                        :: ee, ns, s, n
         character(len=20)              :: char

         do ee = lbound(Tdomain%list_DIRICBC,1),ubound(Tdomain%list_DIRICBC,1)
            ns = Tdomain%list_DIRICBC(ee)
            do s = lbound(Tdomain%nsurfsource(ns)%index,1),ubound(Tdomain%nsurfsource(ns)%index,1)
               write(char,*) Tdomain%nsurfsource(ns)%index(s)
               block: &
               do n = 0,size(Tdomain%sSurfaces)-1
                  if (Tdomain%sSurfaces(n)%name=="surface"//adjustl(char(1:len_trim(char)))) then
                      call init_dirichlet_surface_MT(Tdomain, Tdomain%sSurfaces(n))
                      exit block
                  endif
               enddo block
            enddo
         enddo

    end subroutine init_dirichlet_unstructuredMesh

    subroutine dirichlet_gll_map(dirichlet, dirichlet_out)

        implicit none
        integer, dimension(:), intent( in  )              :: dirichlet
        integer, dimension(:), allocatable, intent(out)   :: dirichlet_out
        integer, dimension(:), allocatable                :: dummy
        integer                                           :: pp, n, m

        allocate(dummy(0:size(dirichlet)-1))
        pp=0
        do n=0,size(dirichlet)-1
           bloc : &
           do m=0,n
              if (dirichlet(m)==dirichlet(n)) exit bloc
           enddo bloc
           if (n==m) then
              dummy(pp) = dirichlet(n)
              pp=pp+1
           endif
        enddo
        allocate(dirichlet_out(0:pp-1))
        dirichlet_out = dummy(0:pp-1)
        deallocate(dummy)

    end subroutine dirichlet_gll_map

    subroutine init_dirichlet_surface_MT(Tdomain, surf)

        implicit none
        type (domain), intent (INOUT)      :: Tdomain
        type (SurfaceT), intent(INOUT)     :: surf
        integer, dimension(:), allocatable :: dummy
        !
        Tdomain%sdom%n_dirich = surf%surf_sl%nbtot
        if (Tdomain%sdom%n_dirich/=0) then
            if (allocated(Tdomain%sdom%dirich)) then
                allocate(dummy(0:size(Tdomain%sdom%dirich)-1))
                dummy=Tdomain%sdom%dirich
                deallocate(Tdomain%sdom%dirich)
                call dirichlet_gll_map((/dummy, surf%surf_sl%map/),Tdomain%sdom%dirich)
                Tdomain%sdom%n_dirich = size(Tdomain%sdom%dirich)
            else
               allocate(Tdomain%sdom%dirich(0:surf%surf_sl%nbtot-1))
               Tdomain%sdom%dirich = surf%surf_sl%map
            endif
        end if
        Tdomain%sdomdg%n_dirich = surf%surf_sldg%nbtot
        if (Tdomain%sdomdg%n_dirich/=0) then
            if (allocated(Tdomain%sdomdg%dirich)) then
                allocate(dummy(0:size(Tdomain%sdomdg%dirich)-1))
                dummy=Tdomain%sdomdg%dirich
                deallocate(Tdomain%sdomdg%dirich)
                call dirichlet_gll_map((/dummy, surf%surf_sldg%map/),Tdomain%sdomdg%dirich)
                Tdomain%sdomdg%n_dirich = size(Tdomain%sdomdg%dirich)
            else
               allocate(Tdomain%sdomdg%dirich(0:surf%surf_sldg%nbtot-1))
               Tdomain%sdomdg%dirich = surf%surf_sldg%map
            endif
        end if
        Tdomain%fdom%n_dirich = surf%surf_fl%nbtot
        if (Tdomain%fdom%n_dirich/=0) then
            if (allocated(Tdomain%fdom%dirich)) then
                allocate(dummy(0:size(Tdomain%fdom%dirich)-1))
                dummy=Tdomain%fdom%dirich
                deallocate(Tdomain%fdom%dirich)
                call dirichlet_gll_map((/dummy, surf%surf_fl%map/),Tdomain%fdom%dirich)
                Tdomain%fdom%n_dirich = size(Tdomain%fdom%dirich)
            else
                allocate(Tdomain%fdom%dirich(0:surf%surf_fl%nbtot-1))
                Tdomain%fdom%dirich = surf%surf_fl%map
            endif
        end if
        Tdomain%spmldom%n_dirich = surf%surf_spml%nbtot
        if (Tdomain%spmldom%n_dirich/=0) then
            if (allocated(Tdomain%spmldom%dirich)) then
                allocate(dummy(0:size(Tdomain%spmldom%dirich)-1))
                dummy=Tdomain%spmldom%dirich
                deallocate(Tdomain%spmldom%dirich)
                call dirichlet_gll_map((/dummy, surf%surf_spml%map/),Tdomain%spmldom%dirich)
                Tdomain%spmldom%n_dirich = size(Tdomain%spmldom%dirich)
            else
                allocate(Tdomain%spmldom%dirich(0:surf%surf_spml%nbtot-1))
                Tdomain%spmldom%dirich = surf%surf_spml%map
            endif
        end if
        Tdomain%fpmldom%n_dirich = surf%surf_fpml%nbtot
        if (Tdomain%fpmldom%n_dirich/=0) then
            if (allocated(Tdomain%fpmldom%dirich)) then
                allocate(dummy(0:size(Tdomain%fpmldom%dirich)-1))
                dummy=Tdomain%fpmldom%dirich
                deallocate(Tdomain%fpmldom%dirich)
                call dirichlet_gll_map((/dummy, surf%surf_fpml%map/),Tdomain%fpmldom%dirich)
                Tdomain%fpmldom%n_dirich = size(Tdomain%fpmldom%dirich)
            else
                allocate(Tdomain%fpmldom%dirich(0:surf%surf_fpml%nbtot-1))
                Tdomain%fpmldom%dirich = surf%surf_fpml%map
            endif
        end if
        if (allocated(dummy)) deallocate(dummy)

    end subroutine init_dirichlet_surface_MT


    subroutine init_dirichlet_surface(Tdomain, surf)
        type (domain), intent (INOUT) :: Tdomain
        type (SurfaceT), intent(INOUT) :: surf
        !
        Tdomain%sdom%n_dirich = surf%surf_sl%nbtot
        if (Tdomain%sdom%n_dirich/=0) then
            allocate(Tdomain%sdom%dirich(0:surf%surf_sl%nbtot-1))
            Tdomain%sdom%dirich = surf%surf_sl%map
        end if
        Tdomain%sdomdg%n_dirich = surf%surf_sldg%nbtot
        if (Tdomain%sdomdg%n_dirich/=0) then
            allocate(Tdomain%sdomdg%dirich(0:surf%surf_sldg%nbtot-1))
            Tdomain%sdomdg%dirich = surf%surf_sldg%map
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
        real(fpp) :: Mass

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

                call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                    Tdomain%Comm_data%Data(n)%IGiveSDG, Tdomain%sdomdg%MassMat, k)

                ! Domain SOLID PML
                if (Tdomain%Comm_data%Data(n)%nsolpml>0) then
#ifndef CPML
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%DumpMass, k)
#else
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%DumpMat, k)
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%MasUMat, k)
#endif
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%MassMat, k)
                end if

                ! Domain FLUID
                call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                    Tdomain%Comm_data%Data(n)%IGiveF, Tdomain%fdom%MassMat, k)

                ! Domain FLUID PML
                if (Tdomain%Comm_data%Data(n)%nflupml>0) then
#ifndef CPML
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%DumpMass, k)
#else
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%DumpMat, k)
                    call comm_give_data(Tdomain%Comm_data%Data(n)%Give, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%MasUMat, k)
#endif
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

                call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                    Tdomain%Comm_data%Data(n)%IGiveSDG, Tdomain%sdomdg%MassMat, k)

                ! Domain SOLID PML
                if (Tdomain%Comm_data%Data(n)%nsolpml>0) then
#ifndef CPML
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%DumpMass, k)
#else
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%DumpMat, k)
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%MasUMat, k)
#endif
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveSPML, Tdomain%spmldom%MassMat, k)
                end if

                ! Domain FLUID
                call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                    Tdomain%Comm_data%Data(n)%IGiveF,  Tdomain%fdom%MassMat, k)

                ! Domain FLUID PML
                if (Tdomain%Comm_data%Data(n)%nflupml>0) then
#ifndef CPML
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%DumpMass, k)
#else
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%DumpMat, k)
                    call comm_take_data(Tdomain%Comm_data%Data(n)%Take, &
                        Tdomain%Comm_data%Data(n)%IGiveFPML, Tdomain%fpmldom%MasUMat, k)
#endif
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

        if (Tdomain%sdom%nglltot    /= 0) Tdomain%sdom%MassMat(:) = 1d0/Tdomain%sdom%MassMat(:)
        if (Tdomain%sdomdg%nglltot  /= 0) Tdomain%sdomdg%MassMat(:) = 1d0/Tdomain%sdomdg%MassMat(:)
        if (Tdomain%fdom%nglltot    /= 0) Tdomain%fdom%MassMat(:) = 1d0/Tdomain%fdom%MassMat(:)
        if (Tdomain%spmldom%nglltot /= 0) Tdomain%spmldom%MassMat(:) = 1d0/Tdomain%spmldom%MassMat(:)
        if (Tdomain%fpmldom%nglltot /= 0) Tdomain%fpmldom%MassMat(:) = 1d0/Tdomain%fpmldom%MassMat(:)

    end subroutine inverse_mass_mat

    subroutine init_domains(Tdomain)
        type (domain), intent (INOUT), target :: Tdomain
        !
        if (Tdomain%sdom%nglltot /= 0) call init_domain_solid(Tdomain, Tdomain%sdom)
        if (Tdomain%sdomdg%nglltot /= 0) call init_domain_solid_dg(Tdomain, Tdomain%sdomdg)
        if (Tdomain%fdom%nglltot /= 0) call init_domain_fluid(Tdomain, Tdomain%fdom)
        if (Tdomain%spmldom%nglltot /= 0) call init_domain_solidpml(Tdomain, Tdomain%spmldom)
        if (Tdomain%fpmldom%nglltot /= 0) call init_domain_fluidpml(Tdomain, Tdomain%fpmldom)
    end subroutine init_domains

    subroutine init_materials(Tdomain)
        use build_prop_files
        type (domain), intent (INOUT), target :: Tdomain
        !
        integer :: mat, n
        logical :: isfile

        do mat = 0, Tdomain%n_mat-1

            ! compute rotation matrix if needed
            if (Tdomain%sSubdomain(mat)%is_sph) then
                call get_rotation_to_pole(Tdomain%sSubdomain(mat)%sph_args%theta_chk, &
                                          Tdomain%sSubdomain(mat)%sph_args%phi_chk, &
                                          Tdomain%sSubdomain(mat)%sph_args%R_to_pole_chk)
                Tdomain%sSubdomain(mat)%sph_args%R_from_pole_chk=transpose(Tdomain%sSubdomain(mat)%sph_args%R_to_pole_chk)
            end if
            isfile = Tdomain%sSubdomain(mat)%material_definition == MATERIAL_FILE
            if (isfile) then
                call init_prop_file(Tdomain%sSubdomain(mat))
            end if
            do n = 0,Tdomain%n_elem-1
                ! on traite tous les elements, materiau par materiau...
                if (Tdomain%specel(n)%mat_index/=mat) cycle

                ! Attribute elastic properties from material !!!
                ! Sets Lambda, Mu, Qmu, ... from mat
                call init_material_properties(Tdomain, Tdomain%specel(n), Tdomain%sSubdomain(mat))
            end do
            if (isfile) then
                call cleanup_prop_file(Tdomain%sSubdomain(mat))
            end if
        end do
    end subroutine init_materials

    subroutine init_material_properties(Tdomain, specel, mat)
        use build_prop_files
        use tensor_util
        type(domain), intent(inout) :: Tdomain
        type(element), intent(inout) :: specel
        type(subdomain), intent(in) :: mat
        !
        real(fpp), dimension(0:mat%NGLL-1,0:mat%NGLL-1,0:mat%NGLL-1) :: v0, v1
        real(fpp), dimension(0:mat%NGLL-1,0:mat%NGLL-1,0:mat%NGLL-1) :: lambda, mu, rho, nu, nlkp
        real(fpp), dimension(0:mat%NGLL-1,0:mat%NGLL-1,0:mat%NGLL-1) :: v0h, v1h, eta, A, C, F, L, M
        real(fpp), dimension(0:mat%NGLL-1,0:mat%NGLL-1,0:mat%NGLL-1) :: Qk,Qm
        real(fpp), dimension(1:6,1:6,0:mat%NGLL-1,0:mat%NGLL-1,0:mat%NGLL-1) :: Cij
        real(fpp), dimension(0:2) :: xx,xxr
        integer :: i,j,k,n,idef
        logical :: aniso
        ! integration de la prise en compte du gradient de proprietes
        aniso = .false.
        Qk = mat%Qpression
        Qm = mat%Qmu
        select case(mat%material_definition)
            case(MATERIAL_CONSTANT)
                !    on copie toujours le materiau de base
                !    si le flag gradient est actif alors on peut changer les proprietes
                rho = mat%Ddensity
                select case(mat%deftype)
                case(MATDEF_VP_VS_RHO,MATDEF_VP_VS_RHO_D)
                    v0 = mat%Pspeed
                    v1 = mat%Sspeed
                case(MATDEF_E_NU_RHO,MATDEF_E_NU_RHO_D)
                    v0 = mat%DE
                    v1 = mat%DNU
                case(MATDEF_LAMBDA_MU_RHO,MATDEF_LAMBDA_MU_RHO_d)
                    v0 = mat%DLambda
                    v1 = mat%DMu
                case(MATDEF_KAPPA_MU_RHO,MATDEF_KAPPA_MU_RHO_D)
                    v0 = mat%DKappa
                    v1 = mat%DMu
                case(MATDEF_HOOKE_RHO,MATDEF_HOOKE_RHO_D)
                    aniso = .true.
                    ! XXX TODO REQUIRED ATTENTION
                case(MATDEF_NLKP_VS_RHO,MATDEF_NLKP_VS_RHO_D)
                    v0 = mat%DNlkp
                    v1 = mat%Sspeed
                case(MATDEF_NU_VS_RHO,MATDEF_NU_VS_RHO_D)
                    v0 = mat%DNu
                    v1 = mat%Sspeed
                end select
            case( MATERIAL_FILE )
                ! XXX interpolate rho/v0/v1 from file
                call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(1), v0)
                call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(2), v1)
                call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(3), rho)
                select case(mat%deftype)
                case(MATDEF_VP_VS_RHO_D,MATDEF_E_NU_RHO_D,MATDEF_LAMBDA_MU_RHO_D,&
                    MATDEF_KAPPA_MU_RHO_D,MATDEF_HOOKE_RHO_D,MATDEF_NLKP_VS_RHO_D,&
                    MATDEF_NU_VS_RHO_D)
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(4), Qk)
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(5), Qm)
                case (MATDEF_VTI_ANISO)
                    ! v0, v1, rho call dans le commun 
                    aniso=.true.
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(4), v0h)
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(5), v1h)
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(6), eta)
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(7), Qk)
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(8), Qm)
                case (MATDEF_HOOKE_ANISO)
                    aniso=.true.
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(1),Cij(1,1,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(2),Cij(2,2,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(3),Cij(3,3,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(4),Cij(4,4,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(5),Cij(5,5,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(6),Cij(6,6,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(7),Cij(1,2,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(8),Cij(1,3,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(9),Cij(1,4,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(10),Cij(1,5,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(11),Cij(1,6,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(12),Cij(2,3,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(13),Cij(2,4,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(14),Cij(2,5,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(15),Cij(2,6,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(16),Cij(3,4,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(17),Cij(3,5,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(18),Cij(3,6,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(19),Cij(4,5,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(20),Cij(4,6,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(21),Cij(5,6,:,:,:))
                    call interpolate_elem_field(Tdomain, specel, mat, mat%prop_field(22), rho)
                end select
        end select

        select case(mat%deftype)
        case(MATDEF_VP_VS_RHO,MATDEF_VP_VS_RHO_D)
            nlkp = HUGE_VAL
            mu = rho * v1**2
            lambda = rho*(v0**2 - 2d0 *v1**2)
        case(MATDEF_E_NU_RHO,MATDEF_E_NU_RHO_D)
            nlkp = HUGE_VAL
            lambda = v0*v1/((1d0+v1)*(1d0-2d0*v1))
            mu = v0/(2d0*(1d0+v1))
        case(MATDEF_LAMBDA_MU_RHO,MATDEF_LAMBDA_MU_RHO_D)
            nlkp = HUGE_VAL
            lambda = v0
            mu = v1
        case(MATDEF_KAPPA_MU_RHO,MATDEF_KAPPA_MU_RHO_D)
            nlkp = HUGE_VAL
            mu = v1
            lambda = v0 - 2d0*mu/3d0
        case(MATDEF_HOOKE_RHO,MATDEF_HOOKE_RHO_D)
            ! XXX TODO
            Cij = 0
        case(MATDEF_NLKP_VS_RHO,MATDEF_NLKP_VS_RHO_D)
            nu = mat%DNu
            mu = rho*v1**2
            lambda = 2.0d0*nu*(rho*v1**2)/(1.0d0-2.0d0*nu)
            nlkp = v0
        case(MATDEF_NU_VS_RHO,MATDEF_NU_VS_RHO_D)
            nu = v0
            mu = rho*v1**2
            lambda = 2.0d0*nu*(rho*v1**2)/(1.0d0-2.0d0*nu)
            nlkp = HUGE_VAL
        case(MATDEF_VTI_ANISO)
            A = rho*v0h**2
            C = rho*v0**2
            L = rho*v1**2
            M = rho*v1h**2
            F = eta*(A-2.d0*L)
            !!! elastic tensor Cij in mandel notation
            Cij = 0.d0
            Cij(1,1,:,:,:) = C
            Cij(2,2,:,:,:) = A
            Cij(3,3,:,:,:) = A
            Cij(1,2,:,:,:) = F
            Cij(1,3,:,:,:) = F
            Cij(2,3,:,:,:) = A-2.d0*M
            Cij(4,4,:,:,:) = 2.d0*M
            Cij(5,5,:,:,:) = 2.d0*L
            Cij(6,6,:,:,:) = 2.d0*L
            do i = 2,6
                do j = 1,i-1
                    Cij(i,j,:,:,:) = Cij(j,i,:,:,:)
                end do
            end do
            do k = 0,mat%NGLL-1
                do j = 0,mat%NGLL-1
                    do i = 0,mat%NGLL-1
                        lambda(i,j,k)=lambda_from_Cij(Cij(:,:,i,j,k))
                        mu(i,j,k)=mu_from_Cij(Cij(:,:,i,j,k))
                    end do
                end do
            end do
            if (mat%is_sph) then
                do k = 0,mat%NGLL-1
                    do j = 0,mat%NGLL-1
                        do i = 0,mat%NGLL-1
                            idef = specel%Iglobnum(i,j,k)
                            do n=0,2
                                xx(n) = Tdomain%GlobCoord(n,idef)
                            end do
                            call cart2sph(xx,xxr,.false.)
                            call c_4tensor(Cij(:,:,i,j,k),xxr(1),xxr(2))
                        end do
                    end do
                end do
            end if
        case(MATDEF_HOOKE_ANISO)
            !Cij definis au moment lecture fichier material.spec
            do i = 2,6
                do j = 1,i-1
                    Cij(i,j,:,:,:) = Cij(j,i,:,:,:)
                end do
            end do
            do k = 0,mat%NGLL-1
                do j = 0,mat%NGLL-1
                    do i = 0,mat%NGLL-1
                        lambda(i,j,k)=lambda_from_Cij(Cij(:,:,i,j,k))
                        mu(i,j,k)=mu_from_Cij(Cij(:,:,i,j,k))
                    end do
                end do
            end do
        end select

        select case (specel%domain)
        case (DM_SOLID_CG)
            if (mat%deftype==MATDEF_VTI_ANISO .or. mat%deftype==MATDEF_HOOKE_ANISO) then
                call init_material_tensor_solid(Tdomain%sdom,specel%lnum,mat,rho,lambda,mu,Qk,Qm,Cij)
            else
                call init_material_properties_solid(Tdomain%sdom,specel%lnum,mat,rho,lambda,mu,Qk,Qm,nlkp,Tdomain%nl_flag)
            end if
        case (DM_SOLID_DG)
            call init_material_properties_solid_dg(Tdomain%sdomdg,specel%lnum,mat,rho,lambda,mu)
        case (DM_FLUID_CG)
            call init_material_properties_fluid(Tdomain%fdom,specel%lnum,mat,rho,lambda)
        case (DM_SOLID_CG_PML)
            call init_material_properties_solidpml(Tdomain%spmldom,specel%lnum,mat,rho,lambda,mu)
        case (DM_FLUID_CG_PML)
            call init_material_properties_fluidpml(Tdomain%fpmldom,specel%lnum,mat,rho,lambda)
        case default
            stop "unknown domain"
        end select
    end subroutine init_material_properties

    subroutine finalize_pml_properties(Tdomain)
        type (domain), intent (INOUT), target :: Tdomain
        !
        if (Tdomain%fpmldom%nbelem>0) call finalize_fluidpml_properties(Tdomain, Tdomain%fpmldom)
        if (Tdomain%spmldom%nbelem>0) call finalize_solidpml_properties(Tdomain, Tdomain%spmldom)
    end subroutine finalize_pml_properties

    subroutine init_local_mass(Tdomain, specel)
        type (domain), intent (INOUT), target :: Tdomain
        type (element), intent(inout) :: specel
        !
        integer :: i, j, k, ind
        real(fpp) :: Whei
        real(fpp), dimension(:), allocatable :: GLLw

        integer ngll

        ngll = domain_ngll(Tdomain, specel%domain)
        call domain_gllw(Tdomain, specel%domain, GLLw)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    Whei = GLLw(i)*GLLw(j)*GLLw(k)
                    ind = specel%Idom(i,j,k)
                    select case (specel%domain)
                        case (DM_SOLID_CG)
                            call init_local_mass_solid(Tdomain%sdom,specel,i,j,k,ind,Whei)
                        case (DM_SOLID_DG)
                            call init_local_mass_solid_dg(Tdomain%sdomdg,specel,i,j,k,ind,Whei)
                        case (DM_SOLID_CG_PML)
                            call init_local_mass_solidpml(Tdomain%spmldom,specel,i,j,k,ind,Whei)
                        case (DM_FLUID_CG)
                            call init_local_mass_fluid(Tdomain%fdom,specel,i,j,k,ind,Whei)
                        case (DM_FLUID_CG_PML)
                            call init_local_mass_fluidpml(Tdomain%fpmldom,specel,i,j,k,ind,Whei)
                        case default
                            stop "unknown domain"
                    end select
                enddo
            enddo
        enddo
        deallocate(GLLw)
    end subroutine init_local_mass

!    !---------------------------------------------------------------------------
!    !---------------------------------------------------------------------------
!    !---------------------------------------------------------------------------
!    !---------------------------------------------------------------------------
!    function materialIsConstant(Tdomain, mat) result(authorization)
!
!        !INPUTS
!        type (domain), intent (in), target :: Tdomain
!        type (subdomain), intent(in) :: mat
!
!        !OUTPUT
!        logical :: authorization
!
!        !LOCAL
!        integer :: assocMat
!
!        assocMat = mat%assocMat
!        authorization = .false.
!
!        if(Tdomain%sSubDomain(assocMat)%material_type == "S" .or. &
!            Tdomain%sSubDomain(assocMat)%material_type == "N" .or. &
!            Tdomain%sSubDomain(assocMat)%material_type == "P" .or. &
!            Tdomain%sSubDomain(assocMat)%material_type == "F" .or. &
!            Tdomain%sSubDomain(assocMat)%material_type == "L" .or. &
!            Tdomain%sSubDomain(assocMat)%material_type == "T") then
!            authorization = .true.
!        end if
!
!    end function materialIsConstant

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
