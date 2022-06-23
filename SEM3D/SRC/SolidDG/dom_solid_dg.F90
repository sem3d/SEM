!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_solid_dg
    use constants
    use champs_solid_dg
    implicit none
#include "index.h"

contains

    subroutine allocate_champs_solid_dg(dom, i)
        type(domain_solid_dg), intent (INOUT) :: dom
        integer, intent(in) :: i
        integer :: ngll
        ngll = dom%ngll
        ! Allocate ONE more gll than necessary to use as a dummy target for
        ! indirections for fake elements.
        allocate(dom%champs(i)%Depla (0:dom%nglltot,0:2))
        allocate(dom%champs(i)%Veloc (0:dom%nglltot,0:2))
        ! Traces (last face==dummy)
        allocate(dom%champs(i)%tr_Veloc(0:ngll-1,0:ngll-1,0:2,0:dom%nbface))

        dom%champs(i)%Depla = 0d0
        dom%champs(i)%Veloc = 0d0
    end subroutine allocate_champs_solid_dg

    subroutine allocate_dom_solid_dg (Tdomain, dom)
        use sdomain
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_solid_dg), intent (INOUT) :: dom
        !
        integer :: nbelem, nblocks, ngll, i
        !

        ngll    = dom%ngll
        nbelem  = dom%nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call init_dombase(dom)

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! Do not allocate if not needed (save allocation/RAM)
            ! Wexo can have glls without elements
            nblocks = dom%nblocks

            allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Mu_     (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            ! 0:9 0: face, 1-3 : I0 4-6: DI 7-9:DJ  face(i,j) = elem(I0+i*DI+j*DJ)m 10: side
            allocate(dom%Itrace_     (0:5,0:10,0:nblocks-1, 0:VCHUNK-1))
        end if
        call compute_trace_numbering(Tdomain, dom)
        ! Allocation et initialisation de champs0 et champs1 pour les solides
        if (dom%nglltot /= 0) then
            do i=0,Tdomain%TimeD%nsubsteps
                call allocate_champs_solid_dg(dom, i)
            end do
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid DG domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_solid_dg

    subroutine compute_trace_numbering (Tdomain, dom)
        use sdomain
        use mindex
        implicit none
        type(domain), intent(INOUT) :: TDomain
        type(domain_solid_dg), intent (INOUT) :: dom
        !
        integer :: n, nf, lnum, eb, ec, side, k, ngll, nnf
        integer, dimension(0:3) :: elface
        integer, dimension(0:2) :: i0, di, dj
        
        !
        !Recollecting at the element level, from faces, edges and vertices.
        do n = 0,Tdomain%n_elem-1
            lnum = Tdomain%specel(n)%lnum
            eb = lnum/VCHUNK
            ec = mod(lnum,VCHUNK)
            !Taking information from faces
            do nf = 0,5
                nnf = Tdomain%specel(n)%Near_Faces(nf)
                do k=0,3
                    elface(k) = Tdomain%specel(n)%Control_nodes(face_def(k,nf))
                end do
                call ind_elem_face(ngll, nf, Tdomain%sFace(nnf)%inodes, elface, i0, di, dj)

                if (Tdomain%sFace(nnf)%elem_0==n) then
                    side = 0
                else if (Tdomain%sFace(nnf)%elem_1==n) then
                    side = 1
                else
                    write(*,*) "Inconsistency in compute_trace_numbering"
                    stop 1
                endif
                dom%Itrace_(nf, 0, ec, eb) = Tdomain%sFace(nnf)%lnum
                dom%Itrace_(nf, 1:3, ec, eb) = i0
                dom%Itrace_(nf, 4:6, ec, eb) = di
                dom%Itrace_(nf, 7:9, ec, eb) = dj
                dom%Itrace_(nf, 10, ec, eb)  = side

                ! How to use i0,di,dj
!                do i=1,Tdomain%sFace(nnf)%ngll-2
!                    do j=1,Tdomain%sFace(nnf)%ngll-2
!                        idxi = i0(0)+i*di(0)+j*dj(0)
!                        idxj = i0(1)+i*di(1)+j*dj(1)
!                        idxk = i0(2)+i*di(2)+j*dj(2)
!                        ind = Tdomain%sFace(nnf)%Iglobnum_Face(i,j)
!                        Tdomain%specel(n)%Iglobnum(idxi,idxj,idxk) = ind
!                        ind = Tdomain%sFace(nnf)%Idom(i,j)
!                        Tdomain%specel(n)%Idom(idxi,idxj,idxk) = ind
!                    end do
!                end do
            end do
        end do
    end subroutine compute_trace_numbering
    
    subroutine deallocate_dom_solid_dg (dom)
        implicit none
        type(domain_solid_dg), intent (INOUT) :: dom
        integer :: i

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_Mu     )) deallocate(dom%m_Mu     )

        do i=0,size(dom%champs)-1
            if(allocated(dom%champs(i)%Depla )) deallocate(dom%champs(i)%Depla )
            if(allocated(dom%champs(i)%Veloc )) deallocate(dom%champs(i)%Veloc )
        end do

        call deallocate_dombase(dom)

        end subroutine deallocate_dom_solid_dg


    subroutine get_solid_dg_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy,&
        eps_vol, eps_dev, sig_dev, dUdX)
        use deriv3d
        implicit none
        !
        type(domain_solid_dg), intent(inout)     :: dom
        integer, intent(in), dimension(0:)    :: out_variables
        integer, intent(in)                   :: lnum
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldU, fieldV, fieldA
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: fieldP
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: P_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: S_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:8) :: dUdX
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: sig_dev
        !
        logical                  :: flag_gradU
        integer                  :: ngll, i, j, k, ind
        real(fpp)                :: DXX, DXY, DXZ
        real(fpp)                :: DYX, DYY, DYZ
        real(fpp)                :: DZX, DZY, DZZ
        real(fpp)                :: divU
        real(fpp)                :: xmu, xlambda, xeps_vol
        real(fpp)                :: x2mu, xlambda2mu
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        !
        integer :: bnum, ee, flag_gradUint

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        xlambda = 0d0
        xmu = 0d0
        x2mu = 0d0
        flag_gradUint = 0
        flag_gradUint = out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS)

        flag_gradUint = flag_gradUint     + &
            out_variables(OUT_PRESSION)   + &
            out_variables(OUT_DUDX)       + &
            out_variables(OUT_EPS_VOL)    + &
            out_variables(OUT_EPS_DEV)    + &
            out_variables(OUT_STRESS_DEV)

        flag_gradU = flag_gradUint /= 0
        ngll = dom%ngll

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        fieldU(i,j,k,:) = dom%champs(0)%Depla(ind,:)
                    end if
                end do
            end do
        end do

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    if (flag_gradU .or. (out_variables(OUT_DUDX) == 1)) then
                        invgrad_ijk = dom%InvGrad_(:,:,i,j,k,bnum,ee) ! cache for performance
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,0),DXX,DYX,DZX)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,1),DXY,DYY,DZY)
                        call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                             invgrad_ijk,fieldU(:,:,:,2),DXZ,DYZ,DZZ)
                        divU = DXX+DYY+DZZ
                        if (out_variables(OUT_DUDX) == 1) then
                            dUdX(i,j,k,0) = DXX
                            dUdX(i,j,k,1) = DXY
                            dUdX(i,j,k,2) = DXZ
                            dUdX(i,j,k,3) = DYX
                            dUdX(i,j,k,4) = DYY
                            dUdX(i,j,k,5) = DYZ
                            dUdX(i,j,k,6) = DZX
                            dUdX(i,j,k,7) = DZY
                            dUdX(i,j,k,8) = DZZ
                        end if
                    end if
                    if (out_variables(OUT_VITESSE) == 1) fieldV(i,j,k,:) = dom%champs(0)%Veloc(ind,:)
                    if (out_variables(OUT_ACCEL) == 1) fieldA(i,j,k,:) = dom%champs(1)%Veloc(ind,:) ! dV/dt
                    if (flag_gradU) then
                        ! PRESSION
                        if (out_variables(OUT_PRESSION) == 1) then
                            fieldP(i,j,k) = -(dom%Lambda_(i,j,k,bnum,ee)+two*M_1_3*dom%Mu_(i,j,k,bnum,ee))*divU
                        endif
                        ! EPSVOL
                        if (out_variables(OUT_EPS_VOL) == 1) then
                            eps_vol(i,j,k) = divU
                        end if
                        ! ELASTIC-VISCOELASTIC MODULA
                        if (out_variables(OUT_ENERGYP) == 1 .or. &
                            out_variables(OUT_ENERGYS) == 1 .or. &
                            out_variables(OUT_STRESS_DEV) == 1) then
                            xmu     = dom%Mu_    (i,j,k,bnum,ee)
                            xlambda = dom%Lambda_(i,j,k,bnum,ee)
                            xeps_vol = DXX + DYY + DZZ
                            x2mu       = two * xmu
                            xlambda2mu = xlambda + x2mu
                        endif
                        ! P-ENERGY
                        if (out_variables(OUT_ENERGYP) == 1) then
                            P_energy(i,j,k) = ((0.5d0*xlambda) + xmu) * xeps_vol**2d0
                        end if
                        ! S-ENERGY
                        if (out_variables(OUT_ENERGYS) == 1) then
                            S_energy(i,j,k) = xmu/2.0d0 * (           &
                                                    (DZY - DYZ)**2d0  &
                                                  + (DXZ - DZX)**2d0  &
                                                  + (DYX - DXY)**2d0  &
                                                  )
                        end if
                        ! DEVIATORIC STRAIN
                        if (out_variables(OUT_EPS_DEV) == 1) then
                            eps_dev(i,j,k,0:5) = zero
                            eps_dev(i,j,k,0) = DXX - M_1_3 * divU
                            eps_dev(i,j,k,1) = DYY - M_1_3 * divU
                            eps_dev(i,j,k,2) = DZZ - M_1_3 * divU
                            eps_dev(i,j,k,3) = (DXY + DYX)
                            eps_dev(i,j,k,4) = (DZX + DXZ)
                            eps_dev(i,j,k,5) = (DZY + DYZ)
                        endif
                        ! DEVIATORIC STRESS
                        if (out_variables(OUT_STRESS_DEV) == 1) then
                            sig_dev(i,j,k,0:5) = zero
                            sig_dev(i,j,k,0) = x2mu * (DXX - M_1_3 * divU)
                            sig_dev(i,j,k,1) = x2mu * (DYY - M_1_3 * divU)
                            sig_dev(i,j,k,2) = x2mu * (DZZ - M_1_3 * divU)
                            sig_dev(i,j,k,3) =  xmu * (DXY + DYX)
                            sig_dev(i,j,k,4) =  xmu * (DXZ + DZX)
                            sig_dev(i,j,k,5) =  xmu * (DYZ + DZY)
                        endif
                    end if
                enddo
            enddo
        enddo
        !
    end subroutine get_solid_dg_dom_var


    subroutine get_solid_dg_dom_elem_energy(dom, lnum, P_energy, S_energy, R_energy, C_energy)
        use deriv3d
        implicit none
        !
        type(domain_solid_dg), intent(inout)          :: dom
        integer, intent(in)                        :: lnum
        real(fpp), dimension(:,:,:), allocatable, intent(inout) :: P_energy, S_energy, R_energy !R_energy = Residual energy (tend to zero as propagation takes place)
        real(fpp), dimension(:,:,:), allocatable, intent(inout) :: C_energy !Cinetic energy
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV

        integer                  :: ngll, i, j, k, ind
        real(fpp)                :: xmu, xlambda, xdensity
        real(fpp)                :: xeps_vol
        real(fpp), dimension(0:2,0:2) :: invgrad_ijk
        real(fpp), dimension(0:2) ::xvel
        real(fpp) :: dUx_dx,dUx_dy,dUx_dz
        real(fpp) :: dUy_dx,dUy_dy,dUy_dz
        real(fpp) :: dUz_dx,dUz_dy,dUz_dz
        !
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)


        ngll = dom%ngll

        !Dellocation
        if(allocated(S_energy)) then
            if(size(S_energy) /= ngll*ngll*ngll) deallocate(S_energy)
        end if

        if(allocated(P_energy)) then
            if(size(P_energy) /= ngll*ngll*ngll) deallocate(P_energy)
        end if

        if(allocated(R_energy)) then
            if(size(R_energy) /= ngll*ngll*ngll) deallocate(R_energy)
        end if

        if(allocated(C_energy)) then
            if(size(C_energy) /= ngll*ngll*ngll) deallocate(C_energy)
        end if

        if(allocated(fieldU)) then
            if(size(fieldU) /= ngll*ngll*ngll) deallocate(fieldU)
        end if

        if(allocated(fieldV)) then
            if(size(fieldV) /= ngll*ngll*ngll) deallocate(fieldV)
        end if

        !Allocation
        if(.not. allocated(S_energy)) allocate(S_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if(.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if(.not. allocated(R_energy)) allocate(R_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        if(.not. allocated(C_energy)) allocate(C_energy(0:ngll-1,0:ngll-1,0:ngll-1))
        S_energy = -1
        P_energy = -1


        if(.not. allocated(fieldU)) allocate(fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        if(.not. allocated(fieldV)) allocate(fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:2))

        ! First, get displacement.
        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    fieldU(i,j,k,:) = dom%champs(0)%Depla(ind,:)
                    fieldV(i,j,k,:) = dom%champs(0)%Veloc(ind,:)
                enddo
            enddo
        enddo

        ! Then, get the energies.
        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ! Compute gradU with displacement.
                    invgrad_ijk = dom%InvGrad_(:,:,i,j,k,bnum,ee) ! cache for performance

                    call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                         invgrad_ijk,fieldU(:,:,:,0),dUx_dx,dUx_dy,dUx_dz)
                    call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                         invgrad_ijk,fieldU(:,:,:,1),dUy_dx,dUy_dy,dUy_dz)
                    call physical_part_deriv_ijk(i,j,k,ngll,dom%hprime,&
                         invgrad_ijk,fieldU(:,:,:,2),dUz_dx,dUz_dy,dUz_dz)


                    ! Get other variables.
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    xeps_vol = dUx_dx + dUy_dy + dUz_dz

                    xmu     = dom%Mu_    (i,j,k,bnum,ee)
                    xlambda = dom%Lambda_(i,j,k,bnum,ee)
                    xdensity = dom%Density_ (i,j,k,bnum,ee)
                    xvel     = fieldV(i,j,k,:)

                    P_energy(i,j,k) = ((0.5d0*xlambda) + xmu) * xeps_vol**2d0
                    S_energy(i,j,k) = xmu/2.0d0 * (                       &
                                                    (dUz_dy - dUy_dz)**2d0  &
                                                  + (dUx_dz - dUz_dx)**2d0  &
                                                  + (dUy_dx - dUx_dy)**2d0  &
                                                  )
                    R_energy(i,j,k) = 2.0d0*xmu*(dUx_dy*dUy_dx + dUx_dz*dUz_dx + dUy_dz*dUz_dy) &
                                     -2.0d0*xmu*(dUx_dx*dUy_dy + dUx_dx*dUz_dz + dUy_dy*dUz_dz)

                    C_energy(i,j,k) = 0.5d0*xdensity*(xvel(0)**2.0d0 + xvel(1)**2.0d0 + xvel(2)**2.0d0)

                    !PAPER: The Energy Partitioning and the Diffusive Character of the Seismic Coda, Shapiro et al, 2000

                enddo
            enddo
        enddo

        if(allocated(fieldU)) deallocate(fieldU)
        if(allocated(fieldV)) deallocate(fieldV)

    end subroutine get_solid_dg_dom_elem_energy

    subroutine init_domain_solid_dg(Tdomain, dom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solid_dg), intent(inout) :: dom

        dom%dt = Tdomain%TimeD%dtmin
        dom%n_dirich = 0
    end subroutine init_domain_solid_dg

    subroutine init_material_properties_solid_dg(dom, lnum, mat, density, lambda, mu)
        use ssubdomains
        type(domain_solid_dg), intent(inout) :: dom
        integer, intent(in) :: lnum
        type (subdomain), intent(in) :: mat
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: density
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: lambda
        real(fpp), intent(in), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: mu
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        dom%Density_(:,:,:,bnum,ee) = density
        dom%Lambda_ (:,:,:,bnum,ee) = lambda
        dom%Mu_     (:,:,:,bnum,ee) = mu

    end subroutine init_material_properties_solid_dg

    subroutine init_local_mass_solid_dg(dom,specel,i,j,k,ind,Whei)
        use selement
        type(domain_solid_dg), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real(fpp) :: Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ! Solid.

        specel%MassMat(i,j,k) = Whei*dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_solid_dg

    subroutine forces_int_solid_dg(dom, var, dvdt, bnum)
        use m_calcul_forces_dg

        type(domain_solid_dg), intent (INOUT) :: dom
        type(champssolid_dg), intent(in)    :: var
        type(champssolid_dg), intent(inout) :: dvdt
        integer, intent(in) :: bnum
        !
        integer :: ngll,i,j,k,i_dir,ee,idx
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: Fox,Foy,Foz
        real(fpp), dimension(0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: Depla

        ngll = dom%ngll

        ! ELASTIC CASE: DISPLACEMENT PREDICTION
        do i_dir = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        do ee = 0, VCHUNK-1
                            idx = dom%Idom_(i,j,k,bnum,ee)
                            Depla(ee,i,j,k,i_dir) = var%Depla(idx,i_dir)
                        enddo
                    enddo
                enddo
            enddo
        enddo

        Fox = 0d0
        Foy = 0d0
        Foz = 0d0

        call calcul_forces_solid_dg_iso(dom,bnum,Fox,Foy,Foz,Depla)

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                        idx = dom%Idom_(i,j,k,bnum,ee)
                        dvdt%Veloc(idx,0) = dvdt%Veloc(idx,0)-Fox(ee,i,j,k)
                        dvdt%Veloc(idx,1) = dvdt%Veloc(idx,1)-Foy(ee,i,j,k)
                        dvdt%Veloc(idx,2) = dvdt%Veloc(idx,2)-Foz(ee,i,j,k)
                    enddo
               enddo
            enddo
        enddo
    end subroutine forces_int_solid_dg

    subroutine newmark_predictor_solid_dg(dom, f0, f1)
        type(domain_solid_dg), intent (INOUT) :: dom
        integer :: f0, f1
        !
        dom%champs(f1)%Veloc = 0d0
    end subroutine newmark_predictor_solid_dg

    subroutine newmark_corrector_solid_dg(dom, dt, f0, f1)
        type(domain_solid_dg), intent (INOUT) :: dom
        real(fpp):: dt
        integer, intent(in) :: f0, f1
        !
        integer :: i_dir, n, indpml
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            dom%champs(f1)%Veloc(indpml,:) = 0.
        enddo
        do i_dir = 0,2
!$omp simd linear(n)
            do n = 0,dom%nglltot-1
                dom%champs(f1)%Veloc(n,i_dir) = dom%champs(f1)%Veloc(n,i_dir) * dom%MassMat(n)
                dom%champs(f0)%Veloc(n,i_dir) = dom%champs(f0)%Veloc(n,i_dir) + dt * dom%champs(f1)%Veloc(n,i_dir)
                dom%champs(f0)%Depla(n,i_dir) = dom%champs(f0)%Depla(n,i_dir) + dt * dom%champs(f0)%Veloc(n,i_dir)
            end do
        enddo
    end subroutine newmark_corrector_solid_dg

    function solid_Pspeed_dg(dom, lnum, i, j, k) result(Pspeed)
        type(domain_solid_dg), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed, M
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        M = dom%Lambda_(i,j,k,bnum,ee) + 2.*dom%Mu_(i,j,k,bnum,ee)
        Pspeed = sqrt(M/dom%Density_(i,j,k,bnum,ee))
    end function solid_Pspeed_dg

    subroutine lddrk_init_solid_dg(dom, f2)
        type(domain_solid_dg), intent (INOUT) :: dom
        integer, intent(in) :: f2

        if (dom%nglltot==0) return
        dom%champs(f2)%Veloc = 0d0
    end subroutine lddrk_init_solid_dg

    subroutine lddrk_update_solid_dg(dom, f0, f1, f2, dt, cb, cg)
        type(domain_solid_dg), intent (INOUT) :: dom
        integer, intent(in) :: f0, f1, f2
        real(fpp), intent(in) :: cb, cg, dt
        integer :: i, n

        ! f2  contains forces computation  ie f2 = dU/dt
        ! f1 : w(n+1) = cb*w(n) + dt*dU/dt
        ! f0 : U(n+1) = U(n) + cg*w(n+1)
        if (dom%nglltot==0) return
        ! Only solid for starters
        do i = 0,2
            do n = 0,dom%nglltot-1
                dom%champs(f2)%Veloc(n,i) = dom%champs(f2)%Veloc(n,i) * dom%MassMat(n)
                dom%champs(f1)%Veloc(n,i) = cb*dom%champs(f1)%Veloc(n,i) + dt*dom%champs(f2)%Veloc(n,i)
                dom%champs(f1)%Depla(n,i) = cb*dom%champs(f1)%Depla(n,i) + dt*dom%champs(f0)%Veloc(n,i)
                dom%champs(f0)%Depla(n,i) = dom%champs(f0)%Depla(n,i) + cg*dom%champs(f1)%Depla(n,i)
                dom%champs(f0)%Veloc(n,i) = dom%champs(f0)%Veloc(n,i) + cg*dom%champs(f1)%Veloc(n,i)
            end do
        end do

!        do n = 0,dom%nglltot-1
!            dom%champs(f2)%Veloc(n,:) = dom%champs(f2)%Veloc(n,:) * dom%MassMat(n)
!            dom%champs(f2)%Eps(n,:)   = dom%champs(f2)%Eps  (n,:) * dom%MassMat(n)
!        end do
!        do n = 0,dom%nglltot-1
!            dom%champs(f1)%Veloc(n,:) = cb*dom%champs(f1)%Veloc(n,:) + dt*dom%champs(f2)%Veloc(n,:)
!            dom%champs(f1)%Eps(n,:)   = cb*dom%champs(f1)%Eps(n,:) + dt*dom%champs(f2)%Eps(n,:)
!        end do
!        do n = 0,dom%nglltot-1
!            dom%champs(f0)%Veloc(n,:) = dom%champs(f0)%Veloc(n,:) + cg*dom%champs(f1)%Veloc(n,:)
!            dom%champs(f0)%Eps(n,:)   = dom%champs(f0)%Eps  (n,:) + cg*dom%champs(f1)%Eps  (n,:)
!        end do
    end subroutine lddrk_update_solid_dg


end module dom_solid_dg
