!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

!! Convolutional Perfectly Matched Layers implemented according 2 references:
!! R1. Improved forward wave propagation and adjoint-based sensitivity kernel calculations using a numerically stable finite-element PML
!!     Zhinan Xie, Dimitri Komatitsch, Roland Martin, Rene Matzen
!!     Geophysical Journal International, 2014, 198, 1714-1747
!! R2. An efficient finite element time-domain formulation for the elastic second-order wave equation: a non-split complex frequency shifted convolutional PML
!!     Rene Matzen
!!     International Journal For Numerical Methods In Engineering, 2011, 88, 951-973

module dom_solidpml
    use constants
    use sdomain
    use champs_solidpml
    implicit none
#include "index.h"

    !! CPML parameters: for the very first implementation, parameters are hard-coded. Later on, they should be valuated according to input.spec
    real(fpp), private, parameter :: c_x = 1., c_y = 1., c_z = 1.
    integer,   private, parameter :: n_x = 2,  n_y = 2,  n_z = 2
    real(fpp), private, parameter :: r_c = 0.001
    integer,   private, parameter :: kappa_0 = 1, kappa_1 = 0
    real(fpp), private, parameter :: L_x = -1., L_y = -1., L_z = -1.

contains

    subroutine allocate_dom_solidpml (Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_solidpml), intent (INOUT) :: dom
        !
        integer nbelem, ngll, nblocks
        !

        ngll   = dom%ngll
        nbelem = dom%nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call compute_gll_data(ngll, dom%gllc, dom%gllw, dom%hprime, dom%htprime)

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)

            nblocks = ((nbelem+VCHUNK-1)/VCHUNK)
            dom%nblocks = nblocks

            allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Mu_     (0:ngll-1, 0:ngll-1, 0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))

            allocate (dom%Jacob_  (        0:ngll-1,0:ngll-1,0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            allocate (dom%InvGrad_(0:2,0:2,0:ngll-1,0:ngll-1,0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))

            allocate(dom%Idom_(0:ngll-1,0:ngll-1,0:ngll-1, 0:nblocks-1, 0:VCHUNK-1))
            dom%m_Idom = 0

            allocate(dom%Diagonal_Stress1_(0:ngll-1,0:ngll-1,0:ngll-1,0:2, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Diagonal_Stress2_(0:ngll-1,0:ngll-1,0:ngll-1,0:2, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Diagonal_Stress3_(0:ngll-1,0:ngll-1,0:ngll-1,0:2, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Residual_Stress1_(0:ngll-1,0:ngll-1,0:ngll-1,0:2, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Residual_Stress2_(0:ngll-1,0:ngll-1,0:ngll-1,0:2, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Residual_Stress3_(0:ngll-1,0:ngll-1,0:ngll-1,0:2, 0:nblocks-1, 0:VCHUNK-1))
            dom%Diagonal_Stress1_(:,:,:,:,:,:) = 0d0
            dom%Diagonal_Stress2_(:,:,:,:,:,:) = 0d0
            dom%Diagonal_Stress3_(:,:,:,:,:,:) = 0d0
            dom%Residual_Stress1_(:,:,:,:,:,:) = 0d0
            dom%Residual_Stress2_(:,:,:,:,:,:) = 0d0
            dom%Residual_Stress3_(:,:,:,:,:,:) = 0d0
            allocate(dom%PMLDumpSx_(0:ngll-1,0:ngll-1,0:ngll-1,0:1, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%PMLDumpSy_(0:ngll-1,0:ngll-1,0:ngll-1,0:1, 0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%PMLDumpSz_(0:ngll-1,0:ngll-1,0:ngll-1,0:1, 0:nblocks-1, 0:VCHUNK-1))
            dom%PMLDumpSx_(:,:,:,:,:,:) = 0d0
            dom%PMLDumpSy_(:,:,:,:,:,:) = 0d0
            dom%PMLDumpSz_(:,:,:,:,:,:) = 0d0
        end if

        ! Allocation et initialisation de champs0 pour les PML solides
        if (dom%nglltot /= 0) then
            allocate(dom%champs0%ForcesPML(0:dom%nglltot-1,0:2,0:2))
            allocate(dom%champs0%VelocPML (0:dom%nglltot-1,0:2,0:2))
            allocate(dom%champs0%DumpV    (0:dom%nglltot-1,0:2,0:2))
            allocate(dom%champs1%ForcesPML(0:dom%nglltot-1,0:2,0:2))
            allocate(dom%champs1%VelocPML (0:dom%nglltot-1,0:2,0:2))
            allocate(dom%champs1%DumpV    (0:dom%nglltot-1,0:2,0:2))
            dom%champs0%VelocPML  = 0d0
            dom%champs0%ForcesPML = 0d0
            dom%champs0%DumpV     = 0d0
            dom%champs1%VelocPML  = 0d0
            dom%champs1%ForcesPML = 0d0
            dom%champs1%DumpV     = 0d0

            ! Allocation de MassMat pour les PML solides
            allocate(dom%MassMat(0:dom%nglltot-1))
            dom%MassMat = 0d0

            allocate(dom%DumpMass(0:dom%nglltot-1,0:2))
            dom%DumpMass = 0d0
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid cpml domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_solidpml

    subroutine deallocate_dom_solidpml (dom)
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_Mu     )) deallocate(dom%m_Mu     )

        if(allocated(dom%m_Jacob  )) deallocate(dom%m_Jacob  )
        if(allocated(dom%m_InvGrad)) deallocate(dom%m_InvGrad)

        if(allocated(dom%m_Idom)) deallocate(dom%m_Idom)

        if(allocated(dom%gllc))    deallocate(dom%gllc)
        if(allocated(dom%gllw))    deallocate(dom%gllw)
        if(allocated(dom%hprime))  deallocate(dom%hprime)
        if(allocated(dom%htprime)) deallocate(dom%htprime)

        if(allocated(dom%m_Diagonal_Stress1)) deallocate(dom%m_Diagonal_Stress1)
        if(allocated(dom%m_Diagonal_Stress2)) deallocate(dom%m_Diagonal_Stress2)
        if(allocated(dom%m_Diagonal_Stress3)) deallocate(dom%m_Diagonal_Stress3)
        if(allocated(dom%m_Residual_Stress1)) deallocate(dom%m_Residual_Stress1)
        if(allocated(dom%m_Residual_Stress2)) deallocate(dom%m_Residual_Stress2)
        if(allocated(dom%m_Residual_Stress3)) deallocate(dom%m_Residual_Stress3)
        if(allocated(dom%m_PMLDumpSx  ))      deallocate(dom%m_PMLDumpSx  )
        if(allocated(dom%m_PMLDumpSy  ))      deallocate(dom%m_PMLDumpSy  )
        if(allocated(dom%m_PMLDumpSz  ))      deallocate(dom%m_PMLDumpSz  )

        if(allocated(dom%champs0%VelocPML )) deallocate(dom%champs0%VelocPML )
        if(allocated(dom%champs0%ForcesPML)) deallocate(dom%champs0%ForcesPML)
        if(allocated(dom%champs0%DumpV    )) deallocate(dom%champs0%DumpV    )
        if(allocated(dom%champs1%VelocPML )) deallocate(dom%champs1%VelocPML )
        if(allocated(dom%champs1%ForcesPML)) deallocate(dom%champs1%ForcesPML)
        if(allocated(dom%champs1%DumpV    )) deallocate(dom%champs1%DumpV    )

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)

        if(allocated(dom%DumpMass)) deallocate(dom%DumpMass)
    end subroutine deallocate_dom_solidpml

    subroutine get_solidpml_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain_solidpml)                      :: dom
        integer, dimension(0:8)                    :: out_variables
        integer                                    :: lnum
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        !
        logical :: flag_gradU
        integer :: ngll, i, j, k, ind
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        flag_gradU = (out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS) + &
            out_variables(OUT_EPS_VOL) + &
            out_variables(OUT_EPS_DEV) + &
            out_variables(OUT_STRESS_DEV)) /= 0

        ngll = dom%ngll

        do k=0,ngll-1
            do j=0,ngll-1
                do i=0,ngll-1
                    ind = dom%Idom_(i,j,k,bnum,ee)

                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        if(.not. allocated(fieldU)) allocate(fieldU(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldU(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_VITESSE) == 1) then
                        if(.not. allocated(fieldV)) allocate(fieldV(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldV(i,j,k,:) = dom%champs0%VelocPml(ind,:,0) + &
                                          dom%champs0%VelocPml(ind,:,1) + &
                                          dom%champs0%VelocPml(ind,:,2)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
                        fieldA(i,j,k,:) = dom%Massmat(ind) * ( dom%champs1%ForcesPml(ind,:,0) + &
                                                               dom%champs1%ForcesPml(ind,:,1) + &
                                                               dom%champs1%ForcesPml(ind,:,2) )
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:ngll-1,0:ngll-1,0:ngll-1))
                        fieldP(i,j,k) = 0d0
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        if(.not. allocated(eps_vol)) allocate(eps_vol(0:ngll-1,0:ngll-1,0:ngll-1))
                        eps_vol(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        if(.not. allocated(P_energy)) allocate(P_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                        P_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        if(.not. allocated(S_energy)) allocate(S_energy(0:ngll-1,0:ngll-1,0:ngll-1))
                        S_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        if(.not. allocated(eps_dev)) allocate(eps_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
                        eps_dev(i,j,k,:) = 0.
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        if(.not. allocated(sig_dev)) allocate(sig_dev(0:ngll-1,0:ngll-1,0:ngll-1,0:5))
                        sig_dev(i,j,k,:) = 0.
                    end if
                enddo
            enddo
        enddo
    end subroutine get_solidpml_dom_var

    subroutine init_material_properties_solidpml(dom, lnum, i, j, k, density, lambda, mu)
        type(domain_solidpml), intent(inout) :: dom
        integer, intent(in) :: lnum
        integer, intent(in) :: i, j, k ! -1 means :
        real(fpp), intent(in) :: density
        real(fpp), intent(in) :: lambda
        real(fpp), intent(in) :: mu
        !
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        if (i==-1 .and. j==-1 .and. k==-1) then
            dom%Density_(:,:,:,bnum,ee) = density
            dom%Lambda_ (:,:,:,bnum,ee) = lambda
            dom%Mu_     (:,:,:,bnum,ee) = mu
        else
            dom%Density_(i,j,k,bnum,ee) = density
            dom%Lambda_ (i,j,k,bnum,ee) = lambda
            dom%Mu_     (i,j,k,bnum,ee) = mu
        end if
    end subroutine init_material_properties_solidpml

    subroutine init_local_mass_solidpml(dom,specel,i,j,k,ind,Whei)
        type(domain_solidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)

        ! Solid.

        specel%MassMat(i,j,k) = Whei*dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_solidpml

    subroutine forces_int_sol_pml(dom, champs1, bnum)
        type(domain_solidpml), intent(inout) :: dom
        type(champssolidpml), intent(inout) :: champs1
        integer :: bnum
        !
        integer :: ngll
        integer :: i, j, k, l, ind, e, ee
        real, dimension(0:VCHUNK-1,0:2,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)  :: Forces1, Forces2, Forces3

        ngll = dom%ngll

        Forces1 = 0d0
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i=0,ngll-1
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    Forces1(ee,0,i,j,k) = Forces1(ee,0,i,j,k) + 0.
                    Forces1(ee,1,i,j,k) = Forces1(ee,1,i,j,k) + 0.
                    Forces1(ee,2,i,j,k) = Forces1(ee,2,i,j,k) + 0.
                    END_SUBELEM_LOOP()
                end do
            end do
        end do

        Forces2 = 0d0
        do k = 0,ngll-1
            do l = 0,ngll-1
                do j = 0,ngll-1
                    do i=0,ngll-1
                        BEGIN_SUBELEM_LOOP(e,ee,bnum)
                        Forces2(ee,0,i,j,k) = Forces2(ee,0,i,j,k) + 0.
                        Forces2(ee,1,i,j,k) = Forces2(ee,1,i,j,k) + 0.
                        Forces2(ee,2,i,j,k) = Forces2(ee,2,i,j,k) + 0.
                        END_SUBELEM_LOOP()
                    end do
                end do
            end do
        end do

        Forces3 = 0
        do l = 0,ngll-1
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i=0,ngll-1
                        BEGIN_SUBELEM_LOOP(e,ee,bnum)
                        Forces3(ee,0,i,j,k) = Forces3(ee,0,i,j,k) + 0.
                        Forces3(ee,1,i,j,k) = Forces3(ee,1,i,j,k) + 0.
                        Forces3(ee,2,i,j,k) = Forces3(ee,2,i,j,k) + 0.
                        END_SUBELEM_LOOP()
                    end do
                end do
            end do
        end do

        ! Assemblage
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    if (e>=dom%nbelem) exit
                    ind = dom%Idom_(i,j,k,bnum,ee)
                    champs1%ForcesPML(ind,:,0) = champs1%ForcesPML(ind,:,0) + Forces1(ee,:,i,j,k)
                    champs1%ForcesPML(ind,:,1) = champs1%ForcesPML(ind,:,1) + Forces2(ee,:,i,j,k)
                    champs1%ForcesPML(ind,:,2) = champs1%ForcesPML(ind,:,2) + Forces3(ee,:,i,j,k)
                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    end subroutine forces_int_sol_pml

    subroutine pred_sol_pml(dom, dt, champs1, bnum)
        implicit none

        type(domain_solidpml), intent(inout) :: dom
        type(champssolidpml), intent(inout) :: champs1
        real(fpp), intent(in) :: dt
        integer :: bnum
        !
        real(fpp), dimension (0:VCHUNK-1,0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: Veloc
        real(fpp) :: dVx_dx, dVx_dy, dVx_dz
        real(fpp) :: dVy_dx, dVy_dy, dVy_dz
        real(fpp) :: dVz_dx, dVz_dy, dVz_dz
        real(fpp) :: dS_dxi, dS_deta, dS_dzeta
        integer :: ngll
        integer :: i, j, k, l, ind, i_dir, e, ee

        ngll = dom%ngll

        do i_dir = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        BEGIN_SUBELEM_LOOP(e,ee,bnum)
                        ind = dom%Idom_(i,j,k,bnum,ee)
                        Veloc(ee,i,j,k,i_dir) = champs1%VelocPML(ind,i_dir,0) + &
                                                champs1%VelocPML(ind,i_dir,1) + &
                                                champs1%VelocPML(ind,i_dir,2)
                        END_SUBELEM_LOOP()
                    enddo
                enddo
            enddo
        enddo

        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
#ifdef SEM_VEC
!$omp simd linear(e,ee)
#endif
                    BEGIN_SUBELEM_LOOP(e,ee,bnum)
                    ! partial of velocity components with respect to xi,eta,zeta
                    part_deriv_ijke(Veloc,0,dS_dxi,dS_deta,dS_dzeta,dVx_dx,dVx_dy,dVx_dz)
                    part_deriv_ijke(Veloc,1,dS_dxi,dS_deta,dS_dzeta,dVy_dx,dVy_dy,dVy_dz)
                    part_deriv_ijke(Veloc,2,dS_dxi,dS_deta,dS_dzeta,dVz_dx,dVz_dy,dVz_dz)

                    END_SUBELEM_LOOP()
                enddo
            enddo
        enddo
    end subroutine pred_sol_pml

    subroutine init_solidpml_properties(Tdomain,specel,mat)
        type (domain), intent (INOUT), target :: Tdomain
        type (element), intent(inout) :: specel
        type (subdomain), intent(in) :: mat
        !
        integer :: ngll
        real(fpp), dimension(:,:,:,:), allocatable :: PMLDumpMass
        integer :: i,j,k,i_dir,ind

        ngll = domain_ngll(Tdomain, specel%domain)

        ! Compute dump mass
        allocate(PMLDumpMass(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        PMLDumpMass = 0d0

        ! Assemble dump mass
        do i_dir = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        ind = specel%Idom(i,j,k)
                        Tdomain%spmldom%DumpMass(ind,i_dir) =   Tdomain%spmldom%DumpMass(ind,i_dir) &
                                                              + PMLDumpMass(i,j,k,i_dir)
                    enddo
                enddo
            enddo
        enddo
        if(allocated(PMLDumpMass)) deallocate(PMLDumpMass)
    end subroutine init_solidpml_properties

    subroutine finalize_solidpml_properties(dom)
      type (domain_solidpml), intent (INOUT), target :: dom
      !
    end subroutine finalize_solidpml_properties

end module dom_solidpml

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
