!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_solidpml
    use sdomain
    use constants
    use champs_solidpml
    use selement
    use ssubdomains
    use sdomain
    use pml
    implicit none
#include "index.h"
#include "loops.h"
#include "gllopt.h"
contains

    subroutine allocate_champs_solidpml(dom, f)
        type(domain_solidpml), intent (INOUT) :: dom
        integer, intent(in) :: f
        !
        ! Allocate ONE more gll than necessary to use as a dummy target for
        ! indirections for fake elements.
        allocate(dom%champs(f)%ForcesPML(0:dom%nglltot,0:2,0:2))
        allocate(dom%champs(f)%VelocPML (0:dom%nglltot,0:2,0:2))
        dom%champs(f)%ForcesPML = 0d0
        dom%champs(f)%VelocPML = 0d0
    end subroutine allocate_champs_solidpml

    subroutine allocate_dom_solidpml (Tdomain, dom)
        use gll3d
        implicit none
        type(domain) :: TDomain
        type(domain_solidpml), intent (INOUT) :: dom
        !
        integer :: nbelem, ngll, nblocks, i
        !

        ngll   = dom%ngll
        nbelem = dom%nbelem
        if (ngll == 0) return ! Domain doesn't exist anywhere
        ! Initialisation poids, points des polynomes de lagranges aux point de GLL
        call init_dombase(dom)

        ! Glls are initialized first, because we can have faces of a domain without elements
        if(nbelem /= 0) then
            ! We can have glls without elements
            ! Do not allocate if not needed (save allocation/RAM)

            nblocks = dom%nblocks

            allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1,0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1,0:VCHUNK-1))
            allocate(dom%Mu_     (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1,0:VCHUNK-1))

            if(Tdomain%TimeD%velocity_scheme)then
                allocate(dom%Stress1_(0:ngll-1,0:ngll-1,0:ngll-1,0:5,0:nblocks-1,0:VCHUNK-1))
                allocate(dom%Stress2_(0:ngll-1,0:ngll-1,0:ngll-1,0:5,0:nblocks-1,0:VCHUNK-1))
                allocate(dom%Stress3_(0:ngll-1,0:ngll-1,0:ngll-1,0:5,0:nblocks-1,0:VCHUNK-1))
                dom%Stress1_(:,:,:,:,:,:) = 0d0
                dom%Stress2_(:,:,:,:,:,:) = 0d0
                dom%Stress3_(:,:,:,:,:,:) = 0d0
                allocate(dom%PMLDumpSx_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nblocks-1,0:VCHUNK-1))
                allocate(dom%PMLDumpSy_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nblocks-1,0:VCHUNK-1))
                allocate(dom%PMLDumpSz_(0:ngll-1,0:ngll-1,0:ngll-1,0:1,0:nblocks-1,0:VCHUNK-1))
                dom%PMLDumpSx_(:,:,:,:,:,:) = 0d0
                dom%PMLDumpSy_(:,:,:,:,:,:) = 0d0
                dom%PMLDumpSz_(:,:,:,:,:,:) = 0d0
                allocate(dom%Veloc(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:2,0:nblocks-1))
            endif
        end if
        ! Allocation et initialisation de champs0 pour les PML solides
        if (dom%nglltot /= 0) then
            do i=0,Tdomain%TimeD%nsubsteps
                call allocate_champs_solidpml(dom, i)
            end do
            allocate(dom%DumpV(0:dom%nglltot,0:1,0:2))
            allocate(dom%DumpMass(0:dom%nglltot,0:2))
            dom%DumpV = 0d0
            dom%DumpMass = 0d0
            dom%DumpV(dom%nglltot,1,:) = 1d0
            dom%DumpMass(dom%nglltot,:) = 1d0
        endif
        if(Tdomain%rank==0) write(*,*) "INFO - solid pml domain : ", dom%nbelem, " elements and ", dom%nglltot, " ngll pts"
    end subroutine allocate_dom_solidpml

    subroutine deallocate_dom_solidpml (dom)
        implicit none
        type(domain_solidpml), intent (INOUT) :: dom
        integer :: f

        if(allocated(dom%m_Density)) deallocate(dom%m_Density)
        if(allocated(dom%m_Lambda )) deallocate(dom%m_Lambda )
        if(allocated(dom%m_Mu     )) deallocate(dom%m_Mu     )

        if(allocated(dom%m_Stress1)) deallocate(dom%m_Stress1)
        if(allocated(dom%m_Stress2)) deallocate(dom%m_Stress2)
        if(allocated(dom%m_Stress3)) deallocate(dom%m_Stress3)
        if(allocated(dom%m_PMLDumpSx  ))      deallocate(dom%m_PMLDumpSx  )
        if(allocated(dom%m_PMLDumpSy  ))      deallocate(dom%m_PMLDumpSy  )
        if(allocated(dom%m_PMLDumpSz  ))      deallocate(dom%m_PMLDumpSz  )

        do f=0,1
            if(allocated(dom%champs(f)%ForcesPML)) deallocate(dom%champs(f)%ForcesPML)
            if(allocated(dom%champs(f)%VelocPML )) deallocate(dom%champs(f)%VelocPML )
        end do
        if(allocated(dom%DumpV    )) deallocate(dom%DumpV    )

        if(allocated(dom%DumpMass)) deallocate(dom%DumpMass)

        call deallocate_dombase(dom)
    end subroutine deallocate_dom_solidpml

    subroutine get_solidpml_dom_var(dom, lnum, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain_solidpml)                      :: dom
        integer, dimension(0:), intent(in)         :: out_variables
        integer                                    :: lnum
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldU, fieldV, fieldA
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: fieldP
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: P_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: S_energy
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1)     :: eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: sig_dev
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
                        fieldU(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_VITESSE) == 1) then
                        fieldV(i,j,k,:) = dom%champs(0)%VelocPml(ind,:,0) + &
                                          dom%champs(0)%VelocPml(ind,:,1) + &
                                          dom%champs(0)%VelocPml(ind,:,2)
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        fieldA(i,j,k,:) = dom%Massmat(ind) * (&
                            dom%champs(1)%ForcesPml(ind,:,0) + &
                            dom%champs(1)%ForcesPml(ind,:,1) + &
                            dom%champs(1)%ForcesPml(ind,:,2) )
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        fieldP(i,j,k) = 0d0
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        eps_vol(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        P_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        S_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        eps_dev(i,j,k,:) = 0.
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        sig_dev(i,j,k,:) = 0.
                    end if
                enddo
            enddo
        enddo
    end subroutine get_solidpml_dom_var

    subroutine init_domain_solidpml(Tdomain, dom)
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solidpml), intent(inout) :: dom
        !
        ! TODO : useless, kill this method, needed for build compatibility SolidPML / SolidCPML
        dom%dt = Tdomain%TimeD%dtmin
    end subroutine init_domain_solidpml

    subroutine start_domain_solidpml(Tdomain, spmldom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solidpml), intent(inout) :: spmldom
        !
        integer :: i

        !$acc enter data copyin(spmldom, spmldom%champs, spmldom%dirich) &
        !$acc&           copyin(spmldom%m_Stress1,spmldom%m_Stress2,spmldom%m_Stress3) &
        !$acc&           copyin(spmldom%hprime,spmldom%gllw, spmldom%m_Idom) &
        !$acc&           copyin(spmldom%m_InvGrad, spmldom%m_Jacob, spmldom%Veloc) &
        !$acc&           copyin(spmldom%m_PMLDumpSx,spmldom%m_PMLDumpSy,spmldom%m_PMLDumpSz) &
        !$acc&           copyin(spmldom%DumpMass,spmldom%DumpV,spmldom%m_Lambda,spmldom%m_Mu) &
        !$acc&           copyin(Tdomain%intSolPml) &
        !$acc&           copyin(Tdomain%intSolPml%surf0, Tdomain%intSolPml%surf0%map) &
        !$acc&           copyin(Tdomain%intSolPml%surf1, Tdomain%intSolPml%surf1%map) &
        !$acc&
        do i = 0,1
            !$acc enter data  copyin(spmldom%champs(i), spmldom%champs(i)%VelocPML, spmldom%champs(i)%ForcesPML)
        end do

    end subroutine start_domain_solidpml

    subroutine stop_domain_solidpml(Tdomain, spmldom)
        use sdomain
        type (domain), intent (INOUT), target :: Tdomain
        type(domain_solidpml), intent(inout) :: spmldom
        !
        integer :: i

        !$acc exit data delete(spmldom, spmldom%champs, spmldom%dirich) &
        !$acc&          delete(spmldom%m_Stress1,spmldom%m_Stress2,spmldom%m_Stress3) &
        !$acc&          delete(spmldom%hprime,spmldom%gllw, spmldom%m_Idom) &
        !$acc&          delete(spmldom%m_InvGrad, spmldom%m_Jacob, spmldom%Veloc) &
        !$acc&          delete(spmldom%m_PMLDumpSx,spmldom%m_PMLDumpSy,spmldom%m_PMLDumpSz) &
        !$acc&          delete(spmldom%DumpMass,spmldom%DumpV,spmldom%m_Lambda,spmldom%m_Mu) &
        !$acc&          delete(Tdomain%intSolPml%surf0, Tdomain%intSolPml%surf0%map) &
        !$acc&          delete(Tdomain%intSolPml%surf1, Tdomain%intSolPml%surf1%map) &
        !$acc&          delete(Tdomain%intSolPml) &
        !$acc&
        do i = 0,1
            !$acc exit data  delete(spmldom%champs(i)%VelocPML, spmldom%champs(i)%ForcesPML)
        end do

    end subroutine stop_domain_solidpml

    subroutine init_material_properties_solidpml(dom, lnum, mat, density, lambda, mu)
        type(domain_solidpml), intent(inout) :: dom
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
    end subroutine init_material_properties_solidpml

    subroutine init_local_mass_solidpml(dom,specel,i,j,k,ind,Whei)
        type(domain_solidpml), intent (INOUT) :: dom
        type (Element), intent (INOUT) :: specel
        integer :: i,j,k,ind
        real(fpp) Whei
        !
        integer :: bnum, ee
        bnum = specel%lnum/VCHUNK
        ee = mod(specel%lnum,VCHUNK)
        ! Solid.

        specel%MassMat(i,j,k) = Whei*dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)
    end subroutine init_local_mass_solidpml

    subroutine forces_int_sol_pml_mainloop(dom, i1)
        use m_calcul_forces_pml
        type(domain_solidpml), intent(inout) :: dom
        integer, intent(in) :: i1
        !
        select case(dom%ngll)
            NGLLDISPATCHCALL_4(calcul_forces_pml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_5(calcul_forces_pml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_6(calcul_forces_pml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_7(calcul_forces_pml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_8(calcul_forces_pml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_9(calcul_forces_pml,,(dom,dom%ngll,dom%champs(i1)))
            NGLLDISPATCHCALL_N(calcul_forces_pml,,(dom,dom%ngll,dom%champs(i1)))
        end select
    end subroutine forces_int_sol_pml_mainloop


    subroutine init_solidpml_properties(Tdomain,specel,mat)
        type (domain), intent (INOUT), target :: Tdomain
        type (element), intent(inout) :: specel
        type (subdomain), intent(in) :: mat
        !
        integer :: ngll, lnum
        real(fpp), dimension(:,:,:), allocatable :: temp_PMLx,temp_PMLy
        real(fpp), dimension(:,:,:), allocatable :: wx,wy,wz
        real(fpp) :: dt
        real(fpp), dimension(:,:,:,:), allocatable :: PMLDumpMass
        integer :: i,j,k,idx,m,ind
        real(fpp), dimension(:,:,:), allocatable   :: Vp
        real(fpp), dimension(:,:,:,:), allocatable :: coords
        !
        integer :: bnum, ee
        lnum = specel%lnum
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)

        dt = Tdomain%TimeD%dtmin

        ngll = domain_ngll(Tdomain, specel%domain)

        allocate(Vp(0:ngll-1,0:ngll-1,0:ngll-1))
        Vp = sqrt((Tdomain%spmldom%Lambda_ (:,:,:,bnum,ee) + &
            2. * Tdomain%spmldom%Mu_(:,:,:,bnum,ee))/Tdomain%spmldom%Density_(:,:,:,bnum,ee))

        !- definition of the attenuation coefficient in PMLs (alpha in the literature)
        allocate(wx(0:ngll-1,0:ngll-1,0:ngll-1))
        allocate(wy(0:ngll-1,0:ngll-1,0:ngll-1))
        allocate(wz(0:ngll-1,0:ngll-1,0:ngll-1))

        allocate(coords(0:ngll-1,0:ngll-1,0:ngll-1,0:2))

        DO K=0,ngll-1
            DO J=0,ngll-1
                DO I=0,ngll-1
                    idx = specel%Iglobnum(I,J,K)
                    coords(I,J,K,:) = Tdomain%GlobCoord(:,idx)
                END DO
            END DO
        END DO
        call define_alpha_PML(coords, 0, ngll, Vp, mat%pml_width, mat%pml_pos, mat%Apow, mat%npow, wx)
        call define_alpha_PML(coords, 1, ngll, Vp, mat%pml_width, mat%pml_pos, mat%Apow, mat%npow, wy)
        call define_alpha_PML(coords, 2, ngll, Vp, mat%pml_width, mat%pml_pos, mat%Apow, mat%npow, wz)

        !- M-PMLs
        if(Tdomain%logicD%MPML)then
            allocate(temp_PMLx(0:ngll-1,0:ngll-1,0:ngll-1))
            allocate(temp_PMLy(0:ngll-1,0:ngll-1,0:ngll-1))
            temp_PMLx(:,:,:) = wx(:,:,:)
            temp_PMLy(:,:,:) = wy(:,:,:)
            wx(:,:,:) = wx(:,:,:)+Tdomain%MPML_coeff*(wy(:,:,:)+wz(:,:,:))
            wy(:,:,:) = wy(:,:,:)+Tdomain%MPML_coeff*(temp_PMLx(:,:,:)+wz(:,:,:))
            wz(:,:,:) = wz(:,:,:)+Tdomain%MPML_coeff*(temp_PMLx(:,:,:)+temp_PMLy(:,:,:))
            deallocate(temp_PMLx,temp_PMLy)
        end if

        allocate(PMLDumpMass(0:ngll-1,0:ngll-1,0:ngll-1,0:2))
        PMLDumpMass = 0d0

        !- strong formulation for stresses. Dumped mass elements, convolutional terms.
        ! Compute DumpS(x,y,z) and DumpMass(0,1,2)
        call define_PML_DumpInit(ngll,dt,wx,specel%MassMat, &
            Tdomain%spmldom%PMLDumpSx_(:,:,:,:,bnum,ee),PMLDumpMass(:,:,:,0))
        call define_PML_DumpInit(ngll,dt,wy,specel%MassMat, &
            Tdomain%spmldom%PMLDumpSy_(:,:,:,:,bnum,ee),PMLDumpMass(:,:,:,1))
        call define_PML_DumpInit(ngll,dt,wz,specel%MassMat, &
            Tdomain%spmldom%PMLDumpSz_(:,:,:,:,bnum,ee),PMLDumpMass(:,:,:,2))
        deallocate(wx,wy,wz)

        ! Assemble dump mass
        do m = 0,2
            do k = 0,ngll-1
                do j = 0,ngll-1
                    do i = 0,ngll-1
                        ind = specel%Idom(i,j,k)
                        Tdomain%spmldom%DumpMass(ind,m) =   Tdomain%spmldom%DumpMass(ind,m) &
                                                          + PMLDumpMass(i,j,k,m)
                    enddo
                enddo
            enddo
        enddo
        if(allocated(PMLDumpMass)) deallocate(PMLDumpMass)

        deallocate(Vp)
    end subroutine init_solidpml_properties

    subroutine finalize_solidpml_properties(Tdomain,dom)
      type (domain), intent (INOUT), target :: Tdomain
      type (domain_solidpml), intent (INOUT), target :: dom
      !
      call define_PML_DumpEnd(dom%nglltot, dom%MassMat, dom%DumpMass, dom%DumpV)
    end subroutine finalize_solidpml_properties

    subroutine newmark_predictor_solidpml(dom, Tdomain, f0, f1)
        type(domain_solidpml), intent (INOUT) :: dom
        type (domain), intent (INOUT) :: Tdomain
        integer, intent(in) :: f0, f1
        !
        integer :: n, i, j, indsol, indpml
        real(fpp) :: bega, dt

        bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
        dt = Tdomain%TimeD%dtmin

!!        !$acc   parallel async(1)  &
!!        !$acc&  present(Tdomain,Tdomain%sdom,dom) &
!!        !$acc&  present(dom%champs(f0),dom%champs(f1)) &
!!        !$acc&  present(dom%champs(f0)%VelocPML) &
!!        !$acc&  present(dom%champs(f1)%VelocPML,dom%champs(f1)%ForcesPML) &
!!        !$acc&  present(Tdomain%sdom%champs(f0)%Veloc) &
!!        !$acc&  present(Tdomain%intSolPml%surf0, Tdomain%intSolPml%surf0%map) &
!!        !$acc&  present(Tdomain%intSolPml%surf1, Tdomain%intSolPml%surf1%map) &
!!        !$acc&

        !$acc   parallel async(1)  &
        !$acc&  present(Tdomain,dom) &
        !$acc&  present(dom%champs(f1)) &
        !$acc&  present(dom%champs(f1)%VelocPML) &
        !$acc&
        !$acc loop collapse(3)
        do n = 0,dom%nglltot
            do i=0,2
                do j=0,2
                    dom%champs(f1)%ForcesPML(n,i,j) = 0.
                end do
            end do
        end do
        !$acc end parallel

        !$acc   parallel async(1)  &
        !$acc&  present(Tdomain,Tdomain%spmldom,dom) &
        !$acc&  present(dom%champs(f0),dom%champs(f1)) &
        !$acc&  present(dom%champs(f0)%VelocPML) &
        !$acc&  present(dom%champs(f1)%VelocPML,dom%champs(f1)%ForcesPML) &
        !$acc&  present(Tdomain%sdom%champs(f0)%Veloc) &
        !$acc&  present(Tdomain%intSolPml%surf0, Tdomain%intSolPml%surf0%map) &
        !$acc&  present(Tdomain%intSolPml%surf1, Tdomain%intSolPml%surf1%map) &
        !$acc&
        !$acc loop private(indsol,indpml)
        do n = 0,Tdomain%intSolPml%surf0%nbtot-1
            ! Couplage à l'interface solide / PML
            indsol = Tdomain%intSolPml%surf0%map(n)
            indpml = Tdomain%intSolPml%surf1%map(n)
            do i=0,2
                !$acc atomic write
                dom%champs(f0)%VelocPML(indpml,i,0) = Tdomain%sdom%champs(f0)%Veloc(indsol,i)
                dom%champs(f0)%VelocPML(indpml,i,1) = 0.
                dom%champs(f0)%VelocPML(indpml,i,2) = 0.
            end do
        enddo
        ! Prediction
        !$acc end parallel

        !$acc   parallel async(1)  &
        !$acc&  present(dom) &
        !$acc&  present(dom%champs(f0),dom%champs(f1)) &
        !$acc&  present(dom%champs(f0)%VelocPML) &
        !$acc&  present(dom%champs(f1)%VelocPML,dom%champs(f1)%ForcesPML) &
        !$acc&
        !$acc loop collapse(3)
        do n = 0,dom%nglltot
            do i=0,2
                do j=0,2
                    dom%champs(f1)%VelocPML(n,i,j) = dom%champs(f0)%VelocPML(n,i,j) + &
                        dt*(0.5-bega)*dom%champs(f1)%ForcesPML(n,i,j)
                end do
            end do
        end do
        !$acc end parallel
    end subroutine newmark_predictor_solidpml

    subroutine newmark_corrector_solidpml(dom, dt, t, f0, f1)
        type(domain_solidpml), intent (INOUT) :: dom
        real(fpp), intent(in) :: dt, t
        integer, intent(in) :: f0, f1
        !
        integer  :: n, id, indpml, jp

        !$acc parallel async(1) &
        !$acc&   present(dom, dom%DumpV, dom%champs) &
        !$acc&   present(dom%champs(f0)%VelocPML, dom%champs(f1)%ForcesPML) &
        !$acc&   present(dom%dirich) &
        !$acc&
        !$acc loop collapse(3)
        do jp = 0,2
            do id = 0,2
                !$omp simd linear(n)
                do n = 0,dom%nglltot-1
                    dom%champs(f0)%VelocPML(n,id,jp) = dom%DumpV(n,0,jp) * dom%champs(f0)%VelocPML(n,id,jp) + &
                        dt * dom%DumpV(n,1,jp) * dom%champs(f1)%ForcesPML(n,id,jp)
                end do
            end do
        end do
        !TODO Eventuellement : DeplaPML(:,:) = DeplaPML(:,:) + dt * VelocPML(:,:)
        !$acc loop
        do n = 0, dom%n_dirich-1
            indpml = dom%dirich(n)
            do jp = 0,2
                do id = 0,2
                    dom%champs(f0)%VelocPML(indpml,id,jp) = 0.
                end do
            end do
        enddo
        !$acc end parallel
    end subroutine newmark_corrector_solidpml

    function solidpml_Pspeed(dom, lnum, i, j, k) result(Pspeed)
        type(domain_solidpml), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        !
        real(fpp) :: Pspeed, M
        integer :: bnum, ee
        bnum = lnum/VCHUNK
        ee = mod(lnum,VCHUNK)
        M = dom%Lambda_(i,j,k,bnum,ee) + 2.*dom%Mu_(i,j,k,bnum,ee)
        Pspeed = sqrt(M/dom%Density_(i,j,k,bnum,ee))
    end function solidpml_Pspeed

    subroutine couplage_pml_solid(Tdomain, sdom, spmldom, i0, i1)
        use dom_solid
        implicit none
        type(domain), intent(inout)  :: Tdomain
        type(domain_solid), intent (inout) :: sdom
        type(domain_solidpml), intent (inout) :: spmldom
        integer, intent(in) :: i0, i1
        !
        integer  :: n, indsol, indpml, i
        real(fpp) :: force
        !
        !$acc  parallel loop async(1) &
        !$acc& present(Tdomain, Tdomain%intSolPml, sdom, spmldom) &
        !$acc& present(Tdomain%intSolPml%surf0, Tdomain%intSolPml%surf0%map) &
        !$acc& present(Tdomain%intSolPml%surf1, Tdomain%intSolPml%surf1%map) &
        !$acc& present(sdom%champs(i1),spmldom%champs(i1)) &
        !$acc& present(sdom%champs(i1)%Veloc) &
        !$acc& present(spmldom%champs(i1)%ForcesPML) &
        !$acc&
        do n = 0,Tdomain%intSolPml%surf0%nbtot-1
            indsol = Tdomain%intSolPml%surf0%map(n)
            indpml = Tdomain%intSolPml%surf1%map(n)
            do i=0,2
                force = spmldom%champs(i1)%ForcesPML(indpml,i,0) + &
                    spmldom%champs(i1)%ForcesPML(indpml,i,1) + &
                    spmldom%champs(i1)%ForcesPML(indpml,i,2)
                !$acc atomic update
                sdom%champs(i1)%Veloc(indsol,i) = sdom%champs(i1)%Veloc(indsol,i) + force
            end do
        enddo

    end subroutine couplage_pml_solid

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
