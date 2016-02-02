!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!

module dom_fluidpml
    use sdomain
    use constants
    use champs_fluidpml
    use selement
    use sdomain
    use ssubdomains
    implicit none

contains

  subroutine allocate_dom_fluidpml (Tdomain, dom)
        implicit none
        type(domain) :: TDomain
        type(domain_fluidpml) :: dom
        !
        integer nbelem, ngllx, nglly, ngllz
        !

        dom%ngllx = Tdomain%specel(0)%ngllx ! Temporaire: ngll* doit passer sur le domaine a terme
        dom%nglly = Tdomain%specel(0)%nglly ! Temporaire: ngll* doit passer sur le domaine a terme
        dom%ngllz = Tdomain%specel(0)%ngllz ! Temporaire: ngll* doit passer sur le domaine a terme

        nbelem  = dom%nbelem
        ngllx   = dom%ngllx
        nglly   = dom%nglly
        ngllz   = dom%ngllz

        if(Tdomain%TimeD%velocity_scheme)then
            allocate(dom%Veloc(0:ngllx-1,0:nglly-1,0:ngllz-1,0:2,0:nbelem-1))
            dom%Veloc = 0d0
            allocate(dom%PMLDumpSx  (0:ngllx-1,0:nglly-1,0:ngllz-1,0:1,0:nbelem-1))
            allocate(dom%PMLDumpSy  (0:ngllx-1,0:nglly-1,0:ngllz-1,0:1,0:nbelem-1))
            allocate(dom%PMLDumpSz  (0:ngllx-1,0:nglly-1,0:ngllz-1,0:1,0:nbelem-1))
            dom%PMLDumpSx   = 0d0
            dom%PMLDumpSy   = 0d0
            dom%PMLDumpSz   = 0d0
        endif

        ! Allocation et initialisation de champs0 pour les PML fluides
        if (dom%ngll /= 0) then
            allocate(dom%champs1%fpml_Forces(0:dom%ngll-1,0:2))
            allocate(dom%champs0%fpml_VelPhi(0:dom%ngll-1,0:2))
            allocate(dom%champs0%fpml_Phi   (0:dom%ngll-1,0:2))
            allocate(dom%champs1%fpml_VelPhi(0:dom%ngll-1,0:2))
            allocate(dom%champs0%fpml_DumpV (0:dom%ngll-1,0:1,0:2))
            dom%champs1%fpml_Forces = 0d0
            dom%champs0%fpml_VelPhi = 0d0
            dom%champs0%fpml_Phi = 0d0
            dom%champs0%fpml_DumpV = 0d0

            ! Allocation de MassMat pour les PML fluides
            allocate(dom%MassMat(0:dom%ngll-1))
            dom%MassMat = 0d0

            allocate(dom%DumpMass(0:dom%ngll-1,0:2))
            dom%DumpMass = 0d0
        endif
    end subroutine allocate_dom_fluidpml

    subroutine deallocate_dom_fluidpml (dom)
        implicit none
        type(domain_fluidpml) :: dom

        if(allocated(dom%Veloc))       deallocate(dom%Veloc      )
        if(allocated(dom%PMLDumpSx  )) deallocate(dom%PMLDumpSx  )
        if(allocated(dom%PMLDumpSy  )) deallocate(dom%PMLDumpSy  )
        if(allocated(dom%PMLDumpSz  )) deallocate(dom%PMLDumpSz  )

        if(allocated(dom%champs1%fpml_Forces)) deallocate(dom%champs1%fpml_Forces)
        if(allocated(dom%champs0%fpml_VelPhi)) deallocate(dom%champs0%fpml_VelPhi)
        if(allocated(dom%champs0%fpml_Phi   )) deallocate(dom%champs0%fpml_Phi   )
        if(allocated(dom%champs1%fpml_VelPhi)) deallocate(dom%champs1%fpml_VelPhi)
        if(allocated(dom%champs0%fpml_DumpV )) deallocate(dom%champs0%fpml_DumpV )

        if(allocated(dom%MassMat)) deallocate(dom%MassMat)

        if(allocated(dom%DumpMass)) deallocate(dom%DumpMass)
    end subroutine deallocate_dom_fluidpml

    subroutine get_fluidpml_dom_var(Tdomain, el, out_variables, &
        fieldU, fieldV, fieldA, fieldP, P_energy, S_energy, eps_vol, eps_dev, sig_dev)
        implicit none
        !
        type(domain)                               :: TDomain
        integer, dimension(0:8)                    :: out_variables
        type(element)                              :: el
        real(fpp), dimension(:,:,:,:), allocatable :: fieldU, fieldV, fieldA
        real(fpp), dimension(:,:,:), allocatable   :: fieldP
        real(fpp), dimension(:,:,:), allocatable   :: P_energy, S_energy, eps_vol
        real(fpp), dimension(:,:,:,:), allocatable :: eps_dev
        real(fpp), dimension(:,:,:,:), allocatable :: sig_dev
        !
        logical :: flag_gradU
        integer :: nx, ny, nz, i, j, k, ind

        flag_gradU = (out_variables(OUT_ENERGYP) + &
            out_variables(OUT_ENERGYS) + &
            out_variables(OUT_EPS_VOL) + &
            out_variables(OUT_EPS_DEV) + &
            out_variables(OUT_STRESS_DEV)) /= 0

        nx = el%ngllx
        ny = el%nglly
        nz = el%ngllz

        do k=0,nz-1
            do j=0,ny-1
                do i=0,nx-1
                    ind = el%Idom(i,j,k)

                    if (flag_gradU .or. (out_variables(OUT_DEPLA) == 1)) then
                        if(.not. allocated(fieldU)) allocate(fieldU(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldU(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_VITESSE) == 1) then
                        if(.not. allocated(fieldV)) allocate(fieldV(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldV(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_ACCEL) == 1) then
                        if(.not. allocated(fieldA)) allocate(fieldA(0:nx-1,0:ny-1,0:nz-1,0:2))
                        fieldA(i,j,k,:) = 0d0
                    end if

                    if (out_variables(OUT_PRESSION) == 1) then
                        if(.not. allocated(fieldP)) allocate(fieldP(0:nx-1,0:ny-1,0:nz-1))
                        fieldP(i,j,k) = 0d0
                    end if

                    if (out_variables(OUT_EPS_VOL) == 1) then
                        if(.not. allocated(eps_vol)) allocate(eps_vol(0:nx-1,0:ny-1,0:nz-1))
                        eps_vol(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYP) == 1) then
                        if(.not. allocated(P_energy)) allocate(P_energy(0:nx-1,0:ny-1,0:nz-1))
                        P_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_ENERGYS) == 1) then
                        if(.not. allocated(S_energy)) allocate(S_energy(0:nx-1,0:ny-1,0:nz-1))
                        S_energy(i,j,k) = 0.
                    end if

                    if (out_variables(OUT_EPS_DEV) == 1) then
                        if(.not. allocated(eps_dev)) allocate(eps_dev(0:nx-1,0:ny-1,0:nz-1,0:5))
                        eps_dev(i,j,k,:) = 0.
                    end if

                    if (out_variables(OUT_STRESS_DEV) == 1) then
                        if(.not. allocated(sig_dev)) allocate(sig_dev(0:nx-1,0:ny-1,0:nz-1,0:5))
                        sig_dev(i,j,k,:) = 0.
                    end if
                enddo
            enddo
        enddo
    end subroutine get_fluidpml_dom_var

  subroutine forces_int_flu_pml(dom, mat, champs1, Elem, lnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        type (subdomain), intent(IN) :: mat
        type(champsfluidpml), intent(inout) :: champs1
        type (Element), intent(INOUT) :: Elem
        integer :: lnum
        !
        integer :: m1, m2, m3
        integer :: i, j, k, l, ind
        real :: sum_vx, sum_vy, sum_vz, acoeff
        real, dimension(0:2,0:Elem%ngllx-1,0:Elem%nglly-1,0:Elem%ngllz-1)  :: ForcesFl

        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

        do k = 0,m3-1
            do j = 0,m2-1
                do i=0,m1-1
                    ind = Elem%Idom(i,j,k)
                    sum_vx = 0d0
                    sum_vy = 0d0
                    sum_vz = 0d0
                    do l = 0,m1-1
                        acoeff = - mat%hprimex(i,l)*mat%GLLwx(l)*mat%GLLwy(j)*mat%GLLwz(k)*Elem%Jacob(l,j,k)
                        sum_vx = sum_vx + acoeff*Elem%InvGrad(0,0,l,j,k)*dom%Veloc(l,j,k,0,lnum)
                        sum_vy = sum_vy + acoeff*Elem%InvGrad(1,0,l,j,k)*dom%Veloc(l,j,k,1,lnum)
                        sum_vz = sum_vz + acoeff*Elem%InvGrad(2,0,l,j,k)*dom%Veloc(l,j,k,2,lnum)
                    end do
                    ForcesFl(0,i,j,k) = sum_vx
                    ForcesFl(1,i,j,k) = sum_vy
                    ForcesFl(2,i,j,k) = sum_vz
                end do
            end do
        end do
        do k = 0,m3-1
            do l = 0,m2-1
                do j = 0,m2-1
                    do i=0,m1-1
                        ind = Elem%Idom(i,j,k)
                        acoeff = - mat%hprimey(j,l)*mat%GLLwx(i)*mat%GLLwy(l)*mat%GLLwz(k)*Elem%Jacob(i,l,k)
                        sum_vx = acoeff*Elem%InvGrad(0,1,i,l,k)*dom%Veloc(i,l,k,0,lnum)
                        sum_vy = acoeff*Elem%InvGrad(1,1,i,l,k)*dom%Veloc(i,l,k,1,lnum)
                        sum_vz = acoeff*Elem%InvGrad(2,1,i,l,k)*dom%Veloc(i,l,k,2,lnum)
                        ForcesFl(0,i,j,k) = ForcesFl(0,i,j,k) + sum_vx
                        ForcesFl(1,i,j,k) = ForcesFl(1,i,j,k) + sum_vy
                        ForcesFl(2,i,j,k) = ForcesFl(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do
        ! TODO reorder loops ?
        do l = 0,m3-1
            do k = 0,m3-1
                do j = 0,m2-1
                    do i=0,m1-1
                        ind = Elem%Idom(i,j,k)
                        acoeff = - mat%hprimez(k,l)*mat%GLLwx(i)*mat%GLLwy(j)*mat%GLLwz(l)*Elem%Jacob(i,j,l)
                        sum_vx = acoeff*Elem%InvGrad(0,2,i,j,l)*dom%Veloc(i,j,l,0,lnum)
                        sum_vy = acoeff*Elem%InvGrad(1,2,i,j,l)*dom%Veloc(i,j,l,1,lnum)
                        sum_vz = acoeff*Elem%InvGrad(2,2,i,j,l)*dom%Veloc(i,j,l,2,lnum)
                        ForcesFl(0,i,j,k) = ForcesFl(0,i,j,k) + sum_vx
                        ForcesFl(1,i,j,k) = ForcesFl(1,i,j,k) + sum_vy
                        ForcesFl(2,i,j,k) = ForcesFl(2,i,j,k) + sum_vz
                    end do
                end do
            end do
        end do
        
        ! Assemblage
        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    ind = Elem%Idom(i,j,k)
                    ! We should have atomic adds with openmp // here
                    champs1%fpml_Forces(ind,0) = champs1%fpml_Forces(ind,0) + ForcesFl(0,i,j,k)
                    champs1%fpml_Forces(ind,1) = champs1%fpml_Forces(ind,1) + ForcesFl(1,i,j,k)
                    champs1%fpml_Forces(ind,2) = champs1%fpml_Forces(ind,2) + ForcesFl(2,i,j,k)
                enddo
            enddo
        enddo

        return
    end subroutine forces_int_flu_pml

    subroutine pred_flu_pml(dom, mat, dt, champs1, Elem, lnum)
        type (domain_fluidpml), intent (INOUT) :: dom
        type (subdomain), intent(IN) :: mat
        real, intent(in) :: dt
        type(champsfluidpml), intent(inout) :: champs1
        type(Element), intent(inout) :: Elem
        integer :: lnum
        !
        real, dimension(0:Elem%ngllx-1, 0:Elem%nglly-1, 0:Elem%ngllz-1) :: dVelPhi_dx, dVelPhi_dy, dVelPhi_dz
        integer :: m1, m2, m3
        integer :: i, j, k, ind
        real, dimension(:,:,:), allocatable :: VelPhi

        m1 = Elem%ngllx ; m2 = Elem%nglly ; m3 = Elem%ngllz

        allocate(VelPhi(0:m1-1,0:m2-1,0:m3-1))
        ! prediction in the element
        ! We do the sum V1+V2+V3 and F1+F2+F3 here
        do k = 0,m3-1
            do j = 0,m2-1
                do i = 0,m1-1
                    ind = Elem%Idom(i,j,k)
                    VelPhi(i,j,k) = champs1%fpml_VelPhi(ind,0) + &
                        champs1%fpml_VelPhi(ind,1) + &
                        champs1%fpml_VelPhi(ind,2)
                enddo
            enddo
        enddo
        ! XXX DumpS{xyz}(:,:,:,1) doit etre multiplie par 1/density
        ! d(rho*Phi)_d(xi,eta,zeta)
        call physical_part_deriv(m1,m2,m3,mat%htprimex,mat%hprimey,mat%hprimez,Elem%InvGrad, &
            VelPhi(:,:,:), dVelPhi_dx, dVelPhi_dy, dVelPhi_dz)

        ! prediction for (physical) velocity (which is the equivalent of a stress, here)
        ! V_x^x
        dom%Veloc(:,:,:,0,lnum) = dom%PMLDumpSx(:,:,:,0,lnum) * dom%Veloc(:,:,:,0,lnum) + &
                                  dom%PMLDumpSx(:,:,:,1,lnum) * Dt * dVelPhi_dx
        ! V_x^y = 0
        ! V_x^z = 0
        ! V_y^x = 0
        ! V_y^y
        dom%Veloc(:,:,:,1,lnum) = dom%PMLDumpSy(:,:,:,0,lnum) * dom%Veloc(:,:,:,1,lnum) + &
                                  dom%PMLDumpSy(:,:,:,1,lnum) * Dt * dVelPhi_dy
        ! V_y^z = 0
        ! V_z^x = 0
        ! V_z^y = 0
        ! V_z^z
        dom%Veloc(:,:,:,2,lnum) = dom%PMLDumpSz(:,:,:,0,lnum) * dom%Veloc(:,:,:,2,lnum) + &
                                  dom%PMLDumpSz(:,:,:,1,lnum) * Dt * dVelPhi_dz
        return
    end subroutine Pred_Flu_Pml

end module dom_fluidpml

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
