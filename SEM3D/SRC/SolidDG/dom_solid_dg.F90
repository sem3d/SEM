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
        integer :: ngll, nf

        ngll = dom%ngll
        nf   = dom%nbelem

        ! Q
        allocate(dom%champs(i)%Q(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:ngll-1,0:8,0:dom%nblocks-1))

        ! Traces (last face==dummy)
        allocate(dom%champs(i)%trace_Q(0:VCHUNK-1,0:ngll-1,0:ngll-1,0:8,0:nf/VCHUNK))

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
            allocate(dom%Qasm    (0:dom%nglltot,0:8))
            allocate(dom%Density_(0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Lambda_ (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            allocate(dom%Mu_     (0:ngll-1, 0:ngll-1, 0:ngll-1,0:nblocks-1, 0:VCHUNK-1))
            ! 0:9 0: face, 1-3 : I0 4-6: DI 7-9:DJ  face(i,j) = elem(I0+i*DI+j*DJ)m 10: side
            allocate(dom%Itrace_ (0:5,0:10,0:nblocks-1, 0:VCHUNK-1))
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
                dom%Itrace_(nf, 0,   ec, eb) = Tdomain%sFace(nnf)%lnum
                dom%Itrace_(nf, 1:3, ec, eb) = i0
                dom%Itrace_(nf, 4:6, ec, eb) = di
                dom%Itrace_(nf, 7:9, ec, eb) = dj
                dom%Itrace_(nf, 10,  ec, eb) = side

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
        if(allocated(dom%Qasm     )) deallocate(dom%Qasm     )

        do i=0,size(dom%champs)-1
            if(allocated(dom%champs(i)%Q       )) deallocate(dom%champs(i)%Q       )
            if(allocated(dom%champs(i)%trace_Q )) deallocate(dom%champs(i)%trace_Q )
        end do

        call deallocate_dombase(dom)

    end subroutine deallocate_dom_solid_dg


    subroutine get_solid_dg_dom_var(dom, lnum, out_variables,  &
         fieldV, fieldA, eps_vol, eps_dev)

        use deriv3d
    
        implicit none

        type(domain_solid_dg), intent(inout) :: dom
        integer, intent(in), dimension(0:)   :: out_variables
        integer, intent(in)                  :: lnum
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldV
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:2) :: fieldA
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1) :: eps_vol
        real(fpp), dimension(0:dom%ngll-1,0:dom%ngll-1,0:dom%ngll-1,0:5) :: eps_dev
        integer :: ngll, i, j, k
        integer :: bnum, ee
        real(fpp) :: exx, eyy, ezz

        bnum = lnum/VCHUNK

        ee = mod(lnum,VCHUNK)

        ngll = dom%ngll
        
        do k = 0,ngll-1
            do j = 0,ngll-1
                do i = 0,ngll-1
                    if (out_variables(OUT_VITESSE) == 1) fieldV(i,j,k,:) = dom%champs(0)%Q(ee,i,j,k,6:8,bnum)
                    if (out_variables(OUT_ACCEL) == 1) fieldA(i,j,k,:) = dom%champs(2)%Q(ee,i,j,k,6:8,bnum)
                    if (out_variables(OUT_EPS_VOL) == 1) then
                       exx = dom%champs(0)%Q(ee,i,j,k,0,bnum)
                       eyy = dom%champs(0)%Q(ee,i,j,k,1,bnum)
                       ezz = dom%champs(0)%Q(ee,i,j,k,2,bnum)
                       eps_vol(i,j,k) = exx+eyy+ezz
                    end if
                    if (out_variables(OUT_EPS_DEV) == 1) eps_dev(i,j,k,:) = dom%champs(0)%Q(ee,i,j,k,0:5,bnum)
                enddo
            enddo
        enddo
    
    end subroutine get_solid_dg_dom_var


    subroutine get_solid_dg_dom_elem_energy(dom, lnum, P_energy, S_energy, R_energy, C_energy)

        use deriv3d

        implicit none

        type(domain_solid_dg), intent(inout) :: dom
        integer, intent(in)                  :: lnum
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
                    fieldU(i,j,k,:) = 0 !dom%champs(0)%Depla(ind,:)
                    fieldV(i,j,k,:) = dom%Qasm(ind,:)
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
        integer :: bnum, ee

        bnum = lnum/VCHUNK
        ee   = mod(lnum,VCHUNK)

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
        integer :: bnum, ee

        bnum = specel%lnum/VCHUNK
        ee   = mod(specel%lnum,VCHUNK)

        specel%MassMat(i,j,k) = Whei*dom%Density_(i,j,k,bnum,ee)*dom%Jacob_(i,j,k,bnum,ee)
        dom%MassMat(ind)      = dom%MassMat(ind) + specel%MassMat(i,j,k)

    end subroutine init_local_mass_solid_dg

    subroutine forces_int_solid_dg(dom, var, dvdt, bnum)

        use m_calcul_forces_dg

        type(domain_solid_dg), intent (INOUT) :: dom
        type(champssolid_dg),  intent(IN)     :: var
        type(champssolid_dg),  intent(INOUT)  :: dvdt
        integer, intent(in) :: bnum
        
        integer :: ngll

        ngll = dom%ngll

        call calcul_forces_solid_dg_iso(dom,bnum,var%Q(:,:,:,:,:,bnum), dvdt%Q(:,:,:,:,:,bnum))

    end subroutine forces_int_solid_dg

    function solid_Pspeed_dg(dom, lnum, i, j, k) result(Pspeed)
        
        type(domain_solid_dg), intent (IN) :: dom
        integer, intent(in) :: lnum, i, j, k
        real(fpp) :: Pspeed, M
        integer :: bnum, ee
        
        bnum   = lnum/VCHUNK
        ee     = mod(lnum,VCHUNK)
        M      = dom%Lambda_(i,j,k,bnum,ee) + 2.*dom%Mu_(i,j,k,bnum,ee)
        Pspeed = sqrt(M/dom%Density_(i,j,k,bnum,ee))
    
    end function solid_Pspeed_dg

    subroutine lddrk_init_solid_dg(dom, f2)

        type(domain_solid_dg), intent (INOUT) :: dom
        integer, intent(in) :: f2
    
        if (dom%nglltot == 0) return

        !dom%champs(f2)%Q(:,:,:,:,:,:) = 0.d0
        !dom%champs(0)%Q(:,:,:,:,:,:) = 0.d0
        dom%champs(0)%Q(:,0,:,:,6,:) = 1.d0
        
      end subroutine lddrk_init_solid_dg

    subroutine lddrk_update_solid_dg(dom, f0, f1, f2, dt, cb, cg)

        type(domain_solid_dg), intent (INOUT) :: dom
        real(fpp), intent(in) :: cb, cg, dt
        integer :: n, i, j, k, ee, ngll, p
        integer :: f0, f1, f2
        real(fpp), dimension(0:8) :: val

        ngll = dom%ngll

        ! f2  contains forces computation  ie f2 = dU/dt
        ! f1 : w(n+1) = cb*w(n) + dt*dU/dt
        ! f0 : U(n+1) = U(n) + cg*w(n+1)
        if (dom%nglltot==0) return
        ! Only solid for starters
        do n = 0,dom%nblocks-1
           do k = 0,ngll-1
              do j = 0,ngll-1
                 do i = 0,ngll-1
                    do ee = 0, VCHUNK-1
                       val = cb*dom%champs(f1)%Q(ee,i,j,k,:,n) + dt*dom%champs(f2)%Q(ee,i,j,k,:,n)
                       dom%champs(f1)%Q(ee,i,j,k,:,n) = val
                       dom%champs(f0)%Q(ee,i,j,k,:,n) = dom%champs(f0)%Q(ee,i,j,k,:,n) + cg*val
                    enddo
                 enddo
              enddo
           enddo
        enddo

    end subroutine lddrk_update_solid_dg


end module dom_solid_dg
