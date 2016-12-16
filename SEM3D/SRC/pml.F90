!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module pml
    use constants
    use champs_solidpml
    use champs_fluidpml
    implicit none

contains

    subroutine define_alpha_PML(Coord, dir, ngll, vp, pml_width, pml_pos, Apow, npow, alpha)
        !- routine determines attenuation profile in an PML layer (see Festa & Vilotte)
        !   dir = attenuation's direction, ldir_attenu = the logical giving the orientation
        integer, intent(in) :: dir, ngll, npow
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1, 0:2), intent(in) :: Coord
        real(fpp), intent(in)  :: Apow
        real(fpp), dimension(0:2), intent(in) :: pml_pos, pml_width
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: vp
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: alpha
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1)  :: ri
        real(fpp) :: invdh, coef
        if (pml_width(dir)==0d0) then
            alpha(:,:,:) = 0d0
            return
        end if

        ! Until we can do better, the plane equation is only parallel to one axis
        ! when more coefficient are != 0 it means we have a corner
        invdh = 1d0/abs(pml_width(dir))
        coef = 1/pml_width(dir)
        ri = coef*(Coord(:,:,:,dir)-pml_pos(dir))
        alpha = Apow * Vp * invdh *  (ri)**npow
    end subroutine define_alpha_PML

    subroutine define_PML_DumpInit(ngll,dt,alpha,&
        MassMat,DumpS,DumpMass)
        !- defining parameters related to stresses and mass matrix elements, in the case of
        !    a PML, along a given splitted direction:
        integer, intent(in)  :: ngll
        real(fpp), intent(in) :: dt
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: alpha
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(in) :: MassMat
        real(fpp), dimension(0:1,0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: DumpS
        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1), intent(out) :: DumpMass

        real(fpp), dimension(0:ngll-1,0:ngll-1,0:ngll-1)  :: Id

        Id = 1d0

        DumpS(1,:,:,:) = Id + 0.5d0*dt*alpha
        DumpS(1,:,:,:) = 1d0/DumpS(1,:,:,:)
        DumpS(0,:,:,:) = (Id - 0.5d0*dt*alpha)*DumpS(1,:,:,:)
        DumpMass(:,:,:) = 0.5d0*MassMat(:,:,:)*alpha(:,:,:)*dt

        return
    end subroutine define_PML_DumpInit

    subroutine define_PML_DumpEnd(nglltot,Massmat,DumpMass,DumpV)
        implicit none
        integer, intent(in)   :: nglltot
        real(fpp), dimension(0:nglltot), intent(in) :: MassMat
        real(fpp), dimension(0:nglltot,0:2), intent(in) :: DumpMass
        real(fpp), dimension(0:nglltot,0:1,0:2), intent(out) :: DumpV

        DumpV(:,1,:) = spread(MassMat(:),2,3) + DumpMass(:,:)
        DumpV(:,1,:) = 1d0/DumpV(:,1,:)
        DumpV(:,0,:) = spread(MassMat(:),2,3) - DumpMass(:,:)
        DumpV(:,0,:) = DumpV(:,0,:) * DumpV(:,1,:)

        return
    end subroutine define_PML_DumpEnd
    !
    subroutine cpml_compute_coefs(m, a0, dt, cf0, cf1, cf2)
        integer, intent(in) :: m
        real(fpp), intent(in) :: a0, dt
        real(fpp), intent(out) :: cf0, cf1, cf2
        select case (m)
        case (CPML_MIDPOINT)
            call cpml_compute_coefs_midpoint_ctr(a0,dt,cf0,cf1,cf2)
        case (CPML_ORDER1)
            call  cpml_compute_coefs_O1(a0,dt,cf0,cf1,cf2)
        case (CPML_ORDER2)
            call cpml_compute_coefs_O2(a0,dt,cf0,cf1,cf2)
        case default
            stop "Unknown CPML integration scheme"
        end select
    end subroutine cpml_compute_coefs
    !
    subroutine cpml_compute_coefs_midpoint(a0, dt, cf0, cf1, cf2)
        real(fpp), intent(in) :: a0, dt
        real(fpp), intent(out) :: cf0, cf1, cf2
        !
        real(fpp) :: c
        ! Update convolution term (implicit midpoint)
        c = (1d0+0.5d0*a0*dt)
        cf0 = (1d0-0.5d0*a0*dt)/c
        cf1 = dt/c
        cf2 = 0d0
    end subroutine cpml_compute_coefs_midpoint
    !
    subroutine cpml_compute_coefs_midpoint_ctr(a0, dt, cf0, cf1, cf2)
        real(fpp), intent(in) :: a0, dt
        real(fpp), intent(out) :: cf0, cf1, cf2
        !
        real(fpp) :: c
        real(fpp), parameter :: theta=0.25d0
        ! Update convolution term (implicit midpoint)
        c = (1d0+0.5d0*a0*dt)
        cf0 = (1d0-0.5d0*a0*dt)/c
        cf1 = theta*dt/c
        cf2 = (1d0-theta)*dt/c
    end subroutine cpml_compute_coefs_midpoint_ctr
    !
    subroutine cpml_compute_coefs_O1(a0, dt, cf0, cf1, cf2)
        real(fpp), intent(in) :: a0, dt
        real(fpp), intent(out) :: cf0, cf1, cf2
        ! First order CPML
        cf0 = exp(-a0*dt)
        if (abs(a0)<1e-8) then
            cf1 = dt
        else
            cf1 = (1d0-cf0)/a0
        endif
        cf2 = 0d0
    end subroutine cpml_compute_coefs_O1
    !
    subroutine cpml_compute_coefs_O2(a0, dt, cf0, cf1, cf2)
        real(fpp), intent(in) :: a0, dt
        real(fpp), intent(out) :: cf0, cf1, cf2
        ! First order CPML
        cf0 = exp(-a0*dt)
        if (abs(a0)<1e-8) then
            cf1 = 0.5d0*dt
            cf2 = cf1*exp(-0.5d0*a0*dt)
        else
            cf1 = (1d0-exp(-0.5d0*a0*dt))/a0
            cf2 = cf1*exp(-0.5d0*a0*dt)
        endif
    end subroutine cpml_compute_coefs_O2
    !
    subroutine copy_cpml_coordinates(Tdomain, dom, dmtype)
        use sdomain
        type(domain), intent (INOUT) :: Tdomain
        class(dombase_cpml), intent(inout) :: dom
        integer, intent(in) :: dmtype
        integer :: i, j, k, n, indL, indG

        ! Handle on node global coords : mandatory to compute
        ! distances in the PML (compute_dxi_alpha_kappa)
        ! TODO precompute usefull coeffs instead of copying coords...
        allocate(dom%GlobCoord(0:2,0:dom%nglltot-1))
        do n=0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%domain/=dmtype) cycle
            do k = 0,dom%ngll-1
                do j = 0,dom%ngll-1
                    do i = 0,dom%ngll-1
                        indG = Tdomain%specel(n)%Iglobnum(i,j,k)
                        indL = Tdomain%specel(n)%Idom(i,j,k)
                        dom%GlobCoord(:,indL) = Tdomain%GlobCoord(:,indG)
                    end do
                end do
            end do
        end do
    end subroutine copy_cpml_coordinates

    subroutine cpml_allocate_multi_dir(Tdomain, dom, dmtype)
        use sdomain
        use gll3d
        implicit none
        type(domain) :: TDomain
        class(dombase_cpml), intent (INOUT) :: dom
        integer, intent(in) :: dmtype
        !
        integer :: dir1_count, dir2_count, ndir, mi, n, ngll
        dir1_count = 0
        dir2_count = 0
        ngll = dom%ngll
        do n=0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%domain/=dmtype) cycle
            ndir = 0
            mi = Tdomain%specel(n)%mat_index
            if (Tdomain%sSubDomain(mi)%pml_width(0)/=0d0) ndir = ndir + 1
            if (Tdomain%sSubDomain(mi)%pml_width(1)/=0d0) ndir = ndir + 1
            if (Tdomain%sSubDomain(mi)%pml_width(2)/=0d0) ndir = ndir + 1
            if (ndir>=2) dir1_count = dir1_count + 1
            if (ndir>=3) dir2_count = dir2_count + 1
        end do
        write(*,*) "Allocating: PML-N1:", dir1_count
        write(*,*) "Allocating: PML-N2:", dir2_count
        if (dir1_count>0) then
            allocate(dom%Alpha_1(0:ngll-1,0:ngll-1,0:ngll-1,0:dir1_count-1))
            allocate(dom%Kappa_1(0:ngll-1,0:ngll-1,0:ngll-1,0:dir1_count-1))
            allocate(dom%dxi_k_1(0:ngll-1,0:ngll-1,0:ngll-1,0:dir1_count-1))

        endif
        if (dir2_count>0) then
            allocate(dom%Alpha_2(0:ngll-1,0:ngll-1,0:ngll-1,0:dir2_count-1))
            allocate(dom%Kappa_2(0:ngll-1,0:ngll-1,0:ngll-1,0:dir2_count-1))
            allocate(dom%dxi_k_2(0:ngll-1,0:ngll-1,0:ngll-1,0:dir2_count-1))
        end if
        dom%dir1_count = dir1_count
        dom%dir2_count = dir2_count
    end subroutine cpml_allocate_multi_dir

    subroutine compute_dxi_alpha_kappa(dom, xi, L, wpml, Pspeed, alpha, kappa, dxi)
        type(domain_fluidpml), intent(inout) :: dom
        real(fpp), intent(in) :: xi, L, wpml, Pspeed
        real(fpp), intent(out) :: alpha, kappa, dxi
        !
        real(fpp) :: xoverl, d0
        integer :: lnum

        xoverl = xi/wpml
        if (xoverl > 1d0) xoverl = 1d0
        if (xoverl < 0d0) xoverl = 0d0

        d0 = -1.*(dom%cpml_n+1)*Pspeed*log(dom%cpml_rc)
        d0 = d0/(2*wpml)

        kappa = dom%cpml_kappa_0 + dom%cpml_kappa_1 * xoverl
        dxi   = dom%cpml_c*d0*(xoverl)**dom%cpml_n / kappa
        alpha = dom%alphamax*(1. - xoverl) ! alpha*: (76) from Ref1
    end subroutine compute_dxi_alpha_kappa

end module pml
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
