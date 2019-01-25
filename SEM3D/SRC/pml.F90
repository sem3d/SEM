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
#ifdef CPML
    !
    subroutine cpml_coefs_midpoint2(a0, dt, cf0, cf1)
        real(fpp), intent(in) :: a0, dt
        real(fpp), intent(out) :: cf0, cf1
        !
        real(fpp) :: c
        ! Update convolution term (implicit midpoint)
        c = (1d0+0.5d0*a0*dt)
        cf1 = dt/c
        cf0 = -a0*cf1
    end subroutine cpml_coefs_midpoint2

    !
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
        allocate(dom%GlobCoord(0:2,0:dom%nglltot))
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

    subroutine cpml_reorder_elements(Tdomain, dom, dmtype)
        use sdomain
        use gll3d
        implicit none
        type(domain) :: TDomain
        class(dombase_cpml), intent (INOUT) :: dom
        integer, intent(in) :: dmtype
        !
        integer, allocatable, dimension(:) :: types
        integer :: dir1_count, dir2_count, ndir, mi, n, dir
        !integer :: itype, num, lnum
        dir1_count = 0
        dir2_count = 0
        !
        allocate(types(0:Tdomain%n_elem-1))
        !
        ! We sort elements according to PML attenuation direction, and associate a type :
        ! Type=0 PML attenuation in X only
        ! Type=1 PML attenuation in Y only
        ! Type=2 PML attenuation in Z only
        ! Type=3 all other
        ! We also count the number of elements with 2 or 3 attenuating directions
        types = -1
        do n=0,Tdomain%n_elem-1
            if (Tdomain%specel(n)%domain/=dmtype) cycle
            ndir = 0
            dir = -1
            mi = Tdomain%specel(n)%mat_index
            if (Tdomain%sSubDomain(mi)%pml_width(0)/=0d0) then
                ndir = ndir + 1
                dir = 0
            endif
            if (Tdomain%sSubDomain(mi)%pml_width(1)/=0d0) then
                ndir = ndir + 1
                dir = 1
            endif
            if (Tdomain%sSubDomain(mi)%pml_width(2)/=0d0) then
                ndir = ndir + 1
                dir = 2
            endif
            if (ndir==1) then
                types(n) = dir
            else
                types(n) = 3
            end if
            if (ndir>=2) dir1_count = dir1_count + 1
            if (ndir>=3) dir2_count = dir2_count + 1
        end do
        ! Now renumber elements according to their type
!        num = 0
!        do itype=0,3
!            ! Skip num to the next VCHUNK
!            num = VCHUNK*((num+VCHUNK-1)/VCHUNK)
!            do n=0,Tdomain%n_elem-1
!                if (types(n)/=itype) cycle
!                write(*,*) Tdomain%specel(n)%lnum, "->", num
!                Tdomain%specel(n)%lnum = num
!                num = num+1
!            end do
!        end do
!        if (dom%nbelem /= num) stop 1
!        ! Cheat dom%nelem so we allocate all chunks correctly
!        dom%nbelem = num
        dom%dir1_count = dir1_count
        dom%dir2_count = dir2_count
        deallocate(types)
    end subroutine cpml_reorder_elements

    subroutine cpml_allocate_multi_dir(Tdomain, dom)
        use sdomain
        use gll3d
        implicit none
        type(domain) :: TDomain
        class(dombase_cpml), intent (INOUT) :: dom
        !
        integer :: dir1_count, dir2_count, ngll
        ngll = dom%ngll
        dir1_count = dom%dir1_count
        dir2_count = dom%dir2_count
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
    end subroutine cpml_allocate_multi_dir

    subroutine compute_dxi_alpha_kappa(dom, xi, L, wpml, Pspeed, alpha, kappa, dxi)
        type(domain_fluidpml), intent(inout) :: dom
        real(fpp), intent(in) :: xi, L, wpml, Pspeed
        real(fpp), intent(out) :: alpha, kappa, dxi
        !
        real(fpp) :: xoverl, d0

        xoverl = xi/wpml
        if (xoverl > 1d0) xoverl = 1d0
        if (xoverl < 0d0) xoverl = 0d0

        d0 = -1.*(dom%cpml_n+1)*Pspeed*log(dom%cpml_rc)
        d0 = d0/(2*wpml)

        kappa = dom%cpml_kappa_0 + dom%cpml_kappa_1 * xoverl
        dxi   = dom%cpml_c*d0*(xoverl)**dom%cpml_n / kappa
        alpha = dom%alphamax*(1. - xoverl) ! alpha*: (76) from Ref1
    end subroutine compute_dxi_alpha_kappa
#endif
    function isclose(a,b)
        real(fpp) :: a,b
        logical :: isclose
        isclose = abs(a-b)<SMALLFPP
    end function isclose

    subroutine get_coefs_L_abc(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: a3b, a4b, a5b
        !
        real(fpp) :: k
        k = k0*k1*k2
        a3b = k*a0*a0*d0*(d1+a1-a0)*(d2+a2-a0)/((a1-a0)*(a2-a0))
        a4b = k*a1*a1*d1*(d0+a0-a1)*(d2+a2-a1)/((a0-a1)*(a2-a1))
        a5b = k*a2*a2*d2*(d0+a0-a2)*(d1+a1-a2)/((a0-a2)*(a1-a2))
    end subroutine get_coefs_L_abc

    subroutine get_coefs_L_aba(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: a3b, a4b, a5b
        !
        real(fpp) :: k
        k = k0*k1*k2
        a3b = k*a0*(a0**3*d0 + a0**3*d2 - 2*a0**2*a1*d0 - 2*a0**2*a1*d2 - a0**2*d0*d1 &
            - 2*a0**2*d0*d2 - a0**2*d1*d2 + a0*a1**2*d0 + a0*a1**2*d2 + a0*a1*d0*d1 &
            + 4*a0*a1*d0*d2 + a0*a1*d1*d2 + a0*d0*d1*d2 - 2*a1**2 *d0*d2 - 2*a1*d0*d1*d2)/(a0 - a1)**2
        a4b = k*a1**2*d1*( a0 - a1 + d0)*(a0 - a1 + d2)/(a0 - a1)**2
        a5b = k*a0**2*d0*d2*(a0 - a1 - d1)/(a0 - a1)
    end subroutine get_coefs_L_aba

    subroutine get_coefs_L_aac(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: a3b, a4b, a5b
        !
        real(fpp) :: k
        k = k0*k1*k2
        a3b = k*a0*(a0**3*d0 + a0**3*d1 - 2*a0**2*a2*d0 - 2*a0**2*a2*d1 - 2*a0**2 *d0*d1 &
            - a0**2*d0*d2 - a0**2*d1*d2 + a0*a2**2*d0 + a0*a2**2*d1 + 4*a0*a2*d0*d1 + a0*a2*d0*d2 &
            + a0*a2*d1*d2 + a0*d0*d1*d2 - 2*a2**2 *d0*d1 - 2*a2*d0*d1*d2)/(a0 - a2)**2
        a4b = k*a0**2*d0*d1*(a0 - a2 - d2)/(a0 - a2)
        a5b = k*a2**2*d2*( a0 - a2 + d0)*(a0 - a2 + d1)/(a0 - a2)**2
    end subroutine get_coefs_L_aac

    subroutine get_coefs_L_abb(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: a3b, a4b, a5b
        !
        real(fpp) :: k
        k = k0*k1*k2
        a3b = k*a0**2*d0*(a0 - a1 - d1)*(a0 - a1 - d2)/(a0 - a1)**2
        a4b = k* a1*(a0**2*a1*d1 + a0**2*a1*d2 - 2*a0**2*d1*d2 - 2*a0*a1**2*d1 &
            - 2*a0*a1**2*d2 + a0*a1*d0*d1 + a0*a1*d0*d2 + 4*a0*a1*d1*d2 &
            - 2*a0*d0*d1*d2 + a1**3*d1 + a1**3*d2 - a1**2*d0*d1 - a1**2 *d0*d2 &
            - 2*a1**2*d1*d2 + a1*d0*d1*d2)/(a0 - a1)**2
        a5b = k*a1**2*d1*d2*(a0 - a1 + d0)/(a0 - a1)
    end subroutine get_coefs_L_abb

    subroutine get_coefs_L_aaa(k0,k1,k2,a0,a1,a2,d0,d1,d2, a3b, a4b, a5b)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: a3b, a4b, a5b
        !
        real(fpp) :: k
        k = k0*k1*k2
        a3b = k*(a0**2*d0 + a0**2*d1 + a0**2*d2 - 2*a0*d0*d1 - 2*a0*d0*d2 - 2*a0*d1*d2 + d0*d1*d2)
        a4b = k*a0*(a0*d0*d1 + a0*d0*d2 + a0*d1*d2 - 2*d0*d1*d2)
        a5b = k*a0**2*d0*d1*d2
    end subroutine get_coefs_L_aaa

    subroutine update_roots_L1_3d(sel, c, b1, b2, b3)
        integer, intent(in) :: sel, c
        real(fpp), intent(inout), dimension(0:14) :: b1, b2, b3
        !
        select case (sel)
        case(CMP_ABC)
        case(CMP_AAC)
            b1(c) = b2(c); b2(c) = 0.
        case(CMP_ABA)
        case(CMP_ABB)
        case(CMP_AAA)
            b1(c) = b2(c); b2(c) = 0.
        end select
    end subroutine update_roots_L1_3d

    subroutine update_roots_L2_3d(sel, c, b1, b2, b3)
        integer, intent(in) :: sel, c
        real(fpp), intent(inout), dimension(0:14) :: b1, b2, b3
        !
        select case (sel)
        case(CMP_ABC)
        case(CMP_AAC)
        case(CMP_ABA)
            b1(c) = b3(c); b3(c) = 0.
        case(CMP_ABB)
            b2(c) = b3(c); b3(c) = 0.
        case(CMP_AAA)
            b1(c) = b3(c); b3(c) = 0.
        end select
    end subroutine update_roots_L2_3d

    function compare_roots(a,b,c) result(sel)
        real(fpp), intent(in) :: a,b,c
        integer :: sel
        sel = CMP_ABC
        if (isclose(a,b)) sel = CMP_AAC
        if (isclose(a,c)) then
            if (sel==0) then
                sel = CMP_ABA
            else
                sel = CMP_AAA
            end if
        end if
        if (isclose(b,c)) then
            if (sel==0) then
                sel = CMP_ABB
            else
                sel = CMP_AAA
            end if
        end if
    end function compare_roots

    subroutine get_coefs_Lijk_2d(dom, ee, bnum, i, j, k, kB120, kB021, kB012, b0, b1, b2, e0, e1)
        class(dombase_cpml), intent (INOUT) :: dom
        integer, intent(in) :: ee, bnum, i, j, k, kB120, kB021, kB012
        real(fpp), dimension(0:), intent (OUT) :: b0, b1, b2 ! CAUTION : dimension is 0:2 for FluidCPML, and, 0:5 for SolidCPML
        real(fpp), dimension(0:), intent (OUT) :: e0, e1 ! CAUTION : dimension is 0:2 for FluidCPML, and, 0:8 for SolidCPML
        !
        real(fpp) :: k0, d0, a0
        real(fpp) :: k1, d1, a1
        integer   :: i1, dim0, dim1, Is0xs1, Is0os1, Is1os0
        real(fpp) :: AA, BB, CC, DD ! common subexpressions

        k0 = dom%Kappa_0(ee,i,j,k,bnum)
        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        i1 = dom%I1(ee,bnum)
        k1 = dom%Kappa_1(i,j,k,i1)
        a1 = dom%Alpha_1(i,j,k,i1)
        d1 = dom%dxi_k_1(i,j,k,i1)
        if (k0<0.or.k1<0) then
            write(*,*) "Erreur:", bnum, ee, i,j,k, k0
            stop 1
        endif

        dim0 = dom%D0(ee,bnum)
        dim1 = dom%D1(ee,bnum)
        ! The 3 terms L120, L021, L012 are one of
        ! F-1(s0.s1) F-1(s0/s1) and F-1(s1/s0) depending
        ! on which coordinates are in slot 0/1
        if (dim0==0.and.dim1==1) then
            ! For X=0,Y=1
            Is1os0 = kB120
            Is0os1 = kB021
            Is0xs1 = kB012
        else if (dim0==0.and.dim1==2) then
            ! For X=0,Z=1
            Is1os0 = kB120
            Is0xs1 = kB021
            Is0os1 = kB012
        else if (dim0==1.and.dim1==2) then
            ! For Y=0,Z=1
            Is0xs1 = kB120
            Is1os0 = kB021
            Is0os1 = kB012
        else
            stop 1
        end if

        b0(Is1os0) = k1/k0
        b0(Is0os1) = k0/k1
        b0(Is0xs1) = k0*k1

        AA = (a0+d0 - (a1+d1))
        BB = (a0+d0 -  a1)
        CC = (a0-a1)
        DD = (a0 - (a1+d1))

        e0(Is1os0) = a0+d0
        e1(Is1os0) = a1
        if (.not. isclose(e0(Is1os0),e1(Is1os0))) then
            b1(Is1os0) = -b0(Is1os0)*d0*AA/BB !b3_120
            b2(Is1os0) =  b0(Is1os0)*d1*CC/BB !b1_120
        else
            b1(Is1os0) = b0(Is1os0)*(d1+a0-a1)
            b2(Is1os0) = b0(Is1os0)*(a0-a1)*d1
        end if

        e0(Is0os1) = a0
        e1(Is0os1) = a1+d1
        if (.not. isclose(e0(Is0os1),e1(Is0os1))) then
            b1(Is0os1) = b0(Is0os1)*d0*CC/DD  !b1_021
            b2(Is0os1) =-b0(Is0os1)*d1*AA/DD  !b3_021
        else
            b1(Is0os1) = b0(Is0os1)*(d0+a1-a0)
            b2(Is0os1) = b0(Is0os1)*(a1-a0)*d0
        end if
        e0(Is0xs1) = a0
        e1(Is0xs1) = a1
        if (.not. isclose(e0(Is0xs1),e1(Is0xs1))) then
            b1(Is0xs1) = b0(Is0xs1)*d0*DD/CC  !b1_012
            b2(Is0xs1) = b0(Is0xs1)*d1*BB/CC  !b2_012
        else
            b1(Is0xs1) = b0(Is0xs1)*(d0+d1)   !b1_012
            b2(Is0xs1) = b0(Is0xs1)*d1*d0     !b2_012
        end if
    end subroutine get_coefs_Lijk_2d

    subroutine get_coefs_Lijk_3d(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3, sel)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: cb0, cb1, cb2, cb3
        integer  , intent(out) :: sel
        !
        real(fpp) :: b2
        !
        b2 = a2+d2
        sel = compare_roots(a0,a1,b2)
        select case(sel)
        case (CMP_ABC)
            call get_coefs_Lijk_abc(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        case (CMP_AAC)
            call get_coefs_Lijk_aac(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        case (CMP_ABA)
            call get_coefs_Lijk_aba(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        case (CMP_ABB)
            call get_coefs_Lijk_abb(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        case (CMP_AAA)
            call get_coefs_Lijk_aaa(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        case default
            stop 1
        end select
    end subroutine get_coefs_Lijk_3d

    subroutine get_coefs_Lijk_abc(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: cb0, cb1, cb2, cb3
        !
        cb0 =  k0*k1/k2
        cb1 =  cb0*(a0-a2)*d0*(a0-a1-d1)/((a0-a1)*(a0-a2-d2))
        cb2 =  cb0*(a1-a2)*(a1-a0-d0)*d1/((a1-a0)*(a1-a2-d2))
        cb3 = -cb0*d2*(a2-a0+d2-d0)*(a2-a1+d2-d1)/((a2+d2-a0)*(a2+d2-a1))
    end subroutine get_coefs_Lijk_abc

    subroutine get_coefs_Lijk_aac(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: cb0, cb1, cb2, cb3
        !
        cb0 =  k0*k1/k2

        cb2 = cb0*d0*d1*(a0 - a2)/(a0 - a2 - d2)
        cb3 = - cb0*d2*(a0 - a2 + d0 - d2)*(a0 - a2 + d1 - d2)/(a0 - a2 - d2)**2
        cb1 = cb0*(a0**2*d0 + a0**2*d1 - 2*a0*a2*d0 - 2*a0*a2*d1 - a0*d0*d2 &
            -a0*d1*d2 + a2**2*d0 + a2**2*d1 + a2*d0*d2 + a2*d1*d2 &
            + d0*d1*d2)/(a0 - a2 - d2)**2
    end subroutine get_coefs_Lijk_aac

    subroutine get_coefs_Lijk_abb(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: cb0, cb1, cb2, cb3
        !
        cb0 =  k0*k1/k2
        cb1 = cb0*d0*(a0 - a2)*(a0 - a1 - d1)/(a0 - a1)**2
        cb3 = -cb0*d1*(a1 - a2)*(a0 - a1 + d0)/(a0 - a1)
        cb2 = -cb0*(a0**2*a1 - a0**2*a2 - a0**2*d1 - 2*a0*a1**2 + 2*a0*a1*a2 &
            + a0*a1*d0 + 2*a0*a1*d1 - a0*a2*d0 - a0*d0*d1 + a1**3 - a1**2*a2 &
            - a1**2*d0 - a1**2*d1 + a1*a2*d0 + a2*d0*d1)/(a0 - a1)**2
    end subroutine get_coefs_Lijk_abb

    subroutine get_coefs_Lijk_aba(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: cb0, cb1, cb2, cb3
        !
        cb0 =  k0*k1/k2
        cb3 = -cb0*d0*(a0 - a2)*(a0 - a1 - d1)/(a0 - a1)
        cb2 = -cb0*d1*(a1 - a2)* (a0 - a1 + d0)/(a0 - a1)**2
        cb1 = -cb0*(a0**3 - 2*a0**2 *a1 - a0**2*a2 - a0**2*d0 - a0**2*d1 &
            + a0*a1**2 + 2*a0*a1*a2 + 2*a0*a1*d0 + a0*a1*d1 + a0*a2*d1 - &
            a1**2*a2 - a1**2*d0 - a1*a2*d1 - a1*d0*d1 + a2*d0*d1)/(a0 - a1)**2
    end subroutine get_coefs_Lijk_aba

    subroutine get_coefs_Lijk_aaa(k0,k1,k2,a0,a1,a2,d0,d1,d2, cb0, cb1, cb2, cb3)
        real(fpp), intent(in)  :: k0,k1,k2,a0,a1,a2,d0,d1,d2
        real(fpp), intent(out) :: cb0, cb1, cb2, cb3
        !
        cb0 =  k0*k1/k2
        cb3 = - cb0*d0*d1*(a0 - a2)
        cb1 = - cb0*(a0 - a2 - d0 - d1)
        cb2 = - cb0*(a0*d0 + a0*d1 - a2*d0 - a2*d1 - d0*d1)
    end subroutine get_coefs_Lijk_aaa

    subroutine cpml_only_one_root(dom, i, j, k, bnum, ee, ndir)
        class(dombase_cpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee, ndir
        integer :: idx

        select case(ndir)
        case(2)
            call cpml_only_one_root_2d(dom, i, j, k, bnum, ee)
        case(3)
            call cpml_only_one_root_3d(dom, i, j, k, bnum, ee)
        end select
        if (ndir>=1) then
            if (dom%Alpha_0(ee,i,j,k,bnum)<0.) dom%Alpha_0(ee,i,j,k,bnum) = 0.
        endif
        if (ndir>=2) then
            idx = dom%I1(ee,bnum)
            if (dom%Alpha_1(i,j,k,idx)<0.) dom%Alpha_1(i,j,k,idx) = 0.
        endif
        if (ndir>=3) then
            idx = dom%I2(ee,bnum)
            if (dom%Alpha_2(i,j,k,idx)<0.) dom%Alpha_2(i,j,k,idx) = 0.
        endif
    end subroutine cpml_only_one_root

    subroutine cpml_only_one_root_2d(dom, i, j, k, bnum, ee)
        class(dombase_cpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        !
        real(fpp) :: a0, a1, d0, d1, cfp
        integer :: i1, i2

        if (dom%cpml_one_root == 0) return

        i1 = dom%I1(ee,bnum)
        i2 = dom%I2(ee,bnum)

        cfp = 1. + dom%cpml_one_root/100.

        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        a1 = dom%Alpha_1(i,j,k,i1)
        d1 = dom%dxi_k_1(i,j,k,i1)

        ! On modifie a1/b1 (en fonction de a0/b0). On ne modifie pas a0/b0.
        if (isclose(a1,a0)) then ! s0s1/s2 + s2=1
            a1 = a1*cfp
        end if
        if (isclose(a1+d1,a0)) then ! s0s2/s1 + s2=1
            d1 = d1*cfp
        end if
        if (isclose(a1,a0+d0)) then ! s1s2/s0 + s2=1
            d0 = d0*cfp
        end if
        dom%dxi_k_0(ee,i,j,k,bnum) = d0
        dom%Alpha_1(i,j,k,i1) = a1
        dom%dxi_k_1(i,j,k,i1) = d1
    end subroutine cpml_only_one_root_2d

    subroutine cpml_only_one_root_3d(dom, i, j, k, bnum, ee)
        class(dombase_cpml), intent (INOUT) :: dom
        integer, intent(in) :: i, j, k, bnum, ee
        !
        real(fpp) :: a0, a1, a2, d0, d1, d2, cfp
        integer :: i1, i2, sel

        if (dom%cpml_one_root == 0) return

        i1 = dom%I1(ee,bnum)
        i2 = dom%I2(ee,bnum)

        cfp = real(dom%cpml_one_root,fpp)/100.

        a0 = dom%Alpha_0(ee,i,j,k,bnum)
        d0 = dom%dxi_k_0(ee,i,j,k,bnum)
        a1 = dom%Alpha_1(i,j,k,i1)
        d1 = dom%dxi_k_1(i,j,k,i1)
        ! La premiere passe assure que a1/b1 > a0/b0
        ! On modifie a2/b2 (en fonction de a0/b0 et a1/b1). On ne modifie pas ni a0/b0, ni a1/b1.
        a2 = dom%Alpha_2(i,j,k,i2)
        d2 = dom%dxi_k_2(i,j,k,i2)

        sel = compare_roots(a0,a1,a2)
        select case(sel)
        case (CMP_ABC)
        case (CMP_ABA)
            write(*,*) "Case L ABA",a0,a2
            a2 = a0 + cfp*abs(a1-a0)
        case (CMP_AAC)
            write(*,*) "Case L AAC",a0,a1
            a1 = a0 + cfp*abs(a2-a1)
        case (CMP_ABB)
            write(*,*) "Case L ABB",a0,a1
            a2 = a1 + cfp*abs(a1-a0)
        case (CMP_AAA)
            write(*,*) "Case L AAA",a0,a1
            if (a0/=0) then
                a1 = a0*(1+cfp)
                a2 = a1*(1+cfp)
            else
                a1 = 1e-6
                a2 = 2e-6
            end if
        case default
            stop 1
        end select

        sel = compare_roots(a0,a1,a2+d2)
        select case(sel)
        case (CMP_ABC)
        case (CMP_ABA)
            write(*,*) "Case L012 ABA",a0,a2+d2
            d2 = d2*(1 + cfp)
        case (CMP_AAC)
            write(*,*) "Can't happen AAC:A0A1B2 ", a0, a1
            stop 1
        case (CMP_ABB)
            write(*,*) "Case L012 ABB",a1,a2+d2
            d2 = d1*(1 + cfp)
        case (CMP_AAA)
            write(*,*) "Can't happen AAA:A0A1B2"
            stop 1
        case default
            stop 1
        end select
        sel = compare_roots(a0,a1+d1,a2)
        select case(sel)
        case (CMP_ABC)
        case (CMP_ABA)
            write(*,*) "Can't happen ABA:A0B1A2"
            stop 1
        case (CMP_AAC)
            write(*,*) "Case L021 AAC",a0,a1+d1
            d1 = d1*(1 + cfp)
        case (CMP_ABB)
            write(*,*) "Case L021 ABB",a1+d1,a2
            d1 = d1*(1 + cfp)
        case (CMP_AAA)
            write(*,*) "Can't happen AAA:A0B1A2"
            stop 1
        case default
            stop 1
        end select
        sel = compare_roots(a0+d0,a1,a2)
        select case(sel)
        case (CMP_ABC)
        case (CMP_ABA)
            write(*,*) "Case L120 ABA",a0+d0,a2
            d0 = d0*(1 + cfp)
        case (CMP_AAC)
            write(*,*) "Case L120 AAC",a0+d0,a1
            d0 = d0*(1 + cfp)
        case (CMP_ABB)
            write(*,*) "Can't happen ABB:B0A1A2"
            stop 1
        case (CMP_AAA)
            write(*,*) "Can't happen"
            stop 1
        case default
            stop 1
        end select

        dom%Alpha_0(ee,i,j,k,bnum) = a0
        dom%dxi_k_0(ee,i,j,k,bnum) = d0
        dom%Alpha_1(i,j,k,i1)      = a1
        dom%dxi_k_1(i,j,k,i1)      = d1
        dom%Alpha_2(i,j,k,i2)      = a2
        dom%dxi_k_2(i,j,k,i2)      = d2

        sel = compare_roots(a0,a1,a2)
        if (sel/=CMP_ABC) stop 1
        sel = compare_roots(a0,a1,a2+d2)
        if (sel/=CMP_ABC) stop 1
        sel = compare_roots(a0,a1+d1,a2)
        if (sel/=CMP_ABC) stop 1
        sel = compare_roots(a0+d0,a1,a2)
        if (sel/=CMP_ABC) stop 1

    end subroutine cpml_only_one_root_3d

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
