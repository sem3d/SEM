!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file define_arr.F90
!!\brief Contient la subroutine define_arrays.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief La routine define_arrays definit des tableaux concernant les proprietes physiques du materiau.
!! les proprietes physiques concernees sont : l'elasticite, le tenseur de contraintes, les conditions limites,...
!! Elle assure egalement les communications entre processeurs de SEM pour diffuser l'information.
!! \param type (domain),intent (INOUT), target Tdomain
!<


subroutine define_arrays(Tdomain)

    use sdomain
    use scompute_coeff_HDG
    use mpi

    implicit none
    type (domain),intent (INOUT), target :: Tdomain

    ! local variables

    integer :: n,mat,ngllx,ngllz,ngll,i,j,idef,n_elem,w_face,nv_aus,nf,npow,powOmc,nface
    integer :: i_send, n_face_pointed, i_proc, nv, i_stock, tag_send, tag_receive, ierr
    integer ,  dimension (MPI_STATUS_SIZE) :: status
    real :: vp,ri,rj,dx,LocMassmat_Vertex,Apow,OmegaCprime,PI,Omega_c,dxdxi,dzdeta, Dt
    real, external :: pow
    real, dimension (:), allocatable :: LocMassMat1D, LocMassMat1D_Down, Send_bt, Receive_Bt
    real, dimension (:,:), allocatable :: xix,etax, xiz,etaz,Jac, Rlam,Rmu,RKmod,Whei,Id,wx, wz
    real, dimension (:,:), allocatable :: LocMassMat,OmegaCutx,OmegaCutz,du_du_x,du_du_z
    real, dimension (:,:), allocatable :: duux,duuz,wx_prime,wz_prime

    ! Gaetano Festa, modified 01/06/2004
    ! Modification (MPI) 13/10/2005
    ! Modifications by S. Terrana for introduction of CPML : may 2014
    ! Modifications by S. Terrana for introduction of ADEPML : june 2014

    PI = 4. * atan(1.)

    ! Attribute elastic properties from material

    do n = 0, Tdomain%n_elem - 1

        mat = Tdomain%specel(n)%mat_index
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz

        do j = 0, ngllz - 1
            do i = 0, ngllx - 1
                Tdomain%specel(n)%Density(i,j) = Tdomain%sSubDomain(mat)%Ddensity
                Tdomain%specel(n)%Lambda(i,j) = Tdomain%sSubDomain(mat)%DLambda
                Tdomain%specel(n)%Mu(i,j) = Tdomain%sSubDomain(mat)%DMu
            enddo
        enddo

        ! For Atmospheric Waveguide :
        !call modify_atmospheric_rho(Tdomain,n)

        allocate (xix (0:ngllx-1,0:ngllz-1))
        allocate (xiz (0:ngllx-1,0:ngllz-1))
        allocate (etax (0:ngllx-1,0:ngllz-1))
        allocate (etaz (0:ngllx-1,0:ngllz-1))
        allocate (RKmod (0:ngllx-1,0:ngllz-1))
        allocate (Whei (0:ngllx-1,0:ngllz-1))
        allocate (wx (0:ngllx-1,0:ngllz-1))
        allocate (wz (0:ngllx-1,0:ngllz-1))
        allocate (wx_prime (0:ngllx-1,0:ngllz-1))
        allocate (wz_prime (0:ngllx-1,0:ngllz-1))
        allocate (OmegaCutx (0:ngllx-1,0:ngllz-1))
        allocate (OmegaCutz (0:ngllx-1,0:ngllz-1))
        allocate (du_du_x (0:ngllx-1,0:ngllz-1))
        allocate (du_du_z (0:ngllx-1,0:ngllz-1))
        allocate (duux (0:ngllx-1,0:ngllz-1))
        allocate (duuz (0:ngllx-1,0:ngllz-1))
        allocate (Id (0:ngllx-1,0:ngllz-1))
        allocate (Jac(0:ngllx-1,0:ngllz-1))
        allocate (Rlam(0:ngllx-1,0:ngllz-1))
        allocate (Rmu(0:ngllx-1,0:ngllz-1))

        do j = 0,ngllz-1
            do i = 0,ngllx-1
                Whei (i,j) = Tdomain%sSubdomain(mat)%GLLwx(i)*  Tdomain%sSubdomain(mat)%GLLwz(j)
            enddo
        enddo

        xix  = Tdomain%specel(n)%InvGrad(:,:,0,0)
        xiz = Tdomain%specel(n)%InvGrad(:,:,1,0)
        etax  = Tdomain%specel(n)%InvGrad(:,:,0,1)
        etaz = Tdomain%specel(n)%InvGrad(:,:,1,1)
        !
        Jac  = Tdomain%specel(n)%Jacob
        Rlam = Tdomain%specel(n)%Lambda
        Rmu  = Tdomain%specel(n)%Mu
        RKmod = Rlam + 2. * Rmu
        Tdomain%specel(n)%MassMat  = Whei*Tdomain%specel(n)%Density*Jac

        if(Tdomain%specel(n)%Type_DG .NE. GALERKIN_CONT) then ! Discontinuous Galerkin
            if(Tdomain%specel(n)%acoustic) then
                Tdomain%specel(n)%Acoeff(:,:,0) = Whei*xix*Jac
                Tdomain%specel(n)%Acoeff(:,:,1) = Whei*etax*Jac
                Tdomain%specel(n)%Acoeff(:,:,2) = Whei*xiz*Jac
                Tdomain%specel(n)%Acoeff(:,:,3) = Whei*etaz*Jac
                Tdomain%specel(n)%Acoeff(:,:,4)= Whei*Jac
            else ! Discontinuous Galerkin, usual
                Tdomain%specel(n)%Acoeff(:,:,0) = Whei*xix*Jac
                Tdomain%specel(n)%Acoeff(:,:,1) = Whei*etax*Jac
                Tdomain%specel(n)%Acoeff(:,:,2) = Whei*xiz*Jac
                Tdomain%specel(n)%Acoeff(:,:,3) = Whei*etaz*Jac
                Tdomain%specel(n)%Acoeff(:,:,4) = Whei*xix*Rlam*Jac
                Tdomain%specel(n)%Acoeff(:,:,5) = 2.*Whei*xix*Rmu*Jac
                Tdomain%specel(n)%Acoeff(:,:,6) = Whei*xiz*Rlam*Jac
                Tdomain%specel(n)%Acoeff(:,:,7) = 2.*Whei*xiz*Rmu*Jac
                Tdomain%specel(n)%Acoeff(:,:,8) = Whei*etax*Rlam*Jac
                Tdomain%specel(n)%Acoeff(:,:,9) = 2.*Whei*etax*Rmu*Jac
                Tdomain%specel(n)%Acoeff(:,:,10)= Whei*etaz*Rlam*Jac
                Tdomain%specel(n)%Acoeff(:,:,11)= 2.*Whei*etaz*Rmu*Jac
                Tdomain%specel(n)%Acoeff(:,:,12)= Whei*Jac
            endif
        else ! Continuous Galerkin (usual SEM) without PML
            if ((.NOT. Tdomain%specel(n)%PML)) then ! .AND. (.NOT. Tdomain%type_timeInteg == TIME_INTEG_RK4)) then
                Tdomain%specel(n)%Acoeff(:,:,0) = -Whei*(RKmod*xix**2+Rmu*xiz**2) *Jac
                Tdomain%specel(n)%Acoeff(:,:,1) = -Whei*(RKmod*xix*etax+Rmu*  xiz*etaz)*Jac
                Tdomain%specel(n)%Acoeff(:,:,2) = -Whei*(Rlam+Rmu)*xix*xiz*Jac
                Tdomain%specel(n)%Acoeff(:,:,3) = -Whei*(Rlam*xix*etaz+Rmu*xiz*etax)*Jac
                Tdomain%specel(n)%Acoeff(:,:,4) = -Whei*(RKmod*etax**2+Rmu*etaz**2)*Jac
                Tdomain%specel(n)%Acoeff(:,:,5) = -Whei*(Rlam*etax*xiz+Rmu*etaz*xix)*Jac
                Tdomain%specel(n)%Acoeff(:,:,6) = -Whei*(Rlam+Rmu)*etaz*etax*Jac
                Tdomain%specel(n)%Acoeff(:,:,7) = -Whei*(RKmod*xiz**2+Rmu*xix**2)*Jac
                Tdomain%specel(n)%Acoeff(:,:,8) = -Whei*(RKmod*xiz*etaz+Rmu*xix*etax)*Jac
                Tdomain%specel(n)%Acoeff(:,:,9) = -Whei*(RKmod*etaz**2+Rmu* etax**2)*Jac
            else
                Tdomain%specel(n)%Acoeff(:,:,0) = Rkmod*xix
                Tdomain%specel(n)%Acoeff(:,:,1) = Rkmod*etax
                Tdomain%specel(n)%Acoeff(:,:,2) = Rlam * xiz
                Tdomain%specel(n)%Acoeff(:,:,3) = Rlam *etaz
                Tdomain%specel(n)%Acoeff(:,:,4) = Rlam*xix
                Tdomain%specel(n)%Acoeff(:,:,5) = Rlam*etax
                Tdomain%specel(n)%Acoeff(:,:,6) = Rkmod * xiz
                Tdomain%specel(n)%Acoeff(:,:,7) = Rkmod *etaz
                Tdomain%specel(n)%Acoeff(:,:,8) = Rmu*xiz
                Tdomain%specel(n)%Acoeff(:,:,9) = Rmu*etaz
                Tdomain%specel(n)%Acoeff(:,:,10) = Rmu * xix
                Tdomain%specel(n)%Acoeff(:,:,11) = Rmu *etax
                Tdomain%specel(n)%Acoeff(:,:,12) = -Whei*xix*Jac
                Tdomain%specel(n)%Acoeff(:,:,13) = -Whei*etax*Jac
                Tdomain%specel(n)%Acoeff(:,:,14) = -Whei * xiz*Jac
                Tdomain%specel(n)%Acoeff(:,:,15) = -Whei *etaz*Jac
                if (Tdomain%specel(n)%CPML) then
                    Tdomain%specel(n)%Acoeff(:,:,16) = Whei * Jac
                    Tdomain%specel(n)%Acoeff(:,:,17) = 1./(Whei*Jac)
                endif
            endif
        endif

        if (Tdomain%specel(n)%PML) then
            ! Defining PML properties
            Apow = Tdomain%sSubdomain(mat)%Apow ; npow = Tdomain%sSubdomain(mat)%npow
            Omega_c = 2. * PI * Tdomain%sSubdomain(mat)%freq ! Frequence de coupure
            ! PowOmc is the exponent of the power law of decreasing Omega_c (pulsation de coupure)
            ! in the PML. Usually, Omega_C obbey to a law Omega_c(x) = 2*pi*freq_c (1-(x/L)^{powOmc})
            ! powOmc is set to 1 because it produces better absorbtion on the cases we have studied.
            powOmc = 1
            if (Tdomain%sSubDomain(mat)%Px) then
                ! Computation of dx : the horizontal length of the PML element
                idef = Tdomain%specel(n)%Iglobnum (0,0); dx = Tdomain%GlobCoord (0,idef)
                idef = Tdomain%specel(n)%Iglobnum (ngllx-1,0); dx = abs (Tdomain%GlobCoord (0,idef) -dx);
                if (Tdomain%sSubDomain(mat)%Left) then
                    do i = 0,ngllx-1
                        ri = 0.5*(1+Tdomain%sSubDomain(mat)%GLLcx(ngllx-1-i))
                        vp = Rkmod(i,0)/Tdomain%specel(n)%Density(i,0)
                        vp = sqrt(vp)
                        dxdxi = -0.5
                        wx(i,0:ngllz-1) = pow (ri,vp,1,dx,Apow,npow)
                        wx_prime(i,0:ngllz-1) = dxdxi * npow * pow (ri,vp,1,dx,Apow,npow-1)
                        OmegaCutx(i,0:ngllz-1) = Omega_c * (1 - ri**PowOmc)
                        OmegaCprime = -Omega_c * PowOmc * ri**(PowOmc-1)
                        if (powOmc == 0) then
                            OmegaCutx(i,0:ngllz-1) = Omega_c
                            OmegaCprime = 0.
                        endif
                        du_du_x(i,0:ngllz-1) = dxdxi* (npow*pow(ri,vp,1,dx,Apow,npow-1) * (OmegaCutx(i,0) + wx(i,0)) &
                                                     - wx(i,0) * (npow*pow(ri,vp,1,dx,Apow,npow-1) + OmegaCprime))
                        duux(i,0:ngllz-1) = dxdxi * wx(i,0) * (npow*pow(ri,vp,1,dx,Apow,npow-1) + OmegaCprime) &
                                          * (wx(i,0) + OmegaCutx(i,0)) * Tdomain%sSubdomain(mat)%Dt
                    enddo
                else
                    do i = 0,ngllx-1
                        ri = 0.5*(1+Tdomain%sSubDomain(mat)%GLLcx(i))
                        vp = Rkmod(i,0)/Tdomain%specel(n)%Density(i,0)
                        vp = sqrt(vp)
                        dxdxi = 0.5
                        wx(i,0:ngllz-1) = pow (ri,vp,1,dx,Apow,npow)
                        wx_prime(i,0:ngllz-1) = dxdxi * npow * pow (ri,vp,1,dx,Apow,npow-1)
                        OmegaCutx(i,0:ngllz-1)  = Omega_c * (1 - ri**PowOmc)
                        OmegaCprime = -Omega_c * PowOmc * ri**(PowOmc-1)
                        if (powOmc == 0) then
                            OmegaCutx(i,0:ngllz-1) = Omega_c
                            OmegaCprime = 0.
                        endif
                        du_du_x(i,0:ngllz-1) = dxdxi* (npow*pow(ri,vp,1,dx,Apow,npow-1) * (OmegaCutx(i,0) + wx(i,0)) &
                                                     - wx(i,0) * (npow*pow(ri,vp,1,dx,Apow,npow-1) + OmegaCprime))
                        duux(i,0:ngllz-1) = dxdxi * wx(i,0) * (npow*pow(ri,vp,1,dx,Apow,npow-1) + OmegaCprime) &
                                          * (wx(i,0) + OmegaCutx(i,0)) * Tdomain%sSubdomain(mat)%Dt
                        !du_du_x(i,0:ngllz-1) = 0.5 * npow*pow(ri,vp,1,dx,Apow,npow-1) * (OmegaCutx(i,0) + wx(i,0)) &
                        !                     - wx(i,0) * (0.5*npow*pow(ri,vp,1,dx,Apow,npow-1) + OmegaCprime)
                        !duux(i,0:ngllz-1) = 0.5 * wx(i,0) * (npow*pow(ri,vp,1,dx,Apow,npow-1) + OmegaCprime) &
                        !                  * (wx(i,0) + OmegaCutx(i,0)) * Tdomain%sSubdomain(mat)%Dt
                    enddo
                endif
            else
                wx = 0.
            endif
            if (Tdomain%sSubDomain(mat)%Pz) then
                ! Computation of dx : the vertical heigth of the PML element
                idef = Tdomain%specel(n)%Iglobnum (0,0); dx = Tdomain%GlobCoord (1,idef)
                idef = Tdomain%specel(n)%Iglobnum (0,ngllz-1); dx = abs (Tdomain%GlobCoord (1,idef) -dx)
                if (Tdomain%sSubDomain(mat)%Down) then
                    do j = 0,ngllz-1
                        rj = 0.5*(1+Tdomain%sSubdomain(mat)%GLLcz(ngllz-1-j))
                        vp = Rkmod(0,j)/Tdomain%specel(n)%Density(0,j)
                        vp = sqrt(vp)
                        dzdeta = -0.5
                        wz(0:ngllx-1,j) = pow (rj,vp,1,dx,Apow,npow)
                        wz_prime(0:ngllx-1,j) = dzdeta * npow * pow (rj,vp,1,dx,Apow,npow-1)
                        OmegaCutz(0:ngllx-1,j) = Omega_c * (1 - rj**PowOmc)
                        OmegaCprime = -Omega_c * PowOmc * rj**(PowOmc-1)
                        if (powOmc == 0) then
                            OmegaCutz(0:ngllx-1,j) = Omega_c
                            OmegaCprime = 0.
                        endif
                        du_du_z(0:ngllx-1,j) = dzdeta* (npow*pow(rj,vp,1,dx,Apow,npow-1) * (OmegaCutz(0,j) + wz(0,j)) &
                                             - wz(0,j) * (npow*pow(rj,vp,1,dx,Apow,npow-1) + OmegaCprime))
                        duuz(0:ngllx-1,j) = dzdeta * wz(0,j) * (npow*pow(rj,vp,1,dx,Apow,npow-1) + OmegaCprime) &
                                          * (wz(0,j) + OmegaCutz(0,j)) * Tdomain%sSubdomain(mat)%Dt
                    enddo
                else
                    do j = 0,ngllz-1
                        rj = 0.5*(1+Tdomain%sSubdomain(mat)%GLLcz(j))
                        vp = Rkmod(0,j)/Tdomain%specel(n)%Density(0,j)
                        vp = sqrt(vp)
                        dzdeta = 0.5
                        wz(0:ngllx-1,j) = pow (rj,vp,1,dx,Apow,npow)
                        wz_prime(0:ngllx-1,j) = dzdeta * npow * pow (rj,vp,1,dx,Apow,npow-1)
                        OmegaCutz(0:ngllx-1,j) = Omega_c * (1 - rj**PowOmc)
                        OmegaCprime = -Omega_c * PowOmc * rj**(PowOmc-1)
                        if (powOmc == 0) then
                            OmegaCutz(0:ngllx-1,j) = Omega_c
                            OmegaCprime = 0.
                        endif
                        du_du_z(0:ngllx-1,j) = dzdeta* (npow*pow(rj,vp,1,dx,Apow,npow-1) * (OmegaCutz(0,j) + wz(0,j)) &
                                             - wz(0,j) * (npow*pow(rj,vp,1,dx,Apow,npow-1) + OmegaCprime))
                        duuz(0:ngllx-1,j) = dzdeta * wz(0,j) * (npow*pow(rj,vp,1,dx,Apow,npow-1) + OmegaCprime) &
                                          * (wz(0,j) + OmegaCutz(0,j)) * Tdomain%sSubdomain(mat)%Dt
                        !du_du_z(0:ngllx-1,j) = 0.5*npow*pow(rj,vp,1,dx,Apow,npow-1) * (OmegaCutz(0,j) + wz(0,j)) &
                        !                     - wz(0,j) * (0.5*npow*pow(rj,vp,1,dx,Apow,npow-1) + OmegaCprime)
                        !duuz(0:ngllx-1,j) = wz(0,j) * (0.5*npow*pow(rj,vp,1,dx,Apow,npow-1) + OmegaCprime) &
                        !                  * (wz(0,j) + OmegaCutz(0,j)) * Tdomain%sSubdomain(mat)%Dt
                    enddo
                endif
            else
                wz = 0.
            endif
            Id = 1

            if (Tdomain%specel(n)%CPML .OR. Tdomain%specel(n)%ADEPML) then
                if (Tdomain%sSubDomain(mat)%Px) then
                    if (Tdomain%specel(n)%CPML) then
                        Tdomain%specel(n)%Bx(:,:) = exp(-(wx(:,:) + OmegaCutx(:,:)) * Tdomain%sSubdomain(mat)%Dt)
                        Tdomain%specel(n)%Ax(:,:) = wx(:,:) *  (Tdomain%specel(n)%Bx(:,:) - Id(:,:)) / (wx(:,:) + OmegaCutx(:,:))
                        if (Tdomain%sSubDomain(mat)%freq == 0.) Tdomain%specel(n)%Ax(:,:) = Tdomain%specel(n)%Bx (:,:) - Id(:,:)
                    else ! ADE-PML for HDG
                        Tdomain%specel(n)%Bx(:,:) = wx(:,:) + OmegaCutx(:,:)
                        Tdomain%specel(n)%Ax(:,:) = wx(:,:)
                    endif
               endif
               if (Tdomain%sSubDomain(mat)%Pz) then
                   if (Tdomain%specel(n)%CPML) then
                       Tdomain%specel(n)%Bz(:,:) = exp(-(wz(:,:) + OmegaCutz(:,:)) * Tdomain%sSubdomain(mat)%Dt)
                       Tdomain%specel(n)%Az(:,:) = wz(:,:) *  (Tdomain%specel(n)%Bz(:,:) - Id(:,:)) / (wz(:,:) + OmegaCutz(:,:))
                       if (Tdomain%sSubDomain(mat)%freq == 0.) Tdomain%specel(n)%Az(:,:) = Tdomain%specel(n)%Bz (:,:) - Id(:,:)
                   else ! ADE-PML for HDG
                       Tdomain%specel(n)%Bz(:,:) = wz(:,:) + OmegaCutz(:,:)
                       Tdomain%specel(n)%Az(:,:) = wz(:,:)
                   endif
               endif
            else ! Usual PML
                Tdomain%specel(n)%DumpSx(:,:,1) = Id(:,:) + 0.5 * Tdomain%sSubdomain(mat)%Dt * wx(:,:)
                Tdomain%specel(n)%DumpSx (:,:,1) = 1./ Tdomain%specel(n)%DumpSx (:,:,1)
                Tdomain%specel(n)%DumpSx (:,:,0) = (Id(:,:) - Tdomain%sSubdomain(mat)%Dt * 0.5 * wx(:,:)) * Tdomain%specel(n)%DumpSx(:,:,1)

                Tdomain%specel(n)%DumpSz(:,:,1) = Id(:,:) + 0.5 * Tdomain%sSubdomain(mat)%Dt * wz(:,:)
                Tdomain%specel(n)%DumpSz (:,:,1) = 1./ Tdomain%specel(n)%DumpSz (:,:,1)
                Tdomain%specel(n)%DumpSz (:,:,0) = (Id(:,:) - Tdomain%sSubdomain(mat)%Dt * 0.5 * wz(:,:)) * Tdomain%specel(n)%DumpSz(:,:,1)

                Tdomain%specel(n)%DumpMass(:,:,0) = 0.5 * Tdomain%specel(n)%Density * Whei * Tdomain%sSubdomain(mat)%Dt * wx * Jac
                Tdomain%specel(n)%DumpMass(:,:,1) = 0.5 * Tdomain%specel(n)%Density * Whei * Tdomain%sSubdomain(mat)%Dt * wz * Jac
            endif

        endif
        deallocate(xix,xiz,etax,etaz,Id,wx,wz,OmegaCutx,OmegaCutz,du_du_x,du_du_z, &
                   duux,duuz,wx_prime,wz_prime,Whei,RKmod,Jac,Rmu,Rlam)

    enddo


    ! Communication inside the processor - Assembling Masses

    do nf = 0, Tdomain%n_face-1
        n_elem = Tdomain%sFace(nf)%Near_element(0)
        w_face = Tdomain%sFace(nf)%Which_face(0)
        ! Assemblage des masses sur les Faces
        if (Tdomain%specel(n_elem)%type_DG == GALERKIN_CONT) &
            call getMass_element2face (Tdomain,n_elem,nf,w_face,.true.)
        if (Tdomain%sFace(nf)%PML .and. (.not.Tdomain%sFace(nf)%CPML) .and. &
           (.not.Tdomain%sFace(nf)%ADEPML)) call getDumpMass_element2face (Tdomain,n_elem,nf,w_face,.true.)
        n_elem = Tdomain%sFace(nf)%Near_element(1)
        if (n_elem > -1) then
            w_face = Tdomain%sFace(nf)%Which_face(1)
            if (Tdomain%specel(n_elem)%type_DG == GALERKIN_CONT) &
                call getMass_element2face (Tdomain,n_elem,nf,w_face,Tdomain%sFace(nf)%coherency)
            if (Tdomain%sFace(nf)%PML .and. (.not.Tdomain%sFace(nf)%CPML) .and. &
               (.not.Tdomain%sFace(nf)%ADEPML)) call getDumpMass_element2face (Tdomain,n_elem,nf,w_face,Tdomain%sFace(nf)%coherency)
        endif
    enddo

    ! Assemblage des masses sur les vertexes
    do n = 0, Tdomain%n_elem -1
      if (Tdomain%specel(n)%type_DG == GALERKIN_CONT) then
        ngllx = Tdomain%specel(n)%ngllx; ngllz = Tdomain%specel(n)%ngllz
        nv_aus = Tdomain%specel(n)%Near_Vertex(0)
        Tdomain%sVertex(nv_aus)%MassMat = Tdomain%sVertex(nv_aus)%MassMat + Tdomain%specel(n)%MassMat(0,0)
        if (Tdomain%sVertex(nv_aus)%PML .and. (.not.Tdomain%sVertex(nv_aus)%CPML) &
                                        .and. (.not.Tdomain%sVertex(nv_aus)%ADEPML)) then
            Tdomain%sVertex(nv_aus)%DumpMass(:) = Tdomain%sVertex(nv_aus)%DumpMass(:) &
                                                + Tdomain%specel(n)%DumpMass(0,0,:)
        endif

        nv_aus = Tdomain%specel(n)%Near_Vertex(1)
        Tdomain%sVertex(nv_aus)%MassMat = Tdomain%sVertex(nv_aus)%MassMat + Tdomain%specel(n)%MassMat(ngllx-1,0)
        if (Tdomain%sVertex(nv_aus)%PML .and. (.not.Tdomain%sVertex(nv_aus)%CPML) &
                                        .and. (.not.Tdomain%sVertex(nv_aus)%ADEPML)) then
            Tdomain%sVertex(nv_aus)%DumpMass(:) = Tdomain%sVertex(nv_aus)%DumpMass(:) &
                                                + Tdomain%specel(n)%DumpMass(ngllx-1,0,:)
        endif

        nv_aus = Tdomain%specel(n)%Near_Vertex(2)
        Tdomain%sVertex(nv_aus)%MassMat = Tdomain%sVertex(nv_aus)%MassMat + Tdomain%specel(n)%MassMat(ngllx-1,ngllz-1)
        if (Tdomain%sVertex(nv_aus)%PML .and. (.not.Tdomain%sVertex(nv_aus)%CPML) &
                                        .and. (.not.Tdomain%sVertex(nv_aus)%ADEPML)) then
            Tdomain%sVertex(nv_aus)%DumpMass(:) = Tdomain%sVertex(nv_aus)%DumpMass(:) &
                                                + Tdomain%specel(n)%DumpMass(ngllx-1,ngllz-1,:)
        endif

        nv_aus = Tdomain%specel(n)%Near_Vertex(3)
        Tdomain%sVertex(nv_aus)%MassMat = Tdomain%sVertex(nv_aus)%MassMat + Tdomain%specel(n)%MassMat(0,ngllz-1)
        if (Tdomain%sVertex(nv_aus)%PML .and. (.not.Tdomain%sVertex(nv_aus)%CPML) &
                                        .and. (.not.Tdomain%sVertex(nv_aus)%ADEPML)) then
            Tdomain%sVertex(nv_aus)%DumpMass(:) = Tdomain%sVertex(nv_aus)%DumpMass(:) &
                                                + Tdomain%specel(n)%DumpMass(0,ngllz-1,:)
        endif
      endif
    enddo

    !################  Communications between processors  ###################

    ! Duplicate Vertex values
    do n = 0, Tdomain%n_communications - 1
        do nv = 0, Tdomain%sWall(n)%n_vertices-1
            nv_aus = Tdomain%sWall(n)%Vertex_List(nv)
            Tdomain%sVertex(nv_aus)%Double_Value(0) = Tdomain%sVertex(nv_aus)%MassMat
        enddo
    enddo

    ! Communications for PML
    if (Tdomain%any_PML) then
        do i_proc = 0, Tdomain%n_communications - 1
            allocate (Tdomain%sWall(i_proc)%Send_data_2(0:Tdomain%sWall(i_proc)%n_points_pml-1,0:1))
            allocate (Tdomain%sWall(i_proc)%Receive_data_2(0:Tdomain%sWall(i_proc)%n_points_pml-1,0:1))
            Tdomain%sWall(i_proc)%Send_data_2 = 0;Tdomain%sWall(i_proc)%Receive_data_2 = 0
            i_send = Tdomain%Communication_list (i_proc)
            i_stock = 0
            do nf = 0, Tdomain%sWall(i_proc)%n_pml_faces - 1
                n_face_pointed = Tdomain%sWall(i_proc)%FacePML_List(nf)
                ngll = Tdomain%sFace(n_face_pointed)%ngll
                if (Tdomain%sWall(i_proc)%FacePML_Coherency(nf)) then
                    Tdomain%sWall(i_proc)%Send_data_2(i_stock:i_stock+ngll-3,0:1) = Tdomain%sFace(n_face_pointed)%DumpMass (1:ngll-2,0:1)
                else
                    do j = 1, ngll-2
                        Tdomain%sWall(i_proc)%Send_data_2(i_stock+j-1,0:1) = Tdomain%sFace(n_face_pointed)%DumpMass (ngll-1-j,0:1)
                    enddo
                endif
                i_stock = i_stock + ngll - 2
            enddo

            tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank + 600
            tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 600
            if (Tdomain%sWall(i_proc)%n_points_pml > 0) then
                call MPI_SEND (Tdomain%sWall(i_proc)%Send_data_2,2*Tdomain%sWall(i_proc)%n_points_pml, &
                    MPI_DOUBLE_PRECISION, i_send, tag_send, Tdomain%communicateur, ierr )
                call MPI_RECV (Tdomain%sWall(i_proc)%Receive_data_2, 2*Tdomain%sWall(i_proc)%n_points_pml, &
                    MPI_DOUBLE_PRECISION, i_send, tag_receive, Tdomain%communicateur, status, ierr )
            endif

            i_stock = 0
            do nf = 0, Tdomain%sWall(i_proc)%n_pml_faces - 1
                n_face_pointed = Tdomain%sWall(i_proc)%FacePML_List(nf)
                ngll = Tdomain%sFace(n_face_pointed)%ngll
                if (Tdomain%sWall(i_proc)%FacePML_Coherency(nf)) then
                    Tdomain%sFace(n_face_pointed)%DumpMass (1:ngll-2,0:1) =  Tdomain%sFace(n_face_pointed)%DumpMass(1:ngll-2,0:1) + &
                        Tdomain%sWall(i_proc)%Receive_data_2(i_stock:i_stock+ngll-3,0:1)
                else
                    do j = 1, ngll-2
                        Tdomain%sFace(n_face_pointed)%DumpMass (ngll-1-j,0:1) =  Tdomain%sFace(n_face_pointed)%DumpMass (ngll-1-j,0:1)+ &
                            Tdomain%sWall(i_proc)%Receive_data_2(i_stock+j-1,0:1)
                    enddo
                endif
                i_stock = i_stock + ngll - 2
            enddo
            deallocate (Tdomain%sWall(i_proc)%Send_data_2)
            deallocate (Tdomain%sWall(i_proc)%Receive_data_2)
        enddo
    endif

    ! Communications for mass assembling
    do i_proc = 0, Tdomain%n_communications - 1
        allocate (Tdomain%sWall(i_proc)%Send_data_1(0:Tdomain%sWall(i_proc)%n_points-1))
        allocate (Tdomain%sWall(i_proc)%Receive_data_1(0:Tdomain%sWall(i_proc)%n_points-1))
        i_send = Tdomain%Communication_list (i_proc)
        i_stock = 0
        do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
            n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
            ngll = Tdomain%sFace(n_face_pointed)%ngll
            if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                Tdomain%sWall(i_proc)%Send_data_1(i_stock:i_stock+ngll-3) = Tdomain%sFace(n_face_pointed)%MassMat (1:ngll-2)
            else
                do j = 1, ngll-2
                    Tdomain%sWall(i_proc)%Send_data_1(i_stock+j-1) = Tdomain%sFace(n_face_pointed)%MassMat (ngll-1-j)
                enddo
            endif
            i_stock = i_stock + ngll - 2
        enddo
        do nv = 0, Tdomain%sWall(i_proc)%n_vertices - 1
            nv_aus =  Tdomain%sWall(i_proc)%Vertex_List(nv)
            Tdomain%sWall(i_proc)%Send_data_1 (i_stock)  = Tdomain%sVertex(nv_aus)%Double_Value(0)
            i_stock = i_stock + 1
        enddo
        tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank + 700
        tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 700

        call MPI_SEND (Tdomain%sWall(i_proc)%Send_data_1, Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION, i_send, &
            tag_send, Tdomain%communicateur, ierr )
        call MPI_RECV (Tdomain%sWall(i_proc)%Receive_data_1, Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION, i_send, &
            tag_receive, Tdomain%communicateur, status, ierr )

        i_stock = 0
        do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
            n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
            ngll = Tdomain%sFace(n_face_pointed)%ngll
            if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                Tdomain%sFace(n_face_pointed)%MassMat (1:ngll-2) =  Tdomain%sFace(n_face_pointed)%MassMat (1:ngll-2) + &
                    Tdomain%sWall(i_proc)%Receive_data_1(i_stock:i_stock+ngll-3)
            else
                do j = 1, ngll-2
                    Tdomain%sFace(n_face_pointed)%MassMat (ngll-1-j) =  Tdomain%sFace(n_face_pointed)%MassMat (ngll-1-j)+ &
                        Tdomain%sWall(i_proc)%Receive_data_1(i_stock+j-1)
                enddo
            endif
            i_stock = i_stock + ngll - 2
        enddo
        do nv = 0, Tdomain%sWall(i_proc)%n_vertices - 1
            nv_aus =  Tdomain%sWall(i_proc)%Vertex_List(nv)
            Tdomain%sVertex(nv_aus)%MassMat = Tdomain%sVertex(nv_aus)%MassMat  + Tdomain%sWall(i_proc)%Receive_data_1 (i_stock)
            i_stock = i_stock + 1
        enddo
        deallocate (Tdomain%sWall(i_proc)%Send_data_1)
        deallocate (Tdomain%sWall(i_proc)%Receive_data_1)
    enddo

    ! Exchanging datas Lambda and Mu for DG faces
    do i_proc = 0, Tdomain%n_communications - 1
        allocate (Tdomain%sWall(i_proc)%Send_data_2(0:Tdomain%sWall(i_proc)%n_points-1,0:1))
        allocate (Tdomain%sWall(i_proc)%Receive_data_2(0:Tdomain%sWall(i_proc)%n_points-1,0:1))
        i_send = Tdomain%Communication_list (i_proc)
        i_stock = 0
        do nf = 0, Tdomain%sWall(i_proc)%n_faces-1
            n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
            ngll = Tdomain%sFace(n_face_pointed)%ngll
            if(Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                Tdomain%sWall(i_proc)%Send_data_2(i_stock:i_stock+ngll-1,0) = Tdomain%sFace(n_face_pointed)%Lambda_m
                Tdomain%sWall(i_proc)%Send_data_2(i_stock:i_stock+ngll-1,1) = Tdomain%sFace(n_face_pointed)%Mu_m
            else
                do j=0,ngll-1
                    Tdomain%sWall(i_proc)%Send_data_2(i_stock+j,0)= Tdomain%sFace(n_face_pointed)%Lambda_m(ngll-1-j)
                    Tdomain%sWall(i_proc)%Send_data_2(i_stock+j,1)= Tdomain%sFace(n_face_pointed)%Mu_m(ngll-1-j)
                enddo
            endif
            i_stock = i_stock + ngll
        enddo
        tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank + 800
        tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 800

        call MPI_SEND (Tdomain%sWall(i_proc)%Send_data_2, 2*Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION, &
            i_send, tag_send, Tdomain%communicateur, ierr )
        call MPI_RECV (Tdomain%sWall(i_proc)%Receive_data_2,2*Tdomain%sWall(i_proc)%n_points, MPI_DOUBLE_PRECISION,&
            i_send, tag_receive, Tdomain%communicateur, status, ierr )

        i_stock = 0
        do nf = 0, Tdomain%sWall(i_proc)%n_faces - 1
            n_face_pointed = Tdomain%sWall(i_proc)%Face_List(nf)
            ngll = Tdomain%sFace(n_face_pointed)%ngll
            if (Tdomain%sWall(i_proc)%Face_Coherency(nf)) then
                Tdomain%sFace(n_face_pointed)%Lambda_p(0:ngll-1) = Tdomain%sWall(i_proc)%Receive_data_2(i_stock:i_stock+ngll-1,0)
                Tdomain%sFace(n_face_pointed)%Mu_p(0:ngll-1)     = Tdomain%sWall(i_proc)%Receive_data_2(i_stock:i_stock+ngll-1,1)
            else
                do j=0,ngll-1
                    Tdomain%sFace(n_face_pointed)%Lambda_p(0:ngll-1-j) = Tdomain%sWall(i_proc)%Receive_data_2(i_stock+j,0)
                    Tdomain%sFace(n_face_pointed)%Mu_p(0:ngll-1-j) = Tdomain%sWall(i_proc)%Receive_data_2(i_stock+j,1)
                enddo
            endif
            i_stock = i_stock + ngll
        enddo
        deallocate (Tdomain%sWall(i_proc)%Send_data_2)
        deallocate (Tdomain%sWall(i_proc)%Receive_data_2)
    enddo


    ! Transmission des coeffs Mu et Lambda aux faces
    do n = 0,Tdomain%n_elem-1
        do nf=0,3
            nface = Tdomain%specel(n)%Near_Face(nf)
            if((Tdomain%Specel(n)%Type_DG == GALERKIN_DG_WEAK) .OR. &
               (Tdomain%Specel(n)%Type_DG == GALERKIN_DG_STRONG)) &
                 call get_MuLambdaRho_el2f(Tdomain,n,nface)
        enddo
    enddo
    ! Calcul des coefficients pour les Fluxs Godunov (DG ou HDG)
    do nf = 0, Tdomain%n_face-1
       if(Tdomain%sFace(nf)%type_Flux .EQ. FLUX_GODUNOV) &
           call coeffs_flux_godunov(Tdomain,nf)
    enddo
    ! Prolongement par continuite des proprietes du milieu pour surface libre ou absorbante
    do nf = 0, Tdomain%n_face-1
        if(Tdomain%sFace(nf)%freesurf .OR. Tdomain%sFace(nf)%abs) then
            if(Tdomain%sFace(nf)%Type_Flux .EQ. FLUX_GODUNOV) &
                call coeffs_freesurf_abs(Tdomain,nf)
        endif
    enddo
    ! Calcul des matrices de Penalisation pour les elements HDG
    do n = 0, Tdomain%n_elem-1
        if(Tdomain%Specel(n)%Type_DG .EQ. GALERKIN_HDG_RP) then
            call compute_MatPen(Tdomain, n)
            call build_K_Expl(Tdomain, n)
        endif
    enddo

    ! Calcul des matrices de coefficients CA^-1 et ED^-1 pour HDG en semi-implicite
    if (Tdomain%Implicitness == TIME_INTEG_SEMI_IMPLICIT) then
        ! Calcul des matrices elementaires pour chaque elem
        do n = 0, Tdomain%n_elem-1
            mat = Tdomain%Specel(n)%mat_index
            Dt  = Tdomain%sSubDomain(mat)%Dt
            call compute_CAinv(Tdomain%Specel(n))
            call compute_EDinv(Tdomain%Specel(n),0.5*Dt)
            !call compute_EJ(Tdomain%Specel(n),0.5*Dt)
            call build_K_on_face(Tdomain,n,0.5*Dt)
            call build_K_on_vertex(Tdomain,n,0.5*Dt)
        enddo
        ! Inversion des matrices K sur les vertexs
        do n = 0, Tdomain%n_vertex-1
            call invert_K_vertex(Tdomain%sVertex(n),n)
        enddo
        do nf = 0, Tdomain%n_face-1
            call invert_K_face(Tdomain%sFace(nf))
        enddo
    endif

    ! Inversion des matrices K sur les faces pour explicit
    do nf = 0, Tdomain%n_face-1
        if (Tdomain%sFace(nf)%type_Flux .EQ. FLUX_HDG) &
            call invert_K_face_Expl(Tdomain%sFace(nf))
    enddo

    ! Coefficients pour visualisation des vertexs (HDG)
    call assemble_coeffs_proj_Vhat (Tdomain)

    ! Preparation des Mortars
    call prepare_mortars(Tdomain)

    ! Preparing and allocating vectors to be exchanged at each time-step
    do i_proc = 0, Tdomain%n_communications-1
        nv_aus = Tdomain%sWall(i_proc)%n_points
        if (Tdomain%any_PML) then
            nv_aus = nv_aus + 2*Tdomain%sWall(i_proc)%n_points_pml
        endif
        ! Modif of the dimensions of the following arrays by S. Terrana the 22/11/2013
        allocate (Tdomain%sWall(i_proc)%Send_data_2(0:nv_aus-1,0:4))
        allocate (Tdomain%sWall(i_proc)%Receive_data_2(0:nv_aus-1,0:4))
        Tdomain%sWall(i_proc)%Send_data_2    = 0.
        Tdomain%sWall(i_proc)%Receive_data_2 = 0.
        Tdomain%sWall(i_proc)%n_points = nv_aus
    enddo

    ! Invert Mass Matrix expression

    do n  = 0, Tdomain%n_elem - 1
        ngllx = Tdomain%specel(n)%ngllx;  ngllz = Tdomain%specel(n)%ngllz
        if (Tdomain%specel(n)%Type_DG == GALERKIN_CONT) then
            ! Continuous Galerkin : MassMatrices need to be resized
            allocate (LocMassMat(1:ngllx-2,1:ngllz-2))
            ! Redefinition of Matrices
            if ((.not. Tdomain%specel(n)%PML) .or. Tdomain%specel(n)%CPML) then
                LocMassMat(:,:) = Tdomain%specel(n)%MassMat(1:ngllx-2,1:ngllz-2)
                LocMassmat = 1./ LocMassMat
                deallocate (Tdomain%specel(n)%MassMat)
                allocate (Tdomain%specel(n)%MassMat(1:ngllx-2,1:ngllz-2) )
                Tdomain%specel(n)%MassMat =  LocMassMat
            else
                if (.not.Tdomain%specel(n)%CPML) then
                    LocMassMat(:,:) = Tdomain%specel(n)%MassMat(1:ngllx-2,1:ngllz-2)
                    Tdomain%specel(n)%DumpVx (:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:ngllz-2,0)
                    Tdomain%specel(n)%DumpVx (:,:,1) = 1./ Tdomain%specel(n)%DumpVx (:,:,1)
                    Tdomain%specel(n)%DumpVx(:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:ngllz-2,0)
                    Tdomain%specel(n)%DumpVx (:,:,0) =    Tdomain%specel(n)%DumpVx (:,:,0) *   Tdomain%specel(n)%DumpVx (:,:,1)

                    Tdomain%specel(n)%DumpVz (:,:,1) = LocMassMat + Tdomain%specel(n)%DumpMass(1:ngllx-2,1:ngllz-2,1)
                    Tdomain%specel(n)%DumpVz (:,:,1) = 1./ Tdomain%specel(n)%DumpVz (:,:,1)
                    Tdomain%specel(n)%DumpVz(:,:,0) = LocMassMat - Tdomain%specel(n)%DumpMass(1:ngllx-2,1:ngllz-2,1)
                    Tdomain%specel(n)%DumpVz (:,:,0) =    Tdomain%specel(n)%DumpVz (:,:,0) *   Tdomain%specel(n)%DumpVz (:,:,1)
                    deallocate (Tdomain%specel(n)%MassMat) ; deallocate (Tdomain%specel(n)%DumpMass)
                    allocate (Tdomain%specel(n)%MassMat(1:ngllx-2,1:ngllz-2) )
                    Tdomain%specel(n)%MassMat =  LocMassMat
                endif
            endif
            deallocate (LocMassMat)
            ! deallocate (Tdomain%specel(n)%Lamda)
            ! deallocate (Tdomain%specel(n)%Mu)
            !ac a cause des sorties eventuelles des capteurs on garde Tdomain%specel(n)%InvGrad
            !ac if (.not. Tdomain%logicD%save_deformation)  deallocate (Tdomain%specel(n)%InvGrad)
        else
            ! Discontinuous Galerkin Case : Mass Mat do NOT need to be resized
            Tdomain%specel(n)%MassMat = 1. / Tdomain%specel(n)%MassMat
            !Tdomain%specel(n)%Acoeff(:,:,12) = 1. / Tdomain%specel(n)%Acoeff(:,:,12)
        endif
    enddo

    do nf  = 0, Tdomain%n_face - 1
        if (Tdomain%sFace(nf)%PML ) then
            if (Tdomain%sFace(nf)%CPML .OR. Tdomain%sFace(nf)%ADEPML) then
               Tdomain%sFace(nf)%MassMat = 1./ Tdomain%sFace(nf)%MassMat
            else ! Normal PML
                Tdomain%sFace(nf)%DumpVx (:,1) =  Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,0)
                Tdomain%sFace(nf)%DumpVx (:,1) = 1./Tdomain%sFace(nf)%DumpVx (:,1)
                Tdomain%sFace(nf)%DumpVx (:,0) =  Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,0)
                Tdomain%sFace(nf)%DumpVx (:,0) = Tdomain%sFace(nf)%DumpVx (:,0) *  Tdomain%sFace(nf)%DumpVx (:,1)

                Tdomain%sFace(nf)%DumpVz (:,1) =  Tdomain%sFace(nf)%MassMat + Tdomain%sFace(nf)%DumpMass(:,1)
                Tdomain%sFace(nf)%DumpVz (:,1) = 1./Tdomain%sFace(nf)%DumpVz (:,1)
                Tdomain%sFace(nf)%DumpVz (:,0) =  Tdomain%sFace(nf)%MassMat - Tdomain%sFace(nf)%DumpMass(:,1)
                Tdomain%sFace(nf)%DumpVz (:,0) = Tdomain%sFace(nf)%DumpVz (:,0) *  Tdomain%sFace(nf)%DumpVz (:,1)
                deallocate (Tdomain%sFace(nf)%DumpMass)
            endif
        else
            Tdomain%sFace(nf)%MassMat = 1./ Tdomain%sFace(nf)%MassMat
        endif
    enddo

    do nv_aus = 0, Tdomain%n_vertex - 1
        if (Tdomain%sVertex(nv_aus)%PML) then
            if (tdomain%sVertex(nv_aus)%CPML) then
                Tdomain%sVertex(nv_aus)%MassMat = 1./Tdomain%sVertex(nv_aus)%MassMat
            else ! usual PML
                Tdomain%sVertex(nv_aus)%DumpVx (1) =  Tdomain%sVertex(nv_aus)%MassMat + Tdomain%sVertex(nv_aus)%DumpMass(0)
                Tdomain%sVertex(nv_aus)%DumpVx (1) = 1./Tdomain%sVertex(nv_aus)%DumpVx (1)
                Tdomain%sVertex(nv_aus)%DumpVx (0) =  Tdomain%sVertex(nv_aus)%MassMat - Tdomain%sVertex(nv_aus)%DumpMass(0)
                Tdomain%sVertex(nv_aus)%DumpVx (0) = Tdomain%sVertex(nv_aus)%DumpVx (0) *  Tdomain%sVertex(nv_aus)%DumpVx (1)

                Tdomain%sVertex(nv_aus)%DumpVz (1) =  Tdomain%sVertex(nv_aus)%MassMat + Tdomain%sVertex(nv_aus)%DumpMass(1)
                Tdomain%sVertex(nv_aus)%DumpVz (1) = 1./Tdomain%sVertex(nv_aus)%DumpVz (1)
                Tdomain%sVertex(nv_aus)%DumpVz (0) =  Tdomain%sVertex(nv_aus)%MassMat - Tdomain%sVertex(nv_aus)%DumpMass(1)
                Tdomain%sVertex(nv_aus)%DumpVz (0) = Tdomain%sVertex(nv_aus)%DumpVz (0) *  Tdomain%sVertex(nv_aus)%DumpVz (1)
                deallocate (Tdomain%sVertex(nv_aus)%DumpMass)
            endif
        else
            Tdomain%sVertex(nv_aus)%MassMat = 1./Tdomain%sVertex(nv_aus)%MassMat
        endif
    enddo


    ! Define Fault properties
    if (Tdomain%logicD%super_object_local_present) then
        do nf = 0, Tdomain%n_fault-1
            do j = 0, Tdomain%sFault(nf)%n_vertex-1
                Tdomain%sFault(nf)%fVertex(j)%Bt  = 0
            enddo

            ! Internal communication of BT
            do j = 0, Tdomain%sFault(nf)%n_face-1
                ngll = Tdomain%sFault(nf)%fFace(j)%ngll
                mat = Tdomain%sFault(nf)%fFace(j)%mat_index
                if (Tdomain%sFault(nf)%fFace(j)%mat_dir ==0 ) then
                    Tdomain%sFault(nf)%fFace(j)%Bt = Tdomain%sFault(nf)%fFace(j)%ds * Tdomain%sSubdomain(mat)%GLLwx
                else
                    Tdomain%sFault(nf)%fFace(j)%Bt = Tdomain%sFault(nf)%fFace(j)%ds * Tdomain%sSubdomain(mat)%GLLwz
                endif
                nv_aus = Tdomain%sFault(nf)%fFace(j)%Face_to_Vertex(0)
                Tdomain%sFault(nf)%fVertex(nv_aus)%Bt = Tdomain%sFault(nf)%fVertex(nv_aus)%Bt + Tdomain%sFault(nf)%fFace(j)%Bt(0)
                nv_aus = Tdomain%sFault(nf)%fFace(j)%Face_to_Vertex(1)
                Tdomain%sFault(nf)%fVertex(nv_aus)%Bt = Tdomain%sFault(nf)%fVertex(nv_aus)%Bt + Tdomain%sFault(nf)%fFace(j)%Bt(ngll-1)
            enddo

            ! External communication of Bt
            do i_proc = 0, Tdomain%n_communications - 1
                i_send = Tdomain%Communication_list (i_proc)
                if (Tdomain%sWall(i_proc)%n_vertex_superobject > 0) then
                    allocate (Send_Bt(0: Tdomain%sWall(i_proc)%n_vertex_superobject-1))
                    allocate (Receive_Bt(0: Tdomain%sWall(i_proc)%n_vertex_superobject-1))

                    do nv = 0, Tdomain%sWall(i_proc)%n_vertex_superobject - 1
                        nv_aus =  Tdomain%sWall(i_proc)%Vertex_SuperObject_List(nv)
                        Send_Bt (nv)  = Tdomain%sFault(nf)%fVertex(nv_aus)%Bt
                    enddo
                    tag_send = i_send * Tdomain%MPI_var%n_proc +Tdomain%MPI_var%my_rank + 500
                    tag_receive = Tdomain%MPI_var%my_rank * Tdomain%MPI_var%n_proc + i_send + 500
                    call MPI_SEND (Send_Bt, Tdomain%sWall(i_proc)%n_vertex_superobject , MPI_DOUBLE_PRECISION, i_send, &
                        tag_send, Tdomain%communicateur, ierr )
                    call MPI_RECV (Receive_Bt, Tdomain%sWall(i_proc)%n_vertex_superobject, MPI_DOUBLE_PRECISION, i_send, &
                        tag_receive, Tdomain%communicateur, status, ierr )


                    do nv = 0, Tdomain%sWall(i_proc)%n_vertex_superobject - 1
                        nv_aus =  Tdomain%sWall(i_proc)%Vertex_SuperObject_List(nv)
                        Tdomain%sFault(nf)%fVertex(nv_aus)%Bt = Tdomain%sFault(nf)%fVertex(nv_aus)%Bt + Receive_Bt (nv)
                    enddo
                    deallocate (Send_Bt)
                    deallocate (Receive_Bt)
                endif
            enddo

            do j = 0, Tdomain%sFault(nf)%n_face-1
                ngll = Tdomain%sFault(nf)%fFace(j)%ngll
                mat = Tdomain%sFault(nf)%fFace(j)%mat_index
                allocate (LocMassMat1D(1:ngll-2))
                allocate (LocMassMat1D_Down(1:ngll-2))
                nv_aus = Tdomain%sFault(nf)%fFace(j)%Face_UP
                LocMassMat1D = Tdomain%sFace(nv_aus)%MassMat
                nv_aus = Tdomain%sFault(nf)%fFace(j)%Face_DOWN
                if (Tdomain%sFault(nf)%fFace(j)%Coherency) then
                    LocMassMat1D_Down =  Tdomain%sFace(nv_aus)%MassMat
                else
                    do i = 0, ngll-1
                        LocMassMat1D_Down(i) =  Tdomain%sFace(nv_aus)%MassMat(ngll-i-1)
                    enddo
                endif
                LocMassMat1D =  LocMassMat1D + LocMassMat1D_Down
                Tdomain%sFault(nf)%fFace(j)%CGamma= LocMassMat1D * Tdomain%sSubDomain(mat)%Dt * Tdomain%sFault(nf)%fFace(j)%Bt(1:ngll-2)
                deallocate (LocMassMat1D, LocMassMat1D_Down)
            enddo

            do j = 0, Tdomain%sFault(nf)%n_vertex-1
                mat = Tdomain%sFault(nf)%fVertex(j)%mat_index
                nv_aus = Tdomain%sFault(nf)%fVertex(j)%Vertex_UP
                LocMassMat_Vertex = Tdomain%sVertex(nv_aus)%MassMat
                nv_aus = Tdomain%sFault(nf)%fVertex(j)%Vertex_DOWN
                LocMassMat_Vertex =  LocMassMat_Vertex+ Tdomain%sVertex(nv_aus)%MassMat
                Tdomain%sFault(nf)%fVertex(j)%CGamma= LocMassMat_Vertex * Tdomain%sSubDomain(mat)%Dt * Tdomain%sFault(nf)%fVertex(j)%Bt
            enddo
        enddo
    endif

    return
end subroutine define_arrays

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
