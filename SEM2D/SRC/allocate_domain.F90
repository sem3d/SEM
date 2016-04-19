!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file allocate_domain.F90
!!\brief Gere l'allocation des domaines.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type(domain), intent (INOUT) Tdomain
!<


subroutine allocate_domain (Tdomain)

  !Modified by Gaetano Festa 01/06/2005

  use sdomain
  use constants
  implicit none

  type(domain), intent (INOUT):: Tdomain

  integer :: n,ngllx,ngllz,ngll,i,j

  do n=0,Tdomain%n_elem-1
     ngllx = Tdomain%specel(n)%ngllx
     ngllz = Tdomain%specel(n)%ngllz
     ngll  = max (ngllx,ngllz)
     allocate (Tdomain%specel(n)%MassMat(0:ngllx-1,0:ngllz-1))
     allocate (Tdomain%specel(n)%Density(0:ngllx-1,0:ngllz-1))
     allocate (Tdomain%specel(n)%Lambda(0:ngllx-1,0:ngllz-1))
     allocate (Tdomain%specel(n)%Mu(0:ngllx-1,0:ngllz-1))
     allocate (Tdomain%specel(n)%Displ(1:ngllx-2,1:ngllz-2,0:1))
     if (Tdomain%specel(n)%type_DG==GALERKIN_CONT) then
         allocate (Tdomain%specel(n)%Veloc (1:ngllx-2, 1:ngllz-2, 0:1))
         allocate (Tdomain%specel(n)%Accel(1:ngllx-2,1:ngllz-2,0:1))
         allocate (Tdomain%specel(n)%V0(1:ngllx-2,1:ngllz-2, 0:1))
         allocate (Tdomain%specel(n)%Forces(0:ngllx-1,0:ngllz-1,0:1))
         Tdomain%specel(n)%Accel = 0
     else ! Discontinuous case
         allocate (Tdomain%specel(n)%Veloc (0:ngllx-1,0:ngllz-1, 0:1))
         allocate (Tdomain%specel(n)%V0(  0:ngllx-1, 0:ngllz-1, 0:1 ) )
         allocate (Tdomain%specel(n)%Forces(0:ngllx-1,0:ngllz-1,0:4))
         allocate (Tdomain%specel(n)%Strain(0:ngllx-1,0:ngllz-1,0:2))
         Tdomain%specel(n)%Strain = 0
         if(Tdomain%specel(n)%Type_DG==GALERKIN_HDG_RP) then
             if (Tdomain%specel(n)%acoustic) then ! Acoustic Element (velocity-pressure)
                 deallocate(Tdomain%specel(n)%Forces,Tdomain%specel(n)%Strain)
                 allocate (Tdomain%specel(n)%Strain(0:ngllx-1,0:ngllz-1,0:0))
                 allocate (Tdomain%specel(n)%Forces(0:ngllx-1,0:ngllz-1,0:2))
                 allocate (Tdomain%specel(n)%MatPen(0:2*(ngllx+ngllz)-1,0:0))
                 allocate (Tdomain%specel(n)%TracFace(0:2*(ngllx+ngllz)-1,0:0))
                 allocate (Tdomain%specel(n)%Vhat(0:2*(ngllx+ngllz)-1,0:0))
                 Tdomain%specel(n)%Strain = 0
             else ! Elastic Element
                 allocate (Tdomain%specel(n)%MatPen(0:2*(ngllx+ngllz)-1,0:2))
                 allocate (Tdomain%specel(n)%TracFace(0:2*(ngllx+ngllz)-1,0:1))
                 allocate (Tdomain%specel(n)%Vhat(0:2*(ngllx+ngllz)-1,0:1))
             endif
                 Tdomain%specel(n)%MatPen   = 0.
                 Tdomain%specel(n)%TracFace = 0.
                 Tdomain%specel(n)%Vhat     = 0.
         endif
     endif
     Tdomain%specel(n)%MassMat = 0
     Tdomain%specel(n)%Density = 0
     Tdomain%specel(n)%Lambda  = 0
     Tdomain%specel(n)%Mu = 0
     Tdomain%specel(n)%Displ = 0.
     Tdomain%specel(n)%Veloc = 0
     Tdomain%specel(n)%Forces = 0
     Tdomain%specel(n)%V0 = 0

     if(Tdomain%specel(n)%CPML .OR. Tdomain%specel(n)%ADEPML) then
         if(Tdomain%specel(n)%ADEPML) then
             allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:12))
         elseif(Tdomain%specel(n)%CPML) then
             allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:17))
         endif
         allocate (Tdomain%specel(n)%Stress (0:ngllx-1, 0:ngllz-1, 0:2))
         allocate (Tdomain%specel(n)%Ax(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Az(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Bx(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Bz(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Ax_prime(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Az_prime(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiVxx (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiVxz (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiVzx (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiVzz (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSxxx(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSzzz(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSxzx(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSxzz(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Psi_store(0:ngllx-1,0:ngllz-1,0:7))
         Tdomain%specel(n)%Stress   = 0.
         Tdomain%specel(n)%Acoeff   = 0.
         Tdomain%specel(n)%Ax       = 0.
         Tdomain%specel(n)%Az       = 0.
         Tdomain%specel(n)%Bx       = 0.
         Tdomain%specel(n)%Bz       = 0.
         Tdomain%specel(n)%Ax_prime = 0.
         Tdomain%specel(n)%Az_prime = 0.
         Tdomain%specel(n)%PsiVxx   = 0.
         Tdomain%specel(n)%PsiVxz   = 0.
         Tdomain%specel(n)%PsiVzx   = 0.
         Tdomain%specel(n)%PsiVzz   = 0.
         Tdomain%specel(n)%PsiSxxx  = 0.
         Tdomain%specel(n)%PsiSzzz  = 0.
         Tdomain%specel(n)%PsiSxzx  = 0.
         Tdomain%specel(n)%PsiSxzz  = 0.
         Tdomain%specel(n)%Psi_store = 0.
     elseif (Tdomain%specel(n)%PML ) then
         allocate (Tdomain%specel(n)%Stress (  0:ngllx-1, 0:ngllz-1, 0:2 ))
         allocate (Tdomain%specel(n)%Forces1 (0:ngllx-1,0:ngllz-1,0:1) )
         allocate (Tdomain%specel(n)%Forces2 (0:ngllx-1,0:ngllz-1,0:1) )
         allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:15))
         allocate (Tdomain%specel(n)%Stress1 (0:ngllx-1,0:ngllz-1,0:2) )
         allocate (Tdomain%specel(n)%Stress2 (0:ngllx-1,0:ngllz-1,0:2) )
         allocate (Tdomain%specel(n)%Veloc1 (1:ngllx-2,1:ngllz-2,0:1) )
         allocate (Tdomain%specel(n)%Veloc2 (1:ngllx-2,1:ngllz-2,0:1) )
         allocate (Tdomain%specel(n)%DumpSx (0:ngllx-1,0:ngllz-1,0:1) )
         allocate (Tdomain%specel(n)%DumpSz (0:ngllx-1,0:ngllz-1,0:1) )
         allocate (Tdomain%specel(n)%DumpVx (1:ngllx-2,1:ngllz-2,0:1) )
         allocate (Tdomain%specel(n)%DumpVz (1:ngllx-2,1:ngllz-2,0:1) )
         allocate (Tdomain%specel(n)%DumpMass(0:ngllx-1,0:ngllz-1,0:1))
         Tdomain%specel(n)%Stress = 0.
         Tdomain%specel(n)%Veloc1 = 0.
         Tdomain%specel(n)%Veloc2 = 0.
         Tdomain%specel(n)%Stress1 = 0.
         Tdomain%specel(n)%Stress2 = 0.
         Tdomain%specel(n)%DumpMass= 0.
     else ! Case Element is not PML
         if(Tdomain%specel(n)%Type_DG==GALERKIN_CONT) then
             allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:15))
         else ! Discontinuous Case
             if (Tdomain%specel(n)%Type_DG==GALERKIN_HDG_RP .AND. Tdomain%specel(n)%acoustic) then
                 allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:4))
             else ! Discontinuous Galerkin, usual
                 allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:12))
             endif
         endif
         Tdomain%specel(n)%Acoeff = 0.
     endif
 enddo

  do n = 0, Tdomain%n_face-1
     Tdomain%sFace(n)%is_computed = .false.
     ngll = Tdomain%sFace(n)%ngll
#ifdef MKA3D
     allocate (Tdomain%sFace(n)%ForcesMka(1:ngll-2,0:1 ) )
     Tdomain%sFace(n)%ForcesMka=0.
#endif
     allocate (Tdomain%sFace(n)%MassMat(1:ngll-2))
     allocate (Tdomain%sFace(n)%Veloc (1:ngll-2, 0:1 ) )
     allocate (Tdomain%sFace(n)%Forces(1:ngll-2,0:1 ) )
     allocate (Tdomain%sFace(n)%Accel(1:ngll-2, 0:1))
     allocate (Tdomain%sFace(n)%V0( 1:ngll-2, 0:1 ) )
     allocate (Tdomain%sFace(n)%Displ( 1:ngll-2, 0:1 ) )
     Tdomain%sFace(n)%MassMat = 0
     Tdomain%sFace(n)%Veloc = 0
     Tdomain%sFace(n)%Accel = 0
     Tdomain%sFace(n)%V0 = 0
     Tdomain%sFace(n)%Displ = 0
     Tdomain%sFace(n)%Forces= 0
     if (Tdomain%sFace(n)%type_Flux .EQ. FLUX_GODUNOV ) then !
         ! Allocation coefficients for Godunov Fluxes
         allocate (Tdomain%sFace(n)%k0(0:ngll-1))
         allocate (Tdomain%sFace(n)%k1(0:ngll-1))
         allocate (Tdomain%sFace(n)%r1(0:ngll-1,0:4))
         allocate (Tdomain%sFace(n)%r2(0:ngll-1,0:4))
         allocate (Tdomain%sFace(n)%r3(0:ngll-1,0:4))
         allocate (Tdomain%sFace(n)%Zp_m(0:ngll-1))
         allocate (Tdomain%sFace(n)%Zp_p(0:ngll-1))
         allocate (Tdomain%sFace(n)%Zs_m(0:ngll-1))
         allocate (Tdomain%sFace(n)%Zs_p(0:ngll-1))
         allocate (Tdomain%sFace(n)%Flux(0:ngll-1, 0:4 ) )
         allocate (Tdomain%sFace(n)%Veloc_p(0:ngll-1, 0:1))
         allocate (Tdomain%sFace(n)%Veloc_m(0:ngll-1, 0:1))
         allocate (Tdomain%sFace(n)%Strain_p(0:ngll-1, 0:2))
         allocate (Tdomain%sFace(n)%Strain_m(0:ngll-1, 0:2))
         allocate (Tdomain%sFace(n)%Mu_p(0:ngll-1))
         allocate (Tdomain%sFace(n)%Mu_m(0:ngll-1))
         allocate (Tdomain%sFace(n)%Lambda_p(0:ngll-1))
         allocate (Tdomain%sFace(n)%Lambda_m(0:ngll-1))
         allocate (Tdomain%sFace(n)%Rho_p(0:ngll-1))
         allocate (Tdomain%sFace(n)%Rho_m(0:ngll-1))
         allocate (Tdomain%sFace(n)%Flux_p(0:ngll-1,0:4))
         Tdomain%sFace(n)%k0 = 0
         Tdomain%sFace(n)%k1 = 0
         Tdomain%sFace(n)%r1 = 0
         Tdomain%sFace(n)%r2 = 0
         Tdomain%sFace(n)%r3 = 0
         Tdomain%sFace(n)%Zp_m = 0
         Tdomain%sFace(n)%Zp_p = 0
         Tdomain%sFace(n)%Zs_m = 0
         Tdomain%sFace(n)%Zs_p = 0
         Tdomain%sFace(n)%Flux  = 0
         Tdomain%sFace(n)%Veloc_p = 0
         Tdomain%sFace(n)%Veloc_m = 0
         Tdomain%sFace(n)%Strain_p = 0
         Tdomain%sFace(n)%Strain_m = 0
         Tdomain%sFace(n)%Mu_p = 0
         Tdomain%sFace(n)%Mu_m = 0
         Tdomain%sFace(n)%Lambda_p = 0
         Tdomain%sFace(n)%Lambda_m = 0
         Tdomain%sFace(n)%Rho_p = 0
         Tdomain%sFace(n)%Rho_m = 0
         Tdomain%sFace(n)%Flux_p = 0
     elseif (Tdomain%sFace(n)%type_Flux .EQ. FLUX_HDG ) then
         deallocate(Tdomain%sFace(n)%Veloc,Tdomain%sFace(n)%V0,Tdomain%sFace(n)%Accel)
         if (Tdomain%sFace(n)%acoustic) then ! ACOUSTIC CASE
             allocate (Tdomain%sFace(n)%Veloc(0:ngll-1,0:0))
             allocate (Tdomain%sFace(n)%KinvExpl(0:ngll-1,0:0))
             allocate (Tdomain%sFace(n)%SmbrTrac(0:ngll-1,0:0))
             allocate (Tdomain%sFace(n)%V0(0:ngll-1,0:0))
         else ! ELASTIC CASE
             allocate (Tdomain%sFace(n)%Veloc(0:ngll-1,0:1))
             allocate (Tdomain%sFace(n)%KinvExpl(0:ngll-1,0:2))
             allocate (Tdomain%sFace(n)%SmbrTrac(0:ngll-1,0:1))
             allocate (Tdomain%sFace(n)%V0(0:ngll-1,0:1))
         endif
         Tdomain%sFace(n)%Veloc = 0
         Tdomain%sFace(n)%KinvExpl = 0
         Tdomain%sFace(n)%SmbrTrac = 0
         Tdomain%sFace(n)%V0 = 0
     endif

     if (Tdomain%sFace(n)%PML .AND. (.NOT.Tdomain%sFace(n)%CPML)) then
         allocate (Tdomain%sFace(n)%Forces1 (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%Forces2 (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%Veloc1 (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%Veloc2 (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%DumpVx (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%DumpVz (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%DumpMass(1:ngll-2,0:1))
         Tdomain%sFace(n)%Veloc1 = 0.
         Tdomain%sFace(n)%Veloc2 = 0.
         Tdomain%sFace(n)%DumpMass = 0.
     endif

     ! Putting a flag for Continuous-Discontinuous interface :
     i = Tdomain%sFace(n)%Near_Element(0)
     j = Tdomain%sFace(n)%Near_Element(1)
     if ( (j.NE.-1) .AND. (Tdomain%specel(i)%type_DG .NE. Tdomain%specel(j)%type_DG)) then
        write(*,*) "Changing CG-HDG for face : ",n
        deallocate(Tdomain%sFace(n)%Veloc,Tdomain%sFace(n)%Forces)
        allocate(Tdomain%sFace(n)%Veloc (0:ngll-1,0:1))
        allocate(Tdomain%sFace(n)%Forces(0:ngll-1,0:1))
        if (.NOT. allocated(Tdomain%sFace(n)%Accel)) allocate(Tdomain%sFace(n)%Accel(0:ngll-1,0:1))
        if (.NOT. allocated(Tdomain%sFace(n)%KinvExpl)) allocate(Tdomain%sFace(n)%KinvExpl(0:ngll-1,0:2))
        if (.NOT. allocated(Tdomain%sFace(n)%SmbrTrac)) allocate(Tdomain%sFace(n)%SmbrTrac(0:ngll-1,0:1))
        Tdomain%type_elem = COUPLE_CG_HDG
        Tdomain%sFace(n)%Veloc = 0.
        Tdomain%sFace(n)%Forces = 0.
        Tdomain%sFace(n)%KinvExpl = 0.
        Tdomain%sFace(n)%SmbrTrac = 0.
        Tdomain%sFace(n)%type_DG = COUPLE_CG_HDG
        Tdomain%sFace(n)%type_Flux = FLUX_HDG
        i = Tdomain%sFace(n)%Near_Vertex(0)
        j = Tdomain%sFace(n)%Near_Vertex(1)
        Tdomain%sVertex(i)%Type_DG = GALERKIN_CONT
        Tdomain%sVertex(j)%Type_DG = GALERKIN_CONT
        ! Coefficients of integration on face
        !j = Tdomain%sFace(n)%Near_Element(0)
        !i = Tdomain%sFace(n)%Which_Face(0)
        !allocate(Tdomain%sFace(n)%Coeff_Integr(0:ngll-1))
        !call get_iminimax(Tdomain%specel(j),i,imin,imax)
        !Tdomain%sFace(n)%Coeff_Integr(:) = Tdomain%specel(j)%Coeff_Integr_Faces(imin:imax)
     endif
  enddo

  do n = 0, Tdomain%n_vertex-1
     Tdomain%sVertex(n)%is_computed = .false.
     allocate (Tdomain%sVertex(n)%Veloc (0:1))
     allocate (Tdomain%sVertex(n)%Forces (0:1))
     allocate (Tdomain%sVertex(n)%Accel (0:1))
     allocate (Tdomain%sVertex(n)%V0 (0:1))
     allocate (Tdomain%sVertex(n)%Displ(0:1))
#ifdef MKA3D
     allocate (Tdomain%sVertex(n)%ForcesMka(0:1 ) )
     Tdomain%sVertex(n)%ForcesMka=0
#endif

     Tdomain%sVertex(n)%MassMat = 0
     Tdomain%sVertex(n)%Veloc = 0
     Tdomain%sVertex(n)%Accel = 0
     Tdomain%sVertex(n)%V0 = 0
     Tdomain%sVertex(n)%Displ = 0

     if (Tdomain%sVertex(n)%PML .AND. (.NOT.Tdomain%sFace(n)%CPML)) then
         allocate (Tdomain%sVertex(n)%Forces1 (0:1) )
         allocate (Tdomain%sVertex(n)%Forces2 (0:1) )
         allocate (Tdomain%sVertex(n)%Veloc1 (0:1) )
         allocate (Tdomain%sVertex(n)%Veloc2 (0:1) )
         allocate (Tdomain%sVertex(n)%DumpVx (0:1) )
         allocate (Tdomain%sVertex(n)%DumpVz (0:1) )
         allocate (Tdomain%sVertex(n)%DumpMass(0:1))
         Tdomain%sVertex(n)%Veloc1 = 0
         Tdomain%sVertex(n)%Veloc2 = 0
         Tdomain%sVertex(n)%DumpMass = 0
     endif
 enddo

  ! Init for Runge-Kutta Low Storage time integration.
  if(Tdomain%type_TimeInteg==TIME_INTEG_RK4) then
     do n = 0, Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        if (Tdomain%specel(n)%Type_DG == GALERKIN_CONT) then
           allocate (Tdomain%specel(n)%Vect_RK(1:ngllx-2,1:ngllz-2,0:3))
        else
           if (Tdomain%specel(n)%acoustic) then
               allocate (Tdomain%specel(n)%Vect_RK(0:ngllx-1,0:ngllz-1,0:2))
           else ! Elastic
               allocate (Tdomain%specel(n)%Vect_RK(0:ngllx-1,0:ngllz-1,0:4))
           endif
        endif
        Tdomain%specel(n)%Vect_RK = 0.
     enddo
     do n = 0, Tdomain%n_face-1
        ngll = Tdomain%sFace(n)%ngll
        if (Tdomain%sface(n)%Type_DG == GALERKIN_CONT .OR. &
            Tdomain%sface(n)%Type_DG == COUPLE_CG_HDG) then
           allocate (Tdomain%sface(n)%Vect_RK(1:ngll-2,0:3))
           Tdomain%sface(n)%Vect_RK = 0.
        endif
     enddo
     do n = 0, Tdomain%n_vertex-1
        if (Tdomain%svertex(n)%Type_DG == GALERKIN_CONT) then
           allocate (Tdomain%sVertex(n)%Vect_RK(0:3))
           Tdomain%svertex(n)%Vect_RK = 0.
        endif
     enddo
  else if (Tdomain%Implicitness==TIME_INTEG_SEMI_IMPLICIT) then
     do n = 0, Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        allocate (Tdomain%specel(n)%Strain0(0:ngllx-1,0:ngllz-1,0:2))
        allocate (Tdomain%specel(n)%CAinv(0:2*(ngllx+ngllz)-1,0:1,0:2))
        allocate (Tdomain%specel(n)%EDinv(0:2*(ngllx+ngllz)-1,0:1,0:1))
        allocate (Tdomain%specel(n)%Dinv (0:2*(ngllx+ngllz)-1,0:2))
     enddo
     do n = 0, Tdomain%n_face-1
         ngll = Tdomain%sFace(n)%ngll
         allocate (Tdomain%sFace(n)%Kinv(0:ngll-1,0:2))
         Tdomain%sFace(n)%Kinv = 0.
     enddo
     do n = 0, Tdomain%n_vertex-1
        i = Tdomain%sVertex(n)%Valence
        allocate (Tdomain%sVertex(n)%Kmat(0:(2*i-1),0:(2*i-1)))
        allocate (Tdomain%sVertex(n)%smbrLambda(0:(2*i-1)))
        allocate (Tdomain%sVertex(n)%Lambda (0:(2*i-1)))
        allocate (Tdomain%sVertex(n)%K_up (1:2*i*(2*i+1)/2))
        Tdomain%sVertex(n)%smbrLambda = 0.
        Tdomain%sVertex(n)%Kmat = 0.
        Tdomain%sVertex(n)%K_up = 0.
     enddo
  else if (Tdomain%type_timeInteg==TIME_INTEG_MIDPOINT .OR. &
           Tdomain%type_timeInteg==TIME_INTEG_MIDPOINT_ITER ) then
     do n = 0, Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        allocate (Tdomain%specel(n)%Strain0(0:ngllx-1,0:ngllz-1,0:2))
    enddo
  endif

  ! Special addition for Lamb test :
  if (Tdomain%LogicD%Lamb_test) then
      do n=0,Tdomain%n_face-1
          i = Tdomain%sFace(n)%Near_Vertex(0)
          j = Tdomain%sFace(n)%Near_Vertex(1)
          i = Tdomain%sVertex(i)%Glob_numbering
          j = Tdomain%sVertex(j)%Glob_numbering
          if (abs(Tdomain%coord_nodes(1,i)) .LT. 1.E-7 .and. abs(Tdomain%coord_nodes(1,j)) .LT. 1.E-7 &
              .and. Tdomain%sFace(n)%Reflex ) then
              Tdomain%sFace(n)%freesurf = .true.
              Tdomain%sFace(n)%Abs      = .false.
              Tdomain%sFace(n)%Reflex   = .false.
              i = Tdomain%sFace(n)%Near_Vertex(0)
              j = Tdomain%sFace(n)%Near_Vertex(1)
              Tdomain%sVertex(i)%Abs      = .false.
              Tdomain%sVertex(i)%Reflex   = .false.
              Tdomain%sVertex(j)%Abs      = .false.
              Tdomain%sVertex(j)%Reflex   = .false.
          endif
      enddo
      do n=0,Tdomain%n_face-1  ! Traitement pour les coins
          if (Tdomain%sFace(n)%Reflex) then
              i = Tdomain%sFace(n)%Near_Vertex(0)
              j = Tdomain%sFace(n)%Near_Vertex(1)
              Tdomain%sVertex(i)%Reflex   = .true.
              Tdomain%sVertex(j)%Reflex   = .true.
          endif
      enddo
  endif

  return
end subroutine allocate_domain

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
