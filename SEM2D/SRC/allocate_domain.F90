!>
!!\file allocate_domain.F90
!!\brief Gère l'allocation des domaines.
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

  integer :: n,ngllx,ngllz,ngll

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
     else ! Discontinuous case
         allocate (Tdomain%specel(n)%Veloc (0:ngllx-1,0:ngllz-1, 0:1))
         allocate (Tdomain%specel(n)%Accel(0:ngllx-1,0:ngllz-1,0:1))
         allocate (Tdomain%specel(n)%V0(  0:ngllx-1, 0:ngllz-1, 0:1 ) )
         allocate (Tdomain%specel(n)%Forces(0:ngllx-1,0:ngllz-1,0:4))
         allocate (Tdomain%specel(n)%Strain(0:ngllx-1,0:ngllz-1,0:2))
         if(Tdomain%specel(n)%Type_DG==GALERKIN_HDG_RP) then
             allocate (Tdomain%specel(n)%MatPen(0:2*(ngllx+ngllz)-1,0:2))
             allocate (Tdomain%specel(n)%TracFace(0:2*(ngllx+ngllz)-1,0:1))
             allocate (Tdomain%specel(n)%Vhat(0:2*(ngllx+ngllz)-1,0:1))
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
     Tdomain%specel(n)%Accel = 0
     Tdomain%specel(n)%V0 = 0
     Tdomain%specel(n)%Strain = 0

     if(Tdomain%specel(n)%CPML .OR. Tdomain%specel(n)%ADEPML) then
         if(Tdomain%specel(n)%ADEPML) then
             allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:12))
         elseif(Tdomain%specel(n)%CPML) then
             allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:17))
         endif
         allocate (Tdomain%specel(n)%Stress (0:ngllx-1, 0:ngllz-1, 0:2))
         allocate (Tdomain%specel(n)%Axi (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Aeta(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Bxi (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Beta(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Axi_prime (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%Aeta_prime(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiVxxi (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiVxeta(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiVzxi (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiVzeta(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSxxxi (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSxxeta(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSzzxi (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSzzeta(0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSxzxi (0:ngllx-1,0:ngllz-1))
         allocate (Tdomain%specel(n)%PsiSxzeta(0:ngllx-1,0:ngllz-1))
         Tdomain%specel(n)%Stress    = 0.
         Tdomain%specel(n)%Acoeff    = 0.
         Tdomain%specel(n)%Axi       = 0.
         Tdomain%specel(n)%Aeta      = 0.
         Tdomain%specel(n)%Bxi       = 0.
         Tdomain%specel(n)%Beta      = 0.
         Tdomain%specel(n)%Axi_prime = 0.
         Tdomain%specel(n)%Aeta_prime= 0.
         Tdomain%specel(n)%PsiVxxi   = 0.
         Tdomain%specel(n)%PsiVxeta  = 0.
         Tdomain%specel(n)%PsiVzxi   = 0.
         Tdomain%specel(n)%PsiVzeta  = 0.
         Tdomain%specel(n)%PsiSxxxi  = 0.
         Tdomain%specel(n)%PsiSxxeta = 0.
         Tdomain%specel(n)%PsiSzzxi  = 0.
         Tdomain%specel(n)%PsiSzzeta = 0.
         Tdomain%specel(n)%PsiSxzxi  = 0.
         Tdomain%specel(n)%PsiSxzeta = 0.
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
         Tdomain%specel(n)%Stress = 0.
         Tdomain%specel(n)%Veloc1 = 0.
         Tdomain%specel(n)%Veloc2 = 0.
         Tdomain%specel(n)%Stress1 = 0.
         Tdomain%specel(n)%Stress2 = 0.

         if (Tdomain%specel(n)%FPML) then
             allocate (Tdomain%specel(n)%Isx(0:ngllx-1,0:ngllz-1))
             allocate (Tdomain%specel(n)%Isz(0:ngllx-1,0:ngllz-1))
             allocate (Tdomain%specel(n)%Ivx(0:ngllx-1,0:ngllz-1))
             allocate (Tdomain%specel(n)%Ivz(0:ngllx-1,0:ngllz-1))
             allocate (Tdomain%specel(n)%IStress1(0:ngllx-1,0:ngllz-1,0:2))
             allocate (Tdomain%specel(n)%IStress2(0:ngllx-1,0:ngllz-1,0:2))
             allocate (Tdomain%specel(n)%IVeloc1(1:ngllx-2,1:ngllz-2,0:1))
             allocate (Tdomain%specel(n)%IVeloc2(1:ngllx-2,1:ngllz-2,0:1))
             allocate (Tdomain%specel(n)%DumpMass(0:ngllx-1,0:ngllz-1,0:3))
             Tdomain%specel(n)%IStress1 = 0.
             Tdomain%specel(n)%IStress2 = 0.
             Tdomain%specel(n)%IVeloc1 = 0.
             Tdomain%specel(n)%IVeloc2 = 0.
         else
             allocate (Tdomain%specel(n)%DumpMass(0:ngllx-1,0:ngllz-1,0:1))
         endif
     else ! Case Element is not PML
         if(Tdomain%specel(n)%Type_DG==GALERKIN_CONT) then
             allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:9))
         else ! Discontinuous Case
             allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:12))
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
     Tdomain%sFace(n)%MassMat = 0
     Tdomain%sFace(n)%Veloc = 0
     Tdomain%sFace(n)%Accel = 0
     Tdomain%sFace(n)%V0 = 0
     Tdomain%sFace(n)%Displ = 0
     Tdomain%sFace(n)%Forces= 0
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
         Tdomain%sFace(n)%k0 = 0
         Tdomain%sFace(n)%k1 = 0
         Tdomain%sFace(n)%r1 = 0
         Tdomain%sFace(n)%r2 = 0
         Tdomain%sFace(n)%r3 = 0
         Tdomain%sFace(n)%Zp_m = 0
         Tdomain%sFace(n)%Zp_p = 0
         Tdomain%sFace(n)%Zs_m = 0
         Tdomain%sFace(n)%Zs_p = 0
     elseif (Tdomain%sFace(n)%type_Flux .EQ. FLUX_HDG ) then
         deallocate(Tdomain%sFace(n)%Veloc)
         allocate (Tdomain%sFace(n)%Veloc(0:ngll-1,0:1))
         allocate (Tdomain%sFace(n)%invMatPen(0:ngll-1,0:2))
         allocate (Tdomain%sFace(n)%Traction(0:ngll-1,0:1))
         Tdomain%sFace(n)%invMatPen = 0
         Tdomain%sFace(n)%Traction  = 0
         Tdomain%sFace(n)%Veloc = 0
     endif

     if (Tdomain%sFace(n)%PML .AND. (.NOT.Tdomain%sFace(n)%CPML)) then
         allocate (Tdomain%sFace(n)%Forces1 (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%Forces2 (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%Veloc1 (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%Veloc2 (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%DumpVx (1:ngll-2,0:1) )
         allocate (Tdomain%sFace(n)%DumpVz (1:ngll-2,0:1) )

         Tdomain%sFace(n)%Veloc1 = 0.
         Tdomain%sFace(n)%Veloc2 = 0.

         if (Tdomain%sFace(n)%FPML) then
             allocate (Tdomain%sFace(n)%Ivx(1:ngll-2))
             allocate (Tdomain%sFace(n)%Ivz(1:ngll-2))
             allocate (Tdomain%sFace(n)%IVeloc1(1:ngll-2,0:1))
             allocate (Tdomain%sFace(n)%IVeloc2(1:ngll-2,0:1))
             allocate (Tdomain%sFace(n)%DumpMass(1:ngll-2,0:3))
             Tdomain%sFace(n)%IVeloc1 = 0.
             Tdomain%sFace(n)%IVeloc2 = 0.
             Tdomain%sFace(n)%Ivx = 0.
             Tdomain%sFace(n)%Ivz = 0.
         else
             allocate (Tdomain%sFace(n)%DumpMass(1:ngll-2,0:1))
         endif
         Tdomain%sFace(n)%DumpMass = 0.
     endif

     ! Flags on absorbing/free surface Faces
     if (Tdomain%sFace(n)%Near_Element(1) == -1) then
         if (Tdomain%type_bc == DG_BC_ABS) then
             Tdomain%sFace(n)%abs      = .true.
             Tdomain%sFace(n)%freesurf = .false.
         elseif (Tdomain%type_bc == DG_BC_FREE) then
             Tdomain%sFace(n)%abs      = .false.
             Tdomain%sFace(n)%freesurf = .true.
         endif
     else
         Tdomain%sFace(n)%abs = .false.
         Tdomain%sFace(n)%freesurf = .false.
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
         Tdomain%sVertex(n)%Veloc1 = 0
         Tdomain%sVertex(n)%Veloc2 = 0

         if (Tdomain%sFace(n)%FPML) then
             allocate (Tdomain%sVertex(n)%Ivx(0:0))
             allocate (Tdomain%sVertex(n)%Ivz(0:0))
             allocate (Tdomain%sVertex(n)%IVeloc1(0:1))
             allocate (Tdomain%sVertex(n)%IVeloc2(0:1))
             allocate (Tdomain%sVertex(n)%DumpMass(0:1))
             Tdomain%sVertex(n)%IVeloc1 = 0.
             Tdomain%sVertex(n)%IVeloc2 = 0.
         else
             allocate (Tdomain%sVertex(n)%DumpMass(0:1))
         endif
         Tdomain%sVertex(n)%DumpMass = 0
     endif
 enddo

  ! Init for Runge-Kutta Low Storage time integration.
  if(Tdomain%type_TimeInteg==TIME_INTEG_RK4) then
     do n = 0, Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        if (Tdomain%specel(n)%Type_DG == GALERKIN_CONT) then
           allocate (Tdomain%specel(n)%Vect_RK(0:ngllx-1,0:ngllz-1,0:1))
        else
           allocate (Tdomain%specel(n)%Vect_RK(0:ngllx-1,0:ngllz-1,0:4))
        endif
        Tdomain%specel(n)%Vect_RK = 0.
        if (Tdomain%specel(n)%ADEPML) then
            allocate (Tdomain%specel(n)%Psi_RK(0:ngllx-1,0:ngllz-1,0:9))
            Tdomain%specel(n)%Psi_RK = 0.
        endif
     enddo
     do n = 0, Tdomain%n_face-1
        ngll = Tdomain%sFace(n)%ngll
        if (Tdomain%sface(n)%Type_DG == GALERKIN_CONT) then
           allocate (Tdomain%sface(n)%Vect_RK(0:ngll-1,0:1))
        else
           allocate (Tdomain%sface(n)%Vect_RK(0:ngll-1,0:4))
        endif
        Tdomain%sface(n)%Vect_RK = 0.
     enddo
     do n = 0, Tdomain%n_vertex-1
        if (Tdomain%svertex(n)%Type_DG == GALERKIN_CONT) then
           allocate (Tdomain%sVertex(n)%Vect_RK(0:1))
        else
           allocate (Tdomain%sVertex(n)%Vect_RK(0:4))
        endif
        Tdomain%svertex(n)%Vect_RK = 0.
     enddo
  endif

  !! Special addition for Lamb test : A SUPPRIMER !!!!!!!!
!  do n=0,Tdomain%n_face-1
!      i = Tdomain%sFace(n)%Near_Vertex(0)
!      j = Tdomain%sFace(n)%Near_Vertex(1)
!      i = Tdomain%sVertex(i)%Glob_numbering
!      j = Tdomain%sVertex(j)%Glob_numbering
!      if (Tdomain%coord_nodes(1,i)==0. .and. Tdomain%coord_nodes(1,j)==0. &
!          .and. Tdomain%sFace(n)%Abs ) then
!          Tdomain%sFace(n)%freesurf = .true.
!          Tdomain%sFace(n)%Abs      = .false.
!      endif
!  enddo
  !!! FIN A SUPPRIMER !!!!!!!!

  return
end subroutine allocate_domain
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
