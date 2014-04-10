!>
!!\file compute_coeff_flux.F90
!!\brief Contient la subroutine compute_coeff_flux.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Assure le calcul des coefficients pour les Flux Godunov associes a chaque face.
!!
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) n_face
!<

subroutine coeffs_flux_godunov (Tdomain, n_face)
  use sdomain
  implicit none
  type (Domain), intent (INOUT) :: Tdomain
  integer, intent (IN) :: n_face

  ! local variables
  integer ::  ngll, ngllx, ngllz,i
  integer ::  n_elem0, n_elem1, w_face

  n_elem0 = Tdomain%sFace(n_face)%Near_Element(0)
  n_elem1 = Tdomain%sFace(n_face)%Near_Element(1)
  ngll = Tdomain%sFace(n_face)%ngll

  w_face = Tdomain%sFace(n_face)%Which_Face(0)
  ngllx  = Tdomain%specel(n_elem0)%ngllx
  ngllz  = Tdomain%specel(n_elem0)%ngllz

  !! Computing coefficients (impedences) on the - side of the face (i.e. on side of Near_Element(0))
  if (w_face == 0 ) then
     Tdomain%sFace(n_face)%Zp_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0:ngllx-1,0) * &
          sqrt((Tdomain%specel(n_elem0)%Lambda(0:ngllx-1,0) &
          + 2.* Tdomain%specel(n_elem0)%Mu(0:ngllx-1,0)) / Tdomain%specel(n_elem0)%Density(0:ngllx-1,0))
     Tdomain%sFace(n_face)%Zs_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0:ngllx-1,0) * &
          sqrt(Tdomain%specel(n_elem0)%Mu(0:ngllx-1,0) / Tdomain%specel(n_elem0)%Density(0:ngllx-1,0))
  else if (w_face == 1 ) then
     Tdomain%sFace(n_face)%Zp_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(ngllx-1,0:ngllz-1) * &
          sqrt((Tdomain%specel(n_elem0)%Lambda(ngllx-1,0:ngllz-1) &
          + 2.* Tdomain%specel(n_elem0)%Mu(ngllx-1,0:ngllz-1)) / Tdomain%specel(n_elem0)%Density(ngllx-1,0:ngllz-1))
     Tdomain%sFace(n_face)%Zs_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(ngllx-1,0:ngllz-1) * &
          sqrt(Tdomain%specel(n_elem0)%Mu(ngllx-1,0:ngllz-1) / Tdomain%specel(n_elem0)%Density(ngllx-1,0:ngllz-1))
  else if (w_face == 2 ) then
     Tdomain%sFace(n_face)%Zp_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0:ngllx-1,ngllz-1) * &
          sqrt((Tdomain%specel(n_elem0)%Lambda(0:ngllx-1,ngllz-1) &
          + 2.* Tdomain%specel(n_elem0)%Mu(0:ngllx-1,ngllz-1)) / Tdomain%specel(n_elem0)%Density(0:ngllx-1,ngllz-1))
     Tdomain%sFace(n_face)%Zs_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0:ngllx-1,ngllz-1) * &
          sqrt(Tdomain%specel(n_elem0)%Mu(0:ngllx-1,ngllz-1) / Tdomain%specel(n_elem0)%Density(0:ngllx-1,ngllz-1))
  else
     Tdomain%sFace(n_face)%Zp_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0,0:ngllz-1) * &
          sqrt((Tdomain%specel(n_elem0)%Lambda(0,0:ngllz-1) &
          + 2.* Tdomain%specel(n_elem0)%Mu(0,0:ngllz-1)) / Tdomain%specel(n_elem0)%Density(0,0:ngllz-1))
     Tdomain%sFace(n_face)%Zs_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0,0:ngllz-1) * &
          sqrt(Tdomain%specel(n_elem0)%Mu(0,0:ngllz-1) / Tdomain%specel(n_elem0)%Density(0,0:ngllz-1))
  endif

  !! Computing coefficients (impedences) on the + side of the face (i.e. on side of Near_Element(1))
  if( Tdomain%sFace(n_face)%Near_Element(1) > -1) then
     w_face = Tdomain%sFace(n_face)%Which_Face(1)
     ngllx  = Tdomain%specel(n_elem1)%ngllx
     ngllz  = Tdomain%specel(n_elem1)%ngllz
     if (Tdomain%sFace(n_face)%Coherency) then  !! Case COHERENT
        if (w_face == 0 ) then
           Tdomain%sFace(n_face)%Zp_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0:ngllx-1,0) * &
                sqrt((Tdomain%specel(n_elem1)%Lambda(0:ngllx-1,0) &
                + 2.* Tdomain%specel(n_elem1)%Mu(0:ngllx-1,0)) / Tdomain%specel(n_elem1)%Density(0:ngllx-1,0))
           Tdomain%sFace(n_face)%Zs_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0:ngllx-1,0) * &
                sqrt(Tdomain%specel(n_elem1)%Mu(0:ngllx-1,0) / Tdomain%specel(n_elem1)%Density(0:ngllx-1,0))
        else if (w_face == 1 ) then
           Tdomain%sFace(n_face)%Zp_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(ngllx-1,0:ngllz-1) * &
                sqrt((Tdomain%specel(n_elem1)%Lambda(ngllx-1,0:ngllz-1) &
                + 2.* Tdomain%specel(n_elem1)%Mu(ngllx-1,0:ngllz-1)) / Tdomain%specel(n_elem1)%Density(ngllx-1,0:ngllz-1))
           Tdomain%sFace(n_face)%Zs_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(ngllx-1,0:ngllz-1) * &
                sqrt(Tdomain%specel(n_elem1)%Mu(ngllx-1,0:ngllz-1) / Tdomain%specel(n_elem1)%Density(ngllx-1,0:ngllz-1))
        else if (w_face == 2 ) then
           Tdomain%sFace(n_face)%Zp_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0:ngllx-1,ngllz-1) * &
                sqrt((Tdomain%specel(n_elem1)%Lambda(0:ngllx-1,ngllz-1) &
                + 2.* Tdomain%specel(n_elem1)%Mu(0:ngllx-1,ngllz-1)) / Tdomain%specel(n_elem1)%Density(0:ngllx-1,ngllz-1))
           Tdomain%sFace(n_face)%Zs_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0:ngllx-1,ngllz-1) * &
                sqrt(Tdomain%specel(n_elem1)%Mu(0:ngllx-1,ngllz-1) / Tdomain%specel(n_elem1)%Density(0:ngllx-1,ngllz-1))
        else
           Tdomain%sFace(n_face)%Zp_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0,0:ngllz-1) * &
                sqrt((Tdomain%specel(n_elem1)%Lambda(0,0:ngllz-1) &
                + 2.* Tdomain%specel(n_elem1)%Mu(0,0:ngllz-1)) / Tdomain%specel(n_elem1)%Density(0,0:ngllz-1))
           Tdomain%sFace(n_face)%Zs_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0,0:ngllz-1) * &
                sqrt(Tdomain%specel(n_elem1)%Mu(0,0:ngllz-1) / Tdomain%specel(n_elem1)%Density(0,0:ngllz-1))
        endif
     else !! Case NOT COHERENT
        if (w_face == 0 ) then
           do i=0,ngll-1
              Tdomain%sFace(n_face)%Zp_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1-i,0) * &
                   sqrt((Tdomain%specel(n_elem1)%Lambda(ngllx-1-i,0) &
                   + 2.* Tdomain%specel(n_elem1)%Mu(ngllx-1-i,0)) / Tdomain%specel(n_elem1)%Density(ngllx-1-i,0))
              Tdomain%sFace(n_face)%Zs_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1-i,0) * &
                   sqrt(Tdomain%specel(n_elem1)%Mu(ngllx-1-i,0) / Tdomain%specel(n_elem1)%Density(ngllx-1-i,0))
           end do
        else if (w_face == 1 ) then
           do i=0,ngll-1
              Tdomain%sFace(n_face)%Zp_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1,ngllz-1-i) * &
                   sqrt((Tdomain%specel(n_elem1)%Lambda(ngllx-1,ngllz-1-i) &
                   + 2.* Tdomain%specel(n_elem1)%Mu(ngllx-1,ngllz-1-i)) / Tdomain%specel(n_elem1)%Density(ngllx-1,ngllz-1-i))
              Tdomain%sFace(n_face)%Zs_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1,ngllz-1-i) * &
                   sqrt(Tdomain%specel(n_elem1)%Mu(ngllx-1,ngllz-1-i) / Tdomain%specel(n_elem1)%Density(ngllx-1,ngllz-1-i))
           end do
        else if (w_face == 2 ) then
           do i=0,ngll-1
              Tdomain%sFace(n_face)%Zp_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1-i,ngllz-1) * &
                   sqrt((Tdomain%specel(n_elem1)%Lambda(ngllx-1-i,ngllz-1) &
                   + 2.* Tdomain%specel(n_elem1)%Mu(ngllx-1-i,ngllz-1)) / Tdomain%specel(n_elem1)%Density(ngllx-1-i,ngllz-1))
              Tdomain%sFace(n_face)%Zs_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1-i,ngllz-1) * &
                   sqrt(Tdomain%specel(n_elem1)%Mu(ngllx-1-i,ngllz-1) / Tdomain%specel(n_elem1)%Density(ngllx-1-i,ngllz-1))
           end do
        else
           do i=0,ngll-1
              Tdomain%sFace(n_face)%Zp_p(i) = Tdomain%specel(n_elem1)%Density(0,ngllz-1-i) * &
                   sqrt((Tdomain%specel(n_elem1)%Lambda(0,ngllz-1-i) &
                   + 2.* Tdomain%specel(n_elem1)%Mu(0,ngllz-1-i)) / Tdomain%specel(n_elem1)%Density(0,ngllz-1-i))
              Tdomain%sFace(n_face)%Zs_p(i) = Tdomain%specel(n_elem1)%Density(0,ngllz-1-i) * &
                   sqrt(Tdomain%specel(n_elem1)%Mu(0,ngllz-1-i) / Tdomain%specel(n_elem1)%Density(0,ngllz-1-i))
           end do
        endif
     endif
  endif
  !! Computation of the coefficients k0, k1, and vector r1 for Godunov fluxs
  do i=0,ngll-1
     Tdomain%sFace(n_face)%k0(i) = 1./(Tdomain%sFace(n_face)%Zp_m(i) + Tdomain%sFace(n_face)%Zp_p(i))
     Tdomain%sFace(n_face)%k1(i) = 1./(Tdomain%sFace(n_face)%Zs_m(i) + Tdomain%sFace(n_face)%Zs_p(i))
     Tdomain%sFace(n_face)%r1(i,0) = Tdomain%sFace(n_face)%normal(0)**2
     Tdomain%sFace(n_face)%r1(i,1) = Tdomain%sFace(n_face)%normal(1)**2
     Tdomain%sFace(n_face)%r1(i,2) = Tdomain%sFace(n_face)%normal(0) * Tdomain%sFace(n_face)%normal(1)
  end do
  return
end subroutine  coeffs_flux_godunov
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!


!>
!! \brief get_MuLambda_el2f sends the Mu and Lambda coefficients from the elements to the faces
!!  Used only for Godunov-type fluxes
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) nelem
!! \param integer, intent (IN) nface
!<
subroutine get_MuLambdaRho_el2f (Tdomain, nelem, nface)
  use sdomain
  implicit none
  type (Domain), intent (INOUT) :: Tdomain
  integer, intent (IN) :: nelem
  integer, intent (IN) :: nface

  ! local variables
  integer ::  ngll, ngllx, ngllz, w_face, i
  logical :: coherency

  ngll  = Tdomain%sFace(nface)%ngll
  ngllx = Tdomain%specel(nelem)%ngllx
  ngllz = Tdomain%specel(nelem)%ngllz
  coherency = Tdomain%sFace(nface)%coherency

  ! Determinates the w_face value
  if (Tdomain%sFace(nface)%Near_Element(0)==nelem) then
      w_face = Tdomain%sFace(nface)%Which_Face(0)
  else
      w_face = Tdomain%sFace(nface)%Which_Face(1)
  endif

  if (Tdomain%sFace(nface)%Near_Element(0)==nelem) then
     if (w_face == 0 ) then
        Tdomain%sFace(nface)%Mu_m  = Tdomain%specel(nelem)%Mu (0:ngll-1,0)
        Tdomain%sFace(nface)%Lambda_m = Tdomain%specel(nelem)%Lambda(0:ngll-1,0)
        Tdomain%sFace(nface)%Rho_m = Tdomain%specel(nelem)%Density(0:ngll-1,0)
     else if (w_face == 1 ) then
        Tdomain%sFace(nface)%Mu_m  = Tdomain%specel(nelem)%Mu (ngllx-1,0:ngll-1)
        Tdomain%sFace(nface)%Lambda_m = Tdomain%specel(nelem)%Lambda(ngllx-1,0:ngll-1)
        Tdomain%sFace(nface)%Rho_m = Tdomain%specel(nelem)%Density(ngllx-1,0:ngll-1)
     else if (w_face == 2 ) then
        Tdomain%sFace(nface)%Mu_m  = Tdomain%specel(nelem)%Mu (0:ngll-1,ngllz-1)
        Tdomain%sFace(nface)%Lambda_m = Tdomain%specel(nelem)%Lambda(0:ngll-1,ngllz-1)
        Tdomain%sFace(nface)%Rho_m = Tdomain%specel(nelem)%Density(0:ngll-1,ngllz-1)
     else
        Tdomain%sFace(nface)%Mu_m  = Tdomain%specel(nelem)%Mu (0,0:ngll-1)
        Tdomain%sFace(nface)%Lambda_m = Tdomain%specel(nelem)%Lambda(0,0:ngll-1)
        Tdomain%sFace(nface)%Rho_m = Tdomain%specel(nelem)%Density(0,0:ngll-1)
     endif
  else if (coherency) then
     if (w_face == 0 ) then
        Tdomain%sFace(nface)%Mu_p  = Tdomain%specel(nelem)%Mu (0:ngll-1,0)
        Tdomain%sFace(nface)%Lambda_p = Tdomain%specel(nelem)%Lambda(0:ngll-1,0)
        Tdomain%sFace(nface)%Rho_p = Tdomain%specel(nelem)%Density(0:ngll-1,0)
     else if (w_face == 1 ) then
        Tdomain%sFace(nface)%Mu_p  = Tdomain%specel(nelem)%Mu (ngllx-1,0:ngll-1)
        Tdomain%sFace(nface)%Lambda_p = Tdomain%specel(nelem)%Lambda(ngllx-1,0:ngll-1)
        Tdomain%sFace(nface)%Rho_p = Tdomain%specel(nelem)%Density(ngllx-1,0:ngll-1)
     else if (w_face == 2 ) then
        Tdomain%sFace(nface)%Mu_p  = Tdomain%specel(nelem)%Mu (0:ngll-1,ngllz-1)
        Tdomain%sFace(nface)%Lambda_p = Tdomain%specel(nelem)%Lambda(0:ngll-1,ngllz-1)
        Tdomain%sFace(nface)%Rho_p = Tdomain%specel(nelem)%Density(0:ngll-1,ngllz-1)
     else
        Tdomain%sFace(nface)%Mu_p  = Tdomain%specel(nelem)%Mu (0,0:ngll-1)
        Tdomain%sFace(nface)%Lambda_p = Tdomain%specel(nelem)%Lambda(0,0:ngll-1)
        Tdomain%sFace(nface)%Rho_p = Tdomain%specel(nelem)%Density(0,0:ngll-1)
     end if
  else
     if (w_face == 0 ) then
        do i=0,ngll-1
           Tdomain%sFace(nface)%Mu_p(i)  = Tdomain%specel(nelem)%Mu (ngll-1-i,0)
           Tdomain%sFace(nface)%Lambda_p(i) = Tdomain%specel(nelem)%Lambda(ngll-1-i,0)
           Tdomain%sFace(nface)%Rho_p(i) = Tdomain%specel(nelem)%Density(ngll-1-i,0)
        end do
     else if (w_face == 1 ) then
        do i=0,ngll-1
           Tdomain%sFace(nface)%Mu_p(i)  = Tdomain%specel(nelem)%Mu (ngllx-1,ngll-1-i)
           Tdomain%sFace(nface)%Lambda_p(i) = Tdomain%specel(nelem)%Lambda(ngllx-1,ngll-1-i)
           Tdomain%sFace(nface)%Rho_p(i) = Tdomain%specel(nelem)%Density(ngllx-1,ngll-1-i)
        end do
     else if (w_face == 2 ) then
        do i=0,ngll-1
           Tdomain%sFace(nface)%Mu_p(i)  = Tdomain%specel(nelem)%Mu (ngll-1-i,ngllz-1)
           Tdomain%sFace(nface)%Lambda_p(i) = Tdomain%specel(nelem)%Lambda(ngll-1-i,ngllz-1)
           Tdomain%sFace(nface)%Rho_p(i) = Tdomain%specel(nelem)%Density(ngll-1-i,ngllz-1)
        end do
     else
        do i=0,ngll-1
           Tdomain%sFace(nface)%Mu_p(i)  = Tdomain%specel(nelem)%Mu (0,ngll-1-i)
           Tdomain%sFace(nface)%Lambda_p(i) = Tdomain%specel(nelem)%Lambda(0,ngll-1-i)
           Tdomain%sFace(nface)%Rho_p(i) = Tdomain%specel(nelem)%Density(0,ngll-1-i)
        end do
     endif
  end if
  return

end subroutine get_MuLambdaRho_el2f


!>
!! \brief coeff_freesurf extends the Mu, Lambda, mass, and impedences coeficients from
!! the interior side of a face to the exterior side of it. Used for the faces where
!! a "free surface" or an "absorbing" boundary condition is applied.
!!  Used principaly for Godunov-type fluxes
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) nface
!<
subroutine coeffs_freesurf_abs(Tdomain,nface)
  use sdomain
  implicit none
  type (Domain), intent (INOUT) :: Tdomain
  integer, intent (IN) :: nface

  ! local variables
  integer ::  ngll,  i

  Tdomain%sface(nface)%Mu_p     = Tdomain%sface(nface)%Mu_m
  Tdomain%sface(nface)%Lambda_p = Tdomain%sface(nface)%Lambda_m
  Tdomain%sface(nface)%Rho_p    = Tdomain%sface(nface)%Rho_m
  Tdomain%sface(nface)%Zp_p     = Tdomain%sface(nface)%Zp_m
  Tdomain%sface(nface)%Zs_p     = Tdomain%sface(nface)%Zs_m
  ngll = Tdomain%sface(nface)%ngll
  do i=0,ngll-1
     Tdomain%sFace(nface)%k0(i) = 1./(Tdomain%sFace(nface)%Zp_m(i) + Tdomain%sFace(nface)%Zp_p(i))
     Tdomain%sFace(nface)%k1(i) = 1./(Tdomain%sFace(nface)%Zs_m(i) + Tdomain%sFace(nface)%Zs_p(i))
     Tdomain%sFace(nface)%r1(i,0) = Tdomain%sFace(nface)%normal(0)**2
     Tdomain%sFace(nface)%r1(i,1) = Tdomain%sFace(nface)%normal(1)**2
     Tdomain%sFace(nface)%r1(i,2) = Tdomain%sFace(nface)%normal(0) * Tdomain%sFace(nface)%normal(1)
  end do

end subroutine coeffs_freesurf_abs

!>
!! \brief coeff_freesurf extends the Mu, Lambda, mass, from
!! the interior side of a face to its exterior side. Used for the faces where
!! a "free surface" or an "absorbing" boundary condition is applied.
!!  Used principaly for HDG.
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) nface
!<
subroutine coeffs_freesurf_abs_HDG(Tdomain,nface)
    use sdomain
    implicit none
    type (Domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: nface

    if (Tdomain%sface(nface)%Abs) then
        Tdomain%sface(nface)%Mu_p     = Tdomain%sface(nface)%Mu_m
        Tdomain%sface(nface)%Lambda_p = Tdomain%sface(nface)%Lambda_m
        Tdomain%sface(nface)%Rho_p    = Tdomain%sface(nface)%Rho_m
        call compute_InvMatPen(Tdomain%sFace(nface))
    elseif(Tdomain%sface(nface)%freesurf) then
        Tdomain%sface(nface)%Mu_p     = 0.
        Tdomain%sface(nface)%Lambda_p = 0.
        Tdomain%sface(nface)%Rho_p    = 0.
        call compute_InvMatPen(Tdomain%sFace(nface))
    endif

end subroutine coeffs_freesurf_abs_HDG

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
