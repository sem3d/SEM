!>
!!\file compute_coeff_flux.F90
!!\brief Contient la subroutine compute_coeff_flux.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Assure le calcul des coefficients pour les Flux Godunov associées à chaque face.
!!
!! \param type (Domain), intent (INOUT) Tdomain
!! \param integer, intent (IN) n_face
!<


subroutine compute_coeff_flux (Tdomain, n_face)
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
                     sqrt(Tdomain%specel(n_elem0)%Mu(0:ngllx-1,0)) / Tdomain%specel(n_elem0)%Density(0:ngllx-1,0))
    else if (w_face == 1 ) then
       Tdomain%sFace(n_face)%Zp_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(ngllx-1,0:ngllz-1) * &
                     sqrt((Tdomain%specel(n_elem0)%Lambda(ngllx-1,0:ngllz-1) &
                     + 2.* Tdomain%specel(n_elem0)%Mu(ngllx-1,0:ngllz-1)) / Tdomain%specel(n_elem0)%Density(ngllx-1,0:ngllz-1))
       Tdomain%sFace(n_face)%Zs_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(ngllx-1,0:ngllz-1) * &
                     sqrt(Tdomain%specel(n_elem0)%Mu(ngllx-1,0:ngllz-1)) / Tdomain%specel(n_elem0)%Density(ngllx-1,0:ngllz-1))
    else if (w_face == 2 ) then
       Tdomain%sFace(n_face)%Zp_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0:ngllx-1,ngllz-1) * &
                     sqrt((Tdomain%specel(n_elem0)%Lambda(0:ngllx-1,ngllz-1) &
                     + 2.* Tdomain%specel(n_elem0)%Mu(0:ngllx-1,ngllz-1)) / Tdomain%specel(n_elem0)%Density(0:ngllx-1,ngllz-1))
       Tdomain%sFace(n_face)%Zs_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0:ngllx-1,ngllz-1) * &
                     sqrt(Tdomain%specel(n_elem0)%Mu(0:ngllx-1,ngllz-1)) / Tdomain%specel(n_elem0)%Density(0:ngllx-1,ngllz-1))
    else
       Tdomain%sFace(n_face)%Zp_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0,0:ngllz-1) * &
                     sqrt((Tdomain%specel(n_elem0)%Lambda(0,0:ngllz-1) &
                     + 2.* Tdomain%specel(n_elem0)%Mu(0,0:ngllz-1)) / Tdomain%specel(n_elem0)%Density(0,0:ngllz-1))
       Tdomain%sFace(n_face)%Zs_m(0:ngll-1) = Tdomain%specel(n_elem0)%Density(0,0:ngllz-1) * &
                     sqrt(Tdomain%specel(n_elem0)%Mu(0,0:ngllz-1)) / Tdomain%specel(n_elem0)%Density(0,0:ngllz-1))
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
                     sqrt(Tdomain%specel(n_elem1)%Mu(0:ngllx-1,0)) / Tdomain%specel(n_elem1)%Density(0:ngllx-1,0))
          else if (w_face == 1 ) then
             Tdomain%sFace(n_face)%Zp_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(ngllx-1,0:ngllz-1) * &
                     sqrt((Tdomain%specel(n_elem1)%Lambda(ngllx-1,0:ngllz-1) &
                     + 2.* Tdomain%specel(n_elem1)%Mu(ngllx-1,0:ngllz-1)) / Tdomain%specel(n_elem1)%Density(ngllx-1,0:ngllz-1))
             Tdomain%sFace(n_face)%Zs_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(ngllx-1,0:ngllz-1) * &
                     sqrt(Tdomain%specel(n_elem1)%Mu(ngllx-1,0:ngllz-1)) / Tdomain%specel(n_elem1)%Density(ngllx-1,0:ngllz-1))
          else if (w_face == 2 ) then
             Tdomain%sFace(n_face)%Zp_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0:ngllx-1,ngllz-1) * &
                     sqrt((Tdomain%specel(n_elem1)%Lambda(0:ngllx-1,ngllz-1) &
                     + 2.* Tdomain%specel(n_elem1)%Mu(0:ngllx-1,ngllz-1)) / Tdomain%specel(n_elem1)%Density(0:ngllx-1,ngllz-1))
             Tdomain%sFace(n_face)%Zs_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0:ngllx-1,ngllz-1) * &
                     sqrt(Tdomain%specel(n_elem1)%Mu(0:ngllx-1,ngllz-1)) / Tdomain%specel(n_elem1)%Density(0:ngllx-1,ngllz-1))
          else
             Tdomain%sFace(n_face)%Zp_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0,0:ngllz-1) * &
                     sqrt((Tdomain%specel(n_elem1)%Lambda(0,0:ngllz-1) &
                     + 2.* Tdomain%specel(n_elem1)%Mu(0,0:ngllz-1)) / Tdomain%specel(n_elem1)%Density(0,0:ngllz-1))
             Tdomain%sFace(n_face)%Zs_p(0:ngll-1) = Tdomain%specel(n_elem1)%Density(0,0:ngllz-1) * &
                     sqrt(Tdomain%specel(n_elem1)%Mu(0,0:ngllz-1)) / Tdomain%specel(n_elem1)%Density(0,0:ngllz-1))
          endif
      else !! Case NOT COHERENT
          if (w_face == 0 ) then
             do i=0,ngll-1
                Tdomain%sFace(n_face)%Zp_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1-i,0) * &
                     sqrt((Tdomain%specel(n_elem1)%Lambda(ngllx-1-i,0) &
                     + 2.* Tdomain%specel(n_elem1)%Mu(ngllx-1-i,0)) / Tdomain%specel(n_elem1)%Density(ngllx-1-i,0))
                Tdomain%sFace(n_face)%Zs_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1-i,0) * &
                     sqrt(Tdomain%specel(n_elem1)%Mu(ngllx-1-i,0)) / Tdomain%specel(n_elem1)%Density(ngllx-1-i,0))
             end do
          else if (w_face == 1 ) then
             do i=0,ngll-1
                Tdomain%sFace(n_face)%Zp_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1,ngllz-1-i) * &
                     sqrt((Tdomain%specel(n_elem1)%Lambda(ngllx-1,ngllz-1-i) &
                     + 2.* Tdomain%specel(n_elem1)%Mu(ngllx-1,ngllz-1-i)) / Tdomain%specel(n_elem1)%Density(ngllx-1,ngllz-1-i))
                Tdomain%sFace(n_face)%Zs_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1,ngllz-1-i) * &
                     sqrt(Tdomain%specel(n_elem1)%Mu(ngllx-1,ngllz-1-i)) / Tdomain%specel(n_elem1)%Density(ngllx-1,ngllz-1-i))
             end do
          else if (w_face == 2 ) then
             do i=0,ngll-1
                Tdomain%sFace(n_face)%Zp_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1-i,ngllz-1) * &
                     sqrt((Tdomain%specel(n_elem1)%Lambda(ngllx-1-i,ngllz-1) &
                     + 2.* Tdomain%specel(n_elem1)%Mu(ngllx-1-i,ngllz-1)) / Tdomain%specel(n_elem1)%Density(ngllx-1-i,ngllz-1))
                Tdomain%sFace(n_face)%Zs_p(i) = Tdomain%specel(n_elem1)%Density(ngllx-1-i,ngllz-1) * &
                     sqrt(Tdomain%specel(n_elem1)%Mu(ngllx-1-i,ngllz-1)) / Tdomain%specel(n_elem1)%Density(ngllx-1-i,ngllz-1))
             end do
          else
             do i=0,ngll-1
                Tdomain%sFace(n_face)%Zp_p(i) = Tdomain%specel(n_elem1)%Density(0,ngllz-1-i) * &
                     sqrt((Tdomain%specel(n_elem1)%Lambda(0,ngllz-1-i) &
                     + 2.* Tdomain%specel(n_elem1)%Mu(0,ngllz-1-i)) / Tdomain%specel(n_elem1)%Density(0,ngllz-1-i))
                Tdomain%sFace(n_face)%Zs_p(i) = Tdomain%specel(n_elem1)%Density(0,ngllz-1-i) * &
                     sqrt(Tdomain%specel(n_elem1)%Mu(0,ngllz-1-i)) / Tdomain%specel(n_elem1)%Density(0,ngllz-1-i))
             end do
          endif
      endif

      !! Computation of the coefficients k0, k1, and vector r1 for Godunov fluxs
      do i=0,ngll
         Tdomain%sFace(n_face)%k0(i) = 1./(Tdomain%sFace(n_face)%Zp_m(i) + Tdomain%sFace(n_face)%Zp_p(i))
         Tdomain%sFace(n_face)%k1(i) = 1./(Tdomain%sFace(n_face)%Zs_m(i) + Tdomain%sFace(n_face)%Zs_p(i))
         Tdomain%sFace(n_face)%r1(i,0) = Tdomain%sFace(n_face)%normal(0)**2
         Tdomain%sFace(n_face)%r1(i,1) = Tdomain%sFace(n_face)%normal(1)**2
         Tdomain%sFace(n_face)%r1(i,2) = Tdomain%sFace(n_face)%normal(0) * Tdomain%sFace(n_face)%normal(1)
      end do
    return
end subroutine  compute_coeff_flux
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
