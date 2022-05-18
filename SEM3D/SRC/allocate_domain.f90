!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file allocate_domain.f90
!!\brief Gére l'allocation des domaines.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief Allocation des attributs de la structure domain
!!
!! \param type(domain), intent (INOUT) Tdomain
!<
module sdomain_alloc
    use sdomain
    use dom_solid
    use dom_solid_dg
    use dom_solidpml
    use dom_fluid
    use dom_fluidpml
    implicit none
contains

subroutine allocate_domain (Tdomain)

    type(domain), intent (INOUT) :: Tdomain
    integer :: n,ngll
    !integer :: mat, randSize, assocMat

    if(.false.) then
        write(*,*) "Tdomain%any_sdom = ", Tdomain%any_sdom
        write(*,*) "Tdomain%any_fdom = ", Tdomain%any_fdom
        write(*,*) "Tdomain%any_spml = ", Tdomain%any_spml
        write(*,*) "Tdomain%any_fpml = ", Tdomain%any_fpml
    end if
    if(Tdomain%any_sdom)   call allocate_dom_solid   (Tdomain, Tdomain%sdom)
    if(Tdomain%any_sdomdg) call allocate_dom_solid_dg(Tdomain, Tdomain%sdomdg)
    if(Tdomain%any_fdom)   call allocate_dom_fluid   (Tdomain, Tdomain%fdom)
    if(Tdomain%any_spml)   call allocate_dom_solidpml(Tdomain, Tdomain%spmldom)
    if(Tdomain%any_fpml)   call allocate_dom_fluidpml(Tdomain, Tdomain%fpmldom)

    do n = 0,Tdomain%n_elem-1
        ngll = domain_ngll(Tdomain, Tdomain%specel(n)%domain)
        allocate(Tdomain%specel(n)%MassMat(0:ngll-1, 0:ngll-1, 0:ngll-1))
        allocate(Tdomain%specel(n)%En_S_avg(0:ngll-2, 0:ngll-2, 0:ngll-2)) 
        allocate(Tdomain%specel(n)%En_P_avg(0:ngll-2, 0:ngll-2, 0:ngll-2))

    enddo
end subroutine allocate_domain

end module sdomain_alloc

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
