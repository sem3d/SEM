!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file allocate_domain.f90
!!\brief G�re l'allocation des domaines.
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
    use dom_solidpml
    use dom_fluid
    use dom_fluidpml
    implicit none
contains

subroutine allocate_domain (Tdomain)

    type(domain), intent (INOUT) :: Tdomain
    integer :: n,ngllx,nglly,ngllz
    integer :: mat, randSize, assocMat

    do mat = 0,Tdomain%n_mat-1
        assocMat = Tdomain%sSubdomain(mat)%assocMat
        if(Tdomain%sSubDomain(assocMat)%material_type == "R") then
            if (Tdomain%subD_exist(mat)) then
                call random_seed(size = randSize)

                if(.not.(allocated(Tdomain%sSubdomain(mat)%chosenSeed))) then
                    allocate(Tdomain%sSubdomain(mat)%chosenSeed(randSize))
                end if
            end if
        end if
    end do

    call allocate_dom_solid   (Tdomain, Tdomain%sdom)
    call allocate_dom_fluid   (Tdomain, Tdomain%fdom)
    call allocate_dom_solidpml(Tdomain, Tdomain%spmldom)
    call allocate_dom_fluidpml(Tdomain, Tdomain%fpmldom)

    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        allocate(Tdomain%specel(n)%Density(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%MassMat(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%Lambda (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%Mu     (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
        allocate(Tdomain%specel(n)%Kappa  (0:ngllx-1, 0:nglly-1, 0:ngllz-1))
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
!! vim: set sw=4 ts=8 et tw=80 smartindent :
