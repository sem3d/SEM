!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file source_select.f90
!! \brief Choix de l''element qui porte une source ponctuelle
!!
!<

module msource_select
    use constants
    implicit none
contains

    !> Vrai pour un domaine fluide (meme classement que les sources : tout ce qui n''est pas solide)
    pure logical function domain_is_fluid(dom)
        integer, intent(in) :: dom
        domain_is_fluid = .not. (dom == DM_SOLID_CG .or. dom == DM_SOLID_DG .or. dom == DM_SOLID_CG_PML)
    end function domain_is_fluid

    !> Indice (1..n) du premier candidat qui contient le point (coordonnees locales dans
    !! [-1-eps, 1+eps]) et dont le domaine convient au type de source ; 0 si aucun.
    !! Un point sur une frontiere a plusieurs candidats : on passe au suivant quand le
    !! domaine ne convient pas (source solide sur une interface fluide-solide par exemple).
    pure integer function select_source_candidate(n, localcoord, doms, fluid_source, eps)
        integer, intent(in) :: n
        real(fpp), dimension(0:2,n), intent(in) :: localcoord
        integer, dimension(n), intent(in) :: doms
        logical, intent(in) :: fluid_source
        real(fpp), intent(in) :: eps
        integer :: i

        select_source_candidate = 0
        do i = 1,n
            if (any(localcoord(:,i) < -1-eps) .or. any(localcoord(:,i) > 1+eps)) cycle
            if (domain_is_fluid(doms(i)) .neqv. fluid_source) cycle
            select_source_candidate = i
            return
        end do
    end function select_source_candidate

end module msource_select
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
