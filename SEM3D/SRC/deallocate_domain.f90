!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file deallocate_domain.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine deallocate_domain (Tdomain)

    ! Modified by Gaetano Festa 25/02/2005
    ! Modified by Paul Cupillard 08/12/2005

    use sdomain
    use dom_solid
    use dom_solidpml
    use dom_fluid
    use dom_fluidpml

    implicit none

    type(domain), intent (INOUT):: Tdomain

    integer :: n

    deallocate (Tdomain%GlobCoord)
    deallocate (Tdomain%Coord_Nodes)
    if(allocated(Tdomain%not_PML_List)) deallocate (Tdomain%not_PML_List)
    if(allocated(Tdomain%subD_exist)) deallocate (Tdomain%subD_exist)

    call deallocate_dom_solid   (Tdomain%sdom)
    call deallocate_dom_fluid   (Tdomain%fdom)
    call deallocate_dom_solidpml(Tdomain%spmldom)
    call deallocate_dom_fluidpml(Tdomain%fpmldom)

    do n = 0,Tdomain%n_elem-1
        deallocate (Tdomain%specel(n)%MassMat)
        deallocate (Tdomain%specel(n)%IglobNum)
        deallocate (Tdomain%specel(n)%Control_Nodes)
    enddo

    !purge -fuites memoire
    deallocate (Tdomain%sComm)

    do n = 0, Tdomain%n_mat-1
        if (associated(Tdomain%sSubdomain(n)%GLLcz, Tdomain%sSubdomain(n)%GLLcx) .or. &
            associated(Tdomain%sSubdomain(n)%GLLcz, Tdomain%sSubdomain(n)%GLLcy )) then
            nullify (Tdomain%sSubdomain(n)%GLLcz)
            nullify(Tdomain%sSubdomain(n)%GLLwz)
            nullify(Tdomain%sSubdomain(n)%hprimez)
            nullify(Tdomain%sSubdomain(n)%hTprimez)
        else
            deallocate (Tdomain%sSubdomain(n)%GLLcz)
            deallocate (Tdomain%sSubdomain(n)%GLLwz)
            deallocate (Tdomain%sSubdomain(n)%hprimez)
            deallocate (Tdomain%sSubdomain(n)%hTprimez)
        endif
        if (associated(Tdomain%sSubdomain(n)%GLLcy, Tdomain%sSubdomain(n)%GLLcx) ) then
            nullify (Tdomain%sSubdomain(n)%GLLcy)
            nullify(Tdomain%sSubdomain(n)%GLLwy)
            nullify(Tdomain%sSubdomain(n)%hprimey)
            nullify(Tdomain%sSubdomain(n)%hTprimey)
        else
            deallocate (Tdomain%sSubdomain(n)%GLLcy)
            deallocate (Tdomain%sSubdomain(n)%GLLwy)
            deallocate (Tdomain%sSubdomain(n)%hprimey)
            deallocate (Tdomain%sSubdomain(n)%hTprimey)
        endif
        deallocate (Tdomain%sSubdomain(n)%GLLcx)
        deallocate (Tdomain%sSubdomain(n)%GLLwx)
        deallocate (Tdomain%sSubdomain(n)%hprimex)
        deallocate (Tdomain%sSubdomain(n)%hTprimex)

    enddo

    deallocate (Tdomain%sSubdomain)

    do n = 0, Tdomain%n_source-1
        if (Tdomain%rank==Tdomain%sSource(n)%proc) then
            if (Tdomain%sSource(n)%i_type_source==2) deallocate (Tdomain%sSource(n)%coeff)
            if (Tdomain%sSource(n)%i_time_function==3) deallocate (Tdomain%sSource(n)%timefunc)
        endif
    enddo

    if (Tdomain%logicD%any_source) then
        deallocate (Tdomain%sSource)
    endif

#ifdef COUPLAGE
    !purge - fuites memoire
    do n = 0, Tdomain%n_face-1
        deallocate (Tdomain%sFace(n)%ForcesExt)
        deallocate (Tdomain%sFace(n)%tsurfsem)
    enddo
    do n = 0, Tdomain%n_edge-1
        deallocate (Tdomain%sEdge(n)%ForcesExt)
        deallocate (Tdomain%sEdge(n)%tsurfsem)
    enddo
    do n = 0, Tdomain%n_vertex-1
        deallocate (Tdomain%sVertex(n)%ForcesExt)
    enddo
#endif

    deallocate (Tdomain%specel)
    deallocate (Tdomain%sFace)
    deallocate (Tdomain%sEdge)
    deallocate (Tdomain%sVertex)

    return
end subroutine deallocate_domain

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
