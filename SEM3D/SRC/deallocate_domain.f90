!>
!! \file deallocate_domain.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine deallocate_domain (Tdomain, rg)

    ! Modified by Gaetano Festa 25/02/2005
    ! Modified by Paul Cupillard 08/12/2005


    use sdomain

    implicit none

    type(domain), intent (INOUT):: Tdomain
    integer, intent (IN) :: rg


    integer :: n


    deallocate (Tdomain%GlobCoord)
    deallocate (Tdomain%Coord_Nodes)

    do n = 0,Tdomain%n_elem-1
        deallocate (Tdomain%specel(n)%Density)
        deallocate (Tdomain%specel(n)%MassMat)
        deallocate (Tdomain%specel(n)%IglobNum)
        deallocate (Tdomain%specel(n)%Control_Nodes)
        deallocate (Tdomain%specel(n)%Jacob)
        if (Tdomain%TimeD%velocity_scheme) then
            deallocate (Tdomain%specel(n)%Veloc )
            !  modif mariotti fevrier 2007 cea capteur displ
            deallocate (Tdomain%specel(n)%Displ)

            deallocate (Tdomain%specel(n)%Accel)
            deallocate (Tdomain%specel(n)%V0 )
            deallocate (Tdomain%specel(n)%Forces)
            if (Tdomain%specel(n)%PML) then
                !  modif mariotti fevrier 2007 cea
                deallocate (Tdomain%specel(n)%Lambda)
                deallocate (Tdomain%specel(n)%Kappa)
                deallocate (Tdomain%specel(n)%Mu)

                deallocate (Tdomain%specel(n)%Acoeff)
                deallocate (Tdomain%specel(n)%spml%Diagonal_Stress)
                deallocate (Tdomain%specel(n)%spml%Diagonal_Stress1)
                deallocate (Tdomain%specel(n)%spml%Diagonal_Stress2)
                deallocate (Tdomain%specel(n)%spml%Diagonal_Stress3)
                deallocate (Tdomain%specel(n)%spml%Residual_Stress)
                deallocate (Tdomain%specel(n)%spml%Residual_Stress1)
                deallocate (Tdomain%specel(n)%spml%Residual_Stress2)
                deallocate (Tdomain%specel(n)%spml%Veloc1)
                deallocate (Tdomain%specel(n)%spml%Veloc2)
                deallocate (Tdomain%specel(n)%spml%Veloc3)
                deallocate (Tdomain%specel(n)%spml%Forces1)
                deallocate (Tdomain%specel(n)%spml%Forces2)
                deallocate (Tdomain%specel(n)%spml%Forces3)
                deallocate (Tdomain%specel(n)%spml%DumpSx)
                deallocate (Tdomain%specel(n)%spml%DumpSy)
                deallocate (Tdomain%specel(n)%spml%DumpSz)
                deallocate (Tdomain%specel(n)%spml%DumpVx)
                deallocate (Tdomain%specel(n)%spml%DumpVy)
                deallocate (Tdomain%specel(n)%spml%DumpVz)
            else
                if (Tdomain%aniso) then
                    deallocate (Tdomain%specel(n)%Cij)
                    if (Tdomain%n_sls>0) then
                        deallocate (Tdomain%specel(n)%Lambda)
                        deallocate (Tdomain%specel(n)%Kappa)
                        deallocate (Tdomain%specel(n)%Mu)
                    endif
                else
                    deallocate (Tdomain%specel(n)%Lambda)
                    deallocate (Tdomain%specel(n)%Kappa)
                    deallocate (Tdomain%specel(n)%Mu)
                endif
                if (Tdomain%n_sls>0) then
                    if (Tdomain%aniso) then
                        deallocate (Tdomain%specel(n)%Q)
                    else
                        !             deallocate (Tdomain%specel(n)%Kappa)
                        deallocate (Tdomain%specel(n)%Qs)
                        deallocate (Tdomain%specel(n)%Qp)
                        deallocate (Tdomain%specel(n)%onemPbeta)
                        deallocate (Tdomain%specel(n)%factor_common_P)
                        deallocate (Tdomain%specel(n)%alphaval_P)
                        deallocate (Tdomain%specel(n)%betaval_P)
                        deallocate (Tdomain%specel(n)%gammaval_P)
                        deallocate (Tdomain%specel(n)%epsilonvol_)
                        deallocate (Tdomain%specel(n)%R_vol_)
                    endif
                    deallocate (Tdomain%specel(n)%onemSbeta)
                    deallocate (Tdomain%specel(n)%factor_common_3)
                    deallocate (Tdomain%specel(n)%alphaval_3)
                    deallocate (Tdomain%specel(n)%betaval_3)
                    deallocate (Tdomain%specel(n)%gammaval_3)
                    deallocate (Tdomain%specel(n)%epsilondev_xx_)
                    deallocate (Tdomain%specel(n)%epsilondev_yy_)
                    deallocate (Tdomain%specel(n)%epsilondev_xy_)
                    deallocate (Tdomain%specel(n)%epsilondev_xz_)
                    deallocate (Tdomain%specel(n)%epsilondev_yz_)
                    deallocate (Tdomain%specel(n)%R_xx_)
                    deallocate (Tdomain%specel(n)%R_yy_)
                    deallocate (Tdomain%specel(n)%R_xy_)
                    deallocate (Tdomain%specel(n)%R_xz_)
                    deallocate (Tdomain%specel(n)%R_yz_)
                endif

            endif
        endif
        if (.not. Tdomain%specel(n)%PML) then
            deallocate (Tdomain%specel(n)%InvGrad)     !purge fuites memoire Gsa
        endif
    enddo

    do n = 0, Tdomain%n_face-1
        deallocate (Tdomain%sFace(n)%MassMat)
        deallocate (Tdomain%sFace(n)%Veloc)

        !  modif mariotti fevrier 2007 cea capteur displ
        deallocate (Tdomain%sFace(n)%Displ)

        deallocate (Tdomain%sFace(n)%Forces)
        deallocate (Tdomain%sFace(n)%Accel)
        deallocate (Tdomain%sFace(n)%V0)
        if (Tdomain%sFace(n)%PML) then
            deallocate (Tdomain%sFace(n)%spml%Forces1)
            deallocate (Tdomain%sFace(n)%spml%Forces2)
            deallocate (Tdomain%sFace(n)%spml%Forces3)
            deallocate (Tdomain%sFace(n)%spml%Veloc1)
            deallocate (Tdomain%sFace(n)%spml%Veloc2)
            deallocate (Tdomain%sFace(n)%spml%Veloc3)
            deallocate (Tdomain%sFace(n)%spml%DumpVx)
            deallocate (Tdomain%sFace(n)%spml%DumpVy)
            deallocate (Tdomain%sFace(n)%spml%DumpVz)
        else
            !  modif mariotti fevrier 2007 cea capteur displ
            !        deallocate (Tdomain%sFace(n)%Displ)
        endif
    enddo

    do n = 0,Tdomain%n_edge-1
        deallocate (Tdomain%sEdge(n)%MassMat)
        deallocate (Tdomain%sEdge(n)%Veloc)
        !  modif mariotti fevrier 2007 cea capteur displ
        deallocate (Tdomain%sEdge(n)%Displ)

        deallocate (Tdomain%sEdge(n)%Forces)
        deallocate (Tdomain%sEdge(n)%Accel)
        deallocate (Tdomain%sEdge(n)%V0)
        if (Tdomain%sEdge(n)%PML) then
            deallocate (Tdomain%sEdge(n)%spml%Forces1)
            deallocate (Tdomain%sEdge(n)%spml%Forces2)
            deallocate (Tdomain%sEdge(n)%spml%Forces3)
            deallocate (Tdomain%sEdge(n)%spml%Veloc1)
            deallocate (Tdomain%sEdge(n)%spml%Veloc2)
            deallocate (Tdomain%sEdge(n)%spml%Veloc3)
            deallocate (Tdomain%sEdge(n)%spml%DumpVx)
            deallocate (Tdomain%sEdge(n)%spml%DumpVy)
            deallocate (Tdomain%sEdge(n)%spml%DumpVz)
        else
            !  modif mariotti fevrier 2007 cea capteur displ
            !        deallocate (Tdomain%sEdge(n)%Displ)
        endif
    enddo

    do n = 0,Tdomain%n_proc-1
        if (Tdomain%sComm(n)%ngll>0) then
            deallocate (Tdomain%sComm(n)%GiveForces)
            deallocate (Tdomain%sComm(n)%TakeForces)
        endif
        if (Tdomain%sComm(n)%ngllPML>0) then
            deallocate (Tdomain%sComm(n)%GiveForcesPML)
            deallocate (Tdomain%sComm(n)%TakeForcesPML)
        endif
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
        if (rg==Tdomain%sSource(n)%proc) then
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
        deallocate (Tdomain%sFace(n)%ForcesMka)
        deallocate (Tdomain%sFace(n)%tsurfsem)
    enddo
    do n = 0, Tdomain%n_edge-1
        deallocate (Tdomain%sEdge(n)%ForcesMka)
        deallocate (Tdomain%sEdge(n)%tsurfsem)
    enddo
    do n = 0, Tdomain%n_vertex-1
        deallocate (Tdomain%sVertex(n)%ForcesMka)
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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
