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

    implicit none

    type(domain), intent (INOUT):: Tdomain

    integer :: n, code

    deallocate (Tdomain%GlobCoord)
    deallocate (Tdomain%Coord_Nodes)
    if(allocated(Tdomain%not_PML_List)) deallocate (Tdomain%not_PML_List)
    if(allocated(Tdomain%subD_exist)) deallocate (Tdomain%subD_exist)

    do n = 0,Tdomain%n_elem-1
        deallocate (Tdomain%specel(n)%Density)
        deallocate (Tdomain%specel(n)%MassMat)
        deallocate (Tdomain%specel(n)%IglobNum)
        deallocate (Tdomain%specel(n)%Control_Nodes)
        deallocate (Tdomain%specel(n)%Jacob)
        if (Tdomain%TimeD%velocity_scheme) then
            if (Tdomain%specel(n)%PML) then
                !  modif mariotti fevrier 2007 cea
                deallocate (Tdomain%specel(n)%Lambda)
                deallocate (Tdomain%specel(n)%Kappa)
                deallocate (Tdomain%specel(n)%Mu)

                if (Tdomain%specel(n)%solid) then
                    deallocate (Tdomain%specel(n)%sl%Acoeff)
                    deallocate (Tdomain%specel(n)%slpml%Diagonal_Stress)
                    deallocate (Tdomain%specel(n)%slpml%Diagonal_Stress1)
                    deallocate (Tdomain%specel(n)%slpml%Diagonal_Stress2)
                    deallocate (Tdomain%specel(n)%slpml%Diagonal_Stress3)
                    deallocate (Tdomain%specel(n)%slpml%Residual_Stress)
                    deallocate (Tdomain%specel(n)%slpml%Residual_Stress1)
                    deallocate (Tdomain%specel(n)%slpml%Residual_Stress2)
                    deallocate (Tdomain%specel(n)%slpml%Forces1)
                    deallocate (Tdomain%specel(n)%slpml%Forces2)
                    deallocate (Tdomain%specel(n)%slpml%Forces3)
                endif
                deallocate (Tdomain%specel(n)%xpml%DumpSx)
                deallocate (Tdomain%specel(n)%xpml%DumpSy)
                deallocate (Tdomain%specel(n)%xpml%DumpSz)
            else
                if (Tdomain%aniso) then
                    if (Tdomain%specel(n)%solid) deallocate (Tdomain%specel(n)%sl%Cij)
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
                if (Tdomain%specel(n)%solid .and. Tdomain%n_sls>0) then
                    if (Tdomain%aniso) then
                        deallocate (Tdomain%specel(n)%sl%Q)
                    else
                        !             deallocate (Tdomain%specel(n)%Kappa)
                        deallocate (Tdomain%specel(n)%sl%Qs)
                        deallocate (Tdomain%specel(n)%sl%Qp)
                        deallocate (Tdomain%specel(n)%sl%onemPbeta)
                        deallocate (Tdomain%specel(n)%sl%factor_common_P)
                        deallocate (Tdomain%specel(n)%sl%alphaval_P)
                        deallocate (Tdomain%specel(n)%sl%betaval_P)
                        deallocate (Tdomain%specel(n)%sl%gammaval_P)
                        deallocate (Tdomain%specel(n)%sl%epsilonvol_)
                        deallocate (Tdomain%specel(n)%sl%R_vol_)
                    endif
                    deallocate (Tdomain%specel(n)%sl%onemSbeta)
                    deallocate (Tdomain%specel(n)%sl%factor_common_3)
                    deallocate (Tdomain%specel(n)%sl%alphaval_3)
                    deallocate (Tdomain%specel(n)%sl%betaval_3)
                    deallocate (Tdomain%specel(n)%sl%gammaval_3)
                    deallocate (Tdomain%specel(n)%sl%epsilondev_xx_)
                    deallocate (Tdomain%specel(n)%sl%epsilondev_yy_)
                    deallocate (Tdomain%specel(n)%sl%epsilondev_xy_)
                    deallocate (Tdomain%specel(n)%sl%epsilondev_xz_)
                    deallocate (Tdomain%specel(n)%sl%epsilondev_yz_)
                    deallocate (Tdomain%specel(n)%sl%R_xx_)
                    deallocate (Tdomain%specel(n)%sl%R_yy_)
                    deallocate (Tdomain%specel(n)%sl%R_xy_)
                    deallocate (Tdomain%specel(n)%sl%R_xz_)
                    deallocate (Tdomain%specel(n)%sl%R_yz_)
                endif

            endif
        endif
        deallocate (Tdomain%specel(n)%InvGrad)     !purge fuites memoire Gsa
        if(Tdomain%specel(n)%solid) then
            deallocate (Tdomain%specel(n)%sl)
            if (Tdomain%specel(n)%PML) then
                deallocate(Tdomain%specel(n)%slpml)
            endif
        end if
    enddo

    do n = 0, Tdomain%n_face-1
        if (Tdomain%sFace(n)%PML) then
        else
            !  modif mariotti fevrier 2007 cea capteur displ
            !        deallocate (Tdomain%sFace(n)%Displ)
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

        if (Tdomain%any_Random) then
            call MPI_COMM_FREE (Tdomain%subDComm(n),code)
        end if
    enddo

    if (Tdomain%any_Random) then
        if(allocated(Tdomain%subDComm)) deallocate (Tdomain%subDComm)
    end if
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
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
