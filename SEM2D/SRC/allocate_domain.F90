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

    implicit none

    type(domain), intent (INOUT):: Tdomain

    integer :: n,ngllx,ngllz,ngll

    do n=0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        ngllz = Tdomain%specel(n)%ngllz
        allocate (Tdomain%specel(n)%MassMat(0:ngllx-1,0:ngllz-1))
        allocate (Tdomain%specel(n)%Density(0:ngllx-1,0:ngllz-1))
        allocate (Tdomain%specel(n)%Lambda(0:ngllx-1,0:ngllz-1))
        allocate (Tdomain%specel(n)%Mu(0:ngllx-1,0:ngllz-1))
        allocate (Tdomain%specel(n)%Veloc (  1:ngllx-2, 1:ngllz-2, 0:1 ) )
        allocate (Tdomain%specel(n)%Forces(  0:ngllx-1, 0:ngllz-1, 0:1 ) )
        allocate (Tdomain%specel(n)%Accel(1:ngllx-2,1:ngllz-2,0:1))
        allocate (Tdomain%specel(n)%V0(  1:ngllx-2, 1:ngllz-2, 0:1 ) )

        Tdomain%specel(n)%Veloc = 0
        Tdomain%specel(n)%Forces = 0
        Tdomain%specel(n)%Accel = 0
        Tdomain%specel(n)%V0 = 0

        if(Tdomain%specel(n)%CPML ) then
            allocate (Tdomain%specel(n)%Stress (0:ngllx-1, 0:ngllz-1, 0:2))
            allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:15))
            allocate (Tdomain%specel(n)%Axi (0:ngllx-1,0:ngllz-1))
            allocate (Tdomain%specel(n)%Aeta(0:ngllx-1,0:ngllz-1))
            allocate (Tdomain%specel(n)%Bxi (0:ngllx-1,0:ngllz-1))
            allocate (Tdomain%specel(n)%Beta(0:ngllx-1,0:ngllz-1))
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
        else ! Case no PML
            allocate (Tdomain%specel(n)%Acoeff(0:ngllx-1,0:ngllz-1,0:9))
            allocate (Tdomain%specel(n)%Displ(1:ngllx-2,1:ngllz-2,0:1))
            Tdomain%specel(n)%Displ = 0
        endif
    enddo

    do n = 0, Tdomain%n_face-1
        ngll = Tdomain%sFace(n)%ngll
        allocate (Tdomain%sFace(n)%MassMat(1:ngll-2))

        allocate (Tdomain%sFace(n)%Veloc (1:ngll-2, 0:1 ) )
        allocate (Tdomain%sFace(n)%Forces(1:ngll-2,0:1 ) )

#ifdef MKA3D
        allocate (Tdomain%sFace(n)%ForcesMka(1:ngll-2,0:1 ) )
        Tdomain%sFace(n)%ForcesMka=0.
#endif


        allocate (Tdomain%sFace(n)%Accel(1:ngll-2, 0:1))
        allocate (Tdomain%sFace(n)%V0( 1:ngll-2, 0:1 ) )
        Tdomain%sFace(n)%MassMat = 0
        Tdomain%sFace(n)%Veloc = 0
        Tdomain%sFace(n)%Accel = 0
        Tdomain%sFace(n)%V0 = 0
        Tdomain%sFace(n)%Forces = 0

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
        else
            allocate (Tdomain%sFace(n)%Displ(1:ngll-2,0:1))
            Tdomain%sFace(n)%Displ = 0
        endif
    enddo

    do n = 0, Tdomain%n_vertex-1
        allocate (Tdomain%sVertex(n)%Veloc (0:1))
        allocate (Tdomain%sVertex(n)%Forces (0:1))
        allocate (Tdomain%sVertex(n)%Accel (0:1))
        allocate (Tdomain%sVertex(n)%V0 (0:1))

#ifdef MKA3D
        allocate (Tdomain%sVertex(n)%ForcesMka(0:1 ) )
        Tdomain%sVertex(n)%ForcesMka=0
#endif


        Tdomain%sVertex(n)%MassMat = 0
        Tdomain%sVertex(n)%Veloc = 0
        Tdomain%sVertex(n)%Accel = 0
        Tdomain%sVertex(n)%V0 = 0

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
        else
            allocate (Tdomain%sVertex(n)%Displ(0:1))
            Tdomain%sVertex(n)%Displ = 0
        endif
    enddo

    return
end subroutine allocate_domain
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
