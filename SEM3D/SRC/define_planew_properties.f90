!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!! \file define_planew_properties.f90
!! \brief
!! \author
!! \version 1.0
!! \date
!!
!<

subroutine define_planew_properties (Tdomain)

    use sdomain
    use semdatafiles

    implicit none

    type(Domain), intent (INOUT) :: Tdomain
    integer                      :: ngll,ngll1,ngll2, nf,ne,nv,mat_index,rg

    rg = Tdomain%rank

    do nf = 0, Tdomain%sPlaneW%n_faces-1
        ngll1 = Tdomain%sPlaneW%pFace(nf)%ngll1
        ngll2 = Tdomain%sPlaneW%pFace(nf)%ngll2
        Tdomain%sPlaneW%pFace(nf)%mat_index = Tdomain%sPlaneW%pParam%mat_index
        allocate (Tdomain%sPlaneW%pFace(nf)%Btn(0:ngll1-1,0:ngll2-1,0:2))
        allocate (Tdomain%sPlaneW%pFace(nf)%MassMat_Up(1:ngll1-2,1:ngll2-2))
        allocate (Tdomain%sPlaneW%pFace(nf)%MassMat_Down(1:ngll1-2,1:ngll2-2))
        allocate (Tdomain%sPlaneW%pFace(nf)%Forces_Up(1:ngll1-2,1:ngll2-2,0:2))
        allocate (Tdomain%sPlaneW%pFace(nf)%Forces_Down(1:ngll1-2,1:ngll2-2,0:2))
    enddo

    do ne = 0, Tdomain%sPlaneW%n_edges-1
        ngll = Tdomain%sPlaneW%pEdge(ne)%ngll
        Tdomain%sPlaneW%pEdge(ne)%mat_index = Tdomain%sPlaneW%pParam%mat_index
        allocate (Tdomain%sPlaneW%pEdge(ne)%Btn(1:ngll-2,0:2))
        allocate (Tdomain%sPlaneW%pEdge(ne)%MassMat_Up(1:ngll-2))
        allocate (Tdomain%sPlaneW%pEdge(ne)%MassMat_Down(1:ngll-2))
        allocate (Tdomain%sPlaneW%pEdge(ne)%Forces_Up(1:ngll-2,0:2))
        allocate (Tdomain%sPlaneW%pEdge(ne)%Forces_Down(1:ngll-2,0:2))
        Tdomain%sPlaneW%pEdge(ne)%Btn = 0
    enddo

    do nv = 0, Tdomain%sPlaneW%n_vertices-1
        Tdomain%sPlaneW%pVertex(nv)%mat_index = Tdomain%sPlaneW%pParam%mat_index
        Tdomain%sPlaneW%pVertex(nv)%Btn = 0
    enddo

    mat_index = Tdomain%sPlaneW%pParam%mat_index
    Tdomain%sPlaneW%pParam%Lambda = Tdomain%sSubDomain(mat_index)%DLambda
    Tdomain%sPlaneW%pParam%Kappa = Tdomain%sSubDomain(mat_index)%DKappa
    Tdomain%sPlaneW%pParam%Mu = Tdomain%sSubDomain(mat_index)%DMu

    if (Tdomain%sPlaneW%pParam%wtype == 'S') then
        Tdomain%sPlaneW%pParam%speed = Tdomain%sSubDomain(mat_index)%Sspeed
    else
        Tdomain%sPlaneW%pParam%speed = Tdomain%sSubDomain(mat_index)%Pspeed
    endif

    return


end subroutine define_planew_properties


function getwtype(wtype)

Implicit none
Integer, intent(in) :: wtype
character           :: getwtype

   select case (wtype)
        case (1)
           getwtype = "P"
        case (2)
           getwtype = "S"
        case (3)
           getwtype = "R"
        case (4)
           getwtype = "A"
   end select

end function getwtype

!subroutine read_planeW_input(Tdomain)


!use sdomain
!use semdatafiles

!implicit none
! 
!type(Domain), intent (INOUT) :: Tdomain 
!integer                      :: i, rg
!character(Len=MAX_FILE_SIZE) :: fnamef
!character(len=800)           :: parametric_var
!character                    :: getwtype
!
!Tdomain%sPlaneW%pParam%mat_index = Tdomain%config%plane_wave_mat
!Tdomain%sPlaneW%pParam%wtype= getwtype(Tdomain%config%plane_wave_type)
!Tdomain%sPlaneW%pParam%lx   = Tdomain%config%plane_wave_L(1)
!Tdomain%sPlaneW%pParam%ly   = Tdomain%config%plane_wave_L(2)
!Tdomain%sPlaneW%pParam%lz   = Tdomain%config%plane_wave_L(3)
!Tdomain%sPlaneW%pParam%xs   = Tdomain%config%plane_wave_C(1)
!Tdomain%sPlaneW%pParam%ys   = Tdomain%config%plane_wave_C(2)
!Tdomain%sPlaneW%pParam%zs   = Tdomain%config%plane_wave_C(3)
!Tdomain%sPlaneW%pParam%f0   = Tdomain%config%plane_wave_f0
!
!if (Tdomain%config%plane_wave_whatbc.eq.1) then
!   Tdomain%logicD%Neumann = .true.
!   Tdomain%Neumann%Neu_Param%lx = Tdomain%sPlaneW%pParam%lx
!   Tdomain%Neumann%Neu_Param%ly = Tdomain%sPlaneW%pParam%ly
!   Tdomain%Neumann%Neu_Param%lz = Tdomain%sPlaneW%pParam%lz
!   Tdomain%Neumann%Neu_Param%xs = Tdomain%sPlaneW%pParam%xs
!   Tdomain%Neumann%Neu_Param%ys = Tdomain%sPlaneW%pParam%ys
!   Tdomain%Neumann%Neu_Param%zs = Tdomain%sPlaneW%pParam%zs 
!   Tdomain%Neumann%Neu_Param%f0 = Tdomain%sPlaneW%pParam%f0
!endif
!
!if (Tdomain%sPlaneW%pParam%wtype=="A") then
!   Tdomain%sPlaneW%pSource%PlaneW_dim   = Tdomain%config%plane_wave_dim
!   Tdomain%sPlaneW%pSource%PlaneW_varia(1:LEN_TRIM(fromcstr(Tdomain%config%plane_wave_varia))) = &
!                                        TRIM(fromcstr(Tdomain%config%plane_wave_varia))
!   Tdomain%sPlaneW%pSource%PlaneW_source(1:LEN_TRIM(fromcstr(Tdomain%config%plane_wave_source)))= &
!                                           TRIM(fromcstr(Tdomain%config%plane_wave_source))
!   Tdomain%sPlaneW%pSource%PlaneW_funcx(1:LEN_TRIM(fromcstr(Tdomain%config%plane_wave_funcx))) = &
!                                            TRIM(fromcstr(Tdomain%config%plane_wave_funcx))
!   if (Tdomain%sPlaneW%pSource%PlaneW_dim.gt.1) then
!       Tdomain%sPlaneW%pSource%PlaneW_funcy(1:LEN_TRIM(fromcstr(Tdomain%config%plane_wave_funcy))) = &
!                                                 TRIM(fromcstr(Tdomain%config%plane_wave_funcy))
!   endif
!   if (Tdomain%sPlaneW%pSource%PlaneW_dim.gt.2) then
!      Tdomain%sPlaneW%pSource%PlaneW_funcz(1:LEN_TRIM(fromcstr(Tdomain%config%plane_wave_funcz))) = &
!                                                  TRIM(fromcstr(Tdomain%config%plane_wave_funcz))
!   endif
!   if (Tdomain%sPlaneW%pSource%PlaneW_source.eq."M") then
!       Tdomain%sPlaneW%pSource%PlaneW_funcxy(1:LEN_TRIM(fromcstr(Tdomain%config%plane_wave_funcxy)))= &
!                                          TRIM(fromcstr(Tdomain%config%plane_wave_funcxy))
!        if (Tdomain%sPlaneW%pSource%PlaneW_dim.eq.3) then
!            Tdomain%sPlaneW%pSource%PlaneW_funcxz(1:LEN_TRIM(fromcstr(Tdomain%config%plane_wave_funcxz)))= &
!                                                TRIM(fromcstr(Tdomain%config%plane_wave_funcxz))
!            Tdomain%sPlaneW%pSource%PlaneW_funcyz(1:LEN_TRIM(fromcstr(Tdomain%config%plane_wave_funcyz)))= &
!                                                 TRIM(fromcstr(Tdomain%config%plane_wave_funcyz))
!         endif
!    endif
!    !!
!    !!
!    Tdomain%sPlaneW%pSource%PlaneW_paramvar=Tdomain%config%plane_wave_paramvar
!    if (Tdomain%sPlaneW%pSource%PlaneW_paramvar==1) then
!        Tdomain%sPlaneW%pSource%PlaneW_nparamvar = Tdomain%config%plane_wave_nparamvar
!       do i=1,Tdomain%config%plane_wave_nparamvar
!            Tdomain%sPlaneW%pSource%PlaneW_paravalue(i) = Tdomain%config%plane_wave_Paravalue(i)
!        enddo
!        parametric_var(1:LEN_TRIM(fromcstr(Tdomain%config%plane_wave_paramname)))=TRIM(fromcstr(Tdomain%config%plane_wave_paramname))
!        call Split(parametric_var, Tdomain%config%plane_wave_nparamvar, &
!                   Tdomain%sPlaneW%pSource%PlaneW_paramname)
!     endif
!     !!
!     !! Ecriture par un seul proc
!     rg = Tdomain%rank
!     if (rg.eq.0) then
!     write(*,*)
!     write(*,*) "Plane wave Source analytical defined "
!     write(*,*) "-------------------------------------"
!     write(*,*) "Va = "//TRIM(Tdomain%sPlaneW%pSource%PlaneW_varia)
!     write(*,*) "F1 = "//TRIM(Tdomain%sPlaneW%pSource%PlaneW_funcx)
!     if (Tdomain%sPlaneW%pSource%PlaneW_dim.gt.1) then
!         write(*,*) "F2 = "//TRIM(Tdomain%sPlaneW%pSource%PlaneW_funcy)
!     endif
!     if (Tdomain%sPlaneW%pSource%PlaneW_dim.gt.2) then
!         write(*,*) "F3 = "//TRIM(Tdomain%sPlaneW%pSource%PlaneW_funcz)
!     endif
!     if (Tdomain%sPlaneW%pSource%PlaneW_source.eq."M") then
!         write(*,*) "F4 = "//TRIM(Tdomain%sPlaneW%pSource%PlaneW_funcxy)
!         if (Tdomain%sPlaneW%pSource%PlaneW_dim.gt.2) then
!             write(*,*) "F5 = "//TRIM(Tdomain%sPlaneW%pSource%PlaneW_funcxz)
!             write(*,*) "F6 = "//TRIM(Tdomain%sPlaneW%pSource%PlaneW_funcyz)
!         endif
!         write(*,*)
!     endif
!     if (Tdomain%sPlaneW%pSource%PlaneW_paramvar==1) then
!         write(*,*) "Param = "//TRIM(parametric_var)
!     endif
!     endif
!endif
!
!end subroutine read_planeW_input

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
