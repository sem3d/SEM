module build_prop_files

    use sdomain
    use mpi
    use constants
    use scomm
    use scommutils

    implicit none

    private :: gcoord, gindex
#include "index.h"

contains

    ! etant donne une grille de n points regulierement espaces entre xmin et xmax
    ! revoie grille(i)
    function gcoord(i, n, xmin, xmax) result(x)
        integer, intent(in) :: i,n
        real(fpp), intent(in) :: xmin, xmax
        real(fpp) :: x

        x = xmin + i*(xmax-xmin)/(n-1)
    end function gcoord

    ! soit une grille de n points reguliers entre xmin et xmax, renvoie
    ! i tel que x(i)<=x<x(i+1)   x(0) = xmin x(n-1)=xmax
    function gindex(x, n, xmin, xmax) result(i)
        integer, intent(in) :: n
        real(fpp), intent(in) :: x, xmin, xmax
        integer :: i

        i = floor( ((n-1)*(x-xmin))/(xmax-xmin) )
    end function gindex

    subroutine init_prop_file(mat)
        use hdf5
        use sem_hdf5
        integer :: i
        type(subdomain), intent(inout) :: mat

        select case(mat%deftype)
        case(MATDEF_VP_VS_RHO)
            mat%prop_field(1)%propName = "Vp"
            mat%prop_field(2)%propName = "Vs"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_E_NU_RHO)
            mat%prop_field(1)%propName = "E"
            mat%prop_field(2)%propName = "Nu"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_LAMBDA_MU_RHO)
            mat%prop_field(1)%propName = "Lambda"
            mat%prop_field(2)%propName = "Mu"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_KAPPA_MU_RHO)
            mat%prop_field(1)%propName = "Kappa"
            mat%prop_field(2)%propName = "Mu"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_HOOKE_RHO)
            stop "Not supported yet"
        case(MATDEF_NLKP_VS_RHO)
            mat%prop_field(1)%propName = "NLkin"
            mat%prop_field(2)%propName = "Vs"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_NU_VS_RHO)
            mat%prop_field(1)%propName = "Nu"
            mat%prop_field(2)%propName = "Vs"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_VTI_ANISO)
            mat%prop_field(1)%propName = "Vpv"
            mat%prop_field(2)%propName = "Vph"
            mat%prop_field(3)%propName = "Vsv"
            mat%prop_field(4)%propName = "Vsh"
            mat%prop_field(5)%propName = "Eta"
            mat%prop_field(6)%propName = "Rho"
            mat%prop_field(7)%propName = "Qkappa"
            mat%prop_field(8)%propName = "Qmu"
        end select

        do i = 1,mat%n_prop
            call init_prop_file_field(mat, mat%prop_field(i))
        end do

    end subroutine init_prop_file

    ! Check for existence of a group in HDF5 file
    function prop_check_var(pid, pname, gid)
        use hdf5
        integer(HID_T), intent(in) :: pid
        integer(HID_T), intent(out) :: gid
        character(len=100) :: pname
        !
        logical :: ok, prop_check_var
        integer :: hdferr
        !
        prop_check_var = .false.
        call H5Lexists_f(pid, trim(adjustl(pname)), ok, hdferr)
        if (ok) then
            call H5Gopen_f(pid, trim(adjustl(pname)), gid, hdferr)
            if (hdferr == 0) prop_check_var = .true.
        endif
    end function prop_check_var

    subroutine init_prop_file_field(mat, pf)
        use hdf5
        use sem_hdf5
        type(subdomain), intent(inout) :: mat
        type(PropertyField), intent(inout) :: pf
        !
        integer :: k, hdferr
        integer(HID_T) :: file_id, grp_id
        integer(HSIZE_T), dimension(:), allocatable :: dims
        logical :: subgrp
        real(fpp) :: minmax_swap
        real(fpp), dimension(0:2) :: min_bound_loc,max_bound_loc,xxr

        if (.not. mat%present) return
        call init_hdf5()
        call h5fopen_f(trim(pf%propFilePath), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
        if(hdferr /= 0) then
            write(*,*) "Could not open file:", trim(pf%propFilePath)
            stop 1
        end if
        subgrp = prop_check_var(file_id, pf%propName, grp_id)
        if (.not. subgrp) then
            grp_id = file_id
        end if
        call read_attr_real_vec(grp_id, "xMinGlob", pf%MinBound)
        call read_attr_real_vec(grp_id, "xMaxGlob", pf%MaxBound)
        call read_dims(grp_id, "samples", dims)

        pf%NN = int(dims)
        ! On va calculer les indices i0,i1 j0,j1,k0,k1 tels que
        ! i0 plus grand entier tel que x(i0)<MinBound_loc(0), 
        ! i1 plus petit entier tel que x(i1)>MaxBound_loc(0), etc... 

        !!! useless since we force full model read, to delete
        min_bound_loc=mat%MinBound_Loc
        max_bound_loc=mat%MaxBound_Loc
        if (mat%is_sph) then
            xxr=matmul(mat%sph_args%R_from_pole_chk,min_bound_loc)
            call cart2sph(xxr,min_bound_loc,.true.)
            xxr=matmul(mat%sph_args%R_from_pole_chk,max_bound_loc)
            call cart2sph(xxr,max_bound_loc,.true.)
            !!! reverse latitude correction
            minmax_swap=min_bound_loc(1)
            min_bound_loc(1)=max_bound_loc(1)
            max_bound_loc(1)=minmax_swap
        end if

        do k = 0,2
            pf%imin(k) = gindex(min_bound_loc(k), pf%NN(k), pf%MinBound(k), pf%MaxBound(k))
            pf%imax(k) = gindex(max_bound_loc(k), pf%NN(k), pf%MinBound(k), pf%MaxBound(k))+1
            pf%step(k) = (pf%MaxBound(k)-pf%MinBound(k))/(pf%NN(k)-1)
            if ((pf%imax(k)-pf%imin(k))<1) pf%imin(k) = pf%imax(k)-1
            if (pf%imin(k)<0) then
                pf%imin(k) = 0
                if (pf%imax(k)<1) pf%imax(k) = 1
            endif
            if (pf%imax(k)>=pf%NN(k)) then
                pf%imax(k) = pf%NN(k)-1
                if (pf%imin(k)>(pf%imax(k)-1)) pf%imin(k) = pf%imax(k)-1
            endif
            if ((pf%imax(k)-pf%imin(k))<1) pf%imax(k) = pf%imin(k)
            !!!if (pf%imin(k)<0) pf%imin(k) = 0
            !!!if (pf%imax(k)<pf%imin(k)) pf%imax(k) = pf%imin(k)
            !!!if (pf%imax(k)>=pf%NN(k)) pf%imax(k) = pf%NN(k)-1
            !!!if (pf%imin(k)>pf%imax(k)) pf%imin(k) = pf%imax(k)
            !!! spherical material :: force full model read for all procs
            if (mat%is_sph) then
                pf%imin(k)=0
                pf%imax(k)=pf%NN(k)-1
            end if
        end do

        call read_subset_3d_real(grp_id, "samples", pf%imin, pf%imax, pf%var)

        if (subgrp) call H5Gclose_f(grp_id, hdferr)
        call H5Fclose_f(file_id, hdferr)
    end subroutine init_prop_file_field

    subroutine cleanup_prop_file(mat)
        integer :: i
        type(subdomain) :: mat

        if (.not. mat%present) return
        do i = 1,mat%n_prop
            deallocate(mat%prop_field(i)%var)
        end do
        deallocate(mat%prop_field)

    end subroutine cleanup_prop_file


    subroutine interpolate_elem_field(Tdomain, specel, mat, pf, field)
        type(domain), intent(inout) :: Tdomain
        type(element), intent(inout) :: specel
        type(subdomain), intent(in) :: mat
        type(PropertyField), intent(in) :: pf
        !
        real(fpp), dimension(0:mat%ngll-1,0:mat%ngll-1,0:mat%ngll-1) :: field
        real(fpp), dimension(0:2) :: xx,xxr ! node coord
        real(fpp), dimension(0:2) :: aa   ! node interpolation coeffs
        real(fpp) :: pos, val
        integer :: i,j,k     ! gll index
        integer :: idef      ! global coord index
        integer :: n
        integer, dimension(0:2) :: ii  ! index of x,y,z inside material domain
        logical :: pml


        pml = is_pml(mat)
        do k = 0,mat%ngll-1
            do j = 0,mat%ngll-1
                do i = 0,mat%ngll-1
                    idef = specel%Iglobnum(i,j,k)
                    do n=0,2
                        xx(n) = Tdomain%GlobCoord(n,idef)
                        ! Traitement PML : on echantillonne au bord de la PML uniquement
                        if (pml) then
                            if (mat%pml_width(n)>0) then
                                if (xx(n)>mat%pml_pos(n)) xx(n)=mat%pml_pos(n)
                            end if
                            if (mat%pml_width(n)<0) then
                                if (xx(n)<mat%pml_pos(n)) xx(n)=mat%pml_pos(n)
                            end if
                        end if
                    end do
                    ! spherical material
                    if (mat%is_sph) then
                        xxr=matmul(mat%sph_args%R_from_pole_chk,xx)
                        call cart2sph(xxr,xx,.true.)
                    end if
                    do n=0,2
                        pos = (xx(n)-pf%MinBound(n))/pf%step(n)
                        ii(n) = floor(pos)
                        if (ii(n)  < pf%imin(n)) ii(n) = pf%imin(n)
                        if (ii(n) >= pf%imax(n)) ii(n) = pf%imax(n)-1
                        aa(n) = pos-ii(n)
                        if (aa(n)<0.) aa(n)=0.
                        if (aa(n)>1.) aa(n)=1.
                    end do
                    ! trilinear interpolation
                    val =       (1.-aa(0))*(1.-aa(1))*(1.-aa(2))*pf%var(ii(0)  ,ii(1)  ,ii(2)  )
                    val = val + (   aa(0))*(1.-aa(1))*(1.-aa(2))*pf%var(ii(0)+1,ii(1)  ,ii(2)  )
                    val = val + (1.-aa(0))*(   aa(1))*(1.-aa(2))*pf%var(ii(0)  ,ii(1)+1,ii(2)  )
                    val = val + (   aa(0))*(   aa(1))*(1.-aa(2))*pf%var(ii(0)+1,ii(1)+1,ii(2)  )
                    val = val + (1.-aa(0))*(1.-aa(1))*(   aa(2))*pf%var(ii(0)  ,ii(1)  ,ii(2)+1)
                    val = val + (   aa(0))*(1.-aa(1))*(   aa(2))*pf%var(ii(0)+1,ii(1)  ,ii(2)+1)
                    val = val + (1.-aa(0))*(   aa(1))*(   aa(2))*pf%var(ii(0)  ,ii(1)+1,ii(2)+1)
                    val = val + (   aa(0))*(   aa(1))*(   aa(2))*pf%var(ii(0)+1,ii(1)+1,ii(2)+1)
                    field(i,j,k) = val
                end do
            end do
        end do
    end subroutine interpolate_elem_field

    subroutine get_rotation_to_pole(theta ,phi, M)
    !!! compute the rotation matrix from the chunk center to the pole
    !!! theta :: chunk center colatitude (radians)
    !!! phi :: chunk center longitude (radians)
        implicit none
        real(fpp), intent(in) :: theta, phi
        real(fpp), intent(out), dimension(0:2,0:2) :: M
        real(fpp) :: ct, cp, st, sp

        ct = cos(theta*M_PI/180.0_fpp)
        cp = cos(phi*M_PI/180.0_fpp)
        st = sin(theta*M_PI/180.0_fpp)
        sp = sin(phi*M_PI/180.0_fpp)

        M(0,0) = ct*cp
        M(0,1) = ct*sp
        M(0,2) = -st
        M(1,0) = -sp
        M(1,1) = cp
        M(1,2) = 0.
        M(2,0) = st*cp
        M(2,1) = st*sp
        M(2,2) = ct

    end subroutine get_rotation_to_pole


    subroutine get_rotation_from_pole(theta ,phi, M)
    !!! compute the rotation matrix from the pole to the chunk center
    !!! theta :: chunk center colatitude (radians)
    !!! phi :: chunk center longitude (radians)
        implicit none
        real(fpp), intent(in) :: theta, phi
        real(fpp), intent(out), dimension(0:2,0:2) :: M
        real(fpp) :: ct, cp, st, sp 

        ct = cos(theta*M_PI/180.0_fpp)
        cp = cos(phi*M_PI/180.0_fpp)
        st = sin(theta*M_PI/180.0_fpp)
        sp = sin(phi*M_PI/180.0_fpp)

        M(0,0) = ct*cp
        M(0,1) = -sp
        M(0,2) = st*cp
        M(1,0) = ct*sp
        M(1,1) = cp
        M(1,2) = st*sp
        M(2,0) = -st
        M(2,1) = 0.
        M(2,2) = ct

    end subroutine get_rotation_from_pole


    subroutine cart2sph_to_check(xyz, rtp, rad2deg)
        implicit none
        logical, intent(in) :: rad2deg
        real(fpp), intent(in), dimension(0:2) :: xyz
        real(fpp), intent(out), dimension(0:2) :: rtp
        real(fpp) :: dx

        rtp(0) = dsqrt(xyz(0)**2+xyz(1)**2+xyz(2)**2)
        if (rtp(0)==0.0_fpp) then
            rtp(1) = 0.0_fpp
            rtp(2) = 0.0_fpp
        else
            dx = xyz(2)/rtp(0)
            if (dx>=1.0_fpp) then
                rtp(1) = 0.0_fpp
            else if (dx<=-1.0_fpp) then
                rtp(1) = M_PI
            else
                rtp(1) = dacos(dx)
            end if
            if ((rtp(1)==0.0_fpp).or.(rtp(1)==M_PI)) then
                rtp(2) = 0.0_fpp
            else
                dx = xyz(0)/(rtp(0)*dsin(rtp(1)))
                if (dx>1.0_fpp) then
                    rtp(2) = 0.0_fpp
                else if (dx<1.0_fpp) then
                    rtp(2) = M_PI
                else
                    rtp(2) = dacos(dx)
                    if (xyz(1)<0.0_fpp) rtp(2) = 2*M_PI-rtp(2)
                end if
            end if
        end if

        if (rad2deg) then
            rtp(1)=90.0_fpp-rtp(1)*180.0_fpp/M_PI
            rtp(2)=rtp(2)*180.0_fpp/M_PI
        end if

    end subroutine cart2sph_to_check

    subroutine cart2sph(xyz, rtp, rad2deg)
        implicit none
        logical, intent(in) :: rad2deg
        real(fpp), intent(in), dimension(0:2) :: xyz
        real(fpp), intent(out), dimension(0:2) :: rtp

        rtp(0) = sqrt(xyz(0)**2+xyz(1)**2+xyz(2)**2)
        rtp(1) = acos(xyz(2)/rtp(0))
        rtp(2) = atan2(xyz(1),xyz(0))

        if (rad2deg) then
            rtp(1)=90.0_fpp-rtp(1)*180.0_fpp/M_PI
            rtp(2)=rtp(2)*180.0_fpp/M_PI
        end if

    end subroutine cart2sph



end module build_prop_files
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
