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
        type(subdomain), intent(inout) :: mat

        select case(mat%deftype)
        case(MATDEF_VP_VS_RHO,MATDEF_VP_VS_RHO_D)
            mat%pf(1)%propName = "Vp"
            mat%pf(2)%propName = "Vs"
        case(MATDEF_E_NU_RHO,MATDEF_E_NU_RHO_D)
            mat%pf(1)%propName = "E"
            mat%pf(2)%propName = "Nu"
        case(MATDEF_LAMBDA_MU_RHO,MATDEF_LAMBDA_MU_RHO_D)
            mat%pf(1)%propName = "Lambda"
            mat%pf(2)%propName = "Mu"
        case(MATDEF_KAPPA_MU_RHO,MATDEF_KAPPA_MU_RHO_D)
            mat%pf(1)%propName = "Kappa"
            mat%pf(2)%propName = "Mu"
        case(MATDEF_HOOKE_RHO,MATDEF_HOOKE_RHO_D)
            stop "Not supported yet"
        case(MATDEF_NLKP_VS_RHO,MATDEF_NLKP_VS_RHO_D)
            mat%pf(1)%propName = "NLkin"
            mat%pf(2)%propName = "Vs"
        case(MATDEF_NU_VS_RHO,MATDEF_NU_VS_RHO_D)
            mat%pf(1)%propName = "Nu"
            mat%pf(2)%propName = "Vs"
        end select
        mat%pf(3)%propName = "Rho"
        call init_prop_file_field(mat, mat%pf(1))
        call init_prop_file_field(mat, mat%pf(2))
        call init_prop_file_field(mat, mat%pf(3))
        select case(mat%deftype)
        case(MATDEF_VP_VS_RHO_D,MATDEF_E_NU_RHO_D,MATDEF_LAMBDA_MU_RHO_D,&
            MATDEF_KAPPA_MU_RHO_D,MATDEF_HOOKE_RHO_D,MATDEF_NLKP_VS_RHO_D,&
            MATDEF_NU_VS_RHO_D)
            mat%pf(4)%propName = "Qp"
            mat%pf(5)%propName = "Qs"
            call init_prop_file_field(mat, mat%pf(4))
            call init_prop_file_field(mat, mat%pf(5))
        end select
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

        do k = 0,2
            pf%imin(k) = gindex(mat%MinBound_Loc(k), pf%NN(k), pf%MinBound(k), pf%MaxBound(k))
            pf%imax(k) = gindex(mat%MaxBound_Loc(k), pf%NN(k), pf%MinBound(k), pf%MaxBound(k))+1
            pf%step(k) = (pf%MaxBound(k)-pf%MinBound(k))/(pf%NN(k)-1)
            if (pf%imin(k)<0) pf%imin(k) = 0
            if (pf%imax(k)<pf%imin(k)) pf%imax(k) = pf%imin(k)
            if (pf%imax(k)>=pf%NN(k)) pf%imax(k) = pf%NN(k)-1
            if (pf%imin(k)>pf%imax(k)) pf%imin(k) = pf%imax(k)
        end do

        call read_subset_3d_real(grp_id, "samples", pf%imin, pf%imax, pf%var)
        if (subgrp) call H5Gclose_f(grp_id, hdferr)
        call H5Fclose_f(file_id, hdferr)
    end subroutine init_prop_file_field

    subroutine cleanup_prop_file(mat)
        type(subdomain) :: mat
        if (.not. mat%present) return
        deallocate(mat%pf(1)%var)
        deallocate(mat%pf(2)%var)
        deallocate(mat%pf(3)%var)
    end subroutine cleanup_prop_file


    subroutine interpolate_elem_field(Tdomain, specel, mat, pf, field)
        type(domain), intent(inout) :: Tdomain
        type(element), intent(inout) :: specel
        type(subdomain), intent(in) :: mat
        type(PropertyField), intent(in) :: pf
        !
        real(fpp), dimension(0:mat%ngll-1,0:mat%ngll-1,0:mat%ngll-1) :: field
        real(fpp), dimension(0:2) :: xx ! node coord
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
                        pos = (xx(n)-pf%MinBound(n))/pf%step(n)
                        ii(n) = floor(pos)
                        if (ii(n)  < pf%imin(n)) ii(n) = pf%imin(n)
                        if (ii(n) >= pf%imax(n)) ii(n) = pf%imax(n)-1
                        aa(n) = pos-ii(n)
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
