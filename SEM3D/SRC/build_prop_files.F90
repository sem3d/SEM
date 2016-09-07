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

        call init_prop_file_field(mat, mat%pf(1))
        call init_prop_file_field(mat, mat%pf(2))
        call init_prop_file_field(mat, mat%pf(3))
    end subroutine init_prop_file

    subroutine init_prop_file_field(mat, pf)
        use hdf5
        use sem_hdf5
        type(subdomain), intent(inout) :: mat
        type(PropertyField), intent(inout) :: pf

        integer, dimension(0:2) :: imin, imax
        integer :: k, hdferr
        integer(HID_T) :: file_id, grp_id
        integer(HSIZE_T), dimension(:), allocatable :: dims
        call init_hdf5()
        call h5fopen_f(trim(pf%propFilePath), H5F_ACC_RDONLY_F, file_id, hdferr) !Open File
        if(hdferr /= 0) then
            write(*,*) "Could not open file:", trim(pf%propFilePath)
            stop 1
        end if

        grp_id = file_id
        call read_attr_real_vec(grp_id, "xMinGlob", pf%MinBound)
        call read_attr_real_vec(grp_id, "xMaxGlob", pf%MaxBound)
        call read_dims(grp_id, "samples", dims)
        pf%NN = dims
        ! On va calculer les indices i0,i1 j0,j1,k0,k1 tels que
        ! i0 plus grand entier tel que x(i0)<MinBound_loc(0), 
        ! i1 plus petit entier tel que x(i1)>MaxBound_loc(0), etc...
        do k = 0,2
            pf%imin(k) = gindex(mat%MinBound_Loc(k), pf%NN(k), pf%MinBound(k), pf%MaxBound(k))
            pf%imax(k) = gindex(mat%MaxBound_Loc(k), pf%NN(k), pf%MinBound(k), pf%MaxBound(k))
            pf%step(k) = (pf%MaxBound(k)-pf%MinBound(k))/pf%NN(k)
        end do
        call read_subset_3d_real(grp_id, "samples", pf%imin, pf%imax, pf%var)

    end subroutine init_prop_file_field
    
    subroutine cleanup_prop_file(mat)
        type(subdomain) :: mat
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
                        ii(n) = int(pos)
                        if (ii(n) < pf%imin(n)) ii(n) = pf%imin(n)
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
                     
