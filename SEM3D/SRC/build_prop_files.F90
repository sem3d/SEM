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
        case(MATDEF_VP_VS_RHO,MATDEF_VP_VS_RHO_D)
            mat%prop_field(1)%propName = "Vp"
            mat%prop_field(2)%propName = "Vs"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_E_NU_RHO,MATDEF_E_NU_RHO_D)
            mat%prop_field(1)%propName = "E"
            mat%prop_field(2)%propName = "Nu"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_LAMBDA_MU_RHO,MATDEF_LAMBDA_MU_RHO_D)
            mat%prop_field(1)%propName = "Lambda"
            mat%prop_field(2)%propName = "Mu"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_KAPPA_MU_RHO,MATDEF_KAPPA_MU_RHO_D)
            mat%prop_field(1)%propName = "Kappa"
            mat%prop_field(2)%propName = "Mu"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_HOOKE_RHO,MATDEF_HOOKE_RHO_D)
            stop "Not supported yet"
        case(MATDEF_NLKP_VS_RHO,MATDEF_NLKP_VS_RHO_D)
            mat%prop_field(1)%propName = "NLkin"
            mat%prop_field(2)%propName = "Vs"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_NU_VS_RHO,MATDEF_NU_VS_RHO_D)
            mat%prop_field(1)%propName = "Nu"
            mat%prop_field(2)%propName = "Vs"
            mat%prop_field(3)%propName = "Rho"
        case(MATDEF_VTI_ANISO)
            mat%prop_field(1)%propName = "Vpv"
            mat%prop_field(2)%propName = "Vsv"
            mat%prop_field(3)%propName = "Rho"
            mat%prop_field(4)%propName = "Vph"
            mat%prop_field(5)%propName = "Vsh"
            mat%prop_field(6)%propName = "Eta"
            mat%prop_field(7)%propName = "Qkappa"
            mat%prop_field(8)%propName = "Qmu"
        case(MATDEF_FLUID_ANISO, CSTAR_FLUID)
            mat%prop_field(1)%propName = "K11"
            mat%prop_field(2)%propName = "K22"
            mat%prop_field(3)%propName = "K33"
            mat%prop_field(4)%propName = "K12"
            mat%prop_field(5)%propName = "K13"
            mat%prop_field(6)%propName = "K23"
            mat%prop_field(7)%propName = "Rho"
        case(MATDEF_HOOKE_ANISO,CSTAR)
            mat%prop_field(1)%propName = "C11"
            mat%prop_field(2)%propName = "C22"
            mat%prop_field(3)%propName = "C33"
            mat%prop_field(4)%propName = "C44"
            mat%prop_field(5)%propName = "C55"
            mat%prop_field(6)%propName = "C66"
            mat%prop_field(7)%propName = "C12"
            mat%prop_field(8)%propName = "C13"
            mat%prop_field(9)%propName = "C14"
            mat%prop_field(10)%propName = "C15"
            mat%prop_field(11)%propName = "C16"
            mat%prop_field(12)%propName = "C23"
            mat%prop_field(13)%propName = "C24"
            mat%prop_field(14)%propName = "C25"
            mat%prop_field(15)%propName = "C26"
            mat%prop_field(16)%propName = "C34"
            mat%prop_field(17)%propName = "C35"
            mat%prop_field(18)%propName = "C36"
            mat%prop_field(19)%propName = "C45"
            mat%prop_field(20)%propName = "C46"
            mat%prop_field(21)%propName = "C56"
            mat%prop_field(22)%propName = "Rho"
        !case(CSTAR)
        !    mat%prop_field(1)%propName  = "C11"
        !    mat%prop_field(2)%propName  = "C12"
        !    mat%prop_field(3)%propName  = "C13"
        !    mat%prop_field(4)%propName  = "C14"
        !    mat%prop_field(5)%propName  = "C15"
        !    mat%prop_field(6)%propName  = "C16"
        !    mat%prop_field(7)%propName  = "C22"
        !    mat%prop_field(8)%propName  = "C23"
        !    mat%prop_field(9)%propName  = "C24"
        !    mat%prop_field(10)%propName = "C25"
        !    mat%prop_field(11)%propName = "C26"
        !    mat%prop_field(12)%propName = "C33"
        !    mat%prop_field(13)%propName = "C34"
        !    mat%prop_field(14)%propName = "C35"
        !    mat%prop_field(15)%propName = "C35"
        !    mat%prop_field(16)%propName = "C44"
        !    mat%prop_field(17)%propName = "C45"
        !    mat%prop_field(18)%propName = "C46"
        !    mat%prop_field(19)%propName = "C55"
        !    mat%prop_field(20)%propName = "C56"
        !    mat%prop_field(21)%propName = "C66"
        !    mat%prop_field(22)%propName = "Rho"
        end select

        select case(mat%deftype)
        case(MATDEF_VP_VS_RHO_D,MATDEF_E_NU_RHO_D,MATDEF_LAMBDA_MU_RHO_D,&
            MATDEF_KAPPA_MU_RHO_D,MATDEF_HOOKE_RHO_D,MATDEF_NLKP_VS_RHO_D,&
            MATDEF_NU_VS_RHO_D)
            mat%prop_field(4)%propName = "Qp"
            mat%prop_field(5)%propName = "Qs"
        end select

        select case (mat%deftype)
            case(CSTAR, CSTAR_FLUID)
                call init_prop_file_field_Cstar(mat)
            case default
                do i = 1,mat%n_prop
                    call init_prop_file_field(mat, mat%prop_field(i))
                end do
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
            write(*,*) "Missing properties in h5file",pf%propName !! MC
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
            write(*,*) "k: ", k, "imin: ",pf%imin(k), "imax: ",pf%imax(k), "step: ",pf%step(k)
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

    subroutine init_prop_file_field_Cstar(mat)
        !use, intrinsic :: iso_c_binding
        type(subdomain), intent(inout) :: mat
        real(fpp), dimension(0:2) :: xxr
 
        integer, parameter :: num_comp = 22
        integer, parameter :: ntriu = 21
        integer :: i, j, k, a, b, kk, comp, ielx, iely, ielz, idx, Nd
        integer :: ncomp
        integer :: nelx, nely, nelz, nelements
        integer :: iheader, reclen, icode, ndeg, len_int,len_real
        real(4) :: rxel, ryel, rzel, xs, ys, zs, xs_whole_domain, ys_whole_domain, zs_whole_domain
        integer :: nx, ny, nz, len
        integer :: triu_idx(2, ntriu)
        integer :: iindo(2, ntriu)
        integer, dimension(2) :: ix, iy, iz
        integer :: cx, cy, cz
        integer :: mapping(22)
         
        real(4), dimension(:,:,:,:)  , allocatable :: buffer
        real(4), dimension(:,:,:,:), allocatable :: elemtmp
        integer, dimension(3) :: list
        integer :: l, m
        integer :: irec, ier, nelem_needed
        integer :: rg

        ! File I/O
        integer :: unit, ios
        character(len=256) :: filename

        if (.not. mat%present) return
        call MPI_Comm_rank(MPI_COMM_WORLD, rg, ier)
        
        inquire(iolength=len_int) i
        inquire(iolength=len_real)rxel
        len=8*len_int+6*len_real

        open(newunit=unit, file= trim(mat%prop_field(1)%propFilePath), form='unformatted', access='direct',status='old', action='read', iostat=ios, recl=len)
        if (ios /= 0) then
          write(*,*) "Could not open file:", trim(mat%prop_field(1)%propFilePath)
          stop 1
        end if
       
        read(unit,rec=1) icode,iheader,len,Nd,ndeg,nelx,nely,nelz,rxel,ryel,rzel,xs_whole_domain,ys_whole_domain,zs_whole_domain
        xs=rxel*nelx
        ys=ryel*nely
        zs=rzel*nelz
        close(unit)

        ! 22 = elastic (Nd=6): C11..C66 + rho
        ! 7  = acoustic aniso (Nd=3): K11,K22,K33,K12,K13,K23 + rho
        if (Nd*(Nd+1)/2+1 /= 22 .and. Nd*(Nd+1)/2+1 /= 7) then
            write(*,*) "Error: unsupported Cstar component count on proc", rg, "value =", Nd*(Nd+1)/2+1
            stop "Unsupported Cstar dimension (expected 7 for acoustic or 22 for elastic)"
        end if
                                          
        if (icode /= -82) stop 'Unsupported icode in reading CStar file'
        ncomp = Nd*(Nd+1)/2+1

        ix(0) = floor(mat%MinBound_Loc(0)/rxel)
        ix(1) = ceiling(mat%MaxBound_Loc(0)/rxel)

        iy(0) = floor(mat%MinBound_Loc(1)/ryel)
        iy(1) = ceiling(mat%MaxBound_Loc(1)/ryel)

        iz(0) = floor(mat%MinBound_Loc(2)/rzel)
        iz(1) = ceiling(mat%MaxBound_Loc(2)/rzel)
        
        cx = ix(1)-ix(0)
        cy = iy(1)-iy(0)
        cz = iz(1)-iz(0)
        nelem_needed = (ix(1)-ix(0))*(iy(1)-iy(0))*(iz(1)-iz(0))
        write(*,*) "Cstar read header on proc", rg
        write(*,*) "icode :", icode
        write(*,*) "iheader :", iheader
        write(*,*) "len :", len
        write(*,*) "Nd :", Nd
        write(*,*) "ndeg :", ndeg
        write(*,*) "n el :", nelx, nely, nelz, cx, cy, cz
        write(*,*) "step in :", rxel, ryel, rzel
        write(*,*) "domain size :", xs, ys, zs, xs_whole_domain,ys_whole_domain, zs_whole_domain
        write(*,*) "finished "
        write(*,*) "BB min from mat strucutre :", mat%MinBound_Loc
        write(*,*) "BB max from mat strucutre :", mat%MaxBound_Loc

        write(*,*) "x y z start in struuctured: ", ix(0), iy(0), iz(0)
        write(*,*) "x y z end in structured: ", ix(1), iy(1), iz(1)
        write(*,*) "Number of elements to be read :", nelem_needed
       
        open(newunit=unit, file=trim(mat%prop_field(1)%propFilePath),status='old',access='direct',form='unformatted',recl=len, iostat=ios)
        if (ios /= 0) then
            write(*,*) "Could not open file:",trim(mat%prop_field(1)%propFilePath)
            stop 1
        end if
                                                 
        allocate(buffer(Nd*(Nd+1)/2+1,ndeg+1,ndeg+1,ndeg+1) )
        !allocate(elemtmp(ndeg+1,ndeg+1,ndeg+1,22,nelx+1,nely+1,nelz+1), stat=ios)
        allocate(elemtmp(ncomp,cx+1,cy+1,cz+1), stat=ios)
        if (ios /= 0) then
            print *, "Allocation failed!"
            stop
        endif
        buffer = -1
        elemtmp = 10
        cz = 0
        do ielz = iz(0), iz(1)-1
            cy = 0
            do iely = iy(0), iy(1)-1
                cx = 0
                do ielx = ix(0), ix(1)-1
                    list(0) = ielx+1
                    list(1) = iely+1
                    list(2) = ielz+1
                    irec=iheader+list(0)+(list(1)-1)*nelx+(list(2)-1)*nelx*nely
                    !write(*,*) irec, list
                    read(unit,rec=irec,iostat=ier) buffer
                    if (ier /= 0) then
                        write(*,*) 'READ error at irec=', irec, ' ier=', ier
                    end if
                    do kk=1,Nd*(Nd+1)/2+1
                        !if (buffer(kk,1,1,1)<0) then
                        !    !buffer(kk,1,1,1) = 0
                        !endif
                        !if (buffer(kk,1,1,1) > 1d13) then
                        !    !buffer(kk,1,1,1) = 1d13
                        !endif
                        elemtmp(kk,cx+1,cy+1,cz+1)=buffer(kk,1,1,1)
                    enddo
                    cx = cx+1
                end do
                cy = cy+1
            end do
            cz = cz+1     
        end do
     
        close(unit)
        
        cx = ix(1)-ix(0)
        cy = iy(1)-iy(0)
        cz = iz(1)-iz(0)
        
        if (ncomp == 22) then
            ! Elastic Nd=6: upper-triangle row-major C11,C12,...,C66,Rho
            ! mapped to prop_field order C11,C22,C33,C44,C55,C66,C12,...,Rho
            mapping(1:22) = [1, 7, 8, 9, 10, 11, 2, 12, 13, 14, 15, 3, 16, 17, 18, 4, 19, 20, 5, 21, 6, 22]
        else
            ! Acoustic Nd=3: upper-triangle row-major K11,K12,K13,K22,K23,K33,Rho
            ! mapped to prop_field order K11,K22,K33,K12,K13,K23,Rho
            mapping(1:7) = [1, 4, 5, 2, 6, 3, 7]
        end if
        ! fill mat struct
        do j = 1,mat%n_prop
            i = mapping(j)
            mat%prop_field(i)%MinBound(0) = ix(0)*rxel
            mat%prop_field(i)%MaxBound(0) = ix(1)*rxel
            mat%prop_field(i)%MinBound(1) = iy(0)*ryel
            mat%prop_field(i)%MaxBound(1) = iy(1)*ryel
            mat%prop_field(i)%MinBound(2) = iz(0)*rzel
            mat%prop_field(i)%MaxBound(2) = iz(1)*rzel

            !mat%prop_field(i)%NN(0) = (ix(1)-ix(0))*(ndeg+1) + 1
            !mat%prop_field(i)%NN(1) = (iy(1)-iy(0))*(ndeg+1) + 1
            !mat%prop_field(i)%NN(2) = (iz(1)-iz(0))*(ndeg+1) + 1
            mat%prop_field(i)%NN(0) = (ix(1)-ix(0))
            mat%prop_field(i)%NN(1) = (iy(1)-iy(0))
            mat%prop_field(i)%NN(2) = (iz(1)-iz(0))
            !write(*,*) mat%prop_field(i)%MinBound
            !write(*,*) mat%prop_field(i)%NN

            do k = 0,2
                mat%prop_field(i)%imin(k) = gindex(mat%MinBound_Loc(k), mat%prop_field(i)%NN(k), mat%prop_field(i)%MinBound(k), mat%prop_field(i)%MaxBound(k))
                mat%prop_field(i)%imax(k) = gindex(mat%Maxbound_Loc(k), mat%prop_field(i)%NN(k), mat%prop_field(i)%MinBound(k), mat%prop_field(i)%MaxBound(k))+1
                mat%prop_field(i)%step(k) = (mat%prop_field(i)%MaxBound(k)-mat%prop_field(i)%MinBound(k))/(mat%prop_field(i)%NN(k)-1)
                !write(*,*) "prop: ", i, "k: ", k, "imin: ",mat%prop_field(i)%imin(k), "imax: ",mat%prop_field(i)%imax(k), "step: ",mat%prop_field(i)%step(k)
                if ((mat%prop_field(i)%imax(k)-mat%prop_field(i)%imin(k))<1) mat%prop_field(i)%imin(k) = mat%prop_field(i)%imax(k)-1
                if (mat%prop_field(i)%imin(k)<0) then
                    mat%prop_field(i)%imin(k) = 0
                    if (mat%prop_field(i)%imax(k)<1) mat%prop_field(i)%imax(k) = 1
                end if
                if (mat%prop_field(i)%imax(k)>=mat%prop_field(i)%NN(k)) then
                    mat%prop_field(i)%imax(k) = mat%prop_field(i)%NN(k)-1
                    if (mat%prop_field(i)%imin(k)>(mat%prop_field(i)%imax(k)-1)) mat%prop_field(i)%imin(k) = mat%prop_field(i)%imax(k)-1
                endif
                if ((mat%prop_field(i)%imax(k)-mat%prop_field(i)%imin(k))<1) mat%prop_field(i)%imax(k) = mat%prop_field(i)%imin(k)
                if (mat%is_sph) then
                    mat%prop_field(i)%imin(k)=0
                    mat%prop_field(i)%imax(k)=mat%prop_field(i)%NN(k)-1
                end if
            end do
            ! --------    
            allocate(mat%prop_field(i)%var(0:cx-1,0:cy-1,0:cz-1)) 
        end do
        
        ! populate the elements for all properties with optimal cache locality
        do ielz = 0, cz -1
            do iely = 0, cy -1
                do ielx = 0, cx -1
                    do j = 1,mat%n_prop
                        i = mapping(j)
                        mat%prop_field(i)%var(ielx,iely,ielz) = elemtmp(j,ielx+1,iely+1,ielz+1)
                    end do
                end do
            end do
        end do
        
        do j = 1,mat%n_prop
            i = mapping(j)
            ! ------- check this thing here...
            !if (i == 4) .OR. (i==5) .OR. (i==6) then
            !    mat%prop_field(i)%var = mat%prop_field(i)%var/2
            !end if 
        end do
        
        !interp to reagular grid
        write(*,*) "Cstar readed routine end on proc", rg
        !stop 1

        deallocate(buffer)
        deallocate(elemtmp)
    end subroutine init_prop_file_field_Cstar

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
