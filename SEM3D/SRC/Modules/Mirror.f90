!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module smirror
  use constants, only : fpp
  implicit none
#include "index.h"
  public :: init_mirror, create_mirror_files, dump_mirror_sl, load_mirror_sl, dump_mirror_fl, &
    load_mirror_fl
  private

  integer, dimension(:,:,:,:), allocatable :: map2glltot_sl, map2glltot_fl
  real(fpp), dimension(:), allocatable :: bspl_tmp
  real(fpp), dimension(:,:,:), allocatable :: displ_sl, force_sl
  real(fpp), dimension(:,:), allocatable :: displ_fl, force_fl
  real(fpp) :: d_t,d_tm
  integer :: rnk, recp, n_spl, n_dcm, n_t, n_tm, n_gll, n_glltot_sl, n_glltot_fl

contains

subroutine init_mirror(Tdomain, surf)
  use sspline
  use sdomain
  implicit none
  type(domain), intent(inout) :: Tdomain
  type(SurfaceT), intent(inout) :: surf

  rnk = Tdomain%rank

  if (Tdomain%mirror_type==0) then
    call create_mirror_h5()
    call dump_param_mirror(Tdomain)
  else
    call load_param_mirror(Tdomain)
  endif

  if (n_dcm>1) then
    allocate(bspl_tmp((n_spl+1)*n_dcm+1))
    call bmn(bspl_tmp, n_dcm, n_spl)
  endif

  if (surf%surf_sl%nbtot/=0) then
    call init_mirror_sl(Tdomain, surf%surf_sl)
  endif
  if (surf%surf_fl%nbtot/=0) then
    call init_mirror_fl(Tdomain, surf%surf_fl)
  endif

end subroutine init_mirror

subroutine dump_param_mirror(Tdomain)
  use sdomain
  implicit none
  type(domain), intent(inout) :: Tdomain
  real(fpp) :: fmax

  n_spl = Tdomain%config%mirror_nspl
  d_t = Tdomain%TimeD%dtmin
  n_t = Tdomain%TimeD%ntimeMax

  fmax = Tdomain%config%mirror_fmax
  if (fmax>0.) then
    n_dcm = int(0.25/(fmax*d_t))
    n_tm = int(n_t/n_dcm)+(n_spl+1)
    if (Tdomain%rank==0) write(*,'("--> SEM : mirror decimation : ",2i6)') n_dcm,n_tm
  else
    n_dcm = 1
    n_tm = n_t
  endif

  d_tm = d_t*n_dcm

  call write_mirror_h5_attr()

end subroutine dump_param_mirror

subroutine load_param_mirror(Tdomain)
  use sdomain
  implicit none
  type(domain), intent(inout) :: Tdomain

  n_spl = Tdomain%config%mirror_nspl
  d_t = Tdomain%TimeD%dtmin
  n_t = Tdomain%TimeD%ntimeMax

  call read_mirror_h5_attr()

  n_dcm = -1
  if (modulo(d_tm,d_t)==0.d0) n_dcm = int(d_tm/d_t)

  if (Tdomain%rank==0) write(*,'("--> SEM : mirror decimation : ",2i6)') n_dcm,n_tm

end subroutine load_param_mirror

subroutine init_mirror_sl(Tdomain, surf)
  use sdomain
  implicit none
  type(domain), intent (inout) :: Tdomain
  type(surf_num), intent(inout) :: surf
  integer :: n_elmtot

  n_elmtot = Tdomain%sdom%nblocks
  n_gll = Tdomain%sdom%ngll
  if (surf%nbtot/=0) then
    allocate(map2glltot_sl(0:n_elmtot-1,0:n_gll-1,0:n_gll-1,0:n_gll-1))
    call map_mirror_sl(Tdomain, surf)
    if (n_glltot_sl>0) then
      write(*,'("--> SEM : mirror nodes_sl : ",i3,i6)') rnk,n_glltot_sl
      allocate(Tdomain%sdom%mirror_sl%map(0:n_elmtot-1,0:n_gll-1,0:n_gll-1,0:n_gll-1))
      allocate(Tdomain%sdom%mirror_sl%fields(6,n_glltot_sl))
      allocate(displ_sl(3,n_glltot_sl,n_spl+1))
      allocate(force_sl(3,n_glltot_sl,n_spl+1))
      Tdomain%sdom%mirror_sl%n_glltot = n_glltot_sl
      Tdomain%sdom%mirror_sl%n_gll = n_gll
      Tdomain%sdom%mirror_sl%map = map2glltot_sl
      displ_sl = 0.
      force_sl = 0.
    else
      Tdomain%sdom%mirror_sl%n_glltot = n_glltot_sl
    endif
    deallocate(map2glltot_sl)
  endif

end subroutine init_mirror_sl

subroutine init_mirror_fl(Tdomain, surf)
  use sdomain
  implicit none
  type(domain), intent (inout) :: Tdomain
  type(surf_num), intent(inout) :: surf
  integer :: n_elmtot

  n_elmtot = Tdomain%fdom%nblocks
  n_gll = Tdomain%fdom%ngll
  if (surf%nbtot/=0) then
    allocate(map2glltot_fl(0:n_elmtot-1,0:n_gll-1,0:n_gll-1,0:n_gll-1))
    call map_mirror_fl(Tdomain, surf)
    if (n_glltot_fl>0) then
      write(*,'("--> SEM : mirror nodes_fl : ",i3,i6)') rnk,n_glltot_fl
      allocate(Tdomain%fdom%mirror_fl%map(0:n_elmtot-1,0:n_gll-1,0:n_gll-1,0:n_gll-1))
      allocate(Tdomain%fdom%mirror_fl%fields(2,n_glltot_fl))
      allocate(displ_fl(n_glltot_fl,n_spl+1))
      allocate(force_fl(n_glltot_fl,n_spl+1))
      Tdomain%fdom%mirror_fl%n_glltot = n_glltot_fl
      Tdomain%fdom%mirror_fl%n_gll = n_gll
      Tdomain%fdom%mirror_fl%map = map2glltot_fl
      displ_fl = 0.
      force_fl = 0.
    else
      Tdomain%fdom%mirror_fl%n_glltot = n_glltot_fl
    endif
    deallocate(map2glltot_fl)
  endif

end subroutine init_mirror_fl

subroutine map_mirror_sl(Tdomain, surf)
  use sdomain
  implicit none
  type(domain), intent(in) :: Tdomain
  type(surf_num), intent(in) :: surf
  character(len=80) :: fname
  integer :: i_surf,i_node,i_dir,f,n,m,e,i,j,k
  real(fpp) :: weight
  integer, dimension(:), allocatable :: ind
  logical, dimension(:,:,:,:), allocatable :: tmp
  real(fpp), dimension(:,:,:,:,:), allocatable :: tmp_n
  real(fpp), dimension(:,:), allocatable :: nodesinfos
  real(fpp), dimension(0:2,0:3) :: cnodes

  allocate(tmp(0:Tdomain%sdom%nblocks-1,0:n_gll-1,0:n_gll-1,0:n_gll-1))
  allocate(tmp_n(0:Tdomain%sdom%nblocks-1,0:n_gll-1,0:n_gll-1,0:n_gll-1,0:2))
  allocate(ind(0:Tdomain%sdom%nblocks-1))
  tmp = .false.
  tmp_n = 0.d0
  do i_surf = 0,surf%n_faces-1
    f = surf%if_faces(i_surf)
    n = surf%if_norm(i_surf)
    if (n>0) then
      m = 0
      e = Tdomain%sFace(f)%elem_0
    else
      m = n_gll-1
      e = Tdomain%sFace(f)%elem_1
    endif
    if (e==-1) cycle
    ind(Tdomain%specel(e)%lnum) = e
    do i_node = 0,3
      cnodes(0:2,i_node) = Tdomain%Coord_Nodes(0:2,Tdomain%sFace(f)%inodes(i_node))
    enddo
    call mirror_face_normal(cnodes,i_dir)
    select case (i_dir)
      case(0)
        tmp(Tdomain%specel(e)%lnum,m,:,:) = .true.
        tmp_n(Tdomain%specel(e)%lnum,m,:,:,0) = m*2.d0/(n_gll-1)-1.d0
      case(1)
        tmp(Tdomain%specel(e)%lnum,:,m,:) = .true.
        tmp_n(Tdomain%specel(e)%lnum,:,m,:,1) = m*2.d0/(n_gll-1)-1.d0
      case(2)
        tmp(Tdomain%specel(e)%lnum,:,:,m) = .true.
        tmp_n(Tdomain%specel(e)%lnum,:,:,m,2) = m*2.d0/(n_gll-1)-1.d0
    end select
  enddo
  map2glltot_sl = -1
  n_glltot_sl = 0
  allocate(nodesinfos(0:(Tdomain%sdom%nblocks*n_gll*n_gll*n_gll-1),0:6))
  do k = 0,n_gll-1
    do j = 0,n_gll-1
      do i = 0,n_gll-1
        do e = 0,Tdomain%sdom%nblocks-1
          if (tmp(e,i,j,k)) then
            n_glltot_sl = n_glltot_sl+1
            map2glltot_sl(e,i,j,k) = n_glltot_sl
            i_node = Tdomain%specel(ind(e))%Iglobnum(i,j,k)
            nodesinfos(n_glltot_sl-1,0:2) = Tdomain%GlobCoord(0:2,i_node)
            nodesinfos(n_glltot_sl-1,3:5) = tmp_n(e,i,j,k,1:3)
            nodesinfos(n_glltot_sl-1,6) = Tdomain%sdom%GLLw(i)*Tdomain%sdom%GLLw(j)*Tdomain%sdom%GLLw(k)
          endif
        enddo
      enddo
    enddo
  enddo
  if (n_glltot_sl>0) call dump_coupling_nodes_sl(n_glltot_sl,nodesinfos)
  deallocate(nodesinfos)
  deallocate(ind)
  deallocate(tmp)
  deallocate(tmp_n)

end subroutine map_mirror_sl

subroutine dump_coupling_nodes_sl(n,nodesinfos)
  use semdatafiles
  implicit none
  integer, intent(in) :: n
  real(fpp), dimension(0:n-1,0:6), intent(in) :: nodesinfos
  character(len=MAX_FILE_SIZE) :: fnamef
  integer :: fu,i

  call semname_mirrorfile_nodes(rnk,fnamef)
  fu = 1111+rnk
  open(fu,file=fnamef)
  do i = 0,n-1
    write(fu,'(3e16.8,f16.8,3f5.1)') nodesinfos(i,0:6)
  enddo
  close(fu)

end subroutine dump_coupling_nodes_sl

subroutine map_mirror_fl(Tdomain, surf)
  use sdomain
  implicit none
  type(domain), intent(in) :: Tdomain
  type(surf_num), intent(in) :: surf
  integer :: i_surf,i_node,i_dir,f,n,m,e,i,j,k
  logical, dimension(:,:,:,:), allocatable :: tmp
  real(fpp), dimension(0:2,0:3) :: cnodes

  allocate(tmp(0:Tdomain%fdom%nblocks-1,0:n_gll-1,0:n_gll-1,0:n_gll-1))
  tmp = .false.
  do i_surf = 0,surf%n_faces-1
    f = surf%if_faces(i_surf)
    n = surf%if_norm(i_surf)
    if (n>0) then
      m = 0
      e = Tdomain%sFace(f)%elem_0
    else
      m = n_gll-1
      e = Tdomain%sFace(f)%elem_1
    endif
    if (e==-1) cycle
    do i_node = 0,3
      cnodes(0:2,i_node) = Tdomain%Coord_Nodes(0:2,Tdomain%sFace(f)%inodes(i_node))
    enddo
    call mirror_face_normal(cnodes, i_dir)
    select case (i_dir)
      case(0)
        tmp(Tdomain%specel(e)%lnum,m,:,:) = .true.
      case(1)
        tmp(Tdomain%specel(e)%lnum,:,m,:) = .true.
      case(2)
        tmp(Tdomain%specel(e)%lnum,:,:,m) = .true.
    end select
  enddo
  map2glltot_fl = -1
  n_glltot_fl = 0
  do k = 0,n_gll-1
    do j = 0,n_gll-1
      do i = 0,n_gll-1
        do e = 0,Tdomain%fdom%nblocks-1
          if (tmp(e,i,j,k)) then
            n_glltot_fl = n_glltot_fl+1
            map2glltot_fl(e,i,j,k) = n_glltot_fl
          endif
        enddo
      enddo
    enddo
  enddo
  deallocate(tmp)

end subroutine map_mirror_fl

subroutine mirror_face_normal(cnodes, fdir)
  use shape_geom_3d
  implicit none
  real(fpp), dimension(0:2, 0:3), intent(in) :: cnodes
  real(fpp), dimension(0:2) :: dfx, dfy,fnorm
  integer, intent(out) :: fdir
  integer :: i

  do i = 0,2
    dfx(i) = -cnodes(i,0)+cnodes(i,1)+cnodes(i,2)-cnodes(i,3)
    dfy(i) = -cnodes(i,0)-cnodes(i,1)+cnodes(i,2)+cnodes(i,3)
  enddo
  call cross_prod(dfx, dfy, fnorm)
  fdir = 0
  if (abs(fnorm(1))>abs(fnorm(fdir))) fdir = 1
  if (abs(fnorm(2))>abs(fnorm(fdir))) fdir = 2

end subroutine mirror_face_normal

subroutine create_mirror_files(Tdomain)
  use sdomain
  implicit none
  type(domain), intent(inout) :: Tdomain

  if (Tdomain%TimeD%NtimeMin==0.and.Tdomain%mirror_type==0) then
    call create_mirror_h5()
  endif

end subroutine create_mirror_files

subroutine create_mirror_h5()
  use semdatafiles
  use sem_hdf5
  use HDF5
  implicit none
  character (len=MAX_FILE_SIZE) :: fnamef
  character (len=40) :: dname
  integer(HID_T) :: fid, dset_id
  integer :: hdferr

  call init_hdf5()
  call semname_mirrorfile_h5(rnk, fnamef)
  call h5fcreate_f(fnamef, H5F_ACC_TRUNC_F, fid, hdferr)
  !! Displacement and force fields in solid
  dname = "Fields_sl"
  call create_dset_2d(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, &
    int(6,HSIZE_T), int(H5S_UNLIMITED_F,HSIZE_T), dset_id)
  call h5dclose_f(dset_id, hdferr)
  !! Displacement and force fields in fluid
  dname = "Fields_fl"
  call create_dset_2d(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, &
    int(2,HSIZE_T), int(H5S_UNLIMITED_F,HSIZE_T), dset_id)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine create_mirror_h5

subroutine dump_mirror_sl(dom, ntime)
  use champs_solid
  use sspline
  implicit none
  type(domain_solid), intent(inout) :: dom
  integer, intent(in) :: ntime
  integer :: i,j
  real(fpp), dimension(:), allocatable :: ddiag
  real(fpp), dimension(:,:), allocatable :: diag

  if (n_dcm==1) then
    call append_mirror_h5_sl(dom)
  else
    do i = 1,n_spl+1
      j = mod(ntime,n_dcm)+(n_spl+1-i)*n_dcm+1
      displ_sl(:,:,i) = displ_sl(:,:,i)+bspl_tmp(j)*dom%mirror_sl%fields(1:3,:)
      force_sl(:,:,i) = force_sl(:,:,i)+bspl_tmp(j)*dom%mirror_sl%fields(4:6,:)
    enddo
    if (mod(ntime,n_dcm)==n_dcm-1) then
      if (rnk==0) write(*,'("--> SEM : dump mirror_sl, iteration : ",i6.6)') ntime
      dom%mirror_sl%fields(1:3,:) = displ_sl(:,:,1)
      dom%mirror_sl%fields(4:6,:) = force_sl(:,:,1)
      call append_mirror_h5_sl(dom)
      do i = 1,n_spl
        displ_sl(:,:,i) = displ_sl(:,:,i+1)
        force_sl(:,:,i) = force_sl(:,:,i+1)
      enddo
      displ_sl(:,:,n_spl+1) = 0.
      force_sl(:,:,n_spl+1) = 0.
    endif
    if (ntime==n_t-1) then
      if (rnk==0) write(*,'("--> SEM : dump mirror_sl, iteration : ",i6.6)') ntime
      do i = 1,n_spl+1
        dom%mirror_sl%fields(1:3,:) = displ_sl(:,:,i)
        dom%mirror_sl%fields(4:6,:) = force_sl(:,:,i)
        call append_mirror_h5_sl(dom)
      enddo
      allocate(diag(n_spl+1,n_tm),ddiag(n_tm))
      call bmdiag(n_dcm,n_spl,n_tm,n_t,diag)
      call bchfac(diag,n_spl+1,n_tm,ddiag)
      call bchslv_sl(dom,diag,n_spl+1,n_tm)
      deallocate(diag,ddiag)
    endif
  endif

end subroutine dump_mirror_sl

subroutine dump_mirror_fl(dom, ntime)
  use champs_fluid
  use sspline
  implicit none
  type(domain_fluid), intent(inout) :: dom
  integer, intent(in) :: ntime
  integer :: i,j
  real(fpp), dimension(:), allocatable :: ddiag
  real(fpp), dimension(:,:), allocatable :: diag

  if (n_dcm==1) then
    call append_mirror_h5_fl(dom)
  else
    do i = 1,n_spl+1
      j = mod(ntime,n_dcm)+(n_spl+1-i)*n_dcm+1
      displ_fl(:,i) = displ_fl(:,i)+bspl_tmp(j)*dom%mirror_fl%fields(1,:)
      force_fl(:,i) = force_fl(:,i)+bspl_tmp(j)*dom%mirror_fl%fields(2,:)
    enddo
    if (mod(ntime,n_dcm)==n_dcm-1) then
      if (rnk==0) write(*,'("--> SEM : dump mirror_fl, iteration : ",i6.6)') ntime
      dom%mirror_fl%fields(1,:) = displ_fl(:,1)
      dom%mirror_fl%fields(2,:) = force_fl(:,1)
      call append_mirror_h5_fl(dom)
      do i = 1,n_spl
        displ_fl(:,i) = displ_fl(:,i+1)
        force_fl(:,i) = force_fl(:,i+1)
      enddo
      displ_fl(:,n_spl+1) = 0.
      force_fl(:,n_spl+1) = 0.
    endif
    if (ntime==n_t-1) then
      if (rnk==0) write(*,'("--> SEM : dump mirror_fl, iteration : ",i6.6)') ntime
      do i = 1,n_spl+1
        dom%mirror_fl%fields(1,:) = displ_fl(:,i)
        dom%mirror_fl%fields(2,:) = force_fl(:,i)
        call append_mirror_h5_fl(dom)
      enddo
      allocate(diag(n_spl+1,n_tm),ddiag(n_tm))
      call bmdiag(n_dcm,n_spl,n_tm,n_t,diag)
      call bchfac(diag,n_spl+1,n_tm,ddiag)
      call bchslv_fl(dom,diag,n_spl+1,n_tm)
      deallocate(diag,ddiag)
    endif
  endif

end subroutine dump_mirror_fl

subroutine load_mirror_sl(dom, ntime)
  use sspline
  use champs_solid
  implicit none
  type(domain_solid), intent(inout) :: dom
  integer, intent(in) :: ntime
  integer :: i,j,j1,j2,ntime_tmp,recn
  real(fpp) :: t,tnt,tmp

  if (dom%mirror_type==1) then
    ntime_tmp = ntime
  elseif (dom%mirror_type==2) then
    ntime_tmp = n_t-ntime-1
  else
    write(*,*) "STOP, unauthorized action",rnk
  endif

  if (n_dcm==1) then
    call read_mirror_h5_sl(dom, ntime_tmp)
  elseif (n_dcm>1) then
    if (mod(ntime_tmp,n_dcm)==0) then
      if (rnk==0) write(*,'("--> SEM : read mirror_sl, iteration : ",i6.6)') ntime_tmp
      do i = 1,n_spl+1
        j = int(ntime_tmp/n_dcm)+1
        call read_mirror_h5_sl(dom, j+i-1-1)
        displ_sl(:,:,i) = dom%mirror_sl%fields(1:3,:)
        force_sl(:,:,i) = dom%mirror_sl%fields(4:6,:)
      enddo
    endif
    dom%mirror_sl%fields(1:3,:) = 0.
    dom%mirror_sl%fields(4:6,:) = 0.
    do i = 1,n_spl+1
      j = mod(ntime_tmp,n_dcm)+(n_spl+1-i)*n_dcm+1
      dom%mirror_sl%fields(1:3,:) = dom%mirror_sl%fields(1:3,:)+bspl_tmp(j)*displ_sl(:,:,i)
      dom%mirror_sl%fields(4:6,:) = dom%mirror_sl%fields(4:6,:)+bspl_tmp(j)*force_sl(:,:,i)
    enddo
  else
    t = (ntime_tmp-1)*d_t
    tnt = t/d_tm
    recn = int(tnt)
    tnt = tnt-dble(recn)
    if (ntime_tmp==0.or.recn/=recp) then
      j1 = 1
      j2 = n_spl+1
      if (ntime_tmp>0) then
        if (recn==recp+1) then
          displ_sl(:,:,1:n_spl) = displ_sl(:,:,2:n_spl+1)
          force_sl(:,:,1:n_spl) = force_sl(:,:,2:n_spl+1)
          j1 = n_spl+1
        elseif (recn==recp-1) then
          displ_sl(:,:,2:n_spl+1) = displ_sl(:,:,1:n_spl)
          force_sl(:,:,2:n_spl+1) = force_sl(:,:,1:n_spl)
          j2 = 1
        endif
      endif
      do j = j1,j2
        if (recn+j<1.or.recn+j>n_tm) then
          displ_sl(:,:,j) = 0.d0
          force_sl(:,:,j) = 0.d0
        else
          call read_mirror_h5_sl(dom, recn+j-1)
          displ_sl(:,:,j) = dom%mirror_sl%fields(1:3,:)
          force_sl(:,:,j) = dom%mirror_sl%fields(4:6,:)
        endif
      enddo
      recp = recn
    endif
    dom%mirror_sl%fields(1:3,:) = 0.d0
    dom%mirror_sl%fields(4:6,:) = 0.d0
    do i = 1,n_spl+1
      tmp = bspln(0,n_spl,dble(tnt+i-1))
      dom%mirror_sl%fields(1:3,:) = dom%mirror_sl%fields(1:3,:)+tmp*displ_sl(:,:,n_spl+2-i)
      dom%mirror_sl%fields(4:6,:) = dom%mirror_sl%fields(4:6,:)+tmp*force_sl(:,:,n_spl+2-i)
    enddo
  endif

end subroutine load_mirror_sl

subroutine load_mirror_fl(dom, ntime)
  use sspline
  use champs_fluid
  implicit none
  type(domain_fluid), intent(inout) :: dom
  integer, intent(in) :: ntime
  integer :: i,j,j1,j2,ntime_tmp,recn
  real(fpp) :: t,tnt,tmp

  if (dom%mirror_type==1) then
    ntime_tmp = ntime
  elseif (dom%mirror_type==2) then
    ntime_tmp = n_t-ntime-1
  else
    write(*,*) "STOP, unauthorized action",rnk
  endif

  if (n_dcm==1) then
    call read_mirror_h5_fl(dom, ntime_tmp)
  elseif (n_dcm>1) then
    if (mod(ntime_tmp,n_dcm)==0) then
      if (rnk==0) write(*,'("--> SEM : read mirror_fl, iteration : ",i6.6)') ntime_tmp
      do i = 1,n_spl+1
        j = int(ntime_tmp/n_dcm)+1
        call read_mirror_h5_fl(dom, j+i-1-1)
        displ_fl(:,i) = dom%mirror_fl%fields(1,:)
        force_fl(:,i) = dom%mirror_fl%fields(2,:)
      enddo
    endif
    dom%mirror_fl%fields(1,:) = 0.
    dom%mirror_fl%fields(2,:) = 0.
    do i = 1,n_spl+1
      j = mod(ntime_tmp,n_dcm)+(n_spl+1-i)*n_dcm+1
      dom%mirror_fl%fields(1,:) = dom%mirror_fl%fields(1,:)+bspl_tmp(j)*displ_fl(:,i)
      dom%mirror_fl%fields(2,:) = dom%mirror_fl%fields(2,:)+bspl_tmp(j)*force_fl(:,i)
    enddo
  else
    t = (ntime_tmp-1)*d_t
    tnt = t/d_tm
    recn = int(tnt)
    tnt = tnt-dble(recn)
    if (ntime_tmp==0.or.recn/=recp) then
      j1 = 1
      j2 = n_spl+1
      if (ntime_tmp>0) then
        if (recn==recp+1) then
          displ_fl(:,1:n_spl) = displ_fl(:,2:n_spl+1)
          force_fl(:,1:n_spl) = force_fl(:,2:n_spl+1)
          j1 = n_spl+1
        elseif (recn==recp-1) then
          displ_fl(:,2:n_spl+1) = displ_fl(:,1:n_spl)
          force_fl(:,2:n_spl+1) = force_fl(:,1:n_spl)
          j2 = 1
        endif
      endif
      do j = j1,j2
        if (recn+j<1.or.recn+j>n_tm) then
          displ_fl(:,j) = 0.d0
          force_fl(:,j) = 0.d0
        else
          call read_mirror_h5_fl(dom, recn+j-1)
          displ_fl(:,j) = dom%mirror_fl%fields(1,:)
          force_fl(:,j) = dom%mirror_fl%fields(2,:)
        endif
      enddo
      recp = recn
    endif
    dom%mirror_fl%fields(1,:) = 0.d0
    dom%mirror_fl%fields(2,:) = 0.d0
    do i = 1,n_spl+1
      tmp = bspln(0,n_spl,dble(tnt+i-1))
      dom%mirror_fl%fields(1,:) = dom%mirror_fl%fields(1,:)+tmp*displ_fl(:,n_spl+2-i)
      dom%mirror_fl%fields(2,:) = dom%mirror_fl%fields(2,:)+tmp*force_fl(:,n_spl+2-i)
    enddo
  endif

end subroutine load_mirror_fl

subroutine bchslv_sl(dom, w, nbands, nrow)
  use champs_solid
  implicit none
  type(domain_solid), intent(inout) :: dom
  integer, intent(in) :: nrow,nbands
  real(fpp), intent(in) :: w(nbands,nrow)
  integer :: j,n

  !! Forward substitution, Solve L*Y = B.
  do n = 1,nrow-nbands+1
    if (n==1) then
      do j = 0,nbands-1
        call read_mirror_h5_sl(dom, j+n-1)
        displ_sl(:,:,j+1) = dom%mirror_sl%fields(1:3,:)
        force_sl(:,:,j+1) = dom%mirror_sl%fields(4:6,:)
      enddo
    else
      do j = 0,nbands-2
        displ_sl(:,:,j+1) = displ_sl(:,:,j+2)
        force_sl(:,:,j+1) = force_sl(:,:,j+2)
      enddo
      call read_mirror_h5_sl(dom, n+nbands-2)
      displ_sl(:,:,nbands) = dom%mirror_sl%fields(1:3,:)
      force_sl(:,:,nbands) = dom%mirror_sl%fields(4:6,:)
    endif
    do j = 1,nbands-1
      displ_sl(:,:,j+1) = displ_sl(:,:,j+1)-w(j+1,n)*displ_sl(:,:,1)
      force_sl(:,:,j+1) = force_sl(:,:,j+1)-w(j+1,n)*force_sl(:,:,1)
    enddo
    dom%mirror_sl%fields(1:3,:) = displ_sl(:,:,1)
    dom%mirror_sl%fields(4:6,:) = force_sl(:,:,1)
    call write_mirror_h5_sl(dom, n-1)
  enddo
  do n = nrow-nbands+2,nrow
    do j = 0,nbands-2
      displ_sl(:,:,j+1) = displ_sl(:,:,j+2)
      force_sl(:,:,j+1) = force_sl(:,:,j+2)
    enddo
    do j = 1,nrow-n
      displ_sl(:,:,j+1) = displ_sl(:,:,j+1)-w(j+1,n)*displ_sl(:,:,1)
      force_sl(:,:,j+1) = force_sl(:,:,j+1)-w(j+1,n)*force_sl(:,:,1)
    enddo
    dom%mirror_sl%fields(1:3,:) = displ_sl(:,:,1)
    dom%mirror_sl%fields(4:6,:) = force_sl(:,:,1)
    call write_mirror_h5_sl(dom, n-1)
  enddo
  !! Back substitution, Solve L'*X = D**(-1)*Y.
  do n = nrow,nrow-nbands+2,-1
    do j = nrow-n,1,-1
      displ_sl(:,:,j+1) = displ_sl(:,:,j)
      force_sl(:,:,j+1) = force_sl(:,:,j)
    enddo
    call read_mirror_h5_sl(dom, n-1)
    displ_sl(:,:,1) = dom%mirror_sl%fields(1:3,:)*w(1,n)
    force_sl(:,:,1) = dom%mirror_sl%fields(4:6,:)*w(1,n)
    do j = 1,nrow-n
      displ_sl(:,:,1) = displ_sl(:,:,1)-w(j+1,n)*displ_sl(:,:,j+1)
      force_sl(:,:,1) = force_sl(:,:,1)-w(j+1,n)*force_sl(:,:,j+1)
    enddo
    dom%mirror_sl%fields(1:3,:) = displ_sl(:,:,1)
    dom%mirror_sl%fields(4:6,:) = force_sl(:,:,1)
    call write_mirror_h5_sl(dom, n-1)
  enddo
  do n = nrow-nbands+1,1,-1
    do j = nbands-1,1,-1
      displ_sl(:,:,j+1) = displ_sl(:,:,j)
      force_sl(:,:,j+1) = force_sl(:,:,j)
    enddo
    call read_mirror_h5_sl(dom, n-1)
    displ_sl(:,:,1) = dom%mirror_sl%fields(1:3,:)*w(1,n)
    force_sl(:,:,1) = dom%mirror_sl%fields(4:6,:)*w(1,n)
    do j = 1,nbands-1
      displ_sl(:,:,1) = displ_sl(:,:,1)-w(j+1,n)*displ_sl(:,:,j+1)
      force_sl(:,:,1) = force_sl(:,:,1)-w(j+1,n)*force_sl(:,:,j+1)
    enddo
    dom%mirror_sl%fields(1:3,:) = displ_sl(:,:,1)
    dom%mirror_sl%fields(4:6,:) = force_sl(:,:,1)
    call write_mirror_h5_sl(dom, n-1)
  enddo
  return

end subroutine bchslv_sl

subroutine bchslv_fl(dom, w, nbands, nrow)
  use champs_fluid
  implicit none
  type(domain_fluid), intent(inout) :: dom
  integer, intent(in) :: nrow,nbands
  real(fpp), intent(in) :: w(nbands,nrow)
  integer :: j,n

  !! Forward substitution, Solve L*Y = B.
  do n = 1,nrow-nbands+1
    if (n==1) then
      do j = 0,nbands-1
        call read_mirror_h5_fl(dom, j+n-1)
        displ_fl(:,j+1) = dom%mirror_fl%fields(1,:)
        force_fl(:,j+1) = dom%mirror_fl%fields(2,:)
      enddo
    else
      do j = 0,nbands-2
        displ_fl(:,j+1) = displ_fl(:,j+2)
        force_fl(:,j+1) = force_fl(:,j+2)
      enddo
      call read_mirror_h5_fl(dom, n+nbands-2)
      displ_fl(:,nbands) = dom%mirror_fl%fields(1,:)
      force_fl(:,nbands) = dom%mirror_fl%fields(2,:)
    endif
    do j = 1,nbands-1
      displ_fl(:,j+1) = displ_fl(:,j+1)-w(j+1,n)*displ_fl(:,1)
      force_fl(:,j+1) = force_fl(:,j+1)-w(j+1,n)*force_fl(:,1)
    enddo
    dom%mirror_fl%fields(1,:) = displ_fl(:,1)
    dom%mirror_fl%fields(2,:) = force_fl(:,1)
    call write_mirror_h5_fl(dom, n-1)
  enddo
  do n = nrow-nbands+2,nrow
    do j = 0,nbands-2
      displ_fl(:,j+1) = displ_fl(:,j+2)
      force_fl(:,j+1) = force_fl(:,j+2)
    enddo
    do j = 1,nrow-n
      displ_fl(:,j+1) = displ_fl(:,j+1)-w(j+1,n)*displ_fl(:,1)
      force_fl(:,j+1) = force_fl(:,j+1)-w(j+1,n)*force_fl(:,1)
    enddo
    dom%mirror_fl%fields(1,:) = displ_fl(:,1)
    dom%mirror_fl%fields(2,:) = force_fl(:,1)
    call write_mirror_h5_fl(dom, n-1)
  enddo
  !! Back substitution, Solve L'*X = D**(-1)*Y.
  do n = nrow,nrow-nbands+2,-1
    do j = nrow-n,1,-1
      displ_fl(:,j+1) = displ_fl(:,j)
      force_fl(:,j+1) = force_fl(:,j)
    enddo
    call read_mirror_h5_fl(dom, n-1)
    displ_fl(:,1) = dom%mirror_fl%fields(1,:)*w(1,n)
    force_fl(:,1) = dom%mirror_fl%fields(2,:)*w(1,n)
    do j = 1,nrow-n
      displ_fl(:,1) = displ_fl(:,1)-w(j+1,n)*displ_fl(:,j+1)
      force_fl(:,1) = force_fl(:,1)-w(j+1,n)*force_fl(:,j+1)
    enddo
    dom%mirror_fl%fields(1,:) = displ_fl(:,1)
    dom%mirror_fl%fields(2,:) = force_fl(:,1)
    call write_mirror_h5_fl(dom, n-1)
  enddo
  do n = nrow-nbands+1,1,-1
    do j = nbands-1,1,-1
      displ_fl(:,j+1) = displ_fl(:,j)
      force_fl(:,j+1) = force_fl(:,j)
    enddo
    call read_mirror_h5_fl(dom, n-1)
    displ_fl(:,1) = dom%mirror_fl%fields(1,:)*w(1,n)
    force_fl(:,1) = dom%mirror_fl%fields(2,:)*w(1,n)
    do j = 1,nbands-1
      displ_fl(:,1) = displ_fl(:,1)-w(j+1,n)*displ_fl(:,j+1)
      force_fl(:,1) = force_fl(:,1)-w(j+1,n)*force_fl(:,j+1)
    enddo
    dom%mirror_fl%fields(1,:) = displ_fl(:,1)
    dom%mirror_fl%fields(2,:) = force_fl(:,1)
    call write_mirror_h5_fl(dom, n-1)
  enddo
  return

end subroutine bchslv_fl

subroutine append_mirror_h5_sl(dom)
  use champs_solid
  use semdatafiles
  use sem_hdf5
  implicit none
  type(domain_solid), intent(inout) :: dom
  character (len=40) :: dname
  character (len=MAX_FILE_SIZE) :: fnamef
  integer(HID_T) :: fid, dset_id
  integer :: hdferr

  call semname_mirrorfile_h5(rnk, fnamef)
  call h5fopen_f(fnamef, H5F_ACC_RDWR_F, fid, hdferr)
  dname = "Fields_sl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call append_dataset_2d(dset_id, dom%mirror_sl%fields, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine append_mirror_h5_sl

subroutine append_mirror_h5_fl(dom)
  use champs_fluid
  use semdatafiles
  use sem_hdf5
  implicit none
  type(domain_fluid), intent(inout) :: dom
  character (len=40) :: dname
  character (len=MAX_FILE_SIZE) :: fnamef
  integer(HID_T) :: fid, dset_id
  integer :: hdferr

  call semname_mirrorfile_h5(rnk, fnamef)
  call h5fopen_f(fnamef, H5F_ACC_RDWR_F, fid, hdferr)
  dname = "Fields_fl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call append_dataset_2d(dset_id, dom%mirror_fl%fields, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine append_mirror_h5_fl

subroutine read_mirror_h5_sl(dom, n)
  use champs_solid
  use semdatafiles
  use sem_hdf5
  implicit none
  type(domain_solid), intent(inout) :: dom
  integer, intent(in) :: n
  character (len=40) :: dname
  character (len=MAX_FILE_SIZE) :: fnamef
  integer(HID_T) :: fid, dset_id
  integer :: hdferr
  integer(HSIZE_T), dimension(2) :: offset, dcount

  offset = (/0,n*dom%mirror_sl%n_glltot/)
  dcount = (/6,dom%mirror_sl%n_glltot/)
  call semname_mirrorfile_h5(rnk, fnamef)
  call h5fopen_f(fnamef, H5F_ACC_RDONLY_F, fid, hdferr)
  dname = "Fields_sl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call read_datasubset_2d(dset_id, offset, dcount, dom%mirror_sl%fields, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine read_mirror_h5_sl

subroutine read_mirror_h5_fl(dom, n)
  use champs_fluid
  use semdatafiles
  use sem_hdf5
  implicit none
  type(domain_fluid), intent(inout) :: dom
  integer, intent(in) :: n
  character (len=40) :: dname
  character (len=MAX_FILE_SIZE) :: fnamef
  integer(HID_T) :: fid, dset_id
  integer :: hdferr
  integer(HSIZE_T), dimension(2) :: offset, dcount

  offset = (/0,n*dom%mirror_fl%n_glltot/)
  dcount = (/2,dom%mirror_fl%n_glltot/)
  call semname_mirrorfile_h5(rnk, fnamef)
  call h5fopen_f(fnamef, H5F_ACC_RDONLY_F, fid, hdferr)
  dname = "Fields_fl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call read_datasubset_2d(dset_id, offset, dcount, dom%mirror_fl%fields, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine read_mirror_h5_fl

subroutine write_mirror_h5_sl(dom, n)
  use champs_solid
  use semdatafiles
  use sem_hdf5
  implicit none
  type(domain_solid), intent(inout) :: dom
  integer, intent(in) :: n
  character (len=40) :: dname
  character (len=MAX_FILE_SIZE) :: fnamef
  integer(HID_T) :: fid, dset_id
  integer :: hdferr
  integer(HSIZE_T), dimension(2) :: offset, dcount

  offset = (/0,n*dom%mirror_sl%n_glltot/)
  dcount = (/6,dom%mirror_sl%n_glltot/)
  call semname_mirrorfile_h5(rnk, fnamef)
  call h5fopen_f(fnamef, H5F_ACC_RDWR_F, fid, hdferr)
  dname = "Fields_sl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call write_datasubset_2d(dset_id, offset, dcount, dom%mirror_sl%fields, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine write_mirror_h5_sl

subroutine write_mirror_h5_fl(dom, n)
  use champs_fluid
  use semdatafiles
  use sem_hdf5
  implicit none
  type(domain_fluid), intent(inout) :: dom
  integer, intent(in) :: n
  character (len=40) :: dname
  character (len=MAX_FILE_SIZE) :: fnamef
  integer(HID_T) :: fid, dset_id
  integer :: hdferr
  integer(HSIZE_T), dimension(2) :: offset, dcount

  offset = (/0,n*dom%mirror_fl%n_glltot/)
  dcount = (/2,dom%mirror_fl%n_glltot/)
  call semname_mirrorfile_h5(rnk, fnamef)
  call h5fopen_f(fnamef, H5F_ACC_RDWR_F, fid, hdferr)
  dname = "Fields_fl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call write_datasubset_2d(dset_id, offset, dcount, dom%mirror_fl%fields, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine write_mirror_h5_fl

subroutine write_mirror_h5_attr()
  use semdatafiles
  use sem_hdf5
  implicit none
  character(len=40) :: aname
  character(len=MAX_FILE_SIZE) :: fnamef
  integer(HID_T) :: fid,attr_id,memspace,space_id
  integer :: hdferr
  integer(HSIZE_T), dimension(1) :: dims

  dims(1) = 1
  call semname_mirrorfile_h5(rnk,fnamef)
  call h5fopen_f(fnamef,H5F_ACC_RDWR_F,fid,hdferr)
  aname = "Delta_t"
  call h5screate_f(H5S_SCALAR_F,memspace,hdferr)
  call h5acreate_f(fid,aname,H5T_NATIVE_DOUBLE,memspace,attr_id,hdferr)
  call h5awrite_f(attr_id,H5T_NATIVE_DOUBLE,d_tm,dims,hdferr)
  call h5aclose_f(attr_id,hdferr)
  call h5sclose_f(memspace,hdferr)

  aname = "Num_t"
  call h5screate_f(H5S_SCALAR_F,memspace,hdferr)
  call h5acreate_f(fid,aname,H5T_NATIVE_INTEGER,memspace,attr_id,hdferr)
  call h5awrite_f(attr_id,H5T_NATIVE_INTEGER,n_tm,dims,hdferr)
  call h5aclose_f(attr_id,hdferr)
  call h5sclose_f(memspace,hdferr)

  call h5fclose_f(fid,hdferr)

end subroutine write_mirror_h5_attr

subroutine read_mirror_h5_attr()
  use semdatafiles
  use sem_hdf5
  implicit none
  character(len=40) :: aname
  character(len=MAX_FILE_SIZE) :: fnamef
  integer(HID_T) :: fid,attr_id,memspace,space_id
  integer :: hdferr
  integer(HSIZE_T), dimension(1) :: dims,max_dims
  real(fpp), dimension(1) :: val_r
  integer, dimension(1) :: val_i

  dims(1) = 1
  call semname_mirrorfile_h5(rnk,fnamef)
  call h5fopen_f(fnamef,H5F_ACC_RDONLY_F,fid,hdferr)

  aname = "Delta_t"
  call h5aopen_f(fid,aname,attr_id,hdferr)
  call h5aget_space_f(attr_id,space_id,hdferr)
  call h5sget_simple_extent_dims_f(space_id,dims,max_dims,hdferr)
  call h5aread_f(attr_id,H5T_NATIVE_DOUBLE,val_r,dims,hdferr)
  d_tm = val_r(1)
  call h5aclose_f(attr_id,hdferr)
  call h5sclose_f(space_id,hdferr)

  aname = "Num_t"
  call h5aopen_f(fid,aname,attr_id,hdferr)
  call h5aget_space_f(attr_id,space_id,hdferr)
  call h5sget_simple_extent_dims_f(space_id,dims,max_dims,hdferr)
  call h5aread_f(attr_id,H5T_NATIVE_INTEGER,val_i,dims,hdferr)
  n_tm = val_i(1)
  call h5aclose_f(attr_id,hdferr)
  call h5sclose_f(space_id,hdferr)

  call h5fclose_f(fid,hdferr)

end subroutine read_mirror_h5_attr

end module smirror

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
