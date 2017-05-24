!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
module smirror
  use constants, only : fpp
  implicit none
  public :: init_mirror, create_mirror_files, dump_mirror_solid, load_mirror_solid
  private

  integer :: n_elm
  integer, dimension(:), allocatable :: map2elem
  real(fpp), dimension(:), allocatable :: wfunc, bspl_tmp
  real(fpp), dimension(:,:,:), allocatable :: displ_tmp, force_tmp
  integer :: rnk, n_spl, n_dcm, n_t, n_tdwn

contains

subroutine init_mirror(Tdomain, surf)
  use sspline
  use sdomain
  implicit none
  type(domain), intent(inout) :: Tdomain
  type(SurfaceT), intent(inout) :: surf
  integer :: n_elmtot, n_glltot, n_gll1, n_gll3

  !! some local variables
  n_elmtot = Tdomain%sdom%nblocks
  n_gll1 = Tdomain%sdom%ngll
  n_gll3 = n_gll1*n_gll1*n_gll1
  if (surf%surf_sl%nbtot/=0) then
    !! find mirror elements
    allocate(map2elem(0:n_elmtot-1))
    call count_mirror_elem(Tdomain, surf)
    n_glltot = n_elm*n_gll3
    if (n_elm>0) then
      write(*,'("--> SEM : mirror elements : ",i3,i6)') Tdomain%rank,n_elm
      !! define window function
      allocate(wfunc(0:n_elm*n_gll3-1))
      call mirror_window_function(Tdomain, surf)
      !! define decimation properties
      n_spl = Tdomain%config%mirror_nspl
      call find_decim_factor(Tdomain)
      allocate(bspl_tmp((n_spl+1)*n_dcm+1))
      call bmn(bspl_tmp, n_dcm, n_spl)
      !! allocate mirror arrays
      allocate(Tdomain%sdom%mirror%map(0:n_elmtot-1))
      allocate(Tdomain%sdom%mirror%winf(0:n_glltot-1))
      allocate(Tdomain%sdom%mirror%displ(0:2,0:n_glltot-1))
      allocate(Tdomain%sdom%mirror%force(0:2,0:n_glltot-1))
      allocate(displ_tmp(0:2,0:n_glltot-1,n_spl+1))
      allocate(force_tmp(0:2,0:n_glltot-1,n_spl+1))
      !! define mirror properties
      Tdomain%sdom%mirror%n_elem = n_elm
      Tdomain%sdom%mirror%n_glltot = n_glltot
      Tdomain%sdom%mirror%n_gll = n_gll1
      Tdomain%sdom%mirror%map = map2elem
      Tdomain%sdom%mirror%winf = wfunc
      displ_tmp = 0.
      force_tmp = 0.
      deallocate(wfunc)
    else
      Tdomain%sdom%mirror%n_glltot = 0
    endif
    deallocate(map2elem)
  endif

end subroutine init_mirror

subroutine count_mirror_elem(Tdomain, surf)
  use sdomain
  implicit none
  type(domain), intent(in) :: Tdomain
  type(SurfaceT), intent(in) :: surf
  integer :: i, ff, ee
  logical, dimension(:), allocatable :: ltmp

  allocate(ltmp(0:Tdomain%sdom%nblocks-1))
  ltmp = .false.
  do i = 0,surf%surf_sl%n_faces-1
    ff = surf%surf_sl%if_faces(i)
    ee = Tdomain%sFace(ff)%elem_0
    if (surf%surf_sl%if_norm(i)<0) ee = Tdomain%sFace(ff)%elem_1
    if (ee>=0) ltmp(Tdomain%specel(ee)%lnum) = .true.
  enddo
  map2elem = -1
  n_elm = 0
  do i = 0,Tdomain%sdom%nblocks-1
    if (ltmp(i)) then
      map2elem(i) = n_elm
      n_elm = n_elm+1
    endif
  enddo
  deallocate(ltmp)

end subroutine count_mirror_elem

subroutine mirror_window_function(Tdomain, surf)
  use sdomain
  implicit none
  type(domain), intent(in) :: Tdomain
  type(SurfaceT), intent(in) :: surf
  integer :: ii, jj, i, j, k, ff, ee, nn
  integer :: fnorm, lnum, ngll, ind
  real(fpp) :: winx, winy, winz
  real(fpp), dimension(0:2,0:4) :: fnodes

  ngll = Tdomain%sdom%ngll
  wfunc = 1.
  do ii = 0,surf%surf_sl%n_faces-1
    ff = surf%surf_sl%if_faces(ii)
    ee = Tdomain%sFace(ff)%elem_0
    nn = surf%surf_sl%if_norm(ii)
    if (nn<0) ee = Tdomain%sFace(ff)%elem_1
    if (ee==-1) cycle
    lnum = Tdomain%specel(ee)%lnum
    do jj = 0,3
      fnodes(0:2,jj) = Tdomain%Coord_Nodes(0:2,Tdomain%sFace(ff)%inodes(jj))
    enddo
    call mirror_face_normal(fnodes, fnorm) ! find direction for window function
    winx = 1.
    winy = 1.
    winz = 1.
    do k = 0,ngll-1
      if (fnorm==2) winz = (0.5*cos(4.*atan(1.)*k*1./(ngll-1))+0.5)
      do j = 0,ngll-1
        if (fnorm==1) winy = (0.5*cos(4.*atan(1.)*j*1./(ngll-1))+0.5)
        do i = 0,ngll-1
          if (fnorm==0) winx = (0.5*cos(4.*atan(1.)*i*1./(ngll-1))+0.5)
          ind = (map2elem(lnum)*ngll*ngll*ngll)+i+(j*ngll)+(k*ngll*ngll)
          if (nn>0) then
            wfunc(ind) = wfunc(ind)*(1.-winx*winy*winz)
          else
            wfunc(ind) = wfunc(ind)*winx*winy*winz
          endif
        enddo
      enddo
    enddo
  enddo

end subroutine mirror_window_function

subroutine mirror_face_normal(fnodes, fdir)
  use shape_geom_3d
  implicit none
  real(fpp), dimension(0:2, 0:3), intent(in) :: fnodes
  real(fpp), dimension(0:2) :: dfx, dfy, fnorm
  integer, intent(out) :: fdir
  integer :: i

  do i = 0,2
    dfx(i) = -fnodes(i,0)+fnodes(i,1)+fnodes(i,2)-fnodes(i,3)
    dfy(i) = -fnodes(i,0)-fnodes(i,1)+fnodes(i,2)+fnodes(i,3)
  enddo
  call cross_prod(dfx, dfy, fnorm)
  fdir = 0
  if (abs(fnorm(1))>abs(fnorm(fdir))) fdir = 1
  if (abs(fnorm(2))>abs(fnorm(fdir))) fdir = 2

end subroutine mirror_face_normal

subroutine find_decim_factor(Tdomain)
  use sdomain
  implicit none
  type(domain), intent(inout) :: Tdomain
  integer :: i
  real(fpp) :: fmax,d_t

  fmax = Tdomain%config%mirror_fmax
  d_t = Tdomain%TimeD%dtmin
  n_t = Tdomain%TimeD%ntimeMax

  if (fmax>0.) then
    n_dcm = int(0.25/(fmax*d_t))
    n_tdwn = int(n_t/n_dcm)+(n_spl+1)
    if (Tdomain%rank==25) write(*,'("--> SEM : mirror decimation : ",2i6)') n_dcm,n_tdwn
  else
    n_dcm = 1
  endif

end subroutine find_decim_factor

subroutine create_mirror_files(Tdomain)
  use sdomain
  implicit none
  type(domain), intent(inout) :: Tdomain

  Tdomain%sdom%mirror%rank = Tdomain%rank
  rnk = Tdomain%rank
  ! Create mirror files for i_type_mirror=record and !recovery_mode
  if (Tdomain%TimeD%NtimeMin==0 .and. Tdomain%sdom%mirror_type==0) then
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
  !! Displacement field in solid
  dname = "Displ_sl"
  call create_dset_2d(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, &
    int(3,HSIZE_T), int(H5S_UNLIMITED_F,HSIZE_T), dset_id)
  call h5dclose_f(dset_id, hdferr)
  !! Displacement field in fluid
  dname = "Displ_fl"
  call create_dset(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, &
    int(H5S_UNLIMITED_F,HSIZE_T), dset_id)
  call h5dclose_f(dset_id, hdferr)
  !! Force field in solid
  dname = "Force_sl"
  call create_dset_2d(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, &
    int(3,HSIZE_T), int(H5S_UNLIMITED_F,HSIZE_T), dset_id)
  call h5dclose_f(dset_id, hdferr)
  !! Force field in fluid
  dname = "Force_fl"
  call create_dset(fid, trim(adjustl(dname)), H5T_IEEE_F64LE, &
    int(H5S_UNLIMITED_F,HSIZE_T), dset_id)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine create_mirror_h5

subroutine dump_mirror_solid(dom, ntime)
  use champs_solid
  use sspline
  implicit none
  type(domain_solid), intent(inout) :: dom
  integer, intent(in) :: ntime
  integer :: i,j
  real(fpp), dimension(:), allocatable :: ddiag
  real(fpp), dimension(:,:), allocatable :: diag

  if (n_dcm==1) then
    call append_mirror_h5_solid(dom)
  else
    do i = 1,n_spl+1
      j = mod(ntime,n_dcm)+(n_spl+1-i)*n_dcm+1
      displ_tmp(:,:,i) = displ_tmp(:,:,i)+bspl_tmp(j)*dom%mirror%displ(:,:)
      force_tmp(:,:,i) = force_tmp(:,:,i)+bspl_tmp(j)*dom%mirror%force(:,:)
    enddo
    if (mod(ntime,n_dcm)==n_dcm-1) then
      if (rnk==0) write(*,'("--> SEM : dump mirror, iteration : ",i6.6)') ntime
      dom%mirror%displ(:,:) = displ_tmp(:,:,1)
      dom%mirror%force(:,:) = force_tmp(:,:,1)
      call append_mirror_h5_solid(dom)
      do i = 1,n_spl
        displ_tmp(:,:,i) = displ_tmp(:,:,i+1)
        force_tmp(:,:,i) = force_tmp(:,:,i+1)
      enddo
      displ_tmp(:,:,n_spl+1) = 0.
      force_tmp(:,:,n_spl+1) = 0.
    endif
    if (ntime==n_t-1) then
      if (rnk==0) write(*,'("--> SEM : dump mirror, iteration : ",i6.6)') ntime
      do i = 1,n_spl+1
        dom%mirror%displ(:,:) = displ_tmp(:,:,i)
        dom%mirror%force(:,:) = force_tmp(:,:,i)
        call append_mirror_h5_solid(dom)
      enddo
      allocate(diag(n_spl+1,n_tdwn),ddiag(n_tdwn))
      call bmdiag(n_dcm,n_spl,n_tdwn,n_t,diag)
      call bchfac(diag,n_spl+1,n_tdwn,ddiag)
      call bchslv(dom,diag,n_spl+1,n_tdwn)
      deallocate(diag,ddiag)
    endif
  endif

end subroutine dump_mirror_solid

subroutine load_mirror_solid(dom, ntime)
  use champs_solid
  implicit none
  type(domain_solid), intent(inout) :: dom
  integer, intent(in) :: ntime
  integer :: i,j,ntime_tmp

  if (dom%mirror_type==1) then
    ntime_tmp = ntime
  elseif (dom%mirror_type==2) then
    ntime_tmp = n_t-ntime-1
  else
    write(*,*) "STOP, unauthorized action",rnk
  endif

  if (n_dcm==1) then
    call read_mirror_h5_solid(dom, ntime_tmp)
  else
    if (mod(ntime_tmp,n_dcm)==0) then
      if (rnk==0) write(*,'("--> SEM : read mirror, iteration : ",i6.6)') ntime_tmp
      do i = 1,n_spl+1
        j = int(ntime_tmp/n_dcm)+1
        call read_mirror_h5_solid(dom, j+i-1-1)
        displ_tmp(:,:,i) = dom%mirror%displ(:,:)
        force_tmp(:,:,i) = dom%mirror%force(:,:)
      enddo
    endif    
    dom%mirror%displ(:,:) = 0.
    dom%mirror%force(:,:) = 0.
    do i = 1,n_spl+1
      j = mod(ntime_tmp,n_dcm)+(n_spl+1-i)*n_dcm+1
      dom%mirror%displ(:,:) = dom%mirror%displ(:,:)+bspl_tmp(j)*displ_tmp(:,:,i)
      dom%mirror%force(:,:) = dom%mirror%force(:,:)+bspl_tmp(j)*force_tmp(:,:,i)
    enddo
  endif

end subroutine load_mirror_solid

subroutine bchslv(dom, w, nbands, nrow)
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
        call read_mirror_h5_solid(dom, j+n-1)
        displ_tmp(:,:,j+1) = dom%mirror%displ(:,:)
        force_tmp(:,:,j+1) = dom%mirror%force(:,:)
      enddo
    else
      do j = 0,nbands-2
        displ_tmp(:,:,j+1) = displ_tmp(:,:,j+2)
        force_tmp(:,:,j+1) = force_tmp(:,:,j+2)
      enddo
      call read_mirror_h5_solid(dom, n+nbands-2)
      displ_tmp(:,:,nbands) = dom%mirror%displ(:,:)
      force_tmp(:,:,nbands) = dom%mirror%force(:,:)
    endif
    do j = 1,nbands-1
      displ_tmp(:,:,j+1) = displ_tmp(:,:,j+1)-w(j+1,n)*displ_tmp(:,:,1)
      force_tmp(:,:,j+1) = force_tmp(:,:,j+1)-w(j+1,n)*force_tmp(:,:,1)
    enddo
    dom%mirror%displ(:,:) = displ_tmp(:,:,1)
    dom%mirror%force(:,:) = force_tmp(:,:,1)
    call write_mirror_h5_solid(dom, n-1)
  enddo
  do n = nrow-nbands+2,nrow
    do j = 0,nbands-2
      displ_tmp(:,:,j+1) = displ_tmp(:,:,j+2)
      force_tmp(:,:,j+1) = force_tmp(:,:,j+2)
    enddo
    do j = 1,nrow-n
      displ_tmp(:,:,j+1) = displ_tmp(:,:,j+1)-w(j+1,n)*displ_tmp(:,:,1)
      force_tmp(:,:,j+1) = force_tmp(:,:,j+1)-w(j+1,n)*force_tmp(:,:,1)
    enddo
    dom%mirror%displ(:,:) = displ_tmp(:,:,1)
    dom%mirror%force(:,:) = force_tmp(:,:,1)
    call write_mirror_h5_solid(dom, n-1)
  enddo
  !! Back substitution, Solve L'*X = D**(-1)*Y.
  do n = nrow,nrow-nbands+2,-1
    do j = nrow-n,1,-1
      displ_tmp(:,:,j+1) = displ_tmp(:,:,j)
      force_tmp(:,:,j+1) = force_tmp(:,:,j)
    enddo
    call read_mirror_h5_solid(dom, n-1)
    displ_tmp(:,:,1) = dom%mirror%displ(:,:)*w(1,n)
    force_tmp(:,:,1) = dom%mirror%force(:,:)*w(1,n)
    do j = 1,nrow-n
      displ_tmp(:,:,1) = displ_tmp(:,:,1)-w(j+1,n)*displ_tmp(:,:,j+1)
      force_tmp(:,:,1) = force_tmp(:,:,1)-w(j+1,n)*force_tmp(:,:,j+1)
    enddo
    dom%mirror%displ(:,:) = displ_tmp(:,:,1)
    dom%mirror%force(:,:) = force_tmp(:,:,1)
    call write_mirror_h5_solid(dom, n-1)
  enddo
  do n = nrow-nbands+1,1,-1
    do j = nbands-1,1,-1
      displ_tmp(:,:,j+1) = displ_tmp(:,:,j)
      force_tmp(:,:,j+1) = force_tmp(:,:,j)
    enddo
    call read_mirror_h5_solid(dom, n-1)
    displ_tmp(:,:,1) = dom%mirror%displ(:,:)*w(1,n)
    force_tmp(:,:,1) = dom%mirror%force(:,:)*w(1,n)
    do j = 1,nbands-1
      displ_tmp(:,:,1) = displ_tmp(:,:,1)-w(j+1,n)*displ_tmp(:,:,j+1)
      force_tmp(:,:,1) = force_tmp(:,:,1)-w(j+1,n)*force_tmp(:,:,j+1)
    enddo
    dom%mirror%displ(:,:) = displ_tmp(:,:,1)
    dom%mirror%force(:,:) = force_tmp(:,:,1)
    call write_mirror_h5_solid(dom, n-1)
  enddo
  return

end subroutine bchslv

subroutine append_mirror_h5_solid(dom)
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
  dname = "Displ_sl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call append_dataset_2d(dset_id, dom%mirror%displ, hdferr)
  call h5dclose_f(dset_id, hdferr)
  dname = "Force_sl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call append_dataset_2d(dset_id, dom%mirror%force, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine append_mirror_h5_solid

subroutine read_mirror_h5_solid(dom, n)
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

  offset = (/0,n*dom%mirror%n_glltot/)
  dcount = (/3,dom%mirror%n_glltot/)
  call semname_mirrorfile_h5(rnk, fnamef)
  call h5fopen_f(fnamef, H5F_ACC_RDONLY_F, fid, hdferr)
  dname = "Displ_sl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call read_datasubset_2d(dset_id, offset, dcount, dom%mirror%displ, hdferr)
  call h5dclose_f(dset_id, hdferr)
  dname = "Force_sl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call read_datasubset_2d(dset_id, offset, dcount, dom%mirror%force, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine read_mirror_h5_solid

subroutine write_mirror_h5_solid(dom, n)
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

  offset = (/0,n*dom%mirror%n_glltot/)
  dcount = (/3,dom%mirror%n_glltot/)
  call semname_mirrorfile_h5(rnk, fnamef)
  call h5fopen_f(fnamef, H5F_ACC_RDWR_F, fid, hdferr)
  dname = "Displ_sl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call write_datasubset_2d(dset_id, offset, dcount, dom%mirror%displ, hdferr)
  call h5dclose_f(dset_id, hdferr)
  dname = "Force_sl"
  call h5dopen_f(fid, trim(dname), dset_id, hdferr)
  call write_datasubset_2d(dset_id, offset, dcount, dom%mirror%force, hdferr)
  call h5dclose_f(dset_id, hdferr)
  call h5fclose_f(fid, hdferr)

end subroutine write_mirror_h5_solid

end module smirror

!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
