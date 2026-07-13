!! This file is part of SEM
!!
!>
!!\file build_prop_files_2d.F90
!!\brief Read the 2D homogenised Cstar.h5 (written by homofft dump_effectiveC2d_h5)
!!       and interpolate the material tensor onto the GLL points of each
!!       anisotropic-from-file element.
!!
!! HDF5 layout (per property): one group named after the property, with a 2D
!! dataset "samples" (NNx,NNy on a uniform grid) and attributes xMinGlob/xMaxGlob
!! (length 2). The homogenisation already did the GLL->uniform interpolation when
!! writing, so here we read the uniform grid and bilinearly interpolate it back to
!! the SEM2D GLL points (mirror of SEM3D build_prop_files::interpolate_elem_field).
!<
! Just to git
module build_prop_files_2d

    use constants
    use selement
    use sdomain
    use sem_hdf5

    implicit none

    type PropertyField2d
        character(len=100) :: propName = ""
        real(fpp), dimension(0:1) :: MinBound, MaxBound, step
        integer, dimension(0:1)   :: NN
        real(fpp), dimension(:,:), allocatable :: var
    end type PropertyField2d

contains

    !-----------------------------------------------------------------------
    !> Read one property group (pf%propName) of a 2D Cstar.h5: full uniform
    !! grid into pf%var, plus xMinGlob/xMaxGlob and the derived step.
    subroutine read_prop_field_2d(filepath, pf)
        use HDF5
        character(len=*), intent(in)        :: filepath
        type(PropertyField2d), intent(inout):: pf
        integer(HID_T) :: file_id, grp_id, dset_id
        integer(HSIZE_T), dimension(:), allocatable :: dims
        integer(HSIZE_T), dimension(2) :: rdims
        integer :: hdferr, k
        logical :: ok
        ! homofft always writes H5T_NATIVE_DOUBLE; read into a double buffer and
        ! convert, so this is correct whether SEM2D is built single or double (fpp).
        real(kind(0.d0)), dimension(:,:), allocatable :: tmp

        call h5open_f(hdferr)   ! safe if the library is already open (ref-counted)
        call h5fopen_f(trim(filepath), H5F_ACC_RDONLY_F, file_id, hdferr)
        if (hdferr /= 0) then
            write(*,*) "build_prop_files_2d: cannot open ", trim(filepath)
            stop 1
        end if
        call h5lexists_f(file_id, trim(pf%propName), ok, hdferr)
        if (.not. ok) then
            write(*,*) "build_prop_files_2d: missing group '", trim(pf%propName), &
                       "' in ", trim(filepath)
            stop 1
        end if
        call h5gopen_f(file_id, trim(pf%propName), grp_id, hdferr)

        call read_attr_real_vec(grp_id, "xMinGlob", pf%MinBound)
        call read_attr_real_vec(grp_id, "xMaxGlob", pf%MaxBound)
        call read_dims(grp_id, "samples", dims)
        pf%NN(0) = int(dims(1))
        pf%NN(1) = int(dims(2))

        if (allocated(pf%var)) deallocate(pf%var)
        allocate(pf%var(0:pf%NN(0)-1, 0:pf%NN(1)-1))
        allocate(tmp(0:pf%NN(0)-1, 0:pf%NN(1)-1))
        rdims(1) = pf%NN(0)
        rdims(2) = pf%NN(1)
        call h5dopen_f(grp_id, "samples", dset_id, hdferr)
        call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp, rdims, hdferr)
        if (hdferr /= 0) stop "build_prop_files_2d: samples read failed"
        pf%var = real(tmp, fpp)
        deallocate(tmp)
        call h5dclose_f(dset_id, hdferr)
        call h5gclose_f(grp_id, hdferr)
        call h5fclose_f(file_id, hdferr)

        do k = 0, 1
            pf%step(k) = (pf%MaxBound(k) - pf%MinBound(k)) / real(pf%NN(k)-1, fpp)
        end do
    end subroutine read_prop_field_2d

    !-----------------------------------------------------------------------
    !> Bilinear interpolation of a uniform-grid property pf onto the GLL points
    !! of element specel (mirror of SEM3D interpolate_elem_field, 2D).
    subroutine interpolate_elem_field_2d(Tdomain, specel, mat, pf, field)
        type(domain), intent(in)            :: Tdomain
        type(element), intent(in)           :: specel
        type(Subdomain), intent(in)         :: mat
        type(PropertyField2d), intent(in)   :: pf
        real(fpp), dimension(0:specel%ngllx-1,0:specel%ngllz-1), intent(out) :: field
        integer :: i, j, n, idef
        integer, dimension(0:1) :: ii
        real(fpp), dimension(0:1) :: xx, aa
        real(fpp) :: pos, val
        logical :: pml

        pml = is_pml_mat(mat)
        do j = 0, specel%ngllz-1
            do i = 0, specel%ngllx-1
                idef  = specel%Iglobnum(i,j)
                xx(0) = Tdomain%GlobCoord(0, idef)
                xx(1) = Tdomain%GlobCoord(1, idef)
                do n = 0, 1
                    ! Traitement PML : on echantillonne au bord de la PML uniquement.
                    ! Clamping the coordinate to the face freezes the field along the PML
                    ! depth axis (it still varies transversally), so each PML GLL inherits
                    ! the value of the interior GLL it was extruded from.
                    if (pml) then
                        if (mat%pml_width(n) > 0._fpp) then
                            if (xx(n) > mat%pml_pos(n)) xx(n) = mat%pml_pos(n)
                        end if
                        if (mat%pml_width(n) < 0._fpp) then
                            if (xx(n) < mat%pml_pos(n)) xx(n) = mat%pml_pos(n)
                        end if
                    end if
                    pos   = (xx(n) - pf%MinBound(n)) / pf%step(n)
                    ii(n) = floor(pos)
                    if (ii(n) < 0)            ii(n) = 0
                    if (ii(n) >= pf%NN(n)-1)  ii(n) = pf%NN(n)-2
                    aa(n) = pos - ii(n)
                    if (aa(n) < 0._fpp) aa(n) = 0._fpp
                    if (aa(n) > 1._fpp) aa(n) = 1._fpp
                end do
                val =       (1._fpp-aa(0))*(1._fpp-aa(1))*pf%var(ii(0)  , ii(1)  )
                val = val + (       aa(0))*(1._fpp-aa(1))*pf%var(ii(0)+1, ii(1)  )
                val = val + (1._fpp-aa(0))*(       aa(1))*pf%var(ii(0)  , ii(1)+1)
                val = val + (       aa(0))*(       aa(1))*pf%var(ii(0)+1, ii(1)+1)
                field(i,j) = val
            end do
        end do
    end subroutine interpolate_elem_field_2d

    !-----------------------------------------------------------------------
    !> Driver: for every anisotropic-from-file subdomain, read the 2D Cstar
    !! properties and fill the GLL-point material -- elastic: Cij2d (+ Density
    !! from Rho); fluid: IDensTensor2d (+ invKappa2d) and the potential fields.
    subroutine read_aniso_material_2d(Tdomain)
        type(domain), intent(inout) :: Tdomain
        type(PropertyField2d) :: pf
        integer :: mat, n, ngllx, ngllz, c
        integer, parameter :: NEL = 7
        integer, parameter :: NFL = 4
        character(len=8), dimension(NEL) :: enames
        character(len=8), dimension(NFL) :: fnames
        integer, dimension(NEL) :: erow, ecol   ! Voigt (row,col); Rho => ecol<0
        real(fpp), dimension(:,:), allocatable :: fld
        logical :: any_aniso

        enames = [character(len=8) :: 'C11','C12','C13','C22','C23','C33','Rho']
        fnames = [character(len=8) :: 'iRho11','iRho12','iRho22','iKappa']
        erow   = [1, 1, 1, 2, 2, 3, 0]
        ecol   = [1, 2, 3, 2, 3, 3, -1]
        any_aniso = .false.

        do mat = 0, Tdomain%n_mat-1
            if (Tdomain%sSubDomain(mat)%material_definition /= 1) cycle  ! not a file material
            select case (Tdomain%sSubDomain(mat)%deftype)
            case (MATDEF_HOOKE_ANISO, CSTAR)
                any_aniso = .true.
                ! allocate Cij2d on every element of this subdomain
                do n = 0, Tdomain%n_elem-1
                    if (Tdomain%specel(n)%mat_index /= mat) cycle
                    ngllx = Tdomain%specel(n)%ngllx
                    ngllz = Tdomain%specel(n)%ngllz
                    if (.not. allocated(Tdomain%specel(n)%Cij2d)) &
                        allocate(Tdomain%specel(n)%Cij2d(3,3,0:ngllx-1,0:ngllz-1))
                    Tdomain%specel(n)%Cij2d = 0._fpp
                end do
                ! read each property once, interpolate into all elements of this subdomain
                do c = 1, NEL
                    pf%propName = enames(c)
                    call read_prop_field_2d(trim(Tdomain%sSubDomain(mat)%prop_file), pf)
                    do n = 0, Tdomain%n_elem-1
                        if (Tdomain%specel(n)%mat_index /= mat) cycle
                        ngllx = Tdomain%specel(n)%ngllx
                        ngllz = Tdomain%specel(n)%ngllz
                        if (allocated(fld)) deallocate(fld)
                        allocate(fld(0:ngllx-1,0:ngllz-1))
                        call interpolate_elem_field_2d(Tdomain, Tdomain%specel(n), Tdomain%sSubDomain(mat), pf, fld)
                        if (ecol(c) < 0) then                 ! Rho -> Density
                            Tdomain%specel(n)%Density(:,:) = fld
                        else
                            Tdomain%specel(n)%Cij2d(erow(c),ecol(c),:,:) = fld
                            Tdomain%specel(n)%Cij2d(ecol(c),erow(c),:,:) = fld  ! symmetric
                        end if
                    end do
                end do
                if (Tdomain%Mpi_var%my_rank == 0) &
                    write(*,'(a,i0,a)') ' [aniso] elastic Cstar tensor read & interpolated to GLL (subdomain ', mat, ')'
            case (MATDEF_FLUID_ANISO, CSTAR_FLUID)
                any_aniso = .true.
                ! allocate fluid-aniso (velocity-potential) fields on this subdomain's elements
                do n = 0, Tdomain%n_elem-1
                    if (Tdomain%specel(n)%mat_index /= mat) cycle
                    ngllx = Tdomain%specel(n)%ngllx
                    ngllz = Tdomain%specel(n)%ngllz
                    if (.not. allocated(Tdomain%specel(n)%IDensTensor2d)) &
                        allocate(Tdomain%specel(n)%IDensTensor2d(2,2,0:ngllx-1,0:ngllz-1))
                    if (.not. allocated(Tdomain%specel(n)%invKappa2d)) &
                        allocate(Tdomain%specel(n)%invKappa2d(0:ngllx-1,0:ngllz-1))
                    if (.not. allocated(Tdomain%specel(n)%Phi)) &
                        allocate(Tdomain%specel(n)%Phi(1:ngllx-2,1:ngllz-2))
                    if (.not. allocated(Tdomain%specel(n)%VelPhi)) &
                        allocate(Tdomain%specel(n)%VelPhi(1:ngllx-2,1:ngllz-2))
                    if (.not. allocated(Tdomain%specel(n)%ForcesFl)) &
                        allocate(Tdomain%specel(n)%ForcesFl(0:ngllx-1,0:ngllz-1))
                    Tdomain%specel(n)%IDensTensor2d = 0._fpp
                    Tdomain%specel(n)%invKappa2d     = 0._fpp
                    Tdomain%specel(n)%Phi            = 0._fpp
                    Tdomain%specel(n)%VelPhi         = 0._fpp
                    Tdomain%specel(n)%ForcesFl       = 0._fpp
                end do
                ! read each property once, interpolate into all elements of this subdomain
                do c = 1, NFL
                    pf%propName = fnames(c)
                    call read_prop_field_2d(trim(Tdomain%sSubDomain(mat)%prop_file), pf)
                    do n = 0, Tdomain%n_elem-1
                        if (Tdomain%specel(n)%mat_index /= mat) cycle
                        ngllx = Tdomain%specel(n)%ngllx
                        ngllz = Tdomain%specel(n)%ngllz
                        if (allocated(fld)) deallocate(fld)
                        allocate(fld(0:ngllx-1,0:ngllz-1))
                        call interpolate_elem_field_2d(Tdomain, Tdomain%specel(n), Tdomain%sSubDomain(mat), pf, fld)
                        select case (c)
                        case (1)  ! iRho11
                            Tdomain%specel(n)%IDensTensor2d(1,1,:,:) = fld
                        case (2)  ! iRho12 (symmetric)
                            Tdomain%specel(n)%IDensTensor2d(1,2,:,:) = fld
                            Tdomain%specel(n)%IDensTensor2d(2,1,:,:) = fld
                        case (3)  ! iRho22
                            Tdomain%specel(n)%IDensTensor2d(2,2,:,:) = fld
                        case (4)  ! iKappa = 1/kappa*
                            Tdomain%specel(n)%invKappa2d(:,:) = fld
                        end select
                    end do
                end do
                if (Tdomain%Mpi_var%my_rank == 0) &
                    write(*,'(a,i0,a)') ' [aniso] fluid-aniso material read & interpolated to GLL (subdomain ', mat, ')'
            end select
        end do

        if (allocated(fld)) deallocate(fld)
        if (allocated(pf%var)) deallocate(pf%var)

        ! NOTE: coefficient/mass builders (build_aniso_acoeff_2d, build_aniso_fluid_coeff_2d)
        ! are called separately from define_arrays, AFTER this read but BEFORE the mass is
        ! assembled to faces/vertices, so the aniso mass (rho from Cstar / 1/kappa) propagates.
        if (any_aniso .and. Tdomain%Mpi_var%my_rank == 0) &
            write(*,*) ' [aniso] Cstar material read & interpolated to GLL.'
    end subroutine read_aniso_material_2d

    !-----------------------------------------------------------------------
    !> Nearest-isotropic (Voigt) projection of a 2D Voigt tensor Cij2d onto Lame
    !! parameters (lambda, mu). Same Voigt average as homofft (homo/src/cut_cstar.f90:231)
    !! so a PML built from this matches the isotropic reference medium the
    !! homogenisation used. NOTE the convention: SEM2D's Cij2d stores the shear entry
    !! WITHOUT the factor 2 (C(3,3)=mu in the isotropic limit -- see build_aniso_acoeff_2d),
    !! whereas homo stores C(3,3)=2*mu; hence the 4*C(3,3) here vs homo's 2*C(3,3).
    !! Voigt indices: 1=xx, 2=zz, 3=xz. Stability-guarded so the resulting PML medium
    !! stays physical (mu>=0, lambda+2mu>0 -> vp real & >0). Reduces EXACTLY to (lambda,mu)
    !! for an isotropic input (C11=C22=lambda+2mu, C12=lambda, C33=mu).
    pure subroutine project_iso_cij2d(C, lambda, mu)
        real(fpp), dimension(3,3), intent(in) :: C
        real(fpp), intent(out)                :: lambda, mu
        real(fpp) :: K
        K  = 0.25_fpp  * (C(1,1) + 2._fpp*C(1,2) + C(2,2))                       ! lambda + mu
        mu = 0.125_fpp * (C(1,1) + C(2,2) - 2._fpp*C(1,2) + 4._fpp*C(3,3))       ! Voigt mu
        if (mu < 0._fpp) mu = 0._fpp             ! ponytail: physical floor, homogenised mu>0
        lambda = K - mu
        if (lambda + 2._fpp*mu <= 0._fpp) lambda = tiny(lambda)  ! ponytail: keep vp real
    end subroutine project_iso_cij2d

    !-----------------------------------------------------------------------
    !> Overwrite the CG (non-PML) Acoeff of every elastic aniso element
    !! with the full-Cij contraction  A = -Whei*Jac * D * C * D^T, where
    !!   D = [ xix 0 xiz ; etax 0 etaz ; 0 xiz xix ; 0 etaz etax ]   (4x3)
    !! and C = Cij2d (3x3 Voigt: 1=xx, 2=zz, 3=xz). This reduces EXACTLY to the
    !! isotropic Acoeff when C is built from (lambda,mu), so the isotropic limit is
    !! preserved. The kernel (compute_InternalForces_Elem) is unchanged: it already
    !! applies the symmetric 4x4 stored in Acoeff(0:9) to (dUx_dxi,dUx_deta,dUz_dxi,dUz_deta).
    subroutine build_aniso_acoeff_2d(Tdomain)
        type(domain), intent(inout) :: Tdomain
        integer :: n, mat, i, j, ngllx, ngllz, nover
        real(fpp) :: xix, xiz, etax, etaz, Jac, w
        real(fpp), dimension(4,3) :: D
        real(fpp), dimension(3,3) :: C
        real(fpp), dimension(4,4) :: A

        nover = 0
        do n = 0, Tdomain%n_elem-1
            if (.not. allocated(Tdomain%specel(n)%Cij2d)) cycle   ! only elastic aniso elements
            mat = Tdomain%specel(n)%mat_index
            ! Scope: continuous Galerkin, non-PML only. PML elements are handled
            ! upstream in define_arrays via the isotropic projection (project_iso_cij2d
            ! -> Lambda/Mu), so their absorbing medium matches the physical interface;
            ! only DG (genuinely unsupported) still warrants a warning.
            if (Tdomain%specel(n)%PML .or. Tdomain%specel(n)%type_DG /= GALERKIN_CONT) then
                if (Tdomain%specel(n)%type_DG /= GALERKIN_CONT .and. Tdomain%Mpi_var%my_rank == 0) &
                    write(*,'(a,i0,a)') ' [aniso] WARNING: aniso element ', n, &
                        ' is DG -- full-tensor kernel not supported, left on the isotropic path.'
                cycle
            end if
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            do j = 0, ngllz-1
                do i = 0, ngllx-1
                    xix  = Tdomain%specel(n)%InvGrad(i,j,0,0)
                    xiz  = Tdomain%specel(n)%InvGrad(i,j,1,0)
                    etax = Tdomain%specel(n)%InvGrad(i,j,0,1)
                    etaz = Tdomain%specel(n)%InvGrad(i,j,1,1)
                    Jac  = Tdomain%specel(n)%Jacob(i,j)
                    w    = Tdomain%sSubDomain(mat)%GLLwx(i) * Tdomain%sSubDomain(mat)%GLLwz(j)
                    D(1,1)=xix ;   D(1,2)=0._fpp; D(1,3)=xiz
                    D(2,1)=etax;   D(2,2)=0._fpp; D(2,3)=etaz
                    D(3,1)=0._fpp; D(3,2)=xiz ;   D(3,3)=xix
                    D(4,1)=0._fpp; D(4,2)=etaz;   D(4,3)=etax
                    C = Tdomain%specel(n)%Cij2d(:,:,i,j)
                    A = -w*Jac * matmul(matmul(D, C), transpose(D))
                    Tdomain%specel(n)%Acoeff(i,j,0) = A(1,1)
                    Tdomain%specel(n)%Acoeff(i,j,1) = A(1,2)
                    Tdomain%specel(n)%Acoeff(i,j,2) = A(1,3)
                    Tdomain%specel(n)%Acoeff(i,j,3) = A(1,4)
                    Tdomain%specel(n)%Acoeff(i,j,4) = A(2,2)
                    Tdomain%specel(n)%Acoeff(i,j,5) = A(2,3)
                    Tdomain%specel(n)%Acoeff(i,j,6) = A(2,4)
                    Tdomain%specel(n)%Acoeff(i,j,7) = A(3,3)
                    Tdomain%specel(n)%Acoeff(i,j,8) = A(3,4)
                    Tdomain%specel(n)%Acoeff(i,j,9) = A(4,4)
                    ! mass uses the density read from Cstar (Rho), not the subdomain constant
                    Tdomain%specel(n)%MassMat(i,j) = w * Tdomain%specel(n)%Density(i,j) * Jac
                end do
            end do
            nover = nover + 1
        end do
        if (nover > 0 .and. Tdomain%Mpi_var%my_rank == 0) &
            write(*,'(a,i0,a)') ' [aniso] CG Acoeff + mass rebuilt from Cstar on ', nover, ' element(s).'
    end subroutine build_aniso_acoeff_2d

    !-----------------------------------------------------------------------
    !> Build the fluid-aniso scalar stiffness coefficients and the diagonal
    !! 1/kappa mass for the velocity-potential domain.
    !!   AcoeffFl(:,:,0:2) = -Whei*Jac * G^T rho^-1 G   (2x2 symmetric: (1,1),(1,2),(2,2))
    !!   MassMatFl(:,:)    =  Whei * (1/kappa) * Jac
    !! with G = [[xix,etax],[xiz,etaz]]. Reduces to the standard acoustic operator
    !! when rho^-1 = (1/rho) I. CG interior uses the full tensor rho^-1; CG PML uses the
    !! isotropic projection (1/Density) I for the split-field velocity-potential PML.
    subroutine build_aniso_fluid_coeff_2d(Tdomain)
        type(domain), intent(inout) :: Tdomain
        integer :: n, mat, i, j, ngllx, ngllz, nover
        real(fpp) :: xix, xiz, etax, etaz, Jac, w, rinv
        real(fpp), dimension(2,2) :: G, R, A

        nover = 0
        do n = 0, Tdomain%n_elem-1
            if (.not. allocated(Tdomain%specel(n)%IDensTensor2d)) cycle   ! fluid-aniso only
            mat = Tdomain%specel(n)%mat_index
            ! DG is genuinely unsupported.
            if (Tdomain%specel(n)%type_DG /= GALERKIN_CONT) then
                if (Tdomain%Mpi_var%my_rank == 0) &
                    write(*,'(a,i0,a)') ' [aniso] WARNING: fluid-aniso element ', n, &
                        ' is DG -- velocity-potential kernel not supported, skipped.'
                cycle
            end if
            ! PML fluid-aniso: ISOTROPIC velocity-potential PML. AcoeffFl uses the projected
            ! isotropic inverse density R = (1/Density) I (Density = rho_iso, set in
            ! define_arrays). Mass (1/kappa) + DumpMass were overridden in the define_arrays
            ! PML block. The split-field kernels (Prediction_Elem_PML_VelPhi /
            ! compute_InternalForcesFl_PML_Elem) consume AcoeffFl. Do NOT touch MassMat here.
            if (Tdomain%specel(n)%PML) then
                ngllx = Tdomain%specel(n)%ngllx
                ngllz = Tdomain%specel(n)%ngllz
                if (.not. allocated(Tdomain%specel(n)%AcoeffFl)) &
                    allocate(Tdomain%specel(n)%AcoeffFl(0:ngllx-1,0:ngllz-1,0:2))
                do j = 0, ngllz-1
                    do i = 0, ngllx-1
                        xix  = Tdomain%specel(n)%InvGrad(i,j,0,0)
                        xiz  = Tdomain%specel(n)%InvGrad(i,j,1,0)
                        etax = Tdomain%specel(n)%InvGrad(i,j,0,1)
                        etaz = Tdomain%specel(n)%InvGrad(i,j,1,1)
                        Jac  = Tdomain%specel(n)%Jacob(i,j)
                        w    = Tdomain%sSubDomain(mat)%GLLwx(i) * Tdomain%sSubDomain(mat)%GLLwz(j)
                        rinv = 1._fpp / Tdomain%specel(n)%Density(i,j)   ! isotropic rho^-1
                        G(1,1)=xix; G(1,2)=etax
                        G(2,1)=xiz; G(2,2)=etaz
                        A = -w*Jac*rinv * matmul(transpose(G), G)
                        Tdomain%specel(n)%AcoeffFl(i,j,0) = A(1,1)
                        Tdomain%specel(n)%AcoeffFl(i,j,1) = A(1,2)
                        Tdomain%specel(n)%AcoeffFl(i,j,2) = A(2,2)
                    end do
                end do
                nover = nover + 1
                cycle
            end if
            ngllx = Tdomain%specel(n)%ngllx
            ngllz = Tdomain%specel(n)%ngllz
            if (.not. allocated(Tdomain%specel(n)%AcoeffFl)) &
                allocate(Tdomain%specel(n)%AcoeffFl(0:ngllx-1,0:ngllz-1,0:2))
            if (.not. allocated(Tdomain%specel(n)%MassMatFl)) &
                allocate(Tdomain%specel(n)%MassMatFl(0:ngllx-1,0:ngllz-1))
            do j = 0, ngllz-1
                do i = 0, ngllx-1
                    xix  = Tdomain%specel(n)%InvGrad(i,j,0,0)
                    xiz  = Tdomain%specel(n)%InvGrad(i,j,1,0)
                    etax = Tdomain%specel(n)%InvGrad(i,j,0,1)
                    etaz = Tdomain%specel(n)%InvGrad(i,j,1,1)
                    Jac  = Tdomain%specel(n)%Jacob(i,j)
                    w    = Tdomain%sSubDomain(mat)%GLLwx(i) * Tdomain%sSubDomain(mat)%GLLwz(j)
                    G(1,1)=xix; G(1,2)=etax
                    G(2,1)=xiz; G(2,2)=etaz
                    R = Tdomain%specel(n)%IDensTensor2d(:,:,i,j)
                    A = -w*Jac * matmul(matmul(transpose(G), R), G)
                    Tdomain%specel(n)%AcoeffFl(i,j,0) = A(1,1)
                    Tdomain%specel(n)%AcoeffFl(i,j,1) = A(1,2)
                    Tdomain%specel(n)%AcoeffFl(i,j,2) = A(2,2)
                    ! phi rides in component 0 of the 2-comp field, so the diagonal
                    ! mass for the potential equation goes into Elem%MassMat (1/kappa).
                    Tdomain%specel(n)%MassMatFl(i,j) = w * Tdomain%specel(n)%invKappa2d(i,j) * Jac
                    Tdomain%specel(n)%MassMat(i,j)   = Tdomain%specel(n)%MassMatFl(i,j)
                end do
            end do
            nover = nover + 1
        end do
        if (nover > 0 .and. Tdomain%Mpi_var%my_rank == 0) &
            write(*,'(a,i0,a)') ' [aniso] fluid-aniso mass(1/kappa) + stiffness built on ', nover, ' element(s).'
    end subroutine build_aniso_fluid_coeff_2d

end module build_prop_files_2d
