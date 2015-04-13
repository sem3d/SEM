!>
!! \file define_arrays.f90
!! \brief Contient la routine Define_Arrays()
!! \author
!! \version 1.0
!! \date
!!
!<
module mdefinitions
    use sdomain
    use mpi
    use scomm
    implicit none

    public :: define_arrays
    private :: define_Acoeff_iso, define_Acoeff_fluid
    private :: define_Acoeff_PML_iso, define_Acoeff_PML_fluid
    private :: define_alpha_PML
    private :: define_PML_DumpInit, define_PML_DumpEnd
    private :: assemble_DumpMass
    private :: define_FPML_DumpInit, define_FPML_DumpEnd

contains

subroutine Define_Arrays(Tdomain)
    use constants
    use model_earthchunk
    use model_prem
    implicit none

    type (domain), intent (INOUT), target :: Tdomain
    integer :: n, mat
    real, dimension(:,:,:), allocatable :: Whei

    if( Tdomain%earthchunk_isInit/=0) then
        call load_earthchunk(Tdomain%earthchunk_file, Tdomain%earthchunk_delta_lon, Tdomain%earthchunk_delta_lat)
    endif
    do n = 0,Tdomain%n_elem-1
       mat = Tdomain%specel(n)%mat_index
       if ( mat < 0 .or. mat >= Tdomain%n_mat ) then
          print*, "ERROR : inconsistent material index = ", mat
          stop
       end if
       !!! Attribute elastic properties from material !!!
       ! Sets Lambda, Mu, Qmu, ... from mat
       call init_material_properties(Tdomain, Tdomain%specel(n), Tdomain%sSubdomain(mat))
       ! Compute MassMat and Whei (with allocation)
       call init_local_mass_mat(Tdomain, Tdomain%specel(n), Tdomain%sSubdomain(mat),Whei)
       ! Compute Acoeff (for PML)
       call init_element_acoeff(Tdomain, Tdomain%specel(n), Tdomain%sSubdomain(mat),Whei)
       ! Computes DumpS, DumpMass (local),and for FPML :  Iv and Is
       if (Tdomain%specel(n)%PML) then
          call init_pml_properties(Tdomain, Tdomain%specel(n), Tdomain%sSubdomain(mat),Whei)
       end if
       deallocate(Whei)
    enddo
    ! Here we have local mass matrix (not assembled) on elements and
    ! each of faces, edges, vertices containing assembled (on local processor only) mass matrix
    if( Tdomain%earthchunk_isInit/=0) then
         call clean_earthchunk()
    endif
    !- defining Neumann properties (Btn: the complete normal term, ponderated
    !      by Gaussian weights)
    if(Tdomain%logicD%neumann_local_present)then
        call define_FEV_Neumann(Tdomain)
    endif

    call init_solid_fluid_interface(Tdomain)
    call assemble_mass_matrices(Tdomain)
    call finalize_pml_properties(Tdomain)
    call inverse_mass_mat(Tdomain)

    return
end subroutine define_arrays

subroutine init_solid_fluid_interface(Tdomain)
    type (domain), intent (INOUT), target :: Tdomain
    !
    integer :: nf, nnf, n
    integer :: dir
    integer :: ngll1, ngll2
    integer :: ngllx, nglly, ngllz
    !- defining Solid/Fluid faces'properties
    if(Tdomain%logicD%SF_local_present)then
        !  Btn: the complete normal term, ponderated by GLL weights
        call define_Face_SF(Tdomain)
        ! density (which links VelPhi and pressure)
        do nf = 0,Tdomain%SF%SF_n_faces-1
            nnf = Tdomain%SF%SF_Face(nf)%Face(0)
            if(nnf < 0) cycle
            n = Tdomain%sFace(nnf)%which_elem
            dir = Tdomain%sFace(nnf)%dir
            ngll1 = Tdomain%SF%SF_Face(nf)%ngll1
            ngll2 = Tdomain%SF%SF_Face(nf)%ngll2
            ngllx = Tdomain%specel(n)%ngllx
            nglly = Tdomain%specel(n)%nglly
            ngllz = Tdomain%specel(n)%ngllz
            if(dir == 0)then
                Tdomain%SF%SF_Face(nf)%density(0:ngll1-1,0:ngll2-1) =    &
                    Tdomain%specel(n)%density(0:ngllx-1,0:nglly-1,0)
            else if(dir == 1)then
                Tdomain%SF%SF_Face(nf)%density(0:ngll1-1,0:ngll2-1) =    &
                    Tdomain%specel(n)%density(0:ngllx-1,0,0:ngllz-1)
            else if(dir == 2)then
                Tdomain%SF%SF_Face(nf)%density(0:ngll1-1,0:ngll2-1) =    &
                    Tdomain%specel(n)%density(ngllx-1,0:nglly-1,0:ngllz-1)
            else if(dir == 3)then
                Tdomain%SF%SF_Face(nf)%density(0:ngll1-1,0:ngll2-1) =    &
                    Tdomain%specel(n)%density(0:ngllx-1,nglly-1,0:ngllz-1)
            else if(dir == 4)then
                Tdomain%SF%SF_Face(nf)%density(0:ngll1-1,0:ngll2-1) =    &
                    Tdomain%specel(n)%density(0,0:nglly-1,0:ngllz-1)
            else if(dir == 5)then
                Tdomain%SF%SF_Face(nf)%density(0:ngll1-1,0:ngll2-1) =    &
                    Tdomain%specel(n)%density(0:ngllx-1,0:nglly-1,ngllz-1)
            end if
        end do
    endif
end subroutine init_solid_fluid_interface

subroutine assemble_mass_matrices(Tdomain)
    implicit none
    type (domain), intent (INOUT), target :: Tdomain

    integer :: n, indsol, indpml, nsol, nsolpml
    integer :: i,j,k,nf,ne,nv,idx
    real :: Mass

    ! Couplage à l'interface solide / PML
    do n = 0,Tdomain%nbInterfSolPml-1
        indsol = Tdomain%InterfSolPml(n,0)
        indpml = Tdomain%InterfSolPml(n,1)
        Mass = Tdomain%MassMatSol(indsol) + Tdomain%MassMatSolPml(indpml)
        Tdomain%MassMatSol(indsol) = Mass
        Tdomain%MassMatSolPml(indpml+0) = Mass
        Tdomain%MassMatSolPml(indpml+1) = Mass
        Tdomain%MassMatSolPml(indpml+2) = Mass
    enddo

    if(Tdomain%Comm_data%ncomm > 0)then
        do n = 0,Tdomain%Comm_data%ncomm-1
            ! Domain SOLID
            k = 0
            do i=0,Tdomain%Comm_data%Data(n)%nsol-1
                idx = Tdomain%Comm_data%Data(n)%IGiveS(i)
                Tdomain%Comm_data%Data(n)%Give(k) = Tdomain%MassMatSol(idx)
                k=k+1
            end do

            ! Domain SOLID PML
            do i=0,Tdomain%Comm_data%Data(n)%nsolpml-1
                idx = Tdomain%Comm_data%Data(n)%IGiveSPML(i)
                Tdomain%Comm_data%Data(n)%Give(k+0) = Tdomain%DumpMass(idx)
                Tdomain%Comm_data%Data(n)%Give(k+1) = Tdomain%DumpMass(idx+1)
                Tdomain%Comm_data%Data(n)%Give(k+2) = Tdomain%DumpMass(idx+2)
                Tdomain%Comm_data%Data(n)%Give(k+3) = Tdomain%MassMatSolPml(idx)
                k=k+4
            end do

            ! Domain FLUID
            do i=0,Tdomain%Comm_data%Data(n)%nflu-1
                idx = Tdomain%Comm_data%Data(n)%IGiveF(i)
                Tdomain%Comm_data%Data(n)%Give(k) = Tdomain%MassMatFlu(idx)
                k=k+1
            end do

            ! Domain FLUID PML
            do i=0,Tdomain%Comm_data%Data(n)%nflupml-1
                idx = Tdomain%Comm_data%Data(n)%IGiveFPML(i)
                !TODO
                Tdomain%Comm_data%Data(n)%Give(k) = Tdomain%MassMatFluPml(idx)
                k=k+1
            end do

            Tdomain%Comm_data%Data(n)%nsend = k
        end do
        
        ! Exchange
        call exchange_sem_var(Tdomain, 103, Tdomain%Comm_data)

        ! Take
        do n = 0,Tdomain%Comm_data%ncomm-1
            ! Domain SOLID
            k = 0
            do i=0,Tdomain%Comm_data%Data(n)%nsol-1
                idx = Tdomain%Comm_data%Data(n)%ITakeS(i)
                Tdomain%MassMatSol(idx) = Tdomain%MassMatSol(idx) + &
                                          Tdomain%Comm_data%Data(n)%Take(k)
                k = k + 1
            end do

            ! Domain SOLID PML
            do i=0,Tdomain%Comm_data%Data(n)%nsolpml-1
                idx = Tdomain%Comm_data%Data(n)%ITakeSPML(i)
                Tdomain%DumpMass(idx+0) = Tdomain%DumpMass(idx+0) + Tdomain%Comm_data%Data(n)%Take(k+0)
                Tdomain%DumpMass(idx+1) = Tdomain%DumpMass(idx+1) + Tdomain%Comm_data%Data(n)%Take(k+1)
                Tdomain%DumpMass(idx+2) = Tdomain%DumpMass(idx+2) + Tdomain%Comm_data%Data(n)%Take(k+2)
                Tdomain%MassMatSolPml(idx+0) = Tdomain%MassMatSolPml(idx+0) + Tdomain%Comm_data%Data(n)%Take(k+3)
                Tdomain%MassMatSolPml(idx+1) = Tdomain%MassMatSolPml(idx+1) + Tdomain%Comm_data%Data(n)%Take(k+3)
                Tdomain%MassMatSolPml(idx+2) = Tdomain%MassMatSolPml(idx+2) + Tdomain%Comm_data%Data(n)%Take(k+3)
                k = k + 4
            end do

            ! Domain FLUID
            do i=0,Tdomain%Comm_data%Data(n)%nflu-1
                idx = Tdomain%Comm_data%Data(n)%ITakeF(i)
                Tdomain%MassMatFlu(idx) = Tdomain%MassMatFlu(idx) + Tdomain%Comm_data%Data(n)%Take(k)
                k=k+1
            end do

            ! Domain FLUID PML
            do i=0,Tdomain%Comm_data%Data(n)%nflupml-1
                idx = Tdomain%Comm_data%Data(n)%ITakeFPML(i)
                !TODO
                Tdomain%MassMatFluPml(idx) = Tdomain%MassMatFluPml(idx) + Tdomain%Comm_data%Data(n)%Take(k)
                k=k+1
            end do
        end do
    end if

    !!! XXXXXXXXX QUID echange PML->SOLIDE inter CPU... ?
    return
end subroutine assemble_mass_matrices
!---------------------------------------------------------------------------------------
subroutine inverse_mass_mat(Tdomain)
  type (domain), intent (INOUT), target :: Tdomain
  !
  integer :: ngllx, nglly, ngllz, n, nf, ne, nv
  real, dimension(:,:,:), allocatable  :: LocMassMat

    if (Tdomain%ngll_s /= 0)    Tdomain%MassMatSol(:) = 1d0/Tdomain%MassMatSol(:)
    if (Tdomain%ngll_f /= 0)    Tdomain%MassMatFlu(:) = 1d0/Tdomain%MassMatFlu(:)
    if (Tdomain%ngll_pmls /= 0) Tdomain%MassMatSolPml(:) = 1d0/Tdomain%MassMatSolPml(:)

end subroutine inverse_mass_mat


subroutine init_material_properties(Tdomain, specel, mat)
  type (domain), intent (INOUT), target :: Tdomain
  type (element), intent(inout) :: specel
  type (subdomain), intent(in) :: mat
  !
  integer :: ngllx, nglly, ngllz
  !
  ngllx = specel%ngllx
  nglly = specel%nglly
  ngllz = specel%ngllz

  !    integration de la prise en compte du gradient de proprietes



  select case( mat%material_definition)
  case( MATERIAL_CONSTANT )
     !    on copie toujours le materiau de base
     specel%Density = mat%Ddensity
     specel%Lambda = mat%DLambda
     specel%Kappa = mat%DKappa
     specel%Mu = mat%DMu
     !    si le flag gradient est actif alors on peut changer les proprietes
  case( MATERIAL_EARTHCHUNK )
     call initialize_material_earthchunk(specel, mat, Tdomain%GlobCoord, size(Tdomain%GlobCoord,2))

 case( MATERIAL_PREM )
     call initialize_material_prem(specel, mat, Tdomain%GlobCoord, size(Tdomain%GlobCoord,2))

  case( MATERIAL_GRADIENT )
     !    on copie toujours le materiau de base
     specel%Density = mat%Ddensity
     specel%Lambda = mat%DLambda
     specel%Kappa = mat%DKappa
     specel%Mu = mat%DMu
     !    si le flag gradient est actif alors on peut changer les proprietes
     if ( Tdomain%logicD%grad_bassin ) then
        call initialize_material_gradient(Tdomain, specel, mat)
     endif
  end select

  if ((.not. specel%PML) .and. (Tdomain%n_sls>0))  then
     if(specel%solid) then
        if (Tdomain%aniso) then
           specel%sl%Q = mat%Qmu
        else
           specel%sl%Qs = mat%Qmu
           specel%sl%Qp = mat%Qpression
        endif
     endif
  endif
end subroutine init_material_properties


subroutine init_element_acoeff(Tdomain,specel,mat,whei)
  type (domain), intent (INOUT), target :: Tdomain
  type (element), intent(inout) :: specel
  type (subdomain), intent(in) :: mat
  real, dimension(:,:,:), allocatable, intent(in) :: Whei
  !
  integer :: ngllx, nglly, ngllz
  real, dimension(:,:,:), allocatable :: xix,xiy,xiz
  real, dimension(:,:,:), allocatable :: etax,etay,etaz
  real, dimension(:,:,:), allocatable :: zetax,zetay,zetaz
  real, dimension(:,:,:), allocatable :: RKmod

  ngllx = specel%ngllx
  nglly = specel%nglly
  ngllz = specel%ngllz

  allocate(xix(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(xiy(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(xiz(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(etax(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(etay(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(etaz(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(zetax(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(zetay(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(zetaz(0:ngllx-1,0:nglly-1,0:ngllz-1))
  allocate(RKmod(0:ngllx-1,0:nglly-1,0:ngllz-1))

  xix = specel%InvGrad(:,:,:,0,0)
  xiy = specel%InvGrad(:,:,:,1,0)
  xiz = specel%InvGrad(:,:,:,2,0)

  etax = specel%InvGrad(:,:,:,0,1)
  etay = specel%InvGrad(:,:,:,1,1)
  etaz = specel%InvGrad(:,:,:,2,1)

  zetax = specel%InvGrad(:,:,:,0,2)
  zetay = specel%InvGrad(:,:,:,1,2)
  zetaz = specel%InvGrad(:,:,:,2,2)

  RKmod = specel%Lambda + 2. * specel%Mu

  !- verif. for fluid part
  if(.not. specel%solid .and. maxval(specel%Mu) > 1.d-5) stop "Fluid element with a non null shear modulus."


  !- parts of the internal forces terms: Acoeff; to be compared to
  !  general expressions = products of material properties and nabla operators

  if(.not. specel%PML .and. .false.)then
      ! On n'utilise pas les Acoef pour l'instant
      if(specel%solid)then
          call define_Acoeff_iso(ngllx,nglly,ngllz,Rkmod,specel%Mu,specel%Lambda,xix,xiy,xiz,    &
              etax,etay,etaz,zetax,zetay,zetaz,Whei,specel%Jacob,specel%sl%Acoeff)
      else   ! fluid case
          call define_Acoeff_fluid(ngllx,nglly,ngllz,specel%Density,xix,xiy,xiz,    &
              etax,etay,etaz,zetax,zetay,zetaz,Whei,specel%Jacob,specel%fl%Acoeff)
      end if
  endif

  if(specel%PML)then   ! PML case: valid for solid and fluid parts
     if(specel%solid)then
        call define_Acoeff_PML_iso(ngllx,nglly,ngllz,Rkmod,specel%Mu,specel%Lambda, &
             xix,xiy,xiz,    &
             etax,etay,etaz, &
             zetax,zetay,zetaz, &
             Whei,specel%Jacob,specel%sl%Acoeff)
     else  ! fluid case
        call define_Acoeff_PML_fluid(ngllx,nglly,ngllz,specel%Density, &
             xix,xiy,xiz,    &
             etax,etay,etaz, &
             zetax,zetay,zetaz,&
             Whei,specel%Jacob,specel%fl%Acoeff)
     end if
  end if
  deallocate(xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,RKmod)
end subroutine init_element_acoeff

subroutine init_pml_properties(Tdomain,specel,mat,Whei)
    type (domain), intent (INOUT), target :: Tdomain
    type (element), intent(inout) :: specel
    type (subdomain), intent(in) :: mat
    real, dimension(:,:,:), allocatable, intent(in) :: Whei
    !
    integer :: ngllx, nglly, ngllz
    real, dimension(:,:,:), allocatable :: temp_PMLx,temp_PMLy
    real, dimension(:,:,:), allocatable :: RKmod
    real, dimension(:,:,:), allocatable :: wx,wy,wz
    real :: freq

    if (.not. specel%PML) stop "init_pml_properties should not be called for non-pml element"

    ngllx = specel%ngllx
    nglly = specel%nglly
    ngllz = specel%ngllz

    allocate(RKmod(0:ngllx-1,0:nglly-1,0:ngllz-1))
    RKmod = specel%Lambda + 2. * specel%Mu


    ! PML case: valid for solid and fluid parts

     !- definition of the attenuation coefficient in PMLs (alpha in the literature)
     allocate(wx(0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate(wy(0:ngllx-1,0:nglly-1,0:ngllz-1))
     allocate(wz(0:ngllx-1,0:nglly-1,0:ngllz-1))

     call define_alpha_PML(mat%Px,0,mat%Left, &
          ngllx,nglly,ngllz,ngllx,Tdomain%n_glob_points,Tdomain%GlobCoord,       &
          mat%GLLcx,RKmod(:,0,0),                            &
          specel%Density(:,0,0),specel%Iglobnum(0,0,0),    &
          specel%Iglobnum(ngllx-1,0,0),mat%Apow,  &
          mat%npow,wx)
     call define_alpha_PML(mat%Py,1,mat%Forward, &
          ngllx,nglly,ngllz,nglly,Tdomain%n_glob_points,Tdomain%GlobCoord,          &
          mat%GLLcy,RKmod(0,:,0),                               &
          specel%Density(0,:,0),specel%Iglobnum(0,0,0),       &
          specel%Iglobnum(0,nglly-1,0),mat%Apow,     &
          mat%npow,wy)
     call define_alpha_PML(mat%Pz,2,mat%Down, &
          ngllx,nglly,ngllz,ngllz,Tdomain%n_glob_points,Tdomain%GlobCoord,       &
          mat%GLLcz,RKmod(0,0,:),                            &
          specel%Density(0,0,:),specel%Iglobnum(0,0,0),    &
          specel%Iglobnum(0,0,ngllz-1),mat%Apow,  &
          mat%npow,wz)


     !- M-PMLs
     if(Tdomain%logicD%MPML)then
        allocate(temp_PMLx(0:ngllx-1,0:nglly-1,0:ngllz-1))
        allocate(temp_PMLy(0:ngllx-1,0:nglly-1,0:ngllz-1))
        temp_PMLx(:,:,:) = wx(:,:,:)
        temp_PMLy(:,:,:) = wy(:,:,:)
        wx(:,:,:) = wx(:,:,:)+Tdomain%MPML_coeff*(wy(:,:,:)+wz(:,:,:))
        wy(:,:,:) = wy(:,:,:)+Tdomain%MPML_coeff*(temp_PMLx(:,:,:)+wz(:,:,:))
        wz(:,:,:) = wz(:,:,:)+Tdomain%MPML_coeff*(temp_PMLx(:,:,:)+temp_PMLy(:,:,:))
        deallocate(temp_PMLx,temp_PMLy)
     end if

     !- strong formulation for stresses. Dumped mass elements, convolutional terms.
     if(specel%FPML)then
        freq = mat%freq
     else
        freq = 1d0
     end if

     ! Compute DumpS(x,y,z) and DumpMass(0,1,2)
     call define_PML_DumpInit(ngllx,nglly,ngllz,mat%Dt,freq,wx,specel%MassMat, &
         specel%xpml%DumpSx,specel%xpml%DumpMass(:,:,:,0))
     call define_PML_DumpInit(ngllx,nglly,ngllz,mat%Dt,freq,wy,specel%MassMat, &
         specel%xpml%DumpSy,specel%xpml%DumpMass(:,:,:,1))
     call define_PML_DumpInit(ngllx,nglly,ngllz,mat%Dt,freq,wz,specel%MassMat, &
         specel%xpml%DumpSz,specel%xpml%DumpMass(:,:,:,2))

     call assemble_DumpMass(Tdomain,specel)

     if(specel%FPML)then
        ! Compute Is(xyz) and Iv(xyz)
        call define_FPML_DumpInit(ngllx,nglly,ngllz,mat%Dt,        &
             mat%freq,wx,specel%MassMat, &
             specel%xpml%DumpSx,specel%slpml%Isx,specel%slpml%Ivx)
        call define_FPML_DumpInit(ngllx,nglly,ngllz,mat%Dt,        &
             mat%freq,wy,specel%MassMat,             &
             specel%xpml%DumpSy,specel%slpml%Isy,specel%slpml%Ivy)
        call define_FPML_DumpInit(ngllx,nglly,ngllz,mat%Dt,        &
             mat%freq,wz,specel%MassMat,             &
             specel%xpml%DumpSz,specel%slpml%Isz,specel%slpml%Ivz)
     endif
     deallocate(wx,wy,wz)


  deallocate(RKmod)

end subroutine init_pml_properties

subroutine finalize_pml_properties(Tdomain)
  type (domain), intent (INOUT), target :: Tdomain
  !
  integer :: n, nf, ne, nv
  integer :: ngllx, nglly, ngllz
  integer :: ngll1, ngll2, ngll

    do n = 0,Tdomain%n_elem-1
        ngllx = Tdomain%specel(n)%ngllx
        nglly = Tdomain%specel(n)%nglly
        ngllz = Tdomain%specel(n)%ngllz

        if(Tdomain%specel(n)%PML)then   ! dumped masses in PML
            ! Compute DumpV
            call define_PML_DumpEnd(Tdomain%ngll_pmls, Tdomain%MassMatSolPml, &
                Tdomain%DumpMass, Tdomain%champs0%DumpV)
            deallocate(Tdomain%specel(n)%xpml%DumpMass)
        end if
     end do
end subroutine finalize_pml_properties

subroutine init_local_mass_mat(Tdomain, specel, mat, Whei)
  type (domain), intent (INOUT), target :: Tdomain
  type (element), intent(inout) :: specel
  type (subdomain), intent(in) :: mat
  real, dimension(:,:,:), allocatable, intent(out) :: Whei
  !
  integer :: i, j, k
  integer :: n, ind, indsol, indpml
  real :: Mass

  allocate(Whei(0:specel%ngllx-1,0:specel%nglly-1,0:specel%ngllz-1))
  !- general (element) weighting: tensorial property..
  do k = 0,specel%ngllz-1
     do j = 0,specel%nglly-1
        do i = 0,specel%ngllx-1
           Whei(i,j,k) = mat%GLLwx(i)*mat%GLLwy(j)*mat%GLLwz(k)
           if(specel%solid)then
              specel%MassMat(i,j,k) = Whei(i,j,k)*specel%Density(i,j,k)*specel%Jacob(i,j,k)
           else   ! fluid case: inertial term ponderation by the inverse of the bulk modulus
              specel%MassMat(i,j,k) = Whei(i,j,k)*specel%Jacob(i,j,k)/specel%Lambda(i,j,k)
           end if
        enddo
     enddo
  enddo

  if(specel%PML)then ! PML part
    do k = 0,specel%ngllz-1
        do j = 0,specel%nglly-1
            do i = 0,specel%ngllx-1
                if (specel%solid) then
                    ind = specel%slpml%ISolPml(i,j,k)
                    Tdomain%MassMatSolPml(ind) = Tdomain%MassMatSolPml(ind) &
                        + specel%MassMat(i,j,k)
                    Tdomain%MassMatSolPml(ind+1) = Tdomain%MassMatSolPml(ind)
                    Tdomain%MassMatSolPml(ind+2) = Tdomain%MassMatSolPml(ind)
                else
                    !todo
                    !Tdomain%MassMatFlu(specel%Iflu(i,j,k)) = Tdomain%MassMatFlu(specel%Iflu(i,j,k)) &
                    !    + specel%MassMat(i,j,k)
                endif
            enddo
        enddo
    enddo
  else ! Element non PML
    do k = 0,specel%ngllz-1
        do j = 0,specel%nglly-1
            do i = 0,specel%ngllx-1
                if (specel%solid) then
                    Tdomain%MassMatSol(specel%Isol(i,j,k)) = Tdomain%MassMatSol(specel%Isol(i,j,k)) &
                        + specel%MassMat(i,j,k)
                else
                    Tdomain%MassMatFlu(specel%Iflu(i,j,k)) = Tdomain%MassMatFlu(specel%Iflu(i,j,k)) &
                        + specel%MassMat(i,j,k)
                endif
            enddo
        enddo
    enddo
  endif
  !- mass matrix elements
end subroutine init_local_mass_mat

subroutine initialize_material_gradient(Tdomain, specel, mat)
  type (domain), intent (INOUT), target :: Tdomain
  type (element), intent(inout) :: specel
  type (subdomain), intent(in) :: mat
  !
  integer :: i,j,k,ipoint,iflag
  real :: Mu,Lambda,Kappa
  integer :: imx,imy,imz
  real :: xp,yp,zp,xfact
  real :: zg1,zd1,zg2,zd2,zz1,zz2,zfact
  real :: xd1,xg1
  real :: zrho,zrho1,zrho2,zCp,zCp1,zCp2,zCs,zCs1,zCs2
  integer :: ngllx, nglly, ngllz
  integer :: icolonne,jlayer
  !
  ngllx = specel%ngllx
  nglly = specel%nglly
  ngllz = specel%ngllz

  !    debut modification des proprietes des couches de materiaux
  !    bassin    voir programme Surface.f90

  !     n_layer nombre de couches
  !     n_colonne nombre de colonnes en x ici uniquement
  !     x_type == 0 on remet des materiaux  homogenes dans chaque bloc
  !     x_type == 1 on met des gradients pour chaque colonne en interpolant
  !     suivant z
  !       integer  :: n_colonne, n_layer, x_type
  !    x_coord correspond aux abscisses des colonnes
  !       real, pointer, dimension(:) :: x_coord
  !      z_layer profondeur de  linterface pour chaque x de colonne
  !      on definit egalement le materiaux par rho, Cp , Cs
  !       real, pointer, dimension(:,:) :: z_layer, z_rho, z_Cp, z_Cs


  !     on cherche tout d abord a localiser la maille a partir d un
  !     point de Gauss interne milieux (imx,imy,imz)
  imx = 1+(ngllx-1)/2
  imy = 1+(nglly-1)/2
  imz = 1+(ngllz-1)/2
  !     on impose qu une maille appartienne a un seul groupe de gradient de
  !     proprietes
  ipoint = specel%Iglobnum(imx,imy,imz)
  xp = Tdomain%GlobCoord(0,ipoint)
  yp = Tdomain%GlobCoord(1,ipoint)
  zp = Tdomain%GlobCoord(2,ipoint)
  iflag = 0
  if ( Tdomain%sBassin%x_type .eq. 2 ) then
     if ( zp .gt. Tdomain%sBassin%zmax) then
        iflag = 1
     endif
     if ( zp .lt. Tdomain%sBassin%zmin) then
        iflag = 1
     endif
  endif
  !  si iflag nul on peut faire les modifications  pour toute la maille
  if ( iflag .eq. 0 ) then
     icolonne = 0
     xfact = 0.D0
     do i = 1, Tdomain%sBassin%n_colonne
        if ( xp .ge. Tdomain%sBassin%x_coord(i-1) .and.  xp .lt. Tdomain%sBassin%x_coord(i) ) then
           icolonne = i-1
           xfact = (xp - Tdomain%sBassin%x_coord(i-1))/(Tdomain%sBassin%x_coord(i)-Tdomain%sBassin%x_coord(i-1))
        endif
     enddo

     jlayer = 0
     zfact = 0.D0
     do j = 1,Tdomain%sBassin%n_layer
        zg1 = Tdomain%sBassin%z_layer(icolonne,j-1)
        zd1 = Tdomain%sBassin%z_layer(icolonne+1,j-1)
        zz1 = zg1 + xfact*(zd1-zg1)
        zg2 = Tdomain%sBassin%z_layer(icolonne,j)
        zd2 = Tdomain%sBassin%z_layer(icolonne+1,j)
        zz2 = zg2 + xfact*(zd2-zg2)
        if ( zp .ge. zz1 .and. zp .lt. zz2 ) then
           jlayer = j-1
           zfact = ( zp -zz1)/(zz2-zz1)
        endif
     enddo
     !        limite du sous-domaine de gradient
     xg1 = Tdomain%sBassin%x_coord(icolonne)
     xd1 = Tdomain%sBassin%x_coord(icolonne+1)
     zg1 = Tdomain%sBassin%z_layer(icolonne,jlayer)
     zd1 = Tdomain%sBassin%z_layer(icolonne+1,jlayer)
     zg2 = Tdomain%sBassin%z_layer(icolonne,jlayer+1)
     zd2 = Tdomain%sBassin%z_layer(icolonne+1,jlayer+1)
     !
     zrho1 = Tdomain%sBassin%z_rho(icolonne,jlayer)
     zrho2 = Tdomain%sBassin%z_rho(icolonne,jlayer+1)
     zCp1 = Tdomain%sBassin%z_Cp(icolonne,jlayer)
     zCp2 = Tdomain%sBassin%z_Cp(icolonne,jlayer+1)
     zCs1 = Tdomain%sBassin%z_Cs(icolonne,jlayer)
     zCs2 = Tdomain%sBassin%z_Cs(icolonne,jlayer+1)

     !     boucle sur les points de Gauss de la maille
     !     xp, yp, zp coordonnees du point de Gauss
     do k = 0, ngllz -1
        do j = 0,nglly-1
           do i = 0,ngllx-1
              ipoint = specel%Iglobnum(i,j,k)
              xp = Tdomain%GlobCoord(0,ipoint)
              yp = Tdomain%GlobCoord(1,ipoint)
              zp = Tdomain%GlobCoord(2,ipoint)
              if ( Tdomain%sBassin%x_type .ge. 1 ) then
                 !    interpolations  pour le calcul du gradient
                 xfact = ( xp - xg1)/(xd1-xg1)
                 zz1 = zg1 + xfact*(zd1-zg1)
                 zz2 = zg2 + xfact*(zd2-zg2)
                 zfact = ( zp - zz1)/(zz2-zz1)
              else
                  !   on met les memes proprietes dans toute la maille
                  zfact = 0.D0
              endif
              zrho   = zrho1 + zfact*(zrho2-zrho1)
              zCp   = zCp1 + zfact*(zCp2-zCp1)
              zCs   = zCs1 + zfact*(zCs2-zCs1)
              !     calcul des coeffcients elastiques
              Mu     = zrho*zCs*zCs
              Lambda = zrho*(zCp*zCp - zCs*zCs)
              Kappa  = Lambda + 2.D0*Mu/3.D0

              specel%Density(i,j,k) = zrho
              specel%Lambda(i,j,k) = Lambda
              specel%Kappa(i,j,k) = Kappa
              specel%Mu(i,j,k) = Mu
           enddo
        enddo
     enddo
     !    fin test iflag nul
  endif
  !    fin modification des proprietes des couches de materiaux
end subroutine initialize_material_gradient

!---------------------------------------------------------------------------------------
subroutine define_Acoeff_iso(ngllx,nglly,ngllz,Rkmod,Rmu,Rlam, &
     xix,xiy,xiz, &
     etax,etay,etaz, &
     zetax,zetay,zetaz, &
     Whei,Jac,Acoeff)
    integer, intent(in)  :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Rkmod,Rmu,Rlam,  &
        xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:44), intent(out) :: Acoeff

    Acoeff(:,:,:,0) = -Whei*(RKmod*xix**2+Rmu*(xiy**2+xiz**2))*Jac
    Acoeff(:,:,:,1) = -Whei*(RKmod*xix*etax+Rmu*(xiy*etay+xiz*etaz))*Jac
    Acoeff(:,:,:,2) = -Whei*(RKmod*xix*zetax+Rmu*(xiy*zetay+xiz*zetaz))*Jac
    Acoeff(:,:,:,3) = -Whei*(Rlam+Rmu)*xix*xiy*Jac
    Acoeff(:,:,:,4) = -Whei*(Rlam*xix*etay+Rmu*xiy*etax)*Jac
    Acoeff(:,:,:,5) = -Whei*(Rlam*xix*zetay+Rmu*xiy*zetax)*Jac
    Acoeff(:,:,:,6) = -Whei*(Rlam+Rmu)*xix*xiz*Jac
    Acoeff(:,:,:,7) = -Whei*(Rlam*xix*etaz+Rmu*xiz*etax)*Jac
    Acoeff(:,:,:,8) = -Whei*(Rlam*xix*zetaz+rmu*xiz*zetax)*Jac
    Acoeff(:,:,:,9) = -Whei*(RKmod*etax**2+Rmu*(etay**2+etaz**2))*Jac
    Acoeff(:,:,:,10) = -Whei*(RKmod*etax*zetax+Rmu*(etay*zetay+etaz*zetaz))*Jac
    Acoeff(:,:,:,11) = -Whei*(Rlam*etax*xiy+Rmu*etay*xix)*Jac
    Acoeff(:,:,:,12) = -Whei*(Rlam+Rmu)*etay*etax*Jac
    Acoeff(:,:,:,13) = -Whei*(Rlam*etax*zetay+Rmu*etay*zetax)*Jac
    Acoeff(:,:,:,14) = -Whei*(Rlam*etax*xiz+Rmu*etaz*xix)*Jac
    Acoeff(:,:,:,15) = -Whei*(Rlam+Rmu)*etaz*etax*Jac
    Acoeff(:,:,:,16) = -Whei*(Rlam*etax*zetaz+Rmu*etaz*zetax)*Jac
    Acoeff(:,:,:,17) = -Whei*(RKmod*zetax**2+Rmu*(zetay**2+zetaz**2))*Jac
    Acoeff(:,:,:,18) = -Whei*(Rlam*zetax*xiy+Rmu*zetay*xix)*Jac
    Acoeff(:,:,:,19) = -Whei*(Rlam*zetax*etay+Rmu*zetay*etax)*Jac
    Acoeff(:,:,:,20) = -Whei*(Rlam+Rmu)*zetax*zetay*Jac
    Acoeff(:,:,:,21) = -Whei*(Rlam*zetax*xiz+Rmu*zetaz*xix)*Jac
    Acoeff(:,:,:,22) = -Whei*(Rlam*zetax*etaz+Rmu*zetaz*etax)*Jac
    Acoeff(:,:,:,23) = -Whei*(Rlam+Rmu)*zetax*zetaz*Jac
    Acoeff(:,:,:,24) = -Whei*(RKmod*xiy**2+Rmu*(xix**2+xiz**2))*Jac
    Acoeff(:,:,:,25) = -Whei*(RKmod*xiy*etay+Rmu*(xix*etax+xiz*etaz))*Jac
    Acoeff(:,:,:,26) = -Whei*(RKmod*xiy*zetay+Rmu*(xix*zetax+xiz*zetaz))*Jac
    Acoeff(:,:,:,27) = -Whei*(Rlam+Rmu)*xiy*xiz*Jac
    Acoeff(:,:,:,28) = -Whei*(Rlam*etaz*xiy+Rmu*etay*xiz)*Jac
    Acoeff(:,:,:,29) = -Whei*(Rlam*zetaz*xiy+Rmu*zetay*xiz)*Jac
    Acoeff(:,:,:,30) = -Whei*(RKmod*etay**2+Rmu*(etax**2+etaz**2))*Jac
    Acoeff(:,:,:,31) = -Whei*(RKmod*zetay*etay+Rmu*(zetax*etax+zetaz*etaz))*Jac
    Acoeff(:,:,:,32) = -Whei*(Rlam*etay*xiz+Rmu*etaz*xiy)*Jac
    Acoeff(:,:,:,33) = -Whei*(Rlam+Rmu)*etay*etaz*Jac
    Acoeff(:,:,:,34) = -Whei*(Rlam*zetaz*etay+Rmu*zetay*etaz)*Jac
    Acoeff(:,:,:,35) = -Whei*(RKmod*zetay**2+Rmu*(zetax**2+zetaz**2))*Jac
    Acoeff(:,:,:,36) = -Whei*(Rlam*xiz*zetay+Rmu*xiy*zetaz)*Jac
    Acoeff(:,:,:,37) = -Whei*(Rlam*zetay*etaz+Rmu*zetaz*etay)*Jac
    Acoeff(:,:,:,38) = -Whei*(Rlam+Rmu)*zetay*zetaz*Jac
    Acoeff(:,:,:,39) = -Whei*(RKmod*xiz**2+Rmu*(xix**2+xiy**2))*Jac
    Acoeff(:,:,:,40) = -Whei*(RKmod*xiz*etaz+Rmu*(xix*etax+xiy*etay))*Jac
    Acoeff(:,:,:,41) = -Whei*(RKmod*xiz*zetaz+Rmu*(xix*zetax+xiy*zetay))*Jac
    Acoeff(:,:,:,42) = -Whei*(RKmod*etaz**2+Rmu*(etax**2+etay**2))*Jac
    Acoeff(:,:,:,43) = -Whei*(RKmod*zetaz*etaz+Rmu*(zetax*etax+zetay*etay))*Jac
    Acoeff(:,:,:,44) = -Whei*(RKmod*zetaz**2+Rmu*(zetax**2+zetay**2))*Jac

    return

end subroutine define_Acoeff_iso
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine define_Acoeff_fluid(ngllx,nglly,ngllz,Density,xix,xiy,xiz,    &
    etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac,Acoeff)
    integer, intent(in)  :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Density,  &
        xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:5), intent(out) :: Acoeff

    Acoeff(:,:,:,0) = -Whei*(xix**2+xiy**2+xiz**2)*Jac/density
    Acoeff(:,:,:,1) = -Whei*(xix*etax+xiy*etay+xiz*etaz)*Jac/density
    Acoeff(:,:,:,2) = -Whei*(xix*zetax+xiy*zetay+xiz*zetaz)*Jac/density
    Acoeff(:,:,:,3) = -Whei*(etax**2+etay**2+etaz**2)*Jac/density
    Acoeff(:,:,:,4) = -Whei*(etax*zetax+etay*zetay+etaz*zetaz)*Jac/density
    Acoeff(:,:,:,5) = -Whei*(zetax**2+zetay**2+zetaz**2)*Jac/density


end subroutine define_Acoeff_fluid
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
subroutine define_Acoeff_PML_iso(ngllx,nglly,ngllz,Rkmod,Rmu,Rlam,xix,xiy,xiz,    &
    etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac,Acoeff)
    integer, intent(in)  :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: Rkmod,Rmu,Rlam,  &
        xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:35), intent(out) :: Acoeff


    Acoeff(:,:,:,0) = RKmod *xix
    Acoeff(:,:,:,1) = RKmod *etax
    Acoeff(:,:,:,2) = RKmod *zetax
    Acoeff(:,:,:,3) = RLam *xiy
    Acoeff(:,:,:,4) = RLam *etay
    Acoeff(:,:,:,5) = RLam *zetay
    Acoeff(:,:,:,6) = RLam *xiz
    Acoeff(:,:,:,7) = RLam *etaz
    Acoeff(:,:,:,8) = RLam *zetaz
    Acoeff(:,:,:,9) = RLam *xix
    Acoeff(:,:,:,10) = RLam *etax
    Acoeff(:,:,:,11) = RLam *zetax
    Acoeff(:,:,:,12) = RKmod *xiy
    Acoeff(:,:,:,13) = RKmod *etay
    Acoeff(:,:,:,14) = RKmod *zetay
    Acoeff(:,:,:,15) = RKmod *xiz
    Acoeff(:,:,:,16) = RKmod *etaz
    Acoeff(:,:,:,17) = RKmod *zetaz
    Acoeff(:,:,:,18) = RMu *xix
    Acoeff(:,:,:,19) = RMu *etax
    Acoeff(:,:,:,20) = RMu *zetax
    Acoeff(:,:,:,21) = RMu *xiy
    Acoeff(:,:,:,22) = RMu *etay
    Acoeff(:,:,:,23) = RMu *zetay
    Acoeff(:,:,:,24) = RMu *xiz
    Acoeff(:,:,:,25) = RMu *etaz
    Acoeff(:,:,:,26) = RMu *zetaz
    Acoeff(:,:,:,27) = -Whei * xix * Jac
    Acoeff(:,:,:,28) = -Whei * xiy * Jac
    Acoeff(:,:,:,29) = -Whei * xiz * Jac
    Acoeff(:,:,:,30) = -Whei * etax * Jac
    Acoeff(:,:,:,31) = -Whei * etay * Jac
    Acoeff(:,:,:,32) = -Whei * etaz * Jac
    Acoeff(:,:,:,33) = -Whei * zetax * Jac
    Acoeff(:,:,:,34) = -Whei * zetay * Jac
    Acoeff(:,:,:,35) = -Whei * zetaz * Jac

    return

end subroutine define_Acoeff_PML_iso
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine define_Acoeff_PML_fluid(ngllx,nglly,ngllz,density,xix,xiy,xiz,    &
    etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac,Acoeff)
    integer, intent(in)  :: ngllx,nglly,ngllz
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: density,  &
        xix,xiy,xiz,etax,etay,etaz,zetax,zetay,zetaz,Whei,Jac
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:17), intent(out) :: Acoeff


    Acoeff(:,:,:,0) = xix(:,:,:)/density(:,:,:)
    Acoeff(:,:,:,1) = etax(:,:,:)/density(:,:,:)
    Acoeff(:,:,:,2) = zetax(:,:,:)/density(:,:,:)
    Acoeff(:,:,:,3) = xiy(:,:,:)/density(:,:,:)
    Acoeff(:,:,:,4) = etay(:,:,:)/density(:,:,:)
    Acoeff(:,:,:,5) = zetay(:,:,:)/density(:,:,:)
    Acoeff(:,:,:,6) = xiz(:,:,:)/density(:,:,:)
    Acoeff(:,:,:,7) = etaz(:,:,:)/density(:,:,:)
    Acoeff(:,:,:,8) = zetaz(:,:,:)/density(:,:,:)
    Acoeff(:,:,:,9) = -whei(:,:,:)*xix(:,:,:)*Jac(:,:,:)
    Acoeff(:,:,:,10) = -whei(:,:,:)*etax(:,:,:)*Jac(:,:,:)
    Acoeff(:,:,:,11) = -whei(:,:,:)*zetax(:,:,:)*Jac(:,:,:)
    Acoeff(:,:,:,12) = -whei(:,:,:)*xiy(:,:,:)*Jac(:,:,:)
    Acoeff(:,:,:,13) = -whei(:,:,:)*etay(:,:,:)*Jac(:,:,:)
    Acoeff(:,:,:,14) = -whei(:,:,:)*zetay(:,:,:)*Jac(:,:,:)
    Acoeff(:,:,:,15) = -whei(:,:,:)*xiz(:,:,:)*Jac(:,:,:)
    Acoeff(:,:,:,16) = -whei(:,:,:)*etaz(:,:,:)*Jac(:,:,:)
    Acoeff(:,:,:,17) = -whei(:,:,:)*zetaz(:,:,:)*Jac(:,:,:)

    return

end subroutine define_Acoeff_PML_fluid
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
subroutine define_alpha_PML(lattenu,dir,ldir_attenu,ngllx,nglly,ngllz,ngll,n_pts,   &
    Coord,GLLc,Rkmod,Density,ind_min,ind_max,Apow,npow,alpha)
    !- routine determines attenuation profile in an PML layer (see Festa & Vilotte)
    !   dir = attenuation's direction, ldir_attenu = the logical giving the orientation
    logical, intent(in)   :: lattenu,ldir_attenu
    integer, intent(in) :: dir,ngllx,nglly,ngllz,ngll,n_pts,ind_min,ind_max,npow
    real, dimension(0:2,0:n_pts-1), intent(in) :: Coord
    real, dimension(0:ngll-1), intent(in) :: GLLc,RKmod,Density
    real, intent(in)  :: Apow
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: alpha
    integer  :: i
    real  :: dh
    real, dimension(0:ngll-1)  :: ri,vp
    real, external  :: pow

    if(.not. lattenu)then   ! no attenuation in the dir-direction
        alpha(:,:,:) = 0d0
    else  ! yes, attenuation in this dir-direction
        dh = Coord(dir,ind_min)
        dh = abs(Coord(dir,ind_max)-dh)
        if(ldir_attenu)then  ! Left in x, Forward in y, Down in z
            ri(:) = 0.5d0*(1d0+GLLc(ngll-1:0:-1))*float(ngll-1)
        else  ! Right in x, Backward in y, Up in z
            ri(:) = 0.5d0*(1d0+GLLc(0:ngll-1))*float(ngll-1)
        end if
        vp(:) = sqrt(Rkmod(:)/Density(:))
        select case(dir)
        case(0)  ! dir = x
            do i = 0,ngll-1
                alpha(i,0:,0:) = pow(ri(i),vp(i),ngll-1,dh,Apow,npow)
            end do
        case(1)  ! dir = y
            do i = 0,ngll-1
                alpha(0:,i,0:) = pow(ri(i),vp(i),ngll-1,dh,Apow,npow)
            end do
        case(2)  ! dir = z
            do i = 0,ngll-1
                alpha(0:,0:,i) = pow(ri(i),vp(i),ngll-1,dh,Apow,npow)
            end do
        end select
    end if

    return

end subroutine define_alpha_PML
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine define_FPML_DumpInit(ngllx,nglly,ngllz,dt,freq,alpha,MassMat,  &
    DumpS,Is,Iv)
    !- defining parameters related to stresses and mass matrix elements, in the case of
    !    a FPML, along a given splitted direction:
    integer, intent(in)  :: ngllx,nglly,ngllz
    real, intent(in) :: dt, freq
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: alpha
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: MassMat
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1), intent(in) :: DumpS
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: Is,Iv

    Is(:,:,:) = dt*alpha(:,:,:)*freq*DumpS(:,:,:,1)
    Iv(:,:,:) = MassMat(:,:,:)*Dt*alpha(:,:,:)*freq
    return

end subroutine define_FPML_DumpInit
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine define_PML_DumpInit(ngllx,nglly,ngllz,dt,freq,alpha,&
                                 MassMat,DumpS,DumpMass)
    !- defining parameters related to stresses and mass matrix elements, in the case of
    !    a PML, along a given splitted direction:
    integer, intent(in)  :: ngllx,nglly,ngllz
    real, intent(in) :: dt, freq
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: alpha
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(in) :: MassMat
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1,0:1), intent(out) :: DumpS
    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1), intent(out) :: DumpMass

    real, dimension(0:ngllx-1,0:nglly-1,0:ngllz-1)  :: Id

    Id = 1d0

    DumpS(:,:,:,1) = Id + 0.5d0*dt*alpha*freq
    DumpS(:,:,:,1) = 1d0/DumpS(:,:,:,1)
    DumpS(:,:,:,0) = (Id - 0.5d0*dt*alpha*freq)*DumpS(:,:,:,1)
    DumpMass(:,:,:) = 0.5d0*MassMat(:,:,:)*alpha(:,:,:)*dt*freq

    return
end subroutine define_PML_DumpInit

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine assemble_DumpMass(Tdomain,specel)
    type(domain), intent(inout) :: Tdomain
    type (element), intent(inout) :: specel

    integer :: i,j,k,m

    do m = 0,2
        do k = 0,specel%ngllz-1
            do j = 0,specel%nglly-1
                do i = 0,specel%ngllx-1
                    if (specel%solid) then
                        Tdomain%DumpMass(specel%slpml%ISolPml(i,j,k)+m) = &
                            Tdomain%DumpMass(specel%slpml%ISolPml(i,j,k)+m) &
                            + specel%xpml%DumpMass(i,j,k,m)
                    else
                        !todo
                        !Tdomain%DumpMass(specel%flpml%IFluPml(i,j,k),:) = &
                        !    Tdomain%DumpMass(specel%flpml%IFluPml(i,j,k),:) &
                        !    + specel%xpml%DumpMass(i,j,k,:)
                    endif
                enddo
            enddo
        enddo
    enddo

    return
end subroutine assemble_DumpMass
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine define_PML_DumpEnd(ngll_pmls,Massmat,DumpMass,DumpV)
    implicit none
    integer, intent(in)   :: ngll_pmls
    real, dimension(0:ngll_pmls-1), intent(in) :: MassMat, DumpMass
    real, dimension(0:ngll_pmls-1,0:1), intent(out) :: DumpV
    
    DumpV(:,1) = MassMat(:) + DumpMass(:)
    DumpV(:,1) = 1d0/DumpV(:,1)
    DumpV(:,0) = MassMat(:) - DumpMass(:)
    DumpV(:,0) = DumpV(:,0) * DumpV(:,1)

    return
end subroutine define_PML_DumpEnd
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine define_FPML_DumpEnd(ngllx,nglly,ngllz,DumpV,Iv)

    implicit none
    integer, intent(in)   :: ngllx,nglly,ngllz
    real, dimension(1:ngllx-2,1:nglly-2,1:ngllz-2,0:1), intent(in) :: DumpV
    real, dimension(:,:,:), allocatable, intent(inout) :: Iv
    real, dimension(1:ngllx-2,1:nglly-2,1:ngllz-2)  :: LocIv

    LocIv = Iv(1:ngllx-2,1:nglly-2,1:ngllz-2)
    deallocate(Iv)
    allocate(Iv(1:ngllx-2,1:nglly-2,1:ngllz-2))
    Iv = LocIv*DumpV(:,:,:,1)
    return

end subroutine define_FPML_DumpEnd
end module mdefinitions
!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
