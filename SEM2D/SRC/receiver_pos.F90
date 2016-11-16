!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file receiver_pos.F90
!!\brief Contient la subroutine ReceiverPosition.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module treceivers
    use sreceivers
    use sdomain
    implicit none
contains

! ###########################################################
!>
!! \brief
!<
subroutine ReceiverPosition(Tdomain)

    use sdomain
    use mlocations2d
    use constants, only : NCAPT_CACHE
    implicit none
    type (domain), intent (INOUT) :: Tdomain

    ! local variables

    integer, parameter :: NMAXEL=20
    integer, dimension(NMAXEL) :: elems
    double precision, dimension(0:1,NMAXEL) :: coordloc
    double precision, parameter :: EPS = 1D-13
    integer :: nrec, nmax, i, j, n_el, ngllx, ngllz, mat
    double precision :: xc, zc, xi, eta, outx, outz
    logical :: inside

    ! Find the nearest point to the real location of the receivers in GLL scheme

    do nrec = 0, Tdomain%n_receivers - 1

        ! If "collapse mode" is choosen for captors, the captors are relocated on the nearest gll node
        if (Tdomain%capt_loc_type == CAPT_NEAREST_NODE) then
            call relocalize_on_nearest_gauss_pt(Tdomain, nrec)
        endif

        xc = Tdomain%sReceiver(nrec)%xrec
        zc = Tdomain%sReceiver(nrec)%zrec
        nmax = NMAXEL

        ! Find receiver location
        call find_location(Tdomain, xc, zc, nmax, elems, coordloc)

        Tdomain%sReceiver(nrec)%located_here = .false.
        do i=1,nmax! When the receiver is in the mesh
            inside = .true.
            n_el = elems(i)
            xi   = coordloc(0,i)
            eta  = coordloc(1,i)
            if (xi<(-1-EPS) .or. eta<(-1-EPS)) inside = .false.
            if (xi>( 1+EPS) .or. eta>( 1+EPS)) inside = .false.
            if (inside) then
                write (*,'(a,i5,a,i4,a,f10.5,a,f10.5,a,i10,a,f20.17,a,f20.17,a)') " Receiver ", nrec, " : found on proc. ", Tdomain%Mpi_var%my_rank,    &
                                                                                  ", (xc, zc) : (", xc, ", ", zc, ") <=> (element, xi, eta) : (", n_el, &
                                                                                  ", ", xi, ", ", eta, ")"
                Tdomain%sReceiver(nrec)%located_here = .true.
                Tdomain%sReceiver(nrec)%nr           = n_el
                Tdomain%sReceiver(nrec)%xi           = xi
                Tdomain%sReceiver(nrec)%eta          = eta

                exit ! Choose the first one to be consistent with previous behavior
            end if
        end do

        ! Compute interpolation coefficient

        if (Tdomain%sReceiver(nrec)%located_here) then
            n_el  = Tdomain%sReceiver(nrec)%nr
            mat   = Tdomain%specel(n_el)%mat_index
            ngllx = Tdomain%specel(n_el)%ngllx
            ngllz = Tdomain%specel(n_el)%ngllz
            allocate (Tdomain%sReceiver(nrec)%Interp_Coeff(0:ngllx-1,0:ngllz-1))
            do j = 0, ngllz-1
                call pol_lagrange(ngllz, Tdomain%sSubdomain(mat)%GLLcz, j, Tdomain%sReceiver(nrec)%eta, outz)
                do i = 0, ngllx -1
                    call pol_lagrange(ngllx, Tdomain%sSubdomain(mat)%GLLcx, i, Tdomain%sReceiver(nrec)%xi, outx)
                    Tdomain%sReceiver(nrec)%Interp_Coeff(i,j) = outx*outz
                enddo
            enddo
        endif
    enddo

    ! Prepare Post-Processing if HDG
    if(Tdomain%type_elem==GALERKIN_HDG_RP) &
        call prepare_HDG_postprocess(Tdomain)

    ! Prepare to store receiver trace
    i = 0
    do nrec = 0, Tdomain%n_receivers-1
        if (Tdomain%sReceiver(nrec)%located_here) i= i + 1
    enddo
    if (i > 0) then
        allocate (Tdomain%Store_Trace(0:1,0:i-1,0:NCAPT_CACHE-1))
    endif
end subroutine ReceiverPosition


! ###########################################################
!>
!! \brief This subroutine saves the velocity trace corresponding to the receiver number "it".
!! It computes an interpolation of the velocity at the receiver location from all the nodal
!! velocity values of the receiver's element.
!<
subroutine save_trace (Tdomain, it)
    use sdomain
    use orientation
    implicit none

    type (domain), intent (INOUT) :: Tdomain
    integer, intent (IN) :: it

    integer :: ir, nr, i,j, ngllx, ngllz, nsta, ncache, ind, mat
    real :: dum0, dum1
    real, dimension (:,:,:), allocatable :: Field

    ncache = mod(it, NCAPT_CACHE)
    nsta = Tdomain%n_receivers
    ind = 0
    do ir = 0, nsta-1
        if (Tdomain%sReceiver(ir)%located_here) then
            nr = Tdomain%sReceiver(ir)%nr
            dum0 = 0; dum1 = 0
            ngllx = Tdomain%specel(nr)%ngllx
            ngllz = Tdomain%specel(nr)%ngllz

            allocate (Field(0:ngllx-1,0:ngllz-1,0:1))

            call gather_elem_veloc(Tdomain, nr, Field, .false.)
            if (Tdomain%type_timeInteg == TIME_INTEG_RK4) &
                Field = Tdomain%specel(nr)%Veloc

            if (Tdomain%specel(nr)%type_DG==GALERKIN_HDG_RP) then
                mat = Tdomain%specel(nr)%mat_index
                ngllx = ngllx +1 ; ngllz = ngllz +1
                call compute_InternalForces_HDG_Weak(Tdomain%specel(nr), &
                    Tdomain%sSubDomain(mat)%hprimex, Tdomain%sSubDomain(mat)%hTprimez)
                call get_Vhat_f2el(Tdomain,nr)
                call Compute_Traces (Tdomain%specel(nr))
                call inversion_massmat(Tdomain%specel(nr))
                call postprocess_HDG(Tdomain,nr,ir,Field)
                call compute_TracFace (Tdomain%specel(nr))
                Tdomain%specel(nr)%TracFace(:,0) = Tdomain%specel(nr)%TracFace(:,0) * Tdomain%specel(nr)%Coeff_Integr_Faces(:)
                Tdomain%specel(nr)%TracFace(:,1) = Tdomain%specel(nr)%TracFace(:,1) * Tdomain%specel(nr)%Coeff_Integr_Faces(:)
            endif

            do j = 0,ngllz-1
                do i =0,ngllx -1
                    dum0 = dum0 + Tdomain%sReceiver(ir)%Interp_Coeff(i,j) * Field(i,j,0)
                    dum1 = dum1 + Tdomain%sReceiver(ir)%Interp_Coeff(i,j) * Field(i,j,1)
                enddo
            enddo

            Tdomain%Store_Trace(0,ind,ncache) = dum0
            Tdomain%Store_Trace(1,ind,ncache) = dum1
            ind = ind + 1
            deallocate (Field)
        endif
    enddo

    if (ncache == NCAPT_CACHE-1 .or. it == Tdomain%TimeD%NtimeMax-1) then
        call dump_trace(Tdomain)
    end if

    return
end subroutine save_trace



! ###########################################################
!>
!! \brief
!<
subroutine dump_trace (Tdomain)

    use sdomain
    use semdatafiles
    implicit none

    type(Domain), intent (IN) :: Tdomain

    integer :: i, it, it0, it1, ind, ncache, i_err
    character(Len=MAX_FILE_SIZE) :: fnamef
    real :: rtime

    !dumping the traces
    it0 = NCAPT_CACHE*(Tdomain%TimeD%ntime/NCAPT_CACHE)
    if (Tdomain%logicD%run_restart .and. (Tdomain%TimeD%ntime - Tdomain%TimeD%iter_reprise) < NCAPT_CACHE) then
        it0 = Tdomain%TimeD%iter_reprise+1 ! if restart, start dump from restart iteration (not from last NCAPT_CACHE block)
    end if
    it1 = Tdomain%TimeD%ntime
    ind = 0

    do i = 0,Tdomain%n_receivers-1
        if (Tdomain%sReceiver(i)%located_here) then
            call semname_capteur_type(Tdomain%sReceiver(i)%name, ".vel", fnamef)
            open (31, file=trim(fnamef), action="write", form="formatted", position="append", iostat=i_err)
            if (i_err .ne. 0) stop "open trace file KO"
            rtime=(it0+1)*Tdomain%TimeD%dtmin
            do it = it0, it1
                ncache=it-it0
                if (Tdomain%logicD%run_restart .and. (Tdomain%TimeD%ntime - Tdomain%TimeD%iter_reprise) < NCAPT_CACHE) then
                    ncache = mod(it,NCAPT_CACHE) ! if restart, start dump from restart iteration (not from last NCAPT_CACHE block)
                end if
                write(31,*) rtime, Tdomain%Store_Trace(0,ind,ncache), Tdomain%Store_Trace(1,ind,ncache)
                rtime = rtime + Tdomain%TimeD%dtmin
            enddo
            close (31)
            ind = ind + 1
        endif
    enddo
    return
end subroutine dump_trace


! ###########################################################
!>
!! \brief
!<
subroutine read_receiver_file(Tdomain)
    use sdomain
    use semdatafiles
    type(domain), intent(inout) :: Tdomain
    real :: xrec, zrec
    character(Len=100) :: recname
    character(Len=MAX_FILE_SIZE) :: fnamef
    integer :: i

    if (.not. Tdomain%logicD%save_trace) then
        return
    end if

    call semname_read_inputmesh_parametrage(Tdomain%station_file,fnamef)
    open(14,file=fnamef, status="old")

    read (14,*) Tdomain%n_receivers
    read (14,*)
    allocate (Tdomain%sReceiver(0:Tdomain%n_receivers-1))
    do i = 0, Tdomain%n_receivers-1
        read(14,*) recname, xrec, zrec
        Tdomain%sReceiver(i)%Xrec = xrec
        Tdomain%sReceiver(i)%Zrec = zrec
        Tdomain%sReceiver(i)%name = recname
    enddo
    close (14)

end subroutine read_receiver_file



! ###########################################################
!>
!! \brief This subroutine computes the postprocessed
!<
subroutine postprocess_HDG(Tdomain,nelem,nrec,field)

    use selement
    implicit none

    type(domain), intent (in) :: Tdomain
    integer, intent(in)       :: nelem
    integer, intent(in)       :: nrec
    real, dimension (:,:,:), allocatable, intent(inout) :: field
    type(element),  pointer   :: Elem
    type(receiver), pointer   :: rec
    real, dimension (:,:,:), allocatable :: p
    real, dimension (:,:),   allocatable :: xix, xiz, etax, etaz, jac, Amat
    real,  dimension (:),     allocatable :: smbr, WORK
    integer :: i, j, k, l, r, s
    integer :: ngx, ngz, nlin, nddl, lwork, info, mat
    real    :: aux1, aux2, res

    Elem => Tdomain%specel(nelem)
    rec  => Tdomain%sReceiver(nrec)

    ngx = Elem%ngllx ; ngz = Elem%ngllz ; nddl = 2*(ngx+1)*(ngz+1)
    mat = Elem%mat_index

    ! Allocate fields for solving linear system
    allocate(p(0:ngx-1,0:ngz-1,0:2))
    allocate(smbr(1:nddl+3))
    deallocate(field)
    allocate(field(0:ngx,0:ngz,0:1))
    p = Elem%Forces(:,:,0:2)

    ! Allocate inverse Jacobian functions
    allocate (xix (0:ngx+1,0:ngz+1)) ; allocate (xiz (0:ngx+1,0:ngz+1))
    allocate (etax(0:ngx+1,0:ngz+1)) ; allocate (etaz(0:ngx+1,0:ngz+1))
    allocate (Jac(0:ngx+1,0:ngz+1))
    xix  = rec%InvGrad(:,:,0,0) ; xiz  = rec%InvGrad(:,:,1,0)
    etax = rec%InvGrad(:,:,0,1) ; etaz = rec%InvGrad(:,:,1,1)
    jac  = rec%JacobWheiN2(:,:)

    ! Computation of the 2nd member of the linear system for post-processing
    ! Here, (i,j) indexes corresponds to test function w_ij which gives line number of MatPostProc
    ! indexes (k,l) correspond to the nodal value of the quantity p_11, p_22, p_12
    ! Finally, indexes (r,s) correspond to the quadrature points.

    ! Loop over test functions
    do i = 0,ngx
        do j = 0,ngz
            nlin = Ind(i,j,0,ngx,ngz) ; res = 0.
            ! Loop over quadrature points
            do r = 0,ngx+1
                do s = 0,ngz+1
                    aux1 = 0. ; aux2 = 0.
                    ! Loop over nodal p values
                    do k = 0,ngx-1
                        do l = 0,ngz-1
                            aux1 = aux1 + (p(k,l,0)*xix(r,s)  + p(k,l,2)*xiz(r,s)) &
                                *rec%ReinterpNX(k,r)*rec%ReinterpNZ(l,s)
                            aux2 = aux2 + (p(k,l,0)*etax(r,s) + p(k,l,2)*etaz(r,s)) &
                                *rec%ReinterpNX(k,r)*rec%ReinterpNZ(l,s)
                        enddo
                    enddo
                    aux1 = aux1 * rec%ReinterpZ(j,s) * rec%hprimex(i,r)
                    aux2 = aux2 * rec%ReinterpX(i,r) * rec%hprimez(j,s)
                    res = res + (aux1 + aux2) * jac(r,s)
                enddo
            enddo
            smbr(nlin) = res

            nlin = Ind(i,j,1,ngx,ngz) ; res = 0.
            ! Loop over quadrature points
            do r = 0,ngx+1
                do s = 0,ngz+1
                    aux1 = 0. ; aux2 = 0.
                    ! Loop over nodal p values
                    do k = 0,ngx-1
                        do l = 0,ngz-1
                            aux1 = aux1 + (p(k,l,1)*xiz(r,s)  + p(k,l,2)*xix(r,s)) &
                                *rec%ReinterpNX(k,r)*rec%ReinterpNZ(l,s)
                            aux2 = aux2 + (p(k,l,1)*etaz(r,s) + p(k,l,2)*etax(r,s)) &
                                *rec%ReinterpNX(k,r)*rec%ReinterpNZ(l,s)
                        enddo
                    enddo
                    aux1 = aux1 * rec%ReinterpZ(j,s) * rec%hprimex(i,r)
                    aux2 = aux2 * rec%ReinterpX(i,r) * rec%hprimez(j,s)
                    res = res + (aux1 + aux2) * jac(r,s)
                enddo
            enddo
            smbr(nlin) = res
        enddo
    enddo

    ! Avant-avant-derniere et avant-derniere lignes du systeme, avec coefficient de preconditionnement :
    smbr(nddl+1) = sum(Elem%Acoeff(:,:,12)*Elem%Veloc(:,:,0))
    smbr(nddl+2) = sum(Elem%Acoeff(:,:,12)*Elem%Veloc(:,:,1))

    ! Derniere ligne du systeme :
    res = 0.
    do r = 0,ngx-1
        do s = 0,ngz-1
            ! Loop over nodal values of the regular velocity field
            do k = 0,ngx-1
                do l = 0,ngz-1
                    res = res + (Tdomain%Ssubdomain(mat)%hprimex(k,r)*Delta(l,s)*Elem%Veloc(k,l,1)*Elem%InvGrad(r,s,0,0) & ! xix
                                +Tdomain%Ssubdomain(mat)%hprimez(l,s)*Delta(k,r)*Elem%Veloc(k,l,1)*Elem%InvGrad(r,s,0,1) & ! etax
                                -Tdomain%Ssubdomain(mat)%hprimex(k,r)*Delta(l,s)*Elem%Veloc(k,l,0)*Elem%InvGrad(r,s,1,0) & ! xiz
                                -Tdomain%Ssubdomain(mat)%hprimez(l,s)*Delta(k,r)*Elem%Veloc(k,l,0)*Elem%InvGrad(r,s,1,1))& ! etaz
                                *Elem%Acoeff(r,s,12)
                enddo
            enddo
        enddo
    enddo
    smbr(nddl+3) = res

    ! Coefficient preconditionnement :
    smbr(nddl+1) = xix(0,0)*xix(0,0) * smbr(nddl+1)
    smbr(nddl+2) = xix(0,0)*xix(0,0) * smbr(nddl+2)
    smbr(nddl+3) = xix(0,0)          * smbr(nddl+3)


    ! Solve rectangular System using least squares method
    allocate (Amat(1:nddl+3,1:nddl))
    Amat = rec%MatPostProc
    allocate (WORK(1:1000))
    ! Query the optimal workspace
    lwork = -1
    CALL DGELS('N', nddl+3, nddl, 1, Amat, nddl+3, smbr, nddl+3, WORK, lwork, info)
    lwork = INT(WORK(1))
    deallocate(WORK) ; allocate (WORK(lwork))
    call dgels('N', nddl+3, nddl, 1, Amat, nddl+3, smbr, nddl+3, WORK, lwork, info)
    IF( INFO.GT.0 ) THEN
        WRITE(*,*)'The diagonal element ',INFO,' of the triangular '
        WRITE(*,*)'factor of A is zero, so that A does not have full '
        WRITE(*,*)'rank; the least squares solution could not be computed.'
        STOP
    END IF

    ! Reorganize the results as a matrix
    do i = 0,ngx
        do j = 0,ngz
            nlin = Ind(i,j,0,ngx,ngz)
            field(i,j,0) = smbr(nlin)
            nlin = Ind(i,j,1,ngx,ngz)
            field(i,j,1) = smbr(nlin)
        enddo
    enddo
    deallocate(p, smbr, xix, xiz, etax, etaz, jac, WORK, Amat)

end subroutine postprocess_HDG


! ###########################################################
!>
!! \brief This subroutine builds the matrix to be inverted for the HDG post-process.
!! It also builds the reinterpolation matrices used for evaluation of N-order polynomials
!! on a N+1 polynomial Lagrange interpolation.
!<
subroutine prepare_HDG_postprocess(Tdomain)

    use sdomain
    use shape_lin
    use shape_quad
    use splib, only : zelegl, welegl, dmlegl
    implicit none

    type(domain), intent(inout) :: Tdomain
    real, dimension (:),   allocatable :: GLLcxN2, GLLwxN2, GLLpolxN2, GLLczN2, GLLwzN2, GLLpolzN2
    real, dimension (:),   allocatable :: GLLcxN1, GLLwxN1, GLLpolxN1, GLLczN1, GLLwzN1, GLLpolzN1
    real, dimension (:,:), allocatable :: coord, hprimex_aux, hprimez_aux
    real, dimension (0:1,0:1) :: LocInvGrad
    type(receiver), pointer   :: rec
    real    :: xi, eta, Jac, res, outx, outz
    integer :: nrec, ngx, ngz, nr, nddl, i_aus, mat, i, j

    do nrec = 0, Tdomain%n_receivers-1
        ! Matrices Allocation used for Post-Process
        rec => Tdomain%sReceiver(nrec)
        nr = rec%nr
        if (Tdomain%specel(nr)%Acoustic) cycle
        ngx = Tdomain%specel(nr)%ngllx
        ngz = Tdomain%specel(nr)%ngllz
        nddl= 2*(ngx+1)*(ngz+1)
        mat = Tdomain%specel(nr)%mat_index

        ! Remark : Element's interpolation order is (ngx-1,ngz-1)
        ! when in post-traitment, it becomes (ngx,ngz).
        allocate(rec%InvGrad    (0:ngx+1,0:ngz+1,0:1,0:1))
        allocate(rec%JacobWheiN1(0:ngx,  0:ngz))
        allocate(rec%JacobWheiN2(0:ngx+1,0:ngz+1))
        allocate(rec%ReinterpNX (0:ngx-1,0:ngx+1))
        allocate(rec%ReinterpNZ (0:ngz-1,0:ngz+1))
        allocate(rec%ReinterpX  (0:ngx,  0:ngx+1))
        allocate(rec%ReinterpZ  (0:ngz,  0:ngz+1))
        allocate(rec%hprimex(0:ngx,0:ngx+1))
        allocate(rec%hprimez(0:ngz,0:ngz+1))
        allocate(hprimex_aux(0:ngx+1,0:ngx))
        allocate(hprimez_aux(0:ngz+1,0:ngz))
        allocate(rec%MatPostProc(1:nddl+3,1:nddl))
        deallocate(rec%Interp_Coeff)
        allocate(rec%Interp_Coeff(0:ngx,0:ngz))

        ! Computation of GLL positions, local wheights, and Lagrange polynomial derivatives
        allocate (GLLcxN2(0:ngx+1)); allocate (GLLwxN2(0:ngx+1)); allocate (GLLpolxN2(0:ngx+1))
        allocate (GLLczN2(0:ngz+1)); allocate (GLLwzN2(0:ngz+1)); allocate (GLLpolzN2(0:ngz+1))
        allocate (GLLcxN1(0:ngx));   allocate (GLLwxN1(0:ngx));   allocate (GLLpolxN1(0:ngx))
        allocate (GLLczN1(0:ngz));   allocate (GLLwzN1(0:ngz));   allocate (GLLpolzN1(0:ngz))

        call zelegl (ngx+1,GLLcxN2,GLLpolxN2)
        call welegl (ngx+1,GLLcxN2,GLLpolxN2,GLLwxN2)
        call zelegl (ngz+1,GLLczN2,GLLpolzN2)
        call welegl (ngz+1,GLLczN2,GLLpolzN2,GLLwzN2)

        call zelegl (ngx,GLLcxN1,GLLpolxN1)
        call welegl (ngx,GLLcxN1,GLLpolxN1,GLLwxN1)
        call zelegl (ngz,GLLczN1,GLLpolzN1)
        call welegl (ngz,GLLczN1,GLLpolzN1,GLLwzN1)

        call mderlag(ngx+1,GLLcxN1,ngx+2,GLLcxN2,rec%hprimex)
        call mderlag(ngz+1,GLLczN1,ngz+2,GLLczN2,rec%hprimez)

        ! Computation of Jacobian and InvGrad for the source element
        if (Tdomain%n_nodes == 4) then
            allocate (coord(0:1,0:3))
            do i=0,3
                i_aus = Tdomain%specel(nr)%Control_Nodes(i)
                coord(0,i) = Tdomain%Coord_Nodes(0,i_aus)
                coord(1,i) = Tdomain%Coord_Nodes(1,i_aus)
            end do
            do j = 0,ngz+1
                eta = GLLczN2(j)
                do i = 0,ngx+1
                    xi = GLLcxN2(i)
                    ! Computation of the derivative matrix
                    call shape4_local2jacob(coord, xi, eta, LocInvGrad)
                    call invert2(LocInvGrad, Jac)
                    rec%InvGrad(i,j,0:1,0:1) = LocInvGrad(0:1,0:1)
                    rec%JacobWheiN2(i,j) = Jac*GLLwxN2(i)*GLLwzN2(j)
                enddo
            enddo
            do j = 0,ngz
                eta = GLLczN1(j)
                do i = 0,ngx
                    xi = GLLcxN1(i)
                    ! Computation of the derivative matrix
                    call shape4_local2jacob(coord, xi, eta, LocInvGrad)
                    call invert2(LocInvGrad, Jac)
                    rec%JacobWheiN1(i,j) = Jac*GLLwxN1(i)*GLLwzN1(j)
                enddo
            enddo

        elseif (Tdomain%n_nodes == 8) then
            allocate (coord(0:1,0:7))
            do i=0,7
                i_aus = Tdomain%specel(nr)%Control_Nodes(i)
                coord(0,i) = Tdomain%Coord_Nodes(0,i_aus)
                coord(1,i) = Tdomain%Coord_Nodes(1,i_aus)
            end do
            do j = 0,ngz+1
                eta =  GLLczN2(j)
                do i = 0,ngx+1
                    xi = GLLcxN2(i)
                    ! Computation of the derivative matrix
                    call shape8_local2jacob(coord, xi, eta, LocInvGrad)
                    call invert2(LocInvGrad, Jac)
                    rec%InvGrad (i,j,0:1,0:1) = LocInvGrad(0:1,0:1)
                    rec%JacobWheiN2(i,j) = Jac*GLLwxN2(i)*GLLwzN2(j)
                enddo
            enddo
            do j = 0,ngz
                eta = GLLczN1(j)
                do i = 0,ngx
                    xi = GLLcxN1(i)
                    ! Computation of the derivative matrix
                    call shape8_local2jacob(coord, xi, eta, LocInvGrad)
                    call invert2(LocInvGrad, Jac)
                    rec%JacobWheiN1(i,j) = Jac*GLLwxN1(i)*GLLwzN1(j)
                enddo
            enddo
        endif

        ! Computation of reinterpolation matrices evaluating (N = ngll-1) th order Lagrange polynomials
        ! on N+2 Gauss Lobatto points. Computes the L_i(\xi_j), with 0 \leq i \leq N+1 and 0 \leq j \leq N+2.
        do i=0,ngx-1
            do j=0,ngx+1 ! N+2 Gauss Lobato points
                xi = GLLcxN2(j)
                call pol_lagrange(ngx,Tdomain%sSubdomain(mat)%GLLcx,i,xi,res)
                rec%ReinterpNX(i,j) = res
            enddo
        enddo
        do i=0,ngz-1
            do j=0,ngz+1 ! N+2 Gauss Lobato points
                eta = GLLczN2(j)
                call pol_lagrange(ngz,Tdomain%sSubdomain(mat)%GLLcz,i,eta,res)
                rec%ReinterpNZ(i,j) = res
            enddo
        enddo

        ! Computation of reinterpolation matrices evaluating (N+1 = ngll) th order Lagrange polynomials
        ! on N+2 Gauss Lobatto points. Computes the L_i(\xi_j), with 0 \leq i \leq N+1 and 0 \leq j \leq N+2.
        do i=0,ngx
            do j=0,ngx+1 ! N+2 Gauss Lobato points
                xi = GLLcxN2(j)
                call pol_lagrange(ngx+1,GLLcxN1,i,xi,res)
                rec%ReinterpX(i,j) = res
            enddo
        enddo
        do i=0,ngz
            do j=0,ngz+1 ! N+2 Gauss Lobato points
                eta = GLLczN2(j)
                call pol_lagrange(ngz+1,GLLczN1,i,eta,res)
                rec%ReinterpZ(i,j) = res
            enddo
        enddo

        ! Coefficients d'interpolation pour le placement du capteur dans l'element
        ! As usual for any receiver...
        do j = 0,ngz
            call pol_lagrange(ngz+1,GLLczN1,j,rec%eta,outz)
            do i = 0,ngx
                call pol_lagrange(ngx+1,GLLcxN1,i,rec%xi,outx)
                rec%Interp_Coeff(i,j) = outx*outz
            enddo
        enddo

        deallocate(GLLcxN1,GLLwxN1,GLLpolxN1,GLLczN1,GLLwzN1,GLLpolzN1)
        deallocate(GLLcxN2,GLLwxN2,GLLpolxN2,GLLczN2,GLLwzN2,GLLpolzN2)
        deallocate(coord,hprimex_aux,hprimez_aux)

        ! Computation of the coefficients of matrix MatPostProc
        call build_MatPostProc(tdomain,nrec,ngx,ngz)

    enddo

end subroutine prepare_HDG_postprocess


! ###########################################################
!>
!! \brief This subroutine builds the matrix MatPostProc for the n-th receiver
!! that will be inverted in order to compute the postprocessed velocities v^*
!!
!<
subroutine build_MatPostProc(Tdomain,nrec,ngx,ngz)

    implicit none
    type(Domain), intent(inout) :: Tdomain
    integer, intent(in)         :: nrec, ngx, ngz
    type(receiver), pointer     :: rec
    real, dimension(:,:), allocatable :: xix, etax, xiz, etaz
    real    :: vx, vz
    integer :: i, j, k, l, r, s, nlin, ncol

    rec => Tdomain%sReceiver(nrec)

    ! Allocate inverse Jacobian functions
    allocate (xix (0:ngx+1,0:ngz+1)) ; allocate (xiz (0:ngx+1,0:ngz+1))
    allocate (etax(0:ngx+1,0:ngz+1)) ; allocate (etaz(0:ngx+1,0:ngz+1))
    xix  = rec%InvGrad(:,:,0,0)
    xiz  = rec%InvGrad(:,:,1,0)
    etax = rec%InvGrad(:,:,0,1)
    etaz = rec%InvGrad(:,:,1,1)
    rec%MatPostProc = 0.

    ! Compute entries of matrix MatPostProc
    ! Here, (i,j) indexes corresponds to test function w_ij which gives line number of MatPostProc
    ! indexes (k,l) corresponds to the nodal value of the postprocessed velocity v^* and gives
    ! the column number of MatPostProc. Finally, indexes (r,s) correspond to the quadrature points.

    ! Loop over test functions
    do i = 0,ngx
        do j = 0,ngz
            nlin = Ind(i,j,0,ngx,ngz)
            ! Loop over quadrature points
            do k = 0,ngx
                do l = 0,ngz
                    vx = 0. ; vz = 0.
                    ! Loop over nodal values of v^*
                    do r = 0,ngx+1
                        do s = 0,ngz+1
                        vx = vx + (rec%hprimex(k,r)*rec%ReinterpZ(l,s)*xix(r,s) + rec%ReinterpX(k,r)*rec%hprimez(l,s)*etax(r,s)) &
                                 *(rec%hprimex(i,r)*rec%ReinterpZ(j,s)*xix(r,s) + rec%ReinterpX(i,r)*rec%hprimez(j,s)*etax(r,s)) &
                                 * rec%JacobWheiN2(r,s) &
                                 +(rec%hprimex(k,r)*rec%ReinterpZ(l,s)*xiz(r,s) + rec%ReinterpX(k,r)*rec%hprimez(l,s)*etaz(r,s)) &
                                 *(rec%hprimex(i,r)*rec%ReinterpZ(j,s)*xiz(r,s) + rec%ReinterpX(i,r)*rec%hprimez(j,s)*etaz(r,s)) &
                                 * 0.5 * rec%JacobWheiN2(r,s)
                        vz = vz + (rec%hprimex(k,r)*rec%ReinterpZ(l,s)*xix(r,s) + rec%ReinterpX(k,r)*rec%hprimez(l,s)*etax(r,s)) &
                                 *(rec%hprimex(i,r)*rec%ReinterpZ(j,s)*xiz(r,s) + rec%ReinterpX(i,r)*rec%hprimez(j,s)*etaz(r,s)) &
                                 * 0.5 * rec%JacobWheiN2(r,s)
                        enddo
                    enddo
                    ncol = Ind(k,l,0,ngx,ngz)
                    rec%MatPostProc(nlin,ncol) = vx
                    ncol = Ind(k,l,1,ngx,ngz)
                    rec%MatPostProc(nlin,ncol) = vz
                enddo
            enddo
            nlin = Ind(i,j,1,ngx,ngz)
            ! Loop over nodal values of v^*
            do k = 0,ngx
                do l = 0,ngz
                    vx = 0. ; vz = 0.
                    ! Loop over quadrature points
                    do r = 0,ngx+1
                        do s = 0,ngz+1
                        vx = vx + (rec%hprimex(k,r)*rec%ReinterpZ(l,s)*xiz(r,s) + rec%ReinterpX(k,r)*rec%hprimez(l,s)*etaz(r,s)) &
                                 *(rec%hprimex(i,r)*rec%ReinterpZ(j,s)*xix(r,s) + rec%ReinterpX(i,r)*rec%hprimez(j,s)*etax(r,s)) &
                                 * 0.5 * rec%JacobWheiN2(r,s)
                        vz = vz + (rec%hprimex(k,r)*rec%ReinterpZ(l,s)*xiz(r,s) + rec%ReinterpX(k,r)*rec%hprimez(l,s)*etaz(r,s)) &
                                 *(rec%hprimex(i,r)*rec%ReinterpZ(j,s)*xiz(r,s) + rec%ReinterpX(i,r)*rec%hprimez(j,s)*etaz(r,s)) &
                                 * rec%JacobWheiN2(r,s) &
                                 +(rec%hprimex(k,r)*rec%ReinterpZ(l,s)*xix(r,s) + rec%ReinterpX(k,r)*rec%hprimez(l,s)*etax(r,s)) &
                                 *(rec%hprimex(i,r)*rec%ReinterpZ(j,s)*xix(r,s) + rec%ReinterpX(i,r)*rec%hprimez(j,s)*etax(r,s)) &
                                 * 0.5 * rec%JacobWheiN2(r,s)
                        enddo
                    enddo
                    ncol = Ind(k,l,0,ngx,ngz)
                    rec%MatPostProc(nlin,ncol) = vx
                    ncol = Ind(k,l,1,ngx,ngz)
                    rec%MatPostProc(nlin,ncol) = vz
                enddo
            enddo
        enddo
    enddo

    ! Building the three last lines of the matrix system
    ! The two following lines corresponds to the two equalities in the mean between
    ! the post-processed velocities and the regular velocities.
    nlin = 2*(ngx+1)*(ngz+1)+1
    do k = 0,ngx
        do l = 0,ngz
            ncol = Ind(k,l,0,ngx,ngz)
            rec%MatPostProc(nlin,ncol)  = rec%JacobWheiN1(k,l)
            ncol = Ind(k,l,1,ngx,ngz)
            rec%MatPostProc(nlin+1,ncol)= rec%JacobWheiN1(k,l)
        enddo
    enddo
    ! The last line of the matrix corresponds to the equality in the mean between
    ! the rotational of the post-processed velocity and the regular velocities.
    nlin = 2*(ngx+1)*(ngz+1)+3
    do k = 0,ngx
        do l = 0,ngz
            vx = 0. ; vz = 0.
            do r = 0,ngx+1
                do s = 0,ngz+1
                    vx = vx - (rec%hprimex(k,r)*rec%ReinterpZ(l,s)*xiz(r,s) + rec%ReinterpX(k,r)*rec%hprimez(l,s)*etaz(r,s)) * rec%JacobWheiN2(r,s)
                    vz = vz + (rec%hprimex(k,r)*rec%ReinterpZ(l,s)*xix(r,s) + rec%ReinterpX(k,r)*rec%hprimez(l,s)*etax(r,s)) * rec%JacobWheiN2(r,s)
                enddo
            enddo
            ncol = Ind(k,l,0,ngx,ngz)
            rec%MatPostProc(nlin,ncol) = vx
            ncol = Ind(k,l,1,ngx,ngz)
            rec%MatPostProc(nlin,ncol) = vz
        enddo
    enddo

    ! Multiplication coefficients pour preconditionnement
    nlin = 2*(ngx+1)*(ngz+1)+1
    rec%MatPostProc(nlin:nlin+1,:) = xix(0,0)*xix(0,0) * rec%MatPostProc(nlin:nlin+1,:)
    rec%MatPostProc(nlin+2,:)      = xix(0,0)          * rec%MatPostProc(nlin+2,:)

    deallocate(xix, xiz, etax, etaz)

end subroutine build_MatPostProc


! ###########################################################
!>
!! \brief This function gives the position (line or column) in the matrix MatPostProc
!! corresponfing to a couple of indexes (i,j)
!<
function Ind(i,j,d,ngx,ngz) result(index)

    implicit none
    integer, intent(in) :: i,j,d,ngx,ngz
    integer             :: index

    index = d*(ngx+1)*(ngz+1) + i*(ngz+1) + j+1

end function Ind


! ###########################################################
!>
!! \brief This function computes a kronecker Delta for two integers
!! corresponfing to a couple of indexes (i,j)
!<
function Delta(i,j) result(d)

    implicit none
    integer, intent(in) :: i,j
    integer             :: d

    if (i==j) then
        d=1
    else
        d=0
    endif

end function Delta


end module treceivers
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
