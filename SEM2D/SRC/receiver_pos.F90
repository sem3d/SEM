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
            call check_receiver_on_vertex(Tdomain,nrec)
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

    integer :: ir, nr, i,j, ngllx, ngllz, nsta, ncache, ind
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

            do j = 0,ngllz-1
                do i =0,ngllx -1
                    dum0 = dum0 + Tdomain%sReceiver(ir)%Interp_Coeff(i,j) * Field(i,j,0)
                    dum1 = dum1 + Tdomain%sReceiver(ir)%Interp_Coeff(i,j) * Field(i,j,1)
                enddo
            enddo

            if (Tdomain%sReceiver(ir)%on_vertex .AND. Tdomain%specel(nr)%type_dg == GALERKIN_HDG_RP) &
                call build_vertex_vhat_for_receiver(Tdomain, ir, dum0, dum1)

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


subroutine check_receiver_on_vertex(Tdomain,nrec)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: nrec
    type(element), pointer      :: Elem
    integer, dimension(:), allocatable :: near_faces_tmp
    integer :: nv, n, i
    real    :: tol
    tol = 1.E-10 ; nv = -1 ; i=0
    Elem=> Tdomain%specel(Tdomain%sReceiver(nrec)%nr)
    Tdomain%sReceiver(nrec)%on_vertex = .false.

    if (abs(Tdomain%sReceiver(nrec)%eta + 1.) .LE. tol) then
        if (abs(Tdomain%sReceiver(nrec)%xi + 1.) .LE. tol) then
            nv = Elem%Near_Vertex(0)
            write (*,*) "Receiver relocated on Vertex : ", nv
        else if (abs(Tdomain%sReceiver(nrec)%xi - 1.) .LE. tol) then
            nv = Elem%Near_Vertex(1)
            write (*,*) "Receiver relocated on Vertex : ", nv
        endif
    else if (abs(Tdomain%sReceiver(nrec)%eta - 1.) .LE. tol) then
        if (abs(Tdomain%sReceiver(nrec)%xi + 1.) .LE. tol) then
            nv = Elem%Near_Vertex(3)
            write (*,*) "Receiver relocated on Vertex : ", nv
        else if (abs(Tdomain%sReceiver(nrec)%xi - 1.) .LE. tol) then
            nv = Elem%Near_Vertex(2)
            write (*,*) "Receiver relocated on Vertex : ", nv
        endif
    endif

    if (nv .GE. 0) then
        Tdomain%sReceiver(nrec)%on_vertex = .true.
        write (*,*) "If HDG : Vertex-Projection mode activated for receiver : ", nrec
        Tdomain%sReceiver(nrec)%Nv = nv
        allocate (near_faces_tmp(0:20))
        near_faces_tmp (:) = -1
        do n = 0,Tdomain%n_face-1
            if (nv == Tdomain%sFace(n)%Near_Vertex(0) .OR. nv == Tdomain%sFace(n)%Near_Vertex(1)) then
                near_faces_tmp(i) = n
                i = i+1
            endif
        enddo
        allocate (Tdomain%sReceiver(nrec)%near_faces(0:i-1))
        Tdomain%sReceiver(nrec)%near_faces(0:i-1) = near_faces_tmp(0:i-1)
        deallocate(near_faces_tmp)
    endif

end subroutine check_receiver_on_vertex


subroutine build_vertex_vhat_for_receiver(Tdomain, nrec, dum0, dum1)

    use sdomain
    implicit none

    type(domain), intent(inout) :: Tdomain
    integer, intent(in)         :: nrec
    real, intent(inout)         :: dum0, dum1
    integer :: n, nv, nface, ngll

    nv = Tdomain%sReceiver(nrec)%Nv
    Tdomain%sVertex(Nv)%V0 = 0.
    dum0 = 0 ; dum1 = 0.

    do n=0,size(Tdomain%sReceiver(nrec)%near_faces)-1
        nface = Tdomain%sReceiver(nrec)%near_faces(n)
        ngll  =  Tdomain%sFace(nface)%ngll
        if (nv == Tdomain%sFace(nface)%Near_Vertex(0)) then
            Tdomain%sVertex(Nv)%V0 = Tdomain%sVertex(Nv)%V0 &
                + Tdomain%sFace(nface)%Veloc(0,:) * Tdomain%sFace(nface)%Coeff_Integr_ends(0)
        else
            Tdomain%sVertex(Nv)%V0 = Tdomain%sVertex(Nv)%V0 &
                + Tdomain%sFace(nface)%Veloc(ngll-1,:) * Tdomain%sFace(nface)%Coeff_Integr_ends(1)
        endif
    enddo
    Tdomain%sVertex(Nv)%V0 = Tdomain%sVertex(Nv)%V0 * Tdomain%sVertex(Nv)%CoeffAssem
    dum0 = Tdomain%sVertex(Nv)%V0(0)
    dum1 = Tdomain%sVertex(Nv)%V0(1)

end subroutine build_vertex_vhat_for_receiver


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
    real, dimension (:),   allocatable :: GLLcx, GLLwx, GLLpolx, GLLcz, GLLwz, GLLpolz
    real, dimension (:,:), allocatable :: coord
    real, dimension (0:1,0:1) :: LocInvGrad
    real    :: xi, eta, Jac
    integer :: nrec, ngx, ngz, nr, i_aus, mat, i, j

    do nrec = 0,Tdomain%n_receivers-1
        ! Matrices Allocation used for Post-Process
        nr = Tdomain%sReceiver(nrec)%nr
        if (Tdomain%specel(nr)%Acoustic) cycle
        ngx = Tdomain%specel(nr)%ngllx
        ngz = Tdomain%specel(nr)%ngllz
        mat = Tdomain%specel(nr)%mat_index

        ! Remark : Element's interpolation order is (ngx-1,ngz-1)
        ! when in post-traitment, it becomes (ngx,ngz).
        allocate(Tdomain%sReceiver(nrec)%InvGrad(0:ngx,0:ngz,0:1,0:1))
        allocate(Tdomain%sReceiver(nrec)%JacobWhei(0:ngx,0:ngz))
        allocate(Tdomain%sReceiver(nrec)%ReinterpXi (0:ngx-1,0:ngx))
        allocate(Tdomain%sReceiver(nrec)%ReinterpEta(0:ngz-1,0:ngz))
        allocate(Tdomain%sReceiver(nrec)%hprimex(0:ngx,0:ngx))
        allocate(Tdomain%sReceiver(nrec)%hprimez(0:ngz,0:ngz))
        allocate(Tdomain%sReceiver(nrec)%MatPostProc(0:(ngx+ngz+1),0:(ngx+ngz-1)))

        ! Computation of GLL positions, local wheights, and Lagrange polynomial derivatives
        allocate (GLLcx(0:ngx)); allocate(GLLwz(0:ngz)); allocate (GLLpolx(0:ngx))
        allocate (GLLcz(0:ngx)); allocate(GLLwz(0:ngz)); allocate (GLLpolz(0:ngz))

        call zelegl (ngx,GLLcx,GLLpolx)
        call welegl (ngx,GLLcx,GLLpolx,GLLwx)
        call dmlegl (ngx,ngx,GLLcx,GLLpolx,Tdomain%sReceiver(nrec)%hprimex)
        Tdomain%sReceiver(nrec)%hprimex = TRANSPOSE(Tdomain%sReceiver(nrec)%hprimex)

        call zelegl (ngz,GLLcz,GLLpolz)
        call welegl (ngz,GLLcz,GLLpolz,GLLwz)
        call dmlegl (ngz,ngz,GLLcz,GLLpolz,Tdomain%sReceiver(nrec)%hprimez)
        Tdomain%sReceiver(nrec)%hprimez = TRANSPOSE(Tdomain%sReceiver(nrec)%hprimez)

        ! Computation of Jacobian and InvGrad for the source element
        if (Tdomain%n_nodes == 4) then
            allocate (coord(0:1,0:3))
            do i=0,3
                i_aus = Tdomain%specel(nr)%Control_Nodes(i)
                coord(0,i) = Tdomain%Coord_Nodes(0,i_aus)
                coord(1,i) = Tdomain%Coord_Nodes(1,i_aus)
            end do
            do j = 0,ngz
                eta = GLLcz(j)
                do i = 0,ngx
                    xi = GLLcx(i)
                    ! Computation of the derivative matrix
                    call shape4_local2jacob(coord, xi, eta, LocInvGrad)
                    call invert2(LocInvGrad, Jac)
                    Tdomain%sReceiver(nrec)%InvGrad(i,j,0:1,0:1) = LocInvGrad(0:1,0:1)
                    Tdomain%sReceiver(nrec)%JacobWhei(i,j) = Jac*GLLwx(i)*GLLwz(j)
                enddo
            enddo

        elseif (Tdomain%n_nodes == 8) then
            allocate (coord(0:1,0:7))
            do i=0,7
                i_aus = Tdomain%specel(nr)%Control_Nodes(i)
                coord(0,i) = Tdomain%Coord_Nodes(0,i_aus)
                coord(1,i) = Tdomain%Coord_Nodes(1,i_aus)
            end do
            do j = 0,ngz
                eta =   Tdomain%sSubdomain(mat)%GLLcz (j)
                do i = 0,ngx
                    xi = Tdomain%sSubdomain(mat)%GLLcx (i)
                    ! Computation of the derivative matrix
                    call shape8_local2jacob(coord, xi, eta, LocInvGrad)
                    call invert2(LocInvGrad, Jac)
                    Tdomain%sReceiver(nrec)%InvGrad (i,j,0:1,0:1) = LocInvGrad(0:1,0:1)
                    Tdomain%sReceiver(nrec)%JacobWhei (i,j) = Jac*GLLwx(i)*GLLwz(j)
                enddo
            enddo
        endif

        deallocate(GLLcx,GLLwx,GLLpolx,GLLcz,GLLwz,GLLpolz,coord)

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
    allocate (xix (0:ngx,0:ngz)) ; allocate (xiz (0:ngx,0:ngz))
    allocate (etax(0:ngx,0:ngz)) ; allocate (etaz(0:ngx,0:ngz))
    xix  = rec%InvGrad(:,:,0,0)
    xiz  = rec%InvGrad(:,:,1,0)
    etax = rec%InvGrad(:,:,0,1)
    etaz = rec%InvGrad(:,:,1,1)

    ! Compute entries of matrix MatPostProc
    ! Here, (i,j) indexes corresponds to test function w_ij which gives line number of MatPostProc
    ! indexes (k,l) corresponds to the nodal value of the postprocessed velocity v^* and gives
    ! the column number of MatPostProc. Finally, indexes (r,s) correspond to the quadrature points.

    ! Loop over test functions
    do i = 0,ngx
        do j = 0,ngz
            nlin = Ind(i,j,0,ngx,ngz)
            ! Loop over nodal values of v^*
            do k = 0,ngx
                do l = 0,ngz
                    vx = 0. ; vz = 0.
                    ! Loop over quadrature points
                    do r = 0,ngx
                        do s = 0,ngz
                        vx = vx +((rec%hprimex(k,r)*Delta(l,s)*xix(r,s) + Delta(k,r)*rec%hprimez(l,s)*etax(r,s)) &
                                 *(rec%hprimex(i,r)*Delta(j,s)*xix(r,s) + Delta(i,r)*rec%hprimez(j,s)*etax(r,s)) &
                                 +(rec%hprimex(k,r)*Delta(l,s)*xiz(r,s) + Delta(k,r)*rec%hprimez(l,s)*etaz(r,s)) &
                                 *(rec%hprimex(i,r)*Delta(j,s)*xiz(r,s) + Delta(i,r)*rec%hprimez(j,s)*etaz(r,s)) * 0.5) &
                                 *rec%JacobWhei(r,s)
                        vz = vz +((rec%hprimex(k,r)*Delta(l,s)*xix(r,s) + Delta(k,r)*rec%hprimez(l,s)*etax(r,s)) &
                                 *(rec%hprimex(i,r)*Delta(j,s)*xiz(r,s) + Delta(i,r)*rec%hprimez(j,s)*etaz(r,s)) * 0.5) &
                                 *rec%JacobWhei(r,s)
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
                    do r = 0,ngx
                        do s = 0,ngz
                        vx = vx +((rec%hprimex(k,r)*Delta(l,s)*xiz(r,s) + Delta(k,r)*rec%hprimez(l,s)*etaz(r,s)) &
                                 *(rec%hprimex(i,r)*Delta(j,s)*xix(r,s) + Delta(i,r)*rec%hprimez(j,s)*etax(r,s)) * 0.5) &
                                 *rec%JacobWhei(r,s)
                        vz = vz +((rec%hprimex(k,r)*Delta(l,s)*xiz(r,s) + Delta(k,r)*rec%hprimez(l,s)*etaz(r,s)) &
                                 *(rec%hprimex(i,r)*Delta(j,s)*xix(r,s) + Delta(i,r)*rec%hprimez(j,s)*etax(r,s)) &
                                 +(rec%hprimex(k,r)*Delta(l,s)*xix(r,s) + Delta(k,r)*rec%hprimez(l,s)*etax(r,s)) &
                                 *(rec%hprimex(i,r)*Delta(j,s)*xix(r,s) + Delta(i,r)*rec%hprimez(j,s)*etax(r,s)) * 0.5) &
                                 *rec%JacobWhei(r,s)
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

    index = d*ngx*ngz + i*ngx + j

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
