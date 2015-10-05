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
        ! Find receiver location

        xc = Tdomain%sReceiver(nrec)%xrec
        zc = Tdomain%sReceiver(nrec)%zrec
        nmax = NMAXEL
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

    ! Prepare to store receiver trace

    i = 0
    do nrec = 0, Tdomain%n_receivers-1
        if (Tdomain%sReceiver(nrec)%located_here) i= i + 1
    enddo
    if (i > 0) then
        allocate (Tdomain%Store_Trace(0:1,0:i-1,0:NCAPT_CACHE-1))
    endif
end subroutine ReceiverPosition

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

            call gather_elem_veloc(Tdomain, nr, Field)

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
