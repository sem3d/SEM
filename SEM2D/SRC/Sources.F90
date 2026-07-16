!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Sources.F90
!!\brief Assure la gestion des sources.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssources

    use constants
    use stf_helpers
    use semdatafiles, only : MAX_FILE_SIZE

    type :: elem_source
       integer :: nr
       real(fpp) :: eta,xi
       real(fpp) :: invE,nu
       real(fpp), dimension (0:1,0:1) :: Scoeff
       real(fpp), dimension (:,:,:), pointer :: ExtForce
    end type elem_source

    type :: Source
       integer :: i_type_source, i_time_function,ine
       real(fpp), dimension(2) :: dir
       real(fpp), dimension (0:1,0:1) :: moment
       real(fpp) :: Xsource,Zsource, tau_b,cutoff_freq,amplitude,sigma
       type(elem_source), dimension(:), pointer :: Elem
       logical :: located_here
       ! Running time integral int(f dt) for the pressure source (i_type_source=7):
       ! p=-VelPhi and the source enters at the phi-acceleration level, so injecting the
       ! integral makes the pressure equal f(t) (cf. SEM3D Newmark.f90 type-7). Unused by
       ! the other source types.
       real(fpp) :: time_integral = 0._fpp
       ! SOURCE FROM EXTERNAL FILE (i_time_function=5): time-amplitude table, 2 columns
       ! ("t ampli") and one sample per line, read by read_source_file.
       character(len=MAX_FILE_SIZE) :: time_file
       integer :: Nt
       real(fpp), dimension(:), pointer :: ampli, time
       real(fpp), dimension(0:3) :: fh
    end type Source

contains

    !>
    !! \fn function CompSource (Sour,time,np)
    !! \brief
    !!
    !! \param type (source) Sour
    !! \param integer np
    !! \param real time
    !<
    real(fpp) function CompSource (Sour,time)

        type (source) :: Sour
        real(fpp) :: time

        CompSource = 0.
        select case (Sour%i_time_function)
        case (1)
            CompSource = Gaussian_2D (time,Sour%tau_b,Sour%cutoff_freq)
        case (2)
            CompSource = Ricker_2D (time,Sour%tau_b,Sour%cutoff_freq)
        case (3)
            CompSource = 1
        case (5)
            CompSource = Source_File (time,Sour)
        case (16)
            CompSource = Ormsby (time,Sour%tau_b,Sour%fh)
        end select
        CompSource = Sour%amplitude*CompSource

        return
    end function CompSource

    !>
    !! \brief Reads the time-amplitude table of a func=file source (i_time_function=5).
    !! Two columns ("t ampli"), one sample per line. Port of SEM3D Modules/Sources.f90.
    !<
    subroutine read_source_file(Sour)

        type(Source), intent(inout) :: Sour
        integer   :: nb_time_step, i
        real(fpp) :: tr, trr

        i = 0 ; nb_time_step = 0

        ! count
        open(10,file=Sour%time_file,action="read",status="old")
        do
            read(10,*,end=100) tr, trr
            i = i + 1
        end do
100     close(10)
        nb_time_step = i

        if (nb_time_step < 2) then
            write(*,*) "ERROR: source time_file '", trim(Sour%time_file), "' yielded ", &
                       nb_time_step, " sample(s)."
            write(*,*) "  Expected 2 columns (t ampli) and one sample per line."
            stop 1
        end if

        allocate(Sour%time(0:nb_time_step-1),Sour%ampli(0:nb_time_step-1))
        Sour%Nt = nb_time_step

        open(10,file=Sour%time_file,action="read",status="old")
        do i=0, Sour%Nt-1
            read(10,*,end=101) Sour%time(i),Sour%ampli(i)
        end do
101     close(10)

    end subroutine read_source_file

    !>
    !! \brief Linear interpolation in the time-amplitude table: zero before the first sample,
    !! held at the last amplitude after the last one. Port of SEM3D Modules/Sources.f90.
    !<
    real(fpp) function Source_File(tt,Sour)

        type(Source), intent(in) :: Sour
        real(fpp), intent(in)    :: tt
        integer :: i

        if (tt < Sour%time(0)) then
            Source_File = 0.
        else if (tt >= Sour%time(Sour%Nt-1)) then
            Source_File = Sour%ampli(Sour%Nt-1)
        else
            i = 1
            do while (tt > Sour%time(i))
                i = i + 1
            end do
            Source_File = Sour%ampli(i-1) + (tt-Sour%time(i-1)) &
                          * (Sour%ampli(i)-Sour%ampli(i-1)) / (Sour%time(i)-Sour%time(i-1))
        endif
        return

    end function Source_File



end module ssources

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
