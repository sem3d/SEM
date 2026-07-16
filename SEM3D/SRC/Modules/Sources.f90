!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
!>
!!\file Sources.f90
!!\brief Assure la gestion des sources.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

module ssources

    use constants
    use stf_helpers
    implicit none

    type :: Source
       ! GENERAL PARAMETERS (see SourcePosition.f90-Source.f90)
       integer                     :: proc                      ! source belonging processor
       integer                     :: elem                      ! source belonging elem
       integer                     :: i_type_source             ! source type (solid pulse-moment-fluid pulse)
       ! SPATIAL PARAMETERS
       real(fpp), dimension(0:2)   :: refcoord                  ! local coordinates (master element)
       real(fpp), dimension(0:2)   :: dir                       ! source direction
       real(fpp)                   :: Xsource, YSource, Zsource ! source coordinates
       real(fpp)                   :: amplitude_factor          ! amplitude factor
       integer                     :: ind_i, ind_j              ! position in extended source gr
       ! TIME PARAMETERS
       real(fpp)                   :: ts                        ! time shift
       integer                     :: i_time_function           ! source type (ricker, gabor, etc.)
       integer                     :: Nt                        ! Nb of time samples (for a source file)
       real(fpp)                   :: time_integral             ! running int(f dt) for pressure source (i_type_source==7)

       ! MOMENT SOURCE
       real(fpp), dimension(0:2,0:2)             :: InvGrad     ! Inverse Jacobian
       real(fpp), dimension(0:2,0:2)             :: Moment      ! Moment tensor
       real(fpp), dimension (:,:,:,:), pointer   :: coeff       ! weight coefficient
       ! SOURCE FROM EXTERNAL FILE
       character(len = 30)                  :: time_file        ! file name of external source

       real(fpp) :: tau_b,cutoff_freq,Q,X,Y,L,v,d,a
       !   ajout de parametres pour definir Gabor signal source de type 4
       !   ajout de gamma et ts
       real(fpp) ::  gamma

       real(fpp), dimension(0:3) :: fh

       real(fpp), dimension(:), pointer :: timefunc
       real(fpp), dimension(:), pointer :: ampli, time
       real(fpp), allocatable, dimension(:,:,:,:) :: ExtForce
    end type Source

contains

    !>
    !! \brief
    !!
    !! \param type (source) Sour
    !! \param real time
    !<
    real(fpp) function CompSource (Sour,time,ntime)
        implicit none
        type (source), intent(in) :: Sour
        integer, intent(in) :: ntime
        real(fpp), intent(in) :: time

        CompSource = 0d0

        select case (Sour%i_time_function)
        case (1)
            CompSource = Gaussian_3D (time, Sour%ts, Sour%tau_b)
        case (2)
            CompSource = Ricker_3D (time,Sour%tau_b,Sour%cutoff_freq)
        case (3)
            CompSource = Sour%timefunc(ntime)
        case (4)
            CompSource = Gabor (time,Sour%tau_b,Sour%cutoff_freq,Sour%gamma,Sour%ts)
        case (5)
            CompSource = Source_File (time,Sour)
            !CompSource = Source_File2 (ntime,Sour)
        case (6)
            CompSource = Source_Spice_Bench_math(time, Sour%ts, Sour%cutoff_freq, Sour%gamma)
        case (7)
            CompSource = Source_sinewave_math(time, Sour%ts, Sour%cutoff_freq)
        case (8)
            CompSource = Source_square_math(time, Sour%ts, Sour%tau_b, Sour%gamma)
        case (9)
            CompSource = Source_tanh_math(time, Sour%ts, Sour%gamma)
        case (10)
#ifdef CPML
            CompSource = Ricker_3D(time, Sour%tau_b, Sour%cutoff_freq)
#else
            CompSource = Ricker_fl(time, Sour%tau_b, Sour%cutoff_freq)
#endif
        case (11)
            CompSource = Triangle(time, Sour%tau_b)
        case (12)
            CompSource = HSF(time, Sour%tau_b)
        case(13)
            CompSource = DM(time,Sour%tau_b,Sour%Q,Sour%X,Sour%Y,Sour%L,Sour%v,Sour%d,Sour%a)
        case (16)
            CompSource = Ormsby(time, Sour%tau_b, Sour%fh)
        case(15)
            !Extended source - modified 22/03/17 by Filippo and Elif
            CompSource = Source_File (time,Sour)
        end select
        CompSource = CompSource*Sour%amplitude_factor
        

        return
    end function CompSource





    real(fpp) function Source_File(tt,Sour)
        implicit none
        type(source), intent(in)  :: Sour
        real(fpp), intent(in)     :: tt
        integer :: i

        if(tt < Sour%time(0)) then
            Source_File = 0.
        else if(tt >= Sour%time(Sour%Nt-1)) then
            Source_File = Sour%ampli(Sour%Nt-1)
        else
            i = 1
            do while(tt > Sour%time(i))
                i = i + 1
            end do
            Source_File = Sour%ampli(i-1) + (tt-Sour%time(i-1))*(Sour%ampli(i)-Sour%ampli(i-1))/(Sour%time(i)-Sour%time(i-1))
        endif
        return

    end function Source_File


    real(fpp) function Source_File2(ntime,Sour)
        implicit none
        type(source), intent(in) :: Sour
        integer, intent(in)      :: ntime

        if(ntime < Sour%Nt )then
            Source_File2 = Sour%ampli(ntime)
        else
            Source_File2 = Sour%ampli(Sour%Nt-1)
        endif

    end function Source_File2


    subroutine read_source_file(Sour)
        implicit none
        !- lecture directe d'un fichier temps-amplitude pour la source
        type(Source), intent(inout)   :: Sour
        integer                       :: nb_time_step
        integer                       :: i
        real(fpp)                          :: tr, trr

        i = 0 ; nb_time_step = 0
        ! count

        open(10,file=Sour%time_file,action="read",status="old")
        do
            read(10,*,end=100) tr, trr
            i = i + 1
        end do
100     close(10)
        ! nombre de donnees en entree pour la source
        nb_time_step = i

        allocate(Sour%time(0:nb_time_step-1),Sour%ampli(0:nb_time_step-1))
        Sour%Nt = nb_time_step

        open(10,file=Sour%time_file,action="read",status="old")
        do i=0, Sour%Nt-1
            read(10,*,end=101) Sour%time(i),Sour%ampli(i)
        end do
101     close(10)

    end subroutine read_source_file

    subroutine read_source_file_h5(Sour)
        use sem_hdf5
        implicit none
        !- lecture directe d'un fichier temps-amplitude pour la source
        type(Source), intent(inout)      :: Sour
        real(fpp), allocatable, dimension(:,:,:) :: data
        real(fpp), allocatable, dimension(:)     :: dataT
        integer, dimension(0:2)          :: imin, imax
        integer                          :: hdferr
        integer(HID_T)                   :: fid

        ! Creation des coordonnes du chunck a lire
        imin(0) = 0 ;         imin(1) = Sour%ind_j ; imin(2) = Sour%ind_i
        imax(0) = Sour%Nt-1 ; imax(1) = Sour%ind_j ; imax(2) = Sour%ind_i

        ! Ouverture ficher slip history HDF5
        call h5fopen_f(Sour%time_file, H5F_ACC_RDONLY_F, fid, hdferr)

        ! Lecture du dataset
        ! MODIFICATION 1
        call read_subset_3d_real(fid, 'moment', imin, imax, data)
        call read_dset_1d_real(fid, 'time', dataT)

        
        ! Historique de l'amplitude du Slip
        allocate(Sour%ampli(0:Sour%Nt-1))
        allocate(Sour%time(0:Sour%Nt-1))
        Sour%ampli(:) = data(:,imin(1),imin(2))
        Sour%time(:)  = dataT(:)

        call h5fclose_f(fid, hdferr)

        deallocate(data, dataT)

    end subroutine read_source_file_h5


    !   modif pour benchmark can2
    !-------------------------------------------------



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
