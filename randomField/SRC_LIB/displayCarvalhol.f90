module displayCarvalhol

    use mpi
    use constants_RF
    use write_Log_File
    !All display routines
    interface DispCarvalhol
       module procedure Disp1Ddble,   &
           Disp2Ddble,   &
           DispScalDble, &
           Disp1Dint,    &
           Disp2Dint,    &
           Disp1Dchar,   &
           Disp2Dchar,   &
           Disp1Dbool,   &
           Disp2Dbool
    end interface DispCarvalhol

contains

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine Ordering_MPI_Start(comm)

        implicit none

        !INPUT
        integer, optional, intent(in) :: comm

        !LOCAL VARIABLES
        integer, dimension( MPI_STATUS_SIZE ) :: statut
        integer :: rang, code, nb_procs, id = 15, passageCount = 0, effectComm
        integer :: sender

        if(present(comm))       effectComm = comm
        if(.not. present(comm)) effectComm = MPI_COMM_WORLD

        call MPI_COMM_RANK(effectComm, rang, code)
        call MPI_COMM_SIZE(effectComm, nb_procs, code)


        if(rang == 0) then
            sender = nb_procs - 1
        else
            sender = rang - 1
        end if

        if (passageCount == 0 .and. rang == 0) then
        else
            call MPI_RECV (passageCount, 1, MPI_INTEGER, &
                sender, id, effectComm ,statut,code)
        end if

        passageCount = passageCount + 1

        write(*,*) "RANG ", rang, "-------"

    end subroutine Ordering_MPI_Start

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine Ordering_MPI_End(comm)

        implicit none

        !INPUT
        integer, optional, intent(in) :: comm

        !LOCAL VARIABLES
        integer :: rang, code, nb_procs, id = 15, passageCount = 0, effectComm
        integer :: receiver

        if(present(comm))       effectComm = comm
        if(.not. present(comm)) effectComm = MPI_COMM_WORLD

        call MPI_COMM_RANK(effectComm, rang, code)
        call MPI_COMM_SIZE(effectComm, nb_procs, code)

        if(rang == nb_procs - 1) then
            receiver = 0
        else
            receiver = rang + 1
        end if

        call MPI_SEND (passageCount, 1, MPI_INTEGER , &
            receiver, id, effectComm ,code)

        passageCount = passageCount + 1

    end subroutine Ordering_MPI_End

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    subroutine DispScalDble(scalar, title, format, nColumns, mpi, comm, unit_in)
        ! Displays Scalar

        implicit none

        !INPUT
        double precision,            intent(in) :: scalar
        character (len=*), optional, intent(in) :: title, format
        integer,           optional, intent(in) :: nColumns
        logical,           optional, intent(in) :: mpi
        integer          , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        double precision, dimension(:,:), allocatable :: matrix2d
        integer :: effectComm
        integer :: unit = 6

        unit = 6 !Screen
        if(present(unit_in)) unit = unit_in

        if(present(comm))       effectComm = comm
        if(.not. present(comm)) effectComm = MPI_COMM_WORLD

        allocate(matrix2D(1,1));
        matrix2d = scalar;

        if(present(mpi)) then
            call Disp2Ddble(matrix2D, title, format, nColumns, mpi, effectComm, unit_in = unit);
        else
            call Disp2Ddble(matrix2D, title, format, nColumns, unit_in = unit);
        end if

        deallocate(matrix2D);


    end subroutine DispScalDble

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine Disp1Dbool(vector, title, format, nColumns, mpi, comm, unit_in)
        ! Displays 1D Matrix (Vector)

        implicit none

        !INPUT
        logical,             dimension(:)          , intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical                          , optional, intent(in) :: mpi
        integer                          , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        logical, dimension(:,:), allocatable :: matrix2d
        integer :: effectComm
        integer :: unit = 6

        unit = SCREEN !Screen
        if(present(unit_in)) unit = unit_in

        if(present(comm))       effectComm = comm
        if(.not. present(comm)) effectComm = MPI_COMM_WORLD

        allocate(matrix2D(size(vector),1));
        matrix2d(:,1) = vector;

        if(present(mpi)) then
            call Disp2Dbool(matrix2D, title, format, nColumns, mpi, effectComm, unit_in = unit);
        else
            call Disp2Dbool(matrix2D, title, format, nColumns, unit_in = unit);
        end if

        deallocate(matrix2D);

    end subroutine Disp1Dbool

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine Disp1Ddble(vector, title, format, nColumns, mpi, comm, unit_in)
        ! Displays 1D Matrix (Vector)

        implicit none

        !INPUT
        double precision,     dimension(:),          intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical                          , optional, intent(in) :: mpi
        integer                          , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        double precision, dimension(:,:), allocatable :: matrix2d
        integer :: effectComm
        integer :: unit

        unit = SCREEN !Screen
        if(present(unit_in)) unit = unit_in

        if(present(comm))       effectComm = comm
        if(.not. present(comm)) effectComm = MPI_COMM_WORLD

        allocate(matrix2D(size(vector),1));
        matrix2d(:,1) = vector;

        if(present(mpi)) then
            call Disp2Ddble(matrix2D, title, format, nColumns, mpi, effectComm, unit_in = unit);
        else
            call Disp2Ddble(matrix2D, title, format, nColumns, unit_in = unit);
        end if

        deallocate(matrix2D);

    end subroutine Disp1Ddble

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine Disp1Dint(vector, title, format, nColumns, mpi, comm, unit_in)
        ! Displays 1D Matrix (Vector)

        implicit none

        !INPUT
        integer,            dimension(:),            intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical                          , optional, intent(in) :: mpi
        integer                          , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        integer, dimension(:,:), allocatable :: matrix2D
        integer :: effectComm
        integer :: unit

        unit = SCREEN !Screen
        if(present(unit_in)) unit = unit_in

        if(present(comm))       effectComm = comm
        if(.not. present(comm)) effectComm = MPI_COMM_WORLD


        allocate(matrix2D(size(vector),1));
        matrix2d(:,1) = vector;

        if(present(mpi)) then
            call Disp2Dint(matrix2D, title, format, nColumns, mpi, effectComm, unit_in = unit);
        else
            call Disp2Dint(matrix2D, title, format, nColumns, unit_in = unit);
        end if

        deallocate(matrix2D);

    end subroutine Disp1Dint

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine Disp1Dchar(vector, title, format, nColumns, mpi, comm, unit_in)
        ! Displays 1D Matrix (Vector)

        implicit none

        !INPUT
        character(len=*),   dimension(:),            intent(in) :: vector
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical                          , optional, intent(in) :: mpi
        integer                          , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        character(len=30), dimension(:,:), allocatable :: matrix2D
        integer :: effectComm
        integer :: unit

        unit = SCREEN !Screen
        if(present(unit_in)) unit = unit_in

        if(present(comm))       effectComm = comm
        if(.not. present(comm)) effectComm = MPI_COMM_WORLD


        allocate(matrix2D(size(vector),1));
        matrix2d(:,1) = vector;

        if(present(mpi)) then
            call Disp2Dchar(matrix2D, title, format, nColumns, mpi, effectComm, unit_in = unit);
        else
            call Disp2Dchar(matrix2D, title, format, nColumns, unit_in = unit);
        end if
        deallocate(matrix2D);

    end subroutine Disp1Dchar

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine Disp2Ddble(matrix2D, title, format, nColumns, mpi, comm, unit_in)
        ! Displays 2D Matrix, "div" columns at a time

        implicit none

        !INPUT
        double precision, dimension(:, :),           intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical                          , optional, intent(in) :: mpi
        integer                          , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
        character (len=40) :: doubleFmt;
        character (len=10) :: tempFmt;
        integer            :: effectComm
        integer :: unit

        unit = SCREEN !Screen
        if(present(unit_in)) unit = unit_in

        if(unit == log_file_RF_ID .or. unit == SCREEN) then

            if(present(comm))       effectComm = comm
            if(.not. present(comm)) effectComm = MPI_COMM_WORLD


            if(present(mpi)) then
                if(mpi) call Ordering_MPI_Start(effectComm)
            end if

            !write(*,*) "Inside Disp2Ddble"

            if(present(format))                     tempFmt = format;
            if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
            if(.not.present(format))    tempFmt = "F10.5";
            if(.not.present(nColumns))  div = 15;
            write(doubleFmt, fmt="(I3, A)") div, tempFmt;

            write(unit,*) ""
            if(present(title)) write(unit,*) "/////// ", title, " ///////";

            tol = 1E-10;

            write(unit,*) ""
            do k=1, size(matrix2D,2)/div
                write(unit,*)
                write(unit,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
                write(unit,*)
                do j= lbound(matrix2D,1), ubound(matrix2D,1)
                    write(unit,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
                enddo
            enddo

            if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
                write(unit,*)
                write(unit,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
                write(unit,*)

                do j= lbound(matrix2D,1), ubound(matrix2D,1)
                    write(unit,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
                enddo
            end if

            if(present(mpi)) then
                if(mpi) call Ordering_MPI_End(effectComm)
            end if

        end if

    end subroutine Disp2Ddble

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine Disp2Dint(matrix2D, title, format, nColumns, mpi, comm, unit_in)
        ! Displays 2D Matrix, "div" columns at a time

        implicit none

        !INPUT
        integer,           dimension(:, :),          intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical                          , optional, intent(in) :: mpi
        integer                          , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
        character (len=40) :: doubleFmt;
        character (len=10) :: tempFmt;
        integer            :: effectComm
        integer :: unit

        unit = SCREEN !Screen
        if(present(unit_in)) unit = unit_in

        if(unit == log_file_RF_ID .or. unit == SCREEN) then

            if(present(comm))       effectComm = comm
            if(.not. present(comm)) effectComm = MPI_COMM_WORLD


            if(present(mpi)) then
                if(mpi) call Ordering_MPI_Start(effectComm)
            end if

            if(present(format))                     tempFmt = format;
            if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
            if(.not.present(format))    tempFmt = "I6";
            if(.not.present(nColumns))  div = 5;
            write(doubleFmt, fmt="(I3, A)") div, tempFmt;

            write(unit,*) ""
            if(present(title)) write(unit,*) "/////// ", title, " ///////";

            tol = 1E-10;

            write(unit,*) ""
            do k=1, size(matrix2D,2)/div
                write(unit,*)
                write(unit,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
                write(unit,*)
                do j= lbound(matrix2D,1), ubound(matrix2D,1)
                    write(unit,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
                enddo
            enddo

            if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
                write(unit,*)
                write(unit,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
                write(unit,*)

                do j= lbound(matrix2D,1), ubound(matrix2D,1)
                    write(unit,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
                enddo
            end if

            !write(*,*)

            if(present(mpi)) then
                if(mpi) call Ordering_MPI_End(effectComm)
            end if

        end if

    end subroutine Disp2Dint

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine Disp2Dchar(matrix2D, title, format, nColumns, mpi, comm, unit_in)
        ! Displays 2D Matrix, "div" columns at a time

        implicit none

        !INPUT
        character (len=*), dimension(:, :),          intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical                          , optional, intent(in) :: mpi
        integer                          , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
        character (len=40) :: charFmt;
        character (len=10) :: tempFmt;
        integer            :: effectComm
        integer :: unit

        unit = SCREEN !Screen
        if(present(unit_in)) unit = unit_in

        if(unit == log_file_RF_ID .or. unit == SCREEN) then

            if(present(comm))       effectComm = comm
            if(.not. present(comm)) effectComm = MPI_COMM_WORLD

            if(present(mpi)) then
                if(mpi) call Ordering_MPI_Start(effectComm)
            end if

            if(present(format))                     tempFmt = format;
            if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
            if(.not.present(format))    tempFmt = "A";
            if(.not.present(nColumns))  div = 5;
            write(charFmt, fmt="(I3, A)") div, tempFmt;

            write(unit,*) ""
            if(present(title)) write(unit,*) "/////// ", title, " ///////";

            tol = 1E-10;

            write(unit,*) ""
            do k=1, size(matrix2D,2)/div
                write(unit,*)
                write(unit,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
                write(unit,*)
                do j= lbound(matrix2D,1), ubound(matrix2D,1)
                    write(unit,fmt="(I4, A, ("//charFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
                enddo
            enddo

            if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
                write(unit,*)
                write(unit,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
                write(unit,*)

                do j= lbound(matrix2D,1), ubound(matrix2D,1)
                    write(unit,fmt="(I4, A, ("//charFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
                enddo
            end if

            !write(*,*)

            if(present(mpi)) then
                if(mpi) call Ordering_MPI_End(effectComm)
            end if
        end if

    end subroutine Disp2Dchar

    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine Disp2Dbool(matrix2D, title, format, nColumns, mpi, comm, unit_in)
        ! Displays 2D Matrix, "div" columns at a time

        implicit none

        !INPUT
        logical, dimension(:, :),          intent(in) :: matrix2D
        character (len=*),                 optional, intent(in) :: title, format
        integer,                           optional, intent(in) :: nColumns
        logical                          , optional, intent(in) :: mpi
        integer                          , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        integer            :: k, j, div;
        double precision   :: tol;
        character (len=40) :: doubleFmt;
        character (len=10) :: tempFmt;
        integer            :: effectComm
        integer :: unit

        unit = SCREEN !Screen
        if(present(unit_in)) unit = unit_in

        if(unit == log_file_RF_ID .or. unit == SCREEN) then

            if(present(comm))       effectComm = comm
            if(.not. present(comm)) effectComm = MPI_COMM_WORLD


            if(present(mpi)) then
                if(mpi) call Ordering_MPI_Start(effectComm)
            end if

            !write(*,*) "Inside Disp2Ddble"

            if(present(format))                     tempFmt = format;
            if(present(nColumns).and.nColumns.gt.0)     div = nColumns;
            if(.not.present(format))    tempFmt = "L";
            if(.not.present(nColumns))  div = 15;
            write(doubleFmt, fmt="(I3, A)") div, tempFmt;

            write(unit,*) ""
            if(present(title)) write(unit,*) "/////// ", title, " ///////";

            tol = 1E-10;

            write(unit,*) ""
            do k=1, size(matrix2D,2)/div
                !write(*,*)
                write(unit,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", k*div ;
                write(unit,*)
                do j= lbound(matrix2D,1), ubound(matrix2D,1)
                    write(unit,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:k*div)
                enddo
            enddo

            if ((DBLE(size(matrix2D,2))/DBLE(div))-size(matrix2D,2)/div > tol) then
                write(unit,*)
                write(unit,fmt="(A,I3,A,I3)") "Columns", (k-1)*div+1, " to ", ubound(matrix2D,2);
                write(unit,*)

                do j= lbound(matrix2D,1), ubound(matrix2D,1)
                    write(unit,fmt="(I4, A, ("//doubleFmt//"))") j,"->",matrix2D(j,(k-1)*div+1:)
                enddo
            end if

            !write(*,*)

            if(present(mpi)) then
                if(mpi) call Ordering_MPI_End(effectComm)
            end if
        end if

    end subroutine Disp2Dbool


    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------------------

    subroutine Disp3Ddble(matrix, title, format, nColumns, mpi, comm, unit_in)
        ! Displays 1D Matrix (Vector)

        implicit none

        !INPUT
        double precision , dimension(:,:,:),          intent(in) :: matrix
        character (len=*),                  optional, intent(in) :: title, format
        integer,                            optional, intent(in) :: nColumns
        logical                           , optional, intent(in) :: mpi
        integer                           , optional, intent(in) :: comm
        integer          , optional, intent(in) :: unit_in

        !LOCAL VARIABLES
        integer :: i
        integer :: effectComm
        integer :: unit

        unit = SCREEN !Screen
        if(present(unit_in)) unit = unit_in

        if(present(comm))       effectComm = comm
        if(.not. present(comm)) effectComm = MPI_COMM_WORLD

        do i = 1, size(matrix, 3)
            write(unit,*) "---------Slice ", i, "----------"
            if(present(mpi)) then
                call Disp2Ddble(matrix(:,:,i), title, format, nColumns, mpi, effectComm, unit_in = unit);
            else
                call Disp2Ddble(matrix(:,:,i), title, format, nColumns, unit_in = unit);
            end if
        end do

    end subroutine Disp3Ddble

end module displayCarvalhol
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
