!>
!!\file dumptrace.F90
!!\brief Contient la subroutine dump_trace.
!!\author
!!\version 1.0
!!\date 10/03/2009
!!
!<

!>
!! \brief
!!
!! \param type(Domain), intent (IN) Tdomain
!<


subroutine dump_trace (Tdomain)

    use sdomain
    use semdatafiles
    implicit none

    type(Domain), intent (IN) :: Tdomain

    integer :: i,i1,it, nt
    character(Len=MAX_FILE_SIZE) :: fnamef
    real :: rtime

    !dumping the traces
    print *, Tdomain%MPI_var%my_rank
    nt = Tdomain%TimeD%ntimeMax-1
    do i = 0,Tdomain%n_receivers-1
        if (Tdomain%sReceiver(i)%located_here) then
            i1 = i+1

            call semname_dumptrace_tracex(i1,fnamef)
            open (31,file=fnamef, form="formatted")

            call semname_dumptrace_tracez(i1,fnamef)
            open (32,file=fnamef, form="formatted")
            rtime=0.
            do it = 0, nt-1
                rtime = rtime + Tdomain%TimeD%dtmin
                write (31,*) rtime,Tdomain%Store_Trace (0,i,it)
                write (32,*) rtime,Tdomain%Store_Trace (1,i,it)
            enddo
            close (31)
            close (32)
        endif
    enddo
    return
end subroutine dump_trace
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
