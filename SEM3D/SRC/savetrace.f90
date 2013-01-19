!>
!! \file savetrace.f90
!! \brief 
!!
!<

subroutine savetrace (Tdomain,rg,k)

    use sdomain
    use semdatafiles

    implicit none

    type (domain), intent (IN):: Tdomain
    integer, intent (IN) :: rg, k

    integer :: n, ntime, nt1, nt2, ndt2
    real :: dt1, dt2
    character (len=MAX_FILE_SIZE) :: fnamef


    nt1 = Tdomain%TimeD%ntrace
    dt1 = Tdomain%TimeD%dtmin

    do n = 0, Tdomain%n_receivers-1
        ndt2 = Tdomain%sReceiver(n)%ndt
        nt2 = Tdomain%TimeD%ntrace/ndt2
        dt2 = Tdomain%TimeD%dtmin*ndt2
        if (rg == Tdomain%sReceiver(n)%proc) then
            call  semname_savetrace_trace(n,"X",k,fnamef)
            open (61,file=fnamef,status="unknown",form="formatted")
            call semname_savetrace_trace(n,"Y",k,fnamef)
            open (62,file=fnamef,status="unknown",form="formatted")
            call semname_savetrace_trace(n,"Z",k,fnamef)
            open (63,file=fnamef,status="unknown",form="formatted")
            if ( Tdomain%sReceiver(n)%flag == 1 ) then
                do ntime = 0, nt1-1
                    write (61,*) dt1*(ntime+1+k*nt1),Tdomain%sReceiver(n)%StoreTrace(ntime,0)
                    write (62,*) dt1*(ntime+1+k*nt1),Tdomain%sReceiver(n)%StoreTrace(ntime,1)
                    write (63,*) dt1*(ntime+1+k*nt1),Tdomain%sReceiver(n)%StoreTrace(ntime,2)
                enddo
            else if ( Tdomain%sReceiver(n)%flag==2 ) then
                do ntime = 0, nt2-1
                    write (61,*) dt2*(ntime+1+k*nt2),Tdomain%sReceiver(n)%StoreTrace(ntime,0)
                    write (62,*) dt2*(ntime+1+k*nt2),Tdomain%sReceiver(n)%StoreTrace(ntime,1)
                    write (63,*) dt2*(ntime+1+k*nt2),Tdomain%sReceiver(n)%StoreTrace(ntime,2)
                enddo
            endif
        endif
        close(61); close (62); close (63)
    enddo

    return
end subroutine savetrace
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent : !!
