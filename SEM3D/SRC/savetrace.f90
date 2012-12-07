subroutine savetrace (Tdomain,rg,k)

use sdomain

implicit none

type (domain), intent (IN):: Tdomain
integer, intent (IN) :: rg, k

integer :: n, ntime, nt1, nt2, ndt2
real :: dt1, dt2
character (len=100) :: fnamef


do n = 0, Tdomain%n_receivers-1

  nt1 = Tdomain%TimeD%ntrace
  dt1 = Tdomain%TimeD%dtmin
  ndt2 = Tdomain%sReceiver(n)%ndt
  nt2 = Tdomain%TimeD%ntrace/ndt2
  dt2 = Tdomain%TimeD%dtmin*ndt2

      if (rg == Tdomain%sReceiver(n)%proc) then
         if ( mod(k,2) == 1 ) then 
            write (fnamef,"(a,I4.4,a,I3.3)") "trace",n,"X",k
         else
            write (fnamef,"(a,I4.4,a,I3.3)") "trace",n,"X",k 
         endif
         open (61,file=fnamef,status="unknown",form="formatted")
         if ( mod(k,2) == 1 ) then 
            write (fnamef,"(a,I4.4,a,I3.3)") "trace",n,"Y",k
         else
            write (fnamef,"(a,I4.4,a,I3.3)") "trace",n,"Y",k
         endif
         open (62,file=fnamef,status="unknown",form="formatted")
         if ( mod(k,2) == 1 ) then 
            write (fnamef,"(a,I4.4,a,I3.3)") "trace",n,"Z",k
         else
            write (fnamef,"(a,I4.4,a,I3.3)") "trace",n,"Z",k
         endif
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
end subroutine
