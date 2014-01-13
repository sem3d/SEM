!>
!!\file Compute_external_forces.F90
!!\brief Computes the external Forces (i.e. the source of the seism)
!!\version 1.0
!!\date 20/11/2013
!! This algorithm contains the part of the previous Newmark routine 
!! which was calculating the external forces.
!<
subroutine capteurs_veloc (Tdomain,timelocal,ntime,nelem)

  use sdomain

  implicit none
  type (domain), intent (INOUT) :: Tdomain
  real, intent (IN)          :: timelocal
  integer, intent (IN)       :: ntime, nelem
  integer                    :: ngllx, ngllz

  if (.NOT.Tdomain%openfilescapt) then
     open (51,file="capteur1Exx",status="replace",form="formatted")
     open (52,file="capteur1Ezz",status="replace",form="formatted")
     open (53,file="capteur1Exz",status="replace",form="formatted")
     open (54,file="capteur1Vx",status="replace",form="formatted")
     open (55,file="capteur1Vz",status="replace",form="formatted")
     open (56,file="capteur2Exx",status="replace",form="formatted")
     open (57,file="capteur2Ezz",status="replace",form="formatted")
     open (58,file="capteur2Exz",status="replace",form="formatted")
     open (59,file="capteur2Vx",status="replace",form="formatted")
     open (60,file="capteur2Vz",status="replace",form="formatted")
     Tdomain%openfilescapt = .true.
  endif
  
  if (Tdomain%openfilescapt) then
     ngllx = Tdomain%specel(nelem)%ngllx
     ngllz = Tdomain%specel(nelem)%ngllz

     write(51,*) timelocal, Tdomain%specel(nelem)%Forces(ngllx-1,2,0)
     write(52,*) timelocal, Tdomain%specel(nelem)%Forces(ngllx-1,2,1)
     write(53,*) timelocal, Tdomain%specel(nelem)%Forces(ngllx-1,2,2)
     write(54,*) timelocal, Tdomain%specel(nelem)%Forces(ngllx-1,2,3)
     write(55,*) timelocal, Tdomain%specel(nelem)%Forces(ngllx-1,2,4)
     write(56,*) timelocal, Tdomain%specel(nelem+1)%Forces(0,2,0)
     write(57,*) timelocal, Tdomain%specel(nelem+1)%Forces(0,2,1)
     write(58,*) timelocal, Tdomain%specel(nelem+1)%Forces(0,2,2)
     write(59,*) timelocal, Tdomain%specel(nelem+1)%Forces(0,2,3)
     write(60,*) timelocal, Tdomain%specel(nelem+1)%Forces(0,2,4)

     !write(51,*) timelocal, Tdomain%specel(nelem)%Veloc(1,1,0)
     !write(52,*) timelocal, Tdomain%specel(nelem)%Veloc(1,1,1)
     !write(53,*) timelocal, Tdomain%specel(nelem)%Veloc(3,3,0)
     !write(54,*) timelocal, Tdomain%specel(nelem)%Veloc(3,3,1)
  endif

end subroutine capteurs_veloc
