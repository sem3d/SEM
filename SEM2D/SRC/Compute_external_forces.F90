!>
!!\file Compute_external_forces.F90
!!\brief Computes the external Forces (i.e. the source of the seism)
!!\version 1.0
!!\date 20/11/2013
!! This algorithm contains the part of the previous Newmark routine 
!! which was calculating the external forces.
!<
subroutine Compute_external_forces (Tdomain)

  use sdomain

  implicit none
  type (domain), intent (INOUT) :: Tdomain
  integer  :: n, ns, ncc, ngllx, ngllz, i, j, np

  do n = 0, Tdomain%n_source-1
     do ns =0, Tdomain%sSource(n)%ine-1
        ncc = Tdomain%sSource(n)%Elem(ns)%nr
        ngllx = Tdomain%specel(ncc)%ngllx; ngllz = Tdomain%specel(ncc)%ngllz
        
        if (Tdomain%sSource(n)%i_type_source == 1) then
           do j = 0,ngllz-1
              do i = 0,ngllx-1
                 do np = 0,1
                    Tdomain%specel(ncc)%Forces(i,j,np) = Tdomain%specel(ncc)%Forces(i,j,np) +   &
                         CompSource (Tdomain%sSource(n),Tdomain%TimeD%rtime,np)*  Tdomain%sSource(n)%Elem(ns)%ExtForce(i,j)
                 enddo
              enddo
           enddo
           
        else if (Tdomain%sSource(n)%i_type_source == 2 ) then
           do j = 0,ngllz-1
              do i = 0,ngllx-1
                 do np = 0,1
                    Tdomain%specel(ncc)%Forces(i,j,np) =Tdomain%specel(ncc)%Forces(i,j,np) +   &
                         CompSource(Tdomain%sSource(n), Tdomain%TimeD%rtime,np)* Tdomain%sSource(n)%Elem(ns)%Explosion(i,j,np)
                 enddo
              enddo
           enddo
        endif
     enddo
  enddo

end subroutine Compute_external_forces
