!>
!!\file Compute_external_forces.F90
!!\brief Computes the external Forces (i.e. the source of the seism)
!!\version 1.0
!!\date 20/11/2013
!! This algorithm contains the part of the previous Newmark routine 
!! which was calculating the external forces.
!<
subroutine Compute_external_forces (Tdomain,timelocal)

  use sdomain

  implicit none
  type (domain), intent (INOUT) :: Tdomain
  real, intent (INOUT)          :: timelocal
  real, dimension(0:1)          :: Fext
  integer  :: n, ns, ncc, ngllx, ngllz, i, j, np, nDG

  do n = 0, Tdomain%n_source-1
     do ns =0, Tdomain%sSource(n)%ine-1
        ncc = Tdomain%sSource(n)%Elem(ns)%nr
        ngllx = Tdomain%specel(ncc)%ngllx; ngllz = Tdomain%specel(ncc)%ngllz
        if ( Tdomain%specel(ncc)%Type_DG == 2) then
            nDG = 0
        else ! In DG case, the forces has to be put in Forces(:,:,3:4)
             ! instead of Forces(:,:,0:1) ... nDG is then computed here to shift Forces
            nDG = 3
        endif
        
        if (Tdomain%sSource(n)%i_type_source == 1) then
           do j = 0,ngllz-1
              do i = 0,ngllx-1
                 do np = 0,1
                    Fext(np) =  CompSource (Tdomain%sSource(n),timelocal,np) * Tdomain%sSource(n)%Elem(ns)%ExtForce(i,j)
                    Tdomain%specel(ncc)%Forces(i,j,np+nDG) = Tdomain%specel(ncc)%Forces(i,j,np+nDG) + Fext(np)
                 enddo
                 ! Sortie pour la force externe
                 if ((Fext(0).NE.0.) .OR. (Fext(1).NE.0.)) write(90,*) timelocal, Fext(0)
              enddo
           enddo
           
        else if (Tdomain%sSource(n)%i_type_source == 2 ) then
           do j = 0,ngllz-1
              do i = 0,ngllx-1
                 do np = 0,1
                    Fext(np) = CompSource(Tdomain%sSource(n),timelocal,np) * Tdomain%sSource(n)%Elem(ns)%Explosion(i,j,np)
                    Tdomain%specel(ncc)%Forces(i,j,np+nDG) =Tdomain%specel(ncc)%Forces(i,j,np+nDG) + Fext(np)
                 enddo
                 ! Sortie pour la force externe
                 if ((Fext(0).NE.0.) .OR. (Fext(1).NE.0.)) write(90,*) timelocal, Fext(:)
              enddo
           enddo
        endif
     enddo
  enddo

end subroutine Compute_external_forces
