!>
!!\file capteurs_veloc.F90
!!\brief 
!!\version 1.0
!!\date 10/01/2014
!! This algorithm was added because the precedent way of treating 
!! captors in 2D was not working properly
!<
subroutine capteurs_veloc (Tdomain,timelocal)

  use sdomain

  implicit none
  type (domain), intent (INOUT) :: Tdomain
  real, intent (IN)          :: timelocal
  integer                    :: nCaptElem, nCaptFace, ngx, ngz, i, nelem, nface
  
  nCaptElem = 0
  nCaptFace = 0

  do i=0,Tdomain%nCapt-1

     if (Tdomain%type_capteurs(i)==1) then
        nelem = Tdomain%elems_capteurs(nCaptElem,0)
        ngx   = Tdomain%elems_capteurs(nCaptElem,1)
        ngz   = Tdomain%elems_capteurs(nCaptElem,2)
        write(50+i,*) timelocal, Tdomain%specel(nelem)%Veloc(ngx,ngz,0), Tdomain%specel(nelem)%Veloc(ngx,ngz,1)
        nCaptElem = nCaptElem+1
     elseif (Tdomain%type_capteurs(i)==2) then
        nface = Tdomain%faces_capteurs(nCaptFace,0)
        ngx   = Tdomain%faces_capteurs(nCaptFace,1)
        write(50+i,*) timelocal, Tdomain%sface(nface)%Veloc_m(ngx,0), Tdomain%sface(nface)%Veloc_m(ngx,1), &
                                 Tdomain%sface(nface)%Veloc_p(ngx,0), Tdomain%sface(nface)%Veloc_p(ngx,1)
        nCaptFace = nCaptFace+1
     endif
  enddo


end subroutine capteurs_veloc


!>
!!\file capteurs_veloc.F90
!!\brief 
!!\version 1.0
!!\date 15/01/2014
!! This subroutine reads the file capteurs.dat and opens the files
!! which will store the velocities of the points of interest.
!! Therefore, this subroutine finds the elements corresponding to
!! the locations of captors given in capteurs.dat, and it will give
!! the velocities of the central node of the element.
!! BE CAREFULL :  THIS METHOD WORKS FOR MESH COMPOSED OF SQUARES ONLY :
!<
subroutine init_capteurs_veloc (Tdomain)

  use sdomain
  implicit none
  
  type (domain), intent (INOUT)        :: Tdomain
  real, dimension(:,:), allocatable    :: positions
  integer, dimension(:),   allocatable :: tmp_capt_elem, tmp_capt_face
  character(len=50)                    :: fnamef,n2string
  real    :: tol, typical_size, Xcapt, Zcapt
  integer :: n, nv1, nv2, ncaptElem, ncaptFace, nelem, nface, i, ngx, ngz
  logical :: is_inside, is_onface

  ncaptElem = 0 
  ncaptFace = 0
  tol = 1.e-10
  typical_size = sqrt( (Tdomain%coord_nodes(0,1)-Tdomain%coord_nodes(0,0))**2 &
                      +(Tdomain%coord_nodes(1,1)-Tdomain%coord_nodes(1,0))**2 )
  tol = tol*typical_size

  open (15,file="capteurs.dat", status="old", form="formatted")
  write(*,*)"lecture fichier : capteurs.dat"

  read (15,*) Tdomain%nCapt
  read (15,*) 

  allocate (tmp_capt_elem(0:Tdomain%nCapt-1))
  allocate (tmp_capt_face(0:Tdomain%nCapt-1))

  ! Vecteur indiquent le type (element ou face) pour chaque capteur
  allocate (Tdomain%type_capteurs(0:Tdomain%nCapt-1))

  ! Getting positions of the differents captors
  allocate (positions(0:Tdomain%nCapt-1,0:1))
  do n=0,Tdomain%nCapt-1
     read (15,*) positions(n,:)
  enddo
  close(15)

  ! Getting the elemenent where lies each captor
  do n=0,Tdomain%nCapt-1
     is_inside = .false.
     is_onface = .false.
     
     do nelem=0,Tdomain%n_elem-1
        do i=0,3
           nv1 = Tdomain%specel(nelem)%Near_Vertex(i)
           nv2 = Tdomain%specel(nelem)%Near_Vertex(mod(i+1,4))
           call IsCrossProductOk(Tdomain, positions(n,0), positions(n,1), tol, nv1, nv2, is_inside, is_onface)
           if (.NOT. is_inside) exit
        enddo
        if(is_inside) then
           ! Capteur de type "Element"
           tmp_capt_elem(ncaptElem) = nelem
           Tdomain%type_capteurs(n) = 1
           write(*,*) "Capteur ", n, " is in element ",nelem
           ncaptElem = ncaptElem+1
           exit
        endif
        if(is_onface) then
           ! Capteur de type "Face"
           nface = Tdomain%specel(nelem)%Near_Face(i)
           if (ANY(tmp_capt_elem .EQ. nface)) then
              is_onface = .NOT. is_onface
           endif
           if (is_onface) then
              tmp_capt_elem(ncaptFace)= nface
              Tdomain%type_capteurs(n) = 2
              write(*,*) "Capteur ", n, " is on face ",nface
              ncaptFace = ncaptFace+1
           endif
           exit
        endif
     enddo
  enddo

  allocate (Tdomain%elems_capteurs(0:nCaptElem-1,0:2))
  allocate (Tdomain%faces_capteurs(0:nCaptFace-1,0:1))
  Tdomain%elems_capteurs(:,0) = tmp_capt_elem(0:nCaptElem-1)
  Tdomain%faces_capteurs(:,0) = tmp_capt_face(0:nCaptFace-1)

  !Position des points de gauss pour chacun des capteurs dans l'element ou la face
  ! Ici nous prendrons les points milieux
  do n=0,nCaptElem-1
     nelem = Tdomain%elems_capteurs(n,0)
     ngx = floor(0.5*Tdomain%specel(nelem)%ngllx)
     ngz = floor(0.5*Tdomain%specel(nelem)%ngllz)
     Tdomain%elems_capteurs(n,1) = ngx
     Tdomain%elems_capteurs(n,2) = ngz
  enddo
  do n=0,nCaptFace-1
     nface = Tdomain%faces_capteurs(n,0)
     ngx = floor(0.5*Tdomain%sface(nface)%ngll)
     Tdomain%faces_capteurs(n,1) = ngx
  enddo

  deallocate(tmp_capt_elem)
  deallocate(tmp_capt_face)

  ! Ouverture des fichiers de Capteurs
  if (Tdomain%nCapt .GT. 15) then
     STOP "To much capteurs (16 maximum)"
  endif

  do n=0,Tdomain%nCapt-1
     write(n2string,'(i10)')  n
     fnamef = "./traces/Capteur"//trim(adjustl(n2string))
     open(50+n,file=fnamef,status="replace", form="formatted")
  enddo

end subroutine init_capteurs_veloc


!>
!!\file capteurs_veloc.F90
!!\brief 
!!\version 1.0
!!\date 15/01/2014
!! This subroutine computes the cross-product between two vectors
!! being formedby the captor coordinates ant two following vertices
!! of an element.
!<
subroutine IsCrossProductOk(Tdomain, Xcapt, Zcapt, tol, nv1, nv2, is_inside, is_onface)

  use sDomain
  implicit none

  type(Domain), Intent(IN) :: Tdomain
  real, Intent(IN)          :: Xcapt, Zcapt, tol
  integer, Intent(IN)       :: nv1, nv2
  logical, Intent(INOUT)    :: is_inside, is_onface

  real, dimension(0:1)      :: vect1, vect2
  real                      :: crossprod
  integer                   :: n1, n2

  n1 = Tdomain%sVertex(nv1)%Glob_Numbering
  n2 = Tdomain%sVertex(nv2)%Glob_Numbering
  
  Vect1(0) = Tdomain%coord_nodes(0,n1) - Xcapt
  Vect1(1) = Tdomain%coord_nodes(1,n1) - Zcapt
  Vect2(0) = Tdomain%coord_nodes(0,n2) - Xcapt
  Vect2(1) = Tdomain%coord_nodes(1,n2) - Zcapt

  crossprod = Vect1(0)*Vect2(1) - Vect1(1)*Vect2(0)

  if(crossprod .LT. -tol) then
     is_inside = .false.
     is_onface = .false.
  else if(crossprod .GT. tol) then
     is_inside = .true.
     is_onface = .false.
  elseif(abs(crossprod) .LT. tol) then
     is_inside = .false.
     is_onface = .true.
  endif

end subroutine IsCrossProductOk

