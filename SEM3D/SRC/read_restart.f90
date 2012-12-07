subroutine read_restart (Tdomain,rg)

use sdomain

implicit none

type (domain), intent (INOUT):: Tdomain
integer, intent (IN) :: rg

! local variables
integer :: n,ngllx,nglly,ngllz,i,j,k,ngll,ngll1,ngll2
character (len=20) :: fnamef

write (fnamef,"(a,I2.2)") "backup",rg
open (61,file=fnamef,status="unknown",form="formatted")
read(61,*) Tdomain%TimeD%rtime, Tdomain%TimeD%NtimeMin

! Save Fields for Elements
do n = 0,Tdomain%n_elem-1
    ngllx = Tdomain%specel(n)%ngllx;  nglly = Tdomain%specel(n)%nglly; ngllz = Tdomain%specel(n)%ngllz
    if ( .not. Tdomain%specel(n)%PML ) then   
      do k = 1,ngllz-2
        do j = 1,nglly-2
          do i = 1,ngllx-2
               read(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
               read(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
               read(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
               read(61,*) Tdomain%specel(n)%Displ(i,j,k,0)
               read(61,*) Tdomain%specel(n)%Displ(i,j,k,1)
               read(61,*) Tdomain%specel(n)%Displ(i,j,k,2)
           enddo
        enddo
      enddo
    else
      do k = 1,ngllz-2
       do j = 1,nglly-2
        do i = 1,ngllx-2
           read(61,*) Tdomain%specel(n)%Veloc(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Veloc(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Veloc(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Veloc1(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Veloc1(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Veloc1(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Veloc2(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Veloc2(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Veloc2(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Veloc3(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Veloc3(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Veloc3(i,j,k,2)
        enddo
      enddo 
     enddo 
      do k = 0,ngllz-1
       do j = 0,nglly-1
        do i = 0,ngllx-1
           read(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress1(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress2(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Diagonal_Stress3(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Residual_Stress1(i,j,k,2)
           read(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,0)
           read(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,1)
           read(61,*) Tdomain%specel(n)%Residual_Stress2(i,j,k,2)
        enddo
      enddo 
     enddo 
   endif 
enddo

! Save Fields for Faces
do n = 0,Tdomain%n_face-1
   ngll1 = Tdomain%sFace(n)%ngll1; ngll2 = Tdomain%sFace(n)%ngll2
   if (.not. Tdomain%sFace(n)%PML ) then
    do j = 1,ngll2-2
     do i = 1,ngll1-2
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
       read(61,*) Tdomain%sFace(n)%Displ(i,j,0)
       read(61,*) Tdomain%sFace(n)%Displ(i,j,1)
       read(61,*) Tdomain%sFace(n)%Displ(i,j,2)
      enddo
     enddo
   else
     do j = 1,ngll2-2
      do i = 1,ngll1-2
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc(i,j,2)
       read(61,*) Tdomain%sFace(n)%Veloc1(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc1(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc1(i,j,2)
       read(61,*) Tdomain%sFace(n)%Veloc2(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc2(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc2(i,j,2)
       read(61,*) Tdomain%sFace(n)%Veloc3(i,j,0)
       read(61,*) Tdomain%sFace(n)%Veloc3(i,j,1)
       read(61,*) Tdomain%sFace(n)%Veloc3(i,j,2)
      enddo
     enddo
   endif
enddo

! Save Fields for Edges
do n = 0,Tdomain%n_edge-1
   ngll = Tdomain%sEdge(n)%ngll
   if (.not. Tdomain%sEdge(n)%PML ) then
     do i = 1,ngll-2
       read(61,*) Tdomain%sEdge(n)%Veloc(i,0)
       read(61,*) Tdomain%sEdge(n)%Veloc(i,1)
       read(61,*) Tdomain%sEdge(n)%Veloc(i,2)
       read(61,*) Tdomain%sEdge(n)%Displ(i,0)
       read(61,*) Tdomain%sEdge(n)%Displ(i,1)
       read(61,*) Tdomain%sEdge(n)%Displ(i,2)
     enddo
   else
     do i = 1,ngll-2
       read(61,*) Tdomain%sEdge(n)%Veloc(i,0)
       read(61,*) Tdomain%sEdge(n)%Veloc(i,1)
       read(61,*) Tdomain%sEdge(n)%Veloc(i,2)
       read(61,*) Tdomain%sEdge(n)%Veloc1(i,0)
       read(61,*) Tdomain%sEdge(n)%Veloc1(i,1)
       read(61,*) Tdomain%sEdge(n)%Veloc1(i,2)
       read(61,*) Tdomain%sEdge(n)%Veloc2(i,0)
       read(61,*) Tdomain%sEdge(n)%Veloc2(i,1)
       read(61,*) Tdomain%sEdge(n)%Veloc2(i,2)
       read(61,*) Tdomain%sEdge(n)%Veloc3(i,0)
       read(61,*) Tdomain%sEdge(n)%Veloc3(i,1)
       read(61,*) Tdomain%sEdge(n)%Veloc3(i,2)
     enddo
   endif
enddo

! Save Fields for Vertices
do n = 0,Tdomain%n_vertex-1
   if (.not. Tdomain%sVertex(n)%PML ) then
     do i = 0,2
       read(61,*) Tdomain%sVertex(n)%Veloc(i)
       read(61,*) Tdomain%sVertex(n)%Displ(i)
     enddo
   else
     do i = 0,2
       read(61,*) Tdomain%sVertex(n)%Veloc(i)
       read(61,*) Tdomain%sVertex(n)%Veloc1(i)
       read(61,*) Tdomain%sVertex(n)%Veloc2(i)
       read(61,*) Tdomain%sVertex(n)%Veloc3(i)
     enddo
   endif
enddo
close(61) 



return
end subroutine read_restart
