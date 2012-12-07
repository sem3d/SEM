subroutine global_numbering (Tdomain)

use sdomain

implicit none

type(domain), target, intent (INOUT) :: Tdomain

integer :: n, icount, i,j,k, ngllx,nglly,ngllz, nf,ne,nv, ngll1,ngll2


icount = 0

do n = 0,Tdomain%n_elem-1
   ngllx = Tdomain%specel(n)%ngllx
   nglly = Tdomain%specel(n)%nglly
   ngllz = Tdomain%specel(n)%ngllz
   allocate (Tdomain%specel(n)%Iglobnum(0:ngllx-1, 0:nglly-1, 0:ngllz-1))
   Tdomain%specel(n)%Iglobnum = -1
   do k = 1,ngllz-2
       do j = 1,nglly-2
           do i = 1,ngllx-2
               Tdomain%specel(n)%Iglobnum(i,j,k) = icount
               icount = icount + 1
           enddo
       enddo
   enddo
enddo

do n = 0,Tdomain%n_face-1
   ngllx = Tdomain%sFace(n)%ngll1
   nglly = Tdomain%sFace(n)%ngll2
   allocate (Tdomain%sFace(n)%Iglobnum_Face(1:ngllx-2,1:nglly-2))
   Tdomain%sFace(n)%Iglobnum_Face = -1
   do j = 1,nglly-2
       do i = 1,ngllx-2
           Tdomain%sFace(n)%Iglobnum_Face(i,j) = icount
           icount = icount + 1
       enddo
   enddo
enddo

do n = 0,Tdomain%n_edge-1
    ngllx = Tdomain%sEdge(n)%ngll
    allocate (Tdomain%sEdge(n)%Iglobnum_Edge(1:ngllx-2))
    Tdomain%sEdge(n)%Iglobnum_Edge = -1
    do i = 1,ngllx-2
        Tdomain%sEdge(n)%Iglobnum_Edge(i) = icount
        icount = icount + 1
    enddo
enddo

do n = 0,Tdomain%n_vertex-1
    Tdomain%sVertex(n)%Iglobnum_Vertex = icount
    icount = icount + 1
enddo

Tdomain%n_glob_points = icount 
do n = 0,Tdomain%n_elem - 1
    ngllx = Tdomain%specel(n)%ngllx
    nglly = Tdomain%specel(n)%nglly
    ngllz = Tdomain%specel(n)%ngllz

    ! Taking information from faces
nf = Tdomain%specel(n)%near_faces(0)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(0) == 0 ) then
  Tdomain%specel(n)%Iglobnum (1:ngll1-2, 1:ngll2-2, 0) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
else if ( Tdomain%specel(n)%Orient_Faces(0) == 1 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll1-1-i, j, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 2 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (i, ngll2-1-j, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 3 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll1-1-i, ngll2-1-j, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 4 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (j, i, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 5 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll2-1-j, i, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 6 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (j, ngll1-1-i, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(0) == 7 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll2-1-j, ngll1-1-i, 0) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
endif

nf = Tdomain%specel(n)%near_faces(1)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(1) == 0 ) then
    Tdomain%specel(n)%Iglobnum (1:ngll1-2, 0, 1:ngll2-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
else if ( Tdomain%specel(n)%Orient_Faces(1) == 1 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
    Tdomain%specel(n)%Iglobnum (ngll1-1-i, 0, j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 2 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
    Tdomain%specel(n)%Iglobnum (i, 0, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 3 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
    Tdomain%specel(n)%Iglobnum (ngll1-1-i, 0, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 4 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
    Tdomain%specel(n)%Iglobnum (j, 0, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 5 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
    Tdomain%specel(n)%Iglobnum (ngll2-1-j, 0, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 6 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
    Tdomain%specel(n)%Iglobnum (j, 0, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(1) == 7 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
    Tdomain%specel(n)%Iglobnum (ngll2-1-j, 0, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
endif

nf = Tdomain%specel(n)%near_faces(2)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(2) == 0 ) then
  Tdomain%specel(n)%Iglobnum (ngllx-1, 1:ngll1-2, 1:ngll2-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2) 
else if ( Tdomain%specel(n)%Orient_Faces(2) == 1 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngllx-1, ngll1-1-i, j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 2 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngllx-1, i, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 3 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngllx-1, ngll1-1-i, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 4 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngllx-1, j, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 5 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngllx-1, ngll2-1-j, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 6 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngllx-1, j, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
    enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(2) == 7 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngllx-1, ngll2-1-j, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
endif

nf = Tdomain%specel(n)%near_faces(3)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(3) == 0 ) then
  Tdomain%specel(n)%Iglobnum (1:ngll1-2, nglly-1, 1:ngll2-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
else if ( Tdomain%specel(n)%Orient_Faces(3) == 1 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll1-1-i, nglly-1, j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 2 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (i, nglly-1, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 3 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll1-1-i, nglly-1, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 4 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (j, nglly-1, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 5 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll2-1-j,nglly-1,i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 6 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (j, nglly-1, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(3) == 7 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll2-1-j, nglly-1, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
endif

nf = Tdomain%specel(n)%near_faces(4)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(4) == 0 ) then
  Tdomain%specel(n)%Iglobnum (0, 1:ngll1-2, 1:ngll2-2) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
else if ( Tdomain%specel(n)%Orient_Faces(4) == 1 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (0, ngll1-1-i, j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 2 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (0, i, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 3 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (0, ngll1-1-i, ngll2-1-j) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 4 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (0, j, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 5 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (0, ngll2-1-j, i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 6 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (0, j, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(4) == 7 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (0, ngll2-1-j, ngll1-1-i) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
endif

nf = Tdomain%specel(n)%near_faces(5)
ngll1 = Tdomain%sFace(nf)%ngll1
ngll2 = Tdomain%sFace(nf)%ngll2
if ( Tdomain%specel(n)%Orient_Faces(5) == 0 ) then
  Tdomain%specel(n)%Iglobnum (1:ngll1-2, 1:ngll2-2, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(1:ngll1-2, 1:ngll2-2)
else if ( Tdomain%specel(n)%Orient_Faces(5) == 1 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll1-1-i, j, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 2 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (i, ngll2-1-j, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 3 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll1-1-i, ngll2-1-j, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 4 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (j, i, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 5 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll2-1-j, i, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j)
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 6 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (j, ngll1-1-i, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
else if ( Tdomain%specel(n)%Orient_Faces(5) == 7 ) then
 do i=1,ngll1-2
   do j=1,ngll2-2
     Tdomain%specel(n)%Iglobnum (ngll2-1-j, ngll1-1-i, ngllz-1) = Tdomain%sFace(nf)%Iglobnum_Face(i, j) 
   enddo
 enddo
endif


    ! Taking information from edges
    ne = Tdomain%specel(n)%Near_Edges(0)
    if ( Tdomain%specel(n)%Orient_Edges(0) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(1:ngllx-2,0,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    else
      do i=1,ngllx-2 
        Tdomain%specel(n)%Iglobnum(i,0,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllx-1-i)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(1)
    if ( Tdomain%specel(n)%Orient_Edges(1) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(ngllx-1,1:nglly-2,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    else
      do j=1,nglly-2
        Tdomain%specel(n)%Iglobnum(ngllx-1,j,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(nglly-1-j)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(2)
    if ( Tdomain%specel(n)%Orient_Edges(2) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(1:ngllx-2,nglly-1,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    else
      do i=1,ngllx-2
        Tdomain%specel(n)%Iglobnum(i,nglly-1,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllx-1-i)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(3)
    if ( Tdomain%specel(n)%Orient_Edges(3) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(0,1:nglly-2,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    else
      do j=1,nglly-2
        Tdomain%specel(n)%Iglobnum(0,j,0) = Tdomain%sEdge(ne)%Iglobnum_Edge(nglly-1-j)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(4)
    if ( Tdomain%specel(n)%Orient_Edges(4) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(ngllx-1,0,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    else
      do k=1,ngllz-2
        Tdomain%specel(n)%Iglobnum(ngllx-1,0,k) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllz-1-k)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(5)
    if ( Tdomain%specel(n)%Orient_Edges(5) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(1:ngllx-2,0,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    else
      do i=1,ngllx-2
        Tdomain%specel(n)%Iglobnum(i,0,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllx-1-i)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(6)
    if ( Tdomain%specel(n)%Orient_Edges(6) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(0,0,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    else
      do k=1,ngllz-2
        Tdomain%specel(n)%Iglobnum(0,0,k) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllz-1-k)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(7)
    if ( Tdomain%specel(n)%Orient_Edges(7) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    else
      do k=1,ngllz-2
        Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,k) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllz-1-k)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(8)
    if ( Tdomain%specel(n)%Orient_Edges(8) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(ngllx-1,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    else
      do j=1,nglly-2
        Tdomain%specel(n)%Iglobnum(ngllx-1,j,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(nglly-1-j)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(9)
    if ( Tdomain%specel(n)%Orient_Edges(9) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(1:ngllx-2,nglly-1,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllx-2)
    else
      do i=1,ngllx-2
        Tdomain%specel(n)%Iglobnum(i,nglly-1,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllx-1-i)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(10)
    if ( Tdomain%specel(n)%Orient_Edges(10) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(0,nglly-1,1:ngllz-2) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:ngllz-2)
    else
      do k=1,ngllz-2
        Tdomain%specel(n)%Iglobnum(0,nglly-1,k) = Tdomain%sEdge(ne)%Iglobnum_Edge(ngllz-1-k)
      enddo
    endif
    ne = Tdomain%specel(n)%Near_Edges(11)
    if ( Tdomain%specel(n)%Orient_Edges(11) == 0 ) then 
      Tdomain%specel(n)%Iglobnum(0,1:nglly-2,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(1:nglly-2)
    else
      do j=1,nglly-2
        Tdomain%specel(n)%Iglobnum(0,j,ngllz-1) = Tdomain%sEdge(ne)%Iglobnum_Edge(nglly-1-j)
      enddo
    endif


    ! Taking information from vertices
    nv = Tdomain%specel(n)%Near_Vertices(0)
    Tdomain%specel(n)%Iglobnum(0,0,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(1)
    Tdomain%specel(n)%Iglobnum(ngllx-1,0,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(2)
    Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(3)
    Tdomain%specel(n)%Iglobnum(0,nglly-1,0) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(4)
    Tdomain%specel(n)%Iglobnum(0,0,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(5)
    Tdomain%specel(n)%Iglobnum(ngllx-1,0,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(6)
    Tdomain%specel(n)%Iglobnum(ngllx-1,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
    nv = Tdomain%specel(n)%Near_Vertices(7)
    Tdomain%specel(n)%Iglobnum(0,nglly-1,ngllz-1) = Tdomain%sVertex(nv)%Iglobnum_Vertex
enddo


return
end subroutine global_numbering
