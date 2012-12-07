subroutine define_neu_properties (Tdomain)

use sdomain
implicit none

type(Domain), intent (INOUT) :: Tdomain

integer :: ngll,ngll1,ngll2,mat_index,nf,ne,nv,nf_aus


open (24,file=Tdomain%neumann_file,status="old", form="formatted")
read (24,*) 
read (24,*) mat_index
read (24,*) Tdomain%sNeu%nParam%wtype
read (24,*) Tdomain%sNeu%nParam%lx
read (24,*) Tdomain%sNeu%nParam%ly
read (24,*) Tdomain%sNeu%nParam%lz
read (24,*) Tdomain%sNeu%nParam%xs
read (24,*) Tdomain%sNeu%nParam%ys
read (24,*) Tdomain%sNeu%nParam%zs
read (24,*) Tdomain%sNeu%nParam%f0
close(24)

do nf = 0, Tdomain%sNeu%n_faces-1
  ngll1 = Tdomain%sNeu%nFace(nf)%ngll1
  ngll2 = Tdomain%sNeu%nFace(nf)%ngll2
  nf_aus = Tdomain%sNeu%nFace(nf)%Face
  Tdomain%sNeu%nFace(nf)%mat_index = Tdomain%specel(Tdomain%sFace(nf_aus)%Which_Elem)%mat_index
  allocate (Tdomain%sNeu%nFace(nf)%Btn(0:ngll1-1,0:ngll2-1,0:2))
  allocate (Tdomain%sNeu%nFace(nf)%Forces(1:ngll1-2,1:ngll2-2,0:2))
enddo

do ne = 0, Tdomain%sNeu%n_edges-1
  ngll = Tdomain%sNeu%nEdge(ne)%ngll
  Tdomain%sNeu%nEdge(ne)%mat_index =  mat_index       
  allocate (Tdomain%sNeu%nEdge(ne)%Btn(1:ngll-2,0:2))
  allocate (Tdomain%sNeu%nEdge(ne)%Forces(1:ngll-2,0:2))
  Tdomain%sNeu%nEdge(ne)%Btn = 0
enddo

do nv = 0, Tdomain%sNeu%n_vertices-1
  Tdomain%sNeu%nVertex(nv)%mat_index =  mat_index    
  Tdomain%sNeu%nVertex(nv)%Btn = 0
enddo


Tdomain%sNeu%nParam%Lambda = Tdomain%sSubDomain(mat_index)%DLambda
Tdomain%sNeu%nParam%Mu = Tdomain%sSubDomain(mat_index)%DMu

if (Tdomain%sNeu%nParam%wtype == 'S') then
   Tdomain%sNeu%nParam%speed = Tdomain%sSubDomain(mat_index)%Sspeed
else
   Tdomain%sNeu%nParam%speed = Tdomain%sSubDomain(mat_index)%Pspeed
endif


return
end subroutine define_neu_properties
