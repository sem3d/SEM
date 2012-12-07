subroutine Newmark (Tdomain,rg,ntime)

use sdomain

implicit none

include 'mpif.h'

type (domain), intent (INOUT) :: Tdomain
integer, intent (IN) :: rg,ntime

logical, dimension(:), allocatable :: L_Face, L_Edge, L_Vertex
integer ::  n, mat, nelem, ngll1,ngll2,ngll3, ngllx,nglly,ngllz, nf,ne,nv, nf_aus,ne_aus,nv_aus,ipoint,ndt2
integer ::  ngll, code, i,j,k, x,y,z, ngllPML
integer :: shift, I_give_to, I_take_from, n_rings, ntimetrace
integer, parameter :: etiquette=100
integer, dimension(mpi_status_size) :: statut
integer :: flag, nummax
real :: alpha, bega, gam1, Dt, max, maximum
real, dimension (0:2) :: V_free_vertex
real, dimension(:,:), allocatable :: VE_free
real, dimension(:,:,:), allocatable :: VF_free
logical :: sortie

! Predictor-MultiCorrector Newmark Velocity Scheme within a 
! Time staggered Stress-Velocity formulation inside PML
! PML needs to be implemented


if (Tdomain%TimeD%velocity_scheme) then


! Predictor Phase   

alpha = Tdomain%TimeD%alpha
bega = Tdomain%TimeD%beta / Tdomain%TimeD%gamma
gam1 = 1. / Tdomain%TimeD%gamma


do n = 0,Tdomain%n_elem-1 
    mat = Tdomain%specel(n)%mat_index 
    Dt = Tdomain%sSubDomain(mat)%dt

    if (.not. Tdomain%specel(n)%PML) then
        call Prediction_Elem_Veloc (Tdomain%specel(n), alpha, bega, gam1, Dt)
    else
        call get_PMLprediction_v2el (Tdomain, n, bega, dt) 
        call get_PMLprediction_e2el (Tdomain, n, bega, dt)
        call get_PMLprediction_f2el (Tdomain, n, bega, dt)

	   if (Tdomain%specel(n)%FPML) then
              call Prediction_Elem_FPML_Veloc (Tdomain%specel(n),bega, dt, Tdomain%sSubDomain(mat)%hTPrimex, &
                                        Tdomain%sSubDomain(mat)%hPrimey, Tdomain%sSubDomain(mat)%hprimez, rg, n, & 
                                              Tdomain%sSubDomain(mat)%freq)
          else
                call Prediction_Elem_PML_Veloc (Tdomain%specel(n),bega, dt, Tdomain%sSubDomain(mat)%hTPrimex, &
                                        Tdomain%sSubDomain(mat)%hPrimey, Tdomain%sSubDomain(mat)%hprimez, rg, n)
          endif
    endif
enddo

allocate (L_Face(0:Tdomain%n_face-1))
L_Face = .true.
allocate (L_Edge(0:Tdomain%n_edge-1))
L_Edge = .true.
allocate (L_Vertex(0:Tdomain%n_vertex-1))
L_Vertex = .true.
do n = 0,Tdomain%n_elem-1
    mat = Tdomain%specel(n)%mat_index
    Dt = Tdomain%sSubDomain(mat)%dt
    do i = 0,5
        nf = Tdomain%specel(n)%Near_Faces(i)
        if (L_Face(nf)) then
            L_Face(nf) = .false.
            if (.not.Tdomain%sface(nf)%PML)  call Prediction_Face_Veloc (Tdomain%sface(nf), alpha, bega, gam1, dt)
        endif
    enddo
    do i = 0,11
        ne = Tdomain%specel(n)%Near_Edges(i)
        if (L_Edge(ne)) then
            L_Edge(ne) = .false.
            if (.not.Tdomain%sedge(ne)%PML)  call Prediction_Edge_Veloc (Tdomain%sedge(ne), alpha, bega, gam1, dt)
        endif
    enddo
    do i = 0,7
        nv = Tdomain%specel(n)%Near_Vertices(i)
        if (L_Vertex(nv)) then
            L_Vertex(nv) = .false.
            if (.not.Tdomain%svertex(nv)%PML)  call Prediction_Vertex_Veloc (Tdomain%svertex(nv), alpha, bega, gam1, dt)
        endif
    enddo
enddo
deallocate (L_Face,L_Edge,L_Vertex)


! Solution phase

do n = 0,Tdomain%n_elem-1
  mat = Tdomain%specel(n)%mat_index
  if (.not. Tdomain%specel(n)%PML ) then
     call get_Displ_Face2Elem (Tdomain,n)
     call get_Displ_Edge2Elem (Tdomain,n)
     call get_Displ_Vertex2Elem (Tdomain,n)

     call compute_InternalForces_Elem (Tdomain%specel(n), Tdomain%sSubDomain(mat)%hprimex, &
                                       Tdomain%sSubDomain(mat)%hTprimex, Tdomain%sSubDomain(mat)%hprimey, &
                                       Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hprimez, &
                                       Tdomain%sSubDomain(mat)%hTprimez)
  else
     call compute_InternalForces_PML_Elem (Tdomain%specel(n), Tdomain%sSubDomain(mat)%hprimex, &
                                           Tdomain%sSubDomain(mat)%hTprimey, Tdomain%sSubDomain(mat)%hTprimez)
  endif
enddo


! External Forces
if (Tdomain%logicD%any_source) then
do n = 0, Tdomain%n_source-1
   if (rg == Tdomain%sSource(n)%proc) then
      i = Tdomain%Ssource(n)%elem
      x = Tdomain%Ssource(n)%gll(0)
      y = Tdomain%Ssource(n)%gll(1)
      z = Tdomain%Ssource(n)%gll(2)
      nf = Tdomain%Ssource(n)%i_dir
      Tdomain%specel(i)%Forces(x,y,z,nf) = Tdomain%specel(i)%Forces(x,y,z,nf) + &
      CompSource(Tdomain%sSource(n),Tdomain%TimeD%rtime,nf)
   endif
enddo
endif


! Communication of Forces
do nf = 0,Tdomain%n_face-1
    Tdomain%sFace(nf)%Forces = 0
    if (Tdomain%sFace(nf)%PML) then
       Tdomain%sFace(nf)%Forces1 = 0
       Tdomain%sFace(nf)%Forces2 = 0
       Tdomain%sFace(nf)%Forces3 = 0
    endif
enddo
do ne = 0,Tdomain%n_edge-1
    Tdomain%sEdge(ne)%Forces = 0
    if (Tdomain%sEdge(ne)%PML) then
        Tdomain%sEdge(ne)%Forces1 = 0
        Tdomain%sEdge(ne)%Forces2 = 0
        Tdomain%sEdge(ne)%Forces3 = 0
    endif
enddo
do nv = 0,Tdomain%n_vertex-1
    Tdomain%sVertex(nv)%Forces = 0
    if (Tdomain%sVertex(nv)%PML) then
        Tdomain%sVertex(nv)%Forces1 = 0
        Tdomain%sVertex(nv)%Forces2 = 0
        Tdomain%sVertex(nv)%Forces3 = 0
    endif
enddo


do n = 0,Tdomain%n_elem-1
  call get_Forces_Elem2Face(Tdomain,n)
  call get_Forces_Elem2Edge(Tdomain,n)
  call get_Forces_Elem2Vertex(Tdomain,n)
enddo


! MPI communications
if ( Tdomain%n_proc > 1 ) then
do n = 0,Tdomain%n_proc-1
    ngll = 0
    ngllPML = 0
    do i = 0,Tdomain%sComm(n)%nb_faces-1
        nf = Tdomain%sComm(n)%faces(i)
        do j = 1,Tdomain%sFace(nf)%ngll2-2
            do k = 1,Tdomain%sFace(nf)%ngll1-2
                Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sFace(nf)%Forces(k,j,0:2)
                ngll = ngll + 1
            enddo
        enddo
        if (Tdomain%sFace(nf)%PML) then
            do j = 1,Tdomain%sFace(nf)%ngll2-2
                do k = 1,Tdomain%sFace(nf)%ngll1-2
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sFace(nf)%Forces1(k,j,0:2)
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sFace(nf)%Forces2(k,j,0:2)
                    Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sFace(nf)%Forces3(k,j,0:2)
                    ngllPML = ngllPML + 1
                enddo
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_edges-1
        ne = Tdomain%sComm(n)%edges(i)
        do j = 1,Tdomain%sEdge(ne)%ngll-2
            Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sEdge(ne)%Forces(j,0:2)
            ngll = ngll + 1
        enddo
        if (Tdomain%sEdge(ne)%PML) then
            do j = 1,Tdomain%sEdge(ne)%ngll-2
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sEdge(ne)%Forces1(j,0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sEdge(ne)%Forces2(j,0:2)
                Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sEdge(ne)%Forces3(j,0:2)
                ngllPML = ngllPML + 1
            enddo
        endif
    enddo
    do i = 0,Tdomain%sComm(n)%nb_vertices-1
        nv =  Tdomain%sComm(n)%vertices(i)
        Tdomain%sComm(n)%GiveForces(ngll,0:2) = Tdomain%sVertex(nv)%Forces(0:2)
        ngll = ngll + 1
        if (Tdomain%sVertex(nv)%PML) then
            Tdomain%sComm(n)%GiveForcesPML(ngllPML,1,0:2) = Tdomain%sVertex(nv)%Forces1(0:2)
            Tdomain%sComm(n)%GiveForcesPML(ngllPML,2,0:2) = Tdomain%sVertex(nv)%Forces2(0:2)
            Tdomain%sComm(n)%GiveForcesPML(ngllPML,3,0:2) = Tdomain%sVertex(nv)%Forces3(0:2)
            ngllPML = ngllPML + 1
        endif
    enddo
enddo

n = Tdomain%n_proc
do shift = 1,n-1
    I_give_to = rg + shift
    if (I_give_to > n-1)   I_give_to = I_give_to - n
    I_take_from = rg - shift
    if (I_take_from < 0)   I_take_from = I_take_from + n
    if (mod(n,shift)==0 .and. shift/=1) then
        n_rings = shift
    else if (mod(n,n-shift)==0 .and. shift/=n-1) then
        n_rings = n-shift
    else if (mod(n,2)==0 .and. mod(shift,2)==0) then
        n_rings = 2
    else
        n_rings = 1
    endif
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngll>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForces, 3*Tdomain%sComm(I_give_to)%ngll, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngll>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForces, 3*Tdomain%sComm(I_take_from)%ngll, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngll>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForces, 3*Tdomain%sComm(I_take_from)%ngll, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngll>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForces, 3*Tdomain%sComm(I_give_to)%ngll, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
    call MPI_BARRIER (MPI_COMM_WORLD, code)
    do i = 0,n_rings-1
        if (rg==i) then
            if (Tdomain%sComm(I_give_to)%ngllPML>0) then
             call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForcesPML, 9*Tdomain%sComm(I_give_to)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
            endif
            if (Tdomain%sComm(I_take_from)%ngllPML>0) then
             call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForcesPML, 9*Tdomain%sComm(I_take_from)%ngllPML, &
                            MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
            endif
        else
            do j = 0,n/n_rings-1
                if (rg == i + j*n_rings) then
                    if (Tdomain%sComm(I_take_from)%ngllPML>0) then
                     call MPI_RECV (Tdomain%sComm(I_take_from)%TakeForcesPML, 9*Tdomain%sComm(I_take_from)%ngllPML, &
                                    MPI_DOUBLE_PRECISION, I_take_from, etiquette, MPI_COMM_WORLD, statut, code)
                    endif
                    if (Tdomain%sComm(I_give_to)%ngllPML>0) then
                     call MPI_SEND (Tdomain%sComm(I_give_to)%GiveForcesPML, 9*Tdomain%sComm(I_give_to)%ngllPML, &
                                    MPI_DOUBLE_PRECISION, I_give_to, etiquette, MPI_COMM_WORLD, code)
                    endif
                endif
            enddo
        endif
    enddo
    call MPI_BARRIER (MPI_COMM_WORLD, code)
enddo

do n = 0,Tdomain%n_proc-1
    ngll = 0
    ngllPML = 0
    call Comm_Forces_Face(Tdomain,n,ngll,ngllPML)
    call Comm_Forces_Edge(Tdomain,n,ngll,ngllPML)
    call Comm_Forces_Vertex(Tdomain,n,ngll,ngllPML)
enddo
endif


! Forces from the Neumann condition
if (Tdomain%logicD%neumann_local_present) then

   do nf = 0, Tdomain%sNeu%n_faces-1
      ngll1 = Tdomain%sNeu%nFace(nf)%ngll1
      ngll2 = Tdomain%sNeu%nFace(nf)%ngll2
      mat = Tdomain%sNeu%nFace(nf)%mat_index
      call compute_nforces_on_face (Tdomain%sNeu%nFace(nf),Tdomain%sNeu%nParam, & 
                                            Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime,nf)
   enddo   

   do ne = 0, Tdomain%sNeu%n_edges-1
      ngll = Tdomain%sNeu%nEdge(ne)%ngll
      mat = Tdomain%sNeu%nEdge(ne)%mat_index
      call compute_nforces_on_edge (Tdomain%sNeu%nEdge(ne),Tdomain%sNeu%nParam, & 
                                            Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime,ne)
   enddo 

   do nv = 0, Tdomain%sNeu%n_vertices-1   
      call compute_nforces_on_vertex (Tdomain%sNeu%nVertex(nv),Tdomain%sNeu%nParam, &
                                              Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime,nv,rg)
   enddo

   do nf = 0, Tdomain%sNeu%n_faces-1
      ngll1 = Tdomain%sNeu%nFace(nf)%ngll1
      ngll2 = Tdomain%sNeu%nFace(nf)%ngll2
      nf_aus = Tdomain%sNeu%nFace(nf)%Face
      if (.not.Tdomain%sface(nf_aus)%PML) then
        do i = 0,2
           Tdomain%sFace(nf_aus)%Forces(1:ngll1-2,1:ngll2-2,i) = Tdomain%sFace(nf_aus)%Forces(1:ngll1-2,1:ngll2-2,i) - & 
                                                                 Tdomain%sNeu%nFace(nf)%Forces(1:ngll1-2,1:ngll2-2,i)
        enddo
      else
        do i = 0,2
           Tdomain%sFace(nf_aus)%Forces3(1:ngll1-2,1:ngll2-2,i) = Tdomain%sFace(nf_aus)%Forces3(1:ngll1-2,1:ngll2-2,i) - & 
                                                                  Tdomain%sNeu%nFace(nf)%Forces(1:ngll1-2,1:ngll2-2,i)
        enddo
     endif
   enddo

   do ne = 0, Tdomain%sNeu%n_edges-1
      ngll = Tdomain%sNeu%nEdge(ne)%ngll
      ne_aus = Tdomain%sNeu%nEdge(ne)%Edge
      if (.not.Tdomain%sedge(ne_aus)%PML) then
        do i = 0,2
           Tdomain%sEdge(ne_aus)%Forces(1:ngll-2,i) = Tdomain%sEdge(ne_aus)%Forces(1:ngll-2,i) - & 
                                                      Tdomain%sNeu%nEdge(ne)%Forces(1:ngll-2,i)
        enddo
      else
        do i = 0,2
           Tdomain%sEdge(ne_aus)%Forces3(1:ngll-2,i) = Tdomain%sEdge(ne_aus)%Forces3(1:ngll-2,i) - & 
                                                       Tdomain%sNeu%nEdge(ne)%Forces(1:ngll-2,i)
        enddo   
      endif
   enddo

   do nv = 0, Tdomain%sNeu%n_vertices-1
      nv_aus = Tdomain%sNeu%nVertex(nv)%Vertex
      if (.not.Tdomain%svertex(nv_aus)%PML) then
        Tdomain%sVertex(nv_aus)%Forces(0:2) = Tdomain%sVertex(nv_aus)%Forces(0:2) - Tdomain%sNeu%nVertex(nv)%Forces(0:2)
      else
        Tdomain%sVertex(nv_aus)%Forces3(0:2) = Tdomain%sVertex(nv_aus)%Forces3(0:2) - Tdomain%sNeu%nVertex(nv)%Forces(0:2)
      endif
   enddo
endif


! Forces from the super object
if (Tdomain%logicD%super_object_local_present) then
  if (Tdomain%super_object_type == "P") then  ! Plane Wave

   do nf = 0, Tdomain%sPlaneW%n_faces-1
      ngll1 = Tdomain%sPlaneW%pFace(nf)%ngll1
      ngll2 = Tdomain%sPlaneW%pFace(nf)%ngll2
      mat = Tdomain%sPlaneW%pFace(nf)%mat_index
      allocate (VF_free(1:ngll1-2,1:ngll2-2,0:2))
      VF_free = 0
      nf_aus = Tdomain%sPlaneW%pFace(nf)%Face_UP
      call get_vel_face(Tdomain%sFace(nf_aus),VF_free,ngll1,ngll2,Tdomain%sSubdomain(mat)%dt,.false.,0,nf)
      nf_aus = Tdomain%sPlaneW%pFace(nf)%Face_DOWN
      call get_vel_face(Tdomain%sFace(nf_aus),VF_free,ngll1,ngll2,Tdomain%sSubdomain(mat)%dt,.true.,Tdomain%sPlaneW%pFace(nf)%Orient,nf)
      call compute_pforces_on_face (Tdomain%sPlaneW%pFace(nf),Tdomain%sPlaneW%pParam,VF_free, & 
                                            Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime,nf,rg)
      deallocate (VF_free)
   enddo   

   do ne = 0, Tdomain%sPlaneW%n_edges-1
      ngll = Tdomain%sPlaneW%pEdge(ne)%ngll
      mat = Tdomain%sPlaneW%pEdge(ne)%mat_index
      allocate (VE_free(1:ngll-2,0:2))
      VE_free = 0
      ne_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_UP
      call get_vel_edge(Tdomain%sEdge(ne_aus),VE_free,ngll,Tdomain%sSubdomain(mat)%dt,.false.,0)
      ne_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_DOWN
      call get_vel_edge(Tdomain%sEdge(ne_aus),VE_free,ngll,Tdomain%sSubdomain(mat)%dt,.true.,Tdomain%sPlaneW%pEdge(ne)%Orient)
      call compute_pforces_on_edge (Tdomain%sPlaneW%pEdge(ne),Tdomain%sPlaneW%pParam,VE_free, & 
                                            Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime,ne,rg)
      deallocate (VE_free)
   enddo   

   do nv = 0, Tdomain%sPlaneW%n_vertices-1
      nv_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP
      V_free_Vertex = 0
      call get_vel_vertex (Tdomain%sVertex(nv_aus), V_free_Vertex, Tdomain%sSubdomain(mat)%dt, .false.)
      nv_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_DOWN
      call get_vel_vertex (Tdomain%sVertex(nv_aus), V_free_Vertex, Tdomain%sSubdomain(mat)%dt, .true.)	      
      call compute_pforces_on_vertex (Tdomain%sPlaneW%pVertex(nv),Tdomain%sPlaneW%pParam,V_free_vertex, &
                                              Tdomain%sSubdomain(mat)%dt,Tdomain%TimeD%rtime,nv,rg)
   enddo


   do nf = 0, Tdomain%sPlaneW%n_faces-1
      ngll1 = Tdomain%sPlaneW%pFace(nf)%ngll1
      ngll2 = Tdomain%sPlaneW%pFace(nf)%ngll2
      nf_aus = Tdomain%sPlaneW%pFace(nf)%Face_UP
       do i = 0,2
         Tdomain%sFace(nf_aus)%Forces(1:ngll1-2,1:ngll2-2,i) = Tdomain%sFace(nf_aus)%Forces(1:ngll1-2,1:ngll2-2,i) - & 
                                                               Tdomain%sPlaneW%pFace(nf)%Forces_Up(1:ngll1-2,1:ngll2-2,i)
       enddo
      nf_aus = Tdomain%sPlaneW%pFace(nf)%Face_DOWN
       do i = 0,2
         Tdomain%sFace(nf_aus)%Forces(1:ngll1-2,1:ngll2-2,i) = Tdomain%sFace(nf_aus)%Forces(1:ngll1-2,1:ngll2-2,i) + & 
                                                               Tdomain%sPlaneW%pFace(nf)%Forces_Down(1:ngll1-2,1:ngll2-2,i)
       enddo
   enddo

   do ne = 0, Tdomain%sPlaneW%n_edges-1
      ngll = Tdomain%sPlaneW%pEdge(ne)%ngll
      ne_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_UP
       do i = 0,2
         Tdomain%sEdge(ne_aus)%Forces(1:ngll-2,i) = Tdomain%sEdge(ne_aus)%Forces(1:ngll-2,i) - & 
                                                    Tdomain%sPlaneW%pEdge(ne)%Forces_Up(1:ngll-2,i)
       enddo
      ne_aus = Tdomain%sPlaneW%pEdge(ne)%Edge_DOWN
       do i = 0,2
         Tdomain%sEdge(ne_aus)%Forces(1:ngll-2,i) = Tdomain%sEdge(ne_aus)%Forces(1:ngll-2,i) + & 
                                                    Tdomain%sPlaneW%pEdge(ne)%Forces_Down(1:ngll-2,i)
       enddo
   enddo

   do nv = 0, Tdomain%sPlaneW%n_vertices-1
      nv_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_UP
       Tdomain%sVertex(nv_aus)%Forces(0:2) = Tdomain%sVertex(nv_aus)%Forces(0:2) - Tdomain%sPlaneW%pVertex(nv)%Forces_Up(0:2)
      nv_aus = Tdomain%sPlaneW%pVertex(nv)%Vertex_DOWN
       Tdomain%sVertex(nv_aus)%Forces(0:2) = Tdomain%sVertex(nv_aus)%Forces(0:2) + Tdomain%sPlaneW%pVertex(nv)%Forces_Down(0:2)
   enddo

  endif
endif


! Correction phase
allocate (L_Face(0:Tdomain%n_face-1))
L_Face = .true.
allocate (L_Edge(0:Tdomain%n_edge-1))
L_Edge = .true.
allocate (L_Vertex(0:Tdomain%n_vertex-1))
L_Vertex = .true.
do n = 0,Tdomain%n_elem-1
    mat = Tdomain%specel(n)%mat_index
    if (.not. Tdomain%specel(n)%PML) then
        call Correction_Elem_Veloc (Tdomain%specel(n),bega, gam1,Tdomain%sSubDomain(mat)%Dt)
    else 
        if (Tdomain%specel(n)%FPML) then
             call Correction_Elem_FPML_Veloc (Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt, Tdomain%sSubdomain(mat)%freq)
        else
            call Correction_Elem_PML_Veloc (Tdomain%specel(n),Tdomain%sSubDomain(mat)%Dt)
        endif
    endif
    do i = 0,5
        nf = Tdomain%specel(n)%Near_Faces(i)
        if (L_Face(nf)) then
            L_Face(nf) = .false.
            if (.not.Tdomain%sface(nf)%PML) then
                call Correction_Face_Veloc (Tdomain%sface(nf), bega, gam1, dt)
            else
              if (Tdomain%sface(nf)%FPML) then
                   call Correction_Face_FPML_Veloc (Tdomain%sface(nf), dt, Tdomain%sSubdomain(mat)%freq)
              else
                   call Correction_Face_PML_Veloc (Tdomain%sface(nf), dt)
              endif
            endif
        endif
    enddo
    do i = 0,11
        ne = Tdomain%specel(n)%Near_Edges(i)
        if (L_Edge(ne)) then
            L_Edge(ne) = .false.
            if (.not.Tdomain%sedge(ne)%PML) then
                call Correction_Edge_Veloc (Tdomain%sedge(ne), bega, gam1, dt,ne)
            else
              if (Tdomain%sedge(ne)%FPML) then
                   call Correction_Edge_FPML_Veloc (Tdomain%sedge(ne), dt, Tdomain%sSubdomain(mat)%freq)
             else
                   call Correction_Edge_PML_Veloc (Tdomain%sedge(ne), dt)
             endif
            endif
        endif
    enddo
    do i = 0,7
        nv = Tdomain%specel(n)%Near_Vertices(i)
        if (L_Vertex(nv)) then
            L_Vertex(nv) = .false.
            if (.not. Tdomain%svertex(nv)%PML) then
                call Correction_Vertex_Veloc (Tdomain%svertex(nv), bega, gam1, dt)
            else
              if (Tdomain%svertex(nv)%FPML) then
                  call Correction_Vertex_FPML_Veloc (Tdomain%svertex(nv), dt, Tdomain%sSubdomain(mat)%freq)
              else
                  call Correction_Vertex_PML_Veloc (Tdomain%svertex(nv), dt)
             endif
            endif
        endif
    enddo
enddo
deallocate (L_Face,L_Edge,L_Vertex)


! Save Trace

if (Tdomain%logicD%save_trace) then

   do n = 0, Tdomain%n_receivers-1

      ndt2 = Tdomain%sReceiver(n)%ndt

      if (rg == Tdomain%sReceiver(n)%proc) then

         if ( mod (ntime,Tdomain%TimeD%ntrace) == 0 ) then
          if ( Tdomain%sReceiver(n)%flag == 1 ) allocate (Tdomain%sReceiver(n)%StoreTrace (0:Tdomain%TimeD%ntrace-1, 0:2) )
          if ( Tdomain%sReceiver(n)%flag == 2 ) allocate (Tdomain%sReceiver(n)%StoreTrace (0:(Tdomain%TimeD%ntrace-1)/ndt2, 0:2) )
          Tdomain%sReceiver(n)%StoreTrace = 0.
        endif

        i = Tdomain%sReceiver(n)%elem
        ngll1 = Tdomain%specel(i)%ngllx
        ngll2 = Tdomain%specel(i)%nglly
        ngll3 = Tdomain%specel(i)%ngllz

        if ( Tdomain%sReceiver(n)%flag == 1 .or. ( (Tdomain%sReceiver(n)%flag == 2) .and. (mod (ntime+1,ndt2)==0) ) ) then    

         if ( Tdomain%sReceiver(n)%flag == 1 ) ntimetrace = mod (ntime,Tdomain%TimeD%ntrace)
         if ( Tdomain%sReceiver(n)%flag == 2 ) ntimetrace = mod ((ntime+1)/ndt2-1,Tdomain%TimeD%ntrace/ndt2)

         do x = 0,ngll1-1
          do y = 0,ngll2-1
           do z = 0,ngll3-1
              if (x==0) then
                 if (y==0) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(0)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(4)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(6)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(z,:)
                    endif
                 else if (y==ngll2-1) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(3)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(7)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(10)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(z,:)
                    endif
                 else if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(3)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(y,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(11)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(y,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(4)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(y,z,:)
                 endif
              else if (x==ngll1-1) then
                 if (y==0) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(1)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(5)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(4)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(z,:)
                    endif
                 else if (y==ngll2-1) then
                    if (z==0) then
                       nv = Tdomain%specel(i)%Near_Vertices(2)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else if (z==ngll3-1) then
                       nv = Tdomain%specel(i)%Near_Vertices(6)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%svertex(nv)%Veloc(:)
                    else
                       ne = Tdomain%specel(i)%Near_Edges(7)
                       Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(z,:)
                    endif
                 else if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(1)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(y,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(8)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(y,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(2)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(y,z,:)
                 endif
              else if (y==0) then
                 if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(0)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(x,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(5)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(x,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(1)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(x,z,:)
                 endif
              else if (y==ngll2-1) then
                 if (z==0) then
                    ne = Tdomain%specel(i)%Near_Edges(2)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(x,:)
                 else if (z==ngll3-1) then
                    ne = Tdomain%specel(i)%Near_Edges(9)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sedge(ne)%Veloc(x,:)
                 else
                    nf = Tdomain%specel(i)%Near_Faces(3)
                    Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(x,z,:)
                 endif
              else if (z==0) then
                 nf = Tdomain%specel(i)%Near_Faces(0)
                 Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(x,y,:)
              else if (z==ngll3-1) then
                 nf = Tdomain%specel(i)%Near_Faces(5)
                 Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%sface(nf)%Veloc(x,y,:)
              else
                 Tdomain%sReceiver(n)%coeff(x,y,z,:) = Tdomain%specel(i)%Veloc(x,y,z,:)
              endif
           enddo
          enddo
         enddo

         do x = 0,ngll1-1
          do y = 0,ngll2-1
           do z = 0,ngll3-1
              Tdomain%sReceiver(n)%StoreTrace(ntimetrace,:) = Tdomain%sReceiver(n)%StoreTrace(ntimetrace,:) + &
                                Tdomain%sReceiver(n)%coeff(x,y,z,:) * Tdomain%sReceiver(n)%pol(x,y,z)
           enddo
          enddo
         enddo

       endif

      endif
   enddo
endif

else

 write(*,*) "Sorry, Acceleration Scheme not implemented yet"

endif   ! if Scheme


if (rg==0) print *,'Num iteration - Time',ntime,Tdomain%TimeD%rtime

return

end subroutine Newmark
